//! B-spline surface (tensor-product, clamped knots). Source: `Geom_BSplineSurface.hxx`
//!
//! Poles are stored 0-based as a `u × v` grid; weights (if any) match the
//! pole grid element-wise. Evaluation uses Cox-de Boor basis functions in each
//! parametric direction with a rational division when weights are present.

use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};

use occt_core::gp::GpPnt;

mod eval;
mod fit;
mod iso;
mod knots;
mod local;
mod periodic;
mod surface;
#[cfg(test)]
mod tests;

pub use self::eval::eval_bspline_surface;
pub use self::fit::fit_surface_grid;

/// Cache of `Geom_BSplineSurface::Resolution`'s derivative bound, i.e. OCCT's
/// `myMaxDerivInvOk` / `myUMaxDerivInv` / `myVMaxDerivInv`
/// (`Geom_BSplineSurface.hxx:1318-1323`, `Geom_BSplineSurface_1.cxx:2197-2221`).
///
/// The bound is a pure function of the geometry, so `Resolution` computes it
/// once and reuses it until a mutator clears the flag. The port caches the
/// bound itself rather than its reciprocal, so [`GeomBSplineSurface::uv_resolution`]
/// keeps dividing by it, exactly as it did before the latch existed.
///
/// Mutators of the owning surface must call [`Self::invalidate`]; the state is
/// copied by `Clone` (a clone has the same geometry, hence the same bound).
#[derive(Debug)]
pub struct ResolutionLatch {
    ok: AtomicBool,
    max_deriv_u: AtomicU64,
    max_deriv_v: AtomicU64,
}

impl ResolutionLatch {
    pub fn new() -> Self {
        Self {
            ok: AtomicBool::new(false),
            max_deriv_u: AtomicU64::new(0.0f64.to_bits()),
            max_deriv_v: AtomicU64::new(0.0f64.to_bits()),
        }
    }

    /// `myMaxDerivInvOk = false` (`Geom_BSplineSurface_1.cxx:2197-2221` and the
    /// `updateUKnots` / `updateVKnots` / pole-and-weight setters that call it):
    /// the cached bound no longer describes the surface.
    pub fn invalidate(&mut self) {
        self.ok.store(false, Ordering::Relaxed);
    }

    /// `Geom_BSplineSurface::Resolution`: run the `BSplSLib::Resolution` sweep
    /// only while the flag is clear, then reuse the bound.
    fn or_compute(&self, compute: impl FnOnce() -> (f64, f64)) -> (f64, f64) {
        if self.ok.load(Ordering::Acquire) {
            return (
                f64::from_bits(self.max_deriv_u.load(Ordering::Relaxed)),
                f64::from_bits(self.max_deriv_v.load(Ordering::Relaxed)),
            );
        }
        let (max_deriv_u, max_deriv_v) = compute();
        self.max_deriv_u.store(max_deriv_u.to_bits(), Ordering::Relaxed);
        self.max_deriv_v.store(max_deriv_v.to_bits(), Ordering::Relaxed);
        self.ok.store(true, Ordering::Release);
        (max_deriv_u, max_deriv_v)
    }
}

impl Default for ResolutionLatch {
    fn default() -> Self {
        Self::new()
    }
}

impl Clone for ResolutionLatch {
    fn clone(&self) -> Self {
        let out = Self::new();
        if self.ok.load(Ordering::Acquire) {
            out.max_deriv_u
                .store(self.max_deriv_u.load(Ordering::Relaxed), Ordering::Relaxed);
            out.max_deriv_v
                .store(self.max_deriv_v.load(Ordering::Relaxed), Ordering::Relaxed);
            out.ok.store(true, Ordering::Release);
        }
        out
    }
}

/// Non-rational or rational tensor-product B-spline surface.
///
/// The head of the representation follows OCCT: `knots_u`/`knots_v` are the
/// flat sequences (`myFlatKnots`), while the distinct knots and multiplicities
/// (`myUKnots`/`myUMults`, `myVKnots`/`myVMults`) are stored beside them and
/// rebuilt by `update_u_knots` / `update_v_knots` exactly where OCCT calls
/// `updateUKnots()` / `updateVKnots()` (`Geom_BSplineSurface.cxx:1148-1192`).
/// The head fields are `pub(crate)` so the geometry crate's own rebuilders can
/// read them, but no other crate can mutate them: every mutation from outside
/// goes through a setter that keeps the stored arrays in step, and read access
/// is through [`Self::knots_u`] and friends.
#[derive(Clone)]
pub struct GeomBSplineSurface {
    pub poles: Vec<Vec<GpPnt>>,
    pub(crate) knots_u: Vec<f64>,
    pub(crate) knots_v: Vec<f64>,
    pub(crate) deg_u: usize,
    pub(crate) deg_v: usize,
    pub weights: Option<Vec<Vec<f64>>>,
    /// `Geom_BSplineSurface` UPeriodic (`StepToGeom::MakeBSplineSurface` cxx:1101-1110).
    pub(crate) u_periodic: bool,
    /// `Geom_BSplineSurface` VPeriodic (`StepToGeom::MakeBSplineSurface` cxx:1120-1129).
    pub(crate) v_periodic: bool,
    /// `myUKnots` — distinct U knots inside the surface's own parameter window.
    pub(crate) uknots: Vec<f64>,
    /// `myUMults`.
    pub(crate) umults: Vec<i32>,
    /// `myVKnots`.
    pub(crate) vknots: Vec<f64>,
    /// `myVMults`.
    pub(crate) vmults: Vec<i32>,
    /// `Geom_BSplineSurface::myMaxDerivInvOk` / `myUMaxDerivInv` / `myVMaxDerivInv`.
    pub res_latch: ResolutionLatch,
}

/// `FirstParameter()..LastParameter()` of one direction from the flat vector,
/// with the same fallback [`GeomBSplineSurface::u_range`] uses for a malformed
/// or empty representation.
fn knot_window(flat: &[f64], degree: usize) -> (f64, f64) {
    if flat.is_empty() || degree >= flat.len() {
        return (0.0, 1.0);
    }
    (flat[degree], flat[flat.len() - 1 - degree])
}

impl GeomBSplineSurface {
    /// Build a non-rational surface. Knot counts must be `#poles + degree + 1`.
    pub fn new(
        poles: Vec<Vec<GpPnt>>,
        knots_u: Vec<f64>,
        knots_v: Vec<f64>,
        deg_u: usize,
        deg_v: usize,
    ) -> Result<Self, String> {
        if deg_u < 1 || deg_v < 1 {
            return Err("GeomBSplineSurface::new: degrees must be >= 1".to_string());
        }
        let nu = poles.len();
        if nu == 0 {
            return Err("GeomBSplineSurface::new: empty pole grid".to_string());
        }
        let nv = poles[0].len();
        if nv == 0 {
            return Err("GeomBSplineSurface::new: zero-width pole grid".to_string());
        }
        for row in &poles {
            if row.len() != nv {
                return Err("GeomBSplineSurface::new: ragged pole grid".to_string());
            }
        }
        if knots_u.len() != nu + deg_u + 1 {
            return Err(format!(
                "GeomBSplineSurface::new: knots_u count {} != poles_u {nu} + deg {deg_u} + 1",
                knots_u.len()
            ));
        }
        if knots_v.len() != nv + deg_v + 1 {
            return Err(format!(
                "GeomBSplineSurface::new: knots_v count {} != poles_v {nv} + deg {deg_v} + 1",
                knots_v.len()
            ));
        }
        if !knots_u.windows(2).all(|w| w[0] <= w[1]) {
            return Err("GeomBSplineSurface::new: knots_u not non-decreasing".to_string());
        }
        if !knots_v.windows(2).all(|w| w[0] <= w[1]) {
            return Err("GeomBSplineSurface::new: knots_v not non-decreasing".to_string());
        }
        let (u_periodic, v_periodic) = Self::periodic_flags_from_flat(
            nu, nv, deg_u, deg_v, &knots_u, &knots_v,
        );
        let mut s = Self {
            poles,
            knots_u,
            knots_v,
            deg_u,
            deg_v,
            weights: None,
            u_periodic,
            v_periodic,
            uknots: Vec::new(),
            umults: Vec::new(),
            vknots: Vec::new(),
            vmults: Vec::new(),
            res_latch: ResolutionLatch::new(),
        };
        s.update_u_knots();
        s.update_v_knots();
        Ok(s)
    }

    // ---- representation head (OCCT `myFlatKnots` / `myUKnots` / `myUMults`) ----

    /// `myFlatKnots` in U.
    pub fn knots_u(&self) -> &[f64] {
        &self.knots_u
    }

    /// `myFlatKnots` in V.
    pub fn knots_v(&self) -> &[f64] {
        &self.knots_v
    }

    /// `myDeg`.
    pub fn deg_u(&self) -> usize {
        self.deg_u
    }

    /// `myDeg` in V.
    pub fn deg_v(&self) -> usize {
        self.deg_v
    }

    /// `IsUPeriodic()`.
    pub fn u_periodic(&self) -> bool {
        self.u_periodic
    }

    /// `IsVPeriodic()`.
    pub fn v_periodic(&self) -> bool {
        self.v_periodic
    }

    /// `myUKnots` — the distinct U knots as stored.
    pub fn uknots(&self) -> &[f64] {
        &self.uknots
    }

    /// `myUMults`.
    pub fn umults(&self) -> &[i32] {
        &self.umults
    }

    /// `myVKnots`.
    pub fn vknots(&self) -> &[f64] {
        &self.vknots
    }

    /// `myVMults`.
    pub fn vmults(&self) -> &[i32] {
        &self.vmults
    }

    /// Replace the flat U knot vector (`updateUKnots()` follows, as it does
    /// after every OCCT mutator that assigns `myFlatKnots`).
    pub fn set_uknots(&mut self, flat: Vec<f64>) {
        self.knots_u = flat;
        self.update_u_knots();
    }

    /// Replace the flat V knot vector.
    pub fn set_vknots(&mut self, flat: Vec<f64>) {
        self.knots_v = flat;
        self.update_v_knots();
    }

    /// Restore `IsUPeriodic()` on a rebuilt copy (OCCT carries the flag on the
    /// object rather than in the flat vector). The stored knot window follows.
    pub fn set_u_periodic_flag(&mut self, periodic: bool) {
        if self.u_periodic == periodic {
            return;
        }
        self.u_periodic = periodic;
        self.update_u_knots();
    }

    /// Restore `IsVPeriodic()` on a rebuilt copy.
    pub fn set_v_periodic_flag(&mut self, periodic: bool) {
        if self.v_periodic == periodic {
            return;
        }
        self.v_periodic = periodic;
        self.update_v_knots();
    }

    /// `Geom_BSplineSurface::updateUKnots()` (`Geom_BSplineSurface.cxx:1148-1192`):
    /// rebuild `myUKnots`/`myUMults` from the flat vector and clear the cached
    /// resolution bound (`myMaxDerivInvOk = false`, `cxx:1150-1152`).
    pub(crate) fn update_u_knots(&mut self) {
        let (k, m) = knots::distinct_km(
            &self.knots_u,
            self.u_periodic,
            knot_window(&self.knots_u, self.deg_u),
        );
        self.uknots = k;
        self.umults = m;
        self.res_latch.invalidate();
    }

    /// `Geom_BSplineSurface::updateVKnots()`.
    pub(crate) fn update_v_knots(&mut self) {
        let (k, m) = knots::distinct_km(
            &self.knots_v,
            self.v_periodic,
            knot_window(&self.knots_v, self.deg_v),
        );
        self.vknots = k;
        self.vmults = m;
        self.res_latch.invalidate();
    }

    /// Assemble the head directly from an already-flat knot sequence and
    /// re-derive `myUKnots`/`myUMults` / `myVKnots`/`myVMults`
    /// (`updateUKnots()` / `updateVKnots()`), for the knot-insertion /
    /// segmentation / removal paths and for readers that already hold the flat
    /// `KnotSequence` forms (OCCT's `SetPoles` + assignment of the flat knots
    /// followed by `updateUKnots()`).
    pub fn from_parts(
        poles: Vec<Vec<GpPnt>>,
        weights: Option<Vec<Vec<f64>>>,
        knots_u: Vec<f64>,
        knots_v: Vec<f64>,
        deg_u: usize,
        deg_v: usize,
        u_periodic: bool,
        v_periodic: bool,
    ) -> Self {
        let mut s = Self {
            poles,
            knots_u,
            knots_v,
            deg_u,
            deg_v,
            weights,
            u_periodic,
            v_periodic,
            uknots: Vec::new(),
            umults: Vec::new(),
            vknots: Vec::new(),
            vmults: Vec::new(),
            res_latch: ResolutionLatch::new(),
        };
        s.update_u_knots();
        s.update_v_knots();
        s
    }

    /// Build a rational surface; `weights` must match the pole grid size.
    pub fn rational(
        poles: Vec<Vec<GpPnt>>,
        weights: Vec<Vec<f64>>,
        knots_u: Vec<f64>,
        knots_v: Vec<f64>,
        deg_u: usize,
        deg_v: usize,
    ) -> Result<Self, String> {
        let mut s = Self::new(poles, knots_u, knots_v, deg_u, deg_v)?;
        if s.poles.len() != weights.len() || s.poles[0].len() != weights[0].len() {
            return Err("GeomBSplineSurface::rational: weight grid mismatch".to_string());
        }
        s.weights = Some(weights);
        Ok(s)
    }

    pub fn nb_poles_u(&self) -> usize { self.poles.len() }
    pub fn nb_poles_v(&self) -> usize { self.poles.first().map_or(0, |r| r.len()) }
    pub fn is_rational(&self) -> bool { self.weights.is_some() }

    /// Distinct knots and multiplicities from a flat knot vector.
    pub fn unique_knots_mults(flat: &[f64]) -> (Vec<f64>, Vec<i32>) {
        let mut knots = Vec::new();
        let mut mults = Vec::new();
        for &k in flat {
            if knots.last() == Some(&k) {
                if let Some(m) = mults.last_mut() {
                    *m += 1;
                }
            } else {
                knots.push(k);
                mults.push(1);
            }
        }
        (knots, mults)
    }

    /// `Geom_BSplineSurface::UKnots()` / `UMultiplicities()`: the distinct U
    /// knots and their multiplicities as stored, i.e. restricted to the
    /// surface's own parameter window (`myUKnots` / `myUMults`). The extension
    /// knots of a periodic flat sequence lie outside that window and are not
    /// stored.
    ///
    /// Returns owned copies; hot paths read [`Self::uknots`] / [`Self::umults`].
    pub fn distinct_knots_and_mults_u(&self) -> (Vec<f64>, Vec<i32>) {
        (self.uknots.clone(), self.umults.clone())
    }

    /// `Geom_BSplineSurface::VKnots()` / `VMultiplicities()`; see
    /// [`Self::distinct_knots_and_mults_u`].
    pub fn distinct_knots_and_mults_v(&self) -> (Vec<f64>, Vec<i32>) {
        (self.vknots.clone(), self.vmults.clone())
    }

    /// `Geom_BSplineSurface::UKnotDistribution()`
    /// (`Geom_BSplineSurface_1.cxx:1851-1855`), i.e. `myUKnotSet` as
    /// `UpdateUKnots` sets it (`Geom_BSplineSurface.cxx:1153`:
    /// `BSplCLib::KnotAnalysis(myUDeg, myUPeriodic, myUKnots, myUMults,
    /// myUKnotSet, MaxKnotMult)`).
    pub fn u_knot_distribution(&self) -> occt_core::kernel::geomabs::BSplKnotDistribution {
        occt_core::bspl::knot_analysis::knot_analysis(
            self.deg_u as i32,
            self.u_periodic,
            &self.uknots,
            &self.umults,
        )
        .0
    }

    /// `Geom_BSplineSurface::VKnotDistribution()`
    /// (`Geom_BSplineSurface_1.cxx:1858-1862`); see
    /// [`Self::u_knot_distribution`] (`Geom_BSplineSurface.cxx:1201`).
    pub fn v_knot_distribution(&self) -> occt_core::kernel::geomabs::BSplKnotDistribution {
        occt_core::bspl::knot_analysis::knot_analysis(
            self.deg_v as i32,
            self.v_periodic,
            &self.vknots,
            &self.vmults,
        )
        .0
    }

    /// `StepToGeom::MakeBSplineSurface` periodic test (`cxx:1101-1129`).
    pub fn should_be_periodic(n_poles: usize, degree: usize, mults: &[i32]) -> bool {
        let sum_mult: i32 = mults.iter().sum();
        let n_unique = mults.len();
        if n_unique == 0 {
            return false;
        }
        if sum_mult == (n_poles + degree + 1) as i32 {
            false
        } else if mults[0] == mults[n_unique - 1] && (sum_mult - mults[0]) == n_poles as i32 {
            true
        } else {
            false
        }
    }

    fn periodic_flags_from_flat(
        nu: usize,
        nv: usize,
        deg_u: usize,
        deg_v: usize,
        knots_u: &[f64],
        knots_v: &[f64],
    ) -> (bool, bool) {
        let (_, um) = Self::unique_knots_mults(knots_u);
        let (_, vm) = Self::unique_knots_mults(knots_v);
        (
            Self::should_be_periodic(nu, deg_u, &um),
            Self::should_be_periodic(nv, deg_v, &vm),
        )
    }

    /// `Geom_BSplineSurface` from unique knots + multiplicities.
    pub fn from_poles_knots_mults(
        poles: Vec<Vec<GpPnt>>,
        u_knots: Vec<f64>,
        v_knots: Vec<f64>,
        u_mults: Vec<i32>,
        v_mults: Vec<i32>,
        deg_u: usize,
        deg_v: usize,
    ) -> Result<Self, String> {
        let flat_u = occt_core::bspl::banded_interp::knot_sequence(&u_knots, &u_mults, deg_u as i32);
        let flat_v = occt_core::bspl::banded_interp::knot_sequence(&v_knots, &v_mults, deg_v as i32);
        let mut s = Self::new(poles, flat_u, flat_v, deg_u, deg_v)?;
        // Prefer the STEP/unique mult vectors (`cxx:1101-1129`) over flat-derived.
        s.u_periodic = Self::should_be_periodic(s.poles.len(), deg_u, &u_mults);
        s.v_periodic = Self::should_be_periodic(s.poles.first().map_or(0, |r| r.len()), deg_v, &v_mults);
        s.update_u_knots();
        s.update_v_knots();
        Ok(s)
    }
}

/// Clamped uniform knot vector for `n` poles of degree `p`
/// (length `n + p + 1`, multiplicity `p + 1` at both ends).


/// Clamped uniform knot vectors for a surface with `nu × nv` poles.
pub use occt_math::spline_surface::bspline_surface_uniform_knots;
