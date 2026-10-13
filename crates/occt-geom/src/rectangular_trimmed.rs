//! Rectangular trimmed surface. Source: `Geom_RectangularTrimmedSurface` (TKG3d).
//!
//! A patch of a basis surface limited by two U bounds and two V bounds. Both
//! constructors and `SetTrim` normalise the window the way
//! `Geom_RectangularTrimmedSurface.cxx` does: a periodic direction folds the
//! window with `ElCLib::AdjustPeriodic`, a non-periodic direction swaps the
//! bounds and inverts the requested sense when `U1 > U2`, and a zero-length or
//! out-of-range window raises `Standard_ConstructionError` (returned here as
//! `Err`, since Rust has no exceptions).

use std::sync::Arc;

use occt_core::elib::clib::adjust_periodic;
use occt_core::gp::{GpCone, GpCylinder, GpPln, GpPnt, GpSphere, GpTorus, GpTrsf, GpVec};
use occt_core::precision::Precision;

use crate::curve::Curve;
use crate::geom_adaptor_local::{try_local_d1, try_local_d2};
use crate::offset_surface::GeomOffsetSurface;
use crate::surface::Surface;
use crate::trimmed::GeomTrimmedCurveBasis;

/// Portion of a surface limited by U and/or V parameter bounds.
/// Source: `Geom_RectangularTrimmedSurface.cxx`.
#[derive(Clone)]
pub struct GeomRectangularTrimmedSurface {
    basis: Arc<dyn Surface>,
    utrim1: f64,
    vtrim1: f64,
    utrim2: f64,
    vtrim2: f64,
    is_utrimmed: bool,
    is_vtrimmed: bool,
}

/// Adaptor domain tolerance used by `UTrim`/`VTrim` (`Adaptor3d_CurveOnSurface`
/// EvalFirstLastSurf passes `Precision::PConfusion`).
const ADAPTOR_TOL: f64 = Precision::PCONFUSION;

/// The `T.IsNull() ? S->Copy() : T->BasisSurface()->Copy()` preamble shared by
/// both `Geom_RectangularTrimmedSurface` constructors (`.cxx:91-100`,
/// `.cxx:123-132`): an incoming patch is replaced by its basis surface.
fn copy_untrimmed(s: &Arc<dyn Surface>) -> Arc<dyn Surface> {
    s.rectangular_trimmed_basis()
        .unwrap_or_else(|| Arc::from(s.clone_dyn()))
}

fn ieee_remainder(x: f64, y: f64) -> f64 {
    if !y.is_finite() || y.abs() < f64::EPSILON {
        return x;
    }
    x - y * (x / y).round()
}

fn closed_when_trimmed(periodic: bool, period: f64, a: f64, b: f64) -> bool {
    if !periodic {
        return false;
    }
    let length = b - a;
    length > Precision::PCONFUSION
        && ieee_remainder(length, period).abs() <= Precision::PCONFUSION
}

impl GeomRectangularTrimmedSurface {
    /// `Geom_RectangularTrimmedSurface(S, U1, U2, V1, V2, USense, VSense)`
    /// (`Geom_RectangularTrimmedSurface.cxx:71-111`) with OCCT's default
    /// `USense = VSense = true` (`.hxx:64-65`).
    pub fn uv(s: Arc<dyn Surface>, u1: f64, u2: f64, v1: f64, v2: f64) -> Result<Self, String> {
        Self::uv_sense(s, u1, u2, v1, v2, true, true)
    }

    /// `Geom_RectangularTrimmedSurface(S, U1, U2, V1, V2, USense, VSense)`
    /// (`.cxx:71-111`); the sense-carrying form used by
    /// `StepToGeom::MakeRectangularTrimmedSurface` (`StepToGeom.cxx:1882`).
    pub fn uv_sense(
        s: Arc<dyn Surface>,
        u1: f64,
        u2: f64,
        v1: f64,
        v2: f64,
        u_sense: bool,
        v_sense: bool,
    ) -> Result<Self, String> {
        let mut out = Self {
            basis: copy_untrimmed(&s),
            utrim1: u1,
            vtrim1: v1,
            utrim2: u2,
            vtrim2: v2,
            is_utrimmed: true,
            is_vtrimmed: true,
        };
        // An offset basis is recut by nesting this patch inside the offset
        // (`.cxx:102-109`).
        if let (Some(ob), Some(d)) = (out.basis.offset_basis_surface(), out.basis.offset_distance())
        {
            let inner = Self::uv_sense(ob, u1, u2, v1, v2, u_sense, v_sense)?;
            out.basis = Arc::new(GeomOffsetSurface::new(Arc::new(inner), d));
        }
        out.set_trim_uv(u1, u2, v1, v2, u_sense, v_sense)?;
        Ok(out)
    }

    /// The windowed `GeomAdaptor_Surface` view used by
    /// `Adaptor3d_CurveOnSurface::EvalFirstLastSurf` (`Adaptor3d_CurveOnSurface.cxx:1778-1779`):
    /// `mySurface->UTrim(U1, U2, Tol)->VTrim(V1, V2, Tol)` returns an adaptor
    /// whose domain is `[U1,U2] x [V1,V2]`, not a
    /// `Geom_RectangularTrimmedSurface`, so no `SetTrim` normalisation happens
    /// here. `GeomAdaptor_Surface::Load(S, UFirst, ULast, VFirst, VLast, TolU,
    /// TolV)` stores the window verbatim, unwrapping an incoming patch to its
    /// basis (`GeomAdaptor_Surface.cxx:411-427`, `.hxx:167-180`), and raises
    /// only for an inverted window (`GeomAdaptor_Surface.hxx:172-175`).
    pub fn window(
        s: Arc<dyn Surface>,
        u1: f64,
        u2: f64,
        v1: f64,
        v2: f64,
    ) -> Result<Self, String> {
        if u1 > u2 || v1 > v2 {
            return Err("GeomAdaptor_Surface::Load".into());
        }
        Ok(Self {
            basis: copy_untrimmed(&s),
            utrim1: u1,
            vtrim1: v1,
            utrim2: u2,
            vtrim2: v2,
            is_utrimmed: true,
            is_vtrimmed: true,
        })
    }

    /// `Geom_RectangularTrimmedSurface(S, Param1, Param2, UTrim, Sense)`
    /// (`Geom_RectangularTrimmedSurface.cxx:113-168`) with OCCT's default
    /// `Sense = true` (`.hxx:83`).
    pub fn one_param(s: Arc<dyn Surface>, param1: f64, param2: f64, u_trim: bool) -> Result<Self, String> {
        Self::one_param_sense(s, param1, param2, u_trim, true)
    }
    /// `Geom_RectangularTrimmedSurface(S, Param1, Param2, UTrim, Sense)`
    /// (`.cxx:113-168`).
    pub fn one_param_sense(
        s: Arc<dyn Surface>,
        param1: f64,
        param2: f64,
        u_trim: bool,
        sense: bool,
    ) -> Result<Self, String> {
        let incoming = s.rectangular_trimmed_spec();
        let basis = copy_untrimmed(&s);
        let (udeb, ufin) = basis.u_range();
        let (vdeb, vfin) = basis.v_range();
        let mut out = Self {
            basis,
            utrim1: if u_trim { param1 } else { udeb },
            utrim2: if u_trim { param2 } else { ufin },
            vtrim1: if u_trim { vdeb } else { param1 },
            vtrim2: if u_trim { vfin } else { param2 },
            is_utrimmed: u_trim,
            is_vtrimmed: !u_trim,
        };
        if let (Some(ob), Some(d)) = (out.basis.offset_basis_surface(), out.basis.offset_distance())
        {
            let inner = Self::one_param_sense(ob, param1, param2, u_trim, sense)?;
            out.basis = Arc::new(GeomOffsetSurface::new(Arc::new(inner), d));
        }
        // Re-trimming an incoming patch keeps the trim already carried in the
        // other direction (`.cxx:142-155`).
        if let Some((t_u1, t_u2, t_v1, t_v2, t_u_trim, t_v_trim)) = incoming {
            if u_trim && t_v_trim {
                out.set_trim_uv(param1, param2, t_v1, t_v2, sense, true)?;
                return Ok(out);
            } else if !u_trim && t_u_trim {
                out.set_trim_uv(t_u1, t_u2, param1, param2, true, sense)?;
                return Ok(out);
            }
        }
        // OCCT reaches `SetTrim(Param1, Param2, UTrim, Sense)` with the trim
        // flags of a fresh patch (both false), so `is_vtrimmed`/`is_utrimmed`
        // select the single-direction branch there (`.cxx:210-247`).
        out.set_trim_one_param(param1, param2, u_trim, sense)?;
        Ok(out)
    }

    /// `Geom_RectangularTrimmedSurface::SetTrim(U1, U2, V1, V2, USense, VSense)`
    /// (`Geom_RectangularTrimmedSurface.cxx:196-205`).
    pub fn set_trim_uv(
        &mut self,
        u1: f64,
        u2: f64,
        v1: f64,
        v2: f64,
        u_sense: bool,
        v_sense: bool,
    ) -> Result<(), String> {
        self.set_trim(u1, u2, v1, v2, true, true, u_sense, v_sense)
    }

    /// `Geom_RectangularTrimmedSurface::SetTrim(Param1, Param2, UTrim, Sense)`
    /// (`Geom_RectangularTrimmedSurface.cxx:210-247`).
    pub fn set_trim_one_param(
        &mut self,
        param1: f64,
        param2: f64,
        u_trim: bool,
        sense: bool,
    ) -> Result<(), String> {
        // Dummy arguments to call the general SetTrim.
        let (dummy_a, dummy_b, dummy_sense) = (0.0, 0.0, true);
        if u_trim {
            if self.is_vtrimmed {
                self.set_trim(
                    param1,
                    param2,
                    self.vtrim1,
                    self.vtrim2,
                    true,
                    true,
                    sense,
                    dummy_sense,
                )
            } else {
                self.set_trim(
                    param1,
                    param2,
                    dummy_a,
                    dummy_b,
                    true,
                    false,
                    sense,
                    dummy_sense,
                )
            }
        } else if self.is_utrimmed {
            self.set_trim(
                self.utrim1,
                self.utrim2,
                param1,
                param2,
                true,
                true,
                dummy_sense,
                sense,
            )
        } else {
            self.set_trim(
                dummy_a,
                dummy_b,
                param1,
                param2,
                false,
                true,
                dummy_sense,
                sense,
            )
        }
    }

    /// `Geom_RectangularTrimmedSurface::SetTrim` (8 arguments,
    /// `Geom_RectangularTrimmedSurface.cxx:249-336`).
    ///
    /// The window is normalised against the basis bounds, then the basis is
    /// reversed for whichever direction ended up with an inverted sense.
    fn set_trim(
        &mut self,
        u1: f64,
        u2: f64,
        v1: f64,
        v2: f64,
        u_trim: bool,
        v_trim: bool,
        u_sense: bool,
        v_sense: bool,
    ) -> Result<(), String> {
        let mut u_same_sense = true;
        let mut v_same_sense = true;
        let (udeb, ufin) = self.basis.u_range();
        let (vdeb, vfin) = self.basis.v_range();

        // Trimming the U-Direction.
        self.is_utrimmed = u_trim;
        if !u_trim {
            self.utrim1 = udeb;
            self.utrim2 = ufin;
        } else {
            if u1 == u2 {
                return Err("Geom_RectangularTrimmedSurface::U1==U2".into());
            }
            if self.basis.is_u_periodic() {
                u_same_sense = u_sense;
                // Set utrim1 in the range [udeb, ufin] and utrim2 in the range
                // [utrim1, utrim1 + period].
                self.utrim1 = u1;
                self.utrim2 = u2;
                let preci = ((self.utrim2 - self.utrim1).abs() / 2.0).min(Precision::PCONFUSION);
                adjust_periodic(udeb, ufin, preci, &mut self.utrim1, &mut self.utrim2);
            } else {
                if u1 < u2 {
                    u_same_sense = u_sense;
                    self.utrim1 = u1;
                    self.utrim2 = u2;
                } else {
                    u_same_sense = !u_sense;
                    self.utrim1 = u2;
                    self.utrim2 = u1;
                }
                if (udeb - self.utrim1 > Precision::PCONFUSION)
                    || (self.utrim2 - ufin > Precision::PCONFUSION)
                {
                    return Err(
                        "Geom_RectangularTrimmedSurface::U parameters out of range".into(),
                    );
                }
            }
        }

        // Trimming the V-Direction.
        self.is_vtrimmed = v_trim;
        if !v_trim {
            self.vtrim1 = vdeb;
            self.vtrim2 = vfin;
        } else {
            if v1 == v2 {
                return Err("Geom_RectangularTrimmedSurface::V1==V2".into());
            }
            if self.basis.is_v_periodic() {
                v_same_sense = v_sense;
                self.vtrim1 = v1;
                self.vtrim2 = v2;
                let preci = ((self.vtrim2 - self.vtrim1).abs() / 2.0).min(Precision::PCONFUSION);
                adjust_periodic(vdeb, vfin, preci, &mut self.vtrim1, &mut self.vtrim2);
            } else {
                if v1 < v2 {
                    v_same_sense = v_sense;
                    self.vtrim1 = v1;
                    self.vtrim2 = v2;
                } else {
                    v_same_sense = !v_sense;
                    self.vtrim1 = v2;
                    self.vtrim2 = v1;
                }
                if (vdeb - self.vtrim1 > Precision::PCONFUSION)
                    || (self.vtrim2 - vfin > Precision::PCONFUSION)
                {
                    return Err(
                        "Geom_RectangularTrimmedSurface::V parameters out of range".into(),
                    );
                }
            }
        }

        if !u_same_sense {
            self.u_reverse()?;
        }
        if !v_same_sense {
            self.v_reverse()?;
        }
        Ok(())
    }

    /// `Geom_RectangularTrimmedSurface::UReverse`
    /// (`Geom_RectangularTrimmedSurface.cxx:342-349`): both window ends are
    /// mapped through the basis's `UReversedParameter`, the basis is U-reversed
    /// and the window is re-applied with `USense = VSense = true`.
    pub fn u_reverse(&mut self) -> Result<(), String> {
        let u1 = self.basis.u_reversed_parameter(self.utrim2);
        let u2 = self.basis.u_reversed_parameter(self.utrim1);
        let Some(reversed) = self.basis.u_reversed() else {
            return Err("Geom_RectangularTrimmedSurface::UReverse".into());
        };
        self.basis = reversed;
        self.set_trim(
            u1,
            u2,
            self.vtrim1,
            self.vtrim2,
            self.is_utrimmed,
            self.is_vtrimmed,
            true,
            true,
        )
    }

    /// `Geom_RectangularTrimmedSurface::VReverse`
    /// (`Geom_RectangularTrimmedSurface.cxx:371-378`).
    pub fn v_reverse(&mut self) -> Result<(), String> {
        let v1 = self.basis.v_reversed_parameter(self.vtrim2);
        let v2 = self.basis.v_reversed_parameter(self.vtrim1);
        let Some(reversed) = self.basis.v_reversed() else {
            return Err("Geom_RectangularTrimmedSurface::VReverse".into());
        };
        self.basis = reversed;
        self.set_trim(
            self.utrim1,
            self.utrim2,
            v1,
            v2,
            self.is_utrimmed,
            self.is_vtrimmed,
            true,
            true,
        )
    }

    /// The basis surface (`Geom_RectangularTrimmedSurface::BasisSurface`).
    pub fn basis_surface(&self) -> &Arc<dyn Surface> {
        &self.basis
    }
}

impl Surface for GeomRectangularTrimmedSurface {
    fn d0(&self, u: f64, v: f64) -> GpPnt {
        self.basis.d0(u, v)
    }

    fn d1(&self, u: f64, v: f64) -> (GpPnt, GpVec, GpVec) {
        // `GeomAdaptor_Surface::EvalD1` BSpline LocalD1 when UV is on the
        // restricted adaptor domain end (`cxx:1129-1195`).
        if let Some(bs) = self.basis.osculating_bspline() {
            if let Some(r) = try_local_d1(
                &bs,
                u,
                v,
                self.utrim1,
                self.utrim2,
                self.vtrim1,
                self.vtrim2,
                ADAPTOR_TOL,
                ADAPTOR_TOL,
            ) {
                return r;
            }
        }
        self.basis.d1(u, v)
    }

    fn d2(&self, u: f64, v: f64) -> (GpPnt, GpVec, GpVec, GpVec, GpVec, GpVec) {
        // `GeomAdaptor_Surface::EvalD2` BSpline LocalD2 (`cxx:1273-1292`).
        if let Some(bs) = self.basis.osculating_bspline() {
            if let Some(r) = try_local_d2(
                &bs,
                u,
                v,
                self.utrim1,
                self.utrim2,
                self.vtrim1,
                self.vtrim2,
                ADAPTOR_TOL,
                ADAPTOR_TOL,
            ) {
                return r;
            }
        }
        self.basis.d2(u, v)
    }

    fn u_range(&self) -> (f64, f64) {
        (self.utrim1, self.utrim2)
    }

    fn v_range(&self) -> (f64, f64) {
        (self.vtrim1, self.vtrim2)
    }

    fn is_u_periodic(&self) -> bool {
        // `Geom_RectangularTrimmedSurface.cxx:509-523`.
        if !self.basis.is_u_periodic() {
            return false;
        }
        if !self.is_utrimmed {
            return true;
        }
        closed_when_trimmed(true, self.basis.u_period(), self.utrim1, self.utrim2)
    }

    fn is_v_periodic(&self) -> bool {
        // `Geom_RectangularTrimmedSurface.cxx:534-548`.
        if !self.basis.is_v_periodic() {
            return false;
        }
        if !self.is_vtrimmed {
            return true;
        }
        closed_when_trimmed(true, self.basis.v_period(), self.vtrim1, self.vtrim2)
    }

    fn is_u_closed(&self) -> bool {
        // `Geom_RectangularTrimmedSurface.cxx:559-576`.
        if !self.is_utrimmed {
            return self.basis.is_u_closed();
        }
        closed_when_trimmed(
            self.basis.is_u_periodic(),
            self.basis.u_period(),
            self.utrim1,
            self.utrim2,
        )
    }

    fn is_v_closed(&self) -> bool {
        // `Geom_RectangularTrimmedSurface.cxx:580-597`.
        if !self.is_vtrimmed {
            return self.basis.is_v_closed();
        }
        closed_when_trimmed(
            self.basis.is_v_periodic(),
            self.basis.v_period(),
            self.vtrim1,
            self.vtrim2,
        )
    }

    fn continuity(&self) -> u8 {
        self.basis.continuity()
    }

    /// `Geom_RectangularTrimmedSurface::IsCNu` (`cxx:482-486`).
    fn is_cnu(&self, n: i32) -> Result<bool, String> {
        if n < 0 {
            return Err("Geom_RectangularTrimmedSurface::IsCNu: Geom_RectangularTrimmedSurface.cxx:485".into());
        }
        self.basis.is_cnu(n)
    }

    /// `Geom_RectangularTrimmedSurface::IsCNv` (`cxx:490-496`).
    fn is_cnv(&self, n: i32) -> Result<bool, String> {
        if n < 0 {
            return Err("Geom_RectangularTrimmedSurface::IsCNv: Geom_RectangularTrimmedSurface.cxx:494".into());
        }
        self.basis.is_cnv(n)
    }

    /// `Geom_RectangularTrimmedSurface::EvalDN` (`cxx:419-428`).
    fn try_eval_dn(&self, u: f64, v: f64, nu: i32, nv: i32) -> Result<GpVec, String> {
        if nu + nv < 1 || nu < 0 || nv < 0 {
            return Err(
                "Geom_RectangularTrimmedSurface::EvalDN: Geom_RectangularTrimmedSurface.cxx:425"
                    .into(),
            );
        }
        self.basis.try_eval_dn(u, v, nu, nv)
    }

    /// `Geom_RectangularTrimmedSurface::EvalD3` (`cxx:412-414`) forwards.
    fn try_eval_d3_uv(&self, u: f64, v: f64, nu: i32, nv: i32) -> Result<GpVec, String> {
        self.basis.try_eval_d3_uv(u, v, nu, nv)
    }

    fn transform(&mut self, t: &GpTrsf) {
        let mut b = self.basis.clone_dyn();
        b.transform(t);
        self.basis = Arc::from(b);
    }

    fn clone_dyn(&self) -> Box<dyn Surface> {
        Box::new(self.clone())
    }

    fn rectangular_trimmed_basis(&self) -> Option<Arc<dyn Surface>> {
        Some(self.basis.clone())
    }

    fn rectangular_trimmed_spec(&self) -> Option<(f64, f64, f64, f64, bool, bool)> {
        Some((
            self.utrim1,
            self.utrim2,
            self.vtrim1,
            self.vtrim2,
            self.is_utrimmed,
            self.is_vtrimmed,
        ))
    }

    /// `Geom_RectangularTrimmedSurface::UReversedParameter` (`cxx:354-357`):
    /// the basis surface's.
    fn u_reversed_parameter(&self, u: f64) -> f64 {
        self.basis.u_reversed_parameter(u)
    }

    /// `Geom_RectangularTrimmedSurface::VReversedParameter` (`cxx:383-386`):
    /// the basis surface's.
    fn v_reversed_parameter(&self, v: f64) -> f64 {
        self.basis.v_reversed_parameter(v)
    }

    /// `Geom_Surface::UReversed` (`Geom_Surface.cxx:33-38`) = `Copy()` +
    /// [`Self::u_reverse`].
    fn u_reversed(&self) -> Option<Arc<dyn Surface>> {
        let mut copy = self.clone();
        copy.u_reverse().ok()?;
        Some(Arc::new(copy))
    }

    /// `Geom_Surface::VReversed` (`Geom_Surface.cxx:43-48`) = `Copy()` +
    /// [`Self::v_reverse`].
    fn v_reversed(&self) -> Option<Arc<dyn Surface>> {
        let mut copy = self.clone();
        copy.v_reverse().ok()?;
        Some(Arc::new(copy))
    }

    // `GeomAdaptor_Surface::load` unwraps a `Geom_RectangularTrimmedSurface` to
    // its basis surface and keeps only the parameter range
    // (`GeomAdaptor_Surface.cxx:423-425`), so the `GetType()`-equivalent queries
    // must report the basis type.
    fn gp_pln(&self) -> Option<GpPln> {
        self.basis.gp_pln()
    }
    fn gp_sphere(&self) -> Option<GpSphere> {
        self.basis.gp_sphere()
    }
    fn gp_cylinder(&self) -> Option<GpCylinder> {
        self.basis.gp_cylinder()
    }
    fn gp_cone(&self) -> Option<GpCone> {
        self.basis.gp_cone()
    }
    fn gp_torus(&self) -> Option<GpTorus> {
        self.basis.gp_torus()
    }

    fn u_iso_curve(&self, u: f64) -> Option<Arc<dyn Curve>> {
        // `Geom_RectangularTrimmedSurface.cxx:444-458`. `Geom_TrimmedCurve`
        // keeps the basis parameter (`Geom_TrimmedCurve.cxx:212-215`).
        let c = self.basis.u_iso_curve(u)?;
        if self.is_vtrimmed {
            Some(Arc::new(GeomTrimmedCurveBasis::new(c, self.vtrim1, self.vtrim2)))
        } else {
            Some(c)
        }
    }

    fn v_iso_curve(&self, v: f64) -> Option<Arc<dyn Curve>> {
        // `Geom_RectangularTrimmedSurface.cxx:463-477`.
        let c = self.basis.v_iso_curve(v)?;
        if self.is_utrimmed {
            Some(Arc::new(GeomTrimmedCurveBasis::new(c, self.utrim1, self.utrim2)))
        } else {
            Some(c)
        }
    }

    fn u_period(&self) -> f64 {
        self.basis.u_period()
    }

    fn v_period(&self) -> f64 {
        self.basis.v_period()
    }

    fn bspline_surface_poles(&self) -> Option<&[Vec<GpPnt>]> {
        self.basis.bspline_surface_poles()
    }

    fn bspline_surface_uknots(&self) -> Option<&[f64]> { self.basis.bspline_surface_uknots() }

    fn bspline_surface_vknots(&self) -> Option<&[f64]> { self.basis.bspline_surface_vknots() }

    fn bspline_surface_weights(&self) -> Option<&[Vec<f64>]> {
        self.basis.bspline_surface_weights()
    }

    fn nb_u_poles(&self) -> i32 {
        self.basis.nb_u_poles()
    }
    fn nb_v_poles(&self) -> i32 {
        self.basis.nb_v_poles()
    }
    fn nb_u_intervals(&self, continuity: u8) -> i32 {
        self.basis.nb_u_intervals(continuity)
    }
    fn nb_v_intervals(&self, continuity: u8) -> i32 {
        self.basis.nb_v_intervals(continuity)
    }
    fn u_intervals(&self, continuity: u8) -> Vec<f64> {
        self.basis.u_intervals(continuity)
    }
    fn v_intervals(&self, continuity: u8) -> Vec<f64> {
        self.basis.v_intervals(continuity)
    }
    fn u_degree(&self) -> i32 {
        self.basis.u_degree()
    }
    fn v_degree(&self) -> i32 {
        self.basis.v_degree()
    }
    fn revolution_basis_curve(&self) -> Option<Arc<dyn Curve>> {
        self.basis.revolution_basis_curve()
    }
    fn extrusion_basis_curve(&self) -> Option<Arc<dyn Curve>> {
        self.basis.extrusion_basis_curve()
    }

    fn uv_resolution(&self, r3d: f64) -> Option<(f64, f64)> {
        self.basis.uv_resolution(r3d)
    }
}
