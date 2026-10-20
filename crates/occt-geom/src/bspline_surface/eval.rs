//! Evaluation of `Geom_BSplineSurface` (`EvalD0..D3`, `EvalDN`,
//! `HomogeneousD1`) plus the free-standing [`eval_bspline_surface`].

use super::knots::distinct_km;
use super::*;
use occt_core::gp::{GpPnt, GpVec};

impl GeomBSplineSurface {
    /// `Geom_BSplineSurface::EvalDN` (`Geom_BSplineSurface_1.cxx:279-312`).
    /// `UIndex`/`VIndex` are 0 so `LocateParameter` always runs. Knots are
    /// already flat (`NoMults`).
    pub fn eval_dn_bspl(&self, u: f64, v: f64, nu: i32, nv: i32) -> GpVec {
        if nu + nv < 1 || nu < 0 || nv < 0 {
            // cxx:281-284 throws `Geom_UndefinedDerivative`.
            return GpVec::new(0.0, 0.0, 0.0);
        }
        self.dn_orders(u, v, nu, nv)
    }

    /// Distinct knots / multiplicities of one direction, restricted to the
    /// surface's own parameter window for a periodic direction (the flat
    /// periodic sequence carries one period of extension knots outside it).
    fn distinct_km_pub(flat: &[f64], periodic: bool, range: (f64, f64)) -> (Vec<f64>, Vec<i32>) {
        distinct_km(flat, periodic, range)
    }

    /// Shared knot/mult args for `BSplSLib::{D2,D3,DN}` on this surface.
    /// Borrowed from the stored representation, so an evaluation allocates
    /// nothing for the knots.
    fn bspl_lib_args(
        &self,
    ) -> (
        Option<(&[f64], &[i32])>,
        Option<(&[f64], &[i32])>,
        bool,
    ) {
        let rat = self.weights.is_some();
        if self.has_periodic_direction() {
            (
                Some((self.uknots.as_slice(), self.umults.as_slice())),
                Some((self.vknots.as_slice(), self.vmults.as_slice())),
                rat,
            )
        } else {
            (None, None, rat)
        }
    }

    /// `BSplSLib::HomogeneousD1` on this surface (`BSplSLib.cxx:974-1059`).
    pub fn homogeneous_d1_bspl(
        &self,
        u: f64,
        v: f64,
    ) -> occt_core::bspl::surface_d1::HomogSurfD1 {
        let (uk, vk, rat) = self.bspl_lib_args();
        match (uk, vk) {
            (Some(uk), Some(vk)) => occt_core::bspl::surface_d1::homogeneous_d1(
                u,
                v,
                0,
                0,
                &self.poles,
                self.weights.as_deref(),
                uk.0,
                vk.0,
                Some(uk.1),
                Some(vk.1),
                self.deg_u as i32,
                self.deg_v as i32,
                rat,
                rat,
                self.u_periodic,
                self.v_periodic,
            ),
            _ => occt_core::bspl::surface_d1::homogeneous_d1(
                u,
                v,
                0,
                0,
                &self.poles,
                self.weights.as_deref(),
                &self.knots_u,
                &self.knots_v,
                None,
                None,
                self.deg_u as i32,
                self.deg_v as i32,
                rat,
                rat,
                self.u_periodic,
                self.v_periodic,
            ),
        }
    }

    /// `Geom_BSplineSurface::EvalD0` (`Geom_BSplineSurface_1.cxx:131-159`) via
    /// `BSplSLib::D0` (`surface_d0::surface_d0`).
    pub fn eval_d0_bspl(&self, u: f64, v: f64) -> GpPnt {
        let (uk, vk, rat) = self.bspl_lib_args();
        match (uk, vk) {
            (Some(uk), Some(vk)) => occt_core::bspl::surface_d0::surface_d0(
                u,
                v,
                0,
                0,
                &self.poles,
                self.weights.as_deref(),
                uk.0,
                vk.0,
                Some(uk.1),
                Some(vk.1),
                self.deg_u as i32,
                self.deg_v as i32,
                rat,
                rat,
                self.u_periodic,
                self.v_periodic,
            ),
            _ => occt_core::bspl::surface_d0::surface_d0(
                u,
                v,
                0,
                0,
                &self.poles,
                self.weights.as_deref(),
                &self.knots_u,
                &self.knots_v,
                None,
                None,
                self.deg_u as i32,
                self.deg_v as i32,
                rat,
                rat,
                self.u_periodic,
                self.v_periodic,
            ),
        }
    }

    /// `Geom_BSplineSurface::EvalD1` (`Geom_BSplineSurface_1.cxx:163-195`) via
    /// `BSplSLib::D1` (`surface_d1::surface_d1`).
    pub fn eval_d1_bspl(&self, u: f64, v: f64) -> occt_core::bspl::surface_d1::SurfD1 {
        let (uk, vk, rat) = self.bspl_lib_args();
        match (uk, vk) {
            (Some(uk), Some(vk)) => occt_core::bspl::surface_d1::surface_d1(
                u,
                v,
                0,
                0,
                &self.poles,
                self.weights.as_deref(),
                uk.0,
                vk.0,
                Some(uk.1),
                Some(vk.1),
                self.deg_u as i32,
                self.deg_v as i32,
                rat,
                rat,
                self.u_periodic,
                self.v_periodic,
            ),
            _ => occt_core::bspl::surface_d1::surface_d1(
                u,
                v,
                0,
                0,
                &self.poles,
                self.weights.as_deref(),
                &self.knots_u,
                &self.knots_v,
                None,
                None,
                self.deg_u as i32,
                self.deg_v as i32,
                rat,
                rat,
                self.u_periodic,
                self.v_periodic,
            ),
        }
    }

    /// `Geom_BSplineSurface::EvalD2` (`Geom_BSplineSurface_1.cxx:199-235`) via
    /// `BSplSLib::D2` (`surface_d2::surface_d2`).
    pub fn eval_d2_bspl(&self, u: f64, v: f64) -> occt_core::bspl::surface_d2::SurfD2 {
        let (uk, vk, rat) = self.bspl_lib_args();
        match (uk, vk) {
            (Some(uk), Some(vk)) => occt_core::bspl::surface_d2::surface_d2(
                u,
                v,
                0,
                0,
                &self.poles,
                self.weights.as_deref(),
                uk.0,
                vk.0,
                Some(uk.1),
                Some(vk.1),
                self.deg_u as i32,
                self.deg_v as i32,
                rat,
                rat,
                self.u_periodic,
                self.v_periodic,
            ),
            _ => occt_core::bspl::surface_d2::surface_d2(
                u,
                v,
                0,
                0,
                &self.poles,
                self.weights.as_deref(),
                &self.knots_u,
                &self.knots_v,
                None,
                None,
                self.deg_u as i32,
                self.deg_v as i32,
                rat,
                rat,
                self.u_periodic,
                self.v_periodic,
            ),
        }
    }

    /// `Geom_BSplineSurface::EvalD3` (`Geom_BSplineSurface_1.cxx:239-274`) via
    /// `BSplSLib::D3` (`surface_d3::surface_d3`).
    pub fn eval_d3_bspl(&self, u: f64, v: f64) -> occt_core::elib::surface_eval_d3::SurfD3 {
        let (uk, vk, rat) = self.bspl_lib_args();
        match (uk, vk) {
            (Some(uk), Some(vk)) => occt_core::bspl::surface_d3::surface_d3(
                u,
                v,
                0,
                0,
                &self.poles,
                self.weights.as_deref(),
                uk.0,
                vk.0,
                Some(uk.1),
                Some(vk.1),
                self.deg_u as i32,
                self.deg_v as i32,
                rat,
                rat,
                self.u_periodic,
                self.v_periodic,
            ),
            _ => occt_core::bspl::surface_d3::surface_d3(
                u,
                v,
                0,
                0,
                &self.poles,
                self.weights.as_deref(),
                &self.knots_u,
                &self.knots_v,
                None,
                None,
                self.deg_u as i32,
                self.deg_v as i32,
                rat,
                rat,
                self.u_periodic,
                self.v_periodic,
            ),
        }
    }

    /// Periodic-surface value / derivative through the ported `BSplSLib::DN`
    /// (`prepare_eval::surface_dn`), whose pole indexing wraps across the period
    /// (`prepare_eval.rs:174-207`). `nu`/`nv` are the derivative orders in U/V.
    fn dn_orders(&self, u: f64, v: f64, nu: i32, nv: i32) -> GpVec {
        let rat = self.weights.is_some();
        if self.has_periodic_direction() {
            // OCCT's EvalDN passes the distinct knots and multiplicities
            // (Geom_BSplineSurface_1.cxx:279-312), whose periodic arm wraps the
            // knot window; the flat NoMults arm cannot.
            return occt_core::bspl::prepare_eval::surface_dn(
                u,
                v,
                nu,
                nv,
                0,
                0,
                &self.poles,
                self.weights.as_deref(),
                &self.uknots,
                &self.vknots,
                Some(&self.umults),
                Some(&self.vmults),
                self.deg_u as i32,
                self.deg_v as i32,
                rat,
                rat,
                self.u_periodic,
                self.v_periodic,
            );
        }
        occt_core::bspl::prepare_eval::surface_dn(
            u,
            v,
            nu,
            nv,
            0,
            0,
            &self.poles,
            self.weights.as_deref(),
            &self.knots_u,
            &self.knots_v,
            None,
            None,
            self.deg_u as i32,
            self.deg_v as i32,
            rat,
            rat,
            self.u_periodic,
            self.v_periodic,
        )
    }

    pub(super) fn has_periodic_direction(&self) -> bool {
        self.u_periodic || self.v_periodic
    }
}

/// All non-zero degree-`degree` basis functions at `u` (NURBS A2.2), returned
/// as a full-length vector of length `#poles`.
///
/// The parameter is deliberately NOT clamped into the knot domain.
/// `BSplSLib::PrepareEval` (`BSplSLib.cxx:313-364`) resolves each parameter
/// through `BSplCLib::LocateParameter` (`BSplCLib.cxx:321-364`), which wraps
/// the parameter into the period only when that direction is periodic
/// (`BSplCLib.cxx:246-249`) and otherwise returns the parameter unchanged,
/// clamping only the knot-span index to `[First, Last - 1]` and skipping
/// zero-width spans (`BSplCLib.cxx:277-312`). `BSplCLib::Eval`
/// (`BSplCLib.cxx:865-1000`) then evaluates the located span's polynomial, so
/// an out-of-domain parameter is EXTRAPOLATED with the end span.
///
/// Clamping `u` to the domain here instead fabricated a deviation of up to one
/// chord of the surface's period for any pcurve whose parameter range straddles
/// the surface's parameter domain (a closed pcurve on a surface whose domain
/// spans exactly one revolution always does), which
/// `ShapeAnalysis_Edge::CheckSameParameter` (`ShapeAnalysis_Edge.cxx:787-795`
/// via `BRepLib_ValidateEdge`) reported as a real deviation and
/// `ShapeFix_Edge::FixSameParameter` wrote back as an inflated edge tolerance.
pub(super) fn basis_funs(knots: &[f64], degree: usize, u: f64, periodic: bool) -> Vec<f64> {
    let n = knots.len() - degree - 2; // last control-point index
    // `BSplCLib::LocateParameter(Degree, Knots, Mults = nullptr, U, Periodic,
    // KnotIndex, NewU)` with incoming `KnotIndex = 0`, so the `[FromK1, ToK2]`
    // form (`BSplCLib.cxx:218-317`) always runs. The index is 1-based there.
    let (span_1, u) = occt_core::bspl::locate::locate_parameter(degree as i32, knots, None, u, periodic, 0);
    let span = (span_1 - 1).clamp(degree as i32, n as i32) as usize;
    let mut nvec = vec![0.0; degree + 1];
    nvec[0] = 1.0;
    let mut left = vec![0.0; degree + 1];
    let mut right = vec![0.0; degree + 1];
    for j in 1..=degree {
        left[j] = u - knots[span + 1 - j];
        right[j] = knots[span + j] - u;
        let mut saved = 0.0;
        for r in 0..j {
            let denom = right[r + 1] + left[j - r];
            let temp = if denom.abs() > 1e-300 { nvec[r] / denom } else { 0.0 };
            nvec[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        nvec[j] = saved;
    }
    let mut out = vec![0.0; n + 1];
    for j in 0..=degree {
        let idx = span - degree + j;
        if idx <= n {
            out[idx] = nvec[j];
        }
    }
    out
}

/// Evaluate the tensor-product B-spline surface at `(u, v)`.
///
/// The rational case computes `Σ N_i(u)N_j(v) w_ij P_ij / Σ N_i(u)N_j(v) w_ij`;
/// the non-rational case is the same with unit weights.
pub fn eval_bspline_surface(s: &GeomBSplineSurface, u: f64, v: f64) -> GpPnt {
    let nu = s.poles.len();
    if nu == 0 {
        return GpPnt::zero();
    }
    let nv = s.poles[0].len();
    if nv == 0 {
        return GpPnt::zero();
    }
    let bu = basis_funs(&s.knots_u, s.deg_u, u, s.u_periodic);
    let bv = basis_funs(&s.knots_v, s.deg_v, v, s.v_periodic);

    let mut acc = [0.0_f64, 0.0, 0.0];
    match &s.weights {
        Some(w) => {
            let mut den = 0.0;
            for i in 0..nu {
                for j in 0..nv {
                    let b = bu[i] * bv[j] * w[i][j];
                    let p = s.poles[i][j];
                    acc[0] += b * p.x();
                    acc[1] += b * p.y();
                    acc[2] += b * p.z();
                    den += b;
                }
            }
            if den.abs() < 1e-300 {
                return GpPnt::zero();
            }
            GpPnt::new(acc[0] / den, acc[1] / den, acc[2] / den)
        }
        None => {
            for i in 0..nu {
                for j in 0..nv {
                    let b = bu[i] * bv[j];
                    let p = s.poles[i][j];
                    acc[0] += b * p.x();
                    acc[1] += b * p.y();
                    acc[2] += b * p.z();
                }
            }
            GpPnt::new(acc[0], acc[1], acc[2])
        }
    }
}
