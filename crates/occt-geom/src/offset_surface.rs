//! Offset 3D surface. Source: `Geom_OffsetSurface.hxx`
//!
//! The offset surface at `(u, v)` is the basis surface point displaced along
//! its unit normal by the signed offset: `d0 = basis.d0 + offset * n`.
//! First partials follow `Geom_OffsetSurfaceUtils::EvaluateD1`.
//! `EvalD3` and `EvalDN` follow `EvaluateD3` (`pxx:1427`) and `EvaluateDN`
//! (`pxx:1637`). `mmapptt_` orders -1, 0, 1, and 2 are ported.
//! `mmfmtb1_` transposes the IsoU tables (`AdvApp2Var_ApproxF2var.cxx:5634`).
//! `MakeApprox` stays on one 3D subspace (`AdvApp2Var_Iso.cxx:122`).
//! The rational Bezier `Segment` arm runs when weights are stored
//! (`Geom_BezierSurface.cxx:961`). `EvalD0`, `EvalD1`, `EvalD2`, and `EvalDN` use the rational B-spline
//! evaluator (`Geom_BezierSurface.cxx:1425`). `UIso` and `VIso` call
//! `BSplSLib::Iso` (`Geom_BezierSurface.cxx:1776`). A trimmed rational
//! Bezier becomes a rational B-spline (`GeomConvert_1.cxx:714`).
//! Curve `Segment` trims weights with the poles (`Geom_BezierCurve.cxx:395`).
//! `SetWeight` turns a curve or surface rational (`Geom_BezierCurve.cxx:455`,
//! `Geom_BezierSurface.cxx:1179`). `SetWeightCol` and `SetWeightRow`
//! rewrite one column or a row span (`Geom_BezierSurface.cxx:1219`,
//! `:1262`). `UReverse` and `VReverse` swap poles and weights
//! (`:1305`, `:1352`). `ExchangeUV` transposes poles and weights
//! (`Geom_BezierSurface.cxx:490`). `Increase` elevates degree through
//! `BSplSLib::IncreaseDegree` (`:532`). Insert and remove pole columns
//! and rows (`Geom_BezierSurface.cxx:652-937`). `SetPole` and column/row
//! setters (`:1047-1175`). Weighted construction (`:416`). Curve
//! `Increase`, insert/remove poles, and `Reverse`
//! (`Geom_BezierCurve.cxx:141`, `:210`, `:298`, `:352`).
use std::sync::Arc;

use occt_core::gp::{GpPnt, GpTrsf, GpVec};
use occt_core::precision::{APPROXIMATION, RESOLUTION};

use crate::adv_approx::ApproxAFunction3d;
use crate::bspline_curve::GeomBSplineCurve;
use crate::curve::Curve;
use crate::offset_surface_eval_dn::{evaluate_d3, evaluate_dn};
use crate::offset_surface_utils::{evaluate_d0, evaluate_d1};
use crate::osculating_surface::{OsculatingSurface, OSCULATING_TOL};
use crate::surface::Surface;

/// `GeomAdaptor_Surface` Offset interval shape (`cxx:664-682`).
fn offset_interval_continuity(s: u8) -> u8 {
    match s {
        0 => 2,
        2 => 4,
        4 => 5,
        _ => 6,
    }
}

/// Surface obtained by offsetting a basis surface along its normals by a
/// constant signed distance. Source: `Geom_OffsetSurface.hxx`.
#[derive(Clone)]
pub struct GeomOffsetSurface {
    basis: Arc<dyn Surface>,
    offset: f64,
    osc: OsculatingSurface,
}

impl GeomOffsetSurface {
    /// Wraps `basis` offset by the signed `distance` along the basis normal.
    pub fn new(basis: Arc<dyn Surface>, distance: f64) -> Self {
        let osc = OsculatingSurface::new(basis.clone(), OSCULATING_TOL);
        Self { basis, offset: distance, osc }
    }

    /// The basis surface.
    pub fn basis_surface(&self) -> &Arc<dyn Surface> {
        &self.basis
    }

    /// The signed offset distance.
    pub fn offset(&self) -> f64 {
        self.offset
    }

    /// Unit normal of the basis surface at `(u, v)` (the offset direction).
    /// Prefer `evaluate_d0` / `d0` for the offset point; this remains for iso.
    fn unit_normal(&self, u: f64, v: f64) -> GpVec {
        let (_, du, dv) = self.basis.d1(u, v);
        let n = du.xyz().crossed(dv.xyz());
        let m = n.modulus();
        if m > 1e-12 {
            GpVec::new(n.x / m, n.y / m, n.z / m)
        } else {
            GpVec::new(0.0, 0.0, 1.0)
        }
    }

    /// `Geom_OffsetSurface::UIso` (`Geom_OffsetSurface.cxx:601-652`).
    pub fn u_iso(&self, uu: f64) -> Option<Arc<dyn Curve>> {
        if self.basis.is_surface_of_linear_extrusion() {
            return self.u_iso_extrusion(uu);
        }
        self.approx_iso(true, uu)
    }

    /// `Geom_OffsetSurface::VIso` (`Geom_OffsetSurface.cxx:657-687`).
    pub fn v_iso(&self, vv: f64) -> Option<Arc<dyn Curve>> {
        self.approx_iso(false, vv)
    }

    /// Extrusion arm of `UIso` (`cxx:607-623`).
    fn u_iso_extrusion(&self, uu: f64) -> Option<Arc<dyn Curve>> {
        let a_l = self.basis.u_iso_curve(uu)?;
        let (_, d1u, d1v) = self.basis.d1(uu, 0.0);
        let mut dir = d1u.crossed(&d1v);
        if dir.square_magnitude() < RESOLUTION {
            return Some(a_l);
        }
        dir.normalize();
        dir = GpVec::new(
            dir.x() * self.offset,
            dir.y() * self.offset,
            dir.z() * self.offset,
        );
        Some(Arc::from(a_l.translated(&dir)))
    }

    fn approx_iso(&self, u_iso: bool, iso_par: f64) -> Option<Arc<dyn Curve>> {
        let (first, last) = if u_iso {
            self.v_range()
        } else {
            self.u_range()
        };
        if !first.is_finite() || !last.is_finite() {
            return None;
        }
        let surf = self.clone();
        let eval = |_a: f64, _b: f64, t: f64, deriv: i32, out: &mut [f64]| -> i32 {
            if deriv == 0 {
                let p = if u_iso {
                    surf.d0(iso_par, t)
                } else {
                    surf.d0(t, iso_par)
                };
                out[0] = p.x();
                out[1] = p.y();
                out[2] = p.z();
            } else {
                let (_, du, dv) = if u_iso {
                    surf.d1(iso_par, t)
                } else {
                    surf.d1(t, iso_par)
                };
                let d = if u_iso { dv } else { du };
                out[0] = d.x();
                out[1] = d.y();
                out[2] = d.z();
            }
            0
        };
        let approx = ApproxAFunction3d::approx_c1(first, last, 14, 100, APPROXIMATION, &eval).ok()?;
        if !approx.done && !approx.has_result {
            return None;
        }
        GeomBSplineCurve::from_poles_knots_mults(
            approx.poles,
            approx.knots,
            approx.mults,
            approx.degree.max(0) as usize,
        )
        .ok()
        .map(|c| Arc::new(c) as Arc<dyn Curve>)
    }
}

impl Surface for GeomOffsetSurface {
    fn d0(&self, u: f64, v: f64) -> GpPnt {
        // `Geom_OffsetSurfaceUtils::EvaluateD0` (`pxx:586-783`).
        let (p, d1u, d1v) = self.basis.d1(u, v);
        match evaluate_d0(
            u,
            v,
            self.basis.as_ref(),
            self.offset,
            Some(&self.osc),
            p,
            d1u,
            d1v,
        ) {
            Some(r) => r,
            None => {
                let basis_p = self.basis.d0(u, v);
                let n = self.unit_normal(u, v);
                GpPnt::from_xyz(&basis_p.coord.added(&n.xyz().multiplied(self.offset)))
            }
        }
    }

    fn d1(&self, u: f64, v: f64) -> (GpPnt, GpVec, GpVec) {
        // `Geom_OffsetSurfaceUtils::EvaluateD1` (`pxx:804-1095`).
        let (p, d1u, d1v, d2u, d2v, d2uv) = self.basis.d2(u, v);
        match evaluate_d1(
            u,
            v,
            self.basis.as_ref(),
            self.offset,
            Some(&self.osc),
            p,
            d1u,
            d1v,
            d2u,
            d2v,
            d2uv,
        ) {
            Some(r) => r,
            None => (self.d0(u, v), d1u, d1v),
        }
    }

    fn u_range(&self) -> (f64, f64) {
        self.basis.u_range()
    }

    fn v_range(&self) -> (f64, f64) {
        self.basis.v_range()
    }

    fn is_u_periodic(&self) -> bool {
        self.basis.is_u_periodic()
    }

    fn is_v_periodic(&self) -> bool {
        self.basis.is_v_periodic()
    }

    fn continuity(&self) -> u8 {
        self.basis.continuity()
    }

    /// `Geom_OffsetSurface::IsCNu` (`cxx:692-696`): the basis at order `N + 1`.
    fn is_cnu(&self, n: i32) -> Result<bool, String> {
        if n < 0 {
            return Err("Geom_OffsetSurface::IsCNu: Geom_OffsetSurface.cxx:694".into());
        }
        self.basis.is_cnu(n.checked_add(1).ok_or_else(|| {
            "Geom_OffsetSurface::IsCNu: Geom_OffsetSurface.cxx:695".to_string()
        })?)
    }

    /// `Geom_OffsetSurface::IsCNv` (`cxx:700-704`).
    fn is_cnv(&self, n: i32) -> Result<bool, String> {
        if n < 0 {
            return Err("Geom_OffsetSurface::IsCNv: Geom_OffsetSurface.cxx:702".into());
        }
        self.basis.is_cnv(n.checked_add(1).ok_or_else(|| {
            "Geom_OffsetSurface::IsCNv: Geom_OffsetSurface.cxx:703".to_string()
        })?)
    }

    fn transform(&mut self, _t: &GpTrsf) {
        // Transform is not carried into the basis (matches the simplified
        // derivative handling above).
    }

    fn clone_dyn(&self) -> Box<dyn Surface> {
        Box::new(self.clone())
    }

    fn is_offset_surface(&self) -> bool {
        true
    }

    fn offset_basis_surface(&self) -> Option<Arc<dyn Surface>> {
        Some(self.basis.clone())
    }

    fn offset_distance(&self) -> Option<f64> {
        Some(self.offset)
    }

    /// `Geom_Surface::UReversed` (`Geom_Surface.cxx:33-38`) =
    /// `Geom_OffsetSurface::UReverse` (`Geom_OffsetSurface.cxx:281-286`):
    /// the basis is U-reversed and the signed offset is negated. `None` when
    /// the basis itself has no `UReversed` branch.
    fn u_reversed(&self) -> Option<Arc<dyn Surface>> {
        let basis = self.basis.u_reversed()?;
        Some(Arc::new(Self::new(basis, -self.offset)))
    }

    /// `Geom_Surface::VReversed` (`Geom_Surface.cxx:43-48`) =
    /// `Geom_OffsetSurface::VReverse` (`Geom_OffsetSurface.cxx:297-302`):
    /// the basis is V-reversed and the signed offset is negated. `None` when
    /// the basis itself has no `VReversed` branch.
    fn v_reversed(&self) -> Option<Arc<dyn Surface>> {
        let basis = self.basis.v_reversed()?;
        Some(Arc::new(Self::new(basis, -self.offset)))
    }

    /// `Geom_OffsetSurface::UReversedParameter`
    /// (`Geom_OffsetSurface.cxx:290-294`): the basis surface's.
    fn u_reversed_parameter(&self, u: f64) -> f64 {
        self.basis.u_reversed_parameter(u)
    }

    /// `Geom_OffsetSurface::VReversedParameter`
    /// (`Geom_OffsetSurface.cxx:306-310`): the basis surface's.
    fn v_reversed_parameter(&self, v: f64) -> f64 {
        self.basis.v_reversed_parameter(v)
    }

    fn u_iso_curve(&self, u: f64) -> Option<Arc<dyn Curve>> {
        self.u_iso(u)
    }

    fn v_iso_curve(&self, v: f64) -> Option<Arc<dyn Curve>> {
        self.v_iso(v)
    }

    fn nb_u_poles(&self) -> i32 {
        self.basis.nb_u_poles()
    }
    fn nb_v_poles(&self) -> i32 {
        self.basis.nb_v_poles()
    }
    fn nb_u_intervals(&self, continuity: u8) -> i32 {
        // `GeomAdaptor_Surface::NbUIntervals` Offset (`cxx:664-684`).
        self.basis.nb_u_intervals(offset_interval_continuity(continuity))
    }
    fn nb_v_intervals(&self, continuity: u8) -> i32 {
        self.basis.nb_v_intervals(offset_interval_continuity(continuity))
    }
    fn u_intervals(&self, continuity: u8) -> Vec<f64> {
        self.basis.u_intervals(offset_interval_continuity(continuity))
    }
    fn v_intervals(&self, continuity: u8) -> Vec<f64> {
        self.basis.v_intervals(offset_interval_continuity(continuity))
    }
    fn extrusion_basis_curve(&self) -> Option<Arc<dyn Curve>> {
        self.basis.extrusion_basis_curve()
    }

    fn uv_resolution(&self, r3d: f64) -> Option<(f64, f64)> {
        // `GeomAdaptor_Surface::UResolution` Offset arm (`cxx:1882-1884`).
        self.basis.uv_resolution(r3d)
    }

    /// `Geom_OffsetSurface::EvalDN` (`cxx:466-496`).
    fn try_eval_dn(&self, u: f64, v: f64, nu: i32, nv: i32) -> Result<GpVec, String> {
        if nu + nv < 1 || nu < 0 || nv < 0 {
            return Err("Geom_OffsetSurface::EvalDN: Geom_OffsetSurface.cxx:470".into());
        }
        evaluate_dn(u, v, nu, nv, self.basis.as_ref(), self.offset, Some(&self.osc)).ok_or_else(
            || "Geom_OffsetSurface::EvalDN: Geom_OffsetSurface.cxx:494".to_string(),
        )
    }

    /// `Geom_OffsetSurface::EvalD3` (`cxx:428-461`): `D3UUV` is `(2,1)`, `D3UVV` is `(1,2)`.
    fn try_eval_d3_uv(&self, u: f64, v: f64, nu: i32, nv: i32) -> Result<GpVec, String> {
        let third = evaluate_d3(u, v, self.basis.as_ref(), self.offset, Some(&self.osc))
            .ok_or_else(|| "Geom_OffsetSurface::EvalD3: Geom_OffsetSurface.cxx:459".to_string())?;
        Ok(match (nu, nv) {
            (2, 1) => third.d3uuv,
            (1, 2) => third.d3uvv,
            (3, 0) => third.d3u,
            (0, 3) => third.d3v,
            _ => {
                return Err(
                    "Geom_OffsetSurface::EvalD3: Geom_OffsetSurface.cxx:459".into(),
                );
            }
        })
    }
}
