//! Offset 3D surface. Source: `Geom_OffsetSurface.hxx`
//!
//! The offset surface at `(u, v)` is the basis surface point displaced along
//! its unit normal by the signed offset: `d0 = basis.d0 + offset · n`.
//! First partials follow `Geom_OffsetSurfaceUtils::EvaluateD1`.
use std::sync::Arc;

use occt_core::gp::{GpPnt, GpTrsf, GpVec};
use occt_core::precision::{APPROXIMATION, RESOLUTION};

use crate::adv_approx::ApproxAFunction3d;
use crate::bspline_curve::GeomBSplineCurve;
use crate::curve::Curve;
use crate::offset_surface_utils::evaluate_d1;
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
        let eval = |t: f64, deriv: i32, out: &mut [f64]| -> i32 {
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
        let p = self.basis.d0(u, v);
        let n = self.unit_normal(u, v);
        GpPnt::from_xyz(&p.coord.added(&n.xyz().multiplied(self.offset)))
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
}
