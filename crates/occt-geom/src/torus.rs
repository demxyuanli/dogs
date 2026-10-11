use std::sync::Arc;
use occt_core::elib::{slib, surface_eval, surface_eval_d3};
use occt_core::gp::{GpTorus, GpPnt, GpVec, GpTrsf};
use crate::circle::GeomCircle;
use crate::curve::Curve;
use crate::surface::Surface;

#[derive(Debug, Clone)]
pub struct GeomTorus { pos: GpTorus }

impl GeomTorus {
    pub fn new(t: GpTorus) -> Self { Self { pos: t } }
    /// `Geom_ToroidalSurface::Torus`.
    pub fn torus(&self) -> &GpTorus { &self.pos }
}

impl Surface for GeomTorus {
    fn d0(&self, u: f64, v: f64) -> GpPnt { slib::torus_value(&self.pos, u, v) }
    fn d1(&self, u: f64, v: f64) -> (GpPnt, GpVec, GpVec) {
        surface_eval::torus_d1(&self.pos, u, v)
    }
    /// `Geom_ToroidalSurface::D2` → `ElSLib::TorusD2` (`ElSLib.cxx:1039-1100`).
    /// Replaces the trait's central-difference default (audit A15).
    fn d2(&self, u: f64, v: f64) -> (GpPnt, GpVec, GpVec, GpVec, GpVec, GpVec) {
        surface_eval::torus_d2(&self.pos, u, v)
    }
    /// `Geom_ToroidalSurface::EvalDN` (`Geom_ToroidalSurface.cxx:283-293`).
    fn try_eval_dn(&self, u: f64, v: f64, nu: i32, nv: i32) -> Result<GpVec, String> {
        if nu + nv < 1 || nu < 0 || nv < 0 {
            return Err("Geom_UndefinedDerivative".into());
        }
        Ok(surface_eval::torus_dn(
            u,
            v,
            self.pos.position(),
            self.pos.major_radius(),
            self.pos.minor_radius(),
            nu,
            nv,
        ))
    }
    /// `Geom_ToroidalSurface::EvalD3` (`Geom_ToroidalSurface.cxx:260-281`).
    fn try_eval_d3_uv(&self, u: f64, v: f64, nu: i32, nv: i32) -> Result<GpVec, String> {
        let d3 = surface_eval_d3::torus_d3(&self.pos, u, v);
        match (nu, nv) {
            (2, 1) => Ok(d3.d3_uuv),
            (1, 2) => Ok(d3.d3_uvv),
            _ => Ok(GpVec::zero()),
        }
    }
    fn u_range(&self) -> (f64, f64) { (0.0, 2.0*std::f64::consts::PI) }
    fn v_range(&self) -> (f64, f64) { (0.0, 2.0*std::f64::consts::PI) }
    fn is_u_periodic(&self) -> bool { true }
    fn is_v_periodic(&self) -> bool { true }
    fn gp_torus(&self) -> Option<GpTorus> { Some(self.pos) }
    fn continuity(&self) -> u8 { 6 }
    fn transform(&mut self, t: &GpTrsf) { self.pos.transform(t); }
    fn clone_dyn(&self) -> Box<dyn Surface> { Box::new(self.clone()) }

    /// `Geom_ToroidalSurface::UReversed` (`Geom_ToroidalSurface.cxx`):
    /// a copy with `UReverse()` applied (`gp_Torus.hxx:123` -> `pos.YReverse()`).
    fn u_reversed(&self) -> Option<Arc<dyn Surface>> {
        let mut pos = self.pos;
        pos.u_reverse();
        Some(Arc::new(GeomTorus::new(pos)))
    }

    /// `Geom_ToroidalSurface::UIso` (`Geom_ToroidalSurface.cxx:305-310`):
    /// `Geom_Circle(ElSLib::TorusUIso(pos, majorRadius, minorRadius, U))`,
    /// the minor circle whose plane contains the axis direction and the
    /// radial direction at `U`, centered `majorRadius` out along that radius.
    fn u_iso_curve(&self, u: f64) -> Option<Arc<dyn Curve>> {
        Some(Arc::new(GeomCircle::new(slib::torus_u_iso(
            &self.pos.pos,
            self.pos.major_radius,
            self.pos.minor_radius,
            u,
        ))))
    }

    /// `Geom_ToroidalSurface::VIso` (`Geom_ToroidalSurface.cxx:314-319`):
    /// `Geom_Circle(ElSLib::TorusVIso(pos, majorRadius, minorRadius, V))`,
    /// the parallel of radius `majorRadius + minorRadius*cos(V)`.
    fn v_iso_curve(&self, v: f64) -> Option<Arc<dyn Curve>> {
        Some(Arc::new(GeomCircle::new(slib::torus_v_iso(
            &self.pos.pos,
            self.pos.major_radius,
            self.pos.minor_radius,
            v,
        ))))
    }
}
