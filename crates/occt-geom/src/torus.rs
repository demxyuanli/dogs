use std::sync::Arc;
use occt_core::elib::{slib, surface_eval};
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
    fn u_range(&self) -> (f64, f64) { (0.0, 2.0*std::f64::consts::PI) }
    fn v_range(&self) -> (f64, f64) { (0.0, 2.0*std::f64::consts::PI) }
    fn is_u_periodic(&self) -> bool { true }
    fn is_v_periodic(&self) -> bool { true }
    fn gp_torus(&self) -> Option<GpTorus> { Some(self.pos) }
    fn continuity(&self) -> u8 { 3 } // T-64: OCCT=GeomAbs_CN(6), blocked by consumer
    fn transform(&mut self, t: &GpTrsf) { self.pos.transform(t); }
    fn clone_dyn(&self) -> Box<dyn Surface> { Box::new(self.clone()) }

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
