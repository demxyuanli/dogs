use std::sync::Arc;
use occt_core::elib::{slib, surface_eval};
use occt_core::gp::{GpCylinder, GpPnt, GpVec, GpTrsf};
use crate::circle::GeomCircle;
use crate::curve::Curve;
use crate::line::GeomLine;
use crate::surface::Surface;

#[derive(Debug, Clone)]
pub struct GeomCylinder { pos: GpCylinder }

impl GeomCylinder {
    pub fn new(c: GpCylinder) -> Self { Self { pos: c } }
    /// `Geom_CylindricalSurface::Cylinder`.
    pub fn cylinder(&self) -> &GpCylinder { &self.pos }
}

impl Surface for GeomCylinder {
    fn d0(&self, u: f64, v: f64) -> GpPnt { slib::cylinder_value(&self.pos, u, v) }
    fn d1(&self, u: f64, v: f64) -> (GpPnt, GpVec, GpVec) {
        surface_eval::cylinder_d1(&self.pos, u, v)
    }
    fn u_range(&self) -> (f64, f64) { (0.0, 2.0*std::f64::consts::PI) }
    fn v_range(&self) -> (f64, f64) { (f64::NEG_INFINITY, f64::INFINITY) }
    fn is_u_periodic(&self) -> bool { true }
    fn gp_cylinder(&self) -> Option<GpCylinder> { Some(self.pos.clone()) }
    fn continuity(&self) -> u8 { 3 }
    fn transform(&mut self, t: &GpTrsf) { self.pos.transform(t); }
    fn clone_dyn(&self) -> Box<dyn Surface> { Box::new(self.clone()) }

    /// `Geom_CylindricalSurface::UIso` (`Geom_CylindricalSurface.cxx:294-298`):
    /// `Geom_Line(ElSLib::CylinderUIso(pos, radius, U))`.
    /// `ElSLib::CylinderUIso` (`ElSLib.cxx:1716-1723`) is the generatrix
    /// through `CylinderValue(U, 0)` directed along the cylinder axis
    /// (`CylinderD1(U, 0).DV`, which is the `pos` main direction).
    fn u_iso_curve(&self, u: f64) -> Option<Arc<dyn Curve>> {
        Some(Arc::new(GeomLine::new(slib::cylinder_u_iso(
            &self.pos.pos,
            self.pos.radius,
            u,
        ))))
    }

    /// `Geom_CylindricalSurface::VIso` (`Geom_CylindricalSurface.cxx:302-306`):
    /// `Geom_Circle(ElSLib::CylinderVIso(pos, radius, V))`.
    /// `ElSLib::CylinderVIso` (`ElSLib.cxx:1781-1788`) is the circle of the
    /// placement plane translated by `V` along the axis.
    fn v_iso_curve(&self, v: f64) -> Option<Arc<dyn Curve>> {
        Some(Arc::new(GeomCircle::new(slib::cylinder_v_iso(
            &self.pos.pos,
            self.pos.radius,
            v,
        ))))
    }
}
