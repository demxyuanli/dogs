use occt_core::elib::{slib, surface_eval};
use occt_core::gp::{GpCylinder, GpPnt, GpVec, GpTrsf};
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
}
