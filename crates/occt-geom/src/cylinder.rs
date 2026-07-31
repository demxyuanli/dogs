use occt_core::gp::{GpCylinder, GpPnt, GpVec, GpTrsf};
use crate::surface::Surface;
use occt_core::elib::slib;

#[derive(Debug, Clone)]
pub struct GeomCylinder { pos: GpCylinder }

impl GeomCylinder {
    pub fn new(c: GpCylinder) -> Self { Self { pos: c } }
}

impl Surface for GeomCylinder {
    fn d0(&self, u: f64, v: f64) -> GpPnt { slib::cylinder_value(&self.pos, u, v) }
    fn d1(&self, u: f64, v: f64) -> (GpPnt, GpVec, GpVec) { (slib::cylinder_value(&self.pos,u,v), GpVec::zero(), GpVec::zero()) }
    fn u_range(&self) -> (f64, f64) { (0.0, 2.0*std::f64::consts::PI) }
    fn v_range(&self) -> (f64, f64) { (f64::NEG_INFINITY, f64::INFINITY) }
    fn is_u_periodic(&self) -> bool { true }
    fn continuity(&self) -> u8 { 3 }
    fn transform(&mut self, t: &GpTrsf) { self.pos.transform(t); }
    fn clone_dyn(&self) -> Box<dyn Surface> { Box::new(self.clone()) }
}
