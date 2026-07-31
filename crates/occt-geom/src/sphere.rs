use occt_core::gp::{GpSphere, GpPnt, GpVec, GpTrsf};
use crate::surface::Surface;
use occt_core::elib::slib;

#[derive(Debug, Clone)]
pub struct GeomSphere { pos: GpSphere }

impl GeomSphere {
    pub fn new(s: GpSphere) -> Self { Self { pos: s } }
}

impl Surface for GeomSphere {
    fn d0(&self, u: f64, v: f64) -> GpPnt { slib::sphere_value(&self.pos, u, v) }
    fn d1(&self, u: f64, v: f64) -> (GpPnt, GpVec, GpVec) { (slib::sphere_value(&self.pos,u,v), GpVec::zero(), GpVec::zero()) }
    fn u_range(&self) -> (f64, f64) { (0.0, 2.0*std::f64::consts::PI) }
    fn v_range(&self) -> (f64, f64) { (-std::f64::consts::FRAC_PI_2, std::f64::consts::FRAC_PI_2) }
    fn is_u_periodic(&self) -> bool { true }
    fn continuity(&self) -> u8 { 3 }
    fn transform(&mut self, t: &GpTrsf) { self.pos.transform(t); }
    fn clone_dyn(&self) -> Box<dyn Surface> { Box::new(self.clone()) }
}
