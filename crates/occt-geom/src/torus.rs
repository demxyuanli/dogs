use occt_core::gp::{GpTorus, GpPnt, GpVec, GpTrsf};
use crate::surface::Surface;
use occt_core::elib::slib;

#[derive(Debug, Clone)]
pub struct GeomTorus { pos: GpTorus }

impl GeomTorus {
    pub fn new(t: GpTorus) -> Self { Self { pos: t } }
}

impl Surface for GeomTorus {
    fn d0(&self, u: f64, v: f64) -> GpPnt { slib::torus_value(&self.pos, u, v) }
    fn d1(&self, u: f64, v: f64) -> (GpPnt, GpVec, GpVec) { (slib::torus_value(&self.pos,u,v), GpVec::zero(), GpVec::zero()) }
    fn u_range(&self) -> (f64, f64) { (0.0, 2.0*std::f64::consts::PI) }
    fn v_range(&self) -> (f64, f64) { (0.0, 2.0*std::f64::consts::PI) }
    fn is_u_periodic(&self) -> bool { true }
    fn is_v_periodic(&self) -> bool { true }
    fn continuity(&self) -> u8 { 3 }
    fn transform(&mut self, t: &GpTrsf) { self.pos.transform(t); }
    fn clone_dyn(&self) -> Box<dyn Surface> { Box::new(self.clone()) }
}
