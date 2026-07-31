use occt_core::gp::{GpPln, GpPnt, GpVec, GpTrsf};
use crate::surface::Surface;
use occt_core::elib::slib;

#[derive(Debug, Clone)]
pub struct GeomPlane { pos: GpPln }

impl GeomPlane {
    pub fn new(pl: GpPln) -> Self { Self { pos: pl } }
    pub fn pln(&self) -> &GpPln { &self.pos }
}

impl Surface for GeomPlane {
    fn d0(&self, u: f64, v: f64) -> GpPnt { slib::plane_value(&self.pos, u, v) }
    fn d1(&self, u: f64, v: f64) -> (GpPnt, GpVec, GpVec) { slib::plane_d1(&self.pos, u, v) }
    fn u_range(&self) -> (f64, f64) { (f64::NEG_INFINITY, f64::INFINITY) }
    fn v_range(&self) -> (f64, f64) { (f64::NEG_INFINITY, f64::INFINITY) }
    fn continuity(&self) -> u8 { 3 }
    fn transform(&mut self, t: &GpTrsf) { self.pos.transform(t); }
    fn clone_dyn(&self) -> Box<dyn Surface> { Box::new(self.clone()) }
}
