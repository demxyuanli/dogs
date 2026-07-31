//! 3D hyperbola curve. Source: `Geom_Hyperbola.hxx`
use occt_core::gp::{GpHypr, GpPnt, GpVec, GpTrsf};
use crate::curve::Curve;
use occt_core::elib::clib;

#[derive(Debug, Clone)]
pub struct GeomHyperbola { pos: GpHypr }

impl GeomHyperbola {
    pub fn new(h: GpHypr) -> Self { Self { pos: h } }
}

impl Curve for GeomHyperbola {
    fn d0(&self, u: f64) -> GpPnt { clib::hyperbola_value(&self.pos, u) }
    fn d1(&self, u: f64) -> (GpPnt, GpVec) { clib::hyperbola_d1(&self.pos, u) }
    fn d2(&self, u: f64) -> (GpPnt, GpVec, GpVec) { (self.d0(u), self.d1(u).1, GpVec::zero()) }
    fn first_parameter(&self) -> f64 { f64::NEG_INFINITY }
    fn last_parameter(&self) -> f64 { f64::INFINITY }
    fn continuity(&self) -> u8 { 3 }
    fn transform(&mut self, t: &GpTrsf) { self.pos.transform(t); }
    fn reverse(&mut self) {}
    fn clone_dyn(&self) -> Box<dyn Curve> { Box::new(self.clone()) }
}
