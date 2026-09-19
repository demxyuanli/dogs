//! 3D circle curve. Source: `Geom_Circle.hxx`
use occt_core::gp::{GpCirc, GpPnt, GpVec, GpTrsf};
use crate::curve::Curve;
use occt_core::elib::clib;

#[derive(Debug, Clone)]
pub struct GeomCircle { pos: GpCirc }

impl GeomCircle {
    pub fn new(c: GpCirc) -> Self { Self { pos: c } }
    pub fn circ(&self) -> &GpCirc { &self.pos }
}

impl Curve for GeomCircle {
    fn d0(&self, u: f64) -> GpPnt { clib::circle_value(&self.pos, u) }
    fn d1(&self, u: f64) -> (GpPnt, GpVec) { clib::circle_d1(&self.pos, u) }
    fn d2(&self, u: f64) -> (GpPnt, GpVec, GpVec) { clib::circle_d2(&self.pos, u) }
    fn first_parameter(&self) -> f64 { 0.0 }
    fn last_parameter(&self) -> f64 { 2.0 * std::f64::consts::PI }
    fn is_periodic(&self) -> bool { true }
    fn period(&self) -> f64 { 2.0 * std::f64::consts::PI }
    fn circle_radius(&self) -> Option<f64> { Some(self.pos.radius.abs()) }
    fn gp_circ(&self) -> Option<GpCirc> { Some(self.pos.clone()) }
    fn continuity(&self) -> u8 { 3 } // T-64: OCCT=GeomAbs_CN(6), blocked by consumer
    fn transform(&mut self, t: &GpTrsf) { self.pos.transform(t); }
    fn reverse(&mut self) { self.pos.radius = -self.pos.radius; }
    fn clone_dyn(&self) -> Box<dyn Curve> { Box::new(self.clone()) }
}
