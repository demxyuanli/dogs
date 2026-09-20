//! 3D parabola curve. Source: `Geom_Parabola.hxx`
use occt_core::gp::{GpParab, GpPnt, GpVec, GpTrsf};
use crate::curve::Curve;
use occt_core::elib::clib;

#[derive(Debug, Clone)]
pub struct GeomParabola { pos: GpParab }

impl GeomParabola {
    pub fn new(p: GpParab) -> Self { Self { pos: p } }
}

impl Curve for GeomParabola {
    fn d0(&self, u: f64) -> GpPnt { clib::parabola_value(&self.pos, u) }
    fn d1(&self, u: f64) -> (GpPnt, GpVec) { clib::parabola_d1(&self.pos, u) }
    fn d2(&self, u: f64) -> (GpPnt, GpVec, GpVec) { (self.d0(u), self.d1(u).1, GpVec::zero()) }
    fn first_parameter(&self) -> f64 { f64::NEG_INFINITY }
    fn last_parameter(&self) -> f64 { f64::INFINITY }
    fn continuity(&self) -> u8 { 3 } // T-64: OCCT=GeomAbs_CN(6), blocked by consumer
    fn gp_parabola(&self) -> Option<GpParab> { Some(self.pos.clone()) }
    fn transform(&mut self, t: &GpTrsf) { self.pos.transform(t); }
    fn reverse(&mut self) { self.pos.focal = -self.pos.focal; }
    fn clone_dyn(&self) -> Box<dyn Curve> { Box::new(self.clone()) }
}
