//! 3D line curve. Source: `Geom_Line.hxx`
use occt_core::gp::{GpLin, GpPnt, GpVec, GpTrsf, GpAx1};
use crate::curve::Curve;
use occt_core::elib::clib;

#[derive(Debug, Clone)]
pub struct GeomLine { pos: GpLin }

impl GeomLine {
    pub fn new(l: GpLin) -> Self { Self { pos: l } }
    pub fn from_pnt_dir(p: GpPnt, d: occt_core::gp::GpDir) -> Self { Self { pos: GpLin::from_pnt_dir(p, d) } }
    pub fn lin(&self) -> &GpLin { &self.pos }
}

impl Curve for GeomLine {
    fn d0(&self, u: f64) -> GpPnt { clib::line_value(&self.pos, u) }
    fn d1(&self, u: f64) -> (GpPnt, GpVec) { clib::line_d1(&self.pos, u) }
    fn d2(&self, u: f64) -> (GpPnt, GpVec, GpVec) { clib::line_d2(&self.pos, u) }
    fn first_parameter(&self) -> f64 { f64::NEG_INFINITY }
    fn last_parameter(&self) -> f64 { f64::INFINITY }
    fn continuity(&self) -> u8 { 3 }
    fn transform(&mut self, t: &GpTrsf) { self.pos.transform(t); }
    fn reverse(&mut self) { self.pos.pos.vdir.reverse(); }
    fn clone_dyn(&self) -> Box<dyn Curve> { Box::new(self.clone()) }
}
