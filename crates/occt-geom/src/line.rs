//! 3D line curve. Source: `Geom_Line.hxx`
use occt_core::gp::{GpLin, GpPnt, GpVec, GpTrsf};
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
    /// `Geom_Line::EvalDN` (`Geom_Line.cxx:207-219`) → `ElCLib::LineDN`
    /// (`ElCLib.cxx:911-918`). The `N < 1` throw becomes a zero vector, the
    /// convention of `Curve::eval_dn` (`curve.rs:12-25`).
    fn eval_dn(&self, _u: f64, n: i32) -> GpVec {
        if n < 1 { return GpVec::zero(); }
        clib::line_dn(&self.pos, n)
    }
    fn first_parameter(&self) -> f64 { f64::NEG_INFINITY }
    fn last_parameter(&self) -> f64 { f64::INFINITY }
    fn continuity(&self) -> u8 { 6 }
    fn transform(&mut self, t: &GpTrsf) { self.pos.transform(t); }
    fn reverse(&mut self) { self.pos.pos.vdir.reverse(); }
    fn clone_dyn(&self) -> Box<dyn Curve> { Box::new(self.clone()) }
    fn is_line(&self) -> bool { true }
    fn gp_line(&self) -> Option<GpLin> { Some(self.pos.clone()) }
}
