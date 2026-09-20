//! 3D ellipse curve. Source: `Geom_Ellipse.hxx`
use occt_core::gp::{GpElips, GpPnt, GpVec, GpTrsf};
use crate::curve::Curve;
use occt_core::elib::clib;

#[derive(Debug, Clone)]
pub struct GeomEllipse { pos: GpElips }

impl GeomEllipse {
    pub fn new(e: GpElips) -> Self { Self { pos: e } }
    pub fn elips(&self) -> &GpElips { &self.pos }
}

impl Curve for GeomEllipse {
    fn d0(&self, u: f64) -> GpPnt { clib::ellipse_value(&self.pos, u) }
    fn d1(&self, u: f64) -> (GpPnt, GpVec) { clib::ellipse_d1(&self.pos, u) }
    fn d2(&self, u: f64) -> (GpPnt, GpVec, GpVec) { clib::ellipse_d2(&self.pos, u) }
    /// `Geom_Ellipse::EvalD3` → `ElCLib::EllipseD3` (`ElCLib.cxx:464-491`).
    /// Needed by `Geom_OffsetCurve`'s `CalculateD2` (`D2Ndir`).
    fn d3(&self, u: f64) -> (GpPnt, GpVec, GpVec, GpVec) { clib::ellipse_d3(&self.pos, u) }
    fn first_parameter(&self) -> f64 { 0.0 }
    fn last_parameter(&self) -> f64 { 2.0 * std::f64::consts::PI }
    fn is_periodic(&self) -> bool { true }
    fn period(&self) -> f64 { 2.0 * std::f64::consts::PI }
    fn gp_ellipse(&self) -> Option<GpElips> { Some(self.pos.clone()) }
    fn continuity(&self) -> u8 { 3 } // T-64: OCCT=GeomAbs_CN(6), blocked by consumer
    fn transform(&mut self, t: &GpTrsf) { self.pos.transform(t); }
    fn reverse(&mut self) { std::mem::swap(&mut self.pos.major_radius, &mut self.pos.minor_radius); }
    fn clone_dyn(&self) -> Box<dyn Curve> { Box::new(self.clone()) }
}
