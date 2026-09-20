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
    fn d2(&self, u: f64) -> (GpPnt, GpVec, GpVec) { clib::hyperbola_d2(&self.pos, u) }
    /// `Geom_Hyperbola::EvalD3` → `ElCLib::HyperbolaD3` (`ElCLib.cxx:494-516`).
    /// Needed by `Geom_OffsetCurve`'s `CalculateD2` (`D2Ndir`).
    fn d3(&self, u: f64) -> (GpPnt, GpVec, GpVec, GpVec) { clib::hyperbola_d3(&self.pos, u) }
    /// `Geom_Hyperbola::EvalDN` (`Geom_Hyperbola.cxx:269-276`) →
    /// `ElCLib::HyperbolaDN` (`ElCLib.cxx:996-1016`). `N < 1` returns a zero
    /// vector (see `curve.rs:12-25`).
    fn eval_dn(&self, u: f64, n: i32) -> GpVec {
        if n < 1 { return GpVec::zero(); }
        clib::hyperbola_dn(&self.pos, u, n)
    }
    fn first_parameter(&self) -> f64 { f64::NEG_INFINITY }
    fn last_parameter(&self) -> f64 { f64::INFINITY }
    fn continuity(&self) -> u8 { 3 } // T-64: OCCT=GeomAbs_CN(6), blocked by consumer
    fn gp_hyperbola(&self) -> Option<GpHypr> { Some(self.pos.clone()) }
    fn transform(&mut self, t: &GpTrsf) { self.pos.transform(t); }
    fn reverse(&mut self) {}
    fn clone_dyn(&self) -> Box<dyn Curve> { Box::new(self.clone()) }
}
