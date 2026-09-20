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
    fn d2(&self, u: f64) -> (GpPnt, GpVec, GpVec) { clib::parabola_d2(&self.pos, u) }
    // `Geom_Parabola::EvalD3` (`Geom_Parabola.cxx:195-200`) sets `V3 = (0,0,0)`
    // — the trait default already returns a zero third derivative.
    /// `Geom_Parabola::EvalDN` (`Geom_Parabola.cxx:205-212`) →
    /// `ElCLib::ParabolaDN` (`ElCLib.cxx:1020-1045`). `N < 1` returns a zero
    /// vector (see `curve.rs:12-25`).
    fn eval_dn(&self, u: f64, n: i32) -> GpVec {
        if n < 1 { return GpVec::zero(); }
        clib::parabola_dn(&self.pos, u, n)
    }
    fn first_parameter(&self) -> f64 { f64::NEG_INFINITY }
    fn last_parameter(&self) -> f64 { f64::INFINITY }
    fn continuity(&self) -> u8 { 6 }
    fn gp_parabola(&self) -> Option<GpParab> { Some(self.pos.clone()) }
    fn transform(&mut self, t: &GpTrsf) { self.pos.transform(t); }
    fn reverse(&mut self) { self.pos.focal = -self.pos.focal; }
    fn clone_dyn(&self) -> Box<dyn Curve> { Box::new(self.clone()) }
}
