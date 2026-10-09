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
    /// `Geom_Circle::EvalD3` → `ElCLib::CircleD3` (`ElCLib.cxx:435-460`).
    /// Needed by `Geom_OffsetCurve`'s `CalculateD2` (`D2Ndir`).
    fn d3(&self, u: f64) -> (GpPnt, GpVec, GpVec, GpVec) { clib::circle_d3(&self.pos, u) }
    /// `Geom_Circle::EvalDN` (`Geom_Circle.cxx:184-191`) → `ElCLib::CircleDN`
    /// (`ElCLib.cxx:922-953`), orders above 3 included. `N < 1` returns a zero
    /// vector instead of the OCCT throw (see `curve.rs:12-25`).
    fn eval_dn(&self, u: f64, n: i32) -> GpVec {
        if n < 1 { return GpVec::zero(); }
        clib::circle_dn(&self.pos, u, n)
    }
    fn first_parameter(&self) -> f64 { 0.0 }
    fn last_parameter(&self) -> f64 { 2.0 * std::f64::consts::PI }
    fn is_periodic(&self) -> bool { true }
    fn period(&self) -> f64 { 2.0 * std::f64::consts::PI }
    fn circle_radius(&self) -> Option<f64> { Some(self.pos.radius.abs()) }
    fn gp_circ(&self) -> Option<GpCirc> { Some(self.pos.clone()) }
    fn continuity(&self) -> u8 { 6 }
    /// `GeomAdaptor_Curve::Resolution`'s `GeomAbs_Circle` arm
    /// (`GeomAdaptor_Curve.cxx:1122-1132`).
    fn resolution(&self, r3d: f64) -> f64 {
        let r = self.pos.radius.abs();
        if r > r3d / 2.0 {
            2.0 * (r3d / (2.0 * r)).asin()
        } else {
            2.0 * std::f64::consts::PI
        }
    }
    fn transform(&mut self, t: &GpTrsf) { self.pos.transform(t); }
    /// `Geom_Circle::Reverse` -> `Geom_Conic::Reverse` (`Geom_Conic.cxx:23-28`):
    /// reverse the frame's main direction. `ReversedParameter` is the affine
    /// form (`reversed_parameter` default), i.e. `2*pi - U` here.
    fn reverse(&mut self) {
        let z = self.pos.pos.direction().reversed();
        self.pos.pos.set_direction(z);
    }
    fn clone_dyn(&self) -> Box<dyn Curve> { Box::new(self.clone()) }
}
