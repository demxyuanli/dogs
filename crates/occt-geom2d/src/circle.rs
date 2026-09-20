//! 2D circle curve. Source: `Geom2d_Circle.hxx`
use occt_core::elib::clib;
use occt_core::gp::{GpCirc2d, GpPnt2d, GpVec2d, GpTrsf2d};
use crate::curve::Curve2d;

#[derive(Debug, Clone)]
pub struct Geom2dCircle { pos: GpCirc2d }

impl Geom2dCircle {
    pub fn new(c: GpCirc2d) -> Self { Self { pos: c } }
    pub fn from_center_radius(center: GpPnt2d, radius: f64) -> Self { Self { pos: GpCirc2d::new(GpCirc2d::default().pos, radius) } }
    pub fn circ(&self) -> &GpCirc2d { &self.pos }
}

impl Curve2d for Geom2dCircle {
    fn d0(&self, u: f64) -> GpPnt2d {
        // `Geom2d_Circle::EvalD0` (`cxx:164-167`): `ElCLib::CircleValue`.
        clib::circle2d_value(&self.pos, u)
    }
    fn d1(&self, u: f64) -> (GpPnt2d, GpVec2d) {
        // `Geom2d_Circle::EvalD1` (`cxx:171-176`): `ElCLib::CircleD1`.
        clib::circle2d_d1(&self.pos, u)
    }
    fn d2(&self, u: f64) -> (GpPnt2d, GpVec2d, GpVec2d) {
        // `Geom2d_Circle::EvalD2` (`cxx:180-185`): `ElCLib::CircleD2`.
        clib::circle2d_d2(&self.pos, u)
    }
    /// `Geom2d_Circle::EvalD3` (`cxx:189-194`): `ElCLib::CircleD3` (`cxx:809-839`).
    /// Needed by `Geom2d_OffsetCurve`'s `CalculateD2` (`D2Ndir`).
    fn d3(&self, u: f64) -> (GpPnt2d, GpVec2d, GpVec2d, GpVec2d) {
        clib::circle2d_d3(&self.pos, u)
    }
    fn first_parameter(&self) -> f64 { 0.0 }
    fn last_parameter(&self) -> f64 { 2.0 * std::f64::consts::PI }
    fn is_periodic(&self) -> bool { true }
    fn period(&self) -> f64 { 2.0 * std::f64::consts::PI }
    fn continuity(&self) -> u8 { 3 } // T-64: OCCT=GeomAbs_CN(6), blocked by consumer
    fn transform(&mut self, t: &GpTrsf2d) { self.pos.transform(t); }
    fn reverse(&mut self) { self.pos.radius = -self.pos.radius; }
    fn clone_dyn(&self) -> Box<dyn Curve2d> { Box::new(self.clone()) }
    fn gp_circ2d(&self) -> Option<occt_core::gp::GpCirc2d> { Some(self.pos) }
}
