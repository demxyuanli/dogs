//! 2D circle curve. Source: `Geom2d_Circle.hxx`
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
        GpPnt2d::new(self.pos.location().x() + self.pos.radius * u.cos(), self.pos.location().y() + self.pos.radius * u.sin())
    }
    fn d1(&self, u: f64) -> (GpPnt2d, GpVec2d) {
        let r = self.pos.radius;
        let p = self.d0(u);
        (p, GpVec2d::new(-r * u.sin(), r * u.cos()))
    }
    fn d2(&self, u: f64) -> (GpPnt2d, GpVec2d, GpVec2d) {
        let r = self.pos.radius;
        let p = self.d0(u);
        (p, GpVec2d::new(-r * u.sin(), r * u.cos()), GpVec2d::new(-r * u.cos(), -r * u.sin()))
    }
    fn first_parameter(&self) -> f64 { 0.0 }
    fn last_parameter(&self) -> f64 { 2.0 * std::f64::consts::PI }
    fn is_periodic(&self) -> bool { true }
    fn period(&self) -> f64 { 2.0 * std::f64::consts::PI }
    fn continuity(&self) -> u8 { 3 }
    fn transform(&mut self, t: &GpTrsf2d) { self.pos.transform(t); }
    fn reverse(&mut self) { self.pos.radius = -self.pos.radius; }
    fn clone_dyn(&self) -> Box<dyn Curve2d> { Box::new(self.clone()) }
}
