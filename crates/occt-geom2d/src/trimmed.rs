//! Trimmed 2D curve. Source: `Geom2d_TrimmedCurve.hxx`
use std::sync::Arc;
use occt_core::gp::{GpPnt2d, GpVec2d, GpTrsf2d};
use crate::curve::Curve2d;

#[derive(Clone)]
pub struct Geom2dTrimmedCurve {
    basis: Arc<dyn Curve2d>,
    first: f64,
    last: f64,
}

impl Geom2dTrimmedCurve {
    pub fn new(curve: Arc<dyn Curve2d>, first: f64, last: f64) -> Self {
        let a = first.min(last); let b = first.max(last);
        Self { basis: curve, first: a, last: b }
    }
    pub fn basis_curve(&self) -> &Arc<dyn Curve2d> { &self.basis }
    pub fn trim_parameter(&self, u: f64) -> f64 { self.first + u * (self.last - self.first) }
}

impl Curve2d for Geom2dTrimmedCurve {
    fn d0(&self, u: f64) -> GpPnt2d { self.basis.d0(self.first + u * (self.last - self.first)) }
    fn d1(&self, u: f64) -> (GpPnt2d, GpVec2d) {
        let t = self.first + u * (self.last - self.first);
        let (p, d) = self.basis.d1(t);
        (p, GpVec2d::new(d.x() * (self.last - self.first), d.y() * (self.last - self.first)))
    }
    fn d2(&self, u: f64) -> (GpPnt2d, GpVec2d, GpVec2d) {
        let t = self.first + u * (self.last - self.first);
        let (p, d1, d2) = self.basis.d2(t);
        let s = self.last - self.first;
        (p, GpVec2d::new(d1.x()*s, d1.y()*s), GpVec2d::new(d2.x()*s*s, d2.y()*s*s))
    }
    fn first_parameter(&self) -> f64 { 0.0 }
    fn last_parameter(&self) -> f64 { 1.0 }
    fn continuity(&self) -> u8 { self.basis.continuity() }
    fn transform(&mut self, t: &GpTrsf2d) { /* basis is shared, cannot mutate */ }
    fn reverse(&mut self) { std::mem::swap(&mut self.first, &mut self.last); }
    fn clone_dyn(&self) -> Box<dyn Curve2d> { Box::new(self.clone()) }
}
