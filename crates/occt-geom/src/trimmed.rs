//! Trimmed 3D curve. Source: `Geom_TrimmedCurve.hxx`
use std::sync::Arc;
use occt_core::gp::{GpPnt, GpVec, GpTrsf};
use crate::curve::Curve;

#[derive(Clone)]
pub struct GeomTrimmedCurve {
    basis: Arc<dyn Curve>,
    first: f64,
    last: f64,
}

impl GeomTrimmedCurve {
    pub fn new(curve: Arc<dyn Curve>, first: f64, last: f64) -> Self {
        let a = first.min(last); let b = first.max(last);
        Self { basis: curve, first: a, last: b }
    }
    pub fn basis_curve(&self) -> &Arc<dyn Curve> { &self.basis }
}

impl Curve for GeomTrimmedCurve {
    fn d0(&self, u: f64) -> GpPnt { self.basis.d0(self.first + u * (self.last - self.first)) }
    fn d1(&self, u: f64) -> (GpPnt, GpVec) {
        let t = self.first + u * (self.last - self.first);
        let (p, d) = self.basis.d1(t);
        let s = self.last - self.first;
        (p, GpVec::new(d.x()*s, d.y()*s, d.z()*s))
    }
    fn d2(&self, u: f64) -> (GpPnt, GpVec, GpVec) {
        let t = self.first + u * (self.last - self.first);
        let (p, d1, d2) = self.basis.d2(t);
        let s = self.last - self.first;
        (p, GpVec::new(d1.x()*s,d1.y()*s,d1.z()*s), GpVec::new(d2.x()*s*s,d2.y()*s*s,d2.z()*s*s))
    }
    fn first_parameter(&self) -> f64 { 0.0 }
    fn last_parameter(&self) -> f64 { 1.0 }
    fn continuity(&self) -> u8 { self.basis.continuity() }
    fn circle_radius(&self) -> Option<f64> { self.basis.circle_radius() }
    fn gp_circ(&self) -> Option<occt_core::gp::GpCirc> { self.basis.gp_circ() }
    fn is_geom_trimmed(&self) -> bool { true }
    fn trimmed_basis_range(&self) -> Option<(f64, f64)> {
        // `GeomAdaptor_Curve::load` (`cxx:252-254`) unwraps nested trims.
        if let Some((bf, bl)) = self.basis.trimmed_basis_range() {
            Some((
                bf + self.first * (bl - bf),
                bf + self.last * (bl - bf),
            ))
        } else {
            Some((self.first, self.last))
        }
    }
    fn is_line(&self) -> bool { self.basis.is_line() }
    fn bspline_poles(&self) -> Option<&[GpPnt]> { self.basis.bspline_poles() }
    fn bezier_poles(&self) -> Option<&[GpPnt]> { self.basis.bezier_poles() }
    fn transform(&mut self, _t: &GpTrsf) {}
    fn reverse(&mut self) { std::mem::swap(&mut self.first, &mut self.last); }
    fn clone_dyn(&self) -> Box<dyn Curve> { Box::new(self.clone()) }
}
