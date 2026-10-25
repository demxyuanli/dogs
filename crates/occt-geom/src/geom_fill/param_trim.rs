//! `Geom_TrimmedCurve` over a basis parameter range, as `BRepFill` builds it.
//!
//! The parameter is the basis parameter. This is not
//! `occt_geom::GeomTrimmedCurve`, which remaps the span onto `[0, 1]`.

use std::sync::Arc;

use crate::curve::Curve;

/// `Geom_TrimmedCurve(C, U1, U2)` as `BRepFill` builds it (`BRepFill.cxx:244`).
pub struct ParameterTrim {
    basis: Arc<dyn Curve>,
    first: f64,
    last: f64,
}

impl ParameterTrim {
    pub fn new(basis: Arc<dyn Curve>, first: f64, last: f64) -> Self {
        Self { basis, first, last }
    }
}

impl Curve for ParameterTrim {
    fn d0(&self, u: f64) -> occt_core::gp::GpPnt {
        self.basis.d0(u)
    }
    fn d1(&self, u: f64) -> (occt_core::gp::GpPnt, occt_core::gp::GpVec) {
        self.basis.d1(u)
    }
    fn d2(&self, u: f64) -> (occt_core::gp::GpPnt, occt_core::gp::GpVec, occt_core::gp::GpVec) {
        self.basis.d2(u)
    }
    fn first_parameter(&self) -> f64 {
        self.first
    }
    fn last_parameter(&self) -> f64 {
        self.last
    }
    fn continuity(&self) -> u8 {
        self.basis.continuity()
    }
    fn is_geom_trimmed(&self) -> bool {
        true
    }
    fn untrimmed_basis(&self) -> Option<(Arc<dyn Curve>, f64, f64)> {
        match self.basis.untrimmed_basis() {
            Some((inner, _, _)) => Some((inner, self.first, self.last)),
            None => Some((self.basis.clone(), self.first, self.last)),
        }
    }
    fn transform(&mut self, t: &occt_core::gp::GpTrsf) {
        let mut basis = self.basis.clone_dyn();
        basis.transform(t);
        self.basis = Arc::from(basis);
    }
    /// `Geom_TrimmedCurve::Reverse` (`Geom_TrimmedCurve.cxx:88-96`).
    fn reverse(&mut self) {
        let first = self.basis.reversed_parameter(self.last);
        let last = self.basis.reversed_parameter(self.first);
        let mut basis = self.basis.clone_dyn();
        basis.reverse();
        self.basis = Arc::from(basis);
        self.first = first;
        self.last = last;
    }
    fn clone_dyn(&self) -> Box<dyn Curve> {
        Box::new(ParameterTrim {
            basis: Arc::from(self.basis.clone_dyn()),
            first: self.first,
            last: self.last,
        })
    }
}
