//! `ChFi3d_SearchSing` — the singularity search function
//! `F(t) = (C1(t) - C2(t)) . (C1'(t) - C2'(t))`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_SearchSing.cxx:27-57` and
//! `.hxx:22-52`.

use std::sync::Arc;

use occt_core::gp::GpVec;
use occt_geom::curve::Curve;
use occt_math::function_set_root::FunctionWithDerivative;

/// `ChFi3d_SearchSing` (`hxx:22-52`): a `math_FunctionWithDerivative` over two
/// curves.
pub struct SearchSing {
    /// `myC1`
    my_c1: Arc<dyn Curve>,
    /// `myC2`
    my_c2: Arc<dyn Curve>,
}

impl SearchSing {
    /// `ChFi3d_SearchSing(C1, C2)` (`cxx:27-32`).
    pub fn new(c1: Arc<dyn Curve>, c2: Arc<dyn Curve>) -> Self {
        Self { my_c1: c1, my_c2: c2 }
    }
}

impl FunctionWithDerivative for SearchSing {
    /// `Value(X, F)` (`cxx:34-43`): `F = (P1P2) . (V2 - V1)`.
    fn value(&mut self, x: f64, f: &mut f64) -> bool {
        let (p1, v1) = self.my_c1.d1(x);
        let (p2, v2) = self.my_c2.d1(x);
        let v = GpVec::from_pnts(&p1, &p2);
        *f = v.dot(&v2.subtracted(&v1));
        true
    }

    /// `Derivative(X, D)` (`cxx:45-53`):
    /// `D = |V2 - V1|^2 + (P1P2) . (W2 - W1)`.
    fn derivative(&mut self, x: f64, d: &mut f64) -> bool {
        let (p1, v1, w1) = self.my_c1.d2(x);
        let (p2, v2, w2) = self.my_c2.d2(x);
        let v = GpVec::from_pnts(&p1, &p2);
        let v_prim = v2.subtracted(&v1);
        *d = v_prim.square_magnitude() + v.dot(&w2.subtracted(&w1));
        true
    }

    /// `Values(X, F, D)` (`cxx:55-57`): both at once, from the second
    /// derivatives.
    fn values(&mut self, x: f64, f: &mut f64, d: &mut f64) -> bool {
        let (p1, v1, w1) = self.my_c1.d2(x);
        let (p2, v2, w2) = self.my_c2.d2(x);
        let v = GpVec::from_pnts(&p1, &p2);
        let v_prim = v2.subtracted(&v1);
        *f = v.dot(&v_prim);
        *d = v_prim.square_magnitude() + v.dot(&w2.subtracted(&w1));
        true
    }
}
