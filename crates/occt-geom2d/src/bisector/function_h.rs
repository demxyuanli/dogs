//! `Bisector_FunctionH` (`Bisector_FunctionH.cxx:19-85`, `.hxx:30-58`).
//!
//! The scalar function whose zero gives the parameter `v` on the second curve
//! of a two-curve bisector:
//!
//! ```text
//! F(v) = P1P2(v) . ( ||T2(v)|| * T1 + T2(v) )
//! ```
//!
//! It is the `math_FunctionWithDerivative` that `Bisector_BisecCC` feeds to
//! `math_BissecNewton` (`cxx:684`) and to `math_FunctionRoot`
//! (`cxx:691`) / `math_FunctionRoots` (`cxx:867`).
//!
//! Two Rust traits model the same OCCT interface: `occt_core`'s
//! `MathFunctionWithDerivative` (used by `FunctionRoots`) and `occt-math`'s
//! `FunctionWithDerivative` (used by `MathFunctionRoot`). Both are implemented
//! here by forwarding to the same inherent helpers.

use std::sync::Arc;

use occt_core::gp::{GpPnt2d, GpVec2d};
use occt_core::math_fn::{MathFunction, MathFunctionWithDerivative};

use crate::curve::Curve2d;

/// `class Bisector_FunctionH` (`Bisector_FunctionH.hxx:30-58`).
#[derive(Clone)]
pub struct BisectorFunctionH {
    /// `curve2` (`.hxx:55`).
    curve2: Arc<dyn Curve2d>,
    /// `p1` (`.hxx:53`).
    p1: GpPnt2d,
    /// `t1` (`.hxx:54`), normalized by the constructor.
    t1: GpVec2d,
}

impl BisectorFunctionH {
    /// `Bisector_FunctionH(C2, P1, T1)` (`cxx:21-29`).
    ///
    /// `t1.Normalize()` raises `Standard_ConstructionError` on a null vector in
    /// OCCT; the port keeps the null vector (the error is dropped) so that the
    /// evaluation stays total, exactly like `gp_Vec2d::Normalize` callers
    /// elsewhere in this crate.
    pub fn new(c2: Arc<dyn Curve2d>, p1: GpPnt2d, t1: GpVec2d) -> Self {
        let mut t1 = t1;
        let _ = t1.normalize();
        Self {
            curve2: c2,
            p1,
            t1,
        }
    }

    /// `Value(X, F)` body (`cxx:35-50`).
    fn eval_value(&self, x: f64) -> f64 {
        let (p2, t2) = self.curve2.d1(x);
        let norm_t2 = t2.magnitude();
        let ax = norm_t2 * self.t1.x() - t2.x();
        let ay = norm_t2 * self.t1.y() - t2.y();
        (self.p1.x() - p2.x()) * ax + (self.p1.y() - p2.y()) * ay
    }

    /// `Values(X, F, D)` body (`cxx:58-84`).
    fn eval_values(&self, x: f64) -> (f64, f64) {
        let (p2, t2, t2v) = self.curve2.d2(x);
        let norm_t2 = t2.magnitude();
        let ax = norm_t2 * self.t1.x() - t2.x();
        let ay = norm_t2 * self.t1.y() - t2.y();
        let f = (self.p1.x() - p2.x()) * ax + (self.p1.y() - p2.y()) * ay;
        let scal = t2.dot(&t2v) / norm_t2;
        let d_ax = scal * self.t1.x() - t2v.x();
        let d_ay = scal * self.t1.y() - t2v.y();
        let d = -t2.x() * ax - t2.y() * ay
            + (self.p1.x() - p2.x()) * d_ax
            + (self.p1.y() - p2.y()) * d_ay;
        (f, d)
    }
}

impl MathFunction for BisectorFunctionH {
    /// `Bisector_FunctionH::Value` (`cxx:35-50`): always returns `true`.
    fn value(&mut self, x: f64, f: &mut f64) -> bool {
        *f = self.eval_value(x);
        true
    }
}

impl MathFunctionWithDerivative for BisectorFunctionH {
    /// `Bisector_FunctionH::Derivative` (`cxx:54-58`): forwards to `Values`.
    fn derivative(&mut self, x: f64, d: &mut f64) -> bool {
        let (_, dv) = self.eval_values(x);
        *d = dv;
        true
    }

    /// `Bisector_FunctionH::Values` (`cxx:62-85`).
    fn values(&mut self, x: f64, f: &mut f64, d: &mut f64) -> bool {
        let (fv, dv) = self.eval_values(x);
        *f = fv;
        *d = dv;
        true
    }
}

impl occt_math::FunctionWithDerivative for BisectorFunctionH {
    fn value(&mut self, x: f64, f: &mut f64) -> bool {
        *f = self.eval_value(x);
        true
    }

    fn derivative(&mut self, x: f64, d: &mut f64) -> bool {
        let (_, dv) = self.eval_values(x);
        *d = dv;
        true
    }

    fn values(&mut self, x: f64, f: &mut f64, d: &mut f64) -> bool {
        let (fv, dv) = self.eval_values(x);
        *f = fv;
        *d = dv;
        true
    }
}
