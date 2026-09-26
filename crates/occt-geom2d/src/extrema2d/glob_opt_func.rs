//! `Extrema_GlobOptFuncCCC0/1/2` -- the two-variable curve/curve distance
//! function that `Extrema_GGenExtCC::Perform` hands to `math_GlobOptMin`.
//!
//! Source: `Extrema_GlobOptFuncCC.cxx` -- the 2D statics `_Value` (`:49-65`),
//! `_Gradient` (`:96-119`) and `_Hessian` (`:156-185`), and the
//! `Extrema_GlobOptFuncCCC2(const Adaptor2d_Curve2d&, const Adaptor2d_Curve2d&)`
//! constructor plus its `Value`/`Gradient`/`Values` forwarders
//! (`:318-326`, `:337-347`, `:351-361`, `:372-388`);
//! `Extrema_GlobOptFuncCC.hxx:71-99`.
//!
//! OCCT declares three classes only so that `math_GlobOptMin::computeLocalExtremum`
//! can `dynamic_cast` its way to the best local engine (Hessian -> Newton,
//! gradient -> BFGS, value -> Powell; `math_GlobOptMin.cxx:266-339`). The maths is
//! the same in all three; the port's `GlobOptMin` takes a value closure, so one
//! struct with all three accessors is provided here and the engine selection
//! stays UNPORTED on the `GlobOptMin` side (see
//! `occt_math::globoptmin::GlobOptMin::set_continuity`).
//!
//! **OCCT quirk kept verbatim**: `_Value` returns the *distance*
//! (`C2.Value(v).Distance(C1.Value(u))`, 2D at `cxx:63`) while `_Gradient` and
//! `_Hessian` are the first/second derivatives of the *squared* distance (2D at
//! `cxx:113-116` and `cxx:174-183`, both scaled by 2). The minimizer is the same
//! (distance is monotone in squared distance), the scaling only affects step
//! sizes.

use super::prelude::*;

/// `Extrema_GlobOptFuncCCC2(C1, C2)` (`Extrema_GlobOptFuncCC.cxx:318-326`).
pub(super) struct GlobOptFuncCCC2<'a> {
    c1: &'a dyn Curve2d,
    c2: &'a dyn Curve2d,
}

impl<'a> GlobOptFuncCCC2<'a> {
    /// `Extrema_GlobOptFuncCCC2(const Adaptor2d_Curve2d&, const Adaptor2d_Curve2d&)`.
    pub(super) fn new(c1: &'a dyn Curve2d, c2: &'a dyn Curve2d) -> Self {
        Self { c1, c2 }
    }

    /// `_NbVariables()` (`cxx:24-27`).
    pub(super) fn nb_variables(&self) -> i32 {
        2
    }

    /// `_Value` (`cxx:49-65`): the distance `|C2(v) - C1(u)|`; `None` when the
    /// parameters leave the curves' ranges (`cxx:57-61`), which is how OCCT
    /// reports `false` to `math_GlobOptMin`.
    pub(super) fn value(&self, u: f64, v: f64) -> Option<f64> {
        if u < self.c1.first_parameter()
            || u > self.c1.last_parameter()
            || v < self.c2.first_parameter()
            || v > self.c2.last_parameter()
        {
            return None;
        }
        Some(self.c2.d0(v).distance(&self.c1.d0(u)))
    }

    /// `_Gradient` (`cxx:96-119`), the gradient of the squared distance:
    /// `G1 = -2 (C2 - C1)·C1'`, `G2 = +2 (C2 - C1)·C2'`.
    pub(super) fn gradient(&self, u: f64, v: f64) -> Option<(f64, f64)> {
        if u < self.c1.first_parameter()
            || u > self.c1.last_parameter()
            || v < self.c2.first_parameter()
            || v > self.c2.last_parameter()
        {
            return None;
        }
        let (p1, d1) = self.c1.d1(u);
        let (p2, d2) = self.c2.d1(v);
        let d = GpVec2d::new(p2.x() - p1.x(), p2.y() - p1.y());
        let g1 = -d.dot(&d1);
        let g2 = d.dot(&d2);
        Some((2.0 * g1, 2.0 * g2))
    }

    /// `Values(X, F, G)` = `Value && Gradient` (`cxx:365-371`).
    pub(super) fn values(&self, u: f64, v: f64) -> Option<(f64, (f64, f64))> {
        Some((self.value(u, v)?, self.gradient(u, v)?))
    }

    /// `_Hessian` (`cxx:156-185`), the Hessian of the squared distance:
    /// `H11 = 2(|C1'|^2 - (C2-C1)·C1'')`, `H12 = H21 = -2 C2'·C1'`,
    /// `H22 = 2(|C2'|^2 + (C2-C1)·C2'')`.
    pub(super) fn hessian(&self, u: f64, v: f64) -> Option<[[f64; 2]; 2]> {
        if u < self.c1.first_parameter()
            || u > self.c1.last_parameter()
            || v < self.c2.first_parameter()
            || v > self.c2.last_parameter()
        {
            return None;
        }
        let (p1, d1, d1b) = self.c1.d2(u);
        let (p2, d2, d2b) = self.c2.d2(v);
        let d = GpVec2d::new(p2.x() - p1.x(), p2.y() - p1.y());
        let h11 = 2.0 * (d1.square_magnitude() - d.dot(&d1b));
        let h12 = -2.0 * d2.dot(&d1);
        let h22 = 2.0 * (d2.square_magnitude() + d.dot(&d2b));
        Some([[h11, h12], [h12, h22]])
    }
}