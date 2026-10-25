//! `GeomFill_SnglrFunc` — the `(ratio * C') x C''` curve whose zeroes are the
//! singular points of the Frenet trihedron.
//!
//! Source: `TKGeomAlgo/GeomFill/GeomFill_SnglrFunc.hxx` and
//! `GeomFill_SnglrFunc.cxx:23-159` (OCCT 8.0.0).
//!
//! Translation notes:
//! - The OCCT `EvalD0`/`EvalD1`/`EvalD2`/`EvalD3`/`EvalDN` split becomes the
//!   [`Curve`] `d0`/`d1`/`d2`/`d3`/`eval_dn` overrides; every one of them
//!   reads the same `D1`/`D2`/`D3` (and `DN(4)`, `DN(5)`) derivatives of
//!   `myHCurve`.
//! - `Standard_RangeError` for `N < 1` and `Standard_NotImplemented` for
//!   `N > 3` (`cxx:132-146`) cannot be raised from [`Curve::eval_dn`], whose
//!   other implementations return a zero vector for the same orders, so this
//!   class follows that convention.
//! - `GetType()` is `GeomAbs_OtherCurve` (`cxx:155-158`), so every type query
//!   keeps the [`Curve`] default.

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpTrsf, GpVec};
use occt_core::kernel::geomabs::Shape;

use crate::curve::Curve;

/// `GeomFill_SnglrFunc` (`GeomFill_SnglrFunc.hxx:33-77`).
#[derive(Clone)]
pub struct SnglrFunc {
    /// `myHCurve` (`hxx:72`).
    curve: Arc<dyn Curve>,
    /// `ratio` (`hxx:73`).
    ratio: f64,
}

impl SnglrFunc {
    /// `GeomFill_SnglrFunc` (`cxx:23-27`): `ratio = 1`.
    pub fn new(curve: Arc<dyn Curve>) -> Self {
        Self { curve, ratio: 1.0 }
    }

    /// `GeomFill_SnglrFunc::SetRatio` (`cxx:39-42`).
    pub fn set_ratio(&mut self, ratio: f64) {
        self.ratio = ratio;
    }

    /// `ratio` (`hxx:73`).
    pub fn ratio(&self) -> f64 {
        self.ratio
    }

    /// `GeomFill_SnglrFunc::ShallowCopy` (`cxx:31-37`).
    pub fn shallow_copy(&self) -> Self {
        self.clone()
    }
}

impl Curve for SnglrFunc {
    /// `GeomFill_SnglrFunc::EvalD0` (`cxx:100-104`):
    /// `(D1 * ratio) ^ D2`.
    fn d0(&self, u: f64) -> GpPnt {
        let (_, d1, d2) = self.curve.d2(u);
        GpPnt::from_xyz(&d1.multiplied_scalar(self.ratio).crossed(&d2).xyz())
    }

    /// `GeomFill_SnglrFunc::EvalD1` (`cxx:106-113`).
    fn d1(&self, u: f64) -> (GpPnt, GpVec) {
        let (_, d1, d2, d3) = self.curve.d3(u);
        let dc = d1.multiplied_scalar(self.ratio);
        (
            GpPnt::from_xyz(&dc.crossed(&d2).xyz()),
            dc.crossed(&d3),
        )
    }

    /// `GeomFill_SnglrFunc::EvalD2` (`cxx:115-121`).
    fn d2(&self, u: f64) -> (GpPnt, GpVec, GpVec) {
        let (_, d1, d2, d3) = self.curve.d3(u);
        let d4 = self.curve.eval_dn(u, 4);
        (
            GpPnt::from_xyz(&d1.crossed(&d2).multiplied_scalar(self.ratio).xyz()),
            d1.crossed(&d3).multiplied_scalar(self.ratio),
            d2.crossed(&d3).added(&d1.crossed(&d4)).multiplied_scalar(self.ratio),
        )
    }

    /// `GeomFill_SnglrFunc::EvalD3` (`cxx:123-131`).
    fn d3(&self, u: f64) -> (GpPnt, GpVec, GpVec, GpVec) {
        let (_, d1, d2, d3) = self.curve.d3(u);
        let d4 = self.curve.eval_dn(u, 4);
        let d5 = self.curve.eval_dn(u, 5);
        (
            GpPnt::from_xyz(&d1.crossed(&d2).multiplied_scalar(self.ratio).xyz()),
            d1.crossed(&d3).multiplied_scalar(self.ratio),
            d2.crossed(&d3).added(&d1.crossed(&d4)).multiplied_scalar(self.ratio),
            d1
                .crossed(&d5)
                .added(&d2.crossed(&d4).multiplied_scalar(2.0))
                .multiplied_scalar(self.ratio),
        )
    }

    /// `GeomFill_SnglrFunc::EvalDN` (`cxx:133-149`).
    fn eval_dn(&self, u: f64, n: i32) -> GpVec {
        match n {
            1 => self.d1(u).1,
            2 => self.d2(u).2,
            3 => self.d3(u).3,
            _ => GpVec::zero(),
        }
    }

    /// `GeomFill_SnglrFunc::FirstParameter` (`cxx:44-47`).
    fn first_parameter(&self) -> f64 {
        self.curve.first_parameter()
    }

    /// `GeomFill_SnglrFunc::LastParameter` (`cxx:49-52`).
    fn last_parameter(&self) -> f64 {
        self.curve.last_parameter()
    }

    /// `GeomFill_SnglrFunc::IsPeriodic` (`cxx:92-95`).
    fn is_periodic(&self) -> bool {
        self.curve.is_periodic()
    }

    /// `GeomFill_SnglrFunc::Period` (`cxx:97-100`).
    fn period(&self) -> f64 {
        self.curve.period()
    }

    /// `GeomFill_SnglrFunc` has no `Continuity`; the OCCT class derives every
    /// interval query from `myHCurve` (`cxx:57-88`), so the mapped shape is
    /// reported here.
    fn continuity(&self) -> u8 {
        self.curve.continuity()
    }

    /// `GeomFill_SnglrFunc::NbIntervals` (`cxx:57-70`): the requested shape is
    /// raised to `C2` / `C3` / `CN` before it reaches `myHCurve`.
    fn nb_intervals(&self, s: u8) -> i32 {
        self.curve.nb_intervals(snglr_shape(s) as u8)
    }

    /// `GeomFill_SnglrFunc::Intervals` (`cxx:72-88`).
    fn parameter_intervals(&self, s: u8) -> Vec<f64> {
        self.curve.parameter_intervals(snglr_shape(s) as u8)
    }

    /// `GeomFill_SnglrFunc::Resolution` (`cxx:151-154`):
    /// `Precision::Parametric(R3D)`.
    fn resolution(&self, r3d: f64) -> f64 {
        occt_core::precision::Precision::parametric(r3d, occt_core::precision::CONFUSION)
    }

    // ponytail: adapter over an immutable Arc, like `GeomFill_SnglrFunc`
    // itself, which has no `Transform` / `Reverse`.
    fn transform(&mut self, _t: &GpTrsf) {}

    fn reverse(&mut self) {}

    fn clone_dyn(&self) -> Box<dyn Curve> {
        Box::new(self.clone())
    }
}

/// The `GeomAbs_Shape` mapping shared by `NbIntervals` and `Intervals`
/// (`GeomFill_SnglrFunc.cxx:59-68`, `cxx:74-85`).
fn snglr_shape(s: u8) -> Shape {
    let s = shape_from_u8(s);
    if s == Shape::C0 {
        Shape::C2
    } else if s == Shape::C1 {
        Shape::C3
    } else if s as u8 >= Shape::C2 as u8 {
        Shape::CN
    } else {
        // `HCS` starts at `GeomAbs_C0` and `G1` / `G2` match no arm
        // (`GeomFill_SnglrFunc.cxx:57-68`).
        Shape::C0
    }
}

/// Decode the `u8` continuity used by [`crate::curve::Curve`] into
/// `GeomAbs_Shape` (the discriminants of `occt_core::kernel::geomabs::Shape`).
pub fn shape_from_u8(s: u8) -> Shape {
    match s {
        0 => Shape::C0,
        1 => Shape::G1,
        2 => Shape::C1,
        3 => Shape::G2,
        4 => Shape::C2,
        5 => Shape::C3,
        _ => Shape::CN,
    }
}
