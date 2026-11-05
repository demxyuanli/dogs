//! `Geom2dConvert_PPoint` — a point of a 2D adaptor curve with its parameter
//! and first derivative.
//!
//! Source: `ModelingData/TKGeomBase/Geom2dConvert/Geom2dConvert_PPoint.hxx:31-79`,
//! `Geom2dConvert_PPoint.cxx:22-47`.

use occt_core::gp::{GpPnt2d, GpVec2d};
use occt_core::precision::PCONFUSION;

use crate::curve::Curve2d;

/// `RealLast()` (`Standard_Real.hxx:179-182`).
const REAL_LAST: f64 = f64::MAX;

/// `Geom2dConvert_PPoint` (`hxx:31-79`): a point on the adaptor curve, with
/// its parameter and first derivative.
#[derive(Debug, Clone, Copy)]
pub struct PPoint {
    param: f64,
    point: GpPnt2d,
    d1: GpVec2d,
}

impl PPoint {
    /// Empty constructor (`hxx:33-39`): `myParameter = RealLast()`, both
    /// coordinates zero.
    pub fn empty() -> Self {
        Self {
            param: REAL_LAST,
            point: GpPnt2d::zero(),
            d1: GpVec2d::zero(),
        }
    }

    /// `Geom2dConvert_PPoint(theParameter, thePoint, theD1)` (`hxx:42-48`).
    pub fn new(param: f64, point: GpPnt2d, d1: GpVec2d) -> Self {
        Self { param, point, d1 }
    }

    /// `Geom2dConvert_PPoint(theParameter, theAdaptor)` (`cxx:22-28`):
    /// `theAdaptor.D1(theParameter, myPoint, myD1)`.
    pub fn on_curve(param: f64, c: &dyn Curve2d) -> Self {
        let (point, d1) = c.d1(param);
        Self { param, point, d1 }
    }

    /// `Dist` (`hxx:51-54`): the distance between the two **points**, not the
    /// parameters.
    pub fn dist(&self, other: &Self) -> f64 {
        self.point.distance(&other.point)
    }

    /// `Parameter()`.
    pub fn parameter(&self) -> f64 {
        self.param
    }

    /// `Point()`.
    pub fn point(&self) -> GpPnt2d {
        self.point
    }

    /// `D1()`.
    pub fn d1(&self) -> GpVec2d {
        self.d1
    }

    /// `SetD1` (`hxx:71`).
    pub fn set_d1(&mut self, d1: GpVec2d) {
        self.d1 = d1;
    }
}

/// `operator==` / `operator!=` (`cxx:38-47`): parameter equality within
/// `Precision::PConfusion()`.
impl PartialEq for PPoint {
    fn eq(&self, other: &Self) -> bool {
        (self.param - other.param).abs() <= PCONFUSION
    }
}
