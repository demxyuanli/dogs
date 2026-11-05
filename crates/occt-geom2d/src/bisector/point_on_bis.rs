//! `Bisector_PointOnBis` (`Bisector_PointOnBis.cxx:18-139`, `.hxx:30-75`): one
//! vertex of the `Bisector_PolyBis` polygon, carrying the parameter on each of
//! the two support curves, the parameter on the bisector itself, the distance
//! and the point.

use occt_core::gp::GpPnt2d;

/// `class Bisector_PointOnBis` (`Bisector_PointOnBis.hxx:30-75`).
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct BisectorPointOnBis {
    /// `param1` (`.hxx:70`).
    param1: f64,
    /// `param2` (`.hxx:71`).
    param2: f64,
    /// `paramBis` (`.hxx:72`).
    param_bis: f64,
    /// `distance` (`.hxx:73`).
    distance: f64,
    /// `point` (`.hxx:68`).
    point: GpPnt2d,
    /// `infinite` (`.hxx:74`).
    infinite: bool,
}

impl BisectorPointOnBis {
    /// `Bisector_PointOnBis()` (`cxx:20-28`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `Bisector_PointOnBis(Param1, Param2, ParamBis, Distance, P)` (`cxx:32-48`).
    pub fn with_params(
        param1: f64,
        param2: f64,
        param_bis: f64,
        distance: f64,
        p: GpPnt2d,
    ) -> Self {
        Self {
            param1,
            param2,
            param_bis,
            distance,
            point: p,
            infinite: false,
        }
    }

    /// `ParamOnC1(Param)` (`cxx:52-56`).
    pub fn set_param_on_c1(&mut self, param: f64) {
        self.param1 = param;
    }

    /// `ParamOnC2(Param)` (`cxx:60-64`).
    pub fn set_param_on_c2(&mut self, param: f64) {
        self.param2 = param;
    }

    /// `ParamOnBis(Param)` (`cxx:68-72`).
    pub fn set_param_on_bis(&mut self, param: f64) {
        self.param_bis = param;
    }

    /// `Distance(Distance)` (`cxx:76-80`).
    pub fn set_distance(&mut self, distance: f64) {
        self.distance = distance;
    }

    /// `Point(P)` (`cxx:84-88`).
    pub fn set_point(&mut self, p: GpPnt2d) {
        self.point = p;
    }

    /// `IsInfinite(Infinite)` (`cxx:92-96`).
    pub fn set_is_infinite(&mut self, infinite: bool) {
        self.infinite = infinite;
    }

    /// `ParamOnC1()` (`cxx:100-104`).
    pub fn param_on_c1(&self) -> f64 {
        self.param1
    }

    /// `ParamOnC2()` (`cxx:108-112`).
    pub fn param_on_c2(&self) -> f64 {
        self.param2
    }

    /// `ParamOnBis()` (`cxx:116-120`).
    pub fn param_on_bis(&self) -> f64 {
        self.param_bis
    }

    /// `Distance()` (`cxx:124-128`).
    pub fn distance(&self) -> f64 {
        self.distance
    }

    /// `Point()` (`cxx:132-137`).
    pub fn point(&self) -> GpPnt2d {
        self.point
    }

    /// `IsInfinite()` (`Bisector_PointOnBis.hxx:63`).
    pub fn is_infinite(&self) -> bool {
        self.infinite
    }
}
