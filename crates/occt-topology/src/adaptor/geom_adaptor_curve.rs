use std::sync::Arc;

use occt_core::gp::{GpCirc, GpLin, GpPnt, GpTrsf, GpVec};
use occt_geom::curve::Curve;

/// `GeomAdaptor_Curve` — a `Geom_Curve` plus the range it is evaluated on.
///
/// Source: `TKG3d/GeomAdaptor/GeomAdaptor_Curve.cxx:262-271` (`Load(C, U1,
/// U2)`), `:243-245` (`myFirst` / `myLast`), `:576-585` (`IsClosed`).
///
/// `Load` keeps the *basis* parameters: `Value(U)` forwards to
/// `myCurve->Value(U)` and only `FirstParameter` / `LastParameter` report the
/// trimmed range. That is what the two `ChFi3d` callers below rely on —
/// `ChFi3d_ComputeArete` loads `GeomAdaptor_Curve(C3d, Pardeb, Parfin)` and
/// hands it to `ChFi3d_ComputePCurv`, whose p-curve is trimmed to the same
/// basis parameters.
#[derive(Clone)]
pub struct GeomAdaptorCurve {
    /// `myCurve`.
    curve: Option<Arc<dyn Curve>>,
    /// `myFirst`.
    first: f64,
    /// `myLast`.
    last: f64,
}

impl GeomAdaptorCurve {
    /// `GeomAdaptor_Curve()` (`cxx:28-33`).
    pub fn new() -> Self {
        Self {
            curve: None,
            first: 0.0,
            last: 0.0,
        }
    }

    /// `Load(C)` (`GeomAdaptor_Curve.cxx:225-231`): the curve's own range.
    pub fn of_curve(curve: Arc<dyn Curve>) -> Self {
        Self {
            first: curve.first_parameter(),
            last: curve.last_parameter(),
            curve: Some(curve),
        }
    }

    /// `Load(C, U1, U2)` (`cxx:262-271`).
    pub fn with_range(curve: Arc<dyn Curve>, first: f64, last: f64) -> Self {
        Self {
            curve: Some(curve),
            first,
            last,
        }
    }

    /// `Curve()` (`hxx`).
    pub fn curve(&self) -> Option<&Arc<dyn Curve>> {
        self.curve.as_ref()
    }

    /// `FirstParameter()` (`cxx:243-245`).
    pub fn first_parameter(&self) -> f64 {
        self.first
    }

    /// `LastParameter()` (`cxx:243-245`).
    pub fn last_parameter(&self) -> f64 {
        self.last
    }

    /// `Value(U)`.
    pub fn value(&self, u: f64) -> GpPnt {
        self.geom_curve().d0(u)
    }

    /// See [`CurveAdaptor::geom_curve`]: a null `myCurve` is
    /// `Standard_NullObject`.
    fn geom_curve(&self) -> &dyn Curve {
        self.curve
            .as_ref()
            .map(|c| c.as_ref())
            .expect("GeomAdaptor_Curve: no curve loaded (Standard_NullObject)")
    }
}

impl Default for GeomAdaptorCurve {
    fn default() -> Self {
        Self::new()
    }
}

/// `Adaptor3d_Curve` view of [`GeomAdaptorCurve`], the same way
/// [`CurveAdaptor`] is the view of `BRepAdaptor_Curve`.
impl Curve for GeomAdaptorCurve {
    fn d0(&self, u: f64) -> GpPnt {
        self.geom_curve().d0(u)
    }

    fn d1(&self, u: f64) -> (GpPnt, GpVec) {
        self.geom_curve().d1(u)
    }

    fn d2(&self, u: f64) -> (GpPnt, GpVec, GpVec) {
        self.geom_curve().d2(u)
    }

    fn d3(&self, u: f64) -> (GpPnt, GpVec, GpVec, GpVec) {
        self.geom_curve().d3(u)
    }

    fn first_parameter(&self) -> f64 {
        self.first
    }

    fn last_parameter(&self) -> f64 {
        self.last
    }

    fn continuity(&self) -> u8 {
        self.geom_curve().continuity()
    }

    fn nb_intervals(&self, continuity: u8) -> i32 {
        match self.curve.as_ref() {
            Some(c) => c.nb_intervals(continuity),
            None => 1,
        }
    }

    fn parameter_intervals(&self, continuity: u8) -> Vec<f64> {
        match self.curve.as_ref() {
            Some(c) => c.parameter_intervals(continuity),
            None => vec![self.first, self.last],
        }
    }

    fn resolution(&self, r3d: f64) -> f64 {
        match self.curve.as_ref() {
            Some(c) => c.resolution(r3d),
            None => r3d * occt_core::precision::PCONFUSION / occt_core::precision::CONFUSION,
        }
    }

    fn is_periodic(&self) -> bool {
        self.curve.as_ref().is_some_and(|c| c.is_periodic())
    }

    fn period(&self) -> f64 {
        self.curve.as_ref().map(|c| c.period()).unwrap_or(0.0)
    }

    fn circle_radius(&self) -> Option<f64> {
        self.curve.as_ref().and_then(|c| c.circle_radius())
    }

    fn is_line(&self) -> bool {
        self.curve.as_ref().is_some_and(|c| c.is_line())
    }

    fn gp_line(&self) -> Option<GpLin> {
        self.curve.as_ref().and_then(|c| c.gp_line())
    }

    fn gp_circ(&self) -> Option<GpCirc> {
        self.curve.as_ref().and_then(|c| c.gp_circ())
    }

    fn nurbs_degree(&self) -> Option<usize> {
        self.curve.as_ref().and_then(|c| c.nurbs_degree())
    }

    fn bspline_poles(&self) -> Option<&[GpPnt]> {
        self.curve.as_ref().and_then(|c| c.bspline_poles())
    }

    fn bspline_weights(&self) -> Option<&[f64]> {
        self.curve.as_ref().and_then(|c| c.bspline_weights())
    }

    fn bezier_poles(&self) -> Option<&[GpPnt]> {
        self.curve.as_ref().and_then(|c| c.bezier_poles())
    }

    /// `Adaptor3d_Curve` has no `Transform`; only `Geom_Curve` does.
    fn transform(&mut self, _t: &GpTrsf) {}

    /// `Adaptor3d_Curve` has no `Reverse`; only `Geom_Curve` does.
    fn reverse(&mut self) {}

    fn clone_dyn(&self) -> Box<dyn Curve> {
        Box::new(self.clone())
    }
}
