use std::sync::Arc;

use occt_core::gp::{GpCirc, GpElips, GpHypr, GpLin, GpParab, GpPnt, GpTrsf, GpVec};
use occt_geom::curve::Curve;

use crate::brep_tool::BRepTool;
use crate::shape::Edge;
use occt_geom::geom_abs_curve_type::{curve_kind, GeomAbsCurveType};

/// `BRepAdaptor_Curve(E)` (`BRepAdaptor_Curve.hxx:36-98`).
#[derive(Clone)]
pub struct CurveAdaptor {
    /// `GeomAdaptor_Curve::Curve()` after `Load(BRep_Tool::Curve(E, L))`.
    curve: Option<Arc<dyn Curve>>,
    /// `myFirst` (`GeomAdaptor_Curve.cxx:243`).
    first: f64,
    /// `myLast`.
    last: f64,
    /// `myTypeCurve` (`GeomAdaptor_Curve.cxx:238-313`).
    kind: GeomAbsCurveType,
    /// `BRepAdaptor_Curve::Tolerance()` = `BRep_Tool::Tolerance(E)`.
    tolerance: f64,
}

/// `RwLock`-free view of `BRepAdaptor_Curve`'s members: the inner `dyn Curve`
/// carries no `Debug`, so it prints as loaded / not loaded.
impl std::fmt::Debug for CurveAdaptor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("CurveAdaptor")
            .field("has_curve", &self.curve.is_some())
            .field("first", &self.first)
            .field("last", &self.last)
            .field("kind", &self.kind)
            .field("tolerance", &self.tolerance)
            .finish()
    }
}

impl CurveAdaptor {
    /// `BRepAdaptor_Curve(const TopoDS_Edge&)` — `Initialize(E)` with no face,
    /// so the 3D curve is used when the edge has one and the pcurve adaptor
    /// otherwise.
    pub fn of_edge(e: &Edge) -> Self {
        let curve = BRepTool::edge_curve_adaptor_world(e);
        let (first, last) = BRepTool::edge_parameters(e);
        let kind = match curve.as_ref() {
            Some(c) => curve_kind(c.as_ref()),
            None => GeomAbsCurveType::OtherCurve,
        };
        let tolerance = BRepTool::edge_tolerance(e);
        Self {
            curve,
            first,
            last,
            kind,
            tolerance,
        }
    }

    /// `Adaptor3d_Curve::Curve()`.
    pub fn curve(&self) -> Option<&Arc<dyn Curve>> {
        self.curve.as_ref()
    }

    /// `Adaptor3d_Curve::FirstParameter()`.
    pub fn first_parameter(&self) -> f64 {
        self.first
    }

    /// `Adaptor3d_Curve::LastParameter()`.
    pub fn last_parameter(&self) -> f64 {
        self.last
    }

    /// `Adaptor3d_Curve::GetType()`.
    pub fn get_type(&self) -> GeomAbsCurveType {
        self.kind
    }

    /// `Adaptor3d_Curve::Continuity()` — needed by
    /// `LProp_CurveUtils::IsTangentDefined` (`LProp_CurveUtils.hxx:332-378`).
    pub fn continuity(&self) -> u8 {
        match self.curve.as_ref() {
            Some(c) => c.continuity(),
            None => 0,
        }
    }

    /// `Adaptor3d_Curve::D0(U, P)`.
    pub fn d0(&self, u: f64) -> Option<GpPnt> {
        self.curve.as_ref().map(|c| c.d0(u))
    }

    /// `Adaptor3d_Curve::D1(U, P, D1)`.
    pub fn d1(&self, u: f64) -> Option<(GpPnt, GpVec)> {
        self.curve.as_ref().map(|c| c.d1(u))
    }

    /// `Adaptor3d_Curve::D2(U, P, D1, D2)`.
    pub fn d2(&self, u: f64) -> Option<(GpPnt, GpVec, GpVec)> {
        self.curve.as_ref().map(|c| c.d2(u))
    }

    /// `Adaptor3d_Curve::D3(U, P, D1, D2, D3)`.
    pub fn d3(&self, u: f64) -> Option<(GpPnt, GpVec, GpVec, GpVec)> {
        self.curve.as_ref().map(|c| c.d3(u))
    }

    /// `Adaptor3d_Curve::Resolution(R3d)` — `GeomAdaptor_Curve::Resolution`
    /// (`GeomAdaptor_Curve.cxx:1116-1148`), i.e. the underlying curve's.
    pub fn resolution(&self, r3d: f64) -> f64 {
        match self.curve.as_ref() {
            Some(c) => c.resolution(r3d),
            None => r3d * occt_core::precision::PCONFUSION / occt_core::precision::CONFUSION,
        }
    }

    /// `BRepAdaptor_Curve::Tolerance()` = `BRep_Tool::Tolerance(E)`.
    pub fn tolerance(&self) -> f64 {
        self.tolerance
    }

    /// `BRepAdaptor_Curve::Trim(First, Last, Tol)` (`BRepAdaptor_Curve.cxx:464-471`)
    /// clones the adaptor and calls
    /// `GeomAdaptor_Curve::Trim(First, Last, Tol)` (`GeomAdaptor_Curve.cxx:302-315`),
    /// which only overwrites `myFirst` / `myLast`: evaluation keeps the basis
    /// parameters, so the inner curve is untouched.
    pub fn trim(&self, first: f64, last: f64) -> CurveAdaptor {
        let mut trimmed = self.clone();
        trimmed.first = first;
        trimmed.last = last;
        trimmed
    }

    /// `Adaptor3d_Curve::IsClosed()` — `GeomAdaptor_Curve::IsClosed`
    /// (`GeomAdaptor_Curve.cxx:576-585`): an infinite range is never closed,
    /// otherwise the two ends must coincide within `Precision::Confusion()`.
    pub fn is_closed(&self) -> bool {
        let Some(curve) = self.curve.as_ref() else {
            return false;
        };
        if occt_core::precision::Precision::is_infinite(self.last)
            || occt_core::precision::Precision::is_infinite(self.first)
        {
            return false;
        }
        curve.d0(self.first).distance(&curve.d0(self.last)) <= occt_core::precision::CONFUSION
    }
}

/// `BRepAdaptor_Curve` *is* an `Adaptor3d_Curve` (`BRepAdaptor_Curve.hxx:36-98`),
/// so this is the trait view for algorithms that take one. `ChFiDS_Spine`
/// needs it for its `GCPnts_AbscissaPoint` calls on `myCurve`
/// (`ChFiDS_Spine.cxx:556-561`, `:609`, `:737`, `:786`, `:836`): those read the
/// *edge* range from the adaptor (`FirstParameter()` / `LastParameter()` store
/// `myFirst` / `myLast`, `GeomAdaptor_Curve.cxx:243-245`), not the inner
/// curve's own range, which is why the range is forwarded from `self.first` /
/// `self.last` here.
///
/// A null inner curve makes every geometric query raise `Standard_NullObject`
/// in OCCT; the trait cannot, so those unwrap with that message.
impl Curve for CurveAdaptor {
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
        CurveAdaptor::resolution(self, r3d)
    }

    fn is_periodic(&self) -> bool {
        match self.curve.as_ref() {
            Some(c) => c.is_periodic(),
            None => false,
        }
    }

    fn period(&self) -> f64 {
        match self.curve.as_ref() {
            Some(c) => c.period(),
            None => 0.0,
        }
    }

    fn circle_radius(&self) -> Option<f64> {
        self.curve.as_ref().and_then(|c| c.circle_radius())
    }

    fn gp_line(&self) -> Option<GpLin> {
        self.curve.as_ref().and_then(|c| c.gp_line())
    }

    fn gp_circ(&self) -> Option<GpCirc> {
        self.curve.as_ref().and_then(|c| c.gp_circ())
    }

    fn gp_ellipse(&self) -> Option<GpElips> {
        self.curve.as_ref().and_then(|c| c.gp_ellipse())
    }

    fn gp_hyperbola(&self) -> Option<GpHypr> {
        self.curve.as_ref().and_then(|c| c.gp_hyperbola())
    }

    fn gp_parabola(&self) -> Option<GpParab> {
        self.curve.as_ref().and_then(|c| c.gp_parabola())
    }

    fn is_line(&self) -> bool {
        self.curve.as_ref().is_some_and(|c| c.is_line())
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

    fn bezier_weights(&self) -> Option<&[f64]> {
        self.curve.as_ref().and_then(|c| c.bezier_weights())
    }

    /// `Adaptor3d_Curve` has no `Transform`; only `Geom_Curve` does.
    fn transform(&mut self, _t: &GpTrsf) {}

    /// `Adaptor3d_Curve` has no `Reverse`; only `Geom_Curve` does.
    fn reverse(&mut self) {}

    fn clone_dyn(&self) -> Box<dyn Curve> {
        Box::new(self.clone())
    }
}

impl CurveAdaptor {
    /// The inner curve, or the `Standard_NullObject` OCCT raises when the
    /// adaptor was never initialized.
    fn geom_curve(&self) -> &dyn Curve {
        self.curve
            .as_ref()
            .map(|c| c.as_ref())
            .expect("BRepAdaptor_Curve: no curve loaded (Standard_NullObject)")
    }
}
