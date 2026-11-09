//! `BiTgte_CurveOnVertex` - the constant curve of a filler rolling on a vertex.
//!
//! Source: `ModelingAlgorithms/TKOffset/BiTgte/BiTgte_CurveOnVertex.cxx:40-253`
//! (OCCT 8.0.0), `.hxx:45-135`.
//!
//! The adaptor only stores `BRep_Tool::Range(EonF)` and the vertex point
//! (`cxx:58-61`): `FirstParameter`/`LastParameter` answer the `EonF` range and
//! `EvalD0` answers the vertex point for every parameter (`cxx:132-135`).
//! `GetType()` is always `GeomAbs_OtherCurve` and `Circle()` always raises
//! `Standard_NoSuchObject` (`cxx:173-190`).
//!
//! Translation notes:
//! - `BRep_Tool::Range(EonF, f, l)` is [`crate::brep_tool::BRepTool::edge_parameters`],
//!   `BRep_Tool::Pnt(V)` is [`crate::brep_tool::BRepTool::vertex_point`].
//! - Every other `Adaptor3d_Curve` method `throw Standard_NotImplemented`
//!   (`cxx:80-253`); the port keeps the same panic.
//! - `Adaptor3d_Curve` has no `Transform`/`Reverse`, so those trait methods are
//!   empty as in [`crate::brep_adaptor_comp_curve`].

use occt_core::gp::{GpCirc, GpLin, GpPnt, GpTrsf, GpVec};
use occt_geom::curve::Curve;

use crate::brep_tool::BRepTool;
use crate::geom_abs_curve_type::GeomAbsCurveType;
use crate::shape::{Edge, Vertex};

/// `BiTgte_CurveOnVertex` (`BiTgte_CurveOnVertex.hxx:45-135`).
#[derive(Clone)]
pub struct CurveOnVertex {
    /// `myFirst` (`hxx:133`).
    first: f64,
    /// `myLast` (`hxx:134`).
    last: f64,
    /// `myPnt` (`hxx:135`).
    pnt: GpPnt,
}

impl Default for CurveOnVertex {
    /// `BiTgte_CurveOnVertex()` (`cxx:42-46`): `myFirst = myLast = 0`.
    fn default() -> Self {
        Self {
            first: 0.0,
            last: 0.0,
            pnt: GpPnt::default(),
        }
    }
}

impl CurveOnVertex {
    /// `BiTgte_CurveOnVertex(const TopoDS_Edge&, const TopoDS_Vertex&)`
    /// (`cxx:49-54`).
    pub fn new(eon_f: &Edge, v: &Vertex) -> Self {
        let mut a = Self::default();
        a.init(eon_f, v);
        a
    }

    /// `Init(EonF, V)` (`cxx:58-62`).
    pub fn init(&mut self, eon_f: &Edge, v: &Vertex) {
        let (f, l) = BRepTool::edge_parameters(eon_f);
        self.first = f;
        self.last = l;
        self.pnt = BRepTool::vertex_point(v);
    }

    /// `Adaptor3d_Curve::FirstParameter()` (`cxx:65-68`).
    pub fn first_parameter(&self) -> f64 {
        self.first
    }

    /// `Adaptor3d_Curve::LastParameter()` (`cxx:72-75`).
    pub fn last_parameter(&self) -> f64 {
        self.last
    }

    /// `Adaptor3d_Curve::GetType()` (`cxx:175-178`): always
    /// `GeomAbs_OtherCurve`.
    pub fn get_type(&self) -> GeomAbsCurveType {
        GeomAbsCurveType::OtherCurve
    }

    /// `Adaptor3d_Curve::D0(U, P)` (`cxx:132-135`): the vertex point.
    pub fn d0(&self, _u: f64) -> GpPnt {
        self.pnt
    }
}

impl Curve for CurveOnVertex {
    fn d0(&self, u: f64) -> GpPnt {
        CurveOnVertex::d0(self, u)
    }

    /// `D1(U, P, D1)` (`cxx:139-142`): `Standard_NotImplemented`.
    fn d1(&self, _u: f64) -> (GpPnt, GpVec) {
        panic!("Standard_NotImplemented: BiTgte_CurveOnVertex");
    }

    /// `D2(U, P, D1, D2)` (`cxx:145-149`): `Standard_NotImplemented`.
    fn d2(&self, _u: f64) -> (GpPnt, GpVec, GpVec) {
        panic!("Standard_NotImplemented: BiTgte_CurveOnVertex");
    }

    /// `D3(U, P, D1, D2, D3)` (`cxx:152-156`): `Standard_NotImplemented`.
    fn d3(&self, _u: f64) -> (GpPnt, GpVec, GpVec, GpVec) {
        panic!("Standard_NotImplemented: BiTgte_CurveOnVertex");
    }

    /// `DN(U, N)` (`cxx:160-163`): `Standard_NotImplemented`.
    fn eval_dn(&self, _u: f64, _n: i32) -> GpVec {
        panic!("Standard_NotImplemented: BiTgte_CurveOnVertex");
    }

    fn first_parameter(&self) -> f64 {
        CurveOnVertex::first_parameter(self)
    }

    fn last_parameter(&self) -> f64 {
        CurveOnVertex::last_parameter(self)
    }

    /// `Continuity()` (`cxx:80-83`): `Standard_NotImplemented`.
    fn continuity(&self) -> u8 {
        panic!("Standard_NotImplemented: BiTgte_CurveOnVertex");
    }

    /// `NbIntervals(S)` (`cxx:87-90`): `Standard_NotImplemented`.
    fn nb_intervals(&self, _continuity: u8) -> i32 {
        panic!("Standard_NotImplemented: BiTgte_CurveOnVertex");
    }

    /// `Intervals(T, S)` (`cxx:94-97`): `Standard_NotImplemented`.
    fn parameter_intervals(&self, _continuity: u8) -> Vec<f64> {
        panic!("Standard_NotImplemented: BiTgte_CurveOnVertex");
    }

    /// `IsClosed()` (`cxx:110-113`): `Standard_NotImplemented`.
    fn is_closed(&self) -> bool {
        panic!("Standard_NotImplemented: BiTgte_CurveOnVertex");
    }

    /// `IsPeriodic()` (`cxx:117-120`): `Standard_NotImplemented`.
    fn is_periodic(&self) -> bool {
        panic!("Standard_NotImplemented: BiTgte_CurveOnVertex");
    }

    /// `Period()` (`cxx:124-127`): `Standard_NotImplemented`.
    fn period(&self) -> f64 {
        panic!("Standard_NotImplemented: BiTgte_CurveOnVertex");
    }

    /// `Resolution(R3d)` (`cxx:166-169`): `Standard_NotImplemented`.
    fn resolution(&self, _r3d: f64) -> f64 {
        panic!("Standard_NotImplemented: BiTgte_CurveOnVertex");
    }

    /// `Line()` (`cxx:180-183`): `Standard_NotImplemented`.
    fn gp_line(&self) -> Option<GpLin> {
        panic!("Standard_NotImplemented: BiTgte_CurveOnVertex");
    }

    /// `Circle()` (`cxx:187-190`): always `Standard_NoSuchObject`.
    fn gp_circ(&self) -> Option<GpCirc> {
        None
    }

    /// `Ellipse()` (`cxx:194-197`): `Standard_NotImplemented`.
    fn gp_ellipse(&self) -> Option<occt_core::gp::GpElips> {
        panic!("Standard_NotImplemented: BiTgte_CurveOnVertex");
    }

    /// `Hyperbola()` (`cxx:201-204`): `Standard_NotImplemented`.
    fn gp_hyperbola(&self) -> Option<occt_core::gp::GpHypr> {
        panic!("Standard_NotImplemented: BiTgte_CurveOnVertex");
    }

    /// `Parabola()` (`cxx:208-211`): `Standard_NotImplemented`.
    fn gp_parabola(&self) -> Option<occt_core::gp::GpParab> {
        panic!("Standard_NotImplemented: BiTgte_CurveOnVertex");
    }

    /// `Degree()` (`cxx:215-218`): `Standard_NotImplemented`.
    fn nurbs_degree(&self) -> Option<usize> {
        panic!("Standard_NotImplemented: BiTgte_CurveOnVertex");
    }

    /// `IsRational()` (`cxx:222-225`): `Standard_NotImplemented`.
    fn bspline_weights(&self) -> Option<&[f64]> {
        panic!("Standard_NotImplemented: BiTgte_CurveOnVertex");
    }

    /// `NbPoles()` (`cxx:229-232`): `Standard_NotImplemented`.
    fn bspline_poles(&self) -> Option<&[GpPnt]> {
        panic!("Standard_NotImplemented: BiTgte_CurveOnVertex");
    }

    /// `NbKnots()` (`cxx:236-239`): `Standard_NotImplemented`.
    fn bspline_knots(&self) -> Option<&[f64]> {
        panic!("Standard_NotImplemented: BiTgte_CurveOnVertex");
    }

    /// `Bezier()` (`cxx:243-246`): `Standard_NotImplemented`.
    fn bezier_poles(&self) -> Option<&[GpPnt]> {
        panic!("Standard_NotImplemented: BiTgte_CurveOnVertex");
    }

    /// `BSpline()` (`cxx:250-253`): `Standard_NotImplemented`.
    fn as_bspline_curve(&self) -> Option<&occt_geom::GeomBSplineCurve> {
        panic!("Standard_NotImplemented: BiTgte_CurveOnVertex");
    }

    /// `Adaptor3d_Curve` has no `Transform`.
    fn transform(&mut self, _t: &GpTrsf) {}

    /// `Adaptor3d_Curve` has no `Reverse`.
    fn reverse(&mut self) {}

    fn clone_dyn(&self) -> Box<dyn Curve> {
        Box::new(self.clone())
    }
}
