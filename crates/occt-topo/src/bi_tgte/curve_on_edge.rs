//! `BiTgte_CurveOnEdge` - the centre curve of a filler rolling on an edge.
//!
//! Source: `ModelingAlgorithms/TKOffset/BiTgte/BiTgte_CurveOnEdge.cxx:40-302`
//! (OCCT 8.0.0), `.hxx:45-141`.
//!
//! The adaptor keeps two trimmed curves:
//! - `myCurv` = `BRep_Tool::Curve(Edge)` trimmed to that edge's own range
//!   (`cxx:80-81`), i.e. the edge the filler rolls *on*;
//! - `myConF` = `BRep_Tool::Curve(EonF)` trimmed the same way (`cxx:84-85`),
//!   i.e. the edge the filler rolls *along*.
//!
//! `EvalD0` (`cxx:170-177`) evaluates `myConF` and projects that point onto
//! `myCurv` with `GeomAPI_ProjectPointOnCurve`; `FirstParameter`/
//! `LastParameter` are `myConF`'s (`cxx:106-117`). `Init` only classifies the
//! pair as a zero-radius `gp_Circ` when the rolling edge is a line coaxial
//! with a circular `EonF` (`cxx:91-101`).
//!
//! Translation notes:
//! - `GeomAdaptor_Curve` is [`curve_kind`] (the type chain of
//!   `GeomAdaptor_Curve.cxx:238-313`), `Curv.Line()` / `ConF.Circle()` are the
//!   `Curve::gp_line()` / `Curve::gp_circ()` queries.
//! - `myConF->Value(U)` is `Curve::d0`; `GeomAPI_ProjectPointOnCurve` is
//!   [`occt_geom::geom_api::project_point_on_curve`], whose `IsDone()` rule is
//!   `IsDone() && NbExt() > 0` (`GeomAPI_ProjectPointOnCurve.cxx:135-155`).
//!   OCCT's `NearestPoint()` on a not-done projector raises
//!   `StdFail_NotDone`; `None` from the port panics with that name.
//! - OCCT's `Continuity`/`NbIntervals`/`Intervals`/`Trim`/`IsClosed`/
//!   `IsPeriodic`/`Period`/`D1`/`D2`/`D3`/`DN`/`Resolution`/`Line`/`Ellipse`/
//!   `Hyperbola`/`Parabola`/`Degree`/`IsRational`/`NbPoles`/`NbKnots`/`Bezier`/
//!   `BSpline` all `throw Standard_NotImplemented` (`cxx:120-302`); the port
//!   keeps the same panic, and `Curve::gp_circ` answers `None` when `Circle()`
//!   would raise `Standard_NoSuchObject` (`cxx:230-237`).
//! - `Adaptor3d_Curve` has no `Transform`/`Reverse`, so the trait methods are
//!   empty as in [`crate::brep_adaptor_comp_curve`].

use std::sync::Arc;

use occt_core::gp::{GpCirc, GpLin, GpPnt, GpTrsf, GpVec};
use occt_geom::curve::Curve;
use occt_geom::geom_api::project_point_on_curve;
use occt_geom::trimmed::GeomTrimmedCurve;

use crate::brep_tool::BRepTool;
use crate::geom_abs_curve_type::{curve_kind, GeomAbsCurveType};
use crate::shape::Edge;

/// `BiTgte_CurveOnEdge` (`BiTgte_CurveOnEdge.hxx:45-141`).
#[derive(Clone)]
pub struct CurveOnEdge {
    /// `myEdge` (`hxx:135`).
    edge: Option<Edge>,
    /// `myEonF` (`hxx:136`).
    eon_f: Option<Edge>,
    /// `myCurv` (`hxx:137`) - `Geom_TrimmedCurve` over the edge curve.
    curv: Option<Arc<dyn Curve>>,
    /// `myConF` (`hxx:138`) - `Geom_TrimmedCurve` over the `EonF` curve.
    con_f: Option<Arc<dyn Curve>>,
    /// `myType` (`hxx:139`).
    curve_type: GeomAbsCurveType,
    /// `myCirc` (`hxx:140`).
    circ: GpCirc,
}

impl Default for CurveOnEdge {
    /// `BiTgte_CurveOnEdge()` (`cxx:42-45`): only `myType` is initialised.
    fn default() -> Self {
        Self {
            edge: None,
            eon_f: None,
            curv: None,
            con_f: None,
            curve_type: GeomAbsCurveType::OtherCurve,
            circ: GpCirc::new(occt_core::gp::GpAx2::default(), 0.0),
        }
    }
}

impl CurveOnEdge {
    /// `BiTgte_CurveOnEdge(const TopoDS_Edge&, const TopoDS_Edge&)`
    /// (`cxx:49-55`).
    pub fn new(eon_f: &Edge, edge: &Edge) -> Self {
        let mut a = Self {
            edge: Some(edge.clone()),
            eon_f: Some(eon_f.clone()),
            ..Self::default()
        };
        a.init(eon_f, edge);
        a
    }

    /// `Init(EonF, Edge)` (`cxx:76-102`).
    pub fn init(&mut self, eon_f: &Edge, edge: &Edge) {
        self.edge = Some(edge.clone());
        self.curv = Some(Self::trimmed_edge_curve(edge));

        self.eon_f = Some(eon_f.clone());
        self.con_f = Some(Self::trimmed_edge_curve(eon_f));

        self.curve_type = GeomAbsCurveType::OtherCurve;
        if let (Some(curv), Some(con_f)) = (self.curv.as_ref(), self.con_f.as_ref()) {
            if curve_kind(curv.as_ref()) == GeomAbsCurveType::Line
                && curve_kind(con_f.as_ref()) == GeomAbsCurveType::Circle
            {
                if let (Some(a1), Some(a2)) = (curv.gp_line(), con_f.gp_circ()) {
                    if a1
                        .position()
                        .is_coaxial_tol(&a2.axis(), occt_core::precision::ANGULAR, occt_core::precision::CONFUSION)
                    {
                        self.curve_type = GeomAbsCurveType::Circle;
                        self.circ = GpCirc::new(a2.position(), 0.0);
                    }
                }
            }
        }
    }

    /// `myCurv = new Geom_TrimmedCurve(BRep_Tool::Curve(myEdge, f, l), f, l)`
    /// (`cxx:80-81` and `cxx:84-85`). A null `BRep_Tool::Curve` makes the OCCT
    /// `Geom_TrimmedCurve` constructor raise `Standard_NullObject`.
    fn trimmed_edge_curve(edge: &Edge) -> Arc<dyn Curve> {
        let curve = BRepTool::edge_curve(edge).expect("Standard_NullObject: BRep_Tool::Curve");
        let (f, l) = BRepTool::edge_parameters(edge);
        Arc::new(GeomTrimmedCurve::new(curve, f, l))
    }

    /// `Adaptor3d_Curve::FirstParameter()` (`cxx:106-109`).
    pub fn first_parameter(&self) -> f64 {
        self.con_f
            .as_ref()
            .expect("Standard_NullObject: BiTgte_CurveOnEdge::myConF")
            .first_parameter()
    }

    /// `Adaptor3d_Curve::LastParameter()` (`cxx:113-116`).
    pub fn last_parameter(&self) -> f64 {
        self.con_f
            .as_ref()
            .expect("Standard_NullObject: BiTgte_CurveOnEdge::myConF")
            .last_parameter()
    }

    /// `Adaptor3d_Curve::GetType()` (`cxx:217-220`).
    pub fn get_type(&self) -> GeomAbsCurveType {
        self.curve_type
    }

    /// `Adaptor3d_Curve::Circle()` (`cxx:230-237`): `Standard_NoSuchObject`
    /// unless `GetType()` is `GeomAbs_Circle`.
    pub fn circle(&self) -> GpCirc {
        if self.curve_type != GeomAbsCurveType::Circle {
            panic!("Standard_NoSuchObject: BiTgte_CurveOnEdge::Circle");
        }
        self.circ.clone()
    }

    /// `Adaptor3d_Curve::D0(U, P)` (`cxx:170-177`): `myConF->Value(U)` projected
    /// onto `myCurv`.
    pub fn d0(&self, u: f64) -> GpPnt {
        let p = self
            .con_f
            .as_ref()
            .expect("Standard_NullObject: BiTgte_CurveOnEdge::myConF")
            .d0(u);
        let curv = self
            .curv
            .as_ref()
            .expect("Standard_NullObject: BiTgte_CurveOnEdge::myCurv");
        match project_point_on_curve(curv.as_ref(), &p, 0.0) {
            Some(proj) => proj.point,
            // `GeomAPI_ProjectPointOnCurve::NearestPoint()` on a failed
            // projector: `StdFail_NotDone`.
            None => panic!("StdFail_NotDone: BiTgte_CurveOnEdge::EvalD0"),
        }
    }
}

impl Curve for CurveOnEdge {
    fn d0(&self, u: f64) -> GpPnt {
        CurveOnEdge::d0(self, u)
    }

    /// `D1(U, P, D1)` (`cxx:181-184`): `Standard_NotImplemented`.
    fn d1(&self, _u: f64) -> (GpPnt, GpVec) {
        panic!("Standard_NotImplemented: BiTgte_CurveOnEdge");
    }

    /// `D2(U, P, D1, D2)` (`cxx:187-191`): `Standard_NotImplemented`.
    fn d2(&self, _u: f64) -> (GpPnt, GpVec, GpVec) {
        panic!("Standard_NotImplemented: BiTgte_CurveOnEdge");
    }

    /// `D3(U, P, D1, D2, D3)` (`cxx:194-198`): `Standard_NotImplemented`.
    fn d3(&self, _u: f64) -> (GpPnt, GpVec, GpVec, GpVec) {
        panic!("Standard_NotImplemented: BiTgte_CurveOnEdge");
    }

    /// `DN(U, N)` (`cxx:202-206`): `Standard_NotImplemented`.
    fn eval_dn(&self, _u: f64, _n: i32) -> GpVec {
        panic!("Standard_NotImplemented: BiTgte_CurveOnEdge");
    }

    fn first_parameter(&self) -> f64 {
        CurveOnEdge::first_parameter(self)
    }

    fn last_parameter(&self) -> f64 {
        CurveOnEdge::last_parameter(self)
    }

    /// `Continuity()` (`cxx:120-123`): `Standard_NotImplemented`.
    fn continuity(&self) -> u8 {
        panic!("Standard_NotImplemented: BiTgte_CurveOnEdge");
    }

    /// `NbIntervals(S)` (`cxx:127-130`): `Standard_NotImplemented`.
    fn nb_intervals(&self, _continuity: u8) -> i32 {
        panic!("Standard_NotImplemented: BiTgte_CurveOnEdge");
    }

    /// `Intervals(T, S)` (`cxx:134-137`): `Standard_NotImplemented`.
    fn parameter_intervals(&self, _continuity: u8) -> Vec<f64> {
        panic!("Standard_NotImplemented: BiTgte_CurveOnEdge");
    }

    /// `IsClosed()` (`cxx:150-153`): `Standard_NotImplemented`.
    fn is_closed(&self) -> bool {
        panic!("Standard_NotImplemented: BiTgte_CurveOnEdge");
    }

    /// `IsPeriodic()` (`cxx:157-160`): `Standard_NotImplemented`.
    fn is_periodic(&self) -> bool {
        panic!("Standard_NotImplemented: BiTgte_CurveOnEdge");
    }

    /// `Period()` (`cxx:164-167`): `Standard_NotImplemented`.
    fn period(&self) -> f64 {
        panic!("Standard_NotImplemented: BiTgte_CurveOnEdge");
    }

    /// `Resolution(R3d)` (`cxx:209-212`): `Standard_NotImplemented`.
    fn resolution(&self, _r3d: f64) -> f64 {
        panic!("Standard_NotImplemented: BiTgte_CurveOnEdge");
    }

    /// `Line()` (`cxx:225-228`): `Standard_NotImplemented`.
    fn gp_line(&self) -> Option<GpLin> {
        panic!("Standard_NotImplemented: BiTgte_CurveOnEdge");
    }

    /// `Circle()` (`cxx:230-237`): the stored zero-radius circle when the type
    /// is `GeomAbs_Circle`, `Standard_NoSuchObject` (here `None`) otherwise.
    fn gp_circ(&self) -> Option<GpCirc> {
        if self.curve_type == GeomAbsCurveType::Circle {
            Some(self.circ.clone())
        } else {
            None
        }
    }

    /// `Ellipse()` (`cxx:242-245`): `Standard_NotImplemented`.
    fn gp_ellipse(&self) -> Option<occt_core::gp::GpElips> {
        panic!("Standard_NotImplemented: BiTgte_CurveOnEdge");
    }

    /// `Hyperbola()` (`cxx:249-252`): `Standard_NotImplemented`.
    fn gp_hyperbola(&self) -> Option<occt_core::gp::GpHypr> {
        panic!("Standard_NotImplemented: BiTgte_CurveOnEdge");
    }

    /// `Parabola()` (`cxx:256-259`): `Standard_NotImplemented`.
    fn gp_parabola(&self) -> Option<occt_core::gp::GpParab> {
        panic!("Standard_NotImplemented: BiTgte_CurveOnEdge");
    }

    /// `Degree()` (`cxx:263-266`): `Standard_NotImplemented`.
    fn nurbs_degree(&self) -> Option<usize> {
        panic!("Standard_NotImplemented: BiTgte_CurveOnEdge");
    }

    /// `IsRational()` (`cxx:270-273`): `Standard_NotImplemented`.
    fn bspline_weights(&self) -> Option<&[f64]> {
        panic!("Standard_NotImplemented: BiTgte_CurveOnEdge");
    }

    /// `NbPoles()` (`cxx:277-280`): `Standard_NotImplemented`.
    fn bspline_poles(&self) -> Option<&[GpPnt]> {
        panic!("Standard_NotImplemented: BiTgte_CurveOnEdge");
    }

    /// `NbKnots()` (`cxx:283-286`): `Standard_NotImplemented`.
    fn bspline_knots(&self) -> Option<&[f64]> {
        panic!("Standard_NotImplemented: BiTgte_CurveOnEdge");
    }

    /// `Bezier()` (`cxx:290-293`): `Standard_NotImplemented`.
    fn bezier_poles(&self) -> Option<&[GpPnt]> {
        panic!("Standard_NotImplemented: BiTgte_CurveOnEdge");
    }

    /// `BSpline()` (`cxx:297-301`): `Standard_NotImplemented`.
    fn as_bspline_curve(&self) -> Option<&occt_geom::GeomBSplineCurve> {
        panic!("Standard_NotImplemented: BiTgte_CurveOnEdge");
    }

    /// `Adaptor3d_Curve` has no `Transform`.
    fn transform(&mut self, _t: &GpTrsf) {}

    /// `Adaptor3d_Curve` has no `Reverse`.
    fn reverse(&mut self) {}

    fn clone_dyn(&self) -> Box<dyn Curve> {
        Box::new(self.clone())
    }
}

