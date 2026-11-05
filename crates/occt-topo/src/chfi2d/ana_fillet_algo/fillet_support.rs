//! Support routines of `ChFi2d_AnaFilletAlgo.cxx` (OCCT 8.0.0).
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | `isCW(AC)` (`cxx:39-69`) | 31 | [`is_cw`] |
//! | `IsEqual(p1, p2)` 3D (`cxx:72-75`) | 4 | [`is_equal`] |
//! | `IsEqual(p1, p2)` 2D (`cxx:77-80`) | 4 | [`is_equal_2d`] |
//! | `Cut(plane, e1, e2)` (`cxx:1050-1101`) | 52 | [`cut`] |
//!
//! Translation notes:
//! - `BRepAdaptor_Curve` (`cxx:41-47`) is
//!   [`crate::brep_lib_make_edge::adaptor_curve`] over the stored curve:
//!   trim unwrapping plus the interval in basis parameters, which is the
//!   parameter space `AC.Value` / `AC.FirstParameter` / `AC.LastParameter`
//!   speak,
//! - `occ::down_cast<Geom_Circle>(AC.Curve().Curve())` (`cxx:44`) answers a
//!   null handle for every other curve type and OCCT then dereferences it
//!   (`cxx:51`); the port reads [`Curve::gp_circ`] and falls back to `false`
//!   on that null arm,
//! - `gp_Ax3 plane = AC.Circle().Position();` (`cxx:48`) is
//!   `GpAx2::to_ax3` over the circle's `gp_Ax2`,
//! - `BRepBuilderAPI_MakeEdge(C, U1, U2)` (`cxx:1081`, `:1087`) is
//!   [`crate::brep_lib_make_edge::init_3d`] plus
//!   [`crate::brep_lib_make_edge::add_endpoint_vertices`]; its `!IsDone()`
//!   arms make `Edge()` throw in OCCT, so they are the `None` here.

use std::f64::consts::PI;
use std::sync::Arc;

use occt_core::gp::{GpPln, GpPnt, GpPnt2d, GpVec};
use occt_core::precision::{CONFUSION, SQUARE_CONFUSION};
use occt_geom::curve::Curve;
use occt_geom::geom_api::ExtremaCurveCurve;
use occt_geom::projlib;

use crate::brep_lib_make_edge;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::Edge;

/// `static bool isCW(const BRepAdaptor_Curve& AC)` (`cxx:40-69`): `true` when
/// the arc runs clockwise in its own plane, decided by comparing the angle to
/// the mid-parameter point against the angle to the end point.
pub(super) fn is_cw(curve: &Arc<dyn Curve>, first: f64, last: f64) -> bool {
    // `occ::down_cast<Geom_Circle>(AC.Curve().Curve())` (`cxx:44`).
    let Some(circle) = curve.gp_circ() else {
        return false;
    };
    // `gp_Pnt start = AC.Value(f); gp_Pnt end = AC.Value(l);` (`cxx:45-46`).
    let start = curve.d0(first);
    let end = curve.d0(last);
    let center = circle.location();
    // `gp_Ax3 plane = AC.Circle().Position();` (`cxx:48`).
    let plane = circle.position().to_ax3();

    // `circle->D0(0.5 * (f + l), m);` (`cxx:51`).
    let middle = curve.d0(0.5 * (first + last));

    let startv = GpVec::from_pnts(&center, &start);
    let endv = GpVec::from_pnts(&center, &end);
    let middlev = GpVec::from_pnts(&center, &middle);
    let normal = GpVec::from_xyz(plane.direction().xyz());

    // `cxx:55-60`.
    let mut middlea = startv.angle_with_ref(&middlev, &normal);
    while middlea < 0.0 {
        middlea += 2.0 * PI;
    }
    // `cxx:61-66`.
    let mut enda = startv.angle_with_ref(&endv, &normal);
    while enda < 0.0 {
        enda += 2.0 * PI;
    }
    // `cxx:68`.
    middlea > enda
}

/// `static bool IsEqual(const gp_Pnt&, const gp_Pnt&)` (`cxx:72-75`).
pub(super) fn is_equal(p1: &GpPnt, p2: &GpPnt) -> bool {
    p1.square_distance(p2) < SQUARE_CONFUSION
}

/// `static bool IsEqual(const gp_Pnt2d&, const gp_Pnt2d&)` (`cxx:77-80`).
pub(super) fn is_equal_2d(p1: &GpPnt2d, p2: &GpPnt2d) -> bool {
    p1.square_distance(p2) < SQUARE_CONFUSION
}

/// `BRep_Tool::Curve(E, f, l)` in the `BRepAdaptor_Curve` parameter space: the
/// untrimmed basis and the interval expressed in that basis's parameters.
pub(super) fn curve_range(edge: &Edge) -> Option<(Arc<dyn Curve>, f64, f64)> {
    let stored = BRepTool::edge_curve(edge)?;
    Some(brep_lib_make_edge::adaptor_curve(&stored))
}

/// `BRepBuilderAPI_MakeEdge(C, U1, U2)`: `BRepLib_MakeEdge::Init` with null
/// vertices (`BRepLib_MakeEdge.cxx:740-775`). `None` is the `!IsDone()` of
/// `ParameterOutOfRange` / `LineThroughIdenticPoints`
/// (`BRepLib_MakeEdge.cxx:653-663`), which makes `BRepBuilderAPI_MakeEdge::Edge()`
/// throw in OCCT.
pub(super) fn make_edge_over(curve: &Arc<dyn Curve>, first: f64, last: f64) -> Option<Edge> {
    let build = brep_lib_make_edge::init_3d(curve, first, last)?;
    let builder = TopoBuilder::new();
    let mut edge = builder.make_edge(curve.clone(), build.first, build.last);
    brep_lib_make_edge::add_endpoint_vertices(&builder, &mut edge, &build);
    Some(edge)
}

/// `BRepBuilderAPI_MakeEdge::Init(C, P1, P2)` (`BRepLib_MakeEdge.cxx:522-541`),
/// the form the two shrunk neighbours are built with
/// (`ChFi2d_AnaFilletAlgo.cxx:519`, `:550`). `None` is
/// `BRepLib_PointProjectionFailed`.
pub(super) fn make_edge_points(
    curve: &Arc<dyn Curve>,
    p1: &GpPnt,
    p2: &GpPnt,
) -> Option<Edge> {
    let build = brep_lib_make_edge::init_curve_points(curve, p1, p2)?;
    let builder = TopoBuilder::new();
    let mut edge = builder.make_edge(curve.clone(), build.first, build.last);
    brep_lib_make_edge::add_endpoint_vertices(&builder, &mut edge, &build);
    Some(edge)
}

/// `bool ChFi2d_AnaFilletAlgo::Cut(const gp_Pln&, TopoDS_Edge&, TopoDS_Edge&)`
/// (`cxx:1051-1101`): trims the two edges at their first intersection, and
/// answers the projected intersection point the caller writes into
/// `(x12, y12)` / `(x21, y21)` (`cxx:1093-1097`).
///
/// `None` is the C++ `false` of both early exits (`cxx:1080`, `cxx:1100`);
/// `theE1` / `theE2` keep their values there, exactly like the C++ code.
pub(super) fn cut(plane: &GpPln, the_e1: &mut Edge, the_e2: &mut Edge) -> Option<GpPnt2d> {
    // `cxx:1057-1058`: `BRep_Tool::Curve(theE1, f1, l1)`.
    let (c1, f1, l1) = curve_range(the_e1)?;
    let (c2, f2, l2) = curve_range(the_e2)?;

    // `GeomAPI_ExtremaCurveCurve extrema(c1, c2, f1, l1, f2, l2);` (`cxx:1059`).
    let extrema = ExtremaCurveCurve::with_ranges(c1.as_ref(), c2.as_ref(), f1, l1, f2, l2);
    let mut found = false;
    let mut point = GpPnt::new(0.0, 0.0, 0.0);
    let mut param1 = 0.0;
    let mut param2 = 0.0;
    // `cxx:1060-1077`.
    let nb = extrema.nb_extrema();
    if nb > 0 {
        for i in 1..=nb {
            let d = extrema.distance(i);
            if d < CONFUSION {
                let (u1, u2) = extrema.parameters(i);
                if (l1 - u1).abs() > CONFUSION && (f2 - u2).abs() > CONFUSION {
                    found = true;
                    param1 = u1;
                    param2 = u2;
                    point = extrema.points(i).0;
                    break;
                }
            }
        }
    }
    if !found {
        return None;
    }

    // `cxx:1081-1091`: `BRepBuilderAPI_MakeEdge(c1, f1, param1)` into `theE1`,
    // then `BRepBuilderAPI_MakeEdge(c2, param2, l2)` into `theE2`.
    let new_e1 = make_edge_over(&c1, f1, param1)?;
    *the_e1 = new_e1;
    let new_e2 = make_edge_over(&c2, param2, l2)?;
    *the_e2 = new_e2;

    // `cxx:1093-1096`: `ProjLib::Project(thePlane, p)`.
    Some(projlib::project_pln_pnt(plane, &point))
}
