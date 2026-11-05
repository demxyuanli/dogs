//! `ChFi2d_FilletAlgo::Result` — the fillet arc and the two cut initial edges.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi2d/ChFi2d_FilletAlgo.cxx:592-759`
//! (OCCT 8.0.0).
//!
//! The body picks the stored parameter whose fillet centre is nearest to
//! `thePoint` (`cxx:606-648`), builds the circle arc between the two points of
//! tangency (`cxx:650-688`), then cuts each initial edge at its own point of
//! tangency (`cxx:690-757`).
//!
//! Translation notes:
//! - `ElSLib::PlaneParameters` / `PlaneValue` are
//!   [`occt_core::elib::slib::plane_parameters`] /
//!   [`occt_core::elib::slib::plane_value`],
//! - `ElCLib::Value(0.0, aCircle->Circ())` is
//!   [`occt_core::elib::clib::circle_value`], `Geom_Circle::D1` is
//!   [`occt_core::elib::clib::circle_d1`],
//! - `GeomAPI_ProjectPointOnCurve(P, aCircle)` is
//!   [`occt_geom::extrema_pc::point_curve_extrema_range`] over the circle's own
//!   `[0, 2*PI]` range; OCCT's `LowerDistanceParameter` throws on
//!   `NbPoints() == 0`, the port answers `None`,
//! - `BRepBuilderAPI_MakeEdge(C, u1, u2)` is
//!   [`brep_lib_make_edge::init_3d`] + [`crate::builder::TopoBuilder::make_edge`]
//!   + [`brep_lib_make_edge::add_endpoint_vertices`], the same three steps
//!   `BRepLib_MakeEdge::Init` takes.

use std::sync::Arc;

use occt_core::elib::{clib, slib};
use occt_core::gp::{GpAx2, GpCirc, GpPnt, GpPnt2d};
use occt_core::precision::CONFUSION;
use occt_geom::circle::GeomCircle;
use occt_geom::curve::Curve;

use crate::brep_lib_make_edge;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::Edge;

use super::fillet_point::FilletPoint;
use super::{plane_normal, FilletAlgo, FilletAlgoResult};

/// `TopoDS_Edge ChFi2d_FilletAlgo::Result(...)` (`cxx:592-759`).
pub(super) fn result(
    algo: &mut FilletAlgo,
    the_point: &GpPnt,
    i_solution: i32,
) -> Option<FilletAlgoResult> {
    // `double aX, aY; gp_Pnt2d aTargetPoint2d; ElSLib::PlaneParameters(...)`
    // (`cxx:602-604`). `aTargetPoint2d` is the only reader of `aX`/`aY`.
    let (a_x, a_y) = slib::plane_parameters(&algo.plane.pos, the_point);
    let a_target_point2d = GpPnt2d::new(a_x, a_y);

    // `choose the nearest circle` (`cxx:606-648`).
    let mut a_distance = 0.0;
    let mut a_nearest: Option<FilletPoint> = None;
    let mut i_sol = 1i32;
    for a in 1..=algo.result_params.len() {
        algo.start_side = algo.result_orientation[a - 1];
        let param = algo.result_params[a - 1];
        let mut a_point = FilletPoint::new(param);
        algo.fill_point(&mut a_point, param + 1.0);
        if a_point.has_solution(algo.radius).is_none() {
            continue;
        }
        let mut a_p = f64::MAX;
        if i_solution == -1 {
            a_p = (a_point.center().distance(&a_target_point2d) - algo.radius).abs();
        } else if i_solution == i_sol {
            a_p = 0.0;
        }
        if a_nearest.is_none() || a_p < a_distance {
            a_nearest = Some(a_point);
            a_distance = a_p;
        }
        if i_solution == i_sol {
            break;
        }
        i_sol += 1;
    }
    let a_nearest = a_nearest?;

    // `create circle edge` (`cxx:650-657`).
    let a_center = slib::plane_value(&algo.plane, a_nearest.center().x(), a_nearest.center().y());
    let a_circ = GpCirc::new(
        GpAx2::from_axis(a_center, plane_normal(&algo.plane)),
        algo.radius,
    );
    let a_circle: Arc<dyn Curve> = Arc::new(GeomCircle::new(a_circ.clone()));

    // `myCurve1->D0(aNearest->getParam(), aPoint2d1)` / same for the second
    // (`cxx:658-661`).
    let a_point2d1 = algo.curve1.as_ref()?.d0(a_nearest.param());
    let a_point2d2 = algo.curve2.as_ref()?.d0(a_nearest.param2());
    let a_point1 = slib::plane_value(&algo.plane, a_point2d1.x(), a_point2d1.y());
    let a_point2 = slib::plane_value(&algo.plane, a_point2d2.x(), a_point2d2.y());

    // `GeomAPI_ProjectPointOnCurve aProj(thePoint, aCircle)` (`cxx:663-672`);
    // the circle's own range is `[0, 2*PI]`.
    let (circle_first, circle_last) =
        (a_circle.first_parameter(), a_circle.last_parameter());
    let target = occt_geom::extrema_pc::point_curve_extrema_range(
        a_circle.as_ref(),
        the_point,
        circle_first,
        circle_last,
    )?;
    let mut a_target_param = target.u1;

    // `There is a bug in Open CASCADE in calculation of nearest point to a
    // circle near the parameter 0.0 Therefore I check this extrema point
    // manually` (`cxx:667-673`). OCCT also assigns `aPointOnCircle = p0`
    // there; nothing reads it afterwards, so the port keeps only the
    // parameter.
    let p0 = clib::circle_value(&a_circ, 0.0);
    if p0.distance(the_point) < target.p2.distance(the_point) {
        a_target_param = 0.0;
    }

    // `aProj.Perform(aPoint1); aParam1 = aProj.LowerDistanceParameter();` and
    // the same for `aPoint2` (`cxx:675-678`).
    let a_param1 = occt_geom::extrema_pc::point_curve_extrema_range(
        a_circle.as_ref(),
        &a_point1,
        circle_first,
        circle_last,
    )?
    .u1;
    let a_param2 = occt_geom::extrema_pc::point_curve_extrema_range(
        a_circle.as_ref(),
        &a_point2,
        circle_first,
        circle_last,
    )?
    .u1;

    // `cxx:679-688`: which of the two parameters is the arc's start.
    let mut is_out = (a_param1 < a_target_param && a_param2 < a_target_param)
        || (a_param1 > a_target_param && a_param2 > a_target_param);
    if a_param1 > a_param2 {
        is_out = !is_out;
    }
    let a_result = make_edge_range(
        &a_circle,
        if is_out { a_param2 } else { a_param1 },
        if is_out { a_param1 } else { a_param2 },
    )?;

    // `divide edges` (`cxx:690-757`).
    let curve_a = BRepTool::edge_curve(algo.edge1.as_ref()?)?;
    let (mut a_start, mut an_end) = BRepTool::edge_parameters(algo.edge1.as_ref()?);
    // `aCurve->D1(aNearest->getParam(), aPoint1, aDir); aCircle->D1(aParam1,
    // aPoint1, aCircleDir);` (`cxx:695-698`).
    let (_, a_dir) = curve_a.d1(a_nearest.param());
    let (_, a_circle_dir) = clib::circle_d1(&a_circ, a_param1);
    let cond = if a_circle_dir.angle(&a_dir) > std::f64::consts::FRAC_PI_2 {
        !is_out
    } else {
        is_out
    };
    if cond {
        a_start = a_nearest.param();
    } else {
        an_end = a_nearest.param();
    }
    // `Check the case when start and end are identical` (`cxx:708-714`).
    if (a_start - an_end).abs() < CONFUSION {
        an_end = a_start + CONFUSION;
    }
    // `BRepBuilderAPI_MakeEdge aDivider1(aCurve, aStart, anEnd);` (`cxx:715`).
    let a_divider1 = make_edge_range(&curve_a, a_start, an_end)?;

    let curve_b = BRepTool::edge_curve(algo.edge2.as_ref()?)?;
    let (mut a_start, mut an_end) = BRepTool::edge_parameters(algo.edge2.as_ref()?);
    let (_, a_dir) = curve_b.d1(a_nearest.param2());
    let (_, a_circle_dir) = clib::circle_d1(&a_circ, a_param2);
    // `(aCircleDir.Angle(aDir) > M_PI / 2.0) ? aIsOut : !aIsOut` (`cxx:732`).
    let cond = if a_circle_dir.angle(&a_dir) > std::f64::consts::FRAC_PI_2 {
        is_out
    } else {
        !is_out
    };
    if cond {
        a_start = a_nearest.param2();
    } else {
        an_end = a_nearest.param2();
    }
    if (a_start - an_end).abs() < CONFUSION {
        an_end = a_start + CONFUSION;
    }
    let a_divider2 = make_edge_range(&curve_b, a_start, an_end)?;

    // `if (myEdgesExchnged) { theEdge1 = aDivider2.Edge(); } else { theEdge2 =
    // aDivider2.Edge(); }` (`cxx:718-756`).
    let (edge1, edge2) = if algo.edges_exchanged {
        (a_divider2, a_divider1)
    } else {
        (a_divider1, a_divider2)
    };
    Some(FilletAlgoResult {
        fillet: a_result,
        edge1,
        edge2,
    })
}

/// `BRepBuilderAPI_MakeEdge(C, u1, u2)` -> `BRepLib_MakeEdge::Init(C, p1, p2)`:
/// the edge over `[first, last]` plus the two endpoint vertices.
fn make_edge_range(curve: &Arc<dyn Curve>, first: f64, last: f64) -> Option<Edge> {
    let build = brep_lib_make_edge::init_3d(curve, first, last)?;
    let builder = TopoBuilder::new();
    let mut edge = builder.make_edge(Arc::clone(curve), build.first, build.last);
    brep_lib_make_edge::add_endpoint_vertices(&builder, &mut edge, &build);
    Some(edge)
}
