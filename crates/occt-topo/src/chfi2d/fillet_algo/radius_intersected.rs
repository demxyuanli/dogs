//! `IsRadiusIntersected` — the static validity test of `ChFi2d_FilletAlgo`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi2d/ChFi2d_FilletAlgo.cxx:162-254`
//! (OCCT 8.0.0).
//!
//! The body intersects the 2D curve with the infinite line through the fillet
//! radius (`cxx:179-180`) and rejects the candidate when an intersection point
//! is an endpoint of that radius (`cxx:191-204`) or when the two radius vectors
//! point in opposite directions (`cxx:204-207`). The `NbSegments` loop
//! (`cxx:208-252`) re-tests the same three conditions on
//! `theCurve->Value(FirstParameter())` and `theCurve->Value(LastParameter())`.
//!
//! Translation notes:
//! - `gp_Vec2d(aPoint, theStart)` is `theStart - aPoint`; `GpVec2d` has no
//!   two-point constructor, so both vectors are built from the difference,
//! - `anInter.Segment(a, aCurve)` is commented out in OCCT itself
//!   (`//not implemented (bug in OCC)`, `cxx:210`), so `aCurve` stays
//!   `theCurve`; the port reads the curve's own first and last points, which is
//!   what that commented-out call leaves behind,
//! - `Geom2dAPI_ProjectPointOnCurve aProjInt(...)` is re-`Init`ed before every
//!   read in the segment loop, so a fresh
//!   [`ProjectPointOnCurve2d`] per read is equivalent.

use std::sync::Arc;

use occt_core::gp::{GpDir2d, GpPnt2d, GpVec2d};
use occt_core::precision::{ANGULAR, CONFUSION, SQUARE_CONFUSION};
use occt_geom2d::curve::Curve2d;
use occt_geom2d::geom2d_api_inter_curve_curve::InterCurveCurve2d;
use occt_geom2d::geom2d_api_project_point_on_curve::ProjectPointOnCurve2d;
use occt_geom2d::line::Geom2dLine;

/// `IsRadiusIntersected(theCurve, theCurveMin, theCurveMax, theStart, theEnd,
/// theStartConnected)` (`cxx:166-254`).
pub(super) fn is_radius_intersected(
    the_curve: &Arc<dyn Curve2d>,
    the_curve_min: f64,
    the_curve_max: f64,
    the_start: &GpPnt2d,
    the_end: &GpPnt2d,
    the_start_connected: bool,
) -> bool {
    // `if (theStart.SquareDistance(theEnd) < Precision::SquareConfusion())`
    // (`cxx:176-179`).
    if the_start.square_distance(the_end) < SQUARE_CONFUSION {
        return false;
    }

    // `Geom2d_Line(theStart, gp_Dir2d(gp_Vec2d(theStart, theEnd)))`
    // (`cxx:179`): a zero-length direction raises `Standard_ConstructionError`,
    // excluded by the squared-distance test above but kept as a `false` here.
    let Ok(dir) = GpDir2d::new(the_end.x() - the_start.x(), the_end.y() - the_start.y()) else {
        return false;
    };
    let line: Arc<dyn Curve2d> = Arc::new(Geom2dLine::from_pnt_dir(*the_start, dir));
    let an_inter = InterCurveCurve2d::with_curves(the_curve, &line, CONFUSION);

    // `for (a = anInter.NbPoints(); a > 0; a--)` (`cxx:183-208`).
    let mut a = an_inter.nb_points();
    while a > 0 {
        let Some(a_point) = an_inter.point(a) else {
            a -= 1;
            continue;
        };
        a -= 1;
        if check_point(
            the_curve,
            the_curve_min,
            the_curve_max,
            the_start,
            the_end,
            the_start_connected,
            &a_point,
            false,
        ) {
            return true;
        }
    }

    // `for (a = anInter.NbSegments(); a > 0; a--)` (`cxx:209-252`).
    let mut a = an_inter.nb_segments();
    while a > 0 {
        // `// anInter.Segment(a, aCurve); //not implemented (bug in OCC)`
        // leaves `aCurve == theCurve` (`cxx:208`, `:210`).
        let first = the_curve.d0(the_curve.first_parameter());
        if check_point(
            the_curve,
            the_curve_min,
            the_curve_max,
            the_start,
            the_end,
            the_start_connected,
            &first,
            true,
        ) {
            return true;
        }
        let last = the_curve.d0(the_curve.last_parameter());
        if check_point(
            the_curve,
            the_curve_min,
            the_curve_max,
            the_start,
            the_end,
            the_start_connected,
            &last,
            true,
        ) {
            return true;
        }
        a -= 1;
    }
    false
}

/// The three tests of `cxx:186-206` / `cxx:214-251`.
///
/// `strict` selects the segment loop's acceptance test
/// (`LowerDistanceParameter() < Precision::Confusion()`, `cxx:216`, `:237`)
/// over the point loop's rejection test
/// (`LowerDistanceParameter() > Precision::Confusion()`, `cxx:188`).
fn check_point(
    the_curve: &Arc<dyn Curve2d>,
    the_curve_min: f64,
    the_curve_max: f64,
    the_start: &GpPnt2d,
    the_end: &GpPnt2d,
    the_start_connected: bool,
    a_point: &GpPnt2d,
    strict: bool,
) -> bool {
    let proj =
        ProjectPointOnCurve2d::with_curve_range(a_point, the_curve, the_curve_min, the_curve_max);
    let Some(parameter) = proj.lower_distance_parameter() else {
        return false; // `NbPoints() < 1`: the point is not on the edge.
    };
    // `aProjInt.LowerDistanceParameter() > Precision::Confusion()` (point
    // loop) / `< Precision::Confusion()` (segment loop, guarded by
    // `NbPoints()`).
    let on_edge = if strict {
        parameter < CONFUSION
    } else {
        parameter <= CONFUSION
    };
    if !on_edge {
        return false;
    }

    if a_point.distance(the_start) < CONFUSION && !the_start_connected {
        return true;
    }
    if a_point.distance(the_end) < CONFUSION {
        return true;
    }
    // `gp_Vec2d(aPoint, theStart).IsOpposite(gp_Vec2d(aPoint, theEnd),
    // Precision::Angular())`.
    let v_start = GpVec2d::new(the_start.x() - a_point.x(), the_start.y() - a_point.y());
    let v_end = GpVec2d::new(the_end.x() - a_point.x(), the_end.y() - a_point.y());
    v_start.is_opposite(&v_end, ANGULAR)
}
