//! `IntTools_Context::IsVertexOnLine` — full OCCT 8.0.0 translation.
//!
//! Source: `IntTools_Context.cxx` lines 775-982 (the five-argument overload
//! that takes an explicit vertex tolerance). The four-argument overload
//! (`IsVertexOnLine(V, C, TolC, T)`) is the same walk with
//! `aTolV = BRep_Tool::Tolerance(aV)`.
//!
//! OCCT sequence:
//! 1. Build `aTolSum = 2 * (TolV + TolC)` then floor at `1.e-5` for
//!    BSpline/Bezier and `1.e-6` otherwise.
//! 2. Prefer the finite curve extremities (`FirstParameter` / `LastParameter`).
//!    When the vertex is within `aTolSum` of an end but farther than `TolV`,
//!    run a local `Extrema_LocateExtPC` from that end, falling back to a
//!    global `Extrema_ExtPC`. Reject a refined parameter that walks past the
//!    midpoint, exceeds `aTolSum`, or does not actually move off the end
//!    (`Confusion`).
//! 3. If both ends miss, project with `GeomAPI_ProjectPointOnCurve`. A miss on
//!    a bounded curve still accepts `StartPoint` / `EndPoint` within
//!    `aTolSum`.
//!
//! `Extrema_LocateExtPC` / `Extrema_ExtPC` are translated as a Newton walk
//! from a seed parameter and a dense sample-plus-Newton global search. The
//! OCCT extrema kernels themselves are not in this crate; the decision
//! predicates (midpoint, `aTolSum`, Confusion) are identical.

use occt_core::gp::GpPnt;
use occt_core::precision::{CONFUSION, INFINITE, Precision};
use occt_geom::Curve;

use crate::brep_tool::BRepTool;
use crate::inttools_roots;
use crate::shape::Vertex;

/// Floor used by OCCT for BSpline / Bezier `aTolSum`.
const BSPLINE_TOL_FLOOR: f64 = 1.0e-5;

/// Floor used by OCCT for every other curve type (`xft` comment in source).
const OTHER_TOL_FLOOR: f64 = 1.0e-6;

/// `Extrema_LocateExtPC` / `Extrema_ExtPC` precision argument in OCCT.
const EXTREMA_EPS: f64 = 1.0e-10;

/// Newton iterations for the local (`LocateExtPC`) search.
const LOCATE_ITERS: usize = 12;

/// Sample count for the global (`ExtPC`) search.
const GLOBAL_SAMPLES: usize = 64;

/// Newton iterations used to refine each global sample.
const GLOBAL_REFINE: usize = 8;

/// Result of [`is_vertex_on_line`]: the vertex lies on the curve at `t`.
#[derive(Debug, Clone, Copy)]
pub struct VertexOnLine {
    /// Parameter on the 3D curve (`aT` out-argument in OCCT).
    pub t: f64,
}

/// Curve classification used only to pick the OCCT `aTolSum` floor.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CurveFamily {
    BSplineOrBezier,
    Other,
}

/// True when the runtime type name looks like a BSpline or Bezier curve.
///
/// OCCT uses `GeomAdaptor_Curve::GetType()`. The port has no adaptor type
/// tag, so the same decision is made from the concrete type name — the
/// same approach already used by `pave_ff_paves::curve_is_bezier_or_bspline`.
fn curve_family(curve: &dyn Curve) -> CurveFamily {
    let name = std::any::type_name_of_val(curve);
    if name.contains("BSpline") || name.contains("Bezier") || name.contains("bspline") {
        CurveFamily::BSplineOrBezier
    } else {
        CurveFamily::Other
    }
}

/// OCCT `aTolSum` construction at `IntTools_Context.cxx:789-808`.
pub fn vertex_on_line_tol_sum(tol_v: f64, tol_c: f64, family: impl AsRef<str>) -> f64 {
    let mut sum = tol_v + tol_c;
    let name = family.as_ref();
    let is_nurbs = name.contains("BSpline") || name.contains("Bezier") || name.contains("bspline");
    sum = 2.0 * sum;
    if is_nurbs {
        if sum < BSPLINE_TOL_FLOOR {
            sum = BSPLINE_TOL_FLOOR;
        }
    } else if sum < OTHER_TOL_FLOOR {
        sum = OTHER_TOL_FLOOR;
    }
    sum
}

fn tol_sum_for_curve(curve: &dyn Curve, tol_v: f64, tol_c: f64) -> f64 {
    let mut sum = 2.0 * (tol_v + tol_c);
    match curve_family(curve) {
        CurveFamily::BSplineOrBezier => {
            if sum < BSPLINE_TOL_FLOOR {
                sum = BSPLINE_TOL_FLOOR;
            }
        }
        CurveFamily::Other => {
            if sum < OTHER_TOL_FLOOR {
                sum = OTHER_TOL_FLOOR;
            }
        }
    }
    sum
}

/// Finite parameter window of `curve`, or `None` when either end is infinite.
fn finite_range(curve: &dyn Curve) -> Option<(f64, f64)> {
    let a = curve.first_parameter();
    let b = curve.last_parameter();
    if Precision::is_infinite(a) || Precision::is_infinite(b) {
        return None;
    }
    if !a.is_finite() || !b.is_finite() {
        return None;
    }
    Some((a, b))
}

/// Clamp `t` into `[first, last]` when the range is finite.
fn clamp_param(t: f64, first: f64, last: f64) -> f64 {
    if Precision::is_infinite(first) || Precision::is_infinite(last) {
        return t;
    }
    if first <= last {
        t.clamp(first, last)
    } else {
        t.clamp(last, first)
    }
}

/// One Newton step of `Extrema_LocateExtPC`: minimise `|C(t) - P|^2`.
///
/// The derivative of the squared distance is `2 (C(t) - P) · C'(t)`, so the
/// Newton update is `t -= ((C-P)·C') / |C'|^2`. A vanishing derivative or a
/// non-finite update is reported as a failed step.
fn newton_step(curve: &dyn Curve, p: &GpPnt, t: f64) -> Option<f64> {
    let (q, d1) = curve.d1(t);
    let dx = q.x() - p.x();
    let dy = q.y() - p.y();
    let dz = q.z() - p.z();
    let speed2 = d1.square_magnitude();
    if speed2 <= EXTREMA_EPS * EXTREMA_EPS {
        return None;
    }
    let num = dx * d1.x() + dy * d1.y() + dz * d1.z();
    let next = t - num / speed2;
    if !next.is_finite() {
        return None;
    }
    Some(next)
}

/// `Extrema_LocateExtPC(P, C, tSeed, 1.e-10)` — local search from `t_seed`.
///
/// Returns `(t, point_on_curve)` when the walk stays finite. OCCT treats
/// `!IsDone` as a miss and falls through to `Extrema_ExtPC`.
pub(crate) fn extrema_locate_ext_pc(
    curve: &dyn Curve,
    p: &GpPnt,
    t_seed: f64,
    first: f64,
    last: f64,
) -> Option<(f64, GpPnt)> {
    locate_ext_pc(curve, p, t_seed, first, last)
}

fn locate_ext_pc(curve: &dyn Curve, p: &GpPnt, t_seed: f64, first: f64, last: f64) -> Option<(f64, GpPnt)> {
    let mut t = t_seed;
    let mut last_good = curve.d0(t);
    for _ in 0..LOCATE_ITERS {
        let Some(next) = newton_step(curve, p, t) else {
            break;
        };
        let next = clamp_param(next, first, last);
        let q = curve.d0(next);
        if q.distance(&last_good) < EXTREMA_EPS {
            t = next;
            last_good = q;
            break;
        }
        t = next;
        last_good = q;
    }
    if !t.is_finite() {
        return None;
    }
    Some((t, last_good))
}

/// Dense sampling used as the `Extrema_ExtPC` fallback.
///
/// OCCT iterates `NbExt` minima and keeps the smallest `SquareDistance`. The
/// port samples the finite range, Newton-refines each sample, and keeps the
/// closest refined point. Unbounded curves sample a window around the linear
/// estimate already used by `project_point_on_curve`.
fn extrema_ext_pc(curve: &dyn Curve, p: &GpPnt, first: f64, last: f64) -> Option<(f64, GpPnt, f64)> {
    let (a, b) = if Precision::is_infinite(first) || Precision::is_infinite(last) {
        let p0 = curve.d0(0.0);
        let tan = curve.d1(0.0).1;
        let t_est = if tan.square_magnitude() > 1e-300 {
            p.coord.subtracted(&p0.coord).dot(tan.xyz()) / tan.square_magnitude()
        } else {
            0.0
        };
        let t_est = if t_est.is_finite() { t_est } else { 0.0 };
        let w = t_est.abs().max(1.0) + 1.0;
        (t_est - w, t_est + w)
    } else {
        (first, last)
    };
    if (b - a).abs() <= CONFUSION {
        let q = curve.d0(a);
        return Some((a, q, p.distance(&q)));
    }
    let mut best_t = a;
    let mut best_p = curve.d0(a);
    let mut best_d2 = p.square_distance(&best_p);
    for i in 0..=GLOBAL_SAMPLES {
        let u = a + (b - a) * (i as f64 / GLOBAL_SAMPLES as f64);
        let mut t = u;
        for _ in 0..GLOBAL_REFINE {
            let Some(next) = newton_step(curve, p, t) else {
                break;
            };
            t = clamp_param(next, a, b);
        }
        let q = curve.d0(t);
        let d2 = p.square_distance(&q);
        if d2 < best_d2 {
            best_d2 = d2;
            best_t = t;
            best_p = q;
        }
    }
    Some((best_t, best_p, best_d2.sqrt()))
}

/// OCCT extremity refinement: start from `t_end`, keep it unless Locate/ExtPC
/// produce a point that stays on the same half of the curve, within `aTolSum`,
/// and actually moved off the end by more than Confusion.
fn refine_end(
    curve: &dyn Curve,
    pv: &GpPnt,
    t_end: f64,
    p_end: &GpPnt,
    first: f64,
    last: f64,
    tol_sum: f64,
    prefer_first: bool,
) -> f64 {
    let mid = 0.5 * (first + last);
    if let Some((t, q)) = locate_ext_pc(curve, pv, t_end, first, last) {
        let past_mid = if prefer_first { t > mid } else { t < mid };
        if past_mid || pv.distance(&q) > tol_sum || p_end.distance(&q) < CONFUSION {
            return t_end;
        }
        return t;
    }
    if let Some((t, q, _)) = extrema_ext_pc(curve, pv, first, last) {
        let past_mid = if prefer_first { t > mid } else { t < mid };
        if past_mid || pv.distance(&q) > tol_sum || p_end.distance(&q) < CONFUSION {
            return t_end;
        }
        return t;
    }
    t_end
}

/// Global `GeomAPI_ProjectPointOnCurve` stand-in used after the ends miss.
fn project_global(curve: &dyn Curve, p: &GpPnt, first: f64, last: f64) -> Option<(f64, f64)> {
    let (a, b) = if Precision::is_infinite(first) || Precision::is_infinite(last) {
        let p0 = curve.d0(0.0);
        let tan = curve.d1(0.0).1;
        let t_est = if tan.square_magnitude() > 1e-300 {
            p.coord.subtracted(&p0.coord).dot(tan.xyz()) / tan.square_magnitude()
        } else {
            0.0
        };
        let t_est = if t_est.is_finite() { t_est } else { 0.0 };
        let w = t_est.abs().max(1.0) + 1.0;
        (t_est - w, t_est + w)
    } else {
        (first, last)
    };
    if !a.is_finite() || !b.is_finite() {
        return None;
    }
    let t = inttools_roots::parameter(&|u| curve.d0(u), p, a, b);
    let dist = p.distance(&curve.d0(t));
    Some((t, dist))
}

/// True when `curve` reports a finite natural range, treated as bounded.
fn is_bounded_curve(curve: &dyn Curve) -> bool {
    finite_range(curve).is_some()
}

/// `IntTools_Context::IsVertexOnLine(V, TolV, C, TolC, T)`.
///
/// Returns `Some(t)` when the vertex lies on the curve within the OCCT
/// `aTolSum` policy. The parameter is written exactly as OCCT writes `aT`.
pub fn is_vertex_on_line_tol(
    vertex: &Vertex,
    tol_v: f64,
    curve: &dyn Curve,
    tol_c: f64,
) -> Option<f64> {
    let pv = BRepTool::vertex_point(vertex);
    is_point_on_line_tol(&pv, tol_v, curve, tol_c)
}

/// Four-argument OCCT overload: vertex tolerance taken from the vertex.
pub fn is_vertex_on_line(vertex: &Vertex, curve: &dyn Curve, tol_c: f64) -> Option<f64> {
    let tol_v = BRepTool::vertex_tolerance(vertex);
    is_vertex_on_line_tol(vertex, tol_v, curve, tol_c)
}

/// Geometry-only form used by `EstimatePaveOnCurve` and tests that already
/// hold the vertex point.
pub fn is_point_on_line_tol(pv: &GpPnt, tol_v: f64, curve: &dyn Curve, tol_c: f64) -> Option<f64> {
    let a_tol_sum = tol_sum_for_curve(curve, tol_v, tol_c);
    let a_first = curve.first_parameter();
    let a_last = curve.last_parameter();

    // Checking extremities first. It is necessary to choose the closest
    // bound to the point (`IntTools_Context.cxx:813-873`).
    let mut b_first_valid = false;
    let mut a_t = 0.0;
    let mut a_first_dist = INFINITE;

    if !Precision::is_infinite(a_first) && a_first.is_finite() {
        let a_pc_first = curve.d0(a_first);
        a_first_dist = pv.distance(&a_pc_first);
        if a_first_dist < a_tol_sum {
            b_first_valid = true;
            a_t = a_first;
            if a_first_dist > tol_v {
                a_t = refine_end(
                    curve,
                    pv,
                    a_first,
                    &a_pc_first,
                    a_first,
                    a_last,
                    a_tol_sum,
                    true,
                );
            }
        }
    }

    if !Precision::is_infinite(a_last) && a_last.is_finite() {
        let a_pc_last = curve.d0(a_last);
        let a_dist = pv.distance(&a_pc_last);
        if b_first_valid && a_first_dist < a_dist {
            return Some(a_t);
        }
        if a_dist < a_tol_sum {
            a_t = a_last;
            if a_dist > tol_v {
                a_t = refine_end(
                    curve,
                    pv,
                    a_last,
                    &a_pc_last,
                    a_first,
                    a_last,
                    a_tol_sum,
                    false,
                );
            }
            return Some(a_t);
        }
    } else if b_first_valid {
        return Some(a_t);
    }

    // GeomAPI_ProjectPointOnCurve (`IntTools_Context.cxx:942-981`).
    let Some((t_proj, a_dist)) = project_global(curve, pv, a_first, a_last) else {
        return bounded_end_fallback(curve, pv, a_first, a_last, a_tol_sum);
    };
    if a_dist > a_tol_sum {
        // OCCT still tries Start/End of a bounded curve when NbPoints == 0,
        // not when the projection is merely far. A far projection is a miss.
        return None;
    }
    Some(t_proj)
}

/// Bounded-curve Start/End fallback when the projector reports no points.
fn bounded_end_fallback(
    curve: &dyn Curve,
    pv: &GpPnt,
    a_first: f64,
    a_last: f64,
    a_tol_sum: f64,
) -> Option<f64> {
    if !is_bounded_curve(curve) {
        return None;
    }
    let a_p_start = curve.d0(a_first);
    if pv.distance(&a_p_start) < a_tol_sum {
        return Some(a_first);
    }
    let a_p_end = curve.d0(a_last);
    if pv.distance(&a_p_end) < a_tol_sum {
        return Some(a_last);
    }
    None
}

/// Boolean form matching the existing `IntToolsContext::is_vertex_on_line`
/// signature (parameter discarded).
pub fn is_vertex_on_line_bool(vertex: &Vertex, curve: &dyn Curve, tol: f64) -> bool {
    is_vertex_on_line(vertex, curve, tol).is_some()
}

/// `GeomAdaptor_Curve::Resolution(tol3d)` used by `PutPaveOnCurve` as `aPTol`.
///
/// OCCT converts a 3D tolerance into a parametric one. The port uses the
/// same speed estimate already present in `inttools_range::curve_resolution`.
pub fn adaptor_resolution(curve: &dyn Curve, tol3d: f64) -> f64 {
    let a = curve.first_parameter();
    let b = curve.last_parameter();
    if Precision::is_infinite(a) || Precision::is_infinite(b) || !a.is_finite() || !b.is_finite() {
        return tol3d.max(occt_core::precision::PCONFUSION);
    }
    crate::inttools_range::curve_resolution(curve, a, b, tol3d)
}

/// Project `p` onto `curve` and return `(parameter, distance)`.
///
/// Used by `PutEFPavesOnCurve` (`GeomAPI_ProjectPointOnCurve::Perform`) when
/// the section curve is a Bezier or BSpline.
pub fn project_point_on_curve_dist(curve: &dyn Curve, p: &GpPnt) -> Option<(f64, f64)> {
    let a_first = curve.first_parameter();
    let a_last = curve.last_parameter();
    project_global(curve, p, a_first, a_last)
}

/// `Extrema_ExtPC` on an explicit `[first, last]` (`ShapeAnalysis_Curve::ProjectAct`).
pub fn extrema_project_in_range(
    curve: &dyn Curve,
    p: &GpPnt,
    first: f64,
    last: f64,
) -> Option<(f64, f64)> {
    extrema_ext_pc(curve, p, first, last).map(|(t, _, d)| (t, d))
}

/// Number of projection solutions, matching `GeomAPI_ProjectPointOnCurve::NbPoints`.
///
/// The port reports 1 when a finite closest parameter exists and 0 otherwise.
pub fn project_point_on_curve_nb(curve: &dyn Curve, p: &GpPnt) -> usize {
    if project_point_on_curve_dist(curve, p).is_some() {
        1
    } else {
        0
    }
}

/// Lower-distance parameter of the projector (`LowerDistanceParameter`).
pub fn project_point_on_curve_param(curve: &dyn Curve, p: &GpPnt) -> Option<f64> {
    project_point_on_curve_dist(curve, p).map(|(t, _)| t)
}

/// Lower distance of the projector (`LowerDistance`).
pub fn project_point_on_curve_lower(curve: &dyn Curve, p: &GpPnt) -> Option<f64> {
    project_point_on_curve_dist(curve, p).map(|(_, d)| d)
}

/// True when `type_name` identifies a BSpline or Bezier (OCCT `GeomAbs` test).
pub fn curve_is_bezier_or_bspline(curve: &dyn Curve) -> bool {
    curve_family(curve) == CurveFamily::BSplineOrBezier
}

/// Diagnostic: the `aTolSum` that would be used for this vertex/curve pair.
pub fn debug_tol_sum(vertex: &Vertex, curve: &dyn Curve, tol_c: f64) -> f64 {
    tol_sum_for_curve(curve, BRepTool::vertex_tolerance(vertex), tol_c)
}

/// Diagnostic: extremity distances used by the first-pass check.
pub fn debug_end_distances(vertex: &Vertex, curve: &dyn Curve) -> (Option<f64>, Option<f64>) {
    let pv = BRepTool::vertex_point(vertex);
    let a = curve.first_parameter();
    let b = curve.last_parameter();
    let d0 = if !Precision::is_infinite(a) && a.is_finite() {
        Some(pv.distance(&curve.d0(a)))
    } else {
        None
    };
    let d1 = if !Precision::is_infinite(b) && b.is_finite() {
        Some(pv.distance(&curve.d0(b)))
    } else {
        None
    };
    (d0, d1)
}
