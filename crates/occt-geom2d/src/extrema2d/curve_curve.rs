use super::prelude::*;
use super::*;
use crate::geom2d_int::Geom2dIntGInter;

/// `Extrema_ExtPC2d` solution turned into the port's `Extrema2d` pair carrier:
/// `u1`/`u2` both hold the curve parameter, `p1` the query point.
fn point_solution2d(p: &GpPnt2d, s: &ExtPc2dSolution) -> Extrema2d {
    Extrema2d { p1: *p, p2: s.point, distance: p.distance(&s.point), u1: s.u, u2: s.u }
}

/// All extremal distances between a point and a 2D curve, in the order
/// `Extrema_ExtPC2d` reports them.
///
/// Faithful route: [`ExtremaExtPC2d`] = `Extrema_GGExtPC<Adaptor2d_Curve2d, ...>`
/// (`Extrema_ExtPC2d.hxx:31-38`, body `Extrema_GGExtPC.hxx`). The previous body
/// was a self-invented `span/0.1` seed grid (`analytic_solvers::build_samples2d`)
/// with bracket sign changes plus Newton (`newton_point_curve_all2d`); both are
/// deleted. OCCT neither sorts nor dedups this list, and it only adds a range
/// endpoint when the query point already coincides with it
/// (`Extrema_GGExtPC.hxx:473-502`), hence the empty list on `!IsDone()`.
pub fn point_curve_extrema2d_all(c: &dyn Curve2d, p: &GpPnt2d) -> Vec<Extrema2d> {
    let ext = ExtremaExtPC2d::new(c, p);
    if !ext.is_done() {
        return Vec::new();
    }
    ext.solutions().iter().map(|s| point_solution2d(p, s)).collect()
}

/// Closest point of `p` on `c` as `Geom2dAPI_ProjectPointOnCurve::Perform`
/// reports it (`Geom2dAPI_ProjectPointOnCurve.cxx:59-86`): `Extrema_ExtPC2d`,
/// done only when `IsDone() && NbExt() > 0`, then the smallest
/// `SquareDistance(i)` with the first index winning ties. `None` is the
/// not-done/no-extremum case.
pub fn point_curve_extrema2d(c: &dyn Curve2d, p: &GpPnt2d) -> Option<Extrema2d> {
    let ext = ExtremaExtPC2d::new(c, p);
    if !ext.is_done() || ext.nb_ext() == 0 {
        return None;
    }
    let mut dist2_min = ext.square_distance(1);
    let mut index = 1;
    for i in 2..=ext.nb_ext() {
        let dist2 = ext.square_distance(i);
        if dist2 < dist2_min {
            dist2_min = dist2;
            index = i;
        }
    }
    Some(point_solution2d(p, &ext.solutions()[index - 1]))
}

/// Minimum distance between two curves (all extrema, first = closest).
///
/// Delegates to `curve_curve_extrema2d_all`. `samples` is kept for signature
/// compatibility but no longer drives the computation.
pub fn curve_curve_extrema2d(c1: &dyn Curve2d, c2: &dyn Curve2d, _samples: usize) -> Vec<Extrema2d> {
    curve_curve_extrema2d_all(c1, c2)
}

// ---------------------------------------------------------------------------
// Curve–curve extrema (port of `Extrema_ExtCC2d`; its analytic
// `Extrema_ExtElC2d` dispatch lives in `ext_el_c2d.rs`).
// ---------------------------------------------------------------------------

pub(super) fn pair2d_cc(p1: GpPnt2d, u1: f64, p2: GpPnt2d, u2: f64) -> Extrema2d {
    Extrema2d { p1, p2, distance: p1.distance(&p2), u1, u2 }
}

// The analytic `Extrema_ExtElC2d` arms (line/line, line/circle, circle/circle and
// the conic pairs) live in `extrema2d/ext_el_c2d.rs` as `ExtremaExtElC2d`, which
// `Extrema_ExtCC2d` dispatches to (`Extrema_ExtCC2d.cxx:104-145`). The
// module-level `line_line_extrema2d` / `line_circle_extrema2d` /
// `circle_circle_extrema2d` copies that used to live here were a second,
// divergent engine: the line/line copy fabricated a `myPoint` for the parallel
// arm (`Extrema_ExtElC2d.cxx:70-76` only sets `mySqDist`), and both circle arms
// appended zero-distance intersection pairs that `Extrema_ExtElC2d` never
// produces (`cxx:150-165` returns exactly the two extremal pairs; `cxx:354-367`
// returns exactly the four collinear pairs), on top of a local dedupe/sort.
// They had no caller other than `extrema2d/tests.rs`, so they are deleted rather
// than kept as a second engine.

/// Map `u` into `[lo, lo+period)` for a periodic curve (OCCT `ElCLib::InPeriod`).
pub(super) fn in_period2d(u: f64, lo: f64, period: f64) -> f64 {
    if period <= 0.0 || !lo.is_finite() {
        u
    } else {
        lo + (u - lo).rem_euclid(period)
    }
}

/// All extremal distances between two planar curves, in the order
/// `Extrema_ExtCC2d` reports them.
///
/// Faithful route: [`ExtremaExtCC2d`] = `Extrema_ExtCC2d` (`Extrema_ExtCC2d.hxx`,
/// `Extrema_ExtCC2d.cxx:104-145` for the `Extrema_ExtElC2d` dispatch and
/// `:435-451` for the general arm) with its engine `Extrema_ECC2d` =
/// `Extrema_GGenExtCC`. The previous body re-implemented that dispatch locally
/// and fell back to the self-invented seed grid plus Newton
/// (`newton_curve_curve_all2d`), which is deleted.
pub fn curve_curve_extrema2d_all(c1: &dyn Curve2d, c2: &dyn Curve2d) -> Vec<Extrema2d> {
    let ext = ExtremaExtCC2d::new(c1, c2);
    (1..=ext.nb_ext())
        .map(|n| {
            let (u1, p1, u2, p2) = ext.points(n);
            pair2d_cc(p1, u1, p2, u2)
        })
        .collect()
}

/// Local extremum of the distance between two planar curves from seed
/// parameters. Port of `Extrema_LocateExtCC2d`
/// (`Extrema_LocateExtCC2d.cxx:26-82`), which runs
/// `Extrema_LocECC2d` = `Extrema_GenLocateExtCC` -> `Extrema_CCLocFOfLocECC2d`
/// (`Extrema_GFuncExtCC` 2D) under `math_FunctionSetRoot`. The previous body was
/// a self-invented Newton loop (`refine_curve_curve2d_newton` /
/// `f_and_jac2d`), both deleted.
pub fn locate_extcc2d(c1: &dyn Curve2d, c2: &dyn Curve2d, u0: f64, v0: f64) -> Option<Extrema2d> {
    let loc = ExtremaLocateExtCC2d::new(c1, c2, u0, v0);
    if !loc.is_done() {
        return None;
    }
    let (p1, p2) = loc.point();
    Some(pair2d_cc(p1.value(), p1.parameter(), p2.value(), p2.parameter()))
}

/// All curve-curve intersection points, as `Geom2dAPI_InterCurveCurve` reports
/// them (`Geom2dAPI_InterCurveCurve.cxx:54-106`): the curves go to
/// `Geom2dInt_GInter(AC1, AC2, Tol, Tol)`, `NbPoints()` is 0 unless `IsDone()`,
/// and each point is `Point(i).Value()`.
///
/// The previous body was a self-invented 33x33 parameter grid that sampled both
/// curves and pushed every pair closer than `tol` with a local `1e-6` dedup; no
/// OCCT control flow does that, so it is deleted.
pub fn curve_curve_intersections2d(c1: &dyn Curve2d, c2: &dyn Curve2d, tol: f64) -> Vec<GpPnt2d> {
    let mut inter = Geom2dIntGInter::new();
    inter.perform_curves(c1, c2, tol, tol);
    if !inter.is_done() {
        return Vec::new();
    }
    (1..=inter.nb_points()).map(|n| *inter.point(n).value()).collect()
}

/// Distance from a point to the closest point on a polyline.
pub fn point_polyline_extrema2d(poly: &[GpPnt2d], p: &GpPnt2d) -> Extrema2d {
    let mut best: Option<Extrema2d> = None;
    for w in poly.windows(2) {
        let (a, b) = (w[0], w[1]);
        let abx = b.x() - a.x();
        let aby = b.y() - a.y();
        let len2 = abx * abx + aby * aby;
        let t = if len2 < 1e-30 {
            0.0
        } else {
            (((p.x() - a.x()) * abx + (p.y() - a.y()) * aby) / len2).clamp(0.0, 1.0)
        };
        let q = GpPnt2d::new(a.x() + t * abx, a.y() + t * aby);
        let d = p.distance(&q);
        if best.as_ref().map_or(true, |e: &Extrema2d| d < e.distance) {
            best = Some(Extrema2d { p1: *p, p2: q, distance: d, u1: t, u2: t });
        }
    }
    best.unwrap()
}

/// Tangent direction of a 2D curve at `u` (normalized).
pub fn extrema_tangent2d(c: &dyn Curve2d, u: f64) -> GpVec2d {
    let (_, d) = c.d1(u);
    let m = d.magnitude();
    if m > 1e-12 {
        GpVec2d::new(d.x() / m, d.y() / m)
    } else {
        GpVec2d::new(1.0, 0.0)
    }
}
