//! `Extrema_ExtPC` public entry points.

use super::*;

/// `Extrema_ExtPC(P, C, Uinf, Usup)` — range-aware Perform.
///
/// Analytic Line/Circle stay in `point_curve_extrema_all`; this path covers
/// BSpline / Bezier / OtherCurve matching `GGExtPC::Perform`.
pub fn extrema_ext_pc_range(
    curve: &dyn Curve,
    point: &GpPnt,
    u_inf: f64,
    u_sup: f64,
) -> Vec<ExtPcSolution> {
    extrema_ext_pc_range_done(curve, point, u_inf, u_sup).0
}

/// `Extrema_ExtPC` (`Extrema_GGExtPC`) keeping its `IsDone` result.
///
/// The free-form `default:` arm sets `mydone = IntExtIsDone`
/// (`Extrema_GGExtPC.hxx:468`), so a **non-empty** solution list can coexist
/// with `IsDone() == false`: an interval boundary point is appended at
/// `hxx:452` (when `s1 * s2 < 0`) before any `IntervalPerform` succeeded.
/// `Extrema_GLocateExtPC`'s elementary arm tests `IsDone()` before scanning
/// the solutions (`Extrema_GLocateExtPC.hxx:283-296`), so callers that mirror
/// that arm need the flag.
pub fn extrema_ext_pc_range_done(
    curve: &dyn Curve,
    point: &GpPnt,
    u_inf: f64,
    u_sup: f64,
) -> (Vec<ExtPcSolution>, bool) {
    let (mut u_inf, mut u_sup) = (u_inf, u_sup);
    if u_inf > u_sup {
        std::mem::swap(&mut u_inf, &mut u_sup);
    }
    let mut state = ExtPcState::new(curve, u_inf, u_sup);
    if Precision::is_infinite(u_inf) {
        state.dist1 = f64::MAX;
    } else {
        state.pf = curve.d0(u_inf);
        state.dist1 = point.square_distance(&state.pf);
    }
    if Precision::is_infinite(u_sup) {
        state.dist2 = f64::MAX;
    } else {
        state.pl = curve.d0(u_sup);
        state.dist2 = point.square_distance(&state.pl);
    }

    // `Extrema_GGExtPC` curve-type switch: elementary curves go to
    // `Extrema_ExtPElC`, BSpline/Bezier/OtherCurve to the `default:` arm
    // (`Extrema_GGExtPC.hxx:390-502`).
    if let Some(pairs) = super::super::point_curve::ext_pelc_all(curve, point, u_inf, u_sup) {
        // `ExtPElC` filters to `[Uinf, Usup]` itself (`Extrema_ExtPElC.cxx:180`)
        // and sets `myIsMin` per solution (`cxx:77`, `:184`); `myDone = true`
        // regardless of the number of solutions (`cxx:189`).
        for (e, is_min) in pairs {
            state.add_sol(e.u1, e.p2, e.distance * e.distance, is_min);
        }
        state.done = true;
    } else if curve.bspline_knots().is_some() {
        perform_bspline(&mut state, curve, point);
        postprocess_ends(&mut state);
    } else {
        perform_general(&mut state, curve, point);
        postprocess_ends(&mut state);
    }
    let done = state.done;
    (state.sols, done)
}

/// Closest minimum among `IsMin` solutions (`ShapeAnalysis_Curve.cxx:280-301`).
pub fn extrema_ext_pc_min_in_range(
    curve: &dyn Curve,
    point: &GpPnt,
    u_inf: f64,
    u_sup: f64,
) -> Option<(f64, GpPnt, f64)> {
    let sols = extrema_ext_pc_range(curve, point, u_inf, u_sup);
    let mut best: Option<ExtPcSolution> = None;
    for s in sols {
        if !s.is_min {
            continue;
        }
        if best.map_or(true, |b| s.sq_dist < b.sq_dist) {
            best = Some(s);
        }
    }
    best.map(|s| (s.u, s.point, s.sq_dist.sqrt()))
}
