//! `Perform` on a general curve.

use super::*;

/// Default / OtherCurve / Offset path (`GGExtPC.hxx:390-469`).
/// Bezier early-return is `hxx:182-188`.
pub(crate) fn perform_general(state: &mut ExtPcState, curve: &dyn Curve, point: &GpPnt) {
    if let Some(poles) = curve.bezier_poles() {
        state.int_u_inf = state.u_inf;
        state.int_u_sup = state.u_sup;
        state.sample = (poles.len() as i32) * 2;
        state.interval_perform(curve, point);
        return;
    }

    // `GGExtPC.hxx:395-405`: C2 intervals when NbIntervals>1, else DeflCurvIntervals.
    let c2 = curve.parameter_intervals(2);
    let n_c2 = c2.len().saturating_sub(1);
    let intervals = if n_c2 > 1 {
        c2
    } else {
        defl_curv_intervals(curve)
    };
    let n = intervals.len().saturating_sub(1).max(1);
    let mut maxint = 0.0;
    for i in 0..n {
        let dt = intervals[i + 1] - intervals[i];
        if maxint < dt {
            maxint = dt;
        }
    }
    if maxint <= 0.0 {
        maxint = (state.u_sup - state.u_inf).abs().max(CONFUSION);
    }
    let is_periodic = curve.is_periodic();
    let mut int_ext_done = false;
    let mut s2 = 0.0;
    for i in 0..n {
        state.int_u_inf = intervals[i];
        state.int_u_sup = intervals[i + 1];
        // `RealToInt` truncates toward zero (`GGExtPC.hxx:424`).
        state.sample = ((MAX_SAMPLE as f64) * (state.int_u_sup - state.int_u_inf) / maxint) as i32;
        state.sample = state.sample.max(3);

        let mut an_inf = state.int_u_inf;
        let mut sup = state.int_u_sup;
        if is_periodic {
            let period = curve.period();
            an_inf = elclib_in_period(state.int_u_inf, state.u_inf, state.u_inf + period);
            sup = state.int_u_sup + (an_inf - state.int_u_inf);
        }
        if state.u_inf > sup || state.u_sup < an_inf {
            continue;
        }
        if state.u_inf >= an_inf {
            an_inf = state.u_inf;
        }
        if state.u_sup <= sup {
            sup = state.u_sup;
        }
        if (sup - an_inf) <= state.tol_u {
            continue;
        }
        state.int_u_inf = an_inf;
        state.int_u_sup = sup;

        if i != 0 {
            let (pp, v1) = curve.d1(state.int_u_inf);
            let s1 = GpVec::from_pnts(point, &pp).dot(&v1);
            if s1 * s2 < 0.0 {
                state.add_sol(state.int_u_inf, pp, pp.square_distance(point), s1 < 0.0);
            }
        }
        if i + 1 != n {
            let (pp, v1) = curve.d1(state.int_u_sup);
            s2 = GpVec::from_pnts(point, &pp).dot(&v1);
        }
        state.interval_perform(curve, point);
        int_ext_done = int_ext_done || state.done;
    }
    state.done = int_ext_done;
}

pub(crate) fn postprocess_ends(state: &mut ExtPcState) {
    // `GGExtPC.hxx:474-502` for BSpline / Offset / OtherCurve.
    if state.dist1 >= CONFUSION * CONFUSION && state.dist2 >= CONFUSION * CONFUSION {
        return;
    }
    let mut is_first = false;
    let mut is_last = false;
    for s in &state.sols {
        if (s.u - state.u_inf).abs() < state.tol_u {
            is_first = true;
        } else if (state.u_sup - s.u).abs() < state.tol_u {
            is_last = true;
        }
    }
    if !is_first && state.dist1 < CONFUSION * CONFUSION {
        state.sols.insert(
            0,
            ExtPcSolution {
                u: state.u_inf,
                point: state.pf,
                sq_dist: state.dist1,
                is_min: true,
            },
        );
    }
    if !is_last && state.dist2 < CONFUSION * CONFUSION {
        state.add_sol(state.u_sup, state.pl, state.dist2, true);
    }
    state.done = true;
}

