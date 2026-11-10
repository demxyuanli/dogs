//! `Extrema_GGExtPC::Perform` on a B-spline curve.

use super::*;

pub(crate) fn perform_bspline(state: &mut ExtPcState, curve: &dyn Curve, point: &GpPnt) {
    let flat = match curve.bspline_knots() {
        Some(k) if !k.is_empty() => k,
        _ => {
            perform_general(state, curve, point);
            return;
        }
    };
    let degree = curve.nurbs_degree().unwrap_or(1) as i32;
    let (knots, mults) = unique_knots_mults(flat);
    if knots.len() < 2 {
        perform_general(state, curve, point);
        return;
    }
    let first_idx = first_u_knot_index(degree, &mults);
    let last_idx = last_u_knot_index(degree, &mults);
    let tol_coeff = (state.u_sup - state.u_inf) * PCONFUSION;
    let mut period_jump = 0.0;
    if curve.is_periodic() {
        let period = curve.period();
        if period.abs() > 0.0 {
            let mut shift = ((state.u_inf - knots[(first_idx - 1) as usize]) / period) as i32;
            if state.u_inf < knots[(first_idx - 1) as usize] - tol_coeff {
                shift -= 1;
            }
            period_jump = period * shift as f64;
        }
    }

    let mut first_used = first_idx;
    let mut last_used = last_idx;
    for an_idx in first_idx..=last_idx {
        let knot = knots[(an_idx - 1) as usize] + period_jump;
        if state.u_inf >= knot - tol_coeff {
            first_used = an_idx;
        } else {
            break;
        }
    }
    for an_idx in (first_idx..=last_idx).rev() {
        let knot = knots[(an_idx - 1) as usize] + period_jump;
        if state.u_sup <= knot + tol_coeff {
            last_used = an_idx;
        } else {
            break;
        }
    }
    if first_used == last_used {
        first_used = first_idx;
        last_used = first_idx + 1;
    }

    state.sample = degree + 1;

    if state.sample == 2 {
        // Degree-1 knot-span arm (`GGExtPC.hxx:233-297`).
        let mut tmin = 0.0;
        let mut distmin = f64::MAX;
        let mut pmin = GpPnt::zero();
        let mut min2 = 0.0;
        for an_idx in first_used..last_used {
            let mut f = knots[(an_idx - 1) as usize] + period_jump;
            let mut l = knots[an_idx as usize] + period_jump;
            if an_idx == first_used {
                f = state.u_inf;
            } else if an_idx == last_used - 1 {
                l = state.u_sup;
            }
            let p1_val = curve.d0(f);
            let p2_val = curve.d0(l);
            let base1 = GpVec::from_pnts(point, &p1_val);
            let base2 = GpVec::from_pnts(point, &p2_val);
            let v = GpVec::from_pnts(&p2_val, &p1_val);
            let val1 = v.dot(&base1);
            let val2 = v.dot(&base2);
            let min1 = if an_idx == first_used {
                point.square_distance(&p1_val)
            } else {
                let prev = min2;
                if distmin > prev {
                    distmin = prev;
                    tmin = f;
                    pmin = p1_val;
                }
                prev
            };
            min2 = point.square_distance(&p2_val);
            let min_sq = min1.min(min2);
            let min_der = val1.abs().min(val2.abs());
            if !(Precision::is_infinite(val1) || Precision::is_infinite(val2))
                && (val1 * val2 <= 0.0
                    || min_sq < 100.0 * CONFUSION * CONFUSION
                    || 2.0 * min_der < CONFUSION)
            {
                state.int_u_inf = f;
                state.int_u_sup = l;
                state.interval_perform(curve, point);
            }
        }
        if distmin.is_finite() {
            let mut is_to_add = true;
            for s in &state.sols {
                if !(distmin < s.sq_dist && (s.u - tmin).abs() > state.tol_u) {
                    is_to_add = false;
                    break;
                }
            }
            if is_to_add {
                state.add_sol(tmin, pmin, distmin, true);
            }
        }
    } else {
        // Degree>=2 sample + IntervalPerform (`GGExtPC.hxx:300-385`).
        let n_spans = (last_used - first_used).max(1) as usize;
        let cap = (state.sample as usize) * n_spans + 1;
        let mut vals = Vec::with_capacity(cap);
        let mut params = Vec::with_capacity(cap);
        for an_idx in first_used..last_used {
            let mut f = knots[(an_idx - 1) as usize] + period_jump;
            let mut l = knots[an_idx as usize] + period_jump;
            if an_idx == first_used {
                f = state.u_inf;
            }
            if an_idx == last_used - 1 {
                l = state.u_sup;
            }
            let step = (l - f) / state.sample as f64;
            for pnt_idx in 0..state.sample {
                let cur = f + step * pnt_idx as f64;
                vals.push(curve.d0(cur).square_distance(point));
                params.push(cur);
            }
        }
        vals.push(curve.d0(state.u_sup).square_distance(point));
        params.push(state.u_sup);

        for an_idx in 1..vals.len().saturating_sub(1) {
            if vals[an_idx] <= CONFUSION * CONFUSION {
                state.add_sol(
                    params[an_idx],
                    curve.d0(params[an_idx]),
                    vals[an_idx],
                    true,
                );
            }
            if (vals[an_idx] >= vals[an_idx + 1] && vals[an_idx] >= vals[an_idx - 1])
                || (vals[an_idx] <= vals[an_idx + 1] && vals[an_idx] <= vals[an_idx - 1])
            {
                state.int_u_inf = params[an_idx - 1];
                state.int_u_sup = params[an_idx + 1];
                state.interval_perform(curve, point);
            }
        }

        if state.dist1 > CONFUSION * CONFUSION && state.dist1.is_finite() && params.len() >= 2 {
            let (p1, v1) = curve.d1(params[0]);
            let (p2, v2) = curve.d1(params[1]);
            let base1 = GpVec::from_pnts(point, &p1);
            let base2 = GpVec::from_pnts(point, &p2);
            let val1 = v1.dot(&base1);
            let val2 = v2.dot(&base2);
            if !(Precision::is_infinite(val1) || Precision::is_infinite(val2))
                && (val1 * val2 <= 0.0
                    || base1.dot(&base2) <= 0.0
                    || 2.0 * val1.abs() < CONFUSION)
            {
                state.int_u_inf = params[0];
                state.int_u_sup = params[1];
                state.interval_perform(curve, point);
            }
        }
        if state.dist2 > CONFUSION * CONFUSION && state.dist2.is_finite() && params.len() >= 2 {
            let n = params.len();
            let (p1, v1) = curve.d1(params[n - 2]);
            let (p2, v2) = curve.d1(params[n - 1]);
            let base1 = GpVec::from_pnts(point, &p1);
            let base2 = GpVec::from_pnts(point, &p2);
            let val1 = v1.dot(&base1);
            let val2 = v2.dot(&base2);
            if !(Precision::is_infinite(val1) || Precision::is_infinite(val2))
                && (val1 * val2 <= 0.0
                    || base1.dot(&base2) <= 0.0
                    || 2.0 * val2.abs() < CONFUSION)
            {
                state.int_u_inf = params[n - 2];
                state.int_u_sup = params[n - 1];
                state.interval_perform(curve, point);
            }
        }
    }
    state.done = true;
}

