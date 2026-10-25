//! Curve reparametrization helpers of `GeomFill_Gordon`.
//!
//! Source: `ModelingAlgorithms/TKGeomAlgo/GeomFill/GeomFill_Gordon.cxx`
//! (OCCT 8.0.0), anonymous namespace helpers: `getKinkParameters`
//! (`cxx:101-126`), `linspaceWithBreaks` (`cxx:138-192`),
//! `invertReparamBisection` (`cxx:200-238`) and `reparamCurve`
//! (`cxx:250-429`), plus the file-local constants (`cxx:33-83`).
//!
//! `reparamCurve` builds its 2D interpolation with `Geom2dAPI_Interpolate`
//! (`cxx:306-310`); the port calls
//! [`occt_geom2d::geom2d_api_interpolate::Interpolate`], which is the port of
//! that class' non-periodic path. OCCT passes `PeriodicFlag = false` there, so
//! the unported periodic branch is not reachable.

use occt_core::gp::GpPnt2d;
use occt_core::precision::{CONFUSION, PCONFUSION};
use occt_geom2d::geom2d_api_interpolate::Interpolate;
use occt_geom2d::Curve2d;

use crate::approx::bspline_approx_interp::BSplineApproxInterp;
use crate::bspline_curve::GeomBSplineCurve;
use crate::curve::Curve;

/// Kink angle threshold in radians (6 degrees) (`cxx:38-39`).
const THE_KINK_ANGLE_RAD: f64 = 6.0 / 180.0 * std::f64::consts::PI;
/// Epsilon offset for evaluating derivatives near knots (`cxx:41-42`).
const THE_KINK_EPS: f64 = 1.0e-8;
/// Relative threshold for replacing uniform samples by break points (`cxx:52-53`).
const THE_REPLACE_EPS: f64 = 0.3;
/// Epsilon for duplicate-parameter compaction (`cxx:55-56`).
const THE_PARAM_UNIQUE_EPS: f64 = 1.0e-14;
/// Maximum iterations for reparametrization inversion by bisection (`cxx:58-59`).
const THE_MAX_BISECT_ITER: usize = 50;
/// Convergence tolerance for bisection-based inversion (`cxx:61-62`).
const THE_BISECT_TOL: f64 = 1.0e-12;
/// Parameter tolerance for duplicate/interior filtering (`cxx:64-65`).
const THE_REPARAM_PARAM_TOL: f64 = 1.0e-10;
/// Offset added to current control-point count when selecting sample count (`cxx:67-68`).
const THE_SAMPLE_COUNT_CP_OFFSET: i32 = 10;
/// Offset added to break count when selecting sample count (`cxx:69-70`).
const THE_SAMPLE_COUNT_BREAK_OFFSET: i32 = 2;
/// Minimum sample count for robust approximation (`cxx:72-73`).
const THE_SAMPLE_COUNT_MIN: i32 = 30;
/// Maximum sample count for robust approximation (`cxx:75-76`).
const THE_SAMPLE_COUNT_MAX: i32 = 200;
/// Target approximation degree during curve reparametrization (`cxx:78-79`).
const THE_APPROX_DEGREE: i32 = 3;
/// Divisor used for rough control-point estimation from sample count (`cxx:81-82`).
const THE_CP_ESTIMATE_DIV: usize = 3;

/// `getKinkParameters` (`cxx:101-126`): parameters where the knot
/// multiplicity equals the degree and the two one-sided tangents differ by
/// more than 6 degrees.
pub fn kink_parameters(curve: &GeomBSplineCurve) -> Vec<f64> {
    let mut kinks = Vec::new();
    let degree = curve.degree();
    let knots = curve.distinct_knots().to_vec();
    let mults = curve.multiplicities().to_vec();
    // 1-based `aKnotIdx` from 2 to `NbKnots() - 1`.
    for knot_idx in 2..knots.len() {
        if mults[knot_idx - 1] as usize != degree {
            continue;
        }
        let knot = knots[knot_idx - 1];
        let tan_left = curve.eval_dn(knot - THE_KINK_EPS, 1);
        let tan_right = curve.eval_dn(knot + THE_KINK_EPS, 1);
        if tan_left.magnitude() > CONFUSION && tan_right.magnitude() > CONFUSION {
            if tan_left.angle(&tan_right) > THE_KINK_ANGLE_RAD {
                kinks.push(knot);
            }
        }
    }
    kinks
}

/// `linspaceWithBreaks` (`cxx:138-192`): uniform grid with the break
/// parameters snapped to nearby grid points (or appended) and compacted.
pub fn linspace_with_breaks(umin: f64, umax: f64, nb_values: i32, breaks: &[f64]) -> Vec<f64> {
    let du = (umax - umin) / (nb_values - 1) as f64;
    let mut result: Vec<f64> = Vec::with_capacity(nb_values.max(1) as usize + breaks.len());
    for value_idx in 0..nb_values {
        result.push(umin + value_idx as f64 * du);
    }
    for &brk in breaks {
        let mut is_replaced = false;
        for grid in result.iter_mut() {
            if (*grid - brk).abs() < du * THE_REPLACE_EPS {
                *grid = brk;
                is_replaced = true;
                break;
            }
        }
        if !is_replaced {
            result.push(brk);
        }
    }
    result.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    // Compact duplicates in place (`cxx:178-191`).
    let mut write_idx = 0usize;
    for unique_idx in 1..result.len() {
        if (result[unique_idx] - result[write_idx]).abs() > THE_PARAM_UNIQUE_EPS {
            write_idx += 1;
            result[write_idx] = result[unique_idx];
        }
    }
    result.truncate(write_idx + 1);
    result
}

/// `invertReparamBisection` (`cxx:200-238`): the new parameter whose
/// reparametrization value equals `old_param`.
pub fn invert_reparam_bisection(reparam_func: &dyn Curve2d, old_param: f64) -> f64 {
    let mut lo = reparam_func.first_parameter();
    let mut hi = reparam_func.last_parameter();
    let f_lo = reparam_func.d0(lo).x();
    let f_hi = reparam_func.d0(hi).x();
    let is_increasing = f_hi > f_lo;
    for _ in 0..THE_MAX_BISECT_ITER {
        let mid = 0.5 * (lo + hi);
        let f_mid = reparam_func.d0(mid).x();
        if (f_mid - old_param).abs() < THE_BISECT_TOL {
            return mid;
        }
        let is_mid_below = if is_increasing {
            f_mid < old_param
        } else {
            f_mid > old_param
        };
        if is_mid_below {
            lo = mid;
        } else {
            hi = mid;
        }
        if (hi - lo).abs() < THE_BISECT_TOL {
            return 0.5 * (lo + hi);
        }
    }
    0.5 * (lo + hi)
}

/// `reparamCurve` (`cxx:250-429`): approximate `theCurve` by a B-spline whose
/// parameters at the old intersection parameters become the target ones.
///
/// `Err` marks OCCT's `return false` paths (`cxx:309-312`, `:421-424`).
pub fn reparam_curve(
    curve: &mut GeomBSplineCurve,
    old_params: &[f64],
    new_params: &[f64],
) -> Result<(), String> {
    let mut old_seq: Vec<f64> = vec![curve.first_parameter()];
    let mut new_seq: Vec<f64> = vec![0.0];
    for idx in 0..old_params.len() {
        let old = old_params[idx];
        let new = new_params[idx];
        let last_old = *old_seq.last().unwrap();
        let last_new = *new_seq.last().unwrap();
        if (old - last_old).abs() < THE_REPARAM_PARAM_TOL
            || (new - last_new).abs() < THE_REPARAM_PARAM_TOL
        {
            continue;
        }
        if (old - curve.last_parameter()).abs() < THE_REPARAM_PARAM_TOL
            || (new - 1.0).abs() < THE_REPARAM_PARAM_TOL
        {
            continue;
        }
        old_seq.push(old);
        new_seq.push(new);
    }
    let last_old = *old_seq.last().unwrap();
    if old_seq.len() < 2 || (last_old - curve.last_parameter()).abs() > THE_REPARAM_PARAM_TOL {
        old_seq.push(curve.last_parameter());
        new_seq.push(1.0);
    }
    if old_seq.len() < 2 {
        return Ok(());
    }

    let points: Vec<GpPnt2d> = old_seq.iter().map(|v| GpPnt2d::new(*v, 0.0)).collect();
    let mut interp = Interpolate::with_parameters(points, new_seq, false, PCONFUSION)
        .map_err(|_| "Geom2dAPI_Interpolate".to_string())?;
    interp.perform();
    let reparam_func = interp
        .curve()
        .ok_or_else(|| "Geom2dAPI_Interpolate::Perform".to_string())?;

    let kinks = kink_parameters(curve);
    let nb_interior_params = new_params.len().saturating_sub(2) as i32;
    let nb_kinks = kinks.len() as i32;
    let nb_breaks = nb_interior_params + nb_kinks;

    let mut kinks_in_new = vec![0.0f64; nb_kinks.max(1) as usize];
    for kink_idx in 0..kinks.len() {
        kinks_in_new[kink_idx] = invert_reparam_bisection(reparam_func, kinks[kink_idx]);
    }

    let max_cp = curve.nb_poles() as i32;
    let mut nb_samples = (max_cp + THE_SAMPLE_COUNT_CP_OFFSET)
        .max(nb_breaks + THE_SAMPLE_COUNT_BREAK_OFFSET);
    nb_samples = THE_SAMPLE_COUNT_MIN.max(THE_SAMPLE_COUNT_MAX.min(nb_samples));

    let sample_params: Vec<f64> = if nb_breaks > 0 {
        let mut breaks: Vec<f64> = new_params[1..new_params.len() - 1].to_vec();
        for kink_idx in 0..kinks.len() {
            breaks.push(kinks_in_new[kink_idx]);
        }
        linspace_with_breaks(0.0, 1.0, nb_samples, &breaks)
    } else {
        (0..nb_samples)
            .map(|i| i as f64 / (nb_samples - 1) as f64)
            .collect()
    };

    let nb_actual = sample_params.len();
    let first_param = curve.first_parameter();
    let last_param = curve.last_parameter();
    let mut sampled_pnts: Vec<occt_core::gp::GpPnt> = Vec::with_capacity(nb_actual);
    let mut sampled_new_params: Vec<f64> = Vec::with_capacity(nb_actual);
    for &new_param in &sample_params {
        let old_pt = reparam_func.d0(new_param);
        let mut old_param = old_pt.x();
        if old_param < first_param {
            old_param = first_param;
        }
        if old_param > last_param {
            old_param = last_param;
        }
        sampled_pnts.push(curve.d0(old_param));
        sampled_new_params.push(new_param);
    }

    let nb_interp = 2 + nb_kinks;
    let mut nb_cp = curve.nb_poles() as i32;
    nb_cp = nb_cp.max((nb_actual as i32) / THE_CP_ESTIMATE_DIV as i32);
    nb_cp = nb_cp.max(nb_interp + THE_APPROX_DEGREE + 1);
    nb_cp = nb_cp.max(THE_APPROX_DEGREE + 1);

    let mut approx = BSplineApproxInterp::new(&sampled_pnts, nb_cp, THE_APPROX_DEGREE, false);
    approx.interpolate_point(0, false);
    approx.interpolate_point(nb_actual as i32 - 1, false);
    for kink_idx in 1..=kinks.len() {
        let kink_new = kinks_in_new[kink_idx - 1];
        let mut best_idx = 0usize;
        let mut best_dist = f64::MAX;
        for (sample_idx, &sample_param) in sampled_new_params.iter().enumerate() {
            let dist = (sample_param - kink_new).abs();
            if dist < best_dist {
                best_dist = dist;
                best_idx = sample_idx;
            }
        }
        if best_idx > 0 && best_idx < nb_actual - 1 {
            approx.interpolate_point(best_idx as i32, true);
        }
    }
    approx.perform_with_params(&sampled_new_params);
    if !approx.is_done() {
        return Err("Approx_BSplineApproxInterp::Perform".to_string());
    }
    if let Some(new_curve) = approx.curve() {
        *curve = new_curve.clone();
    }
    Ok(())
}
