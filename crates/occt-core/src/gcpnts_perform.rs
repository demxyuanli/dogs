//! `GCPnts_TangentialDeflection::PerformCurve` / `PerformLinear` / `EvaluateDu`.
//!
//! Source: `GCPnts_TangentialDeflection.cxx:321-348, 302-322, 522-955`.

use crate::gp::{GpPnt, GpVec};
use crate::precision::CONFUSION;

#[path = "gcpnts_estim.rs"]
mod estim;
use estim::estim_defl;

const US3: f64 = 1.0 / 3.0;

/// Point + second derivative at a parameter (`Adaptor3d_Curve::D0` / `D2`).
pub trait CurveSecondDeriv {
    fn point(&self, u: f64) -> GpPnt;
    fn d2(&self, u: f64) -> (GpPnt, GpVec, GpVec);
}

/// `GCPnts_TangentialDeflection::PerformLinear`.
pub fn perform_linear<C: CurveSecondDeriv>(
    curve: &C,
    first: f64,
    last: f64,
    min_nb: usize,
) -> (Vec<f64>, Vec<GpPnt>) {
    let mut params = Vec::new();
    let mut points = Vec::new();
    params.push(first);
    points.push(curve.point(first));
    if min_nb > 2 {
        let du = (last - first) / min_nb as f64;
        let mut u = first + du;
        for _ in 2..min_nb {
            params.push(u);
            points.push(curve.point(u));
            u += du;
        }
    }
    params.push(last);
    points.push(curve.point(last));
    (params, points)
}

/// `GCPnts_TangentialDeflection::PerformCurve` (`cxx:522-916`).
///
/// `intervals` are CN breakpoints including the range ends (OCCT
/// `NbIntervals(GeomAbs_CN)` + `Intervals`). A single `[first, last]` pair is
/// valid when the curve exposes no knot structure.
///
/// `degree_min_nb` is `max(Degree+1, min_nb)` for BSpline/Bezier (`cxx:568-578`).
pub fn perform_tangential_curve<C: CurveSecondDeriv>(
    curve: &C,
    first: f64,
    last: f64,
    angular_deflection: f64,
    curvature_deflection: f64,
    min_nb: usize,
    u_tol: f64,
    min_len: f64,
    intervals: &[f64],
    degree_min_nb: usize,
) -> (Vec<f64>, Vec<GpPnt>) {
    let (first, last) = if first <= last { (first, last) } else { (last, first) };
    let min_nb = min_nb.max(2);
    let min_len = min_len.max(CONFUSION);
    let u_tol = if u_tol > 0.0 { u_tol } else { 1.0e-9 };

    let last_point = curve.point(last);
    let mut params = Vec::new();
    let mut points = Vec::new();

    let mut dusave = (last - first) * US3;
    let mut du = dusave;
    let mut current_point = GpPnt::zero();
    let mut not_done = true;
    evaluate_du(
        curve,
        first,
        curvature_deflection,
        min_len,
        &mut current_point,
        &mut du,
        &mut not_done,
    );
    params.push(first);
    points.push(current_point);

    let intervs = if intervals.len() >= 2 {
        intervals.to_vec()
    } else {
        vec![first, last]
    };

    if not_done || du > 5.0 * dusave {
        let v1 = GpVec::from_pnts(&current_point, &last_point);
        let l1 = v1.magnitude();
        if l1 > CONFUSION {
            let mut is_line = true;
            let nb_points = degree_min_nb.max(min_nb).max(3);
            let mut sample_u = 0.0;
            let mut i = 0usize;
            while i + 1 < intervs.len() && is_line {
                let mut u0 = intervs[i];
                let mut u1 = intervs[i + 1];
                if u1 < first || u0 > last {
                    i += 1;
                    continue;
                }
                if u0 < first && u1 > first {
                    u0 = first;
                }
                if u0 < last && u1 > last {
                    u1 = last;
                }
                let delta = (u1 - u0) / nb_points as f64;
                let mut j = 1usize;
                while j <= nb_points && is_line {
                    sample_u = u0 + j as f64 * delta;
                    let mid = curve.point(sample_u);
                    let v2 = GpVec::from_pnts(&current_point, &mid);
                    let l2 = v2.magnitude();
                    if l2 > CONFUSION {
                        let angle = v2.cross_magnitude(&v1) / (l1 * l2);
                        let mut a_tol = 1.0e-2 * angular_deflection;
                        if a_tol > 1.0e-2 {
                            a_tol = 1.0e-2;
                        } else if a_tol < 1.0e-7 {
                            a_tol = 1.0e-7;
                        }
                        is_line = angle < a_tol;
                    }
                    j += 1;
                }
                i += 1;
            }
            if is_line {
                return perform_linear(curve, first, last, min_nb);
            }
            let mut dummy = GpPnt::zero();
            evaluate_du(
                curve,
                sample_u,
                curvature_deflection,
                min_len,
                &mut dummy,
                &mut du,
                &mut not_done,
            );
        } else {
            du = (last - first) / 2.1;
            let middle_u = first + du;
            let mid = curve.point(middle_u);
            let l1 = current_point.distance(&mid);
            if l1 < CONFUSION {
                params.push(last);
                points.push(last_point);
                return fill_min_points(curve, params, points, min_nb);
            }
        }
    }

    if du > dusave {
        du = dusave;
    } else {
        dusave = du;
    }
    if du < u_tol {
        du = last - first;
        if du < u_tol {
            params.push(last);
            points.push(last_point);
            return fill_min_points(curve, params, points, min_nb);
        }
    }

    let angle_max = angular_deflection * 0.5;
    let mut more_points = true;
    let mut u1 = first;
    let mut u2 = first;
    let mut idx0 = 0usize;
    let mut need_check = false;
    let mut prev_point = points[points.len() - 1];
    let mut guard = 0u32;

    while more_points {
        guard += 1;
        if guard > 1_000_000 {
            break;
        }
        idx0 = interval_idx(u1, &intervs, idx0);
        u2 += du;
        if u2 >= last {
            u2 = last;
            current_point = last_point;
            du = u2 - u1;
            dusave = du;
        } else {
            current_point = curve.point(u2);
        }

        let mut coef = 0.0;
        let mut a_coef = 0.0;
        let mut f_coef = 0.0;
        let mut too_large = false;
        let mut correction = true;
        let mut too_small = false;
        let mut corr_guard = 0u32;

        while correction {
            corr_guard += 1;
            if corr_guard > 10_000 {
                break;
            }
            if need_check {
                let idx1 = interval_idx(u2, &intervs, idx0);
                if idx1 > idx0 && idx0 + 1 < intervs.len() {
                    let span = (intervs[idx0 + 1] - intervs[idx0]) * US3;
                    if du > span {
                        du = span;
                        u2 = u1 + du;
                        if u2 > last {
                            u2 = last;
                        }
                        current_point = curve.point(u2);
                    }
                }
            }
            let middle_u = (u1 + u2) * 0.5;
            let mid = curve.point(middle_u);

            let v1 = GpVec::from_pnts(&prev_point, &current_point);
            let v2 = GpVec::from_pnts(&prev_point, &mid);
            let l1 = v1.magnitude();
            f_coef = if l1 > min_len {
                v1.cross_magnitude(&v2) / (l1 * curvature_deflection)
            } else {
                0.0
            };

            let v1b = GpVec::from_pnts(&mid, &current_point);
            let l1b = v1b.magnitude();
            let l2 = v2.magnitude();
            a_coef = if l1b > min_len && l2 > min_len {
                (v1b.cross_magnitude(&v2) / (l1b * l2)) / angle_max
            } else {
                0.0
            };
            coef = a_coef.max(f_coef);

            if need_check && coef < 0.55 {
                need_check = false;
                du = dusave;
                u2 = u1 + du;
                if u2 > last {
                    u2 = last;
                }
                current_point = curve.point(u2);
                continue;
            }

            if coef <= 1.0 {
                if (last - u2).abs() < u_tol {
                    params.push(last);
                    points.push(last_point);
                    more_points = false;
                    correction = false;
                } else if coef >= 0.55 || too_large {
                    params.push(u2);
                    points.push(current_point);
                    prev_point = current_point;
                    correction = false;
                    need_check = true;
                } else if too_small {
                    correction = false;
                    prev_point = current_point;
                } else {
                    too_small = true;
                    du += ((u2 - u1) * (1.0 - coef)).min(du * US3);
                    u2 = u1 + du;
                    if u2 > last {
                        u2 = last;
                    }
                    current_point = curve.point(u2);
                }
            } else if coef >= 1.5 {
                let last_stored = points[points.len() - 1];
                if prev_point.distance(&last_stored) > CONFUSION {
                    params.push(u1);
                    points.push(prev_point);
                }
                u2 = middle_u;
                du = u2 - u1;
                current_point = mid;
            } else {
                du *= 0.9;
                u2 = u1 + du;
                current_point = curve.point(u2);
                too_large = true;
            }
        }

        du = u2 - u1;
        if more_points {
            if u1 > first {
                if f_coef > a_coef {
                    let mut nd = true;
                    evaluate_du(
                        curve,
                        u2,
                        curvature_deflection,
                        min_len,
                        &mut current_point,
                        &mut du,
                        &mut nd,
                    );
                    if nd {
                        du += (du - dusave) * (du / dusave);
                        if du > 1.5 * dusave {
                            du = 1.5 * dusave;
                        }
                        if du < 0.75 * dusave {
                            du = 0.75 * dusave;
                        }
                    }
                } else {
                    du += (du - dusave) * (du / dusave);
                    if du > 1.5 * dusave {
                        du = 1.5 * dusave;
                    }
                    if du < 0.75 * dusave {
                        du = 0.75 * dusave;
                    }
                }
            }
            if du < u_tol {
                du = last - u2;
                if du < u_tol {
                    params.push(last);
                    points.push(last_point);
                    more_points = false;
                } else if du * US3 > u_tol {
                    du *= US3;
                }
            }
            u1 = u2;
            dusave = du;
        }
    }

    // Second-to-last midpoint snap (`cxx:880-896`).
    let n = points.len();
    if n >= 3 {
        let middle_u = 0.5 * (last + params[n - 3]);
        params[n - 2] = middle_u;
        points[n - 2] = curve.point(middle_u);
    }

    // `cxx:905-915` fill min points, then `cxx:916-955` EstimDefl splits.
    let (mut params, mut points) = fill_min_points(curve, params, points, min_nb);
    split_by_interval_defl(
        curve,
        &mut params,
        &mut points,
        curvature_deflection,
        angle_max,
        u_tol,
        min_len,
        last - first,
    );
    (params, points)
}

/// `EvaluateDu` (`cxx:302-322`).
fn evaluate_du<C: CurveSecondDeriv>(
    curve: &C,
    u: f64,
    curvature_deflection: f64,
    min_len: f64,
    p: &mut GpPnt,
    du: &mut f64,
    not_done: &mut bool,
) {
    let (pt, t, n) = curve.d2(u);
    *p = pt;
    let lt = t.magnitude();
    if lt > CONFUSION && n.magnitude() > CONFUSION {
        let lc = n.cross_magnitude(&t);
        let ln = lc / lt;
        if ln > CONFUSION {
            *du = (8.0 * curvature_deflection.max(min_len) / ln).sqrt();
            *not_done = false;
        }
    }
}

fn interval_idx(param: f64, intervs: &[f64], previous: usize) -> usize {
    let mut idx = previous;
    while idx + 1 < intervs.len() {
        if param >= intervs[idx] && param <= intervs[idx + 1] {
            break;
        }
        idx += 1;
    }
    idx
}

/// `EstimAngl` (`cxx:69-81`).
fn estim_angl(p1: &GpPnt, pm: &GpPnt, p2: &GpPnt) -> f64 {
    let v1 = GpVec::from_pnts(p1, pm);
    let v2 = GpVec::from_pnts(pm, p2);
    let l = v1.magnitude() * v2.magnitude();
    if l > 1e-12 {
        v1.cross_magnitude(&v2) / l
    } else {
        0.0
    }
}

/// Interval split when `EstimDefl` / `EstimAngl` exceed the tolerances
/// (`cxx:916-955`).
fn split_by_interval_defl<C: CurveSecondDeriv>(
    curve: &C,
    params: &mut Vec<f64>,
    points: &mut Vec<GpPnt>,
    curvature_deflection: f64,
    angle_max: f64,
    u_tol: f64,
    min_len: f64,
    curve_span: f64,
) {
    let min_len2 = min_len * min_len;
    let mut i = 0usize;
    let mut nbp = points.len();
    let max_nbp = 10 * nbp.max(1);
    while i + 1 < points.len() {
        let u1 = params[i];
        let u2 = params[i + 1];
        if u2 - u1 <= u_tol {
            i += 1;
            continue;
        }
        let (dmax, umax) = estim_defl(curve, u1, u2, u_tol, curve_span);
        let p1 = points[i];
        let p2 = points[i + 1];
        let mid = curve.point(umax);
        let amax = estim_angl(&p1, &mid, &p2);
        if dmax > curvature_deflection || amax > angle_max {
            if umax - u1 > u_tol && u2 - umax > u_tol {
                if p1.square_distance(&mid) > min_len2 && p2.square_distance(&mid) > min_len2 {
                    params.insert(i + 1, umax);
                    points.insert(i + 1, mid);
                    nbp += 1;
                    if nbp > max_nbp {
                        break;
                    }
                    continue;
                }
            }
        }
        i += 1;
    }
}

fn fill_min_points<C: CurveSecondDeriv>(
    curve: &C,
    mut params: Vec<f64>,
    mut points: Vec<GpPnt>,
    min_nb: usize,
) -> (Vec<f64>, Vec<GpPnt>) {
    while points.len() < min_nb {
        let nbp = points.len();
        if nbp < 2 {
            break;
        }
        let mut i = 1usize;
        while i < points.len() {
            let middle_u = 0.5 * (params[i - 1] + params[i]);
            let mid = curve.point(middle_u);
            params.insert(i, middle_u);
            points.insert(i, mid);
            i += 2;
        }
        if points.len() == nbp {
            break;
        }
    }
    (params, points)
}
