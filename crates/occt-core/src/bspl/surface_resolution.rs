//! `BSplSLib::Resolution` (`BSplSLib.cxx:3519-3858`).
//!
//! Estimates the parametric steps `(UTolerance, VTolerance)` that correspond
//! to a 3D displacement of `tolerance_3d`, from the max first-derivative
//! bound on the pole grid. `flat_u` / `flat_v` are already-expanded knot
//! sequences (`BSplCLib::KnotSequence`).

use crate::gp::GpPnt;

fn pole(poles: &[Vec<GpPnt>], u: i32, v: i32) -> GpPnt {
    poles[(u - 1) as usize][(v - 1) as usize]
}

fn weight(weights: &[Vec<f64>], u: i32, v: i32) -> f64 {
    weights[(u - 1) as usize][(v - 1) as usize]
}

fn wrap_index(i: i32, len: i32) -> i32 {
    if len <= 0 {
        return 1;
    }
    (i - 1).rem_euclid(len) + 1
}

/// `BSplSLib::Resolution`. Returns `(UTolerance, VTolerance)`.
pub fn bspline_surface_resolution(
    poles: &[Vec<GpPnt>],
    weights: Option<&[Vec<f64>]>,
    flat_u: &[f64],
    flat_v: &[f64],
    deg_u: i32,
    deg_v: i32,
    u_rational: bool,
    v_rational: bool,
    tolerance_3d: f64,
) -> (f64, f64) {
    if poles.is_empty() || poles[0].is_empty() || flat_u.is_empty() || flat_v.is_empty() {
        return (0.0, 0.0);
    }

    let mut max_derivative = [0.0f64, 0.0];
    let poles_len_u = poles.len() as i32;
    let poles_len_v = poles[0].len() as i32;
    let num_poles_u = flat_u.len() as i32 - (deg_u + 1);
    let num_poles_v = flat_v.len() as i32 - (deg_v + 1);
    if num_poles_u < 2 || num_poles_v < 1 {
        return (0.0, 0.0);
    }

    let mut min_weights = 0.0;
    if u_rational || v_rational {
        if let Some(w) = weights {
            min_weights = w.iter().flatten().copied().fold(f64::INFINITY, f64::min);
        }
    }

    let ud1 = deg_u + 1;
    let vd1 = deg_v + 1;

    if u_rational {
        if let Some(w) = weights {
            let ud2 = deg_u << 1;
            let vd2 = deg_v << 1;
            for ii in 2..=num_poles_u {
                let ii_index = wrap_index(ii, poles_len_u);
                let ii_minus = wrap_index(ii - 1, poles_len_u);
                let span = flat_u[(ii + deg_u - 1) as usize] - flat_u[(ii - 1) as usize];
                let inverse = 1.0 / span;
                let mut lower_u = ii - ud1;
                if lower_u < 1 {
                    lower_u = 1;
                }
                let mut upper_u = ii + ud2 + 1;
                if upper_u > num_poles_u {
                    upper_u = num_poles_u;
                }
                for jj in 1..=num_poles_v {
                    let jj_index = wrap_index(jj, poles_len_v);
                    let mut lower_v = jj - vd1;
                    if lower_v < 1 {
                        lower_v = 1;
                    }
                    let mut upper_v = jj + vd2 + 1;
                    if upper_v > num_poles_v {
                        upper_v = num_poles_v;
                    }
                    let pij = pole(poles, ii_index, jj_index);
                    let wij = weight(w, ii_index, jj_index);
                    let pmj = pole(poles, ii_minus, jj_index);
                    let wmj = weight(w, ii_minus, jj_index);
                    for pp in lower_u..=upper_u {
                        let pp_index = wrap_index(pp, poles_len_u);
                        for qq in lower_v..=upper_v {
                            let qq_index = wrap_index(qq, poles_len_v);
                            let ppq = pole(poles, pp_index, qq_index);
                            let mut value = 0.0;
                            let mut factor = (ppq.x() - pij.x()) * wij - (ppq.x() - pmj.x()) * wmj;
                            if factor < 0.0 {
                                factor = -factor;
                            }
                            value += factor;
                            factor = (ppq.y() - pij.y()) * wij - (ppq.y() - pmj.y()) * wmj;
                            if factor < 0.0 {
                                factor = -factor;
                            }
                            value += factor;
                            factor = (ppq.z() - pij.z()) * wij - (ppq.z() - pmj.z()) * wmj;
                            if factor < 0.0 {
                                factor = -factor;
                            }
                            value += factor;
                            value *= inverse;
                            if max_derivative[0] < value {
                                max_derivative[0] = value;
                            }
                        }
                    }
                }
            }
            max_derivative[0] /= min_weights;
        }
    } else {
        for ii in 2..=num_poles_u {
            let ii_index = wrap_index(ii, poles_len_u);
            let ii_minus = wrap_index(ii - 1, poles_len_u);
            let span = flat_u[(ii + deg_u - 1) as usize] - flat_u[(ii - 1) as usize];
            let inverse = 1.0 / span;
            for jj in 1..=num_poles_v {
                let jj_index = wrap_index(jj, poles_len_v);
                let pij = pole(poles, ii_index, jj_index);
                let pmj = pole(poles, ii_minus, jj_index);
                let mut value = 0.0;
                let mut factor = pij.x() - pmj.x();
                if factor < 0.0 {
                    factor = -factor;
                }
                value += factor;
                factor = pij.y() - pmj.y();
                if factor < 0.0 {
                    factor = -factor;
                }
                value += factor;
                factor = pij.z() - pmj.z();
                if factor < 0.0 {
                    factor = -factor;
                }
                value += factor;
                value *= inverse;
                if max_derivative[0] < value {
                    max_derivative[0] = value;
                }
            }
        }
    }
    max_derivative[0] *= deg_u as f64;

    if v_rational {
        if let Some(w) = weights {
            let ud2 = deg_u << 1;
            let vd2 = deg_v << 1;
            for ii in 2..=num_poles_v {
                let ii_index = wrap_index(ii, poles_len_v);
                let ii_minus = wrap_index(ii - 1, poles_len_v);
                let span = flat_v[(ii + deg_v - 1) as usize] - flat_v[(ii - 1) as usize];
                let inverse = 1.0 / span;
                let mut lower_v = ii - vd1;
                if lower_v < 1 {
                    lower_v = 1;
                }
                let mut upper_v = ii + vd2 + 1;
                if upper_v > num_poles_v {
                    upper_v = num_poles_v;
                }
                for jj in 1..=num_poles_u {
                    let jj_index = wrap_index(jj, poles_len_u);
                    let mut lower_u = jj - ud1;
                    if lower_u < 1 {
                        lower_u = 1;
                    }
                    let mut upper_u = jj + ud2 + 1;
                    if upper_u > num_poles_u {
                        upper_u = num_poles_u;
                    }
                    let pji = pole(poles, jj_index, ii_index);
                    let wji = weight(w, jj_index, ii_index);
                    let pjm = pole(poles, jj_index, ii_minus);
                    let wjm = weight(w, jj_index, ii_minus);
                    for pp in lower_v..=upper_v {
                        let pp_index = wrap_index(pp, poles_len_v);
                        for qq in lower_u..=upper_u {
                            let qq_index = wrap_index(qq, poles_len_u);
                            let pqp = pole(poles, qq_index, pp_index);
                            let mut value = 0.0;
                            let mut factor = (pqp.x() - pji.x()) * wji - (pqp.x() - pjm.x()) * wjm;
                            if factor < 0.0 {
                                factor = -factor;
                            }
                            value += factor;
                            factor = (pqp.y() - pji.y()) * wji - (pqp.y() - pjm.y()) * wjm;
                            if factor < 0.0 {
                                factor = -factor;
                            }
                            value += factor;
                            factor = (pqp.z() - pji.z()) * wji - (pqp.z() - pjm.z()) * wjm;
                            if factor < 0.0 {
                                factor = -factor;
                            }
                            value += factor;
                            value *= inverse;
                            if max_derivative[1] < value {
                                max_derivative[1] = value;
                            }
                        }
                    }
                }
            }
            max_derivative[1] /= min_weights;
        }
    } else {
        if num_poles_v < 2 {
            return (0.0, 0.0);
        }
        for ii in 2..=num_poles_v {
            let ii_index = wrap_index(ii, poles_len_v);
            let ii_minus = wrap_index(ii - 1, poles_len_v);
            let span = flat_v[(ii + deg_v - 1) as usize] - flat_v[(ii - 1) as usize];
            let inverse = 1.0 / span;
            for jj in 1..=num_poles_u {
                let jj_index = wrap_index(jj, poles_len_u);
                let pji = pole(poles, jj_index, ii_index);
                let pjm = pole(poles, jj_index, ii_minus);
                let mut value = 0.0;
                let mut factor = pji.x() - pjm.x();
                if factor < 0.0 {
                    factor = -factor;
                }
                value += factor;
                factor = pji.y() - pjm.y();
                if factor < 0.0 {
                    factor = -factor;
                }
                value += factor;
                factor = pji.z() - pjm.z();
                if factor < 0.0 {
                    factor = -factor;
                }
                value += factor;
                value *= inverse;
                if max_derivative[1] < value {
                    max_derivative[1] = value;
                }
            }
        }
    }
    max_derivative[1] *= deg_v as f64;
    max_derivative[0] *= std::f64::consts::SQRT_2;
    max_derivative[1] *= std::f64::consts::SQRT_2;
    if max_derivative[0] != 0.0 && max_derivative[1] != 0.0 {
        (
            tolerance_3d / max_derivative[0],
            tolerance_3d / max_derivative[1],
        )
    } else {
        (0.0, 0.0)
    }
}
