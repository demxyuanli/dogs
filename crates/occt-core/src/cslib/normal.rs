//! `CSLib::Normal` second-derivative and MaxOrder overloads.
//! Source: `CSLib.cxx:84-151` (D2) and `CSLib.cxx:183-387` (MaxOrder).

use crate::cslib::poly_def::NormalPolyDef;
use crate::gp::{GpDir, GpVec};
use crate::math_fn::MathFunction;
use crate::math_function_roots::FunctionRoots;
use crate::precision::{COMPUTATIONAL, CONFUSION, PCONFUSION};

/// Angular tolerance for MaxOrder lambda ratios (`CSLib.cxx:34`).
const PARALLEL_ANGULAR_TOL: f64 = 1e-6;
/// `THE_MAX_ROOT_ITERATIONS` (`CSLib.cxx:37`).
const MAX_ROOT_ITERATIONS: i32 = 200;
/// `THE_ROOT_FINDING_TOL` (`CSLib.cxx:40`).
const ROOT_FINDING_TOL: f64 = 1e-5;

/// `CSLib_NormalStatus` (`CSLib_NormalStatus.hxx:24-35`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CSLibNormalStatus {
    Singular,
    Defined,
    InfinityOfSolutions,
    D1NuIsNull,
    D1NvIsNull,
    D1NIsNull,
    D1NuNvRatioIsNull,
    D1NvNuRatioIsNull,
    D1NuIsParallelD1Nv,
}

/// `CSLib::Normal(D1U, D1V, MagTol, Status, Normal)` (`CSLib.cxx:156-179`).
pub fn normal_d1_mag(d1u: &GpVec, d1v: &GpVec, mag_tol: f64) -> (CSLibNormalStatus, Option<GpDir>) {
    let d1u_mag = d1u.magnitude();
    let d1v_mag = d1v.magnitude();
    let cross = d1u.crossed(d1v);
    let n_mag = cross.magnitude();
    if n_mag <= mag_tol || d1u_mag <= mag_tol || d1v_mag <= mag_tol {
        return (CSLibNormalStatus::Singular, None);
    }
    match (GpDir::from_vec(d1u), GpDir::from_vec(d1v)) {
        (Ok(du), Ok(dv)) => match du.crossed(&dv) {
            Ok(n) => (CSLibNormalStatus::Defined, Some(n)),
            Err(_) => (CSLibNormalStatus::Singular, None),
        },
        _ => (CSLibNormalStatus::Singular, None),
    }
}

/// `CSLib::Normal(D1U, D1V, D2U, D2V, D2UV, SinTol, Done, Status, Normal)`.
/// Source: `CSLib.cxx:84-151`.
pub fn normal_d2(
    d1u: &GpVec,
    d1v: &GpVec,
    d2u: &GpVec,
    d2v: &GpVec,
    d2uv: &GpVec,
    sin_tol: f64,
) -> (bool, CSLibNormalStatus, Option<GpDir>) {
    let d1_nu = d2u.crossed(d1v).added(&d1u.crossed(d2uv));
    let d1_nv = d2uv.crossed(d1v).added(&d1u.crossed(d2v));
    let ld1_nu = d1_nu.square_magnitude();
    let ld1_nv = d1_nv.square_magnitude();

    if ld1_nu <= COMPUTATIONAL && ld1_nv <= COMPUTATIONAL {
        return (false, CSLibNormalStatus::D1NIsNull, None);
    }
    if ld1_nu < COMPUTATIONAL {
        return match GpDir::from_vec(&d1_nv) {
            Ok(n) => (true, CSLibNormalStatus::D1NuIsNull, Some(n)),
            Err(_) => (false, CSLibNormalStatus::D1NIsNull, None),
        };
    }
    if ld1_nv < COMPUTATIONAL {
        return match GpDir::from_vec(&d1_nu) {
            Ok(n) => (true, CSLibNormalStatus::D1NvIsNull, Some(n)),
            Err(_) => (false, CSLibNormalStatus::D1NIsNull, None),
        };
    }
    if (ld1_nv / ld1_nu) <= COMPUTATIONAL {
        return (false, CSLibNormalStatus::D1NvNuRatioIsNull, None);
    }
    if (ld1_nu / ld1_nv) <= COMPUTATIONAL {
        return (false, CSLibNormalStatus::D1NuNvRatioIsNull, None);
    }

    let d1n_cross = d1_nu.crossed(&d1_nv);
    let sin2 = d1n_cross.square_magnitude() / (ld1_nu * ld1_nv);
    if sin2 < sin_tol * sin_tol {
        match GpDir::from_vec(&d1_nu) {
            Ok(n) => (true, CSLibNormalStatus::D1NuIsParallelD1Nv, Some(n)),
            Err(_) => (false, CSLibNormalStatus::D1NIsNull, None),
        }
    } else {
        (false, CSLibNormalStatus::InfinityOfSolutions, None)
    }
}

/// `CSLib::Normal(MaxOrder, DerNUV, ...)` (`CSLib.cxx:183-387`).
/// `der_nuv[u][v]` is `theDerNUV(u, v)`.
pub fn normal_max_order(
    max_order: i32,
    der_nuv: &[Vec<GpVec>],
    sin_tol: f64,
    u: f64,
    v: f64,
    u_min: f64,
    u_max: f64,
    v_min: f64,
    v_max: f64,
) -> (CSLibNormalStatus, Option<GpDir>, i32, i32) {
    let mut an_order = -1;
    let mut found_u_idx = 0;
    let mut found = false;
    let mut a_d = GpVec::new(0.0, 0.0, 0.0);
    while !found && an_order < max_order {
        an_order += 1;
        found_u_idx = an_order;
        while found_u_idx >= 0 && !found {
            let v_idx = an_order - found_u_idx;
            a_d = der_nuv[found_u_idx as usize][v_idx as usize];
            let norme = a_d.magnitude();
            found = norme >= sin_tol;
            found_u_idx -= 1;
        }
    }
    let order_u = found_u_idx + 1;
    let order_v = an_order - order_u;
    if !found {
        return (CSLibNormalStatus::Singular, None, order_u, order_v);
    }
    if an_order == 0 {
        return match GpDir::from_vec(&a_d) {
            Ok(n) => (CSLibNormalStatus::Defined, Some(n), order_u, order_v),
            Err(_) => (CSLibNormalStatus::Singular, None, order_u, order_v),
        };
    }
    let vk0 = der_nuv[order_u as usize][order_v as usize];
    let mut ratio = vec![0.0; (an_order as usize) + 1];
    let mut ratio_idx = 0;
    let mut is_defined = false;
    while ratio_idx <= an_order && !is_defined {
        let der_vec = der_nuv[ratio_idx as usize][(an_order - ratio_idx) as usize];
        if der_vec.magnitude() <= sin_tol {
            ratio[ratio_idx as usize] = 0.0;
        } else if der_vec.is_parallel_ang(&vk0, PARALLEL_ANGULAR_TOL) {
            let mut mag_ratio = der_vec.magnitude() / vk0.magnitude();
            if der_vec.is_opposite_ang(&vk0, PARALLEL_ANGULAR_TOL) {
                mag_ratio = -mag_ratio;
            }
            ratio[ratio_idx as usize] = mag_ratio;
        } else {
            is_defined = true;
        }
        ratio_idx += 1;
    }
    if is_defined {
        return match GpDir::from_vec(&a_d) {
            Ok(n) => (CSLibNormalStatus::Defined, Some(n), order_u, order_v),
            Err(_) => (CSLibNormalStatus::Singular, None, order_u, order_v),
        };
    }
    let mut a_inf = -std::f64::consts::PI;
    let mut a_sup = std::f64::consts::PI;
    let is_fu = (u - u_min).abs() < PCONFUSION;
    let is_lu = (u - u_max).abs() < PCONFUSION;
    let is_fv = (v - v_min).abs() < PCONFUSION;
    let is_lv = (v - v_max).abs() < PCONFUSION;
    if is_lu {
        a_inf = std::f64::consts::FRAC_PI_2;
        a_sup = 3.0 * std::f64::consts::FRAC_PI_2;
        if is_lv {
            a_inf = std::f64::consts::PI;
        }
        if is_fv {
            a_sup = std::f64::consts::PI;
        }
    } else if is_fu {
        a_sup = std::f64::consts::FRAC_PI_2;
        a_inf = -std::f64::consts::FRAC_PI_2;
        if is_lv {
            a_sup = 0.0;
        }
        if is_fv {
            a_inf = 0.0;
        }
    } else if is_lv {
        a_inf = -std::f64::consts::PI;
        a_sup = 0.0;
    } else if is_fv {
        a_inf = 0.0;
        a_sup = std::f64::consts::PI;
    }
    let mut changes_sign = false;
    let mut v_prec = 0.0;
    let mut v_suiv = 0.0;
    let mut poly = NormalPolyDef::new(an_order, &ratio);
    let find_roots = FunctionRoots::new(
        &mut poly,
        a_inf,
        a_sup,
        MAX_ROOT_ITERATIONS,
        ROOT_FINDING_TOL,
        CONFUSION,
        CONFUSION,
        0.0,
    );
    if find_roots.is_done() && find_roots.nb_solutions() > 0 {
        let nb_sol = find_roots.nb_solutions();
        let mut sol = vec![0.0; (nb_sol as usize) + 2];
        for root_idx in 1..=nb_sol {
            sol[root_idx as usize] = find_roots.value(root_idx);
        }
        sol[1..=nb_sol as usize].sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        sol[0] = a_inf;
        sol[(nb_sol as usize) + 1] = a_sup;
        let mut first = 0;
        for interval_idx in 0..=nb_sol {
            if (sol[(interval_idx as usize) + 1] - sol[interval_idx as usize]).abs() > PCONFUSION {
                poly.value(
                    (sol[interval_idx as usize] + sol[(interval_idx as usize) + 1]) / 2.0,
                    &mut v_suiv,
                );
                if first == 0 {
                    first = interval_idx;
                    changes_sign = false;
                    v_prec = v_suiv;
                } else {
                    changes_sign = changes_sign || (v_prec * v_suiv) < 0.0;
                    v_prec = v_suiv;
                }
            }
        }
    } else {
        changes_sign = false;
        poly.value(a_inf, &mut v_suiv);
    }
    if changes_sign {
        (CSLibNormalStatus::InfinityOfSolutions, None, order_u, order_v)
    } else {
        let a_sign = if v_suiv > 0.0 { 1.0 } else { -1.0 };
        match GpDir::from_vec(&vk0) {
            Ok(n) => {
                let n = if a_sign < 0.0 { n.reversed() } else { n };
                (CSLibNormalStatus::Defined, Some(n), order_u, order_v)
            }
            Err(_) => (CSLibNormalStatus::Singular, None, order_u, order_v),
        }
    }
}
