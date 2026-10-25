//! `BSplCLib::Eval` with an extrapolation mode (`BSplCLib.cxx:3456-3637`
//! weighted, `:3640-3790` plain).
//!
//! Both overloads sum `Poles * Weights * BsplineBasis` over the span for each
//! requested derivative order and, when the parameter falls outside the flat
//! knot window with `ExtrapMode < Degree`, rerun the sum at the window edge and
//! evaluate the Taylor series of the basis there. `LocalRequest` / `Index`
//! bookkeeping follows OCCT line for line; only the 1-valued `ArrayDimension`
//! of `Law_BSpline` is carried here (`BSplCLib::Eval` itself is templated on
//! that dimension, `BSplCLib.hxx:1360`).

use crate::bspl::eval_basis::eval_bspline_basis;
use crate::bspl::plib_eval::eval_polynomial;

/// Shared preparation of `BSplCLib::Eval` (`cxx:3467-3536` / `:3650-3718`).
struct EvalPlan {
    local_request: i32,
    delta: f64,
    modulus: i32,
    extrapolating: [bool; 2],
    first_non_zero: i32,
    basis: Vec<Vec<f64>>,
}

fn prepare(
    parameter: f64,
    is_periodic: bool,
    derivative_request: i32,
    extrap_mode: [i32; 2],
    degree: i32,
    flat_knots: &[f64],
) -> Option<EvalPlan> {
    let order = degree + 1;
    let mut local_request = derivative_request;
    let mut local_parameter = parameter;
    let mut extrapolating = [false; 2];
    // `FlatKnots(2)` and `FlatKnots(Upper() - 1)` in the 1-based array.
    let k2 = *flat_knots.get(1)?;
    let k_upper_minus_1 = *flat_knots.get(flat_knots.len().checked_sub(2)?)?;
    if is_periodic {
        let period = k_upper_minus_1 - k2;
        while local_parameter > k_upper_minus_1 {
            local_parameter -= period;
        }
        while local_parameter < k2 {
            local_parameter += period;
        }
    }
    if parameter < k2 && local_request < extrap_mode[0] && extrap_mode[0] < degree {
        local_request = extrap_mode[0];
        local_parameter = k2;
        extrapolating[0] = true;
    }
    if parameter > k_upper_minus_1 && local_request < extrap_mode[1] && extrap_mode[1] < degree {
        local_request = extrap_mode[1];
        local_parameter = k_upper_minus_1;
        extrapolating[1] = true;
    }
    let delta = parameter - local_parameter;
    if local_request >= order {
        local_request = degree;
    }
    let modulus = if is_periodic {
        flat_knots.len() as i32 - degree - 1
    } else {
        flat_knots.len() as i32 - degree
    };
    let (first_non_zero, basis) =
        eval_bspline_basis(local_request, order, flat_knots, local_parameter, false).ok()?;
    Some(EvalPlan {
        local_request,
        delta,
        modulus,
        extrapolating,
        first_non_zero,
        basis,
    })
}

/// One homogeneous sum over the span, `cxx:3545-3562` / `:3591-3601`.
fn sum_over_span(
    start: i32,
    order: i32,
    modulus: i32,
    basis_row: &[f64],
    poles: &[f64],
    weights: Option<&[f64]>,
) -> f64 {
    let mut index1 = start;
    let mut value = 0.0;
    for &b in basis_row.iter().take(order as usize) {
        let p = poles.get((index1 - 1) as usize).copied().unwrap_or(0.0);
        let w = weights
            .map(|w| w.get((index1 - 1) as usize).copied().unwrap_or(0.0))
            .unwrap_or(1.0);
        value += p * w * b;
        index1 = index1 % modulus + 1;
    }
    value
}

fn sum_weights_over_span(
    start: i32,
    order: i32,
    modulus: i32,
    basis_row: &[f64],
    weights: &[f64],
) -> f64 {
    let mut index1 = start;
    let mut value = 0.0;
    for &b in basis_row.iter().take(order as usize) {
        let w = weights.get((index1 - 1) as usize).copied().unwrap_or(0.0);
        value += w * b;
        index1 = index1 % modulus + 1;
    }
    value
}

/// `BSplCLib::Eval(Parameter, PeriodicFlag, DerivativeRequest, ExtrapMode,
/// Degree, FlatKnots, ArrayDimension, Poles, Weights, PolesResults,
/// WeightsResults)` (`BSplCLib.cxx:3456-3637`) for `ArrayDimension == 1`.
///
/// Returns `(poles_results, weights_results)`, each `LocalRequest + 1` long,
/// or `None` when `EvalBsplineBasis` fails (OCCT returns without writing).
pub fn eval_1d_weighted(
    parameter: f64,
    is_periodic: bool,
    derivative_request: i32,
    extrap_mode: [i32; 2],
    degree: i32,
    flat_knots: &[f64],
    poles: &[f64],
    weights: &[f64],
) -> Option<(Vec<f64>, Vec<f64>)> {
    let plan = prepare(
        parameter,
        is_periodic,
        derivative_request,
        extrap_mode,
        degree,
        flat_knots,
    )?;
    let order = degree + 1;
    let n_orders = (plan.local_request + 1) as usize;
    let mut results = vec![0.0; n_orders];
    let mut weight_results = vec![0.0; n_orders];
    if !plan.extrapolating[0] && !plan.extrapolating[1] {
        for ii in 0..n_orders {
            results[ii] = sum_over_span(
                plan.first_non_zero,
                order,
                plan.modulus,
                &plan.basis[ii],
                poles,
                Some(weights),
            );
            weight_results[ii] = sum_weights_over_span(
                plan.first_non_zero,
                order,
                plan.modulus,
                &plan.basis[ii],
                weights,
            );
        }
    } else {
        // `cxx:3569-3636`: Taylor coefficient arrays at the clamped parameter,
        // each divided by `ii!` (`Inverse /= ii`), then `PLib::EvalPolynomial`.
        let mut new_request = derivative_request;
        if new_request > degree {
            new_request = degree;
        }
        let mut taylor = vec![0.0; n_orders];
        let mut inverse = 1.0;
        for ii in 1..=n_orders {
            taylor[ii - 1] = sum_over_span(
                plan.first_non_zero,
                order,
                plan.modulus,
                &plan.basis[ii - 1],
                poles,
                Some(weights),
            ) * inverse;
            inverse /= ii as f64;
        }
        eval_polynomial(plan.delta, new_request, degree, 1, &taylor, &mut results);

        let mut weight_taylor = vec![0.0; n_orders];
        inverse = 1.0;
        for ii in 1..=n_orders {
            weight_taylor[ii - 1] = sum_weights_over_span(
                plan.first_non_zero,
                order,
                plan.modulus,
                &plan.basis[ii - 1],
                weights,
            ) * inverse;
            inverse /= ii as f64;
        }
        eval_polynomial(
            plan.delta,
            new_request,
            degree,
            1,
            &weight_taylor,
            &mut weight_results,
        );
    }
    Some((results, weight_results))
}

/// `BSplCLib::Eval(Parameter, PeriodicFlag, DerivativeRequest, ExtrapMode,
/// Degree, FlatKnots, ArrayDimension, Poles, Results)` (`BSplCLib.cxx:3640-3790`)
/// for `ArrayDimension == 1`. Returns `LocalRequest + 1` values.
pub fn eval_1d(    parameter: f64,
    is_periodic: bool,
    derivative_request: i32,
    extrap_mode: [i32; 2],
    degree: i32,
    flat_knots: &[f64],
    poles: &[f64],
) -> Option<Vec<f64>> {
    let plan = prepare(
        parameter,
        is_periodic,
        derivative_request,
        extrap_mode,
        degree,
        flat_knots,
    )?;
    let order = degree + 1;
    let n_orders = (plan.local_request + 1) as usize;
    let mut results = vec![0.0; n_orders];
    if !plan.extrapolating[0] && !plan.extrapolating[1] {
        for ii in 0..n_orders {
            results[ii] = sum_over_span(
                plan.first_non_zero,
                order,
                plan.modulus,
                &plan.basis[ii],
                poles,
                None,
            );
        }
    } else {
        let mut new_request = derivative_request;
        if new_request > degree {
            new_request = degree;
        }
        let mut taylor = vec![0.0; n_orders];
        let mut inverse = 1.0;
        for ii in 1..=n_orders {
            taylor[ii - 1] = sum_over_span(
                plan.first_non_zero,
                order,
                plan.modulus,
                &plan.basis[ii - 1],
                poles,
                None,
            ) * inverse;
            inverse /= ii as f64;
        }
        eval_polynomial(plan.delta, new_request, degree, 1, &taylor, &mut results);
    }
    Some(results)
}

/// `BSplCLib::Eval(Parameter, PeriodicFlag, DerivativeRequest, ExtrapMode,
/// Degree, FlatKnots, ArrayDimension, Poles, Results)` (`BSplCLib.cxx:3640-3790`)
/// for a general `ArrayDimension`: the pole block of one index is `dimension`
/// long and the result is derivative-major, `(LocalRequest + 1) * dimension`
/// values.
pub fn eval_dim(
    parameter: f64,
    is_periodic: bool,
    derivative_request: i32,
    extrap_mode: [i32; 2],
    degree: i32,
    flat_knots: &[f64],
    dimension: usize,
    poles: &[f64],
) -> Option<Vec<f64>> {
    let plan = prepare(
        parameter,
        is_periodic,
        derivative_request,
        extrap_mode,
        degree,
        flat_knots,
    )?;
    let dim = dimension.max(1);
    let order = degree + 1;
    let n_orders = (plan.local_request + 1) as usize;
    let mut results = vec![0.0; n_orders * dim];
    let span_sum = |basis_row: &[f64]| -> Vec<f64> {
        let mut index1 = plan.first_non_zero;
        let mut value = vec![0.0; dim];
        for &b in basis_row.iter().take(order as usize) {
            let base = index1.max(1) as usize - 1;
            for (k, slot) in value.iter_mut().enumerate() {
                *slot += poles.get(base * dim + k).copied().unwrap_or(0.0) * b;
            }
            index1 = index1 % plan.modulus + 1;
        }
        value
    };
    if !plan.extrapolating[0] && !plan.extrapolating[1] {
        for ii in 0..n_orders {
            let value = span_sum(&plan.basis[ii]);
            results[ii * dim..(ii + 1) * dim].copy_from_slice(&value);
        }
    } else {
        let mut new_request = derivative_request;
        if new_request > degree {
            new_request = degree;
        }
        let mut taylor = vec![0.0; n_orders * dim];
        let mut inverse = 1.0;
        for ii in 1..=n_orders {
            let value = span_sum(&plan.basis[ii - 1]);
            for k in 0..dim {
                taylor[(ii - 1) * dim + k] = value[k] * inverse;
            }
            inverse /= ii as f64;
        }
        eval_polynomial(plan.delta, new_request, degree, dim, &taylor, &mut results);
    }
    Some(results)
}
