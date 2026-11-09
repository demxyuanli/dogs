//! `BSplCLib::MergeBSplineKnots`, `FunctionMultiply`, and the non-extrapolating
//! arm of `BSplCLib::Eval` used by them.
//! Source: `BSplCLib_2.cxx:861-934`, `:1001-1028`, `:1063-1173`;
//! `BSplCLib.cxx:3640-3743` (in-range path only).

use crate::bspl::banded_interp::{eval_bspline_basis_flat, interpolate_contact, knot_sequence, schoenberg_points};
use crate::bspl::knots::{nb_poles, reparametrize};

/// `BSplCLib::MergeBSplineKnots` (`BSplCLib_2.cxx:1063-1173`).
///
/// Returns `(num_poles, new_knots, new_mults)`.
pub fn merge_bspline_knots(
    tolerance: f64,
    start_value: f64,
    end_value: f64,
    degree1: i32,
    knots1: &[f64],
    mults1: &[i32],
    degree2: i32,
    knots2: &[f64],
    mults2: &[i32],
) -> (i32, Vec<f64>, Vec<i32>) {
    if start_value < end_value - tolerance {
        let mut k1 = knots1.to_vec();
        let mut k2 = knots2.to_vec();
        let degree = degree1 + degree2;
        reparametrize(start_value, end_value, &mut k1);
        reparametrize(start_value, end_value, &mut k2);

        let mut num_knots = 0i32;
        let mut jj = 0usize;
        for ii in 0..k1.len() {
            while jj < k2.len() && k2[jj] <= k1[ii] - tolerance {
                jj += 1;
                num_knots += 1;
            }
            while jj < k2.len() && k2[jj] <= k1[ii] + tolerance {
                jj += 1;
            }
            num_knots += 1;
        }
        let mut new_knots = vec![0.0; num_knots as usize];
        let mut new_mults = vec![0i32; num_knots as usize];
        let mut nk = 0usize;
        jj = 0;
        for ii in 0..k1.len() {
            while jj < k2.len() && k2[jj] <= k1[ii] - tolerance {
                new_knots[nk] = k2[jj];
                new_mults[nk] = mults2[jj] + degree1;
                jj += 1;
                nk += 1;
            }
            let mut set_mults_flag = false;
            while jj < k2.len() && k2[jj] <= k1[ii] + tolerance {
                let continuity = (degree1 - mults1[ii]).min(degree2 - mults2[jj]);
                set_mults_flag = true;
                new_mults[nk] = degree - continuity;
                jj += 1;
            }
            new_knots[nk] = k1[ii];
            if !set_mults_flag {
                new_mults[nk] = mults1[ii] + degree2;
            }
            nk += 1;
        }
        let num_knots = nk;
        new_knots.truncate(num_knots);
        new_mults.truncate(num_knots);
        if !new_mults.is_empty() {
            new_mults[0] = degree + 1;
            new_mults[num_knots - 1] = degree + 1;
        }
        let mut index = 0i32;
        for m in &new_mults {
            index += *m;
        }
        let num_poles = index - degree - 1;
        (num_poles, new_knots, new_mults)
    } else {
        let degree = degree1 + degree2;
        let new_knots = vec![start_value, end_value];
        let new_mults = vec![degree + 1, degree + 1];
        let num_poles = nb_poles(degree, false, &new_mults);
        (num_poles, new_knots, new_mults)
    }
}

/// Non-extrapolating value evaluation of a flat-knot B-spline
/// (`BSplCLib.cxx:3718-3743` when both extrapolating flags are clear).
fn eval_poles_at(u: f64, degree: i32, flat: &[f64], poles: &[f64], dim: usize, out: &mut [f64]) {
    let order = (degree + 1) as usize;
    if order == 0 || flat.len() < order || out.len() < dim {
        for o in out.iter_mut().take(dim) {
            *o = 0.0;
        }
        return;
    }
    let n_poles = flat.len() - order;
    let Ok((first_nz, basis)) = eval_bspline_basis_flat(0, order, flat, u) else {
        for o in out.iter_mut().take(dim) {
            *o = 0.0;
        }
        return;
    };
    for k in 0..dim {
        out[k] = 0.0;
    }
    for j in 0..order {
        let pi = first_nz + j;
        if pi >= n_poles {
            break;
        }
        let b = basis[j];
        for k in 0..dim {
            out[k] += poles[pi * dim + k] * b;
        }
    }
}

/// Evaluator callback: `(derivative_request, start_end, parameter) -> Result<f64, ()>`.
/// `start_end` is `[u_first, u_last]` of the target flat knot vector.
pub type EvaluatorFn<'a> = dyn Fn(i32, [f64; 2], f64) -> Result<f64, ()> + 'a;

/// `BSplCLib::FunctionMultiply` dimension form (`BSplCLib_2.cxx:861-934`).
///
/// `poles` / `new_poles` are row-major `[n][dimension]`. Returns `Err(1)` when
/// the evaluator fails (`theStatus = 1`); interpolation pivot failure returns
/// its status code.
pub fn function_multiply(
    function: &EvaluatorFn<'_>,
    bspline_degree: i32,
    bspline_flat_knots: &[f64],
    poles: &[f64],
    dimension: usize,
    flat_knots: &[f64],
    new_degree: i32,
    new_poles: &mut [f64],
) -> Result<(), i32> {
    let num_new_poles = flat_knots.len() as i32 - new_degree - 1;
    if num_new_poles <= 0 || dimension == 0 {
        return Err(1);
    }
    let num_new = num_new_poles as usize;
    let start_end = [
        flat_knots[new_degree as usize],
        flat_knots[num_new],
    ];
    let mut parameters = schoenberg_points(new_degree as usize, flat_knots, num_new);
    if let Some(p0) = parameters.first_mut() {
        if *p0 < start_end[0] {
            *p0 = start_end[0];
        }
    }
    if let Some(pn) = parameters.last_mut() {
        if *pn > start_end[1] {
            *pn = start_end[1];
        }
    }
    let mut samples = vec![0.0; num_new * dimension];
    for (i, &param) in parameters.iter().enumerate() {
        let result = function(0, start_end, param).map_err(|_| 1)?;
        let mut buf = vec![0.0; dimension];
        eval_poles_at(
            param,
            bspline_degree,
            bspline_flat_knots,
            poles,
            dimension,
            &mut buf,
        );
        for k in 0..dimension {
            samples[i * dimension + k] = buf[k] * result;
        }
    }
    let contact = vec![0i32; num_new];
    interpolate_contact(
        new_degree as usize,
        flat_knots,
        &parameters,
        &contact,
        &mut samples,
        dimension,
    )?;
    if new_poles.len() < num_new * dimension {
        return Err(1);
    }
    new_poles[..num_new * dimension].copy_from_slice(&samples);
    Ok(())
}

/// Convenience: build a flat knot sequence from distinct knots/mults.
pub fn flat_knot_sequence(knots: &[f64], mults: &[i32], degree: i32) -> Vec<f64> {
    knot_sequence(knots, mults, degree)
}

/// `BSplCLib::FunctionReparameterise` dimension form (`BSplCLib_2.cxx:938-997`).
///
/// Samples the evaluator at Schoenberg parameters of `flat_knots`, evaluates
/// the source B-spline at those mapped parameters, then interpolates.
pub fn function_reparameterise(
    function: &EvaluatorFn<'_>,
    bspline_degree: i32,
    bspline_flat_knots: &[f64],
    poles: &[f64],
    dimension: usize,
    flat_knots: &[f64],
    new_degree: i32,
    new_poles: &mut [f64],
) -> Result<(), i32> {
    let num_new_poles = flat_knots.len() as i32 - new_degree - 1;
    if num_new_poles <= 0 || dimension == 0 {
        return Err(1);
    }
    let num_new = num_new_poles as usize;
    let start_end = [
        flat_knots[new_degree as usize],
        flat_knots[num_new],
    ];
    let parameters = schoenberg_points(new_degree as usize, flat_knots, num_new);
    let mut samples = vec![0.0; num_new * dimension];
    for (i, &param) in parameters.iter().enumerate() {
        let mapped = function(0, start_end, param).map_err(|_| 1)?;
        let mut buf = vec![0.0; dimension];
        eval_poles_at(
            mapped,
            bspline_degree,
            bspline_flat_knots,
            poles,
            dimension,
            &mut buf,
        );
        samples[i * dimension..(i + 1) * dimension].copy_from_slice(&buf);
    }
    let contact = vec![0i32; num_new];
    interpolate_contact(
        new_degree as usize,
        flat_knots,
        &parameters,
        &contact,
        &mut samples,
        dimension,
    )?;
    if new_poles.len() < num_new * dimension {
        return Err(1);
    }
    new_poles[..num_new * dimension].copy_from_slice(&samples);
    Ok(())
}
