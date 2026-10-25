//! `BSplCLib::MovePointAndTangent` (`BSplCLib_2.cxx:567-857`) for a 1D curve.
//!
//! Builds the two cubic blending functions that vanish on the first/last
//! `StartingCondition + 2` / `EndingCondition + 2` Schoenberg points, evaluates
//! their value and first derivative at `U`, inverts the 2x2 system to obtain
//! the pair of coefficients that moves the curve by `Delta` with derivative
//! `DeltaDerivatives`, and rewrites every pole as
//! `Pole + f1 * a + f2 * b` (`cxx:834-846`).

use crate::bspl::banded_interp::schoenberg_points;
use crate::bspl::eval_extrap::{eval_1d, eval_1d_weighted};
use crate::bspl::locate::locate_parameter_range;

/// Result of [`move_point_and_tangent_1d`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MovePointError {
    /// OCCT's `ErrorStatus` out parameter: `1` = not enough freedom left at the
    /// ends for the requested conditions, `2` = the input conditions/degree are
    /// out of range (`BSplCLib_2.cxx:848-856`).
    Status(i32),
    /// `gp_Mat2d::Invert` raised `Standard_ConstructionError`
    /// (`gp_Mat2d.cxx:132-134`): the blending functions are linearly dependent
    /// at `U`.
    SingularMatrix,
}

/// `BSplCLib::MovePointAndTangent` (`BSplCLib_2.cxx:567-857`), 1D.
///
/// `poles` holds the curve's own poles; `weights` is the (rational) weight
/// array of the same length, used only for the homogenised evaluation
/// (`cxx:769-805`). On success the moved pole array is returned; on
/// `Err(Status(_))` OCCT leaves `NewPoles` untouched.
pub fn move_point_and_tangent_1d(
    u: f64,
    delta: f64,
    delta_derivative: f64,
    tolerance: f64,
    degree: i32,
    starting_condition: i32,
    ending_condition: i32,
    poles: &[f64],
    weights: Option<&[f64]>,
    flat_knots: &[f64],
) -> Result<Vec<f64>, MovePointError> {
    let order = degree + 1;
    let num_knots = flat_knots.len() as i32;
    let num_poles = num_knots - order;
    let conditions = starting_condition + ending_condition + 4;
    if !(starting_condition >= -1
        && starting_condition <= degree
        && ending_condition >= -1
        && ending_condition <= degree
        && conditions <= num_poles)
    {
        return Err(MovePointError::Status(2));
    }
    // `cxx:611-635`: the parameter must leave room for the requested end
    // conditions.
    let start_index = 1 + degree;
    let end_index = num_knots - degree;
    let fk = |i: i32| flat_knots[(i - 1) as usize];
    let mut free = true;
    if starting_condition == -1 {
        free = free && fk(start_index) <= u;
    } else {
        free = free && fk(start_index) + tolerance < u;
    }
    if ending_condition == -1 {
        free = free && fk(end_index) >= u;
    } else {
        free = free && fk(end_index) - tolerance > u;
    }
    if !free {
        return Err(MovePointError::Status(1));
    }
    // Rust-only bounds guard: the pole loop below (`BSplCLib_2.cxx:836`) is
    // bounded by `FlatKnots.Length() - Order`, which for a periodic curve is
    // `poles.len() + Degree`; OCCT would then read and write past the end of
    // `Poles`/`NewPoles`.
    if poles.len() as i32 != num_poles {
        return Err(MovePointError::Status(2));
    }

    // `cxx:639-658`.
    let sp = schoenberg_points(degree as usize, flat_knots, num_poles as usize);
    let sp_at = |i: i32| sp[(i - 1) as usize];
    let start_num_poles = starting_condition + 2;
    let end_num_poles = num_poles - ending_condition - 1;
    let (index, _new_parameter) = locate_parameter_range(
        &sp,
        u,
        false,
        start_num_poles,
        end_num_poles,
        0.0,
        1.0,
    );
    let other_index = if index == start_num_poles {
        index + 1
    } else if index == end_num_poles {
        index - 1
    } else if u - fk(index) < fk(index + 1) - u {
        index - 1
    } else {
        index + 1
    };
    let power = 3;

    // `cxx:673-691`.
    let start_value = if start_num_poles == 1 {
        sp_at(1) - (sp_at(num_poles) - sp_at(1))
    } else {
        sp_at(start_num_poles - 1)
    };
    let end_value = if end_num_poles == num_poles {
        sp_at(num_poles) + (sp_at(num_poles) - sp_at(1))
    } else {
        sp_at(end_num_poles + 1)
    };

    let n = num_poles as usize;
    let mut first = vec![0.0; n];
    let mut second = vec![0.0; n];
    // Arrays are 1-based in OCCT; `*_at(i)` below takes the same 1-based index.
    let mut divide = 1.0 / (sp_at(index) - start_value);
    for ii in start_num_poles..=index {
        let value = (sp_at(ii) - start_value) * divide;
        first[(ii - 1) as usize] = value.powi(power);
    }
    divide = 1.0 / (end_value - sp_at(index));
    for ii in index..=end_num_poles {
        let value = (end_value - sp_at(ii)) * divide;
        first[(ii - 1) as usize] = value.powi(power);
    }
    divide = 1.0 / (sp_at(other_index) - start_value);
    for ii in start_num_poles..=other_index {
        let value = (sp_at(ii) - start_value) * divide;
        second[(ii - 1) as usize] = value.powi(power);
    }
    divide = 1.0 / (end_value - sp_at(other_index));
    for ii in other_index..=end_num_poles {
        let value = (end_value - sp_at(ii)) * divide;
        second[(ii - 1) as usize] = value.powi(power);
    }

    // `cxx:761-810`: value and first derivative of both blending functions.
    let extrap_mode = [degree, degree];
    let mut results = [[0.0f64; 2]; 2];
    for (k, func) in [&first, &second].into_iter().enumerate() {
        if let Some(w) = weights {
            let (p, wd) = eval_1d_weighted(u, false, 1, extrap_mode, degree, flat_knots, func, w)
                .ok_or(MovePointError::Status(2))?;
            let rational = crate::bspl::plib_rational::rational_derivatives_split(1, 1, &p, &wd);
            results[k][0] = rational[0];
            results[k][1] = rational[1];
        } else {
            let r = eval_1d(u, false, 1, extrap_mode, degree, flat_knots, func)
                .ok_or(MovePointError::Status(2))?;
            results[k][0] = r[0];
            results[k][1] = r[1];
        }
    }

    // `gp_Mat2d a_matrix` built as `SetValue(ii + 1, jj + 1, results[ii][jj])`,
    // then `a_matrix.Invert()` (`cxx:811-821`, `gp_Mat2d.cxx:125-140`).
    let m = [
        [results[0][0], results[0][1]],
        [results[1][0], results[1][1]],
    ];
    let new_mat = [[m[1][1], -m[0][1]], [-m[1][0], m[0][0]]];
    let det = new_mat[0][0] * new_mat[1][1] - new_mat[0][1] * new_mat[1][0];
    if det.abs() <= crate::precision::REAL_SMALL {
        return Err(MovePointError::SingularMatrix);
    }
    let inverse = 1.0 / det;
    let inv = [
        [new_mat[0][0] * inverse, new_mat[1][0] * inverse],
        [new_mat[0][1] * inverse, new_mat[1][1] * inverse],
    ];
    // `cxx:824-832`, `ArrayDimension == 1`.
    let the_a = inv[0][0] * delta + inv[1][0] * delta_derivative;
    let the_b = inv[0][1] * delta + inv[1][1] * delta_derivative;

    let mut new_poles = Vec::with_capacity(n);
    for ii in 0..n {
        let mut value = poles[ii];
        value += first[ii] * the_a;
        value += second[ii] * the_b;
        new_poles.push(value);
    }
    Ok(new_poles)
}
