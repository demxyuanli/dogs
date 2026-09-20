//! `BSplCLib::EvalBsplineBasis` — B-spline basis values (and derivatives) at a
//! parameter, in OCCT's own layout.
//!
//! Source: `BSplCLib_2.cxx:429-563` (the function lives in the `_2` split of
//! `BSplCLib`). **`BSplCLib::BasisFuns` does not exist in OCCT 8.0.0** (whole
//! tree grep: no hit in `BSplCLib.cxx`/`.hxx`/`_2.cxx`); this is the routine the
//! 8.0.0 `BSplCLib::Eval` calls (`BSplCLib.cxx:3526`) and the one the audit's
//! A15/T-51 item mistook for `BasisFuns`.
//!
//! The result matrix is OCCT's: `DerivativeRequest + 1` rows × `Order` columns,
//! row 1 the values `B_i(t) … B_{i+k-1}(t)`, row `r+1` the `r`-th derivative,
//! with `FirstNonZeroBsplineIndex` (`i`, 1-based) saying where the block starts.
//! `ErrorCode` follows OCCT: `1` = matrix too small (unreachable here, the port
//! allocates the same `(DerivativeRequest+1) × Order` matrix as
//! `BSplCLib_LocalMatrix`, so the guard of `cxx:466-470` always passes), `2` =
//! vanishing knot span (`cxx:509-512`/`:538-541`, compared against
//! **`gp::Resolution()`**, i.e. [`REAL_SMALL`] — *not* an invented `1e-15`).

use super::locate::locate_parameter_range;
use crate::precision::REAL_SMALL;

/// `BSplCLib::EvalBsplineBasis` (`BSplCLib_2.cxx:429-563`).
///
/// `flat_knots` is the flattened knot vector (0-based slice holding the same
/// values as OCCT's 1-based `NCollection_Array1`), `order = degree + 1`.
/// Returns `(FirstNonZeroBsplineIndex, basis)` with `basis[r][c]` = OCCT
/// `BsplineBasis(r + 1, c + 1)`, or the OCCT error code.
pub fn eval_bspline_basis(
    derivative_request: i32,
    order: i32,
    flat_knots: &[f64],
    parameter: f64,
    is_periodic: bool,
) -> Result<(i32, Vec<Vec<f64>>), i32> {
    let mut local_request = derivative_request;
    if derivative_request >= order {
        local_request = order - 1;
    }
    if order <= 0 || local_request < 0 {
        return Err(1);
    }

    let n_cols = order as usize;
    let n_rows = (derivative_request + 1) as usize;
    // `BSplCLib_LocalMatrix BsplineBasis(LocalRequest, Order)` (`cxx:3525`):
    // rows 1..=DerivativeRequest+1, cols 1..=Order, row-major (`aPos =
    // (row - LowerRow) * nCols + (col - LowerCol)`, `NCollection_Array2.hxx:317`).
    let n_poles = flat_knots.len() as i32 - order;

    // `BSplCLib::LocateParameter(Degree, Knots, U, Periodic, FromK1, ToK2, …)`
    // (`BSplCLib.cxx:189-214`): `Degree = Order - 1`, `FromK1 = Order`,
    // `ToK2 = aNumPoles + 1`; periodic uses the trimmed knot window as the
    // period bounds, non-periodic passes `(0., 1.)`.
    let (u_first, u_last) = if is_periodic {
        (
            flat_knots[(order - 1) as usize],
            flat_knots[flat_knots.len() - order as usize],
        )
    } else {
        (0.0, 1.0)
    };
    let (knot_index, new_param) = locate_parameter_range(
        flat_knots,
        parameter,
        is_periodic,
        order,
        n_poles + 1,
        u_first,
        u_last,
    );

    let first_non_zero = knot_index - order + 1;

    // Raw pointer view of OCCT, translated: `ii` is 1-based into `flat_knots`,
    // the port's slice is 0-based, so rebase by `FlatKnots.Lower() == 1`.
    let ii = (knot_index - 1) as usize;
    let knots = |k: i32| flat_knots[(ii as i32 + k) as usize];

    let mut data = vec![0.0f64; n_rows * n_cols];
    data[0] = 1.0;
    let a_local_request = local_request;

    for qq in 2..=order - a_local_request {
        data[(qq - 1) as usize] = 0.0;
        for pp in 1..=qq - 1 {
            let scale = knots(pp) - knots(1 - qq + pp);
            if scale.abs() < REAL_SMALL {
                return Err(2);
            }
            let factor = (new_param - knots(1 - qq + pp)) / scale;
            let saved = factor * data[(pp - 1) as usize];
            data[(pp - 1) as usize] *= 1.0 - factor;
            data[(pp - 1) as usize] += data[(qq - 1) as usize];
            data[(qq - 1) as usize] = saved;
        }
    }

    for qq in order - a_local_request + 1..=order {
        for pp in 1..=qq - 1 {
            data[((order - qq + 1) * order + (pp - 1)) as usize] = data[(pp - 1) as usize];
        }
        data[(qq - 1) as usize] = 0.0;

        for ss in order - a_local_request + 1..=qq {
            data[((order - ss + 1) * order + (qq - 1)) as usize] = 0.0;
        }

        for pp in 1..=qq - 1 {
            let scale = knots(pp) - knots(1 - qq + pp);
            if scale.abs() < REAL_SMALL {
                return Err(2);
            }
            let inverse = 1.0 / scale;
            let factor = (new_param - knots(1 - qq + pp)) * inverse;
            let mut saved = factor * data[(pp - 1) as usize];
            data[(pp - 1) as usize] *= 1.0 - factor;
            data[(pp - 1) as usize] += data[(qq - 1) as usize];
            data[(qq - 1) as usize] = saved;
            let local_inverse = f64::from(qq - 1) * inverse;

            for ss in order - a_local_request + 1..=qq {
                let row = ((order - ss + 1) * order) as usize;
                saved = local_inverse * data[row + (pp - 1) as usize];
                data[row + (pp - 1) as usize] *= -local_inverse;
                data[row + (pp - 1) as usize] += data[row + (qq - 1) as usize];
                data[row + (qq - 1) as usize] = saved;
            }
        }
    }

    let basis: Vec<Vec<f64>> = (0..n_rows)
        .map(|r| data[r * n_cols..(r + 1) * n_cols].to_vec())
        .collect();
    Ok((first_non_zero, basis))
}
