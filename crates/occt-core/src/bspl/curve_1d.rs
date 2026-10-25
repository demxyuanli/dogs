//! 1D B-spline evaluation: `BSplCLib::BuildEval` / `PrepareEval` /
//! `D0` / `D1` / `D2` / `D3` / `DN` for the `Dimension == 1` overloads.
//!
//! Source: `BSplCLib_2.cxx:47-90` (`BuildEval`, the 1D specialization),
//! `:91-139` (`PrepareEval`), `:142-166` (`D0`), `:171-197` (`D1`),
//! `:201-236` (`D2`), `:240-285` (`D3`), `:288-320` (`DN`), and
//! `BSplCLib.cxx:865-900` (`BSplCLib::Eval`, the `Dimension == 1` arm).
//!
//! Translation notes:
//! - OCCT indexes every `NCollection_Array1` 1-based; the port keeps 0-based
//!   slices and rebases the 1-based arithmetic where it matters.
//! - The homogeneous pole buffer of the rational arm is `(value, weight)`
//!   pairs, so `Dimension` is 2 there and 1 otherwise (`cxx:128-138`).
//! - `PLib::RationalDerivative` is called with `Dimension == 1`
//!   (`cxx:189`, `:222`, `:258`, `:307`), which is the *generic* arm of
//!   `PLib.cxx:437-557` (not the `Dimension == 3` arm).

use crate::bspl::bohm::bohm;
use crate::bspl::build_knots::build_knots;
use crate::bspl::knots::pole_index;
use crate::bspl::locate::locate_parameter;
use crate::bspl::plib_rational::rational_derivative_1d;

/// `BSplCLib::BuildEval` for 1D (`BSplCLib_2.cxx:47-90`).
///
/// The local pole buffer is `Degree + 1` values in the non-rational arm and
/// `2 * (Degree + 1)` values `(value, weight)` in the rational arm.
fn build_eval(degree: i32, index: i32, poles: &[f64], weights: Option<&[f64]>) -> Vec<f64> {
    let p_lower = 1i32;
    let p_upper = poles.len() as i32;
    let mut ip = p_lower + index - 1;
    let mut lp = Vec::new();
    match weights {
        None => {
            for _ in 0..=degree {
                ip += 1;
                if ip > p_upper {
                    ip = p_lower;
                }
                lp.push(poles[(ip - 1) as usize]);
            }
        }
        Some(wts) => {
            for _ in 0..=degree {
                ip += 1;
                if ip > p_upper {
                    ip = p_lower;
                }
                let w = wts[(ip - 1) as usize];
                lp.push(poles[(ip - 1) as usize] * w);
                lp.push(w);
            }
        }
    }
    lp
}

/// `BSplCLib::IsRational(Weights, I1, I2)` (`BSplCLib.cxx:842-861`).
///
/// OCCT walks `i` from `I1 - f` to `I3 - 1` with `f == Weights.Lower() == 1`
/// and `I3 == I2 - f`, comparing `Weights(1 + (i % l))` against
/// `Weights(1 + ((i + 1) % l))`. In 0-based terms that is the pair
/// `weights[i % l]` / `weights[(i + 1) % l]`.
fn is_rational_span(weights: &[f64], i1: i32, i2: i32) -> bool {
    let l = weights.len() as i32;
    if l == 0 {
        return false;
    }
    let i3 = i2 - 1; // `I2 - f` with `f == 1`.
    let mut i = i1 - 1;
    while i < i3 {
        let a = weights[(i.rem_euclid(l)) as usize];
        let b = weights[((i + 1).rem_euclid(l)) as usize];
        if a != b {
            return true;
        }
        i += 1;
    }
    false
}

/// The 1D `PrepareEval` result (`BSplCLib_2.cxx:91-139`).
pub struct PreparedEval1D {
    /// The parameter after `LocateParameter` (period-wrapped when periodic).
    pub u: f64,
    /// 1 for the polynomial arm, 2 for the rational arm.
    pub dim: i32,
    pub rational: bool,
    /// The local knot window (`BuildKnots`), length `2 * Degree`.
    pub knots: Vec<f64>,
    /// The local homogeneous pole buffer.
    pub poles: Vec<f64>,
}

/// `BSplCLib::PrepareEval` for 1D (`BSplCLib_2.cxx:91-139`).
pub fn prepare_eval_1d(
    u: f64,
    index: i32,
    degree: i32,
    periodic: bool,
    poles: &[f64],
    weights: Option<&[f64]>,
    knots: &[f64],
    mults: Option<&[i32]>,
) -> PreparedEval1D {
    let (mut knot_index, uu) = locate_parameter(degree, knots, mults, u, periodic, index);
    let local_knots = build_knots(degree, knot_index, periodic, knots, mults);
    if mults.is_none() {
        knot_index -= 1 + degree;
    } else if let Some(m) = mults {
        knot_index = pole_index(degree, knot_index, periodic, m);
    }
    let mut rational = weights.is_some();
    if let Some(w) = weights {
        let w_lower = 1 + knot_index;
        rational = is_rational_span(w, w_lower, w_lower + degree);
    }
    let (dim, lp) = if rational {
        (2, build_eval(degree, knot_index, poles, weights))
    } else {
        (1, build_eval(degree, knot_index, poles, None))
    };
    PreparedEval1D {
        u: uu,
        dim,
        rational,
        knots: local_knots,
        poles: lp,
    }
}

/// `BSplCLib::Eval` for `Dimension == 1` (`BSplCLib.cxx:878-899`).
///
/// The de Boor step over the local pole buffer; the port uses the shared
/// [`crate::bspl::bohm::eval`] whose `Dimension == 1` walk is the same
/// pointer arithmetic.
fn eval_1d(u: f64, degree: i32, knots: &[f64], poles: &mut [f64]) {
    crate::bspl::bohm::eval(u, degree, knots, 1, poles);
}

/// `BSplCLib::D0` for 1D (`BSplCLib_2.cxx:142-166`).
pub fn curve_1d_d0(
    u: f64,
    index: i32,
    degree: i32,
    periodic: bool,
    poles: &[f64],
    weights: Option<&[f64]>,
    knots: &[f64],
    mults: Option<&[i32]>,
) -> f64 {
    let mut p = prepare_eval_1d(u, index, degree, periodic, poles, weights, knots, mults);
    eval_1d(p.u, degree, &p.knots, &mut p.poles);
    if p.rational {
        p.poles[0] / p.poles[1]
    } else {
        p.poles[0]
    }
}

/// `BSplCLib::D1` for 1D (`BSplCLib_2.cxx:171-197`).
pub fn curve_1d_d1(
    u: f64,
    index: i32,
    degree: i32,
    periodic: bool,
    poles: &[f64],
    weights: Option<&[f64]>,
    knots: &[f64],
    mults: Option<&[i32]>,
) -> (f64, f64) {
    let mut p = prepare_eval_1d(u, index, degree, periodic, poles, weights, knots, mults);
    bohm(p.u, degree, 1, &p.knots, p.dim, &mut p.poles);
    if p.rational {
        let ders = rational_derivative_1d(degree, 1, &p.poles, true);
        (ders[0], ders[1])
    } else {
        (p.poles[0], p.poles[1])
    }
}

/// `BSplCLib::D2` for 1D (`BSplCLib_2.cxx:201-236`).
pub fn curve_1d_d2(
    u: f64,
    index: i32,
    degree: i32,
    periodic: bool,
    poles: &[f64],
    weights: Option<&[f64]>,
    knots: &[f64],
    mults: Option<&[i32]>,
) -> (f64, f64, f64) {
    let mut p = prepare_eval_1d(u, index, degree, periodic, poles, weights, knots, mults);
    bohm(p.u, degree, 2, &p.knots, p.dim, &mut p.poles);
    let (value, d1, d2) = if p.rational {
        let ders = rational_derivative_1d(degree, 2, &p.poles, true);
        (ders[0], ders[1], ders[2])
    } else if degree < 2 {
        // `cxx:229-231`: `V2 = 0.` when the polynomial degree is below 2.
        (p.poles[0], p.poles[1], 0.0)
    } else {
        (p.poles[0], p.poles[1], p.poles[2])
    };
    (value, d1, d2)
}

/// `BSplCLib::D3` for 1D (`BSplCLib_2.cxx:240-285`).
pub fn curve_1d_d3(
    u: f64,
    index: i32,
    degree: i32,
    periodic: bool,
    poles: &[f64],
    weights: Option<&[f64]>,
    knots: &[f64],
    mults: Option<&[i32]>,
) -> (f64, f64, f64, f64) {
    let mut p = prepare_eval_1d(u, index, degree, periodic, poles, weights, knots, mults);
    bohm(p.u, degree, 3, &p.knots, p.dim, &mut p.poles);
    let (value, d1, d2, d3) = if p.rational {
        let ders = rational_derivative_1d(degree, 3, &p.poles, true);
        (ders[0], ders[1], ders[2], ders[3])
    } else {
        let v2 = if degree < 2 { 0.0 } else { p.poles[2] };
        let v3 = if degree < 3 { 0.0 } else { p.poles[3] };
        (p.poles[0], p.poles[1], v2, v3)
    };
    (value, d1, d2, d3)
}

/// `BSplCLib::DN` for 1D (`BSplCLib_2.cxx:288-320`).
pub fn curve_1d_dn(
    u: f64,
    n: i32,
    index: i32,
    degree: i32,
    periodic: bool,
    poles: &[f64],
    weights: Option<&[f64]>,
    knots: &[f64],
    mults: Option<&[i32]>,
) -> f64 {
    let mut p = prepare_eval_1d(u, index, degree, periodic, poles, weights, knots, mults);
    bohm(p.u, degree, n, &p.knots, p.dim, &mut p.poles);
    if p.rational {
        // `BSplCLib::DN` calls `PLib::RationalDerivative` with `All == false`
        // (`BSplCLib_2.cxx:307-310`).
        rational_derivative_1d(degree, n, &p.poles, false)[0]
    } else if n > degree {
        // `cxx:313-316`: `VN = 0.` above the degree.
        0.0
    } else {
        p.poles[n as usize]
    }
}
