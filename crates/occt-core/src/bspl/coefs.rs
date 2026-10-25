//! `BSplCLib::PolesCoefficients` / `CoefsD0` / `CoefsD1` — the "cache" syntax
//! used to evaluate Bezier curves as power-basis polynomials.
//!
//! Sources (OCCT 8.0.0):
//! - `BSplCLib::PolesCoefficients` (`BSplCLib_BzSyntaxes.cxx:62-84`), which is
//!   `BSplCLib_BuildCache(0, 1, false, Deg, FlatBezierKnots(Deg), ...)`
//!   (`BSplCLib_CurveComputation.pxx:2052-2070`),
//! - `BSplCLib_BuildCache` (`BSplCLib_CurveComputation.pxx:1527-1600`),
//! - `BSplCLib::CoefsD0/D1` (`BSplCLib.lxx:47-84`), which forward to
//!   `CacheD0`/`CacheD1` (`BSplCLib_CurveComputation.pxx:1266-1344`) with
//!   `Degree = Length - 1`, `CacheParameter = 0`, `SpanLength = 1`,
//! - `PLib::EvalPolynomial` / `NoDerivativeEvalPolynomial`
//!   (`PLib.cxx:945-1100`).
//!
//! Translation note: `BuildCache` runs `BSplCLib::Bohm(U, Deg, Deg, localKnots,
//! ...)` first, then divides the `k`-th slot by `(k-1)!` (`pxx:1584-1589`,
//! `SpanDomain == 1`); the result is the polynomial's Taylor coefficient
//! vector, which `EvalPolynomial` then evaluates directly.

use crate::gp::{GpPnt, GpPnt2d, GpVec, GpVec2d};

use super::bohm::bohm;
use super::build_knots::build_knots;
use super::locate::locate_parameter;

/// `BSplCLib::FlatBezierKnots(Degree)` (`BSplCLib.cxx:4971-4977`): `Degree + 1`
/// zeros followed by `Degree + 1` ones (`cxx:4940-4966` builds the table).
pub fn flat_bezier_knots(degree: i32) -> Vec<f64> {
    let d = degree.max(1) as usize;
    let mut k = Vec::with_capacity(2 * (d + 1));
    for _ in 0..=d {
        k.push(0.0);
    }
    for _ in 0..=d {
        k.push(1.0);
    }
    k
}

/// Local knot window `BuildCache` hands to `Bohm`: `PrepareEval_T`
/// (`pxx:789-793`) locates `U` in the flat knots and calls `BuildKnots`.
fn cache_knots(u: f64, degree: i32, flat: &[f64]) -> Vec<f64> {
    let (index, _) = locate_parameter(degree, flat, None, u, false, 0);
    build_knots(degree, index, false, flat, None)
}

/// `BSplCLib_BuildCache` (`pxx:1527-1600`) for a non-rational Bezier curve:
/// `poles` -> the `(Degree + 1) x Dimension` coefficient buffer.
fn build_cache(poles: &[f64], degree: i32, dimension: i32) -> Vec<f64> {
    let dim = dimension.max(1) as usize;
    let deg = degree.max(0) as usize;
    if poles.len() < (deg + 1) * dim {
        return vec![0.0; (deg + 1) * dim];
    }
    let flat = flat_bezier_knots(degree);
    let local = cache_knots(0.0, degree, &flat);
    let mut buf = poles[..(deg + 1) * dim].to_vec();
    // `Bohm(u, Degree, Degree, *dc.knots, LocalDimension, *dc.poles)`
    // (`pxx:1557`).
    bohm(0.0, degree, degree, &local, dimension, &mut buf);
    // `Traits::CoordsToPointMultiplied(CachePoles(ii), dc.poles + LocalIndex,
    // LocalValue)` with `LocalValue *= SpanDomain / ii`, `SpanDomain == 1`
    // (`pxx:1584-1589`).
    let mut local_value = 1.0f64;
    for ii in 1..=deg + 1 {
        let base = (ii - 1) * dim;
        for k in 0..dim {
            buf[base + k] *= local_value;
        }
        local_value *= 1.0 / ii as f64;
    }
    buf
}

/// `BSplCLib::PolesCoefficients(Poles, NoWeights, CachePoles, NoWeights)`
/// (`BSplCLib_BzSyntaxes.cxx:62-73`).
pub fn poles_coefficients(poles: &[GpPnt]) -> Vec<GpPnt> {
    let mut flat = Vec::with_capacity(poles.len() * 3);
    for p in poles {
        flat.push(p.x());
        flat.push(p.y());
        flat.push(p.z());
    }
    let degree = poles.len() as i32 - 1;
    let c = build_cache(&flat, degree, 3);
    (0..poles.len())
        .map(|i| GpPnt::new(c[3 * i], c[3 * i + 1], c[3 * i + 2]))
        .collect()
}

/// `BSplCLib::PolesCoefficients` for `gp_Pnt2d`
/// (`BSplCLib_BzSyntaxes.cxx:75-84`).
pub fn poles_coefficients2d(poles: &[GpPnt2d]) -> Vec<GpPnt2d> {
    let mut flat = Vec::with_capacity(poles.len() * 2);
    for p in poles {
        flat.push(p.x());
        flat.push(p.y());
    }
    let degree = poles.len() as i32 - 1;
    let c = build_cache(&flat, degree, 2);
    (0..poles.len())
        .map(|i| GpPnt2d::new(c[2 * i], c[2 * i + 1]))
        .collect()
}

/// `PLib::EvalPolynomial(Par, DerivativeRequest, Degree, Dimension, coeffs,
/// results)` (`PLib.cxx:945-1025`): `results` holds `DerivativeRequest + 1`
/// `Dimension`-long blocks, `[f, f', ...]`.
pub fn eval_polynomial(
    par: f64,
    derivative_request: i32,
    degree: i32,
    dimension: i32,
    coeffs: &[f64],
) -> Vec<f64> {
    let dim = dimension.max(0) as usize;
    let deg = degree.max(0) as usize;
    let req = derivative_request.max(0) as usize;
    let res_size = (1 + req) * dim;
    let mut res = vec![0.0f64; res_size];
    if dim == 0 || coeffs.len() < (deg + 1) * dim {
        return res;
    }
    // `cxx:1003-1023`: walk the coefficients from the highest power down,
    // shifting the derivative chain one block towards the front each round.
    for a_deg in (0..=deg).rev() {
        let mut ptr = res_size - dim;
        for a_deriv in (1..=req).rev() {
            let original = ptr - dim;
            for ind in 0..dim {
                res[ptr + ind] = res[ptr + ind] * par + res[original + ind] * a_deriv as f64;
            }
            ptr = original;
        }
        for ind in 0..dim {
            res[ptr + ind] = res[ptr + ind] * par + coeffs[a_deg * dim + ind];
        }
    }
    res
}

/// `PLib::NoDerivativeEvalPolynomial` (`PLib.cxx:1030-1100`): the value only,
/// same coefficient layout as `eval_polynomial`.
pub fn no_derivative_eval_polynomial(
    par: f64,
    degree: i32,
    dimension: i32,
    coeffs: &[f64],
) -> Vec<f64> {
    let dim = dimension.max(0) as usize;
    let deg = degree.max(0) as usize;
    let mut res = vec![0.0f64; dim];
    if dim == 0 || coeffs.len() < (deg + 1) * dim {
        return res;
    }
    for a_deg in (0..=deg).rev() {
        for ind in 0..dim {
            res[ind] = res[ind] * par + coeffs[a_deg * dim + ind];
        }
    }
    res
}

/// `BSplCLib::CoefsD0(U, Poles, NoWeights, P)` (`BSplCLib.lxx:47-51`).
pub fn coefs_d0(u: f64, coeffs: &[GpPnt]) -> GpPnt {
    let degree = coeffs.len() as i32 - 1;
    let c = flatten3(coeffs);
    let r = no_derivative_eval_polynomial(u, degree, 3, &c);
    GpPnt::new(r[0], r[1], r[2])
}

/// `BSplCLib::CoefsD1(U, Poles, NoWeights, P, Vec)` (`BSplCLib.lxx:67-73`).
///
/// `CacheD1` divides the derivative by `SpanLength` (`pxx:1327`), which is `1`
/// in the `CoefsD1` forwarding call (`lxx:72`).
pub fn coefs_d1(u: f64, coeffs: &[GpPnt]) -> (GpPnt, GpVec) {
    let degree = coeffs.len() as i32 - 1;
    let c = flatten3(coeffs);
    let r = eval_polynomial(u, 1, degree, 3, &c);
    (
        GpPnt::new(r[0], r[1], r[2]),
        GpVec::new(r[3], r[4], r[5]),
    )
}

/// `BSplCLib::CoefsD0` for `gp_Pnt2d` (`BSplCLib.lxx:55-61`).
pub fn coefs_d0_2d(u: f64, coeffs: &[GpPnt2d]) -> GpPnt2d {
    let degree = coeffs.len() as i32 - 1;
    let c = flatten2(coeffs);
    let r = no_derivative_eval_polynomial(u, degree, 2, &c);
    GpPnt2d::new(r[0], r[1])
}

/// `BSplCLib::CoefsD1` for `gp_Pnt2d` (`BSplCLib.lxx:76-84`).
pub fn coefs_d1_2d(u: f64, coeffs: &[GpPnt2d]) -> (GpPnt2d, GpVec2d) {
    let degree = coeffs.len() as i32 - 1;
    let c = flatten2(coeffs);
    let r = eval_polynomial(u, 1, degree, 2, &c);
    (GpPnt2d::new(r[0], r[1]), GpVec2d::new(r[2], r[3]))
}

fn flatten3(poles: &[GpPnt]) -> Vec<f64> {
    let mut out = Vec::with_capacity(poles.len() * 3);
    for p in poles {
        out.push(p.x());
        out.push(p.y());
        out.push(p.z());
    }
    out
}

fn flatten2(poles: &[GpPnt2d]) -> Vec<f64> {
    let mut out = Vec::with_capacity(poles.len() * 2);
    for p in poles {
        out.push(p.x());
        out.push(p.y());
    }
    out
}
