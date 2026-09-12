//! Rational (NURBS) surface operations: evaluation and weight utilities.

use crate::gp::{GpPnt, GpVec};
use super::eval;

/// Evaluate a rational tensor-product B-spline surface at (u, v).
///
/// Two-stage scheme: for every v-column evaluate a rational u-curve (yielding
/// intermediate poles and weights), then evaluate a rational v-curve through
/// the intermediate poles using the intermediate weights. Zero-weight spans
/// fall back to the unweighted point via `eval_curve_rational`.
pub fn eval_rational(
    poles: &[GpPnt],
    weights: &[f64],
    n_u: usize,
    n_v: usize,
    ku: &[f64],
    kv: &[f64],
    du: usize,
    dv: usize,
    u: f64,
    v: f64,
) -> GpPnt {
    if n_u == 0 || n_v == 0 || poles.is_empty() || weights.is_empty() {
        return GpPnt::zero();
    }
    let mut temp_poles = Vec::with_capacity(n_v);
    let mut temp_weights = Vec::with_capacity(n_v);
    for j in 0..n_v {
        let col: Vec<GpPnt> = (0..n_u).map(|i| poles[i * n_v + j]).collect();
        let wcol: Vec<f64> = (0..n_u).map(|i| weights[i * n_v + j]).collect();
        temp_poles.push(eval::eval_curve_rational(&col, &wcol, ku, du, u));
        temp_weights.push(eval_weight_curve(&wcol, ku, du, u));
    }
    eval::eval_curve_rational(&temp_poles, &temp_weights, kv, dv, v)
}

/// Rational surface D1 via homogeneous poles `A = w P`.
/// Source: `BSplSLib::D1` + `BSplSLib::RationalDerivative` for `(N,M)=(1,1)`
/// (`BSplSLib.cxx:855-970`, `cxx:87-120`): `S = A/w`,
/// `dS = (dA * w - A * dw) / w^2`.
pub fn eval_surface_rational_d1(
    poles: &[GpPnt],
    weights: &[f64],
    n_u: usize,
    n_v: usize,
    knots_u: &[f64],
    knots_v: &[f64],
    degree_u: usize,
    degree_v: usize,
    u: f64,
    v: f64,
) -> (GpPnt, GpVec, GpVec) {
    if poles.len() != weights.len() || poles.is_empty() || n_u == 0 || n_v == 0 {
        return (GpPnt::zero(), GpVec::zero(), GpVec::zero());
    }
    let a_poles: Vec<GpPnt> = poles
        .iter()
        .zip(weights.iter())
        .map(|(p, w)| GpPnt::new(p.x() * w, p.y() * w, p.z() * w))
        .collect();
    let w_poles: Vec<GpPnt> = weights.iter().map(|&w| GpPnt::new(w, 0.0, 0.0)).collect();
    let (a, da_u, da_v) =
        eval::eval_surface_d1(&a_poles, n_u, n_v, knots_u, knots_v, degree_u, degree_v, u, v);
    let (wpt, dw_u, dw_v) =
        eval::eval_surface_d1(&w_poles, n_u, n_v, knots_u, knots_v, degree_u, degree_v, u, v);
    let w = wpt.x();
    if w.abs() < 1e-30 {
        return (a, da_u, da_v);
    }
    let ww = w * w;
    let p = GpPnt::new(a.x() / w, a.y() / w, a.z() / w);
    let du = GpVec::new(
        (da_u.x() * w - a.x() * dw_u.x()) / ww,
        (da_u.y() * w - a.y() * dw_u.x()) / ww,
        (da_u.z() * w - a.z() * dw_u.x()) / ww,
    );
    let dv = GpVec::new(
        (da_v.x() * w - a.x() * dw_v.x()) / ww,
        (da_v.y() * w - a.y() * dw_v.x()) / ww,
        (da_v.z() * w - a.z() * dw_v.x()) / ww,
    );
    (p, du, dv)
}

/// Rational surface D2 via homogeneous poles `A = w P`.
/// Source: `BSplSLib::D2` + `BSplSLib::RationalDerivative` for `(N,M)=(2,2)`
/// (`BSplSLib.cxx:1063-1247`, `cxx:87-120`):
/// `S = A/w`, `dS = (dA * w - A * dw) / w^2`,
/// `Suu = (Auu - 2 Su wu - S wuu) / w`,
/// `Svv = (Avv - 2 Sv wv - S wvv) / w`,
/// `Suv = (Auv - Su wv - Sv wu - S wuv) / w`.
pub fn eval_surface_rational_d2(
    poles: &[GpPnt],
    weights: &[f64],
    n_u: usize,
    n_v: usize,
    knots_u: &[f64],
    knots_v: &[f64],
    degree_u: usize,
    degree_v: usize,
    u: f64,
    v: f64,
) -> (GpPnt, GpVec, GpVec, GpVec, GpVec, GpVec) {
    if poles.len() != weights.len() || poles.is_empty() || n_u == 0 || n_v == 0 {
        return (
            GpPnt::zero(),
            GpVec::zero(),
            GpVec::zero(),
            GpVec::zero(),
            GpVec::zero(),
            GpVec::zero(),
        );
    }
    let a_poles: Vec<GpPnt> = poles
        .iter()
        .zip(weights.iter())
        .map(|(p, w)| GpPnt::new(p.x() * w, p.y() * w, p.z() * w))
        .collect();
    let w_poles: Vec<GpPnt> = weights.iter().map(|&w| GpPnt::new(w, 0.0, 0.0)).collect();
    let (a, da_u, da_v, d2a_u, d2a_v, d2a_uv) = eval::eval_surface_d2(
        &a_poles, n_u, n_v, knots_u, knots_v, degree_u, degree_v, u, v,
    );
    let (wpt, dw_u, dw_v, d2w_u, d2w_v, d2w_uv) = eval::eval_surface_d2(
        &w_poles, n_u, n_v, knots_u, knots_v, degree_u, degree_v, u, v,
    );
    let w = wpt.x();
    if w.abs() < 1e-30 {
        return (a, da_u, da_v, d2a_u, d2a_v, d2a_uv);
    }
    let ww = w * w;
    let p = GpPnt::new(a.x() / w, a.y() / w, a.z() / w);
    let wu = dw_u.x();
    let wv = dw_v.x();
    let wuu = d2w_u.x();
    let wvv = d2w_v.x();
    let wuv = d2w_uv.x();
    let du = GpVec::new(
        (da_u.x() * w - a.x() * wu) / ww,
        (da_u.y() * w - a.y() * wu) / ww,
        (da_u.z() * w - a.z() * wu) / ww,
    );
    let dv = GpVec::new(
        (da_v.x() * w - a.x() * wv) / ww,
        (da_v.y() * w - a.y() * wv) / ww,
        (da_v.z() * w - a.z() * wv) / ww,
    );
    let d2u = GpVec::new(
        (d2a_u.x() - 2.0 * du.x() * wu - p.x() * wuu) / w,
        (d2a_u.y() - 2.0 * du.y() * wu - p.y() * wuu) / w,
        (d2a_u.z() - 2.0 * du.z() * wu - p.z() * wuu) / w,
    );
    let d2v = GpVec::new(
        (d2a_v.x() - 2.0 * dv.x() * wv - p.x() * wvv) / w,
        (d2a_v.y() - 2.0 * dv.y() * wv - p.y() * wvv) / w,
        (d2a_v.z() - 2.0 * dv.z() * wv - p.z() * wvv) / w,
    );
    let d2uv = GpVec::new(
        (d2a_uv.x() - du.x() * wv - dv.x() * wu - p.x() * wuv) / w,
        (d2a_uv.y() - du.y() * wv - dv.y() * wu - p.y() * wuv) / w,
        (d2a_uv.z() - du.z() * wv - dv.z() * wu - p.z() * wuv) / w,
    );
    (p, du, dv, d2u, d2v, d2uv)
}

/// Scalar B-spline evaluation of a weight vector (the intermediate weight of
/// the rational-surface scheme). Reuses the curve evaluator on x = weight
/// points.
fn eval_weight_curve(weights: &[f64], knots: &[f64], degree: usize, u: f64) -> f64 {
    let pts: Vec<GpPnt> = weights.iter().map(|&w| GpPnt::new(w, 0.0, 0.0)).collect();
    eval::eval_curve(&pts, knots, degree, u).x()
}

/// True when every weight is positive enough to be usable (> 1e-12).
pub fn weights_grid_ok(weights: &[f64]) -> bool {
    weights.iter().all(|&w| w > 1e-12)
}

/// Normalize every u-row of the weight grid by its row maximum.
pub fn normalize_weights_grid(weights: &mut [f64], n_u: usize, n_v: usize) {
    for i in 0..n_u {
        let row_max = weights[i * n_v..(i + 1) * n_v].iter().cloned().fold(f64::MIN, f64::max);
        if row_max > 1e-12 {
            for w in &mut weights[i * n_v..(i + 1) * n_v] {
                *w /= row_max;
            }
        }
    }
}

/// Bounds check for a pole index in an n_poles_u x n_poles_v grid.
pub fn surface_contains_pole(n_poles_u: usize, n_poles_v: usize, i_u: usize, i_v: usize) -> bool {
    i_u < n_poles_u && i_v < n_poles_v
}

/// Clamp every weight to the usable range [1e-12, 1e12].
pub fn project_weights(weights: &[f64]) -> Vec<f64> {
    weights.iter().map(|&w| w.clamp(1e-12, 1e12)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn weights_grid_ok_accepts_positive() {
        assert!(weights_grid_ok(&[1.0, 2.0, 0.5]));
        assert!(!weights_grid_ok(&[1.0, 0.0, 0.5]));
        assert!(!weights_grid_ok(&[1e-13, 2.0]));
    }

    #[test]
    fn normalize_weights_grid_rows_max_to_one() {
        let mut w = vec![0.5, 2.0, 1.0, 4.0];
        normalize_weights_grid(&mut w, 2, 2);
        assert!((w[0] - 0.25).abs() < 1e-12);
        assert!((w[1] - 1.0).abs() < 1e-12);
        assert!((w[2] - 0.25).abs() < 1e-12);
        assert!((w[3] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn surface_contains_pole_bounds() {
        assert!(surface_contains_pole(4, 5, 3, 4));
        assert!(!surface_contains_pole(4, 5, 4, 0));
        assert!(!surface_contains_pole(4, 5, 0, 5));
    }

    #[test]
    fn eval_rational_bilinear() {
        let poles = vec![
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(0.0, 1.0, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
            GpPnt::new(1.0, 1.0, 0.0),
        ];
        let weights = vec![1.0; 4];
        let ku = vec![0.0, 0.0, 1.0, 1.0];
        let kv = vec![0.0, 0.0, 1.0, 1.0];
        let p = eval_rational(&poles, &weights, 2, 2, &ku, &kv, 1, 1, 0.5, 0.5);
        assert!((p.x() - 0.5).abs() < 1e-12);
        assert!((p.y() - 0.5).abs() < 1e-12);
    }
}
