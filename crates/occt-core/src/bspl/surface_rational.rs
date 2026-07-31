//! Rational (NURBS) surface operations: evaluation and weight utilities.

use crate::gp::GpPnt;
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
