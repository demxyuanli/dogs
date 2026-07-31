//! Bezier extraction from B-spline. Source: `BSplCLib.cxx` — Boehm, FlatBezier
//! Converts B-spline curves to piecewise Bezier representation via knot insertion.

use crate::gp::GpPnt;

/// Flat Bezier coefficients: extract Bezier segments from B-spline curve.
/// Returns vector of (bezier_poles_per_segment, start_parameter, end_parameter).
/// Source: BSplCLib::FlatBezierCoefficients
pub fn flat_bezier_coefficients(poles: &[GpPnt], weights: Option<&[f64]>,
                                 knots: &[f64], degree: usize,
                                 tolerance: f64) -> Vec<Vec<GpPnt>> {
    let mut result = Vec::new();
    let n_segments = count_segments(knots, degree);

    for seg in 0..n_segments {
        let u0 = knots[degree + seg];
        let u1 = knots[degree + seg + 1];
        let mult_start = super::knots::multiplicity(knots, u0);
        let mult_end = super::knots::multiplicity(knots, u1);

        let insertions_needed = (degree + 1 - mult_start).max(degree + 1 - mult_end);
        if insertions_needed == 0 {
            // Already Bezier — extract poles for this segment
            let bez_poles: Vec<GpPnt> = poles[seg..=seg + degree].to_vec();
            result.push(bez_poles);
        } else {
            // Insert knots to make Bezier
            let mut local_poles = poles[seg..].to_vec();
            let mut local_knots = knots[seg..].to_vec();

            // Insert interior knots until degree+1 multiplicity
            for _ in 0..insertions_needed {
                let mid = (u0 + u1) * 0.5;
                let idx = super::knots::hunt(&local_knots, mid);
                boehm_insert(&mut local_poles, &local_knots, idx, mid, degree, None);
                local_knots = super::knots::insert_knot(&local_knots, mid, 1);
            }
            let bez_poles = local_poles[..=degree].to_vec();
            result.push(bez_poles);
        }
    }
    result
}

/// Boehm knot insertion algorithm. Inserts knot u (which lies in the span
/// [knots[k], knots[k+1]) with k = hunt(knots, u)) into the B-spline curve,
/// adding one control point. Modifies poles and optionally weights in-place.
/// Source: BSplCLib::Boehm
pub fn boehm_insert(poles: &mut Vec<GpPnt>, knots: &[f64], _idx: usize, u: f64,
                     degree: usize, weights: Option<&mut Vec<f64>>) {
    let n = poles.len();
    if n == 0 { return; }
    let k = super::knots::hunt(knots, u).min(n - 1);
    let m = super::knots::multiplicity(knots, u);
    let p = degree;

    // New control points Q_0..Q_n:
    //   Q_i = P_i                      for i <= k-p
    //   Q_i = (1-a_i) P_{i-1} + a_i P_i  for k-p+1 <= i <= k-m
    //   Q_i = P_{i-1}                  for i >= k-m+1
    // where a_i = (u - U_i) / (U_{i+p} - U_i).
    let mut new_poles = Vec::with_capacity(n + 1);
    for i in 0..=n {
        if i <= k.saturating_sub(p) {
            new_poles.push(poles[i]);
        } else if i <= k.saturating_sub(m) {
            let alpha = knot_alpha(knots, u, i, p);
            new_poles.push(lerp_point(poles[i - 1], poles[i], alpha));
        } else {
            new_poles.push(poles[i - 1]);
        }
    }
    *poles = new_poles;

    if let Some(w) = weights {
        let mut new_w = Vec::with_capacity(w.len() + 1);
        for i in 0..=n {
            if i <= k.saturating_sub(p) {
                new_w.push(w[i]);
            } else if i <= k.saturating_sub(m) {
                let alpha = knot_alpha(knots, u, i, p);
                new_w.push((1.0 - alpha) * w[i - 1] + alpha * w[i]);
            } else {
                new_w.push(w[i - 1]);
            }
        }
        *w = new_w;
    }
}

/// Interpolation factor a_i = (u - U_i) / (U_{i+p} - U_i), 0 when degenerate.
fn knot_alpha(knots: &[f64], u: f64, i: usize, p: usize) -> f64 {
    if i + p < knots.len() && (knots[i + p] - knots[i]).abs() > 1e-30 {
        ((u - knots[i]) / (knots[i + p] - knots[i])).clamp(0.0, 1.0)
    } else {
        0.0
    }
}

/// Linear interpolation a + t * (b - a).
fn lerp_point(a: GpPnt, b: GpPnt, t: f64) -> GpPnt {
    GpPnt::new(
        a.x() + t * (b.x() - a.x()),
        a.y() + t * (b.y() - a.y()),
        a.z() + t * (b.z() - a.z()),
    )
}

/// Count segments in B-spline between repeating interior knots.
fn count_segments(knots: &[f64], degree: usize) -> usize {
    let mut count = 0usize;
    for i in degree..knots.len() - degree - 1 {
        if knots[i] < knots[i + 1] { count += 1; }
    }
    count.max(1)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gp::GpPnt;

    #[test]
    fn flat_bezier_linear() {
        let poles = vec![GpPnt::new(0.,0.,0.), GpPnt::new(1.,0.,0.)];
        let knots = vec![0.,0.,1.,1.];
        let segments = flat_bezier_coefficients(&poles, None, &knots, 1, 0.001);
        assert_eq!(segments.len(), 1);
        assert_eq!(segments[0].len(), 2);
    }
}
