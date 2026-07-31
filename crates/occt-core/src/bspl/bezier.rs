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

/// Boehm knot insertion algorithm. Inserts knot u at position idx in the B-spline curve.
/// Modifies poles and optionally weights in-place.
/// Source: BSplCLib::Boehm
pub fn boehm_insert(poles: &mut Vec<GpPnt>, knots: &[f64], idx: usize, u: f64,
                     degree: usize, mut weights: Option<&mut Vec<f64>>) {
    let n = poles.len();
    let new_idx = super::knots::hunt(knots, u);
    let mult = super::knots::multiplicity(knots, u);

    // Compute new poles using de Boor insertion
    poles.insert(new_idx + 1, GpPnt::zero());
    if let Some(ref mut w) = weights { w.insert(new_idx + 1, 1.0); }

    for i in (new_idx - degree + 1)..=new_idx {
        let i0 = i.max(0);
        let ki = idx + i - new_idx + degree;
        let alpha = if (knots[ki] - knots[i0]).abs() > 1e-30 {
            (u - knots[i0]) / (knots[ki] - knots[i0])
        } else { 0.0 };

        if alpha > 0.0 && alpha < 1.0 {
            let pi = poles[i];
            let pi1 = poles[i + 1];
            poles[i] = GpPnt::new(
                (1.0 - alpha) * pi.x() + alpha * pi1.x(),
                (1.0 - alpha) * pi.y() + alpha * pi1.y(),
                (1.0 - alpha) * pi.z() + alpha * pi1.z(),
            );
            poles[i + 1] = GpPnt::new(
                alpha * pi.x() + (1.0 - alpha) * pi1.x(),
                alpha * pi.y() + (1.0 - alpha) * pi1.y(),
                alpha * pi.z() + (1.0 - alpha) * pi1.z(),
            );
            if let Some(ref mut w) = weights {
                let wi = w[i];
                let wi1 = w[i + 1];
                w[i] = (1.0 - alpha) * wi + alpha * wi1;
                w[i + 1] = alpha * wi + (1.0 - alpha) * wi1;
            }
        }
    }
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
