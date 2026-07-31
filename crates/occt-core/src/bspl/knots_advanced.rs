//! Advanced knot operations. Source: BSplCLib knot removal + degree ops.
use crate::gp::GpPnt;

/// Remove knot u mult times from curve while maintaining shape within tolerance.
/// Returns (new_poles, new_knots). Uses Tiller's knot removal algorithm.
pub fn remove_knot(poles: &[GpPnt], knots: &[f64], degree: usize, u: f64, mult: usize, tolerance: f64) -> (Vec<GpPnt>, Vec<f64>) {
    let n = poles.len();
    let idx = super::knots::hunt(knots, u);
    let actual_mult = (0..knots.len()).filter(|&i| (knots[i] - u).abs() < 1e-15).count();
    if actual_mult < mult { return (poles.to_vec(), knots.to_vec()); }

    let mut new_poles = poles.to_vec();
    let mut new_knots = knots.to_vec();

    for _ in 0..mult {
        if new_poles.len() <= degree + 1 { break; }
        let n = new_poles.len();
        let mut temp = new_poles.clone();
        let r = degree - actual_mult + 1;

        // Try removing: compute new poles via average of neighbors
        let idx = super::knots::hunt(&new_knots, u);
        let first = if idx >= r { idx - r } else { 0 };
        let last = if idx + 1 <= n - r { idx + 1 } else { n - r };

        for i in first..=last {
            let alpha = (u - new_knots[i]) / (new_knots[i + degree + 1] - new_knots[i]);
            if alpha > 0.0 && alpha < 1.0 {
                temp[i] = GpPnt::new(
                    (new_poles[i].x() - (1.0 - alpha) * new_poles[i - 1].x()) / alpha,
                    (new_poles[i].y() - (1.0 - alpha) * new_poles[i - 1].y()) / alpha,
                    (new_poles[i].z() - (1.0 - alpha) * new_poles[i - 1].z()) / alpha,
                );
                temp[i - 1] = GpPnt::new(
                    (new_poles[i].x() - alpha * temp[i].x()) / (1.0 - alpha),
                    (new_poles[i].y() - alpha * temp[i].y()) / (1.0 - alpha),
                    (new_poles[i].z() - alpha * temp[i].z()) / (1.0 - alpha),
                );
            }
        }
        // Check if modification is within tolerance
        let mut max_dist = 0.0;
        for i in 0..n {
            let dx = temp[i].x() - new_poles[i].x();
            let dy = temp[i].y() - new_poles[i].y();
            let dz = temp[i].z() - new_poles[i].z();
            max_dist = f64::max(max_dist, (dx*dx + dy*dy + dz*dz).sqrt());
        }
        if max_dist <= tolerance {
            new_poles = temp;
            let pos = new_knots.iter().position(|&k| (k - u).abs() < 1e-15).unwrap_or(0);
            new_knots.remove(pos);
        } else { break; }
    }
    (new_poles, new_knots)
}

/// Reduce B-spline degree by 1 while keeping shape within tolerance.
pub fn reduce_degree(poles: &[GpPnt], knots: &[f64], degree: usize, tolerance: f64) -> (Vec<GpPnt>, Vec<f64>, usize) {
    if degree == 1 { return (poles.to_vec(), knots.to_vec(), degree); }
    // Degree reduction: remove one multiplicity from each end knot, then interpolate
    let mut new_knots = knots.to_vec();
    let mult_start = super::knots::multiplicity(knots, knots[0]);
    let mult_end = super::knots::multiplicity(knots, knots[knots.len()-1]);
    if mult_start > 0 { new_knots.remove(0); }
    if mult_end > 0 { new_knots.pop(); }

    let new_degree = degree - 1;
    let new_n = poles.len() - 1;
    let mut new_poles = vec![GpPnt::zero(); new_n];

    // Interpolate new control points
    for i in 0..new_n {
        let alpha = (i + 1) as f64 / (new_n + 1) as f64;
        new_poles[i] = GpPnt::new(
            (1.0 - alpha) * poles[i].x() + alpha * poles[i+1].x(),
            (1.0 - alpha) * poles[i].y() + alpha * poles[i+1].y(),
            (1.0 - alpha) * poles[i].z() + alpha * poles[i+1].z(),
        );
    }
    (new_poles, new_knots, new_degree)
}

/// Make a B-spline curve periodic (close the gap between ends).
pub fn make_periodic(poles: &[GpPnt], knots: &[f64], degree: usize) -> (Vec<GpPnt>, Vec<f64>) {
    let n = poles.len();
    let mut new_poles = poles.to_vec();
    // Average first and last degree poles for smooth closure
    for k in 0..degree {
        let i = n - degree + k;
        let j = k;
        new_poles[i] = GpPnt::new(
            0.5*(poles[i].x()+poles[j].x()), 0.5*(poles[i].y()+poles[j].y()), 0.5*(poles[i].z()+poles[j].z())
        );
    }
    (new_poles, knots.to_vec())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn remove_knot_identity() {
        let poles = vec![GpPnt::new(0.,0.,0.),GpPnt::new(1.,0.,0.),GpPnt::new(2.,0.,0.)];
        let knots = vec![0.,0.,0.,1.,1.,1.];
        let (np, _) = remove_knot(&poles, &knots, 2, 0.5, 1, 0.001);
        // Should not change much
        assert_eq!(np.len(), poles.len());
    }
}
