//! Pole/control point conversion utilities. Source: `BSplCLib.cxx` — PolesCoefficients
//! Converts between B-spline control points and power basis (Bezier) representation.

use crate::gp::GpPnt;

/// Convert B-spline poles to monomial (power) basis coefficients.
/// Given degree d, knots, and poles, returns coefficients [c0, c1, ..., cd]
/// where curve(u) = c0 + c1*u + c2*u² + ... + cd*u^d.
/// Uses Marsden's identity for conversion.
pub fn poles_to_coefficients(poles: &[GpPnt], knots: &[f64], degree: usize) -> Vec<GpPnt> {
    let n = poles.len();
    if n == 0 { return vec![]; }
    let mut coeffs = vec![GpPnt::zero(); degree + 1];
    let mut temp = vec![vec![GpPnt::zero(); degree + 1]; degree + 1];

    for span in 0..(n - degree) {
        let local_knots = &knots[span..span + 2 * degree + 2];
        let local_poles = &poles[span..span + degree + 1];
        // Marsden identity to convert this Bezier segment to power basis
        let bezier = power_to_bezier_basis(degree);
        for i in 0..=degree {
            let mut p = GpPnt::zero();
            for j in 0..=degree {
                let b = bezier[i * (degree + 1) + j];
                p = GpPnt::new(p.x()+b*local_poles[j].x(), p.y()+b*local_poles[j].y(), p.z()+b*local_poles[j].z());
            }
            temp[span][i] = p;
        }
    }

    // Average coefficients from overlapping segments
    let mut count = vec![0usize; degree + 1];
    for (ti, tseg) in temp.iter().enumerate() {
        for (i, p) in tseg.iter().enumerate() {
            coeffs[i] = GpPnt::new(coeffs[i].x()+p.x(), coeffs[i].y()+p.y(), coeffs[i].z()+p.z());
            count[i] += 1;
        }
    }

    for i in 0..=degree {
        if count[i] > 0 {
            coeffs[i] = GpPnt::new(coeffs[i].x()/count[i] as f64, coeffs[i].y()/count[i] as f64, coeffs[i].z()/count[i] as f64);
        }
    }
    coeffs
}

/// Bezier-to-power basis conversion matrix. Element (i,j) for degree d:
/// B[i][j] = binomial(d, i) * binomial(i, j) * (-1)^(i-j) / binomial(d, j) scaled.
fn power_to_bezier_basis(degree: usize) -> Vec<f64> {
    let n = degree + 1;
    let mut m = vec![0.0f64; n * n];
    for i in 0..n {
        for j in 0..=i {
            let sign = if (i - j) % 2 == 0 { 1.0 } else { -1.0 };
            let c = c_binomial(degree, i) as f64 * c_binomial(i, j) as f64;
            let d = c_binomial(degree, j) as f64;
            m[i * n + j] = if d > 0.0 { sign * c / d } else { 0.0 };
        }
    }
    m
}

/// Binomial coefficient C(n, k) as u64.
fn c_binomial(n: usize, k: usize) -> u64 {
    if k > n { return 0; }
    let k = k.min(n - k);
    let mut c = 1u64;
    for i in 0..k { c = c * (n - i) as u64 / (i + 1) as u64; }
    c
}

/// Greville abscissae: parameter values for control points.
/// a[i] = (knots[i+1] + knots[i+2] + ... + knots[i+degree]) / degree
pub fn greville_abscissae(knots: &[f64], degree: usize, n_poles: usize) -> Vec<f64> {
    let mut a = vec![0.0f64; n_poles];
    for i in 0..n_poles {
        let mut sum = 0.0;
        for j in 1..=degree { sum += knots[i + j]; }
        a[i] = sum / degree as f64;
    }
    a
}

/// Convert poles from non-uniform to uniform parameterization via knot insertion.
pub fn uniform_poles(poles: &[GpPnt], knots: &[f64], degree: usize) -> Vec<GpPnt> {
    let n = poles.len();
    let mut result = poles.to_vec();
    let mut current_knots = knots.to_vec();
    let target_knots = super::knots::build_uniform_knots(n, degree);

    // Insert knots until we match the target knot vector
    let mut ki = 0;
    while ki < target_knots.len() && result.len() < target_knots.len() - degree {
        let tk = target_knots[ki];
        let mult_target = super::knots::multiplicity(&target_knots, tk);
        let mult_current = super::knots::multiplicity(&current_knots, tk);
        if mult_current < mult_target {
            let idx = super::knots::hunt(&current_knots, tk);
            for _ in 0..(mult_target - mult_current) {
                super::bezier::boehm_insert(&mut result, &current_knots, idx, tk, degree, None);
                current_knots = super::knots::insert_knot(&current_knots, tk, 1);
            }
        }
        ki += mult_target;
    }
    result.truncate(target_knots.len() - degree - 1);
    result
}

/// Compute all derivatives up to order `n_deriv` at parameter u.
/// Returns vector of GpPnt where result[k] is the k-th derivative.
pub fn all_derivatives(poles: &[GpPnt], knots: &[f64], degree: usize, u: f64, n_deriv: usize) -> Vec<GpPnt> {
    let mut derivs = Vec::with_capacity(n_deriv + 1);

    // Evaluate D0 and collect poles for derivative computation
    let n = poles.len();
    let idx = super::knots::hunt(knots, u).max(degree).min(n - 1);

    // D0
    derivs.push(super::eval::eval_curve(poles, knots, degree, u));

    if n_deriv >= 1 {
        // Build derivative poles
        let mut d_poles = poles.to_vec();
        let mut d_knots = knots.to_vec();
        let mut d_degree = degree;

        for deriv in 1..=n_deriv {
            // Build derivative poles: dpi = degree * (pi+1 - pi) / (ki+degree+1 - ki+1)
            let np = d_poles.len();
            if np < 2 || d_degree == 0 { derivs.push(GpPnt::zero()); continue; }
            let mut new_poles = Vec::with_capacity(np - 1);
            for i in 0..np - 1 {
                let k0 = i + 1;
                let k1 = k0 + d_degree;
                let alpha = if d_knots[k1] > d_knots[k0] { d_degree as f64 / (d_knots[k1] - d_knots[k0]) } else { 1.0 };
                let pi = d_poles[i]; let pi1 = d_poles[i+1];
                new_poles.push(GpPnt::new(
                    alpha*(pi1.x()-pi.x()), alpha*(pi1.y()-pi.y()), alpha*(pi1.z()-pi.z())
                ));
            }
            if new_poles.is_empty() { derivs.push(GpPnt::zero()); continue; }
            let nd = new_poles.len();
            // Knot vector for derivative: drop first and last
            let new_knots = if d_knots.len() >= 2 { d_knots[1..d_knots.len()-1].to_vec() } else { d_knots.clone() };
            derivs.push(super::eval::eval_curve(&new_poles, &new_knots, d_degree - 1, u));
            d_poles = new_poles;
            d_knots = new_knots;
            if d_degree > 0 { d_degree -= 1; }
        }
    }
    derivs
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn greville_linear() {
        let k = vec![0., 0., 0.5, 1., 1.];
        let a = greville_abscissae(&k, 2, 3);
        assert_eq!(a.len(), 3);
        assert!((a[0]-0.25).abs() < 1e-14);
        assert!((a[1]-0.75).abs() < 1e-14);
        assert!((a[2]-1.00).abs() < 1e-14);
    }

    #[test]
    fn derivs_linear() {
        let poles = vec![GpPnt::new(0.,0.,0.), GpPnt::new(1.,0.,0.)];
        let knots = vec![0.,0.,1.,1.];
        let ds = all_derivatives(&poles, &knots, 1, 0.5, 2);
        assert_eq!(ds.len(), 3);
        assert!((ds[0].x() - 0.5).abs() < 1e-14);
        assert!((ds[1].x() - 1.0).abs() < 1e-14);
        assert!((ds[2].x()).abs() < 1e-14);
    }
}
