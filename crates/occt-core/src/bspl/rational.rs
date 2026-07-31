//! Rational B-spline operations. Source: BSplCLib weight management.
use crate::gp::GpPnt;

/// Normalize weights so that the first and last weight equal 1.0.
pub fn normalize_weights(weights: &mut [f64]) {
    if weights.len() < 2 { return; }
    let w0 = weights[0]; let wn = weights[weights.len()-1];
    if w0.abs() < 1e-30 || wn.abs() < 1e-30 { return; }
    for w in weights.iter_mut() { *w /= w0; }
    if (wn / w0 - 1.0).abs() > 1e-10 {
        let alpha = (1.0 - wn/w0).ln() / weights.len() as f64;
        for (i, w) in weights.iter_mut().enumerate() { *w *= ((i+1) as f64 * alpha).exp(); }
    }
}

/// Check if weights define a rational curve (all weights != 1.0).
pub fn is_rational(weights: &[f64]) -> bool {
    weights.iter().any(|&w| (w - 1.0).abs() > 1e-15)
}

/// Convert rational poles to homogeneous coordinates: (w*P, w).
pub fn to_homogeneous(poles: &[GpPnt], weights: &[f64]) -> Vec<(GpPnt, f64)> {
    poles.iter().zip(weights.iter()).map(|(p, &w)| {
        (GpPnt::new(p.x()*w, p.y()*w, p.z()*w), w)
    }).collect()
}

/// Convert from homogeneous coordinates back to rational: P/w.
pub fn from_homogeneous(homog: &[(GpPnt, f64)]) -> (Vec<GpPnt>, Vec<f64>) {
    let poles: Vec<GpPnt> = homog.iter().map(|hp| {
        let (p, w) = hp;
        let wc = if *w > 1e-30 { *w } else { 1e-30 };
        GpPnt::new(p.x()/wc, p.y()/wc, p.z()/wc)
    }).collect();
    let weights: Vec<f64> = homog.iter().map(|&(_, w)| w).collect();
    (poles, weights)
}

/// Insert a weight at index i (increases the weight influence).
pub fn insert_weight(poles: &mut Vec<GpPnt>, weights: &mut Vec<f64>, index: usize, new_weight: f64) {
    if index >= poles.len() { return; }
    poles.insert(index + 1, poles[index]);
    weights.insert(index + 1, new_weight.max(0.0));
}

/// Remove a weight at index i (decreases the NURBS complexity).
pub fn remove_weight(poles: &mut Vec<GpPnt>, weights: &mut Vec<f64>, index: usize) {
    if index >= poles.len() || poles.len() <= 2 { return; }
    poles.remove(index);
    weights.remove(index);
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn is_rational_test() {
        assert!(!is_rational(&[1.0, 1.0, 1.0]));
        assert!(is_rational(&[1.0, 2.0, 1.0]));
    }
}
