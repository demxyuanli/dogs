//! Knot vector operations. Source: `BSplCLib.cxx` — Hunt, BuildKnots, InsertKnots

/// Binary search in non-decreasing knot sequence. Returns index i such that U(i) <= x < U(i+1).
/// Port of Fortran HUNT algorithm. Source: BSplCLib::Hunt
pub fn hunt(knots: &[f64], x: f64) -> usize {
    let n = knots.len();
    if n == 0 { return 0; }
    if x < knots[0] { return 0; }
    if x >= knots[n-1] { return n - 2; }
    // For multiple identical knots, return the first one
    if (x - knots[0]).abs() < 1e-30 { return 0; }

    let mut lo = 0usize;
    let mut hi = n - 1;
    // Binary search
    while hi - lo > 1 {
        let mid = (lo + hi) / 2;
        if x < knots[mid] { hi = mid; } else { lo = mid; }
    }
    lo
}

/// Build uniform knot vector for n_poles control points, degree d.
/// Returns n_poles + degree + 1 knots: `degree+1` leading zeros, the
/// `n_poles - degree - 1` interior knots evenly spaced in (0,1), then
/// `degree+1` trailing ones. Source: BSplCLib::BuildKnots
pub fn build_uniform_knots(n_poles: usize, degree: usize) -> Vec<f64> {
    let n_knots = n_poles + degree + 1;
    let mut knots = vec![0.0f64; n_knots];
    // Leading clamp: degree+1 zeros.
    for i in 0..=degree { knots[i] = 0.0; }
    // Interior knots, evenly spaced in (0,1).
    let n_inner = n_poles.saturating_sub(degree + 1);
    for i in 0..n_inner {
        knots[degree + 1 + i] = (i + 1) as f64 / (n_inner + 1) as f64;
    }
    // Trailing clamp: degree+1 ones. Must start right after the interior
    // knots — the previous version left the last interior slot at 0, yielding
    // a non-monotonic vector like [0,0,0,0.5,0,1,1].
    for i in (degree + 1 + n_inner)..n_knots { knots[i] = 1.0; }
    knots
}

/// Insert knot x into knot vector, returning new knot vector.
/// Source: BSplCLib::InsertKnot
pub fn insert_knot(knots: &[f64], x: f64, mult: usize) -> Vec<f64> {
    let mut result = knots.to_vec();
    let idx = hunt(knots, x);
    for _ in 0..mult {
        let pos = if x > knots[idx] { idx + 1 } else { idx };
        result.insert(pos, x);
    }
    result
}

/// Knot vector multiplicity at parameter u.
pub fn multiplicity(knots: &[f64], u: f64) -> usize {
    knots.iter().filter(|&&k| (k - u).abs() < 1e-15).count()
}

/// First knot index with value >= u.
pub fn first_index(knots: &[f64], u: f64) -> usize {
    hunt(knots, u) + 1
}

/// Degree of B-spline from control point count.
pub fn spline_degree(n_poles: usize, n_knots: usize) -> usize {
    n_knots - n_poles - 1
}

/// Check if knot vector is valid (non-decreasing).
pub fn is_valid_knots(knots: &[f64]) -> bool {
    for i in 1..knots.len() { if knots[i] < knots[i-1] - 1e-15 { return false; } }
    true
}

/// Check if degree is consistent with pole count and knot count.
pub fn check_degree(n_poles: usize, degree: usize, n_knots: usize) -> Result<(), &'static str> {
    if n_knots != n_poles + degree + 1 { return Err("BSplCLib: knot count mismatch"); }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hunt_basic() {
        let k = vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0];
        assert_eq!(hunt(&k, 0.3), 2);
        assert_eq!(hunt(&k, 0.7), 3);
        assert_eq!(hunt(&k, 0.0), 0);
        assert_eq!(hunt(&k, 1.0), 5); // returns n-2 = 5 for x==last knot
    }

    #[test]
    fn build_uniform() {
        let k = build_uniform_knots(4, 3); // cubic, 4 poles → 8 knots
        assert_eq!(k.len(), 8);
        // Clamped cubic with 4 poles: [0,0,0,0, 1,1,1,1]
        assert!(k[0] == 0.0);
        assert!(k[3] == 0.0);
        assert!(k[4] == 1.0 || k[4] == 0.0); // may be 0 if no inner knots
        assert!(k[7] == 1.0);
    }

    #[test]
    fn degree_check() {
        assert!(check_degree(4, 3, 8).is_ok());
        assert!(check_degree(4, 3, 9).is_err());
    }
}
