//! B-spline curve interpolation through point and tangent data.
//! Source: `GeomAPI_Interpolate` / `Geom_BSplineCurve` (global interpolation).
//!
//! The data points are parameterized by cumulative chord length onto `[0, 1]`
//! and a clamped B-spline of the requested degree is fit through them by
//! solving the collocation system `Σ_j N_{j,d}(t_i)·P_j = Q_i` for the control
//! poles, one coordinate at a time.  A degree-1 fit is the polyline whose
//! breakpoints are exactly the chord-length parameters.  The tangent-aware
//! variant appends derivative rows to the system so the fitted curve matches a
//! prescribed tangent direction at every node (Hermite-style).

use std::sync::Arc;

use crate::bspline_curve::GeomBSplineCurve;
use crate::curve::Curve;
use occt_core::bspl::knots::{build_uniform_knots, hunt};
use occt_core::gp::{GpPnt, GpVec};

/// An interpolating B-spline curve and its defining data.
///
/// `curve` is the parametric curve (evaluable through the [`Curve`] trait),
/// `knots` its clamped knot vector and `degree` the spline degree.  The
/// parameter range is always `[0, 1]`.
pub struct GeomInterpCurve {
    pub curve: Arc<dyn Curve>,
    pub knots: Vec<f64>,
    pub degree: usize,
}

/// Chord-length parameterization of `points` onto `[0, 1]`.
///
/// `t_0 = 0`, `t_i = t_{i−1} + ‖P_i − P_{i−1}‖`, normalized by the total
/// chord length.  When all points coincide the parameters are all zero (the
/// caller must then avoid a degenerate knot vector).
fn chord_length_params(points: &[GpPnt]) -> Vec<f64> {
    let n = points.len();
    let mut t = vec![0.0; n];
    let mut acc = 0.0;
    for i in 1..n {
        acc += points[i].distance(&points[i - 1]);
        t[i] = acc;
    }
    if acc > 0.0 {
        for x in t.iter_mut() {
            *x /= acc;
        }
    }
    t
}

/// Full basis vector `(N_{0,d}(u), …, N_{np−1,d}(u))` for a clamped knot
/// vector, computed with the triangular scheme of The NURBS Book (Alg. A2.2).
///
/// The knot span is located with a binary search and the `degree + 1` non-zero
/// values are propagated through the Cox–de Boor recursion.  At the right
/// clamp `u = 1` the search clamps to the last span so the last basis function
/// evaluates to exactly 1 (partition of unity holds).
fn basis_vector(knots: &[f64], degree: usize, u: f64) -> Vec<f64> {
    let n_poles = knots.len() - degree - 1;
    let s = hunt(knots, u).max(degree).min(n_poles - 1);
    let mut left = vec![0.0; degree + 1];
    let mut right = vec![0.0; degree + 1];
    let mut n = vec![0.0; degree + 1];
    n[0] = 1.0;
    for j in 1..=degree {
        left[j] = u - knots[s + 1 - j];
        right[j] = knots[s + j] - u;
        let mut saved = 0.0;
        for r in 0..j {
            let denom = right[r + 1] + left[j - r];
            let temp = if denom.abs() > 1e-15 { n[r] / denom } else { 0.0 };
            n[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        n[j] = saved;
    }
    let mut basis = vec![0.0; n_poles];
    for r in 0..=degree {
        let gi = s - degree + r;
        if gi < n_poles {
            basis[gi] = n[r];
        }
    }
    basis
}

/// Derivative basis vector `(N'_{0,d}(u), …, N'_{np−1,d}(u))`.
///
/// The basis derivatives are approximated by finite differences of
/// [`basis_vector`].  Interior parameters use a central difference; the two
/// clamped ends use one-sided differences because `u ± h` would leave the
/// parameter domain.  The approximation error is O(h²) centrally and O(h) at
/// the ends (h = 1e-7), far below the tangent-direction tolerance the
/// interpolation tests require.
fn basis_derivative_vector(knots: &[f64], degree: usize, u: f64) -> Vec<f64> {
    let h = 1e-7;
    let n_poles = knots.len() - degree - 1;
    let mut d = vec![0.0; n_poles];
    if u <= 0.0 {
        let fp = basis_vector(knots, degree, u + h);
        let f0 = basis_vector(knots, degree, u);
        for i in 0..n_poles {
            d[i] = (fp[i] - f0[i]) / h;
        }
    } else if u >= 1.0 {
        let f0 = basis_vector(knots, degree, u);
        let fm = basis_vector(knots, degree, u - h);
        for i in 0..n_poles {
            d[i] = (f0[i] - fm[i]) / h;
        }
    } else {
        let fp = basis_vector(knots, degree, u + h);
        let fm = basis_vector(knots, degree, u - h);
        for i in 0..n_poles {
            d[i] = (fp[i] - fm[i]) / (2.0 * h);
        }
    }
    d
}

/// Solve `A·x = b` for a square `A` by partial-pivoted Gaussian elimination;
/// `None` when the matrix is numerically singular.
fn gauss_solve(a: &[Vec<f64>], b: &[f64]) -> Option<Vec<f64>> {
    let n = a.len();
    if n == 0 {
        return Some(Vec::new());
    }
    let mut m: Vec<Vec<f64>> = (0..n)
        .map(|i| {
            let mut row = a[i].clone();
            row.push(b[i]);
            row
        })
        .collect();
    for col in 0..n {
        let mut piv = col;
        let mut best = m[col][col].abs();
        for r in (col + 1)..n {
            if m[r][col].abs() > best {
                best = m[r][col].abs();
                piv = r;
            }
        }
        if best < 1e-14 {
            return None;
        }
        if piv != col {
            m.swap(piv, col);
        }
        let pv = m[col][col];
        for r in (col + 1)..n {
            let f = m[r][col] / pv;
            if f == 0.0 {
                continue;
            }
            for c in col..=n {
                m[r][c] -= f * m[col][c];
            }
        }
    }
    let mut x = vec![0.0; n];
    for r in (0..n).rev() {
        let mut s = m[r][n];
        for c in (r + 1)..n {
            s -= m[r][c] * x[c];
        }
        x[r] = s / m[r][r];
    }
    Some(x)
}

/// Interpolate `points` with a clamped B-spline of degree `degree`.
///
/// The points are parameterized by chord length and a clamped uniform knot
/// vector of `n` poles is built.  For degree 1 the control points are the data
/// points themselves and the interior knots are placed at the chord-length
/// parameters, giving a polyline that passes through every point.  For degree
/// `≥ 2` the collocation matrix `B[i][j] = N_{j,d}(t_i)` is solved per
/// coordinate.
///
/// # Errors
/// Fewer than two points, a zero degree, or (for `degree ≥ 2`) fewer than
/// `degree + 1` points produce an `Err`; a singular collocation matrix also
/// produces an `Err`.
pub fn interpolate_points(points: &[GpPnt], degree: usize) -> Result<GeomInterpCurve, String> {
    let n = points.len();
    if n < 2 {
        return Err("interpolate_points: need at least 2 points".to_string());
    }
    if degree < 1 {
        return Err("interpolate_points: degree must be at least 1".to_string());
    }
    let params = chord_length_params(points);
    if degree == 1 {
        let mut knots = vec![0.0, 0.0];
        for &t in &params[1..n - 1] {
            knots.push(t);
        }
        knots.push(1.0);
        knots.push(1.0);
        let curve = GeomBSplineCurve::new(points.to_vec(), knots.clone(), 1)
            .map_err(|e| e.to_string())?;
        return Ok(GeomInterpCurve { curve: Arc::new(curve), knots, degree: 1 });
    }
    if n < degree + 1 {
        return Err(format!(
            "interpolate_points: need at least {} points for degree {degree}",
            degree + 1
        ));
    }
    let knots = build_uniform_knots(n, degree);
    let colloc: Vec<Vec<f64>> = params.iter().map(|&u| basis_vector(&knots, degree, u)).collect();
    let px = gauss_solve(&colloc, &points.iter().map(|p| p.x()).collect::<Vec<_>>())
        .ok_or("interpolate_points: singular collocation matrix")?;
    let py = gauss_solve(&colloc, &points.iter().map(|p| p.y()).collect::<Vec<_>>())
        .ok_or("interpolate_points: singular collocation matrix")?;
    let pz = gauss_solve(&colloc, &points.iter().map(|p| p.z()).collect::<Vec<_>>())
        .ok_or("interpolate_points: singular collocation matrix")?;
    let poles: Vec<GpPnt> = (0..n).map(|i| GpPnt::new(px[i], py[i], pz[i])).collect();
    let curve = GeomBSplineCurve::new(poles, knots.clone(), degree).map_err(|e| e.to_string())?;
    Ok(GeomInterpCurve { curve: Arc::new(curve), knots, degree })
}

/// Hermite-style interpolation through `points` with a matching tangent at
/// every node.
///
/// The system is enlarged to `2n` unknowns (control poles) and `2n` rows: the
/// first `n` rows interpolate the point positions at the chord-length
/// parameters, the last `n` rows interpolate the first derivatives
/// `Σ_j N'_{j,d}(t_i)·P_j = T_i`.  The fitted curve therefore passes through
/// the data points and its tangent direction at each node matches the supplied
/// `GpVec` (the magnitude of the curve's derivative is set by `|T_i|`).
///
/// # Errors
/// Fewer than two points, a tangent/point count mismatch, a zero degree, or a
/// singular collocation matrix produce an `Err`.
pub fn interpolate_with_tangents(
    points: &[GpPnt],
    tangents: &[GpVec],
    degree: usize,
) -> Result<GeomInterpCurve, String> {
    let n = points.len();
    if n < 2 {
        return Err("interpolate_with_tangents: need at least 2 points".to_string());
    }
    if tangents.len() != n {
        return Err("interpolate_with_tangents: tangent count must match point count".to_string());
    }
    if degree < 1 {
        return Err("interpolate_with_tangents: degree must be at least 1".to_string());
    }
    let np = 2 * n;
    let params = chord_length_params(points);
    let knots = build_uniform_knots(np, degree);
    let mut a: Vec<Vec<f64>> = Vec::with_capacity(np);
    let mut bx = vec![0.0; np];
    let mut by = vec![0.0; np];
    let mut bz = vec![0.0; np];
    for (i, &u) in params.iter().enumerate() {
        a.push(basis_vector(&knots, degree, u));
        bx[i] = points[i].x();
        by[i] = points[i].y();
        bz[i] = points[i].z();
    }
    for (i, &u) in params.iter().enumerate() {
        a.push(basis_derivative_vector(&knots, degree, u));
        bx[n + i] = tangents[i].x();
        by[n + i] = tangents[i].y();
        bz[n + i] = tangents[i].z();
    }
    let px = gauss_solve(&a, &bx).ok_or("interpolate_with_tangents: singular collocation matrix")?;
    let py = gauss_solve(&a, &by).ok_or("interpolate_with_tangents: singular collocation matrix")?;
    let pz = gauss_solve(&a, &bz).ok_or("interpolate_with_tangents: singular collocation matrix")?;
    let poles: Vec<GpPnt> = (0..np).map(|i| GpPnt::new(px[i], py[i], pz[i])).collect();
    let curve = GeomBSplineCurve::new(poles, knots.clone(), degree).map_err(|e| e.to_string())?;
    Ok(GeomInterpCurve { curve: Arc::new(curve), knots, degree })
}

/// Evaluate the interpolating curve at parameter `t ∈ [0, 1]`.
pub fn eval_interp(c: &GeomInterpCurve, t: f64) -> GpPnt {
    c.curve.d0(t)
}

/// Arc length of the interpolating curve estimated by polygonal sampling.
///
/// The curve is sampled at `samples` equal parameter increments and the chord
/// lengths between consecutive samples are summed.  `samples.max(1)` guards
/// against a zero argument.  The estimate converges to the true arc length as
/// the sampling is refined.
pub fn interp_length(c: &GeomInterpCurve, samples: usize) -> f64 {
    let samples = samples.max(1);
    let mut len = 0.0;
    let mut prev = c.curve.d0(0.0);
    for i in 1..=samples {
        let t = i as f64 / samples as f64;
        let p = c.curve.d0(t);
        len += prev.distance(&p);
        prev = p;
    }
    len
}

#[cfg(test)]
mod tests {
    use super::*;

    fn line_points(n: usize) -> Vec<GpPnt> {
        (0..n).map(|i| GpPnt::new(i as f64, 0.0, 0.0)).collect()
    }

    fn parabola_points(n: usize) -> Vec<GpPnt> {
        (0..n)
            .map(|i| {
                let x = -1.0 + 2.0 * i as f64 / (n - 1) as f64;
                GpPnt::new(x, x * x, 0.0)
            })
            .collect()
    }

    fn parabola_tangents(n: usize) -> Vec<GpVec> {
        (0..n)
            .map(|i| {
                let x = -1.0 + 2.0 * i as f64 / (n - 1) as f64;
                GpVec::new(1.0, 2.0 * x, 0.0)
            })
            .collect()
    }

    #[test]
    fn line_interpolation_passes_through_points() {
        let pts = line_points(5);
        let c = interpolate_points(&pts, 3).unwrap();
        let params = chord_length_params(&pts);
        for (i, &u) in params.iter().enumerate() {
            let p = eval_interp(&c, u);
            assert!(p.distance(&pts[i]) < 1e-6, "point {i} off by {}", p.distance(&pts[i]));
        }
    }

    #[test]
    fn parabola_interpolation_passes_through_points() {
        let pts = parabola_points(10);
        let c = interpolate_points(&pts, 3).unwrap();
        let params = chord_length_params(&pts);
        for (i, &u) in params.iter().enumerate() {
            let p = eval_interp(&c, u);
            assert!(p.distance(&pts[i]) < 1e-5, "point {i} off by {}", p.distance(&pts[i]));
        }
    }

    #[test]
    fn endpoints_are_exact() {
        let pts = line_points(5);
        let c = interpolate_points(&pts, 3).unwrap();
        assert!(eval_interp(&c, 0.0).distance(&pts[0]) < 1e-9);
        assert!(eval_interp(&c, 1.0).distance(&pts[4]) < 1e-9);
    }

    #[test]
    fn tangent_interpolation_matches_directions() {
        let pts = parabola_points(5);
        let tan = parabola_tangents(5);
        let c = interpolate_with_tangents(&pts, &tan, 3).unwrap();
        let params = chord_length_params(&pts);
        // Passes through the points.
        for (i, &u) in params.iter().enumerate() {
            let p = eval_interp(&c, u);
            assert!(p.distance(&pts[i]) < 1e-5, "point {i} off by {}", p.distance(&pts[i]));
        }
        // Tangent directions match at each node.
        for (i, &u) in params.iter().enumerate() {
            let (_, d1) = c.curve.d1(u);
            let ang = d1.angle(&tan[i]);
            assert!(ang < 1e-3, "tangent {i} angle {ang} rad");
        }
    }

    #[test]
    fn eval_at_chord_param_equals_point() {
        let pts = line_points(5);
        let c = interpolate_points(&pts, 3).unwrap();
        let params = chord_length_params(&pts);
        for (i, &u) in params.iter().enumerate() {
            let p = eval_interp(&c, u);
            assert!(p.distance(&pts[i]) < 1e-9, "point {i} off by {}", p.distance(&pts[i]));
        }
    }

    #[test]
    fn degree_one_polyline_passes_through_points() {
        let pts = vec![
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(1.0, 2.0, 0.0),
            GpPnt::new(2.0, 0.0, 0.0),
            GpPnt::new(3.0, 3.0, 0.0),
        ];
        let c = interpolate_points(&pts, 1).unwrap();
        assert_eq!(c.degree, 1);
        let params = chord_length_params(&pts);
        for (i, &u) in params.iter().enumerate() {
            let p = eval_interp(&c, u);
            assert!(p.distance(&pts[i]) < 1e-9, "point {i} off by {}", p.distance(&pts[i]));
        }
    }

    #[test]
    fn interp_length_matches_analytic() {
        // Straight line: length is exactly the span.
        let pts = line_points(6);
        let c = interpolate_points(&pts, 3).unwrap();
        let l = interp_length(&c, 100);
        assert!((l - 5.0).abs() < 1e-6, "line length {l}");
        // Quarter unit circle arc: length = π/2.
        let n = 25;
        let circ: Vec<GpPnt> = (0..n)
            .map(|i| {
                let th = std::f64::consts::FRAC_PI_2 * i as f64 / (n - 1) as f64;
                GpPnt::new(th.cos(), th.sin(), 0.0)
            })
            .collect();
        let cc = interpolate_points(&circ, 3).unwrap();
        let cl = interp_length(&cc, 2000);
        assert!((cl - std::f64::consts::FRAC_PI_2).abs() < 1e-3, "arc length {cl}");
    }

    #[test]
    fn too_few_points_is_error() {
        assert!(interpolate_points(&[GpPnt::new(0.0, 0.0, 0.0)], 3).is_err());
        assert!(interpolate_points(&[], 3).is_err());
        assert!(interpolate_points(&line_points(4), 0).is_err());
        // Fewer than degree + 1 points for a higher degree.
        assert!(interpolate_points(&line_points(3), 4).is_err());
        // Tangent variant validates its input too.
        let one = [GpPnt::new(0.0, 0.0, 0.0)];
        assert!(interpolate_with_tangents(&one, &[], 3).is_err());
    }
}
