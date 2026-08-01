//! Tensor-product B-spline surface — clamped-knot construction, Cox–de Boor
//! evaluation and grid interpolation. Source: `Geom_BSplineSurface`.
//!
//! A degree `pu × pv` tensor-product B-spline surface is defined by two knot
//! vectors and a control-pole grid. Evaluation uses the Cox–de Boor recurrence
//! in each parameter (de Boor's algorithm on the u rows and v columns). When a
//! weight grid is present the surface is a rational NURBS surface.

use crate::{MathMatrix, MathVector};

/// Tensor-product B-spline surface.
///
/// `poles[i][j]` is the control point with u-index `i` and v-index `j`; the
/// grid has `poles.len()` rows (u direction) and `poles[0].len()` columns
/// (v direction). `weights[i][j]` optionally stores the homogeneous weight for
/// each pole; when `Some`, evaluation divides by the summed weights (NURBS).
///
/// Knot vectors are clamped (the first and last `deg + 1` entries coincide
/// with the parameter bounds), matching OCCT's `Geom_BSplineSurface` for
/// non-periodic surfaces.
#[derive(Debug, Clone)]
pub struct BSplineSurface {
    pub knots_u: Vec<f64>,
    pub knots_v: Vec<f64>,
    pub deg_u: usize,
    pub deg_v: usize,
    pub poles: Vec<Vec<[f64; 3]>>,
    pub weights: Option<Vec<Vec<f64>>>,
}

/// Total knot count across both knot vectors for a surface with `nu`×`nv`
/// poles and degrees `deg_u`, `deg_v`.
///
/// Each direction satisfies `#poles = #knots − deg − 1`, so the total is
/// `(nu + deg_u + 1) + (nv + deg_v + 1)`. This is a knot-count sanity helper:
/// a caller can check a surface built elsewhere before using it.
pub fn bspline_surface_degree(deg_u: usize, deg_v: usize, nu: usize, nv: usize) -> usize {
    (nu + deg_u + 1) + (nv + deg_v + 1)
}

/// Clamped (open-uniform) knot vectors for an `nu`×`nv` pole grid with degrees
/// `deg_u`, `deg_v`. Returns `(knots_u, knots_v)`.
///
/// The first `deg + 1` entries of a direction are 0.0 and the last `deg + 1`
/// are 1.0; the interior knots are uniformly spaced. The result parameterizes
/// the unit square `[0, 1]²`.
pub fn bspline_surface_uniform_knots(
    nu: usize,
    nv: usize,
    deg_u: usize,
    deg_v: usize,
) -> (Vec<f64>, Vec<f64>) {
    (uniform_clamped_knots(nu, deg_u), uniform_clamped_knots(nv, deg_v))
}

/// Clamped uniform knot vector for `n` poles of degree `deg`.
///
/// Requires `n > deg` for a non-degenerate B-spline; the caller validates this
/// before building surfaces/curves.
fn uniform_clamped_knots(n: usize, deg: usize) -> Vec<f64> {
    let interior = n.saturating_sub(deg + 1);
    let mut k = Vec::with_capacity(n + deg + 1);
    for _ in 0..=deg {
        k.push(0.0);
    }
    for i in 1..=interior {
        k.push(i as f64 / (interior + 1) as f64);
    }
    for _ in 0..=deg {
        k.push(1.0);
    }
    k
}

/// Find the knot span `[knots[s], knots[s+1])` containing `u`.
///
/// `s` is returned in `deg..=n` where `n = #knots − deg − 2` (the last pole
/// index). Out-of-range parameters are clamped to the surface ends, so the
/// returned span is always valid for the basis recursion.
fn find_span(knots: &[f64], deg: usize, u: f64) -> usize {
    let n = knots.len() - deg - 2; // number of poles minus one
    if u <= knots[deg] {
        return deg;
    }
    if u >= knots[n + 1] {
        return n;
    }
    // Binary search: knots[lo] <= u < knots[lo+1].
    let mut lo = deg;
    let mut hi = n + 1;
    while hi - lo > 1 {
        let mid = (lo + hi) / 2;
        if u < knots[mid] {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    lo
}

/// Non-zero B-spline basis values on the span containing `u`.
///
/// Returns a vector of length `deg + 1` where `out[r] = N_{span−deg+r, deg}(u)`.
/// The de Boor/Cox–de Boor recursion is used; repeated-knot denominators that
/// vanish contribute zero (the limit of the basis function there).
fn basis_vals(span: usize, deg: usize, u: f64, knots: &[f64]) -> Vec<f64> {
    let mut left = vec![0.0; deg + 1];
    let mut right = vec![0.0; deg + 1];
    let mut n = vec![0.0; deg + 1];
    n[0] = 1.0;
    for j in 1..=deg {
        left[j] = u - knots[span + 1 - j];
        right[j] = knots[span + j] - u;
        let mut saved = 0.0;
        for r in 0..j {
            let den = right[r + 1] + left[j - r];
            let temp = if den.abs() < 1e-300 { 0.0 } else { n[r] / den };
            n[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        n[j] = saved;
    }
    n
}

/// Full basis vector `(N_{0,deg}(u), …, N_{np−1,deg}(u))` for a clamped knot
/// vector with `np` poles. Zeros outside the active support. Useful for
/// partition-of-unity checks and for building collocation matrices.
pub fn bspline_basis_1d(knots: &[f64], deg: usize, u: f64) -> Vec<f64> {
    let np = knots.len() - deg - 1;
    let mut out = vec![0.0; np];
    let span = find_span(knots, deg, u);
    let vals = basis_vals(span, deg, u, knots);
    for r in 0..=deg {
        let idx = span - deg + r;
        if idx < np {
            out[idx] = vals[r];
        }
    }
    out
}

/// Evaluate a 1-D B-spline curve (optionally rational) at parameter `t`.
///
/// This is the building block used by [`eval_bspline_surface`]: the surface is
/// a tensor product, so evaluating along v for each active u basis function and
/// combining with the u basis gives the surface point.
pub fn eval_bspline_curve_1d(
    poles: &[[f64; 3]],
    weights: Option<&[f64]>,
    knots: &[f64],
    deg: usize,
    t: f64,
) -> [f64; 3] {
    let np = poles.len();
    if np == 0 {
        return [0.0; 3];
    }
    let span = find_span(knots, deg, t);
    let vals = basis_vals(span, deg, t, knots);
    let mut num = [0.0; 3];
    let mut denom = 0.0;
    for r in 0..=deg {
        let idx = span - deg + r;
        if idx >= np {
            continue;
        }
        let w = weights.map_or(1.0, |ws| ws[idx]);
        let b = vals[r] * w;
        num[0] += b * poles[idx][0];
        num[1] += b * poles[idx][1];
        num[2] += b * poles[idx][2];
        denom += b;
    }
    if denom.abs() < 1e-300 {
        return [0.0; 3];
    }
    [num[0] / denom, num[1] / denom, num[2] / denom]
}

/// Evaluate the tensor-product B-spline surface at parameter `(u, v)`.
///
/// The Cox–de Boor basis is computed in each direction and the double sum over
/// the active control-pole grid is formed. For a rational surface (`weights`
/// present) the weighted sum is divided by the weight sum, which reproduces
/// standard NURBS evaluation.
pub fn eval_bspline_surface(s: &BSplineSurface, u: f64, v: f64) -> [f64; 3] {
    let nu_p = s.poles.len();
    if nu_p == 0 || s.poles[0].len() == 0 {
        return [0.0; 3];
    }
    let nv_p = s.poles[0].len();
    let span_u = find_span(&s.knots_u, s.deg_u, u);
    let vals_u = basis_vals(span_u, s.deg_u, u, &s.knots_u);
    let span_v = find_span(&s.knots_v, s.deg_v, v);
    let vals_v = basis_vals(span_v, s.deg_v, v, &s.knots_v);

    let mut num = [0.0; 3];
    let mut denom = 0.0;
    for ru in 0..=s.deg_u {
        let i = span_u - s.deg_u + ru;
        if i >= nu_p {
            continue;
        }
        for rv in 0..=s.deg_v {
            let j = span_v - s.deg_v + rv;
            if j >= nv_p {
                continue;
            }
            let w = match &s.weights {
                Some(ws) => ws[i][j],
                None => 1.0,
            };
            let b = vals_u[ru] * vals_v[rv] * w;
            num[0] += b * s.poles[i][j][0];
            num[1] += b * s.poles[i][j][1];
            num[2] += b * s.poles[i][j][2];
            denom += b;
        }
    }
    if denom.abs() < 1e-300 {
        return [0.0; 3];
    }
    [num[0] / denom, num[1] / denom, num[2] / denom]
}

/// Axis-aligned bounding box of the control poles: `(min, max)`.
///
/// Useful for fast rejection tests before sampling the surface; the surface
/// itself lies inside the convex hull of its poles, so this is a conservative
/// bound.
pub fn surface_bbox(s: &BSplineSurface) -> ([f64; 3], [f64; 3]) {
    let mut mn = [f64::INFINITY; 3];
    let mut mx = [f64::NEG_INFINITY; 3];
    for row in &s.poles {
        for p in row {
            for c in 0..3 {
                mn[c] = mn[c].min(p[c]);
                mx[c] = mx[c].max(p[c]);
            }
        }
    }
    (mn, mx)
}

/// Validate the internal consistency of a surface: knot counts must match the
/// pole counts (`#knots = #poles + deg + 1` in each direction), knots must be
/// non-decreasing, and the weight grid (when present) must match the pole grid.
pub fn validate_bspline_surface(s: &BSplineSurface) -> Result<(), String> {
    let nu = s.poles.len();
    if nu == 0 {
        return Err("BSplineSurface: no u-poles".to_string());
    }
    let nv = s.poles[0].len();
    if nv == 0 {
        return Err("BSplineSurface: no v-poles".to_string());
    }
    if s.poles.iter().any(|row| row.len() != nv) {
        return Err("BSplineSurface: ragged pole grid".to_string());
    }
    if s.knots_u.len() != nu + s.deg_u + 1 {
        return Err("BSplineSurface: u knot count does not match poles".to_string());
    }
    if s.knots_v.len() != nv + s.deg_v + 1 {
        return Err("BSplineSurface: v knot count does not match poles".to_string());
    }
    if s.knots_u.windows(2).any(|w| w[1] < w[0]) || s.knots_v.windows(2).any(|w| w[1] < w[0]) {
        return Err("BSplineSurface: knot vector not non-decreasing".to_string());
    }
    if let Some(ws) = &s.weights {
        if ws.len() != nu || ws.iter().any(|row| row.len() != nv) {
            return Err("BSplineSurface: weight grid does not match poles".to_string());
        }
    }
    Ok(())
}

/// Interpolate a `nu`×`nv` grid of 3-D points with a tensor-product B-spline
/// surface of degrees `deg_u`, `deg_v`.
///
/// The interpolation is performed in two passes of 1-D cubic B-spline
/// interpolation: first each u-row is interpolated along v to obtain
/// intermediate control poles, then each intermediate v-column is interpolated
/// along u. The resulting surface passes through every grid point exactly at
/// the parameters `u = i/(nu−1)`, `v = j/(nv−1)` (matching
/// [`eval_grid_surface`]). For `deg_u = deg_v = 3` this is the standard cubic
/// B-spline interpolation; higher degrees use the same clamped collocation
/// solve.
///
/// # Errors
/// Returns an error when the grid dimensions do not match `nu`×`nv`, when a
/// direction has fewer poles than `deg + 1`, or when a collocation system is
/// singular.
pub fn interpolate_grid(
    points: &[&[[f64; 3]]],
    nu: usize,
    nv: usize,
    deg_u: usize,
    deg_v: usize,
) -> Result<BSplineSurface, String> {
    if points.len() != nu {
        return Err(format!("interpolate_grid: expected {nu} rows, got {}", points.len()));
    }
    for (i, row) in points.iter().enumerate() {
        if row.len() != nv {
            return Err(format!("interpolate_grid: row {i} expected {nv} points, got {}", row.len()));
        }
    }
    if nu < deg_u + 1 || nv < deg_v + 1 {
        return Err(format!(
            "interpolate_grid: grid {nu}×{nv} must be at least {}×{}",
            deg_u + 1,
            deg_v + 1
        ));
    }

    // Pass 1: interpolate each u-row along v.
    let mut intermediate = vec![vec![[0.0; 3]; nv]; nu];
    for i in 0..nu {
        let column: Vec<[f64; 3]> = (0..nv).map(|j| points[i][j]).collect();
        let poles_v = interpolate_curve_1d(&column, deg_v)?;
        for j in 0..nv {
            intermediate[i][j] = poles_v[j];
        }
    }

    // Pass 2: interpolate each intermediate v-column along u.
    let mut poles = vec![vec![[0.0; 3]; nv]; nu];
    for j in 0..nv {
        let row: Vec<[f64; 3]> = (0..nu).map(|i| intermediate[i][j]).collect();
        let poles_u = interpolate_curve_1d(&row, deg_u)?;
        for i in 0..nu {
            poles[i][j] = poles_u[i];
        }
    }

    let (knots_u, knots_v) = bspline_surface_uniform_knots(nu, nv, deg_u, deg_v);
    Ok(BSplineSurface {
        knots_u,
        knots_v,
        deg_u,
        deg_v,
        poles,
        weights: None,
    })
}

/// 1-D clamped B-spline interpolation of `points` (returning the control
/// poles), using uniform parameters `t_i = i/(n−1)` on `[0, 1]`.
fn interpolate_curve_1d(points: &[[f64; 3]], deg: usize) -> Result<Vec<[f64; 3]>, String> {
    let n = points.len();
    if n < deg + 1 {
        return Err(format!("interpolate_curve_1d: need at least {} points", deg + 1));
    }
    let knots = uniform_clamped_knots(n, deg);
    let mut colloc = MathMatrix::new(1, n, 1, n);
    for i in 0..n {
        let t = if n == 1 { 0.0 } else { i as f64 / (n - 1) as f64 };
        let span = find_span(&knots, deg, t);
        let vals = basis_vals(span, deg, t, &knots);
        for r in 0..=deg {
            let j = span - deg + r;
            colloc.set_value(i + 1, j + 1, vals[r]);
        }
    }
    let mut poles = vec![[0.0; 3]; n];
    for coord in 0..3 {
        let rhs = MathVector::from_slice(&points.iter().map(|p| p[coord]).collect::<Vec<_>>());
        let sol = colloc.solve(&rhs).map_err(|e| e.to_string())?;
        for i in 0..n {
            poles[i][coord] = sol.value(i + 1);
        }
    }
    Ok(poles)
}

/// Sample the surface at `nu×nv` parameters spanning `[0, 1]` in each
/// direction. The sample parameters are the grid nodes `i/(nu−1)` and
/// `j/(nv−1)`, so for a surface built by [`interpolate_grid`] the sampled
/// points reproduce the input grid exactly.
pub fn eval_grid_surface(s: &BSplineSurface, nu: usize, nv: usize) -> Vec<[f64; 3]> {
    let mut out = Vec::with_capacity(nu * nv);
    for i in 0..nu {
        let u = if nu == 1 { 0.5 } else { i as f64 / (nu - 1) as f64 };
        for j in 0..nv {
            let v = if nv == 1 { 0.5 } else { j as f64 / (nv - 1) as f64 };
            out.push(eval_bspline_surface(s, u, v));
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build an `nu×nv` grid of 3-D points `(x, y, f(x, y))` on `[0, 1]²`.
    fn analytic_grid(nu: usize, nv: usize, f: &dyn Fn(f64, f64) -> f64) -> Vec<Vec<[f64; 3]>> {
        let mut g = vec![vec![[0.0; 3]; nv]; nu];
        for i in 0..nu {
            let x = i as f64 / (nu - 1) as f64;
            for j in 0..nv {
                let y = j as f64 / (nv - 1) as f64;
                g[i][j] = [x, y, f(x, y)];
            }
        }
        g
    }

    fn rows_of(g: &[Vec<[f64; 3]>]) -> Vec<&[[f64; 3]]> {
        g.iter().map(|r| r.as_slice()).collect()
    }

    #[test]
    fn uniform_knot_sanity() {
        let (ku, kv) = bspline_surface_uniform_knots(5, 6, 3, 2);
        // #knots = #poles + deg + 1 in each direction.
        assert_eq!(ku.len(), 5 + 3 + 1);
        assert_eq!(kv.len(), 6 + 2 + 1);
        assert_eq!(bspline_surface_degree(3, 2, 5, 6), 9 + 9);
        // Non-decreasing and clamped at both ends.
        assert!(ku.windows(2).all(|w| w[0] <= w[1]));
        assert!(kv.windows(2).all(|w| w[0] <= w[1]));
        assert_eq!(ku[0], 0.0);
        assert_eq!(*ku.last().unwrap(), 1.0);
        assert_eq!(kv[0], 0.0);
        assert_eq!(*kv.last().unwrap(), 1.0);
    }

    #[test]
    fn eval_matches_pole_at_clamped_knot() {
        let s = BSplineSurface {
            knots_u: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
            deg_u: 3,
            deg_v: 3,
            poles: vec![
                vec![[0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 2.0, 0.0], [0.0, 3.0, 0.0]],
                vec![[1.0, 0.0, 1.0], [1.0, 1.0, 1.0], [1.0, 2.0, 1.0], [1.0, 3.0, 1.0]],
                vec![[2.0, 0.0, 2.0], [2.0, 1.0, 2.0], [2.0, 2.0, 2.0], [2.0, 3.0, 2.0]],
                vec![[3.0, 0.0, 3.0], [3.0, 1.0, 3.0], [3.0, 2.0, 3.0], [3.0, 3.0, 3.0]],
            ],
            weights: None,
        };
        // u = 0 selects the first pole row, v = 0 the first column.
        let p = eval_bspline_surface(&s, 0.0, 0.0);
        assert!((p[0] - 0.0).abs() < 1e-12);
        assert!((p[1] - 0.0).abs() < 1e-12);
        assert!((p[2] - 0.0).abs() < 1e-12);
        let p = eval_bspline_surface(&s, 0.0, 1.0);
        assert!((p[1] - 3.0).abs() < 1e-12);
        let p = eval_bspline_surface(&s, 1.0, 0.0);
        assert!((p[0] - 3.0).abs() < 1e-12);
        // Interior evaluation is finite.
        let p = eval_bspline_surface(&s, 0.5, 0.5);
        assert!(p.iter().all(|x| x.is_finite()));
    }

    #[test]
    fn interpolate_grid_analytic() {
        // z = x·y on [0,1]² — a bilinear surface reproduced exactly by cubic
        // tensor-product splines, so the grid samples are within 1e-6.
        let g = analytic_grid(4, 4, &|x, y| x * y);
        let rows = rows_of(&g);
        let s = interpolate_grid(&rows, 4, 4, 3, 3).unwrap();
        assert!(validate_bspline_surface(&s).is_ok());
        let sampled = eval_grid_surface(&s, 4, 4);
        let mut idx = 0;
        for i in 0..4 {
            let x = i as f64 / 3.0;
            for j in 0..4 {
                let y = j as f64 / 3.0;
                let got = sampled[idx];
                let want = [x, y, x * y];
                for c in 0..3 {
                    assert!(
                        (got[c] - want[c]).abs() < 1e-6,
                        "at ({i},{j}) c{c}: got {} want {}",
                        got[c],
                        want[c]
                    );
                }
                idx += 1;
            }
        }
    }

    #[test]
    fn finite_difference_derivative_is_finite() {
        let g = analytic_grid(5, 5, &|x, y| x * x + y * y);
        let rows = rows_of(&g);
        let s = interpolate_grid(&rows, 5, 5, 3, 3).unwrap();
        let h = 1e-6;
        let (u, v) = (0.37, 0.61);
        let p1 = eval_bspline_surface(&s, u + h, v);
        let p2 = eval_bspline_surface(&s, u - h, v);
        let du = [(p1[0] - p2[0]) / (2.0 * h), (p1[1] - p2[1]) / (2.0 * h), (p1[2] - p2[2]) / (2.0 * h)];
        assert!(du.iter().all(|x| x.is_finite()), "du = {du:?}");
        let q1 = eval_bspline_surface(&s, u, v + h);
        let q2 = eval_bspline_surface(&s, u, v - h);
        let dv = [(q1[0] - q2[0]) / (2.0 * h), (q1[1] - q2[1]) / (2.0 * h), (q1[2] - q2[2]) / (2.0 * h)];
        assert!(dv.iter().all(|x| x.is_finite()), "dv = {dv:?}");
    }

    #[test]
    fn corner_interpolation_exact() {
        let g = analytic_grid(4, 4, &|x, y| x * y);
        let rows = rows_of(&g);
        let s = interpolate_grid(&rows, 4, 4, 3, 3).unwrap();
        for (u, v, want) in [
            (0.0, 0.0, [0.0, 0.0, 0.0]),
            (1.0, 0.0, [1.0, 0.0, 0.0]),
            (0.0, 1.0, [0.0, 1.0, 0.0]),
            (1.0, 1.0, [1.0, 1.0, 1.0]),
        ] {
            let got = eval_bspline_surface(&s, u, v);
            for c in 0..3 {
                assert!(
                    (got[c] - want[c]).abs() < 1e-10,
                    "corner ({u},{v}) c{c}: got {} want {}",
                    got[c],
                    want[c]
                );
            }
        }
    }

    #[test]
    fn bump_grid_5x5() {
        // A smooth bump z = sin(x·y) on a 5×5 grid.
        let g = analytic_grid(5, 5, &|x, y| (x * y).sin());
        let rows = rows_of(&g);
        let s = interpolate_grid(&rows, 5, 5, 3, 3).unwrap();
        let sampled = eval_grid_surface(&s, 5, 5);
        let mut idx = 0;
        for i in 0..5 {
            for j in 0..5 {
                let x = i as f64 / 4.0;
                let y = j as f64 / 4.0;
                let got = sampled[idx];
                let want = [x, y, (x * y).sin()];
                for c in 0..3 {
                    assert!(
                        (got[c] - want[c]).abs() < 1e-8,
                        "bump ({i},{j}) c{c}: got {} want {}",
                        got[c],
                        want[c]
                    );
                }
                idx += 1;
            }
        }
    }

    #[test]
    fn rational_weights_match_unweighted_when_one() {
        let g = analytic_grid(4, 4, &|x, y| x * y);
        let rows = rows_of(&g);
        let mut s = interpolate_grid(&rows, 4, 4, 3, 3).unwrap();
        let nu = s.poles.len();
        let nv = s.poles[0].len();
        let all_ones = Some(vec![vec![1.0; nv]; nu]);
        s.weights = all_ones.clone();
        for (u, v) in [(0.0, 0.0), (0.3, 0.7), (1.0, 1.0), (0.5, 0.25)] {
            let a = eval_bspline_surface(&s, u, v);
            s.weights = None;
            let b = eval_bspline_surface(&s, u, v);
            s.weights = all_ones.clone();
            for c in 0..3 {
                assert!((a[c] - b[c]).abs() < 1e-12, "weight=1 mismatch at ({u},{v})");
            }
        }
    }
}
