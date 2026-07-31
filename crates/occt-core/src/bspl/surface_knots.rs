//! B-spline surface knot operations: seam closing, periodic wrapping,
//! trimming, reversing, and transposing the control-point grid.

use crate::gp::GpPnt;
use super::knots;
use super::surface;

/// Average the first and last v-pole of each u-row so a periodic v-surface
/// closes smoothly at its seam. Returns the modified pole grid.
pub fn compute_seam_poles(poles: &[GpPnt], n_u: usize, n_v: usize) -> Vec<GpPnt> {
    let mut result = poles.to_vec();
    if n_u == 0 || n_v == 0 {
        return result;
    }
    for i in 0..n_u {
        let a = i * n_v;
        let b = i * n_v + n_v - 1;
        let mid = GpPnt::new(
            0.5 * (poles[a].x() + poles[b].x()),
            0.5 * (poles[a].y() + poles[b].y()),
            0.5 * (poles[a].z() + poles[b].z()),
        );
        result[a] = mid;
        result[b] = mid;
    }
    result
}

/// Make the surface periodic in u by wrapping the last `degree_u + 1` rows onto
/// the front of the grid and extending the u-knot vector by one period.
/// Returns (poles, knots_u).
pub fn make_periodic_u(
    poles: &[GpPnt],
    n_u: usize,
    n_v: usize,
    degree_u: usize,
    knots_u: &[f64],
) -> (Vec<GpPnt>, Vec<f64>) {
    let wrap = degree_u + 1;
    if n_u == 0 || knots_u.is_empty() {
        return (poles.to_vec(), knots_u.to_vec());
    }
    let period = knots_u[knots_u.len() - 1] - knots_u[0];
    let mut new_knots = Vec::with_capacity(knots_u.len() + wrap);
    for &k in knots_u.iter().take(wrap.min(knots_u.len())) {
        new_knots.push(k - period);
    }
    new_knots.extend_from_slice(knots_u);

    let start = n_u.saturating_sub(wrap);
    let mut new_poles = Vec::with_capacity((n_u + wrap) * n_v);
    for i in start..n_u {
        new_poles.extend_from_slice(&poles[i * n_v..(i + 1) * n_v]);
    }
    new_poles.extend_from_slice(poles);
    (new_poles, new_knots)
}

/// Make the surface periodic in v by wrapping the last `degree_v + 1` columns
/// onto the front of each u-row and extending the v-knot vector by one period.
/// Returns (poles, knots_v).
pub fn make_periodic_v(
    poles: &[GpPnt],
    n_u: usize,
    n_v: usize,
    degree_v: usize,
    knots_v: &[f64],
) -> (Vec<GpPnt>, Vec<f64>) {
    let wrap = degree_v + 1;
    if n_v == 0 || knots_v.is_empty() {
        return (poles.to_vec(), knots_v.to_vec());
    }
    let period = knots_v[knots_v.len() - 1] - knots_v[0];
    let mut new_knots = Vec::with_capacity(knots_v.len() + wrap);
    for &k in knots_v.iter().take(wrap.min(knots_v.len())) {
        new_knots.push(k - period);
    }
    new_knots.extend_from_slice(knots_v);

    let start = n_v.saturating_sub(wrap);
    let new_n_v = n_v + wrap;
    let mut new_poles = Vec::with_capacity(n_u * new_n_v);
    for i in 0..n_u {
        for j in start..n_v {
            new_poles.push(poles[i * n_v + j]);
        }
        for j in 0..n_v {
            new_poles.push(poles[i * n_v + j]);
        }
    }
    (new_poles, new_knots)
}

/// Last index of `u` in the knot vector `ks`, or 0 when absent.
fn last_index(ks: &[f64], u: f64) -> usize {
    ks.iter().rposition(|&k| (k - u).abs() < 1e-15).unwrap_or(0)
}

/// Trim the surface to the sub-rectangle [u1,u2] x [v1,v2] by inserting knots
/// at u1, u2 (u-direction) and v1, v2 (v-direction) up to full multiplicity and
/// extracting the enclosed control grid. Returns (poles, knots_u, knots_v).
pub fn trim_surface(
    poles: &[GpPnt],
    knots_u: &[f64],
    knots_v: &[f64],
    degree_u: usize,
    degree_v: usize,
    u1: f64,
    u2: f64,
    v1: f64,
    v2: f64,
) -> (Vec<GpPnt>, Vec<f64>, Vec<f64>) {
    let mut n_u = knots_u.len().saturating_sub(degree_u + 1);
    let mut n_v = knots_v.len().saturating_sub(degree_v + 1);
    let mut p = poles.to_vec();
    let mut ku = knots_u.to_vec();
    let mut kv = knots_v.to_vec();

    // Bring u1 and u2 up to full multiplicity (degree_u + 1), one knot at a
    // time so the knot vector stays consistent between insertions.
    for u in [u1, u2] {
        let extra = (degree_u + 1).saturating_sub(knots::multiplicity(&ku, u));
        for _ in 0..extra {
            let (np, nk) = surface::insert_knot_u(&p, n_u, n_v, &ku, &kv, degree_u, degree_v, u, 1);
            p = np;
            ku = nk;
            n_u += 1;
        }
    }
    // Same for v1 and v2 in the v-direction.
    for v in [v1, v2] {
        let extra = (degree_v + 1).saturating_sub(knots::multiplicity(&kv, v));
        for _ in 0..extra {
            let (np, nk) = surface::insert_knot_v(&p, n_u, n_v, &ku, &kv, degree_u, degree_v, v, 1);
            p = np;
            kv = nk;
            n_v += 1;
        }
    }

    // The sub-surface over [u1,u2] x [v1,v2] is bounded by the full-multiplicity
    // knots: control rows/cols a..=b and knots a..=b + degree + 1.
    let lu1 = last_index(&ku, u1);
    let lu2 = last_index(&ku, u2);
    let lv1 = last_index(&kv, v1);
    let lv2 = last_index(&kv, v2);

    let a_u = lu1.saturating_sub(degree_u);
    let b_u = lu2.saturating_sub(degree_u + 1);
    let a_v = lv1.saturating_sub(degree_v);
    let b_v = lv2.saturating_sub(degree_v + 1);

    let eu = (b_u + degree_u + 1).min(ku.len().saturating_sub(1));
    let ev = (b_v + degree_v + 1).min(kv.len().saturating_sub(1));
    let new_knots_u = if a_u <= eu { ku[a_u..=eu].to_vec() } else { Vec::new() };
    let new_knots_v = if a_v <= ev { kv[a_v..=ev].to_vec() } else { Vec::new() };

    let mut new_poles = Vec::new();
    for i in a_u..=b_u {
        for j in a_v..=b_v {
            new_poles.push(p[i * n_v + j]);
        }
    }
    (new_poles, new_knots_u, new_knots_v)
}

/// Reverse the u-rows of the pole grid (row-major: n_u rows x n_v columns).
pub fn reverse_u(poles: &mut [GpPnt], n_u: usize, n_v: usize) {
    for i in 0..n_u / 2 {
        for j in 0..n_v {
            poles.swap(i * n_v + j, (n_u - 1 - i) * n_v + j);
        }
    }
}

/// Reverse the v-columns of the pole grid (within every u-row).
pub fn reverse_v(poles: &mut [GpPnt], n_u: usize, n_v: usize) {
    for i in 0..n_u {
        for j in 0..n_v / 2 {
            poles.swap(i * n_v + j, i * n_v + (n_v - 1 - j));
        }
    }
}

/// Transpose the pole grid (swap u and v directions). Row-major in, row-major
/// out with n_v rows x n_u columns.
pub fn exchange_u_v(poles: &[GpPnt], n_u: usize, n_v: usize) -> Vec<GpPnt> {
    let mut result = vec![GpPnt::zero(); n_u * n_v];
    for i in 0..n_u {
        for j in 0..n_v {
            result[j * n_u + i] = poles[i * n_v + j];
        }
    }
    result
}

/// Guess a (n_u, n_v) pole-grid split from the total pole count using sqrt.
pub fn pole_grid_dimensions(n_total: usize) -> (usize, usize) {
    if n_total == 0 {
        return (0, 0);
    }
    let n_u = (n_total as f64).sqrt().ceil() as usize;
    let n_v = (n_total + n_u - 1) / n_u;
    (n_u, n_v)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn grid(rows: usize, cols: usize) -> Vec<GpPnt> {
        let mut v = Vec::with_capacity(rows * cols);
        for i in 0..rows {
            for j in 0..cols {
                v.push(GpPnt::new(i as f64, j as f64, 0.0));
            }
        }
        v
    }

    #[test]
    fn reverse_u_reverses_rows() {
        let mut p = grid(3, 2);
        reverse_u(&mut p, 3, 2);
        assert_eq!(p[0].x(), 2.0);
        assert_eq!(p[1].y(), 1.0);
        assert_eq!(p[2].x(), 1.0);
        assert_eq!(p[4].x(), 0.0);
        assert_eq!(p[5].y(), 1.0);
    }

    #[test]
    fn exchange_u_v_transposes() {
        let p = grid(2, 3);
        let t = exchange_u_v(&p, 2, 3);
        assert_eq!(t.len(), 6);
        assert_eq!(t[0].x(), 0.0);
        assert_eq!(t[0].y(), 0.0);
        assert_eq!(t[1].x(), 1.0);
        assert_eq!(t[1].y(), 0.0);
        assert_eq!(t[2].x(), 0.0);
        assert_eq!(t[2].y(), 1.0);
        assert_eq!(t[3].x(), 1.0);
        assert_eq!(t[3].y(), 1.0);
        assert_eq!(t[4].x(), 0.0);
        assert_eq!(t[4].y(), 2.0);
        assert_eq!(t[5].x(), 1.0);
        assert_eq!(t[5].y(), 2.0);
    }

    #[test]
    fn pole_grid_dimensions_guesses_sqrt() {
        assert_eq!(pole_grid_dimensions(16), (4, 4));
        assert_eq!(pole_grid_dimensions(12), (4, 3));
        assert_eq!(pole_grid_dimensions(9), (3, 3));
    }

    #[test]
    fn trim_surface_bilinear() {
        let poles = vec![
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(0.0, 1.0, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
            GpPnt::new(1.0, 1.0, 0.0),
        ];
        let ku = vec![0.0, 0.0, 1.0, 1.0];
        let kv = vec![0.0, 0.0, 1.0, 1.0];
        let (tp, nku, nkv) = trim_surface(&poles, &ku, &kv, 1, 1, 0.25, 0.75, 0.25, 0.75);
        assert_eq!(nku, vec![0.25, 0.25, 0.75, 0.75]);
        assert_eq!(nkv, vec![0.25, 0.25, 0.75, 0.75]);
        assert_eq!(tp.len(), 4);
        assert!((tp[0].x() - 0.25).abs() < 1e-12);
        assert!((tp[0].y() - 0.25).abs() < 1e-12);
        assert!((tp[1].y() - 0.75).abs() < 1e-12);
        assert!((tp[2].x() - 0.75).abs() < 1e-12);
        assert!((tp[3].x() - 0.75).abs() < 1e-12);
        assert!((tp[3].y() - 0.75).abs() < 1e-12);
    }
}
