//! Grid interpolation of `Geom_BSplineSurface` (`fit_surface_grid`).

use super::eval::basis_funs;
use super::*;
use occt_core::gp::GpPnt;

/// Collocation matrix `A[i][k] = N_k(u_i)` for square interpolation solves.
fn collocation_matrix(knots: &[f64], degree: usize, params: &[f64]) -> Vec<Vec<f64>> {
    let n = params.len();
    let mut a = vec![vec![0.0; n]; n];
    for (i, &u) in params.iter().enumerate() {
        // Interpolation nodes are inside the knot domain, so the periodic
        // wrap of `BSplCLib::LocateParameter` is a no-op here.
        let b = basis_funs(knots, degree, u, false);
        for k in 0..n {
            a[i][k] = b[k];
        }
    }
    a
}

/// Small square solve via Gaussian elimination with partial pivoting.
fn solve_linear_small(a: &[Vec<f64>], b: &[f64]) -> Result<Vec<f64>, String> {
    let n = b.len();
    let mut m = a.to_vec();
    let mut rhs = b.to_vec();
    for col in 0..n {
        let mut piv = col;
        for r in (col + 1)..n {
            if m[r][col].abs() > m[piv][col].abs() {
                piv = r;
            }
        }
        if m[piv][col].abs() < 1e-30 {
            return Err(format!("fit_surface_grid: singular collocation at column {col}"));
        }
        m.swap(col, piv);
        rhs.swap(col, piv);
        let d = m[col][col];
        for r in (col + 1)..n {
            let f = m[r][col] / d;
            for c in col..n {
                m[r][c] -= f * m[col][c];
            }
            rhs[r] -= f * rhs[col];
        }
    }
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let mut s = rhs[i];
        for j in (i + 1)..n {
            s -= m[i][j] * x[j];
        }
        x[i] = s / m[i][i];
    }
    Ok(x)
}

/// Interpolate an `nu × nv` grid of points with a clamped uniform B-spline
/// surface of degrees `deg_u × deg_v`, passing exactly through every grid node.
///
/// Degree 1 uses the grid points directly as poles. Higher degrees use the
/// two-pass 1-D B-spline collocation solve (per column, then per row), so the
/// result is an interpolant through all nodes.
pub fn fit_surface_grid(
    points: &[Vec<GpPnt>],
    deg_u: usize,
    deg_v: usize,
) -> Result<GeomBSplineSurface, String> {
    let nu = points.len();
    if nu == 0 {
        return Err("fit_surface_grid: empty grid".to_string());
    }
    let nv = points[0].len();
    if nv == 0 {
        return Err("fit_surface_grid: zero-width grid".to_string());
    }
    for row in points {
        if row.len() != nv {
            return Err("fit_surface_grid: ragged grid".to_string());
        }
    }
    if nu < deg_u + 1 || nv < deg_v + 1 {
        return Err("fit_surface_grid: grid too small for the requested degree".to_string());
    }
    let (knots_u, knots_v) = bspline_surface_uniform_knots(nu, nv, deg_u, deg_v);
    if deg_u == 1 && deg_v == 1 {
        return Ok(GeomBSplineSurface::from_parts(
            points.to_vec(),
            None,
            knots_u,
            knots_v,
            deg_u,
            deg_v,
            false,
            false,
        ));
    }

    // Uniform grid parameters in [0, 1].
    let us: Vec<f64> = (0..nu).map(|i| i as f64 / (nu - 1) as f64).collect();
    let vs: Vec<f64> = (0..nv).map(|j| j as f64 / (nv - 1) as f64).collect();
    let au = collocation_matrix(&knots_u, deg_u, &us);
    let av = collocation_matrix(&knots_v, deg_v, &vs);

    // Pass 1: interpolate each v-column in u -> intermediate lattice C.
    let mut c = vec![vec![GpPnt::zero(); nv]; nu];
    for j in 0..nv {
        let mut bx = vec![0.0; nu];
        let mut by = vec![0.0; nu];
        let mut bz = vec![0.0; nu];
        for i in 0..nu {
            bx[i] = points[i][j].x();
            by[i] = points[i][j].y();
            bz[i] = points[i][j].z();
        }
        let cx = solve_linear_small(&au, &bx)?;
        let cy = solve_linear_small(&au, &by)?;
        let cz = solve_linear_small(&au, &bz)?;
        for i in 0..nu {
            c[i][j] = GpPnt::new(cx[i], cy[i], cz[i]);
        }
    }
    // Pass 2: interpolate each u-row in v -> final poles.
    let mut poles = vec![vec![GpPnt::zero(); nv]; nu];
    for i in 0..nu {
        let mut bx = vec![0.0; nv];
        let mut by = vec![0.0; nv];
        let mut bz = vec![0.0; nv];
        for j in 0..nv {
            bx[j] = c[i][j].x();
            by[j] = c[i][j].y();
            bz[j] = c[i][j].z();
        }
        let px = solve_linear_small(&av, &bx)?;
        let py = solve_linear_small(&av, &by)?;
        let pz = solve_linear_small(&av, &bz)?;
        for j in 0..nv {
            poles[i][j] = GpPnt::new(px[j], py[j], pz[j]);
        }
    }
    Ok(GeomBSplineSurface::from_parts(
        poles,
        None,
        knots_u,
        knots_v,
        deg_u,
        deg_v,
        false,
        false,
    ))
}
