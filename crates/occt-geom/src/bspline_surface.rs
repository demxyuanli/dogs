//! B-spline surface (tensor-product, clamped knots). Source: `Geom_BSplineSurface.hxx`
//!
//! Poles are stored 0-based as a `u × v` grid; weights (if any) match the
//! pole grid element-wise. Evaluation uses Cox-de Boor basis functions in each
//! parametric direction with a rational division when weights are present.

use occt_core::gp::{GpPnt, GpTrsf, GpVec};

use crate::surface::Surface;

/// Non-rational or rational tensor-product B-spline surface.
#[derive(Clone)]
pub struct GeomBSplineSurface {
    pub poles: Vec<Vec<GpPnt>>,
    pub knots_u: Vec<f64>,
    pub knots_v: Vec<f64>,
    pub deg_u: usize,
    pub deg_v: usize,
    pub weights: Option<Vec<Vec<f64>>>,
}

impl GeomBSplineSurface {
    /// Build a non-rational surface. Knot counts must be `#poles + degree + 1`.
    pub fn new(
        poles: Vec<Vec<GpPnt>>,
        knots_u: Vec<f64>,
        knots_v: Vec<f64>,
        deg_u: usize,
        deg_v: usize,
    ) -> Result<Self, String> {
        if deg_u < 1 || deg_v < 1 {
            return Err("GeomBSplineSurface::new: degrees must be >= 1".to_string());
        }
        let nu = poles.len();
        if nu == 0 {
            return Err("GeomBSplineSurface::new: empty pole grid".to_string());
        }
        let nv = poles[0].len();
        if nv == 0 {
            return Err("GeomBSplineSurface::new: zero-width pole grid".to_string());
        }
        for row in &poles {
            if row.len() != nv {
                return Err("GeomBSplineSurface::new: ragged pole grid".to_string());
            }
        }
        if knots_u.len() != nu + deg_u + 1 {
            return Err(format!(
                "GeomBSplineSurface::new: knots_u count {} != poles_u {nu} + deg {deg_u} + 1",
                knots_u.len()
            ));
        }
        if knots_v.len() != nv + deg_v + 1 {
            return Err(format!(
                "GeomBSplineSurface::new: knots_v count {} != poles_v {nv} + deg {deg_v} + 1",
                knots_v.len()
            ));
        }
        if !knots_u.windows(2).all(|w| w[0] <= w[1]) {
            return Err("GeomBSplineSurface::new: knots_u not non-decreasing".to_string());
        }
        if !knots_v.windows(2).all(|w| w[0] <= w[1]) {
            return Err("GeomBSplineSurface::new: knots_v not non-decreasing".to_string());
        }
        Ok(Self { poles, knots_u, knots_v, deg_u, deg_v, weights: None })
    }

    /// Build a rational surface; `weights` must match the pole grid size.
    pub fn rational(
        poles: Vec<Vec<GpPnt>>,
        weights: Vec<Vec<f64>>,
        knots_u: Vec<f64>,
        knots_v: Vec<f64>,
        deg_u: usize,
        deg_v: usize,
    ) -> Result<Self, String> {
        let mut s = Self::new(poles, knots_u, knots_v, deg_u, deg_v)?;
        if s.poles.len() != weights.len() || s.poles[0].len() != weights[0].len() {
            return Err("GeomBSplineSurface::rational: weight grid mismatch".to_string());
        }
        s.weights = Some(weights);
        Ok(s)
    }

    pub fn nb_poles_u(&self) -> usize { self.poles.len() }
    pub fn nb_poles_v(&self) -> usize { self.poles.first().map_or(0, |r| r.len()) }
    pub fn is_rational(&self) -> bool { self.weights.is_some() }
}

/// Clamped uniform knot vector for `n` poles of degree `p`
/// (length `n + p + 1`, multiplicity `p + 1` at both ends).
fn clamped_uniform_knots(n: usize, p: usize) -> Vec<f64> {
    let len = n + p + 1;
    let interior = if n >= p + 1 { n - p - 1 } else { 0 };
    let mut k = Vec::with_capacity(len);
    for _ in 0..=p {
        k.push(0.0);
    }
    for i in 1..=interior {
        k.push(i as f64 / (interior + 1) as f64);
    }
    while k.len() < len {
        k.push(1.0);
    }
    k
}

/// Clamped uniform knot vectors for a surface with `nu × nv` poles.
pub fn bspline_surface_uniform_knots(
    nu: usize,
    nv: usize,
    deg_u: usize,
    deg_v: usize,
) -> (Vec<f64>, Vec<f64>) {
    (clamped_uniform_knots(nu, deg_u), clamped_uniform_knots(nv, deg_v))
}

/// Find the knot span index `s` with `knots[s] <= u <= knots[s + 1]`
/// (The NURBS Book algorithm A2.1, clamped to the valid range).
fn find_span(knots: &[f64], degree: usize, u: f64) -> usize {
    let n = knots.len() - degree - 2; // last control-point index
    if u >= knots[n + 1] {
        return n;
    }
    if u <= knots[degree] {
        return degree;
    }
    let mut low = degree;
    let mut high = n + 1;
    let mut mid = (low + high) / 2;
    while u < knots[mid] || u >= knots[mid + 1] {
        if u < knots[mid] {
            high = mid;
        } else {
            low = mid;
        }
        mid = (low + high) / 2;
    }
    mid
}

/// All non-zero degree-`degree` basis functions at `u` (NURBS A2.2), returned
/// as a full-length vector of length `#poles`.
fn basis_funs(knots: &[f64], degree: usize, u: f64) -> Vec<f64> {
    let n = knots.len() - degree - 2; // last control-point index
    let u = if u < knots[degree] {
        knots[degree]
    } else if u > knots[n + 1] {
        knots[n + 1]
    } else {
        u
    };
    let span = find_span(knots, degree, u);
    let mut nvec = vec![0.0; degree + 1];
    nvec[0] = 1.0;
    let mut left = vec![0.0; degree + 1];
    let mut right = vec![0.0; degree + 1];
    for j in 1..=degree {
        left[j] = u - knots[span + 1 - j];
        right[j] = knots[span + j] - u;
        let mut saved = 0.0;
        for r in 0..j {
            let denom = right[r + 1] + left[j - r];
            let temp = if denom.abs() > 1e-300 { nvec[r] / denom } else { 0.0 };
            nvec[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        nvec[j] = saved;
    }
    let mut out = vec![0.0; n + 1];
    for j in 0..=degree {
        let idx = span - degree + j;
        if idx <= n {
            out[idx] = nvec[j];
        }
    }
    out
}

/// Evaluate the tensor-product B-spline surface at `(u, v)`.
///
/// The rational case computes `Σ N_i(u)N_j(v) w_ij P_ij / Σ N_i(u)N_j(v) w_ij`;
/// the non-rational case is the same with unit weights.
pub fn eval_bspline_surface(s: &GeomBSplineSurface, u: f64, v: f64) -> GpPnt {
    let nu = s.poles.len();
    if nu == 0 {
        return GpPnt::zero();
    }
    let nv = s.poles[0].len();
    if nv == 0 {
        return GpPnt::zero();
    }
    let bu = basis_funs(&s.knots_u, s.deg_u, u);
    let bv = basis_funs(&s.knots_v, s.deg_v, v);

    let mut acc = [0.0_f64, 0.0, 0.0];
    match &s.weights {
        Some(w) => {
            let mut den = 0.0;
            for i in 0..nu {
                for j in 0..nv {
                    let b = bu[i] * bv[j] * w[i][j];
                    let p = s.poles[i][j];
                    acc[0] += b * p.x();
                    acc[1] += b * p.y();
                    acc[2] += b * p.z();
                    den += b;
                }
            }
            if den.abs() < 1e-300 {
                return GpPnt::zero();
            }
            GpPnt::new(acc[0] / den, acc[1] / den, acc[2] / den)
        }
        None => {
            for i in 0..nu {
                for j in 0..nv {
                    let b = bu[i] * bv[j];
                    let p = s.poles[i][j];
                    acc[0] += b * p.x();
                    acc[1] += b * p.y();
                    acc[2] += b * p.z();
                }
            }
            GpPnt::new(acc[0], acc[1], acc[2])
        }
    }
}

/// Collocation matrix `A[i][k] = N_k(u_i)` for square interpolation solves.
fn collocation_matrix(knots: &[f64], degree: usize, params: &[f64]) -> Vec<Vec<f64>> {
    let n = params.len();
    let mut a = vec![vec![0.0; n]; n];
    for (i, &u) in params.iter().enumerate() {
        let b = basis_funs(knots, degree, u);
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
        return Ok(GeomBSplineSurface {
            poles: points.to_vec(),
            knots_u,
            knots_v,
            deg_u,
            deg_v,
            weights: None,
        });
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
    Ok(GeomBSplineSurface { poles, knots_u, knots_v, deg_u, deg_v, weights: None })
}

impl Surface for GeomBSplineSurface {
    fn d0(&self, u: f64, v: f64) -> GpPnt {
        eval_bspline_surface(self, u, v)
    }

    fn d1(&self, u: f64, v: f64) -> (GpPnt, GpVec, GpVec) {
        let p = self.d0(u, v);
        let (u0, u1) = self.u_range();
        let (v0, v1) = self.v_range();
        let hu = if u1 > u0 { (u1 - u0) * 1e-4 } else { 1e-6 };
        let hv = if v1 > v0 { (v1 - v0) * 1e-4 } else { 1e-6 };
        let pu = self.d0(u + hu, v);
        let pm = self.d0(u - hu, v);
        let pv = self.d0(u, v + hv);
        let pn = self.d0(u, v - hv);
        let du = GpVec::new(
            (pu.x() - pm.x()) / (2.0 * hu),
            (pu.y() - pm.y()) / (2.0 * hu),
            (pu.z() - pm.z()) / (2.0 * hu),
        );
        let dv = GpVec::new(
            (pv.x() - pn.x()) / (2.0 * hv),
            (pv.y() - pn.y()) / (2.0 * hv),
            (pv.z() - pn.z()) / (2.0 * hv),
        );
        (p, du, dv)
    }

    fn u_range(&self) -> (f64, f64) {
        if self.knots_u.is_empty() || self.deg_u >= self.knots_u.len() {
            return (0.0, 1.0);
        }
        let lo = self.knots_u[self.deg_u];
        let hi = self.knots_u[self.knots_u.len() - 1 - self.deg_u];
        (lo, hi)
    }

    fn v_range(&self) -> (f64, f64) {
        if self.knots_v.is_empty() || self.deg_v >= self.knots_v.len() {
            return (0.0, 1.0);
        }
        let lo = self.knots_v[self.deg_v];
        let hi = self.knots_v[self.knots_v.len() - 1 - self.deg_v];
        (lo, hi)
    }

    fn continuity(&self) -> u8 {
        self.deg_u.min(self.deg_v).min(u8::MAX as usize) as u8
    }

    fn transform(&mut self, t: &GpTrsf) {
        for row in self.poles.iter_mut() {
            for p in row.iter_mut() {
                *p = p.transformed(t);
            }
        }
    }

    fn clone_dyn(&self) -> Box<dyn Surface> {
        Box::new(self.clone())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn uniform_knot_sanity() {
        let (ku, kv) = bspline_surface_uniform_knots(4, 3, 2, 1);
        assert_eq!(ku.len(), 4 + 2 + 1);
        assert_eq!(kv.len(), 3 + 1 + 1);
        assert_eq!(ku[0], 0.0);
        assert_eq!(*ku.last().unwrap(), 1.0);
        assert!(ku.windows(2).all(|w| w[0] <= w[1]));
        assert!(kv.windows(2).all(|w| w[0] <= w[1]));
        // End multiplicity degree + 1.
        assert_eq!(ku[..3].iter().filter(|&&x| x == 0.0).count(), 3);
        assert_eq!(ku[ku.len() - 3..].iter().filter(|&&x| x == 1.0).count(), 3);
        assert_eq!(kv[..2].iter().filter(|&&x| x == 0.0).count(), 2);
        // Interior knot is midway for a 4-pole degree-2 vector.
        assert!((ku[3] - 0.5).abs() < 1e-12);
    }

    #[test]
    fn eval_at_clamped_corner_is_first_pole() {
        let poles = vec![
            vec![GpPnt::new(1.0, 2.0, 3.0), GpPnt::new(4.0, 5.0, 6.0)],
            vec![GpPnt::new(7.0, 8.0, 9.0), GpPnt::new(10.0, 11.0, 12.0)],
        ];
        let (ku, kv) = bspline_surface_uniform_knots(2, 2, 1, 1);
        let s = GeomBSplineSurface::new(poles, ku, kv, 1, 1).unwrap();
        let p = s.d0(0.0, 0.0);
        assert!((p.x() - 1.0).abs() < 1e-12);
        assert!((p.y() - 2.0).abs() < 1e-12);
        assert!((p.z() - 3.0).abs() < 1e-12);
        // Opposite corner: last pole.
        let q = s.d0(1.0, 1.0);
        assert!((q.x() - 10.0).abs() < 1e-12);
        assert!((q.y() - 11.0).abs() < 1e-12);
    }

    #[test]
    fn degree1_plane_grid() {
        // z = 2x + y + 1 sampled at the four corners.
        let poles = vec![
            vec![GpPnt::new(0.0, 0.0, 1.0), GpPnt::new(0.0, 1.0, 2.0)],
            vec![GpPnt::new(1.0, 0.0, 3.0), GpPnt::new(1.0, 1.0, 4.0)],
        ];
        let (ku, kv) = bspline_surface_uniform_knots(2, 2, 1, 1);
        let s = GeomBSplineSurface::new(poles, ku, kv, 1, 1).unwrap();
        let p = s.d0(0.5, 0.5);
        assert!((p.x() - 0.5).abs() < 1e-12);
        assert!((p.y() - 0.5).abs() < 1e-12);
        assert!((p.z() - 2.5).abs() < 1e-12, "z = {}", p.z());
    }

    #[test]
    fn xy_grid_fit_matches_analytic() {
        // Grid points on z = u·v; degree-1 surface reproduces the bilinear field.
        let (nu, nv) = (4, 4);
        let points: Vec<Vec<GpPnt>> = (0..nu)
            .map(|i| {
                (0..nv)
                    .map(|j| {
                        let u = i as f64 / (nu - 1) as f64;
                        let v = j as f64 / (nv - 1) as f64;
                        GpPnt::new(u, v, u * v)
                    })
                    .collect()
            })
            .collect();
        let s = fit_surface_grid(&points, 1, 1).unwrap();
        let p = s.d0(0.3, 0.7);
        assert!((p.x() - 0.3).abs() < 1e-6);
        assert!((p.y() - 0.7).abs() < 1e-6);
        assert!((p.z() - 0.21).abs() < 1e-6, "z = {}", p.z());
    }

    #[test]
    fn transform_translates_all_poles() {
        let poles = vec![
            vec![GpPnt::new(0.0, 0.0, 0.0), GpPnt::new(0.0, 1.0, 0.0)],
            vec![GpPnt::new(1.0, 0.0, 0.0), GpPnt::new(1.0, 1.0, 0.0)],
        ];
        let (ku, kv) = bspline_surface_uniform_knots(2, 2, 1, 1);
        let mut s = GeomBSplineSurface::new(poles, ku, kv, 1, 1).unwrap();
        let before = s.d0(0.4, 0.6);
        let mut trsf = GpTrsf::identity();
        trsf.set_translation_vec(&GpVec::new(1.0, 2.0, 3.0));
        s.transform(&trsf);
        let after = s.d0(0.4, 0.6);
        assert!((after.x() - before.x() - 1.0).abs() < 1e-12);
        assert!((after.y() - before.y() - 2.0).abs() < 1e-12);
        assert!((after.z() - before.z() - 3.0).abs() < 1e-12);
    }

    #[test]
    fn rational_weights_matches_weighted_formula() {
        // 2×2 degree-1 surface; at (0.5, 0.5) every basis product is 0.25.
        let poles = vec![
            vec![GpPnt::new(0.0, 0.0, 0.0), GpPnt::new(0.0, 1.0, 0.0)],
            vec![GpPnt::new(1.0, 0.0, 0.0), GpPnt::new(1.0, 1.0, 0.0)],
        ];
        let weights = vec![vec![1.0, 2.0], vec![3.0, 4.0]];
        let (ku, kv) = bspline_surface_uniform_knots(2, 2, 1, 1);
        let s = GeomBSplineSurface::rational(poles, weights, ku, kv, 1, 1).unwrap();
        let p = s.d0(0.5, 0.5);
        // Σ w P = (1·0+2·0+3·1+4·1, 1·0+2·1+3·0+4·1, 0) = (7, 6, 0), Σ w = 10.
        assert!((p.x() - 0.7).abs() < 1e-12, "x = {}", p.x());
        assert!((p.y() - 0.6).abs() < 1e-12, "y = {}", p.y());
        assert!(p.z().abs() < 1e-12);
        // Non-rational version differs (weights pull the point toward heavy poles).
        let s_nr = GeomBSplineSurface::new(
            vec![
                vec![GpPnt::new(0.0, 0.0, 0.0), GpPnt::new(0.0, 1.0, 0.0)],
                vec![GpPnt::new(1.0, 0.0, 0.0), GpPnt::new(1.0, 1.0, 0.0)],
            ],
            bspline_surface_uniform_knots(2, 2, 1, 1).0,
            bspline_surface_uniform_knots(2, 2, 1, 1).1,
            1,
            1,
        )
        .unwrap();
        let p_nr = s_nr.d0(0.5, 0.5);
        assert!((p_nr.x() - 0.5).abs() < 1e-12);
    }

    #[test]
    fn fit_surface_grid_passes_through_nodes() {
        // Cubic interpolation of z = u² + v³ over a 5×5 grid.
        let (nu, nv) = (5, 5);
        let points: Vec<Vec<GpPnt>> = (0..nu)
            .map(|i| {
                (0..nv)
                    .map(|j| {
                        let u = i as f64 / (nu - 1) as f64;
                        let v = j as f64 / (nv - 1) as f64;
                        GpPnt::new(u, v, u * u + v * v * v)
                    })
                    .collect()
            })
            .collect();
        let s = fit_surface_grid(&points, 3, 3).unwrap();
        for i in 0..nu {
            for j in 0..nv {
                let u = i as f64 / (nu - 1) as f64;
                let v = j as f64 / (nv - 1) as f64;
                let p = s.d0(u, v);
                assert!(
                    (p.z() - points[i][j].z()).abs() < 1e-6,
                    "node ({i},{j}) z = {} vs {}",
                    p.z(),
                    points[i][j].z()
                );
            }
        }
    }

    #[test]
    fn d1_nonzero_for_curved_surface() {
        let (nu, nv) = (5, 5);
        let points: Vec<Vec<GpPnt>> = (0..nu)
            .map(|i| {
                (0..nv)
                    .map(|j| {
                        let u = i as f64 / (nu - 1) as f64;
                        let v = j as f64 / (nv - 1) as f64;
                        GpPnt::new(u, v, u * u * v)
                    })
                    .collect()
            })
            .collect();
        let s = fit_surface_grid(&points, 3, 3).unwrap();
        let (p, du, dv) = s.d1(0.5, 0.5);
        assert!((p.z() - 0.125).abs() < 1e-4, "z = {}", p.z());
        assert!(du.magnitude() > 1e-6, "du should be non-zero");
        assert!(dv.magnitude() > 1e-6, "dv should be non-zero");
    }
}
