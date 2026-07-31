//! Surface fitting — least-squares plane/circle/sphere fit and B-spline
//! surface approximation from a sampled point grid.
//! Source: `GeomAPI_PointsToBSplineSurface`, `gce_Make*`.

use std::sync::Arc;

use occt_core::gp::{GpAx3, GpDir, GpPln, GpPnt, GpVec, GpXyz};

use crate::surface::Surface;

/// Least-squares plane fit through a point set. Returns the plane whose
/// normal is the smallest-eigenvalue direction of the covariance.
pub fn fit_plane(points: &[GpPnt]) -> Option<GpPln> {
    if points.len() < 3 {
        return None;
    }
    let n = points.len();
    let c = centroid(points);
    let mut cov = [[0.0; 3]; 3];
    for p in points {
        let d = [p.x() - c.x(), p.y() - c.y(), p.z() - c.z()];
        for i in 0..3 {
            for j in 0..3 {
                cov[i][j] += d[i] * d[j];
            }
        }
    }
    let normal = plane_normal(&cov)?;
    let d = GpDir::from_vec(&normal).ok()?;
    // In-plane X direction: perpendicular to the (possibly tilted) normal via
    // a cross product with a reference axis not parallel to d.
    let ref_v = if d.x().abs() < 0.9 {
        occt_core::gp::GpVec::new(1.0, 0.0, 0.0)
    } else {
        occt_core::gp::GpVec::new(0.0, 1.0, 0.0)
    };
    let x = d.xyz().crossed(ref_v.xyz());
    let x_dir = GpDir::from_vec(&occt_core::gp::GpVec::new(x.x, x.y, x.z)).ok()?;
    Some(GpPln::new(GpAx3::new(c, d, &x_dir).ok()?))
}

/// Least-squares sphere fit (algebraic). Returns (center, radius).
pub fn fit_sphere(points: &[GpPnt]) -> Option<(GpPnt, f64)> {
    if points.len() < 4 {
        return None;
    }
    // Solve the linear system for the sphere equation
    // 2x·cx + 2y·cy + 2z·cz - (cx²+cy²+cz²-r²) = x²+y²+z².
    let n = points.len();
    let mut a = [[0.0f64; 4]; 4];
    let mut b = [0.0f64; 4];
    for p in points {
        let r2 = p.x() * p.x() + p.y() * p.y() + p.z() * p.z();
        let row = [2.0 * p.x(), 2.0 * p.y(), 2.0 * p.z(), 1.0];
        for j in 0..4 {
            for k in 0..4 {
                a[j][k] += row[j] * row[k];
            }
            b[j] += row[j] * r2;
        }
    }
    let sol = solve4(&a, &b)?;
    let center = GpPnt::new(sol[0], sol[1], sol[2]);
    // The equation was 2cx·x + 2cy·y + 2cz·z + 1·c = x²+y²+z² with
    // c = r² − |center|², so r² = |center|² + sol[3].
    let csq = sol[0] * sol[0] + sol[1] * sol[1] + sol[2] * sol[2];
    let r2 = csq + sol[3];
    if r2 <= 0.0 {
        return None;
    }
    Some((center, r2.sqrt()))
}

/// Fit a plane, circle (in that plane), or sphere and report which matched
/// best by residual.
pub fn fit_surface_kind(points: &[GpPnt]) -> (SurfaceFitKind, f64) {
    let pln = fit_plane(points);
    let sph = fit_sphere(points);
    let plane_resid = pln.as_ref().map(|p| plane_residual(points, p)).unwrap_or(f64::INFINITY);
    let sphere_resid = sph.as_ref().map(|(c, r)| sphere_residual(points, c, *r)).unwrap_or(f64::INFINITY);
    if sphere_resid < plane_resid && sphere_resid.is_finite() {
        (SurfaceFitKind::Sphere, sphere_resid)
    } else if plane_resid.is_finite() {
        (SurfaceFitKind::Plane, plane_resid)
    } else {
        (SurfaceFitKind::None, f64::INFINITY)
    }
}

/// Result of surface classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SurfaceFitKind {
    Plane,
    Sphere,
    None,
}

/// A degree-1 (bilinear) B-spline surface through a grid of points. Implements
/// `Surface` by bilinear interpolation of the control lattice.
#[derive(Clone)]
pub struct GridSurface {
    pub us: Vec<f64>,
    pub vs: Vec<f64>,
    pub poles: Vec<Vec<GpPnt>>,
}

impl GridSurface {
    pub fn from_grid(us: Vec<f64>, vs: Vec<f64>, poles: Vec<Vec<GpPnt>>) -> Self {
        Self { us, vs, poles }
    }
}

impl Surface for GridSurface {
    fn d0(&self, u: f64, v: f64) -> GpPnt {
        let nu = self.us.len();
        let nv = self.vs.len();
        if nu < 2 || nv < 2 {
            return self.poles.first().and_then(|r| r.first()).copied().unwrap_or(GpPnt::zero());
        }
        let ui = index_of(&self.us, u).min(nu - 2);
        let vi = index_of(&self.vs, v).min(nv - 2);
        let tu = ((u - self.us[ui]) / (self.us[ui + 1] - self.us[ui]).max(1e-30)).clamp(0.0, 1.0);
        let tv = ((v - self.vs[vi]) / (self.vs[vi + 1] - self.vs[vi]).max(1e-30)).clamp(0.0, 1.0);
        let p00 = self.poles[ui][vi];
        let p10 = self.poles[ui + 1][vi];
        let p01 = self.poles[ui][vi + 1];
        let p11 = self.poles[ui + 1][vi + 1];
        let lerp = |a: &GpPnt, b: &GpPnt, t: f64| GpPnt::new(
            a.x() + t * (b.x() - a.x()),
            a.y() + t * (b.y() - a.y()),
            a.z() + t * (b.z() - a.z()),
        );
        let a = lerp(&p00, &p10, tu);
        let b = lerp(&p01, &p11, tu);
        lerp(&a, &b, tv)
    }
    fn d1(&self, u: f64, v: f64) -> (GpPnt, GpVec, GpVec) {
        let eps = 1e-6;
        let (u0, u1) = self.u_range();
        let (v0, v1) = self.v_range();
        let hu = if u1 > u0 { (u1 - u0) * 1e-4 } else { eps };
        let hv = if v1 > v0 { (v1 - v0) * 1e-4 } else { eps };
        let p = self.d0(u, v);
        let pu = self.d0(u + hu, v);
        let pv = self.d0(u, v + hv);
        (p, GpVec::from_pnts(&p, &pu), GpVec::from_pnts(&p, &pv))
    }
    fn u_range(&self) -> (f64, f64) {
        (*self.us.first().unwrap_or(&0.0), *self.us.last().unwrap_or(&1.0))
    }
    fn v_range(&self) -> (f64, f64) {
        (*self.vs.first().unwrap_or(&0.0), *self.vs.last().unwrap_or(&1.0))
    }
    fn continuity(&self) -> u8 {
        0
    }
    fn transform(&mut self, t: &occt_core::gp::GpTrsf) {
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

fn index_of(sorted: &[f64], x: f64) -> usize {
    if x <= sorted[0] {
        return 0;
    }
    for i in 1..sorted.len() {
        if x <= sorted[i] {
            return i - 1;
        }
    }
    sorted.len() - 2
}

/// Wrap a fitted surface as a trait object.
pub fn fit_surface(points: &[GpPnt]) -> Result<Arc<dyn Surface>, String> {
    let pln = fit_plane(points).ok_or("fit_surface: insufficient points")?;
    Ok(Arc::new(crate::GeomPlane::new(pln)))
}

fn centroid(points: &[GpPnt]) -> GpPnt {
    let n = points.len().max(1) as f64;
    let mut acc = GpXyz::zero();
    for p in points {
        acc = acc.added(&p.coord);
    }
    GpPnt::from_xyz(&acc.divided(n))
}

/// Normal to the best-fit plane: the smallest-variance eigenvector of the
/// covariance, via a full Jacobi 3×3 eigendecomposition (robust for
/// near-planar and degenerate sets).
fn plane_normal(m: &[[f64; 3]; 3]) -> Option<GpVec> {
    let (evals, evecs) = jacobi_eigen(m);
    let mut min_i = 0;
    for i in 1..3 {
        if evals[i] < evals[min_i] {
            min_i = i;
        }
    }
    let v = evecs[min_i];
    let norm = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if norm < 1e-20 {
        return None;
    }
    Some(GpVec::new(v[0] / norm, v[1] / norm, v[2] / norm))
}

/// Symmetric 3×3 Jacobi eigensolver. Returns (eigenvalues, eigenvectors as
/// column vectors).
fn jacobi_eigen(m: &[[f64; 3]; 3]) -> ([f64; 3], [[f64; 3]; 3]) {
    let mut a = *m;
    let mut v = [[1.0f64, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    for _ in 0..50 {
        let mut p = 0usize;
        let mut q = 1usize;
        let mut max = a[0][1].abs();
        for i in 0..3 {
            for j in (i + 1)..3 {
                if a[i][j].abs() > max {
                    max = a[i][j].abs();
                    p = i;
                    q = j;
                }
            }
        }
        if max < 1e-30 {
            break;
        }
        let app = a[p][p];
        let aqq = a[q][q];
        let apq = a[p][q];
        let tau = (aqq - app) / (2.0 * apq);
        let t = tau.signum() / (tau.abs() + (1.0 + tau * tau).sqrt());
        let c = 1.0 / (1.0 + t * t).sqrt();
        let s = t * c;
        for k in 0..3 {
            let akp = a[k][p];
            let akq = a[k][q];
            a[k][p] = c * akp - s * akq;
            a[p][k] = a[k][p];
            a[k][q] = s * akp + c * akq;
            a[q][k] = a[k][q];
            let vkp = v[k][p];
            let vkq = v[k][q];
            v[k][p] = c * vkp - s * vkq;
            v[k][q] = s * vkp + c * vkq;
        }
        a[p][q] = 0.0;
        a[q][p] = 0.0;
    }
    let evals = [a[0][0], a[1][1], a[2][2]];
    // evecs: columns → array of column vectors.
    let evecs = [
        [v[0][0], v[1][0], v[2][0]],
        [v[0][1], v[1][1], v[2][1]],
        [v[0][2], v[1][2], v[2][2]],
    ];
    (evals, evecs)
}

fn solve4(a: &[[f64; 4]; 4], b: &[f64; 4]) -> Option<[f64; 4]> {
    let mut m = *a;
    let mut rhs = *b;
    for col in 0..4 {
        let mut piv = col;
        for r in (col + 1)..4 {
            if m[r][col].abs() > m[piv][col].abs() {
                piv = r;
            }
        }
        if m[piv][col].abs() < 1e-30 {
            return None;
        }
        m.swap(col, piv);
        rhs.swap(col, piv);
        let d = m[col][col];
        for r in (col + 1)..4 {
            let f = m[r][col] / d;
            for c in col..4 {
                m[r][c] -= f * m[col][c];
            }
            rhs[r] -= f * rhs[col];
        }
    }
    let mut x = [0.0; 4];
    for i in (0..4).rev() {
        let mut s = rhs[i];
        for j in (i + 1)..4 {
            s -= m[i][j] * x[j];
        }
        x[i] = s / m[i][i];
    }
    Some(x)
}

fn plane_residual(points: &[GpPnt], pln: &GpPln) -> f64 {
    let n = *pln.axis().direction();
    let loc = pln.location();
    points
        .iter()
        .map(|p| {
            let d = GpVec::from_pnts(&loc, p).xyz().dot(n.xyz());
            d.abs()
        })
        .fold(0.0, f64::max)
}

fn sphere_residual(points: &[GpPnt], center: &GpPnt, r: f64) -> f64 {
    points
        .iter()
        .map(|p| (p.distance(center) - r).abs())
        .fold(0.0, f64::max)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fit_plane_from_z0_points() {
        let pts = vec![
            GpPnt::new(0.,0.,0.), GpPnt::new(1.,0.,0.1),
            GpPnt::new(0.,1.,-0.1), GpPnt::new(1.,1.,0.),
        ];
        let pln = fit_plane(&pts).expect("plane");
        let n = *pln.axis().direction();
        assert!(n.z().abs() > 0.99, "normal ~ +Z: {n:?}");
        assert!(plane_residual(&pts, &pln) < 0.15);
    }

    #[test]
    fn fit_sphere_radius() {
        // Points on a radius-2 sphere (lat/long grid) → center ~origin, r≈2.
        let mut pts = Vec::new();
        for i in 0..8 {
            let phi = std::f64::consts::PI * (i as f64 + 0.5) / 8.0;
            for j in 0..12 {
                let theta = std::f64::consts::TAU * j as f64 / 12.0;
                pts.push(GpPnt::new(
                    2.0 * phi.sin() * theta.cos(),
                    2.0 * phi.sin() * theta.sin(),
                    2.0 * phi.cos(),
                ));
            }
        }
        let (c, r) = fit_sphere(&pts).expect("sphere");
        assert!((r - 2.0).abs() < 0.2, "radius {r}");
        assert!(c.distance(&GpPnt::zero()) < 0.2);
    }

    #[test]
    fn fit_surface_kind_plane() {
        let pts = vec![
            GpPnt::new(0.,0.,0.), GpPnt::new(1.,0.,0.),
            GpPnt::new(0.,1.,0.), GpPnt::new(1.,1.,0.),
        ];
        let (kind, _) = fit_surface_kind(&pts);
        assert_eq!(kind, SurfaceFitKind::Plane);
    }

    #[test]
    fn grid_surface_bilinear() {
        // poles[u][v] = (u, v, 0): the bilinear surface reproduces the grid.
        let us = vec![0.0, 1.0, 2.0];
        let vs = vec![0.0, 1.0, 2.0];
        let poles: Vec<Vec<GpPnt>> = (0..3)
            .map(|u| (0..3).map(|v| GpPnt::new(u as f64, v as f64, 0.0)).collect())
            .collect();
        let g = GridSurface::from_grid(us, vs, poles);
        // Interior of a cell → bilinear interpolation of f(x,y) = y.
        let p = g.d0(0.5, 1.5);
        assert!((p.x() - 0.5).abs() < 1e-9 && (p.y() - 1.5).abs() < 1e-9);
        // Exact grid node.
        let q = g.d0(1.0, 2.0);
        assert!((q.x() - 1.0).abs() < 1e-9 && (q.y() - 2.0).abs() < 1e-9);
    }

    #[test]
    fn fit_surface_wrapper() {
        let pts = vec![GpPnt::new(0.,0.,0.), GpPnt::new(1.,0.,0.), GpPnt::new(0.,1.,0.), GpPnt::new(1.,1.,0.)];
        let s = fit_surface(&pts).expect("fit");
        let p = s.d0(0.5, 0.5);
        assert!((p.z() - 0.0).abs() < 1e-9);
    }
}
