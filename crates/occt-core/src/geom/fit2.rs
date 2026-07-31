//! Least-squares geometric fitting (lines, circles, planes).
//! Source: `gce_MakeLin`/`gce_MakeCirc`, `math_LeastSquare`.
use crate::gp::{GpAx1, GpAx3, GpDir, GpPln, GpPnt, GpPnt2d, GpVec, GpVec2d, GpXyz};

/// Fit a least-squares line through ≥2 points.
/// Returns (centroid, unit direction along the principal axis).
pub fn fit_line(points: &[GpPnt]) -> Option<(GpPnt, GpVec)> {
    if points.len() < 2 { return None; }
    let c = centroid(points);
    let cov = covariance3(points, &c);
    let (eig, vecs) = eigen_3x3(cov);
    let mut imax = 0;
    for i in 1..3 { if eig[i] > eig[imax] { imax = i; } }
    if eig[imax] <= crate::precision::RESOLUTION { return None; }
    let d = GpVec::new(vecs[0][imax], vecs[1][imax], vecs[2][imax]);
    Some((c, d.normalized()))
}

/// Kasa (algebraic) circle fit: center + radius. Points must be near-coplanar.
/// Returns None if degenerate (collinear/coincident points, non-finite radius).
pub fn fit_circle(points: &[GpPnt]) -> Option<(GpPnt, f64)> {
    let n = points.len();
    if n < 3 { return None; }
    // Fit a plane, then fit the circle in the plane's local 2D frame.
    let pln = fit_plane(points)?;
    let z = pln.pos.direction();
    let x = perpendicular(&z);
    let y = z.crossed(&x).ok()?;
    let c3 = centroid(points);
    let mut pts2 = Vec::with_capacity(n);
    for p in points {
        let v = GpXyz::new(p.x() - c3.x(), p.y() - c3.y(), p.z() - c3.z());
        pts2.push(GpPnt2d::new(v.dot(x.xyz()), v.dot(y.xyz())));
    }
    let (c2, r) = fit_circle2d(&pts2)?;
    if !r.is_finite() || r <= 0.0 { return None; }
    // Map the 2D center back to 3D.
    let center = GpPnt::new(
        c3.x() + c2.x() * x.x() + c2.y() * y.x(),
        c3.y() + c2.x() * x.y() + c2.y() * y.y(),
        c3.z() + c2.x() * x.z() + c2.y() * y.z(),
    );
    Some((center, r))
}

/// Least-squares plane: normal = smallest-eigenvalue eigenvector of the
/// covariance matrix, offset through the centroid.
pub fn fit_plane(points: &[GpPnt]) -> Option<GpPln> {
    if points.len() < 3 { return None; }
    let c = centroid(points);
    let cov = covariance3(points, &c);
    let (eig, vecs) = eigen_3x3(cov);
    let mut imax = 0;
    let mut imin = 0;
    for i in 1..3 {
        if eig[i] > eig[imax] { imax = i; }
        if eig[i] < eig[imin] { imin = i; }
    }
    if eig[imax] <= crate::precision::RESOLUTION { return None; } // all coincident
    let normal = GpDir::from_xyz(&GpXyz::new(vecs[0][imin], vecs[1][imin], vecs[2][imin])).ok()?;
    Some(GpPln::new(GpAx3::from_ax1(&GpAx1::new(c, normal))))
}

/// Fit a least-squares line through ≥2 2D points.
/// Returns (centroid, unit direction along the principal axis).
pub fn fit_line2d(points: &[GpPnt2d]) -> Option<(GpPnt2d, GpVec2d)> {
    let n = points.len();
    if n < 2 { return None; }
    let mut cx = 0.0; let mut cy = 0.0;
    for p in points { cx += p.x(); cy += p.y(); }
    cx /= n as f64; cy /= n as f64;
    let mut sxx = 0.0; let mut sxy = 0.0; let mut syy = 0.0;
    for p in points {
        let dx = p.x() - cx; let dy = p.y() - cy;
        sxx += dx * dx; sxy += dx * dy; syy += dy * dy;
    }
    if sxx + syy <= crate::precision::RESOLUTION { return None; } // all coincident
    let theta = 0.5 * (2.0 * sxy).atan2(sxx - syy);
    Some((GpPnt2d::new(cx, cy), GpVec2d::new(theta.cos(), theta.sin())))
}

/// Kasa (algebraic) 2D circle fit: center + radius.
/// Returns None if degenerate (collinear/coincident points, non-finite radius).
pub fn fit_circle2d(points: &[GpPnt2d]) -> Option<(GpPnt2d, f64)> {
    let n = points.len();
    if n < 3 { return None; }
    // Solve the linearized problem: x²+y² = 2cx·x + 2cy·y + (r² - cx² - cy²).
    let mut sx = 0.0; let mut sy = 0.0;
    let mut sxx = 0.0; let mut sxy = 0.0; let mut syy = 0.0;
    let mut sxu = 0.0; let mut syu = 0.0; let mut su = 0.0;
    for p in points {
        let x = p.x(); let y = p.y();
        let u = x * x + y * y;
        sx += x; sy += y;
        sxx += x * x; sxy += x * y; syy += y * y;
        sxu += x * u; syu += y * u; su += u;
    }
    let m = [[sxx, sxy, sx], [sxy, syy, sy], [sx, sy, n as f64]];
    let b = [sxu, syu, su];
    let sol = solve3x3(&m, &b)?;
    let cx = sol[0] * 0.5;
    let cy = sol[1] * 0.5;
    let r2 = sol[2] + cx * cx + cy * cy;
    if !r2.is_finite() || r2 <= 0.0 { return None; }
    Some((GpPnt2d::new(cx, cy), r2.sqrt()))
}

/// Σ (distance² − r²)² — algebraic residual of a circle fit, for diagnostics.
pub fn residual_square(points: &[GpPnt], p: &GpPnt, r: f64) -> f64 {
    let mut sum = 0.0;
    for q in points {
        let e = q.square_distance(p) - r * r;
        sum += e * e;
    }
    sum
}

// --- internals ---

fn centroid(points: &[GpPnt]) -> GpPnt {
    let n = points.len() as f64;
    let mut cx = 0.0; let mut cy = 0.0; let mut cz = 0.0;
    for p in points { cx += p.x(); cy += p.y(); cz += p.z(); }
    GpPnt::new(cx / n, cy / n, cz / n)
}

fn covariance3(points: &[GpPnt], c: &GpPnt) -> [[f64; 3]; 3] {
    let mut cov = [[0.0f64; 3]; 3];
    for p in points {
        let dx = p.x() - c.x(); let dy = p.y() - c.y(); let dz = p.z() - c.z();
        cov[0][0] += dx * dx; cov[0][1] += dx * dy; cov[0][2] += dx * dz;
        cov[1][0] += dy * dx; cov[1][1] += dy * dy; cov[1][2] += dy * dz;
        cov[2][0] += dz * dx; cov[2][1] += dz * dy; cov[2][2] += dz * dz;
    }
    let n = points.len() as f64;
    for i in 0..3 { for j in 0..3 { cov[i][j] /= n; } }
    cov
}

/// Unit vector perpendicular to `d` (uses the coordinate axis least aligned with d).
fn perpendicular(d: &GpDir) -> GpDir {
    let (ax, ay, az) = (d.x().abs(), d.y().abs(), d.z().abs());
    let v = if ax <= ay && ax <= az {
        GpXyz::new(0.0, d.z(), -d.y())
    } else if ay <= az {
        GpXyz::new(-d.z(), 0.0, d.x())
    } else {
        GpXyz::new(d.y(), -d.x(), 0.0)
    };
    GpDir::from_xyz(&v).unwrap_or_default()
}

/// Jacobi eigenvalue decomposition for a symmetric 3×3 matrix.
/// Returns (eigenvalues, eigenvectors as columns).
fn eigen_3x3(mut a: [[f64; 3]; 3]) -> ([f64; 3], [[f64; 3]; 3]) {
    let mut v = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let mut eig = [a[0][0], a[1][1], a[2][2]];
    for _ in 0..50 {
        let mut max_off = 0.0;
        let mut p = 0usize; let mut q = 1usize;
        for i in 0..3 {
            for j in (i + 1)..3 {
                if a[i][j].abs() > max_off { max_off = a[i][j].abs(); p = i; q = j; }
            }
        }
        if max_off < 1e-15 { break; }
        let theta = 0.5 * (eig[q] - eig[p]) / a[p][q];
        let t = 1.0 / (theta.abs() + (1.0 + theta * theta).sqrt());
        let t = if theta < 0.0 { -t } else { t };
        let c = 1.0 / (1.0 + t * t).sqrt();
        let s = t * c;
        let tau = s / (1.0 + c);
        let h = t * a[p][q];
        eig[p] -= h; eig[q] += h;
        a[p][q] = 0.0;
        for j in 0..3 {
            if j != p && j != q {
                let g = a[if j < p { j } else { p }][if j < p { p } else { j }];
                let hh = a[if j < q { j } else { q }][if j < q { q } else { j }];
                a[p][j] = g - s * (hh + g * tau);
                a[q][j] = hh + s * (g - hh * tau);
            }
        }
        for j in 0..3 {
            let g = v[j][p]; let h = v[j][q];
            v[j][p] = g - s * (h + g * tau);
            v[j][q] = h + s * (g - h * tau);
        }
    }
    (eig, v)
}

/// Solve a 3×3 linear system via Cramer's rule. None if singular.
fn solve3x3(a: &[[f64; 3]; 3], b: &[f64; 3]) -> Option<[f64; 3]> {
    let det = |m: &[[f64; 3]; 3]| {
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    };
    let d0 = det(a);
    if d0.abs() < 1e-300 { return None; }
    let mut m1 = *a; m1[0][0] = b[0]; m1[1][0] = b[1]; m1[2][0] = b[2];
    let mut m2 = *a; m2[0][1] = b[0]; m2[1][1] = b[1]; m2[2][1] = b[2];
    let mut m3 = *a; m3[0][2] = b[0]; m3[1][2] = b[1]; m3[2][2] = b[2];
    Some([det(&m1) / d0, det(&m2) / d0, det(&m3) / d0])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn line_collinear() {
        let pts = vec![
            GpPnt::new(0.0, 0.0, 0.0), GpPnt::new(1.0, 1.0, 1.0),
            GpPnt::new(2.0, 2.0, 2.0), GpPnt::new(3.0, 3.0, 3.0),
        ];
        let (c, d) = fit_line(&pts).unwrap();
        assert!(c.distance(&GpPnt::new(1.5, 1.5, 1.5)) < 1e-9);
        // direction parallel to (1,1,1): cross product ~ 0
        let cross = d.crossed(&GpVec::new(1.0, 1.0, 1.0));
        assert!(cross.magnitude() < 1e-6);
        assert!((d.magnitude() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn circle_unit() {
        let s = std::f64::consts::FRAC_1_SQRT_2;
        let pts = vec![
            GpPnt::new(1.0, 0.0, 0.0), GpPnt::new(0.0, 1.0, 0.0),
            GpPnt::new(-1.0, 0.0, 0.0), GpPnt::new(0.0, -1.0, 0.0),
            GpPnt::new(s, s, 0.0), GpPnt::new(-s, s, 0.0),
        ];
        let (c, r) = fit_circle(&pts).unwrap();
        assert!(c.distance(&GpPnt::new(0.0, 0.0, 0.0)) < 1e-6);
        assert!((r - 1.0).abs() < 1e-6);
    }

    #[test]
    fn plane_z0() {
        let pts = vec![
            GpPnt::new(0.0, 0.0, 0.0), GpPnt::new(1.0, 0.0, 0.0),
            GpPnt::new(0.0, 1.0, 0.0), GpPnt::new(1.0, 1.0, 0.0),
        ];
        let pln = fit_plane(&pts).unwrap();
        let n = pln.pos.direction();
        assert!(n.z().abs() > 0.999);
        assert!(n.x().abs() < 1e-6);
        assert!(n.y().abs() < 1e-6);
    }

    #[test]
    fn circle2d_unit() {
        let s = std::f64::consts::FRAC_1_SQRT_2;
        let pts = vec![
            GpPnt2d::new(1.0, 0.0), GpPnt2d::new(0.0, 1.0),
            GpPnt2d::new(-1.0, 0.0), GpPnt2d::new(0.0, -1.0),
            GpPnt2d::new(s, s),
        ];
        let (c, r) = fit_circle2d(&pts).unwrap();
        assert!(c.distance(&GpPnt2d::new(0.0, 0.0)) < 1e-9);
        assert!((r - 1.0).abs() < 1e-9);
    }

    #[test]
    fn line2d_axis() {
        let pts = vec![GpPnt2d::new(0.0, 0.0), GpPnt2d::new(1.0, 0.0), GpPnt2d::new(2.0, 0.0)];
        let (c, d) = fit_line2d(&pts).unwrap();
        assert!(c.distance(&GpPnt2d::new(1.0, 0.0)) < 1e-9);
        assert!(d.y().abs() < 1e-12);
        assert!(d.x().abs() > 0.999);
    }

    #[test]
    fn fit_degenerate() {
        let pts = vec![GpPnt::new(1.0, 1.0, 1.0); 4];
        assert!(fit_line(&pts).is_none());
        assert!(fit_plane(&pts).is_none());
        assert!(fit_circle(&pts).is_none());
    }
}
