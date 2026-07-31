//! PCA-based OBB computation. Source: Bnd_OBB + Jacobi eigenvalue decomposition.
use crate::gp::{GpPnt, GpDir, GpXyz};
use crate::bnd::BndOBB;

/// Compute optimal OBB for point cloud using PCA (principal component analysis).
/// Returns OBB with center and axes aligned to principal directions.
pub fn compute_obb_pca(points: &[GpPnt]) -> BndOBB {
    let n = points.len();
    if n == 0 { return BndOBB::new(); }
    if n == 1 { return BndOBB::from_params(&points[0], &GpDir::from_axis(crate::gp::dir::DirAxis::X), &GpDir::from_axis(crate::gp::dir::DirAxis::Y), &GpDir::from_axis(crate::gp::dir::DirAxis::Z), 0.0, 0.0, 0.0); }

    // Compute centroid
    let mut cx = 0.0f64; let mut cy = 0.0f64; let mut cz = 0.0f64;
    for p in points { cx += p.x(); cy += p.y(); cz += p.z(); }
    cx /= n as f64; cy /= n as f64; cz /= n as f64;
    let centroid = GpPnt::new(cx, cy, cz);

    // Compute 3x3 covariance matrix
    let mut cov = [[0.0f64; 3]; 3];
    for p in points {
        let dx = p.x() - cx; let dy = p.y() - cy; let dz = p.z() - cz;
        cov[0][0] += dx*dx; cov[0][1] += dx*dy; cov[0][2] += dx*dz;
        cov[1][0] += dy*dx; cov[1][1] += dy*dy; cov[1][2] += dy*dz;
        cov[2][0] += dz*dx; cov[2][1] += dz*dy; cov[2][2] += dz*dz;
    }
    for i in 0..3 { for j in 0..3 { cov[i][j] /= n as f64; } }

    // Jacobi eigenvalue decomposition on 3x3 symmetric matrix
    let (vals, vecs) = jacobi_3x3(cov);

    // Principal axes = eigenvectors
    let xdir = GpDir::from_xyz(&GpXyz::new(vecs[0][0], vecs[1][0], vecs[2][0])).unwrap_or_default();
    let ydir = GpDir::from_xyz(&GpXyz::new(vecs[0][1], vecs[1][1], vecs[2][1])).unwrap_or_default();
    let zdir = GpDir::from_xyz(&GpXyz::new(vecs[0][2], vecs[1][2], vecs[2][2])).unwrap_or_default();

    // Compute half-sizes via projection onto each axis
    let axes = [xdir.xyz(), ydir.xyz(), zdir.xyz()];
    let mut hdims = [0.0f64; 3];
    for p in points {
        let d = GpXyz::new(p.x()-cx, p.y()-cy, p.z()-cz);
        for k in 0..3 { hdims[k] = f64::max(hdims[k], d.dot(axes[k]).abs()); }
    }

    BndOBB::from_params(&centroid, &xdir, &ydir, &zdir, hdims[0], hdims[1], hdims[2])
}

/// Jacobi eigenvalue decomposition for 3x3 symmetric matrix.
/// Returns (eigenvalues [3], eigenvectors [3][3] as column vectors).
fn jacobi_3x3(mut a: [[f64; 3]; 3]) -> ([f64; 3], [[f64; 3]; 3]) {
    let mut v = [[1.0,0.,0.],[0.,1.,0.],[0.,0.,1.]];
    let mut eig = [a[0][0], a[1][1], a[2][2]];
    let max_iter = 50;

    for _iter in 0..max_iter {
        // Find max off-diagonal
        let mut max_off: f64 = 0.0; let mut p = 0usize; let mut q = 1usize;
        for i in 0..3 { for j in (i+1)..3 { if a[i][j].abs() > max_off { max_off = a[i][j].abs(); p = i; q = j; }}}
        if max_off < 1e-15 { break; }

        let theta = 0.5 * (eig[q] - eig[p]) / a[p][q];
        let t = 1.0 / (theta.abs() + (1.0 + theta*theta).sqrt());
        if theta < 0.0 { let t = -t; }
        let c = 1.0/(1.0+t*t).sqrt();
        let s = t * c;
        let tau = s/(1.0+c);

        // Update eigenvalues
        let h = t * a[p][q];
        eig[p] -= h; eig[q] += h;
        a[p][q] = 0.0;

        // Rotate rows/cols
        for j in 0..3 {
            if j != p && j != q {
                let g = a[if j<p{j}else{p}][if j<p{p}else{j}];
                let h = a[if j<q{j}else{q}][if j<q{q}else{j}];
                a[p][j] = g - s*(h+g*tau);
                a[q][j] = h + s*(g-h*tau);
            }
        }
        // Rotate eigenvectors
        for j in 0..3 {
            let g = v[j][p]; let h = v[j][q];
            v[j][p] = g - s*(h+g*tau);
            v[j][q] = h + s*(g-h*tau);
        }
    }
    (eig, v)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn obb_pca_cube() {
        let pts = vec![
            GpPnt::new(0.,0.,0.),GpPnt::new(1.,0.,0.),GpPnt::new(0.,1.,0.),GpPnt::new(1.,1.,0.),
            GpPnt::new(0.,0.,1.),GpPnt::new(1.,0.,1.),GpPnt::new(0.,1.,1.),GpPnt::new(1.,1.,1.),
        ];
        let obb = compute_obb_pca(&pts);
        assert!(!obb.is_void());
        let (hx, hy, hz) = obb.half_sizes();
        assert!(hx > 0.4 && hx < 0.7);
        assert!(hy > 0.4 && hy < 0.7);
        assert!(hz > 0.4 && hz < 0.7);
    }
}
