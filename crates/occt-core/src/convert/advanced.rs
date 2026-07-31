//! Advanced coordinate conversions: spherical/cylindrical coordinates,
//! rotations about arbitrary axes, normal transformation, barycentric
//! coordinates and scalar interpolation helpers.

use crate::gp::{GpDir, GpPnt, GpTrsf, GpXyz};

/// Cartesian point to spherical coordinates `(r, theta, phi)` relative to
/// `origin`. `theta` is the polar angle from the +Z axis, `phi` the azimuth
/// measured from the +X axis in the XY plane.
pub fn cartesian_to_spherical_point(p: &GpPnt, origin: &GpPnt) -> (f64, f64, f64) {
    let dx = p.x() - origin.x();
    let dy = p.y() - origin.y();
    let dz = p.z() - origin.z();
    let r = (dx * dx + dy * dy + dz * dz).sqrt();
    if r < f64::EPSILON {
        return (0.0, 0.0, 0.0);
    }
    let theta = (dz / r).clamp(-1.0, 1.0).acos();
    let phi = dy.atan2(dx);
    (r, theta, phi)
}

/// Inverse of [`cartesian_to_spherical_point`].
pub fn spherical_to_cartesian_point(r: f64, theta: f64, phi: f64, origin: &GpPnt) -> GpPnt {
    let st = theta.sin();
    GpPnt::new(
        origin.x() + r * st * phi.cos(),
        origin.y() + r * st * phi.sin(),
        origin.z() + r * theta.cos(),
    )
}

/// Cartesian point to cylindrical coordinates `(rho, theta, axial)` about the
/// `axis` line through `origin`. `axial` is the projection onto the axis,
/// `rho` the perpendicular distance, `theta` the azimuth of the perpendicular
/// component measured from +X.
pub fn cartesian_to_cylindrical_point(p: &GpPnt, origin: &GpPnt, axis: &GpDir) -> (f64, f64, f64) {
    let dx = p.x() - origin.x();
    let dy = p.y() - origin.y();
    let dz = p.z() - origin.z();
    let (ax, ay, az) = (axis.x(), axis.y(), axis.z());
    let axial = dx * ax + dy * ay + dz * az;
    let rx = dx - axial * ax;
    let ry = dy - axial * ay;
    let rho = (rx * rx + ry * ry).sqrt();
    let theta = ry.atan2(rx);
    (rho, theta, axial)
}

/// Rotate `p` about the `axis_dir` line through `axis_origin` by `angle`
/// radians using Rodrigues' rotation formula.
pub fn rotate_point_about_axis(
    p: &GpPnt,
    axis_origin: &GpPnt,
    axis_dir: &GpDir,
    angle: f64,
) -> GpPnt {
    let vx = p.x() - axis_origin.x();
    let vy = p.y() - axis_origin.y();
    let vz = p.z() - axis_origin.z();
    let (kx, ky, kz) = (axis_dir.x(), axis_dir.y(), axis_dir.z());
    let c = angle.cos();
    let s = angle.sin();
    let dot = vx * kx + vy * ky + vz * kz;
    // k x v
    let (cx, cy, cz) = (
        ky * vz - kz * vy,
        kz * vx - kx * vz,
        kx * vy - ky * vx,
    );
    GpPnt::new(
        axis_origin.x() + vx * c + cx * s + kx * dot * (1.0 - c),
        axis_origin.y() + vy * c + cy * s + ky * dot * (1.0 - c),
        axis_origin.z() + vz * c + cz * s + kz * dot * (1.0 - c),
    )
}

/// Transform a normal vector by the inverse-transpose of `t`'s linear part,
/// so normals stay perpendicular under non-uniform scaling. Translation is
/// ignored.
pub fn transform_normal(t: &GpTrsf, normal: &GpXyz) -> GpXyz {
    let m = t.inverted().vectorial_part();
    let x = m.value(0, 0) * normal.x() + m.value(1, 0) * normal.y() + m.value(2, 0) * normal.z();
    let y = m.value(0, 1) * normal.x() + m.value(1, 1) * normal.y() + m.value(2, 1) * normal.z();
    let z = m.value(0, 2) * normal.x() + m.value(1, 2) * normal.y() + m.value(2, 2) * normal.z();
    GpXyz::new(x, y, z)
}

/// Barycentric coordinates `(u, v, w)` of `p` in triangle `abc`, computed in
/// the XY plane. Returns `(0, 0, 0)` for a degenerate triangle.
pub fn barycentric(a: &GpPnt, b: &GpPnt, c: &GpPnt, p: &GpPnt) -> (f64, f64, f64) {
    let v0x = b.x() - a.x();
    let v0y = b.y() - a.y();
    let v1x = c.x() - a.x();
    let v1y = c.y() - a.y();
    let v2x = p.x() - a.x();
    let v2y = p.y() - a.y();
    let d00 = v0x * v0x + v0y * v0y;
    let d01 = v0x * v1x + v0y * v1y;
    let d11 = v1x * v1x + v1y * v1y;
    let d20 = v2x * v0x + v2y * v0y;
    let d21 = v2x * v1x + v2y * v1y;
    let denom = d00 * d11 - d01 * d01;
    if denom.abs() < 1e-12 {
        return (0.0, 0.0, 0.0);
    }
    let v = (d11 * d20 - d01 * d21) / denom;
    let w = (d00 * d21 - d01 * d20) / denom;
    (1.0 - v - w, v, w)
}

/// Bilinear interpolation of the four corner values at `(u, v)` in `[0, 1]^2`.
pub fn bilinear_interp(p00: f64, p10: f64, p01: f64, p11: f64, u: f64, v: f64) -> f64 {
    p00 * (1.0 - u) * (1.0 - v) + p10 * u * (1.0 - v) + p01 * (1.0 - u) * v + p11 * u * v
}

/// Linear interpolation from `a` to `b` at parameter `t`.
pub fn lerp(a: f64, b: f64, t: f64) -> f64 {
    a + (b - a) * t
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gp::GpPnt;

    fn assert_approx(a: f64, b: f64, tol: f64) {
        assert!((a - b).abs() < tol, "{a} != {b}");
    }

    #[test]
    fn lerp_half() {
        assert_approx(lerp(0.0, 1.0, 0.5), 0.5, 1e-12);
    }

    #[test]
    fn barycentric_of_centroid() {
        let a = GpPnt::new(0.0, 0.0, 0.0);
        let b = GpPnt::new(1.0, 0.0, 0.0);
        let c = GpPnt::new(0.0, 1.0, 0.0);
        let p = GpPnt::new(1.0 / 3.0, 1.0 / 3.0, 0.0);
        let (u, v, w) = barycentric(&a, &b, &c, &p);
        assert_approx(u, 1.0 / 3.0, 1e-12);
        assert_approx(v, 1.0 / 3.0, 1e-12);
        assert_approx(w, 1.0 / 3.0, 1e-12);
    }

    #[test]
    fn bilinear_center() {
        assert_approx(bilinear_interp(0.0, 1.0, 1.0, 0.0, 0.5, 0.5), 0.5, 1e-12);
    }
}
