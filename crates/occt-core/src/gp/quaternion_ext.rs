//! Extended quaternion operations — slerp, axis/angle conversion, Euler
//! angles, point rotation.
//! Source: `gp_Quaternion.hxx` (SetVectorAndAngle, slerp, Euler).

use crate::gp::dir::GpDir;
use crate::gp::pnt::GpPnt;
use crate::gp::quaternion::GpQuaternion;
use crate::gp::vec::GpVec;

/// Conjugate quaternion (negate the vector part).
pub fn conjugate(q: &GpQuaternion) -> GpQuaternion {
    GpQuaternion::new(-q.x, -q.y, -q.z, q.w)
}

/// Inverse: conjugate / |q|² (unit quaternions: just the conjugate).
pub fn inverse(q: &GpQuaternion) -> GpQuaternion {
    let n2 = q.x * q.x + q.y * q.y + q.z * q.z + q.w * q.w;
    if n2 < 1e-30 {
        return GpQuaternion::identity();
    }
    GpQuaternion::new(-q.x / n2, -q.y / n2, -q.z / n2, q.w / n2)
}

/// Hamilton product `a * b` (apply `b` first, then `a`, when used to rotate).
pub fn multiply(a: &GpQuaternion, b: &GpQuaternion) -> GpQuaternion {
    GpQuaternion::new(
        a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
        a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x,
        a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w,
        a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z,
    )
}

/// Dot product of two quaternions (cos of half the relative angle when both
/// are unit).
pub fn dot(a: &GpQuaternion, b: &GpQuaternion) -> f64 {
    a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w
}

/// Quaternion from an axis and angle (radians).
pub fn from_axis_angle(axis: &GpDir, angle: f64) -> GpQuaternion {
    let s = (angle * 0.5).sin();
    GpQuaternion::new(axis.x() * s, axis.y() * s, axis.z() * s, (angle * 0.5).cos())
}

/// Quaternion from a rotation axis as a vector (normalized internally).
pub fn from_axis_angle_vec(axis: &GpVec, angle: f64) -> GpQuaternion {
    match GpDir::from_vec(axis) {
        Ok(d) => from_axis_angle(&d, angle),
        Err(_) => GpQuaternion::identity(),
    }
}

/// Extract (axis, angle) from a quaternion; angle in [0, π].
pub fn to_axis_angle(q: &GpQuaternion) -> (GpDir, f64) {
    let n = (q.x * q.x + q.y * q.y + q.z * q.z).sqrt();
    if n < 1e-20 {
        return (GpDir::new(1.0, 0.0, 0.0).unwrap(), 0.0);
    }
    let angle = 2.0 * q.w.clamp(-1.0, 1.0).acos();
    let d = GpDir::new(q.x / n, q.y / n, q.z / n).unwrap_or_else(|_| GpDir::new(1.0, 0.0, 0.0).unwrap());
    (d, angle)
}

/// Rotate a point by the quaternion: p' = q · p · q⁻¹.
pub fn rotate_point(q: &GpQuaternion, p: &GpPnt) -> GpPnt {
    let v = GpVec::new(p.x(), p.y(), p.z());
    let qv = GpQuaternion::new(v.x(), v.y(), v.z(), 0.0);
    let r = multiply(&multiply(q, &qv), &conjugate(q));
    GpPnt::new(r.x, r.y, r.z)
}

/// Spherical linear interpolation between two unit quaternions (`t` in [0,1]).
/// Takes the short way (negates `b` if the dot product is negative).
pub fn slerp(a: &GpQuaternion, b: &GpQuaternion, t: f64) -> GpQuaternion {
    let mut b = *b;
    let mut cos_theta = dot(a, &b);
    if cos_theta < 0.0 {
        b.x = -b.x;
        b.y = -b.y;
        b.z = -b.z;
        b.w = -b.w;
        cos_theta = -cos_theta;
    }
    let t = t.clamp(0.0, 1.0);
    let (q0, q1) = (*a, b);
    if cos_theta > 0.9995 {
        // Nearly parallel: linear interpolation + normalize.
        let x = q0.x + t * (q1.x - q0.x);
        let y = q0.y + t * (q1.y - q0.y);
        let z = q0.z + t * (q1.z - q0.z);
        let w = q0.w + t * (q1.w - q0.w);
        let n = (x * x + y * y + z * z + w * w).sqrt();
        return if n > 1e-30 { GpQuaternion::new(x / n, y / n, z / n, w / n) } else { q0 };
    }
    let theta = cos_theta.clamp(-1.0, 1.0).acos();
    let sin_theta = theta.sin();
    let w0 = ((1.0 - t) * theta).sin() / sin_theta;
    let w1 = (t * theta).sin() / sin_theta;
    GpQuaternion::new(
        w0 * q0.x + w1 * q1.x,
        w0 * q0.y + w1 * q1.y,
        w0 * q0.z + w1 * q1.z,
        w0 * q0.w + w1 * q1.w,
    )
}

/// Quaternion from intrinsic Z-Y-X (yaw, pitch, roll) Euler angles (radians).
pub fn from_euler_zyx(yaw: f64, pitch: f64, roll: f64) -> GpQuaternion {
    // sin_cos() returns (sin, cos); bind accordingly.
    let (sy, cy) = (yaw * 0.5).sin_cos();
    let (sp, cp) = (pitch * 0.5).sin_cos();
    let (sr, cr) = (roll * 0.5).sin_cos();
    GpQuaternion::new(
        sr * cp * cy - cr * sp * sy,
        cr * sp * cy + sr * cp * sy,
        cr * cp * sy - sr * sp * cy,
        cr * cp * cy + sr * sp * sy,
    )
}

/// Extract intrinsic Z-Y-X Euler angles from a quaternion (radians), matching
/// the composition in `from_euler_zyx` (yaw about Z, pitch about Y, roll
/// about X).
pub fn to_euler_zyx(q: &GpQuaternion) -> (f64, f64, f64) {
    let (x, y, z, w) = (q.x, q.y, q.z, q.w);
    let sin_pitch = 2.0 * (w * y - z * x);
    let pitch = if sin_pitch.abs() >= 1.0 {
        std::f64::consts::FRAC_PI_2 * sin_pitch.signum()
    } else {
        sin_pitch.asin()
    };
    let yaw = (2.0 * (w * z + x * y)).atan2(1.0 - 2.0 * (y * y + z * z));
    let roll = (2.0 * (w * x + y * z)).atan2(1.0 - 2.0 * (x * x + y * y));
    (yaw, pitch, roll)
}

/// Angle (radians) between two rotation quaternions, in [0, π].
pub fn angle_between(a: &GpQuaternion, b: &GpQuaternion) -> f64 {
    2.0 * dot(a, b).clamp(-1.0, 1.0).abs().acos()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::precision::ANGULAR;

    fn q_unit() -> GpQuaternion {
        GpQuaternion::identity()
    }

    #[test]
    fn axis_angle_roundtrip() {
        let q = from_axis_angle(&GpDir::new(0.0, 0.0, 1.0).unwrap(), std::f64::consts::FRAC_PI_2);
        let (d, a) = to_axis_angle(&q);
        assert!((a - std::f64::consts::FRAC_PI_2).abs() < 1e-12);
        assert!(d.z() > 0.999);
    }

    #[test]
    fn rotate_point_quarter_turn_z() {
        // 90° about Z maps (1,0,0) → (0,1,0).
        let q = from_axis_angle(&GpDir::new(0.0, 0.0, 1.0).unwrap(), std::f64::consts::FRAC_PI_2);
        let p = rotate_point(&q, &GpPnt::new(1.0, 0.0, 0.0));
        assert!((p.x() - 0.0).abs() < 1e-12);
        assert!((p.y() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn slerp_halfway_is_midpoint() {
        let id = q_unit();
        let q90 = from_axis_angle(&GpDir::new(0.0, 0.0, 1.0).unwrap(), std::f64::consts::FRAC_PI_2);
        let mid = slerp(&id, &q90, 0.5);
        let (_, a) = to_axis_angle(&mid);
        assert!((a - std::f64::consts::FRAC_PI_4).abs() < 1e-9);
        // Endpoints.
        let e0 = slerp(&id, &q90, 0.0);
        let e1 = slerp(&id, &q90, 1.0);
        assert!(angle_between(&e0, &id) < 1e-9);
        assert!(angle_between(&e1, &q90) < 1e-9);
    }

    #[test]
    fn euler_roundtrip() {
        let q = from_euler_zyx(0.3, 0.4, 0.5);
        let (y, p, r) = to_euler_zyx(&q);
        assert!((y - 0.3).abs() < 1e-9, "yaw {y}");
        assert!((p - 0.4).abs() < 1e-9, "pitch {p}");
        assert!((r - 0.5).abs() < 1e-9, "roll {r}");
    }

    #[test]
    fn inverse_composes_to_identity() {
        let q = from_axis_angle(&GpDir::new(0.0, 1.0, 0.0).unwrap(), 0.7);
        let qinv = inverse(&q);
        let prod = multiply(&q, &qinv);
        assert!(angle_between(&prod, &q_unit()) < 1e-9);
    }

    #[test]
    fn angle_between_zero_for_same() {
        let a = from_axis_angle(&GpDir::new(1.0, 0.0, 0.0).unwrap(), 1.0);
        assert!(angle_between(&a, &a) < 1e-12);
        let _ = ANGULAR;
    }
}
