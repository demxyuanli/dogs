//! Extended `GpTrsf` operations — decomposition into translation/rotation/
//! scale, interpolation, vector rotation, affine inverse, frame construction.
//! Source: OCCT `gp_Trsf` and `math_Recipes` (trsf decomposition).

use crate::gp::dir::GpDir;
use crate::gp::mat::GpMat;
use crate::gp::pnt::GpPnt;
use crate::gp::quaternion::GpQuaternion;
use crate::gp::quaternion_ext::slerp;
use crate::gp::trsf::GpTrsf;
use crate::gp::trsf_form::TrsfForm;
use crate::gp::vec::GpVec;

/// The translation part of `t` as a vector.
pub fn translation_of(t: &GpTrsf) -> GpVec {
    GpVec::from_xyz(&t.loc)
}

/// The linear part of `t` (rotation + uniform scale) as a 3×3 matrix.
pub fn rotation_matrix_of(t: &GpTrsf) -> GpMat {
    t.matrix.multiply_scalar(t.scale)
}

/// The uniform scale factor stored in `t`.
pub fn scale_factor(t: &GpTrsf) -> f64 {
    t.scale
}

/// Decompose `t` into (translation, rotation matrix, scale).
pub fn decompose(t: &GpTrsf) -> (GpVec, GpMat, f64) {
    (GpVec::from_xyz(&t.loc), t.matrix.clone(), t.scale)
}

/// Interpolate between two transforms: lerp translation and scale, slerp the
/// rotation (converted to quaternions). `t` is clamped to [0, 1].
pub fn interpolate_transforms(a: &GpTrsf, b: &GpTrsf, t: f64) -> GpTrsf {
    let t = t.clamp(0.0, 1.0);
    let loc = a.loc.multiply(1.0 - t).add(&b.loc.multiply(t));
    let scale = a.scale * (1.0 - t) + b.scale * t;

    let mut qa = GpQuaternion::identity();
    let mut qb = GpQuaternion::identity();
    qa.set_matrix(&a.matrix);
    qb.set_matrix(&b.matrix);
    let q = slerp(&qa, &qb, t);

    GpTrsf {
        scale,
        shape: TrsfForm::CompoundTrsf,
        matrix: q.get_matrix(),
        loc,
    }
}

/// Apply only the rotation part (including uniform scale) of `t` to `v`,
/// ignoring the translation.
pub fn rotate_vector(t: &GpTrsf, v: &GpVec) -> GpVec {
    GpVec::from_xyz(&rotation_matrix_of(t).multiplied(v.xyz()))
}

/// Affine inverse of `t`: inverse rotation (transpose) over the scale, and
/// translation `-(R⁻¹ · trans) / scale`.
pub fn inverse_affine(t: &GpTrsf) -> GpTrsf {
    let inv_scale = 1.0 / t.scale;
    let inv_rot = t.matrix.transpose();
    let inv_loc = inv_rot.multiplied(&t.loc.multiply(-1.0)).multiply(inv_scale);
    GpTrsf {
        scale: inv_scale,
        shape: TrsfForm::CompoundTrsf,
        matrix: inv_rot,
        loc: inv_loc,
    }
}

/// Build the transform that maps the standard frame onto the frame at
/// `origin` with the given X and Y directions (Z = X × Y).
pub fn trsf_from_basis(origin: &GpPnt, x_dir: &GpDir, y_dir: &GpDir) -> GpTrsf {
    let z = x_dir
        .crossed(y_dir)
        .unwrap_or_else(|_| GpDir::new(0.0, 0.0, 1.0).unwrap());
    GpTrsf {
        scale: 1.0,
        shape: TrsfForm::Rotation,
        matrix: GpMat::from_cols(x_dir.xyz(), y_dir.xyz(), z.xyz()),
        loc: origin.coord,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gp::ax1::GpAx1;
    use crate::gp::dir::DirAxis;
    use std::f64::consts::FRAC_PI_2;

    fn assert_trsf_close(a: &GpTrsf, b: &GpTrsf, tol: f64) {
        assert!((a.scale - b.scale).abs() < tol, "scale {} vs {}", a.scale, b.scale);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (a.matrix.m[i][j] - b.matrix.m[i][j]).abs() < tol,
                    "matrix[{i}][{j}] {} vs {}",
                    a.matrix.m[i][j],
                    b.matrix.m[i][j]
                );
            }
        }
        assert!((a.loc.x - b.loc.x).abs() < tol, "loc.x");
        assert!((a.loc.y - b.loc.y).abs() < tol, "loc.y");
        assert!((a.loc.z - b.loc.z).abs() < tol, "loc.z");
    }

    #[test]
    fn translation_of_translated() {
        let mut t = GpTrsf::identity();
        t.set_translation_vec(&GpVec::new(3.0, -2.0, 7.0));
        let tr = translation_of(&t);
        assert!((tr.x() - 3.0).abs() < 1e-12);
        assert!((tr.y() + 2.0).abs() < 1e-12);
        assert!((tr.z() - 7.0).abs() < 1e-12);
    }

    #[test]
    fn scale_factor_of_scaled() {
        let mut t = GpTrsf::identity();
        t.set_scale(&GpPnt::new(1.0, 2.0, 3.0), 2.0).unwrap();
        assert!((scale_factor(&t) - 2.0).abs() < 1e-12);
        let (_, _, s) = decompose(&t);
        assert!((s - 2.0).abs() < 1e-12);
        // The rotation matrix (with scale) has determinant scale^3 = 8.
        assert!((rotation_matrix_of(&t).determinant() - 8.0).abs() < 1e-9);
    }

    #[test]
    fn interpolate_translation_halfway() {
        let a = GpTrsf::identity();
        let mut b = GpTrsf::identity();
        b.set_translation_vec(&GpVec::new(10.0, 0.0, 0.0));
        let mid = interpolate_transforms(&a, &b, 0.5);
        let tr = translation_of(&mid);
        assert!((tr.x() - 5.0).abs() < 1e-12, "x = {}", tr.x());
        assert!((tr.y() - 0.0).abs() < 1e-12);
        assert!((tr.z() - 0.0).abs() < 1e-12);
    }

    #[test]
    fn rotate_vector_quarter_turn_z() {
        let mut t = GpTrsf::identity();
        t.set_rotation_ax1(&GpAx1::from_axis(DirAxis::Z), FRAC_PI_2).unwrap();
        let v = rotate_vector(&t, &GpVec::new(1.0, 0.0, 0.0));
        assert!((v.x() - 0.0).abs() < 1e-12, "x = {}", v.x());
        assert!((v.y() - 1.0).abs() < 1e-12, "y = {}", v.y());
        assert!((v.z() - 0.0).abs() < 1e-12);
    }

    #[test]
    fn inverse_affine_composes_to_identity() {
        let mut t = GpTrsf::identity();
        t.set_rotation_ax1(&GpAx1::from_axis(DirAxis::Z), 0.7).unwrap();
        let mut tr = GpTrsf::identity();
        tr.set_translation_vec(&GpVec::new(4.0, 5.0, 6.0));
        t.multiply(&tr); // t = R * T
        let inv = inverse_affine(&t);
        let prod = &t * &inv;
        assert_trsf_close(&prod, &GpTrsf::identity(), 1e-9);
    }

    #[test]
    fn trsf_from_basis_roundtrip() {
        let origin = GpPnt::new(1.0, 2.0, 3.0);
        let x = GpDir::new(0.0, 1.0, 0.0).unwrap();
        let y = GpDir::new(0.0, 0.0, 1.0).unwrap();
        let t = trsf_from_basis(&origin, &x, &y);

        let p0 = GpPnt::new(0.0, 0.0, 0.0).transformed(&t);
        assert!((p0.x() - 1.0).abs() < 1e-12);
        assert!((p0.y() - 2.0).abs() < 1e-12);
        assert!((p0.z() - 3.0).abs() < 1e-12);

        let px = GpPnt::new(1.0, 0.0, 0.0).transformed(&t);
        assert!((px.x() - 1.0).abs() < 1e-12);
        assert!((px.y() - 3.0).abs() < 1e-12);
        assert!((px.z() - 3.0).abs() < 1e-12);

        let inv = t.inverted().unwrap();
        let back = px.transformed(&inv);
        assert!((back.x() - 1.0).abs() < 1e-12);
        assert!((back.y() - 0.0).abs() < 1e-12);
        assert!((back.z() - 0.0).abs() < 1e-12);
    }
}
