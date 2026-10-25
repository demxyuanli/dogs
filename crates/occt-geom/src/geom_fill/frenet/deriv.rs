//! The free functions of `GeomFill_Frenet.cxx` and
//! `GeomFill_ConstantBiNormal.cxx`: `FDeriv`, `DDeriv`, `CosAngle` and
//! `RotateTrihedron`.
//!
//! Source: `GeomFill_Frenet.cxx:42-97` and `cxx:397-470` (OCCT 8.0.0);
//! `FDeriv` / `DDeriv` are duplicated verbatim in
//! `GeomFill_ConstantBiNormal.cxx:33-54`.

use occt_core::gp::GpVec;
use occt_core::precision::{ANGULAR, RESOLUTION};

/// `FDeriv` (`GeomFill_Frenet.cxx:45-52`): computes `(F/|F|)'`.
pub(crate) fn fderiv(f: &GpVec, df: &GpVec) -> GpVec {
    let norma = f.magnitude();
    df.subtracted(&f.multiplied_scalar(f.dot(df) / (norma * norma)))
        .divided(norma)
}

/// `DDeriv` (`GeomFill_Frenet.cxx:57-67`): computes `(F/|F|)''`.
pub(crate) fn dderiv(f: &GpVec, df: &GpVec, d2f: &GpVec) -> GpVec {
    let norma = f.magnitude();
    let fd = f.dot(df);
    d2f.subtracted(&df.multiplied_scalar(2.0 * fd / (norma * norma)))
        .divided(norma)
        .subtracted(&f.multiplied_scalar(
            (df.square_magnitude() + f.dot(d2f) - 3.0 * fd * fd / (norma * norma))
                / (norma * norma * norma),
        ))
}

/// `CosAngle` (`GeomFill_Frenet.cxx:73-96`): cosine between two vectors, `1`
/// when either is degenerate, clamped to `[-1, 1]`.
pub(crate) fn cos_angle(v1: &GpVec, v2: &GpVec) -> f64 {
    let tol = RESOLUTION;
    let m1 = v1.magnitude();
    let m2 = v2.magnitude();
    if m1 <= tol || m2 <= tol {
        return 1.0;
    }
    let mut ca = v1.dot(v2) / (m1 * m2);
    if ca > 1.0 {
        ca = 1.0;
    }
    if ca < -1.0 {
        ca = -1.0;
    }
    ca
}

/// `GeomFill_Frenet::RotateTrihedron` (`cxx:402-470`): revolves the trihedron
/// so that `Tangent` coincides with `NewTangent`.
///
/// The two candidate rotation senses (`aV*1` / `aV*2`) are applied to the
/// tangent and the one closer to `NewTangent` wins, exactly as in OCCT.
pub(crate) fn rotate_trihedron(
    tangent: &mut GpVec,
    normal: &mut GpVec,
    bi_normal: &mut GpVec,
    new_tangent: &GpVec,
) -> bool {
    let an_inf_cos = ANGULAR.cos();
    let tol = RESOLUTION;

    let mut axis = tangent.crossed(new_tangent);
    let nt = axis.magnitude();
    if nt <= tol {
        return true;
    }
    axis = axis.divided(nt);

    let (px, py, pz) = (axis.x(), axis.y(), axis.z());
    let ca = cos_angle(tangent, new_tangent);
    let add_ca = 1.0 - ca;
    let sa = (1.0 - ca * ca).sqrt();

    let v11 = GpVec::new(
        add_ca * px * px + ca,
        add_ca * px * py - pz * sa,
        add_ca * px * pz + py * sa,
    );
    let v12 = GpVec::new(
        add_ca * px * px + ca,
        add_ca * px * py + pz * sa,
        add_ca * px * pz - py * sa,
    );
    let v21 = GpVec::new(
        add_ca * px * py + pz * sa,
        add_ca * py * py + ca,
        add_ca * py * pz - px * sa,
    );
    let v22 = GpVec::new(
        add_ca * px * py - pz * sa,
        add_ca * py * py + ca,
        add_ca * py * pz + px * sa,
    );
    let v31 = GpVec::new(
        add_ca * px * pz - py * sa,
        add_ca * py * pz + px * sa,
        add_ca * pz * pz + ca,
    );
    let v32 = GpVec::new(
        add_ca * px * pz + py * sa,
        add_ca * py * pz - px * sa,
        add_ca * pz * pz + ca,
    );

    let t1 = GpVec::new(tangent.dot(&v11), tangent.dot(&v21), tangent.dot(&v31));
    let t2 = GpVec::new(tangent.dot(&v12), tangent.dot(&v22), tangent.dot(&v32));

    if cos_angle(&t1, new_tangent) >= cos_angle(&t2, new_tangent) {
        *tangent = t1;
        *normal = GpVec::new(normal.dot(&v11), normal.dot(&v21), normal.dot(&v31));
        *bi_normal = GpVec::new(
            bi_normal.dot(&v11),
            bi_normal.dot(&v21),
            bi_normal.dot(&v31),
        );
    } else {
        *tangent = t2;
        *normal = GpVec::new(normal.dot(&v12), normal.dot(&v22), normal.dot(&v32));
        *bi_normal = GpVec::new(
            bi_normal.dot(&v12),
            bi_normal.dot(&v22),
            bi_normal.dot(&v32),
        );
    }

    cos_angle(tangent, new_tangent) >= an_inf_cos
}
