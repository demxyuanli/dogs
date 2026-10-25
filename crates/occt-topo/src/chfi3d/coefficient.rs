//! `ChFi3d_Coefficient` — the `(du, dv)` of a 3D vector in a surface frame.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_0.cxx:1451-1463`,
//! declared in `ChFi3d_Builder_0.hxx`.

use occt_core::gp::GpVec;

/// `ChFi3d_Coefficient(V3d, D1u, D1v, DU, DV)` (`cxx:1451-1463`): solves
/// `V3d = DU * D1u + DV * D1v` through the first fundamental form. No
/// degeneracy test — a null Gram determinant divides by zero exactly as OCCT
/// does.
pub fn chfi3d_coefficient(v3d: &GpVec, d1u: &GpVec, d1v: &GpVec) -> (f64, f64) {
    let aa = d1u.square_magnitude();
    let bb = d1u.dot(d1v);
    let cc = d1v.square_magnitude();
    let dd = d1u.dot(v3d);
    let ee = d1v.dot(v3d);
    let delta = aa * cc - bb * bb;
    let du = (dd * cc - ee * bb) / delta;
    let dv = (aa * ee - bb * dd) / delta;
    (du, dv)
}
