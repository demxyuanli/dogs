//! CSLib ports plus a **port-internal** surface classifier.
//!
//! **Ported from `CSLib/`**: `class2d` (`CSLib_Class2d`), `normal` /
//! `dn_normal` (`CSLib::Normal` / `CSLib::DNNormal`), `poly_def`
//! (`CSLib_NormalPolyDef`).
//!
//! **Provenance (audit A9)**: `CSLibResult`, `classify_point` and
//! `sphere_normal` below are **not** OCCT translations — the `CSLib` package
//! contains only `CSLib.cxx/hxx`, `CSLib_Class2d`, `CSLib_NormalPolyDef` and
//! the status enums (`DerivativeStatus`, `NormalStatus`), i.e. nothing that
//! classifies a point against a surface. OCCT classifies points against faces
//! through `BRepClass_FaceClassifier` / `IntTools_FClass2d` on the pcurve
//! (ported in `occt-topo/fclass2d`).

pub mod class2d;
pub mod dn_normal;
pub mod normal;
pub mod poly_def;
pub use class2d::{Class2d, Class2dResult};
pub use dn_normal::{dn_normal, dnnuv, dnnuv2};
pub use normal::{normal_d1_mag, normal_d2, normal_max_order, CSLibNormalStatus};
pub use poly_def::NormalPolyDef;

/// Result of a surface classification test.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CSLibResult { Inside, Outside, On }

/// Classify point relative to closed surface using normal at nearest point.
/// normal: outward-pointing surface normal at the reference point.
/// point: the test point.
/// tolerance: distance tolerance for "on" classification.
pub fn classify_point(point: &[f64; 3], surface_point: &[f64; 3], normal: &[f64; 3], tolerance: f64) -> CSLibResult {
    let dx = point[0] - surface_point[0];
    let dy = point[1] - surface_point[1];
    let dz = point[2] - surface_point[2];
    let dist = dx * normal[0] + dy * normal[1] + dz * normal[2];

    if dist.abs() <= tolerance { CSLibResult::On }
    else if dist > 0.0 { CSLibResult::Outside }
    else { CSLibResult::Inside }
}

/// Compute the normal at a point on a sphere (center at origin).
pub fn sphere_normal(center: &[f64; 3], point: &[f64; 3]) -> [f64; 3] {
    let mut n = [point[0]-center[0], point[1]-center[1], point[2]-center[2]];
    let len = (n[0]*n[0] + n[1]*n[1] + n[2]*n[2]).sqrt();
    if len > 1e-30 { n[0] /= len; n[1] /= len; n[2] /= len; }
    n
}
