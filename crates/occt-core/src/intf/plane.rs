//! Port of the static helpers of `Intf`
//! (`src/ModelingAlgorithms/TKGeomAlgo/Intf/Intf.hxx`, `Intf.cxx`).
//!
//! * `Intf::PlaneEquation`  `Intf.cxx:24-44`
//! * `Intf::Contain`        `Intf.cxx:48-54`

use crate::gp::{GpPnt, GpXyz};
use crate::precision::REAL_SMALL;

/// `Intf::PlaneEquation(P1, P2, P3, NormalVector, PolarDistance)`
/// (`Intf.cxx:24-44`).
///
/// Returns the unit normal of the triangle `(P1, P2, P3)` and the polar
/// distance `Normal . P1`. When the raw normal is shorter than
/// `gp::Resolution()` the normal is returned *unnormalized* and the polar
/// distance is `0.`, exactly as OCCT leaves it (`Intf.cxx:35-38`).
pub fn plane_equation(p1: &GpPnt, p2: &GpPnt, p3: &GpPnt) -> (GpXyz, f64) {
    let v1 = p2.xyz().subtracted(p1.xyz());
    let v2 = p3.xyz().subtracted(p2.xyz());
    let v3 = p1.xyz().subtracted(p3.xyz());
    let mut normal = v1
        .crossed(&v2)
        .add(&v2.crossed(&v3))
        .add(&v3.crossed(&v1));
    let norm_len = normal.modulus();
    if norm_len < REAL_SMALL {
        (normal, 0.0)
    } else {
        normal = normal.divided(norm_len);
        let polar_distance = normal.dot(p1.xyz());
        (normal, polar_distance)
    }
}

/// `Intf::Contain(P1, P2, P3, ThePnt)` (`Intf.cxx:48-54`).
///
/// True when `ThePnt` lies inside the triangle `(P1, P2, P3)` as
/// seen from the three oriented edge cross products.
pub fn contain(p1: &GpPnt, p2: &GpPnt, p3: &GpPnt, the_pnt: &GpPnt) -> bool {
    let v1 = p2.xyz().subtracted(p1.xyz()).crossed(&the_pnt.xyz().subtracted(p1.xyz()));
    let v2 = p3.xyz().subtracted(p2.xyz()).crossed(&the_pnt.xyz().subtracted(p2.xyz()));
    let v3 = p1.xyz().subtracted(p3.xyz()).crossed(&the_pnt.xyz().subtracted(p3.xyz()));
    v1.dot(&v2) >= 0.0 && v2.dot(&v3) >= 0.0 && v3.dot(&v1) >= 0.0
}
