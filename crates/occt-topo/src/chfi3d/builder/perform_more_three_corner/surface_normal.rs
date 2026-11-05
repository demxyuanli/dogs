//! `SummarizeNormal` — one face normal added to the running sum at the corner
//! vertex.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_CnCrn.cxx:1113-1141`
//! (OCCT 8.0.0).
//!
//! The sum is a `gp_Vec&` out-parameter in OCCT that the caller normalizes after
//! every addition; the port takes it as `&mut GpVec`. The early `return` when
//! the normal is degenerate leaves the sum as it was — it is *not* renormalized,
//! which is what the OCCT body does.

use occt_core::gp::GpVec;
use occt_core::precision::PCONFUSION;

use crate::abs::Orientation;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::{edge_vertices, is_same};

/// `SummarizeNormal(V1, Fcur, Ecur, SumFaceNormalAtV1)` (`cxx:1113-1141`).
///
/// `Err` carries the `BRep_Tool::Surface(Fcur)` dereference OCCT performs
/// unconditionally (`cxx:1127`).
pub fn summarize_normal(
    sum: &mut GpVec,
    v1: &TopoShape,
    fcur: &Face,
    ecur: &Edge,
) -> Result<(), String> {
    // `BRep_Tool::UVPoints(Ecur, Fcur, uv1, uv2)` (`cxx:1119`).
    let (uv1, uv2) = BRepTool::uv_points(ecur, fcur);
    // `if (!V1.IsSame(TopExp::FirstVertex(Ecur))) uv1 = uv2;` (`cxx:1120-1123`).
    let first = edge_vertices(ecur).0;
    let uv = match &first {
        Some(f) if is_same(&f.0, v1) => uv1,
        _ => uv2,
    };

    // `BRep_Tool::Surface(Fcur)->D1(uv1.X(), uv1.Y(), P, d1U, d1V)`
    // (`cxx:1125-1127`).
    let surface = BRepTool::face_surface(fcur)
        .ok_or_else(|| "SummarizeNormal: the face has no surface".to_string())?;
    let (_, d1u, d1v) = surface.d1(uv.x(), uv.y());

    // `gp_Vec N = d1U.Crossed(d1V);` then the `TopAbs_REVERSED` reversal
    // (`cxx:1128-1132`).
    let mut n = d1u.crossed(&d1v);
    if fcur.0.orientation() == Orientation::Reversed {
        n.reverse();
    }

    // `if (N.SquareMagnitude() <= Precision::PConfusion()) return;`
    // (`cxx:1134-1137`).
    if n.square_magnitude() <= PCONFUSION {
        return Ok(());
    }

    // `SumFaceNormalAtV1 += N.Normalized(); SumFaceNormalAtV1.Normalize();`
    // (`cxx:1139-1140`).
    *sum = sum.added(&n.normalized());
    sum.normalize();
    Ok(())
}
