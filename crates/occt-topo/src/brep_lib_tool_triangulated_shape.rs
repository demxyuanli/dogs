//! `BRepLib_ToolTriangulatedShape`.
//! Source: `BRepLib_ToolTriangulatedShape.cxx:27-80`.

use occt_core::gp::{GpDir, GpPnt, GpPnt2d, GpXyz};
use occt_core::poly::connect::PolyConnect;
use occt_core::poly::triangulation::Triangle;
use occt_core::poly::poly_triangulation::PolyTriangulation;
use occt_geom::Surface;

/// `Precision::Confusion()`: `constexpr double aTol = Precision::Confusion();`
/// (`BRepLib_ToolTriangulatedShape.cxx:46`).
pub const CONFUSION: f64 = 1.0e-7;

/// `BRepLib_ToolTriangulatedShape::ComputeNormals(theFace, theTris,
/// thePolyConnect)` (`BRepLib_ToolTriangulatedShape.cxx:27-80`): one normal per
/// node of a triangulation that carries none.
///
/// `surface` is `BRep_Tool::Surface(aZeroFace)` (`cxx:37`, the forward
/// oriented, unlocated face surface) and `uv_nodes` is
/// `Poly_Triangulation::HasUVNodes` (`cxx:39`); `nodes` / `triangles` are the
/// triangulation itself. With either of the first two missing OCCT falls into
/// `Poly::ComputeNormals` (`cxx:40-44`).
pub fn compute_normals(
    surface: Option<&dyn Surface>,
    nodes: &[GpPnt],
    uv_nodes: &[GpPnt2d],
    triangles: &[Triangle],
) -> Vec<GpPnt> {
    let mut poly =
        PolyTriangulation::from_parts_uv(nodes.to_vec(), uv_nodes.to_vec(), triangles.to_vec());
    let surf = match surface {
        Some(s) if !uv_nodes.is_empty() => s,
        _ => {
            // `Poly::ComputeNormals(theTris)` (`cxx:42`) ->
            // `Poly_Triangulation::ComputeNormals` (`Poly_Triangulation.cxx:436-470`).
            poly.compute_normals();
            return (0..poly.nb_nodes()).map(|i| poly.normal(i)).collect();
        }
    };

    // `Poly_Connect& thePolyConnect` (`cxx:28`): OCCT reloads it when it does
    // not point at this triangulation (`cxx:54-57`); here it is built once.
    let mut connect = PolyConnect::from_triangulation(&poly);
    let mut normals = Vec::with_capacity(poly.nb_nodes());
    for i in 0..poly.nb_nodes() {
        // `if (GeomLib::NormEstim(aSurf, theTris->UVNode(aNodeIter), aTol,
        // aNorm) > 1) { ... }` (`cxx:53`): status 0/1 leaves the estimate in
        // `aNorm`, 2/3 takes the flat fallback below. `theTris->SetNormal`
        // (`cxx:78`) then stores it on the node; `RWMesh_FaceIterator::normal`
        // (`cxx:59-67`) reads exactly those stored values back for export.
        let uv = uv_nodes[i];
        let (status, estimate) =
            crate::meshing::geomlib_norm::norm_estim(surf, uv.x(), uv.y(), CONFUSION);
        let normal = match estimate {
            Some(dir) if status <= 1 => GpPnt::new(dir.x(), dir.y(), dir.z()),
            _ => flat_normal(&poly, &mut connect, i),
        };
        normals.push(normal);
    }
    normals
}

/// The `cxx:59-78` fallback: the normalised sum of the unit normals of the
/// triangles sharing the node (`Poly_Connect` fan, `cxx:61-68`), `gp::DZ()`
/// when that sum is null (`cxx:75-76`).
fn flat_normal(poly: &PolyTriangulation, connect: &mut PolyConnect, node: usize) -> GpPnt {
    let mut eq_plan = GpXyz::zero();
    connect.initialize(node);
    while connect.more() {
        let tri = poly.triangle(connect.value());
        let p0 = poly.node(tri.n0).coord;
        let p1 = poly.node(tri.n1).coord;
        let p2 = poly.node(tri.n2).coord;
        let v1 = p1.subtracted(&p0);
        let v2 = p2.subtracted(&p1);
        let vv = v1.crossed(&v2);
        let modulus = vv.modulus();
        if modulus >= CONFUSION {
            eq_plan = eq_plan.add(&vv.divided(modulus));
        }
        connect.next();
    }
    if eq_plan.modulus() > CONFUSION {
        if let Ok(dir) = GpDir::new(eq_plan.x, eq_plan.y, eq_plan.z) {
            return GpPnt::new(dir.x(), dir.y(), dir.z());
        }
    }
    // `gp::DZ()`.
    GpPnt::new(0.0, 0.0, 1.0)
}
