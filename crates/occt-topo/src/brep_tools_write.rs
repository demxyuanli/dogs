//! `BRepTools` triangulation dump -- `BRepTools::Write`'s triangle output.
//!
//! OCCT writes the stored `Poly_Triangulation` of each face; the port has no
//! stored triangulation, so it runs the one mesher (`BRepMesh_IncrementalMesh`,
//! `meshing::incremental_mesh`) and prints the result. That makes this an IO
//! concern above the `BRep_Tool` data queries, which is why it does not live in
//! `occt-topology::brep_tools`.

use crate::shape::TopoShape;

/// Dump the incremental-mesh triangulation of `shape` as ASCII
/// (`BRepTools` triangulation write). Each triangle is `v0 v1 v2` after a
/// vertex list. `deflection` is the mesher linear deflection.
pub fn write_triangulation(shape: &TopoShape, deflection: f64) -> String {
    // One mesher only (`BRepMesh_IncrementalMesh`, via `meshing::incremental_mesh`).
    // A shape whose triangulation fails contributes no triangles; OCCT has no
    // second mesher to fall back to (audit A17/A18/A19).
    let mesh = crate::meshing::incremental_mesh::incremental_mesh_to_shape_mesh(shape, deflection)
        .unwrap_or_else(|_| crate::mesh::ShapeMesh {
            vertices: Vec::new(),
            triangles: Vec::new(),
            source_shape: crate::abs::ShapeType::Shape,
        });
    let mut out = String::new();
    out.push_str(&format!("Triangulation {}\n", mesh.vertices.len()));
    for p in &mesh.vertices {
        out.push_str(&format!("{} {} {}\n", p.x(), p.y(), p.z()));
    }
    out.push_str(&format!("Triangles {}\n", mesh.triangles.len()));
    for t in &mesh.triangles {
        out.push_str(&format!("{} {} {}\n", t.n0, t.n1, t.n2));
    }
    out
}
