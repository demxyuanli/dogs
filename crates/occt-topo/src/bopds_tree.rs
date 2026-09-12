//! Sub-shape extraction helpers used while filling BopdsDS.
use std::collections::HashSet;
use std::sync::Arc;

use crate::abs::ShapeType;
use crate::shape::{Edge, Face, Shell, Solid, TopoShape, Wire};
use crate::topo_tools_full::{edges_of, edges_of_wire, faces_of, shapes_of, vertices_of, wires_of_face};

// ---------------------------------------------------------------------------
// Sub-shape extraction helpers
// ---------------------------------------------------------------------------

/// Whether two parameter ranges overlap within `tol`.
pub(crate) fn ranges_overlap(a1: f64, a2: f64, b1: f64, b2: f64, tol: f64) -> bool {
    let (a1, a2) = if a1 <= a2 { (a1, a2) } else { (a2, a1) };
    let (b1, b2) = if b1 <= b2 { (b1, b2) } else { (b2, b1) };
    a1 <= b2 + tol && b1 <= a2 + tol
}

/// Direct vertex children of a shape (used for edges).
pub(crate) fn direct_vertex_children(shape: &TopoShape) -> Vec<TopoShape> {
    shape
        .tshape
        .read()
        .unwrap()
        .children
        .iter()
        .filter(|h| h.shape_type() == ShapeType::Vertex)
        .cloned()
        .collect()
}

/// All direct children of a shape.
pub(crate) fn direct_children(shape: &TopoShape) -> Vec<TopoShape> {
    shape
        .tshape
        .read()
        .unwrap()
        .children
        .clone()
}

/// The two boundary vertices of an edge.
pub(crate) fn edge_vertex_shapes(edge: &Edge) -> Vec<TopoShape> {
    direct_vertex_children(&edge.0)
}

/// The edges of a wire.
pub(crate) fn wire_edge_shapes(wire: &Wire) -> Vec<TopoShape> {
    edges_of_wire(wire).into_iter().map(|e| e.0).collect()
}

/// The boundary sub-shapes of a face: the distinct edges of its wires plus
/// its distinct vertices (wires themselves are not boundary sub-shapes).
pub(crate) fn face_sub_shapes(face: &Face) -> Vec<TopoShape> {
    let mut seen = HashSet::new();
    let mut out = Vec::new();
    for w in wires_of_face(face) {
        for e in edges_of_wire(&w) {
            if seen.insert(Arc::as_ptr(&e.0.tshape) as usize) {
                out.push(e.0);
            }
        }
    }
    for v in vertices_of(&face.0) {
        if seen.insert(Arc::as_ptr(&v.0.tshape) as usize) {
            out.push(v.0);
        }
    }
    out
}

/// The faces of a shell.
pub(crate) fn shell_face_shapes(shell: &Shell) -> Vec<TopoShape> {
    faces_of(&shell.0).into_iter().map(|f| f.0).collect()
}

/// The shells of a solid (structural children).
pub(crate) fn solid_shell_shapes(solid: &Solid) -> Vec<TopoShape> {
    shapes_of(&solid.0, ShapeType::Shell)
}

/// The prepared boundary sub-shapes of a solid: its faces and edges.
pub(crate) fn solid_face_edge_shapes(solid: &Solid) -> Vec<TopoShape> {
    let mut seen = HashSet::new();
    let mut out = Vec::new();
    for f in faces_of(&solid.0) {
        if seen.insert(Arc::as_ptr(&f.0.tshape) as usize) {
            out.push(f.0);
        }
    }
    for e in edges_of(&solid.0) {
        if seen.insert(Arc::as_ptr(&e.0.tshape) as usize) {
            out.push(e.0);
        }
    }
    out
}

/// The structural children used to walk the whole subtree of a shape while
/// appending it (a solid walks its shells, a face its edges and vertices).
pub(crate) fn structural_children(shape: &TopoShape) -> Vec<TopoShape> {
    match shape.shape_type() {
        ShapeType::Vertex => Vec::new(),
        ShapeType::Edge => edge_vertex_shapes(&Edge(shape.clone())),
        ShapeType::Wire => wire_edge_shapes(&Wire(shape.clone())),
        ShapeType::Face => face_sub_shapes(&Face(shape.clone())),
        ShapeType::Shell => shell_face_shapes(&Shell(shape.clone())),
        ShapeType::Solid => solid_shell_shapes(&Solid(shape.clone())),
        _ => direct_children(shape),
    }
}

/// The prepared boundary sub-shapes of a shape, in the form the DS stores
/// (matching OCCT after `BOPDS_DS::Init` + `prepareFaces`/`prepareSolids`):
/// a solid lists faces and edges, a shell its faces, a face its edges and
/// vertices, an edge its vertices.
pub(crate) fn prepared_sub_shapes(shape: &TopoShape) -> Vec<TopoShape> {
    match shape.shape_type() {
        ShapeType::Vertex => Vec::new(),
        ShapeType::Edge => edge_vertex_shapes(&Edge(shape.clone())),
        ShapeType::Wire => wire_edge_shapes(&Wire(shape.clone())),
        ShapeType::Face => face_sub_shapes(&Face(shape.clone())),
        ShapeType::Shell => shell_face_shapes(&Shell(shape.clone())),
        ShapeType::Solid => solid_face_edge_shapes(&Solid(shape.clone())),
        _ => direct_children(shape),
    }
}
