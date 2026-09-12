//! `BRepTools` remainder: UV bounds, outer wire, triangulation dump.
//!
//! Source: `ModelingData/TKBRep/BRepTools`. [`add_uv_bounds`] delegates to
//! [`crate::brep_uv_bounds`] (including the periodic B-spline extra probe).
//! [`write_triangulation`] dumps the incremental-mesh triangles of a shape
//! (`BRepTools::Write` triangulation block).

use occt_core::bnd::BndBox2d;

use crate::brep_tool::BRepTool;
use crate::brep_uv_bounds::{add_uv_bounds_edge, add_uv_bounds_face, add_uv_bounds_wire};
use crate::shape::{Edge, Face, TopoShape, Wire};
use crate::topo_tools_full::wires_of_face;

/// `BRepTools::AddUVBounds(Face, Box2d)`.
pub fn add_uv_bounds(face: &Face, box2d: &mut BndBox2d) {
    add_uv_bounds_face(face, box2d);
}

/// `BRepTools::AddUVBounds(Face, Wire, Box2d)`.
pub fn add_uv_bounds_on_wire(face: &Face, wire: &Wire, box2d: &mut BndBox2d) {
    add_uv_bounds_wire(face, wire, box2d);
}

/// `BRepTools::AddUVBounds(Face, Edge, Box2d)`.
pub fn add_uv_bounds_on_edge(face: &Face, edge: &Edge, box2d: &mut BndBox2d) {
    add_uv_bounds_edge(face, edge, box2d);
}

/// `BRepTools::UVBounds(Face)` → `(u_min, u_max, v_min, v_max)`.
pub fn uv_bounds(face: &Face) -> (f64, f64, f64, f64) {
    let mut b = BndBox2d::new();
    add_uv_bounds_face(face, &mut b);
    match b.get() {
        Some((u0, v0, u1, v1)) => (u0, u1, v0, v1),
        None => (0.0, 1.0, 0.0, 1.0),
    }
}

/// `BRepTools::OuterWire`. The wire whose UV box contains every other wire
/// of the face; if none, the first wire.
pub fn outer_wire(face: &Face) -> Option<Wire> {
    let wires = wires_of_face(face);
    if wires.is_empty() {
        return None;
    }
    let mut best = 0usize;
    let mut best_span = -1.0;
    for (i, w) in wires.iter().enumerate() {
        let mut b = BndBox2d::new();
        add_uv_bounds_wire(face, w, &mut b);
        let span = match b.get() {
            Some((u0, v0, u1, v1)) => (u1 - u0).abs() * (v1 - v0).abs(),
            None => 0.0,
        };
        if span > best_span {
            best_span = span;
            best = i;
        }
    }
    Some(wires[best].clone())
}

/// `BRepTools::Map3DEdges` — unique 3-D edges of `shape`.
pub fn map_3d_edges(shape: &TopoShape) -> Vec<Edge> {
    crate::topo_tools_full::edges_of(shape)
}

/// Dump the incremental-mesh triangulation of `shape` as ASCII
/// (`BRepTools` triangulation write). Each triangle is `v0 v1 v2` after a
/// vertex list. `deflection` is the mesher linear deflection.
pub fn write_triangulation(shape: &TopoShape, deflection: f64) -> String {
    let mesh = crate::meshing::incremental_mesh::incremental_mesh_to_shape_mesh(shape, deflection)
        .or_else(|_| crate::brepmesh::incremental_mesh(shape, deflection).map(|im| im.mesh))
        .unwrap_or_else(|_| crate::shape_mesh::mesh_shape(shape, deflection));
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

/// `BRepTools::Compare` vertices: same point within the larger tolerance.
pub fn compare_vertices(a: &crate::shape::Vertex, b: &crate::shape::Vertex) -> bool {
    let pa = BRepTool::vertex_point(a);
    let pb = BRepTool::vertex_point(b);
    let tol = BRepTool::vertex_tolerance(a).max(BRepTool::vertex_tolerance(b));
    pa.distance(&pb) <= tol
}
