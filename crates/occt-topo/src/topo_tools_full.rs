//! Complete topology traversal and classification — a fuller port of
//! `TopExp.hxx` / `TopTools.hxx` than the lightweight `topexp.rs` module.
//!
//! These helpers walk the real children tree of a `TShape`, uniquifying by
//! `TShape` identity (shared edges/vertices appear once), and mirror the
//! OCCT `TopExp::MapShapes`, `TopExp::FirstVertex/LastVertex`,
//! `TopTools::MapOfShape` idioms.

use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use occt_core::gp::GpPnt;

use crate::abs::{Orientation, ShapeType};
use crate::iterator::cumulated_children;
use crate::shape::{Edge, Face, TopoShape, Vertex, Wire};
use crate::tgeometry::GeometryRegistry;

/// Map every distinct sub-shape of the given types, keyed by type.
/// Mirrors `TopExp::MapShapes` (which also accepts an `Allocator`).
pub fn map_shapes(shape: &TopoShape, types: &[ShapeType]) -> HashMap<ShapeType, Vec<TopoShape>> {
    let mut out: HashMap<ShapeType, Vec<TopoShape>> = HashMap::new();
    for &t in types {
        out.entry(t).or_default();
    }
    let mut seen = HashSet::new();
    let mut stack = vec![shape.clone()];
    while let Some(s) = stack.pop() {
        // Key on TShape identity only — OCCT's TopTools_ShapeMapHasher (used by
        // TopTools_IndexedMapOfShape / TopExp::MapShapes) keys with IsSame, so a
        // shared edge used Forward in one face and Reversed in another (e.g. the
        // box's 12 edges) counts once. The orientation is carried on the returned
        // shape view; a seam's two opposite-oriented occurrences stay reachable
        // through `edges_of_wire` (direct children), which does not deduplicate.
        let key = Arc::as_ptr(&s.tshape) as usize;
        if !seen.insert(key) {
            continue;
        }
        if types.contains(&s.shape_type()) {
            out.get_mut(&s.shape_type()).unwrap().push(s.clone());
        }
        for k in cumulated_children(&s).into_iter().rev() {
            stack.push(k);
        }
    }
    out
}

/// Every distinct sub-shape of `shape` in pre-order (root excluded, like
/// `TopExp_Explorer`). Includes all topological types.
pub fn all_subshapes(shape: &TopoShape) -> Vec<TopoShape> {
    let mut out = Vec::new();
    let mut seen = HashSet::new();
    let mut stack = vec![shape.clone()];
    while let Some(s) = stack.pop() {
        let key = Arc::as_ptr(&s.tshape) as usize;
        if !seen.insert(key) {
            continue;
        }
        out.push(s.clone());
        for k in cumulated_children(&s).into_iter().rev() {
            stack.push(k);
        }
    }
    out
}

/// Per-type occurrence counts (distinct shapes), root included.
pub fn shape_counts(shape: &TopoShape) -> HashMap<ShapeType, usize> {
    let mut out = HashMap::new();
    for s in all_subshapes(shape) {
        *out.entry(s.shape_type()).or_insert(0) += 1;
    }
    out
}

/// Distinct sub-shapes of one type.
pub fn shapes_of(shape: &TopoShape, t: ShapeType) -> Vec<TopoShape> {
    map_shapes(shape, &[t]).remove(&t).unwrap_or_default()
}

/// Distinct vertices of `shape`.
pub fn vertices_of(shape: &TopoShape) -> Vec<Vertex> {
    shapes_of(shape, ShapeType::Vertex).into_iter().map(Vertex).collect()
}

/// Distinct edges of `shape`.
pub fn edges_of(shape: &TopoShape) -> Vec<Edge> {
    shapes_of(shape, ShapeType::Edge).into_iter().map(Edge).collect()
}

/// Distinct wires of `shape`.
pub fn wires_of(shape: &TopoShape) -> Vec<Wire> {
    shapes_of(shape, ShapeType::Wire).into_iter().map(Wire).collect()
}

/// Distinct faces of `shape`.
pub fn faces_of(shape: &TopoShape) -> Vec<Face> {
    shapes_of(shape, ShapeType::Face).into_iter().map(Face).collect()
}

/// Direct child wires of a face (the face's boundary wires).
pub fn wires_of_face(face: &Face) -> Vec<Wire> {
    cumulated_children(&face.0)
        .into_iter()
        .filter(|s| s.shape_type() == ShapeType::Wire)
        .map(Wire)
        .collect()
}

/// Direct child edges of a wire, with `TopoDS_Iterator` Compose so a reversed
/// parent wire yields reversed edge views. Does not uniquify: a seam edge
/// stored twice (Forward and Reversed) is returned twice.
pub fn edges_of_wire(wire: &Wire) -> Vec<Edge> {
    cumulated_children(&wire.0)
        .into_iter()
        .filter(|s| s.shape_type() == ShapeType::Edge)
        .map(Edge)
        .collect()
}

/// The two endpoint vertices of an edge (`TopExp::FirstVertex` / `LastVertex`,
/// `CumOri = false`). First is the child stored `FORWARD`; last is `REVERSED`.
pub fn edge_vertices(edge: &Edge) -> (Option<Vertex>, Option<Vertex>) {
    let stored = edge
        .0
        .tshape
        .read()
        .expect("poisoned TShape lock")
        .children
        .clone();
    let verts: Vec<TopoShape> = stored
        .into_iter()
        .filter(|s| s.shape_type() == ShapeType::Vertex)
        .collect();
    let first = verts
        .iter()
        .find(|s| s.orientation() == Orientation::Forward)
        .cloned()
        .map(Vertex);
    let last = verts
        .iter()
        .find(|s| s.orientation() == Orientation::Reversed)
        .cloned()
        .map(Vertex);
    if first.is_some() || last.is_some() {
        return (first, last);
    }
    // Vertices added without FORWARD/REVERSED storage (pre-MakeEdge write).
    let first = verts.first().cloned().map(Vertex);
    let last = verts.get(1).or_else(|| verts.last()).cloned().map(Vertex);
    (first, last)
}

/// Whether two shapes reference the same `TShape` data (identical, not equal).
pub fn is_same(a: &TopoShape, b: &TopoShape) -> bool {
    Arc::ptr_eq(&a.tshape, &b.tshape)
}

/// Location of a vertex from the geometry side-table (`BRep_Tool::Pnt`).
pub fn vertex_position(v: &Vertex) -> GpPnt {
    GeometryRegistry::global().vertex_point(&v.0)
}

/// Number of distinct shapes of each type under `shape` (verbose form).
pub fn count_types(shape: &TopoShape) -> Vec<(ShapeType, usize)> {
    let mut v: Vec<(ShapeType, usize)> = shape_counts(shape).into_iter().collect();
    v.sort_by_key(|(t, _)| *t as u8);
    v
}

/// A wire is closed if its edges chain end-to-end (each edge's end equals the
/// next edge's start, in child order) and the last edge's end returns to the
/// first edge's start. Mirrors the ordered `TopExp` convention; edges stored
/// in reversed order relative to the wire traversal will read as open.
pub fn wire_is_closed(wire: &Wire) -> bool {
    let edges = edges_of_wire(wire);
    if edges.len() < 3 {
        return false;
    }
    let mut prev_end: Option<GpPnt> = None;
    for e in &edges {
        let (a, b) = edge_vertices(e);
        let (Some(a), Some(b)) = (a, b) else { return false };
        let pa = vertex_position(&a);
        let pb = vertex_position(&b);
        if let Some(pe) = prev_end {
            if pa.distance(&pe) > 1e-9 {
                return false;
            }
        }
        prev_end = Some(pb);
    }
    if let (Some(first), Some(pe)) = (edges.first(), prev_end) {
        let (a, _) = edge_vertices(first);
        return a.map_or(false, |a| vertex_position(&a).distance(&pe) < 1e-9);
    }
    false
}

/// Structural sanity: every vertex child sits on its parent edge, every edge
/// child on its parent wire, every wire on a face, every face on a shell,
/// every shell on a solid. Mirrors `BRepCheck_Analyser`'s basic pass.
pub fn structure_is_valid(shape: &TopoShape) -> bool {
    fn valid_child(parent: ShapeType, child: ShapeType) -> bool {
        matches!(
            (parent, child),
            (ShapeType::Compound, _)
                | (ShapeType::Wire, ShapeType::Edge)
                | (ShapeType::Face, ShapeType::Wire)
                | (ShapeType::Shell, ShapeType::Face)
                | (ShapeType::Solid, ShapeType::Shell)
                | (ShapeType::Edge, ShapeType::Vertex)
        )
    }
    let mut stack = vec![shape.clone()];
    let mut seen = HashSet::new();
    while let Some(s) = stack.pop() {
        if !seen.insert(Arc::as_ptr(&s.tshape)) {
            continue;
        }
        let parent = s.shape_type();
        for k in cumulated_children(&s) {
            let child = k.shape_type();
            if !valid_child(parent, child) {
                return false;
            }
            stack.push(k);
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::builder::TopoBuilder;
    use crate::primitives::BRepPrimBox;
    use occt_core::gp::GpPnt;

    #[test]
    fn box_counts_and_types() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let c = shape_counts(&b.solid.0);
        assert_eq!(c[&ShapeType::Vertex], 8);
        assert_eq!(c[&ShapeType::Edge], 12);
        assert_eq!(c[&ShapeType::Face], 6);
        assert_eq!(c[&ShapeType::Wire], 6);
        assert_eq!(c[&ShapeType::Shell], 1);
        assert_eq!(c[&ShapeType::Solid], 1);
    }

    #[test]
    fn map_shapes_and_specific_collectors() {
        let b = BRepPrimBox::make_box(1.0, 2.0, 3.0);
        let map = map_shapes(&b.solid.0, &[ShapeType::Vertex, ShapeType::Face]);
        assert_eq!(map[&ShapeType::Vertex].len(), 8);
        assert_eq!(map[&ShapeType::Face].len(), 6);
        assert!(vertices_of(&b.solid.0).len() == 8);
        assert!(edges_of(&b.solid.0).len() == 12);
        assert!(faces_of(&b.solid.0).len() == 6);
    }

    #[test]
    fn wire_closed_and_edge_vertices() {
        let b = TopoBuilder::new();
        // A closed square wire.
        let p = [GpPnt::new(0.,0.,0.), GpPnt::new(1.,0.,0.),
                 GpPnt::new(1.,1.,0.), GpPnt::new(0.,1.,0.)];
        let e1 = b.make_edge_segment(&p[0], &p[1]);
        let e2 = b.make_edge_segment(&p[1], &p[2]);
        let e3 = b.make_edge_segment(&p[2], &p[3]);
        let e4 = b.make_edge_segment(&p[3], &p[0]);
        let wire = b.make_wire(&[e1.clone(), e2, e3, e4]);
        assert!(wire_is_closed(&wire));

        let (a, z) = edge_vertices(&e1);
        assert!(vertex_position(&a.unwrap()).is_equal(&p[0]));
        assert!(vertex_position(&z.unwrap()).is_equal(&p[1]));
    }

    #[test]
    fn structure_valid_box() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        assert!(structure_is_valid(&b.solid.0));
    }

    #[test]
    fn is_same_and_wires_of_face() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let faces = faces_of(&b.solid.0);
        assert!(!faces.is_empty());
        let w = wires_of_face(&faces[0]);
        assert_eq!(w.len(), 1);
        let ws = edges_of_wire(&w[0]);
        assert_eq!(ws.len(), 4);
        assert!(is_same(&faces[0].0, &faces[0].0));
    }
}
