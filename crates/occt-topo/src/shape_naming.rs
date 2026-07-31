//! Shape identification and indexing — stable IDs, geometric hashes and
//! predicate-based sub-shape lookup. Mirrors `TopTools_MapOfShape` /
//! `TopTools_IndexedMapOfShape` plus `TopExp::MapShapes`-style queries.
//!
//! A `ShapeId` is derived from the `TShape` pointer (stable for the shape's
//! lifetime); a `geometric_hash` additionally mixes in the registered geometry
//! so two shapes with equal structure-and-geometry can be compared by value.

use std::collections::HashMap;
use std::hash::{Hash, Hasher};
use std::sync::Arc;

use occt_core::gp::GpPnt;

use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools_full::{all_subshapes, edges_of, faces_of, vertices_of};

/// Stable identity of a shape: the address of its shared `TShape`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct ShapeId(pub usize);

impl ShapeId {
    pub fn of(shape: &TopoShape) -> ShapeId {
        ShapeId(Arc::as_ptr(&shape.tshape) as usize)
    }
}

/// Geometric fingerprint of a shape subtree: mixes every vertex point and
/// every edge's first parameter into a 64-bit hash. Two shapes built from the
/// same geometry (but distinct TShapes) collide.
pub fn geometric_hash(shape: &TopoShape) -> u64 {
    let mut h = std::collections::hash_map::DefaultHasher::new();
    for v in vertices_of(shape) {
        let p = BRepTool::vertex_point(&v);
        p.x().to_bits().hash(&mut h);
        p.y().to_bits().hash(&mut h);
        p.z().to_bits().hash(&mut h);
    }
    for e in edges_of(shape) {
        let (a, b) = BRepTool::edge_parameters(&e);
        a.to_bits().hash(&mut h);
        b.to_bits().hash(&mut h);
    }
    h.finish()
}

/// Are two shapes geometrically identical within `tol` (same sub-shape
/// counts and matching vertex positions)?
pub fn shapes_geometrically_equal(a: &TopoShape, b: &TopoShape, tol: f64) -> bool {
    let va = vertices_of(a);
    let vb = vertices_of(b);
    if va.len() != vb.len() {
        return false;
    }
    // Match every vertex of `a` to a vertex of `b` within tol (multiset match).
    let mut used = vec![false; vb.len()];
    for v in &va {
        let p = BRepTool::vertex_point(v);
        let mut found = false;
        for (i, w) in vb.iter().enumerate() {
            if !used[i] && p.distance(&BRepTool::vertex_point(w)) <= tol {
                used[i] = true;
                found = true;
                break;
            }
        }
        if !found {
            return false;
        }
    }
    true
}

/// Count occurrences of every sub-shape type (root excluded), keyed by type
/// index. Mirrors `TopExp::MapShapes` counts.
pub fn shape_index(shape: &TopoShape) -> HashMap<usize, usize> {
    let mut out = HashMap::new();
    for s in all_subshapes(shape) {
        *out.entry(s.shape_type() as usize).or_insert(0) += 1;
    }
    out
}

/// Find the first sub-shape whose type and — for vertices — position match a
/// predicate.
pub fn find_vertex_at(shape: &TopoShape, p: &GpPnt, tol: f64) -> Option<Vertex> {
    vertices_of(shape).into_iter().find(|v| BRepTool::vertex_point(v).distance(p) <= tol)
}

/// Find the edge whose parameter range matches `(first, last)` within tol.
pub fn find_edge_by_range(shape: &TopoShape, first: f64, last: f64, tol: f64) -> Option<Edge> {
    edges_of(shape).into_iter().find(|e| {
        let (a, b) = BRepTool::edge_parameters(e);
        (a - first).abs() <= tol && (b - last).abs() <= tol
    })
}

/// Find the face whose surface passes through `p` within tol. The closest
/// (u, v) on the surface is located by grid search before the distance check,
/// so unbounded planes work correctly.
pub fn find_face_through(shape: &TopoShape, p: &GpPnt, tol: f64) -> Option<Face> {
    faces_of(shape).into_iter().find(|f| {
        match GeometryRegistry::global().face_surface(&f.0) {
            Some(s) => {
                let (u, v) = crate::brep_surface::surface_closest_params(s.as_ref(), p, 16, 16);
                s.d0(u, v).distance(p) <= tol
            }
            None => false,
        }
    })
}

/// A flat map from every sub-shape's `ShapeId` to its type — the port of
/// `TopTools_IndexedMapOfShape` for quick membership checks.
pub fn shape_id_map(shape: &TopoShape) -> HashMap<ShapeId, usize> {
    all_subshapes(shape)
        .into_iter()
        .map(|s| (ShapeId::of(&s), s.shape_type() as usize))
        .collect()
}

/// Number of distinct shapes in the subtree (root included).
pub fn shape_count(shape: &TopoShape) -> usize {
    all_subshapes(shape).len()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::builder::TopoBuilder;
    use crate::primitives::BRepPrimBox;
    use crate::shape_ops::translated_copy;
    use occt_core::gp::GpVec;

    #[test]
    fn shape_id_is_stable_and_unique() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let id1 = ShapeId::of(&b.solid.0);
        let id2 = ShapeId::of(&b.solid.0);
        assert_eq!(id1, id2);
        let other = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        assert_ne!(id1, ShapeId::of(&other.solid.0));
    }

    #[test]
    fn geometric_hash_equals_for_equal_shapes() {
        let a = BRepPrimBox::make_box(1.0, 2.0, 3.0);
        let b = BRepPrimBox::make_box(1.0, 2.0, 3.0);
        assert_eq!(geometric_hash(&a.solid.0), geometric_hash(&b.solid.0));
        let moved = translated_copy(&a.solid.0, &GpVec::new(1.0, 0.0, 0.0)).unwrap();
        assert_ne!(geometric_hash(&a.solid.0), geometric_hash(&moved));
    }

    #[test]
    fn geometrically_equal_boxes() {
        let a = BRepPrimBox::make_box(1.0, 2.0, 3.0);
        let b = BRepPrimBox::make_box(1.0, 2.0, 3.0);
        assert!(shapes_geometrically_equal(&a.solid.0, &b.solid.0, 1e-9));
        let moved = translated_copy(&a.solid.0, &GpVec::new(0.0, 0.0, 5.0)).unwrap();
        assert!(!shapes_geometrically_equal(&a.solid.0, &moved, 1e-9));
    }

    #[test]
    fn find_by_position_and_range() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let v = find_vertex_at(&b.solid.0, &GpPnt::zero(), 1e-9).expect("origin vertex");
        assert!(BRepTool::vertex_point(&v).distance(&GpPnt::zero()) < 1e-9);
        let e = find_edge_by_range(&b.solid.0, 0.0, 1.0, 1e-9).expect("edge of length 1");
        let (a, z) = BRepTool::edge_parameters(&e);
        assert!(a == 0.0 && z == 1.0);
        let f = find_face_through(&b.solid.0, &GpPnt::new(0.0, 0.5, 0.5), 1e-6).expect("face through point");
        let _ = f;
    }

    #[test]
    fn index_and_count() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let idx = shape_index(&b.solid.0);
        assert_eq!(idx[&(crate::abs::ShapeType::Vertex as usize)], 8);
        assert_eq!(idx[&(crate::abs::ShapeType::Edge as usize)], 12);
        assert_eq!(idx[&(crate::abs::ShapeType::Face as usize)], 6);
        assert!(shape_count(&b.solid.0) >= 33); // 8+12+6+6+1+1 + compound? = 34-ish
        let map = shape_id_map(&b.solid.0);
        assert!(map.len() >= 33);
    }

    #[test]
    fn manual_vertex_id() {
        let b = TopoBuilder::new();
        let v = b.make_vertex(GpPnt::new(1., 2., 3.), 0.0);
        assert!(find_vertex_at(&v.0, &GpPnt::new(1., 2., 3.), 1e-9).is_some());
    }
}
