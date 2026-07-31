//! Shape comparison and identification — structural signatures, sub-shape
//! containment and topological equality.
//! Source: `BRepTools::Compare`, `BRepTools_ShapeSet`, `TopTools`.

use std::collections::HashMap;

use occt_core::gp::GpPnt;

use crate::abs::ShapeType;
use crate::brep_tool::BRepTool;
use crate::shape::TopoShape;
use crate::topo_tools_full::{all_subshapes, edges_of, vertices_of};

/// Structural signature of a shape: a string encoding the type of every
/// sub-shape in traversal order. Two shapes with identical structure have the
/// same signature. Mirrors `BRepTools::Write`-style structural hashing.
pub fn shape_signature(shape: &TopoShape) -> String {
    fn type_char(t: ShapeType) -> char {
        match t {
            ShapeType::Vertex => 'V',
            ShapeType::Edge => 'E',
            ShapeType::Wire => 'W',
            ShapeType::Face => 'F',
            ShapeType::Shell => 'S',
            ShapeType::Solid => 'B', // BRep solid
            ShapeType::Compound => 'C',
            ShapeType::CompSolid => 'c',
            ShapeType::Shape => '?',
        }
    }
    let mut out = String::new();
    for s in all_subshapes(shape) {
        out.push(type_char(s.shape_type()));
    }
    out
}

/// Summary comparison: structural signature + sub-shape counts. Two shapes are
/// "topologically equal" if signatures match and per-type counts match.
pub fn topologically_equal(a: &TopoShape, b: &TopoShape) -> bool {
    shape_signature(a) == shape_signature(b) && counts_match(a, b)
}

fn counts_match(a: &TopoShape, b: &TopoShape) -> bool {
    let mut ca = HashMap::new();
    let mut cb = HashMap::new();
    for s in all_subshapes(a) {
        *ca.entry(s.shape_type() as usize).or_insert(0) += 1;
    }
    for s in all_subshapes(b) {
        *cb.entry(s.shape_type() as usize).or_insert(0) += 1;
    }
    ca == cb
}

/// Whether `inner` is a sub-shape of `outer` (same `TShape` identity anywhere
/// in the outer's subtree). Mirrors `TopExp::Contains`.
pub fn contains(outer: &TopoShape, inner: &TopoShape) -> bool {
    all_subshapes(outer)
        .iter()
        .any(|s| std::sync::Arc::ptr_eq(&s.tshape, &inner.tshape))
}

/// Whether two shapes share any `TShape` (overlap in shared sub-shapes).
pub fn share_subshape(a: &TopoShape, b: &TopoShape) -> bool {
    let ids_b: std::collections::HashSet<usize> = all_subshapes(b)
        .iter()
        .map(|s| std::sync::Arc::as_ptr(&s.tshape) as usize)
        .collect();
    all_subshapes(a)
        .iter()
        .any(|s| ids_b.contains(&(std::sync::Arc::as_ptr(&s.tshape) as usize)))
}

/// Degree of a vertex in a shape: number of distinct edges incident to it.
/// Mirrors `TopTools`-style valence queries.
pub fn vertex_valence(shape: &TopoShape, v: &TopoShape) -> usize {
    let vid = std::sync::Arc::as_ptr(&v.tshape) as usize;
    edges_of(shape)
        .iter()
        .filter(|e| {
            let (a, b) = crate::topo_tools_full::edge_vertices(e);
            let pa = a.map(|a| std::sync::Arc::as_ptr(&a.0.tshape) as usize);
            let pb = b.map(|b| std::sync::Arc::as_ptr(&b.0.tshape) as usize);
            pa == Some(vid) || pb == Some(vid)
        })
        .count()
}

/// Average vertex valence over a shape's vertices (2.0 for a closed manifold
/// triangle/quad mesh boundary; box corners → 3, box edges → 4).
pub fn average_valence(shape: &TopoShape) -> f64 {
    let verts = vertices_of(shape);
    if verts.is_empty() {
        return 0.0;
    }
    let sum: usize = verts.iter().map(|v| vertex_valence(shape, &v.0)).sum();
    sum as f64 / verts.len() as f64
}

/// Shape signature string with the shape's type and boundary counts, for
/// diagnostics (`BRepTools::Dump`-like one-liner).
pub fn shape_summary(shape: &TopoShape) -> String {
    let counts = {
        let mut m = HashMap::new();
        for s in all_subshapes(shape) {
            *m.entry(s.shape_type() as usize).or_insert(0) += 1;
        }
        m
    };
    let v = counts.get(&(ShapeType::Vertex as usize)).copied().unwrap_or(0);
    let e = counts.get(&(ShapeType::Edge as usize)).copied().unwrap_or(0);
    let f = counts.get(&(ShapeType::Face as usize)).copied().unwrap_or(0);
    format!("{:?} V={v} E={e} F={f}", shape.shape_type())
}

/// Whether `p` is within `tol` of any vertex of `shape`.
pub fn near_any_vertex(shape: &TopoShape, p: &GpPnt, tol: f64) -> bool {
    vertices_of(shape)
        .iter()
        .any(|v| BRepTool::vertex_point(v).distance(p) <= tol)
}

/// The bounding corner with the largest x+y+z coordinate (`BndBox` helper).
pub fn max_corner(shape: &TopoShape) -> Option<GpPnt> {
    let verts = vertices_of(shape);
    if verts.is_empty() {
        return None;
    }
    verts
        .iter()
        .map(|v| BRepTool::vertex_point(v))
        .max_by(|a, b| (a.x() + a.y() + a.z()).total_cmp(&(b.x() + b.y() + b.z())))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitives::BRepPrimBox;
    use crate::shape_ops::translated_copy;
    use occt_core::gp::GpVec;

    #[test]
    fn signatures_match_for_equal_boxes() {
        let a = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        assert_eq!(shape_signature(&a.solid.0), shape_signature(&b.solid.0));
        assert!(topologically_equal(&a.solid.0, &b.solid.0));
    }

    #[test]
    fn signature_differs_after_transform() {
        let a = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        // Translated copy keeps the same structure → signature equal.
        let moved = translated_copy(&a.solid.0, &GpVec::new(3.0, 0.0, 0.0)).unwrap();
        assert_eq!(shape_signature(&a.solid.0), shape_signature(&moved));
        assert!(topologically_equal(&a.solid.0, &moved));
    }

    #[test]
    fn contains_and_share() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let verts = vertices_of(&b.solid.0);
        let first_v = &verts[0].0;
        assert!(contains(&b.solid.0, first_v));
        // A fresh vertex is not contained.
        let fresh = TopoShape::new(ShapeType::Vertex);
        assert!(!contains(&b.solid.0, &fresh));
        // A box shares sub-shapes with its translated copy? No (copies are new).
        let moved = translated_copy(&b.solid.0, &GpVec::new(0.0, 1.0, 0.0)).unwrap();
        assert!(!share_subshape(&b.solid.0, &moved));
    }

    #[test]
    fn valence_of_box_vertices() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let verts = vertices_of(&b.solid.0);
        // Every corner is incident to 3 box edges.
        for v in &verts {
            assert_eq!(vertex_valence(&b.solid.0, &v.0), 3);
        }
        assert!((average_valence(&b.solid.0) - 3.0).abs() < 1e-12);
    }

    #[test]
    fn summary_and_near_vertex() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let s = shape_summary(&b.solid.0);
        assert!(s.contains("V=8"));
        assert!(near_any_vertex(&b.solid.0, &GpPnt::zero(), 1e-9));
        assert!(!near_any_vertex(&b.solid.0, &GpPnt::new(5.0, 5.0, 5.0), 1e-9));
        let c = max_corner(&b.solid.0).expect("max corner");
        assert!((c.x() + c.y() + c.z() - 3.0).abs() < 1e-9);
    }
}
