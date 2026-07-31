//! Geometric access to topology — a port of `BRep_Tool`.
//!
//! Real geometry (curves, surfaces, tolerances) is not yet stored on shapes in
//! this port, so the accessors return documented placeholders while keeping
//! the `BRep_Tool` call shapes. `BRepTool` is a unit struct used purely as a
//! namespace for static-like queries.

use crate::abs::ShapeType;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use occt_core::gp::{GpPnt, GpVec};

/// Namespace for BRep_Tool-style static queries.
pub struct BRepTool;

impl BRepTool {
    /// The location of a vertex. Placeholder: returns the origin.
    pub fn vertex_point(v: &Vertex) -> GpPnt {
        let _ = v;
        GpPnt::zero()
    }

    /// Vertex tolerance. Placeholder.
    pub fn vertex_tolerance(v: &Vertex) -> f64 {
        let _ = v;
        0.0
    }

    /// The curve underlying an edge. Placeholder: would return Handle(Geom_Curve).
    pub fn edge_curve(e: &Edge) -> Option<GpPnt> {
        let _ = e;
        None
    }

    /// Parameter range of an edge. Placeholder: unbounded.
    pub fn edge_parameters(e: &Edge) -> (f64, f64) {
        let _ = e;
        (f64::NEG_INFINITY, f64::INFINITY)
    }

    /// Edge tolerance. Placeholder.
    pub fn edge_tolerance(e: &Edge) -> f64 {
        let _ = e;
        0.0
    }

    /// The surface underlying a face. Placeholder: would return Handle(Geom_Surface).
    pub fn face_surface(f: &Face) -> Option<GpPnt> {
        let _ = f;
        None
    }

    /// Face tolerance. Placeholder.
    pub fn face_tolerance(f: &Face) -> f64 {
        let _ = f;
        0.0
    }

    /// Endpoint vertices of an edge as curve points. Placeholder.
    pub fn edge_vertices(e: &Edge) -> Option<(GpPnt, GpPnt)> {
        let _ = e;
        None
    }

    /// Whether the edge is same-parameter. Placeholder.
    pub fn same_parameter(e: &Edge) -> bool {
        let _ = e;
        false
    }

    /// Whether the edge is degenerated. Placeholder.
    pub fn is_degenerated(e: &Edge) -> bool {
        let _ = e;
        false
    }

    /// The first (`i == 0`) or last parameter of the edge.
    pub fn parameter_on_edge(e: &Edge, i: usize) -> f64 {
        let (first, last) = Self::edge_parameters(e);
        if i == 0 {
            first
        } else {
            last
        }
    }

    /// A closed edge has its first parameter equal to its last within tolerance.
    pub fn is_closed_edge(e: &Edge) -> bool {
        let (first, last) = Self::edge_parameters(e);
        (first - last).abs() <= Self::edge_tolerance(e)
    }

    /// Structural integrity of a shape list. Always true in this port:
    /// TShape types are fixed at construction.
    pub fn check_shape_integrity(shapes: &[TopoShape]) -> bool {
        let _ = shapes;
        true
    }

    /// Human-readable type name.
    pub fn shape_type_str(t: ShapeType) -> &'static str {
        t.to_str()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn edge() -> Edge {
        Edge::from(TopoShape::new(ShapeType::Edge))
    }

    #[test]
    fn edge_parameters_are_unbounded() {
        let (first, last) = BRepTool::edge_parameters(&edge());
        assert_eq!(first, f64::NEG_INFINITY);
        assert_eq!(last, f64::INFINITY);
    }

    #[test]
    fn unbounded_edge_is_not_closed() {
        assert!(!BRepTool::is_closed_edge(&edge()));
    }

    #[test]
    fn shape_type_str_names_vertex() {
        assert_eq!(BRepTool::shape_type_str(ShapeType::Vertex), "Vertex");
    }

    #[test]
    fn empty_shape_list_is_integral() {
        let shapes: Vec<TopoShape> = Vec::new();
        assert!(BRepTool::check_shape_integrity(&shapes));
    }
}
