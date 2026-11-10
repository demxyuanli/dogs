//! `TopExp::Vertices` / face-cast helpers local to the candidate search.


use crate::abs::ShapeType;
use crate::shape::{Edge, Face, TopoShape};



/// `TopExp::Vertices(edge)` helper kept local so the module does not depend on
/// the `Shape` wrappers of the caller.
#[allow(dead_code)]
fn edge_end_points(edge: &TopoShape) -> Option<(TopoShape, TopoShape)> {
    match crate::topo_tools::edge_vertices(&Edge(edge.clone())) {
        (Some(a), Some(b)) => Some((a.0, b.0)),
        _ => None,
    }
}

#[allow(dead_code)]
fn face_of(shape: &TopoShape) -> Option<Face> {
    if shape.shape_type() != ShapeType::Face {
        return None;
    }
    Face::wrap(shape.clone())
}
