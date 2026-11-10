use std::sync::Arc;


use crate::abs::ShapeType;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::topexp::Explorer;


pub(super) fn collect_vertices(shape: &TopoShape) -> Vec<Vertex> {
    let mut out = Vec::new();
    let mut seen = std::collections::HashSet::new();
    let mut ex = Explorer::new(shape, ShapeType::Vertex);
    while ex.more() {
        let s = ex.current().clone();
        ex.next();
        if seen.insert(Arc::as_ptr(&s.tshape) as usize) {
            if let Some(v) = Vertex::wrap(s) {
                out.push(v);
            }
        }
    }
    out
}

pub(super) fn collect_edges(shape: &TopoShape) -> Vec<Edge> {
    let mut out = Vec::new();
    let mut seen = std::collections::HashSet::new();
    let mut ex = Explorer::new(shape, ShapeType::Edge);
    while ex.more() {
        let s = ex.current().clone();
        ex.next();
        if seen.insert(Arc::as_ptr(&s.tshape) as usize) {
            if let Some(e) = Edge::wrap(s) {
                out.push(e);
            }
        }
    }
    out
}

pub(super) fn collect_faces(shape: &TopoShape) -> Vec<Face> {
    let mut out = Vec::new();
    let mut seen = std::collections::HashSet::new();
    let mut ex = Explorer::new(shape, ShapeType::Face);
    while ex.more() {
        let s = ex.current().clone();
        ex.next();
        if seen.insert(Arc::as_ptr(&s.tshape) as usize) {
            if let Some(f) = Face::wrap(s) {
                out.push(f);
            }
        }
    }
    out
}
