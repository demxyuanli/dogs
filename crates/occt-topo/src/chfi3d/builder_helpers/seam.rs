//! Seam detection: `ChFi3d_Couture`, `ChFi3d_CoutureOnVertex` and
//! `ChFi3d_IsPseudoSeam`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_0.cxx:5750-5832`.

use crate::brep_tool::BRepTool;
use crate::brep_tools::is_really_closed;
use crate::shape::{Edge, Face, Vertex};
use crate::topo_tools::{edge_vertices, explorer_shapes, is_same};
use crate::abs::ShapeType;

/// `ChFi3d_Couture(F, couture, edgecouture)` (`cxx:5750-5769`): the first seam
/// edge of the face.
///
/// OCCT asks `BRep_Tool::IsClosed(Ecur, Surf, Loc)` with the face's own surface
/// and location, which is [`BRepTool::is_closed_edge_face`] here.
pub fn couture(f: &Face) -> Option<Edge> {
    for value in explorer_shapes(&f.0, ShapeType::Edge) {
        let Some(cur) = Edge::wrap(value) else {
            continue;
        };
        if BRepTool::is_closed_edge_face(&cur, f) {
            return Some(cur);
        }
    }
    None
}

/// `ChFi3d_CoutureOnVertex(F, V, couture, edgecouture)` (`cxx:5773-5801`): the
/// first seam edge of the face that carries the vertex `v` at either end.
pub fn couture_on_vertex(f: &Face, v: &Vertex) -> Option<Edge> {
    for value in explorer_shapes(&f.0, ShapeType::Edge) {
        let Some(cur) = Edge::wrap(value) else {
            continue;
        };
        if !BRepTool::is_closed_edge_face(&cur, f) {
            continue;
        }
        let (vf, vl) = edge_vertices(&cur);
        let carries = |e: Option<Vertex>| e.is_some_and(|e| is_same(&e.0, &v.0));
        if carries(vf) || carries(vl) {
            return Some(cur);
        }
    }
    None
}

/// `ChFi3d_IsPseudoSeam(E, F)` (`cxx:5805-5832`): true when the seam edge `E`
/// shares a vertex with another really-closed edge of `F`.
///
/// The shared vertex makes the two seams a single pseudoseam pair: `E` is a
/// regular seam of the face rather than one of the two sides of a degenerate
/// vertex.
pub fn is_pseudo_seam(e: &Edge, f: &Face) -> bool {
    if !BRepTool::is_closed_edge_face(e, f) {
        return false;
    }
    let (vf, vl) = edge_vertices(e);
    for value in explorer_shapes(&f.0, ShapeType::Edge) {
        let Some(cur) = Edge::wrap(value) else {
            continue;
        };
        if is_same(&cur.0, &e.0) {
            continue;
        }
        let (v1, v2) = edge_vertices(&cur);
        let touches = |a: Option<Vertex>, b: Option<Vertex>| {
            [a, b].into_iter().flatten().any(|x| {
                [vf.clone(), vl.clone()]
                    .into_iter()
                    .flatten()
                    .any(|y| is_same(&x.0, &y.0))
            })
        };
        if touches(v1, v2) && is_really_closed(&cur, f) {
            return true;
        }
    }
    false
}
