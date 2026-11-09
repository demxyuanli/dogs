//! `LocOpe::Closed` (`LocOpe.cxx:43-126`).
//!
//! Decides whether a wire closes up on the face: its first and last vertex must
//! be the same, and the wire's first edge must start at the same 2D point the
//! last edge ends at (within `Precision::PConfusion(Precision::Confusion())`),
//! which is what distinguishes a genuine closed loop from a wire whose ends
//! merely touch at a vertex.

use occt_core::precision::CONFUSION;

use crate::abs::{Orientation, ShapeType};
use crate::boptools_2d;
use crate::shape::{Edge, Face, TopoShape, Wire};
use crate::topo_tools::{explorer_shapes, explorer_shapes_avoiding, wire_vertices};

/// `LocOpe::Closed(const TopoDS_Wire& W, const TopoDS_Face& F)`
/// (`LocOpe.cxx:43-115`).
pub fn closed_wire(w: &Wire, f: &Face) -> bool {
    // `cxx:45-50`: `TopExp::Vertices(W, Vf, Vl)` with the default `CumOri`.
    let (vf, vl) = wire_vertices(w);
    let same = match (&vf, &vl) {
        (None, None) => true,
        (Some(a), Some(b)) => crate::topo_tools::is_same(&a.0, &b.0),
        _ => false,
    };
    if !same {
        return false;
    }
    let Some(vf) = vf else {
        // Two null vertices are `IsSame` in OCCT, then `TopoDS::Edge(exp.Current())`
        // on an exhausted explorer (`cxx:69`) is what the reference does - an
        // empty wire has no first edge. Answer `false` rather than fabricate one.
        return false;
    };
    let Some(vl) = vl else {
        return false;
    };

    // `cxx:54-69`: the edge whose child vertex at `Vf` is FORWARD, i.e. the edge
    // `Vf` starts. `TopExp_Explorer(W.Oriented(FORWARD), TopAbs_EDGE)` then
    // `TopExp_Explorer(edge, TopAbs_VERTEX)`; both compose through
    // `TopoDS_Iterator`, so the child orientations below are the composed ones.
    let edges = explorer_shapes_avoiding_forward(&w.0);
    let Some(ef) = find_edge_at_vertex(&edges, &vf, Orientation::Forward) else {
        return false;
    };
    // `cxx:71-87`: the edge whose child vertex at `Vl` is REVERSED.
    let Some(el) = find_edge_at_vertex(&edges, &vl, Orientation::Reversed) else {
        return false;
    };

    // `cxx:89-108`: the 2D endpoints in the edge direction.
    let (pf, pl) = match (
        endpoint_2d(&ef, f, ef.0.orientation() == Orientation::Forward),
        endpoint_2d(&el, f, el.0.orientation() != Orientation::Forward),
    ) {
        (Some(pf), Some(pl)) => (pf, pl),
        _ => return false,
    };

    // `cxx:110-114`: `Precision::PConfusion(Precision::Confusion())`, which is
    // `Parametric(P, T) = P / T` with both arguments `Confusion()`.
    let tol2d = CONFUSION / CONFUSION;
    pf.distance(&pl) <= tol2d
}

/// `LocOpe::Closed(const TopoDS_Edge& E, const TopoDS_Face& F)`
/// (`LocOpe.cxx:119-126`): the same test on a one-edge wire.
pub fn closed_edge(e: &Edge, f: &Face) -> bool {
    let mut w = crate::shape::Wire::new();
    let mut child = e.0.clone();
    child.set_orientation(Orientation::Forward);
    crate::builder::TopoBuilder::new().add(&mut w.0, &child);
    closed_wire(&w, f)
}

/// `TopExp_Explorer(W.Oriented(TopAbs_FORWARD), TopAbs_EDGE)`: the wire oriented
/// forward, so the stored child orientations reach the explorer unchanged.
fn explorer_shapes_avoiding_forward(w: &TopoShape) -> Vec<Edge> {
    let mut forward = w.clone();
    forward.set_orientation(Orientation::Forward);
    explorer_shapes(&forward, ShapeType::Edge)
        .into_iter()
        .filter_map(Edge::wrap)
        .collect()
}

/// The first edge of `edges` whose child vertex `v` carries `ori`.
fn find_edge_at_vertex(edges: &[Edge], v: &crate::shape::Vertex, ori: Orientation) -> Option<Edge> {
    for e in edges {
        for child in explorer_shapes_avoiding(&e.0, ShapeType::Vertex, ShapeType::Edge) {
            if child.orientation() == ori && crate::topo_tools::is_same(&child, &v.0) {
                return Some(e.clone());
            }
        }
    }
    None
}

/// The p-curve endpoint the reference reads: `C2d->Value(f)` when `take_first`.
fn endpoint_2d(e: &Edge, f: &Face, take_first: bool) -> Option<occt_core::gp::GpPnt2d> {
    let (c2d, first, last) = boptools_2d::curve_on_surface_range(e, f)?;
    Some(c2d.d0(if take_first { first } else { last }))
}
