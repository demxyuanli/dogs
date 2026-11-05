//! The vertex helpers of `BRepFill_OffsetWire`: `VertexFromNode`
//! (`cxx:2310-2345`), `StoreInMap` (`cxx:2348-2376`), `DoubleOrNotInside`
//! (`cxx:2604-2628`), `IsSmallClosedEdge` (`cxx:2631-2657`),
//! `CheckSmallParamOnEdge` (`cxx:2922-2940`) and `IsInnerEdge`
//! (`cxx:2575-2601`).

use std::collections::HashMap;

use occt_core::gp::{GpPnt, GpPnt2d};
use occt_core::precision::{CONFUSION, PCONFUSION};

use crate::abs::ShapeType;
use crate::brep_fill_offset_wire::maps::ShapeToShapeMap;

use crate::builder::TopoBuilder;
use crate::mat::NodeId;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{is_same, map_shapes_and_ancestors};

/// `VertexFromNode(aNode, Offset, PN, MapNodeVertex, VN)` (`cxx:2310-2345`):
/// `Some` when the node sits exactly `Offset` away, i.e. it is a vertex of the
/// offset. The node's vertex is shared through `map_node_vertex`.
pub(crate) fn vertex_from_node(
    node: NodeId,
    node_distance: f64,
    node_infinite: bool,
    offset: f64,
    pn: &GpPnt2d,
    map_node_vertex: &mut HashMap<NodeId, TopoShape>,
) -> Option<Vertex> {
    // `cxx:2320`: `!Infinite() && |Distance() - Offset| < Confusion`.
    if node_infinite || (node_distance - offset).abs() >= CONFUSION {
        return None;
    }
    // `cxx:2326-2336`.
    if let Some(v) = map_node_vertex.get(&node) {
        return Some(Vertex(v.clone()));
    }
    let builder = TopoBuilder::new();
    let p = GpPnt::new(pn.x(), pn.y(), 0.0);
    let v = builder.make_vertex(p, CONFUSION);
    map_node_vertex.insert(node, v.0.clone());
    Some(v)
}

/// `StoreInMap(V1, V2, MapVV)` (`cxx:2348-2376`): make `V1` point at the same
/// offset vertex as `V2`, and repoint every value that already was `V1`.
pub(crate) fn store_in_map(v1: &TopoShape, v2: &TopoShape, map_vv: &mut ShapeToShapeMap) {
    let mut new_v = v2.clone();
    if let Some(found) = map_vv.find_from_key(v2) {
        new_v = found.clone();
    }
    if let Some(slot) = map_vv.change_from_key(v1) {
        *slot = new_v.clone();
    }
    for i in 1..=map_vv.extent() {
        if let Some(slot) = map_vv.value_mut(i) {
            if is_same(slot, v1) {
                *slot = new_v.clone();
            }
        }
    }
    map_vv.add(v1, new_v);
}

/// `DoubleOrNotInside(LV, V)` (`cxx:2604-2628`): `true` when `V` appears twice
/// in `detromp`'s list or does not appear at all.
pub(crate) fn double_or_not_inside(detromp: &[TopoShape], v: &Vertex) -> bool {
    let mut seen = false;
    for s in detromp {
        if is_same(s, &v.0) {
            if seen {
                return true;
            }
            seen = true;
        }
    }
    !seen
}

/// `IsSmallClosedEdge(anEdge, aVertex)` (`cxx:2631-2657`): the p-curve's first,
/// last and middle points all sit on the vertex.
///
/// Deviation: OCCT reads the first `BRep_CurveRepresentation` of the edge, which
/// may be the 3D curve; the port reads the first p-curve representation, the one
/// that exists for every edge this class builds (`BRepLib_MakeEdge(G2d, Plane)`).
pub(crate) fn is_small_closed_edge(edge: &Edge, v: &Vertex) -> bool {
    let reg = GeometryRegistry::global();
    let reps = reg.edge_pcurve_reps(&edge.0);
    let Some((_, pc, f, l)) = reps.first() else {
        return false;
    };
    let pv = reg.vertex_point(&v.0);
    let pv2d = GpPnt2d::new(pv.x(), pv.y());
    let p_first = pc.d0(*f);
    let p_last = pc.d0(*l);
    let p_mid = pc.d0((f + l) * 0.5);
    let tol = reg.vertex_tolerance(&v.0) * 1.5;
    p_first.distance(&pv2d) <= tol
        && p_last.distance(&pv2d) <= tol
        && p_mid.distance(&pv2d) <= tol
}

/// `CheckSmallParamOnEdge(anEdge)` (`cxx:2922-2940`): `false` for a null-length
/// range.
pub(crate) fn check_small_param_on_edge(edge: &Edge) -> bool {
    let (f, l) = GeometryRegistry::global().edge_parameters(&edge.0);
    (l - f).abs() >= PCONFUSION
}

/// `IsInnerEdge(ProE, AllSpine, TrPar1, TrPar2)` (`cxx:2575-2601`): `Some(range)`
/// when `pro_e` is an edge that shares both its vertices with other edges of the
/// spine, i.e. an interior edge whose ends must not become offset end vertices.
pub(crate) fn is_inner_edge(pro_e: &TopoShape, all_spine: &Face) -> Option<(f64, f64)> {
    // `cxx:2580-2584`.
    if pro_e.shape_type() != ShapeType::Edge {
        return None;
    }
    let an_edge = Edge(pro_e.clone());

    // `cxx:2587-2598`.
    let vemap = map_shapes_and_ancestors(&all_spine.0, ShapeType::Vertex, ShapeType::Edge);
    for list in vemap.values() {
        if list.len() == 1 && is_same(&an_edge.0, &list[0]) {
            return None;
        }
    }

    // `cxx:2600-2601`.
    Some(GeometryRegistry::global().edge_parameters(&an_edge.0))
}

/// `BRep_Builder::UpdateVertex(V, Tol)` (`BRep_Builder.cxx:1113-1125`): the
/// stored tolerance becomes the larger of `tol` and the current one.
pub(crate) fn update_vertex_tolerance(v: &Vertex, tol: f64) {
    let reg = GeometryRegistry::global();
    let Some(mut geom) = reg.vertex_geom(&v.0) else {
        return;
    };
    geom.tolerance = geom.tolerance.max(tol);
    reg.set_vertex(&v.0, geom);
}

/// `BRep_Builder::UpdateVertex(V, P, Tol)` (`BRep_Builder.cxx:1096-1111`): the
/// point is replaced and the tolerance raised.
pub(crate) fn update_vertex_point(v: &Vertex, p: &GpPnt, tol: f64) {
    let reg = GeometryRegistry::global();
    let current = reg.vertex_tolerance(&v.0);
    reg.set_vertex(
        &v.0,
        crate::tgeometry::VertexGeom {
            point: *p,
            tolerance: current.max(tol),
        },
    );
}

/// `BRep_TVertex::Pnt(P)` (`BRep_TVertex.hxx:63-66`): the point is replaced and
/// the tolerance is left untouched — the direct member write
/// `PerformWithBiLo` uses at `cxx:1159` / `:1167`, as opposed to
/// `BRep_Builder::UpdateVertex`.
pub(crate) fn set_vertex_point(v: &Vertex, p: &GpPnt) {
    let reg = GeometryRegistry::global();
    let tolerance = reg.vertex_tolerance(&v.0);
    reg.set_vertex(&v.0, crate::tgeometry::VertexGeom { point: *p, tolerance });
}
