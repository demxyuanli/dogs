//! `BRepFill_Filling::BuildWires` / `FindExtremitiesOfHoles` and the two file
//! static helpers of `BRepFill_Filling.cxx`.
//!
//! Source: `TKBool/BRepFill/BRepFill_Filling.cxx:64-140` (`MakeFinVec`,
//! `WireFromList`), `:392-479` (`BuildWires`), `:485-567`
//! (`FindExtremitiesOfHoles`) (OCCT 8.0.0).
//!
//! Translation notes:
//! - `NCollection_List<TopoDS_Shape>` becomes `Vec<TopoShape>`; `Remove(itl)`
//!   and `RemoveFirst()` are index removals, so the "first matching edge" walk
//!   is a plain forward loop with `remove(i)`.
//! - `TopExp::Vertices(E, V1, V2)` (`TopExp.cxx:182-233`, `CumOri = false`) is
//!   [`crate::topo_tools::edge_vertices`]; the `CumOri = true` form used by
//!   `WireFromList` is [`oriented_edge_vertices`], which swaps the stored pair
//!   when the edge itself is `REVERSED` — exactly what `TopoDS_Iterator(E, true)`
//!   reports as the cumulated orientation.
//! - A null `TopoDS_Vertex` is `None`; `IsSame` between a null and a non-null
//!   vertex is false in OCCT, so a missing vertex never matches.
//! - `BRepLib_MakeWire::Add` failures are discarded, as the reference discards
//!   `MW.Error()`.

use occt_core::gp::GpVec;

use crate::abs::Orientation;
use crate::brep_lib_make_wire::MakeWire;
use crate::brep_tool::BRepTool;
use crate::brep_tools_wire_explorer::WireExplorer;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, TopoShape, Vertex, Wire};
use crate::topo_tools::{edge_vertices, is_same, shape_map_key, wire_vertices};

/// The `myOldNewMap` of `BRepFill_Filling` (`hxx:216`):
/// `NCollection_DataMap<TopoDS_Shape, TopoDS_Shape, TopTools_ShapeMapHasher>`.
pub type OldNewMap = std::collections::HashMap<crate::topo_tools::ShapeMapKey, TopoShape>;

/// `TopExp::Vertices(E, V1, V2, /*CumOri=*/true)` (`TopExp.cxx:182-233`): the
/// first child reads `FORWARD` and the second `REVERSED` after the edge's own
/// orientation is applied, so a `REVERSED` edge yields the stored pair swapped.
pub fn oriented_edge_vertices(edge: &Edge) -> (Option<Vertex>, Option<Vertex>) {
    let (first, last) = edge_vertices(edge);
    if edge.0.orientation() == Orientation::Reversed {
        (last, first)
    } else {
        (first, last)
    }
}

/// `static gp_Vec MakeFinVec(aWire, aVertex)` (`cxx:64-83`): the vector from the
/// far end of `aWire` to `aVertex`.
///
/// `None` is the reference's null `Origin` — `BRep_Tool::Pnt` on a null vertex,
/// a crash in OCCT. It cannot happen for the call sites of
/// [`find_extremities_of_holes`], which always pass a vertex of the wire.
fn make_fin_vec(a_wire: &Wire, a_vertex: &Vertex) -> Option<GpVec> {
    let mut origin: Option<Vertex> = None;
    let mut explo = WireExplorer::of_wire(a_wire);
    while explo.more() {
        if let Some(current) = explo.current_edge() {
            let (vfirst, vlast) = edge_vertices(&current);
            if let Some(v) = vfirst.as_ref() {
                if is_same(&v.0, &a_vertex.0) {
                    origin = vlast;
                    break;
                }
            }
            if let Some(v) = vlast.as_ref() {
                if is_same(&v.0, &a_vertex.0) {
                    origin = vfirst;
                    break;
                }
            }
        }
        explo.next();
    }
    origin.map(|o| {
        GpVec::from_pnts(
            &BRepTool::vertex_point_world(&o),
            &BRepTool::vertex_point_world(a_vertex),
        )
    })
}

/// Both handles name the same `TShape` at the same location
/// (`TopoDS_Shape::IsSame`).
fn same_vertex(a: &Option<Vertex>, b: &Option<Vertex>) -> bool {
    match (a, b) {
        (Some(x), Some(y)) => is_same(&x.0, &y.0),
        _ => false,
    }
}

/// `static TopoDS_Wire WireFromList(NCollection_List<TopoDS_Shape>& Edges)`
/// (`cxx:85-140`): chain the edges into one wire, reversing each edge that only
/// shares its last vertex with the current chain.
///
/// The `Message::SendWarning` arm (`cxx:128-133`) is the `break`: the wire is
/// returned incomplete. `Edges.First()` on an empty list raises in OCCT and is
/// unreachable from `Build`, which passes the edges it just built.
pub fn wire_from_list(edges: &[TopoShape]) -> Wire {
    let b = TopoBuilder::new();
    let mut a_wire = Wire::new();
    let mut edges: Vec<TopoShape> = edges.to_vec();
    if edges.is_empty() {
        a_wire.0.set_closed(true);
        return a_wire;
    }
    let mut an_edge = edges.remove(0);
    b.add(&mut a_wire.0, &an_edge);
    let (mut v1, mut v2) = oriented_edge_vertices(&Edge(an_edge.clone()));

    while !edges.is_empty() {
        let mut found: Option<usize> = None;
        for (i, candidate) in edges.iter().enumerate() {
            let mut edge = Edge(candidate.clone());
            let (v3, v4) = oriented_edge_vertices(&edge);
            if same_vertex(&v1, &v3)
                || same_vertex(&v1, &v4)
                || same_vertex(&v2, &v3)
                || same_vertex(&v2, &v4)
            {
                if same_vertex(&v1, &v3) {
                    edge.0.reverse();
                    v1 = v4;
                } else if same_vertex(&v1, &v4) {
                    v1 = v3;
                } else if same_vertex(&v2, &v3) {
                    v2 = v4;
                } else {
                    edge.0.reverse();
                    v2 = v3;
                }
                an_edge = edge.0;
                found = Some(i);
                break;
            }
        }
        match found {
            // `if (!itl.More())` (`cxx:128-133`): the walk ran off the list.
            None => break,
            Some(i) => {
                b.add(&mut a_wire.0, &an_edge);
                edges.remove(i);
            }
        }
    }

    a_wire.0.set_closed(true);
    a_wire
}

/// `void BRepFill_Filling::BuildWires(EdgeList, WireList)` (`cxx:392-479`).
///
/// `edge_list` is consumed; each geometric-proximity join records the
/// old-to-new edge mapping `myOldNewMap.Bind(CurEdge.Oriented(FORWARD),
/// NewEdge.Oriented(FORWARD))` (`cxx:457-458`) — the map key ignores the
/// orientation (`TopTools_ShapeMapHasher`), so the port keys on
/// [`shape_map_key`].
pub fn build_wires(edge_list: &mut Vec<TopoShape>, old_new_map: &mut OldNewMap) -> Vec<TopoShape> {
    let mut wire_list: Vec<TopoShape> = Vec::new();
    while !edge_list.is_empty() {
        let mut mw = MakeWire::new();
        let first_edge = Edge(edge_list.remove(0));
        let _ = mw.add(&first_edge);

        loop {
            let cur_wire = mw.wire();
            let (v_wire_0, v_wire_1) = wire_vertices(&cur_wire);
            let v_wire = [v_wire_0, v_wire_1];
            let mut found = false;

            'search: for i in 0..edge_list.len() {
                let cur_edge = Edge(edge_list[i].clone());
                let (v_edge_0, v_edge_1) = edge_vertices(&cur_edge);
                let v_edge = [v_edge_0, v_edge_1];
                for wi in 0..2 {
                    for ei in 0..2 {
                        if same_vertex(&v_wire[wi], &v_edge[ei]) {
                            let _ = mw.add(&cur_edge);
                            edge_list.remove(i);
                            found = true;
                            break 'search;
                        }
                    }
                }
            }

            if !found {
                // `cxx:438-473`: geometric coincidence, gated by *both* vertex
                // tolerances. A null `V_wire[i]` never takes part in OCCT either
                // (`BRep_Tool::Pnt(null)` is a crash there).
                let p_wire = [
                    v_wire[0].as_ref().map(BRepTool::vertex_point_world),
                    v_wire[1].as_ref().map(BRepTool::vertex_point_world),
                ];
                'search2: for i in 0..edge_list.len() {
                    let cur_edge = Edge(edge_list[i].clone());
                    let (v_edge_0, v_edge_1) = edge_vertices(&cur_edge);
                    let v_edge = [v_edge_0, v_edge_1];
                    for wi in 0..2 {
                        for ei in 0..2 {
                            let (Some(pw), Some(vw), Some(ve)) =
                                (&p_wire[wi], &v_wire[wi], &v_edge[ei])
                            else {
                                continue;
                            };
                            let ve_pnt = BRepTool::vertex_point_world(ve);
                            let a_dist = pw.distance(&ve_pnt);
                            if a_dist < BRepTool::vertex_tolerance(vw)
                                && a_dist < BRepTool::vertex_tolerance(ve)
                            {
                                let _ = mw.add(&cur_edge);
                                if let Some(new_edge) = mw.last_edge().cloned() {
                                    old_new_map.insert(
                                        shape_map_key(&cur_edge.0),
                                        new_edge.0.oriented(Orientation::Forward),
                                    );
                                }
                                edge_list.remove(i);
                                found = true;
                                break 'search2;
                            }
                        }
                    }
                }
            }

            if !found {
                // `cxx:475-479`: end of the current wire, begin the next one.
                wire_list.push(mw.wire().0);
                break;
            }
        }
    }
    wire_list
}

/// `void BRepFill_Filling::FindExtremitiesOfHoles(WireList, VerSeq)`
/// (`cxx:485-567`): the ordered extremity pairs of the open wires left over by
/// [`build_wires`], walking the wires by minimal turning angle.
///
/// The reference appends to `VerSeq` without clearing it; the port does the
/// same.
pub fn find_extremities_of_holes(wire_list: &[TopoShape], ver_seq: &mut Vec<TopoShape>) {
    let mut wire_seq: Vec<TopoShape> = wire_list.to_vec();
    if wire_seq.is_empty() {
        // `WireSeq(1)` on an empty sequence raises in OCCT; `Build` always has
        // at least one boundary edge, so this cannot be reached.
        return;
    }
    let the_wire = Wire(wire_seq.remove(0));

    // `BRep_Tool::IsClosed(theWire)` (`cxx:500-503`).
    if BRepTool::is_closed(&the_wire.0) {
        return;
    }

    let (mut vfirst, mut vlast) = wire_vertices(&the_wire);
    if same_vertex(&vfirst, &vlast) {
        // `cxx:508-513`: closed indeed, despite the earlier test.
        return;
    }

    let Some(mut last_vertex) = vlast.clone() else {
        return;
    };
    let Some(mut fin_vec) = make_fin_vec(&the_wire, &last_vertex) else {
        return;
    };
    ver_seq.push(last_vertex.0.clone());

    while !wire_seq.is_empty() {
        let mut min_vtx: Option<Vertex> = None;
        let mut min_ind = 0usize;
        let mut is_last = false;
        let mut min_angle = std::f64::consts::PI;
        let p_vertex = BRepTool::vertex_point_world(&last_vertex);

        for i in 0..wire_seq.len() {
            let cur_wire = Wire(wire_seq[i].clone());
            (vfirst, vlast) = wire_vertices(&cur_wire);
            if let Some(first) = vfirst.as_ref() {
                let dir = GpVec::from_pnts(&p_vertex, &BRepTool::vertex_point_world(first));
                let angle = fin_vec.angle(&dir);
                if angle < min_angle {
                    min_angle = angle;
                    min_vtx = Some(first.clone());
                    min_ind = i;
                    is_last = true;
                }
            }
            if let Some(last) = vlast.as_ref() {
                let dir = GpVec::from_pnts(&p_vertex, &BRepTool::vertex_point_world(last));
                let angle = fin_vec.angle(&dir);
                if angle < min_angle {
                    min_angle = angle;
                    min_vtx = Some(last.clone());
                    min_ind = i;
                    is_last = false;
                }
            }
        }

        // `MinVtx` stays null when no wire was examined; `Append(null)` is a
        // null `TopoDS_Shape` entry, which the port keeps as a bare shape with
        // no `TShape` (it can only happen on a non-convex contour, where the
        // reference then dereferences the null vertex downstream anyway).
        match min_vtx {
            Some(v) => ver_seq.push(v.0),
            None => ver_seq.push(TopoShape::new(crate::abs::ShapeType::Vertex)),
        }

        let min_wire = Wire(wire_seq[min_ind].clone());
        let (min_first, min_last) = wire_vertices(&min_wire);
        if is_last {
            if let Some(l) = min_last {
                if let Some(f) = make_fin_vec(&min_wire, &l) {
                    fin_vec = f;
                }
                last_vertex = l;
            }
        } else if let Some(f) = min_first {
            if let Some(v) = make_fin_vec(&min_wire, &f) {
                fin_vec = v;
            }
            last_vertex = f;
        }
        ver_seq.push(last_vertex.0.clone());
        wire_seq.remove(min_ind);
    }

    let (vfirst, _) = wire_vertices(&the_wire);
    if let Some(v) = vfirst {
        ver_seq.push(v.0);
    }
}
