//! `TopOpeBRepBuild_Tools2d` — the wire path walker of `CorrectFace2d`.
//!
//! Source: `TKBool/TopOpeBRepBuild/TopOpeBRepBuild_Tools2d.cxx` (OCCT 8.0.0):
//! `Path` `:39-66`, the static `BuildPath` `:26-30` / `:70-158`,
//! `MakeMapOfShapeVertexInfo` `:162-210`, `DumpMapOfShapeVertexInfo` `:212-243`.
//!
//! The only caller is `TopOpeBRepBuild_CorrectFace2d::MakeRightWire`
//! (`TopOpeBRepBuild_CorrectFace2d.cxx:664-665`), which needs the edges of one
//! wire in connection order (`NCollection_List`, prepend order).
//!
//! `Path` is a depth-first walk: at every vertex it takes the next outgoing
//! edge (`VertexInfo::CurrentOut`), pushes it on the result list and recurses on
//! the far vertex. The walk succeeds only when it gets back to the first vertex
//! after passing exactly as many edges as the wire stores.

use crate::abs::{Orientation, ShapeType};
use crate::shape::{Edge, TopoShape};
use crate::topo_tools::{edge_vertices, explorer_shapes, is_same};
use crate::topopebrep::build::vertex_info::VertexInfo;
use crate::topopebrep::ds::data_structure::IndexedShapeData;

/// `TopExp::Vertices(E, V1, V2, CumOri = true)`.
///
/// `TopExp.cxx:Vertices` fills `V1` from the `TopAbs_FORWARD` vertex of the
/// edge and `V2` from the `TopAbs_REVERSED` one, then swaps them when the edge
/// itself is `TopAbs_REVERSED`.
fn edge_vertices_cumulative(e: &TopoShape) -> (Option<TopoShape>, Option<TopoShape>) {
    let Some(edge) = Edge::wrap(e.clone()) else {
        return (None, None);
    };
    let (v1, v2) = edge_vertices(&edge);
    let (v1, v2) = (v1.map(|v| v.0), v2.map(|v| v.0));
    if e.orientation() == Orientation::Reversed {
        (v2, v1)
    } else {
        (v1, v2)
    }
}

/// `TopExp::FirstVertex(E)`.
pub fn first_vertex(e: &TopoShape) -> Option<TopoShape> {
    edge_vertices_cumulative(e).0
}

/// `TopExp::LastVertex(E)`.
pub fn last_vertex(e: &TopoShape) -> Option<TopoShape> {
    edge_vertices_cumulative(e).1
}

/// `TopOpeBRepBuild_Tools2d::Path(aWire, aResList)` (`cxx:39-66`).
///
/// Returns the wire's edges in the order `BuildPath` collected them, or an
/// empty list when the walk fails (`aResList` is cleared unconditionally in
/// OCCT, `cxx:64-65`, so a failed walk also answers an empty list there).
pub fn path(a_wire: &TopoShape) -> Vec<TopoShape> {
    let a_nb_edges = explorer_shapes(a_wire, ShapeType::Edge).len() as i32;

    let mut m: IndexedShapeData<VertexInfo> = IndexedShapeData::new();
    make_map_of_shape_vertex_info(a_wire, &mut m);

    let Some(my_vertex0) = m.find_key(1).cloned() else {
        // `M.FindKey(1)` on an empty map: OCCT reads out of range (`cxx:59`).
        return Vec::new();
    };

    let mut my_res_list: Vec<TopoShape> = Vec::new();
    let mut a_break_flag = 0;
    build_path(
        &my_vertex0,
        None,
        &my_vertex0,
        a_nb_edges,
        &mut m,
        0,
        &mut a_break_flag,
        &mut my_res_list,
    );
    my_res_list
}

/// `BuildPath(myVertex0, myInputEdge, myInputVertex, aNbEdges, M, anEdgesCount,
/// aBreakFlag, myResList)` (`cxx:70-158`).
///
/// `anEdgesCount` is passed by value (C++ `int`), so the count of a failed
/// branch does not leak into the caller; `aBreakFlag` is `1` (path closed) /
/// `2` (dead end) and stops the whole walk.
///
/// The `anEdgesCount--` right before the `return` at `cxx:148` writes a
/// by-value copy that is never read again; the port keeps the statement for
/// control-flow parity and allows the resulting dead-assignment lint.
#[allow(unused_assignments)]
fn build_path(
    my_vertex0: &TopoShape,
    my_input_edge: Option<&TopoShape>,
    my_input_vertex: &TopoShape,
    a_nb_edges: i32,
    m: &mut IndexedShapeData<VertexInfo>,
    an_edges_count: i32,
    a_break_flag: &mut i32,
    my_res_list: &mut Vec<TopoShape>,
) {
    let mut stop_flag = 0;

    if *a_break_flag == 1 {
        return;
    }

    let Some(v_info) = m.change_from_key(my_input_vertex) else {
        // `M.ChangeFromKey` on a missing vertex throws in OCCT (`cxx:89`); the
        // map is built from the same wire, so this arm is unreachable.
        return;
    };
    v_info.set_current_in(my_input_edge);
    v_info.prepare(my_res_list.as_slice());
    let a_nb_cases = v_info.nb_cases();
    if a_nb_cases == 0 {
        *a_break_flag = 2;
    }

    let mut an_edges_count = an_edges_count;

    for j in 1..=a_nb_cases {
        let Some(v_info) = m.change_from_key(my_input_vertex) else {
            return;
        };
        let my_edge = v_info.current_out();
        let a_found_out = v_info.found_out();
        if a_found_out == 0 {
            // `FoundOut == 0`: no outgoing edge left at this vertex.
            *a_break_flag = 2;
            return;
        }
        let Some(my_edge) = my_edge else {
            continue;
        };

        if stop_flag != 0 {
            // The previous sub-path was wrong: drop its last edge.
            v_info.remove_passed();
            if !my_res_list.is_empty() {
                my_res_list.remove(0);
            }
            stop_flag = 0;
            an_edges_count -= 1;
        }

        v_info.append_passed(&my_edge);
        my_res_list.insert(0, my_edge.clone());
        an_edges_count += 1;

        let my_vertex = if my_edge.orientation() == Orientation::Forward {
            last_vertex(&my_edge)
        } else {
            first_vertex(&my_edge)
        };
        let Some(my_vertex) = my_vertex else {
            continue;
        };

        if is_same(&my_vertex, my_vertex0) && an_edges_count == a_nb_edges {
            *a_break_flag = 1;
            return;
        }

        build_path(
            my_vertex0,
            Some(&my_edge),
            &my_vertex,
            a_nb_edges,
            m,
            an_edges_count,
            a_break_flag,
            my_res_list,
        );

        if *a_break_flag == 1 {
            return;
        }

        if *a_break_flag == 2 {
            // Come back to this vertex and go on with the next outgoing edge.
            let Some(v_info) = m.change_from_key(my_input_vertex) else {
                return;
            };
            if j == a_nb_cases {
                v_info.remove_passed();
                if !my_res_list.is_empty() {
                    my_res_list.remove(0);
                }
                an_edges_count -= 1;
                return;
            }
            stop_flag = 1;
            *a_break_flag = 0;
        }
    }
}

/// `MakeMapOfShapeVertexInfo(aWire, M)` (`cxx:162-210`).
///
/// Every vertex of every edge of the wire maps to a [`VertexInfo`] holding the
/// edge as `EdgesOut` (vertex `TopAbs_FORWARD`) or `EdgesIn` (vertex
/// `TopAbs_REVERSED`); a vertex with neither a single in nor a single out edge
/// is marked `Smart`.
pub fn make_map_of_shape_vertex_info(a_wire: &TopoShape, m: &mut IndexedShapeData<VertexInfo>) {
    for an_edge in explorer_shapes(a_wire, ShapeType::Edge) {
        for a_vertex in explorer_shapes(&an_edge, ShapeType::Vertex) {
            if m.find_index(&a_vertex) == 0 {
                m.add(&a_vertex, VertexInfo::new());
            }
            let Some(a_v_info) = m.change_from_key(&a_vertex) else {
                continue;
            };
            a_v_info.set_vertex(&a_vertex);
            let an_or = a_vertex.orientation();
            if an_or == Orientation::Forward {
                a_v_info.add_out(&an_edge);
            } else if an_or == Orientation::Reversed {
                a_v_info.add_in(&an_edge);
            }
        }
    }

    let a_nb = m.extent();
    for i in 1..=a_nb {
        let Some(a_v_info) = m.change_from_index(i) else {
            continue;
        };
        let a_nb_edges_in = a_v_info.edges_in().extent();
        let a_nb_edges_out = a_v_info.edges_out().extent();
        if a_nb_edges_in != 1 && a_nb_edges_out != 1 {
            a_v_info.set_smart(true);
        }
    }
}

/// `DumpMapOfShapeVertexInfo(M)` (`cxx:212-243`): `printf` output only, no
/// port; kept so the caller's name resolution does not change.
pub fn dump_map_of_shape_vertex_info(_m: &IndexedShapeData<VertexInfo>) {}
