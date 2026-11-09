//! `SplitUEdges` and the bounded-edge rebuild it needs.
//!
//! Moved here verbatim from the single-file port; the `BRepFill_TrimShellCorner.cxx`
//! branch structure is unchanged.

use super::*;
use super::choose_section::*;

/// `SplitUEdges(theUEdges, theDS, theCrossDirection, theHistMap)`
/// (`cxx:1289-1392`).
///
/// For every row (the pair of consecutive `u`-iso edges of one face) that is
/// not already replaced in the map: look for a common vertex on the two edges;
/// when the `EE` search fails, fall back to the `VV` interference of the row's
/// `LastVertex(e1)` / `FirstVertex(e2)`. On success the two edges are re-cut at
/// the new vertex and bound in the map.
pub(super) fn split_u_edges(
    the_u_edges: &ShapeTable,
    the_ds: &BopdsDS,
    the_cross_direction: &GpVec,
    the_hist_map: &mut HashMap<ShapeMapKey, Vec<TopoShape>>,
) -> bool {
    let low_col = 1usize;
    let up_col = the_u_edges.nb_cols();
    let up_row = the_u_edges.nb_rows();

    for ueit in 1..=up_row {
        let Some(a_e1) = the_u_edges.value(ueit, low_col).cloned() else {
            continue;
        };
        let Some(a_e2) = the_u_edges.value(ueit, up_col).cloned() else {
            continue;
        };
        if the_hist_map.contains_key(&shape_map_key(&a_e1))
            || the_hist_map.contains_key(&shape_map_key(&a_e2))
        {
            continue;
        }
        let (Some(an_e_index1), Some(an_e_index2)) = (the_ds.index(&a_e1), the_ds.index(&a_e2))
        else {
            continue;
        };

        let mut a_common_vertex: Option<TopoShape> = None;
        let mut apar1 = 0.0f64;
        let mut apar2 = 0.0f64;
        let mut b_vertex_found = find_common_vertex(
            the_ds,
            an_e_index1,
            an_e_index2,
            the_cross_direction,
            &mut a_common_vertex,
            &mut apar1,
            &mut apar2,
        );

        if !b_vertex_found {
            // `cxx:1327-1353`: `VV` fallback on the row's end/start vertices.
            let v1 = edge_vertices(&Edge(a_e1.clone())).1;
            let v2 = edge_vertices(&Edge(a_e2.clone())).0;
            if let (Some(v1), Some(v2)) = (v1, v2) {
                let vindex1 = the_ds.index(&v1.0).unwrap_or(usize::MAX);
                let vindex2 = the_ds.index(&v2.0).unwrap_or(usize::MAX);
                for vv in the_ds.interf_vv() {
                    let matches = (vindex1 == vv.index1() && vindex2 == vv.index2())
                        || (vindex1 == vv.index2() && vindex2 == vv.index1());
                    if !matches {
                        continue;
                    }
                    let Some(n_new) = vv.new_index() else {
                        continue;
                    };
                    a_common_vertex = the_ds.shape(n_new).cloned();
                    // `cxx:1349-1350`: `BRep_Tool::Parameter(V1, Edge(aE1))`
                    // and the same on `aE2`.
                    apar1 = BRepTool::parameter_of_vertex(&a_e1, &v1.0);
                    apar2 = BRepTool::parameter_of_vertex(&a_e2, &v2.0);
                    b_vertex_found = true;
                    break;
                }
            }
        }

        if !b_vertex_found {
            continue;
        }
        let Some(a_common_vertex) = a_common_vertex else {
            continue;
        };

        // `cxx:1357-1385`: `EmptyCopied` + `Add(E, V1)` + `Add(E, V2)` +
        // `Range(E, f, l)`. The copy keeps the geometry of the source edge, so
        // the port rebuilds it on the same curve with the new bounds.
        let (Some(curve1), Some(curve2)) = (
            BRepTool::edge_curve(&Edge(a_e1.clone())),
            BRepTool::edge_curve(&Edge(a_e2.clone())),
        ) else {
            continue;
        };
        let (f1, _) = BRepTool::edge_parameters(&Edge(a_e1.clone()));
        let (_, l2) = BRepTool::edge_parameters(&Edge(a_e2.clone()));
        let verts1 = edge_vertices(&Edge(a_e1.clone()));
        let verts2 = edge_vertices(&Edge(a_e2.clone()));
        let (Some(a_v1), Some(a_v2)) = (verts1.0, verts2.1) else {
            continue;
        };

        let Some(new_e1) = rebuild_bounded_edge(
            &curve1,
            &a_v1,
            &Vertex(a_common_vertex.clone()),
            f1,
            apar1,
        ) else {
            continue;
        };
        let Some(new_e2) = rebuild_bounded_edge(
            &curve2,
            &Vertex(a_common_vertex.clone()),
            &a_v2,
            apar2,
            l2,
        ) else {
            continue;
        };

        the_hist_map
            .insert(shape_map_key(&a_e1), vec![new_e1.0.clone()]);
        the_hist_map
            .insert(shape_map_key(&a_e2), vec![new_e2.0.clone()]);
    }
    true
}

/// The `B.Begin/Add(E, V)/Range(E, f, l)` sequence of
/// `BRepFill_TrimShellCorner.cxx:1360-1385` on the curve of an existing edge.
pub(super) fn rebuild_bounded_edge(
    curve: &std::sync::Arc<dyn occt_geom::curve::Curve>,
    v_first: &Vertex,
    v_last: &Vertex,
    first: f64,
    last: f64,
) -> Option<Edge> {
    let b = TopoBuilder::new();
    let build = init_3d_vertices(curve, Some(v_first), Some(v_last), first, last)?;
    let mut e = b.make_edge(curve.clone(), build.first, build.last);
    add_endpoint_vertices(&b, &mut e, &build);
    Some(e)
}

