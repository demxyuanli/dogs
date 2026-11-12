//! `Build` phase 1: the two wire walks and `myPaths`.
//!
//! Source: `ModelingAlgorithms/TKOffset/BRepOffsetAPI/BRepOffsetAPI_MiddlePath.cxx`
//! (OCCT 8.0.0): wire walks `cxx:339-386`, `myPaths` initialization
//! `cxx:390-421`, `myPaths` filling `cxx:424-497`.

use crate::abs::ShapeType;
use crate::brep_offset_api_middle_path::{BRepOffsetApiMiddlePath, Ctx};
use crate::brep_tools_wire_explorer::{edge_vertices_cum_ori, last_vertex_cum_ori, WireExplorer};
use crate::iterator::cumulated_children;
use crate::shape::{Edge, TopoShape, Vertex};
use crate::topo_tools::{map_shapes_and_ancestors, shape_map_key, AncestorMap};

/// `Build` up to and including the filling of `myPaths` (`cxx:339-497`).
///
/// Answers `false` for the `return` arm of `cxx:416-421` (a start vertex with
/// no free starting edge), which leaves the object not done.
pub(super) fn collect(mp: &mut BRepOffsetApiMiddlePath, ctx: &mut Ctx) -> Result<bool, String> {
    // `cxx:339-342`.
    let mut start_vertices: Vec<TopoShape> = Vec::new();
    let mut edge_seq: Vec<Option<TopoShape>> = Vec::new();

    // `cxx:344-352`.
    let mut wexp = WireExplorer::of_wire(&mp.my_start_wire);
    while wexp.more() {
        if let Some(v) = wexp.current_vertex() {
            start_vertices.push(v.0.clone());
        }
        if let Some(e) = wexp.current() {
            edge_seq.push(Some(e.0.clone()));
        }
        wexp.next();
    }
    // `cxx:351-353`.
    if !mp.my_closed_section {
        if let Some(v) = wexp.current_vertex() {
            start_vertices.push(v.0.clone());
        }
    }
    // `cxx:354`.
    ctx.edge_seq = edge_seq;
    ctx.sections_edges.push(ctx.edge_seq.clone());

    // `cxx:356-365`.
    let mut wexp = WireExplorer::of_wire(&mp.my_end_wire);
    while wexp.more() {
        if let Some(v) = wexp.current_vertex() {
            ctx.end_vertices.insert(shape_map_key(&v.0));
        }
        if let Some(e) = wexp.current() {
            ctx.end_edges.insert(shape_map_key(&e.0));
        }
        wexp.next();
    }
    if !mp.my_closed_section {
        if let Some(v) = wexp.current_vertex() {
            ctx.end_vertices.insert(shape_map_key(&v.0));
        }
    }

    // `cxx:367-375`: `TopoDS_Iterator` over both wires.
    for child in cumulated_children(&mp.my_start_wire.0) {
        if child.shape_type() == ShapeType::Edge {
            mp.my_start_wire_edges.insert(shape_map_key(&child));
        }
    }
    for child in cumulated_children(&mp.my_end_wire.0) {
        if child.shape_type() == ShapeType::Edge {
            mp.my_end_wire_edges.insert(shape_map_key(&child));
        }
    }

    // `cxx:377-382`.
    ctx.vemap = AncestorMap::map_shapes_and_ancestors(
        &mp.my_initial_shape,
        ShapeType::Vertex,
        ShapeType::Edge,
    );
    ctx.efmap = map_shapes_and_ancestors(&mp.my_initial_shape, ShapeType::Edge, ShapeType::Face);

    // `cxx:390-421`: the first free edge at every start vertex.
    for i in 0..start_vertices.len() {
        let mut edges: Vec<TopoShape> = Vec::new();
        let le: Vec<TopoShape> = ctx
            .vemap
            .find_from_key(&start_vertices[i])
            .map(|s| s.to_vec())
            .unwrap_or_default();
        for it in le {
            let Some(mut an_edge) = Edge::wrap(it) else {
                return Err("Standard_TypeMismatch: TopoDS::Edge".into());
            };
            // `cxx:397`.
            if mp.my_start_wire_edges.contains(&shape_map_key(&an_edge.0)) {
                continue;
            }
            // `cxx:399-411`.
            let (v1, v2) = edge_vertices_cum_ori(&an_edge);
            if let Some(v1) = &v1 {
                if crate::topo_tools::is_same(&v1.0, &start_vertices[i]) {
                    if let Some(v2) = &v2 {
                        ctx.cur_vertices.insert(shape_map_key(&v2.0));
                    }
                } else {
                    an_edge.0.reverse();
                    if let Some(v2) = &v2 {
                        ctx.cur_vertices.insert(shape_map_key(&v2.0));
                    }
                }
            }
            edges.push(an_edge.0.clone());
            break;
        }
        // `cxx:413-421`.
        if edges.is_empty() {
            return Ok(false);
        }
        mp.my_paths.push(edges);
    }

    // `cxx:424-497`.
    fill_paths(mp, ctx);
    Ok(true)
}

/// `Filling of "myPaths"` (`cxx:424-497`).
fn fill_paths(mp: &mut BRepOffsetApiMiddlePath, ctx: &mut Ctx) {
    let mut next_vertices: Vec<TopoShape> = Vec::new();
    loop {
        // `cxx:427-488`.
        for i in 0..mp.my_paths.len() {
            let the_shape = mp.my_paths[i].last().cloned().expect("path");
            let (the_edge, the_vertex): (Option<Edge>, Option<Vertex>) =
                if the_shape.shape_type() == ShapeType::Edge {
                    // `cxx:433-436`.
                    let e = Edge::wrap(the_shape).expect("path element cast to edge");
                    let v = last_vertex_cum_ori(&e);
                    (Some(e), v)
                } else {
                    // `cxx:437-441`: the last segment was punctual.
                    let n = mp.my_paths[i].len();
                    let e = Edge::wrap(mp.my_paths[i][n - 2].clone())
                        .expect("previous path element cast to edge");
                    (Some(e), Vertex::wrap(the_shape))
                };

            // `cxx:443-446`.
            if let Some(v) = &the_vertex {
                if ctx.end_vertices.contains(&shape_map_key(&v.0)) {
                    continue;
                }
            }
            let Some(the_vertex) = the_vertex else {
                continue;
            };

            // `cxx:448-470`.
            let le: Vec<TopoShape> = ctx
                .vemap
                .find_from_key(&the_vertex.0)
                .map(|s| s.to_vec())
                .unwrap_or_default();
            let mut next_edge_candidates: Vec<Edge> = Vec::new();
            for it in le {
                let Some(mut an_edge) = Edge::wrap(it) else {
                    continue;
                };
                // `cxx:452-455`.
                if let Some(e) = &the_edge {
                    if crate::topo_tools::is_same(&an_edge.0, &e.0) {
                        continue;
                    }
                }
                let (v1, v2) = edge_vertices_cum_ori(&an_edge);
                // `cxx:456-466`: `if (V1.IsSame(theVertex)) NextVertex = V2;
                // else { anEdge.Reverse(); NextVertex = V1; }`.
                let v1_is_the_vertex = match &v1 {
                    Some(a) => crate::topo_tools::is_same(&a.0, &the_vertex.0),
                    None => false,
                };
                let next_vertex = if v1_is_the_vertex {
                    v2.clone()
                } else {
                    an_edge.0.reverse();
                    v1.clone()
                };
                let Some(next_vertex) = next_vertex else {
                    continue;
                };
                // `cxx:467-469`.
                if !ctx.cur_vertices.contains(&shape_map_key(&next_vertex.0)) {
                    next_edge_candidates.push(an_edge);
                }
            }

            // `cxx:471-486`.
            if next_edge_candidates.is_empty() {
                continue;
            }
            if next_edge_candidates.len() > 1 {
                // `cxx:474-476`: a punctual segment of the path.
                mp.my_paths[i].push(the_vertex.0.clone());
            } else {
                let an_edge = next_edge_candidates.remove(0);
                mp.my_paths[i].push(an_edge.0.clone());
                // `cxx:483`: `NextVertex = TopExp::LastVertex(anEdge, true)`.
                if let Some(nv) = last_vertex_cum_ori(&an_edge) {
                    next_vertices.push(nv.0.clone());
                }
            }
        }

        // `cxx:488-490`.
        if next_vertices.is_empty() {
            break;
        }
        // `cxx:492-496`.
        for it in &next_vertices {
            ctx.cur_vertices.insert(shape_map_key(it));
        }
        next_vertices.clear();
    }
}
