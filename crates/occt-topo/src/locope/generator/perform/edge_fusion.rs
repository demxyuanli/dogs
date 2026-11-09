//! `LocOpe_Generator.cxx:442-522`: the fusion of the twin edges.
//!
//! Every edge bound in `theEEMap` is rebuilt from an empty copy: for each of its
//! vertices the reference either takes the vertex the twin edge carries, with
//! the parameter of that twin vertex re-expressed on the edge (`NewParameter`),
//! or keeps the original vertex. A rebuilt edge whose vertices moved goes to
//! `myModShapes` and its twin and itself go to `toRemove`.

use crate::brep_builder_update_vertex::update_vertex_on_edge;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::locope::shape_maps::ShapeValueMap;
use crate::shape::{Edge, TopoShape, Vertex};
use crate::topo_tools::is_same;

use super::super::{empty_copied_edge, statics};
use super::Stage;

/// `cxx:454-522`.
pub(super) fn fuse_edges(stage: &mut Stage, mod_shapes: &mut ShapeValueMap<Vec<TopoShape>>) {
    let entries: Vec<(TopoShape, Vec<TopoShape>)> = stage
        .ee_map
        .iter()
        .map(|(k, v)| (k.clone(), v.clone()))
        .collect();
    let builder = TopoBuilder::new();
    for (edg, twins) in entries {
        let Some(edge) = Edge::wrap(edg.clone()) else {
            continue;
        };
        // `cxx:461`: `BRep_Tool::Range(edg, f, l)`.
        let (f, l) = BRepTool::edge_parameters(&edge);
        // `cxx:462-467`.
        let mut new_edge = empty_copied_edge(&edge);
        let mut keep_new_edge = false;
        let mut last_twin: Option<TopoShape> = None;
        for vtx in crate::topo_tools::explorer_shapes(
            &edg.oriented(crate::abs::Orientation::Forward),
            crate::abs::ShapeType::Vertex,
        ) {
            // `cxx:469`.
            let mut prm = BRepTool::parameter_of_vertex(&edg, &vtx);
            let mut new_vtx: Option<TopoShape> = None;
            // `cxx:471-494`: the first twin that carries `vtx` and another
            // vertex hands that other vertex over.
            for twin in &twins {
                last_twin = Some(twin.clone());
                let twin_vertices =
                    crate::topo_tools::explorer_shapes(twin, crate::abs::ShapeType::Vertex);
                if !twin_vertices.iter().any(|x| is_same(x, &vtx)) {
                    continue;
                }
                for other in &twin_vertices {
                    if is_same(other, &vtx) {
                        continue;
                    }
                    new_vtx = Some(other.clone());
                    prm = statics::new_parameter(&edg, &vtx, &new_edge.0, other);
                    break;
                }
                break;
            }
            // `cxx:497-511`.
            if stage.to_remove.contains(&vtx) || (prm < l && prm > f) {
                // The reference adds and updates `newvtx` even when the twin
                // scan left it null (`cxx:499-501`); the port keeps a null out
                // of the edge instead.
                let Some(new_vtx) = new_vtx else {
                    continue;
                };
                builder.add(&mut new_edge.0, &new_vtx.oriented(vtx.orientation()));
                let tol = Vertex::wrap(new_vtx.clone())
                    .map(|v| BRepTool::vertex_tolerance(&v))
                    .unwrap_or(0.0);
                // `cxx:501`: `B.UpdateVertex(NewVtx, prm, newedg, tol)`. The
                // only failure the reference raises for is an infinite
                // parameter (`BRep_Builder.cxx:1225-1228`).
                let _ = update_vertex_on_edge(&new_vtx, prm, &new_edge.0, tol);
                // `cxx:502`: the twin edge goes.
                if let Some(twin) = last_twin.clone() {
                    stage.to_remove.add(&twin);
                }
                keep_new_edge = true;
            } else {
                builder.add(&mut new_edge.0, &vtx.oriented(vtx.orientation()));
                let tol = Vertex::wrap(vtx.clone())
                    .map(|v| BRepTool::vertex_tolerance(&v))
                    .unwrap_or(0.0);
                let _ = update_vertex_on_edge(&vtx, prm, &new_edge.0, tol);
            }
        }
        // `cxx:513-521`.
        if keep_new_edge {
            if !mod_shapes.is_bound(&edg) {
                mod_shapes.bind(&edg, Vec::new());
            }
            if let Some(list) = mod_shapes.change(&edg) {
                list.push(new_edge.0.clone());
            }
            stage.to_remove.add(&edg);
        }
    }
}
