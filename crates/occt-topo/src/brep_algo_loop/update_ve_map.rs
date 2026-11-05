//! `BRepAlgo_Loop::UpdateVEmap` (`BRepAlgo_Loop.cxx:1018-1151`): merge the
//! vertex -> edges map of the edges that carry a single image, then rebuild the
//! map from the surviving edges.

use std::collections::HashSet;

use occt_core::precision::CONFUSION;

use crate::abs::{Orientation, ShapeType};
use crate::brep_algo_image::BRepAlgoImage;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::iterator::cumulated_children;
use crate::shape::{Edge, TopoShape, Vertex};
use crate::topo_tools::{edge_vertices, is_same, shape_map_key, ShapeMapKey};

use super::mve::VertexEdgeMap;

/// `UpdateVEmap(theVEmap)` (`cxx:1018-1151`).
pub fn update_ve_map(
    the_ve_map: &mut VertexEdgeMap,
    image_vv: &mut BRepAlgoImage,
    const_edges: &[TopoShape],
) {
    // `cxx:1021-1046`: group the single-edge vertices by their `myImageVV`
    // image. `NCollection_IndexedDataMap` keeps insertion order.
    let mut ver_lver: Vec<(TopoShape, Vec<TopoShape>)> = Vec::new();
    let mut ver_lver_index: Vec<ShapeMapKey> = Vec::new();

    for ii in 1..=the_ve_map.extent() {
        let Some(a_vertex) = the_ve_map.find_key(ii).cloned() else {
            continue;
        };
        let a_e_list = the_ve_map.find_from_index(ii).cloned().unwrap_or_default();
        if a_e_list.len() == 1 && image_vv.is_image(&a_vertex) {
            let Ok(pro) = image_vv.image_from(&a_vertex) else {
                continue;
            };
            let key = shape_map_key(&pro);
            let entry = if let Some(i) = ver_lver_index.iter().position(|k| *k == key) {
                Some(i)
            } else {
                ver_lver.push((pro, Vec::new()));
                ver_lver_index.push(key);
                Some(ver_lver.len() - 1)
            };
            if let Some(i) = entry {
                ver_lver[i].1.push(a_vertex.oriented(Orientation::Forward));
            }
        }
    }

    if ver_lver.is_empty() {
        return;
    }

    let builder = TopoBuilder::new();

    // `cxx:1052-1135`.
    for i in 0..ver_lver.len() {
        let a_v_list = ver_lver[i].1.clone();
        if a_v_list.len() == 1 {
            continue;
        }

        // `cxx:1060-1077`: the barycentre axis of the vertex cloud.
        let mut a_max_tol = 0.0f64;
        let mut points = Vec::with_capacity(a_v_list.len());
        for a_vertex in &a_v_list {
            let a_tol = BRepTool::vertex_tolerance(&Vertex(a_vertex.clone()));
            if a_tol > a_max_tol {
                a_max_tol = a_tol;
            }
            points.push(BRepTool::vertex_point(&Vertex(a_vertex.clone())));
        }

        let (an_axis, _is_singular) =
            occt_geom::geom_lib::geom_lib_axe_of_inertia(&points, CONFUSION);
        let a_centre = an_axis.location();
        let mut a_max_dist = 0.0f64;
        for p in &points {
            let sq = a_centre.square_distance(p);
            if sq > a_max_dist {
                a_max_dist = sq;
            }
        }
        let a_max_dist = a_max_dist.sqrt();
        if a_max_dist > a_max_tol {
            a_max_tol = a_max_dist;
        }

        // `cxx:1088-1111`: the constant vertex, i.e. the one on a const edge.
        let mut a_const_vertex: Option<TopoShape> = None;
        for a_vertex in &a_v_list {
            let Some(a_e_list) = the_ve_map.find_from_key(a_vertex).cloned() else {
                continue;
            };
            let Some(an_edge) = a_e_list.first() else {
                continue;
            };
            if const_edges.iter().any(|c| is_same(c, an_edge)) {
                a_const_vertex = Some(a_vertex.clone());
                break;
            }
        }
        if a_const_vertex.is_none() {
            a_const_vertex = a_v_list.first().cloned();
        }
        let Some(a_const_vertex) = a_const_vertex else {
            continue;
        };

        // `cxx:1112`: `aBB.UpdateVertex(aConstVertex, aCentre, aMaxTol)`.
        let cv = Vertex(a_const_vertex.clone());
        cv.set_point(a_centre);
        cv.set_tolerance(a_max_tol);

        // `cxx:1114-1132`: re-attach every other vertex of the group to it.
        for a_vertex in &a_v_list {
            if is_same(a_vertex, &a_const_vertex) {
                continue;
            }
            let Some(a_e_list) = the_ve_map.find_from_key(a_vertex).cloned() else {
                continue;
            };
            let Some(an_edge) = a_e_list.first().cloned() else {
                continue;
            };
            let mut an_edge = an_edge.oriented(Orientation::Forward);
            let (a_v1, a_v2) = edge_vertices(&Edge(an_edge.clone()));
            let vertex_to_remove = match (&a_v1, &a_v2) {
                (Some(v1), _) if is_same(&v1.0, a_vertex) => v1.0.clone(),
                (_, Some(v2)) => v2.0.clone(),
                _ => continue,
            };
            an_edge.set_free(true);
            builder.remove(&mut an_edge, &vertex_to_remove);
            builder.add(
                &mut an_edge,
                &a_const_vertex.oriented(vertex_to_remove.orientation()),
            );
        }
    }

    // `cxx:1137-1150`: rebuild the map from the distinct edges.
    let mut emap: Vec<TopoShape> = Vec::new();
    let mut seen: HashSet<ShapeMapKey> = HashSet::new();
    for ii in 1..=the_ve_map.extent() {
        let a_e_list = the_ve_map.find_from_index(ii).cloned().unwrap_or_default();
        for e in a_e_list {
            if seen.insert(shape_map_key(&e)) {
                emap.push(e);
            }
        }
    }

    *the_ve_map = VertexEdgeMap::new();
    for e in emap {
        for v in cumulated_children(&e) {
            if v.shape_type() != ShapeType::Vertex {
                continue;
            }
            if !the_ve_map.contains(&v) {
                the_ve_map.add(&v, Vec::new());
            }
            if let Some(list) = the_ve_map.change_from_key(&v) {
                list.push(e.clone());
            }
        }
    }
}
