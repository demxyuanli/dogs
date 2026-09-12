//! `BOPAlgo_PaveFiller::PerformVV` and `MakeSDVertices`
//! (`BOPAlgo_PaveFiller_1.cxx:45/136`).
//!
//! Interfering vertex pairs are clustered (`FillMap` / `MakeBlocks`). Each
//! connected component is fused into one same-domain vertex. Vertices of the
//! same argument that land in one component raise a self-interference warning.
//! After fusion, `InitPaveBlocksForVertex` is called for every SD key.

use std::collections::HashMap;

use occt_core::precision::CONFUSION;

use crate::abs::ShapeType;
use crate::algo_tools::AlgoTools;
use crate::algo_tools_range::make_vertex_from_list;
use crate::bopalgo_tools::{fill_map_pair, make_blocks_int};
use crate::bbox_from_geometry::shape_bbox;
use crate::brep_tool::BRepTool;
use crate::pave_force_ee::init_pave_blocks_for_vertex;
use crate::pave_filler::PaveFiller;
use crate::pave_intersect::collect_pairs;
use crate::shape::{TopoShape, Vertex};

/// `BOPAlgo_PaveFiller::MakeSDVertices`.
pub fn make_sd_vertices(
    f: &mut PaveFiller,
    the_vert_indices: &[usize],
    the_add_interfs: bool,
) -> Result<usize, String> {
    if the_vert_indices.is_empty() {
        return Err("make_sd_vertices: empty vertex list".into());
    }
    let mut n_sd: Option<usize> = None;
    let mut a_lv: Vec<TopoShape> = Vec::new();
    for &n_x in the_vert_indices {
        if let Some(n_sd1) = f.ds().has_shape_sd(n_x) {
            let Some(a_vsd1) = f.ds().shape(n_sd1).cloned() else {
                continue;
            };
            if n_sd.is_none() {
                n_sd = Some(n_sd1);
            } else {
                a_lv.push(a_vsd1);
            }
        }
        if let Some(a_v) = f.ds().shape(n_x).cloned() {
            a_lv.push(a_v);
        }
    }
    let a_vn = make_vertex_from_list(&a_lv)?;
    let n_v = if let Some(n) = n_sd {
        if let Some(shape) = f.ds().shape(n).cloned() {
            let vtx = Vertex(shape);
            let p = BRepTool::vertex_point(&Vertex(a_vn.clone()));
            let tol = BRepTool::vertex_tolerance(&Vertex(a_vn.clone()));
            vtx.set_point(p);
            vtx.set_tolerance(tol);
        }
        n
    } else {
        f.ds_mut().append(a_vn.clone())?
    };
    if let Some(shape) = f.ds().shape(n_v).cloned() {
        let mut a_box = shape_bbox(&shape);
        let tol = BRepTool::vertex_tolerance(&Vertex(shape));
        a_box.set_gap(tol + CONFUSION);
        let _ = a_box;
    }
    for (i, &n1) in the_vert_indices.iter().enumerate() {
        f.ds_mut().add_shape_sd(n1, n_v);
        let i_r1 = f.ds().rank(n1);
        for &n2 in &the_vert_indices[i + 1..] {
            if i_r1 == f.ds().rank(n2) {
                f.add_warning(format!(
                    "self-interfering shape: vertices {n1} and {n2} of one argument"
                ));
            }
            if the_add_interfs {
                f.ds_mut().add_interf_vv(n1, n2, Some(n_v));
            }
        }
    }
    Ok(n_v)
}

/// `BOPAlgo_PaveFiller::PerformVV`.
pub fn perform_vv(f: &mut PaveFiller) -> Result<(), String> {
    let pairs = collect_pairs(f.ds(), ShapeType::Vertex, ShapeType::Vertex);
    if pairs.is_empty() {
        return Ok(());
    }
    let mut a_mili: HashMap<usize, Vec<usize>> = HashMap::new();
    let fuzzy = f.fuzzy_value();
    for (n1, n2) in pairs {
        if f.ds().has_interf_pair(n1, n2) {
            fill_map_pair(n1, n2, &mut a_mili);
            continue;
        }
        let n1sd = f.ds().has_shape_sd(n1).unwrap_or(n1);
        let n2sd = f.ds().has_shape_sd(n2).unwrap_or(n2);
        let Some(v1s) = f.ds().shape(n1sd).cloned() else {
            continue;
        };
        let Some(v2s) = f.ds().shape(n2sd).cloned() else {
            continue;
        };
        let p2 = BRepTool::vertex_point(&Vertex(v2s));
        if AlgoTools::compute_vv(&v1s, &p2, fuzzy) == 1 {
            fill_map_pair(n1, n2, &mut a_mili);
        }
    }
    let blocks = make_blocks_int(&a_mili);
    for a_li in blocks {
        make_sd_vertices(f, &a_li, true)?;
    }
    let keys: Vec<usize> = f.ds().shapes_sd().keys().copied().collect();
    for n1 in keys {
        init_pave_blocks_for_vertex(f, n1);
    }
    Ok(())
}

/// Restricted-pair VV used by `RepeatIntersection`.
pub fn perform_vv_pairs(f: &mut PaveFiller, pairs: &[(usize, usize)]) -> Result<(), String> {
    if pairs.is_empty() {
        return Ok(());
    }
    let mut a_mili: HashMap<usize, Vec<usize>> = HashMap::new();
    let fuzzy = f.fuzzy_value();
    for &(n1, n2) in pairs {
        if f.ds().has_interf_pair(n1, n2) {
            fill_map_pair(n1, n2, &mut a_mili);
            continue;
        }
        let n1sd = f.ds().has_shape_sd(n1).unwrap_or(n1);
        let n2sd = f.ds().has_shape_sd(n2).unwrap_or(n2);
        let Some(v1s) = f.ds().shape(n1sd).cloned() else {
            continue;
        };
        let Some(v2s) = f.ds().shape(n2sd).cloned() else {
            continue;
        };
        let p2 = BRepTool::vertex_point(&Vertex(v2s));
        if AlgoTools::compute_vv(&v1s, &p2, fuzzy) == 1 {
            fill_map_pair(n1, n2, &mut a_mili);
        }
    }
    for a_li in make_blocks_int(&a_mili) {
        make_sd_vertices(f, &a_li, true)?;
    }
    Ok(())
}
