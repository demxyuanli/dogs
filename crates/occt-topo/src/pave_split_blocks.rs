//! `BOPAlgo_PaveFiller::SplitPaveBlocks` (`BOPAlgo_PaveFiller_2.cxx:419`).
//!
//! For every edge in `theMEdges` whose pave blocks carry extra paves, replace
//! each block by the elementary blocks produced by `BOPDS_PaveBlock::Update`.
//! After the split:
//!
//! - `UpdatePaveBlockWithSDVertices` + `FillShrunkData` run on each new block;
//! - a block with no valid shrunk range (or a non-splittable range whose
//!   bounding vertices interfere) has those vertices unified by
//!   `MakeSDVertices`, and `InitPaveBlocksForVertex` is queued for them;
//! - common-block members are regrouped: open blocks with the same bound
//!   vertices inherit the original common-block faces; closed blocks
//!   (`nV1 == nV2`) are clustered by `ComputePE` coincidence before a new
//!   common block is made.

use std::collections::{HashMap, HashSet};

use crate::abs::ShapeType;
use crate::algo_tools::AlgoTools;
use crate::algo_tools_range::point_on_edge;
use crate::bopds::{BopdsCommonBlock, BopdsPaveBlock};
use crate::brep_tool::BRepTool;
use crate::int_tools_full::IntToolsContext;
use crate::pave_common::fill_shrunk_data_for_block;
use crate::pave_ff_exist::pb_key;
use crate::pave_filler::PaveFiller;
use crate::pave_force_ee::init_pave_blocks_for_vertex;
use crate::pave_intersect::make_sd_vertices;
use crate::pave_ve::make_new_common_block;
use crate::shape::{Edge, Vertex};

fn cb_identity(cb: &BopdsCommonBlock) -> (Vec<usize>, Vec<usize>) {
    let mut edges = cb.indices().to_vec();
    edges.sort_unstable();
    let mut faces = cb.faces().to_vec();
    faces.sort_unstable();
    (edges, faces)
}

fn update_pave_block_with_sd_vertices(f: &PaveFiller, pb: &mut BopdsPaveBlock) {
    let (n1, n2) = pb.indices();
    pb.set_indices(
        f.ds().get_same_domain_index(n1),
        f.ds().get_same_domain_index(n2),
    );
}

fn fill_shrunk_on_block(f: &mut PaveFiller, pb: &mut BopdsPaveBlock) {
    let n_e = if pb.original_edge() != usize::MAX {
        pb.original_edge()
    } else {
        pb.edge()
    };
    let Some(es) = f.ds().shape(n_e).cloned() else {
        pb.clear_shrunk_data();
        return;
    };
    let edge = Edge(es);
    let tol = f.fuzzy_value();
    fill_shrunk_data_for_block(f, &edge, tol, pb);
}

fn vertices_interfere(f: &PaveFiller, n_v1: usize, n_v2: usize) -> bool {
    let Some(s1) = f.ds().shape(n_v1).cloned() else {
        return false;
    };
    let Some(s2) = f.ds().shape(n_v2).cloned() else {
        return false;
    };
    let p2 = BRepTool::vertex_point(&Vertex(s2));
    AlgoTools::compute_vv(&s1, &p2, f.fuzzy_value()) == 1
}

fn max_vertex_tolerance(f: &PaveFiller, n_e: usize) -> f64 {
    let mut tol = 0.0_f64;
    if let Some(si) = f.ds().shape_info(n_e) {
        for &sub in si.sub_shapes() {
            if f.ds().shape_info(sub).map(|s| s.shape_type()) != Some(ShapeType::Vertex) {
                continue;
            }
            if let Some(vs) = f.ds().shape(sub) {
                tol = tol.max(BRepTool::vertex_tolerance(&Vertex(vs.clone())));
            }
        }
    }
    if let Some(es) = f.ds().shape(n_e) {
        tol = tol.max(BRepTool::edge_tolerance(&Edge(es.clone())));
    }
    tol
}

/// `BOPAlgo_PaveFiller::SplitPaveBlocks`.
pub fn split_pave_blocks(
    f: &mut PaveFiller,
    the_m_edges: &[usize],
    the_add_interfs: bool,
) -> Result<(), String> {
    let mut edges: Vec<usize> = if the_m_edges.is_empty() {
        (0..f.ds().nb_source_shapes())
            .filter(|&i| {
                f.ds().shape_info(i).map(|s| s.shape_type()) == Some(ShapeType::Edge)
                    && f.ds().pave_blocks(i).iter().any(|pb| pb.is_to_update())
            })
            .collect()
    } else {
        the_m_edges.to_vec()
    };
    edges.sort_unstable();
    edges.dedup();
    if edges.is_empty() {
        return Ok(());
    }

    let mut pair_fence: HashSet<(usize, usize)> = HashSet::new();
    let mut vertices_to_init: HashSet<usize> = HashSet::new();
    let mut mcb_new_pb: HashMap<(Vec<usize>, Vec<usize>), Vec<BopdsPaveBlock>> = HashMap::new();
    let mut mcb_faces: HashMap<(Vec<usize>, Vec<usize>), Vec<usize>> = HashMap::new();
    let mut mcb_closed: HashMap<(Vec<usize>, Vec<usize>), bool> = HashMap::new();

    for n_e in edges {
        if !f.ds().has_pave_blocks(n_e) {
            continue;
        }
        let old = f.ds().pave_blocks(n_e).to_vec();
        let mut new_list: Vec<BopdsPaveBlock> = Vec::new();
        for mut pb in old {
            if !pb.is_to_update() {
                new_list.push(pb);
                continue;
            }
            let cb = f.ds().common_block(&pb).cloned();
            let mut lpbn: Vec<BopdsPaveBlock> = Vec::new();
            pb.update(&mut lpbn, true);
            for mut pbn in lpbn {
                update_pave_block_with_sd_vertices(f, &mut pbn);
                fill_shrunk_on_block(f, &mut pbn);
                let b_has_valid_range = pbn.has_shrunk_data();
                let b_check_dist = b_has_valid_range && !pbn.is_splittable();
                if !b_has_valid_range || b_check_dist {
                    let (n_v1, n_v2) = pbn.indices();
                    if n_v1 == n_v2 {
                        continue;
                    }
                    let mut unify = !b_has_valid_range;
                    if b_check_dist && vertices_interfere(f, n_v1, n_v2) {
                        unify = true;
                    }
                    if unify {
                        let key = (n_v1.min(n_v2), n_v1.max(n_v2));
                        if pair_fence.insert(key) {
                            if let Err(e) = make_sd_vertices(f.ds_mut(), &[n_v1, n_v2], the_add_interfs)
                            {
                                f.add_error(e);
                            }
                            vertices_to_init.insert(n_v1);
                            vertices_to_init.insert(n_v2);
                        }
                        continue;
                    }
                }
                new_list.push(pbn.clone());
                if let Some(ref a_cb) = cb {
                    let id = cb_identity(a_cb);
                    mcb_new_pb.entry(id.clone()).or_default().push(pbn);
                    mcb_faces.entry(id.clone()).or_insert_with(|| a_cb.faces().to_vec());
                    let closed = a_cb
                        .pave_block1()
                        .map(|p| {
                            let (a, b) = p.indices();
                            a == b
                        })
                        .unwrap_or(false);
                    mcb_closed.entry(id).or_insert(closed);
                }
            }
        }
        *f.ds_mut().change_pave_blocks_mut(n_e) = new_list;
    }

    rebuild_common_blocks(f, mcb_new_pb, mcb_faces, mcb_closed);

    for n_v in vertices_to_init {
        init_pave_blocks_for_vertex(f, n_v);
    }
    Ok(())
}

fn rebuild_common_blocks(
    f: &mut PaveFiller,
    mcb_new_pb: HashMap<(Vec<usize>, Vec<usize>), Vec<BopdsPaveBlock>>,
    mcb_faces: HashMap<(Vec<usize>, Vec<usize>), Vec<usize>>,
    mcb_closed: HashMap<(Vec<usize>, Vec<usize>), bool>,
) {
    for (id, lpbn) in mcb_new_pb {
        let faces = mcb_faces.get(&id).cloned().unwrap_or_default();
        let b_is_closed = mcb_closed.get(&id).copied().unwrap_or(false);
        let mut by_bounds: HashMap<(usize, usize), Vec<BopdsPaveBlock>> = HashMap::new();
        for pb in lpbn {
            let (n1, n2) = pb.indices();
            let key = (n1.min(n2), n1.max(n2));
            by_bounds.entry(key).or_default().push(pb);
        }
        for (_key, mut a_lpb) in by_bounds {
            if !b_is_closed {
                make_new_common_block(f, &a_lpb, &faces);
                continue;
            }
            cluster_closed_common_blocks(f, &mut a_lpb, &faces);
        }
    }
}

/// Closed common-block regrouping (`_2.cxx:572-616`).
///
/// Pave blocks of a closed original common block are clustered by projecting
/// the mid-point of the first member onto the remaining edges (`ComputePE`).
/// Each coincidence group becomes a new common block with the original faces.
fn cluster_closed_common_blocks(
    f: &mut PaveFiller,
    a_lpb: &mut Vec<BopdsPaveBlock>,
    faces: &[usize],
) {
    let ctx = IntToolsContext::new();
    let fuzzy = f.fuzzy_value();
    while !a_lpb.is_empty() {
        let mut a_lpb_cb: Vec<BopdsPaveBlock> = Vec::new();
        let mut a_pm_first = occt_core::gp::GpPnt::new(0.0, 0.0, 0.0);
        let mut a_tol_e_first = 0.0_f64;
        let mut rest: Vec<BopdsPaveBlock> = Vec::new();
        for a_pb in a_lpb.drain(..) {
            if a_lpb_cb.is_empty() {
                a_lpb_cb.push(a_pb.clone());
                let n_e = a_pb.original_edge();
                a_tol_e_first = max_vertex_tolerance(f, n_e);
                if let Some(es) = f.ds().shape(n_e).cloned() {
                    let e_first = Edge(es);
                    let a_tm = 0.5 * (a_pb.pave1().parameter() + a_pb.pave2().parameter());
                    if let Some(p) = point_on_edge(&e_first, a_tm) {
                        a_pm_first = p;
                    } else if let Ok(p) = AlgoTools::point_on_edge(&e_first, a_tm) {
                        a_pm_first = p;
                    }
                }
                continue;
            }
            let n_e = a_pb.original_edge();
            let a_tol_e = max_vertex_tolerance(f, n_e);
            let Some(es) = f.ds().shape(n_e).cloned() else {
                rest.push(a_pb);
                continue;
            };
            let a_e = Edge(es);
            let (i_err, a_t_out, _dist) =
                ctx.compute_pe_pnt(&a_pm_first, a_tol_e_first + a_tol_e + fuzzy, &a_e);
            let t1 = a_pb.pave1().parameter();
            let t2 = a_pb.pave2().parameter();
            if i_err == 0 && a_t_out > t1 && a_t_out < t2 {
                a_lpb_cb.push(a_pb);
            } else {
                rest.push(a_pb);
            }
        }
        make_new_common_block(f, &a_lpb_cb, faces);
        *a_lpb = rest;
    }
}

/// Fence helper used when comparing newly written common-block members.
#[allow(dead_code)]
fn pb_seen(seen: &mut HashSet<crate::pave_ff_exist::PbKey>, pb: &BopdsPaveBlock) -> bool {
    !seen.insert(pb_key(pb))
}
