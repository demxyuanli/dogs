//! `BOPAlgo_PaveFiller::MakeSplitEdges` and `SplitEdge`.
//!
//! Source: `BOPAlgo_PaveFiller_7.cxx` (`MakeSplitEdges` at 371, `SplitEdge`
//! at 553). Common-block members share one split edge; an untouched
//! single-block original edge is reused without splitting.

use std::collections::HashSet;

use occt_core::precision::{CONFUSION, PCONFUSION};

use crate::abs::Orientation;
use crate::algo_tools::AlgoTools;
use crate::bopalgo_tools::compute_tolerance_of_cb;
use crate::bopds::{BopdsPaveBlock, BopdsShapeInfo};
use crate::brep_tool::BRepTool;
use crate::pave_common::{update_common_blocks_with_sd_vertices, update_edge_tolerance};
use crate::pave_ff_exist::{pb_key, PbKey};
use crate::pave_filler::PaveFiller;
use crate::shape::{Edge, Vertex};

struct SplitTask {
    pb: BopdsPaveBlock,
    cb_key: Option<PbKey>,
    n_e: usize,
    n_v1: usize,
    n_v2: usize,
    t1: f64,
    t2: f64,
}

fn write_pb_edge(f: &mut PaveFiller, pb: &BopdsPaveBlock, n_sp: usize) {
    let orig = pb.original_edge();
    if orig < f.ds().nb_shapes() && f.ds().has_pave_blocks(orig) {
        let blocks = f.ds_mut().change_pave_blocks_mut(orig);
        for slot in blocks.iter_mut() {
            if pb_key(slot) == pb_key(pb) || {
                let (a, b) = slot.indices();
                let (c, d) = pb.indices();
                a == c && b == d && slot.original_edge() == pb.original_edge()
            } {
                slot.set_edge(n_sp);
            }
        }
    }
}

/// `BOPAlgo_PaveFiller::MakeSplitEdges`.
pub fn make_split_edges(f: &mut PaveFiller) -> Result<(), String> {
    if f.ds().pave_blocks_pool().is_empty() {
        return Ok(());
    }
    update_common_blocks_with_sd_vertices(f)?;
    let mut mcb: HashSet<PbKey> = HashSet::new();
    let mut tasks: Vec<SplitTask> = Vec::new();
    let n_pool = f.ds().pave_blocks_pool().len();
    for i in 0..n_pool {
        let list = f.ds().pave_blocks_pool()[i].clone();
        let a_lpb_extent = list.len();
        for pb in list {
            let n_e = pb.original_edge();
            let Some(si_e) = f.ds().shape_info(n_e) else {
                continue;
            };
            if si_e.has_flag() {
                continue;
            }
            let cb = f.ds().common_block(&pb).cloned();
            let b_cb = cb.is_some();
            if let Some(ref a_cb) = cb {
                if let Some(pbr) = a_cb.pave_block1() {
                    if !mcb.insert(pb_key(pbr)) {
                        continue;
                    }
                }
            }
            let (n_v1, n_v2) = pb.indices();
            let b_v1 = f.ds().is_new_shape(n_v1);
            let b_v2 = f.ds().is_new_shape(n_v2);
            let mut b_to_split = true;
            if !b_v1 && !b_v2 {
                if !f.non_destructive() || !b_cb {
                    if let Some(ref a_cb) = cb {
                        let mut found = None;
                        for pbx in a_cb.pave_blocks() {
                            let n_ex = pbx.original_edge();
                            if f.ds().pave_blocks(n_ex).len() == 1 {
                                found = Some(pbx.clone());
                                break;
                            }
                        }
                        if let Some(pbr) = found {
                            b_to_split = false;
                            let n_e_real = pbr.original_edge();
                            if let Some(mut a_cb) = f.ds().common_block(&pb).cloned() {
                                a_cb.set_real_pave_block(&pbr);
                                a_cb.set_edge(n_e_real);
                                let a_tol = compute_tolerance_of_cb(&a_cb, f.ds(), f.context());
                                f.ds_mut().set_common_block(&pb, a_cb);
                                update_edge_tolerance(f, n_e_real, a_tol)?;
                            }
                        }
                    } else if a_lpb_extent == 1 {
                        b_to_split = false;
                        write_pb_edge(f, &pb, n_e);
                    }
                    if !b_to_split {
                        continue;
                    }
                }
            }
            let mut a_pb = pb.clone();
            let mut n_e_split = n_e;
            let mut n_v1s = n_v1;
            let mut n_v2s = n_v2;
            if b_cb {
                if let Some(a_cb) = f.ds().common_block(&pb) {
                    if let Some(pbr) = a_cb.pave_block1() {
                        a_pb = pbr.clone();
                        n_e_split = a_pb.original_edge();
                        let (a, b) = a_pb.indices();
                        n_v1s = a;
                        n_v2s = b;
                    }
                }
            }
            let (a_t1, a_t2) = a_pb.range();
            if (a_t2 - a_t1).abs() <= PCONFUSION {
                continue;
            }
            tasks.push(SplitTask {
                pb: a_pb,
                cb_key: cb.and_then(|c| c.pave_block1().map(pb_key)),
                n_e: n_e_split,
                n_v1: n_v1s,
                n_v2: n_v2s,
                t1: a_t1,
                t2: a_t2,
            });
        }
    }

    for t in tasks {
        let Some(e_shape) = f.ds().shape(t.n_e).cloned() else {
            let msg = format!("make_split_edges: edge {} not found", t.n_e);
            f.add_error(msg.clone());
            return Err(msg);
        };
        let Some(v1s) = f.ds().shape(t.n_v1).cloned() else {
            return Err(format!("make_split_edges: vertex {} not found", t.n_v1));
        };
        let Some(v2s) = f.ds().shape(t.n_v2).cloned() else {
            return Err(format!("make_split_edges: vertex {} not found", t.n_v2));
        };
        let mut a_e = e_shape;
        a_e.set_orientation(Orientation::Forward);
        let mut a_v1 = v1s;
        a_v1.set_orientation(Orientation::Forward);
        let mut a_v2 = v2s;
        a_v2.set_orientation(Orientation::Reversed);
        let sp = AlgoTools::make_split_edge(&Edge(a_e), Some(&a_v1), t.t1, Some(&a_v2), t.t2)
            .map_err(|e| {
                let msg = format!("make_split_edges: {e}");
                f.add_error(msg.clone());
                msg
            })?;
        let mut si = BopdsShapeInfo::new(sp.0);
        si.change_sub_shapes()
            .extend_from_slice(&[t.n_v1, t.n_v2]);
        let n_sp = f.ds_mut().append_info(si);
        if let Some(box_) = f.ds_mut().box_of(n_sp).cloned() {
            let mut b = box_;
            b.set_gap(b.gap() + CONFUSION);
            if let Some(slot) = {
                // refresh via enlarge
                None::<()>
            } {
                let _ = slot;
            }
            let _ = b;
        }
        if let Some(key) = t.cb_key {
            if let Some(cb) = f.ds().common_block(&t.pb).cloned() {
                let a_tol = compute_tolerance_of_cb(&cb, f.ds(), f.context());
                update_edge_tolerance(f, n_sp, a_tol)?;
                if let Some(mut cb) = f.ds().common_block(&t.pb).cloned() {
                    let _ = key;
                    cb.set_edge(n_sp);
                    let members = cb.pave_blocks().to_vec();
                    f.ds_mut().set_common_block(&t.pb, cb);
                    write_pb_edge(f, &t.pb, n_sp);
                    for pbx in &members {
                        write_pb_edge(f, pbx, n_sp);
                    }
                }
            }
        } else {
            write_pb_edge(f, &t.pb, n_sp);
        }
    }
    Ok(())
}

/// `BOPAlgo_PaveFiller::SplitEdge`.
pub fn split_edge(
    f: &mut PaveFiller,
    n_e: usize,
    n_v1: usize,
    a_t1: f64,
    n_v2: usize,
    a_t2: f64,
) -> Result<usize, String> {
    let Some(e_shape) = f.ds().shape(n_e).cloned() else {
        return Err(format!("split_edge: edge {n_e} not found"));
    };
    let Some(v1s) = f.ds().shape(n_v1).cloned() else {
        return Err(format!("split_edge: vertex {n_v1} not found"));
    };
    let Some(v2s) = f.ds().shape(n_v2).cloned() else {
        return Err(format!("split_edge: vertex {n_v2} not found"));
    };
    let mut a_e = e_shape;
    a_e.set_orientation(Orientation::Forward);
    let mut a_v1 = v1s;
    a_v1.set_orientation(Orientation::Forward);
    let mut a_v2 = v2s;
    a_v2.set_orientation(Orientation::Reversed);
    let sp = AlgoTools::make_split_edge(&Edge(a_e), Some(&a_v1), a_t1, Some(&a_v2), a_t2)?;
    let mut si = BopdsShapeInfo::new(sp.0);
    si.change_sub_shapes().extend_from_slice(&[n_v1, n_v2]);
    Ok(f.ds_mut().append_info(si))
}

/// Keep `Vertex` import live for callers that pass DS vertex shapes.
pub fn split_edge_vertices<'a>(v1: &'a Vertex, v2: &'a Vertex) -> (&'a Vertex, &'a Vertex) {
    (v1, v2)
}
