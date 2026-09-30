//! `UpdatePaveBlocks` remainder and `RemovePaveBlocks`.
//!
//! Source: `BOPAlgo_PaveFiller_6.cxx` (`UpdatePaveBlocks` at 3679,
//! `RemovePaveBlocks` at 3815).
//!
//! After SD substitution, a pave block whose two vertices collapse onto one
//! and that has no shrunk data is a micro-edge and is collected for
//! `RemovePaveBlocks`. Removal walks the pave-block pool, every F/F section
//! curve, and FaceInfo On/In/Sc tuples.
//! T-97: items below are faithful ports of the named OCCT source, but their
//! OCCT-side consumers are not all ported yet, so parts are not called from this
//! crate. The `dead_code` allowance is deliberate: **pending wiring**, not dead
//! code. Do not delete them to silence warnings (see
//! specs/_a3n00_gap_analysis.md §9.309/§9.310); wire the consumer instead.
#![allow(dead_code)]

use std::collections::HashSet;


use crate::bopds::{BopdsDS, BopdsPave, BopdsPaveBlock};
use crate::brep_tool::BRepTool;
use crate::pave_blocks::make_split_edge;
use crate::pave_common::fill_shrunk_data_for_block;
use crate::pave_ff_exist::{pb_key, PbKey};
use crate::pave_filler::PaveFiller;
use crate::shape::Edge;

fn pb_fence_key(pb: &BopdsPaveBlock) -> PbKey {
    pb_key(pb)
}

/// Collect every pave block on F/F section curves and in the pool.
fn collect_all_pave_blocks(ds: &BopdsDS) -> Vec<BopdsPaveBlock> {
    let mut out = Vec::new();
    for ff in ds.interf_ff() {
        for nc in ff.curves() {
            out.extend(nc.pave_blocks().iter().cloned());
        }
    }
    for list in ds.pave_blocks_pool() {
        out.extend(list.iter().cloned());
    }
    out
}

fn write_edge_on_pb(ds: &mut BopdsDS, pb: &BopdsPaveBlock, n_sp: usize) {
    let orig = pb.original_edge();
    if orig < ds.nb_shapes() && ds.has_pave_blocks(orig) {
        let blocks = ds.change_pave_blocks_mut(orig);
        for slot in blocks.iter_mut() {
            if pb_key(slot) == pb_key(pb) {
                slot.set_edge(n_sp);
            }
        }
    }
    for ff in ds.interf_ff_mut() {
        for nc in ff.change_curves() {
            for slot in nc.change_pave_blocks() {
                if pb_key(slot) == pb_key(pb) {
                    slot.set_edge(n_sp);
                }
            }
        }
    }
}

/// `BOPAlgo_PaveFiller::RemovePaveBlocks` (`_6.cxx:3815`).
pub fn remove_pave_blocks(ds: &mut BopdsDS, the_edges: &HashSet<usize>) {
    ds.remove_pave_blocks(the_edges);
}

/// `BOPAlgo_PaveFiller::UpdatePaveBlocks` (`_6.cxx:3679`).
pub fn update_pave_blocks(f: &mut PaveFiller, dm_new_sd: &std::collections::HashMap<usize, usize>) -> Result<(), String> {
    if dm_new_sd.is_empty() {
        return Ok(());
    }
    let all = collect_all_pave_blocks(f.ds());
    let mut fence: HashSet<PbKey> = HashSet::new();
    let mut micro: HashSet<usize> = HashSet::new();
    let mut rebuilt: Vec<(BopdsPaveBlock, usize, bool)> = Vec::new();
    for mut a_pb in all {
        let b_cb = f.ds().is_common_block(&a_pb);
        if b_cb {
            if let Some(cb) = f.ds().common_block(&a_pb) {
                if let Some(first) = cb.pave_block1() {
                    a_pb = first.clone();
                }
            }
        }
        if !fence.insert(pb_fence_key(&a_pb)) {
            continue;
        }
        let (mut n_v0, mut n_v1) = a_pb.indices();
        let (a_t0, a_t1) = a_pb.range();
        let was_regular = n_v0 != n_v1;
        let mut b_rebuild = false;
        if let Some(&n) = dm_new_sd.get(&n_v0) {
            n_v0 = n;
            a_pb.set_pave1(BopdsPave::new(n_v0, a_t0));
            b_rebuild = true;
        }
        if let Some(&n) = dm_new_sd.get(&n_v1) {
            n_v1 = n;
            a_pb.set_pave2(BopdsPave::new(n_v1, a_t1));
            b_rebuild = true;
        }
        if !b_rebuild {
            continue;
        }
        let mut n_e = a_pb.edge();
        if n_e == 0 {
            n_e = a_pb.original_edge();
        }
        let is_deg = f.ds().shape_info(n_e).map(|s| s.has_flag()).unwrap_or(false)
            || f
                .ds()
                .shape(n_e)
                .map(|s| BRepTool::is_degenerated(&Edge(s.clone())))
                .unwrap_or(false);
        if was_regular && !is_deg && n_v0 == n_v1 {
            if let Some(e_shape) = f.ds().shape(n_e).cloned() {
                let edge = Edge(e_shape);
                let tol = BRepTool::edge_tolerance(&edge);
                fill_shrunk_data_for_block(f, &edge, tol, &mut a_pb);
            }
            if !a_pb.has_shrunk_data() {
                micro.insert(n_e);
                continue;
            }
        }
        rebuilt.push((a_pb, n_e, b_cb));
    }
    for (a_pb, n_e, b_cb) in rebuilt {
        let (n_v0, n_v1) = a_pb.indices();
        let (a_t0, a_t1) = a_pb.range();
        let n_sp = make_split_edge(f, n_e, n_v0, a_t0, n_v1, a_t1)?;
        if b_cb {
            if let Some(mut cb) = f.ds().common_block(&a_pb).cloned() {
                cb.set_edge(n_sp);
                f.ds_mut().set_common_block(&a_pb, cb);
            }
        } else {
            write_edge_on_pb(f.ds_mut(), &a_pb, n_sp);
        }
    }
    if !micro.is_empty() {
        remove_pave_blocks(f.ds_mut(), &micro);
    }
    Ok(())
}
