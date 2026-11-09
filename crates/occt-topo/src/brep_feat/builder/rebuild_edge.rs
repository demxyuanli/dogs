//! `BRepFeat_Builder::RebuildEdge` (`BRepFeat_Builder.cxx:562-722`).
//!
//! A boundary edge of a tool face has images that are partly retained and
//! partly dropped. When dropping some of them splits the edge differently, the
//! kept pieces no longer cover the face boundary; this rebuilds the split edges
//! from the paves of the removed pieces plus the original bound vertices.
//!
//! Every step keeps the numbering of the reference:
//! 1. the original bound vertices (`aMV`, `aMVOr`),
//! 2. the paves of the blocks whose split edge was dropped (`aMI`),
//! 3. the vertices that split the source edge,
//! 4. the ext paves handed to `BOPDS_PaveBlock::Update`,
//! 5. the re-split and `BOPTools_AlgoTools::MakeSplitEdge` for the new pieces.

use std::collections::HashSet;

use crate::abs::Orientation;
use crate::algo_tools::AlgoTools;
use crate::bop_builder2::BopBuilder;
use crate::bopds_pave::BopdsPaveBlock;
use crate::boptools_2d::build_pcurve_for_edge_on_face;
use crate::shape::{Edge, Face, TopoShape};

use super::shape_set::ShapeSet;

/// Rebuilds the images of the edge `the_e` on the face `the_f`, keeping only
/// the parts not covered by `a_me`, and appends them to `a_lim`
/// (`BRepFeat_Builder::RebuildEdge`).
pub fn rebuild_edge(
    b: &mut BopBuilder,
    the_e: &TopoShape,
    the_f: &Face,
    a_me: &ShapeSet,
    a_lim: &mut Vec<TopoShape>,
) -> Result<(), String> {
    let Some(n_e) = b.ds().index(the_e) else {
        return Ok(());
    };
    // 1. collect origin vertices to aMV map.
    let a_ls: Vec<usize> =
        b.ds().shape_info(n_e).map(|si| si.sub_shapes().to_vec()).unwrap_or_default();
    let sd = b.shapes_sd().clone();
    let mut mv: HashSet<usize> = HashSet::new();
    let mut mv_or: HashSet<usize> = HashSet::new();
    for n_v in &a_ls {
        let n_vx = sd.get(n_v).copied().unwrap_or(*n_v);
        mv.insert(n_vx);
        mv_or.insert(n_vx);
    }

    // The pave blocks of the source edge (`myDS->ChangePaveBlocks(nE)`).
    let lpb: Vec<BopdsPaveBlock> = b.ds().pave_blocks(n_e).to_vec();

    // 2. collect vertices that should be removed to aMI map.
    let mut mi: HashSet<usize> = HashSet::new();
    let mut mpb: Vec<BopdsPaveBlock> = Vec::new();
    for pb in &lpb {
        let a_e1 = b.ds().shape(pb.edge()).cloned();
        let in_me = a_e1.as_ref().map(|s| a_me.contains(s)).unwrap_or(false);
        if in_me {
            let (n_v1, n_v2) = pb.indices();
            mi.insert(n_v1);
            mi.insert(n_v2);
        } else {
            mpb.push(pb.clone());
        }
    }
    // 3. collect vertices that split the source shape.
    for pb in &lpb {
        let (n_v1, n_v2) = pb.indices();
        if !mi.contains(&n_v1) {
            mv.insert(n_v1);
        }
        if !mi.contains(&n_v2) {
            mv.insert(n_v2);
        }
    }
    // 4. collect ext paves.
    let mut pb_new = BopdsPaveBlock::new();
    let mut m_add: HashSet<usize> = HashSet::new();
    for pb in &lpb {
        let (n_v1, n_v2) = pb.indices();
        if mv.contains(&n_v1) && (m_add.insert(n_v1) || mv_or.contains(&n_v1)) {
            pb_new.append_ext_pave(pb.pave1());
        }
        if mv.contains(&n_v2) && (m_add.insert(n_v2) || mv_or.contains(&n_v2)) {
            pb_new.append_ext_pave(pb.pave2());
        }
    }

    // 5. split edge by the new set of vertices. `aLPB` is the data-structure
    // list of the source edge: OCCT clears it and refills it with the
    // elementary blocks of the re-split (`BRepFeat_Builder.cxx:670-672`).
    let mut e_fwd = the_e.clone();
    e_fwd.set_orientation(Orientation::Forward);
    a_lim.clear();
    pb_new.set_original_edge(n_e);
    {
        let list = b.ds_mut().change_pave_blocks_mut(n_e);
        list.clear();
        pb_new.update(list, false);
    }
    let lpb_new: Vec<BopdsPaveBlock> = b.ds().pave_blocks(n_e).to_vec();

    // `aPB->SetEdge(nSp)` rewrites the data-structure blocks, so the new edge
    // indices are applied after the walk (the walk itself reads the DS).
    let mut new_edges: Vec<(usize, usize)> = Vec::new();
    for (pos, pb) in lpb_new.iter().enumerate() {
        let (n_v1, t1) = pb.pave1().contents();
        let (n_v2, t2) = pb.pave2().contents();
        // check if it is the old pave block.
        let mut b_old = false;
        for pb1 in &mpb {
            let (n_v11, n_v21) = pb1.indices();
            let (t11, t21) = pb1.range();
            if n_v1 == n_v11 && n_v2 == n_v21 && t1 == t11 && t2 == t21 {
                if let Some(e_im) = b.ds().shape(pb1.edge()).cloned() {
                    a_lim.push(e_im);
                }
                b_old = true;
                break;
            }
        }
        if b_old {
            continue;
        }
        let Some(v1) = b.ds().shape(n_v1).cloned() else { continue };
        let Some(v2) = b.ds().shape(n_v2).cloned() else { continue };
        let mut v1 = v1;
        v1.set_orientation(Orientation::Forward);
        let mut v2 = v2;
        v2.set_orientation(Orientation::Reversed);
        let a_sp = AlgoTools::make_split_edge(
            &Edge(e_fwd.clone()),
            Some(&v1),
            t1,
            Some(&v2),
            t2,
        )?;
        let _ = build_pcurve_for_edge_on_face(&a_sp, the_f);
        let n_sp = b.ds_mut().append(a_sp.0.clone())?;
        new_edges.push((pos, n_sp));
        a_lim.push(a_sp.0.clone());
    }
    for (pos, n_sp) in new_edges {
        if let Some(pb) = b.ds_mut().change_pave_blocks_mut(n_e).get_mut(pos) {
            pb.set_edge(n_sp);
        }
    }
    Ok(())
}
