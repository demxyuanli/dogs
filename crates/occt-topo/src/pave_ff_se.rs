//! `PutSEInOtherFaces`, `RemoveMicroSectionEdges`, `RemoveMicroEdges`.
//!
//! Source: `BOPAlgo_PaveFiller_6.cxx` (`PutSEInOtherFaces` at 4277,
//! `RemoveMicroSectionEdges` at 4308, `RemoveMicroEdges` at 4388) and
//! `BOPTools_AlgoTools::IsMicroEdge` (`BOPTools_AlgoTools.cxx:2067`).

use std::collections::HashSet;

use occt_core::precision::CONFUSION;

use crate::abs::ShapeType;
use crate::bopds::BopdsPaveBlock;
use crate::brep_tool::BRepTool;
use crate::inttools_range::ShrunkRange;
use crate::pave_common::fill_shrunk_data_for_block;
use crate::pave_ff::CoupleOfPaveBlocks;
use crate::pave_ff_exist::pb_key;
use crate::pave_filler::PaveFiller;
use crate::pave_force_ef::force_interf_ef_on;
use crate::shape::{Edge, Face, TopoShape};

fn pb_has_edge(pb: &BopdsPaveBlock) -> bool {
    pb.edge() != 0
}

/// `BOPTools_AlgoTools::IsMicroEdge`.
///
/// Degenerated or non-geometric edges are micro. Otherwise a
/// `IntTools_ShrunkRange` is computed on the edge's own vertex parameters;
/// the edge is micro when that computation is not done. When
/// `check_splittable` is set, a done-but-not-splittable range is also micro.
pub fn is_micro_edge(edge: &Edge, check_splittable: bool) -> bool {
    if BRepTool::is_degenerated(edge) {
        return true;
    }
    if BRepTool::edge_curve(edge).is_none() {
        return true;
    }
    let mut sr = ShrunkRange::new();
    let _ = sr.set_shrunk_range(edge, &Face::new(), CONFUSION);
    if !sr.is_done() {
        return true;
    }
    if check_splittable && !sr.is_splittable() {
        return true;
    }
    false
}

/// `BOPAlgo_PaveFiller::PutSEInOtherFaces` (`_6.cxx:4277`).
///
/// Collect every section pave block from F/F curves and force an edge/face
/// intersection against faces that did not produce that section
/// (`ForceInterfEF(aMPBScAll, false)`).
pub fn put_se_in_other_faces(f: &mut PaveFiller) -> Result<(), String> {
    let mut mpb_sc: Vec<BopdsPaveBlock> = Vec::new();
    let mut fence: HashSet<crate::pave_ff_exist::PbKey> = HashSet::new();
    let a_nb_ff = f.ds().interf_ff().len();
    for i in 0..a_nb_ff {
        let a_nb_c = f.ds().interf_ff()[i].curves().len();
        for j in 0..a_nb_c {
            for pb in f.ds().interf_ff()[i].curves()[j].pave_blocks() {
                if fence.insert(pb_key(pb)) {
                    mpb_sc.push(pb.clone());
                }
            }
        }
    }
    force_interf_ef_on(f, &mpb_sc, false)
}

/// `BOPAlgo_PaveFiller::RemoveMicroSectionEdges` (`_6.cxx:4308`).
///
/// Walks the post-treat map of section shapes. Non-edge entries and section
/// pave blocks that already have a real DS edge are kept. Remaining edges are
/// classified with `IsMicroEdge(..., false)`: a micro edge is dropped from
/// both the map and the originating F/F curve's pave-block list, and its pave
/// block is collected so `PostTreatFF` can unify the bounding vertices.
pub fn remove_micro_section_edges(
    f: &mut PaveFiller,
    the_mscpb: &mut Vec<(TopoShape, CoupleOfPaveBlocks)>,
    the_micro_pb: &mut Vec<BopdsPaveBlock>,
) {
    if the_mscpb.is_empty() {
        return;
    }
    let mut kept: Vec<(TopoShape, CoupleOfPaveBlocks)> = Vec::new();
    let mut drop: Vec<(usize, usize, BopdsPaveBlock)> = Vec::new();
    for (shape, cpb) in the_mscpb.drain(..) {
        if shape.shape_type() != ShapeType::Edge {
            kept.push((shape, cpb));
            continue;
        }
        let Some(pb) = cpb.pb.as_ref() else {
            kept.push((shape, cpb));
            continue;
        };
        if pb_has_edge(pb) {
            kept.push((shape, cpb));
            continue;
        }
        let edge = Edge(shape.clone());
        if !is_micro_edge(&edge, false) {
            kept.push((shape, cpb));
            continue;
        }
        drop.push((cpb.index_interf, cpb.index, pb.clone()));
        the_micro_pb.push(pb.clone());
    }
    for (i_ff, i_c, pb) in drop {
        if let Some(lpbc) = f
            .ds_mut()
            .interf_ff_mut()
            .get_mut(i_ff)
            .and_then(|ff| ff.change_curves().get_mut(i_c))
            .map(|nc| nc.change_pave_blocks())
        {
            lpbc.retain(|p| {
                p.original_edge() != pb.original_edge()
                    || (p.pave1().parameter() - pb.pave1().parameter()).abs() > 1e-12
                    || (p.pave2().parameter() - pb.pave2().parameter()).abs() > 1e-12
            });
        }
    }
    *the_mscpb = kept;
}

/// `BOPAlgo_PaveFiller::RemoveMicroEdges` (`_6.cxx:4388`).
///
/// Every pool list with at least two pave blocks is examined. Flagged
/// (degenerated) source edges are skipped. For each real pave block whose
/// bounding vertices coincide, `FillShrunkData` is run; a block without
/// shrunk data is a micro edge and is later removed from the DS.
pub fn remove_micro_edges(f: &mut PaveFiller) {
    let mut fence: HashSet<crate::pave_ff_exist::PbKey> = HashSet::new();
    let mut micro: HashSet<usize> = HashSet::new();
    let mut jobs: Vec<(usize, BopdsPaveBlock)> = Vec::new();
    {
        let ds = f.ds();
        for list in ds.pave_blocks_pool() {
            if list.len() < 2 {
                continue;
            }
            let orig = list[0].original_edge();
            if ds.shape_info(orig).map(|si| si.has_flag()).unwrap_or(false) {
                continue;
            }
            for pb in list {
                let real = ds.real_pave_block(pb);
                if !fence.insert(pb_key(&real)) {
                    continue;
                }
                let (n1, n2) = real.indices();
                if n1 != n2 {
                    continue;
                }
                jobs.push((if real.edge() != 0 { real.edge() } else { real.original_edge() }, real));
            }
        }
    }
    let tol = f.fuzzy_value();
    for (n_e, mut real) in jobs {
        let Some(shape) = f.ds().shape(n_e).cloned() else {
            continue;
        };
        fill_shrunk_data_for_block(f, &Edge(shape), tol, &mut real);
        if !real.has_shrunk_data() {
            micro.insert(n_e);
        }
    }
    if !micro.is_empty() {
        f.ds_mut().remove_pave_blocks(&micro);
    }
}

/// Raise the two bounding vertices of each micro section pave block so they
/// cover each other, matching the `theMicroPB` loop in `PostTreatFF`
/// (`_6.cxx:1318-1359`). The vertices are then appended to `a_ls` for the
/// nested filler.
pub fn collect_micro_section_vertices(
    f: &mut PaveFiller,
    the_micro_pb: &[BopdsPaveBlock],
    dm_new_sd: &std::collections::HashMap<usize, usize>,
    added_sd: &mut HashSet<usize>,
    a_ls: &mut Vec<TopoShape>,
) {
    for pb in the_micro_pb {
        let (n0, n1) = pb.indices();
        let mut verts: [Option<TopoShape>; 2] = [None, None];
        for (i, n_raw) in [n0, n1].into_iter().enumerate() {
            let n = dm_new_sd.get(&n_raw).copied().unwrap_or(n_raw);
            let Some(s) = f.ds().shape(n).cloned() else {
                continue;
            };
            verts[i] = Some(s.clone());
            if added_sd.insert(n) {
                a_ls.push(s);
            }
        }
        let (Some(v0), Some(v1)) = (verts[0].take(), verts[1].take()) else {
            continue;
        };
        if v0.same_tshape(&v1) {
            continue;
        }
        let p1 = BRepTool::vertex_point(&crate::shape::Vertex(v0.clone()));
        let p2 = BRepTool::vertex_point(&crate::shape::Vertex(v1.clone()));
        let mut a_tol_v1 = BRepTool::vertex_tolerance(&crate::shape::Vertex(v0.clone()));
        let mut a_tol_v2 = BRepTool::vertex_tolerance(&crate::shape::Vertex(v1.clone()));
        let mut a_dist = p1.distance(&p2);
        a_dist -= a_tol_v1 + a_tol_v2;
        if a_dist > 0.0 {
            a_dist *= 0.5;
            a_tol_v1 += a_dist;
            a_tol_v2 += a_dist;
            crate::shape::Vertex(v0).set_tolerance(a_tol_v1);
            crate::shape::Vertex(v1).set_tolerance(a_tol_v2);
        }
    }
}
