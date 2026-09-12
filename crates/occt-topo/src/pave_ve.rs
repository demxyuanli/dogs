//! `BOPAlgo_PaveFiller::PerformVE` and `IntersectVE`
//! (`BOPAlgo_PaveFiller_2.cxx:141/212`).
//!
//! Vertices are grouped per pave block, SD vertices are fenced so the same
//! vertex/edge pair is projected once, then `IntTools_Context::ComputeVE`
//! (`compute_pe_pnt`) writes extra paves onto the block that contains the
//! parameter. `SplitPaveBlocks` follows for the modified edges.

use std::collections::{HashMap, HashSet};

use crate::abs::ShapeType;
use crate::bopds::{BopdsPave, BopdsPaveBlock};
use crate::brep_tool::BRepTool;
use crate::pave_common::{fill_shrunk_data, update_vertex_sd};
use crate::pave_ef::has_interf_shape_sub_shapes;
use crate::pave_filler::PaveFiller;
use crate::pave_ff_exist::pb_key;
use crate::pave_intersect::{add_ext_pave, collect_pairs};
use crate::shape::{Edge, Vertex};

struct VertexEdgeTask {
    n_vsd: usize,
    n_e: usize,
    members: Vec<usize>,
    pb: BopdsPaveBlock,
}

/// `BOPAlgo_PaveFiller::IntersectVE`.
pub fn intersect_ve(
    f: &mut PaveFiller,
    the_ve_pairs: &[(BopdsPaveBlock, Vec<usize>)],
    the_add_interfs: bool,
) -> Result<(), String> {
    if the_ve_pairs.is_empty() {
        return Ok(());
    }
    let mut a_dmvsd: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
    let mut tasks: Vec<VertexEdgeTask> = Vec::new();
    for (a_pb, a_lv) in the_ve_pairs {
        let n_e = a_pb.original_edge();
        let mut a_mvpb: HashSet<usize> = HashSet::new();
        for pb in f.ds().pave_blocks(n_e) {
            let (n1, n2) = pb.indices();
            a_mvpb.insert(n1);
            a_mvpb.insert(n2);
        }
        for &n_v in a_lv {
            let n_vsd = f.ds().has_shape_sd(n_v).unwrap_or(n_v);
            if a_mvpb.contains(&n_vsd) {
                continue;
            }
            if let Some(li) = a_dmvsd.get_mut(&(n_vsd, n_e)) {
                li.push(n_v);
                continue;
            }
            a_dmvsd.insert((n_vsd, n_e), vec![n_v]);
            tasks.push(VertexEdgeTask {
                n_vsd,
                n_e,
                members: vec![n_v],
                pb: a_pb.clone(),
            });
        }
    }
    for t in &mut tasks {
        if let Some(m) = a_dmvsd.get(&(t.n_vsd, t.n_e)) {
            t.members = m.clone();
        }
    }

    let fuzzy = f.fuzzy_value();
    let mut a_m_edges: Vec<usize> = Vec::new();
    for t in tasks {
        let Some(vs) = f.ds().shape(t.n_vsd).cloned() else {
            continue;
        };
        let Some(es) = f.ds().shape(t.n_e).cloned() else {
            continue;
        };
        let v = Vertex(vs);
        let e = Edge(es);
        let (flag, a_t, dist) = f.context().compute_pe_pnt(
            &BRepTool::vertex_point(&v),
            BRepTool::vertex_tolerance(&v) + fuzzy,
            &e,
        );
        if flag != 0 {
            if flag < 0 {
                f.add_warning(format!(
                    "intersect_ve: ComputeVE failed for vertex {} / edge {}",
                    t.n_vsd, t.n_e
                ));
            }
            continue;
        }
        let a_tol_v_new = BRepTool::vertex_tolerance(&v).max(dist);
        let mut n_last = t.n_vsd;
        for n_v in t.members {
            let n_vx = update_vertex_sd(f, n_v, a_tol_v_new)?;
            n_last = n_vx;
            let pbs = f.ds().pave_blocks(t.n_e).to_vec();
            let mut found = None;
            for pb in &pbs {
                let (a_t1, a_t2) = pb.range();
                if a_t > a_t1 && a_t < a_t2 {
                    found = Some(pb.clone());
                    break;
                }
            }
            let Some(mut a_pb) = found else {
                continue;
            };
            a_pb.append_ext_pave(BopdsPave::new(n_vx, a_t));
            crate::pave_ee_aux::write_pave_block(f, &a_pb);
            add_ext_pave(f.ds_mut(), t.n_e, a_t, n_vx);
            if the_add_interfs {
                if f.ds().is_new_shape(n_vx) {
                    f.ds_mut().add_interf_ve(n_v, t.n_e, Some(n_vx));
                } else {
                    f.ds_mut().add_interf_ve(n_v, t.n_e, None);
                }
            }
            if !a_m_edges.contains(&t.n_e) {
                a_m_edges.push(t.n_e);
            }
            let _ = n_last;
            let _ = pb_key(&t.pb);
        }
    }
    if !a_m_edges.is_empty() {
        crate::pave_split_blocks::split_pave_blocks(f, &a_m_edges, the_add_interfs)?;
    }
    Ok(())
}

/// `BOPAlgo_PaveFiller::PerformVE`.
pub fn perform_ve(f: &mut PaveFiller) -> Result<(), String> {
    crate::pave_shrunk::fill_shrunk_data_ve(f)?;
    let pairs = collect_pairs(f.ds(), ShapeType::Vertex, ShapeType::Edge);
    if pairs.is_empty() {
        return Ok(());
    }
    let mut a_mve: HashMap<(usize, u64, u64), (BopdsPaveBlock, Vec<usize>)> = HashMap::new();
    for (n_v, n_e) in pairs {
        let Some(si_e) = f.ds().shape_info(n_e).cloned() else {
            continue;
        };
        if si_e.has_subshape(n_v) {
            continue;
        }
        if si_e.has_flag() {
            continue;
        }
        if f.ds().has_interf_pair(n_v, n_e) {
            continue;
        }
        if has_interf_shape_sub_shapes(f.ds(), n_v, n_e) {
            continue;
        }
        let pbs = f.ds().pave_blocks(n_e);
        if pbs.is_empty() {
            continue;
        }
        let a_pb = pbs[0].clone();
        if !a_pb.is_splittable() {
            continue;
        }
        a_mve
            .entry(pb_key(&a_pb))
            .or_insert_with(|| (a_pb, Vec::new()))
            .1
            .push(n_v);
    }
    let grouped: Vec<(BopdsPaveBlock, Vec<usize>)> =
        a_mve.into_values().collect();
    intersect_ve(f, &grouped, true)
}

/// Restricted-pair VE used by `RepeatIntersection`.
pub fn perform_ve_pairs(f: &mut PaveFiller, pairs: &[(usize, usize)]) -> Result<(), String> {
    fill_shrunk_data(f)?;
    let mut a_mve: HashMap<(usize, u64, u64), (BopdsPaveBlock, Vec<usize>)> = HashMap::new();
    for &(n_v, n_e) in pairs {
        let Some(si_e) = f.ds().shape_info(n_e).cloned() else {
            continue;
        };
        if si_e.has_subshape(n_v) || si_e.has_flag() {
            continue;
        }
        if f.ds().has_interf_pair(n_v, n_e) {
            continue;
        }
        let pbs = f.ds().pave_blocks(n_e);
        if pbs.is_empty() || !pbs[0].is_splittable() {
            continue;
        }
        let a_pb = pbs[0].clone();
        a_mve
            .entry(pb_key(&a_pb))
            .or_insert_with(|| (a_pb, Vec::new()))
            .1
            .push(n_v);
    }
    let grouped: Vec<(BopdsPaveBlock, Vec<usize>)> = a_mve.into_values().collect();
    intersect_ve(f, &grouped, true)
}

/// `MakeNewCommonBlock` (`BOPAlgo_PaveFiller_2.cxx:401`).
pub fn make_new_common_block(
    f: &mut PaveFiller,
    the_lpb: &[BopdsPaveBlock],
    the_l_faces: &[usize],
) {
    if the_lpb.is_empty() {
        return;
    }
    let mut a_cb_new = crate::bopds::BopdsCommonBlock::new();
    a_cb_new.set_pave_blocks(the_lpb.to_vec());
    a_cb_new.set_faces(the_l_faces.to_vec());
    for a_pb in the_lpb {
        f.ds_mut().set_common_block(a_pb, a_cb_new.clone());
    }
}

