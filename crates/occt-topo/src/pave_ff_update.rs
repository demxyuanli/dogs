//! `UpdateExistingPaveBlocks`, `UpdateFaceInfo`, and `UpdatePaveBlocks`.
//!
//! Source: `BOPAlgo_PaveFiller_6.cxx` (`UpdateFaceInfo` at 1673,
//! `UpdateExistingPaveBlocks` at 3278, `UpdatePaveBlocks` at 3679).

use std::collections::{HashMap, HashSet};

use crate::abs::ShapeType;
use crate::bopds::{BopdsCommonBlock, BopdsPave, BopdsPaveBlock};
use crate::brep_tool::BRepTool;
use crate::edge_face::EdgeFace;
use crate::inttools_data::CommonPartType;
use crate::pave_blocks::make_split_edge;
use crate::pave_common::fill_shrunk_data_for_block;
use crate::pave_ff_exist::{pb_key, PbKey};
use crate::pave_filler::PaveFiller;
use crate::pave_intersect::vertex_on_edge;
use crate::shape::{Edge, Face, Vertex};

/// Map of an existing pave block to the post-treat replacements.
pub(crate) type DmExEdges = HashMap<PbKey, Vec<BopdsPaveBlock>>;

/// Map of a pave block to faces it should be projected onto as IN.
pub(crate) type PbFacesMap = HashMap<PbKey, Vec<usize>>;

fn pb_eq(a: &BopdsPaveBlock, b: &BopdsPaveBlock) -> bool {
    pb_key(a) == pb_key(b) && a.original_edge() == b.original_edge()
}

fn tuple_of(pb: &BopdsPaveBlock) -> (usize, f64, f64) {
    (pb.edge(), pb.first, pb.last)
}

fn tuples_eq(a: (usize, f64, f64), b: (usize, f64, f64)) -> bool {
    a.0 == b.0 && (a.1 - b.1).abs() <= 1e-7 && (a.2 - b.2).abs() <= 1e-7
}

fn replace_vert(list: &mut Vec<usize>, old: usize, new: usize) {
    let mut hit = false;
    list.retain(|&v| {
        if v == old {
            hit = true;
            false
        } else {
            true
        }
    });
    if hit && !list.contains(&new) {
        list.push(new);
    }
}

/// `BOPAlgo_PaveFiller::UpdateExistingPaveBlocks` (`_6.cxx:3278`).
pub(crate) fn update_existing_pave_blocks(
    f: &mut PaveFiller,
    pbf: &BopdsPaveBlock,
    lpb: &mut Vec<BopdsPaveBlock>,
    pb_faces: &PbFacesMap,
) {
    if lpb.is_empty() {
        return;
    }
    let b_cb = f.ds().is_common_block(pbf);
    let mut lpb1: Vec<BopdsPaveBlock> = Vec::new();
    let mut a_faces: Vec<usize> = Vec::new();
    if b_cb {
        if let Some(cb) = f.ds().common_block(pbf) {
            a_faces = cb.faces().to_vec();
            if !cb.pave_blocks().is_empty() {
                lpb1 = cb.pave_blocks().to_vec();
            } else {
                lpb1.push(pbf.clone());
            }
        }
    } else {
        lpb1.push(pbf.clone());
    }

    for pb1 in &lpb1 {
        let n_e = pb1.original_edge();
        if n_e == 0 {
            continue;
        }
        let list = f.ds_mut().change_pave_blocks_mut(n_e);
        list.retain(|pb2| !pb_eq(pb1, pb2));
    }

    if b_cb {
        let mut lpb_new: Vec<BopdsPaveBlock> = Vec::new();
        let fuzzy = f.fuzzy_value();
        for pb_value in lpb.iter() {
            let pb_value_paves = [pb_value.pave1(), pb_value.pave2()];
            let mut a_cb = BopdsCommonBlock::new();
            for pb2 in &lpb1 {
                let n_e = pb2.original_edge();
                let mut pb2n = BopdsPaveBlock::new();
                if pb_value.original_edge() == n_e {
                    pb2n.set_pave1(pb_value_paves[0]);
                    pb2n.set_pave2(pb_value_paves[1]);
                } else {
                    let mut pave = [BopdsPave::new(0, 0.0), BopdsPave::new(0, 0.0)];
                    if pb_value_paves[0].index() == pb_value_paves[1].index()
                        && pb2.pave1().index() == pb2.pave2().index()
                    {
                        pave[0] = BopdsPave::new(
                            pb_value_paves[0].index(),
                            pb2.pave1().parameter(),
                        );
                        pave[1] = BopdsPave::new(
                            pb_value_paves[1].index(),
                            pb2.pave2().parameter(),
                        );
                    } else {
                        for i in 0..2 {
                            let n_v = pb_value_paves[i].index();
                            pave[i].set_index(n_v);
                            if n_v == pb2.pave1().index() {
                                pave[i].set_parameter(pb2.pave1().parameter());
                            } else if n_v == pb2.pave2().index() {
                                pave[i].set_parameter(pb2.pave2().parameter());
                            } else {
                                let param = project_vertex_on_original(f, n_v, n_e, pb2, fuzzy);
                                pave[i].set_parameter(param);
                            }
                        }
                    }
                    if pave[1].parameter() < pave[0].parameter() {
                        pave.swap(0, 1);
                    }
                    pb2n.set_pave1(pave[0]);
                    pb2n.set_pave2(pave[1]);
                }
                pb2n.set_edge(pb_value.edge());
                pb2n.set_original_edge(n_e);
                a_cb.add_pave_block(pb2n.clone());
                f.ds_mut().change_pave_blocks_mut(n_e).push(pb2n.clone());
                f.ds_mut().set_common_block(&pb2n, a_cb.clone());
            }
            a_cb.set_faces(a_faces.clone());
            if let Some(pb_new) = a_cb.pave_block1().cloned() {
                lpb_new.push(pb_new);
            }
            f.ds_mut().push_common_block(a_cb);
        }
        *lpb = lpb_new;
    } else {
        let n_e = pbf.original_edge();
        if n_e != 0 {
            for pb in lpb.iter() {
                f.ds_mut().change_pave_blocks_mut(n_e).push(pb.clone());
            }
        }
    }

    let Some(l_faces) = pb_faces.get(&pb_key(pbf)).cloned() else {
        return;
    };
    let fuzzy = f.fuzzy_value();
    for n_f in l_faces {
        f.ds_mut().ensure_face_info(n_f);
        let Some(f_shape) = f.ds().shape(n_f).cloned() else {
            continue;
        };
        let a_f = Face(f_shape);
        for pb in lpb.iter() {
            let on_or_in = f
                .ds()
                .face_info(n_f)
                .map(|fi| {
                    let t = tuple_of(pb);
                    fi.paves_on().iter().copied().any(|x| tuples_eq(x, t))
                        || fi.paves_in().iter().copied().any(|x| tuples_eq(x, t))
                })
                .unwrap_or(false);
            if on_or_in {
                continue;
            }
            let Some(e_shape) = f.ds().shape(pb.edge()).cloned() else {
                continue;
            };
            let mut an_ef = EdgeFace::new();
            an_ef.set_edge(Edge(e_shape));
            an_ef.set_face(a_f.clone());
            an_ef.set_fuzzy_value(fuzzy);
            an_ef.set_range(pb.pave1().parameter(), pb.pave2().parameter());
            if an_ef.perform().is_err() {
                continue;
            }
            let cps = an_ef.common_parts();
            let b_coincide = cps.len() == 1 && cps[0].part_type() == CommonPartType::Edge;
            if !b_coincide {
                continue;
            }
            if f.ds().is_common_block(pb) {
                if let Some(mut cb) = f.ds().common_block(pb).cloned() {
                    cb.add_face(n_f);
                    f.ds_mut().set_common_block(pb, cb);
                }
            } else {
                let mut cb = BopdsCommonBlock::new();
                cb.add_pave_block(pb.clone());
                cb.add_face(n_f);
                f.ds_mut().set_common_block(pb, cb);
            }
            if let Some(fi) = f.ds_mut().face_info_mut(n_f) {
                fi.add_pave_in(pb.edge(), pb.first, pb.last);
            }
        }
    }
}

fn project_vertex_on_original(
    f: &PaveFiller,
    n_v: usize,
    n_e: usize,
    pb2: &BopdsPaveBlock,
    fuzzy: f64,
) -> f64 {
    let Some(v_shape) = f.ds().shape(n_v).cloned() else {
        return pb2.pave1().parameter();
    };
    let Some(e_shape) = f.ds().shape(n_e).cloned() else {
        return pb2.pave1().parameter();
    };
    if let Ok(Some((t, _))) = vertex_on_edge(&f.context(), &Vertex(v_shape.clone()), &Edge(e_shape), fuzzy)
    {
        return t;
    }
    let Some(v1s) = f.ds().shape(pb2.pave1().index()) else {
        return pb2.pave1().parameter();
    };
    let Some(v2s) = f.ds().shape(pb2.pave2().index()) else {
        return pb2.pave2().parameter();
    };
    let p = BRepTool::vertex_point(&Vertex(v_shape));
    let p1 = BRepTool::vertex_point(&Vertex(v1s.clone()));
    let p2 = BRepTool::vertex_point(&Vertex(v2s.clone()));
    if p.square_distance(&p1) < p.square_distance(&p2) {
        pb2.pave1().parameter()
    } else {
        pb2.pave2().parameter()
    }
}

/// `BOPAlgo_PaveFiller::UpdateFaceInfo` (`_6.cxx:1673`).
pub(crate) fn update_face_info_ff(
    f: &mut PaveFiller,
    dm_ex: &mut DmExEdges,
    dm_new_sd: &HashMap<usize, usize>,
    pb_faces: &PbFacesMap,
) {
    let nb_ff = f.ds().interf_ff().len();
    let mut mf: HashSet<usize> = HashSet::new();
    let mut an_edge_lpb: HashMap<usize, Vec<BopdsPaveBlock>> = HashMap::new();

    for i in 0..nb_ff {
        let (n_f1, n_f2) = f.ds().interf_ff()[i].indices();
        f.ds_mut().ensure_face_info(n_f1);
        f.ds_mut().ensure_face_info(n_f2);
        let nb_c = f.ds().interf_ff()[i].curves().len();
        for j in 0..nb_c {
            let mut k = 0;
            while k < f.ds().interf_ff()[i].curves()[j].pave_blocks().len() {
                let pb = f.ds().interf_ff()[i].curves()[j].pave_blocks()[k].clone();
                if let Some(lpb) = dm_ex.get(&pb_key(&pb)).cloned() {
                    let mut lpb = lpb;
                    update_existing_pave_blocks(f, &pb, &mut lpb, pb_faces);
                    for it in &lpb {
                        let n_e = it.edge();
                        an_edge_lpb.entry(n_e).or_default().push(it.clone());
                    }
                    f.ds_mut()
                        .interf_ff_mut()[i]
                        .change_curves()[j]
                        .change_pave_blocks()
                        .remove(k);
                    continue;
                }
                if pb.edge() != 0 {
                    if let Some(fi) = f.ds_mut().face_info_mut(n_f1) {
                        fi.add_pave(pb.edge(), pb.first, pb.last);
                    }
                    if let Some(fi) = f.ds_mut().face_info_mut(n_f2) {
                        fi.add_pave(pb.edge(), pb.first, pb.last);
                    }
                    an_edge_lpb
                        .entry(pb.edge())
                        .or_default()
                        .push(pb.clone());
                }
                k += 1;
            }
        }
        let nb_p = f.ds().interf_ff()[i].points().len();
        for j in 0..nb_p {
            let Some(n_v1) = f.ds().interf_ff()[i].points()[j].index() else {
                continue;
            };
            if let Some(fi) = f.ds_mut().face_info_mut(n_f1) {
                fi.add_vert_sc(n_v1);
                fi.add_vert(n_v1, 0.0, 0.0);
            }
            if let Some(fi) = f.ds_mut().face_info_mut(n_f2) {
                fi.add_vert_sc(n_v1);
                fi.add_vert(n_v1, 0.0, 0.0);
            }
        }
        mf.insert(n_f1);
        mf.insert(n_f2);
    }

    let mut b_new_cb = false;
    let edge_keys: Vec<usize> = an_edge_lpb.keys().copied().collect();
    for n_e in edge_keys {
        let lpb = match an_edge_lpb.get(&n_e).cloned() {
            Some(l) if l.len() > 1 => l,
            _ => continue,
        };
        b_new_cb = true;
        let mut m_faces: HashSet<usize> = HashSet::new();
        let mut m_pbs: Vec<BopdsPaveBlock> = Vec::new();
        let mut a_cb: Option<BopdsCommonBlock> = None;
        for pb in &lpb {
            if !m_pbs.iter().any(|p| pb_eq(p, pb)) {
                m_pbs.push(pb.clone());
            }
            if let Some(cb) = f.ds().common_block(pb) {
                for p in cb.pave_blocks() {
                    if !m_pbs.iter().any(|x| pb_eq(x, p)) {
                        m_pbs.push(p.clone());
                    }
                }
                for &face in cb.faces() {
                    m_faces.insert(face);
                }
                if a_cb.is_none() {
                    a_cb = Some(cb.clone());
                }
            }
        }
        if a_cb.is_none() {
            let mut cb = BopdsCommonBlock::new();
            cb.set_pave_blocks(lpb.clone());
            for pb in &lpb {
                f.ds_mut().set_common_block(pb, cb.clone());
            }
            f.ds_mut().push_common_block(cb);
        } else if let Some(mut cb) = a_cb {
            let mut lpb_new = Vec::new();
            for pb in &m_pbs {
                f.ds_mut().set_common_block(pb, cb.clone());
                lpb_new.push(pb.clone());
            }
            cb.set_pave_blocks(lpb_new);
            cb.set_faces(m_faces.into_iter().collect());
            f.ds_mut().set_common_block(&lpb[0], cb);
        }
    }

    let b_verts = !dm_new_sd.is_empty();
    let b_edges = !dm_ex.is_empty() || b_new_cb;
    if !b_verts && !b_edges {
        return;
    }

    let faces: Vec<usize> = mf.iter().copied().collect();
    for n_f1 in faces {
        if b_verts {
            if let Some(fi) = f.ds_mut().face_info_mut(n_f1) {
                for (&n_v1, &n_v2) in dm_new_sd {
                    replace_vert(&mut fi.verts_on, n_v1, n_v2);
                    replace_vert(&mut fi.verts_in, n_v1, n_v2);
                }
            }
        }
        if b_edges {
            replace_face_paves(f, n_f1, dm_ex);
        }
    }
}

fn replace_face_paves(f: &mut PaveFiller, n_f: usize, dm_ex: &DmExEdges) {
    let Some(fi) = f.ds().face_info(n_f).cloned() else {
        return;
    };
    let mut fence: HashSet<PbKey> = HashSet::new();
    let groups = [fi.paves_on().to_vec(), fi.paves_in().to_vec(), fi.paves().to_vec()];
    let mut rebuilt = [Vec::new(), Vec::new(), Vec::new()];
    for (g, src) in groups.iter().enumerate() {
        for &t in src {
            let replacements: Vec<(usize, f64, f64)> = dm_ex
                .iter()
                .find(|(k, _)| k.0 == t.0 && k.1 == t.1.to_bits() && k.2 == t.2.to_bits())
                .map(|(_, lpb)| lpb.iter().map(tuple_of).collect())
                .unwrap_or_default();
            if !replacements.is_empty() {
                for r in replacements {
                    let key = (r.0, r.1.to_bits(), r.2.to_bits());
                    if fence.insert(key) {
                        rebuilt[g].push(r);
                    }
                }
            } else {
                let key = (t.0, t.1.to_bits(), t.2.to_bits());
                if fence.insert(key) {
                    rebuilt[g].push(t);
                }
            }
        }
    }
    if let Some(fi) = f.ds_mut().face_info_mut(n_f) {
        fi.paves_on.clear();
        fi.paves_in.clear();
        fi.paves.clear();
        for t in &rebuilt[0] {
            fi.add_pave_on(t.0, t.1, t.2);
        }
        for t in &rebuilt[1] {
            fi.add_pave_in(t.0, t.1, t.2);
        }
        for t in &rebuilt[2] {
            fi.add_pave(t.0, t.1, t.2);
        }
    }
}

/// `BOPAlgo_PaveFiller::UpdatePaveBlocks` (`_6.cxx:3679`).
pub(crate) fn update_pave_blocks_ff(
    f: &mut PaveFiller,
    dm_new_sd: &HashMap<usize, usize>,
) -> Result<(), String> {
    if dm_new_sd.is_empty() {
        return Ok(());
    }
    let mut an_all: Vec<BopdsPaveBlock> = Vec::new();
    for ff in f.ds().interf_ff() {
        for nc in ff.curves() {
            an_all.extend(nc.pave_blocks().iter().cloned());
        }
    }
    for list in f.ds().pave_blocks_pool() {
        an_all.extend(list.iter().cloned());
    }

    let mut mpb: HashSet<PbKey> = HashSet::new();
    let mut micro: HashSet<usize> = HashSet::new();
    for mut a_pb in an_all {
        if f.ds().is_common_block(&a_pb) {
            if let Some(cb) = f.ds().common_block(&a_pb) {
                if let Some(first) = cb.pave_block1() {
                    a_pb = first.clone();
                }
            }
        }
        if !mpb.insert(pb_key(&a_pb)) {
            continue;
        }
        let (mut n_v1, mut n_v2) = a_pb.indices();
        let (t1, t2) = a_pb.range();
        let was_regular = n_v1 != n_v2;
        let mut b_rebuild = false;
        if let Some(&n) = dm_new_sd.get(&n_v1) {
            n_v1 = n;
            a_pb.set_pave1(BopdsPave::new(n_v1, t1));
            b_rebuild = true;
        }
        if let Some(&n) = dm_new_sd.get(&n_v2) {
            n_v2 = n;
            a_pb.set_pave2(BopdsPave::new(n_v2, t2));
            b_rebuild = true;
        }
        if !b_rebuild {
            continue;
        }
        let mut n_e = a_pb.edge();
        if n_e == 0 {
            n_e = a_pb.original_edge();
        }
        let is_deg = f
            .ds()
            .shape(n_e)
            .map(|s| BRepTool::is_degenerated(&Edge(s.clone())))
            .unwrap_or(false)
            || f.ds().shape_info(n_e).and_then(|s| s.flag()).is_some()
                && f
                    .ds()
                    .shape_info(n_e)
                    .map(|s| s.shape_type() == ShapeType::Edge)
                    .unwrap_or(false)
                    && f
                        .ds()
                        .shape(n_e)
                        .map(|s| BRepTool::is_degenerated(&Edge(s.clone())))
                        .unwrap_or(false);
        if was_regular && !is_deg && n_v1 == n_v2 {
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
        let n_sp = make_split_edge(f, n_e, n_v1, t1, n_v2, t2)?;
        if f.ds().is_common_block(&a_pb) {
            if let Some(mut cb) = f.ds().common_block(&a_pb).cloned() {
                cb.set_edge(n_sp);
                f.ds_mut().set_common_block(&a_pb, cb);
            }
        } else {
            a_pb.set_edge(n_sp);
            let orig = a_pb.original_edge();
            if orig != 0 {
                let list = f.ds_mut().change_pave_blocks_mut(orig);
                if let Some(slot) = list.iter_mut().find(|p| pb_eq(p, &a_pb) || p.edge() == n_e) {
                    slot.set_edge(n_sp);
                    slot.set_pave1(BopdsPave::new(n_v1, t1));
                    slot.set_pave2(BopdsPave::new(n_v2, t2));
                }
            }
            for ff in f.ds_mut().interf_ff_mut() {
                for nc in ff.change_curves() {
                    for pb in nc.change_pave_blocks() {
                        if pb_eq(pb, &a_pb) || (pb.edge() == n_e && (pb.first - t1).abs() <= 1e-7) {
                            pb.set_edge(n_sp);
                            pb.set_pave1(BopdsPave::new(n_v1, t1));
                            pb.set_pave2(BopdsPave::new(n_v2, t2));
                        }
                    }
                }
            }
        }
    }
    if !micro.is_empty() {
        f.ds_mut().remove_pave_blocks(&micro);
    }
    Ok(())
}
