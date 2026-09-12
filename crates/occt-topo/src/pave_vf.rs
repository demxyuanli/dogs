//! `BOPAlgo_PaveFiller::PerformVF` and `TreatVerticesEE`
//! (`BOPAlgo_PaveFiller_4.cxx:139/305`).
//!
//! Glue-full mode only initializes FaceInfo. Otherwise each interfering
//! vertex/face pair is classified with `IntTools_Context::ComputeVF`. Same-
//! domain vertices that share a face are grouped so the expensive projection
//! runs once. `TreatVerticesEE` then classifies new EE vertices against faces
//! they do not already lie ON.

use std::collections::{HashMap, HashSet};

use crate::abs::ShapeType;
use crate::pave_common::update_vertex_sd;
use crate::pave_ef::{compute_vf_uv, has_interf_shape_sub_shapes};
use crate::pave_filler::{GlueEnum, PaveFiller};
use crate::pave_intersect::collect_pairs;
use crate::shape::{Face, Vertex};

/// `BOPAlgo_VertexFace` solver record.
struct VertexFaceTask {
    n_vx: usize,
    n_f: usize,
    members: Vec<usize>,
}

/// `BOPAlgo_PaveFiller::TreatVerticesEE`.
pub fn treat_vertices_ee(f: &mut PaveFiller) -> Result<(), String> {
    let mut a_mi: HashSet<usize> = HashSet::new();
    let mut a_liv: Vec<usize> = Vec::new();
    for it in f.ds().interf_ee() {
        if let Some(n_v) = it.get_index_new() {
            if a_mi.insert(n_v) {
                a_liv.push(n_v);
            }
        }
    }
    if a_liv.is_empty() {
        return Ok(());
    }
    let a_nb_s = f.ds().nb_source_shapes();
    let mut a_lif: Vec<usize> = Vec::new();
    for n_f in 0..a_nb_s {
        let Some(si) = f.ds().shape_info(n_f) else {
            continue;
        };
        if si.shape_type() == ShapeType::Face {
            a_lif.push(n_f);
        }
    }
    if a_lif.is_empty() {
        return Ok(());
    }
    let fuzzy = f.fuzzy_value();
    for &n_f in &a_lif {
        f.ds_mut().ensure_face_info(n_f);
        let on: HashSet<usize> = f
            .ds()
            .face_info(n_f)
            .map(|fi| fi.verts_on().iter().copied().collect())
            .unwrap_or_default();
        for &n_v in &a_liv {
            if on.contains(&n_v) {
                continue;
            }
            let Some(vs) = f.ds().shape(n_v).cloned() else {
                continue;
            };
            let Some(fs) = f.ds().shape(n_f).cloned() else {
                continue;
            };
            let (flag, u, vv, _dummy) = {
                let ctx = f.context_mut();
                compute_vf_uv(ctx, &Vertex(vs), &Face(fs), fuzzy)
            };
            if flag != 0 {
                continue;
            }
            f.ds_mut().add_interf_vf(n_v, n_f, None);
            if let Some(fi) = f.ds_mut().face_info_mut(n_f) {
                fi.add_vert_in(n_v);
                fi.add_vert(n_v, u, vv);
            }
        }
    }
    Ok(())
}

/// `BOPAlgo_PaveFiller::PerformVF`.
pub fn perform_vf(f: &mut PaveFiller) -> Result<(), String> {
    let pairs = collect_pairs(f.ds(), ShapeType::Vertex, ShapeType::Face);
    if f.glue() == GlueEnum::Full {
        for &(_n_v, n_f) in &pairs {
            let is_sub = f
                .ds()
                .shape_info(n_f)
                .map(|si| si.has_subshape(_n_v))
                .unwrap_or(false);
            if !is_sub {
                f.ds_mut().ensure_face_info(n_f);
            }
        }
        return Ok(());
    }
    if pairs.is_empty() {
        return treat_vertices_ee(f);
    }

    let mut a_mvf_pairs: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
    let mut tasks: Vec<VertexFaceTask> = Vec::new();
    for (n_v, n_f) in pairs {
        if f.ds().shape_info(n_f).map(|si| si.has_subshape(n_v)).unwrap_or(false) {
            continue;
        }
        if f.ds().has_interf_pair(n_v, n_f) {
            continue;
        }
        f.ds_mut().ensure_face_info(n_f);
        if has_interf_shape_sub_shapes(f.ds(), n_v, n_f) {
            continue;
        }
        let n_vx = f.ds().has_shape_sd(n_v).unwrap_or(n_v);
        if let Some(members) = a_mvf_pairs.get_mut(&(n_vx, n_f)) {
            members.push(n_v);
            continue;
        }
        a_mvf_pairs.insert((n_vx, n_f), vec![n_v]);
        tasks.push(VertexFaceTask {
            n_vx,
            n_f,
            members: Vec::new(),
        });
    }
    for t in &mut tasks {
        if let Some(m) = a_mvf_pairs.get(&(t.n_vx, t.n_f)) {
            t.members = m.clone();
        }
    }

    let fuzzy = f.fuzzy_value();
    for t in tasks {
        let Some(vs) = f.ds().shape(t.n_vx).cloned() else {
            continue;
        };
        let Some(fs) = f.ds().shape(t.n_f).cloned() else {
            continue;
        };
        let (flag, a_t1, a_t2, a_tol_v_new) = {
            let ctx = f.context_mut();
            compute_vf_uv(ctx, &Vertex(vs), &Face(fs), fuzzy)
        };
        if flag != 0 {
            if flag < 0 {
                f.add_warning(format!(
                    "perform_vf: ComputeVF failed for vertex {} / face {}",
                    t.n_vx, t.n_f
                ));
            }
            continue;
        }
        let mut n_last = t.n_vx;
        for n_v in t.members {
            let n_vx = update_vertex_sd(f, n_v, a_tol_v_new)?;
            n_last = n_vx;
            if f.ds().is_new_shape(n_vx) {
                f.ds_mut().add_interf_vf(n_v, t.n_f, Some(n_vx));
            } else {
                f.ds_mut().add_interf_vf(n_v, t.n_f, None);
            }
            if let Some(fi) = f.ds_mut().face_info_mut(t.n_f) {
                fi.add_vert(n_vx, a_t1, a_t2);
            }
        }
        if let Some(fi) = f.ds_mut().face_info_mut(t.n_f) {
            fi.add_vert_in(n_last);
        }
    }
    treat_vertices_ee(f)
}

/// Restricted-pair VF used by `RepeatIntersection`.
pub fn perform_vf_pairs(f: &mut PaveFiller, pairs: &[(usize, usize)]) -> Result<(), String> {
    if f.glue() == GlueEnum::Full {
        for &(_, n_f) in pairs {
            f.ds_mut().ensure_face_info(n_f);
        }
        return Ok(());
    }
    let fuzzy = f.fuzzy_value();
    let mut hits: Vec<(usize, usize, f64, f64, f64)> = Vec::new();
    for &(n_v, n_f) in pairs {
        if f.ds().shape_info(n_f).map(|si| si.has_subshape(n_v)).unwrap_or(false) {
            continue;
        }
        if f.ds().has_interf_pair(n_v, n_f) {
            continue;
        }
        let n_vsd = f.ds().has_shape_sd(n_v).unwrap_or(n_v);
        let Some(vs) = f.ds().shape(n_vsd).cloned() else {
            continue;
        };
        let Some(fs) = f.ds().shape(n_f).cloned() else {
            continue;
        };
        let (flag, u, vv, tol_new) = {
            let ctx = f.context_mut();
            compute_vf_uv(ctx, &Vertex(vs), &Face(fs), fuzzy)
        };
        if flag == 0 {
            hits.push((n_v, n_f, u, vv, tol_new));
        }
    }
    for (n_v, n_f, u, vv, tol_new) in hits {
        let n_vx = update_vertex_sd(f, n_v, tol_new)?;
        if f.ds().is_new_shape(n_vx) {
            f.ds_mut().add_interf_vf(n_v, n_f, Some(n_vx));
        } else {
            f.ds_mut().add_interf_vf(n_v, n_f, None);
        }
        f.ds_mut().ensure_face_info(n_f);
        if let Some(fi) = f.ds_mut().face_info_mut(n_f) {
            fi.add_vert_in(n_vx);
            fi.add_vert(n_vx, u, vv);
        }
    }
    Ok(())
}
