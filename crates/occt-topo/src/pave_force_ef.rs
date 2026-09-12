//! `BOPAlgo_PaveFiller::ForceInterfEF` (`BOPAlgo_PaveFiller_5.cxx:772/831`).
//!
//! After vertex unification, look for additional EDGE-type common parts
//! between pave blocks and faces that already share both bounding vertices.
//! Matching pairs are intersected with `IntTools_EdgeFace` using a fuzzy
//! value derived from the vertex/face distances, then recorded as IN pave
//! blocks (`FillMap` / `PerformCommonBlocks`).

use std::collections::{HashMap, HashSet};

use occt_core::precision::{CONFUSION, RESOLUTION};

use crate::abs::ShapeType;
use crate::bopalgo_tools::{fill_map_pb_face, perform_common_blocks_faces, PbFaceListMap};
use crate::bbox_from_geometry::shape_bbox;
use crate::bopds::BopdsPaveBlock;
use crate::boptools_2d::intermediate_point;
use crate::brep_surface::face_is_planar;
use crate::brep_tool::BRepTool;
use crate::edge_face::EdgeFace;
use crate::fclass2d::FaceState;
use crate::inttools_data::CommonPartType;
use crate::pave_common::fill_shrunk_data_for_block;
use crate::pave_ee_aux::get_pb_box;
use crate::pave_ff_exist::{pb_key, PbKey};
use crate::pave_filler::PaveFiller;
use crate::shape::{Edge, Face, Vertex};

fn vertex_tol(f: &PaveFiller, n_v: usize) -> f64 {
    f.ds()
        .shape(n_v)
        .map(|s| BRepTool::vertex_tolerance(&Vertex(s.clone())))
        .unwrap_or(CONFUSION)
}

fn curve_is_line_edge(edge: &Edge) -> bool {
    let Some(c) = BRepTool::edge_curve(edge) else {
        return false;
    };
    let (a, b) = (c.first_parameter(), c.last_parameter());
    if !a.is_finite() || !b.is_finite() || (b - a).abs() <= 1e-15 {
        return false;
    }
    let p0 = c.d0(a);
    let p1 = c.d0(0.5 * (a + b));
    let p2 = c.d0(b);
    let v1 = occt_core::gp::GpVec::from_pnts(&p0, &p1);
    let v2 = occt_core::gp::GpVec::from_pnts(&p0, &p2);
    let m = v1.crossed(&v2).magnitude();
    let scale = v1.magnitude() * v2.magnitude();
    m < 1e-6 * scale.max(1e-12)
}

fn pb_on_face_lists(f: &PaveFiller, n_f: usize, pb: &BopdsPaveBlock) -> bool {
    let Some(fi) = f.ds().face_info(n_f) else {
        return false;
    };
    let e = pb.original_edge();
    let (t1, t2) = pb.range();
    let hit = |list: &[(usize, f64, f64)]| {
        list.iter().any(|&(ee, a, b)| {
            ee == e && (a - t1).abs() <= 1e-7 && (b - t2).abs() <= 1e-7
        })
    };
    hit(fi.paves_on()) || hit(fi.paves_in()) || hit(fi.paves())
}

fn collect_face_verts(f: &PaveFiller, n_f: usize) -> HashSet<usize> {
    let mut a_mvf: HashSet<usize> = HashSet::new();
    let Some(fi) = f.ds().face_info(n_f) else {
        return a_mvf;
    };
    for &v in fi.verts_on().iter().chain(fi.verts_in().iter()).chain(fi.verts_sc().iter()) {
        a_mvf.insert(v);
    }
    for &(e, t1, t2) in fi
        .paves_on()
        .iter()
        .chain(fi.paves_in().iter())
        .chain(fi.paves().iter())
    {
        for pb in f.ds().pave_blocks(e) {
            let (a, b) = pb.range();
            if (a - t1).abs() <= 1e-7 && (b - t2).abs() <= 1e-7 {
                let (n1, n2) = pb.indices();
                a_mvf.insert(n1);
                a_mvf.insert(n2);
            }
        }
    }
    a_mvf
}

struct EdgeFaceForceTask {
    n_e: usize,
    n_f: usize,
    pb: BopdsPaveBlock,
    t1: f64,
    t2: f64,
    fuzzy: f64,
}

/// `BOPAlgo_PaveFiller::ForceInterfEF` (no-argument overload).
pub fn force_interf_ef(f: &mut PaveFiller) -> Result<(), String> {
    if !f.is_primary() {
        return Ok(());
    }
    let a_nb_s = f.ds().nb_source_shapes();
    let mut mpb: Vec<BopdsPaveBlock> = Vec::new();
    let mut fence: HashSet<PbKey> = HashSet::new();
    for n_e in 0..a_nb_s {
        let Some(si) = f.ds().shape_info(n_e) else {
            continue;
        };
        if si.shape_type() != ShapeType::Edge {
            continue;
        }
        if !f.ds().has_pave_blocks(n_e) {
            continue;
        }
        if si.has_flag() {
            continue;
        }
        for pb in f.ds().pave_blocks(n_e) {
            let a_pbr = f.ds().real_pave_block(pb);
            if fence.insert(pb_key(&a_pbr)) {
                mpb.push(a_pbr);
            }
        }
    }
    force_interf_ef_on(f, &mpb, true)
}

/// `BOPAlgo_PaveFiller::ForceInterfEF` (pave-block map overload).
pub fn force_interf_ef_on(
    f: &mut PaveFiller,
    the_mpb: &[BopdsPaveBlock],
    the_add_interf: bool,
) -> Result<(), String> {
    if the_mpb.is_empty() {
        return Ok(());
    }
    let mut pb_ready: Vec<BopdsPaveBlock> = Vec::new();
    for pb in the_mpb {
        let mut a_pb = pb.clone();
        if !a_pb.has_shrunk_data() {
            let n_e = a_pb.original_edge();
            if let Some(es) = f.ds().shape(n_e).cloned() {
                let e = Edge(es);
                let tol = BRepTool::edge_tolerance(&e);
                fill_shrunk_data_for_block(f, &e, tol, &mut a_pb);
            }
        }
        if a_pb.has_shrunk_data() {
            pb_ready.push(a_pb);
        }
    }
    if pb_ready.is_empty() {
        return Ok(());
    }

    let b_si_check_mode = f.arguments().len() == 1;
    let a_nb_s = f.ds().nb_source_shapes();
    let mut pb_box: HashMap<(usize, u64, u64), occt_core::bnd::BndBox> = HashMap::new();
    let mut tasks: Vec<EdgeFaceForceTask> = Vec::new();
    for n_f in 0..a_nb_s {
        let Some(si) = f.ds().shape_info(n_f).cloned() else {
            continue;
        };
        if si.shape_type() != ShapeType::Face {
            continue;
        }
        if f.ds().face_info(n_f).is_none() {
            continue;
        }
        let Some(fs) = f.ds().shape(n_f).cloned() else {
            continue;
        };
        let a_f = Face(fs);
        let box_f = f
            .ds()
            .box_of(n_f)
            .cloned()
            .unwrap_or_else(|| shape_bbox(&a_f.0));
        let a_mvf = collect_face_verts(f, n_f);
        for a_pb in &pb_ready {
            let n_e = if a_pb.has_edge() {
                a_pb.edge()
            } else {
                a_pb.original_edge()
            };
            if n_e >= f.ds().nb_shapes() {
                continue;
            }
            if pb_on_face_lists(f, n_f, a_pb) {
                continue;
            }
            let (n_v1, n_v2) = a_pb.indices();
            if !a_mvf.contains(&n_v1) || !a_mvf.contains(&n_v2) {
                continue;
            }
            if !a_pb.has_edge() && f.ds().rank(n_f) == f.ds().rank(a_pb.original_edge()) {
                continue;
            }
            let Some(es) = f.ds().shape(n_e).cloned() else {
                continue;
            };
            let a_e = Edge(es);
            // `_5.cxx:903-910`: BVH Select of the shrunk PB box against the
            // face box. Whole-edge boxes would keep pairs that OCCT rejects.
            if let Some((_, _, _, _, shrunk_box)) = get_pb_box(&a_e, a_pb, &mut pb_box) {
                if shrunk_box.is_out_box(&box_f) {
                    continue;
                }
            }
            let Some(c) = BRepTool::edge_curve(&a_e) else {
                continue;
            };
            let (ts0, ts1, _) = a_pb.shrunk_data();
            let tm = intermediate_point(ts0, ts1);
            let (p_on_e, ve_tgt) = c.d1(tm);
            if ve_tgt.square_magnitude() < RESOLUTION {
                continue;
            }
            let Ok((u, v)) = f.context_mut().project_point_on_face(&a_f, &p_on_e) else {
                continue;
            };
            let Some(surf) = BRepTool::face_surface(&a_f) else {
                continue;
            };
            let p_on_s = surf.d0(u, v);
            let dist_mid = p_on_e.distance(&p_on_s);
            let a_tol_check = if b_si_check_mode {
                f.fuzzy_value()
            } else {
                2.0 * vertex_tol(f, n_v1).max(vertex_tol(f, n_v2))
            };
            if dist_mid > a_tol_check + f.fuzzy_value() {
                continue;
            }
            // `_5.cxx:1033`: `IsPointInFace(face, Pnt2d)` — In only, not On.
            let in_face = f
                .context_mut()
                .state_point_face(&a_f, (u, v), 0.0)
                .map(|st| st == FaceState::In)
                .unwrap_or(false);
            if !in_face {
                continue;
            }
            let mut b_use_add_tol = true;
            if !face_is_planar(&a_f) || !curve_is_line_edge(&a_e) {
                let vf_norm = occt_core::gp::GpVec::from_pnts(&p_on_s, &p_on_e);
                if vf_norm.square_magnitude() > RESOLUTION {
                    let a_cos = vf_norm.normalized().dot(&ve_tgt.normalized());
                    if a_cos.abs() > 0.4226 {
                        b_use_add_tol = false;
                    }
                }
            }
            let mut a_tol_add = 0.0;
            if b_use_add_tol {
                for t in [ts0, ts1] {
                    let p = c.d0(t);
                    if let Ok((uu, vv)) = f.context_mut().project_point_on_face(&a_f, &p) {
                        let d = p.distance(&surf.d0(uu, vv));
                        if d < a_tol_check && d > a_tol_add {
                            a_tol_add = d;
                        }
                    }
                }
                if a_tol_add > 0.0 {
                    a_tol_add -= BRepTool::edge_tolerance(&a_e) + BRepTool::face_tolerance(&a_f);
                    if a_tol_add < 0.0 {
                        a_tol_add = 0.0;
                    }
                }
            }
            let mut b_intersect = a_tol_add > 0.0;
            if !b_intersect {
                let done = f
                    .fpb_done()
                    .get(&n_f)
                    .map(|s| s.contains(&pb_key(a_pb)))
                    .unwrap_or(false);
                b_intersect = !done;
            }
            if !b_intersect {
                continue;
            }
            let (t1, t2) = a_pb.range();
            tasks.push(EdgeFaceForceTask {
                n_e,
                n_f,
                pb: a_pb.clone(),
                t1,
                t2,
                fuzzy: f.fuzzy_value() + a_tol_add,
            });
        }
    }
    if tasks.is_empty() {
        return Ok(());
    }

    let mut mpbli: PbFaceListMap = HashMap::new();
    for t in tasks {
        let Some(es) = f.ds().shape(t.n_e).cloned() else {
            continue;
        };
        let Some(fs) = f.ds().shape(t.n_f).cloned() else {
            continue;
        };
        let mut ef = EdgeFace::new();
        ef.set_edge(Edge(es));
        ef.set_face(Face(fs));
        ef.set_range(t.t1, t.t2);
        ef.set_fuzzy_value(t.fuzzy);
        ef.set_quick_coincidence_check(true);
        if let Err(msg) = ef.perform() {
            f.add_warning(format!(
                "force_interf_ef: intersection failed for edge {} / face {}: {msg}",
                t.n_e, t.n_f
            ));
            continue;
        }
        if !ef.is_done() {
            continue;
        }
        let cps = ef.common_parts();
        if cps.len() != 1 || cps[0].part_type() != CommonPartType::Edge {
            continue;
        }
        if the_add_interf {
            // `_5.cxx:1177-1181`: always Appended, then AddInterf.
            let ix = f.ds_mut().add_interf_ef(t.n_e, t.n_f, None);
            f.ds_mut().set_ef_common_range_at(
                ix,
                cps[0].range.first,
                cps[0].range.last,
            );
        }
        f.ds_mut().ensure_face_info(t.n_f);
        if let Some(fi) = f.ds_mut().face_info_mut(t.n_f) {
            fi.add_pave_block_in(&t.pb);
        }
        if the_add_interf {
            fill_map_pb_face(&t.pb, t.n_f, &mut mpbli);
        }
    }
    if !mpbli.is_empty() {
        let ctx = f.context().clone();
        perform_common_blocks_faces(&mpbli, f.ds_mut(), &ctx);
    }
    Ok(())
}
