//! `BOPAlgo_PaveFiller::PerformEF` (`BOPAlgo_PaveFiller_5.cxx:165`).
//!
//! Collects edge/face solver tasks over every pave block of every interfering
//! pair, runs `IntTools_EdgeFace`, then post-treats VERTEX and EDGE common
//! parts with `CheckFacePaves` / `ForceInterfVF` / `ReduceIntersectionRange`.

use std::collections::{HashMap, HashSet};

use occt_core::precision::INTERSECTION;
use occt_core::toploc::TopLocLocation;

use crate::abs::ShapeType;
use crate::algo_tools_range::{
    correct_range_ef, is_in_range, make_new_vertex_ef, point_on_edge, vertex_parameter,
};
use crate::bopalgo_tools::{fill_map_pb_face, perform_common_blocks_faces, trsf_to_point, PbFaceListMap};
use crate::bbox_from_geometry::shape_bbox;
use crate::bopds::{BopdsDS, BopdsPaveBlock};
use crate::brep_surface::face_is_planar;
use crate::brep_tool::BRepTool;
use crate::edge_face::EdgeFace;
use crate::int_tools_full::IntToolsContext;
use crate::inttools_data::{CommonPartType, IntRange};
use crate::pave_common::update_vertex_sd;
use crate::pave_ee_aux::{get_pb_box, update_vertices_of_cb};
use crate::pave_ef::{
    check_face_paves_index, check_face_paves_vertex, face_paves_on_in, force_interf_vf,
    record_ef_distance, reduce_intersection_range,
};
use crate::pave_ff_exist::pb_key;
use crate::pave_filler::{EdgeRangeDistance, GlueEnum, PaveFiller};
use crate::pave_intersect::collect_pairs;
use crate::pave_new::{perform_new_vertices, NewVertexCpb};
use crate::shape::{Edge, Face, Vertex};

/// One `BOPAlgo_EdgeFace` solver task (`_5.cxx:278-305`).
struct EdgeFaceTask {
    n_e: usize,
    n_f: usize,
    pb: BopdsPaveBlock,
    t1: f64,
    t2: f64,
    ts1: f64,
    ts2: f64,
    new_sr: IntRange,
    pb_range: IntRange,
    express: bool,
}

fn pb_on_face(ds: &BopdsDS, n_f: usize, pbr: &BopdsPaveBlock) -> bool {
    let Some(fi) = ds.face_info(n_f) else {
        return false;
    };
    let e = pbr.original_edge();
    let (t1, t2) = pbr.range();
    fi.paves_on().iter().any(|&(ee, f, l)| {
        ee == e && (f - t1).abs() <= 1e-7 && (l - t2).abs() <= 1e-7
    })
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

/// `BOPAlgo_PaveFiller::PerformEF`.
pub fn perform_ef(f: &mut PaveFiller) -> Result<(), String> {
    crate::pave_shrunk::fill_shrunk_data_ef(f)?;
    let pairs = collect_pairs(f.ds(), ShapeType::Edge, ShapeType::Face);
    if pairs.is_empty() {
        return Ok(());
    }
    if f.glue() == GlueEnum::Full {
        for &(_, n_f) in &pairs {
            f.ds_mut().ensure_face_info(n_f);
        }
        return Ok(());
    }

    let mut pb_box: HashMap<(usize, u64, u64), occt_core::bnd::BndBox> = HashMap::new();
    let mut tasks: Vec<EdgeFaceTask> = Vec::new();
    for (n_e, n_f) in pairs {
        let Some(si_e) = f.ds().shape_info(n_e).cloned() else {
            continue;
        };
        if si_e.has_flag() {
            continue;
        }
        let Some(e_shape) = f.ds().shape(n_e).cloned() else {
            continue;
        };
        let Some(f_shape) = f.ds().shape(n_f).cloned() else {
            continue;
        };
        let a_e = Edge(e_shape);
        let a_f = Face(f_shape);
        f.ds_mut().ensure_face_info(n_f);
        let blocks = f.ds().pave_blocks(n_e).to_vec();
        let (on, inn) = face_paves_on_in(f.ds(), n_f);
        for pb in blocks {
            let pbr = f.ds().real_pave_block(&pb);
            if pb_on_face(f.ds(), n_f, &pbr) {
                continue;
            }
            let Some((a_t1, a_t2, a_ts1, a_ts2, a_bbe)) = get_pb_box(&a_e, &pb, &mut pb_box) else {
                continue;
            };
            if let Some(a_bbf) = f.ds().box_of(n_f) {
                if a_bbf.is_out_box(&a_bbe) {
                    continue;
                }
            }
            let (n_v1, n_v2) = pbr.indices();
            let b_v1 = on.contains(&n_v1) || inn.contains(&n_v1);
            let b_v2 = on.contains(&n_v2) || inn.contains(&n_v2);
            let b_express = b_v1 && b_v2;
            let a_sr = IntRange::new_unchecked(a_ts1, a_ts2);
            let anew_sr = correct_range_ef(&a_e, &a_f, a_sr);
            let a_pb_range0 = IntRange::new_unchecked(a_t1, a_t2);
            let a_pb_range = correct_range_ef(&a_e, &a_f, a_pb_range0);
            f.fpb_done_mut()
                .entry(n_f)
                .or_default()
                .insert(pb_key(&pb));
            tasks.push(EdgeFaceTask {
                n_e,
                n_f,
                pb,
                t1: a_t1,
                t2: a_t2,
                ts1: a_ts1,
                ts2: a_ts2,
                new_sr: anew_sr,
                pb_range: a_pb_range,
                express: b_express,
            });
        }
    }

    let mut mpbl: PbFaceListMap = HashMap::new();
    let mut mvcpb: Vec<NewVertexCpb> = Vec::new();
    let mut miefc: HashSet<usize> = HashSet::new();
    let mut distances: HashMap<(usize, usize), Vec<EdgeRangeDistance>> =
        std::mem::take(f.distances_mut());

    for task in tasks {
        let Some(e_shape) = f.ds().shape(task.n_e).cloned() else {
            continue;
        };
        let Some(f_shape) = f.ds().shape(task.n_f).cloned() else {
            continue;
        };
        let box_e = f
            .ds()
            .box_of(task.n_e)
            .copied()
            .unwrap_or_else(|| shape_bbox(&e_shape));
        let box_f = f
            .ds()
            .box_of(task.n_f)
            .copied()
            .unwrap_or_else(|| shape_bbox(&f_shape));
        let mut a_e = Edge(e_shape);
        let mut a_f = Face(f_shape);
        if let Some(trsf) = trsf_to_point(&box_e, &box_f, None, None) {
            let loc = TopLocLocation::composed(&trsf, TopLocLocation::identity());
            a_e.0.move_location(&loc);
            a_f.0.move_location(&loc);
        }
        let a_tol_e = BRepTool::edge_tolerance(&a_e);
        let a_tol_f = BRepTool::face_tolerance(&a_f);
        let mut ef = EdgeFace::new();
        ef.set_edge(a_e.clone());
        ef.set_face(a_f.clone());
        ef.set_fuzzy_value(f.fuzzy_value());
        ef.set_quick_coincidence_check(task.express);
        ef.set_range(task.pb_range.first, task.pb_range.last);
        if let Err(msg) = ef.perform() {
            f.add_warning(format!("perform_ef: intersection failed: {msg}"));
            continue;
        }
        if !ef.is_done() || ef.error_status() != 0 {
            f.add_warning("perform_ef: intersection not done".into());
            continue;
        }
        let cps = ef.common_parts().to_vec();
        if cps.is_empty() {
            record_ef_distance(
                &mut distances,
                task.n_e,
                task.n_f,
                &task.pb,
                ef.minimal_distance(),
                a_tol_e,
                a_tol_f,
            );
            continue;
        }
        let n_v = task.pb.indices();
        let mut a_ts1 = task.new_sr.first;
        let mut a_ts2 = task.new_sr.last;
        if cps[0].part_type() == CommonPartType::Vertex {
            reduce_intersection_range(
                f.ds(),
                n_v.0,
                n_v.1,
                task.n_e,
                task.n_f,
                &mut a_ts1,
                &mut a_ts2,
            );
        }
        let a_r1 = IntRange::new_unchecked(task.t1, a_ts1);
        let a_r2 = IntRange::new_unchecked(a_ts2, task.t2);
        let (a_mif_on, a_mif_in) = face_paves_on_in(f.ds(), task.n_f);
        let b_line_plane = curve_is_line_edge(&a_e) && face_is_planar(&a_f);
        let b_splittable = task.pb.is_splittable();

        for cp in cps.iter() {
            match cp.part_type() {
                CommonPartType::Vertex => {
                    let a_t = vertex_parameter(cp, cp.vertex_parameter1());
                let Ok(a_vnew) = make_new_vertex_ef(&a_e, a_t, &a_f) else {
                    continue;
                };
                let a_r = cp.range;
                let a_tol_to_decide = 5.0e-8;
                let mut b_is_on_pave = [
                    is_in_range(a_r1, a_r, a_tol_to_decide),
                    is_in_range(a_r2, a_r, a_tol_to_decide),
                ];
                if (b_is_on_pave[0] && b_is_on_pave[1])
                    || (b_line_plane && (b_is_on_pave[0] || b_is_on_pave[1]))
                {
                    let b_v0 = check_face_paves_index(n_v.0, &a_mif_on, &a_mif_in);
                    let b_v1 = check_face_paves_index(n_v.1, &a_mif_on, &a_mif_in);
                    if b_v0 && b_v1 {
                        // `_5.cxx:427-438`: Appended InterfEF, type promoted to
                        // EDGE, FillMap. A later VERTEX record of the same pair
                        // must not collapse this EDGE slot.
                        let ix = f.ds_mut().add_interf_ef(task.n_e, task.n_f, None);
                        f.ds_mut().set_ef_common_range_at(
                            ix,
                            cp.range.first,
                            cp.range.last,
                        );
                        miefc.insert(task.n_f);
                        fill_map_pb_face(&task.pb, task.n_f, &mut mpbl);
                        break;
                    }
                }
                if !b_splittable {
                    continue;
                }
                for j in 0..2 {
                    if b_is_on_pave[j] {
                        let n_vj = if j == 0 { n_v.0 } else { n_v.1 };
                        if !check_face_paves_index(n_vj, &a_mif_on, &a_mif_in)
                            && force_interf_vf(f, n_vj, task.n_f)?
                        {
                            b_is_on_pave[j] = true;
                        } else if !check_face_paves_index(n_vj, &a_mif_on, &a_mif_in) {
                            b_is_on_pave[j] = false;
                        }
                    }
                }
                if b_is_on_pave[0] || b_is_on_pave[1] {
                    let a_pnew = BRepTool::vertex_point(&Vertex(a_vnew.clone()));
                    let a_min_dist = {
                        let ctx = f.context();
                        match ctx.project_point_on_face(&a_f, &a_pnew) {
                            Ok((u, v)) => BRepTool::face_surface(&a_f)
                                .map(|s| s.d0(u, v).distance(&a_pnew))
                                .unwrap_or(f64::MAX),
                            Err(_) => f64::MAX,
                        }
                    };
                    if a_min_dist >= INTERSECTION {
                        continue;
                    }
                    for j in 0..2 {
                        if !b_is_on_pave[j] {
                            continue;
                        }
                        let n_vj = if j == 0 { n_v.0 } else { n_v.1 };
                        let Some(vs) = f.ds().shape(n_vj).cloned() else {
                            continue;
                        };
                        let a_v = Vertex(vs);
                        let a_p = BRepTool::vertex_point(&a_v);
                        let a_dist_pp = a_p.distance(&a_pnew);
                        let a_tol = BRepTool::vertex_tolerance(&a_v);
                        let mut a_max_dist = 1.0e4 * a_tol;
                        if a_tol < 0.01 {
                            a_max_dist = a_max_dist.min(0.1);
                        }
                        if a_dist_pp < a_max_dist {
                            update_vertex_sd(f, n_vj, a_dist_pp)?;
                            f.verts_to_avoid_extension_mut().insert(n_vj);
                        }
                    }
                    continue;
                }
                if check_face_paves_vertex(f.ds(), &Vertex(a_vnew.clone()), &a_mif_on) {
                    continue;
                }
                let mut a_tol_vnew = BRepTool::vertex_tolerance(&Vertex(a_vnew.clone()));
                a_tol_vnew = a_tol_vnew.max(a_tol_e.max(a_tol_f));
                Vertex(a_vnew.clone()).set_tolerance(a_tol_vnew);
                if b_line_plane {
                    a_tol_vnew = a_tol_vnew.max((cp.range.last - cp.range.first) / 2.0);
                }
                let a_pnew = BRepTool::vertex_point(&Vertex(a_vnew.clone()));
                let in_face = f
                    .context_mut()
                    .is_point_in_face(&a_f, &a_pnew, None, a_tol_vnew)
                    .unwrap_or(false);
                if !in_face {
                    continue;
                }
                miefc.insert(task.n_f);
                // `_5.cxx:530-542`: Appended InterfEF, `SetIndexInterf(iX)`.
                let ix = f.ds_mut().add_interf_ef(task.n_e, task.n_f, None);
                f.ds_mut().set_ef_common_range_at(
                    ix,
                    cp.range.first,
                    cp.range.last,
                );
                mvcpb.push(NewVertexCpb {
                    shape: a_vnew,
                    tol: a_tol_vnew,
                    pb1: task.pb.clone(),
                    pb2: task.pb.clone(),
                    n_e1: task.n_e,
                    n_e2: task.n_f,
                    t1: a_t,
                    t2: a_t,
                    range1_first: cp.range.first,
                    range1_last: cp.range.last,
                    is_ee: false,
                    index_interf: Some(ix),
                });
                }
                CommonPartType::Edge => {
                    miefc.insert(task.n_f);
                    // `_5.cxx:549-564`: always Appended; AddInterf even when a
                    // bounding vertex is missing; FillMap only when both vertices
                    // are already On or In.
                    let ix = f.ds_mut().add_interf_ef(task.n_e, task.n_f, None);
                    let b_v0 = check_face_paves_index(n_v.0, &a_mif_on, &a_mif_in);
                    let b_v1 = check_face_paves_index(n_v.1, &a_mif_on, &a_mif_in);
                    if !b_v0 || !b_v1 {
                        break;
                    }
                    f.ds_mut().set_ef_common_range_at(
                        ix,
                        cp.range.first,
                        cp.range.last,
                    );
                    fill_map_pb_face(&task.pb, task.n_f, &mut mpbl);
                }
                _ => {}
            }
        }
        let _ = (task.ts1, task.ts2);
        let _ = point_on_edge;
    }

    *f.distances_mut() = distances;
    if !mpbl.is_empty() {
        let ictx = IntToolsContext::new();
        perform_common_blocks_faces(&mpbl, f.ds_mut(), &ictx);
    }
    update_vertices_of_cb(f)?;
    perform_new_vertices(f, &mvcpb, false)?;
    if f.has_errors() {
        return Err(f
            .errors()
            .first()
            .cloned()
            .unwrap_or_else(|| "perform_ef: errors".into()));
    }
    f.ds_mut().update_face_info_in_faces(&miefc);
    Ok(())
}
