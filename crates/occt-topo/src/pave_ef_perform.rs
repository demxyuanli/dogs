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
    correct_range_ef, algo_tools_is_in_range, make_new_vertex_ef, point_on_edge, algo_tools_vertex_parameter,
};
use crate::bopalgo_tools::{fill_map_pb_face, perform_common_blocks_faces, trsf_to_point, PbFaceListMap};
use crate::bbox_from_geometry::shape_bbox;
use crate::bopds::{BopdsDS, BopdsPaveBlock};
use crate::brep_surface::face_is_plane_type;
use crate::brep_tool::BRepTool;
use crate::edge_face::EdgeFace;
use crate::int_tools_full::IntToolsContext;
use crate::inttools_data::{CommonPartType, IntRange};
use crate::pave_common::update_vertex_sd;
use crate::pave_force_ve::{pave_block_box, update_vertices_of_cb};
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
        let e_val = Edge(e_shape);
        let f_val = Face(f_shape);
        f.ds_mut().ensure_face_info(n_f);
        let blocks = f.ds().pave_blocks(n_e).to_vec();
        let (on, inn) = face_paves_on_in(f.ds(), n_f);
        for pb in blocks {
            let pbr = f.ds().real_pave_block(&pb);
            if pb_on_face(f.ds(), n_f, &pbr) {
                continue;
            }
            let Some((t1_val, t2_val, ts1_val, ts2_val, bbe)) = pave_block_box(&e_val, &pb, &mut pb_box) else {
                continue;
            };
            if let Some(bbf) = f.ds().box_of(n_f) {
                if bbf.is_out_box(&bbe) {
                    continue;
                }
            }
            let (n_v1, n_v2) = pbr.indices();
            let b_v1 = on.contains(&n_v1) || inn.contains(&n_v1);
            let b_v2 = on.contains(&n_v2) || inn.contains(&n_v2);
            let b_express = b_v1 && b_v2;
            let sr = IntRange::new_unchecked(ts1_val, ts2_val);
            let anew_sr = correct_range_ef(&e_val, &f_val, sr);
            let pb_range0 = IntRange::new_unchecked(t1_val, t2_val);
            let pb_range_val = correct_range_ef(&e_val, &f_val, pb_range0);
            f.fpb_done_mut()
                .entry(n_f)
                .or_default()
                .insert(pb_key(&pb));
            tasks.push(EdgeFaceTask {
                n_e,
                n_f,
                pb,
                t1: t1_val,
                t2: t2_val,
                ts1: ts1_val,
                ts2: ts2_val,
                new_sr: anew_sr,
                pb_range: pb_range_val,
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
        let mut e_val = Edge(e_shape);
        let mut f_val = Face(f_shape);
        if let Some(trsf) = trsf_to_point(&box_e, &box_f, None, None) {
            let loc = TopLocLocation::composed(&trsf, TopLocLocation::identity());
            e_val.0.move_location(&loc);
            f_val.0.move_location(&loc);
        }
        let tol_e = BRepTool::edge_tolerance(&e_val);
        let tol_f = BRepTool::face_tolerance(&f_val);
        let mut ef = EdgeFace::new();
        ef.set_edge(e_val.clone());
        ef.set_face(f_val.clone());
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
                tol_e,
                tol_f,
            );
            continue;
        }
        let n_v = task.pb.indices();
        let mut ts1_val = task.new_sr.first;
        let mut ts2_val = task.new_sr.last;
        if cps[0].part_type() == CommonPartType::Vertex {
            reduce_intersection_range(
                f.ds(),
                n_v.0,
                n_v.1,
                task.n_e,
                task.n_f,
                &mut ts1_val,
                &mut ts2_val,
            );
        }
        let r1 = IntRange::new_unchecked(task.t1, ts1_val);
        let r2 = IntRange::new_unchecked(ts2_val, task.t2);
        let (mif_on, mif_in) = face_paves_on_in(f.ds(), task.n_f);
        let b_line_plane = curve_is_line_edge(&e_val) && face_is_plane_type(&f_val);
        let b_splittable = task.pb.is_splittable();

        for cp in cps.iter() {
            match cp.part_type() {
                CommonPartType::Vertex => {
                    let t = algo_tools_vertex_parameter(cp, cp.vertex_parameter1());
                let Ok(vnew) = make_new_vertex_ef(&e_val, t, &f_val) else {
                    continue;
                };
                let r = cp.range;
                let tol_to_decide = 5.0e-8;
                let mut b_is_on_pave = [
                    algo_tools_is_in_range(r1, r, tol_to_decide),
                    algo_tools_is_in_range(r2, r, tol_to_decide),
                ];
                if (b_is_on_pave[0] && b_is_on_pave[1])
                    || (b_line_plane && (b_is_on_pave[0] || b_is_on_pave[1]))
                {
                    let b_v0 = check_face_paves_index(n_v.0, &mif_on, &mif_in);
                    let b_v1 = check_face_paves_index(n_v.1, &mif_on, &mif_in);
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
                        if !check_face_paves_index(n_vj, &mif_on, &mif_in)
                            && force_interf_vf(f, n_vj, task.n_f)?
                        {
                            b_is_on_pave[j] = true;
                        } else if !check_face_paves_index(n_vj, &mif_on, &mif_in) {
                            b_is_on_pave[j] = false;
                        }
                    }
                }
                if b_is_on_pave[0] || b_is_on_pave[1] {
                    let pnew = BRepTool::vertex_point(&Vertex(vnew.clone()));
                    let min_dist = {
                        let ctx = f.context();
                        match ctx.project_point_on_face(&f_val, &pnew) {
                            Ok((u, v)) => BRepTool::face_surface(&f_val)
                                .map(|s| s.d0(u, v).distance(&pnew))
                                .unwrap_or(f64::MAX),
                            Err(_) => f64::MAX,
                        }
                    };
                    if min_dist >= INTERSECTION {
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
                        let v_val = Vertex(vs);
                        let p = BRepTool::vertex_point(&v_val);
                        let dist_pp = p.distance(&pnew);
                        let tol_val = BRepTool::vertex_tolerance(&v_val);
                        let mut max_dist = 1.0e4 * tol_val;
                        if tol_val < 0.01 {
                            max_dist = max_dist.min(0.1);
                        }
                        if dist_pp < max_dist {
                            update_vertex_sd(f, n_vj, dist_pp)?;
                            f.verts_to_avoid_extension_mut().insert(n_vj);
                        }
                    }
                    continue;
                }
                if check_face_paves_vertex(f.ds(), &Vertex(vnew.clone()), &mif_on) {
                    continue;
                }
                let mut tol_vnew = BRepTool::vertex_tolerance(&Vertex(vnew.clone()));
                tol_vnew = tol_vnew.max(tol_e.max(tol_f));
                Vertex(vnew.clone()).set_tolerance(tol_vnew);
                if b_line_plane {
                    tol_vnew = tol_vnew.max((cp.range.last - cp.range.first) / 2.0);
                }
                let pnew = BRepTool::vertex_point(&Vertex(vnew.clone()));
                let in_face = f
                    .context_mut()
                    .is_point_in_face(&f_val, &pnew, None, tol_vnew)
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
                    shape: vnew,
                    tol: tol_vnew,
                    pb1: task.pb.clone(),
                    pb2: task.pb.clone(),
                    n_e1: task.n_e,
                    n_e2: task.n_f,
                    t1: t,
                    t2: t,
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
                    let b_v0 = check_face_paves_index(n_v.0, &mif_on, &mif_in);
                    let b_v1 = check_face_paves_index(n_v.1, &mif_on, &mif_in);
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
