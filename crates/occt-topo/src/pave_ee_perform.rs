//! `BOPAlgo_PaveFiller::PerformEE` (`BOPAlgo_PaveFiller_3.cxx:140`).
//!
//! Builds `BOPAlgo_EdgeEdge` tasks per pave-block pair, runs
//! `IntTools_EdgeEdge`, then post-treats VERTEX hits with `ForceInterfVE` /
//! `UpdateVertex` / `myVertsToAvoidExtension` and EDGE coincidence with
//! `FillMap` + `PerformCommonBlocks`.

use std::collections::{HashMap, HashSet};

use occt_core::precision::{CONFUSION, INTERSECTION};
use occt_core::toploc::TopLocLocation;

use crate::abs::ShapeType;
use crate::algo_tools_range::{
    algo_tools_is_on_pave1, make_new_vertex_ee, point_on_edge,
};
use crate::bopalgo_tools::{
    fill_map_pb_pb, perform_common_blocks_pb, trsf_to_point, PbPbListMap,
};
use crate::bbox_from_geometry::shape_bbox;
use crate::bopds::BopdsPaveBlock;
use crate::brep_tool::BRepTool;
use crate::edge_edge::EdgeEdge;
use crate::int_tools_full::IntToolsContext;
use crate::inttools_data::{CommonPartType, IntRange};
use crate::pave_common::update_vertex_sd;
use crate::pave_force_ve::{force_interf_ve, pave_block_box, update_vertices_of_cb};
use crate::pave_filler::PaveFiller;
use crate::pave_intersect::collect_pairs;
use crate::pave_split_blocks;
use crate::pave_new::{perform_new_vertices, NewVertexCpb};
use crate::shape::{Edge, Vertex};

struct EdgeEdgeTask {
    n_e1: usize,
    n_e2: usize,
    pb1: BopdsPaveBlock,
    pb2: BopdsPaveBlock,
    t11: f64,
    t12: f64,
    t21: f64,
    t22: f64,
    express: bool,
}

fn curve_kind_line_or_circle(edge: &Edge) -> (bool, bool) {
    let Some(c) = BRepTool::edge_curve(edge) else {
        return (false, false);
    };
    let (a, b) = (c.first_parameter(), c.last_parameter());
    if !a.is_finite() || !b.is_finite() || (b - a).abs() <= 1e-15 {
        return (false, false);
    }
    let p0 = c.d0(a);
    let p1 = c.d0(0.5 * (a + b));
    let p2 = c.d0(b);
    let v1 = occt_core::gp::GpVec::from_pnts(&p0, &p1);
    let v2 = occt_core::gp::GpVec::from_pnts(&p0, &p2);
    let m = v1.crossed(&v2).magnitude();
    let scale = v1.magnitude() * v2.magnitude();
    let is_line = m < 1e-6 * scale.max(1e-12);
    let is_circle = !is_line && {
        let d01 = p0.distance(&p1);
        let d12 = p1.distance(&p2);
        let d20 = p2.distance(&p0);
        (d01 - d12).abs() < 0.05 * d01.max(d12).max(1e-12) && d20 > 1e-12
    };
    (is_line, is_circle)
}

/// `BOPAlgo_PaveFiller::PerformEE`.
pub fn perform_ee(f: &mut PaveFiller) -> Result<(), String> {
    crate::pave_shrunk::fill_shrunk_data_ee(f)?;
    let pairs = collect_pairs(f.ds(), ShapeType::Edge, ShapeType::Edge);
    if pairs.is_empty() {
        return Ok(());
    }
    let mut pb_box = HashMap::new();
    let mut tasks: Vec<EdgeEdgeTask> = Vec::new();
    for (n_e1, n_e2) in pairs {
        let Some(si1) = f.ds().shape_info(n_e1).cloned() else {
            continue;
        };
        let Some(si2) = f.ds().shape_info(n_e2).cloned() else {
            continue;
        };
        if si1.has_flag() || si2.has_flag() {
            continue;
        }
        let pbs1 = f.ds().pave_blocks(n_e1).to_vec();
        let pbs2 = f.ds().pave_blocks(n_e2).to_vec();
        if pbs1.is_empty() || pbs2.is_empty() {
            continue;
        }
        let Some(e1s) = f.ds().shape(n_e1).cloned() else {
            continue;
        };
        let Some(e2s) = f.ds().shape(n_e2).cloned() else {
            continue;
        };
        let e1 = Edge(e1s);
        let e2 = Edge(e2s);
        for pb1 in &pbs1 {
            let Some((t11, t12, _ts11, _ts12, bb1)) = pave_block_box(&e1, pb1, &mut pb_box) else {
                continue;
            };
            let (n_v11, n_v12) = pb1.indices();
            for pb2 in &pbs2 {
                let Some((t21, t22, _ts21, _ts22, bb2)) = pave_block_box(&e2, pb2, &mut pb_box)
                else {
                    continue;
                };
                if bb1.is_out_box(&bb2) {
                    continue;
                }
                let (n_v21, n_v22) = pb2.indices();
                let express = (n_v11 == n_v21 && n_v12 == n_v22)
                    || (n_v12 == n_v21 && n_v11 == n_v22);
                tasks.push(EdgeEdgeTask {
                    n_e1,
                    n_e2,
                    pb1: pb1.clone(),
                    pb2: pb2.clone(),
                    t11,
                    t12,
                    t21,
                    t22,
                    express,
                });
            }
        }
    }

    let mut mpblpb: PbPbListMap = HashMap::new();
    let mut mvcpb: Vec<NewVertexCpb> = Vec::new();
    let mut m_edges: Vec<usize> = Vec::new();

    for task in tasks {
        let Some(e1s) = f.ds().shape(task.n_e1).cloned() else {
            continue;
        };
        let Some(e2s) = f.ds().shape(task.n_e2).cloned() else {
            continue;
        };
        let box1 = f
            .ds()
            .box_of(task.n_e1)
            .copied()
            .unwrap_or_else(|| shape_bbox(&e1s));
        let box2 = f
            .ds()
            .box_of(task.n_e2)
            .copied()
            .unwrap_or_else(|| shape_bbox(&e2s));
        let mut e1 = Edge(e1s);
        let mut e2 = Edge(e2s);
        if let Some(trsf) = trsf_to_point(&box1, &box2, None, None) {
            let loc = TopLocLocation::composed(&trsf, TopLocLocation::identity());
            e1.0.move_location(&loc);
            e2.0.move_location(&loc);
        }
        let mut ee = EdgeEdge::with_edges(e1.clone(), e2.clone());
        ee.set_range1(IntRange::new_unchecked(task.t11, task.t12));
        ee.set_range2(IntRange::new_unchecked(task.t21, task.t22));
        ee.set_fuzzy_value(f.fuzzy_value());
        ee.set_quick_coincidence_check(task.express);
        if let Err(msg) = ee.perform() {
            f.add_warning(format!("perform_ee: intersection failed: {msg}"));
            continue;
        }
        if !ee.is_done() {
            continue;
        }
        let mut pb1 = task.pb1.clone();
        let mut pb2 = task.pb2.clone();
        let (t11_val, t12_val) = pb1.range();
        let (t21_val, t22_val) = pb2.range();
        let (ts11, ts12, b_split1) = if pb1.has_shrunk_data() {
            pb1.shrunk_data()
        } else {
            (t11_val, t12_val, false)
        };
        let (ts21, ts22, b_split2) = if pb2.has_shrunk_data() {
            pb2.shrunk_data()
        } else {
            (t21_val, t22_val, false)
        };
        let r11 = IntRange::new_unchecked(t11_val, ts11);
        let r12 = IntRange::new_unchecked(ts12, t12_val);
        let r21 = IntRange::new_unchecked(t21_val, ts21);
        let r22 = IntRange::new_unchecked(ts22, t22_val);
        let (line1, circ1) = curve_kind_line_or_circle(&e1);
        let (line2, circ2) = curve_kind_line_or_circle(&e2);
        let b_analytical = (line1 && circ2) || (circ1 && line2);

        let points = ee.points().to_vec();
        let cps = ee.common_parts().to_vec();
        let nb = points.len() + cps.len();
        if nb == 0 {
            continue;
        }

        for pt in &points {
            if !b_split1 || !b_split2 {
                continue;
            }
            let t1_val = pt.uv1().0;
            let t2_val = pt.uv2().0;
            let tol_val = CONFUSION;
            let cr1 = IntRange::new_unchecked(t1_val, t1_val);
            let cr2 = IntRange::new_unchecked(t2_val, t2_val);
            let mut b_is_on_pave = [
                algo_tools_is_on_pave1(t1_val, r11.first, r11.last, tol_val)
                    || algo_tools_is_on_pave1(r11.first, cr1.first, cr1.last, tol_val),
                algo_tools_is_on_pave1(t1_val, r12.first, r12.last, tol_val)
                    || algo_tools_is_on_pave1(r12.last, cr1.first, cr1.last, tol_val),
                algo_tools_is_on_pave1(t2_val, r21.first, r21.last, tol_val)
                    || algo_tools_is_on_pave1(r21.first, cr2.first, cr2.last, tol_val),
                algo_tools_is_on_pave1(t2_val, r22.first, r22.last, tol_val)
                    || algo_tools_is_on_pave1(r22.last, cr2.first, cr2.last, tol_val),
            ];
            let mut n_v = [0usize; 4];
            let (a, b) = pb1.indices();
            n_v[0] = a;
            n_v[1] = b;
            let (c, d) = pb2.indices();
            n_v[2] = c;
            n_v[3] = d;
            if (b_is_on_pave[0] && b_is_on_pave[2])
                || (b_is_on_pave[0] && b_is_on_pave[3])
                || (b_is_on_pave[1] && b_is_on_pave[2])
                || (b_is_on_pave[1] && b_is_on_pave[3])
            {
                continue;
            }
            let mut is_v_exists = false;
            for j in 0..4 {
                if b_is_on_pave[j] {
                    let pb = if j < 2 { &mut pb2 } else { &mut pb1 };
                    b_is_on_pave[j] = force_interf_ve(f, n_v[j], pb, &mut m_edges)?;
                    if b_is_on_pave[j] {
                        is_v_exists = true;
                    }
                }
            }
            crate::pave_force_ve::write_pave_block(f, &pb1);
            crate::pave_force_ve::write_pave_block(f, &pb2);
            let Ok(vnew) = make_new_vertex_ee(&e1, t1_val, &e2, t2_val) else {
                continue;
            };
            let pnew = BRepTool::vertex_point(&Vertex(vnew.clone()));
            if is_v_exists {
                let p_on_e1 = point_on_edge(&e1, t1_val).unwrap_or(pnew);
                let p_on_e2 = point_on_edge(&e2, t2_val).unwrap_or(pnew);
                if p_on_e1.distance(&p_on_e2) > INTERSECTION {
                    continue;
                }
                for j in 0..4 {
                    if !b_is_on_pave[j] {
                        continue;
                    }
                    let Some(vs) = f.ds().shape(n_v[j]).cloned() else {
                        continue;
                    };
                    let v = Vertex(vs);
                    let dist_pp = BRepTool::vertex_point(&v).distance(&pnew);
                    update_vertex_sd(f, n_v[j], dist_pp)?;
                    f.verts_to_avoid_extension_mut().insert(n_v[j]);
                }
            }
            let mut tol_vnew = BRepTool::vertex_tolerance(&Vertex(vnew.clone()));
            if b_analytical {
                let tol_min = if line1 {
                    0.0
                } else {
                    0.0
                };
                tol_vnew = tol_vnew.max(tol_min);
            }
            let mut i_found = false;
            let mut mv = HashSet::new();
            mv.insert(n_v[0]);
            mv.insert(n_v[1]);
            let mut n_vs: Vec<usize> = Vec::new();
            if mv.contains(&n_v[2]) {
                n_vs.push(n_v[2]);
            }
            if mv.contains(&n_v[3]) {
                n_vs.push(n_v[3]);
            }
            for &n_vx in &n_vs {
                let Some(vs) = f.ds().shape(n_vx).cloned() else {
                    continue;
                };
                let vx = Vertex(vs);
                let tol_vx = BRepTool::vertex_tolerance(&vx);
                let px = BRepTool::vertex_point(&vx);
                let d2 = pnew.square_distance(&px);
                let dt = 100.0 * (tol_vnew + tol_vx) * (tol_vnew + tol_vx);
                if d2 < dt {
                    i_found = true;
                    break;
                }
            }
            if i_found {
                continue;
            }
            f.ds_mut().add_interf_ee(task.n_e1, task.n_e2, None);
            mvcpb.push(NewVertexCpb {
                shape: vnew,
                tol: tol_vnew,
                pb1: pb1.clone(),
                pb2: pb2.clone(),
                n_e1: task.n_e1,
                n_e2: task.n_e2,
                t1: t1_val,
                t2: t2_val,
                range1_first: t1_val,
                range1_last: t1_val,
                is_ee: true,
                index_interf: None,
            });
        }

        for cp in &cps {
            if cp.part_type() != CommonPartType::Edge {
                continue;
            }
            // OCCT: EDGE type with aNbCPrts > 1 is ignored (`_3.cxx:529-533`).
            if cps.len() > 1 {
                continue;
            }
            if !pb1.has_same_bounds(&pb2) {
                continue;
            }
            f.ds_mut().add_interf_ee(task.n_e1, task.n_e2, None);
            fill_map_pb_pb(&pb1, &pb2, &mut mpblpb);
        }
    }

    if !mpblpb.is_empty() {
        let ictx = IntToolsContext::new();
        perform_common_blocks_pb(&mpblpb, f.ds_mut(), &ictx);
    }
    update_vertices_of_cb(f)?;
    perform_new_vertices(f, &mvcpb, true)?;
    if f.has_errors() {
        return Err(f
            .errors()
            .first()
            .cloned()
            .unwrap_or_else(|| "perform_ee: errors".into()));
    }
    if !m_edges.is_empty() {
        for cpb in &mvcpb {
            m_edges.retain(|e| *e != cpb.pb1.original_edge() && *e != cpb.pb2.original_edge());
        }
        if !m_edges.is_empty() {
            pave_split_blocks::split_pave_blocks(f, &m_edges, false)?;
        }
    }
    Ok(())
}
