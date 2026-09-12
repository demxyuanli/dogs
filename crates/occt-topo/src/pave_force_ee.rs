//! `BOPAlgo_PaveFiller::ForceInterfEE` (`BOPAlgo_PaveFiller_3.cxx:997`).
//!
//! After vertex unification, look for additional EDGE-type common parts among
//! pave blocks that share the same bounding vertices. Real VERTEX
//! intersections are ignored here; only coincident EDGE common parts feed
//! `FillMap` / `PerformCommonBlocks`.

use std::collections::{HashMap, HashSet};

use occt_core::precision::{CONFUSION, RESOLUTION};
use occt_geom::Curve;

use crate::abs::ShapeType;
use crate::bopalgo_tools::{fill_map_pb_pb, perform_common_blocks_pb, PbPbListMap};
use crate::bopds::BopdsPaveBlock;
use crate::brep_tool::BRepTool;
use crate::edge_edge::EdgeEdge;
use crate::int_tools_full::IntToolsContext;
use crate::inttools_data::{CommonPartType, IntRange};
use crate::pave_ff_exist::pb_key;
use crate::pave_filler::PaveFiller;
use crate::shape::{Edge, Vertex};

fn curve_is_line(c: &dyn Curve) -> bool {
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

fn vertex_tol(f: &PaveFiller, n_v: usize) -> f64 {
    f.ds()
        .shape(n_v)
        .map(|s| BRepTool::vertex_tolerance(&Vertex(s.clone())))
        .unwrap_or(CONFUSION)
}

/// `BOPDS_DS::InitPaveBlocksForVertex`.
///
/// OCCT walks `myMapVE`; without that map the same vertices are recovered by
/// scanning source edges whose sub-shape list contains the vertex.
pub fn init_pave_blocks_for_vertex(f: &mut PaveFiller, n_v: usize) {
    let n = f.ds().nb_source_shapes();
    let mut edges: Vec<usize> = Vec::new();
    for n_e in 0..n {
        let Some(si) = f.ds().shape_info(n_e) else {
            continue;
        };
        if si.shape_type() != ShapeType::Edge {
            continue;
        }
        if si.has_subshape(n_v) {
            edges.push(n_e);
        }
    }
    for n_e in edges {
        f.ds_mut().init_pave_blocks_for_edge(n_e);
    }
}

fn use_add_tol(
    ctx: &mut IntToolsContext,
    c1: &dyn Curve,
    c2: &dyn Curve,
    e2: &Edge,
    mid: f64,
) -> Option<bool> {
    if curve_is_line(c1) && curve_is_line(c2) {
        return Some(true);
    }
    let (_, vt1) = c1.d1(mid);
    if vt1.square_magnitude() < RESOLUTION {
        return None;
    }
    let vt1 = vt1.normalized();
    let pm = c1.d0(mid);
    let Some(t2) = ctx.project_point_on_edge(e2, &pm) else {
        return None;
    };
    let (_, vt2) = c2.d1(t2);
    if vt2.square_magnitude() < RESOLUTION {
        return None;
    }
    let a_cos = vt1.dot(&vt2.normalized());
    Some(a_cos.abs() >= 0.9063)
}

struct EdgeEdgeForceTask {
    pb1: BopdsPaveBlock,
    pb2: BopdsPaveBlock,
    n_e1: usize,
    n_e2: usize,
    t11: f64,
    t12: f64,
    t21: f64,
    t22: f64,
    fuzzy: f64,
}

/// `BOPAlgo_PaveFiller::ForceInterfEE`.
pub fn force_interf_ee(f: &mut PaveFiller) -> Result<(), String> {
    let a_nb_s = f.ds().nb_source_shapes();
    let mut verts: Vec<usize> = Vec::new();
    for i in 0..a_nb_s {
        let Some(si) = f.ds().shape_info(i) else {
            continue;
        };
        if si.shape_type() == ShapeType::Vertex && f.ds().has_interf(i) {
            verts.push(i);
        }
    }
    for n_v in verts {
        init_pave_blocks_for_vertex(f, n_v);
    }

    let mut pb_map: HashMap<(usize, usize), Vec<BopdsPaveBlock>> = HashMap::new();
    let mut fence: HashSet<(usize, u64, u64)> = HashSet::new();
    for i in 0..a_nb_s {
        let Some(si) = f.ds().shape_info(i) else {
            continue;
        };
        if si.shape_type() != ShapeType::Edge {
            continue;
        };
        if !f.ds().has_pave_blocks(i) {
            continue;
        }
        if si.has_flag() {
            continue;
        }
        let pbs = f.ds().pave_blocks(i).to_vec();
        for pb in pbs {
            let a_pbr = f.ds().real_pave_block(&pb);
            if !fence.insert(pb_key(&a_pbr)) {
                continue;
            }
            let (n_v1, n_v2) = a_pbr.indices();
            pb_map
                .entry((n_v1.min(n_v2), n_v1.max(n_v2)))
                .or_default()
                .push(a_pbr);
        }
    }
    if pb_map.is_empty() {
        return Ok(());
    }

    let b_si_check_mode = f.arguments().len() == 1;
    let mut tasks: Vec<EdgeEdgeForceTask> = Vec::new();
    for ((n_v1, n_v2), group) in &pb_map {
        if group.len() < 2 {
            continue;
        }
        let a_tol_add = if b_si_check_mode {
            f.fuzzy_value()
        } else {
            2.0 * vertex_tol(f, *n_v1).max(vertex_tol(f, *n_v2))
        };
        for (k1, pb1) in group.iter().enumerate() {
            let n_e1 = pb1.original_edge();
            let i_r1 = f.ds().rank(n_e1);
            let (t11, t12) = pb1.range();
            let Some(e1s) = f.ds().shape(n_e1).cloned() else {
                continue;
            };
            let a_e1 = Edge(e1s);
            let Some(c1) = BRepTool::edge_curve(&a_e1) else {
                continue;
            };
            let mid = 0.5 * (t11 + t12);
            let (_, vt1) = c1.d1(mid);
            if vt1.square_magnitude() < RESOLUTION {
                continue;
            }
            let cb1 = f.ds().is_common_block(pb1);
            for pb2 in group.iter().skip(k1 + 1) {
                let n_e2 = pb2.original_edge();
                let i_r2 = f.ds().rank(n_e2);
                if i_r1 == i_r2 {
                    let v1_orig = !f.ds().is_new_shape(*n_v1) && f.ds().rank(*n_v1) == i_r1;
                    let v2_orig = !f.ds().is_new_shape(*n_v2) && f.ds().rank(*n_v2) == i_r2;
                    if v1_orig || v2_orig {
                        continue;
                    }
                }
                if cb1 && f.ds().is_common_block(pb2) {
                    if let (Some(c1b), Some(c2b)) =
                        (f.ds().common_block(pb1), f.ds().common_block(pb2))
                    {
                        if c1b.pave_block1().map(pb_key) == c2b.pave_block1().map(pb_key) {
                            continue;
                        }
                    }
                }
                let Some(e2s) = f.ds().shape(n_e2).cloned() else {
                    continue;
                };
                let a_e2 = Edge(e2s);
                let Some(c2) = BRepTool::edge_curve(&a_e2) else {
                    continue;
                };
                let (t21, t22) = pb2.range();
                let mut ctx = IntToolsContext::new();
                let Some(b_use_add) = use_add_tol(&mut ctx, c1.as_ref(), c2.as_ref(), &a_e2, mid)
                else {
                    continue;
                };
                let fuzzy = if b_use_add {
                    f.fuzzy_value() + a_tol_add
                } else {
                    f.fuzzy_value()
                };
                tasks.push(EdgeEdgeForceTask {
                    pb1: pb1.clone(),
                    pb2: pb2.clone(),
                    n_e1,
                    n_e2,
                    t11,
                    t12,
                    t21,
                    t22,
                    fuzzy,
                });
            }
        }
    }
    if tasks.is_empty() {
        return Ok(());
    }

    let mut mpblpb: PbPbListMap = HashMap::new();
    for t in tasks {
        let Some(e1s) = f.ds().shape(t.n_e1).cloned() else {
            continue;
        };
        let Some(e2s) = f.ds().shape(t.n_e2).cloned() else {
            continue;
        };
        let mut ee = EdgeEdge::with_edges(Edge(e1s), Edge(e2s));
        ee.set_range1(IntRange::new_unchecked(t.t11, t.t12));
        ee.set_range2(IntRange::new_unchecked(t.t21, t.t22));
        ee.set_fuzzy_value(t.fuzzy);
        if let Err(msg) = ee.perform() {
            f.add_warning(format!(
                "force_interf_ee: intersection failed for edges {} / {}: {msg}",
                t.n_e1, t.n_e2
            ));
            continue;
        }
        if !ee.is_done() {
            continue;
        }
        let cps = ee.common_parts();
        if cps.len() != 1 || cps[0].part_type() != CommonPartType::Edge {
            continue;
        }
        if f.ds().rank(t.n_e1) == f.ds().rank(t.n_e2) {
            f.add_warning(format!(
                "acquired self-interference: edges {} and {} of one argument",
                t.n_e1, t.n_e2
            ));
        }
        f.ds_mut().add_interf_ee(t.n_e1, t.n_e2, None);
        if f.ds().is_common_block(&t.pb1) {
            if let Some(cb) = f.ds().common_block(&t.pb1).cloned() {
                for pbx in cb.pave_blocks() {
                    fill_map_pb_pb(&t.pb1, pbx, &mut mpblpb);
                }
            }
        }
        if f.ds().is_common_block(&t.pb2) {
            if let Some(cb) = f.ds().common_block(&t.pb2).cloned() {
                for pbx in cb.pave_blocks() {
                    fill_map_pb_pb(&t.pb2, pbx, &mut mpblpb);
                }
            }
        }
        fill_map_pb_pb(&t.pb1, &t.pb2, &mut mpblpb);
    }
    if !mpblpb.is_empty() {
        let ctx = f.context().clone();
        perform_common_blocks_pb(&mpblpb, f.ds_mut(), &ctx);
    }
    Ok(())
}
