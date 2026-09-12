//! `BOPAlgo_Tools` — common-block assembly, vertex clustering, CB tolerance.
//!
//! Source: `BOPAlgo_Tools.cxx` (`FillMap` at 91, `PerformCommonBlocks` at 107
//! and 191, `ComputeToleranceOfCB` at 248, `IntersectVertices` at 1119).
//! Connected components replace OCCT `MakeBlocks`; vertex pairs use
//! [`crate::bop_aabb_faces::AabbTree`] (`BOPTools_BoxTree`) with the same
//! `Bnd_Box::IsOut` predicate.

use std::collections::{HashMap, HashSet};

use occt_core::bnd::BndBox;
use occt_core::gp::{GpPnt, GpTrsf};
use occt_core::precision::CONFUSION;

use crate::bopds::{BopdsCommonBlock, BopdsDS, BopdsPaveBlock};
use crate::brep_tool::BRepTool;
use crate::int_tools_full::IntToolsContext;
use crate::pave_ff_exist::{pb_key, PbKey};
use crate::shape::{Edge, Face, TopoShape, Vertex};

/// Map of a pave block to the faces it lies on (`FillMap` PB/face overload).
pub type PbFaceListMap = HashMap<PbKey, (BopdsPaveBlock, Vec<usize>)>;

/// Map of pave block to sibling pave blocks (`PerformCommonBlocks` PB/PB).
pub type PbPbListMap = HashMap<PbKey, (BopdsPaveBlock, Vec<BopdsPaveBlock>)>;

/// `BOPAlgo_Tools::FillMap` (pave block, pave block).
pub fn fill_map_pb_pb(a_pb1: &BopdsPaveBlock, a_pb2: &BopdsPaveBlock, a_mpblpb: &mut PbPbListMap) {
    a_mpblpb
        .entry(pb_key(a_pb1))
        .or_insert_with(|| (a_pb1.clone(), Vec::new()))
        .1
        .push(a_pb2.clone());
    a_mpblpb
        .entry(pb_key(a_pb2))
        .or_insert_with(|| (a_pb2.clone(), Vec::new()))
        .1
        .push(a_pb1.clone());
}

/// `BOPAlgo_Tools::FillMap` (pave block, face).
pub fn fill_map_pb_face(a_pb: &BopdsPaveBlock, n_f: usize, a_mpbli: &mut PbFaceListMap) {
    let e = a_mpbli
        .entry(pb_key(a_pb))
        .or_insert_with(|| (a_pb.clone(), Vec::new()));
    if !e.1.contains(&n_f) {
        e.1.push(n_f);
    }
}

/// `BOPAlgo_Tools::FillMap` (integer pair) used by `IntersectVertices`.
pub fn fill_map_pair(n1: usize, n2: usize, mili: &mut HashMap<usize, Vec<usize>>) {
    mili.entry(n1).or_default().push(n2);
    mili.entry(n2).or_default().push(n1);
}

/// Connected components of an undirected adjacency map (`BOPAlgo_Tools::MakeBlocks`).
pub fn make_blocks_int(mili: &HashMap<usize, Vec<usize>>) -> Vec<Vec<usize>> {
    let mut seen: HashSet<usize> = HashSet::new();
    let mut blocks = Vec::new();
    for &start in mili.keys() {
        if !seen.insert(start) {
            continue;
        }
        let mut stack = vec![start];
        let mut block = Vec::new();
        while let Some(n) = stack.pop() {
            block.push(n);
            if let Some(nbrs) = mili.get(&n) {
                for &q in nbrs {
                    if seen.insert(q) {
                        stack.push(q);
                    }
                }
            }
        }
        blocks.push(block);
    }
    blocks
}

/// Connected components of a pave-block adjacency map.
pub fn make_blocks_pb(mpblpb: &PbPbListMap) -> Vec<Vec<BopdsPaveBlock>> {
    let mut seen: HashSet<PbKey> = HashSet::new();
    let mut blocks = Vec::new();
    for (start, (pb0, _)) in mpblpb {
        if !seen.insert(*start) {
            continue;
        }
        let mut stack = vec![*start];
        let mut block: Vec<BopdsPaveBlock> = Vec::new();
        while let Some(k) = stack.pop() {
            if let Some((pb, nbrs)) = mpblpb.get(&k) {
                block.push(pb.clone());
                for n in nbrs {
                    let nk = pb_key(n);
                    if seen.insert(nk) {
                        stack.push(nk);
                    }
                }
            }
        }
        let _ = pb0;
        if block.len() >= 2 {
            blocks.push(block);
        }
    }
    blocks
}

/// `BOPAlgo_Tools::PerformCommonBlocks` (pave-block / pave-block overload).
pub fn perform_common_blocks_pb(
    mpblpb: &PbPbListMap,
    ds: &mut BopdsDS,
    ctx: &IntToolsContext,
) {
    let blocks = make_blocks_pb(mpblpb);
    for a_lpb in blocks {
        if a_lpb.len() < 2 {
            continue;
        }
        let mut a_m_faces: HashSet<usize> = HashSet::new();
        let mut a_l_faces: Vec<usize> = Vec::new();
        let mut a_cb: Option<BopdsCommonBlock> = None;
        for a_pb in &a_lpb {
            if ds.is_common_block(a_pb) {
                if let Some(a_cbx) = ds.common_block(a_pb) {
                    for &n_f in a_cbx.faces() {
                        if a_m_faces.insert(n_f) {
                            a_l_faces.push(n_f);
                        }
                    }
                    if a_cb.is_none() {
                        a_cb = Some(a_cbx.clone());
                    }
                }
            }
        }
        let mut a_cb = a_cb.unwrap_or_else(BopdsCommonBlock::new);
        a_cb.set_pave_blocks(a_lpb.clone());
        a_cb.set_faces(a_l_faces);
        let a_tol = compute_tolerance_of_cb(&a_cb, ds, ctx);
        a_cb.set_tolerance(a_tol);
        for a_pb in &a_lpb {
            ds.set_common_block(a_pb, a_cb.clone());
        }
    }
}

/// `BOPAlgo_Tools::PerformCommonBlocks` (pave-block / face-list overload).
pub fn perform_common_blocks_faces(
    mpbli: &PbFaceListMap,
    ds: &mut BopdsDS,
    ctx: &IntToolsContext,
) {
    for (_k, (a_pb, a_li)) in mpbli {
        let mut a_cb = if ds.is_common_block(a_pb) {
            ds.common_block(a_pb)
                .cloned()
                .unwrap_or_else(BopdsCommonBlock::new)
        } else {
            let mut cb = BopdsCommonBlock::new();
            cb.add_pave_block(a_pb.clone());
            cb
        };
        let mut a_new_faces: Vec<usize> = Vec::new();
        let old = a_cb.faces().to_vec();
        for &n_f in a_li {
            if !old.contains(&n_f) {
                a_new_faces.push(n_f);
            }
        }
        a_cb.append_faces(&a_new_faces);
        ds.set_common_block(a_pb, a_cb.clone());
        if let Some(cb) = ds.common_block(a_pb).cloned() {
            let a_tol = compute_tolerance_of_cb(&cb, ds, ctx);
            if let Some(cbm) = ds.common_block(a_pb).cloned() {
                let mut cbm = cbm;
                cbm.set_tolerance(a_tol);
                ds.set_common_block(a_pb, cbm);
            }
        }
    }
}

/// `BOPAlgo_Tools::ComputeToleranceOfCB` (`BOPAlgo_Tools.cxx:248`).
pub fn compute_tolerance_of_cb(
    the_cb: &BopdsCommonBlock,
    the_ds: &BopdsDS,
    the_ctx: &IntToolsContext,
) -> f64 {
    let Some(a_pbr) = the_cb.pave_block1() else {
        return 0.0;
    };
    let n_e = a_pbr.original_edge();
    let Some(e_or_shape) = the_ds.shape(n_e).cloned() else {
        return 0.0;
    };
    let a_e_or = Edge(e_or_shape);
    let mut a_tol_max = BRepTool::edge_tolerance(&a_e_or);
    if the_cb.pave_blocks().len() < 2 && the_cb.faces().is_empty() {
        return a_tol_max;
    }
    let Some(a_c3d) = BRepTool::edge_curve(&a_e_or) else {
        return a_tol_max;
    };
    let (a_t1, a_t2) = a_pbr.range();
    const A_NB_PNT: i32 = 11;
    let a_dt = (a_t2 - a_t1) / (A_NB_PNT as f64 + 1.0);
    if the_cb.pave_blocks().len() > 1 {
        for a_pb in the_cb.pave_blocks() {
            if pb_key(a_pb) == pb_key(a_pbr) && a_pb.original_edge() == a_pbr.original_edge() {
                continue;
            }
            let n_e2 = a_pb.original_edge();
            let Some(e_shape) = the_ds.shape(n_e2).cloned() else {
                continue;
            };
            let a_e = Edge(e_shape);
            let a_tol = BRepTool::edge_tolerance(&a_e);
            let mut a_t = a_t1;
            for _i in 1..=A_NB_PNT {
                a_t += a_dt;
                let a_p = a_c3d.d0(a_t);
                if let Some(t_proj) = the_ctx.project_point_on_edge(&a_e, &a_p) {
                    if let Some(c2) = BRepTool::edge_curve(&a_e) {
                        let a_tol_new = a_tol + a_p.distance(&c2.d0(t_proj));
                        if a_tol_new > a_tol_max {
                            a_tol_max = a_tol_new;
                        }
                    }
                }
            }
        }
    }
    for &n_f in the_cb.faces() {
        let Some(f_shape) = the_ds.shape(n_f).cloned() else {
            continue;
        };
        let a_f = Face(f_shape);
        let a_tol = BRepTool::face_tolerance(&a_f);
        let mut a_t = a_t1;
        for _i in 1..=A_NB_PNT {
            a_t += a_dt;
            let a_p = a_c3d.d0(a_t);
            if let Ok((u, v)) = the_ctx.project_point_on_face(&a_f, &a_p) {
                if let Some(surf) = BRepTool::face_surface(&a_f) {
                    let a_tol_new = a_tol + a_p.distance(&surf.d0(u, v));
                    if a_tol_new > a_tol_max {
                        a_tol_max = a_tol_new;
                    }
                }
            }
        }
    }
    a_tol_max
}

/// One vertex of the `IntersectVertices` input map, with its extra tolerance.
#[derive(Clone)]
pub struct VertexTol {
    pub shape: TopoShape,
    pub extra_tol: f64,
}

/// `BOPAlgo_Tools::IntersectVertices` (`BOPAlgo_Tools.cxx:1119`).
///
/// The BVH pair selector is a linear AABB + distance scan. Chains are the
/// connected components of interfering vertices; isolated vertices form
/// single-element chains.
pub fn intersect_vertices(the_vertices: &[VertexTol], the_fuzzy: f64) -> Vec<Vec<TopoShape>> {
    let a_nb_v = the_vertices.len();
    if a_nb_v <= 1 {
        return the_vertices
            .iter()
            .map(|v| vec![v.shape.clone()])
            .collect();
    }
    let a_tol_add = the_fuzzy / 2.0;
    let mut boxes: Vec<BndBox> = Vec::with_capacity(a_nb_v);
    let mut tols: Vec<f64> = Vec::with_capacity(a_nb_v);
    for v in the_vertices {
        let a_v = Vertex(v.shape.clone());
        let mut a_tol = BRepTool::vertex_tolerance(&a_v);
        if a_tol < v.extra_tol {
            a_tol = v.extra_tol;
        }
        let mut a_box = BndBox::new();
        a_box.add_point(&BRepTool::vertex_point(&a_v));
        a_box.enlarge(a_tol + a_tol_add);
        boxes.push(a_box);
        tols.push(a_tol);
    }
    let mut mili: HashMap<usize, Vec<usize>> = HashMap::new();
    let mut tree = crate::bop_aabb_faces::AabbTree::new();
    tree.set_size(a_nb_v);
    for (i, box_) in boxes.iter().enumerate() {
        tree.add(i, crate::bnd_tools::bnd2bvh3d(box_));
    }
    tree.build();
    for i in 0..a_nb_v {
        for j in tree.select(&boxes[i]) {
            if j <= i {
                continue;
            }
            let a_v1 = Vertex(the_vertices[i].shape.clone());
            let a_v2 = Vertex(the_vertices[j].shape.clone());
            let a_p1 = BRepTool::vertex_point(&a_v1);
            let a_p2 = BRepTool::vertex_point(&a_v2);
            let a_sum = tols[i] + tols[j] + the_fuzzy;
            if a_p1.distance(&a_p2) <= a_sum.max(CONFUSION) {
                fill_map_pair(i, j, &mut mili);
            }
        }
    }
    let mut chains: Vec<Vec<TopoShape>> = Vec::new();
    let blocks = make_blocks_int(&mili);
    let mut used: HashSet<usize> = HashSet::new();
    for a_li in blocks {
        let mut a_chain = Vec::new();
        for i in &a_li {
            used.insert(*i);
            a_chain.push(the_vertices[*i].shape.clone());
        }
        chains.push(a_chain);
    }
    for i in 0..a_nb_v {
        if !used.contains(&i) && !mili.contains_key(&i) {
            chains.push(vec![the_vertices[i].shape.clone()]);
        }
    }
    chains
}

/// `BOPAlgo_Tools::TrsfToPoint` (`BOPAlgo_Tools.cxx:1912`).
///
/// Returns a translation that moves the unified AABB of `box1`/`box2` toward
/// `point` (default origin) when the objects are far from that point relative
/// to `criteria` (default `1e5`).
pub fn trsf_to_point(
    box1: &BndBox,
    box2: &BndBox,
    point: Option<GpPnt>,
    criteria: Option<f64>,
) -> Option<GpTrsf> {
    let the_point = point.unwrap_or_else(GpPnt::zero);
    let the_criteria = criteria.unwrap_or(1.0e5);
    let mut a_box = *box1;
    a_box.add_box(box2);
    let cmin = a_box.corner_min();
    let cmax = a_box.corner_max();
    let cx = (cmin.x() + cmax.x()) * 0.5;
    let cy = (cmin.y() + cmax.y()) * 0.5;
    let cz = (cmin.z() + cmax.z()) * 0.5;
    let dx = the_point.x() - cx;
    let dy = the_point.y() - cy;
    let dz = the_point.z() - cz;
    let a_pb_dist = (dx * dx + dy * dy + dz * dz).sqrt();
    if a_pb_dist < the_criteria {
        return None;
    }
    let a_b_size = a_box.square_extent().sqrt();
    if a_b_size / a_pb_dist > 1.0 / the_criteria {
        return None;
    }
    let mut trsf = GpTrsf::identity();
    trsf.set_translation_pnts(&cmin, &the_point);
    Some(trsf)
}
