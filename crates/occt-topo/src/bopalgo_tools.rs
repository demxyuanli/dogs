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
pub fn fill_map_pb_pb(pb1: &BopdsPaveBlock, pb2: &BopdsPaveBlock, mpblpb: &mut PbPbListMap) {
    mpblpb
        .entry(pb_key(pb1))
        .or_insert_with(|| (pb1.clone(), Vec::new()))
        .1
        .push(pb2.clone());
    mpblpb
        .entry(pb_key(pb2))
        .or_insert_with(|| (pb2.clone(), Vec::new()))
        .1
        .push(pb1.clone());
}

/// `BOPAlgo_Tools::FillMap` (pave block, face).
pub fn fill_map_pb_face(pb: &BopdsPaveBlock, n_f: usize, mpbli: &mut PbFaceListMap) {
    let e = mpbli
        .entry(pb_key(pb))
        .or_insert_with(|| (pb.clone(), Vec::new()));
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
    for lpb in blocks {
        if lpb.len() < 2 {
            continue;
        }
        let mut m_faces: HashSet<usize> = HashSet::new();
        let mut l_faces: Vec<usize> = Vec::new();
        let mut cb: Option<BopdsCommonBlock> = None;
        for pb in &lpb {
            if ds.is_common_block(pb) {
                if let Some(cbx) = ds.common_block(pb) {
                    for &n_f in cbx.faces() {
                        if m_faces.insert(n_f) {
                            l_faces.push(n_f);
                        }
                    }
                    if cb.is_none() {
                        cb = Some(cbx.clone());
                    }
                }
            }
        }
        let mut cb = cb.unwrap_or_else(BopdsCommonBlock::new);
        cb.set_pave_blocks(lpb.clone());
        cb.set_faces(l_faces);
        let tol = compute_tolerance_of_cb(&cb, ds, ctx);
        cb.set_tolerance(tol);
        for pb in &lpb {
            ds.set_common_block(pb, cb.clone());
        }
    }
}

/// `BOPAlgo_Tools::PerformCommonBlocks` (pave-block / face-list overload).
pub fn perform_common_blocks_faces(
    mpbli: &PbFaceListMap,
    ds: &mut BopdsDS,
    ctx: &IntToolsContext,
) {
    for (_k, (pb, li)) in mpbli {
        let mut a_cb = if ds.is_common_block(pb) {
            ds.common_block(pb)
                .cloned()
                .unwrap_or_else(BopdsCommonBlock::new)
        } else {
            let mut cb = BopdsCommonBlock::new();
            cb.add_pave_block(pb.clone());
            cb
        };
        let mut new_faces: Vec<usize> = Vec::new();
        let old = a_cb.faces().to_vec();
        for &n_f in li {
            if !old.contains(&n_f) {
                new_faces.push(n_f);
            }
        }
        a_cb.append_faces(&new_faces);
        ds.set_common_block(pb, a_cb.clone());
        if let Some(cb) = ds.common_block(pb).cloned() {
            let tol = compute_tolerance_of_cb(&cb, ds, ctx);
            if let Some(cbm) = ds.common_block(pb).cloned() {
                let mut cbm = cbm;
                cbm.set_tolerance(tol);
                ds.set_common_block(pb, cbm);
            }
        }
    }
}

/// `BOPAlgo_Tools::ComputeToleranceOfCB` (`BOPAlgo_Tools.cxx:248`).
pub fn compute_tolerance_of_cb(
    cb: &BopdsCommonBlock,
    ds: &BopdsDS,
    ctx: &IntToolsContext,
) -> f64 {
    let Some(pbr) = cb.pave_block1() else {
        return 0.0;
    };
    let n_e = pbr.original_edge();
    let Some(e_or_shape) = ds.shape(n_e).cloned() else {
        return 0.0;
    };
    let e_or = Edge(e_or_shape);
    let mut tol_max = BRepTool::edge_tolerance(&e_or);
    if cb.pave_blocks().len() < 2 && cb.faces().is_empty() {
        return tol_max;
    }
    let Some(c3d) = BRepTool::edge_curve(&e_or) else {
        return tol_max;
    };
    let (t1, t2) = pbr.range();
    const A_NB_PNT: i32 = 11;
    let dt = (t2 - t1) / (A_NB_PNT as f64 + 1.0);
    if cb.pave_blocks().len() > 1 {
        for pb in cb.pave_blocks() {
            if pb_key(pb) == pb_key(pbr) && pb.original_edge() == pbr.original_edge() {
                continue;
            }
            let n_e2 = pb.original_edge();
            let Some(e_shape) = ds.shape(n_e2).cloned() else {
                continue;
            };
            let e = Edge(e_shape);
            let tol = BRepTool::edge_tolerance(&e);
            let mut t = t1;
            for _i in 1..=A_NB_PNT {
                t += dt;
                let p = c3d.d0(t);
                if let Some(t_proj) = ctx.project_point_on_edge(&e, &p) {
                    if let Some(c2) = BRepTool::edge_curve(&e) {
                        let tol_new = tol + p.distance(&c2.d0(t_proj));
                        if tol_new > tol_max {
                            tol_max = tol_new;
                        }
                    }
                }
            }
        }
    }
    for &n_f in cb.faces() {
        let Some(f_shape) = ds.shape(n_f).cloned() else {
            continue;
        };
        let f = Face(f_shape);
        let tol = BRepTool::face_tolerance(&f);
        let mut t = t1;
        for _i in 1..=A_NB_PNT {
            t += dt;
            let p = c3d.d0(t);
            if let Ok((u, v)) = ctx.project_point_on_face(&f, &p) {
                if let Some(surf) = BRepTool::face_surface(&f) {
                    let tol_new = tol + p.distance(&surf.d0(u, v));
                    if tol_new > tol_max {
                        tol_max = tol_new;
                    }
                }
            }
        }
    }
    tol_max
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
pub fn intersect_vertices(vertices: &[VertexTol], fuzzy: f64) -> Vec<Vec<TopoShape>> {
    let nb_v = vertices.len();
    if nb_v <= 1 {
        return vertices
            .iter()
            .map(|v| vec![v.shape.clone()])
            .collect();
    }
    let tol_add = fuzzy / 2.0;
    let mut boxes: Vec<BndBox> = Vec::with_capacity(nb_v);
    let mut tols: Vec<f64> = Vec::with_capacity(nb_v);
    for v in vertices {
        let a_v = Vertex(v.shape.clone());
        let mut tol = BRepTool::vertex_tolerance(&a_v);
        if tol < v.extra_tol {
            tol = v.extra_tol;
        }
        let mut bnd_box = BndBox::new();
        bnd_box.add_point(&BRepTool::vertex_point(&a_v));
        bnd_box.enlarge(tol + tol_add);
        boxes.push(bnd_box);
        tols.push(tol);
    }
    let mut mili: HashMap<usize, Vec<usize>> = HashMap::new();
    let mut tree = crate::bop_aabb_faces::AabbTree::new();
    tree.set_size(nb_v);
    for (i, box_) in boxes.iter().enumerate() {
        tree.add(i, crate::bnd_tools::bnd2bvh3d(box_));
    }
    tree.build();
    for i in 0..nb_v {
        for j in tree.select(&boxes[i]) {
            if j <= i {
                continue;
            }
            let v1 = Vertex(vertices[i].shape.clone());
            let v2 = Vertex(vertices[j].shape.clone());
            let p1 = BRepTool::vertex_point(&v1);
            let p2 = BRepTool::vertex_point(&v2);
            let sum = tols[i] + tols[j] + fuzzy;
            if p1.distance(&p2) <= sum.max(CONFUSION) {
                fill_map_pair(i, j, &mut mili);
            }
        }
    }
    let mut chains: Vec<Vec<TopoShape>> = Vec::new();
    let blocks = make_blocks_int(&mili);
    let mut used: HashSet<usize> = HashSet::new();
    for li in blocks {
        let mut chain = Vec::new();
        for i in &li {
            used.insert(*i);
            chain.push(vertices[*i].shape.clone());
        }
        chains.push(chain);
    }
    for i in 0..nb_v {
        if !used.contains(&i) && !mili.contains_key(&i) {
            chains.push(vec![vertices[i].shape.clone()]);
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
    let mut bnd_box = *box1;
    bnd_box.add_box(box2);
    let cmin = bnd_box.corner_min();
    let cmax = bnd_box.corner_max();
    let cx = (cmin.x() + cmax.x()) * 0.5;
    let cy = (cmin.y() + cmax.y()) * 0.5;
    let cz = (cmin.z() + cmax.z()) * 0.5;
    let dx = the_point.x() - cx;
    let dy = the_point.y() - cy;
    let dz = the_point.z() - cz;
    let pb_dist = (dx * dx + dy * dy + dz * dz).sqrt();
    if pb_dist < the_criteria {
        return None;
    }
    let b_size = bnd_box.square_extent().sqrt();
    if b_size / pb_dist > 1.0 / the_criteria {
        return None;
    }
    let mut trsf = GpTrsf::identity();
    trsf.set_translation_pnts(&cmin, &the_point);
    Some(trsf)
}
