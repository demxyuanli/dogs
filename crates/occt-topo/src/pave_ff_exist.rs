//! ON/IN `IsExistingPaveBlock`, both `ProcessExistingPaveBlocks` overloads,
//! and `CorrectToleranceOfSE`.
//!
//! Source: `BOPAlgo_PaveFiller_6.cxx` (`IsExistingPaveBlock` at 2047,
//! `ProcessExistingPaveBlocks` at 3072 and 3171, `PreparePostTreatFF` at 3609,
//! `CorrectToleranceOfSE` at 4072).
//! The ON/IN tree is a linear AABB scan (`BoxTree` without BVH).
//! T-97: items below are faithful ports of the named OCCT source, but their
//! OCCT-side consumers are not all ported yet, so parts are not called from this
//! crate. The `dead_code` allowance is deliberate: **pending wiring**, not dead
//! code. Do not delete them to silence warnings (see
//! specs/_a3n00_gap_analysis.md §9.309/§9.310); wire the consumer instead.
#![allow(dead_code)]

use std::collections::{HashMap, HashSet};

use occt_core::bnd::BndBox;
use occt_core::precision::RESOLUTION;
use occt_geom::Curve;

use crate::bbox_from_geometry::shape_bbox;
use crate::bopds::{BopdsDS, BopdsPave, BopdsPaveBlock};

use crate::brep_tool::BRepTool;
use crate::int_tools_full::IntToolsContext;
use crate::pave_ff::CoupleOfPaveBlocks;
use crate::pave_filler::{EdgeRangeDistance, PaveFiller};
use crate::pave_intersect::vertex_on_edge;
use crate::shape::{Edge, TopoShape, Vertex};
use crate::tgeometry::GeometryRegistry;

pub(crate) type PbKey = (usize, u64, u64);

pub(crate) fn pb_key(pb: &BopdsPaveBlock) -> PbKey {
    (pb.edge(), pb.first.to_bits(), pb.last.to_bits())
}

fn tuples_eq(a: (usize, f64, f64), b: (usize, f64, f64)) -> bool {
    a.0 == b.0 && (a.1 - b.1).abs() <= 1e-7 && (a.2 - b.2).abs() <= 1e-7
}

pub(crate) fn pb_in_face(
    on: &[(usize, f64, f64)],
    inn: &[(usize, f64, f64)],
    pb: &BopdsPaveBlock,
) -> bool {
    let t = (pb.edge(), pb.first, pb.last);
    on.iter().chain(inn.iter()).copied().any(|x| tuples_eq(x, t))
}

pub(crate) fn common_on_in(
    on1: &[(usize, f64, f64)],
    in1: &[(usize, f64, f64)],
    on2: &[(usize, f64, f64)],
    in2: &[(usize, f64, f64)],
) -> Vec<(usize, f64, f64)> {
    let in_f2 = |t: (usize, f64, f64)| {
        on2.iter()
            .chain(in2.iter())
            .copied()
            .any(|x| tuples_eq(x, t))
    };
    on1.iter()
        .chain(in1.iter())
        .copied()
        .filter(|&t| in_f2(t))
        .collect()
}

fn pb_in_common(common: &[(usize, f64, f64)], pb: &BopdsPaveBlock) -> bool {
    let t = (pb.edge(), pb.first, pb.last);
    common.iter().copied().any(|x| tuples_eq(x, t))
}

fn original_edge_degenerated(ds: &BopdsDS, pb: &BopdsPaveBlock) -> bool {
    let n = if pb.original_edge() != 0 {
        pb.original_edge()
    } else {
        pb.edge()
    };
    ds.shape(n)
        .map(|s| BRepTool::is_degenerated(&Edge(s.clone())))
        .unwrap_or(false)
}

/// Resolve FaceInfo ON/IN tuples to DS pave blocks (`theMPBOnIn`).
pub(crate) fn resolve_on_in_pbs(
    ds: &BopdsDS,
    tuples: &[(usize, f64, f64)],
) -> Vec<BopdsPaveBlock> {
    let mut out = Vec::new();
    for &(e, first, last) in tuples {
        if e == 0 {
            continue;
        }
        let mut found = false;
        for pb in ds.pave_blocks(e) {
            if (pb.first - first).abs() > 1e-7 || (pb.last - last).abs() > 1e-7 {
                continue;
            }
            if pb.edge() == 0 {
                continue;
            }
            if original_edge_degenerated(ds, pb) {
                continue;
            }
            out.push(pb.clone());
            found = true;
        }
        if found {
            continue;
        }
        if ds.shape(e)
            .map(|s| BRepTool::is_degenerated(&Edge(s.clone())))
            .unwrap_or(false)
        {
            continue;
        }
        let mut pb = BopdsPaveBlock::new();
        pb.set_edge(e);
        pb.set_original_edge(e);
        pb.set_range(first, last);
        if let Some(si) = ds.shape_info(e) {
            let verts: Vec<usize> = si
                .sub_shapes()
                .iter()
                .copied()
                .filter(|&s| {
                    ds.shape_info(s)
                        .map(|x| x.shape_type() == crate::abs::ShapeType::Vertex)
                        .unwrap_or(false)
                })
                .collect();
            if verts.len() >= 2 {
                pb.set_pave1(BopdsPave::new(verts[0], first));
                pb.set_pave2(BopdsPave::new(verts[1], last));
            }
        }
        out.push(pb);
    }
    out
}

fn is_geom_line(c: &dyn Curve) -> bool {
    let (a, b) = (c.first_parameter(), c.last_parameter());
    let samples: [f64; 3] = if a.is_finite() && b.is_finite() && b > a {
        [a, 0.5 * (a + b), b]
    } else {
        [-1.0, 0.0, 1.0]
    };
    let t0 = c.d1(samples[0]).1;
    if t0.square_magnitude() < RESOLUTION * RESOLUTION {
        return false;
    }
    for &u in &samples[1..] {
        let t = c.d1(u).1;
        if t0.cross_magnitude(&t) > 1e-7 * t0.magnitude() * t.magnitude() {
            return false;
        }
    }
    true
}

fn edge_box_out(ds: &BopdsDS, n_e: usize, other: &BndBox) -> bool {
    ds.box_of(n_e).map(|b| b.is_out_box(other)).unwrap_or(false)
}

/// `BOPAlgo_PaveFiller::IsExistingPaveBlock` ON/IN tree overload (`_6.cxx:2047`).
pub(crate) use crate::pave_ff_exist_onin::is_existing_pave_block_on_in;

/// `BOPAlgo_PaveFiller::PreparePostTreatFF` (`_6.cxx:3609`).
pub(crate) fn prepare_post_treat_ff(
    ds: &BopdsDS,
    int: usize,
    cur: usize,
    pb: &BopdsPaveBlock,
    mscpb: &mut Vec<(TopoShape, CoupleOfPaveBlocks)>,
    lpbc: &mut Vec<BopdsPaveBlock>,
) {
    lpbc.push(pb.clone());
    let Some(e) = ds.shape(pb.edge()) else { return };
    mscpb.push((
        e.clone(),
        CoupleOfPaveBlocks {
            index_interf: int,
            index: cur,
            pb: Some(pb.clone()),
        },
    ));
}

/// First `ProcessExistingPaveBlocks` overload (`_6.cxx:3072`).
pub(crate) fn process_existing_pave_blocks(
    ds: &BopdsDS,
    int: usize,
    cur: usize,
    n_f1: usize,
    n_f2: usize,
    es: &TopoShape,
    on_in: &[BopdsPaveBlock],
    on1: &[(usize, f64, f64)],
    in1: &[(usize, f64, f64)],
    on2: &[(usize, f64, f64)],
    in2: &[(usize, f64, f64)],
    distances: &HashMap<(usize, usize), Vec<EdgeRangeDistance>>,
    mscpb: &mut Vec<(TopoShape, CoupleOfPaveBlocks)>,
    lpbc: &mut Vec<BopdsPaveBlock>,
    mpb_add: &mut HashSet<PbKey>,
    pb_faces: &mut HashMap<PbKey, Vec<usize>>,
) {
    let box_es = shape_bbox(es);
    let tol_es = BRepTool::edge_tolerance(&Edge(es.clone()));
    for pbf in on_in {
        if pbf.edge() == 0 {
            continue;
        }
        if edge_box_out(ds, pbf.edge(), &box_es) {
            continue;
        }
        if mpb_add.contains(&pb_key(pbf)) {
            continue;
        }
        let b_in_f1 = pb_in_face(on1, in1, pbf);
        let b_in_f2 = pb_in_face(on2, in2, pbf);
        if b_in_f1 && b_in_f2 {
            mpb_add.insert(pb_key(pbf));
            prepare_post_treat_ff(ds, int, cur, pbf, mscpb, lpbc);
            continue;
        }
        let n_f = if b_in_f1 { n_f2 } else { n_f1 };
        let Some(p_list) = distances.get(&(pbf.original_edge(), n_f)) else {
            continue;
        };
        let (t1, t2) = pbf.range();
        let mut dist = f64::MAX;
        for range_val in p_list {
            if (t1 <= range_val.first && range_val.first <= t2)
                || (t1 <= range_val.last && range_val.last <= t2)
                || (range_val.first <= t1 && t1 <= range_val.last)
                || (range_val.first <= t2 && t2 <= range_val.last)
            {
                dist = range_val.distance;
                break;
            }
        }
        if dist >= f64::MAX {
            continue;
        }
        let tol_ef = ds
            .shape(pbf.edge())
            .map(|s| BRepTool::edge_tolerance(&Edge(s.clone())))
            .unwrap_or(0.0);
        if dist <= tol_es + tol_ef {
            mpb_add.insert(pb_key(pbf));
            prepare_post_treat_ff(ds, int, cur, pbf, mscpb, lpbc);
            let faces = pb_faces.entry(pb_key(pbf)).or_default();
            if !faces.contains(&n_f) {
                faces.push(n_f);
            }
        }
    }
}

/// Bound-vertex `ProcessExistingPaveBlocks` overload (`_6.cxx:3171`).
pub(crate) fn process_existing_pave_blocks_bound(
    f: &mut PaveFiller,
    int: usize,
    n_f1: usize,
    n_f2: usize,
    on_in: &[BopdsPaveBlock],
    dmbv: &HashMap<usize, Vec<usize>>,
    mscpb: &mut Vec<(TopoShape, CoupleOfPaveBlocks)>,
    mpb_add: &mut HashSet<PbKey>,
    pb_faces: &mut HashMap<PbKey, Vec<usize>>,
) {
    if dmbv.is_empty() {
        return;
    }
    let fuzzy = f.fuzzy_value();
    let (on1, in1) = match f.ds().face_info(n_f1) {
        Some(fi) => (fi.paves_on().to_vec(), fi.paves_in().to_vec()),
        None => (Vec::new(), Vec::new()),
    };
    let (on2, in2) = match f.ds().face_info(n_f2) {
        Some(fi) => (fi.paves_on().to_vec(), fi.paves_in().to_vec()),
        None => (Vec::new(), Vec::new()),
    };
    let ctx = f.context().clone();
    for (&i_c, lbv) in dmbv {
        let mut kept_extra: Vec<BopdsPaveBlock> = Vec::new();
        for &n_v in lbv {
            let Some(v_shape) = f.ds().shape(n_v).cloned() else {
                continue;
            };
            if let Some(bv) = f.ds().box_of(n_v) {
                let _ = bv;
            }
            for pb in on_in {
                if pb.pave1().index() == n_v || pb.pave2().index() == n_v {
                    continue;
                }
                if mpb_add.contains(&pb_key(pb)) {
                    continue;
                }
                let n_e = pb.edge();
                if n_e == 0 {
                    continue;
                }
                let Some(e_shape) = f.ds().shape(n_e).cloned() else {
                    continue;
                };
                let Ok(Some((_, _))) =
                    vertex_on_edge(&ctx, &Vertex(v_shape.clone()), &Edge(e_shape), fuzzy)
                else {
                    continue;
                };
                mpb_add.insert(pb_key(pb));
                prepare_post_treat_ff(f.ds(), int, i_c, pb, mscpb, &mut kept_extra);
                let b_in_f1 = pb_in_face(&on1, &in1, pb);
                let b_in_f2 = pb_in_face(&on2, &in2, pb);
                if !b_in_f1 || !b_in_f2 {
                    let n_f = if b_in_f1 { n_f2 } else { n_f1 };
                    let faces = pb_faces.entry(pb_key(pb)).or_default();
                    if !faces.contains(&n_f) {
                        faces.push(n_f);
                    }
                }
            }
        }
        if !kept_extra.is_empty() {
            if let Some(lpbc) = f
                .ds_mut()
                .interf_ff_mut()
                .get_mut(int)
                .and_then(|ff| ff.change_curves().get_mut(i_c))
                .map(|nc| nc.change_pave_blocks())
            {
                lpbc.extend(kept_extra);
            }
        }
    }
}

/// `BOPAlgo_PaveFiller::CorrectToleranceOfSE` (`BOPAlgo_PaveFiller_6.cxx:4072`).
/// Drops section-edge tolerance back to the curve's valid tolerance when
/// PostTreatFF had raised it to cover a large tangential zone, then tries to
/// shrink the connected vertices. Edge write is `BRep_TEdge::Tolerance` only
/// (does not lift vertices).
pub(crate) fn correct_tolerance_of_se(f: &mut PaveFiller) {
    let mut mvi_pbs: HashMap<usize, Vec<BopdsPaveBlock>> = HashMap::new();
    let mut mvi_to_reduce: HashSet<usize> = HashSet::new();
    let mut mpb: HashSet<PbKey> = HashSet::new();

    let nb_ff = f.ds().interf_ff().len();
    for i in 0..nb_ff {
        let nb_c = f.ds().interf_ff()[i].curves().len();
        for k in 0..nb_c {
            let tol_c = f.ds().interf_ff()[i].curves()[k].tolerance();
            let tol_tang = f.ds().interf_ff()[i].curves()[k].tangential_tolerance();
            let mut j = 0;
            while j < f.ds().interf_ff()[i].curves()[k].pave_blocks().len() {
                let pb = f.ds().interf_ff()[i].curves()[k].pave_blocks()[j].clone();
                if pb.edge() == 0 {
                    f.ds_mut()
                        .interf_ff_mut()[i]
                        .change_curves()[k]
                        .change_pave_blocks()
                        .remove(j);
                    continue;
                }
                if !mpb.insert(pb_key(&pb)) {
                    j += 1;
                    continue;
                }
                let mut b_is_reduced = false;
                if pb.original_edge() == 0 && tol_c < tol_tang {
                    let n_e = pb.edge();
                    if let Some(e_shape) = f.ds().shape(n_e).cloned() {
                        let tol_e = BRepTool::edge_tolerance(&Edge(e_shape.clone()));
                        if tol_c < tol_e {
                            if let Some(mut g) = GeometryRegistry::global().edge_geom(&e_shape) {
                                g.tolerance = tol_c;
                                GeometryRegistry::global().set_edge(&e_shape, g);
                            }
                            b_is_reduced = true;
                        }
                    }
                }
                for n_v_raw in [pb.pave1().index(), pb.pave2().index()] {
                    let n_v = f.ds().same_domain_index(n_v_raw);
                    mvi_pbs.entry(n_v).or_default().push(pb.clone());
                    if b_is_reduced {
                        mvi_to_reduce.insert(n_v);
                    }
                }
                j += 1;
            }
        }
    }

    if mvi_to_reduce.is_empty() {
        return;
    }

    let mut mvi_tol: HashMap<usize, f64> = HashMap::new();
    let pool: Vec<Vec<BopdsPaveBlock>> = f.ds().pave_blocks_pool().to_vec();
    for lpb in &pool {
        for pb in lpb {
            if pb.edge() == 0 {
                continue;
            }
            let Some(e_shape) = f.ds().shape(pb.edge()).cloned() else {
                continue;
            };
            let tol_e = BRepTool::edge_tolerance(&Edge(e_shape));
            let (n_v1, n_v2) = pb.indices();
            for n_v in [n_v1, n_v2] {
                if !mvi_to_reduce.contains(&n_v) {
                    continue;
                }
                mvi_tol
                    .entry(n_v)
                    .and_modify(|t| {
                        if tol_e > *t {
                            *t = tol_e;
                        }
                    })
                    .or_insert(tol_e);
                mvi_pbs.entry(n_v).or_default().push(pb.clone());
            }
        }
    }

    let keys: Vec<usize> = mvi_pbs.keys().copied().collect();
    for n_v in keys {
        if !mvi_to_reduce.contains(&n_v) {
            continue;
        }
        let Some(v_shape) = f.ds().shape(n_v).cloned() else {
            continue;
        };
        let v = Vertex(v_shape);
        let tol_v = BRepTool::vertex_tolerance(&v);
        let mut max_tol = mvi_tol.get(&n_v).copied().unwrap_or(0.0);
        if tol_v - max_tol < 0.001 * tol_v {
            continue;
        }
        let p = BRepTool::vertex_point(&v);
        let mut mpb_fence: HashSet<PbKey> = HashSet::new();
        let lpb = mvi_pbs.get(&n_v).cloned().unwrap_or_default();
        for pb in &lpb {
            if !mpb_fence.insert(pb_key(pb)) {
                continue;
            }
            let n_e = pb.edge();
            let Some(e_shape) = f.ds().shape(n_e).cloned() else {
                continue;
            };
            let e = Edge(e_shape);
            let Some(c) = BRepTool::edge_curve(&e) else {
                continue;
            };
            let tol_e = BRepTool::edge_tolerance(&e);
            for pave in [pb.pave1(), pb.pave2()] {
                let n_vsd = f.ds().same_domain_index(pave.index());
                if n_vsd != n_v {
                    continue;
                }
                let pon_e = c.d0(pave.parameter());
                let dist = p.distance(&pon_e) + tol_e;
                if dist > max_tol {
                    max_tol = dist;
                }
            }
        }
        if max_tol < tol_v {
            v.set_tolerance(max_tol);
        }
    }
}
