//! `BOPAlgo_PaveFiller::PostTreatFF` (`_6.cxx:1165-1669`).
//!
//! Nested fuse of the new section vertices/edges, unused stick vertices,
//! vertices of rejected section blocks, and micro-section bounding vertices.
//! Existing (already-DS) section edges are packed into a compound so they
//! do not intersect each other. The nested `PaveFiller` is not primary.
//!
//! After the fuse:
//! - vertex images become DS indices and, when they are not F/F points,
//!   same-domain substitutions;
//! - edge images become either a reused pave block (`aMEPB`) or a new
//!   section pave block on the originating F/F curve;
//! - common-block tolerance is cached (`aMCBTol`) and may raise the
//!   section-curve tolerance;
//! - a one-hop chain update of `aDMNewSD` is applied at the end.

use std::collections::{HashMap, HashSet};

use occt_core::precision::PCONFUSION;

use crate::abs::ShapeType;
use crate::bopalgo_tools::compute_tolerance_of_cb;
use crate::bopds::{BopdsDS, BopdsPave, BopdsPaveBlock, BopdsShapeInfo};
use crate::builder::TopoBuilder;
use crate::iterator::ShapeIterator;
use crate::pave_ff::CoupleOfPaveBlocks;
use crate::pave_ff_exist::pb_key;
use crate::pave_ff_misc::{get_stick_vertices, remove_used_vertices};
use crate::pave_ff_se::collect_micro_section_vertices;
use crate::pave_ff_update::DmExEdges;
use crate::pave_filler::PaveFiller;
use crate::shape::TopoShape;
use crate::tgeometry::GeometryRegistry;

fn pb_has_edge(pb: &BopdsPaveBlock) -> bool {
    pb.edge() != 0
}

fn shape_key(s: &TopoShape) -> usize {
    GeometryRegistry::shape_key(s)
}

/// Unused stick vertices (`_6.cxx:1203-1231`).
///
/// A vertex that survives `RemoveUsedVertices` on exactly one F/F pair is
/// kept. A vertex that appears unused on two or more pairs is dropped
/// (`IndMap.Add` then `VertsUnused.RemoveKey`).
pub fn verts_unused_post_treat(f: &PaveFiller) -> Vec<TopoShape> {
    let nb_ff = f.ds().interf_ff().len();
    let mut ind_map: HashSet<usize> = HashSet::new();
    let mut unused: HashMap<usize, TopoShape> = HashMap::new();
    for i in 0..nb_ff {
        let (n_f1, n_f2) = f.ds().interf_ff()[i].indices();
        let (mut mv, _mv_ef, _mi) = get_stick_vertices(f.ds(), n_f1, n_f2);
        remove_used_vertices(f.ds().interf_ff()[i].curves(), &mut mv);
        for n_v in mv {
            let Some(s) = f.ds().shape(n_v).cloned() else {
                continue;
            };
            if ind_map.insert(n_v) {
                unused.insert(n_v, s);
            } else {
                unused.remove(&n_v);
            }
        }
    }
    unused.into_values().collect()
}

fn register_one_vertex(f: &mut PaveFiller, shape: &TopoShape, cpb: &CoupleOfPaveBlocks) -> Result<(), String> {
    let i_v = f.ds_mut().append(shape.clone())?;
    let i_x = cpb.index_interf;
    let i_p = cpb.index;
    if let Some(np) = f
        .ds_mut()
        .interf_ff_mut()
        .get_mut(i_x)
        .and_then(|ff| ff.change_points().get_mut(i_p))
    {
        np.set_index(i_v);
    }
    Ok(())
}

fn register_one_edge(f: &mut PaveFiller, shape: &TopoShape, cpb: &CoupleOfPaveBlocks) -> Result<(), String> {
    let Some(pb_ref) = cpb.pb.as_ref() else {
        return Ok(());
    };
    if pb_has_edge(pb_ref) {
        return Ok(());
    }
    let (n_v1, n_v2) = pb_ref.indices();
    let (t1, t2) = pb_ref.range();
    let mut si = BopdsShapeInfo::new(shape.clone());
    si.change_sub_shapes().extend_from_slice(&[n_v1, n_v2]);
    let i_e = f.ds_mut().append_info(si);
    let i_x = cpb.index_interf;
    let i_c = cpb.index;
    if let Some(lpbc) = f
        .ds_mut()
        .interf_ff_mut()
        .get_mut(i_x)
        .and_then(|ff| ff.change_curves().get_mut(i_c))
        .map(|nc| nc.change_pave_blocks())
    {
        for pb in lpbc.iter_mut() {
            if pb.indices() == (n_v1, n_v2)
                && (pb.first - t1).abs() <= PCONFUSION
                && (pb.last - t2).abs() <= PCONFUSION
            {
                pb.set_edge(i_e);
            }
        }
    }
    Ok(())
}

fn register_one(f: &mut PaveFiller, shape: &TopoShape, cpb: &CoupleOfPaveBlocks) -> Result<(), String> {
    match shape.shape_type() {
        ShapeType::Vertex => register_one_vertex(f, shape, cpb),
        ShapeType::Edge => register_one_edge(f, shape, cpb),
        _ => Ok(()),
    }
}

fn fast_path_one(
    f: &mut PaveFiller,
    shape: &TopoShape,
    cpb: &CoupleOfPaveBlocks,
    dm_ex: &mut DmExEdges,
) -> Result<(), String> {
    match shape.shape_type() {
        ShapeType::Vertex => register_one_vertex(f, shape, cpb),
        ShapeType::Edge => {
            if let Some(pb1) = cpb.pb.as_ref() {
                if pb_has_edge(pb1) {
                    dm_ex.entry(pb_key(pb1)).or_default().push(pb1.clone());
                    return Ok(());
                }
            }
            register_one_edge(f, shape, cpb)
        }
        _ => Ok(()),
    }
}

fn import_index(main: &mut BopdsDS, nested: &BopdsDS, n: usize) -> Result<usize, String> {
    let n = nested.get_same_domain_index(n);
    let Some(s) = nested.shape(n).cloned() else {
        return Err(format!("post_treat_ff: nested shape {n} missing"));
    };
    if let Some(i) = main.index(&s) {
        return Ok(main.get_same_domain_index(i));
    }
    match s.shape_type() {
        ShapeType::Vertex => main.append(s),
        ShapeType::Edge => {
            let mut si = BopdsShapeInfo::new(s);
            if let Some(info) = nested.shape_info(n) {
                for &sub in info.sub_shapes() {
                    if nested.shape_info(sub).map(|x| x.shape_type()) == Some(ShapeType::Vertex) {
                        si.change_sub_shapes()
                            .push(import_index(main, nested, sub)?);
                    }
                }
            }
            Ok(main.append_info(si))
        }
        _ => main.append(s),
    }
}

fn remove_curve_pb(f: &mut PaveFiller, i_x: usize, i_c: usize, pb1: &BopdsPaveBlock) {
    let (n1, n2) = pb1.indices();
    let (t1, t2) = pb1.range();
    if let Some(lpbc) = f
        .ds_mut()
        .interf_ff_mut()
        .get_mut(i_x)
        .and_then(|ff| ff.change_curves().get_mut(i_c))
        .map(|nc| nc.change_pave_blocks())
    {
        lpbc.retain(|pb| {
            pb.indices() != (n1, n2)
                || (pb.first - t1).abs() > PCONFUSION
                || (pb.last - t2).abs() > PCONFUSION
        });
    }
}

fn bind_sd_from_pave(
    f: &mut PaveFiller,
    nested: &BopdsDS,
    pb1: &BopdsPaveBlock,
    pave_j: &BopdsPave,
    j: usize,
    i_v: usize,
    n_v_nested: usize,
    dm_new_sd: &mut HashMap<usize, usize>,
) {
    let a_p1 = if j == 0 { pb1.pave1() } else { pb1.pave2() };
    if a_p1.index() == i_v {
        return;
    }
    if (a_p1.parameter() - pave_j.parameter()).abs() <= PCONFUSION {
        dm_new_sd.insert(a_p1.index(), i_v);
        f.ds_mut().add_shape_sd(a_p1.index(), i_v);
        return;
    }
    let Some(a_v_pave) = f.ds().shape(a_p1.index()).cloned() else {
        return;
    };
    let Some(n_v_new) = nested.index(&a_v_pave) else {
        return;
    };
    if let Some(n_v_new_sd) = nested.has_shape_sd(n_v_new) {
        if n_v_new_sd == n_v_nested {
            dm_new_sd.insert(a_p1.index(), i_v);
            f.ds_mut().add_shape_sd(a_p1.index(), i_v);
        }
    }
}

fn apply_cb_tolerance(
    f: &mut PaveFiller,
    nested: &PaveFiller,
    pbr: &BopdsPaveBlock,
    i_x: usize,
    i_c: usize,
    mcb_tol: &mut HashMap<(usize, u64, u64), f64>,
) {
    if !nested.ds().is_common_block(pbr) {
        return;
    }
    let Some(cb) = nested.ds().common_block(pbr) else {
        return;
    };
    let key = (pbr.edge(), pbr.first.to_bits(), pbr.last.to_bits());
    let a_tol = *mcb_tol.entry(key).or_insert_with(|| {
        compute_tolerance_of_cb(&cb, nested.ds(), nested.context())
    });
    if let Some(nc) = f
        .ds_mut()
        .interf_ff_mut()
        .get_mut(i_x)
        .and_then(|ff| ff.change_curves().get_mut(i_c))
    {
        if nc.tolerance() < a_tol {
            nc.set_tolerance(a_tol);
        }
    }
}

fn reuse_or_make_pb(
    nested: &BopdsDS,
    pbx: &BopdsPaveBlock,
    pbr: &BopdsPaveBlock,
    i0: usize,
    i1: usize,
    i_e: usize,
    me_pb: &mut HashMap<usize, BopdsPaveBlock>,
) -> BopdsPaveBlock {
    if let Some(existing) = me_pb.get(&i_e) {
        return existing.clone();
    }
    let mut a_pave_r1 = pbr.pave1();
    let mut a_pave_r2 = pbr.pave2();
    let i_r1 = nested
        .shape(a_pave_r1.index())
        .and_then(|s| {
            let _ = s;
            Some(i0)
        })
        .unwrap_or(i0);
    let i_r2 = nested
        .shape(a_pave_r2.index())
        .and_then(|s| {
            let _ = s;
            Some(i1)
        })
        .unwrap_or(i1);
    let _ = pbx;
    a_pave_r1.set_index(i_r1);
    a_pave_r2.set_index(i_r2);
    let mut pb_new = BopdsPaveBlock::new();
    pb_new.set_pave1(a_pave_r1);
    pb_new.set_pave2(a_pave_r2);
    pb_new.set_edge(i_e);
    me_pb.insert(i_e, pb_new.clone());
    pb_new
}

fn process_nested_vertex(
    f: &mut PaveFiller,
    nested: &BopdsDS,
    a_sx: &TopoShape,
    n_sx: usize,
    couples: &HashMap<usize, usize>,
    mscpb: &[(TopoShape, CoupleOfPaveBlocks)],
    dm_new_sd: &mut HashMap<usize, usize>,
) -> Result<(), String> {
    let b_intersection_point = couples.contains_key(&shape_key(a_sx));
    let a_v = if let Some(sd) = nested.has_shape_sd(n_sx) {
        nested
            .shape(sd)
            .cloned()
            .ok_or_else(|| format!("post_treat_ff: nested SD {sd} missing"))?
    } else {
        a_sx.clone()
    };
    let i_v = if let Some(i) = f.ds().index(&a_v) {
        i
    } else {
        f.ds_mut().append_info(BopdsShapeInfo::new(a_v.clone()))
    };
    if !b_intersection_point {
        if let Some(n_main) = f.ds().index(a_sx) {
            if n_main != i_v {
                dm_new_sd.insert(n_main, i_v);
                f.ds_mut().add_shape_sd(n_main, i_v);
            }
        }
    } else if let Some(&ci) = couples.get(&shape_key(a_sx)) {
        let i_x = mscpb[ci].1.index_interf;
        let i_p = mscpb[ci].1.index;
        if let Some(np) = f
            .ds_mut()
            .interf_ff_mut()
            .get_mut(i_x)
            .and_then(|ff| ff.change_points().get_mut(i_p))
        {
            np.set_index(i_v);
        }
    }
    Ok(())
}

fn process_nested_edge(
    f: &mut PaveFiller,
    nested: &PaveFiller,
    a_sx: &TopoShape,
    n_sx: usize,
    couples: &HashMap<usize, usize>,
    mscpb: &[(TopoShape, CoupleOfPaveBlocks)],
    dm_new_sd: &mut HashMap<usize, usize>,
    dm_ex: &mut DmExEdges,
    me_pb: &mut HashMap<usize, BopdsPaveBlock>,
    mcb_tol: &mut HashMap<(usize, u64, u64), f64>,
    a_ls: &mut Vec<TopoShape>,
) -> Result<(), String> {
    let Some(&ci) = couples.get(&shape_key(a_sx)) else {
        return Ok(());
    };
    let cpb = &mscpb[ci].1;
    let Some(pb1) = cpb.pb.as_ref() else {
        return Ok(());
    };
    let b_old = pb_has_edge(pb1);
    let i_x = cpb.index_interf;
    let i_c = cpb.index;
    if b_old {
        dm_ex.entry(pb_key(pb1)).or_default();
    }
    let b_has_pave_blocks = nested.ds().has_pave_blocks(n_sx);
    if !b_has_pave_blocks {
        if b_old {
            dm_ex.entry(pb_key(pb1)).or_default().push(pb1.clone());
        } else {
            register_one_edge(f, a_sx, cpb)?;
        }
        return Ok(());
    }
    let lpbx = nested.ds().pave_blocks(n_sx).to_vec();
    let a_nb_lpbx = lpbx.len();
    let micro = a_nb_lpbx == 0 || (a_nb_lpbx == 1 && !lpbx[0].has_shrunk_data());
    if micro {
        if !b_old {
            remove_curve_pb(f, i_x, i_c, pb1);
        }
        for child in ShapeIterator::of_shape(a_sx) {
            if child.shape_type() == ShapeType::Vertex {
                a_ls.push(child);
            }
        }
        return Ok(());
    }
    if !b_old {
        remove_curve_pb(f, i_x, i_c, pb1);
    }
    for pbx in &lpbx {
        let pbr = nested.ds().real_pave_block(pbx);
        let (n0, n1) = pbx.indices();
        let i0 = import_index(f.ds_mut(), nested.ds(), n0)?;
        let i1 = import_index(f.ds_mut(), nested.ds(), n1)?;
        bind_sd_from_pave(f, nested.ds(), pb1, &pbx.pave1(), 0, i0, n0, dm_new_sd);
        bind_sd_from_pave(f, nested.ds(), pb1, &pbx.pave2(), 1, i1, n1, dm_new_sd);
        let i_e = import_index(f.ds_mut(), nested.ds(), pbr.edge())?;
        apply_cb_tolerance(f, nested, &pbr, i_x, i_c, mcb_tol);
        let mut pb_new = reuse_or_make_pb(nested.ds(), pbx, &pbr, i0, i1, i_e, me_pb);
        if b_old {
            pb_new.set_original_edge(pb1.original_edge());
            me_pb.insert(i_e, pb_new.clone());
            dm_ex.entry(pb_key(pb1)).or_default().push(pb_new);
        } else if let Some(lpbc) = f
            .ds_mut()
            .interf_ff_mut()
            .get_mut(i_x)
            .and_then(|ff| ff.change_curves().get_mut(i_c))
            .map(|nc| nc.change_pave_blocks())
        {
            me_pb.insert(i_e, pb_new.clone());
            lpbc.push(pb_new);
        }
    }
    Ok(())
}

fn chain_update_sd(f: &mut PaveFiller, dm_new_sd: &mut HashMap<usize, usize>) {
    let keys: Vec<usize> = dm_new_sd.keys().copied().collect();
    for k in keys {
        let Some(&v) = dm_new_sd.get(&k) else {
            continue;
        };
        if let Some(&p_sd) = dm_new_sd.get(&v) {
            dm_new_sd.insert(k, p_sd);
            f.ds_mut().add_shape_sd(k, p_sd);
        }
    }
}

fn prepare_arguments(
    f: &PaveFiller,
    mscpb: &[(TopoShape, CoupleOfPaveBlocks)],
    dm_new_sd: &HashMap<usize, usize>,
    _micro_pb: &[BopdsPaveBlock],
    _verts_on_rejected: &[TopoShape],
    _unused: &[TopoShape],
) -> (Vec<TopoShape>, HashSet<usize>) {
    let mut existing: Vec<TopoShape> = Vec::new();
    let mut a_ls: Vec<TopoShape> = Vec::new();
    let mut added_sd: HashSet<usize> = HashSet::new();
    for (shape, cpb) in mscpb.iter().rev() {
        if let Some(pb) = cpb.pb.as_ref() {
            if pb_has_edge(pb) {
                existing.push(shape.clone());
            } else {
                a_ls.push(shape.clone());
            }
        } else {
            a_ls.push(shape.clone());
        }
        for child in ShapeIterator::of_shape(shape) {
            if child.shape_type() != ShapeType::Vertex {
                continue;
            }
            let Some(i_ver) = f.ds().index(&child) else {
                continue;
            };
            let Some(&n_sd) = dm_new_sd.get(&i_ver) else {
                continue;
            };
            if added_sd.insert(n_sd) {
                if let Some(v_sd) = f.ds().shape(n_sd).cloned() {
                    a_ls.push(v_sd);
                }
            }
        }
    }
    if !existing.is_empty() {
        a_ls.push(TopoBuilder::new().make_compound_of(&existing).0);
    }
    (a_ls, added_sd)
}

fn append_rejected_and_unused(
    f: &PaveFiller,
    dm_new_sd: &HashMap<usize, usize>,
    added_sd: &mut HashSet<usize>,
    a_ls: &mut Vec<TopoShape>,
    verts_on_rejected: &[TopoShape],
    unused: &[TopoShape],
) {
    for src in [verts_on_rejected, unused] {
        for a_ver in src {
            let mut a_ver = a_ver.clone();
            if let Some(i_ver) = f.ds().index(&a_ver) {
                if let Some(&p_sd) = dm_new_sd.get(&i_ver) {
                    if let Some(s) = f.ds().shape(p_sd).cloned() {
                        a_ver = s;
                    }
                }
                let n = dm_new_sd.get(&i_ver).copied().unwrap_or(i_ver);
                if !added_sd.insert(n) {
                    continue;
                }
            } else if !added_sd.insert(shape_key(&a_ver)) {
                continue;
            }
            a_ls.push(a_ver);
        }
    }
}

/// `BOPAlgo_PaveFiller::PostTreatFF` (`_6.cxx:1165`).
pub fn post_treat_ff(
    f: &mut PaveFiller,
    mscpb: &[(TopoShape, CoupleOfPaveBlocks)],
    dm_new_sd: &mut HashMap<usize, usize>,
    dm_ex: &mut DmExEdges,
    micro_pb: &[BopdsPaveBlock],
    verts_on_rejected: &[TopoShape],
) -> Result<(), String> {
    if mscpb.is_empty() {
        return Ok(());
    }
    let unused = verts_unused_post_treat(f);
    if mscpb.len() == 1 && unused.is_empty() && micro_pb.is_empty() && verts_on_rejected.is_empty()
    {
        let (shape, cpb) = &mscpb[0];
        return fast_path_one(f, shape, cpb, dm_ex);
    }

    let (mut a_ls, mut added_sd) =
        prepare_arguments(f, mscpb, dm_new_sd, micro_pb, verts_on_rejected, &unused);

    collect_micro_section_vertices(f, micro_pb, dm_new_sd, &mut added_sd, &mut a_ls);
    append_rejected_and_unused(
        f,
        dm_new_sd,
        &mut added_sd,
        &mut a_ls,
        verts_on_rejected,
        &unused,
    );

    if a_ls.is_empty() {
        return Ok(());
    }

    let mut nested = PaveFiller::new();
    nested.set_is_primary(false);
    nested.set_non_destructive(f.non_destructive());
    nested.set_fuzzy_value(f.fuzzy_value());
    nested.set_arguments(&a_ls);
    if let Err(e) = nested.perform() {
        f.add_error("BOPAlgo_AlertPostTreatFF".into());
        return Err(e);
    }
    if nested.has_errors() {
        f.add_error("BOPAlgo_AlertPostTreatFF".into());
        return Err("BOPAlgo_AlertPostTreatFF".into());
    }

    let couples: HashMap<usize, usize> = mscpb
        .iter()
        .enumerate()
        .map(|(i, (s, _))| (shape_key(s), i))
        .collect();

    let mut me_pb: HashMap<usize, BopdsPaveBlock> = HashMap::new();
    let mut mcb_tol: HashMap<(usize, u64, u64), f64> = HashMap::new();

    let mut k = 0;
    while k < a_ls.len() {
        let sx = a_ls[k].clone();
        k += 1;
        if sx.shape_type() == ShapeType::Compound {
            for c in ShapeIterator::of_shape(&sx) {
                a_ls.push(c);
            }
            continue;
        }
        let Some(n_sx) = nested.ds().index(&sx) else {
            continue;
        };
        match sx.shape_type() {
            ShapeType::Vertex => {
                process_nested_vertex(f, nested.ds(), &sx, n_sx, &couples, mscpb, dm_new_sd)?;
            }
            ShapeType::Edge => {
                process_nested_edge(
                    f,
                    &nested,
                    &sx,
                    n_sx,
                    &couples,
                    mscpb,
                    dm_new_sd,
                    dm_ex,
                    &mut me_pb,
                    &mut mcb_tol,
                    &mut a_ls,
                )?;
            }
            _ => {}
        }
    }
    chain_update_sd(f, dm_new_sd);
    Ok(())
}
