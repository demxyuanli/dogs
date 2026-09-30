//! Face/face section-edge construction — `MakeBlocks` + `PostTreatFF`.
//!
//! Source: `BOPAlgo_PaveFiller_6.cxx`. `PerformFF` only writes `InterfFF`
//! curves; this module reads those curves, puts On/In/bound paves, builds
//! section edges, registers them on the DS, and fills `FaceInfo` Sc.
//! T-97: items below are faithful ports of the named OCCT source, but their
//! OCCT-side consumers are not all ported yet, so parts are not called from this
//! crate. The `dead_code` allowance is deliberate: **pending wiring**, not dead
//! code. Do not delete them to silence warnings (see
//! specs/_a3n00_gap_analysis.md §9.309/§9.310); wire the consumer instead.
#![allow(dead_code)]

use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use occt_core::gp::GpPnt;
use occt_core::precision::PCONFUSION;
use occt_geom::Curve;

use crate::abs::ShapeType;
use crate::algo_tools::AlgoTools;
use crate::bopds::{BopdsDS, BopdsPave, BopdsPaveBlock, BopdsShapeInfo};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::int_tools_full::IntToolsContext;
use crate::iterator::ShapeIterator;
use crate::pave_filler::{PaveFiller};
use crate::pave_intersect::make_sd_vertices;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::tgeometry::GeometryRegistry;

/// Couple of a post-treat shape with its originating F/F curve or point.
/// Source: `BOPDS_CoupleOfPaveBlocks`.
pub(crate) struct CoupleOfPaveBlocks {
    pub(crate) index_interf: usize,
    pub(crate) index: usize,
    pub(crate) pb: Option<BopdsPaveBlock>,
}

/// True when the pave block has a real DS edge (`BOPDS_PaveBlock::HasEdge`).
/// The empty constructor stores `edge_index == 0` as "unset".
fn pb_has_edge(pb: &BopdsPaveBlock) -> bool {
    pb.edge() != 0
}

/// `BOPAlgo_PaveFiller::MakeBlocks` (`BOPAlgo_PaveFiller_6.cxx`).
pub fn make_blocks_ff(f: &mut PaveFiller) -> Result<(), String> {
    crate::pave_ff_make::make_blocks_ff(f)
}

fn is_existing_vertex(ds: &BopdsDS, p: &GpPnt, tol: f64, mv_on_in: &HashSet<usize>) -> bool {
    let mut box_p = occt_core::bnd::BndBox::new();
    box_p.add_point(p);
    box_p.enlarge(tol);
    for &n_v in mv_on_in {
        if let Some(bv) = ds.box_of(n_v) {
            if box_p.is_out_box(bv) {
                continue;
            }
        }
        let Some(v) = ds.shape(n_v) else { continue };
        if AlgoTools::compute_vv(v, p, tol) != 0 {
            return true;
        }
    }
    false
}

fn stick_vertices(
    ds: &BopdsDS,
    n_f1: usize,
    n_f2: usize,
) -> (HashSet<usize>, HashSet<usize>, HashSet<usize>) {
    let mut mi = HashSet::new();
    full_shape_map(ds, n_f1, &mut mi);
    full_shape_map(ds, n_f2, &mut mi);
    let mut stick = HashSet::new();
    let mut ef = HashSet::new();
    let take = |arr: &[crate::bopds::BopdsInterf], stick: &mut HashSet<usize>| {
        for it in arr {
            if !it.has_index_new() {
                continue;
            }
            let (s1, s2) = it.indices();
            if mi.contains(&s1) && mi.contains(&s2) {
                let mut n = it.get_index_new().unwrap();
                n = ds.has_shape_sd(n).unwrap_or(n);
                stick.insert(n);
            }
        }
    };
    take(ds.interf_vv(), &mut stick);
    take(ds.interf_ve(), &mut stick);
    take(ds.interf_ee(), &mut stick);
    take(ds.interf_vf(), &mut stick);
    for it in ds.interf_ef() {
        if !it.has_index_new() {
            continue;
        }
        let (s1, s2) = it.indices();
        if mi.contains(&s1) && mi.contains(&s2) {
            let mut n = it.get_index_new().unwrap();
            n = ds.has_shape_sd(n).unwrap_or(n);
            stick.insert(n);
            ef.insert(n);
        }
    }
    (stick, ef, mi)
}

fn full_shape_map(ds: &BopdsDS, n_f: usize, mi: &mut HashSet<usize>) {
    mi.insert(n_f);
    if let Some(si) = ds.shape_info(n_f) {
        for &s in si.sub_shapes() {
            mi.insert(s);
        }
    }
}

fn is_valid_block_for_faces(
    ctx: &mut IntToolsContext,
    curve: &dyn Curve,
    t1: f64,
    t2: f64,
    f1: &Face,
    f2: &Face,
    tol: f64,
) -> bool {
    let samples = [t1, 0.5 * (t1 + t2), t2];
    for t in samples {
        let p = curve.d0(t);
        let ok1 = ctx.is_point_in_on_face(f1, &p, None, tol).unwrap_or(false);
        let ok2 = ctx.is_point_in_on_face(f2, &p, None, tol).unwrap_or(false);
        if !ok1 || !ok2 {
            return false;
        }
    }
    true
}

pub(crate) fn is_existing_pave_block_on_shared(
    ds: &BopdsDS,
    ctx: &IntToolsContext,
    pb: &BopdsPaveBlock,
    curve: &dyn Curve,
    lse: &[usize],
    fuzzy: f64,
) -> Option<(usize, f64)> {
    if lse.is_empty() {
        return None;
    }
    let (t1, t2) = pb.range();
    let tm = crate::boptools_2d::intermediate_point(t1, t2);
    let pm = curve.d0(tm);
    let (n_v1, n_v2) = pb.indices();
    let mut tol: f64 = 0.0;
    if let Some(v) = ds.shape(n_v1) {
        tol = tol.max(BRepTool::vertex_tolerance(&Vertex(v.clone())));
    }
    if let Some(v) = ds.shape(n_v2) {
        tol = tol.max(BRepTool::vertex_tolerance(&Vertex(v.clone())));
    }
    let mut box_pm = occt_core::bnd::BndBox::new();
    box_pm.add_point(&pm);
    box_pm.enlarge(tol);
    for &n_e in lse {
        if n_e == 0 {
            continue;
        }
        if let Some(be) = ds.box_of(n_e) {
            if be.is_out_box(&box_pm) {
                continue;
            }
        }
        let Some(e_shape) = ds.shape(n_e) else { continue };
        let e = Edge(e_shape.clone());
        let a_tol_e = BRepTool::edge_tolerance(&e);
        let a_tol_check = a_tol_e.max(tol) + fuzzy;
        let (flag, _, dist) = ctx.compute_pe_pnt(&pm, a_tol_check, &e);
        if flag == 0 {
            return Some((n_e, dist));
        }
    }
    None
}

fn is_micro_block(
    v1: &TopoShape,
    v2: &TopoShape,
    curve: &dyn Curve,
    t1: f64,
    t2: f64,
    tol_r3d: f64,
) -> bool {
    let p1 = curve.d0(t1);
    let p2 = curve.d0(t2);
    let tv1 = BRepTool::vertex_tolerance(&Vertex(v1.clone())).max(tol_r3d);
    let tv2 = BRepTool::vertex_tolerance(&Vertex(v2.clone())).max(tol_r3d);
    p1.distance(&p2) <= tv1 + tv2
}

fn make_pcurve_on_faces(
    edge: &TopoShape,
    f1: &Face,
    f2: &Face,
    pc1: Option<Arc<dyn occt_geom2d::curve::Curve2d>>,
    pc2: Option<Arc<dyn occt_geom2d::curve::Curve2d>>,
    on_s1: bool,
    on_s2: bool,
) {
    let e = Edge(edge.clone());
    let reg = GeometryRegistry::global();
    if on_s1 {
        if let Some(pc) = pc1 {
            let key = GeometryRegistry::shape_key(&f1.0);
            reg.set_edge_pcurve(edge, key, pc);
        } else if let Ok(pc) = AlgoTools::make_pcurve(&e, f1) {
            let key = GeometryRegistry::shape_key(&f1.0);
            reg.set_edge_pcurve(edge, key, pc);
        }
    }
    if on_s2 {
        if let Some(pc) = pc2 {
            let key = GeometryRegistry::shape_key(&f2.0);
            reg.set_edge_pcurve(edge, key, pc);
        } else if let Ok(pc) = AlgoTools::make_pcurve(&e, f2) {
            let key = GeometryRegistry::shape_key(&f2.0);
            reg.set_edge_pcurve(edge, key, pc);
        }
    }
}

fn make_sd_vertices_ff(
    f: &mut PaveFiller,
    dmvlv: &HashMap<usize, Vec<usize>>,
    dm_new_sd: &mut HashMap<usize, usize>,
) -> Result<(), String> {
    for list in dmvlv.values() {
        if list.len() < 2 {
            continue;
        }
        let n_sd = make_sd_vertices(f.ds_mut(), list, false)?;
        for &n_v in list {
            dm_new_sd.insert(n_v, n_sd);
        }
    }
    Ok(())
}

/// `BOPAlgo_PaveFiller::PostTreatFF` (`BOPAlgo_PaveFiller_6.cxx`).
pub(crate) fn post_treat_ff(
    f: &mut PaveFiller,
    mscpb: &[(TopoShape, CoupleOfPaveBlocks)],
    dm_new_sd: &mut HashMap<usize, usize>,
    dm_ex: &mut crate::pave_ff_update::DmExEdges,
    micro_pb: &[BopdsPaveBlock],
    verts_on_rejected: &[TopoShape],
) -> Result<(), String> {
    crate::pave_ff_post::post_treat_ff(f, mscpb, dm_new_sd, dm_ex, micro_pb, verts_on_rejected)
}

#[allow(dead_code)]
fn post_treat_ff_legacy(
    f: &mut PaveFiller,
    mscpb: &[(TopoShape, CoupleOfPaveBlocks)],
    dm_new_sd: &mut HashMap<usize, usize>,
    dm_ex: &mut crate::pave_ff_update::DmExEdges,
    micro_pb: &[BopdsPaveBlock],
) -> Result<(), String> {
    if mscpb.is_empty() {
        return Ok(());
    }
    let unused = unused_section_vertices(f);
    if mscpb.len() == 1 && unused.is_empty() && micro_pb.is_empty() {
        let (shape, cpb) = &mscpb[0];
        if shape.shape_type() == ShapeType::Edge {
            if let Some(pb1) = cpb.pb.as_ref() {
                if pb_has_edge(pb1) {
                    dm_ex.entry(crate::pave_ff_exist::pb_key(pb1)).or_default().push(pb1.clone());
                    return Ok(());
                }
            }
        }
        register_one(f, shape, cpb)?;
        return Ok(());
    }

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
            let Some(i_ver) = f.ds().index(&child) else { continue };
            let Some(&n_sd) = dm_new_sd.get(&i_ver) else { continue };
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
    for v in unused {
        let Some(i_ver) = f.ds().index(&v) else {
            a_ls.push(v);
            continue;
        };
        let n = dm_new_sd.get(&i_ver).copied().unwrap_or(i_ver);
        if added_sd.insert(n) {
            if let Some(s) = f.ds().shape(n).cloned() {
                a_ls.push(s);
            }
        }
    }
    crate::pave_ff_se::collect_micro_section_vertices(
        f,
        micro_pb,
        dm_new_sd,
        &mut added_sd,
        &mut a_ls,
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
        .map(|(i, (s, _))| (GeometryRegistry::shape_key(s), i))
        .collect();

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
        let Some(n_sx) = nested.ds().index(&sx) else { continue };
        match sx.shape_type() {
            ShapeType::Vertex => {
                let n_import = if let Some(sd) = nested.ds().has_shape_sd(n_sx) {
                    import_index(f.ds_mut(), nested.ds(), sd)?
                } else {
                    import_index(f.ds_mut(), nested.ds(), n_sx)?
                };
                let key = GeometryRegistry::shape_key(&sx);
                if let Some(&ci) = couples.get(&key) {
                    let i_x = mscpb[ci].1.index_interf;
                    let i_p = mscpb[ci].1.index;
                    if let Some(np) = f
                        .ds_mut()
                        .interf_ff_mut()
                        .get_mut(i_x)
                        .and_then(|ff| ff.change_points().get_mut(i_p))
                    {
                        np.set_index(n_import);
                    }
                } else if let Some(n_main) = f.ds().index(&sx) {
                    if n_main != n_import {
                        dm_new_sd.insert(n_main, n_import);
                        f.ds_mut().add_shape_sd(n_main, n_import);
                    }
                }
            }
            ShapeType::Edge => {
                let key = GeometryRegistry::shape_key(&sx);
                let Some(&ci) = couples.get(&key) else { continue };
                let cpb = &mscpb[ci].1;
                let Some(pb1) = cpb.pb.as_ref() else { continue };
                let b_old = pb_has_edge(pb1);
                let i_x = cpb.index_interf;
                let i_c = cpb.index;
                let lpbx = nested.ds().pave_blocks(n_sx).to_vec();
                if lpbx.is_empty() {
                    if !b_old {
                        register_one(f, &sx, cpb)?;
                    }
                    continue;
                }
                let micro = lpbx.len() == 1 && !lpbx[0].has_shrunk_data();
                if micro {
                    if !b_old {
                        remove_curve_pb(f, i_x, i_c, pb1);
                    }
                    for child in ShapeIterator::of_shape(&sx) {
                        if child.shape_type() == ShapeType::Vertex {
                            a_ls.push(child);
                        }
                    }
                    continue;
                }
                if !b_old {
                    remove_curve_pb(f, i_x, i_c, pb1);
                }
                for pbx in &lpbx {
                    let pbr = nested.ds().real_pave_block(pbx);
                    let (n0, n1) = pbx.indices();
                    let i0 = import_index(f.ds_mut(), nested.ds(), n0)?;
                    let i1 = import_index(f.ds_mut(), nested.ds(), n1)?;
                    let i_e = import_index(f.ds_mut(), nested.ds(), pbr.edge())?;
                    let mut pb_new = BopdsPaveBlock::new();
                    pb_new.set_pave1(BopdsPave::new(i0, pbx.pave1().param));
                    pb_new.set_pave2(BopdsPave::new(i1, pbx.pave2().param));
                    pb_new.set_edge(i_e);
                    if nested.ds().is_common_block(pbx) {
                        if let Some(cb) = nested.ds().common_block(pbx) {
                            let a_tol = crate::bopalgo_tools::compute_tolerance_of_cb(
                                cb,
                                nested.ds(),
                                nested.context(),
                            );
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
                    }
                    if b_old {
                        pb_new.set_original_edge(pb1.original_edge());
                        dm_ex
                            .entry(crate::pave_ff_exist::pb_key(pb1))
                            .or_default()
                            .push(pb_new);
                    } else if let Some(lpbc) = f
                        .ds_mut()
                        .interf_ff_mut()
                        .get_mut(i_x)
                        .and_then(|ff| ff.change_curves().get_mut(i_c))
                        .map(|nc| nc.change_pave_blocks())
                    {
                        lpbc.push(pb_new);
                    }
                }
            }
            _ => {}
        }
    }
    Ok(())
}

fn unused_section_vertices(f: &PaveFiller) -> Vec<TopoShape> {
    crate::pave_ff_unused::unused_stick_vertices(f)
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
                    if nested.shape_info(sub).map(|s| s.shape_type()) == Some(ShapeType::Vertex) {
                        si.change_sub_shapes().push(import_index(main, nested, sub)?);
                    }
                }
            }
            Ok(main.append_info(si))
        }
        _ => main.append(s),
    }
}

fn register_one(f: &mut PaveFiller, shape: &TopoShape, cpb: &CoupleOfPaveBlocks) -> Result<(), String> {
    match shape.shape_type() {
        ShapeType::Vertex => {
            let i_v = f.ds_mut().append(shape.clone())?;
            let i_x = cpb.index_interf;
            let i_p = cpb.index;
            if let Some(np) = f.ds_mut().interf_ff_mut().get_mut(i_x).and_then(|ff| ff.change_points().get_mut(i_p))
            {
                np.set_index(i_v);
            }
        }
        ShapeType::Edge => {
            let Some(pb_ref) = cpb.pb.as_ref() else { return Ok(()) };
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
        }
        _ => {}
    }
    Ok(())
}
