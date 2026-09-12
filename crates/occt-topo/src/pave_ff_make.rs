//! `BOPAlgo_PaveFiller::MakeBlocks` — OCCT 8.0.0 translation of the F/F loop.
//!
//! Source: `BOPAlgo_PaveFiller_6.cxx:649-1137`.
//!
//! Gaps filled relative to the previous `pave_ff::make_blocks_ff`:
//! - Glue Off early-out (unchanged).
//! - FF recheck: a pair whose curves produced no elementary pave blocks is
//!   queued and processed again after the first pass (`aFFToRecheck`).
//! - `BRepLib::FindValidRange` micro-section pave blocks go to `aMicroPB`
//!   (and their vertices to `aMVI`) instead of being dropped silently, unless
//!   either bound is a technological bound vertex.
//! - `UpdateSavedTolerance` when an existing shared edge absorbs a block.
//! - Unused vertices in `aMVTol` revert to the saved tolerance and drop their
//!   SD group (`aDMVLV.UnBind`).
//! - Stick paves receive the two faces (crease test).
//! - Closing paves use `FindValidRange`.
//!
//! PostTreat / UpdateFaceInfo / PutSEInOtherFaces stay on the existing
//! `pave_ff` / `pave_ff_update` / `pave_ff_se` implementations.

use std::collections::{HashMap, HashSet};

use occt_core::precision::{CONFUSION, PCONFUSION};

use crate::algo_tools::AlgoTools;
use crate::bopds::{BopdsDS, BopdsPaveBlock};
use crate::brep_tool::BRepTool;
use crate::int_tools_full::IntToolsContext;
use crate::inttools_range;
use crate::pave_ff::{self, CoupleOfPaveBlocks};
use crate::pave_ff_exist;
use crate::pave_ff_pave_bound;
use crate::pave_ff_pave_put;
use crate::pave_ff_pave_stick;
use crate::pave_ff_se;
use crate::pave_ff_update;
use crate::pave_filler::{GlueEnum, PaveFiller};
use crate::pave_intersect::make_sd_vertices;
use crate::shape::{Edge, Face, TopoShape, Vertex};

/// `UpdateSavedTolerance` (`_6.cxx:629`).
///
/// If a vertex of `n_e` is recorded in `mv_tol` with a smaller value than
/// the new edge tolerance, raise the saved value so the later revert does
/// not shrink it below the edge.
pub fn update_saved_tolerance(
    ds: &BopdsDS,
    n_e: usize,
    tol_new: f64,
    mv_tol: &mut HashMap<usize, f64>,
) {
    let Some(si) = ds.shape_info(n_e) else {
        return;
    };
    for &n_v in si.sub_shapes() {
        if let Some(saved) = mv_tol.get_mut(&n_v) {
            if *saved < tol_new {
                *saved = tol_new;
            }
        }
    }
}

/// Revert unused vertices to the tolerance saved in `mv_tol` (`_6.cxx:1074`).
fn revert_unused_vertices(
    f: &mut PaveFiller,
    mv_tol: &HashMap<usize, f64>,
    dmvlv: &mut HashMap<usize, Vec<usize>>,
) {
    for (&n_v, &a_tol) in mv_tol {
        let Some(s) = f.ds().shape(n_v).cloned() else {
            continue;
        };
        Vertex(s.clone()).set_tolerance(a_tol);
        f.ds_mut().refresh_vertex_box(n_v, a_tol + CONFUSION);
        dmvlv.remove(&n_v);
    }
}

#[allow(dead_code)]
fn _is_existing_vertex_local(ds: &BopdsDS, p: &occt_core::gp::GpPnt, tol: f64, mv_on_in: &HashSet<usize>) -> bool {
    let mut box_p = occt_core::bnd::BndBox::new();
    box_p.add_point(p);
    box_p.enlarge(tol);
    for &n_v in mv_on_in {
        if let Some(bv) = ds.box_of(n_v) {
            if box_p.is_out_box(bv) {
                continue;
            }
        }
        let Some(v) = ds.shape(n_v) else {
            continue;
        };
        if AlgoTools::compute_vv(v, p, tol) != 0 {
            return true;
        }
    }
    false
}

fn is_valid_block_for_faces(
    ctx: &mut IntToolsContext,
    curve: &dyn occt_geom::Curve,
    t1: f64,
    t2: f64,
    f1: &Face,
    f2: &Face,
    tol: f64,
    pc1: Option<&dyn occt_geom2d::curve::Curve2d>,
    pc2: Option<&dyn occt_geom2d::curve::Curve2d>,
) -> bool {
    // `IntTools_Context::IsValidBlockForFaces` (`IntTools_Context.cxx:717`):
    // classify only the IntermediatePoint. With a pcurve, sample it at the
    // 3D parameter and run `IsPointInOnFace` (2D). Without a pcurve, run
    // `IsValidPointForFace` (ProjPS + distance + 2D InOn).
    let t = crate::boptools_2d::intermediate_point(t1, t2);
    let p = curve.d0(t);
    valid_block_point_for_face(ctx, f1, &p, t, pc1, tol)
        && valid_block_point_for_face(ctx, f2, &p, t, pc2, tol)
}

/// One face of `IntTools_Context::IsValidBlockForFaces`.
///
/// FaceFace pcurves from `make_pcurve_full` may be UV-speed Geom2dLines that
/// do not share the 3D parameter. Use the 2D branch only when the pcurve at
/// `t` reconstructs to the 3D midpoint within `tol`; otherwise fall through
/// to 3D `IsValidPointForFace`.
fn valid_block_point_for_face(
    ctx: &mut IntToolsContext,
    face: &Face,
    p: &occt_core::gp::GpPnt,
    t: f64,
    pc: Option<&dyn occt_geom2d::curve::Curve2d>,
    tol: f64,
) -> bool {
    if let Some(pc) = pc {
        let uv = pc.d0(t);
        let reconstructed = BRepTool::face_surface(face)
            .map(|s| s.d0(uv.x(), uv.y()).distance(p) <= tol.max(CONFUSION))
            .unwrap_or(false);
        if reconstructed {
            return ctx
                .is_valid_point_for_face((uv.x(), uv.y()), face)
                .unwrap_or(false);
        }
    }
    ctx.is_point_in_on_face(face, p, None, tol).unwrap_or(false)
}

fn make_pcurve_on_faces(
    edge: &TopoShape,
    f1: &Face,
    f2: &Face,
    pc1: Option<std::sync::Arc<dyn occt_geom2d::curve::Curve2d>>,
    pc2: Option<std::sync::Arc<dyn occt_geom2d::curve::Curve2d>>,
    on_s1: bool,
    on_s2: bool,
) {
    let e = Edge(edge.clone());
    let reg = crate::tgeometry::GeometryRegistry::global();
    if on_s1 {
        if let Some(pc) = pc1 {
            let key = crate::tgeometry::GeometryRegistry::shape_key(&f1.0);
            reg.set_edge_pcurve(edge, key, pc);
        } else if let Ok(pc) = AlgoTools::make_pcurve(&e, f1) {
            let key = crate::tgeometry::GeometryRegistry::shape_key(&f1.0);
            reg.set_edge_pcurve(edge, key, pc);
        }
    }
    if on_s2 {
        if let Some(pc) = pc2 {
            let key = crate::tgeometry::GeometryRegistry::shape_key(&f2.0);
            reg.set_edge_pcurve(edge, key, pc);
        } else if let Ok(pc) = AlgoTools::make_pcurve(&e, f2) {
            let key = crate::tgeometry::GeometryRegistry::shape_key(&f2.0);
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

/// `BOPAlgo_PaveFiller::MakeBlocks` (`_6.cxx:649`).
pub fn make_blocks_ff(f: &mut PaveFiller) -> Result<(), String> {
    if f.glue() != GlueEnum::None {
        return Ok(());
    }
    let nb_ff_orig = f.ds().interf_ff().len();
    if nb_ff_orig == 0 {
        return Ok(());
    }

    let mut mscpb: Vec<(TopoShape, CoupleOfPaveBlocks)> = Vec::new();
    let mut dmvlv: HashMap<usize, Vec<usize>> = HashMap::new();
    let mut dm_ex: pave_ff_update::DmExEdges = HashMap::new();
    let mut pb_faces: pave_ff_update::PbFacesMap = HashMap::new();
    let mut mpb_add: HashSet<pave_ff_exist::PbKey> = HashSet::new();
    let mut micro_pb: Vec<BopdsPaveBlock> = Vec::new();
    let mut verts_on_rejected: Vec<TopoShape> = Vec::new();
    let mut ctx_tools = IntToolsContext::new();
    let fuzzy = f.fuzzy_value();
    let distances = f.distances().clone();
    let (pc_on_s1, pc_on_s2) = {
        use crate::pave_blocks::PaveFillerLike;
        f.section_pcurve_on()
    };

    let mut ff_to_recheck: Vec<usize> = Vec::new();
    let mut nb_ff = nb_ff_orig;
    let mut i = 0usize;
    while i < nb_ff {
        let a_cur = if i < nb_ff_orig {
            i
        } else {
            ff_to_recheck[i - nb_ff_orig]
        };
        i += 1;

        let mut mv_tol: HashMap<usize, f64> = HashMap::new();
        let (n_f1, n_f2, nb_c, nb_p) = {
            let ff = &f.ds().interf_ff()[a_cur];
            let (a, b) = ff.indices();
            (a, b, ff.curves().len(), ff.points().len())
        };
        if nb_c == 0 && nb_p == 0 {
            continue;
        }
        let Some(f1_shape) = f.ds().shape(n_f1).cloned() else {
            continue;
        };
        let Some(f2_shape) = f.ds().shape(n_f2).cloned() else {
            continue;
        };
        let face1 = Face(f1_shape);
        let face2 = Face(f2_shape);
        let tol_ff = BRepTool::face_tolerance(&face1).max(BRepTool::face_tolerance(&face2));

        f.ds_mut().ensure_face_info(n_f1);
        f.ds_mut().ensure_face_info(n_f2);

        let (mv_on_in, mv_common, pbs_on_in) = f.ds().sub_shapes_on_in(n_f1, n_f2);
        let lse = f.ds().shared_edges(n_f1, n_f2);
        let (mv_stick, mv_ef, mi) = crate::pave_ff_misc::get_stick_vertices(f.ds(), n_f1, n_f2);
        let (on1, in1) = match f.ds().face_info(n_f1) {
            Some(fi) => (fi.paves_on().to_vec(), fi.paves_in().to_vec()),
            None => (Vec::new(), Vec::new()),
        };
        let (on2, in2) = match f.ds().face_info(n_f2) {
            Some(fi) => (fi.paves_on().to_vec(), fi.paves_in().to_vec()),
            None => (Vec::new(), Vec::new()),
        };
        let on_in_pbs = pave_ff_exist::resolve_on_in_pbs(f.ds(), &pbs_on_in);
        let pb_common = pave_ff_exist::common_on_in(&on1, &in1, &on2, &in2);

        for j in 0..nb_p {
            let p = f.ds().interf_ff()[a_cur].points()[j].pnt().clone();
            if crate::pave_ff_is_exist::is_existing_vertex(
                f.ds(),
                &p,
                tol_ff,
                &mv_on_in,
                fuzzy,
            ) {
                continue;
            }
            let v = AlgoTools::make_new_vertex(&p, tol_ff)?;
            mscpb.push((
                v,
                CoupleOfPaveBlocks {
                    index_interf: a_cur,
                    index: j,
                    pb: None,
                },
            ));
        }

        for j in 0..nb_c {
            f.ds_mut().interf_ff_mut()[a_cur].change_curves()[j].init_pave_block1();
            pave_ff_pave_put::put_paves_on_curve(
                f,
                a_cur,
                j,
                &mv_on_in,
                &mv_common,
                &mi,
                &mv_ef,
                &mut mv_tol,
                &mut dmvlv,
            );
        }
        pave_ff_pave_put::filter_paves_on_curves(f.ds_mut(), a_cur, &mut mv_tol);

        let mut dmbv: HashMap<usize, Vec<usize>> = HashMap::new();
        let mut mv_bounds: HashSet<usize> = HashSet::new();
        for j in 0..nb_c {
            pave_ff_pave_stick::put_stick_paves_on_curve(
                f,
                a_cur,
                j,
                &face1,
                &face2,
                &mi,
                &mv_stick,
                &mut mv_tol,
                &mut dmvlv,
            );
            if nb_c == 1 {
                pave_ff_pave_stick::put_ef_paves_on_curve(
                    f, a_cur, j, &mi, &mv_ef, &mut mv_tol, &mut dmvlv,
                );
            }
            let has_bounds = {
                let nc = &f.ds().interf_ff()[a_cur].curves()[j];
                let (a, b) = nc.range();
                a.is_finite() && b.is_finite()
            };
            if has_bounds {
                let lbv = pave_ff_pave_bound::put_bound_pave_on_curve(
                    f, a_cur, j, &face1, &face2, &mut ctx_tools,
                )?;
                if !lbv.is_empty() {
                    for &n in &lbv {
                        mv_bounds.insert(n);
                    }
                    dmbv.insert(j, lbv);
                }
            }
        }
        for j in 0..nb_c {
            pave_ff_pave_bound::put_closing_pave_on_curve(f, a_cur, j);
        }

        let mut is_to_recheck = nb_c > 0 && i <= nb_ff_orig;
        for j in 0..nb_c {
            let mut elementary: Vec<BopdsPaveBlock> = Vec::new();
            {
                let pb1 = f.ds_mut().interf_ff_mut()[a_cur]
                    .change_curves()[j]
                    .change_pave_block1();
                pb1.update(&mut elementary, false);
            }
            if !elementary.is_empty() {
                is_to_recheck = false;
            }
            let (curve, tol_r3d, pc1, pc2) = {
                let nc = &f.ds().interf_ff()[a_cur].curves()[j];
                let Some(c) = nc.curve().cloned() else {
                    continue;
                };
                (
                    c,
                    nc.tolerance().max(nc.tangential_tolerance()),
                    nc.pcurve1().cloned(),
                    nc.pcurve2().cloned(),
                )
            };
            let mut kept: Vec<BopdsPaveBlock> = Vec::new();
            for pb in elementary {
                let (n_v1, n_v2) = pb.indices();
                let (t1, t2) = pb.range();
                if (t1 - t2).abs() < PCONFUSION {
                    continue;
                }
                if !is_valid_block_for_faces(
                    &mut ctx_tools,
                    curve.as_ref(),
                    t1,
                    t2,
                    &face1,
                    &face2,
                    tol_r3d,
                    pc1.as_deref(),
                    pc2.as_deref(),
                ) {
                    continue;
                }
                if let Some((n_e_out, tol_new)) = crate::pave_ff_is_exist::is_existing_pave_block_on_shared(
                    f.ds(),
                    &ctx_tools,
                    &pb,
                    curve.as_ref(),
                    &lse,
                    fuzzy,
                ) {
                    let _ = crate::pave_common::update_edge_tolerance(f, n_e_out, tol_new);
                    update_saved_tolerance(f.ds(), n_e_out, tol_new, &mut mv_tol);
                    continue;
                }
                let Some(v1s) = f.ds().shape(n_v1).cloned() else {
                    continue;
                };
                let Some(v2s) = f.ds().shape(n_v2).cloned() else {
                    continue;
                };
                let tv1 = BRepTool::vertex_tolerance(&Vertex(v1s.clone())).max(tol_r3d);
                let tv2 = BRepTool::vertex_tolerance(&Vertex(v2s.clone())).max(tol_r3d);
                if inttools_range::find_valid_range(
                    curve.as_ref(),
                    t1,
                    t2,
                    &curve.d0(t1),
                    tv1,
                    &curve.d0(t2),
                    tv2,
                )
                .is_none()
                {
                    if !mv_bounds.contains(&n_v1) && !mv_bounds.contains(&n_v2) {
                        micro_pb.push(pb.clone());
                    }
                    continue;
                }
                if let Some((pb_out, mut tol_new)) = crate::pave_ff_exist_onin::is_existing_pave_block_on_in(
                    f.ds(),
                    &ctx_tools,
                    &pb,
                    curve.as_ref(),
                    tol_r3d,
                    &on_in_pbs,
                    &pb_common,
                    fuzzy,
                ) {
                    let b_in_f1 = pave_ff_exist::pb_in_face(&on1, &in1, &pb_out);
                    let b_in_f2 = pave_ff_exist::pb_in_face(&on2, &in2, &pb_out);
                    if !b_in_f1 || !b_in_f2 {
                        let n_e = pb_out.edge();
                        let nc_tol = f.ds().interf_ff()[a_cur].curves()[j].tolerance();
                        if tol_new < nc_tol {
                            tol_new = nc_tol;
                        }
                        let tol_e = f
                            .ds()
                            .shape(n_e)
                            .map(|s| BRepTool::edge_tolerance(&Edge(s.clone())))
                            .unwrap_or(0.0);
                        if tol_new > tol_e {
                            let _ = crate::pave_common::update_edge_tolerance(f, n_e, tol_new);
                            update_saved_tolerance(f.ds(), n_e, tol_new, &mut mv_tol);
                        }
                        let n_f = if b_in_f1 { n_f2 } else { n_f1 };
                        let key = pave_ff_exist::pb_key(&pb_out);
                        let faces = pb_faces.entry(key).or_default();
                        if !faces.contains(&n_f) {
                            faces.push(n_f);
                        }
                        let (n_out1, n_out2) = pb_out.indices();
                        if n_v1 != n_out1 && n_v1 != n_out2 && !mv_bounds.contains(&n_v1) {
                            if let Some(s) = f.ds().shape(n_v1).cloned() {
                                verts_on_rejected.push(s);
                            }
                        }
                        if n_v2 != n_out1 && n_v2 != n_out2 && !mv_bounds.contains(&n_v2) {
                            if let Some(s) = f.ds().shape(n_v2).cloned() {
                                verts_on_rejected.push(s);
                            }
                        }
                        if mpb_add.insert(key) {
                            pave_ff_exist::prepare_post_treat_ff(
                                f.ds(),
                                a_cur,
                                j,
                                &pb_out,
                                &mut mscpb,
                                &mut kept,
                            );
                        }
                    }
                    continue;
                }
                let edge = AlgoTools::make_edge(
                    curve.clone(),
                    Some(&v1s),
                    t1,
                    Some(&v2s),
                    t2,
                    tol_r3d,
                )?;
                make_pcurve_on_faces(
                    &edge,
                    &face1,
                    &face2,
                    pc1.clone(),
                    pc2.clone(),
                    pc_on_s1,
                    pc_on_s2,
                );
                let cpb = CoupleOfPaveBlocks {
                    index_interf: a_cur,
                    index: j,
                    pb: Some(pb.clone()),
                };
                mscpb.push((edge.clone(), cpb));
                kept.push(pb);
                mv_tol.remove(&n_v1);
                mv_tol.remove(&n_v2);
                crate::pave_ff_exist_es::process_existing_pave_blocks_es(
                    f,
                    a_cur,
                    j,
                    n_f1,
                    n_f2,
                    &edge,
                    &on_in_pbs,
                    &on1,
                    &in1,
                    &on2,
                    &in2,
                    &distances,
                    &mut mscpb,
                    &mut kept,
                    &mut mpb_add,
                    &mut pb_faces,
                );
            }
            {
                let lpbc = f.ds_mut().interf_ff_mut()[a_cur]
                    .change_curves()[j]
                    .change_pave_blocks();
                lpbc.clear();
                lpbc.extend(kept);
            }
        }
        if is_to_recheck {
            ff_to_recheck.push(a_cur);
            nb_ff += 1;
        }
        revert_unused_vertices(f, &mv_tol, &mut dmvlv);
        pave_ff_exist::process_existing_pave_blocks_bound(
            f,
            a_cur,
            n_f1,
            n_f2,
            &on_in_pbs,
            &dmbv,
            &mut mscpb,
            &mut mpb_add,
            &mut pb_faces,
        );
    }

    pave_ff_se::remove_micro_section_edges(f, &mut mscpb, &mut micro_pb);
    let mut dm_new_sd: HashMap<usize, usize> = HashMap::new();
    make_sd_vertices_ff(f, &dmvlv, &mut dm_new_sd)?;
    pave_ff::post_treat_ff(
        f,
        &mscpb,
        &mut dm_new_sd,
        &mut dm_ex,
        &micro_pb,
        &verts_on_rejected,
    )?;
    pave_ff_exist::correct_tolerance_of_se(f);
    pave_ff_update::update_face_info_ff(f, &mut dm_ex, &dm_new_sd, &pb_faces);
    pave_ff_update::update_pave_blocks_ff(f, &dm_new_sd)?;
    pave_ff_se::put_se_in_other_faces(f)?;
    Ok(())
}
