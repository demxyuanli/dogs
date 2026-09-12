//! PutPavesOnCurve, ExtendedTolerance, GetEFPnts, FilterPavesOnCurves,
//! PutBoundPaveOnCurve, PutClosingPaveOnCurve.
//!
//! Source: `BOPAlgo_PaveFiller_6.cxx` (`PutPavesOnCurve` at 2372,
//! `ExtendedTolerance` at 2542, `GetEFPnts` at 2608,
//! `PutEFPavesOnCurve` at 2692, `PutStickPavesOnCurve` at 2748,
//! `FilterPavesOnCurves` at 2437, `PutPaveOnCurve` at 2959).

use std::collections::{HashMap, HashSet};

use occt_core::gp::GpPnt;
use occt_core::precision::{CONFUSION, PCONFUSION, RESOLUTION};
use occt_geom::Curve;

use crate::algo_tools::{AlgoTools, D_TOLERANCE};
use crate::bopds::{BopdsDS, BopdsPave};
use crate::bopds_ff::BopdsCurve;
use crate::brep_tool::BRepTool;
use crate::int_tools_full::IntToolsContext;
use crate::inttools_roots;
use crate::pave_filler::PaveFiller;
use crate::shape::{Edge, Face, Vertex};

/// 2-D point on both faces of an F/F pair (`IntSurf_PntOn2S`).
#[derive(Debug, Clone, Copy)]
pub(crate) struct PntOn2s {
    pub u1: f64,
    pub v1: f64,
    pub u2: f64,
    pub v2: f64,
}

impl PntOn2s {
    pub(crate) fn set_value(&mut self, u_on_f1: f64, v_on_f1: f64, u_on_f2: f64, v_on_f2: f64) {
        self.u1 = u_on_f1;
        self.v1 = v_on_f1;
        self.u2 = u_on_f2;
        self.v2 = v_on_f2;
    }
}

fn curve_is_bezier_or_bspline(c: &dyn Curve) -> bool {
    let name = std::any::type_name_of_val(c);
    name.contains("BSpline") || name.contains("Bezier")
}

fn full_shape_map(ds: &BopdsDS, n_f: usize, mi: &mut HashSet<usize>) {
    mi.insert(n_f);
    if let Some(si) = ds.shape_info(n_f) {
        for &s in si.sub_shapes() {
            mi.insert(s);
        }
    }
}

/// `BOPAlgo_PaveFiller::ExtendedTolerance` (`_6.cxx:2542`).
pub(crate) fn extended_tolerance(
    f: &PaveFiller,
    n_v: usize,
    mi: &HashSet<usize>,
    a_tol_v_ext: &mut f64,
    a_type: i32,
) -> bool {
    if !f.ds().is_new_shape(n_v) {
        return false;
    }
    let mut k = 0usize;
    let mut a_nb_int = 2usize;
    if a_type == 1 {
        a_nb_int = 1;
    } else if a_type == 2 {
        k = 1;
    }
    let Some(v_shape) = f.ds().shape(n_v).cloned() else {
        return false;
    };
    let a_pv = BRepTool::vertex_point(&Vertex(v_shape));
    while k < a_nb_int {
        let lines: Vec<(usize, usize, Option<usize>, f64, f64)> = if k == 0 {
            f.ds()
                .interf_ee()
                .iter()
                .map(|it| {
                    (
                        it.index1(),
                        it.index2(),
                        it.get_index_new(),
                        it.common_first,
                        it.common_last,
                    )
                })
                .collect()
        } else {
            f.ds()
                .interf_ef()
                .iter()
                .map(|it| {
                    (
                        it.index1(),
                        it.index2(),
                        it.get_index_new(),
                        it.common_first,
                        it.common_last,
                    )
                })
                .collect()
        };
        for (i1, i2, index_new, t11, t12) in lines {
            if index_new != Some(n_v) {
                continue;
            }
            if !mi.contains(&i1) || !mi.contains(&i2) {
                continue;
            }
            let Some(e1_shape) = f.ds().shape(i1).cloned() else {
                continue;
            };
            let Ok(a_p11) = AlgoTools::point_on_edge(&Edge(e1_shape.clone()), t11) else {
                continue;
            };
            let Ok(a_p12) = AlgoTools::point_on_edge(&Edge(e1_shape), t12) else {
                continue;
            };
            let a_d1 = a_pv.distance(&a_p11);
            let a_d2 = a_pv.distance(&a_p12);
            let a_d = a_d1.max(a_d2);
            if a_d > *a_tol_v_ext {
                *a_tol_v_ext = a_d;
            }
            return true;
        }
        k += 1;
    }
    false
}

/// `BOPAlgo_PaveFiller::GetEFPnts` (`_6.cxx:2608`).
pub(crate) fn get_ef_pnts(f: &PaveFiller, n_f1: usize, n_f2: usize) -> Vec<PntOn2s> {
    crate::pave_ff_ef_pnts::get_ef_pnts(f, n_f1, n_f2)
}

#[allow(dead_code)]
fn get_ef_pnts_legacy(f: &PaveFiller, n_f1: usize, n_f2: usize) -> Vec<PntOn2s> {
    let mut mi = HashSet::new();
    full_shape_map(f.ds(), n_f1, &mut mi);
    full_shape_map(f.ds(), n_f2, &mut mi);
    let mut out = Vec::new();
    let efs: Vec<(usize, usize, Option<usize>, f64)> = f
        .ds()
        .interf_ef()
        .iter()
        .map(|it| {
            (
                it.index1(),
                it.index2(),
                it.get_index_new(),
                0.5 * (it.common_first + it.common_last),
            )
        })
        .collect();
    for (n_e, n_f_opp, index_new, a_par) in efs {
        if index_new.is_none() {
            continue;
        }
        if !mi.contains(&n_e) || !mi.contains(&n_f_opp) {
            continue;
        }
        let Some(e_shape) = f.ds().shape(n_e).cloned() else {
            continue;
        };
        let Some(f_opp_shape) = f.ds().shape(n_f_opp).cloned() else {
            continue;
        };
        let a_e = Edge(e_shape);
        let a_f_opp = Face(f_opp_shape);
        let Some(a_curve) = BRepTool::edge_curve(&a_e) else {
            continue;
        };
        let n_f = if n_f_opp == n_f1 { n_f2 } else { n_f1 };
        let Some(f_shape) = f.ds().shape(n_f).cloned() else {
            continue;
        };
        let a_f = Face(f_shape);
        let a_point = a_curve.d0(a_par);
        let ctx = f.context();
        let pc = crate::boptools_2d::curve_on_surface(&a_e, &a_f);
        let mut a_pnt = PntOn2s {
            u1: 0.0,
            v1: 0.0,
            u2: 0.0,
            v2: 0.0,
        };
        if let Some(a_pcurve) = pc {
            let a_p2d = a_pcurve.d0(a_par);
            if let Ok((u1, v1)) = ctx.project_point_on_face(&a_f_opp, &a_point) {
                if n_f == n_f1 {
                    a_pnt.set_value(a_p2d.x(), a_p2d.y(), u1, v1);
                } else {
                    a_pnt.set_value(u1, v1, a_p2d.x(), a_p2d.y());
                }
                out.push(a_pnt);
            }
        } else if let (Ok((u1, v1)), Ok((u2, v2))) = (
            ctx.project_point_on_face(&a_f, &a_point),
            ctx.project_point_on_face(&a_f_opp, &a_point),
        ) {
            if n_f == n_f1 {
                a_pnt.set_value(u1, v1, u2, v2);
            } else {
                a_pnt.set_value(u2, v2, u1, v1);
            }
            out.push(a_pnt);
        }
    }
    out
}

/// `BOPAlgo_PaveFiller::PutPavesOnCurve` (`_6.cxx:2372`).
pub(crate) fn put_paves_on_curve(
    f: &mut PaveFiller,
    i: usize,
    j: usize,
    mv_on_in: &HashSet<usize>,
    mv_common: &HashSet<usize>,
    mi: &HashSet<usize>,
    mv_ef: &HashSet<usize>,
    mv_tol: &mut HashMap<usize, f64>,
    dmvlv: &mut HashMap<usize, Vec<usize>>,
) {
    let (tol_r3d, box_c) = {
        let nc = &f.ds().interf_ff()[i].curves()[j];
        (
            nc.tolerance().max(nc.tangential_tolerance()),
            nc.bounding_box().clone(),
        )
    };
    for &n_v in mv_ef {
        put_pave_on_curve(f, i, j, n_v, tol_r3d, mi, mv_tol, dmvlv, 2);
    }
    let others: Vec<usize> = mv_on_in
        .iter()
        .copied()
        .filter(|n| !mv_ef.contains(n))
        .collect();
    for n_v in others {
        if !mv_common.contains(&n_v) {
            if let Some(bv) = f.ds().box_of(n_v) {
                if box_c.is_out_box(bv) {
                    continue;
                }
            }
            if !f.ds().is_new_shape(n_v) {
                continue;
            }
        }
        put_pave_on_curve(f, i, j, n_v, tol_r3d, mi, mv_tol, dmvlv, 1);
    }
}

/// `BOPAlgo_PaveFiller::PutStickPavesOnCurve` (`_6.cxx:2748`).
pub(crate) fn put_stick_paves_on_curve(
    f: &mut PaveFiller,
    i: usize,
    j: usize,
    mi: &HashSet<usize>,
    mv_stick: &HashSet<usize>,
    mv_tol: &mut HashMap<usize, f64>,
    dmvlv: &mut HashMap<usize, Vec<usize>>,
) {
    let bnd = get_bound_paves(f.ds(), &f.ds().interf_ff()[i].curves()[j]);
    if bnd[0].is_some() && bnd[1].is_some() {
        return;
    }
    let already: HashSet<usize> = f.ds().interf_ff()[i].curves()[j]
        .pave_blocks()
        .first()
        .map(|pb| pb.ext_paves().iter().map(|p| p.index).collect())
        .unwrap_or_default();
    let verts: Vec<usize> = mv_stick
        .iter()
        .copied()
        .filter(|n| !already.contains(n))
        .collect();
    if verts.is_empty() {
        return;
    }
    let tol_r3d = {
        let nc = &f.ds().interf_ff()[i].curves()[j];
        nc.tolerance().max(nc.tangential_tolerance())
    };
    for n_v in verts {
        put_pave_on_curve(f, i, j, n_v, tol_r3d, mi, mv_tol, dmvlv, 1);
    }
}

/// `BOPAlgo_PaveFiller::PutEFPavesOnCurve` (`_6.cxx:2692`).
pub(crate) fn put_ef_paves_on_curve(
    f: &mut PaveFiller,
    i: usize,
    j: usize,
    mi: &HashSet<usize>,
    mv_ef: &HashSet<usize>,
    mv_tol: &mut HashMap<usize, f64>,
    dmvlv: &mut HashMap<usize, Vec<usize>>,
) {
    if mv_ef.is_empty() {
        return;
    }
    let Some(curve) = f.ds().interf_ff()[i].curves()[j].curve().cloned() else {
        return;
    };
    if !curve_is_bezier_or_bspline(curve.as_ref()) {
        return;
    }
    let already: HashSet<usize> = f.ds().interf_ff()[i].curves()[j]
        .pave_blocks()
        .first()
        .map(|pb| pb.ext_paves().iter().map(|p| p.index).collect())
        .unwrap_or_default();
    let verts: Vec<usize> = mv_ef
        .iter()
        .copied()
        .filter(|n| !already.contains(n))
        .collect();
    let (a, b) = {
        let nc = &f.ds().interf_ff()[i].curves()[j];
        let (t0, t1) = nc.range();
        if t0.is_finite() && t1.is_finite() {
            (t0, t1)
        } else {
            (curve.first_parameter(), curve.last_parameter())
        }
    };
    for n_v in verts {
        let Some(v_shape) = f.ds().shape(n_v).cloned() else {
            continue;
        };
        let a_pv = BRepTool::vertex_point(&Vertex(v_shape));
        let t = inttools_roots::parameter(&|u| curve.d0(u), &a_pv, a, b);
        let a_dist = a_pv.distance(&curve.d0(t));
        put_pave_on_curve(f, i, j, n_v, a_dist, mi, mv_tol, dmvlv, 2);
    }
}

/// `BOPAlgo_PaveFiller::PutPaveOnCurve` (`_6.cxx:2959`).
pub(crate) fn put_pave_on_curve(
    f: &mut PaveFiller,
    i: usize,
    j: usize,
    n_v: usize,
    tol_r3d: f64,
    mi: &HashSet<usize>,
    mv_tol: &mut HashMap<usize, f64>,
    dmvlv: &mut HashMap<usize, Vec<usize>>,
    i_check_extend: i32,
) {
    let Some(v_shape) = f.ds().shape(n_v).cloned() else {
        return;
    };
    let Some(curve) = f.ds().interf_ff()[i].curves()[j].curve().cloned() else {
        return;
    };
    let range = f.ds().interf_ff()[i].curves()[j].range();
    let fuzzy = f.fuzzy_value();
    let mut a_tol_v = mv_tol
        .get(&n_v)
        .copied()
        .unwrap_or_else(|| BRepTool::vertex_tolerance(&Vertex(v_shape.clone())));
    let mut t_on = vertex_on_curve_param(
        &Vertex(v_shape.clone()),
        curve.as_ref(),
        range,
        a_tol_v,
        tol_r3d,
        fuzzy,
    );
    if t_on.is_none() && i_check_extend != 0 && !f.verts_to_avoid_extension().contains(&n_v) {
        let mut an_extra = a_tol_v;
        if extended_tolerance(f, n_v, mi, &mut an_extra, i_check_extend) {
            if let Some(t) = vertex_on_curve_param(
                &Vertex(v_shape.clone()),
                curve.as_ref(),
                range,
                an_extra,
                tol_r3d,
                fuzzy,
            ) {
                let a_p_on_c = curve.d0(t);
                a_tol_v = a_p_on_c.distance(&BRepTool::vertex_point(&Vertex(v_shape.clone())));
                t_on = Some(t);
            }
        }
    }
    let Some(t) = t_on else {
        return;
    };
    let a_d_tol = D_TOLERANCE;
    let a_p_tol = PCONFUSION.max(tol_r3d.max(a_tol_v));
    let existing = f.ds().interf_ff()[i].curves()[j]
        .pave_blocks()
        .first()
        .and_then(|pb| pb.contains_parameter(t, a_p_tol));
    if let Some(n_used) = existing {
        let list = dmvlv.entry(n_used).or_insert_with(|| vec![n_used]);
        if !list.contains(&n_v) {
            list.push(n_v);
        }
        if !mv_tol.contains_key(&n_used) {
            if let Some(vu) = f.ds().shape(n_used) {
                mv_tol.insert(n_used, BRepTool::vertex_tolerance(&Vertex(vu.clone())));
            }
        }
        if !mv_tol.contains_key(&n_v) {
            mv_tol.insert(n_v, BRepTool::vertex_tolerance(&Vertex(v_shape)));
        }
        return;
    }
    let pb = f.ds_mut().interf_ff_mut()[i].change_curves()[j].change_pave_block1();
    pb.append_ext_pave(BopdsPave::new(n_v, t));
    let p1 = curve.d0(t);
    let p2 = BRepTool::vertex_point(&Vertex(v_shape.clone()));
    let dist = p1.distance(&p2);
    let tol_v = BRepTool::vertex_tolerance(&Vertex(v_shape.clone()));
    if tol_v < dist + a_d_tol {
        if !mv_tol.contains_key(&n_v) {
            mv_tol.insert(n_v, tol_v);
        }
        Vertex(v_shape).set_tolerance(dist + a_d_tol);
        f.ds_mut().refresh_vertex_box(n_v, dist + a_d_tol);
    }
}

fn vertex_on_curve_param(
    v: &Vertex,
    curve: &dyn Curve,
    range: (f64, f64),
    tol_v: f64,
    tol_r3d: f64,
    fuzzy: f64,
) -> Option<f64> {
    let p = BRepTool::vertex_point(v);
    let (mut a, mut b) = range;
    if !a.is_finite() || !b.is_finite() || b <= a {
        a = curve.first_parameter();
        b = curve.last_parameter();
        if !a.is_finite() || !b.is_finite() || b <= a {
            return None;
        }
    }
    let t = inttools_roots::parameter(&|u| curve.d0(u), &p, a, b);
    let dist = p.distance(&curve.d0(t));
    let tol_sum = (2.0 * (tol_v + (tol_r3d + fuzzy).max(0.0))).max(1e-6);
    if dist <= tol_sum {
        Some(t)
    } else {
        None
    }
}

struct PaveBlockDist {
    curve_idx: usize,
    square_dist: f64,
    sin_angle: f64,
    tolerance: f64,
}

/// `BOPAlgo_PaveFiller::FilterPavesOnCurves` (`_6.cxx:2437`).
pub(crate) fn filter_paves_on_curves(
    ds: &mut BopdsDS,
    i: usize,
    mv_tol: &mut HashMap<usize, f64>,
) {
    let nb_c = ds.interf_ff()[i].curves().len();
    let mut by_vert: HashMap<usize, Vec<PaveBlockDist>> = HashMap::new();
    for j in 0..nb_c {
        let nc = &ds.interf_ff()[i].curves()[j];
        let Some(curve) = nc.curve() else { continue };
        let tol_r3d = nc.tolerance().max(nc.tangential_tolerance());
        let Some(pb) = nc.pave_blocks().first() else { continue };
        for pave in pb.ext_paves() {
            let n_v = pave.index;
            let Some(v_shape) = ds.shape(n_v) else { continue };
            let pv = BRepTool::vertex_point(&Vertex(v_shape.clone()));
            let par = pave.param;
            let (ponc, d1) = curve.d1(par);
            let dx = ponc.x() - pv.x();
            let dy = ponc.y() - pv.y();
            let dz = ponc.z() - pv.z();
            let sq_dist = dx * dx + dy * dy + dz * dz;
            let sq_d1 = d1.square_magnitude();
            let mut sin = 0.0;
            if sq_dist > RESOLUTION && sq_d1 > RESOLUTION {
                let cx = dy * d1.z() - dz * d1.y();
                let cy = dz * d1.x() - dx * d1.z();
                let cz = dx * d1.y() - dy * d1.x();
                let cross_sq = cx * cx + cy * cy + cz * cz;
                sin = (cross_sq / (sq_dist * sq_d1)).sqrt();
            }
            by_vert.entry(n_v).or_default().push(PaveBlockDist {
                curve_idx: j,
                square_dist: sq_dist,
                sin_angle: sin,
                tolerance: tol_r3d,
            });
        }
    }
    const SIN_ANGLE_MIN: f64 = 0.5;
    let mut to_remove: Vec<(usize, usize)> = Vec::new();
    let mut restore: Vec<(usize, f64)> = Vec::new();
    for (n_v, list) in &by_vert {
        let min_dist = list
            .iter()
            .map(|d| d.square_dist)
            .fold(f64::MAX, f64::min);
        for d in list {
            let check = 100.0 * (d.tolerance * d.tolerance).max(min_dist);
            if d.square_dist > check && d.sin_angle < SIN_ANGLE_MIN {
                to_remove.push((d.curve_idx, *n_v));
                if let Some(&orig) = mv_tol.get(n_v) {
                    restore.push((*n_v, orig));
                }
            }
        }
    }
    for (j, n_v) in to_remove {
        if let Some(pb) = ds.interf_ff_mut()[i]
            .change_curves()[j]
            .change_pave_blocks()
            .first_mut()
        {
            pb.remove_ext_pave(n_v);
        }
    }
    for (n_v, orig) in restore {
        if let Some(s) = ds.shape(n_v).cloned() {
            Vertex(s).set_tolerance(orig);
        }
    }
}

pub(crate) fn get_bound_paves(ds: &BopdsDS, nc: &BopdsCurve) -> [Option<usize>; 2] {
    let Some(pb) = nc.pave_blocks().first() else {
        return [None, None];
    };
    if pb.ext_paves().is_empty() {
        return [None, None];
    }
    let mut tmin = f64::MAX;
    let mut tmax = f64::MIN;
    let mut nv = [None, None];
    for p in pb.ext_paves() {
        if p.param < tmin {
            nv[0] = Some(p.index);
            tmin = p.param;
        }
        if p.param > tmax {
            nv[1] = Some(p.index);
            tmax = p.param;
        }
    }
    let (t0, t1) = nc.range();
    let Some(curve) = nc.curve() else { return nv };
    let pts = [curve.d0(t0), curve.d0(t1)];
    let tol = nc.tolerance() + CONFUSION;
    for j in 0..2 {
        let Some(n) = nv[j] else { continue };
        let Some(v) = ds.shape(n) else {
            nv[j] = None;
            continue;
        };
        if AlgoTools::compute_vv(v, &pts[j], tol) == 0 {
            nv[j] = None;
        }
    }
    nv
}

/// `BOPAlgo_PaveFiller::PutBoundPaveOnCurve`.
pub(crate) fn put_bound_pave_on_curve(
    f: &mut PaveFiller,
    i: usize,
    j: usize,
    face1: &Face,
    face2: &Face,
    ctx: &mut IntToolsContext,
) -> Result<Vec<usize>, String> {
    let (curve, range, tol_r3d) = {
        let nc = &f.ds().interf_ff()[i].curves()[j];
        let Some(c) = nc.curve().cloned() else {
            return Ok(Vec::new());
        };
        (c, nc.range(), nc.tolerance())
    };
    let (t0, t1) = range;
    let p0 = curve.d0(t0);
    let p1 = curve.d0(t1);
    let bnd = get_bound_paves(f.ds(), &f.ds().interf_ff()[i].curves()[j]);
    let is_closed = p1.distance(&p0) <= CONFUSION;
    if is_closed && (bnd[0].is_some() || bnd[1].is_some()) {
        return Ok(Vec::new());
    }
    let mut lbv = Vec::new();
    let ts = [t0, t1];
    let ps = [p0, p1];
    for j_end in 0..2 {
        if bnd[j_end].is_some() {
            continue;
        }
        if j_end == 1 && is_closed {
            continue;
        }
        let ok_f = ctx
            .is_point_in_on_face(face1, &ps[j_end], None, tol_r3d)
            .unwrap_or(false);
        let ok_g = ctx
            .is_point_in_on_face(face2, &ps[j_end], None, tol_r3d)
            .unwrap_or(false);
        if !ok_f || !ok_g {
            continue;
        }
        let mut vn = AlgoTools::make_new_vertex(&ps[j_end], tol_r3d)?;
        let dist = BRepTool::vertex_point(&Vertex(vn.clone())).distance(&ps[j_end]);
        let need = dist + D_TOLERANCE;
        if BRepTool::vertex_tolerance(&Vertex(vn.clone())) < need {
            Vertex(vn.clone()).set_tolerance(need);
        }
        let n_vn = f.ds_mut().append(vn)?;
        f.ds_mut().refresh_vertex_box(n_vn, need);
        let pb = f.ds_mut().interf_ff_mut()[i].change_curves()[j].change_pave_block1();
        pb.append_ext_pave(BopdsPave::new(n_vn, ts[j_end]));
        lbv.push(n_vn);
    }
    Ok(lbv)
}

/// `BOPAlgo_PaveFiller::PutClosingPaveOnCurve`.
pub(crate) fn put_closing_pave_on_curve(ds: &mut BopdsDS, i: usize, j: usize) {
    let (curve, range, tol_r3d) = {
        let nc = &ds.interf_ff()[i].curves()[j];
        let Some(c) = nc.curve().cloned() else { return };
        (c, nc.range(), nc.tolerance())
    };
    let (t0, t1) = range;
    if !t0.is_finite() || !t1.is_finite() {
        return;
    }
    let p0 = curve.d0(t0);
    let p1 = curve.d0(t1);
    let ts = [t0, t1];
    let ps = [p0, p1];
    let ext: Vec<BopdsPave> = ds.interf_ff()[i].curves()[j]
        .pave_blocks()
        .first()
        .map(|pb| pb.ext_paves().to_vec())
        .unwrap_or_default();
    let mut found: Option<(usize, f64, GpPnt)> = None;
    for pave in &ext {
        for k in 0..2 {
            if (pave.param - ts[k]).abs() < PCONFUSION {
                found = Some((pave.index, ts[1 - k], ps[1 - k].clone()));
                break;
            }
        }
        if found.is_some() {
            break;
        }
    }
    let Some((n_v, t_op, p_op)) = found else {
        return;
    };
    let Some(v_shape) = ds.shape(n_v).cloned() else {
        return;
    };
    let vtx = Vertex(v_shape);
    let tol_v = BRepTool::vertex_tolerance(&vtx);
    let pv = BRepTool::vertex_point(&vtx);
    let tol_p = tol_r3d + CONFUSION;
    let dist = pv.distance(&p_op);
    if dist > tol_v + tol_p {
        return;
    }
    let pb = ds.interf_ff_mut()[i].change_curves()[j].change_pave_block1();
    pb.append_ext_pave1(BopdsPave::new(n_v, t_op));
}
