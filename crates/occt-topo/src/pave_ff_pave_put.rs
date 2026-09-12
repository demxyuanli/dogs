//! PutPavesOnCurve, FilterPavesOnCurves, PutPaveOnCurve.
//!
//! Source: `BOPAlgo_PaveFiller_6.cxx` (`PutPavesOnCurve` at 2372,
//! `FilterPavesOnCurves` at 2437, `PutPaveOnCurve` at 2959).
//!
//! `PutPaveOnCurve` uses `IntTools_Context::IsVertexOnLine` (see
//! `int_tools_vertex_line`) instead of a homemade projector. Filter keeps
//! `max(saved, sqrt(maxDistKept) + Confusion)` on vertices that remain on at
//! least one curve after a grazing-projection drop.

use std::collections::{HashMap, HashSet};

use occt_core::precision::{CONFUSION, RESOLUTION};

use crate::algo_tools::D_TOLERANCE;
use crate::bopds::{BopdsDS, BopdsPave};
use crate::brep_tool::BRepTool;
use crate::int_tools_vertex_line;
use crate::pave_ff_paves::extended_tolerance;
use crate::pave_filler::PaveFiller;
use crate::shape::Vertex;

/// Sinus below which a projection is treated as grazing (`FilterPavesOnCurves`).
const SIN_ANGLE_MIN: f64 = 0.5;

struct PaveBlockDist {
    curve_idx: usize,
    square_dist: f64,
    sin_angle: f64,
    tolerance: f64,
}

/// `BOPAlgo_PaveFiller::PutPavesOnCurve` (`_6.cxx:2372`).
pub fn put_paves_on_curve(
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

/// `BOPAlgo_PaveFiller::PutPaveOnCurve` (`_6.cxx:2959`).
pub fn put_pave_on_curve(
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
    let fuzzy = f.fuzzy_value();
    let mut a_tol_v = mv_tol
        .get(&n_v)
        .copied()
        .unwrap_or_else(|| BRepTool::vertex_tolerance(&Vertex(v_shape.clone())));
    let a_ic_tol = tol_r3d + fuzzy;

    let mut t_on = int_tools_vertex_line::is_vertex_on_line_tol(
        &Vertex(v_shape.clone()),
        a_tol_v,
        curve.as_ref(),
        a_ic_tol,
    );
    if t_on.is_none() && i_check_extend != 0 && !f.verts_to_avoid_extension().contains(&n_v) {
        let mut an_extra = a_tol_v;
        if extended_tolerance(f, n_v, mi, &mut an_extra, i_check_extend) {
            if let Some(t) = int_tools_vertex_line::is_vertex_on_line_tol(
                &Vertex(v_shape.clone()),
                an_extra,
                curve.as_ref(),
                a_ic_tol,
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
    let a_p_tol = int_tools_vertex_line::adaptor_resolution(curve.as_ref(), tol_r3d.max(a_tol_v));
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

/// `BOPAlgo_PaveFiller::FilterPavesOnCurves` (`_6.cxx:2437`).
pub fn filter_paves_on_curves(ds: &mut BopdsDS, i: usize, mv_tol: &mut HashMap<usize, f64>) {
    let nb_c = ds.interf_ff()[i].curves().len();
    let mut by_vert: HashMap<usize, Vec<PaveBlockDist>> = HashMap::new();
    for j in 0..nb_c {
        let nc = &ds.interf_ff()[i].curves()[j];
        let Some(curve) = nc.curve() else {
            continue;
        };
        let tol_r3d = nc.tolerance().max(nc.tangential_tolerance());
        let Some(pb) = nc.pave_blocks().first() else {
            continue;
        };
        for pave in pb.ext_paves() {
            let n_v = pave.index;
            let Some(v_shape) = ds.shape(n_v) else {
                continue;
            };
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

    let mut to_remove: Vec<(usize, usize)> = Vec::new();
    let mut keep_tol: Vec<(usize, f64)> = Vec::new();
    for (n_v, list) in &by_vert {
        let min_dist = list.iter().map(|d| d.square_dist).fold(f64::MAX, f64::min);
        let mut max_dist_kept = -1.0f64;
        let mut is_removed = false;
        for d in list {
            let check = 100.0 * (d.tolerance * d.tolerance).max(min_dist);
            if d.square_dist > check && d.sin_angle < SIN_ANGLE_MIN {
                to_remove.push((d.curve_idx, *n_v));
                is_removed = true;
            } else if d.square_dist > max_dist_kept {
                max_dist_kept = d.square_dist;
            }
        }
        if is_removed && max_dist_kept > 0.0 {
            if let Some(&orig) = mv_tol.get(n_v) {
                let real = orig.max(max_dist_kept.sqrt() + CONFUSION);
                keep_tol.push((*n_v, real));
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
    for (n_v, real) in keep_tol {
        if let Some(s) = ds.shape(n_v).cloned() {
            Vertex(s).set_tolerance(real);
        }
    }
}
