//! `SeveralWlinesProcessing` (`IntPatch_PrmPrmIntersection.cxx:65-271`).
//!
//! Splices a freshly walked line onto an already collected one when the two
//! meet at a vertex: the vertex is replaced by the old line's point and the
//! neighbouring point of the old line is dropped, so the parameter stays
//! monotone along the merged line.

use occt_core::precision::PCONFUSION;
use occt_geom::Surface;

use crate::geom_int::GeomIntLine;
use crate::int_tools_wline::{u_resolution, v_resolution, PatchPoint, WLine, WLineWay};
use crate::intpatch_line::IntPatchLine;
use crate::intsurf_transition::IntSurfTypeTrans;

use super::periodic::make_new_point;

/// `SeveralWlinesProcessing(Surf1, Surf2, SLin, PeriodsArr, Trans1, Trans2,
/// Tol, MaxStepS1, MaxStepS2, WLline)` (`cxx:65-271`).
///
/// `wline` is both read and replaced, as in OCCT's
/// `occ::handle<IntPatch_WLine>&` argument.
pub(crate) fn several_wlines_processing(
    surf1: &dyn Surface,
    surf2: &dyn Surface,
    slin: &[GeomIntLine],
    periods: &[f64; 4],
    trans1: IntSurfTypeTrans,
    trans2: IntSurfTypeTrans,
    tol: f64,
    max_step_s1: f64,
    max_step_s2: f64,
    wline: &mut WLine,
) {
    if slin.is_empty() {
        return;
    }

    // `cnbV` is captured once, before any line is replaced (`cxx:87`).
    let c_nbv = wline.nb_vertex();
    for ci_v in 1..=c_nbv {
        let mut pnt_d_min = 1.0e100;
        let mut vd_min = 0i32;
        let mut wl_d_min = 0i32;

        let c_pv = wline.vertex(ci_v).p;
        let (a_u1, a_v1, a_u2, a_v2) = wline.vertex(ci_v).parameters();
        let a_pcs1 = (a_u1, a_v1);
        let a_pcs2 = (a_u2, a_v2);

        for (i_l, s_line) in slin.iter().enumerate() {
            // `aType != IntPatch_Walking` -> `continue` (`cxx:103-107`): the
            // downcast below is the same test on this port's line enum.
            let Some(a_wline) = s_line.as_wline() else {
                continue;
            };
            let t_nbv = a_wline.nb_vertex();
            // `aU1` / `aV1` / ... keep being overwritten by the inner loop
            // (`cxx:112`); the copies `aPCS1` / `aPCS2` taken before it are the
            // ones compared below.
            for ti_v in 1..=t_nbv {
                let t_pv = a_wline.vertex(ti_v).p;
                let (b_u1, b_v1, b_u2, b_v2) = a_wline.vertex(ti_v).parameters();
                let a_pts1 = (b_u1, b_v1);
                let a_pts2 = (b_u2, b_v2);

                let t_distance = c_pv.distance(&t_pv);
                let u_rs1 = u_resolution(surf1, t_distance);
                let v_rs1 = v_resolution(surf1, t_distance);
                let u_rs2 = u_resolution(surf2, t_distance);
                let v_rs2 = v_resolution(surf2, t_distance);
                let r_max_s1 = u_rs1.max(v_rs1);
                let r_max_s2 = u_rs2.max(v_rs2);

                if r_max_s1 < max_step_s1 && r_max_s2 < max_step_s2 {
                    if pnt_d_min > t_distance && t_distance > PCONFUSION {
                        let sq_dist1 = (a_pcs1.0 - a_pts1.0).powi(2)
                            + (a_pcs1.1 - a_pts1.1).powi(2);
                        let sq_dist2 = (a_pcs2.0 - a_pts2.0).powi(2)
                            + (a_pcs2.1 - a_pts2.1).powi(2);
                        if sq_dist1 < r_max_s1 * r_max_s1 && sq_dist2 < r_max_s2 * r_max_s2 {
                            pnt_d_min = t_distance;
                            vd_min = ti_v;
                            wl_d_min = (i_l + 1) as i32;
                        }
                    }
                }
            }
        }

        if vd_min == 0 {
            continue;
        }

        let a_wline = match slin[(wl_d_min - 1) as usize].as_wline() {
            Some(w) => w,
            None => continue,
        };

        let ti_vpar = a_wline.vertex(vd_min).parameter_on_line() as i32;
        let ci_vpar = wline.vertex(ci_v).parameter_on_line() as i32;
        // OCCT reads the two points here (`cxx:150-153`) but only reuses the
        // `(u21, v21, u22, v22)` pair, for the replacement vertex below.
        let (u21, v21, u22, v22) = a_wline
            .point(ti_vpar)
            .parameters();

        let replace_pnt = *a_wline.point(ti_vpar);
        let c_nb_p = wline.nb_pnts();

        let mut vp_old: Vec<i32> = Vec::with_capacity(c_nbv as usize);
        for i_po in 1..=c_nbv {
            let po = wline.vertex(i_po).parameter_on_line();
            vp_old.push(po as i32);
        }

        let mut remove_next = false;
        let mut remove_prev = false;
        if ci_v == 1 {
            let d_par = (vp_old[(ci_v - 1) as usize] - vp_old[ci_v as usize]).abs();
            if d_par > 10 {
                remove_next = true;
                for i_po in (ci_v + 1)..=c_nbv {
                    let i = (i_po - 1) as usize;
                    vp_old[i] -= 1;
                }
            }
        } else if ci_v == c_nbv {
            let d_par = (vp_old[(ci_v - 1) as usize] - vp_old[(ci_v - 2) as usize]).abs();
            if d_par > 10 {
                remove_prev = true;
                let i = (ci_v - 1) as usize;
                vp_old[i] -= 1;
            }
        } else {
            let d_par_mi = (vp_old[(ci_v - 1) as usize] - vp_old[(ci_v - 2) as usize]).abs();
            let d_par_ma = (vp_old[(ci_v - 1) as usize] - vp_old[ci_v as usize]).abs();
            if d_par_mi > 10 {
                remove_prev = true;
                let i = (ci_v - 1) as usize;
                vp_old[i] -= 1;
            }
            if d_par_ma > 10 {
                remove_next = true;
                for i_po in (ci_v + 1)..=c_nbv {
                    let i = (i_po - 1) as usize;
                    if d_par_mi > 10 {
                        vp_old[i] -= 2;
                    } else {
                        vp_old[i] -= 1;
                    }
                }
            } else if d_par_mi > 10 {
                for i_po in (ci_v + 1)..=c_nbv {
                    let i = (i_po - 1) as usize;
                    vp_old[i] -= 1;
                }
            }
        }

        let p_i = ci_vpar;

        let mut new_line = WLine::new();
        new_line.base = IntPatchLine::with_transitions(false, trans1, trans2);
        new_line.set_creating_way(WLineWay::PrmPrm);
        for i_p in 1..=c_nb_p {
            if p_i == i_p {
                let new_pnt = make_new_point(&replace_pnt, wline.point(i_p), periods);
                new_line.add(new_pnt);
            } else if remove_next && i_p == p_i + 1 {
                continue;
            } else if remove_prev && i_p == p_i - 1 {
                continue;
            } else {
                new_line.add(*wline.point(i_p));
            }
        }

        let pnt3d_v = a_wline.vertex(vd_min).p;
        let mut new_vtx = PatchPoint::new(pnt3d_v, 0.0, 0.0, 0.0, 0.0, 0.0);
        new_vtx.set_value(pnt3d_v, tol, false);
        new_vtx.set_parameters(u21, v21, u22, v22);
        new_vtx.set_parameter(vp_old[(ci_v - 1) as usize] as f64);

        for i_v in 1..=c_nbv {
            if i_v == ci_v {
                new_line.vertices.push(new_vtx.clone());
            } else {
                let mut the_vtx = wline.vertex(i_v).clone();
                the_vtx.set_parameter(vp_old[(i_v - 1) as usize] as f64);
                new_line.vertices.push(the_vtx);
            }
        }

        *wline = new_line;
    }
}
