//! `IntPatch_PrmPrmIntersection::Perform(..., NCollection_List<IntSurf_PntOn2S>&
//! LOfPnts)` (`IntPatch_PrmPrmIntersection.cxx:1827-2172`).
//!
//! Walks a line from every given start point and appends the accepted walking
//! lines to the intersection sequence.

use occt_core::precision::PCONFUSION;
use occt_geom::Surface;

use crate::geom_int::{GeomIntLine, TopolTool};
use crate::int_tools_wline::{PntOn2S, WLine, WLineWay};
use crate::intpatch_line::IntPatchLine;
use crate::intpatch_rst_int::put_vertex_on_line;
use crate::intsurf_transition::IntSurfTypeTrans;
use crate::intwalk::PWalking;

use super::several::several_wlines_processing;
use super::wline_helpers::{add_wline, duplicate_of_lines_processing, is_point_on_line};

/// `Perform(Surf1, D1, Surf2, D2, TolTangency, Epsilon, Deflection,
/// Increment, LOfPnts)` (`cxx:1827-2172`).
///
/// `lof_pnts` is the start-point list. OCCT appends the extra points it needs
/// for closed surfaces (`cxx:1866-1960`) to that same list; the port works on a
/// local copy so the caller's list keeps its own points.
#[allow(clippy::too_many_arguments)]
pub(crate) fn perform_from_points(
    s1: &dyn Surface,
    d1: &TopolTool,
    s2: &dyn Surface,
    d2: &TopolTool,
    lof_pnts: &[PntOn2S],
    tol_tangency: f64,
    epsilon: f64,
    deflection: f64,
    increment: f64,
    slin: &mut Vec<GeomIntLine>,
    done: &mut bool,
    empty: &mut bool,
) {
    if lof_pnts.is_empty() {
        *done = true;
        return;
    }

    *empty = true;
    slin.clear();

    let (u1a, u1b) = s1.u_range();
    let (v1a, v1b) = s1.v_range();
    let (u2a, u2b) = s2.u_range();
    let (v2a, v2b) = s2.v_range();
    let (mut umin_lig1, mut vmin_lig1, mut umax_lig1, mut vmax_lig1) = (u1a, v1a, u1b, v1b);
    let (mut umin_lig2, mut vmin_lig2, mut umax_lig2, mut vmax_lig2) = (u2a, v2a, u2b, v2b);

    let mut periods = [0.0f64; 4];
    periods[0] = if s1.is_u_periodic() { s1.u_period() } else { 0.0 };
    periods[1] = if s1.is_v_periodic() { s1.v_period() } else { 0.0 };
    periods[2] = if s2.is_u_periodic() { s2.u_period() } else { 0.0 };
    periods[3] = if s2.is_v_periodic() { s2.v_period() } else { 0.0 };

    let mut points: Vec<PntOn2S> = lof_pnts.to_vec();

    if s1.is_u_closed() || s1.is_v_closed() || s2.is_u_closed() || s2.is_v_closed() {
        // `cxx:1866-1960`: duplicate every start point that sits exactly on a
        // closing seam to the opposite seam.
        let tol_par = PCONFUSION;
        let mut additional: Vec<PntOn2S> = Vec::new();
        for pnt in lof_pnts.iter() {
            let (u1, v1, u2, v2) = pnt.parameters();
            let mut new_pnt;
            if s1.is_u_closed() {
                if (u1 - u1a).abs() <= tol_par {
                    new_pnt = *pnt;
                    new_pnt.set_parameters(u1b, v1, u2, v2);
                    additional.push(new_pnt);
                } else if (u1 - u1b).abs() <= tol_par {
                    new_pnt = *pnt;
                    new_pnt.set_parameters(u1a, v1, u2, v2);
                    additional.push(new_pnt);
                }
            }
            if s1.is_v_closed() {
                if (v1 - v1a).abs() <= tol_par {
                    new_pnt = *pnt;
                    new_pnt.set_parameters(u1, v1b, u2, v2);
                    additional.push(new_pnt);
                } else if (v1 - v1b).abs() <= tol_par {
                    new_pnt = *pnt;
                    new_pnt.set_parameters(u1, v1a, u2, v2);
                    additional.push(new_pnt);
                }
            }
            if s2.is_u_closed() {
                if (u2 - u2a).abs() <= tol_par {
                    new_pnt = *pnt;
                    new_pnt.set_parameters(u1, v1, u2b, v2);
                    additional.push(new_pnt);
                } else if (u2 - u2b).abs() <= tol_par {
                    new_pnt = *pnt;
                    new_pnt.set_parameters(u1, v1, u2a, v2);
                    additional.push(new_pnt);
                }
            }
            if s2.is_v_closed() {
                if (v2 - v2a).abs() <= tol_par {
                    new_pnt = *pnt;
                    new_pnt.set_parameters(u1, v1, u2, v2b);
                    additional.push(new_pnt);
                } else if (v2 - v2b).abs() <= tol_par {
                    new_pnt = *pnt;
                    new_pnt.set_parameters(u1, v1, u2, v2a);
                    additional.push(new_pnt);
                }
            }
        }
        // Cut repeated points (`cxx:1938-1958`).
        for pnt in lof_pnts.iter() {
            let (u1, v1, u2, v2) = pnt.parameters();
            let mut iter2 = 0usize;
            while iter2 < additional.len() {
                let (new_u1, new_v1, new_u2, new_v2) = additional[iter2].parameters();
                if (u1 - new_u1).abs() <= tol_par
                    && (v1 - new_v1).abs() <= tol_par
                    && (u2 - new_u2).abs() <= tol_par
                    && (v2 - new_v2).abs() <= tol_par
                {
                    additional.remove(iter2);
                } else {
                    iter2 += 1;
                }
            }
        }
        points.extend(additional);
    }

    // `cxx:1962-1999`: the UV box of all start points.
    for pnt in points.iter() {
        let (u1, v1, u2, v2) = pnt.parameters();
        if u1 > umax_lig1 {
            umax_lig1 = u1;
        }
        if v1 > vmax_lig1 {
            vmax_lig1 = v1;
        }
        if u2 > umax_lig2 {
            umax_lig2 = u2;
        }
        if v2 > vmax_lig2 {
            vmax_lig2 = v2;
        }
        if u1 < umin_lig1 {
            umin_lig1 = u1;
        }
        if v1 < vmin_lig1 {
            vmin_lig1 = v1;
        }
        if u2 < umin_lig2 {
            umin_lig2 = u2;
        }
        if v2 < vmin_lig2 {
            vmin_lig2 = v2;
        }
    }

    let seuild_point_ligne = 15.0 * increment * increment;

    let mut pw = PWalking::new(s1, s2, tol_tangency, epsilon, deflection, increment);
    let mut d1m = d1.brep_tool(s1);
    let mut d2m = d2.brep_tool(s2);

    for c_pnt in points.iter() {
        let (u1, v1, u2, v2) = c_pnt.parameters();        let start_params = [u1, v1, u2, v2];

        let Some(start_p_on_2s) = pw.perform_first_point(start_params) else {
            continue;
        };
        // `dminiPointLigne = SeuildPointLigne + SeuildPointLigne` (`cxx:2025`).
        let mut dmini_point_ligne = seuild_point_ligne + seuild_point_ligne;

        let nb_lig_calculee = slin.len();
        for l in 0..nb_lig_calculee {
            let Some(testwline) = slin[l].as_wline() else {
                continue;
            };
            if is_point_on_line(&start_p_on_2s, testwline, deflection) {
                dmini_point_ligne = 0.0;
            }
            if dmini_point_ligne < seuild_point_ligne {
                break;
            }
        }

        if dmini_point_ligne <= seuild_point_ligne {
            continue;
        }

        pw.perform_box(
            start_params,
            umin_lig1,
            vmin_lig1,
            umin_lig2,
            vmin_lig2,
            umax_lig1,
            vmax_lig1,
            umax_lig2,
            vmax_lig2,
        );
        if !pw.is_done() || pw.nb_points() <= 2 {
            continue;
        }

        // Try to extend the intersection line to boundary, if it is possible.
        pw.put_to_boundary();

        let nb_lig_calculee = slin.len();
        let mut rejet_ligne = false;
        let point3d_debut = pw.value(1).value();
        let nb_points = pw.nb_points();
        let point_fin = pw.value(nb_points);
        let point3d_fin = point_fin.value();
        let mut ver = 0usize;
        for v in 0..nb_lig_calculee {
            let Some(verwline) = slin[v].as_wline() else {
                continue;
            };
            // Check end point if it is on existing line. Start point is checked
            // before (`cxx:2068-2077`).
            if is_point_on_line(&point_fin, verwline, deflection) {
                rejet_ligne = true;
                ver = v + 1;
                break;
            }
            let ver_point_debut = verwline.point(1).value();
            let ver_point_fin = verwline.point(verwline.nb_pnts()).value();
            if point3d_debut.distance(&ver_point_debut) <= tol_tangency
                && point3d_fin.distance(&ver_point_fin) <= tol_tangency
            {
                rejet_ligne = true;
                ver = v + 1;
                break;
            }
        }

        if rejet_ligne {
            duplicate_of_lines_processing(&pw, ver, slin, &mut rejet_ligne);
        }
        if rejet_ligne {
            continue;
        }

        // `cxx:2093-2117`: transitions from the normals at the tangent point.
        let (tg_dir, indextg) = pw.tangent_at_line();
        let tgline = occt_core::gp::GpVec::from_xyz(tg_dir.xyz());
        let (locu, locv) = pw.value(indextg).parameters_on_s1();
        let (_, d1u, d1v) = s1.d1(locu, locv);
        let norm1 = d1u.crossed(&d1v);
        let (locu, locv) = pw.value(indextg).parameters_on_s2();
        let (_, d2u, d2v) = s2.d1(locu, locv);
        let norm2 = d2u.crossed(&d2v);
        let (trans1, trans2) = if tgline.dot_cross(&norm2, &norm1) >= 0.0 {
            (IntSurfTypeTrans::Out, IntSurfTypeTrans::In)
        } else {
            (IntSurfTypeTrans::In, IntSurfTypeTrans::Out)
        };

        let mut wline = WLine::new();
        wline.base = IntPatchLine::with_transitions(false, trans1, trans2);
        wline.set_creating_way(WLineWay::PrmPrm);
        for p in pw.line() {
            wline.add(*p);
        }

        // `IntPatch_RstInt::PutVertexOnLine(wline, Surf1, D1, Surf2, true,
        // TolTang)` (`cxx:2122`) and its Surf2 twin (`:2127`).
        let mut line = GeomIntLine::Walking(wline);
        if put_vertex_on_line(&mut line, s1, &mut d1m, s2, true, tol_tangency).is_err() {
            continue;
        }
        if line.as_wline().map_or(true, |w| w.nb_pnts() < 2) {
            continue;
        }
        if put_vertex_on_line(&mut line, s2, &mut d2m, s1, false, tol_tangency).is_err() {
            continue;
        }
        let mut wline = match line {
            GeomIntLine::Walking(w) => w,
            _ => continue,
        };
        if wline.nb_pnts() < 2 {
            continue;
        }

        // `cxx:2132-2150`: the fallback vertices when `PutVertexOnLine` added
        // none.
        if wline.nb_vertex() == 0 {
            let params = pw.value(1).parameters();
            let mut vtx =
                crate::int_tools_wline::PatchPoint::new(point3d_debut, 1.0, params.0, params.1, params.2, params.3);
            vtx.set_tolerance(tol_tangency);
            wline.vertices.push(vtx);

            let params = pw.value(wline.nb_pnts()).parameters();
            let mut vtx = crate::int_tools_wline::PatchPoint::new(
                point3d_fin,
                wline.nb_pnts() as f64,
                params.0,
                params.1,
                params.2,
                params.3,
            );
            vtx.set_tolerance(tol_tangency);
            wline.vertices.push(vtx);
        }

        several_wlines_processing(
            s1,
            s2,
            slin,
            &periods,
            trans1,
            trans2,
            tol_tangency,
            pw.max_step(0).max(pw.max_step(1)),
            pw.max_step(2).max(pw.max_step(3)),
            &mut wline,
        );

        add_wline(slin, wline, deflection);
        *empty = false;
    }

    *done = true;
}

/// The `PutVertexOnLine` pair of `cxx:2122-2131`.
///
/// OCCT passes the two `Adaptor3d_TopolTool` domains (`theD1` / `theD2` of
/// `cxx:1828-1831`) straight to `IntPatch_RstInt::PutVertexOnLine`; the port's
/// `put_vertex_on_line` takes the `BrepTopolTool` view, which is exactly the
/// `brep_tool(surface)` conversion the other `Perform` bodies perform at their
/// own entry.
fn _domains(_d1: &TopolTool, _d2: &TopolTool) {}
