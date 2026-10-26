//! `Perform(..., ClearFlag)` fast branch, tangent zones
//! (`IntPatch_PrmPrmIntersection.cxx:2950-3188`).
//!
//! The zone pass differs from the section-line pass in three ways the `.cxx`
//! spells out: it computes one global UV box for all zones (`:2952-3001`), it
//! calls `SeekAdditionalPoints` but not `PutToBoundary` (`:3065-3074`), and it
//! neither splices with `SeveralWlinesProcessing` nor drops duplicates through
//! `DublicateOfLinesProcessing`.

use occt_geom::Surface;

use crate::geom_int::{GeomIntLine, TopolTool};
use crate::intpolyh::IntPolyhIntersection;
use crate::intwalk::PWalking;

use super::wline_helpers::{add_wline, finish_walking_line, is_point_on_line, line_transitions};

/// `aMinNbPoints` (`cxx:3068`).
const MIN_NB_POINTS: i32 = 40;

#[allow(clippy::too_many_arguments)]
pub(crate) fn walk_tangent_zones(
    inter: &IntPolyhIntersection,
    s1: &dyn Surface,
    d1: &TopolTool,
    s2: &dyn Surface,
    d2: &TopolTool,
    tol_tangency: f64,
    epsilon: f64,
    deflection: f64,
    increment: f64,
    slin: &mut Vec<GeomIntLine>,
    empty: &mut bool,
) {
    let nb_tan_zon = inter.nb_tangent_zones();
    if nb_tan_zon < 1 {
        return;
    }
    let seuil = 15.0 * increment * increment;

    let mut pw = PWalking::new(s1, s2, tol_tangency, epsilon, deflection, increment);
    let mut d1m = d1.brep_tool(s1);
    let mut d2m = d2.brep_tool(s2);

    // `cxx:2952-3001`: UV box of all tangent-zone points. `NbPointsInTangentZone`
    // is always 1, so only `pz = 1` is visited (`:2963`).
    let (mut umin1, mut vmin1, mut umin2, mut vmin2) = (f64::MAX, f64::MAX, f64::MAX, f64::MAX);
    let (mut umax1, mut vmax1, mut umax2, mut vmax2) =
        (-f64::MAX, -f64::MAX, -f64::MAX, -f64::MAX);
    for z in 1..=nb_tan_zon {
        let p = inter.get_tangent_zone_point(z);
        if p.u1 > umax1 {
            umax1 = p.u1;
        }
        if p.v1 > vmax1 {
            vmax1 = p.v1;
        }
        if p.u2 > umax2 {
            umax2 = p.u2;
        }
        if p.v2 > vmax2 {
            vmax2 = p.v2;
        }
        if p.u1 < umin1 {
            umin1 = p.u1;
        }
        if p.v1 < vmin1 {
            vmin1 = p.v1;
        }
        if p.u2 < umin2 {
            umin2 = p.u2;
        }
        if p.v2 < vmin2 {
            vmin2 = p.v2;
        }
    }

    // `cxx:3003-3186`.
    for z in 1..=nb_tan_zon {
        let p = inter.get_tangent_zone_point(z);
        let start = [p.u1, p.v1, p.u2, p.v2];

        let Some(start_p) = pw.perform_first_point(start) else {
            continue;
        };

        let mut dmini = seuil + seuil;
        for l in 0..slin.len() {
            let Some(wl) = slin[l].as_wline() else {
                continue;
            };
            if is_point_on_line(&start_p, wl, deflection) {
                dmini = 0.0;
            }
            if dmini < seuil {
                break;
            }
        }
        if dmini <= seuil {
            continue;
        }

        pw.perform_box(start, umin1, vmin1, umin2, vmin2, umax1, vmax1, umax2, vmax2);
        if !pw.is_done() || pw.nb_points() <= 2 {
            continue;
        }

        // `cxx:3065-3074`: no `PutToBoundary` in this pass.
        let has_been_added = if pw.nb_points() < MIN_NB_POINTS {
            pw.seek_additional_points(MIN_NB_POINTS)
        } else {
            false
        };

        let nbp_walk = pw.nb_points();
        let p_debut = pw.value(1).value();
        let p_fin = pw.value(nbp_walk).value();

        // `cxx:3079-3114`.
        let mut rejet = false;
        for v in 0..slin.len() {
            let Some(vwl) = slin[v].as_wline() else {
                continue;
            };
            // `IsPointOnLine(PointFin, verwline, Deflection)` (`:3096-3100`).
            if is_point_on_line(&pw.value(nbp_walk), vwl, deflection) {
                rejet = true;
                break;
            }
            let v_debut = vwl.point(1).value();
            let v_fin = vwl.point(vwl.nb_pnts()).value();
            if p_debut.distance(&v_debut) < tol_tangency
                || p_fin.distance(&v_fin) < tol_tangency
            {
                rejet = true;
                break;
            }
        }
        if rejet {
            continue;
        }

        // `>` at `cxx:3128`.
        let (trans1, trans2) = line_transitions(&pw, s1, s2, true);
        if let Some(mut wline) = finish_walking_line(
            &pw,
            s1,
            &mut d1m,
            s2,
            &mut d2m,
            trans1,
            trans2,
            p_debut,
            p_fin,
            tol_tangency,
        ) {
            // `wline->EnablePurging(!hasBeenAdded)` (`cxx:3143`), same order
            // argument as in `intpatch_prmprm_clearflag_lines.rs`.
            wline.enable_purging(!has_been_added);
            add_wline(slin, wline, deflection);
            *empty = false;
        }
    }
}
