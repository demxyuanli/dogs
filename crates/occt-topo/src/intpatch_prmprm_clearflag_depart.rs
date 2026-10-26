//! `Perform(..., ClearFlag)` fallback branch
//! (`IntPatch_PrmPrmIntersection.cxx:3199-3402`).
//!
//! Taken when either domain's sample grid exceeds `Limit = 2500` (`:2524`).
//! `PointDepart` (`:3492-3991`) produces approximate start points on the
//! `15 * NbSamples` lattice; every one of them is walked with
//! `PW.Perform(StartParams)` (no UV box, `:3283`), and the result is only
//! checked against the already collected lines (`:3299-3321`) before being
//! appended (`:3391`). There is no `PutToBoundary`, no `SeekAdditionalPoints`
//! and no `SeveralWlinesProcessing` in this branch.

use occt_geom::Surface;

use crate::geom_int::{GeomIntLine, TopolTool};
use crate::int_tools_wline::PntOn2S;
use crate::intwalk::PWalking;

use super::wline_helpers::{
    add_wline, duplicate_of_lines_processing, finish_walking_line, is_point_on_line,
    line_transitions,
};

#[allow(clippy::too_many_arguments)]
pub(crate) fn walk_depart(
    points: &[PntOn2S],
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
    let seuil = 15.0 * increment * increment;

    let mut pw = PWalking::new(s1, s2, tol_tangency, epsilon, deflection, increment);
    let mut d1m = d1.brep_tool(s1);
    let mut d2m = d2.brep_tool(s2);

    for pnt in points {
        let (u1, v1, u2, v2) = pnt.parameters();
        let start = [u1, v1, u2, v2];

        let Some(start_p) = pw.perform_first_point(start) else {
            continue;
        };

        // `cxx:3262-3273`.
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

        // `PW.Perform(StartParams)` (`cxx:3283`): the natural UV box.
        pw.perform(start);
        if !pw.is_done() || pw.nb_points() <= 2 {
            continue;
        }

        let nbp_walk = pw.nb_points();
        let p_debut = pw.value(1).value();
        let p_fin = pw.value(nbp_walk).value();

        // `cxx:3299-3321`.
        let mut rejet = false;
        let mut ver = 0usize;
        for v in 0..slin.len() {
            let Some(vwl) = slin[v].as_wline() else {
                continue;
            };
            if is_point_on_line(&pw.value(nbp_walk), vwl, deflection) {
                rejet = true;
                ver = v + 1;
                break;
            }
            let v_debut = vwl.point(1).value();
            let v_fin = vwl.point(vwl.nb_pnts()).value();
            if p_debut.distance(&v_debut) < tol_tangency
                || p_fin.distance(&v_fin) < tol_tangency
            {
                rejet = true;
                ver = v + 1;
                break;
            }
        }

        // `DublicateOfLinesProcessing(PW, ver, SLin, RejetLigne)` (`cxx:3325`).
        if rejet {
            duplicate_of_lines_processing(&pw, ver, slin, &mut rejet);
        }
        if rejet {
            continue;
        }

        // `>` at `cxx:3343`.
        let (trans1, trans2) = line_transitions(&pw, s1, s2, true);
        if let Some(wline) = finish_walking_line(
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
            add_wline(slin, wline, deflection);
            *empty = false;
        }
    }
}
