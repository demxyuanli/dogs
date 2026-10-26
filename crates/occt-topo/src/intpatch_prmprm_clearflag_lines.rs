//! `Perform(..., ClearFlag)` fast branch, section lines
//! (`IntPatch_PrmPrmIntersection.cxx:2553-2940`).
//!
//! Walks one `IntWalk_PWalking` line from each point of each
//! `IntPolyh_Intersection` section line, trying up to five start points per
//! line (`:2675-2937`), then extends and thickens the result with
//! `PutToBoundary` / `SeekAdditionalPoints` (`:2764-2779`) and splices it into
//! the already collected lines with `SeveralWlinesProcessing` (`:2913-2926`).

use occt_core::gp::GpPnt;
use occt_geom::Surface;

use crate::geom_int::{GeomIntLine, TopolTool};
use crate::intpolyh::IntPolyhIntersection;
use crate::intwalk::PWalking;

use super::several::several_wlines_processing;
use super::wline_helpers::{
    add_wline, duplicate_of_lines_processing, finish_walking_line, is_point_on_line,
    line_transitions,
};

/// `NbDePointsDeDepartDuChmLimit` (`cxx:2676`).
const DEPART_LIMIT: i32 = 5;

/// `aMinNbPoints` (`cxx:2775`, `:3068`).
const MIN_NB_POINTS: i32 = 40;

/// `double aDx <= 2. * Epsilon` (`cxx:2833`).
const LINE_TOL_FACTOR: f64 = 2.0;

#[allow(clippy::too_many_arguments)]
pub(crate) fn walk_section_lines(
    inter: &IntPolyhIntersection,
    s1: &dyn Surface,
    d1: &TopolTool,
    s2: &dyn Surface,
    d2: &TopolTool,
    tol_tangency: f64,
    epsilon: f64,
    deflection: f64,
    increment: f64,
    periods: &[f64; 4],
    slin: &mut Vec<GeomIntLine>,
    empty: &mut bool,
) {
    let nb_lig_sec = inter.nb_section_lines();
    if nb_lig_sec < 1 {
        return;
    }
    let seuil = 15.0 * increment * increment;

    // `IntWalk_PWalking PW(Surf1, Surf2, TolTangency, Epsilon, Deflection,
    // Increment)` (`cxx:2571`).
    let mut pw = PWalking::new(s1, s2, tol_tangency, epsilon, deflection, increment);
    let mut d1m = d1.brep_tool(s1);
    let mut d2m = d2.brep_tool(s2);

    // 1.1 descending sort on `NbPointsInLine` (`cxx:2575-2603`).
    let mut tab_l: Vec<i32> = (1..=nb_lig_sec).collect();
    loop {
        let mut triok = true;
        for b in 1..tab_l.len() {
            let nb_b = inter.nb_points_in_line(tab_l[b]);
            let nb_a = inter.nb_points_in_line(tab_l[b - 1]);
            if nb_b > nb_a {
                tab_l.swap(b, b - 1);
                triok = false;
            }
        }
        if triok {
            break;
        }
    }

    for &line_index in &tab_l {
        let nbp = inter.nb_points_in_line(line_index);
        if nbp == 0 {
            continue;
        }
        let mut tab_pt_dep = vec![0i32; (nbp + 1) as usize];

        // 1.2 UV bounds of the section line (`cxx:2617-2672`).
        let p0 = inter.get_line_point(line_index, 1);
        let (mut umin1, mut vmin1, mut umin2, mut vmin2) = (p0.u1, p0.v1, p0.u2, p0.v2);
        let (mut umax1, mut vmax1, mut umax2, mut vmax2) = (p0.u1, p0.v1, p0.u2, p0.v2);
        for ilig in 2..=nbp {
            let p = inter.get_line_point(line_index, ilig);
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

        // 1.3 start-point search (`cxx:2675-2937`).
        let mut n_dep = 0i32;
        let mut lignetrouvee = false;
        loop {
            n_dep += 1;
            let nbps2 = start_index(n_dep, nbp, &mut n_dep);
            if nbps2 < 1 || nbps2 > nbp {
                // Cannot happen while the `do` condition holds `nbp > 5`; the
                // `.cxx` has no guard and would write out of range.
                break;
            }
            let slot = nbps2 as usize;

            if tab_pt_dep[slot] == 0 {
                tab_pt_dep[slot] = 1;
                let p = inter.get_line_point(line_index, nbps2);
                let start = [p.u1, p.v1, p.u2, p.v2];

                let Some(start_p) = pw.perform_first_point(start) else {
                    if !keep_going(nbp, n_dep, lignetrouvee) {
                        break;
                    }
                    continue;
                };

                // `dminiPointLigne` loop over the lines already collected
                // (`cxx:2726-2742`).
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

                if dmini > seuil {
                    // `PW.Perform(StartParams, UminLig1, ...)` (`cxx:2745-2752`).
                    pw.perform_box(
                        start, umin1, vmin1, umin2, vmin2, umax1, vmax1, umax2, vmax2,
                    );
                    if pw.is_done() && pw.nb_points() > 2 {
                        // `PW.PutToBoundary(Surf1, Surf2)` (`cxx:2767`).
                        pw.put_to_boundary();
                        if pw.nb_points() < 3 {
                            if !keep_going(nbp, n_dep, lignetrouvee) {
                                break;
                            }
                            continue;
                        }
                        // `hasBeenAdded = PW.SeekAdditionalPoints(Surf1, Surf2,
                        // aMinNbPoints)` (`cxx:2777`).
                        let _has_been_added = if pw.nb_points() < MIN_NB_POINTS {
                            pw.seek_additional_points(MIN_NB_POINTS)
                        } else {
                            false
                        };

                        let nbp_walk = pw.nb_points();
                        let p_debut = pw.value(1).value();
                        let p_fin = pw.value(nbp_walk).value();

                        // A-posteriori check against the collected lines
                        // (`cxx:2781-2847`).
                        let mut rejet = false;
                        let mut ver = 0usize;
                        for v in 0..slin.len() {
                            let Some(vwl) = slin[v].as_wline() else {
                                continue;
                            };
                            let v_nbp = vwl.nb_pnts();
                            if v_nbp < 3 {
                                continue;
                            }
                            let v_debut = vwl.point(1).value();
                            let v_fin = vwl.point(v_nbp).value();
                            let d11 = p_debut.distance(&v_debut);
                            let d12 = p_debut.distance(&v_fin);
                            let d21 = p_fin.distance(&v_debut);
                            let d22 = p_fin.distance(&v_fin);
                            if (d11 <= tol_tangency && d22 <= tol_tangency)
                                || (d12 <= tol_tangency && d21 <= tol_tangency)
                            {
                                let mut mx = v_nbp / 2;
                                if v_nbp % 2 != 0 {
                                    mx += 1;
                                }
                                let p_mid = vwl.point(mx).value();
                                for m in 1..nbp_walk {
                                    let a1 = pw.value(m).value();
                                    let a2 = pw.value(m + 1).value();
                                    if let Some(dx) = line_distance(&a1, &a2, &p_mid) {
                                        if dx <= LINE_TOL_FACTOR * epsilon {
                                            rejet = true;
                                            ver = v;
                                            break;
                                        }
                                    }
                                }
                            }
                            if rejet {
                                break;
                            }
                        }

                        if rejet {
                            duplicate_of_lines_processing(&pw, ver, slin, &mut rejet);
                        }
                        if !rejet {
                            // `>=` at `cxx:2867`.
                            let (trans1, trans2) = line_transitions(&pw, s1, s2, false);
                            // **UNPORTED**: `wline->EnablePurging(!hasBeenAdded)`
                            // (`cxx:2882`) sets `myIsPurgerAllowed`, which is
                            // only read by `IntPatch_Intersection.cxx:1356`,
                            // `:1639` and `:2083`; the port does not implement
                            // those three purge passes, so the flag has no
                            // consumer yet.
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
                                lignetrouvee = true;
                                several_wlines_processing(
                                    s1,
                                    s2,
                                    slin,
                                    periods,
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
                        }
                    }
                }
            }

            // `do { ... } while (nbp > 5 && ...)` (`cxx:2931-2934`).
            if !keep_going(nbp, n_dep, lignetrouvee) {
                break;
            }
        }
    }
}

/// `NombreDePointsDeDepartDuCheminement - 3 < nbp || lignetrouvee`
/// (`cxx:2933`), combined with the `nbp > 5` test.
fn keep_going(nbp: i32, n_dep: i32, lignetrouvee: bool) -> bool {
    nbp > 5 && (n_dep < DEPART_LIMIT || !lignetrouvee) && (n_dep - 3 < nbp || lignetrouvee)
}

/// The `switch (NombreDePointsDeDepartDuCheminement)` of `cxx:2681-2708`.
///
/// The `default` arm increments the counter itself (`:2708`), and the `case 1`
/// arm jumps to the limit when the line has fewer than three points
/// (`:2687-2690`); both mutations are reproduced through `n_dep`.
fn start_index(n_dep: i32, nbp: i32, n_dep_mut: &mut i32) -> i32 {
    match n_dep {
        1 => {
            let v = if nbp > 1 { nbp / 2 } else { 1 };
            if nbp < 3 {
                *n_dep_mut = DEPART_LIMIT;
            }
            v
        }
        2 => 1,
        3 => nbp - 1,
        4 => 3 * nbp / 4,
        5 => nbp / 4,
        _ => {
            let v = n_dep - 3;
            *n_dep_mut += 1;
            v
        }
    }
}

/// `gp_Lin(aP1, gp_Dir(gp_Vec(aP1, aP2))).Distance(aPx)` (`cxx:2829-2832`),
/// i.e. the distance from `p` to the infinite line through `a1` along
/// `a1 -> a2`. `None` when the segment is degenerate, which the `.cxx` skips
/// with `if (aVec12.SquareMagnitude() < 1.e-20) continue;` (`:2824-2827`).
fn line_distance(a1: &GpPnt, a2: &GpPnt, p: &GpPnt) -> Option<f64> {
    let vx = a2.x() - a1.x();
    let vy = a2.y() - a1.y();
    let vz = a2.z() - a1.z();
    let sq = vx * vx + vy * vy + vz * vz;
    if sq < 1.0e-20 {
        return None;
    }
    let m = sq.sqrt();
    let dx = vx / m;
    let dy = vy / m;
    let dz = vz / m;
    let wx = p.x() - a1.x();
    let wy = p.y() - a1.y();
    let wz = p.z() - a1.z();
    let t = wx * dx + wy * dy + wz * dz;
    let hx = wx - t * dx;
    let hy = wy - t * dy;
    let hz = wz - t * dz;
    Some((hx * hx + hy * hy + hz * hz).sqrt())
}
