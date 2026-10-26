//! `IntPatch_PrmPrmIntersection::Perform(Surf1, Poly1, D1, TolTangency,
//! Epsilon, Deflection, Increment)` (`IntPatch_PrmPrmIntersection.cxx:932-1466`).
//!
//! One surface intersected with itself through a single `IntPatch_Polyhedron`:
//! `IntPatch_InterferencePolyhedron Interference(Poly1)` (`cxx:940`), then a
//! `IntWalk_PWalking PW(Surf1, Surf1, ...)` (`cxx:954`).
//!
//! Differences against the two-surface `Perform` (`cxx:369-928`), all copied
//! verbatim from the source:
//!
//! * the start point of a section line is the one that **maximises**
//!   `|U1 - U2| + |V1 - V2|` (`cxx:987-1004`), not `LineSec.GetPoint(i)
//!   .Incidence()`;
//! * there is no `TabPtDep` guard: every start point is retried;
//! * `PW.Perform` runs without the UV box clamp (`cxx:1107`, `:1320`);
//! * a start point is only kept when it is not on the diagonal of the
//!   parameter square (`cxx:1053`, `:1264`);
//! * the `dminiPointLigne` probe is the inline 2D box test of
//!   `cxx:1057-1098`, not `IsPointOnLine`;
//! * duplicate detection compares `Point(1)` / `Point(last)` only
//!   (`cxx:1114-1125`);
//! * the degenerate head / tail where both parameter pairs coincide is trimmed
//!   before the line is built (`cxx:1156-1187`);
//! * `PutVertexOnLine(wline, Surf1, D1, Surf1, true, TolTang)` runs once
//!   (`cxx:1217`);
//! * the line is appended directly (`cxx:1238`), without `AddWLine`.

use occt_core::gp::{GpPnt2d, GpVec};
use occt_geom::Surface;

use crate::brep_top_adaptor_topol_tool::TopolTool as BrepTopolTool;
use crate::geom_int::GeomIntLine;
use crate::int_tools_wline::{PatchPoint, WLine, WLineWay};
use crate::intpatch_interference_polyhedron::IntPatchInterferencePolyhedron;
use crate::intpatch_line::IntPatchLine;
use crate::intpatch_polyhedron::IntPatchPolyhedron;
use crate::intpatch_rst_int::put_vertex_on_line;
use crate::intsurf_transition::IntSurfTypeTrans;
use crate::intwalk::PWalking;

use super::section_point::section_point_to_parameters;
use super::wline_helpers::duplicate_of_lines_processing;

/// Tolerance of the diagonal test (`cxx:1053`, `cxx:1264`).
const DIAGONAL_TOL: f64 = 1.0e-7;
/// Tolerance of the degenerate head / tail trim (`cxx:1168`, `cxx:1185`).
const COINCIDENT_TOL: f64 = 1.0e-6;
/// Slack of the per-segment 2D box (`cxx:1083-1086`).
const SEGMENT_SLACK: f64 = 1.0e-7;

/// `Perform(Surf1, Poly1, D1, TolTangency, Epsilon, Deflection, Increment)`
/// (`cxx:932-1466`). The lines are appended to `slin`.
pub(crate) fn perform_self(
    slin: &mut Vec<GeomIntLine>,
    surf: &dyn Surface,
    poly: &IntPatchPolyhedron,
    d1: &mut BrepTopolTool,
    tol_tangency: f64,
    epsilon: f64,
    deflection: f64,
    increment: f64,
) {
    let inter = IntPatchInterferencePolyhedron::with_self(poly);
    let nb_lig_sec = inter.base.nb_section_lines();
    let nb_tan_zon = inter.base.nb_tangent_zones();

    // `IntWalk_PWalking PW(Surf1, Surf1, ...)` (`cxx:954`).
    let mut pw = PWalking::new(surf, surf, tol_tangency, epsilon, deflection, increment);
    // `SeuildPointLigne = 15.0 * Increment * Increment` (`cxx:957`).
    let seuil = 15.0 * increment * increment;

    // `cxx:969-1248`: one start point per section line.
    for ls in 1..=nb_lig_sec {
        let line_sec = inter.base.line_value(ls - 1);
        let nbp = line_sec.number_of_points();
        let mut nbps2: usize = if nbp > 3 { nbp / 2 } else { 1 };
        let mut nb_dep = 0i32;
        let mut ind1 = 0usize;
        let mut ind2 = 0usize;

        loop {
            nb_dep += 1;
            if nb_dep == 1 {
                let mut incidence = 0.0f64;
                // `for (nsp1 = nbp/2; nsp1 >= 1; nsp1--)`
                let mut nsp1 = nbp / 2;
                while nsp1 >= 1 {
                    let (u1, v1, u2, v2) = point_parameters(line_sec.point(nsp1), poly);
                    let cur = (u1 - u2).abs() + (v1 - v2).abs();
                    if cur > incidence {
                        nbps2 = nsp1;
                        incidence = cur;
                    }
                    if nsp1 == 1 {
                        break;
                    }
                    nsp1 -= 1;
                }
                // `for (nsp1 = nbp/2; nsp1 <= nbp; nsp1++)`
                let mut nsp1 = nbp / 2;
                while nsp1 <= nbp {
                    let (u1, v1, u2, v2) = point_parameters(line_sec.point(nsp1), poly);
                    let cur = (u1 - u2).abs() + (v1 - v2).abs();
                    if cur > incidence {
                        nbps2 = nsp1;
                        incidence = cur;
                    }
                    nsp1 += 1;
                }
                if nbp < 3 {
                    nb_dep = 3;
                }
                ind1 = nbps2;
            } else if nb_dep == 2 {
                if ind1 == 1 {
                    nbps2 = nbp / 2;
                    ind2 = nbps2;
                } else {
                    nbps2 = 1;
                    ind2 = 1;
                }
            } else if ind1 == nbp {
                nbps2 = (ind1 + ind2) / 2;
            } else {
                nbps2 = nbp;
            }

            let (u1, v1, u2, v2) = point_parameters(line_sec.point(nbps2), poly);
            attempt(
                slin,
                &mut pw,
                [u1, v1, u2, v2],
                surf,
                d1,
                tol_tangency,
                seuil,
            );

            // `while (nbp > 5 && NombreDePointsDeDepartDuCheminement < 3)`.
            if !(nbp > 5 && nb_dep < 3) {
                break;
            }
        }
    }

    // `cxx:1250-1465`: every point of every tangent zone.
    for z in 0..nb_tan_zon {
        let zone = inter.base.zone_value(z);
        for pz in 1..=zone.number_of_points() {
            let (u1, v1, u2, v2) = point_parameters(zone.point(pz), poly);
            attempt(
                slin,
                &mut pw,
                [u1, v1, u2, v2],
                surf,
                d1,
                tol_tangency,
                seuil,
            );
        }
    }
}

/// `SectionPointToParameters(Sp, Poly1, Poly1, ...)` with fresh outputs.
fn point_parameters(sp: &occt_core::intf::IntfSectionPoint, poly: &IntPatchPolyhedron) -> (f64, f64, f64, f64) {
    let (mut u1, mut v1, mut u2, mut v2) = (0.0, 0.0, 0.0, 0.0);
    section_point_to_parameters(sp, poly, poly, &mut u1, &mut v1, &mut u2, &mut v2);
    (u1, v1, u2, v2)
}

/// One start point: `cxx:1046-1246` (section lines) / `:1256-1463` (zones).
#[allow(clippy::too_many_arguments)]
fn attempt(
    slin: &mut Vec<GeomIntLine>,
    pw: &mut PWalking,
    start: [f64; 4],
    surf: &dyn Surface,
    d1: &mut BrepTopolTool,
    tol_tangency: f64,
    seuil: f64,
) {
    let Some(start_p) = pw.perform_first_point(start) else {
        return;
    };
    let (pu1, pv1, pu2, pv2) = start_p.parameters();
    // `if (Abs(pu1-pu2) > 1e-7 || Abs(pv1-pv2) > 1e-7)` (`cxx:1053`).
    if (pu1 - pu2).abs() <= DIAGONAL_TOL && (pv1 - pv2).abs() <= DIAGONAL_TOL {
        return;
    }

    // `dminiPointLigne` probe (`cxx:1056-1098` / `:1271-1315`).
    let nb_lig = slin.len();
    let mut dmini = seuil + seuil;
    let p1 = GpPnt2d::new(pu1, pv1);
    let p2 = GpPnt2d::new(pu2, pv2);
    let p3 = start_p.value();
    for l in 0..nb_lig {
        if dmini < seuil {
            break;
        }
        let Some(wl) = slin[l].as_wline() else {
            continue;
        };
        if wl.is_out_surf1_box(&p1) || wl.is_out_surf2_box(&p2) || wl.is_out_box(&p3) {
            continue;
        }
        let nb_pnts = wl.nb_pnts();
        let mut ll = 1;
        while ll < nb_pnts && dmini >= seuil {
            let a = *wl.point(ll);
            let b = *wl.point(ll + 1);
            let (a_u1, a_v1, _, _) = a.parameters();
            let (b_u1, b_v1, _, _) = b.parameters();
            // `Test segment (Au1,Av1)-(Bu1,Bv1)` with the 1e-7 slack
            // (`cxx:1070-1087`). Both probes below use the S1 box, as the
            // source does.
            let au1 = a_u1.min(b_u1) - SEGMENT_SLACK;
            let av1 = a_v1.min(b_v1) - SEGMENT_SLACK;
            let bu1 = a_u1.max(b_u1) + SEGMENT_SLACK;
            let bv1 = a_v1.max(b_v1) + SEGMENT_SLACK;
            let inside = |u: f64, v: f64| u >= au1 && u <= bu1 && v >= av1 && v <= bv1;
            if inside(pu1, pv1) || inside(pu2, pv2) {
                dmini = 0.0;
            }
            ll += 1;
        }
    }
    if dmini > seuil {
        // `PW.Perform(StartParams)` without the UV box (`cxx:1107`).
        pw.perform(start);
        if pw.is_done() && pw.nb_points() > 2 {
            build_line(slin, pw, surf, d1, tol_tangency);
        }
    }
}

/// `cxx:1108-1240` / `:1321-1453`: duplicate test, trim, `PutVertexOnLine`,
/// fallback vertices, append.
fn build_line(
    slin: &mut Vec<GeomIntLine>,
    pw: &PWalking,
    surf: &dyn Surface,
    d1: &mut BrepTopolTool,
    tol_tangency: f64,
) {
    let nb_points = pw.nb_points();
    let p_debut = pw.value(1).value();
    let p_fin = pw.value(nb_points).value();

    // `for (ver = 1; ver <= NbLigCalculee; ver++)` (`cxx:1112-1126`).
    let nb_lig = slin.len();
    let mut rejet = false;
    let mut ver = 0usize;
    for v in 0..nb_lig {
        let Some(vwl) = slin[v].as_wline() else {
            continue;
        };
        let v_debut = vwl.point(1).value();
        let v_fin = vwl.point(vwl.nb_pnts()).value();
        if p_debut.distance(&v_debut) < tol_tangency || p_fin.distance(&v_fin) < tol_tangency {
            rejet = true;
            ver = v + 1;
            break;
        }
    }
    if rejet {
        duplicate_of_lines_processing(pw, ver, slin, &mut rejet);
    }
    if rejet {
        return;
    }

    // Transitions (`cxx:1131-1153`).
    let (tg_dir, indextg) = pw.tangent_at_line();
    let tgline = GpVec::from_xyz(tg_dir.xyz());
    let (locu, locv) = pw.value(indextg).parameters_on_s1();
    let (_, d1u, d1v) = surf.d1(locu, locv);
    let norm1 = d1u.crossed(&d1v);
    let (locu, locv) = pw.value(indextg).parameters_on_s2();
    let (_, d2u, d2v) = surf.d1(locu, locv);
    let norm2 = d2u.crossed(&d2v);
    let (trans1, trans2) = if tgline.dot_cross(&norm2, &norm1) > 0.0 {
        (IntSurfTypeTrans::Out, IntSurfTypeTrans::In)
    } else {
        (IntSurfTypeTrans::In, IntSurfTypeTrans::Out)
    };

    // `cxx:1156-1187`: the `imin` / `imax` trim of the degenerate head / tail.
    let line = pw.line();
    let mut i = 1i32;
    let imin = loop {
        let (u1, v1, u2, v2) = line[(i - 1) as usize].parameters();
        let same = (u1 - u2).abs() <= COINCIDENT_TOL && (v1 - v2).abs() <= COINCIDENT_TOL;
        if !(i < nb_points && same) {
            break i;
        }
        i += 1;
    };
    let imin = if imin > 2 { imin - 1 } else { imin };

    let mut i = nb_points + 1;
    let imax = loop {
        i -= 1;
        let (u1, v1, u2, v2) = line[(i - 1) as usize].parameters();
        let same = (u1 - u2).abs() <= COINCIDENT_TOL && (v1 - v2).abs() <= COINCIDENT_TOL;
        if !(i > 2 && same) {
            break i;
        }
    };
    let imax = if imax < nb_points { imax + 1 } else { imax };
    if imin >= imax {
        return;
    }

    let mut wl = WLine::new();
    wl.base = IntPatchLine::with_transitions(false, trans1, trans2);
    wl.set_creating_way(WLineWay::PrmPrm);
    for k in imin..=imax {
        wl.add(line[(k - 1) as usize]);
    }

    // `cxx:1203-1213`: close the loop with a vertex on the last point.
    let p_on_2s_deb = *wl.point(1);
    let p_on_2s_fin = *wl.point(wl.nb_pnts());
    if p_on_2s_deb.value().distance(&p_on_2s_fin.value()) <= tol_tangency {
        let (u1t, v1t, u2t, v2t) = p_on_2s_deb.parameters();
        let mut vtx = PatchPoint::new(p_on_2s_deb.value(), wl.nb_pnts() as f64, u1t, v1t, u2t, v2t);
        vtx.set_value(p_on_2s_deb.value(), tol_tangency, false);
        // `vtx.SetParameters(u2t, v2t, u1t, v1t)` (`cxx:1210`): the source swaps
        // the two parameter pairs, copied verbatim.
        vtx.set_parameters(u2t, v2t, u1t, v1t);
        vtx.set_parameter(wl.nb_pnts() as f64);
        let n = wl.nb_pnts();
        *wl.point_mut(n) = vtx.pnt_on_2s();
    }

    // `PutVertexOnLine(wline, Surf1, D1, Surf1, true, TolTang)` (`cxx:1217`).
    let mut gline = GeomIntLine::Walking(wl);
    if put_vertex_on_line(&mut gline, surf, d1, surf, true, tol_tangency).is_err() {
        return;
    }
    let mut wl = match gline {
        GeomIntLine::Walking(w) => w,
        _ => return,
    };
    if wl.nb_pnts() < 2 {
        return;
    }

    // `cxx:1221-1237`: fallback vertices when none were attached.
    if wl.nb_vertex() == 0 {
        let params = line[0].parameters();
        let mut vtx = PatchPoint::new(p_debut, 1.0, params.0, params.1, params.2, params.3);
        vtx.set_tolerance(tol_tangency);
        wl.vertices.push(vtx);

        let n = wl.nb_pnts();
        let params = line[(n - 1) as usize].parameters();
        let mut vtx = PatchPoint::new(p_fin, n as f64, params.0, params.1, params.2, params.3);
        vtx.set_tolerance(tol_tangency);
        wl.vertices.push(vtx);
    }

    slin.push(GeomIntLine::Walking(wl));
}
