//! `IntPatch_PrmPrmIntersection` - two bi-parametric surfaces.
//!
//! Source: `ModelingAlgorithms/TKGeomAlgo/IntPatch/IntPatch_PrmPrmIntersection.cxx`.
//!
//! The main entry point is the polyhedron-polyhedron `Perform`
//! (`cxx:324-331` / `:337-350` / `:353-367` / `:369-928`):
//!
//! 1. build an `IntPatch_Polyhedron` on each surface (`:326`, `:343`),
//! 2. run `IntPatch_InterferencePolyhedron` on the pair (`:380`),
//! 3. for every section line and every tangent zone, turn the
//!    `Intf_SectionPoint` back into UV parameters (`SectionPointToParameters`,
//!    `:1584-1687`), walk an `IntWalk_PWalking` from there (`:604-614`),
//! 4. reject duplicates (`DublicateOfLinesProcessing`, `:278-312`), compute the
//!    3D transitions from the normals at the tangent point (`:651-676`), write
//!    the restriction vertices (`IntPatch_RstInt::PutVertexOnLine`, `:683-690`)
//!    and append the line (`AddWLine`, `:4085-4143`).
//!
//! The other `Perform` overloads of the `.cxx` are not wired here:
//! `ImplementPolyhedron` / `Perform(Surf1, Poly1, D1, ...)` self-interference
//! (`:932-1583`), the `LOfPnts` start-point list (`:1827-2175`), the single
//! start point (`:2176-2328`) and the grille-based variant (`:2488-3490`) are
//! still unported; `intpatch_prmprm_depart.rs` (the `PointDepart` helper,
//! `:3492-3991`) is kept for the `LOfPnts` overload.

use occt_core::gp::GpVec;
use occt_geom::Surface;

use crate::brep_top_adaptor_topol_tool::TopolTool as BrepTopolTool;
use crate::geom_int::{GeomIntLine, TopolTool};
use crate::int_tools_wline::{PatchPoint, WLine, WLineWay};
use crate::intpatch_interference_polyhedron::IntPatchInterferencePolyhedron;
use crate::intpatch_line::IntPatchLine;
use crate::intpatch_polyhedron::IntPatchPolyhedron;
use crate::intpatch_rst_int::put_vertex_on_line;
use crate::intsurf_transition::IntSurfTypeTrans;
use crate::intwalk::PWalking;

#[path = "intpatch_prmprm_t3bits.rs"]
mod prmprm_t3bits;
#[path = "intpatch_prmprm_depart.rs"]
#[allow(dead_code)]
mod prmprm_depart;
#[path = "intpatch_prmprm_section_point.rs"]
mod section_point;
#[path = "intpatch_prmprm_wline.rs"]
mod wline_helpers;

use section_point::section_point_to_parameters;
use wline_helpers::{add_wline, duplicate_of_lines_processing, is_point_on_line};

/// Parametric x parametric intersection.
#[derive(Clone)]
pub struct PrmPrmIntersection {
    done: bool,
    empty: bool,
    slin: Vec<GeomIntLine>,
}

impl PrmPrmIntersection {
    pub fn new() -> Self {
        Self {
            done: true,
            empty: true,
            slin: Vec::new(),
        }
    }

    /// `Perform(S1, D1, S2, D2, TolTangency, Epsilon, Deflection, Increment)`.
    ///
    /// `IntPatch_PrmPrmIntersection.cxx:324-331` (Poly1 from `D1`'s sample
    /// counts, `IntPatch_Polyhedron(Surf1, D1->NbSamplesU(), D1->NbSamplesV())`)
    /// and `:343` (`IntPatch_Polyhedron Poly2(Surf2)`), then `:369-928`.
    pub fn perform(
        &mut self,
        s1: &dyn Surface,
        d1: &TopolTool,
        s2: &dyn Surface,
        d2: &TopolTool,
        tol_tangency: f64,
        epsilon: f64,
        deflection: f64,
        increment: f64,
    ) {
        self.done = true;
        self.empty = true;
        self.slin.clear();

        let poly1 = IntPatchPolyhedron::with_sizes(s1, d1.nb_samples_u(s1), d1.nb_samples_v(s1));
        let poly2 = IntPatchPolyhedron::new(s2);
        let inter = IntPatchInterferencePolyhedron::with_polyhedra(&poly1, &poly2);

        let nb_lig_sec = inter.base.nb_section_lines();
        let nb_tan_zon = inter.base.nb_tangent_zones();

        // `IntWalk_PWalking PW(Surf1, Surf2, TolTangency, Epsilon, Deflection,
        // Increment)` (`cxx:395`).
        let mut pw = PWalking::new(s1, s2, tol_tangency, epsilon, deflection, increment);
        // `SeuildPointLigne = 15.0 * Increment * Increment` (`cxx:396`).
        let seuil = 15.0 * increment * increment;

        // `PutVertexOnLine` walks the domain's restriction arcs and needs the
        // `Adaptor3d_TopolTool` view OCCT hands to it (`cxx:683`, `:688`).
        let mut d1m = d1.brep_tool(s1);
        let mut d2m = d2.brep_tool(s2);

        if nb_lig_sec >= 1 {
            // Descending sort on `NumberOfPoints` (`cxx:410-433`).
            let mut tab_l: Vec<usize> = (1..=nb_lig_sec).collect();
            loop {
                let mut triok = true;
                for b in 1..tab_l.len() {
                    let nb_b = inter.base.line_value(tab_l[b] - 1).number_of_points();
                    let nb_a = inter.base.line_value(tab_l[b - 1] - 1).number_of_points();
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
                let line_sec = inter.base.line_value(line_index - 1);
                let nbp = line_sec.number_of_points();
                let nbp_i = nbp as i32;
                let mut tab_pt_dep = vec![0i32; nbp + 1];

                // `cxx:450-501`: UV box of the section line.
                let (mut u1, mut v1, mut u2, mut v2) = (0.0, 0.0, 0.0, 0.0);
                section_point_to_parameters(
                    line_sec.point(1),
                    &poly1,
                    &poly2,
                    &mut u1,
                    &mut v1,
                    &mut u2,
                    &mut v2,
                );
                let (mut umin1, mut vmin1, mut umin2, mut vmin2) = (u1, v1, u2, v2);
                let (mut umax1, mut vmax1, mut umax2, mut vmax2) = (u1, v1, u2, v2);
                for ilig in 2..=nbp {
                    section_point_to_parameters(
                        line_sec.point(ilig),
                        &poly1,
                        &poly2,
                        &mut u1,
                        &mut v1,
                        &mut u2,
                        &mut v2,
                    );
                    if u1 > umax1 {
                        umax1 = u1;
                    }
                    if v1 > vmax1 {
                        vmax1 = v1;
                    }
                    if u2 > umax2 {
                        umax2 = u2;
                    }
                    if v2 > vmax2 {
                        vmax2 = v2;
                    }
                    if u1 < umin1 {
                        umin1 = u1;
                    }
                    if v1 < vmin1 {
                        vmin1 = v1;
                    }
                    if u2 < umin2 {
                        umin2 = u2;
                    }
                    if v2 < vmin2 {
                        vmin2 = v2;
                    }
                }
                let bounds = [
                    umin1, vmin1, umin2, vmin2, umax1, vmax1, umax2, vmax2,
                ];

                // `cxx:503-576`: which point of the section line to start from.
                let mut nbps2: i32 = if nbp > 3 { (nbp / 2) as i32 } else { 1 };
                let mut n_dep = 0i32;
                let mut indice1 = 0i32;
                let mut indice2 = 0i32;
                let mut lignetrouvee = false;

                loop {
                    n_dep += 1;
                    if n_dep == 1 {
                        let mut incidence = 3.0;
                        let nbp1_4 = (nbp / 4) as i32;
                        let nbp3_4 = nbp_i - nbp1_4;

                        let mut nsp = (nbp / 2) as i32;
                        while nsp < nbp3_4 {
                            let current = line_sec.point(nsp.max(1) as usize).incidence();
                            if current < incidence {
                                nbps2 = nsp;
                                incidence = 0.9 * current;
                            }
                            nsp += 1;
                        }
                        nsp = (nbp / 2) as i32;
                        while nsp > nbp1_4 {
                            let current = line_sec.point(nsp.max(1) as usize).incidence();
                            if current < incidence {
                                nbps2 = nsp;
                                incidence = 0.9 * current;
                            }
                            nsp -= 1;
                        }
                        if nbp < 3 {
                            n_dep = 3;
                        }
                        indice1 = nbps2;
                    } else if n_dep == 2 {
                        if indice1 == 1 {
                            nbps2 = (nbp / 2) as i32;
                            indice2 = nbps2;
                        } else {
                            nbps2 = 1;
                            indice2 = 1;
                        }
                    } else if n_dep == 3 {
                        if indice1 == nbp_i {
                            nbps2 = (indice1 + indice2) / 2;
                        } else {
                            nbps2 = nbp_i;
                        }
                    } else {
                        nbps2 = n_dep - 3;
                        n_dep += 1;
                    }

                    // OCCT's `if (TabPtDep[nbps2] == 0) { ... }` (`cxx:578-720`).
                    let slot = nbps2 as usize;
                    if slot < tab_pt_dep.len() && tab_pt_dep[slot] == 0 {
                        tab_pt_dep[slot] = 1;
                        section_point_to_parameters(
                            line_sec.point(slot),
                            &poly1,
                            &poly2,
                            &mut u1,
                            &mut v1,
                            &mut u2,
                            &mut v2,
                        );
                        let start = [u1, v1, u2, v2];
                        if try_walking(
                            &mut pw,
                            start,
                            bounds,
                            EndMatch::Both,
                            tol_tangency,
                            deflection,
                            seuil,
                            s1,
                            s2,
                            &mut d1m,
                            &mut d2m,
                            &mut self.slin,
                        ) {
                            lignetrouvee = true;
                        }
                    }

                    // `do { ... } while (nbp > 5 && (((NombreDePoints... < 3) ||
                    // !lignetrouvee) && ((NombreDePoints... - 3 < nbp) ||
                    // lignetrouvee)))` (`cxx:721-723`).
                    let cont = nbp_i > 5
                        && ((n_dep < 3 || !lignetrouvee)
                            && (n_dep - 3 < nbp_i || lignetrouvee));
                    if !cont {
                        break;
                    }
                }
            }
        }

        // `cxx:730-778`: UV box of all tangent zones.
        let mut umin1 = f64::MAX;
        let mut vmin1 = f64::MAX;
        let mut umin2 = f64::MAX;
        let mut vmin2 = f64::MAX;
        let mut umax1 = -umin1;
        let mut vmax1 = -vmin1;
        let mut umax2 = -umin2;
        let mut vmax2 = -vmin2;
        let (mut u1, mut v1, mut u2, mut v2) = (0.0, 0.0, 0.0, 0.0);

        for z in 0..nb_tan_zon {
            let zone = inter.base.zone_value(z);
            for pz in 1..=zone.number_of_points() {
                section_point_to_parameters(
                    zone.point(pz),
                    &poly1,
                    &poly2,
                    &mut u1,
                    &mut v1,
                    &mut u2,
                    &mut v2,
                );
                if u1 > umax1 {
                    umax1 = u1;
                }
                if v1 > vmax1 {
                    vmax1 = v1;
                }
                if u2 > umax2 {
                    umax2 = u2;
                }
                if v2 > vmax2 {
                    vmax2 = v2;
                }
                if u1 < umin1 {
                    umin1 = u1;
                }
                if v1 < vmin1 {
                    vmin1 = v1;
                }
                if u2 < umin2 {
                    umin2 = u2;
                }
                if v2 < vmin2 {
                    vmin2 = v2;
                }
            }
        }
        let bounds = [
            umin1, vmin1, umin2, vmin2, umax1, vmax1, umax2, vmax2,
        ];

        // `cxx:781-926`.
        for z in 0..nb_tan_zon {
            let zone = inter.base.zone_value(z);
            for pz in 1..=zone.number_of_points() {
                section_point_to_parameters(
                    zone.point(pz),
                    &poly1,
                    &poly2,
                    &mut u1,
                    &mut v1,
                    &mut u2,
                    &mut v2,
                );
                let start = [u1, v1, u2, v2];
                // The tangent-zone pass tests the two end points with `||`
                // (`cxx:846-849`), the section-line pass with `&&` (`:637-640`).
                let _ = try_walking(
                    &mut pw,
                    start,
                    bounds,
                    EndMatch::Either,
                    tol_tangency,
                    deflection,
                    seuil,
                    s1,
                    s2,
                    &mut d1m,
                    &mut d2m,
                    &mut self.slin,
                );
            }
        }

        self.empty = self.slin.is_empty();
    }

    pub fn is_done(&self) -> bool {
        self.done
    }

    pub fn is_empty(&self) -> bool {
        self.empty
    }

    pub fn lines(&self) -> &[GeomIntLine] {
        &self.slin
    }
}

impl Default for PrmPrmIntersection {
    fn default() -> Self {
        Self::new()
    }
}

/// How the two end points of the walked line are compared with an existing
/// line (`cxx:637-640` vs `:846-849`).
#[derive(Clone, Copy, PartialEq, Eq)]
enum EndMatch {
    Both,
    Either,
}

/// `PW.PerformFirstPoint` + `PW.Perform` + duplicate test + `PutVertexOnLine` +
/// `AddWLine` (`cxx:585-718` and `:794-921`).
///
/// Returns `true` when a walking line was appended, i.e. OCCT's
/// `lignetrouvee = true`.
#[allow(clippy::too_many_arguments)]
fn try_walking(
    pw: &mut PWalking,
    start: [f64; 4],
    bounds: [f64; 8],
    end_match: EndMatch,
    tol_tangency: f64,
    deflection: f64,
    seuil: f64,
    s1: &dyn Surface,
    s2: &dyn Surface,
    d1: &mut BrepTopolTool,
    d2: &mut BrepTopolTool,
    slin: &mut Vec<GeomIntLine>,
) -> bool {
    // `HasStartPoint = PW.PerformFirstPoint(StartParams, StartPOn2S)`
    // (`cxx:585`).
    let Some(start_p) = pw.perform_first_point(start) else {
        return false;
    };

    // `dminiPointLigne` loop (`cxx:588-602`).
    let nb_lig = slin.len();
    let mut dmini = seuil + seuil;
    for l in 0..nb_lig {
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
        return false;
    }

    pw.perform_box(
        start, bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5], bounds[6],
        bounds[7],
    );
    if !pw.is_done() || pw.nb_points() <= 2 {
        return false;
    }

    let nb_points = pw.nb_points();
    let p_debut = pw.value(1).value();
    let p_fin = pw.value(nb_points).value();

    // `cxx:625-647` / `:831-853`.
    let mut rejet = false;
    let mut ver = 0usize;
    for v in 0..nb_lig {
        let Some(vwl) = slin[v].as_wline() else {
            continue;
        };
        // The end point is checked against existing lines; the start point was
        // checked above.
        if is_point_on_line(&pw.value(nb_points), vwl, deflection) {
            rejet = true;
            ver = v + 1;
            break;
        }
        let v_debut = vwl.point(1).value();
        let v_fin = vwl.point(vwl.nb_pnts()).value();
        let hit = match end_match {
            EndMatch::Both => {
                p_debut.distance(&v_debut) <= tol_tangency
                    && p_fin.distance(&v_fin) <= tol_tangency
            }
            EndMatch::Either => {
                p_debut.distance(&v_debut) < tol_tangency
                    || p_fin.distance(&v_fin) < tol_tangency
            }
        };
        if hit {
            rejet = true;
            ver = v + 1;
            break;
        }
    }

    if rejet {
        duplicate_of_lines_processing(pw, ver, slin, &mut rejet);
    }
    if rejet {
        return false;
    }

    // `cxx:651-676`: transitions from the normals at the tangent point.
    let (tg_dir, indextg) = pw.tangent_at_line();
    let tgline = GpVec::from_xyz(tg_dir.xyz());
    let (locu, locv) = pw.value(indextg).parameters_on_s1();
    let (_, d1u, d1v) = s1.d1(locu, locv);
    let norm1 = d1u.crossed(&d1v);
    let (locu, locv) = pw.value(indextg).parameters_on_s2();
    let (_, d2u, d2v) = s2.d1(locu, locv);
    let norm2 = d2u.crossed(&d2v);
    let (trans1, trans2) = if tgline.dot_cross(&norm2, &norm1) > 0.0 {
        (IntSurfTypeTrans::Out, IntSurfTypeTrans::In)
    } else {
        (IntSurfTypeTrans::In, IntSurfTypeTrans::Out)
    };

    // `new IntPatch_WLine(PW.Line(), false, trans1, trans2)` (`cxx:677-679`).
    let mut wl = WLine::new();
    wl.base = IntPatchLine::with_transitions(false, trans1, trans2);
    wl.set_creating_way(WLineWay::PrmPrm);
    for p in pw.line() {
        wl.add(*p);
    }

    let mut line = GeomIntLine::Walking(wl);
    // `PutVertexOnLine(wline, Surf1, D1, Surf2, true, TolTang)` (`cxx:683`).
    if put_vertex_on_line(&mut line, s1, d1, s2, true, tol_tangency).is_err() {
        return false;
    }
    if line.as_wline().map_or(true, |w| w.nb_pnts() < 2) {
        return false;
    }
    // `PutVertexOnLine(wline, Surf2, D2, Surf1, false, TolTang)` (`cxx:688`).
    if put_vertex_on_line(&mut line, s2, d2, s1, false, tol_tangency).is_err() {
        return false;
    }
    let mut wl = match line {
        GeomIntLine::Walking(w) => w,
        _ => return false,
    };
    if wl.nb_pnts() < 2 {
        return false;
    }

    // `cxx:692-711`: the fallback vertices when `PutVertexOnLine` added none.
    if wl.nb_vertex() == 0 {
        let params = pw.value(1).parameters();
        let mut vtx = PatchPoint::new(p_debut, 1.0, params.0, params.1, params.2, params.3);
        vtx.set_tolerance(tol_tangency);
        wl.vertices.push(vtx);

        let nbp = wl.nb_pnts();
        let params = pw.value(nbp).parameters();
        let mut vtx = PatchPoint::new(p_fin, nbp as f64, params.0, params.1, params.2, params.3);
        vtx.set_tolerance(tol_tangency);
        wl.vertices.push(vtx);
    }

    add_wline(slin, wl, deflection);
    true
}
