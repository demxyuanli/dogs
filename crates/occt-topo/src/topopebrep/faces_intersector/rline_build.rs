//! The `TopOpeBRep_FacesIntersector` walking-line to restriction-line helpers.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRep/TopOpeBRep_FacesIntersector.cxx:650-1058`.
//!
//! Translation notes:
//! - `theDomain` is the port's [`TopolTool`] (`BRepTopAdaptor_TopolTool`);
//!   `theArc` is the `Adaptor2d_Curve2d` it hands out as `Value()`.
//! - The arc identity stored in a vertex / line is [`crate::int_tools_wline::arc_id`],
//!   the address of the `Arc<Curve2dAdaptor>` (`IntPatch_Point.hxx:229` compares
//!   handles by pointer).
//! - `theDomain->Initialize(theDomain->Value())` (`cxx:899`) only primes the
//!   vertex iterator (`BRepTopAdaptor_TopolTool.cxx:98-106`); `Value()` itself
//!   keeps returning the same table entry, so the arc comparand is unchanged.
//! - `theSurface->D1(u, v, ...)` in `TestWLineAlongRestriction` (`cxx:678`)
//!   feeds only the commented-out tolerance block (`cxx:679-684`) and has no
//!   observable effect; the surface argument is kept for signature fidelity.

use occt_geom::Surface;
use occt_geom2d::geom2d_int::find_parameter;

use crate::brep_adaptor_curve2d::Curve2dAdaptor;
use crate::brep_top_adaptor_topol_tool::TopolTool;
use crate::geom_int::RLine;
use crate::int_tools_wline::arc_id;
use crate::int_tools_wline::{PatchPoint, PntOn2S, TransType, WLine};
use crate::intpatch_line::IntPatchLine;
use crate::intpatch_vparams::rline::compute_vertex_parameters_rline;
use crate::intsurf_transition::IntSurfTypeTrans;

/// `TestWLineAlongRestriction(WLine, Rank, Surface, Domain, TolArc)` (`cxx:650-702`):
/// true when every sample of the walking line lands on a restriction of the
/// domain of side `rank`.
pub fn test_wline_along_restriction(
    wline: &WLine,
    rank: i32,
    _surface: Option<&dyn Surface>,
    domain: &mut TopolTool,
    tol_arc: f64,
) -> bool {
    let nb_pnts = wline.nb_pnts();
    let mut along = 0;
    for i in 1..=nb_pnts {
        let pmid: &PntOn2S = wline.point(i);
        let (u, v) = if rank == 1 {
            pmid.parameters_on_s1()
        } else {
            pmid.parameters_on_s2()
        };
        // `cxx:678-684`: the surface derivative call only feeds the
        // commented-out parametric tolerance.
        if domain.is_the_point_on(occt_core::gp::GpPnt2d::new(u, v), tol_arc, true) {
            along += 1;
        }
    }
    along == nb_pnts
}

/// `BuildRLineBasedOnWLine(WLine, Arc, Rank)` (`cxx:704-868`): the walking line
/// re-encoded as a restriction line along `arc`, or `None` when the two ends
/// collapse onto one arc parameter (a closed / degenerated edge) or when the
/// rank is not 1 / 2.
pub fn build_rline_based_on_wline(
    wline: &WLine,
    arc: &std::sync::Arc<Curve2dAdaptor>,
    rank: i32,
) -> Option<RLine> {
    if rank != 1 && rank != 2 {
        return None;
    }

    let (u, v) = if rank == 1 {
        wline.vertex(1).parameters_on_s1()
    } else {
        wline.vertex(1).parameters_on_s2()
    };
    let a_p_on_line = occt_core::gp::GpPnt2d::new(u, v);
    let par1 = find_parameter(
        arc.curve().expect("RLine arc without curve").as_ref(),
        &a_p_on_line,
        1.0e-7,
    );

    let nbvtx = wline.nb_vertex();
    let (u, v) = if rank == 1 {
        wline.vertex(nbvtx).parameters_on_s1()
    } else {
        wline.vertex(nbvtx).parameters_on_s2()
    };
    let a_p_on_line = occt_core::gp::GpPnt2d::new(u, v);
    let par2 = find_parameter(
        arc.curve().expect("RLine arc without curve").as_ref(),
        &a_p_on_line,
        1.0e-7,
    );

    let vtx1 = wline.vertex(1);
    let vtx2 = wline.vertex(nbvtx);
    let tol = if vtx1.tolerance > vtx2.tolerance {
        vtx1.tolerance
    } else {
        vtx2.tolerance
    };

    if let Some(res) = arc.resolution(tol) {
        if (par1 - par2).abs() < res {
            return None;
        }
    }

    let is_on_first = rank == 1;
    let token = arc_id(arc);

    // `new IntPatch_RLine(false, theWLine->TransitionOnS1(), theWLine->TransitionOnS2())`
    let mut an_rline = RLine {
        base: IntPatchLine::with_transitions(
            false,
            wline.base.transition_on_s1(),
            wline.base.transition_on_s2(),
        ),
        points: Vec::new(),
        vertices: Vec::new(),
        has_first_point: false,
        has_last_point: false,
        first_index: 0,
        last_index: 0,
        arc_on_s1: is_on_first,
        arc_on_s2: !is_on_first,
        c2d: arc.curve().cloned(),
        param_f: arc.first_parameter(),
        param_l: arc.last_parameter(),
        uv_arc: None,
        arc_edge: Some(arc.edge().clone()),
    };

    let undecided = TransType::Undecided;

    if par1 < par2 {
        let points: Vec<PntOn2S> = (1..=wline.nb_pnts()).map(|k| wline.point(k).clone()).collect();
        an_rline.add_points(points);

        let mut vtx_first = vtx1.clone();
        vtx_first.set_arc(is_on_first, token, par1, undecided);
        vtx_first.set_parameter(par1);
        an_rline.add_vertex(vtx_first);

        for k in 2..nbvtx {
            let mut vtx = wline.vertex(k).clone();
            let (u, v) = if rank == 1 {
                vtx.parameters_on_s1()
            } else {
                vtx.parameters_on_s2()
            };
            let atmp_point = occt_core::gp::GpPnt2d::new(u, v);
            let apar = find_parameter(
                arc.curve().expect("RLine arc without curve").as_ref(),
                &atmp_point,
                1.0e-7,
            );
            vtx.set_parameter(apar);
            an_rline.add_vertex(vtx);
        }

        let mut vtx_last = vtx2.clone();
        vtx_last.set_arc(is_on_first, token, par2, undecided);
        vtx_last.set_parameter(par2);
        an_rline.add_vertex(vtx_last);
        an_rline.set_first_point(1);
        an_rline.set_last_point(nbvtx);
        compute_vertex_parameters_rline(&mut an_rline);
    } else {
        let points: Vec<PntOn2S> = (1..=wline.nb_pnts())
            .rev()
            .map(|k| wline.point(k).clone())
            .collect();
        an_rline.add_points(points);

        let mut vtx_first = vtx2.clone();
        vtx_first.set_arc(is_on_first, token, par2, undecided);
        vtx_first.set_parameter(par2);
        an_rline.add_vertex(vtx_first);

        for k in (2..nbvtx).rev() {
            let mut vtx: PatchPoint = wline.vertex(k).clone();
            vtx.reverse_transition();
            let (u, v) = if rank == 1 {
                vtx.parameters_on_s1()
            } else {
                vtx.parameters_on_s2()
            };
            let atmp_point = occt_core::gp::GpPnt2d::new(u, v);
            let apar = find_parameter(
                arc.curve().expect("RLine arc without curve").as_ref(),
                &atmp_point,
                1.0e-7,
            );
            vtx.set_parameter(apar);
            an_rline.add_vertex(vtx);
        }

        let mut vtx_last = vtx1.clone();
        vtx_last.set_arc(is_on_first, token, par1, undecided);
        vtx_last.set_parameter(par1);
        an_rline.add_vertex(vtx_last);
        an_rline.set_first_point(1);
        an_rline.set_last_point(nbvtx);
        compute_vertex_parameters_rline(&mut an_rline);
    }

    Some(an_rline)
}

/// `BuildRLine(SeqOfWLine, Rank, Surface, Domain, TolArc)` (`cxx:874-1058`):
/// build one restriction line out of a chain of walking lines that all run
/// along the same domain arc. Returns `None` when no arc carries both chain end
/// vertices, when the ends are the same vertex (a closed edge), or when a
/// walking line strays off the arc.
pub fn build_rline(
    seq_of_wline: &[WLine],
    rank: i32,
    surface: Option<&dyn Surface>,
    domain: &mut TopolTool,
    tol_arc: f64,
) -> Option<RLine> {
    if seq_of_wline.is_empty() {
        return None;
    }
    let a_wline1 = &seq_of_wline[0];
    let a_wline2 = &seq_of_wline[seq_of_wline.len() - 1];
    let p1 = a_wline1.vertex(1);
    let p2 = a_wline2.vertex(a_wline2.nb_vertex());
    let a_v1 = p1.vertex_handle(rank == 1)?;
    let a_v2 = p2.vertex_handle(rank == 1)?;

    // avoid closed and degenerated edges
    if a_v1.is_same(a_v2) {
        return None;
    }

    let arcs: Vec<std::sync::Arc<Curve2dAdaptor>> = {
        domain.init();
        let mut out = Vec::new();
        while domain.more() {
            if let Some(c) = domain.value() {
                out.push(std::sync::Arc::new(c.clone()));
            }
            domain.next();
        }
        out
    };

    for arc_entry in arcs {
        domain.initialize_curve(&arc_entry);

        let mut found_vertex1 = false;
        let mut found_vertex2 = false;
        domain.init_vertex_iterator();
        while (!found_vertex1 || !found_vertex2) && domain.more_vertex() {
            if let Some(v) = domain.vertex() {
                if !found_vertex1 && a_v1.is_same(&v) {
                    found_vertex1 = true;
                }
                if !found_vertex2 && a_v2.is_same(&v) {
                    found_vertex2 = true;
                }
            }
            domain.next_vertex();
        }

        if !(found_vertex1 && found_vertex2) {
            continue;
        }

        let mut build_rline_ok = !seq_of_wline.is_empty();

        for wline in seq_of_wline.iter() {
            if !build_rline_ok {
                break;
            }
            let indexpnt = wline.nb_pnts() / 2;
            if indexpnt < 1 {
                build_rline_ok = false;
                break;
            }
            let p_on2s = wline.point(indexpnt);
            let (u, v) = if rank == 1 {
                p_on2s.parameters_on_s1()
            } else {
                p_on2s.parameters_on_s2()
            };
            let a_p_on_line = occt_core::gp::GpPnt2d::new(u, v);
            let par = find_parameter(
                arc_entry.curve().expect("RLine arc without curve").as_ref(),
                &a_p_on_line,
                1.0e-7,
            );
            let Some(a_p_on_arc) = arc_entry.d0(par) else {
                build_rline_ok = false;
                break;
            };
            let (nad1u, nad1v) = match surface.and_then(|s| {
                let (_, du, dv) = s.d1(u, v);
                Some((du.magnitude(), dv.magnitude()))
            }) {
                Some(v) => v,
                None => (0.0, 0.0),
            };
            let tolu = if nad1u > 1.0e-12 { tol_arc / nad1u } else { 0.1 };
            let tolv = if nad1v > 1.0e-12 { tol_arc / nad1v } else { 0.1 };
            let a_tolerance = if tolu > tolv { tolv } else { tolu };
            if a_p_on_arc.distance(&a_p_on_line) > a_tolerance {
                build_rline_ok = false;
            }
        }

        if !build_rline_ok {
            continue;
        }

        let mut trans1 = IntSurfTypeTrans::Undecided;
        let mut trans2 = IntSurfTypeTrans::Undecided;
        let mut a_line_on2_s: Vec<PntOn2S> = Vec::new();

        for wline in seq_of_wline.iter() {
            if wline.base.transition_on_s1() != IntSurfTypeTrans::Undecided
                && wline.base.transition_on_s1() != IntSurfTypeTrans::Touch
            {
                trans1 = wline.base.transition_on_s1();
            }
            if wline.base.transition_on_s2() != IntSurfTypeTrans::Undecided
                && wline.base.transition_on_s2() != IntSurfTypeTrans::Touch
            {
                trans2 = wline.base.transition_on_s2();
            }

            let param_min_on_line = wline.vertex(1).parameter_on_line() as i32;
            let param_max_on_line = wline.vertex(wline.nb_vertex()).parameter_on_line() as i32;
            for k in param_min_on_line..=param_max_on_line {
                a_line_on2_s.push(wline.point(k).clone());
            }
        }

        let mut emulated_wline = WLine::new();
        emulated_wline.base = IntPatchLine::with_transitions(false, trans1, trans2);
        emulated_wline.points = a_line_on2_s;

        let mut a_first_vertex = p1.clone();
        let mut a_last_vertex = p2.clone();
        a_first_vertex.set_parameter(1.0);
        a_last_vertex.set_parameter(emulated_wline.nb_pnts() as f64);

        emulated_wline.vertices.push(a_first_vertex);

        for apointindex in 2..=a_wline1.nb_vertex() {
            let mut a_point = a_wline1.vertex(apointindex).clone();
            let a_tolerance = if a_point.tolerance > p1.tolerance {
                a_point.tolerance
            } else {
                p1.tolerance
            };
            if a_point.p.is_equal_with_tolerance(&p1.p, a_tolerance) {
                a_point.set_parameter(1.0);
                emulated_wline.vertices.push(a_point);
            }
        }

        let nb_points = emulated_wline.nb_pnts() as f64;
        for apointindex in 1..a_wline2.nb_vertex() {
            let mut a_point = a_wline2.vertex(apointindex).clone();
            let a_tolerance = if a_point.tolerance > p2.tolerance {
                a_point.tolerance
            } else {
                p2.tolerance
            };
            if a_point.p.is_equal_with_tolerance(&p2.p, a_tolerance) {
                a_point.set_parameter(nb_points);
                emulated_wline.vertices.push(a_point);
            }
        }

        emulated_wline.vertices.push(a_last_vertex);

        return build_rline_based_on_wline(&emulated_wline, &arc_entry, rank);
    }

    None
}
