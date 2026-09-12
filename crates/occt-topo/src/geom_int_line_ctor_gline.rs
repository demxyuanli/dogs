//! GLine / Circle / Ellipse branch of `GeomInt_LineConstructor`.

use std::f64::consts::PI;

use occt_core::elib::clib;
use occt_core::gp::{GpPnt, GpPnt2d};
use occt_core::precision::PCONFUSION;

use crate::fclass2d::FaceState;
use crate::geom_int::line_tool;
use crate::geom_int::quadric::{adjust_periodic_uv, surface_parameters};
use crate::geom_int::types::{GLineKind, GeomIntLine, IntPatchIType};
use crate::geom_int::LineConstructor;
use crate::int_tools_wline::PatchPoint;

const TWO_PI: f64 = PI + PI;

pub(crate) fn perform_gline(ctor: &mut LineConstructor, line: &GeomIntLine, tol: f64) {
    ctor.seqp_mut().clear();
    let typl = line.arc_type();
    if typl == IntPatchIType::Circle || typl == IntPatchIType::Ellipse {
        treat_circle(ctor, line, tol);
        ctor.set_done(true);
        return;
    }
    let nbvtx = line_tool::nb_vertex(line);
    let mut intrvtested = false;
    let Some(gline) = line.as_gline() else {
        ctor.set_done(true);
        return;
    };
    let Some((s1, s2)) = ctor.surface_arcs() else {
        ctor.set_done(false);
        return;
    };
    let Some((d1, d2)) = ctor.domain_clones() else {
        ctor.set_done(false);
        return;
    };
    let s1 = s1.as_ref();
    let s2 = s2.as_ref();
    for i in 1..nbvtx {
        let firstp = line_tool::vertex(line, i).parameter_on_line();
        let lastp = line_tool::vertex(line, i + 1).parameter_on_line();
        if (firstp - lastp).abs() > PCONFUSION {
            intrvtested = true;
            let pmid = (firstp + lastp) * 0.5;
            if let Some(p) = gline_point(typl, gline, pmid) {
                if let Some((mut u1, mut v1)) = surface_parameters(s1, &p) {
                    if let Some((mut u2, mut v2)) = surface_parameters(s2, &p) {
                        adjust_periodic_uv(s1, s2, &mut u1, &mut v1, &mut u2, &mut v2);
                        if d1.classify(GpPnt2d::new(u1, v1), tol) != FaceState::Out
                            && d2.classify(GpPnt2d::new(u2, v2), tol) != FaceState::Out
                        {
                            ctor.seqp_mut().push(firstp);
                            ctor.seqp_mut().push(lastp);
                        }
                    }
                }
            }
        }
    }
    if !intrvtested {
        ctor.seqp_mut().push(line_tool::first_parameter(line));
        ctor.seqp_mut().push(line_tool::last_parameter(line));
    }
    ctor.set_done(true);
}

fn treat_circle(ctor: &mut LineConstructor, line: &GeomIntLine, the_tol: f64) {
    let a_type = line.arc_type();
    let Some(a_gline) = line.as_gline() else {
        return;
    };
    if reject_micro_circle(a_gline, a_type, the_tol) {
        return;
    }
    let a_nb_vtx = a_gline.nb_vertex();
    if a_nb_vtx <= 0 {
        ctor.seqp_mut().push(line_tool::first_parameter(line));
        ctor.seqp_mut().push(line_tool::last_parameter(line));
        return;
    }
    let mut a_vtx: Vec<PatchPoint> = (1..=a_nb_vtx).map(|i| *a_gline.vertex(i)).collect();
    a_vtx.sort_by(|a, b| a.param_on_line.total_cmp(&b.param_on_line));
    let a_min_prm = a_vtx.first().map(|v| v.param_on_line + TWO_PI).unwrap_or(TWO_PI);
    a_vtx.push({
        let mut last = a_vtx[0];
        last.param_on_line = a_min_prm;
        last
    });
    reject_duplicates(&mut a_vtx);
    a_vtx.sort_by(|a, b| a.param_on_line.total_cmp(&b.param_on_line));
    let Some((s1, s2)) = ctor.surface_arcs() else {
        return;
    };
    let Some((d1, d2)) = ctor.domain_clones() else {
        return;
    };
    let s1 = s1.as_ref();
    let s2 = s2.as_ref();
    for i in 0..a_vtx.len().saturating_sub(1) {
        let a_t1 = a_vtx[i].param_on_line;
        let a_t2 = a_vtx[i + 1].param_on_line;
        if a_t2 == f64::MAX {
            break;
        }
        let a_tmid = (a_t1 + a_t2) * 0.5;
        let Some(a_pmid) = gline_point(a_type, a_gline, a_tmid) else {
            continue;
        };
        let Some((mut a_u1, mut a_v1)) = surface_parameters(s1, &a_pmid) else {
            continue;
        };
        let Some((mut a_u2, mut a_v2)) = surface_parameters(s2, &a_pmid) else {
            continue;
        };
        adjust_periodic_uv(s1, s2, &mut a_u1, &mut a_v1, &mut a_u2, &mut a_v2);
        if d1.classify(GpPnt2d::new(a_u1, a_v1), the_tol) != FaceState::Out
            && d2.classify(GpPnt2d::new(a_u2, a_v2), the_tol) != FaceState::Out
        {
            ctor.seqp_mut().push(a_t1);
            ctor.seqp_mut().push(a_t2);
        }
    }
}

pub(crate) fn gline_point(
    typl: IntPatchIType,
    gline: &crate::geom_int::types::GLine,
    a_t: f64,
) -> Option<GpPnt> {
    Some(match (&gline.kind, typl) {
        (GLineKind::Lin(l), IntPatchIType::Lin) => clib::line_value(l, a_t),
        (GLineKind::Circ(c), IntPatchIType::Circle) => clib::circle_value(c, a_t),
        (GLineKind::Elips(e), IntPatchIType::Ellipse) => clib::ellipse_value(e, a_t),
        (GLineKind::Hypr(h), IntPatchIType::Hyperbola) => clib::hyperbola_value(h, a_t),
        (GLineKind::Parab(p), IntPatchIType::Parabola) => clib::parabola_value(p, a_t),
        _ => return None,
    })
}

fn reject_micro_circle(
    a_gline: &crate::geom_int::types::GLine,
    a_type: IntPatchIType,
    a_tol3d: f64,
) -> bool {
    match (a_type, &a_gline.kind) {
        (IntPatchIType::Circle, GLineKind::Circ(c)) => c.radius() < a_tol3d,
        (IntPatchIType::Ellipse, GLineKind::Elips(e)) => e.major_radius() < a_tol3d,
        _ => false,
    }
}

fn reject_duplicates(the_vtx: &mut [PatchPoint]) {
    let a_tol_pc = 1000.0 * PCONFUSION;
    if the_vtx.len() < 2 {
        return;
    }
    let upper = the_vtx.len();
    for i in 0..upper.saturating_sub(2) {
        let a_prmi = the_vtx[i].param_on_line;
        if a_prmi == f64::MAX {
            continue;
        }
        for j in (i + 1)..(upper - 1) {
            let a_prmj = the_vtx[j].param_on_line;
            if a_prmj - a_prmi < a_tol_pc {
                the_vtx[j].param_on_line = f64::MAX;
            } else {
                break;
            }
        }
    }
    let a_max_prm = the_vtx[upper - 1].param_on_line;
    for i in (1..upper - 1).rev() {
        let a_prmi = the_vtx[i].param_on_line;
        if a_prmi == f64::MAX {
            continue;
        }
        if (a_max_prm - a_prmi) < a_tol_pc {
            the_vtx[i].param_on_line = f64::MAX;
        } else {
            break;
        }
    }
}
