//! `GeomInt_LineTool`. Source: `GeomInt_LineTool.cxx`.

use occt_core::precision::Precision;

use crate::geom_int::types::{GeomIntLine, IntPatchIType};
use crate::geom_int::LineConstructor;
use crate::int_tools_wline::{PatchPoint, WLine};
use occt_geom::Surface;
use crate::shape::Face;

/// `GeomInt_LineTool::NbVertex`.
pub fn nb_vertex(l: &GeomIntLine) -> i32 {
    match l {
        GeomIntLine::Analytic(a) => a.nb_vertex(),
        GeomIntLine::Restriction(r) => r.nb_vertex(),
        GeomIntLine::Walking(w) => {
            if w.nb_vertex() > 0 {
                w.nb_vertex()
            } else if w.nb_pnts() <= 1 {
                w.nb_pnts()
            } else {
                2
            }
        }
        GeomIntLine::Geometric(g) => g.nb_vertex(),
    }
}

/// `GeomInt_LineTool::Vertex`. Walking lines without stored vertices
/// synthesize the first/last sample (`HasFirstPoint` false path still
/// reports parameters 1 and `NbPnts` as ends).
pub fn vertex(l: &GeomIntLine, i: i32) -> PatchPoint {
    match l {
        GeomIntLine::Analytic(a) => *a.vertex(i),
        GeomIntLine::Restriction(r) => *r.vertex(i),
        GeomIntLine::Walking(w) => {
            if w.nb_vertex() > 0 {
                *w.vertex(i)
            } else {
                let n = w.nb_pnts();
                let idx = if i <= 1 { 1 } else { n };
                let p = w.point(idx);
                PatchPoint::new(p.value(), idx as f64, p.u1, p.v1, p.u2, p.v2)
            }
        }
        GeomIntLine::Geometric(g) => *g.vertex(i),
    }
}

/// `GeomInt_LineTool::FirstParameter`.
pub fn first_parameter(l: &GeomIntLine) -> f64 {
    match l.arc_type() {
        IntPatchIType::Analytic => {
            let a = l.as_aline().expect("analytic");
            if a.has_first_point {
                a.vertex(1).parameter_on_line()
            } else {
                let mut firstp = a.curve.first_parameter();
                firstp += firstp.abs() * f64::EPSILON;
                firstp
            }
        }
        IntPatchIType::Restriction => {
            let r = l.as_rline().expect("rline");
            if r.has_first_point {
                r.vertex(1).parameter_on_line()
            } else {
                -Precision::INFINITE
            }
        }
        IntPatchIType::Walking => {
            let w = l.as_wline().expect("wline");
            if w.has_first_point && !w.vertices.is_empty() {
                w.vertex(1).parameter_on_line()
            } else {
                1.0
            }
        }
        IntPatchIType::Lin | IntPatchIType::Parabola | IntPatchIType::Hyperbola => {
            let g = l.as_gline().expect("gline");
            if g.has_first_point && !g.vertices.is_empty() {
                g.vertex(1).parameter_on_line()
            } else {
                -Precision::INFINITE
            }
        }
        IntPatchIType::Circle | IntPatchIType::Ellipse => {
            let g = l.as_gline().expect("gline");
            if g.has_first_point && !g.vertices.is_empty() {
                g.vertex(1).parameter_on_line()
            } else {
                0.0
            }
        }
    }
}

/// `GeomInt_LineTool::LastParameter`.
pub fn last_parameter(l: &GeomIntLine) -> f64 {
    match l.arc_type() {
        IntPatchIType::Analytic => {
            let a = l.as_aline().expect("analytic");
            if a.has_last_point && a.nb_vertex() > 0 {
                a.vertex(a.nb_vertex()).parameter_on_line()
            } else {
                let mut lastp = a.curve.last_parameter();
                lastp -= lastp.abs() * f64::EPSILON;
                lastp
            }
        }
        IntPatchIType::Restriction => {
            let r = l.as_rline().expect("rline");
            if r.has_last_point && r.nb_vertex() > 0 {
                r.vertex(r.nb_vertex()).parameter_on_line()
            } else {
                Precision::INFINITE
            }
        }
        IntPatchIType::Walking => {
            let w = l.as_wline().expect("wline");
            if w.has_last_point && w.nb_vertex() > 0 {
                w.vertex(w.nb_vertex()).parameter_on_line()
            } else {
                w.nb_pnts() as f64
            }
        }
        IntPatchIType::Lin | IntPatchIType::Parabola | IntPatchIType::Hyperbola => {
            let g = l.as_gline().expect("gline");
            if g.has_last_point && g.nb_vertex() > 0 {
                g.vertex(g.nb_vertex()).parameter_on_line()
            } else {
                Precision::INFINITE
            }
        }
        IntPatchIType::Circle | IntPatchIType::Ellipse => {
            let g = l.as_gline().expect("gline");
            if g.has_last_point && g.nb_vertex() > 0 {
                g.vertex(g.nb_vertex()).parameter_on_line()
            } else {
                2.0 * std::f64::consts::PI
            }
        }
    }
}

/// `GeomInt_LineTool::DecompositionOfWLine` — constructor parts become the
/// `line_parts` argument of `IntTools_WLineTool::DecompositionOfWLine`.
pub fn decomposition_of_wline(
    the_wline: &WLine,
    the_surface1: &dyn Surface,
    the_surface2: &dyn Surface,
    a_tol_sum: f64,
    the_l_constructor: &LineConstructor,
    the_new_lines: &mut Vec<WLine>,
) -> bool {
    let a_nb_parts = the_l_constructor.nb_parts();
    let a_nb_pnts = the_wline.nb_pnts();
    if a_nb_pnts == 0 || a_nb_parts <= 0 {
        return false;
    }
    let mut line_parts = Vec::with_capacity(a_nb_parts as usize);
    for i in 1..=a_nb_parts {
        if let Some((f, l)) = the_l_constructor.part(i) {
            line_parts.push((f as i32, l as i32));
        }
    }
    let dummy = Face::new();
    crate::int_tools_wline::decomposition_of_wline(
        the_wline,
        the_surface1,
        the_surface2,
        &dummy,
        &dummy,
        &line_parts,
        false,
        a_tol_sum,
        the_new_lines,
    )
}
