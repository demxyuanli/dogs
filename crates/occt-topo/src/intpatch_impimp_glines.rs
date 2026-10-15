//! Convert `IntAna_QuadQuadGeo` results into `IntPatch_GLine`s.
//! Source: IntPP / IntPCy / IntPSp / IntPCo GLine construction.

use occt_core::elib::clib;
use occt_core::gp::{GpCirc, GpLin, GpPnt};
use occt_geom::intana::QuadricIntersection;

use crate::geom_int::{GLine, GLineKind, GeomIntLine};
use crate::int_tools_wline::PatchPoint;

use super::quad::{adjust_circ_to_seam, adjust_sphere_circ, ImplicitQuad};

pub(crate) enum PairOutcome {
    Fail,
    Empty,
    Same,
    /// `IntStatus_InfiniteSectionCurve` (CyCyNoGeometric V-extent cap).
    Infinite,
    Result {
        lines: Vec<GeomIntLine>,
        points: Vec<PatchPoint>,
    },
}

pub(crate) fn isolated_point(p: GpPnt) -> PatchPoint {
    PatchPoint::new(p, 0.0, 0.0, 0.0, 0.0, 0.0)
}

fn push_gline(lines: &mut Vec<GeomIntLine>, kind: GLineKind) {
    lines.push(GeomIntLine::Geometric(GLine::new(kind)));
}

fn seam_adjust_circle(circ: &mut GpCirc, seam: Option<&ImplicitQuad>) {
    let Some(q) = seam else {
        return;
    };
    match q {
        ImplicitQuad::Cylinder(cy) => adjust_circ_to_seam(circ, &cy.position()),
        ImplicitQuad::Cone(co) => adjust_circ_to_seam(circ, &co.position()),
        ImplicitQuad::Sphere(sp) => adjust_sphere_circ(circ, sp),
        ImplicitQuad::Torus(to) => adjust_circ_to_seam(circ, to.position()),
        ImplicitQuad::Plane(_) => {}
    }
}

/// Map a closed-form quadric intersection onto GLines / isolated points.
/// `seam` is the surface of revolution used by `AdjustToSeam` (cylinder,
/// sphere, cone, or torus). Circles are left as-is when `seam` is None.
pub(crate) fn from_quadric_intersection(
    qi: QuadricIntersection,
    seam: Option<&ImplicitQuad>,
) -> PairOutcome {
    match qi {
        QuadricIntersection::Same => PairOutcome::Same,
        QuadricIntersection::None => PairOutcome::Empty,
        QuadricIntersection::Point(p) => PairOutcome::Result {
            lines: Vec::new(),
            points: vec![isolated_point(p)],
        },
        QuadricIntersection::Line(l) => {
            let mut lines = Vec::new();
            push_gline(&mut lines, GLineKind::Lin(l));
            PairOutcome::Result {
                lines,
                points: Vec::new(),
            }
        }
        QuadricIntersection::TwoLines(a, b) => {
            let mut lines = Vec::new();
            push_gline(&mut lines, GLineKind::Lin(a));
            push_gline(&mut lines, GLineKind::Lin(b));
            PairOutcome::Result {
                lines,
                points: Vec::new(),
            }
        }
        QuadricIntersection::Circle(mut c) => {
            seam_adjust_circle(&mut c, seam);
            let mut lines = Vec::new();
            push_gline(&mut lines, GLineKind::Circ(c));
            PairOutcome::Result {
                lines,
                points: Vec::new(),
            }
        }
        QuadricIntersection::TwoCircles(mut a, mut b) => {
            seam_adjust_circle(&mut a, seam);
            seam_adjust_circle(&mut b, seam);
            let mut lines = Vec::new();
            push_gline(&mut lines, GLineKind::Circ(a));
            push_gline(&mut lines, GLineKind::Circ(b));
            PairOutcome::Result {
                lines,
                points: Vec::new(),
            }
        }
        QuadricIntersection::Ellipse(e) => {
            let mut lines = Vec::new();
            push_gline(&mut lines, GLineKind::Elips(e));
            PairOutcome::Result {
                lines,
                points: Vec::new(),
            }
        }
        QuadricIntersection::TwoEllipses(a, b) => {
            let mut lines = Vec::new();
            push_gline(&mut lines, GLineKind::Elips(a));
            push_gline(&mut lines, GLineKind::Elips(b));
            PairOutcome::Result {
                lines,
                points: Vec::new(),
            }
        }
        QuadricIntersection::Parabola(p) => {
            let mut lines = Vec::new();
            push_gline(&mut lines, GLineKind::Parab(p));
            PairOutcome::Result {
                lines,
                points: Vec::new(),
            }
        }
        QuadricIntersection::Hyperbola(h) => {
            let mut lines = Vec::new();
            push_gline(&mut lines, GLineKind::Hypr(h));
            PairOutcome::Result {
                lines,
                points: Vec::new(),
            }
        }
    }
}

/// Circles recovered by `plane_torus_circles` (IntPTo).
pub(crate) fn from_circles(circs: Vec<GpCirc>, seam: Option<&ImplicitQuad>) -> PairOutcome {
    if circs.is_empty() {
        return PairOutcome::Empty;
    }
    let mut lines = Vec::new();
    for mut c in circs {
        seam_adjust_circle(&mut c, seam);
        push_gline(&mut lines, GLineKind::Circ(c));
    }
    PairOutcome::Result {
        lines,
        points: Vec::new(),
    }
}

/// `IntCoCo`, `IntPatch_ImpImpIntersection.cxx:8869-8941` / `:8945-8981` /
/// `:8986-9017` / `:9149-9171`: when the two cones share the generatrix through
/// their apexes, that generatrix is part of the solution and the characteristic
/// point becomes a vertex on every conic of the arm
/// (`aPtsol.SetParameter(0.)`), while the generatrix line carries it at
/// `ElCLib::Parameter(linsol, aPChar)`.
pub(crate) fn with_cone_common_gen(
    out: PairOutcome,
    gen: GpLin,
    p_char: GpPnt,
    q1: &ImplicitQuad,
    q2: &ImplicitQuad,
    tol: f64,
) -> PairOutcome {
    match out {
        PairOutcome::Result { mut lines, points } => {
            let (u1, v1) = q1.parameters(&p_char);
            let (u2, v2) = q2.parameters(&p_char);
            // `SetValue(aPChar, Tol, false)` - tolerance set, not a mult-point.
            let vertex = |param: f64| {
                let mut vtx = PatchPoint::new(p_char, param, u1, v1, u2, v2);
                vtx.tolerance = tol;
                vtx
            };
            for line in &mut lines {
                if let GeomIntLine::Geometric(g) = line {
                    g.add_vertex(vertex(0.0));
                }
            }
            let mut glig = GLine::new(GLineKind::Lin(gen));
            glig.add_vertex(vertex(clib::line_parameter(&gen, &p_char)));
            lines.push(GeomIntLine::Geometric(glig));
            PairOutcome::Result { lines, points }
        }
        other => other,
    }
}
