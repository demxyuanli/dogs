//! Convert `IntAna_QuadQuadGeo` results into `IntPatch_GLine`s.
//! Source: IntPP / IntPCy / IntPSp / IntPCo GLine construction.

use occt_core::gp::{GpCirc, GpPnt};
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
