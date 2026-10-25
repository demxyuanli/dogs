//! `TopOpeBRepTool_FuseEdges::SameSupport` (`TopOpeBRepTool_FuseEdges.cxx:624-880`).
//!
//! Two edges are fusible only when their support curves are the same
//! geometric curve, tested per curve type:
//!
//! | type | test | OCCT |
//! |---|---|---|
//! | `Geom_Line` | directions parallel, plus a non-degeneracy guard on the end points | `cxx:694-713` |
//! | `Geom_Circle` | equal radius, coincident centre, parallel axis | `cxx:715-722` |
//! | `Geom_Ellipse` | equal major/minor radii, coincident centre, parallel axis | `cxx:723-733` |
//! | `Geom_BSplineCurve` | contiguous ranges, equal poles / knots / multiplicities / rationality / weights | `cxx:734-818` |
//! | `Geom_BezierCurve` | contiguous ranges, equal poles / rationality / weights | `cxx:819-877` |
//!
//! Any other type answers `false` (`cxx:684-691`).
//!
//! Fidelity notes:
//! - the `Knots()` loop compares `K1(k) - K2(k) > tollin` **without** `abs`
//!   (`cxx:779-782`): the port keeps that asymmetry;
//! - the BSpline pole loop uses `BRep_Tool::Tolerance(E1)` as the tolerance
//!   (`cxx:762`) while the Bezier one uses `Precision::Confusion()`
//!   (`cxx:845-848`);
//! - `Geom_TrimmedCurve` is unwrapped exactly **one** level (`cxx:667-676`), not
//!   recursively.

use std::sync::Arc;

use occt_core::gp::{GpDir, GpPnt};
use occt_core::precision;
use occt_geom::Curve;

use crate::brep_tool::BRepTool;
use crate::shape::{Edge, TopoShape};

/// The five support types `SameSupport` accepts (`cxx:683-691`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SupportKind {
    Line,
    Circle,
    Ellipse,
    BSpline,
    Bezier,
}

/// `C1->DynamicType()` after the `Geom_TrimmedCurve` unwrap (`cxx:667-691`).
fn support_kind(c: &dyn Curve) -> Option<SupportKind> {
    if c.gp_line().is_some() {
        Some(SupportKind::Line)
    } else if c.gp_circ().is_some() {
        Some(SupportKind::Circle)
    } else if c.gp_ellipse().is_some() {
        Some(SupportKind::Ellipse)
    } else if c.bspline_poles().is_some() {
        Some(SupportKind::BSpline)
    } else if c.bezier_poles().is_some() {
        Some(SupportKind::Bezier)
    } else {
        None
    }
}

/// `C = BRep_Tool::Curve(E, loc, f, l)`, transformed by `loc` (`cxx:636-661`)
/// and unwrapped one `Geom_TrimmedCurve` level (`cxx:667-676`).
///
/// `None` is the `C1.IsNull()` / `C2.IsNull()` early answer (`cxx:639-642`,
/// `cxx:652-655`): a degenerated edge carries no 3D curve.
fn support_curve(e: &Edge) -> Option<(Arc<dyn Curve>, f64, f64)> {
    let c = BRepTool::edge_curve(e)?;
    let (f, l) = BRepTool::edge_parameters(e);
    let loc = e.0.location();
    let mut c: Arc<dyn Curve> = if loc.is_identity() {
        c
    } else {
        Arc::from(c.transformed(&loc.transformation()))
    };
    if let Some((basis, _, _)) = c.untrimmed_basis() {
        c = basis;
    }
    Some((c, f, l))
}

/// `BRep_Tool::Pnt(TopExp::FirstVertex(E, true))` / `LastVertex` (`cxx:707-710`).
fn vertex_point(e: &TopoShape, last: bool) -> Option<GpPnt> {
    let v = super::cumulated_vertex(e, last)?;
    let v = crate::shape::Vertex::wrap(v)?;
    Some(BRepTool::vertex_point_world(&v))
}

/// Distinct knots and their multiplicities from a flat knot array
/// (`Geom_BSplineCurve::Knots()` / `Multiplicities()`, `cxx:772-776`).
fn distinct_knots(flat: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let mut knots: Vec<f64> = Vec::new();
    let mut mults: Vec<f64> = Vec::new();
    for &k in flat {
        if let Some(last) = knots.last() {
            if (k - *last).abs() <= f64::EPSILON {
                *mults.last_mut().unwrap() += 1.0;
                continue;
            }
        }
        knots.push(k);
        mults.push(1.0);
    }
    (knots, mults)
}

/// `SameSupport(E1, E2)` (`cxx:624-880`).
pub(super) fn same_support(e1: &Edge, e2: &Edge) -> bool {
    let Some((c1, f1, l1)) = support_curve(e1) else {
        return false;
    };
    let Some((c2, f2, l2)) = support_curve(e2) else {
        return false;
    };
    let (k1, k2) = (support_kind(c1.as_ref()), support_kind(c2.as_ref()));
    // `cxx:678-681`: `typC1 != typC2` -> false.
    if k1.is_none() || k1 != k2 {
        return false;
    }

    // `cxx:693-695`: `tollin = Precision::Confusion()`, `tolang = Precision::Angular()`.
    let tollin = precision::CONFUSION;
    let tolang = precision::ANGULAR;

    match k1.unwrap() {
        // `cxx:696-714`
        SupportKind::Line => {
            let li1 = c1.gp_line().unwrap();
            let li2 = c2.gp_line().unwrap();
            let dir1: GpDir = li1.pos.direction().clone();
            let dir2: GpDir = li2.pos.direction().clone();
            if dir1.is_parallel_tol(&dir2, tolang) {
                // `cxx:706-711`: reject two connex edges that loop onto each
                // other around one point.
                let (pf1, pl1) = (vertex_point(&e1.0, false), vertex_point(&e1.0, true));
                let (pf2, pl2) = (vertex_point(&e2.0, false), vertex_point(&e2.0, true));
                return match (pl1, pf2, pl2, pf1) {
                    (Some(pl1), Some(pf2), Some(pl2), Some(pf1)) => {
                        pl1.distance(&pf2) >= tollin || pl2.distance(&pf1) >= tollin
                    }
                    // A null vertex reaches `BRep_Tool::Pnt` and OCCT raises
                    // `Standard_NullObject`; the port answers "not fusible".
                    _ => false,
                };
            }
            false
        }
        // `cxx:715-722`
        SupportKind::Circle => {
            let ci1 = c1.gp_circ().unwrap();
            let ci2 = c2.gp_circ().unwrap();
            (ci1.radius - ci2.radius).abs() <= tollin
                && ci1.location().square_distance(&ci2.location()) <= tollin * tollin
                && ci1.axis().direction().is_parallel_tol(ci2.axis().direction(), tolang)
        }
        // `cxx:723-733`
        SupportKind::Ellipse => {
            let ci1 = c1.gp_ellipse().unwrap();
            let ci2 = c2.gp_ellipse().unwrap();
            (ci1.major_radius - ci2.major_radius).abs() <= tollin
                && (ci1.minor_radius - ci2.minor_radius).abs() <= tollin
                && ci1.location().square_distance(&ci2.location()) <= tollin * tollin
                && ci1.axis().direction().is_parallel_tol(ci2.axis().direction(), tolang)
        }
        // `cxx:734-818`
        SupportKind::BSpline => {
            // `cxx:737-741`: the end of one curve must not be the beginning of
            // the other.
            if (f2 - l1).abs() > tollin && (f1 - l2).abs() > tollin {
                return false;
            }
            let p1 = c1.bspline_poles().unwrap();
            let p2 = c2.bspline_poles().unwrap();
            if p1.len() != p2.len() {
                return false;
            }
            let (kn1, kn2) = (c1.bspline_knots().unwrap(), c2.bspline_knots().unwrap());
            let (d1, m1) = distinct_knots(kn1);
            let (d2, m2) = distinct_knots(kn2);
            if d1.len() != d2.len() {
                return false;
            }
            // `cxx:762`: the pole tolerance is the *edge* tolerance of `E1`.
            let tol3d = BRepTool::edge_tolerance(e1);
            for (a, b) in p1.iter().zip(p2.iter()) {
                if a.distance(b) > tol3d {
                    return false;
                }
            }
            for i in 0..d1.len() {
                // `cxx:779-782`: `(K1(k) - K2(k)) > tollin`, without `abs`.
                if d1[i] - d2[i] > tollin {
                    return false;
                }
                // `cxx:783-786`
                if (m1[i] - m2[i]).abs() > tollin {
                    return false;
                }
            }
            let r1 = c1.bspline_weights();
            let r2 = c2.bspline_weights();
            // `cxx:789-802`
            if r1.is_some() != r2.is_some() {
                return false;
            }
            // `cxx:804-816`
            if let (Some(w1), Some(w2)) = (r1, r2) {
                for (a, b) in w1.iter().zip(w2.iter()) {
                    if (a - b).abs() > tollin {
                        return false;
                    }
                }
            }
            true
        }
        // `cxx:819-877`
        SupportKind::Bezier => {
            // `cxx:825-828`
            if (f2 - l1).abs() > tollin && (f1 - l2).abs() > tollin {
                return false;
            }
            let p1 = c1.bezier_poles().unwrap();
            let p2 = c2.bezier_poles().unwrap();
            if p1.len() != p2.len() {
                return false;
            }
            // `cxx:842-848`: `Precision::Confusion()` here, not the edge tolerance.
            for (a, b) in p1.iter().zip(p2.iter()) {
                if a.distance(b) > tollin {
                    return false;
                }
            }
            let r1 = c1.bezier_weights();
            let r2 = c2.bezier_weights();
            // `cxx:850-863`
            if r1.is_some() != r2.is_some() {
                return false;
            }
            // `cxx:865-874`
            if let (Some(w1), Some(w2)) = (r1, r2) {
                for (a, b) in w1.iter().zip(w2.iter()) {
                    if (a - b).abs() > tollin {
                        return false;
                    }
                }
            }
            true
        }
    }
}
