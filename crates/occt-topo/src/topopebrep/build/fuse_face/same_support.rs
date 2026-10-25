//! `SameSupport` (`TopOpeBRepBuild_FuseFace.cxx:1181-1398`) — do two edges lie
//! on the same geometric support?
//!
//! The test is the OCCT one, type by type, including its asymmetries:
//! `BRep_Tool::Curve`'s location is folded into the curve (`cxx:1196-1208`),
//! one `Geom_TrimmedCurve` level is stripped (`cxx:1212-1224`), the two
//! `DynamicType`s must match (`cxx:1226-1229`), and only `Geom_Line`,
//! `Geom_Circle`, `Geom_Ellipse`, `Geom_BSplineCurve` and `Geom_BezierCurve`
//! are compared at all (`cxx:1230-1238`); any other support answers false even
//! when both edges have it.
//!
//! The knot comparison is one-sided in OCCT — `(K1(k) - K2(k)) > tollin` at
//! `cxx:1306`, with no `abs` — and is kept that way; the multiplicity test at
//! `cxx:1309` compares two `int`s against `tollin`, which makes it an
//! exact-equality test, also kept.

use std::sync::Arc;

use occt_core::precision::{ANGULAR, CONFUSION};
use occt_geom::Curve;

use crate::brep_tool::BRepTool;
use crate::shape::{Edge, TopoShape};

/// `Geom_Curve::DynamicType()` grouped into the five types `SameSupport`
/// compares (`cxx:1230-1237`); [`CurveKind::Other`] is every other support.
#[derive(PartialEq, Eq, Clone, Copy)]
enum CurveKind {
    Line,
    Circle,
    Ellipse,
    BSpline,
    Bezier,
    Other,
}

fn curve_kind(c: &dyn Curve) -> CurveKind {
    if c.gp_line().is_some() || c.is_line() {
        CurveKind::Line
    } else if c.gp_circ().is_some() {
        CurveKind::Circle
    } else if c.gp_ellipse().is_some() {
        CurveKind::Ellipse
    } else if c.bspline_poles().is_some() {
        CurveKind::BSpline
    } else if c.bezier_poles().is_some() {
        CurveKind::Bezier
    } else {
        CurveKind::Other
    }
}

/// `cxx:1209-1224`: `if (typ == Geom_TrimmedCurve) C = BasisCurve();` — one
/// level, applied after the location was folded in.
fn untrim(c: Arc<dyn Curve>) -> Arc<dyn Curve> {
    if c.is_geom_trimmed() {
        c.trim_basis().unwrap_or(c)
    } else {
        c
    }
}

/// Distinct knots and their multiplicities from the flat knot sequence
/// (`Geom_BSplineCurve::Knots()` + `Multiplicities()`).
fn split_knots(flat: &[f64]) -> (Vec<f64>, Vec<i32>) {
    let mut knots: Vec<f64> = Vec::new();
    let mut mults: Vec<i32> = Vec::new();
    for k in flat {
        match knots.last_mut() {
            Some(last) if *last == *k => *mults.last_mut().expect("parallel arrays") += 1,
            _ => {
                knots.push(*k);
                mults.push(1);
            }
        }
    }
    (knots, mults)
}

/// `SameSupport(E1, E2)` (`cxx:1181-1398`).
pub fn same_support(e1: &TopoShape, e2: &TopoShape) -> bool {
    // `cxx:1188-1191`: null edges never share a support.
    let (Some(edge1), Some(edge2)) = (Edge::wrap(e1.clone()), Edge::wrap(e2.clone())) else {
        return false;
    };

    // `cxx:1196-1208`: `BRep_Tool::Curve(E, loc, f, l)` with the location
    // folded into the curve when it is not the identity.
    let (Some(c1), Some(c2)) = (
        BRepTool::edge_curve_world(&edge1),
        BRepTool::edge_curve_world(&edge2),
    ) else {
        return false;
    };
    let c1 = untrim(c1);
    let c2 = untrim(c2);

    // `cxx:1210-1229`.
    let typ1 = curve_kind(c1.as_ref());
    let typ2 = curve_kind(c2.as_ref());
    if typ1 != typ2 {
        return false;
    }

    // `cxx:1230-1238`: any other support is "not handled" -> false.
    if typ1 == CurveKind::Other {
        return false;
    }

    // `cxx:1241-1243`.
    let tollin = CONFUSION;
    let tolang = ANGULAR;

    match typ1 {
        // `cxx:1244-1252`.
        CurveKind::Line => {
            let (Some(li1), Some(li2)) = (c1.gp_line(), c2.gp_line()) else {
                return false;
            };
            li1.angle(&li2).abs() <= tolang
                && li1.location().square_distance(&li2.location()) <= tollin * tollin
        }
        // `cxx:1253-1258`.
        CurveKind::Circle => {
            let (Some(ci1), Some(ci2)) = (c1.gp_circ(), c2.gp_circ()) else {
                return false;
            };
            (ci1.radius() - ci2.radius()).abs() <= tollin
                && ci1.location().square_distance(&ci2.location()) <= tollin * tollin
        }
        // `cxx:1259-1268`.
        CurveKind::Ellipse => {
            let (Some(ci1), Some(ci2)) = (c1.gp_ellipse(), c2.gp_ellipse()) else {
                return false;
            };
            (ci1.major_radius() - ci2.major_radius()).abs() <= tollin
                && (ci1.minor_radius() - ci2.minor_radius()).abs() <= tollin
                && ci1.location().square_distance(&ci2.location()) <= tollin * tollin
        }
        // `cxx:1269-1344`.
        CurveKind::BSpline => {
            let (Some(p1), Some(p2)) = (c1.bspline_poles(), c2.bspline_poles()) else {
                return false;
            };
            if p1.len() != p2.len() {
                return false;
            }
            let nbpoles = p1.len();
            let (k1, m1) = split_knots(c1.bspline_knots().unwrap_or(&[]));
            let (k2, m2) = split_knots(c2.bspline_knots().unwrap_or(&[]));
            if k1.len() != k2.len() {
                return false;
            }
            // `cxx:1288-1295`: `tol3d = BRep_Tool::Tolerance(E1)`.
            let tol3d = BRepTool::edge_tolerance(&edge1);
            for p in 0..nbpoles {
                if p1[p].distance(&p2[p]) > tol3d {
                    return false;
                }
            }
            // `cxx:1297-1313`: one-sided knot test, exact multiplicity test.
            for k in 0..k1.len() {
                if (k1[k] - k2[k]) > tollin {
                    return false;
                }
                if (m1[k] as f64 - m2[k] as f64).abs() > tollin {
                    return false;
                }
            }
            // `cxx:1315-1328`.
            let r1 = c1.bspline_weights().is_some();
            let r2 = c2.bspline_weights().is_some();
            if r1 != r2 {
                return false;
            }
            // `cxx:1330-1343`.
            if r1 {
                let (Some(w1), Some(w2)) = (c1.bspline_weights(), c2.bspline_weights()) else {
                    return false;
                };
                for w in 0..nbpoles {
                    if (w1[w] - w2[w]).abs() > tollin {
                        return false;
                    }
                }
            }
            true
        }
        // `cxx:1345-1396`.
        CurveKind::Bezier => {
            let (Some(p1), Some(p2)) = (c1.bezier_poles(), c2.bezier_poles()) else {
                return false;
            };
            if p1.len() != p2.len() {
                return false;
            }
            let nbpoles = p1.len();
            // `cxx:1359-1364`: the Bezier branch compares with `tollin`.
            for p in 0..nbpoles {
                if p1[p].distance(&p2[p]) > tollin {
                    return false;
                }
            }
            // `cxx:1366-1381`.
            let r1 = c1.bezier_weights().is_some();
            let r2 = c2.bezier_weights().is_some();
            if r1 != r2 {
                return false;
            }
            // `cxx:1383-1394`.
            if r1 {
                let (Some(w1), Some(w2)) = (c1.bezier_weights(), c2.bezier_weights()) else {
                    return false;
                };
                for w in 0..nbpoles {
                    if (w1[w] - w2[w]).abs() > tollin {
                        return false;
                    }
                }
            }
            true
        }
        CurveKind::Other => false,
    }
}
