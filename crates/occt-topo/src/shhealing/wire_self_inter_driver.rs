//! `ShapeAnalysis_Wire` self-intersection driver and the pair overload of
//! `CheckIntersectingEdges`.
//!
//! | OCCT | lines |
//! |---|---|
//! | `CheckSelfIntersection()` | `ShapeAnalysis_Wire.cxx:372-449` |
//! | `CheckClosed(prec)` | `ShapeAnalysis_Wire.cxx:470-499` |
//! | `CheckIntersectingEdges(num1, num2, ...)` | `ShapeAnalysis_Wire.cxx:1557-1695` |
//!
//! The two per-edge checks the driver loops over are
//! [`check_self_intersecting_edge`] / [`check_intersecting_edges`]
//! (`wire_self_inter.rs`, `cxx:1269-1542`); the driver adds the non-adjacent
//! edge-pair loop behind the `Bnd_Box2d` filter (`cxx:404-447`), which needs
//! the `(num1, num2)` overload this file ports.
//!
//! `StatusSelfIntersection(DONE)` (`ShapeAnalysis_Wire.lxx`) is "any DONE bit",
//! because `ShapeExtend::EncodeStatus(ShapeExtend_DONE)` is `0x00ff`
//! (`ShapeExtend.cxx:76-77`); the port therefore returns
//! `done1 || done2 || done3`.
//!
//! Divergence (not a heuristic): the four endpoint vertices of the pair
//! overload are read with `first_vertex` / `last_vertex`, which answer `None`
//! on an edge without vertices (`TopExp::Vertices` on an infinite edge). OCCT
//! dereferences the null `TopoDS_Vertex` there (`cxx:1603-1608`) and crashes;
//! the port stops the check instead (`cxx:1601-1608` has no other arm).

use super::prelude::*;

use occt_core::bnd::BndBox2d;
use occt_core::gp::GpPnt;
use occt_core::intres2d::{IntRes2dDomain, IntRes2dPosition};
use occt_core::precision::{CONFUSION, PCONFUSION};
use occt_geom2d::geom2d_int::Geom2dIntGInter;

use crate::bnd_lib_add2d;
use crate::boptools_2d::curve_on_surface_oriented;
use crate::brep_tool::BRepTool;

use super::wire_fix::{
    check_connected, check_degenerated, first_vertex, last_vertex, ConnectedCheck,
    DegeneratedCheck,
};
use super::wire_self_inter::{
    check_intersecting_edges, check_self_intersecting_edge, point_on_edge, IntersectingEdgesCheck,
    TOLINT,
};

/// `ShapeAnalysis_Wire::CheckIntersectingEdges(num1, num2, points2d, points3d,
/// errors)` (`ShapeAnalysis_Wire.cxx:1557-1695`): the arbitrary edge pair used
/// by the driver's `Bnd_Box2d` loop.
///
/// Unlike the adjacent pair (`(num - 1, num)`, `cxx:1360-1541`) this overload
/// has no junction vertex to unify and no `isLacking` probe: it only drops the
/// intersection points that land on one of the four endpoint vertices
/// (`cxx:1661-1682`).
pub fn check_intersecting_edges_pair(
    wire: &Wire,
    face: &Face,
    num1: usize,
    num2: usize,
) -> IntersectingEdgesCheck {
    let mut out = IntersectingEdgesCheck::default();
    let edges = edges_of_wire(wire);
    let nb = edges.len();
    // `cxx:1565-1568`: `!IsReady()`.
    if nb == 0 {
        return out;
    }
    // `cxx:1570-1571`.
    let n2 = if num2 > 0 { num2 } else { nb };
    let n1 = if num1 > 0 { num1 } else { nb };
    if n1 > nb || n2 > nb {
        return out;
    }
    // `cxx:1573-1574`.
    let e1 = edges[n1 - 1].clone();
    let e2 = edges[n2 - 1].clone();

    // `cxx:1579-1588`: the two pcurves with `orient = false`.
    let Some((crv1, a1, b1)) = curve_on_surface_oriented(&e1, face, false) else {
        out.fail3 = true;
        return out;
    };
    let Some((crv2, a2, b2)) = curve_on_surface_oriented(&e2, face, false) else {
        out.fail3 = true;
        return out;
    };
    // `cxx:1590-1594`.
    if (a1 - b1).abs() <= PCONFUSION || (a2 - b2).abs() <= PCONFUSION {
        return out;
    }

    // `cxx:1611-1617`.
    let d1 = IntRes2dDomain::bounded(&crv1.d0(a1), a1, TOLINT, &crv1.d0(b1), b1, TOLINT);
    let d2 = IntRes2dDomain::bounded(&crv2.d0(a2), a2, TOLINT, &crv2.d0(b2), b2, TOLINT);
    let mut inter = Geom2dIntGInter::new();
    inter.perform(crv1.as_ref(), &d1, crv2.as_ref(), &d2, TOLINT, TOLINT);
    // `cxx:1618-1622`.
    if !inter.is_done() {
        return out;
    }

    // `cxx:1601-1608`: the four endpoint vertices and their tolerances.
    let (Some(vf1), Some(vl1), Some(vf2), Some(vl2)) = (
        first_vertex(&e1),
        last_vertex(&e1),
        first_vertex(&e2),
        last_vertex(&e2),
    ) else {
        return out;
    };
    let vpts = [
        BRepTool::vertex_point(&vf1),
        BRepTool::vertex_point(&vl1),
        BRepTool::vertex_point(&vf2),
        BRepTool::vertex_point(&vl2),
    ];
    let vto = [
        BRepTool::vertex_tolerance(&vf1),
        BRepTool::vertex_tolerance(&vl1),
        BRepTool::vertex_tolerance(&vf2),
        BRepTool::vertex_tolerance(&vl2),
    ];

    let Some(surf) = BRepTool::face_surface(face) else {
        return out;
    };

    // `cxx:1625-1690`: points first, then the segments' first or last point.
    let nb_points = inter.nb_points();
    let nb_segments = inter.nb_segments();
    for i in 1..=(nb_points + nb_segments) {
        let ip = if i <= nb_points {
            inter.point(i)
        } else {
            let seg = inter.segment(i - nb_points);
            if !seg.has_first_point() || !seg.has_last_point() {
                continue;
            }
            let first = *seg.first_point();
            let t1 = first.transition_of_first();
            let t2 = first.transition_of_second();
            if t1.position_on_curve() == IntRes2dPosition::Middle
                || t2.position_on_curve() == IntRes2dPosition::Middle
            {
                *seg.last_point()
            } else {
                first
            }
        };
        // `cxx:1650-1656`.
        let tr1 = ip.transition_of_first();
        let tr2 = ip.transition_of_second();
        if tr1.position_on_curve() != IntRes2dPosition::Middle
            && tr2.position_on_curve() != IntRes2dPosition::Middle
        {
            continue;
        }
        // `cxx:1657-1660`.
        let param1 = ip.param_on_first();
        let param2 = ip.param_on_second();
        let pi1 = point_on_edge(&e1, surf.as_ref(), crv1.as_ref(), param1);
        let pi2 = point_on_edge(&e2, surf.as_ref(), crv2.as_ref(), param2);

        // `cxx:1661-1682`: a point that sits on an endpoint of *its own* edge
        // is the junction itself, not a self-intersection.
        let ok1 = (0..2).any(|j| {
            let d = pi1.square_distance(&vpts[j]);
            d < vto[j] * vto[j]
        });
        let ok2 = (2..4).any(|j| {
            let d = pi2.square_distance(&vpts[j]);
            d < vto[j] * vto[j]
        });
        if !ok1 || !ok2 {
            // `cxx:1684-1689`.
            let pint = GpPnt::new(
                0.5 * (pi1.x() + pi2.x()),
                0.5 * (pi1.y() + pi2.y()),
                0.5 * (pi1.z() + pi2.z()),
            );
            out.points2d.push(ip);
            out.points3d.push(pint);
            out.errors.push(0.5 * pi1.distance(&pi2));
            out.done1 = true;
        }
    }
    out
}

/// `ShapeAnalysis_Wire::CheckClosed(prec)` (`ShapeAnalysis_Wire.cxx:470-499`):
/// the wire's last edge is connected (or degenerated) to its first one.
///
/// `StatusClosed(DONE)` accepts `DONE1` (from `CheckConnected`) and `DONE2`
/// (from `CheckDegenerated`) alike, and both checks run regardless of the
/// first one's answer (`cxx:478-495`).
pub fn check_closed(wire: &Wire, face: &Face, prec: f64) -> bool {
    // `cxx:473-476`.
    if edges_of_wire(wire).is_empty() {
        return false;
    }
    // `cxx:478-486`.
    let connected = !matches!(check_connected(wire, prec, prec, 1), ConnectedCheck::None);
    // `cxx:488-496`.
    let degenerated = matches!(
        check_degenerated(wire, face, prec, 1),
        DegeneratedCheck::Found { .. }
    );
    connected || degenerated
}

/// `ShapeAnalysis_Wire::CheckSelfIntersection()` (`ShapeAnalysis_Wire.cxx:372-449`).
///
/// The analyzer is constructed with `Precision::Confusion()` wherever the port
/// has a caller (`ChFi2d_AnaFilletAlgo::Perform`,
/// `ChFi2d_AnaFilletAlgo.cxx:358`), so both `myPrecision` and the `prec`
/// argument of the two per-edge checks are that same value.
pub fn check_self_intersection(wire: &Wire, face: &Face) -> bool {
    let mut done1 = false;
    let mut done2 = false;
    let mut done3 = false;
    let nb = edges_of_wire(wire).len();
    // `cxx:376-379`: `!IsReady()`.
    if nb == 0 {
        return false;
    }
    // `cxx:380-402`: the per-edge and adjacent-pair loops.
    for i in 1..=nb {
        if check_self_intersecting_edge(wire, face, i).is_done() {
            done1 = true;
        }
        if check_intersecting_edges(wire, face, i, CONFUSION).is_done() {
            done2 = true;
        }
    }

    // `cxx:404-420`: one `Bnd_Box2d` per edge that has a pcurve; a missing
    // pcurve leaves the box Void, and `Bnd_Box2d::IsOut` reports a Void box as
    // out of everything (`Bnd_Box2d.cxx:232-242`).
    let edges = edges_of_wire(wire);
    let mut boxes: Vec<BndBox2d> = vec![BndBox2d::default(); nb];
    for (i, e) in edges.iter().enumerate() {
        if let Some((c2d, cf, cl)) = curve_on_surface_oriented(e, face, false) {
            bnd_lib_add2d::add_geom2d_range(c2d.as_ref(), cf, cl, CONFUSION, &mut boxes[i]);
        }
    }

    // `cxx:422-441`.
    let closed = check_closed(wire, face, CONFUSION);
    let mut is_fail = false;
    for num1 in 1..nb.saturating_sub(1) {
        // `cxx:427-431`: a closed wire skips the (first, last) pair.
        let fin = if closed && num1 == 1 { nb - 1 } else { nb };
        for num2 in (num1 + 2)..=fin {
            if !boxes[num1 - 1].is_out_box(&boxes[num2 - 1]) {
                let r = check_intersecting_edges_pair(wire, face, num1, num2);
                // `cxx:438-439`: `LastCheckStatus(FAIL1)` / `(DONE1)`.
                is_fail |= r.fail1;
                done3 |= r.done1;
            }
        }
    }
    // `cxx:442-447` marks FAIL3 / DONE3; FAIL3 has no consumer here, because
    // the driver only reports `StatusSelfIntersection(DONE)`.
    let _ = is_fail;

    // `cxx:449`: `StatusSelfIntersection(ShapeExtend_DONE)`.
    done1 || done2 || done3
}
