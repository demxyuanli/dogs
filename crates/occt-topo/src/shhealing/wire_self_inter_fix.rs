//! `ShapeFix_Wire` self-intersection actions: `TryNewPCurve`, `RemoveLoop`,
//! `FixSelfIntersectingEdge` and the `FixSelfIntersection` driver
//! (`ShapeFix_Wire.cxx:1083-2913`).
//!
//! The 2D checks these actions run on live in [`super::wire_self_inter`]
//! (`ShapeAnalysis_Wire::CheckSelfIntersectingEdge` /
//! `CheckIntersectingEdges`), together with the adjacent-edge action
//! `FixIntersectingEdges`. This file holds the single-edge loop removal
//! (`RemoveLoop`, both its internal variants) and the two drivers.

use super::prelude::*;

use occt_core::convert::ParameterisationType;
use occt_core::gp::{GpAx3, GpPln};
use occt_core::intres2d::IntRes2dIntersectionPoint;
use occt_core::precision::{CONFUSION, PCONFUSION};
use occt_geom::convert_bspl::CompCurveToBSplineCurve;
use occt_geom::geom_api::to_2d;
use occt_geom::geom_lib::to_3d;
use occt_geom::trimmed::GeomTrimmedCurveBasis;
use occt_geom::Curve;
use occt_geom2d::curve::Curve2d;

use crate::boptools_2d::curve_on_surface_oriented;
use crate::brep_tool::BRepTool;
use crate::shape_fix_compose_shell::IdentityReShape;
use crate::shhealing::intersection_tool::{
    fix_self_intersect_wire, SHAPE_FIX_INTERSECTION_MAX_TOL,
};
use crate::shhealing::transfer_params::TransferParametersProj;
use crate::shhealing::wire_fix::{
    first_vertex, fix_same_parameter, last_vertex, wire_remove_edge, SHAPE_FIX_MAX_TOLERANCE,
};
use crate::shhealing::wire_self_inter::{check_self_intersecting_edge, fix_intersecting_edges};

/// `BRep_Builder::UpdateVertex(V, Tol)` (`BRep_Builder.cxx:1442-1452`):
/// `BRep_TVertex::UpdateTolerance` is a max (`BRep_TVertex.lxx:33-37`).
fn update_vertex_tolerance(v: &Vertex, tol: f64) {
    v.set_tolerance(BRepTool::vertex_tolerance(v).max(tol));
}

/// `TryNewPCurve(E, face, c2d, first, last, tol)` (`ShapeFix_Wire.cxx:2215-2251`).
///
/// Builds a temporary edge from `E`'s 3D curve over `[f, l]`, attaches the
/// candidate pcurve `c2d` over `[first, last]`, runs
/// `ShapeFix_Edge::FixSameParameter` on it and returns the modified pcurve, its
/// range and the resulting edge tolerance. `None` is the OCCT `return false`.
fn try_new_pcurve(
    e: &Edge,
    face: &Face,
    c2d: &mut Arc<dyn Curve2d>,
    first: &mut f64,
    last: &mut f64,
) -> Option<f64> {
    // `cxx:2224-2228`.
    let crv = BRepTool::edge_curve(e)?;
    let (f, l) = BRepTool::edge_parameters(e);

    // `cxx:2230-2237`: `BRepBuilderAPI_MakeEdge mkedge(crv, f, l)` and
    // `SBE.SetRange3d(mkedge, f, l)`. The temp edge keeps the 3D range `[f, l]`
    // (`make_edge` already sets it); its tolerance is `BRepLib::Precision()`
    // (`BRepLib_MakeEdge.cxx:782`) and it carries the two endpoint vertices
    // (`BRepLib_MakeEdge.cxx:696`), both at `Precision::Confusion()`.
    let builder = TopoBuilder::new();
    let mut edge = builder.make_edge(crv.clone(), f, l);
    let reg = GeometryRegistry::global();
    reg.set_edge_tolerance(&edge.0, CONFUSION);
    let v1 = builder.make_vertex(crv.d0(f), CONFUSION);
    let v2 = builder.make_vertex(crv.d0(l), CONFUSION);
    builder.add_edge_vertices(&mut edge, &v1, &v2);

    let face_key = GeometryRegistry::shape_key(&face.0);
    // `cxx:2242-2245`.
    reg.set_edge_pcurve(&edge.0, face_key, c2d.clone());
    reg.set_pcurve_range(&edge.0, face_key, *first, *last);
    reg.set_same_range(&edge.0, false);

    // `cxx:2247-2250`: `sfe->FixSameParameter(edge, face)`, then read the
    // pcurve, its range and the tolerance back.
    fix_same_parameter(&edge, face);
    let (nc2d, nf, nl) = curve_on_surface_oriented(&edge, face, false)?;
    *c2d = nc2d;
    *first = nf;
    *last = nl;
    Some(BRepTool::edge_tolerance(&edge))
}

/// `RemoveLoop(E, face, IP, tolfact, prec, RemoveLoop3d)`
/// (`ShapeFix_Wire.cxx:2283-2524`).
///
/// Cuts the pcurve loop between the two intersection parameters and stitches
/// the remaining pieces with a null-length patch. `false` when the edge is
/// closed on the face, the pcurve is missing, the loop touches an edge end, or
/// the stitched pcurve is not same-parameter within `prec`.
fn remove_loop(
    e: &Edge,
    face: &Face,
    ip: &IntRes2dIntersectionPoint,
    tolfact: f64,
    prec: f64,
    remove_loop_3d: bool,
) -> bool {
    // `cxx:2290-2293`.
    if BRepTool::is_closed_edge_face(e, face) {
        return false;
    }
    // `cxx:2296-2297`.
    let Some(crv) = BRepTool::edge_curve(e) else {
        return false;
    };
    let (f, l) = BRepTool::edge_parameters(e);

    // `cxx:2299-2306`.
    let (mut t1, mut t2) = (ip.param_on_first(), ip.param_on_second());
    if t1 > t2 {
        std::mem::swap(&mut t1, &mut t2);
    }
    // `cxx:2308-2315`: `sae.PCurve(E, face, c2d, a, b, false)`.
    let Some((c2d, mut a, mut b)) = curve_on_surface_oriented(e, face, false) else {
        return false;
    };
    // `cxx:2320-2322`: `dt = tolfact * GAC.Resolution(prec)`.
    let dt = tolfact * crv.resolution(prec);
    t1 -= dt;
    t2 += dt;
    // `cxx:2327-2330`.
    if t1 <= a || t2 >= b {
        return false;
    }

    // `cxx:2332-2348`: default plane `gp_Pln(gp_Ax3(P(0,0,0), D(Z)))`, replaced
    // by the face plane when the face is planar, otherwise the 3D curve is the
    // pcurve lifted onto that default plane.
    let mut pln = GpPln::new(GpAx3::default());
    let a_plane = BRepTool::face_surface(face).and_then(|s| s.gp_pln());
    let pcurve3d: Arc<dyn Curve> = match &a_plane {
        Some(p) => {
            pln = p.clone();
            crv.clone()
        }
        None => {
            let Some(c) = to_3d(&pln.position().ax2(), c2d.as_ref()) else {
                return false;
            };
            c
        }
    };

    // `cxx:2350-2379`: first segment `[a, t1]`, null-length patch, last segment
    // `[t2, b]`.
    let trim1: Arc<dyn Curve> = Arc::new(GeomTrimmedCurveBasis::new(pcurve3d.clone(), a, t1));
    let Some(mut connect) =
        CompCurveToBSplineCurve::from_curve(trim1.as_ref(), ParameterisationType::TgtThetaOver2)
            .ok()
    else {
        return false;
    };
    // `cxx:2357-2369`: `Geom_BSplineCurve(Poles, Knots, Mults, 1)`.
    let Some(patch) = occt_geom::GeomBSplineCurve::from_poles_knots_mults(
        vec![pcurve3d.d0(t1), pcurve3d.d0(t2)],
        vec![t1, t2],
        vec![2, 2],
        1,
    )
    .ok() else {
        return false;
    };
    if !connect.add(&patch, PCONFUSION, true, false, 0).unwrap_or(false) {
        return false;
    }
    let trim2: Arc<dyn Curve> = Arc::new(GeomTrimmedCurveBasis::new(pcurve3d.clone(), t2, b));
    if !connect.add(trim2.as_ref(), PCONFUSION, true, false, 0).unwrap_or(false) {
        return false;
    }
    // `cxx:2382-2387`.
    let Some(a_new3d) = connect.into_curve() else {
        return false;
    };
    let Some(mut bs) = to_2d(&a_new3d, &pln) else {
        return false;
    };

    let reg = GeometryRegistry::global();
    let face_key = GeometryRegistry::shape_key(&face.0);

    // `cxx:2392-2436`: the old variant, which leaves the 3D loop in place.
    if !remove_loop_3d {
        let mut newtol = 0.0f64;
        // `cxx:2397-2405`: `howMuchPCurves`.
        let nb_c2d = reg.pcurve_rep_count(&e.0);
        if nb_c2d <= 1 && a_plane.is_some() {
            // `cxx:2399-2401`: `B.UpdateEdge(E, bs, face, 0)`.
            reg.set_edge_pcurve(&e.0, face_key, bs.clone());
        } else {
            let Some(t) = try_new_pcurve(e, face, &mut bs, &mut a, &mut b) else {
                return false;
            };
            newtol = t;
        }

        // `cxx:2407-2415`.
        let tol = BRepTool::edge_tolerance(e);
        if newtol > prec.max(tol) {
            return false;
        }
        // `cxx:2417-2421`.
        if (a - f).abs() > PCONFUSION || (b - l).abs() > PCONFUSION {
            return false;
        }
        // `cxx:2422-2431`: `B.UpdateEdge(E, aNew3dCrv, max(newtol, tol))` then
        // re-run `TryNewPCurve`.
        if a_plane.is_some() {
            reg.set_edge_curve3d(&e.0, Arc::new(a_new3d.clone()), newtol.max(tol));
            let Some(t) = try_new_pcurve(e, face, &mut bs, &mut a, &mut b) else {
                return false;
            };
            newtol = t;
        }
        // `cxx:2432-2434`.
        reg.set_edge_pcurve(&e.0, face_key, bs.clone());
        reg.raise_edge_tolerance(&e.0, newtol);
        if let Some(v) = first_vertex(e) {
            update_vertex_tolerance(&v, newtol);
        }
        if let Some(v) = last_vertex(e) {
            update_vertex_tolerance(&v, newtol);
        }
        return true;
    }

    // `cxx:2438-2522`: the `RemoveLoop3d` variant. `ACS.Value(t)` is the pcurve
    // point evaluated on the surface.
    let Some(surf) = BRepTool::face_surface(face) else {
        return false;
    };
    let uv1 = c2d.d0(t1);
    let uv2 = c2d.d0(t2);
    let p1 = surf.d0(uv1.x(), uv1.y());
    let p2 = surf.d0(uv2.x(), uv2.y());
    let pcur_pnt = GpPnt::new(
        0.5 * (p1.x() + p2.x()),
        0.5 * (p1.y() + p2.y()),
        0.5 * (p1.z() + p2.z()),
    );

    // `cxx:2449-2457`: `ShapeAnalysis_TransferParametersProj SFTP(E, face)`.
    let sftp = TransferParametersProj::with_edge_face(e, face);
    let seq3d = sftp.perform(&[t1, t2, 0.5 * (t1 + t2)], false);
    if seq3d.len() < 3 {
        return false;
    }
    // `cxx:2458-2462`.
    let dist1 = pcur_pnt.distance(&crv.d0(seq3d[0]));
    let dist2 = pcur_pnt.distance(&crv.d0(seq3d[1]));
    let dist3 = pcur_pnt.distance(&crv.d0(seq3d[2]));
    let loop_removed_3d = dist3 <= dist1.max(dist2);

    let mut bs1: Option<occt_geom::GeomBSplineCurve> = None;
    if !loop_removed_3d {
        // `cxx:2468-2505`: stitch the 3D curve with a null-length patch too.
        let mut ftrim = seq3d[0];
        let mut ltrim = seq3d[1];
        ftrim -= dt;
        ltrim += dt;
        let trim3: Arc<dyn Curve> =
            Arc::new(GeomTrimmedCurveBasis::new(crv.clone(), f, ftrim));
        let Some(mut connect1) =
            CompCurveToBSplineCurve::from_curve(trim3.as_ref(), ParameterisationType::TgtThetaOver2)
                .ok()
        else {
            return false;
        };
        // `cxx:2471-2489`: `Geom_BSplineCurve(Poles1, Knots1, Mults1, 1)`.
        let Some(patch1) = occt_geom::GeomBSplineCurve::from_poles_knots_mults(
            vec![crv.d0(ftrim), crv.d0(ltrim)],
            vec![ftrim, ltrim],
            vec![2, 2],
            1,
        )
        .ok() else {
            return false;
        };
        if !connect1.add(&patch1, PCONFUSION, true, false, 0).unwrap_or(false) {
            return false;
        }
        let trim4: Arc<dyn Curve> =
            Arc::new(GeomTrimmedCurveBasis::new(crv.clone(), ltrim, l));
        if !connect1.add(trim4.as_ref(), PCONFUSION, true, false, 0).unwrap_or(false) {
            return false;
        }
        let Some(b) = connect1.into_curve() else {
            return false;
        };
        bs1 = Some(b);
    }

    // `cxx:2510-2520`.
    if let Some(b1) = bs1 {
        reg.set_edge_curve3d(&e.0, Arc::new(b1), 0.0);
    }
    reg.set_edge_pcurve(&e.0, face_key, bs.clone());
    reg.set_pcurve_range(&e.0, face_key, f, l);
    reg.set_same_range(&e.0, false);
    fix_same_parameter(e, face);
    true
}

/// `ShapeFix_Wire::FixSelfIntersectingEdge(const int num)`
/// (`ShapeFix_Wire.cxx:2708-2913`).
///
/// Only the reachable `myRemoveLoopMode < 1` arm is ported. The
/// `myRemoveLoopMode == 1` arm (`cxx:2844-2904`) is not modelled: the port has
/// no `myRemoveLoopMode` setting and the STEP reader leaves the default `-1`
/// (`ShapeFix_Wire.cxx:168`), so that arm is unreachable here.
pub fn fix_self_intersecting_edge(wire: &mut Wire, face: &Face, num: usize, preci: f64) -> bool {
    // `cxx:2710-2731`: analysis. The `FAIL1` bit of the analyzer is not
    // modelled as a status here (no `myLastFixStatus` accumulator).
    let check = check_self_intersecting_edge(wire, face, num);
    if !check.done() {
        return false;
    }

    // `cxx:2735-2744`.
    let edges = edges_of_wire(wire);
    let nb = edges.len();
    if nb == 0 || num > nb {
        return false;
    }
    let e = edges[if num > 0 { num } else { nb } - 1].clone();
    let (Some(v1), Some(v2)) = (first_vertex(&e), last_vertex(&e)) else {
        return false;
    };
    let mut tol1 = BRepTool::vertex_tolerance(&v1);
    let mut tol2 = BRepTool::vertex_tolerance(&v2);
    let pnt1 = vertex_position(&v1);
    let pnt2 = vertex_position(&v2);

    // `cxx:2747-2751`.
    let tolfact = 0.1;
    let mut f2d = 0.0f64;
    let mut l2d = 0.0f64;
    let mut c2d: Option<Arc<dyn Curve2d>> = None;
    let mut newtol;

    let mut points2d = check.points2d;
    let mut points3d = check.points3d;
    let mut done = false;
    let mut fail2 = false;

    // `cxx:2753-2841`: `myRemoveLoopMode < 1`, up to 30 iterations.
    for _iter in 0..30 {
        let mut loop_removed = false;
        let (mut prev_first, mut prev_last) = (0.0f64, 0.0f64);
        for i in 0..points2d.len() {
            // `cxx:2760-2767`.
            let pint = points3d[i];
            let dist21 = pnt1.square_distance(&pint);
            let dist22 = pnt2.square_distance(&pint);
            if dist21 < tol1 * tol1 || dist22 < tol2 * tol2 {
                continue;
            }
            newtol = 1.001 * (dist21.min(dist22)).sqrt();

            // `cxx:2770-2795`: `myGeomMode` is true.
            if c2d.is_none() {
                if let Some((c, fa, la)) = curve_on_surface_oriented(&e, face, false) {
                    c2d = Some(c);
                    f2d = fa;
                    l2d = la;
                }
            }
            let firstpar = points2d[i].param_on_first();
            let lastpar = points2d[i].param_on_second();
            if firstpar > prev_first && lastpar < prev_last {
                continue;
            }
            // `std::min(MaxTolerance(), std::max(newtol, Precision()))`.
            let tol_loop = newtol.max(preci).min(SHAPE_FIX_MAX_TOLERANCE);
            if remove_loop(&e, face, &points2d[i], tolfact, tol_loop, false) {
                done = true; // DONE4
                loop_removed = true;
                prev_first = firstpar;
                prev_last = lastpar;
                continue;
            }

            // `cxx:2797-2813`.
            if newtol < SHAPE_FIX_MAX_TOLERANCE {
                done = true; // DONE1
                if dist21 < dist22 {
                    tol1 = newtol;
                    update_vertex_tolerance(&v1, newtol);
                } else {
                    tol2 = newtol;
                    update_vertex_tolerance(&v2, newtol);
                }
            } else {
                fail2 = true; // FAIL2
            }
        }

        // `cxx:2815-2838`: after a loop removal, re-check and continue until
        // the self-intersection disappears (or the iteration cap is hit).
        if loop_removed {
            let pnts = check_self_intersecting_edge(wire, face, num);
            if !pnts.done() {
                break;
            }
            points2d = pnts.points2d;
            points3d = pnts.points3d;
            // `cxx:2829-2831`.
            let reg = GeometryRegistry::global();
            let face_key = GeometryRegistry::shape_key(&face.0);
            if let Some(c) = &c2d {
                reg.set_edge_pcurve(&e.0, face_key, c.clone());
            }
            reg.set_pcurve_range(&e.0, face_key, f2d, l2d);
        } else {
            break;
        }
    }
    let _ = fail2;
    done
}

/// `ShapeFix_Wire::FixSelfIntersection()` (`ShapeFix_Wire.cxx:1083-1280`).
///
/// The two head halves and the non-adjacent tail are wired:
/// `myFixSelfIntersectingEdgeMode` and `myFixIntersectingEdgesMode` are both
/// `-1` (`ShapeFix_Wire.cxx:166-167`), and `NeedFix(flag, need)` is
/// `flag < 0 ? need : flag > 0` (`ShapeFix_Root.lxx:101-104`), so
/// `NeedFix(-1, myClosedMode = true)` requests each. The tail
/// (`cxx:1199-1223`) is `NeedFix(myFixNonAdjacentIntersectingEdgesMode = -1)`
/// with the default `need = true`, so it always runs, and calls
/// `ShapeFix_IntersectionTool::FixSelfIntersectWire`. The `myContext` arms do
/// not exist in the port: this read-in pass has no `ShapeBuild_ReShape`.
///
/// The tail's `cxx:1217` `myAnalyzer->Load(sbwd)` needs no code here -- the
/// port re-derives the analyzer state from the wire on every check
/// ([`check_self_intersecting_edge`]). `cxx:1222` `myShape.Nullify()` and
/// `cxx:1213` `myStatusRemovedSegment` have no counterpart either: the former
/// belongs to the OCCT class-level shape cache that this pass never reads, and
/// the latter is written but never read anywhere in `ShapeFix_Wire.cxx`
/// (`:125`, `:138`, `:1213` are the only occurrences).
pub fn fix_self_intersection(wire: &mut Wire, face: &Face, preci: f64) -> bool {
    // `cxx:1090-1091`.
    let mut nb = edges_of_wire(wire).len();
    let mut done_any = false;

    // `cxx:1096-1119`: `myRemoveLoopMode < 1`.
    let mut num = 1usize;
    while num <= nb {
        if fix_self_intersecting_edge(wire, face, num, preci) {
            done_any = true;
        }
        num += 1;
    }

    // `cxx:1121-1197`: adjacent-edge self-intersections.
    let closed_mode = true; // `ShapeFix_Wire.cxx:171`
    let mut num = if closed_mode { 1usize } else { 2 };
    while nb > 1 && num <= nb {
        let st = fix_intersecting_edges(wire, face, num, preci);
        if !st.done() {
            num += 1;
            continue;
        }
        done_any = true;
        // `cxx:1154-1164`.
        if nb < 3 {
            if st.done7 {
                let _ = fix_intersecting_edges(wire, face, num, preci);
            }
            num += 1;
            continue;
        }
        // `cxx:1166-1191`: drop the removed edge(s) and restart the scan.
        if st.done4 {
            wire_remove_edge(wire, num);
        }
        if st.done3 {
            wire_remove_edge(wire, if num > 1 { num - 1 } else { nb + num - 1 });
        }
        if st.done4 || st.done3 {
            num = if closed_mode { 1 } else { 2 };
            nb = edges_of_wire(wire).len();
        } else {
            let _ = fix_intersecting_edges(wire, face, num, preci);
            num += 1;
        }
    }

    // `cxx:1199-1223`: the non-adjacent tail. The gate
    // `NeedFix(myFixNonAdjacentIntersectingEdgesMode = -1)` takes the default
    // `need = true` (`ShapeFix_Root.lxx:101-104`), so it always runs.
    //
    // `ShapeFix_IntersectionTool ITool(Context(), Precision())` (`cxx:1203`)
    // keeps the `maxtol = 1.0` default (`ShapeFix_IntersectionTool.hxx:46`).
    // That tool only ever *writes* its context -- every use in
    // `ShapeFix_IntersectionTool.cxx` is a `Replace` (it never calls
    // `Apply`) -- and the substitutions die with `ITool` here, so the
    // substitution-free context is faithful.
    let mut itool_ctx = IdentityReShape;
    let st = fix_self_intersect_wire(wire, face, SHAPE_FIX_INTERSECTION_MAX_TOL, &mut itool_ctx);
    if st.done {
        done_any = true; // `cxx:1207`: `myStatusSelfIntersection |= DONE5`
    }
    // `cxx:1209-1223` (`if (NbSplit > 0 || NbRemoved > 0)`) reloads the
    // analyzer from the reworked wire and, with a context, refreshes the
    // class-level wire/shape caches. The port's analyzer is derived from the
    // wire on each check and this pass keeps no shape cache, so the block
    // reduces to the two status no-ops documented above.

    done_any
}
