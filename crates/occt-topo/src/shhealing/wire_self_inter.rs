//! `ShapeAnalysis_Wire` self-intersection checks
//! (`ShapeAnalysis_Wire.cxx:1245-1542`).
//!
//! Both checks feed `ShapeFix_Wire::FixSelfIntersectingEdge`
//! (`ShapeFix_Wire.cxx:2708`) and `ShapeFix_Wire::FixIntersectingEdges`
//! (`ShapeFix_Wire.cxx:2966`, `:3271`): they intersect the edge pcurves with
//! `Geom2dInt_GInter` and hand back the `IntRes2d_IntersectionPoint` sequence
//! plus, per point, the 3D point the fix must move or split on.
//!
//! The two-curve checks reach the dedicated conic-conic overloads ported in
//! `occt_geom2d::geom2d_int` (see specs/_a3n00_gap_analysis.md 9.628);
//! `CheckSelfIntersectingEdge` uses the single-curve `Perform` overload, whose
//! conic arm (`IntCurve_IntCurveCurveGen.gxx:97-106`) intersects nothing by
//! construction and whose non-conic arm runs the polygon intersector ported and
//! wired into `Geom2dInt_GInter::perform_curve` (see
//! specs/_a3n00_gap_analysis.md 9.629 and 9.631).

use super::prelude::*;

use occt_core::gp::{GpDir, GpLin, GpPnt, GpVec};
use occt_core::intres2d::{IntRes2dDomain, IntRes2dIntersectionPoint, IntRes2dPosition};
use occt_core::precision::PCONFUSION;
use occt_geom::Surface;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::geom2d_int::Geom2dIntGInter;

use crate::abs::Orientation;
use crate::boptools_2d::curve_on_surface_oriented;
use crate::brep_tool::BRepTool;
use crate::shhealing::wire_fix::{
    first_vertex, fix_same_parameter, last_vertex, vertices_coincide, SHAPE_FIX_MAX_TOLERANCE,
};
use super::split_tool::cut_edge;

/// `double tolint = 1.0e-10;` (`ShapeAnalysis_Wire.cxx:1296`, `:1426`).
const TOLINT: f64 = 1.0e-10;

/// Out-parameters of `ShapeAnalysis_Wire::CheckSelfIntersectingEdge`
/// (`ShapeAnalysis_Wire.cxx:1269-1342`).
#[derive(Clone, Default)]
pub struct SelfIntersectingEdgeCheck {
    /// `points2d` (`cxx:1270`).
    pub points2d: Vec<IntRes2dIntersectionPoint>,
    /// `points3d` (`cxx:1271`).
    pub points3d: Vec<GpPnt>,
    /// `ShapeExtend_FAIL1` (`cxx:1289`): no pcurve on the face.
    pub fail1: bool,
    /// `ShapeExtend_FAIL2` (`cxx:1312`): a vertex of the edge is null.
    pub fail2: bool,
    /// `ShapeExtend_DONE1` (`cxx:1337`): at least one accepted point.
    pub done1: bool,
}

impl SelfIntersectingEdgeCheck {
    /// `LastCheckStatus(ShapeExtend_FAIL)` (`cxx:2726`, `:2728`).
    pub fn failed(&self) -> bool {
        self.fail1 || self.fail2
    }

    /// `LastCheckStatus(ShapeExtend_DONE)` (`cxx:2729`, `cxx:2821`).
    pub fn done(&self) -> bool {
        self.done1
    }
}

/// `GetPointOnEdge(edge, surf, Crv2d, param)` (`ShapeAnalysis_Wire.cxx:1247-1264`).
///
/// The 3D point comes from the edge's own curve *only* when the edge is
/// `SameParameter`; otherwise `surf->Adaptor3d()->Value(p2d.X(), p2d.Y())`.
pub(super) fn point_on_edge(edge: &Edge, surf: &dyn Surface, c2d: &dyn Curve2d, param: f64) -> GpPnt {
    if BRepTool::same_parameter(edge) {
        // `cxx:1253-1261`: `BRep_Tool::Curve(edge, L, f, l)` then
        // `ConS->Value(param).Transformed(L.Transformation())`.
        if let Some(curve) = BRepTool::edge_curve_world(edge) {
            return curve.d0(param);
        }
    }
    let p2d = c2d.d0(param);
    surf.d0(p2d.x(), p2d.y())
}

/// `ShapeAnalysis_Wire::CheckSelfIntersectingEdge(num, points2d, points3d)`
/// (`ShapeAnalysis_Wire.cxx:1269-1342`).
///
/// `num` is the 1-based wire position (`cxx:1283`: `0` means the last edge).
/// `IsReady()` is satisfied by construction here - the port's callers always
/// hold both the wire and the face - so the guard at `cxx:1277-1281` is not
/// repeated.
pub fn check_self_intersecting_edge(
    wire: &Wire,
    face: &Face,
    num: usize,
) -> SelfIntersectingEdgeCheck {
    let mut out = SelfIntersectingEdgeCheck::default();
    let edges = edges_of_wire(wire);
    let nb = edges.len();
    if nb == 0 || num > nb {
        return out;
    }
    // `cxx:1283`.
    let edge = edges[if num > 0 { num } else { nb } - 1].clone();

    // `cxx:1286-1292`: `sae.PCurve(edge, myFace, Crv, a, b, false)`.
    let Some((crv, a, b)) = curve_on_surface_oriented(&edge, face, false) else {
        out.fail1 = true;
        return out;
    };
    if (a - b).abs() <= PCONFUSION {
        return out;
    }

    // `cxx:1296-1302`: the domain spans `[a, b]` with `tolint` end tolerances,
    // and the intersector is the *single-curve* overload.
    let domain = IntRes2dDomain::bounded(&crv.d0(a), a, TOLINT, &crv.d0(b), b, TOLINT);
    let mut inter = Geom2dIntGInter::new();
    inter.perform_curve(crv.as_ref(), &domain, TOLINT, TOLINT);
    if !inter.is_done() {
        return out;
    }

    // `cxx:1308-1313`.
    let (Some(v1), Some(v2)) = (first_vertex(&edge), last_vertex(&edge)) else {
        out.fail2 = true;
        return out;
    };
    let tol1 = BRepTool::vertex_tolerance(&v1);
    let tol2 = BRepTool::vertex_tolerance(&v2);
    let pnt1 = vertex_position(&v1);
    let pnt2 = vertex_position(&v2);

    // `cxx:1320-1338`. The surface of `mySurf` is the face's surface
    // (`ShapeAnalysis_Wire::Load(face, ...)`; `ShapeAnalysis_Surface::Adaptor3d`
    // is the underlying `Geom_Surface`).
    let Some(surf) = BRepTool::face_surface(face) else {
        return out;
    };
    for i in 1..=inter.nb_points() {
        let ip = inter.point(i);
        let tr1 = ip.transition_of_first();
        let tr2 = ip.transition_of_second();
        if tr1.position_on_curve() != IntRes2dPosition::Middle
            && tr2.position_on_curve() != IntRes2dPosition::Middle
        {
            continue;
        }
        let pint = point_on_edge(&edge, surf.as_ref(), crv.as_ref(), ip.param_on_first());
        let dist21 = pnt1.square_distance(&pint);
        let dist22 = pnt2.square_distance(&pint);
        if dist21 > tol1 * tol1 && dist22 > tol2 * tol2 {
            out.points2d.push(ip);
            out.points3d.push(pint);
            out.done1 = true;
        }
    }
    out
}

/// Out-parameters of `ShapeAnalysis_Wire::CheckIntersectingEdges`
/// (`ShapeAnalysis_Wire.cxx:1360-1541`).
#[derive(Clone, Default)]
pub struct IntersectingEdgesCheck {
    /// `points2d` (`cxx:1362`).
    pub points2d: Vec<IntRes2dIntersectionPoint>,
    /// `points3d` (`cxx:1363`).
    pub points3d: Vec<GpPnt>,
    /// `errors` (`cxx:1364`): `0.5 * pi1.Distance(pi2)` per point.
    pub errors: Vec<f64>,
    /// `ShapeExtend_FAIL1` (`cxx:1387`): a junction vertex is null.
    pub fail1: bool,
    /// `ShapeExtend_FAIL2` (`cxx:1391`): the two junction vertices differ.
    pub fail2: bool,
    /// `ShapeExtend_FAIL3` (`cxx:1403`, `:1407`): a pcurve is missing.
    pub fail3: bool,
    /// `ShapeExtend_DONE1` (`cxx:1538`): at least one accepted point.
    pub done1: bool,
}

impl IntersectingEdgesCheck {
    /// `LastCheckStatus(ShapeExtend_FAIL)` (`cxx:1126-1132`).
    pub fn failed(&self) -> bool {
        self.fail1 || self.fail2 || self.fail3
    }

    /// `LastCheckStatus(ShapeExtend_DONE)` (`cxx:1136`).
    pub fn done(&self) -> bool {
        self.done1
    }
}

/// `ShapeAnalysis_Wire::CheckIntersectingEdges(num, points2d, points3d, errors)`
/// (`ShapeAnalysis_Wire.cxx:1360-1541`).
///
/// `num` is the 1-based position of the *second* edge; the pair checked is
/// `(num - 1, num)` (`cxx:1375-1379`). `preci` is
/// `ShapeFix_Root::Precision()` (`cxx:1458`), `ShapeFix_Root.cxx:26`.
pub fn check_intersecting_edges(
    wire: &Wire,
    face: &Face,
    num: usize,
    preci: f64,
) -> IntersectingEdgesCheck {
    let mut out = IntersectingEdgesCheck::default();
    let edges = edges_of_wire(wire);
    let nb = edges.len();
    // `cxx:1370-1373`: `!IsReady() || NbEdges() < 2`.
    if nb < 2 || num > nb {
        return out;
    }

    // `cxx:1375-1379`. `myWire->Edge(n)` is the stored order, not the
    // `WireData` order (`ShapeAnalysis_Wire.cxx:1377-1378`).
    let n2 = if num > 0 { num } else { nb };
    let n1 = if n2 > 1 { n2 - 1 } else { nb };
    let e1 = edges[n1 - 1].clone();
    let e2 = edges[n2 - 1].clone();

    // `cxx:1381-1392`.
    let (Some(v1), Some(v2)) = (last_vertex(&e1), first_vertex(&e2)) else {
        out.fail1 = true;
        return out;
    };
    if !vertices_coincide(&v1, &v2) {
        out.fail2 = true;
        return out;
    }

    // `cxx:1394-1395`.
    let vp = first_vertex(&e1);
    let vn = last_vertex(&e2);

    // `cxx:1399-1408`.
    let Some((crv1, a1, b1)) = curve_on_surface_oriented(&e1, face, false) else {
        out.fail3 = true;
        return out;
    };
    let Some((crv2, a2, b2)) = curve_on_surface_oriented(&e2, face, false) else {
        out.fail3 = true;
        return out;
    };
    if (a1 - b1).abs() <= PCONFUSION || (a2 - b2).abs() <= PCONFUSION {
        return out;
    }

    // `cxx:1416-1420`.
    let is_forward1 = e1.0.orientation() != Orientation::Reversed;
    let is_forward2 = e2.0.orientation() != Orientation::Reversed;
    let tol0 = BRepTool::vertex_tolerance(&v1).max(BRepTool::vertex_tolerance(&v2));
    let tol = tol0;
    let pnt = vertex_position(&v1);

    // `cxx:1428-1437`.
    let d1 = IntRes2dDomain::bounded(&crv1.d0(a1), a1, TOLINT, &crv1.d0(b1), b1, TOLINT);
    let d2 = IntRes2dDomain::bounded(&crv2.d0(a2), a2, TOLINT, &crv2.d0(b2), b2, TOLINT);

    // `cxx:1438-1448`: the intersection algorithm is not symmetrical, so the
    // edge with the lower number is always the *first* argument - for the
    // `(last, first)` junction (`num == 1`) the roles are swapped.
    let mut inter = Geom2dIntGInter::new();
    if num == 1 {
        inter.perform(crv2.as_ref(), &d2, crv1.as_ref(), &d1, TOLINT, TOLINT);
    } else {
        inter.perform(crv1.as_ref(), &d1, crv2.as_ref(), &d2, TOLINT, TOLINT);
    }
    if !inter.is_done() {
        return out;
    }

    // `cxx:1455-1458`.
    let tole = if BRepTool::same_parameter(&e1) {
        BRepTool::edge_tolerance(&e1)
    } else {
        tol0
    }
    .max(if BRepTool::same_parameter(&e2) {
        BRepTool::edge_tolerance(&e2)
    } else {
        tol0
    });
    let tolt = tol.min(tole.max(preci));

    // `cxx:1460-1461`: `isLacking` is computed once, at the first iteration.
    let mut is_lacking: Option<bool> = None;

    let Some(surf) = BRepTool::face_surface(face) else {
        return out;
    };

    let nb_points = inter.nb_points();
    let nb_segments = inter.nb_segments();
    for i in 1..=(nb_points + nb_segments) {
        // `cxx:1466-1486`: points first, then the segments' first or last point.
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
        let tr1 = ip.transition_of_first();
        let tr2 = ip.transition_of_second();
        if tr1.position_on_curve() != IntRes2dPosition::Middle
            && tr2.position_on_curve() != IntRes2dPosition::Middle
        {
            continue;
        }

        // `cxx:1493-1494`.
        let param1 = if num == 1 { ip.param_on_second() } else { ip.param_on_first() };
        let param2 = if num == 1 { ip.param_on_first() } else { ip.param_on_second() };

        // `cxx:1498-1503`.
        if a1 - param1 > PCONFUSION
            || param1 - b1 > PCONFUSION
            || a2 - param2 > PCONFUSION
            || param2 - b2 > PCONFUSION
        {
            continue;
        }

        // `cxx:1508-1510`.
        let pi1 = point_on_edge(&e1, surf.as_ref(), crv1.as_ref(), param1);
        let pi2 = point_on_edge(&e2, surf.as_ref(), crv2.as_ref(), param2);
        let pint = GpPnt::new(
            0.5 * (pi1.x() + pi2.x()),
            0.5 * (pi1.y() + pi2.y()),
            0.5 * (pi1.z() + pi2.z()),
        );
        let di1 = pi1.square_distance(&pnt);
        let di2 = pi2.square_distance(&pnt);
        let dist2 = di1.max(di2);

        // `cxx:1517-1524`: the once-only `isLacking` test, done like BRepCheck.
        if is_lacking.is_none() {
            let end1 = crv1.d0(if is_forward1 { b1 } else { a1 });
            let end2 = crv2.d0(if is_forward2 { a2 } else { b2 });
            let tol2d = 2.0
                * occt_geom::approx_same_parameter::u_resolution(surf.as_ref(), tol)
                    .max(occt_geom::approx_same_parameter::v_resolution(surf.as_ref(), tol));
            is_lacking = Some(end1.square_distance(&end2) >= tol2d * tol2d);
        }
        let is_lacking = is_lacking.unwrap_or(false);

        // `cxx:1526-1539`.
        let accept = (dist2 > tolt * tolt || is_lacking)
            && match (&vp, &vn) {
                (Some(vp), Some(vn)) => {
                    !vertices_coincide(vp, vn)
                        || dist2 < pint.square_distance(&vertex_position(vp))
                }
                // `BRepTools::Compare(null, null)` - not reachable in practice
                // because `vp`/`vn` come from the same edges as `v1`/`v2`.
                _ => true,
            };
        if accept {
            out.points2d.push(ip);
            out.points3d.push(pint);
            out.errors.push(0.5 * pi1.distance(&pi2));
            out.done1 = true;
        }
    }
    out
}

// ---------------------------------------------------------------------------
// ShapeFix_Wire::FixSelfIntersectingEdge / FixIntersectingEdges actions
// ---------------------------------------------------------------------------

/// `ComputeLocalDeviation(edge, pint, pnt, f, l, face)`
/// (`ShapeFix_Wire.cxx:2917-2964`).
///
/// Samples `NSEG = 10` points of the edge's 3D curve on `[f, l]` and returns
/// the maximum distance to the line through `pint` and `pnt`. `RealLast()`
/// (`f64::MAX`) when the edge carries no 3D curve (`cxx:2925-2931`).
fn compute_local_deviation(
    edge: &Edge,
    pint: &GpPnt,
    pnt: &GpPnt,
    f_in: f64,
    l_in: f64,
    face: &Face,
) -> f64 {
    // `cxx:2922-2931`: `sae.Curve3d(edge, c3d, a, b, false)` - the world curve
    // `BRep_Tool::Curve(edge, L, a, b)` returns, as in `point_on_edge` above.
    let Some(c3d) = BRepTool::edge_curve_world(edge) else {
        return f64::MAX;
    };
    let (a, b) = BRepTool::edge_parameters(edge);

    // `cxx:2932`: `gp_Lin line(pint, gp_Vec(pint, pnt))`.
    let Ok(dir) = GpDir::from_vec(&GpVec::from_pnts(pint, pnt)) else {
        return f64::MAX;
    };
    let line = GpLin::from_pnt_dir(*pint, dir);

    // `cxx:2935-2946`: a pcurve that is a trimmed line remaps the parameters
    // from the pcurve's own range onto the edge's 3D range.
    let mut f = f_in;
    let mut l = l_in;
    if let Some((crv, fp, lp)) = curve_on_surface_oriented(edge, face, false) {
        if crv.trimmed_basis().map(|basis| basis.is_line()).unwrap_or(false) {
            f = a + (f - fp) * (b - a) / (lp - fp);
            l = a + (l - fp) * (b - a) / (lp - fp);
        }
    }

    // `cxx:2948-2962`.
    const NSEG: usize = 10;
    let step = (l - f) / NSEG as f64;
    let mut dev = 0.0f64;
    for i in 1..NSEG {
        let p = c3d.d0(f + i as f64 * step);
        let d = line.distance(&p);
        if dev < d {
            dev = d;
        }
    }
    dev
}

/// `BRep_Builder::UpdateVertex(V, P, Tol)` (`BRep_Builder.cxx:1203-1212`):
/// `BRep_TVertex::Pnt` plus `UpdateTolerance`, which is a max
/// (`BRep_TVertex.hxx`). Locations are identity throughout this port, so the
/// `V.Location().Inverted()` transform of the C++ call is a no-op.
fn update_vertex(v: &Vertex, p: GpPnt, tol: f64) {
    v.set_point(p);
    v.set_tolerance(BRepTool::vertex_tolerance(v).max(tol));
}

/// `BRep_Builder::UpdateVertex(V, Tol)` (`BRep_Builder.cxx:1442-1452`).
fn update_vertex_tolerance(v: &Vertex, tol: f64) {
    v.set_tolerance(BRepTool::vertex_tolerance(v).max(tol));
}

/// `BRep_Builder::UpdateEdge(E, Tol)` (`BRep_Builder.cxx:999-1007`):
/// `BRep_TEdge::UpdateTolerance` is a max (`BRep_TEdge.hxx:50-52`).
fn update_edge_tolerance(e: &Edge, tol: f64) {
    let reg = GeometryRegistry::global();
    reg.set_edge_tolerance(&e.0, BRepTool::edge_tolerance(e).max(tol));
}

/// `ShapeExtend` status bits of `ShapeFix_Wire::FixIntersectingEdges(const int)`
/// (`ShapeFix_Wire.cxx:2966-3269`).
#[derive(Clone, Copy, Default)]
pub struct FixIntersectingEdgesStatus {
    /// `ShapeExtend_FAIL1` (`cxx:2985`): the check reported a failure.
    pub fail1: bool,
    /// `ShapeExtend_FAIL2` (`cxx:3202`): the intersection could not be fixed.
    pub fail2: bool,
    /// `ShapeExtend_DONE1` (`cxx:3175`, `:3186`, `:3195`).
    pub done1: bool,
    /// `ShapeExtend_DONE2` (`cxx:3170`, `:3182`).
    pub done2: bool,
    /// `ShapeExtend_DONE3` (`cxx:3138`, edge 1 removed).
    pub done3: bool,
    /// `ShapeExtend_DONE4` (`cxx:3153`, edge 2 removed).
    pub done4: bool,
    /// `ShapeExtend_DONE6` (`cxx:3117`, edge tolerances increased).
    pub done6: bool,
    /// `ShapeExtend_DONE7` (`cxx:3254`, a cut edge was made same-parameter).
    pub done7: bool,
}

impl FixIntersectingEdgesStatus {
    /// `LastFixStatus(ShapeExtend_DONE)`.
    pub fn done(&self) -> bool {
        self.done1
            || self.done2
            || self.done3
            || self.done4
            || self.done6
            || self.done7
    }
}

/// `ShapeFix_Wire::FixIntersectingEdges(const int num)`
/// (`ShapeFix_Wire.cxx:2966-3269`).
///
/// `num` is the 1-based wire position of the second edge of the pair
/// `(num - 1, num)`. The `myContext` (`ShapeBuild_ReShape`) arms of the cxx
/// (`cxx:3005-3010`, `:3073-3108`, `:3216-3228`, `:3241-3243`, `:3248-3250`)
/// are not modelled - the port always runs with a null context - so the edge
/// copies they make are absent and `myFixEdge->FixSameParameter` runs on the
/// call-site face only (see `fix_same_parameter`).
pub fn fix_intersecting_edges(
    wire: &Wire,
    face: &Face,
    num: usize,
    preci: f64,
) -> FixIntersectingEdgesStatus {
    let mut st = FixIntersectingEdgesStatus::default();
    let edges = edges_of_wire(wire);
    let nb_edges = edges.len();
    // `cxx:2966-2971`.
    if nb_edges < 2 || num > nb_edges {
        return st;
    }

    // `cxx:2972-2989`.
    let check = check_intersecting_edges(wire, face, num, preci);
    if check.failed() {
        st.fail1 = true;
    }
    if !check.done() {
        return st;
    }

    // `cxx:2997-3001`.
    let n2 = if num > 0 { num } else { nb_edges };
    let n1 = if n2 > 1 { n2 - 1 } else { nb_edges };
    let e1 = edges[n1 - 1].clone();
    let e2 = edges[n2 - 1].clone();

    // `cxx:3009-3010`.
    let is_forward1 = e1.0.orientation() == Orientation::Forward;
    let is_forward2 = e2.0.orientation() == Orientation::Forward;

    // `cxx:3011-3013`: `BRep_Tool::Range(E, Face(), a, b)`, the pcurve range.
    let (a1, b1) = curve_on_surface_oriented(&e1, face, false)
        .map(|(_, a, b)| (a, b))
        .unwrap_or((0.0, 0.0));
    let (a2, b2) = curve_on_surface_oriented(&e2, face, false)
        .map(|(_, a, b)| (a, b))
        .unwrap_or((0.0, 0.0));

    // `cxx:3015-3021`.
    let (Some(vp), Some(v1), Some(v2), Some(vn)) = (
        first_vertex(&e1),
        last_vertex(&e1),
        first_vertex(&e2),
        last_vertex(&e2),
    ) else {
        return st;
    };
    let mut tol = BRepTool::vertex_tolerance(&v1);
    let mut pnt = vertex_position(&v1);

    // `cxx:3023-3027`.
    let mut prev_range1 = f64::MAX;
    let mut prev_range2 = f64::MAX;
    let mut cut_edge1 = false;
    let mut cut_edge2 = false;
    let mut is_cut_line = false;

    // `cxx:3032-3204`.
    let nb = check.points3d.len();
    for i in 0..nb {
        let ip = check.points2d[i].clone();
        // `cxx:3035-3036`.
        let param1 = if num == 1 { ip.param_on_second() } else { ip.param_on_first() };
        let param2 = if num == 1 { ip.param_on_first() } else { ip.param_on_second() };

        // `cxx:3038-3042`.
        let new_range1 = ((if is_forward1 { a1 } else { b1 }) - param1).abs();
        let new_range2 = ((if is_forward2 { b2 } else { a2 }) - param2).abs();
        if new_range1 > prev_range1 && new_range2 > prev_range2 {
            continue;
        }

        // `cxx:3044-3046`.
        let pint = check.points3d[i];
        let rad = check.errors[i];
        let mut newtol = 1.0001 * (pnt.distance(&pint) + rad);

        // `cxx:3052-3055`: `myTopoMode` is forced true by `ShapeFix_Shape`
        // (`ShapeFix_Shape.cxx:200`) for the face path.
        let mut loc_may_edit = true;
        if newtol > tol {
            // `cxx:3056-3060`.
            let te1 = rad
                + compute_local_deviation(
                    &e1,
                    &pint,
                    &pnt,
                    param1,
                    if is_forward1 { b1 } else { a1 },
                    face,
                );
            let te2 = rad
                + compute_local_deviation(
                    &e2,
                    &pint,
                    &pnt,
                    if is_forward2 { a2 } else { b2 },
                    param2,
                    face,
                );
            let maxte = te1.max(te2);
            // `cxx:3062-3121`.
            if maxte < SHAPE_FIX_MAX_TOLERANCE && maxte < newtol {
                if BRepTool::edge_tolerance(&e1) < te1 || BRepTool::edge_tolerance(&e2) < te2 {
                    // `cxx:3073-3108`: the context copy arm is absent; OCCT
                    // would `Context()->CopyVertex` and rebuild E1/E2, which is
                    // a no-op with a null context.
                    // `cxx:3109-3114`.
                    let tol1 = 1.000001 * te1;
                    let tol2 = 1.000001 * te2;
                    update_edge_tolerance(&e1, tol1);
                    update_vertex_tolerance(&vp, tol1);
                    update_vertex_tolerance(&v1, tol1);
                    update_edge_tolerance(&e2, tol2);
                    update_vertex_tolerance(&v2, tol2);
                    update_vertex_tolerance(&vn, tol2);
                    st.done6 = true;
                    loc_may_edit = false;
                }
                newtol = 1.000001 * maxte;
            }
        }

        // `cxx:3122-3203`.
        if loc_may_edit || newtol <= SHAPE_FIX_MAX_TOLERANCE {
            prev_range1 = new_range1;
            prev_range2 = new_range2;
            if loc_may_edit {
                // `cxx:3129`.
                newtol = 1.0001 * (pnt.distance(&pint) + rad);
                if std::env::var("PROJDIAG").is_ok() { if let (Some((c1, _, _)), Some((c2, _, _))) = (curve_on_surface_oriented(&e1, face, false), curve_on_surface_oriented(&e2, face, false)) { let v1p = if is_forward1 { b1 } else { a1 }; let v2p = if is_forward2 { a2 } else { b2 }; let q1 = c1.d0(param1); let q1v = c1.d0(v1p); let q2 = c2.d0(param2); let q2v = c2.d0(v2p); eprintln!("PC1 param1={:.6} vtxparam1={:.6} uv1=({:.6},{:.6}) uvV1=({:.6},{:.6}) PC2 param2={:.6} vtxparam2={:.6} uv2=({:.6},{:.6}) uvV2=({:.6},{:.6}) |uv1-uv2|={:.6e} |V1-V2|uv={:.6e} pint=({:.3},{:.3},{:.3}) rad={:.3e} newtol={:.3e} pnt=({:.3},{:.3},{:.3})", param1, v1p, q1.x(), q1.y(), q1v.x(), q1v.y(), param2, v2p, q2.x(), q2.y(), q2v.x(), q2v.y(), ((q1.x()-q2.x()).powi(2)+(q1.y()-q2.y()).powi(2)).sqrt(), ((q1v.x()-q2v.x()).powi(2)+(q1v.y()-q2v.y()).powi(2)).sqrt(), pint.x(), pint.y(), pint.z(), rad, newtol, pnt.x(), pnt.y(), pnt.z()); } }                // `cxx:3132-3146`: `ShapeFix_SplitTool aTool; aTool.CutEdge(...)`.
                if !cut_edge(&e1, if is_forward1 { a1 } else { b1 }, param1, face, &mut is_cut_line)
                {
                    if is_same(&v1.0, &vp.0) {
                        st.done3 = true;
                    } else {
                        loc_may_edit = false;
                    }
                } else {
                    cut_edge1 = true;
                }
                // `cxx:3150-3163`.
                if !cut_edge(&e2, if is_forward2 { b2 } else { a2 }, param2, face, &mut is_cut_line)
                {
                    if is_same(&v2.0, &vn.0) {
                        st.done4 = true;
                    } else {
                        loc_may_edit = false;
                    }
                } else {
                    cut_edge2 = true;
                }
            }

            // `cxx:3167-3198`.
            if loc_may_edit
                && new_range1 <= prev_range1
                && new_range2 <= prev_range2
                && BRepTool::same_parameter(&e1)
                && BRepTool::same_parameter(&e2)
            {
                st.done2 = true;
                pnt = pint;
                if tol <= rad {
                    st.done1 = true;
                    tol = 1.001 * rad;
                }
            } else if is_cut_line {
                st.done2 = true;
                pnt = pint;
                if tol <= rad {
                    st.done1 = true;
                    tol = 1.001 * rad;
                }
            } else if tol < newtol {
                // `cxx:3194-3198`: CCI60005-brep.igs.
                st.done1 = true;
                tol = newtol;
            }
        } else {
            // `cxx:3202`.
            st.fail2 = true;
        }
    }

    // `cxx:3207-3209`.
    if !st.done() {
        return st;
    }

    // `cxx:3211-3228`: the `isChangedEdge` arm only triggers through a
    // non-null context, so both arms reduce to the same two updates here.
    update_vertex(&v1, pnt, tol);
    update_vertex(&v2, pnt, tol);

    // `cxx:3237-3257`: `myFixEdge->FixSameParameter(E1/E2)`.
    if cut_edge1 {
        fix_same_parameter(&e1, face);
    }
    if cut_edge2 && !is_cut_line {
        fix_same_parameter(&e2, face);
    }
    if cut_edge1 || cut_edge2 {
        st.done7 = true;
    }

    st
}
