//! IntTools range/shrunk-range utilities and a lightweight intersection
//! context.
//!
//! Ports the TKBO `IntTools_ShrunkRange`, the range/vertex helpers of
//! `IntTools_Tools`, and a small `IntTools_Context`:
//!
//! - [`ShrunkRange`] computes the *working (shrunk) range* `[t1, t2]` of an
//!   edge's 3D curve — the part of the curve that is not covered by the
//!   tolerance spheres of its endpoint vertices. This is the OCCT
//!   `BRepLib::FindValidRange` walk (step along the curve until the point
//!   leaves the vertex sphere, then bisect to the exact exit), followed by an
//!   arc-length computation of the shrunk span (`GCPnts_AbscissaPoint::Length`
//!   here approximated with `occt_core::elib::measure::curve_arc_length`).
//!   Degenerated edges, edges without a registered 3D curve and unbounded
//!   parameter ranges return `Err` rather than panicking.
//! - [`IntToolsTools`] is a namespace of small static helpers used by the
//!   boolean pipeline (`is_in_range`, `middle_point`, `is_vertex`,
//!   `compute_tolerance`).
//! - [`IntContext`] caches a working edge/face pair and answers point-on-face
//!   and surface-projection queries. The face/solid classifiers (FClass2d,
//!   BRepClass3d_SolidClassifier) arrive in wave B; this wave leaves them as
//!   `None` shells.
//!
//! Note: `occt_core::gcpnts` was checked for a `GCPnts_AbscissaPoint` port but
//! only carries the point-distribution samplers; the arc-length *measurement*
//! is taken from `occt_core::elib::measure::curve_arc_length` instead.

use occt_core::elib::measure::curve_arc_length;
use occt_core::gp::{GpPnt, GpPnt2d};
use occt_core::precision::{CONFUSION, PCONFUSION};
use occt_geom::Curve;

use crate::abs::ShapeType;
use crate::brep_surface::surface_closest_params;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::topo_tools_full::edge_vertices;

// ---------------------------------------------------------------------------
// ShrunkRange
// ---------------------------------------------------------------------------

/// Port of `IntTools_ShrunkRange`.
///
/// The class provides the computation of a working (shrunk) range `[t1, t2]`
/// for the 3D curve of the edge: the range of the curve that is *not* covered
/// by the tolerance spheres of the edge's endpoint vertices. Any later edge
/// splitting/intersection work should be restricted to this range so that the
/// vertex neighbourhoods are left untouched.
#[derive(Debug, Clone)]
pub struct ShrunkRange {
    edge: Option<Edge>,
    face: Option<Face>,
    first: f64,
    last: f64,
    shrunk_first: f64,
    shrunk_last: f64,
    is_done: bool,
    is_splittable: bool,
    length: f64,
    error: Option<String>,
}

impl ShrunkRange {
    /// Create an empty, not-yet-computed shrunk range.
    pub fn new() -> Self {
        Self {
            edge: None,
            face: None,
            first: -99.0,
            last: -99.0,
            shrunk_first: -99.0,
            shrunk_last: -99.0,
            is_done: false,
            is_splittable: false,
            length: 0.0,
            error: None,
        }
    }

    /// Set the working edge/face pair and the edge tolerance, then compute the
    /// shrunk range.
    ///
    /// `tol` plays the role of `BRep_Tool::Tolerance(edge)` in the OCCT
    /// `Perform()`; vertex tolerances are derived from the edge's endpoint
    /// vertices and floored at `tol`. Returns `Err` (leaving `is_done() ==
    /// false`) for degenerated edges, edges without a registered 3D curve,
    /// unbounded parameter ranges, and micro-edges whose vertex spheres cover
    /// the whole range.
    pub fn set_shrunk_range(&mut self, edge: &Edge, face: &Face, tol: f64) -> Result<(), String> {
        self.is_done = false;
        self.is_splittable = false;
        self.shrunk_first = 0.0;
        self.shrunk_last = 0.0;
        self.length = 0.0;
        self.error = None;
        self.edge = Some(edge.clone());
        self.face = Some(face.clone());

        if BRepTool::is_degenerated(edge) {
            return self.fail("edge is degenerated");
        }
        let Some(curve) = BRepTool::edge_curve(edge) else {
            return self.fail("edge has no registered 3D curve");
        };
        let (first, last) = BRepTool::edge_parameters(edge);
        if !first.is_finite() || !last.is_finite() {
            return self.fail("edge has an unbounded parameter range");
        }
        self.first = first;
        self.last = last;
        if last - first < PCONFUSION {
            return self.fail("edge parameter range is too short");
        }

        // Endpoint vertex points and tolerances. OCCT increases the vertex
        // tolerances on Precision::Confusion() to keep correspondence with the
        // intersection precision, and floors them at the edge tolerance.
        let (v1, v2) = edge_vertices(edge);
        let p1 = v1.as_ref().map(|v| BRepTool::vertex_point(v)).unwrap_or_else(|| curve.d0(first));
        let p2 = v2.as_ref().map(|v| BRepTool::vertex_point(v)).unwrap_or_else(|| curve.d0(last));
        let a_tol_e = tol.max(0.0);
        let mut a_tol_v1 = v1.as_ref().map(|v| BRepTool::vertex_tolerance(v)).unwrap_or(0.0);
        let mut a_tol_v2 = v2.as_ref().map(|v| BRepTool::vertex_tolerance(v)).unwrap_or(0.0);
        if a_tol_v1 < a_tol_e {
            a_tol_v1 = a_tol_e;
        }
        if a_tol_v2 < a_tol_e {
            a_tol_v2 = a_tol_e;
        }
        a_tol_v1 += CONFUSION;
        a_tol_v2 += CONFUSION;

        let Some((s1, s2)) = find_valid_range(
            curve.as_ref(),
            first,
            last,
            &p1,
            a_tol_v1,
            &p2,
            a_tol_v2,
        ) else {
            return self.fail("no valid range: vertex spheres cover the whole edge");
        };
        if s2 - s1 < PCONFUSION {
            return self.fail("shrunk range is too short (micro edge)");
        }

        // Parametric tolerance for the arc-length measurement, capped at 1% of
        // the edge range so large-tolerance edges do not over-refine.
        let mut a_ptol_e = curve_resolution(curve.as_ref(), first, last, a_tol_e);
        let a_ptol_e_min = (last - first) / 100.0;
        if a_ptol_e > a_ptol_e_min {
            a_ptol_e = a_ptol_e_min;
        }
        let span = s2 - s1;
        let n = ((span.abs() / a_ptol_e.max(PCONFUSION)).ceil() as usize).clamp(8, 4096);
        let length = curve_arc_length(&|u| curve.d0(u), s1, s2, n);
        if length < CONFUSION {
            return self.fail("shrunk length is too short (micro edge)");
        }

        self.shrunk_first = s1;
        self.shrunk_last = s2;
        self.length = length;
        self.is_done = true;
        // The range is splittable if its length leaves room for a splitting
        // vertex (2*TolE) plus two new micro-edges (2*Confusion).
        self.is_splittable = length > 2.0 * a_tol_e + 2.0 * CONFUSION;
        Ok(())
    }

    /// Whether the shrunk range has been computed.
    pub fn is_done(&self) -> bool {
        self.is_done
    }

    /// The computed shrunk range `(t1, t2)`, or `None` when not done.
    pub fn shrunk_range(&self) -> Option<(f64, f64)> {
        if self.is_done {
            Some((self.shrunk_first, self.shrunk_last))
        } else {
            None
        }
    }

    /// The edge this shrunk range was computed for.
    pub fn edge(&self) -> Option<&Edge> {
        self.edge.as_ref()
    }

    /// The face this shrunk range was computed for.
    pub fn face(&self) -> Option<&Face> {
        self.face.as_ref()
    }

    /// The arc length of the shrunk range, when computed.
    pub fn length(&self) -> f64 {
        self.length
    }

    /// Whether the shrunk range is long enough that the edge can be split at
    /// least once.
    pub fn is_splittable(&self) -> bool {
        self.is_splittable
    }

    /// The error message of the last failed computation, if any.
    pub fn error(&self) -> Option<&str> {
        self.error.as_deref()
    }

    fn fail(&mut self, msg: &str) -> Result<(), String> {
        self.is_done = false;
        self.is_splittable = false;
        self.error = Some(msg.to_string());
        Err(msg.to_string())
    }
}

impl Default for ShrunkRange {
    fn default() -> Self {
        Self::new()
    }
}

/// Port of `BRepLib::FindValidRange` (the `Adaptor3d_Curve` overload).
///
/// Finds the range of `curve` over `[first, last]` that is not covered by the
/// tolerance spheres of the two endpoint points `p1`/`p2`. Returns `None` when
/// no such range exists (the spheres overlap / cover the whole curve).
fn find_valid_range(
    curve: &dyn Curve,
    first: f64,
    last: f64,
    p1: &GpPnt,
    tol1: f64,
    p2: &GpPnt,
    tol2: f64,
) -> Option<(f64, f64)> {
    if last - first < PCONFUSION {
        return None;
    }
    // The bisection target: curve resolution scaled down, floored at the
    // parametric confusion and the relative machine epsilon.
    let max_par = first.abs().max(last.abs());
    let an_eps = (curve_resolution(curve, first, last, tol1) * 0.1)
        .max(f64::EPSILON * max_par.max(1.0))
        .max(PCONFUSION);

    let s1 = find_nearest_valid_point(curve, first, last, true, p1, tol1, an_eps)?;
    if last - s1 < an_eps {
        return None;
    }
    let s2 = find_nearest_valid_point(curve, first, last, false, p2, tol2, an_eps)?;
    if s2 - first < an_eps {
        return None;
    }
    if s1 > s2 {
        // Overlapping tolerance spheres — no valid range.
        return None;
    }
    Some((s1, s2))
}

/// Port of the static `findNearestValidPoint` helper in `BRepLib_1.cxx`.
///
/// Starting from the appointed end of the curve, find the nearest parameter
/// whose curve point is *outside* the tolerance sphere of radius `tol` centred
/// at `vert_pnt`. The walk steps along the curve by a resolution-sized step and
/// then refines the exit by bisection.
///
/// ponytail: the OCCT `aD1Mag` fast-exit for Bézier/B-spline local
/// singularities is omitted — the sphere has the same scale as the step, so the
/// walk exits within a constant number of steps and the bisection dominates.
fn find_nearest_valid_point(
    curve: &dyn Curve,
    first: f64,
    last: f64,
    is_first: bool,
    vert_pnt: &GpPnt,
    tol: f64,
    eps: f64,
) -> Option<f64> {
    let (start_u, end_u) = if is_first { (first, last) } else { (last, first) };
    let sq_tol = tol * tol;
    if curve.d0(start_u).square_distance(vert_pnt) > sq_tol {
        // The vertex does not cover its end of the curve.
        return None;
    }

    let mut step = curve_resolution(curve, first, last, tol) * 1.01;
    if step < eps {
        step = eps;
    }
    if !is_first {
        step = -step;
    }

    let mut is_out = false;
    let mut u_in = start_u;
    let mut u_out = u_in;
    while !is_out {
        u_in = u_out;
        u_out += step;
        if (is_first && u_out > end_u) || (!is_first && u_out < end_u) {
            // The step overshoots the opposite bound: the whole range is inside
            // the sphere unless the bound itself is outside.
            is_out = curve.d0(end_u).square_distance(vert_pnt) > sq_tol;
            if !is_out {
                return None;
            }
            u_out = end_u;
            break;
        }
        is_out = curve.d0(u_out).square_distance(vert_pnt) > sq_tol;
    }

    // Precise solution with binary search on the exit interval.
    let mut a = u_in;
    let mut b = u_out;
    while (b - a).abs() > eps {
        let mid = (a + b) * 0.5;
        if curve.d0(mid).square_distance(vert_pnt) > sq_tol {
            b = mid;
        } else {
            a = mid;
        }
    }
    Some((a + b) * 0.5)
}

/// Approximate `BRepAdaptor_Curve::Resolution(tol)`: the parameter increment
/// over which the curve deviates from a straight chord by roughly `tol`.
///
/// The maximum derivative magnitude `|dp/du|` is estimated by central
/// differences over the curve's range; `Resolution = tol / max_speed`. For an
/// arc-length parametrised line (`max_speed = 1`) this is exactly `tol`.
fn curve_resolution(curve: &dyn Curve, first: f64, last: f64, tol: f64) -> f64 {
    let span = last - first;
    if span.abs() <= 1e-30 || tol <= 0.0 {
        return PCONFUSION;
    }
    let n = 8;
    let h = span / n as f64 * 0.5;
    let mut max_speed = 0.0f64;
    for i in 0..=n {
        let u = first + span * i as f64 / n as f64;
        let speed = curve.d0(u - h).distance(&curve.d0(u + h)) / (2.0 * h);
        if speed > max_speed {
            max_speed = speed;
        }
    }
    if max_speed <= 1e-30 {
        return tol.max(PCONFUSION);
    }
    (tol / max_speed).max(PCONFUSION)
}

// ---------------------------------------------------------------------------
// IntToolsTools — static helpers
// ---------------------------------------------------------------------------

/// Static helpers ported from `IntTools_Tools` (the subset used by the boolean
/// pipeline's range/vertex logic).
pub struct IntToolsTools;

impl IntToolsTools {
    /// Whether `param` lies inside `[first, last]` expanded by `tol`.
    ///
    /// Equivalent to the OCCT `IntTools_Tools::IsOnPave1`: inside the range,
    /// or within `tol` of either endpoint.
    pub fn is_in_range(param: f64, first: f64, last: f64, tol: f64) -> bool {
        let (lo, hi) = (first.min(last), first.max(last));
        param >= lo - tol && param <= hi + tol
    }

    /// A point ~43.2% of the way from `a` to `b` (OCCT
    /// `IntTools_Tools::IntermediatePoint`).
    ///
    /// Deliberately *not* the geometric midpoint: the asymmetry keeps
    /// sampling/intersection results away from symmetric (ambiguous) positions.
    pub fn middle_point(a: f64, b: f64) -> f64 {
        const PAR_T: f64 = 0.43213918;
        (1.0 - PAR_T) * a + PAR_T * b
    }

    /// Whether the point `pnt` on `edge` at `param` falls within `tol` of one
    /// of the edge's vertices (OCCT `IntTools_Tools::IsVertex(Edge, t)`).
    ///
    /// The point is supplied by the caller (it may already have been evaluated
    /// from `param`); the parameter is accepted for signature parity with OCCT
    /// and currently unused.
    pub fn is_vertex(edge: &Edge, _param: f64, pnt: &GpPnt, tol: f64) -> bool {
        let t = tol.max(0.0);
        let tol2 = t * t;
        let kids = edge.0.tshape.read().unwrap().children.clone();
        for h in kids {
            if h.shape_type() != ShapeType::Vertex {
                continue;
            }
            let v = Vertex(h);
            if BRepTool::vertex_point(&v).square_distance(pnt) <= tol2 {
                return true;
            }
        }
        false
    }

    /// The tolerance of a point set: the radius of the smallest sphere centred
    /// at the point-cloud centroid that contains every point.
    ///
    /// An empty or single-point set yields `0.0`.
    pub fn compute_tolerance(pnts: &[GpPnt]) -> f64 {
        let n = pnts.len();
        if n == 0 {
            return 0.0;
        }
        let cx = pnts.iter().map(|p| p.x()).sum::<f64>() / n as f64;
        let cy = pnts.iter().map(|p| p.y()).sum::<f64>() / n as f64;
        let cz = pnts.iter().map(|p| p.z()).sum::<f64>() / n as f64;
        let c = GpPnt::new(cx, cy, cz);
        pnts.iter().map(|p| p.distance(&c)).fold(0.0, f64::max)
    }
}

// ---------------------------------------------------------------------------
// IntContext
// ---------------------------------------------------------------------------

/// Lightweight port of `IntTools_Context`.
///
/// Caches a working edge/face pair and answers point-on-face / surface
/// projection queries through the face's registered surface. The face and
/// solid classifiers (FClass2d, BRepClass3d_SolidClassifier) are intentionally
/// left as `None` shells — wave B of Phase 16 fills them in.
#[derive(Debug, Clone, Default)]
pub struct IntContext {
    edge: Option<Edge>,
    face: Option<Face>,
}

impl IntContext {
    /// Create an empty context.
    pub fn new() -> Self {
        Self { edge: None, face: None }
    }

    /// Cache a working edge/face pair.
    pub fn set_edge(&mut self, edge: &Edge, face: &Face) {
        self.edge = Some(edge.clone());
        self.face = Some(face.clone());
    }

    /// The cached edge, if any.
    pub fn edge(&self) -> Option<&Edge> {
        self.edge.as_ref()
    }

    /// The cached face, if any.
    pub fn face(&self) -> Option<&Face> {
        self.face.as_ref()
    }

    /// Whether the point `p` lies on `face`: its surface projection is within
    /// `tol`, and the projection is inside the face's boundary.
    ///
    /// The closest `(u, v)` parameters are written to `uv` regardless of the
    /// outcome. Faces without a registered surface reject every point.
    pub fn is_point_on_face(&self, face: &Face, p: &GpPnt, uv: &mut GpPnt2d, tol: f64) -> bool {
        let Some(surf) = BRepTool::face_surface(face) else {
            return false;
        };
        let (u, v) = surface_closest_params(surf.as_ref(), p, 32, 32);
        uv.set_coord(u, v);
        if surf.d0(u, v).distance(p) > tol.max(0.0) {
            return false;
        }
        // Inside the face's 2D boundary (planar faces) — the 3D check is
        // self-consistent with the reconstructed face plane used by the polygon.
        crate::inttools::point_on_face(face, p, tol.max(0.0))
    }

    /// Closest `(u, v)` parameters of `p` on the face's surface
    /// (`GeomAPI_ProjectPointOnSurf` equivalent).
    pub fn project_point_on_face(&self, face: &Face, p: &GpPnt) -> Option<(f64, f64)> {
        let surf = BRepTool::face_surface(face)?;
        Some(surface_closest_params(surf.as_ref(), p, 32, 32))
    }

    /// Face classifier — FClass2d arrives in wave B; a shell for now.
    pub fn face_classifier(&self) -> Option<()> {
        None
    }

    /// Solid classifier — BRepClass3d_SolidClassifier arrives in wave B; a
    /// shell for now.
    pub fn solid_classifier(&self) -> Option<()> {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;

    use occt_core::gp::{GpDir, GpLin};
    use occt_geom::GeomLine;

    use crate::brep_extrema::test_box::unit_box;
    use crate::builder::TopoBuilder;
    use crate::shape::{Face, TopoShape};
    use crate::tgeometry::{EdgeGeom, GeometryRegistry};

    fn clear_tree(s: &TopoShape) {
        GeometryRegistry::global().clear_shape(s);
        let children = s.tshape.read().unwrap().children.clone();
        for c in children {
            clear_tree(&c);
        }
    }

    fn dir(x: f64, y: f64, z: f64) -> GpDir {
        GpDir::new(x, y, z).expect("dir")
    }

    // -----------------------------------------------------------------------
    // ShrunkRange
    // -----------------------------------------------------------------------

    #[test]
    fn shrunk_range_unit_box_edge_shrinks_interior() {
        let bx = unit_box();
        // Edge 0 runs from (0,0,0) to (1,0,0) — arc length 1.
        let edge = &bx.edges[0];
        let face = &bx.faces[0];
        let mut sr = ShrunkRange::new();
        sr.set_shrunk_range(edge, face, 1e-7).expect("shrunk range computes");
        assert!(sr.is_done());
        let (s1, s2) = sr.shrunk_range().expect("range");
        assert!(s1 > 0.0, "first shrunk parameter {s1} must be > 0");
        assert!(s2 < 1.0, "last shrunk parameter {s2} must be < 1");
        assert!(s2 > s1, "shrunk range must be non-empty: {s1}..{s2}");
        assert!((sr.length() - 1.0).abs() < 1e-3, "length={}", sr.length());
        assert!(sr.is_splittable());
        assert!(sr.edge().is_some() && sr.face().is_some());
        clear_tree(&bx.solid.0);
    }

    #[test]
    fn shrunk_range_missing_curve_is_err() {
        // A bare edge with no registered geometry is treated as un-shrinkable.
        let e = Edge::new();
        let f = Face::new();
        let mut sr = ShrunkRange::new();
        assert!(sr.set_shrunk_range(&e, &f, 1e-7).is_err());
        assert!(!sr.is_done());
        assert!(sr.shrunk_range().is_none());
        assert!(sr.error().is_some());
    }

    #[test]
    fn shrunk_range_degenerated_edge_is_err() {
        let b = TopoBuilder::new();
        let curve = Arc::new(GeomLine::new(GpLin::from_pnt_dir(GpPnt::zero(), dir(1.0, 0.0, 0.0))));
        let e = b.make_edge(curve, 0.0, 1.0);
        GeometryRegistry::global().set_edge(
            &e.0,
            EdgeGeom {
                curve: BRepTool::edge_curve(&e).expect("curve"),
                first: 0.0,
                last: 1.0,
                tolerance: 0.0,
                same_parameter: true,
                same_range: true,
                degenerated: true,
                pcurves: std::collections::HashMap::new(),
            },
        );
        let mut sr = ShrunkRange::new();
        assert!(sr.set_shrunk_range(&e, &Face::new(), 1e-7).is_err());
        assert!(sr.error().unwrap().contains("degenerated"));
        clear_tree(&e.0);
    }

    #[test]
    fn shrunk_range_micro_edge_is_err() {
        // Vertex spheres (2e-7 each) exceed the edge length — no valid range.
        let b = TopoBuilder::new();
        let curve = Arc::new(GeomLine::new(GpLin::from_pnt_dir(GpPnt::zero(), dir(1.0, 0.0, 0.0))));
        let e = b.make_edge(curve, 0.0, 1e-8);
        let mut sr = ShrunkRange::new();
        assert!(sr.set_shrunk_range(&e, &Face::new(), 1e-7).is_err());
        clear_tree(&e.0);
    }

    // -----------------------------------------------------------------------
    // IntToolsTools
    // -----------------------------------------------------------------------

    #[test]
    fn is_in_range_expands_by_tolerance() {
        assert!(IntToolsTools::is_in_range(0.5, 0.0, 1.0, 1e-9));
        assert!(IntToolsTools::is_in_range(0.0, 0.0, 1.0, 0.0));
        assert!(IntToolsTools::is_in_range(1.0, 0.0, 1.0, 0.0));
        assert!(!IntToolsTools::is_in_range(1.5, 0.0, 1.0, 1e-9));
        // Expanded by tol: slightly-outside parameters are accepted.
        assert!(IntToolsTools::is_in_range(1.0 + 1e-6, 0.0, 1.0, 1e-6));
        assert!(!IntToolsTools::is_in_range(1.0 + 2e-6, 0.0, 1.0, 1e-6));
        // Reversed ranges are handled.
        assert!(IntToolsTools::is_in_range(0.5, 1.0, 0.0, 0.0));
    }

    #[test]
    fn middle_point_matches_occt_intermediate_point() {
        let m = IntToolsTools::middle_point(0.0, 1.0);
        let expected = 0.43213918;
        assert!((m - expected).abs() < 1e-12, "middle={m}, expected {expected}");
        // It is deliberately not the geometric midpoint.
        assert!((m - 0.5).abs() > 1e-3);
        let m2 = IntToolsTools::middle_point(10.0, 20.0);
        assert!((m2 - (10.0 + 0.43213918 * 10.0)).abs() < 1e-12);
    }

    #[test]
    fn is_vertex_detects_endpoint_hits() {
        let bx = unit_box();
        let edge = &bx.edges[0]; // (0,0,0) -> (1,0,0), params [0, 1]
        assert!(IntToolsTools::is_vertex(edge, 0.0, &GpPnt::new(0.0, 0.0, 0.0), 1e-6));
        assert!(IntToolsTools::is_vertex(edge, 1.0, &GpPnt::new(1.0, 0.0, 0.0), 1e-6));
        // A point slightly off the corner but within tolerance counts.
        assert!(IntToolsTools::is_vertex(edge, 0.0, &GpPnt::new(1e-7, 0.0, 0.0), 1e-6));
        // The interior point is not a vertex.
        assert!(!IntToolsTools::is_vertex(edge, 0.5, &GpPnt::new(0.5, 0.0, 0.0), 1e-6));
        clear_tree(&bx.solid.0);
    }

    #[test]
    fn compute_tolerance_centroid_radius() {
        assert_eq!(IntToolsTools::compute_tolerance(&[]), 0.0);
        assert_eq!(IntToolsTools::compute_tolerance(&[GpPnt::new(1.0, 2.0, 3.0)]), 0.0);
        let pts = vec![GpPnt::new(0.0, 0.0, 0.0), GpPnt::new(1.0, 0.0, 0.0)];
        assert!((IntToolsTools::compute_tolerance(&pts) - 0.5).abs() < 1e-12);
        // Centroid of {0, 2} is 1; radius = 1.
        let pts2 = vec![GpPnt::new(0.0, 0.0, 0.0), GpPnt::new(0.0, 0.0, 2.0)];
        assert!((IntToolsTools::compute_tolerance(&pts2) - 1.0).abs() < 1e-12);
    }

    // -----------------------------------------------------------------------
    // IntContext
    // -----------------------------------------------------------------------

    #[test]
    fn context_caches_edge_face() {
        let mut ctx = IntContext::new();
        assert!(ctx.edge().is_none() && ctx.face().is_none());
        let bx = unit_box();
        ctx.set_edge(&bx.edges[0], &bx.faces[0]);
        assert!(ctx.edge().is_some() && ctx.face().is_some());
        // Classifier shells.
        assert!(ctx.face_classifier().is_none());
        assert!(ctx.solid_classifier().is_none());
        clear_tree(&bx.solid.0);
    }

    #[test]
    fn context_point_on_face_inside_outside() {
        let bx = unit_box();
        let face = &bx.faces[0]; // bottom face (z = 0)
        let ctx = IntContext::new();
        let mut uv = GpPnt2d::zero();

        // Inside the face, on the surface.
        assert!(ctx.is_point_on_face(face, &GpPnt::new(0.25, 0.75, 0.0), &mut uv, 1e-9));
        // Projected uv evaluates back to the point.
        let surf = BRepTool::face_surface(face).expect("surface");
        let q = surf.d0(uv.x(), uv.y());
        assert!(q.distance(&GpPnt::new(0.25, 0.75, 0.0)) < 1e-9);

        // In the face's plane but outside its boundary.
        assert!(!ctx.is_point_on_face(face, &GpPnt::new(-0.5, 0.5, 0.0), &mut uv, 1e-9));
        // Above the plane (distance exceeds tolerance).
        assert!(!ctx.is_point_on_face(face, &GpPnt::new(0.5, 0.5, 1.0), &mut uv, 1e-9));
        // On the surface but far outside the search window.
        assert!(!ctx.is_point_on_face(face, &GpPnt::new(5.0, 0.5, 0.0), &mut uv, 1e-9));
        clear_tree(&bx.solid.0);
    }

    #[test]
    fn context_project_point_on_face_uv() {
        let bx = unit_box();
        let face = &bx.faces[0]; // bottom face (z = 0)
        let ctx = IntContext::new();
        let surf = BRepTool::face_surface(face).expect("surface");

        let p = GpPnt::new(0.25, 0.75, 0.0);
        let (u, v) = ctx.project_point_on_face(face, &p).expect("projection");
        let q = surf.d0(u, v);
        assert!(q.distance(&p) < 1e-9, "projection {u},{v} maps to {q:?}, want {p:?}");

        // A point offset from the surface projects onto its closest uv.
        let p2 = GpPnt::new(0.25, 0.75, 0.5);
        let (u2, v2) = ctx.project_point_on_face(face, &p2).expect("projection");
        let q2 = surf.d0(u2, v2);
        assert!((q2.distance(&p2) - 0.5).abs() < 1e-9, "dist={}", q2.distance(&p2));
        clear_tree(&bx.solid.0);
    }
}
