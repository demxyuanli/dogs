//! Full port of `IntTools_Context` (TKBO) — a cached intersection context.
//!
//! The context bundles the geometry/topology toolkit used by the boolean
//! intersection pipeline and caches the reusable tools so repeated queries on
//! the same shape do not rebuild them:
//!
//! - a per-face 2D classifier ([`FClass2d`]) cache, so point-in-face queries on
//!   the same face sample the UV boundary only once;
//! - point/vertex projection helpers onto edges and surfaces;
//! - the block/point validity predicates used to decide whether an
//!   intersection sub-range lies inside a face.
//!
//! Source: `IntTools_Context.hxx/.cxx` (TKBO, `ModelingAlgorithms/TKBO/IntTools`).
//! This is the *complete* context; the lighter [`crate::inttools_range::IntContext`]
//! (Phase 16) is a shell that fills the classifiers in later.
//!
//! The `is_valid_block_for_*` predicates are self-contained: with no
//! intersection curve handle in [`IntRange`], the 1D range is interpreted over
//! the face's first (`u`) parameter at the middle of its `v` range (the
//! u-midline of the face's UV domain).

use std::collections::HashMap;

use occt_core::gp::{GpPnt, GpPnt2d};
use occt_core::precision::{CONFUSION, PCONFUSION};
use occt_geom::Curve;

use crate::brep_surface::{face_uv_bounds, surface_closest_params};
use crate::brep_tool::BRepTool;
use crate::fclass2d::{FaceState, FClass2d};
use crate::inttools_data::IntRange;
use crate::inttools_roots;
use crate::shape::{Edge, Face, Vertex};
use crate::shape_naming::ShapeId;

/// Port of `IntTools_Context`.
///
/// Owns the cached per-face 2D classifiers and answers point-in-face,
/// projection and validity queries. The OCCT class also caches surface/curve
/// projectors, solid classifiers and hatchers; this port caches the
/// [`FClass2d`] classifiers (the hot path) and computes the projections on
/// demand.
#[derive(Debug, Clone, Default)]
pub struct IntToolsContext {
    /// Per-face 2D classifiers, keyed by the face's stable [`ShapeId`]
    /// (`myFClass2dMap`).
    fclass2d_cache: HashMap<ShapeId, FClass2d>,
    /// Point-on-surface projection tolerance (`myPOnSTolerance`).
    pon_s_tolerance: f64,
}

impl IntToolsContext {
    /// Create an empty context.
    pub fn new() -> Self {
        Self {
            fclass2d_cache: HashMap::new(),
            pon_s_tolerance: 1e-12,
        }
    }

    /// Set the point-on-surface projection tolerance and drop cached surface
    /// projectors so the new value takes effect on the next query
    /// (`SetPOnSProjectionTolerance`).
    ///
    /// This port has no cached surface projectors to clear — the tolerance is
    /// stored for parity and future use.
    pub fn set_pon_s_projection_tolerance(&mut self, value: f64) {
        self.pon_s_tolerance = value;
    }

    /// The point-on-surface projection tolerance.
    pub fn pon_s_projection_tolerance(&self) -> f64 {
        self.pon_s_tolerance
    }

    // -----------------------------------------------------------------------
    // 2D face classifier
    // -----------------------------------------------------------------------

    /// The per-face 2D classifier for `face`, built and cached on first use
    /// (`IntTools_Context::FClass2d`).
    ///
    /// The classifier's UV tolerance is the face tolerance floored at
    /// `Precision::PConfusion`, mirroring OCCT (the classifier is built from a
    /// forward-oriented copy of the face and its `BRep_Tool::Tolerance`).
    pub fn fclass2d(&mut self, face: &Face) -> Result<FClass2d, String> {
        let id = ShapeId::of(&face.0);
        if let Some(cl) = self.fclass2d_cache.get(&id) {
            return Ok(cl.clone());
        }
        let tol = BRepTool::face_tolerance(face).max(PCONFUSION);
        let cl = FClass2d::new(face, tol)?;
        self.fclass2d_cache.insert(id, cl.clone());
        Ok(cl)
    }

    /// Number of cached 2D classifiers (diagnostics/tests).
    pub fn fclass2d_cache_len(&self) -> usize {
        self.fclass2d_cache.len()
    }

    /// Drop all cached tools (`IntTools_Context::~IntTools_Context` release of
    /// the classifier map; [`clear_cached`](Self::clear_cached)).
    pub fn clear_cached(&mut self) {
        self.fclass2d_cache.clear();
    }

    // -----------------------------------------------------------------------
    // Point-in-face classification
    // -----------------------------------------------------------------------

    /// State of the 2D point `uv` relative to `face` (`StatePointFace`).
    ///
    /// Delegates to the cached [`FClass2d::perform`]. When the classifier has
    /// not been built yet it is created with `tol` (floored at
    /// `Precision::PConfusion`).
    pub fn state_point_face(&mut self, face: &Face, uv: (f64, f64), tol: f64) -> Result<FaceState, String> {
        let id = ShapeId::of(&face.0);
        if !self.fclass2d_cache.contains_key(&id) {
            let cl = FClass2d::new(face, tol.max(PCONFUSION))?;
            self.fclass2d_cache.insert(id, cl);
        }
        Ok(self.fclass2d_cache[&id].perform(GpPnt2d::new(uv.0, uv.1)))
    }

    /// Whether the 3D point `p` lies *strictly inside* the face `face`
    /// (`IntTools_Context::IsPointInFace`).
    ///
    /// The point is projected onto the face surface (unless `uv` is supplied),
    /// the 3D distance must be within `tol`, and the projection's UV must
    /// classify `In` (a boundary `On` or exterior `Out` point is rejected).
    pub fn is_point_in_face(&mut self, face: &Face, p: &GpPnt, uv: Option<(f64, f64)>, tol: f64) -> Result<bool, String> {
        Ok(self.point_face_state(face, p, uv, tol)? == Some(FaceState::In))
    }

    /// Whether the 3D point `p` lies in *or on* the face `face`
    /// (`IntTools_Context::IsPointInOnFace`).
    ///
    /// Like [`is_point_in_face`](Self::is_point_in_face) but the projection's
    /// UV need only be not-`Out`.
    pub fn is_point_in_on_face(&mut self, face: &Face, p: &GpPnt, uv: Option<(f64, f64)>, tol: f64) -> Result<bool, String> {
        match self.point_face_state(face, p, uv, tol)? {
            Some(s) => Ok(s != FaceState::Out),
            None => Ok(false),
        }
    }

    /// Whether the 2D point `uv` is inside (in or on) the face `face`
    /// (`IntTools_Context::IsValidPointForFace` — the projection-UV check).
    pub fn is_valid_point_for_face(&mut self, uv: (f64, f64), face: &Face) -> Result<bool, String> {
        let cl = self.fclass2d(face)?;
        Ok(cl.perform(GpPnt2d::new(uv.0, uv.1)) != FaceState::Out)
    }

    /// Whether the 2D point `uv` is inside both faces `f1` and `f2`
    /// (`IntTools_Context::IsValidPointForFaces`).
    pub fn is_valid_point_for_faces(&mut self, uv: (f64, f64), f1: &Face, f2: &Face) -> Result<bool, String> {
        Ok(self.is_valid_point_for_face(uv, f1)? && self.is_valid_point_for_face(uv, f2)?)
    }

    /// Whether every point of the 1D parameter block `range` is inside the face
    /// `face` (`IntTools_Context::IsValidBlockForFace`).
    ///
    /// The midpoint and both endpoints of `range` are classified; all three
    /// must be not-`Out`. The range is interpreted over the face's `u`
    /// parameter at the middle of its `v` range (see the module docs).
    pub fn is_valid_block_for_face(&mut self, range: IntRange, face: &Face) -> Result<bool, String> {
        let cl = self.fclass2d(face)?;
        let (_, _, vmin, vmax) = face_block_bounds(&cl, face);
        if !vmin.is_finite() || !vmax.is_finite() {
            // No boundary restriction (open/closed periodic face): the
            // classifier treats every point as In.
            return Ok(true);
        }
        let v_mid = 0.5 * (vmin + vmax);
        let samples = [range.first, 0.5 * (range.first + range.last), range.last];
        for u in samples {
            if cl.perform(GpPnt2d::new(u, v_mid)) == FaceState::Out {
                return Ok(false);
            }
        }
        Ok(true)
    }

    /// Whether every point of the 1D parameter block `range` is inside both
    /// faces `f1` and `f2` (`IntTools_Context::IsValidBlockForFaces`).
    pub fn is_valid_block_for_faces(&mut self, range: IntRange, f1: &Face, f2: &Face) -> Result<bool, String> {
        Ok(self.is_valid_block_for_face(range, f1)? && self.is_valid_block_for_face(range, f2)?)
    }

    // -----------------------------------------------------------------------
    // Projection helpers
    // -----------------------------------------------------------------------

    /// Closest `(u, v)` parameters of `p` on the face's surface
    /// (`IntTools_Context::ProjPS` + `LowerDistanceParameters`).
    pub fn project_point_on_face(&self, face: &Face, p: &GpPnt) -> Result<(f64, f64), String> {
        let Some(surf) = BRepTool::face_surface(face) else {
            return Err("IntToolsContext::project_point_on_face: face has no surface".into());
        };
        Ok(surface_closest_params(surf.as_ref(), p, 32, 32))
    }

    /// Parameter of the closest point of `p` on the edge's 3D curve
    /// (`IntTools_Context::ProjectPointOnEdge`).
    ///
    /// `None` when the edge has no registered curve or a degenerate/unbounded
    /// parameter range.
    pub fn project_point_on_edge(&self, edge: &Edge, p: &GpPnt) -> Option<f64> {
        let curve = BRepTool::edge_curve(edge)?;
        let (a, b) = BRepTool::edge_parameters(edge);
        if !a.is_finite() || !b.is_finite() || b - a <= PCONFUSION {
            return None;
        }
        Some(inttools_roots::parameter(&|u| curve.d0(u), p, a, b))
    }

    // -----------------------------------------------------------------------
    // Geometric classification
    // -----------------------------------------------------------------------

    /// Classify the vertex `vertex` against the edge `edge`
    /// (`IntTools_Context::ComputePE`).
    ///
    /// Returns `0` when the vertex point projects onto the edge within
    /// `tol + edge tolerance + Confusion`, and a negative error code otherwise:
    /// `-1` degenerated edge, `-2` edge has no curve, `-3` projection failed,
    /// `-4` vertex is farther than the tolerance sum.
    pub fn compute_pe(&self, vertex: &Vertex, edge: &Edge, tol: f64) -> i32 {
        if BRepTool::is_degenerated(edge) {
            return -1;
        }
        if BRepTool::edge_curve(edge).is_none() {
            return -2;
        }
        let p = BRepTool::vertex_point(vertex);
        let Some(t) = self.project_point_on_edge(edge, &p) else {
            return -3;
        };
        let Some(curve) = BRepTool::edge_curve(edge) else {
            return -3;
        };
        let dist = p.distance(&curve.d0(t));
        let tol_sum = tol.max(0.0) + BRepTool::edge_tolerance(edge) + CONFUSION;
        if dist > tol_sum {
            return -4;
        }
        0
    }

    /// Classify the two edges `edge1`/`edge2` against each other
    /// (simplified `IntTools_EdgeEdge`).
    ///
    /// Returns `0` when the edges touch/intersect (minimum distance within
    /// `tol + both edge tolerances + Confusion`), `-4` when they are separated,
    /// and negative error codes for degenerate input.
    ///
    /// ponytail: the full `IntTools_EdgeEdge` root-finding is not run here — a
    /// sampled closest-approach distance is sufficient for the "basic
    /// classification" the boolean pipeline needs; upgrade to `EdgeEdge` if
    /// exact common parts are required.
    pub fn compute_ve(&self, edge1: &Edge, edge2: &Edge, tol: f64) -> i32 {
        let Some(c1) = BRepTool::edge_curve(edge1) else { return -2; };
        let Some(c2) = BRepTool::edge_curve(edge2) else { return -2; };
        let (a1, b1) = BRepTool::edge_parameters(edge1);
        let (a2, b2) = BRepTool::edge_parameters(edge2);
        if !a1.is_finite() || !b1.is_finite() || !a2.is_finite() || !b2.is_finite() {
            return -3;
        }
        let dist = min_curve_distance(c1.as_ref(), a1, b1, c2.as_ref(), a2, b2);
        let tol_sum = tol.max(0.0)
            + BRepTool::edge_tolerance(edge1)
            + BRepTool::edge_tolerance(edge2)
            + CONFUSION;
        if dist > tol_sum {
            return -4;
        }
        0
    }

    /// Classify the vertex `vertex` against the face `face`
    /// (`IntTools_Context::ComputeVF`).
    ///
    /// Returns `0` when the vertex projects onto the surface within
    /// `vertex tol + face tol + max(tol, Confusion)` and the projection lies
    /// strictly inside the face; `-1` projection failed, `-2` distance too
    /// large, `-3` projection is out of or on the face boundary.
    pub fn compute_vf(&mut self, vertex: &Vertex, face: &Face, tol: f64) -> i32 {
        let p = BRepTool::vertex_point(vertex);
        let Some(surf) = BRepTool::face_surface(face) else {
            return -1;
        };
        let (u, v) = surface_closest_params(surf.as_ref(), &p, 32, 32);
        let dist = surf.d0(u, v).distance(&p);
        let tol_sum = BRepTool::vertex_tolerance(vertex)
            + BRepTool::face_tolerance(face)
            + tol.max(CONFUSION);
        if dist > tol_sum {
            return -2;
        }
        let state = match self.state_point_face(face, (u, v), tol) {
            Ok(s) => s,
            Err(_) => return -1,
        };
        if state != FaceState::In {
            return -3;
        }
        0
    }

    /// Whether the vertex `vertex` lies on the curve `curve` within tolerance
    /// (`IntTools_Context::IsVertexOnLine`).
    ///
    /// The vertex point is projected onto the curve's natural parameter range
    /// and the distance compared to `2 * (vertex tol + tol)` floored at `1e-6`,
    /// mirroring the OCCT `aTolSum` policy.
    pub fn is_vertex_on_line(&self, vertex: &Vertex, curve: &dyn Curve, tol: f64) -> bool {
        let p = BRepTool::vertex_point(vertex);
        let Some((_, dist)) = project_point_on_curve(&p, curve) else {
            return false;
        };
        let tol_sum = (2.0 * (BRepTool::vertex_tolerance(vertex) + tol.max(0.0))).max(1e-6);
        dist <= tol_sum
    }

    // -----------------------------------------------------------------------
    // Misc
    // -----------------------------------------------------------------------

    /// UV bounds of the face's surface (`IntTools_Context::UVBounds`).
    ///
    /// For an unbounded surface (e.g. a plane) the bounds are infinite.
    pub fn uv_bounds(&self, face: &Face) -> (f64, f64, f64, f64) {
        face_uv_bounds(face)
    }

    // -----------------------------------------------------------------------
    // Helpers
    // -----------------------------------------------------------------------

    /// Surface-projection + distance + 2D-classification of `p` against `face`.
    ///
    /// Returns `Ok(None)` when the 3D distance from `p` to the surface at the
    /// resolved `(u, v)` exceeds `tol`; `Ok(Some(state))` otherwise.
    fn point_face_state(
        &mut self,
        face: &Face,
        p: &GpPnt,
        uv: Option<(f64, f64)>,
        tol: f64,
    ) -> Result<Option<FaceState>, String> {
        let Some(surf) = BRepTool::face_surface(face) else {
            return Err("IntToolsContext::point_face_state: face has no surface".into());
        };
        let (u, v) = match uv {
            Some((u, v)) => (u, v),
            None => surface_closest_params(surf.as_ref(), p, 32, 32),
        };
        if surf.d0(u, v).distance(p) > tol.max(0.0) {
            return Ok(None);
        }
        let state = self.state_point_face(face, (u, v), tol)?;
        Ok(Some(state))
    }
}

// ---------------------------------------------------------------------------
// Free helpers
// ---------------------------------------------------------------------------

/// The finite UV bounds of the face used by the block predicates.
///
/// Prefers the finite surface UV bounds; for unbounded surfaces (planes) falls
/// back to the sampled boundary-ring bounds of the classifier.
fn face_block_bounds(cl: &FClass2d, face: &Face) -> (f64, f64, f64, f64) {
    let (u0, u1, v0, v1) = face_uv_bounds(face);
    if u0.is_finite() && u1.is_finite() && v0.is_finite() && v1.is_finite() {
        return (u0, u1, v0, v1);
    }
    ring_bounds(cl)
}

/// Bounding box of the classifier's sampled UV boundary rings.
fn ring_bounds(cl: &FClass2d) -> (f64, f64, f64, f64) {
    let mut umin = f64::INFINITY;
    let mut umax = f64::NEG_INFINITY;
    let mut vmin = f64::INFINITY;
    let mut vmax = f64::NEG_INFINITY;
    let mut any = false;
    let mut rings: Vec<&[GpPnt2d]> = Vec::new();
    if let Some(r) = cl.outer_ring() {
        rings.push(r);
    }
    rings.extend(cl.hole_rings().iter().map(|r| r.as_slice()));
    for ring in rings {
        for p in ring {
            any = true;
            umin = umin.min(p.x());
            umax = umax.max(p.x());
            vmin = vmin.min(p.y());
            vmax = vmax.max(p.y());
        }
    }
    if any {
        (umin, umax, vmin, vmax)
    } else {
        (f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY)
    }
}

/// Closest `(parameter, distance)` of `p` on `curve` over its natural range.
fn project_point_on_curve(p: &GpPnt, curve: &dyn Curve) -> Option<(f64, f64)> {
    let (a, b) = curve_param_window(curve, p);
    let t = inttools_roots::parameter(&|u| curve.d0(u), p, a, b);
    Some((t, p.distance(&curve.d0(t))))
}

/// A finite parameter window over which to project `p` onto `curve`.
///
/// Bounded curves use their natural range; unbounded curves (infinite lines)
/// project `p` onto the tangent at a reference point and return a window around
/// the estimated closest parameter (exact for straight lines).
fn curve_param_window(curve: &dyn Curve, near: &GpPnt) -> (f64, f64) {
    let a = curve.first_parameter();
    let b = curve.last_parameter();
    if a.is_finite() && b.is_finite() && b > a {
        return (a, b);
    }
    let p0 = curve.d0(0.0);
    let tan = curve.d1(0.0).1;
    let t_est = if tan.square_magnitude() > 1e-300 {
        near.coord.subtracted(&p0.coord).dot(tan.xyz()) / tan.square_magnitude()
    } else {
        0.0
    };
    let t_est = if t_est.is_finite() { t_est } else { 0.0 };
    let w = t_est.abs().max(1.0) + 1.0;
    (t_est - w, t_est + w)
}

/// Minimum distance between two curves over their parameter ranges.
///
/// Samples each curve and projects every sample onto the other curve, keeping
/// the smallest 3D distance. Exact for straight-line edges (the box-test
/// geometry); for curved edges the sampling resolves the closest approach to
/// the projection tolerance.
fn min_curve_distance(
    c1: &dyn Curve,
    a1: f64,
    b1: f64,
    c2: &dyn Curve,
    a2: f64,
    b2: f64,
) -> f64 {
    let n = 32;
    let mut min = f64::INFINITY;
    for i in 0..=n {
        let u = a1 + (b1 - a1) * (i as f64 / n as f64);
        let p = c1.d0(u);
        let t = inttools_roots::parameter(&|v| c2.d0(v), &p, a2, b2);
        min = min.min(p.distance(&c2.d0(t)));
    }
    for i in 0..=n {
        let u = a2 + (b2 - a2) * (i as f64 / n as f64);
        let p = c2.d0(u);
        let t = inttools_roots::parameter(&|v| c1.d0(v), &p, a1, b1);
        min = min.min(p.distance(&c1.d0(t)));
    }
    min
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;

    use occt_core::gp::{GpDir, GpLin};
    use occt_geom::GeomLine;

    use crate::brep_extrema::test_box::unit_box;
    use crate::builder::TopoBuilder;
    use crate::shape::TopoShape;
    use crate::tgeometry::GeometryRegistry;

    /// Drop the global geometry entries owned by a shape subtree.
    fn clear_tree(s: &TopoShape) {
        GeometryRegistry::global().clear_shape(s);
        let children = s.tshape.read().unwrap().children.clone();
        for c in children {
            clear_tree(&c);
        }
    }

    fn dir(x: f64, y: f64, z: f64) -> GpDir {
        GpDir::new(x, y, z).expect("unit direction")
    }

    // -----------------------------------------------------------------------
    // Classifier cache
    // -----------------------------------------------------------------------

    #[test]
    fn classifier_cache_reused_per_face() {
        let bx = unit_box();
        let mut ctx = IntToolsContext::new();
        assert_eq!(ctx.fclass2d_cache_len(), 0);

        // First query on the bottom face builds the classifier.
        assert!(ctx
            .is_point_in_face(&bx.faces[0], &GpPnt::new(0.25, 0.75, 0.0), None, 1e-6)
            .unwrap());
        assert_eq!(ctx.fclass2d_cache_len(), 1);

        // A second query on the same face reuses it.
        assert!(ctx
            .is_point_in_face(&bx.faces[0], &GpPnt::new(0.5, 0.5, 0.0), None, 1e-6)
            .unwrap());
        assert_eq!(ctx.fclass2d_cache_len(), 1);

        // A different face builds a second classifier.
        let _ = ctx
            .is_point_in_face(&bx.faces[1], &GpPnt::new(0.25, 0.25, 1.0), None, 1e-6)
            .unwrap();
        assert_eq!(ctx.fclass2d_cache_len(), 2);

        ctx.clear_cached();
        assert_eq!(ctx.fclass2d_cache_len(), 0);
        clear_tree(&bx.solid.0);
    }

    // -----------------------------------------------------------------------
    // Point-in-face
    // -----------------------------------------------------------------------

    #[test]
    fn is_point_in_face_inside_outside_boundary() {
        let bx = unit_box();
        let face = &bx.faces[0]; // bottom (z = 0)
        let mut ctx = IntToolsContext::new();
        let tol = 1e-6;

        // Face centre -> strictly In.
        assert!(ctx.is_point_in_face(face, &GpPnt::new(0.5, 0.5, 0.0), None, tol).unwrap());
        // In the plane but outside the boundary -> false.
        assert!(!ctx.is_point_in_face(face, &GpPnt::new(-0.5, 0.5, 0.0), None, tol).unwrap());
        // Above the surface (distance > tol) -> false.
        assert!(!ctx.is_point_in_face(face, &GpPnt::new(0.5, 0.5, 1.0), None, tol).unwrap());

        // A boundary point: strict In is false, in-or-on is true, state is On.
        let bnd = GpPnt::new(0.5, 0.0, 0.0);
        assert!(!ctx.is_point_in_face(face, &bnd, None, tol).unwrap());
        assert!(ctx.is_point_in_on_face(face, &bnd, None, tol).unwrap());
        let (u, v) = ctx.project_point_on_face(face, &bnd).unwrap();
        assert_eq!(ctx.state_point_face(face, (u, v), tol).unwrap(), FaceState::On);
        clear_tree(&bx.solid.0);
    }

    #[test]
    fn is_point_in_face_explicit_uv_skips_projection() {
        let bx = unit_box();
        let face = &bx.faces[0];
        let mut ctx = IntToolsContext::new();
        // Supplying the exact UV (0.5, 0.5) classifies In without projecting.
        let p = GpPnt::new(0.5, 0.5, 0.0);
        assert!(ctx.is_point_in_face(face, &p, Some((0.5, 0.5)), 1e-9).unwrap());
        // A supplied UV that is far from the 3D point fails the distance check.
        assert!(!ctx.is_point_in_face(face, &p, Some((0.5, -0.5)), 1e-9).unwrap());
        clear_tree(&bx.solid.0);
    }

    #[test]
    fn state_point_face_in_on_out() {
        let bx = unit_box();
        let face = &bx.faces[0];
        let mut ctx = IntToolsContext::new();
        let tol = 1e-6;
        assert_eq!(ctx.state_point_face(face, (0.5, 0.5), tol).unwrap(), FaceState::In);
        assert_eq!(ctx.state_point_face(face, (0.5, 0.0), tol).unwrap(), FaceState::On);
        assert_eq!(ctx.state_point_face(face, (-0.5, 0.5), tol).unwrap(), FaceState::Out);
        clear_tree(&bx.solid.0);
    }

    #[test]
    fn is_valid_point_for_face_in_on_out() {
        let bx = unit_box();
        let face = &bx.faces[0];
        let mut ctx = IntToolsContext::new();
        assert!(ctx.is_valid_point_for_face((0.5, 0.5), face).unwrap());
        // Boundary counts as valid (in-or-on).
        assert!(ctx.is_valid_point_for_face((0.5, 0.0), face).unwrap());
        assert!(!ctx.is_valid_point_for_face((-0.5, 0.5), face).unwrap());
        assert!(ctx
            .is_valid_point_for_faces((0.5, 0.5), &bx.faces[0], &bx.faces[1])
            .unwrap());
        clear_tree(&bx.solid.0);
    }

    // -----------------------------------------------------------------------
    // Blocks
    // -----------------------------------------------------------------------

    #[test]
    fn is_valid_block_for_face_interior_and_exterior() {
        let bx = unit_box();
        let face = &bx.faces[0];
        let mut ctx = IntToolsContext::new();

        // A u-range inside the face's u-domain, sampled at v = 0.5.
        let inside = IntRange::new(0.2, 0.8).unwrap();
        assert!(ctx.is_valid_block_for_face(inside, face).unwrap());
        // The whole u-domain (endpoints On the boundary) is still valid.
        let whole = IntRange::new(0.0, 1.0).unwrap();
        assert!(ctx.is_valid_block_for_face(whole, face).unwrap());
        // A u-range that leaves the domain is rejected.
        let outside = IntRange::new(1.5, 2.5).unwrap();
        assert!(!ctx.is_valid_block_for_face(outside, face).unwrap());

        // Both faces: valid in both, invalid when one rejects.
        assert!(ctx
            .is_valid_block_for_faces(inside, &bx.faces[0], &bx.faces[1])
            .unwrap());
        assert!(!ctx
            .is_valid_block_for_faces(outside, &bx.faces[0], &bx.faces[1])
            .unwrap());
        clear_tree(&bx.solid.0);
    }

    // -----------------------------------------------------------------------
    // Projection
    // -----------------------------------------------------------------------

    #[test]
    fn project_point_on_edge_returns_parameter() {
        let bx = unit_box();
        let ctx = IntToolsContext::new();
        let edge = &bx.edges[0]; // (0,0,0) -> (1,0,0), params [0, 1]
        let t = ctx.project_point_on_edge(edge, &GpPnt::new(0.5, 0.0, 0.0)).expect("projection");
        assert!((t - 0.5).abs() < 1e-6, "param {t}");
        // A point off the line still projects to its closest parameter.
        let t2 = ctx.project_point_on_edge(edge, &GpPnt::new(0.25, 3.0, 0.0)).expect("projection");
        assert!((t2 - 0.25).abs() < 1e-6, "param {t2}");
        // An edge without geometry has no projection.
        assert!(ctx.project_point_on_edge(&Edge::new(), &GpPnt::zero()).is_none());
        clear_tree(&bx.solid.0);
    }

    #[test]
    fn project_point_on_face_returns_uv() {
        let bx = unit_box();
        let ctx = IntToolsContext::new();
        let face = &bx.faces[0];
        let (u, v) = ctx.project_point_on_face(face, &GpPnt::new(0.25, 0.75, 0.0)).unwrap();
        let surf = BRepTool::face_surface(face).expect("surface");
        let q = surf.d0(u, v);
        assert!(q.distance(&GpPnt::new(0.25, 0.75, 0.0)) < 1e-6, "uv {u},{v}");
        clear_tree(&bx.solid.0);
    }

    // -----------------------------------------------------------------------
    // Geometric classification
    // -----------------------------------------------------------------------

    #[test]
    fn compute_pe_classifies_vertex_on_edge() {
        let bx = unit_box();
        let ctx = IntToolsContext::new();
        // Vertex 0 (0,0,0) lies on edge 0 ((0,0,0) -> (1,0,0)).
        assert_eq!(ctx.compute_pe(&bx.vertices[0], &bx.edges[0], 1e-7), 0);
        // Vertex 6 (1,1,1) is farther than the tolerance from edge 0.
        assert_eq!(ctx.compute_pe(&bx.vertices[6], &bx.edges[0], 1e-7), -4);
        // A bare edge (no geometry) -> -2.
        assert_eq!(ctx.compute_pe(&bx.vertices[0], &Edge::new(), 1e-7), -2);
        clear_tree(&bx.solid.0);
    }

    #[test]
    fn compute_ve_classifies_edges() {
        let bx = unit_box();
        let ctx = IntToolsContext::new();
        // edge 0 and edge 1 share the vertex (1,0,0) -> touch.
        assert_eq!(ctx.compute_ve(&bx.edges[0], &bx.edges[1], 1e-7), 0);
        // edge 0 and edge 4 are parallel disjoint (z=0 vs z=1) -> separated.
        assert_eq!(ctx.compute_ve(&bx.edges[0], &bx.edges[4], 1e-7), -4);
        clear_tree(&bx.solid.0);
    }

    #[test]
    fn compute_vf_classifies_vertex_on_face() {
        let bx = unit_box();
        let b = TopoBuilder::new();
        let interior = b.make_vertex(GpPnt::new(0.5, 0.5, 0.0), 1e-7);
        let above = b.make_vertex(GpPnt::new(0.5, 0.5, 2.0), 1e-7);
        let mut ctx = IntToolsContext::new();

        // Projection strictly inside the bottom face -> on.
        assert_eq!(ctx.compute_vf(&interior, &bx.faces[0], 1e-7), 0);
        // Distance too large -> -2.
        assert_eq!(ctx.compute_vf(&above, &bx.faces[0], 1e-7), -2);
        // Box corner vertex projects onto the face boundary -> -3.
        assert_eq!(ctx.compute_vf(&bx.vertices[0], &bx.faces[0], 1e-7), -3);
        clear_tree(&bx.solid.0);
        clear_tree(&interior.0);
        clear_tree(&above.0);
    }

    #[test]
    fn is_vertex_on_line_detects_hits() {
        let bx = unit_box();
        // A line along edge 0: through (0,0,0) direction (1,0,0).
        let lin = GpLin::from_pnt_dir(GpPnt::zero(), dir(1.0, 0.0, 0.0));
        let curve: Arc<dyn Curve> = Arc::new(GeomLine::new(lin));
        let ctx = IntToolsContext::new();
        // Vertex 0 is on the line.
        assert!(ctx.is_vertex_on_line(&bx.vertices[0], curve.as_ref(), 1e-7));
        // Vertex 6 (1,1,1) is off the line.
        assert!(!ctx.is_vertex_on_line(&bx.vertices[6], curve.as_ref(), 1e-7));
        clear_tree(&bx.solid.0);
    }

    // -----------------------------------------------------------------------
    // Misc
    // -----------------------------------------------------------------------

    #[test]
    fn uv_bounds_delegates_to_surface() {
        let bx = unit_box();
        let ctx = IntToolsContext::new();
        let (u0, u1, v0, v1) = ctx.uv_bounds(&bx.faces[0]);
        // The bottom face is an (unbounded) plane.
        assert!(!u0.is_finite() && !u1.is_finite() && !v0.is_finite() && !v1.is_finite());
        let (a0, a1, b0, b1) = crate::brep_surface::face_uv_bounds(&bx.faces[0]);
        assert_eq!((u0, u1, v0, v1), (a0, a1, b0, b1));
        clear_tree(&bx.solid.0);
    }

    #[test]
    fn projection_tolerance_setter() {
        let mut ctx = IntToolsContext::new();
        assert_eq!(ctx.pon_s_projection_tolerance(), 1e-12);
        ctx.set_pon_s_projection_tolerance(1e-8);
        assert_eq!(ctx.pon_s_projection_tolerance(), 1e-8);
    }
}
