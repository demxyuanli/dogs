use super::prelude::*;
use super::*;

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
    pub(super) fclass2d_cache: HashMap<ShapeId, FClass2d>,
    /// Point-on-surface projection tolerance (`myPOnSTolerance`).
    pub(super) pon_s_tolerance: f64,
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

    /// Drop the cached classifier so the next [`fclass2d`] rebuilds it
    /// (`IntTools_FClass2d::Init` after holes are added).
    pub fn fclass2d_invalidate(&mut self, face: &Face) {
        self.fclass2d_cache.remove(&ShapeId::of(&face.0));
    }

    /// `IntTools_Context::IsInfiniteFace` — true when the 3D box of the face
    /// is open on any side.
    pub fn is_infinite_face(&self, face: &Face) -> bool {
        let a_box = crate::bbox_from_geometry::shape_bbox(&face.0);
        a_box.is_open_xmin()
            || a_box.is_open_xmax()
            || a_box.is_open_ymin()
            || a_box.is_open_ymax()
            || a_box.is_open_zmin()
            || a_box.is_open_zmax()
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
    ///
    /// Planes use the analytic projector (OCCT `GeomAPI_ProjectPointOnSurf` on
    /// `Geom_Plane`). A grid over `sample_bounds` clamps unbounded plane UV to
    /// `[-1, 1]` and can report a ~1 unit miss for a point that lies on the
    /// plane.
    pub fn project_point_on_face(&self, face: &Face, p: &GpPnt) -> Result<(f64, f64), String> {
        let Some(surf) = BRepTool::face_surface(face) else {
            return Err("IntToolsContext::project_point_on_face: face has no surface".into());
        };
        if is_planar(surf.as_ref(), 6, 6, 1e-6) {
            let (u, v, _) = plane_projection(surf.as_ref(), p);
            return Ok((u, v));
        }
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

    /// Point-on-edge classification (`IntTools_Context::ComputePE(gp_Pnt, ...)`).
    ///
    /// Returns `(status, parameter, distance)`. Status `0` is a hit; negative
    /// codes match the OCCT overload: `-2` no curve, `-3` projection missed
    /// the vertices, `-4` farther than `tol_p + edge tol + Confusion`.
    pub fn compute_pe_pnt(&self, p: &GpPnt, tol_p: f64, edge: &Edge) -> (i32, f64, f64) {
        if BRepTool::edge_curve(edge).is_none() {
            return (-2, 0.0, 0.0);
        }
        if let Some(t) = self.project_point_on_edge(edge, p) {
            let Some(curve) = BRepTool::edge_curve(edge) else {
                return (-2, 0.0, 0.0);
            };
            let dist = p.distance(&curve.d0(t));
            let tol_sum = tol_p.max(0.0) + BRepTool::edge_tolerance(edge) + CONFUSION;
            if dist > tol_sum {
                return (-4, t, dist);
            }
            return (0, t, dist);
        }
        let mut dist = f64::MAX;
        let mut t = 0.0;
        for child in ShapeIterator::of_shape(&edge.0) {
            if child.shape_type() != ShapeType::Vertex {
                continue;
            }
            let ori = child.orientation();
            if ori != Orientation::Forward && ori != Orientation::Reversed {
                continue;
            }
            let v = Vertex(child);
            let pv = BRepTool::vertex_point(&v);
            let tol_sum = tol_p.max(0.0) + BRepTool::vertex_tolerance(&v) + CONFUSION;
            let d = p.distance(&pv);
            if d < dist && d < tol_sum {
                dist = d;
                t = if ori == Orientation::Forward {
                    BRepTool::parameter_on_edge(edge, 0)
                } else {
                    BRepTool::parameter_on_edge(edge, 1)
                };
            }
        }
        if occt_core::precision::Precision::is_infinite(dist) {
            return (-3, t, dist);
        }
        (0, t, dist)
    }

    /// Classify the two edges `edge1`/`edge2` against each other
    /// (simplified `IntTools_EdgeEdge`).
    ///
    /// No live PaveFiller caller uses this method: VE projection is
    /// `ComputeVE` (vertex/edge) in [`crate::pave_ve`], and EE common parts
    /// go through [`crate::edge_edge::EdgeEdge`]. Sampled closest-approach
    /// is kept for the `IntTools_Context` surface; do not wire `EdgeEdge` here
    /// until a live caller needs exact common parts.
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
        let Ok((u, v)) = self.project_point_on_face(face, &p) else {
            return -1;
        };
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
        crate::int_tools_vertex_line::is_vertex_on_line_bool(vertex, curve, tol)
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
    pub(super) fn point_face_state(
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
            // `IsValidPointForFace` uses `ProjPS` (`GeomAPI_ProjectPointOnSurf`).
            // The 32x32 grid in `surface_closest_params` clamps an unbounded
            // plane to `[-1, 1]` and can report a ~1 unit miss for a point
            // that lies on the plane (e.g. overlapping-box faces at x=1.5).
            None => self.project_point_on_face(face, p)?,
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
pub(super) fn face_block_bounds(cl: &FClass2d, face: &Face) -> (f64, f64, f64, f64) {
    let (u0, u1, v0, v1) = face_uv_bounds(face);
    if u0.is_finite() && u1.is_finite() && v0.is_finite() && v1.is_finite() {
        return (u0, u1, v0, v1);
    }
    ring_bounds(cl)
}

/// Bounding box of the classifier's sampled UV boundary rings.
pub(super) fn ring_bounds(cl: &FClass2d) -> (f64, f64, f64, f64) {
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
pub(super) fn project_point_on_curve(p: &GpPnt, curve: &dyn Curve) -> Option<(f64, f64)> {
    let (a, b) = curve_param_window(curve, p);
    let t = inttools_roots::parameter(&|u| curve.d0(u), p, a, b);
    Some((t, p.distance(&curve.d0(t))))
}

/// A finite parameter window over which to project `p` onto `curve`.
///
/// Bounded curves use their natural range; unbounded curves (infinite lines)
/// project `p` onto the tangent at a reference point and return a window around
/// the estimated closest parameter (exact for straight lines).
pub(super) fn curve_param_window(curve: &dyn Curve, near: &GpPnt) -> (f64, f64) {
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
pub(super) fn min_curve_distance(
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
