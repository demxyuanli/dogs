//! Port of OCCT edge_discret — Wave 1 BRepMesh.
//!
//! `EdgeDiscret` turns a 3D edge into a deflection-bounded polyline
//! (`BRepMesh_EdgeDiscret` + `BRepMesh_CurveTessellator` + `GCPnts_*Deflection`),
//! `EdgeParameterProvider` maps stored edge parameters to the actual pcurve
//! parameters (`BRepMesh_EdgeParameterProvider`), and `CurveTessellator`
//! flattens an `Arc<dyn Curve>` with adaptive chord-deviation control
//! (`BRepMesh_CurveTessellator`).
//!
//! `MeshParameters` is imported from the sibling `parameters.rs`; `MeshEdge`,
//! `MeshFace` and `MeshModel` are local UV-polygon stand-ins (see below) until
//! the `data_model.rs` topology-indexed types are stable.

use std::sync::Arc;

use occt_core::elib::{clib, slib};
use occt_core::gcpnts::{perform_linear, perform_tangential_curve, CurveSecondDeriv};
use occt_core::gp::{GpAx1, GpAx2, GpCirc, GpDir, GpDir2d, GpLin, GpPnt, GpPnt2d, GpVec};
use occt_core::precision::{ANGULAR, CONFUSION, PCONFUSION};
use occt_geom::{geom_api, Curve, Surface};
use occt_geom2d::curve::Curve2d;

use crate::abs::{Orientation, ShapeType};
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Vertex};

// Meshing parameters — owned by the sibling `parameters.rs` stub.
pub use super::parameters::MeshParameters;

// ---------------------------------------------------------------------------
// Sibling-contract stand-ins
// ---------------------------------------------------------------------------
// ponytail: temporary UV-polygon stand-ins mirroring the sibling-stub contract.
// The real `data_model::{MeshEdge,MeshFace}` are topology-index based and not
// stable yet; swap to `use super::data_model::*;` once the discretizers are
// rewritten against the model's wire/edge/pcurve collections.

/// Discrete edge carrying a 3D curve + parameter range + discretization data.
#[derive(Clone)]
pub struct MeshEdge {
    pub curve: Arc<dyn Curve>,
    pub first: f64,
    pub last: f64,
    pub deflection: f64,
    pub angular_deflection: f64,
    pub same_param: bool,
    pub same_range: bool,
    pub degenerated: bool,
    /// Discretized 2D pcurves (one entry per face the edge bounds).
    pub pcurves: Vec<Vec<GpPnt2d>>,
}

impl MeshEdge {
    pub fn new(curve: Arc<dyn Curve>, first: f64, last: f64) -> Self {
        Self {
            curve,
            first,
            last,
            deflection: 0.1,
            angular_deflection: 0.5,
            same_param: true,
            same_range: true,
            degenerated: false,
            pcurves: Vec::new(),
        }
    }
}

/// Discrete face — boundary wires as UV polygons + mesh control values.
#[derive(Debug, Clone)]
pub struct MeshFace {
    /// Outer wire as a closed UV polygon (first == last point).
    pub outer_wire: Vec<GpPnt2d>,
    /// Inner (hole) wires as closed UV polygons.
    pub inner_wires: Vec<Vec<GpPnt2d>>,
    pub deflection: f64,
    pub id: usize,
}

impl MeshFace {
    pub fn new(outer_wire: Vec<GpPnt2d>) -> Self {
        Self { outer_wire, inner_wires: Vec::new(), deflection: 0.1, id: 0 }
    }
}

/// Discrete model — the mesh data shared between the two discretizers.
#[derive(Default, Clone)]
pub struct MeshModel {
    pub edges: Vec<MeshEdge>,
    pub faces: Vec<MeshFace>,
    pub max_size: f64,
}

// ---------------------------------------------------------------------------
// EdgeParameterProvider
// ---------------------------------------------------------------------------

/// Maps stored polygon parameters to the actual curve parameters.
///
/// Port of `BRepMesh_EdgeParameterProvider`. When an edge is *not* SameParameter
/// the stored parameters are rescaled linearly onto `[first, last]`; for
/// SameParameter edges stored parameters are used verbatim. Also provides the
/// uniform parameter grids used to seed edge discretization.
#[derive(Debug, Clone)]
pub struct EdgeParameterProvider {
    first: f64,
    last: f64,
    is_same_param: bool,
    old_first: f64,
    old_last: f64,
    scale: f64,
    /// `myCurParam` / `myFoundParam` (`hxx:86, 119-136`).
    cur_param: f64,
    found_param: f64,
}

impl EdgeParameterProvider {
    /// Provider over `[first, last]` assuming the edge is SameParameter.
    pub fn new(first: f64, last: f64) -> Self {
        Self {
            first,
            last,
            is_same_param: true,
            old_first: first,
            old_last: last,
            scale: 1.0,
            cur_param: first,
            found_param: first,
        }
    }

    /// Provider aware of the stored parameter range `(old_first, old_last)`.
    /// The linear scale factor `(last - first) / (old_last - old_first)` maps a
    /// stored parameter onto the actual range, mirroring the OCCT constructor.
    pub fn with_stored(first: f64, last: f64, old_first: f64, old_last: f64) -> Self {
        let denom = old_last - old_first;
        let scale = if denom.abs() > f64::EPSILON && first.is_finite() && last.is_finite() {
            (last - first) / denom
        } else {
            1.0
        };
        Self {
            first,
            last,
            is_same_param: false,
            old_first,
            old_last,
            scale,
            cur_param: first,
            found_param: first,
        }
    }

    /// First parameter of the actual range.
    pub fn first(&self) -> f64 {
        self.first
    }

    /// Last parameter of the actual range.
    pub fn last(&self) -> f64 {
        self.last
    }

    /// `(first, last)` of the actual range.
    pub fn range(&self) -> (f64, f64) {
        (self.first, self.last)
    }

    /// Whether the edge is SameParameter.
    pub fn is_same_param(&self) -> bool {
        self.is_same_param
    }

    /// Linear scale factor between the stored and the actual parameter ranges.
    pub fn scale(&self) -> f64 {
        self.scale
    }

    /// Length of the actual parameter range.
    pub fn length(&self) -> f64 {
        self.last - self.first
    }

    /// Clamp `u` into `[first, last]`.
    pub fn clamp(&self, u: f64) -> f64 {
        if self.last >= self.first {
            u.clamp(self.first, self.last)
        } else {
            u
        }
    }

    /// Map a stored parameter onto the actual range. SameParameter edges pass
    /// the value through unchanged.
    ///
    /// Projection-free linear rescale: the `Parameter` result before the
    /// `Extrema_LocateExtPC` refinement. Kept for callers that carry only 2D
    /// data (`tessellate_2d`, `parameter`). Callers that do have the 3D point
    /// and the curve-on-surface must use `parameter_of` below, which implements
    /// the full `BRepMesh_EdgeParameterProvider::Parameter` refinement and IS
    /// the path taken by the live pipeline (`incremental_mesh/discret_root.rs:921,1229`).
    pub fn remap(&self, stored: f64) -> f64 {
        if self.is_same_param {
            return stored;
        }
        self.first + self.scale * (stored - self.old_first)
    }

    /// `BRepMesh_EdgeParameterProvider::Parameter` (`hxx:109-139`).
    /// SameParameter returns `stored`. Otherwise scale, then
    /// `Extrema_LocateExtPC` on the face `BRepAdaptor_Curve` (curve-on-surface)
    /// with the period/regression guard.
    pub fn parameter_of(&mut self, stored: f64, point3d: &GpPnt, cos: &dyn Curve) -> f64 {
        if self.is_same_param {
            return stored;
        }
        let prev_param = self.cur_param;
        self.cur_param = self.first + self.scale * (stored - self.old_first);
        let prev_found = self.found_param;
        self.found_param += self.cur_param - prev_param;
        if let Some((found, _)) = crate::int_tools_vertex_line::extrema_locate_ext_pc(
            cos,
            point3d,
            self.found_param,
            self.first,
            self.last,
        ) {
            if (prev_found < self.found_param && prev_found < found)
                || (prev_found > self.found_param && prev_found > found)
            {
                self.found_param = found;
            }
        }
        self.found_param
    }

    /// `n` parameters uniformly spaced across the actual range, inclusive of
    /// both endpoints. `n` is clamped to at least 2.
    pub fn uniform_parameters(&self, n: usize) -> Vec<f64> {
        let n = n.max(2);
        let span = self.length();
        (0..n)
            .map(|i| self.first + span * i as f64 / (n - 1) as f64)
            .collect()
    }

    /// `n - 1` parameter segments spanning `[first, last]`.
    pub fn segments(&self, n: usize) -> Vec<(f64, f64)> {
        self.uniform_parameters(n).windows(2).map(|w| (w[0], w[1])).collect()
    }

    /// Parameter of the `index`-th point of an `nb`-point grid. When the edge is
    /// not SameParameter the stored grid is remapped onto the actual range.
    pub fn parameter(&self, index: usize, nb: usize) -> f64 {
        let n = nb.max(2);
        let i = index.min(n - 1);
        let stored = self.old_first + (self.old_last - self.old_first) * i as f64 / (n - 1) as f64;
        self.remap(stored)
    }
}

// ---------------------------------------------------------------------------
// CurveTessellator
// ---------------------------------------------------------------------------

/// Flattens a parametric curve into a deflection-bounded polyline.
///
/// Port of `BRepMesh_CurveTessellator` / `GCPnts_TangentialDeflection`: the curve
/// is seeded with at least `min_points` uniform samples and every segment is
/// then refined while its chord (evaluated at the parameter midpoint) deviates
/// from the true curve by more than `deflection`, or while the angle between the
/// segment end tangents exceeds `angular_deflection` (disabled when infinite).
/// Circles use `PerformCircular`. Other types follow
/// `GCPnts_TangentialDeflection::initialize` (`cxx:415-453`).
#[derive(Clone)]
pub struct CurveTessellator {
    curve: Arc<dyn Curve>,
    deflection: f64,
    angular_deflection: f64,
    min_points: usize,
    min_size: f64,
    points: Vec<GpPnt>,
    params: Vec<f64>,
}

struct TessCurve<'a>(&'a dyn Curve);

impl CurveSecondDeriv for TessCurve<'_> {
    fn point(&self, u: f64) -> GpPnt {
        self.0.d0(u)
    }
    fn d2(&self, u: f64) -> (GpPnt, GpVec, GpVec) {
        self.0.d2(u)
    }
}

impl CurveTessellator {
    /// Tessellate `curve` over its natural parameter range.
    pub fn new(curve: Arc<dyn Curve>, deflection: f64, min_points: usize) -> Self {
        let (a, b) = (curve.first_parameter(), curve.last_parameter());
        Self::from_range(curve, a, b, deflection, min_points)
    }

    /// Tessellate `curve` over `[first, last]` so every chord deviates from the
    /// curve by at most `deflection`. At least `min_points` points are emitted.
    pub fn from_range(
        curve: Arc<dyn Curve>,
        first: f64,
        last: f64,
        deflection: f64,
        min_points: usize,
    ) -> Self {
        Self::from_range_angular(curve, first, last, deflection, f64::INFINITY, min_points)
    }

    /// Tessellate `curve` over `[first, last]` bounding both the chord deviation
    /// (`deflection`) and the angular deviation between consecutive end tangents
    /// (`angular_deflection`, `GCPnts_TangentialDeflection`'s angular term).
    pub fn from_range_angular(
        curve: Arc<dyn Curve>,
        first: f64,
        last: f64,
        deflection: f64,
        angular_deflection: f64,
        min_points: usize,
    ) -> Self {
        Self::from_range_angular_min(
            curve,
            first,
            last,
            deflection,
            angular_deflection,
            min_points,
            CONFUSION,
        )
    }

    /// Same as `from_range_angular` with OCCT `Initialize(..., theMinLen)`.
    pub fn from_range_angular_min(
        curve: Arc<dyn Curve>,
        first: f64,
        last: f64,
        deflection: f64,
        angular_deflection: f64,
        min_points: usize,
        min_size: f64,
    ) -> Self {
        let mut t = Self {
            curve,
            deflection: deflection.max(1e-12),
            angular_deflection,
            min_points: min_points.max(2),
            min_size: min_size.max(CONFUSION),
            points: Vec::new(),
            params: Vec::new(),
        };
        t.build(first, last);
        t
    }

    /// Number of tessellation points.
    pub fn points_nb(&self) -> usize {
        self.points.len()
    }

    /// The `index`-th point (0-based). `None` out of range.
    pub fn point(&self, index: usize) -> Option<GpPnt> {
        self.points.get(index).copied()
    }

    /// The parameter of the `index`-th point (0-based). `None` out of range.
    pub fn parameter(&self, index: usize) -> Option<f64> {
        self.params.get(index).copied()
    }

    /// `(point, parameter)` of the `index`-th solution.
    pub fn value(&self, index: usize) -> Option<(GpPnt, f64)> {
        match (self.points.get(index), self.params.get(index)) {
            (Some(&p), Some(&u)) => Some((p, u)),
            _ => None,
        }
    }

    /// All points.
    pub fn points(&self) -> &[GpPnt] {
        &self.points
    }

    /// All parameters.
    pub fn params(&self) -> &[f64] {
        &self.params
    }

    /// `GCPnts_TangentialDeflection::AddPoint` (`cxx:458-491`).
    pub fn add_point(&mut self, pnt: GpPnt, param: f64, is_replace: bool) -> usize {
        let tol = PCONFUSION;
        let nb = self.params.len();
        for i in 0..nb {
            let dist = self.params[i] - param;
            if dist.abs() <= tol {
                if is_replace {
                    if i < self.points.len() {
                        self.points[i] = pnt;
                    }
                    self.params[i] = param;
                }
                return i;
            } else if dist > tol {
                let pi = i.min(self.points.len());
                self.points.insert(pi, pnt);
                self.params.insert(i, param);
                return i;
            }
        }
        self.points.push(pnt);
        self.params.push(param);
        self.params.len() - 1
    }

    /// `BRepMesh_CurveTessellator::addInternalVertices` (`cxx:193-208`).
    ///
    /// `BRep_Tool::Parameter(V,E)` walks `BRep_TVertex` PointOnCurve, which this
    /// port does not store. Fallback: project `BRep_Tool::Pnt` onto the 3D curve.
    pub fn add_internal_vertices(&mut self, edge: &Edge) {
        let stored = edge
            .0
            .tshape
            .read()
            .expect("poisoned TShape lock")
            .children
            .clone();
        for child in stored {
            if child.shape_type() != ShapeType::Vertex {
                continue;
            }
            if child.orientation() != Orientation::Internal {
                continue;
            }
            let vertex = Vertex(child);
            let pnt = BRepTool::vertex_point(&vertex);
            let Some(param) = vertex_parameter_on_edge(&vertex, edge) else {
                continue;
            };
            self.add_point(pnt, param, true);
        }
    }

    /// Measured worst chord deviation. Each segment's chord is sampled at
    /// `samples` interior curve parameters; the furthest distance to the chord
    /// is returned. `samples == 0` samples only the segment midpoint.
    pub fn max_chord_error(&self, samples: usize) -> f64 {
        let mut max_dev = 0.0f64;
        let k = samples.max(1);
        for i in 0..self.points.len().saturating_sub(1) {
            let pa = &self.points[i];
            let pb = &self.points[i + 1];
            let (ua, ub) = if self.params.len() == self.points.len() {
                (self.params[i], self.params[i + 1])
            } else {
                (0.0, 1.0)
            };
            for j in 1..=k {
                let u = ua + (ub - ua) * j as f64 / k as f64;
                let d = point_segment_dist(&self.curve.d0(u), pa, pb);
                if d > max_dev {
                    max_dev = d;
                }
            }
        }
        max_dev
    }

    fn build(&mut self, first: f64, last: f64) {
        let (a, b) = if first.is_finite() && last.is_finite() && last >= first {
            (first, last)
        } else {
            // Unbounded range (e.g. a bare line): fall back to the natural [0,1].
            (0.0, 1.0)
        };
        if b - a < 1e-12 {
            self.points = vec![self.curve.d0(a)];
            self.params = vec![a];
            return;
        }

        // `GCPnts_TangentialDeflection::initialize` dispatches `GeomAbs_Circle`
        // to `PerformCircular` (uniform `ArcAngularStep`). Binary-splitting a
        // 4-point seed on a full circle yields 48 samples at the default 0.5 rad
        // angle; OCCT walks `ceil(span / Du)` instead.
        if let Some(radius) = self.curve.circle_radius() {
            self.build_circular(a, b, radius);
            return;
        }

        // `GCPnts_TangentialDeflection::initialize` (`cxx:415-453`).
        // BSpline intervals are the adaptor-trimmed `NbIntervals(CN)` set
        // (`GeomAdaptor_Curve.cxx:371-413`), not the untrimmed unique knots.
        let adaptor = TessCurve(self.curve.as_ref());
        let two_poles = self
            .curve
            .bspline_poles()
            .or_else(|| self.curve.bezier_poles())
            .is_some_and(|p| p.len() == 2);
        if self.curve.is_line() || two_poles {
            let (params, points) = perform_linear(&adaptor, a, b, self.min_points);
            self.params = params;
            self.points = points;
            return;
        }

        let ang = if self.angular_deflection.is_finite() && self.angular_deflection > 0.0 {
            self.angular_deflection.max(ANGULAR)
        } else {
            std::f64::consts::PI
        };
        let intervals = if let (Some(knots), Some(deg)) =
            (self.curve.bspline_knots(), self.curve.nurbs_degree())
        {
            occt_core::bspl::adaptor_intervals(
                knots,
                deg,
                self.curve.is_periodic(),
                6,
                a,
                b,
                self.curve.resolution(CONFUSION).min(PCONFUSION),
            )
        } else {
            self.curve.parameter_intervals(6)
        };
        let degree_min_nb = self
            .curve
            .nurbs_degree()
            .map(|d| (d + 1).max(self.min_points))
            .unwrap_or(self.min_points);
        let (params, points) = perform_tangential_curve(
            &adaptor,
            a,
            b,
            ang,
            self.deflection,
            self.min_points,
            PCONFUSION,
            self.min_size,
            &intervals,
            degree_min_nb,
        );
        self.params = params;
        self.points = points;
    }

    /// `GCPnts_TangentialDeflection::PerformCircular` + `ArcAngularStep`.
    /// `Initialize` passes `myMinLen` (`cxx:357-360`, `cxx:413`).
    ///
    /// `GeomAdaptor_Curve::load` (`cxx:252-254`) unwraps `Geom_TrimmedCurve` to
    /// the basis, so Circle `U` is radians. Our STEP trim remaps that interval
    /// onto `[0, 1]`; `aDiff` for `ceil(aDiff / Du)` is the basis angle span,
    /// while stored parameters stay in `[first, last]` for SameParameter
    /// Tessellate2d / KeepParam.
    fn build_circular(&mut self, first: f64, last: f64, radius: f64) {
        let ang = if self.angular_deflection.is_finite() && self.angular_deflection > 0.0 {
            self.angular_deflection
        } else {
            std::f64::consts::PI
        };
        let mut du =
            super::range_splitter::arc_angular_step(radius, self.deflection, ang, self.min_size);
        let param_span = last - first;
        let angle_span = if let Some((bf, bl)) = self.curve.trimmed_basis_range() {
            let full = self.curve.last_parameter() - self.curve.first_parameter();
            if full.abs() > 1e-16 {
                (bl - bf).abs() * param_span / full
            } else {
                param_span
            }
        } else {
            param_span
        };
        if du <= 1e-12 {
            du = angle_span.abs();
        }
        let diff = angle_span.abs();
        let mut nb = (diff / du).ceil().min(1.0e6) as i32;
        nb = nb.max(self.min_points as i32 - 1).max(1);
        let du_param = param_span / nb as f64;
        self.params.clear();
        self.points.clear();
        let mut u = first;
        for _ in 0..nb {
            self.params.push(u);
            self.points.push(self.curve.d0(u));
            u += du_param;
        }
        self.params.push(last);
        self.points.push(self.curve.d0(last));
    }
}

// ---------------------------------------------------------------------------
// Deflection
// ---------------------------------------------------------------------------

/// Deflection computation helpers (`BRepMesh_Deflection`).
pub struct Deflection;

impl Deflection {
    /// Absolute deflection for a shape of `max_shape_size` given a relative
    /// deflection. Mirrors `ComputeAbsoluteDeflection` (clamped 0.5×..2×).
    pub fn compute_absolute_deflection(relative: f64, max_shape_size: f64) -> f64 {
        if relative <= 0.0 {
            return 0.0;
        }
        if max_shape_size <= 0.0 {
            return relative;
        }
        let shape_size = max_shape_size.max(relative);
        let coeff = (max_shape_size / (2.0 * shape_size)).clamp(0.5, 2.0);
        coeff * shape_size * relative
    }

    /// Whether an existing polygon deflection fits the required one.
    /// Mirrors `IsConsistent`: current must not exceed required by `ratio`, and
    /// (when quality decrease is allowed) must not undercut it either.
    pub fn is_consistent(
        current: f64,
        required: f64,
        allow_quality_decrease: bool,
        ratio: f64,
    ) -> bool {
        current < (1.0 + ratio) * required
            && (!allow_quality_decrease || current > (1.0 - ratio) * required)
    }
}

// ---------------------------------------------------------------------------
// EdgeDiscret
// ---------------------------------------------------------------------------

/// Discretizes 3D edges into deflection-bounded point sequences.
///
/// Port of `BRepMesh_EdgeDiscret`: computes a per-edge deflection, creates a
/// [`CurveTessellator`] and emits the 3D polyline (optionally mirrored into 2D
/// pcurve points by [`EdgeParameterProvider`]).
#[derive(Debug, Clone)]
pub struct EdgeDiscret {
    params: MeshParameters,
}

impl EdgeDiscret {
    pub fn new(params: MeshParameters) -> Self {
        Self { params }
    }

    /// The parameters in effect.
    pub fn parameters(&self) -> &MeshParameters {
        &self.params
    }

    /// Set the parameters.
    pub fn set_parameters(&mut self, params: MeshParameters) {
        self.params = params;
    }

    /// Deflection used for the given edge: the parameters' deflection, adjusted
    /// for a relative deflection against `max_shape_size`.
    pub fn edge_deflection(&self, max_shape_size: f64) -> f64 {
        if self.params.relative {
            Deflection::compute_absolute_deflection(self.params.deflection, max_shape_size)
        } else {
            self.params.deflection
        }
    }

    /// Create a tessellator for `curve` over `[first, last]`, seeding it with
    /// the minimum point count in effect. Port of `CreateEdgeTessellator`.
    pub fn create_edge_tessellator(
        &self,
        curve: &Arc<dyn Curve>,
        first: f64,
        last: f64,
    ) -> CurveTessellator {
        CurveTessellator::from_range(
            curve.clone(),
            first,
            last,
            self.edge_deflection(0.0),
            2,
        )
    }

    /// Discretize a 3D curve into a deflection-bounded polyline.
    ///
    /// Arc-length/parameter adaptive: the curve is seeded uniformly and refined
    /// while any chord deviates from the curve by more than `deflection`.
    /// Endpoints are always present.
    pub fn discretize_edge(
        curve: &Arc<dyn Curve>,
        first: f64,
        last: f64,
        deflection: f64,
    ) -> Vec<GpPnt> {
        CurveTessellator::from_range(curve.clone(), first, last, deflection, 2).points
    }

    /// Discretize an edge of the mesh model into a 3D polyline using the edge's
    /// own deflection. Degenerated edges collapse to their endpoints.
    pub fn discretize_edge_mesh(&self, edge: &MeshEdge) -> Vec<GpPnt> {
        let def = if edge.deflection > 0.0 {
            edge.deflection
        } else {
            self.edge_deflection(0.0)
        };
        let tess = CurveTessellator::from_range(
            edge.curve.clone(),
            edge.first,
            edge.last,
            def,
            2,
        );
        if edge.degenerated {
            // Degenerated edges contribute their start point only.
            return tess.points.first().copied().map(|p| vec![p]).unwrap_or_default();
        }
        tess.points
    }

    /// Port of `Tessellate3d(theUpdateEnds=true)`: the first and last points are
    /// the edge's vertex points (`BRep_Tool::Pnt(firstVertex)` /
    /// `Pnt(lastVertex)`), the interior points come from the tessellator
    /// (OCCT's indices 2..PointsNb−1). Degenerated edges keep only the start
    /// vertex.
    pub fn tessellate_3d(
        edge: &MeshEdge,
        tessellator: &CurveTessellator,
        first_vertex: GpPnt,
        last_vertex: GpPnt,
    ) -> Vec<GpPnt> {
        let n = tessellator.points_nb();
        if n == 0 {
            return Vec::new();
        }
        let mut pts = Vec::with_capacity(n);
        pts.push(first_vertex);
        if !edge.degenerated {
            // OCCT: `for (i = 2; i < PointsNb(); ++i)` — interior points only.
            for i in 1..n.saturating_sub(1) {
                if let Some(p) = tessellator.point(i) {
                    if pts.last().is_none_or(|q: &GpPnt| q.distance(&p) > 1e-12) {
                        pts.push(p);
                    }
                }
            }
        }
        pts.push(last_vertex);
        pts
    }

    /// Port of `Tessellate2d(theUpdateEnds=true)`: evaluate the 2D pcurve at the
    /// 3D tessellation parameters (mapped through the provider). The provider
    /// passes a SameParameter edge's parameters through verbatim and linearly
    /// rescales a non-SameParameter edge's stored parameters onto the pcurve
    /// range (`BRepMesh_EdgeParameterProvider::Parameter` without the 3D-point
    /// projection — see `EdgeParameterProvider`).
    pub fn tessellate_2d(
        provider: &EdgeParameterProvider,
        tessellator: &CurveTessellator,
        pcurve: &dyn Fn(f64) -> GpPnt2d,
    ) -> Vec<GpPnt2d> {
        let mut out = Vec::with_capacity(tessellator.points_nb());
        for &u in tessellator.params() {
            let actual = provider.remap(u);
            let p2 = pcurve(actual);
            if out.last().is_none_or(|q: &GpPnt2d| q.distance(&p2) > 1e-12) {
                out.push(p2);
            }
        }
        out
    }
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// Distance from `p` to the segment `a..b`.
/// `BRepMesh_CurveTessellator::splitByDeflection2d` (`cxx:157-188`).
///
/// `aNodesNb` is captured before the pcurve loop. Each pcurve rebuilds
/// `aParamArray` from the current discretizer's first `aNodesNb` parameters
/// (`cxx:175-186`), so splits from an earlier pcurve shift which original
/// spans the later pcurves still see.
pub fn split_by_deflection2d(
    curve: &dyn Curve,
    params: &mut Vec<f64>,
    pcurves: &[(Arc<dyn Curve2d>, Arc<dyn Surface>)],
    lin_def: f64,
    min_size: f64,
) {
    if params.len() < 2 {
        return;
    }
    let sq_def = lin_def * lin_def;
    let sq_min = sq_def.max(min_size * min_size);
    let a_nodes_nb = params.len();
    for (c2d, surf) in pcurves {
        // `BRepMesh_CurveTessellator.cxx:168-171`: skip only
        // `STANDARD_TYPE(Geom_Plane)`, not a sampled is_planar classify.
        // Offset / BSpline / Rev stay in the loop (unlike Value cxx:241-246).
        if surf.gp_pln().is_some() {
            continue;
        }
        params.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        params.dedup_by(|a, b| (*a - *b).abs() < PCONFUSION);
        let n = a_nodes_nb.min(params.len());
        let arr = params[..n].to_vec();
        for w in arr.windows(2) {
            split_segment_2d(
                curve,
                c2d.as_ref(),
                surf.as_ref(),
                w[0],
                w[1],
                1,
                sq_def,
                sq_min,
                params,
            );
        }
    }
    params.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    params.dedup_by(|a, b| (*a - *b).abs() < PCONFUSION);
}

/// `BRepMesh_CurveTessellator::splitSegment` (`cxx:274-340`).
fn split_segment_2d(
    c3d: &dyn Curve,
    c2d: &dyn Curve2d,
    surf: &dyn Surface,
    first: f64,
    last: f64,
    iter: i32,
    sq_def: f64,
    sq_min: f64,
    params: &mut Vec<f64>,
) {
    if iter > 10 {
        return;
    }
    if (last - first).abs() < 2.0 * PCONFUSION {
        return;
    }
    let (cf, cl) = (c2d.first_parameter(), c2d.last_parameter());
    if cf.is_finite() && first - cf < -PCONFUSION {
        return;
    }
    if cl.is_finite() && last - cl > PCONFUSION {
        return;
    }
    let uvf = c2d.d0(first);
    let uvl = c2d.d0(last);
    let p3df = surf.d0(uvf.x(), uvf.y());
    let p3dl = surf.d0(uvl.x(), uvl.y());
    if p3df.square_distance(&p3dl) < sq_min {
        return;
    }
    let uvm = GpPnt2d::new(0.5 * (uvf.x() + uvl.x()), 0.5 * (uvf.y() + uvl.y()));
    let mid_surf = surf.d0(uvm.x(), uvm.y());
    let v1 = GpVec::from_pnts(&p3df, &mid_surf);
    if v1.square_magnitude() < sq_min {
        return;
    }
    let mut a_vec = GpVec::from_pnts(&p3df, &p3dl);
    let mag = a_vec.magnitude();
    if mag < 1e-30 {
        return;
    }
    a_vec = a_vec.multiplied_scalar(1.0 / mag);
    let proj = a_vec.multiplied_scalar(v1.dot(&a_vec));
    let dist = v1.subtracted(&proj);
    if dist.square_magnitude() < sq_def {
        return;
    }
    let midpar = 0.5 * (first + last);
    // `cxx:335-336` AddPoint(myCurve.D0(midpar), midpar, false); params-only
    // path records midpar — Tessellate3d evaluates the SameParam 3D curve later.
    params.push(midpar);
    split_segment_2d(
        c3d, c2d, surf, first, midpar, iter + 1, sq_def, sq_min, params,
    );
    split_segment_2d(
        c3d, c2d, surf, midpar, last, iter + 1, sq_def, sq_min, params,
    );
}

/// `BRepMesh_CurveTessellator::Value` (`cxx:223-270`) for CurveOnSurface:
/// drop an interior sample whose pcurve UV is outside the surface range
/// (padded by `U/VResolution(Confusion)`) when the 3D point is farther than
/// the edge tolerance from `surface(UV)`. Analytic and periodic surfaces keep
/// every sample.
///
/// Type gate matches `cxx:241-246`: only `BSpline` / `Bezier` / `OtherSurface`
/// run the UV out-of-range check. `OffsetSurface` is excluded even when the
/// basis exposes poles (`GetType() == GeomAbs_OffsetSurface`).
pub fn curve_tessellator_value_ok(
    pcurve: &dyn Curve2d,
    surface: &dyn Surface,
    parameter: f64,
    point: &GpPnt,
    edge_tol: f64,
) -> bool {
    use super::range_splitter::SurfaceType;
    let ty = super::range_splitter::classify_surface(surface);
    match ty {
        SurfaceType::BSplineSurface
        | SurfaceType::BezierSurface
        | SurfaceType::OtherSurface => {}
        _ => return true,
    }
    if surface.is_u_periodic() || surface.is_v_periodic() {
        return true;
    }
    let uv = pcurve.d0(parameter);
    let (u0, u1) = surface.u_range();
    let (v0, v1) = surface.v_range();
    let (du, dv) = surface
        .uv_resolution(CONFUSION)
        .unwrap_or((CONFUSION, CONFUSION));
    if uv.x() > u0 - du && uv.x() < u1 + du && uv.y() > v0 - dv && uv.y() < v1 + dv {
        return true;
    }
    let on_surf = surface.d0(uv.x(), uv.y());
    point.square_distance(&on_surf) < edge_tol * edge_tol
}

/// Parameter of an INTERNAL vertex on `edge`. PointOnCurve is unported;
/// `GeomAPI_ProjectPointOnCurve` on the 3D curve stands in for
/// `BRep_Tool::Parameter(V,E)` (`BRep_Tool.cxx:1503-1510`).
fn vertex_parameter_on_edge(vertex: &Vertex, edge: &Edge) -> Option<f64> {
    let pnt = BRepTool::vertex_point(vertex);
    let curve = BRepTool::edge_curve(edge)?;
    let tol = BRepTool::vertex_tolerance(vertex);
    geom_api::project_point_on_curve(curve.as_ref(), &pnt, tol).map(|p| p.parameter)
}

/// `BRepMesh_CurveTessellator.cxx:100-114` min-point threshold from `GetType`.
/// Circle uses `circle_radius`. Ellipse is periodic and not a NURBS. Periodic
/// BSpline/Bezier stay at 2 (`GeomAbs_BSplineCurve`). Parabola/hyperbola
/// GetType is unported without a downcast (comment: cxx:103-106).
pub fn tessellator_min_points(curve: &dyn Curve) -> usize {
    if curve.circle_radius().is_some() {
        return 4;
    }
    if curve.is_periodic() && curve.bspline_poles().is_none() && curve.bezier_poles().is_none() {
        return 4;
    }
    2
}

/// `Adaptor3d_CurveOnSurface::EvalKPart` reverse of a classified circle
/// (`cxx:1597-1602` and the same SetDirection pattern on every iso Circle).
fn reverse_kpart_circ(circ: &mut GpCirc) {
    let mut ax = circ.position();
    ax.set_direction(ax.direction().reversed());
    circ.set_position(&ax);
}

/// `to3d(Pl, Circ2d)` (`Adaptor3d_CurveOnSurface.cxx:57-83`).
fn circ2d_to3d_on_plane(pl: &occt_core::gp::GpPln, c: &occt_core::gp::GpCirc2d) -> Option<GpCirc> {
    let loc = c.position().location();
    let p = slib::plane_value(pl, loc.x(), loc.y());
    let vx2 = *c.position().x_direction();
    let vy2 = *c.position().y_direction();
    let vx = GpVec::from_xyz(
        &pl.pos
            .x_direction()
            .xyz()
            .multiplied(vx2.x)
            .added(&pl.pos.y_direction().xyz().multiplied(vx2.y)),
    );
    let vy = GpVec::from_xyz(
        &pl.pos
            .x_direction()
            .xyz()
            .multiplied(vy2.x)
            .added(&pl.pos.y_direction().xyz().multiplied(vy2.y)),
    );
    let n = vx.crossed(&vy);
    let xd = GpDir::from_vec(&vx).ok()?;
    let nd = GpDir::from_vec(&n).ok()?;
    let ax = GpAx2::new(p, nd, xd).ok()?;
    Some(GpCirc::new(ax, c.radius()))
}

/// Rotate a V-iso circle by the 2d line U (`cxx:1594-1596`).
fn rotate_iso_v(circ: &mut GpCirc, axis: &occt_core::gp::GpAx3, u: f64) {
    let drev = axis
        .x_direction()
        .crossed(axis.y_direction())
        .unwrap_or(*axis.x_direction());
    let axe = GpAx1::new(axis.location(), drev);
    circ.rotate(&axe, u);
}

/// `Adaptor3d_CurveOnSurface::EvalKPart` (`cxx:1552-1732`).
/// `GeomAdaptor_Surface` of a rectangular trim reports the basis GetType.
fn eval_k_part(pc: &dyn Curve2d, surf: &dyn Surface) -> (Option<GpCirc>, Option<GpLin>) {
    if let Some(basis) = surf.rectangular_trimmed_basis() {
        return eval_k_part(pc, basis.as_ref());
    }
    let mut circ = None;
    let mut lin = None;
    if let Some(pl) = surf.gp_pln() {
        if let Some(c2) = pc.gp_circ2d() {
            circ = circ2d_to3d_on_plane(&pl, &c2);
        } else if pc.is_line() {
            let (uv, duv) = pc.d1(0.0);
            let (p, d1u, d1v) = surf.d1(uv.x(), uv.y());
            let v = d1u
                .multiplied_scalar(duv.x())
                .added(&d1v.multiplied_scalar(duv.y()));
            if let Ok(dir) = GpDir::from_vec(&v) {
                lin = Some(GpLin::from_pnt_dir(p, dir));
            }
        }
        return (circ, lin);
    }
    let Some(l2) = pc.gp_lin2d() else {
        return (None, None);
    };
    let d = *l2.direction();
    let loc = l2.location();
    const DX2D: GpDir2d = GpDir2d { x: 1.0, y: 0.0 };
    const DY2D: GpDir2d = GpDir2d { x: 0.0, y: 1.0 };
    if d.is_parallel(&DX2D, ANGULAR) {
        if let Some(sph) = surf.gp_sphere() {
            if (loc.y().abs() - std::f64::consts::FRAC_PI_2).abs() >= PCONFUSION {
                let axis = sph.position();
                let mut c = slib::sphere_v_iso(axis, sph.radius(), loc.y());
                rotate_iso_v(&mut c, axis, loc.x());
                if d.is_opposite(&DX2D, ANGULAR) {
                    reverse_kpart_circ(&mut c);
                }
                circ = Some(c);
            }
        } else if let Some(cyl) = surf.gp_cylinder() {
            let axis = cyl.position();
            let mut c = slib::cylinder_v_iso(&axis, cyl.radius(), loc.y());
            rotate_iso_v(&mut c, &axis, loc.x());
            if d.is_opposite(&DX2D, ANGULAR) {
                reverse_kpart_circ(&mut c);
            }
            circ = Some(c);
        } else if let Some(cone) = surf.gp_cone() {
            let axis = cone.position();
            let mut c = slib::cone_v_iso(&axis, cone.radius(), cone.semi_angle(), loc.y());
            rotate_iso_v(&mut c, &axis, loc.x());
            if d.is_opposite(&DX2D, ANGULAR) {
                reverse_kpart_circ(&mut c);
            }
            circ = Some(c);
        } else if let Some(tor) = surf.gp_torus() {
            let axis = tor.position();
            let mut c = slib::torus_v_iso(axis, tor.major_radius(), tor.minor_radius(), loc.y());
            rotate_iso_v(&mut c, axis, loc.x());
            if d.is_opposite(&DX2D, ANGULAR) {
                reverse_kpart_circ(&mut c);
            }
            circ = Some(c);
        }
    } else if d.is_parallel(&DY2D, ANGULAR) {
        if let Some(sph) = surf.gp_sphere() {
            let axis = sph.position();
            let mut c = slib::sphere_u_iso(axis, sph.radius(), 0.0);
            let drev = axis
                .x_direction()
                .crossed(&axis.direction())
                .unwrap_or(*axis.x_direction());
            let axe_y = GpAx1::new(axis.location(), drev);
            c.rotate(&axe_y, loc.y());
            rotate_iso_v(&mut c, axis, loc.x());
            if d.is_opposite(&DY2D, ANGULAR) {
                reverse_kpart_circ(&mut c);
            }
            circ = Some(c);
        } else if let Some(cyl) = surf.gp_cylinder() {
            let mut l = slib::cylinder_u_iso(&cyl.position(), cyl.radius(), loc.x());
            let tr = GpVec::from_xyz(l.direction().xyz()).multiplied_scalar(loc.y());
            l = l.translated_vec(&tr);
            if d.is_opposite(&DY2D, ANGULAR) {
                l.set_direction(l.direction().reversed());
            }
            lin = Some(l);
        } else if let Some(cone) = surf.gp_cone() {
            let mut l = slib::cone_u_iso(
                &cone.position(),
                cone.radius(),
                cone.semi_angle(),
                loc.x(),
            );
            let tr = GpVec::from_xyz(l.direction().xyz()).multiplied_scalar(loc.y());
            l = l.translated_vec(&tr);
            if d.is_opposite(&DY2D, ANGULAR) {
                l.set_direction(l.direction().reversed());
            }
            lin = Some(l);
        } else if let Some(tor) = surf.gp_torus() {
            let mut c = slib::torus_u_iso(
                tor.position(),
                tor.major_radius(),
                tor.minor_radius(),
                loc.x(),
            );
            let axe = *c.position().axis();
            c.rotate(&axe, loc.y());
            if d.is_opposite(&DY2D, ANGULAR) {
                reverse_kpart_circ(&mut c);
            }
            circ = Some(c);
        }
    }
    (circ, lin)
}

/// `Adaptor3d_CurveOnSurface` used by `BRepAdaptor_Curve(edge, face)`.
pub struct CurveOnSurface {
    pcurve: Arc<dyn Curve2d>,
    surface: Arc<dyn Surface>,
    first: f64,
    last: f64,
    kpart_circ: Option<GpCirc>,
    kpart_lin: Option<GpLin>,
    /// `myFirstSurf` from `EvalFirstLastSurf` (`cxx:1778-1783`).
    first_surf: Option<Arc<dyn Surface>>,
    /// `myLastSurf` from `EvalFirstLastSurf` (`cxx:1822-1827`).
    last_surf: Option<Arc<dyn Surface>>,
}

impl CurveOnSurface {
    /// Evaluate `surface(pcurve(t))` on `[first, last]`, with `EvalKPart` and
    /// `EvalFirstLastSurf` (`Adaptor3d_CurveOnSurface::Load` cxx:951-963).
    ///
    /// PORTED: Offset unwrap + `LocatePart_Offset` (BSpline + RevExt) +
    /// `LocatePart_RevExt` + `LocatePart` / `Locate1Coord` (surface+curve) /
    /// `Locate2Coord` (Arr+param) + UTrim/VTrim end patches for D1/D2
    /// (`cxx:1212-1266`, `1833-1866`) + `GeomAdaptor` LocalD1/IfUVBound/Span via
    /// `GeomRectangularTrimmedSurface` (`GeomAdaptor_Surface.cxx:1193-1195`).
    pub fn new(
        pcurve: Arc<dyn Curve2d>,
        surface: Arc<dyn Surface>,
        first: f64,
        last: f64,
    ) -> Self {
        let (kpart_circ, kpart_lin) = eval_k_part(pcurve.as_ref(), surface.as_ref());
        let (first_surf, last_surf) =
            super::cos_locate::eval_first_last_surf(pcurve.as_ref(), &surface, first, last);
        Self {
            pcurve,
            surface,
            first,
            last,
            kpart_circ,
            kpart_lin,
            first_surf,
            last_surf,
        }
    }

    fn end_surf(&self, u: f64) -> Option<&Arc<dyn Surface>> {
        // `Adaptor3d_CurveOnSurface::EvalD1` cxx:1208-1222.
        let tol = PCONFUSION / 10.0;
        if (u - self.first).abs() < tol {
            self.first_surf.as_ref()
        } else if (u - self.last).abs() < tol {
            self.last_surf.as_ref()
        } else {
            None
        }
    }
}

impl Curve for CurveOnSurface {
    fn d0(&self, u: f64) -> GpPnt {
        if let Some(ref c) = self.kpart_circ {
            return clib::circle_value(c, u);
        }
        if let Some(ref l) = self.kpart_lin {
            return clib::line_value(l, u);
        }
        let uv = self.pcurve.d0(u);
        self.surface.d0(uv.x(), uv.y())
    }

    fn d1(&self, u: f64) -> (GpPnt, GpVec) {
        if let Some(ref c) = self.kpart_circ {
            return clib::circle_d1(c, u);
        }
        if let Some(ref l) = self.kpart_lin {
            return clib::line_d1(l, u);
        }
        let (uv, duv) = self.pcurve.d1(u);
        let surf = self.end_surf(u).unwrap_or(&self.surface);
        let (p, su, sv) = surf.d1(uv.x(), uv.y());
        let tan = su
            .multiplied_scalar(duv.x())
            .added(&sv.multiplied_scalar(duv.y()));
        (p, tan)
    }

    fn d2(&self, u: f64) -> (GpPnt, GpVec, GpVec) {
        if let Some(ref c) = self.kpart_circ {
            return clib::circle_d2(c, u);
        }
        if let Some(ref l) = self.kpart_lin {
            return clib::line_d2(l, u);
        }
        let (uv, duv, d2uv) = self.pcurve.d2(u);
        let surf = self.end_surf(u).unwrap_or(&self.surface);
        let (p, su, sv, suu, suv, svv) = surf.d2(uv.x(), uv.y());
        let d1 = su
            .multiplied_scalar(duv.x())
            .added(&sv.multiplied_scalar(duv.y()));
        let d2 = suu
            .multiplied_scalar(duv.x() * duv.x())
            .added(&suv.multiplied_scalar(2.0 * duv.x() * duv.y()))
            .added(&svv.multiplied_scalar(duv.y() * duv.y()))
            .added(&su.multiplied_scalar(d2uv.x()))
            .added(&sv.multiplied_scalar(d2uv.y()));
        (p, d1, d2)
    }

    fn first_parameter(&self) -> f64 {
        self.first
    }
    fn last_parameter(&self) -> f64 {
        self.last
    }
    fn is_periodic(&self) -> bool {
        if self.kpart_circ.is_some() {
            return true;
        }
        self.pcurve.is_periodic()
    }
    fn period(&self) -> f64 {
        if self.kpart_circ.is_some() {
            return 2.0 * std::f64::consts::PI;
        }
        self.pcurve.period()
    }
    fn circle_radius(&self) -> Option<f64> {
        self.kpart_circ.as_ref().map(|c| c.radius())
    }
    fn gp_circ(&self) -> Option<GpCirc> {
        self.kpart_circ.clone()
    }
    fn is_line(&self) -> bool {
        self.kpart_lin.is_some()
    }
    fn continuity(&self) -> u8 {
        self.pcurve.continuity().min(self.surface.continuity())
    }
    fn transform(&mut self, _t: &occt_core::gp::GpTrsf) {}
    fn reverse(&mut self) {
        let tmp = self.first;
        self.first = self.last;
        self.last = tmp;
    }
    fn clone_dyn(&self) -> Box<dyn Curve> {
        Box::new(Self {
            pcurve: Arc::from(self.pcurve.clone_dyn()),
            surface: self.surface.clone(),
            first: self.first,
            last: self.last,
            kpart_circ: self.kpart_circ.clone(),
            kpart_lin: self.kpart_lin,
            first_surf: self.first_surf.clone(),
            last_surf: self.last_surf.clone(),
        })
    }
}

fn point_segment_dist(p: &GpPnt, a: &GpPnt, b: &GpPnt) -> f64 {
    let ab = b.coord.subtracted(&a.coord);
    let len2 = ab.square_modulus();
    if len2 <= f64::EPSILON {
        return p.distance(a);
    }
    let ap = p.coord.subtracted(&a.coord);
    let t = (ap.dot(&ab) / len2).clamp(0.0, 1.0);
    let proj = a.coord.added(&ab.multiplied(t));
    p.coord.subtracted(&proj).modulus()
}

#[cfg(test)]
mod tests {
    use super::*;
    use occt_core::gp::{GpAx2, GpCirc, GpDir, GpLin};
    use occt_geom::{GeomCircle, GeomLine, GeomTrimmedCurve};

    fn segment_curve(len: f64) -> Arc<dyn Curve> {
        let lin = GpLin::from_pnt_dir(GpPnt::zero(), GpDir::new(1.0, 0.0, 0.0).unwrap());
        Arc::new(GeomTrimmedCurve::new(Arc::new(GeomLine::new(lin)), 0.0, len))
    }

    fn arc_curve(radius: f64, from: f64, to: f64) -> Arc<dyn Curve> {
        let circ = GeomCircle::new(GpCirc::new(GpAx2::standard(), radius));
        Arc::new(GeomTrimmedCurve::new(Arc::new(circ), from, to))
    }

    #[test]
    fn straight_edge_endpoints_match() {
        let curve = segment_curve(2.0);
        let pts = EdgeDiscret::discretize_edge(&curve, 0.0, 1.0, 0.05);
        assert!(pts.len() >= 2, "polyline length {}", pts.len());
        assert!(
            pts.first().unwrap().distance(&GpPnt::new(0.0, 0.0, 0.0)) < 1e-9,
            "first endpoint {:?}",
            pts.first()
        );
        assert!(
            pts.last().unwrap().distance(&GpPnt::new(2.0, 0.0, 0.0)) < 1e-9,
            "last endpoint {:?}",
            pts.last()
        );
        // A straight chord deviates nothing: no refinement past the seed.
        assert!(pts.len() <= 64, "unexpected refinement for a line: {}", pts.len());
    }

    #[test]
    fn circle_refines_to_chord_deviation_target() {
        let curve = arc_curve(2.0, 0.0, std::f64::consts::PI);
        let deflection = 0.05;
        let tess = CurveTessellator::from_range(curve, 0.0, 1.0, deflection, 2);
        assert!(tess.points_nb() >= 2);
        // Endpoints: u=0 → (2,0,0), u=1 → (-2,0,0).
        assert!(tess.point(0).unwrap().distance(&GpPnt::new(2.0, 0.0, 0.0)) < 1e-6);
        assert!(tess.point(tess.points_nb() - 1).unwrap().distance(&GpPnt::new(-2.0, 0.0, 0.0)) < 1e-6);
        // The chord deviation of the emitted polyline must meet the target.
        let dev = tess.max_chord_error(16);
        assert!(dev < deflection, "chord deviation {dev} >= deflection {deflection}");
    }

    #[test]
    fn angular_deflection_refines_arc() {
        let curve = arc_curve(1.0, 0.0, std::f64::consts::PI);
        // Loose linear deflection (chord deviation only) vs. tight angular bound:
        // the angular term must add points even when the chord is already short.
        let lin = CurveTessellator::from_range(curve.clone(), 0.0, 1.0, 10.0, 2);
        let ang = CurveTessellator::from_range_angular(curve.clone(), 0.0, 1.0, 10.0, 0.05, 2);
        assert!(
            ang.points_nb() > lin.points_nb(),
            "angular {} should exceed linear {}",
            ang.points_nb(),
            lin.points_nb()
        );
    }

    #[test]
    fn min_points_respected_and_deflection_drives_density() {
        let curve = arc_curve(1.0, 0.0, std::f64::consts::PI);
        let coarse = CurveTessellator::from_range(curve.clone(), 0.0, 1.0, 0.2, 4);
        let fine = CurveTessellator::from_range(curve.clone(), 0.0, 1.0, 0.01, 4);
        assert!(coarse.points_nb() >= 4, "min points violated: {}", coarse.points_nb());
        assert!(fine.points_nb() >= coarse.points_nb(), "finer deflection must refine more");
        assert!(fine.max_chord_error(16) < 0.01);
    }

    #[test]
    fn parameter_provider_remaps_stored_range() {
        let p = EdgeParameterProvider::with_stored(1.0, 3.0, 0.0, 2.0);
        assert!(!p.is_same_param());
        assert!((p.scale() - 1.0).abs() < 1e-12);
        // Stored 1.0 (midpoint of [0,2]) maps to 2.0 (midpoint of [1,3]).
        assert!((p.remap(1.0) - 2.0).abs() < 1e-12);
        // Uniform grid is inclusive of both endpoints.
        let u = p.uniform_parameters(3);
        assert_eq!(u, vec![1.0, 2.0, 3.0]);
        // Same-parameter provider passes values through.
        let same = EdgeParameterProvider::new(0.0, 4.0);
        assert!(same.is_same_param());
        assert_eq!(same.remap(1.5), 1.5);
    }

    #[test]
    fn tessellate_3d_replaces_endpoints_with_vertices() {
        let edge = MeshEdge::new(segment_curve(1.0), 0.0, 1.0);
        let t = CurveTessellator::from_range(edge.curve.clone(), 0.0, 1.0, 0.1, 4);
        let a = GpPnt::new(0.0, 0.0, 0.0);
        let b = GpPnt::new(1.0, 0.0, 0.0);
        let pts = EdgeDiscret::tessellate_3d(&edge, &t, a, b);
        // Endpoints are the vertex points (BRep_Tool::Pnt), not curve values.
        assert_eq!(pts.first(), Some(&a));
        assert_eq!(pts.last(), Some(&b));
        // Interior points come from the tessellator; total ≤ tessellator count.
        assert!(pts.len() <= t.points_nb(), "len {}", pts.len());
        assert!(pts.len() >= 2);
    }

    #[test]
    fn deflection_absolute_and_consistency() {
        // Relative deflection against a 2-unit shape: clamped coefficient ≥ 0.5.
        let abs = Deflection::compute_absolute_deflection(0.01, 2.0);
        assert!(abs > 0.0 && abs <= 0.02);
        assert_eq!(Deflection::compute_absolute_deflection(0.01, 0.0), 0.01);
        assert!(Deflection::is_consistent(0.05, 0.05, false, 0.1));
        assert!(!Deflection::is_consistent(0.5, 0.05, false, 0.1));
        // Quality decrease allowed: too-small current is inconsistent.
        assert!(!Deflection::is_consistent(0.001, 0.05, true, 0.1));
    }
}
