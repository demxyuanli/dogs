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

use occt_core::gp::{GpPnt, GpPnt2d, GpVec};
use occt_geom::Curve;

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
    /// ponytail: OCCT `BRepMesh_EdgeParameterProvider::Parameter` additionally
    /// projects the corresponding 3D point onto the pcurve (`Extrema_LocateExtPC`)
    /// to refine the non-SameParameter result and guard against parameter
    /// regressions on periodic surfaces. That projection needs a pcurve adaptor
    /// + a 3D point, neither carried by this stand-in provider (edges produced
    /// by `make_pcurve_full` are same-parameter, so the path is unreachable
    /// today); the linear rescale below is the faithful degenerate case.
    pub fn remap(&self, stored: f64) -> f64 {
        if self.is_same_param {
            return stored;
        }
        self.first + self.scale * (stored - self.old_first)
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

/// Max recursion depth while splitting a curve segment by deflection.
const MAX_SPLIT_DEPTH: usize = 24;

/// Flattens a parametric curve into a deflection-bounded polyline.
///
/// Port of `BRepMesh_CurveTessellator` / `GCPnts_TangentialDeflection`: the curve
/// is seeded with at least `min_points` uniform samples and every segment is
/// then refined while its chord (evaluated at the parameter midpoint) deviates
/// from the true curve by more than `deflection`, or while the angle between the
/// segment end tangents exceeds `angular_deflection` (disabled when infinite).
#[derive(Clone)]
pub struct CurveTessellator {
    curve: Arc<dyn Curve>,
    deflection: f64,
    angular_deflection: f64,
    min_points: usize,
    points: Vec<GpPnt>,
    params: Vec<f64>,
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
        let mut t = Self {
            curve,
            deflection: deflection.max(1e-12),
            angular_deflection,
            min_points: min_points.max(2),
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

        let n_seed = self.min_points.max(2);
        let span = b - a;
        let mut stack: Vec<(f64, f64, usize)> = Vec::with_capacity(n_seed);
        // Push seed segments last-to-first so they pop in ascending order.
        for i in (0..n_seed - 1).rev() {
            let u0 = a + span * i as f64 / (n_seed - 1) as f64;
            let u1 = a + span * (i + 1) as f64 / (n_seed - 1) as f64;
            stack.push((u0, u1, 0));
        }

        let def = self.deflection;
        let ang = self.angular_deflection;
        let mut final_params: Vec<f64> = Vec::new();
        let mut last: Option<f64> = None;
        while let Some((lo, hi, depth)) = stack.pop() {
            let mid = 0.5 * (lo + hi);
            let (pa, pm, pb) = (self.curve.d0(lo), self.curve.d0(mid), self.curve.d0(hi));
            let linear = point_segment_dist(&pm, &pa, &pb) > def;
            let angular = ang.is_finite() && {
                let (_, ta) = self.curve.d1(lo);
                let (_, tb) = self.curve.d1(hi);
                angle_between(&ta, &tb) > ang
            };
            if (linear || angular) && depth < MAX_SPLIT_DEPTH {
                stack.push((mid, hi, depth + 1));
                stack.push((lo, mid, depth + 1));
            } else {
                if last != Some(lo) {
                    final_params.push(lo);
                    last = Some(lo);
                }
                if last != Some(hi) {
                    final_params.push(hi);
                    last = Some(hi);
                }
            }
        }

        if final_params.first() != Some(&a) {
            final_params.insert(0, a);
        }
        if final_params.last() != Some(&b) {
            final_params.push(b);
        }

        self.params = final_params;
        self.points = self.params.iter().map(|&u| self.curve.d0(u)).collect();
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

/// Angle (radians) between two vectors.
fn angle_between(a: &GpVec, b: &GpVec) -> f64 {
    let denom = a.magnitude() * b.magnitude();
    if denom < 1e-30 {
        return 0.0;
    }
    (a.dot(b) / denom).clamp(-1.0, 1.0).acos()
}

/// Distance from `p` to the segment `a..b`.
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
