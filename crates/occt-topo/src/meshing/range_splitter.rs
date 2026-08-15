//! Port of OCCT BRepMesh range splitters — Wave 4 BRepMesh.
//!
//! Sources (`src/ModelingAlgorithms/TKMesh/BRepMesh/`):
//! - `BRepMesh_DefaultRangeSplitter.{hxx,cxx}` — base splitter: computes the
//!   discrete UV range of a face from its geometric range and the registered
//!   boundary points, plus a face-basis scale (`delta`).
//! - `BRepMesh_{Cylinder,Cone,Sphere,Torus,NURBS}RangeSplitter.{hxx,cxx}` —
//!   analytical splitters generating interior (surface) nodes.
//! - `BRepMesh_UVParamRangeSplitter.hxx`, `BRepMesh_UndefinedRangeSplitter`,
//!   `BRepMesh_BoundaryParamsRangeSplitter`, `BRepMesh_ExtrusionRangeSplitter`.
//!
//! The OCCT classes form an inheritance tree rooted at `BRepMesh_DefaultRangeSplitter`.
//! Rust models that with a [`RangeSplitter`] trait; every concrete splitter embeds
//! its parent as an `inner` field and exposes it through [`RangeSplitter::base`] /
//! [`RangeSplitter::base_mut`], while the virtual hooks (`compute_delta`,
//! `generate_surface_nodes`, `get_undefined_interval_nb`, …) are trait methods the
//! derived types override.
//!
//! `ponytail:` surfaces arrive as `Arc<dyn Surface>` which cannot be downcast, so
//! the analytic type is classified from the periodic flags + `d0` sampling, and
//! radii (`gp_Cylinder::Radius()`, `gp_Torus::MajorRadius()`, …) are measured by
//! sampling diametrically opposite iso-parameter points instead of reading the
//! concrete `gp_*` handles. The NURBS interval machinery reads knots/poles only
//! through `dyn Surface`, so it is approximated by a continuity-based uniform grid.

use std::collections::BTreeSet;
use std::f64::consts::PI;
use std::sync::Arc;

use occt_core::gp::{GpPnt, GpPnt2d, GpVec};
use occt_core::precision::{CONFUSION, PCONFUSION, RESOLUTION, SQUARE_CONFUSION};
use occt_geom::Surface;

use super::data_model::*;
use super::parameters::MeshParameters;

use crate::brep_tool::BRepTool;

/// Sorted set of real parameters used to track boundary U/V values. Mirrors
/// OCCT's `IMeshData::IMapOfReal`; `f64` does not implement `Ord` (NaN), so the
/// values are kept in a sorted `Vec` ordered by `partial_cmp`.
#[derive(Debug, Clone, Default)]
pub struct ParamSet {
    values: Vec<f64>,
}

impl ParamSet {
    /// Creates an empty parameter set.
    pub fn new() -> Self {
        Self { values: Vec::new() }
    }

    /// Inserts `v` keeping the set sorted; duplicates are ignored.
    pub fn insert(&mut self, v: f64) {
        if let Err(pos) = self
            .values
            .binary_search_by(|x| x.partial_cmp(&v).unwrap_or(std::cmp::Ordering::Equal))
        {
            self.values.insert(pos, v);
        }
    }

    /// True when `v` is present.
    pub fn contains(&self, v: &f64) -> bool {
        self.values
            .binary_search_by(|x| x.partial_cmp(v).unwrap_or(std::cmp::Ordering::Equal))
            .is_ok()
    }

    /// Number of stored values.
    pub fn len(&self) -> usize {
        self.values.len()
    }

    /// True when the set is empty.
    pub fn is_empty(&self) -> bool {
        self.values.is_empty()
    }

    /// Clears the set.
    pub fn clear(&mut self) {
        self.values.clear();
    }

    /// Iterator over the sorted values.
    pub fn iter(&self) -> std::slice::Iter<'_, f64> {
        self.values.iter()
    }

    /// The sorted values as a slice.
    pub fn values(&self) -> &[f64] {
        &self.values
    }
}

/// Analytic surface type used to pick a range splitter. Source:
/// `GeomAbs_SurfaceType`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SurfaceType {
    Plane,
    Cylinder,
    Cone,
    Sphere,
    Torus,
    SurfaceOfRevolution,
    SurfaceOfExtrusion,
    BezierSurface,
    BSplineSurface,
    OffsetSurface,
    OtherSurface,
}

/// Classifies a `dyn Surface` into a [`SurfaceType`].
///
/// `ponytail:` trait objects cannot be downcast, so the classification relies on
/// the periodic flags, the finiteness of the parametric ranges, and a `d0`
/// sampling test that separates a cylinder (constant radius along V) from a cone
/// (linearly varying radius along V). Finite non-periodic surfaces (Bezier,
/// BSpline, Offset, …) are indistinguishable here and are routed to the NURBS
/// splitter, which is the safest general choice.
pub fn classify_surface(s: &dyn Surface) -> SurfaceType {
    let up = s.is_u_periodic();
    let vp = s.is_v_periodic();
    let (u0, u1) = s.u_range();
    let (v0, v1) = s.v_range();
    let u_fin = u0.is_finite() && u1.is_finite();
    let v_fin = v0.is_finite() && v1.is_finite();

    match (up, vp) {
        (true, true) => SurfaceType::Torus,
        (true, false) if v_fin => SurfaceType::Sphere,
        (true, false) => {
            if is_cylinder_like(s) {
                SurfaceType::Cylinder
            } else {
                SurfaceType::Cone
            }
        }
        (false, true) => SurfaceType::SurfaceOfRevolution,
        (false, false) if u_fin && v_fin => SurfaceType::BSplineSurface,
        (false, false) => SurfaceType::Plane,
    }
}

/// True when the (u-periodic, unbounded-V) surface keeps a constant radius along
/// V — the cylinder/cone discriminator.
fn is_cylinder_like(s: &dyn Surface) -> bool {
    let r_at = |v: f64| s.d0(0.0, v).distance(&s.d0(PI, v)) * 0.5;
    let r0 = r_at(0.0);
    let r1 = r_at(1.0);
    (r0 - r1).abs() <= 1e-7 * r0.abs().max(r1.abs()).max(1e-7)
}

/// Radius of a cylinder measured from the surface: half the distance between the
/// diametrically opposite points `(0, v)` and `(π, v)`.
fn cylinder_radius(s: &dyn Surface) -> f64 {
    s.d0(0.0, 0.0).distance(&s.d0(PI, 0.0)) * 0.5
}

/// Radius of a sphere measured from the surface.
fn sphere_radius(s: &dyn Surface) -> f64 {
    s.d0(0.0, 0.0).distance(&s.d0(PI, 0.0)) * 0.5
}

/// `(major, minor)` radii of a torus measured from the surface.
fn torus_radii(s: &dyn Surface) -> (f64, f64) {
    let minor = s.d0(0.0, 0.0).distance(&s.d0(0.0, PI)) * 0.5;
    let major_plus_minor = s.d0(0.0, 0.0).distance(&s.d0(PI, 0.0)) * 0.5;
    (major_plus_minor - minor, minor)
}

/// Radius of the U-circle of a cone at parameter `v` measured from the surface.
fn cone_radius_at(s: &dyn Surface, v: f64) -> f64 {
    s.d0(0.0, v).distance(&s.d0(PI, v)) * 0.5
}

/// Angular step for a circle of `radius` under linear (`lin_deflection`) and
/// angular (`ang_deflection`) deflections, with an optional `min_length` floor.
/// Source: `GCPnts_TangentialDeflection::ArcAngularStep`.
pub fn arc_angular_step(
    radius: f64,
    lin_deflection: f64,
    ang_deflection: f64,
    min_length: f64,
) -> f64 {
    let mut du = 0.0;
    let mut min_size_ang = 0.0;
    if radius > CONFUSION {
        du = (1.0 - lin_deflection / radius).max(0.0);
        if min_length > CONFUSION {
            min_size_ang = (min_length / radius).min(PI * 0.5);
        }
    }
    du = 2.0 * du.acos();
    du.min(ang_deflection).max(min_size_ang)
}

/// Parametric step along U or V that corresponds to a 3D displacement of `tol`
/// (finite-difference estimate). Source: `Adaptor3d_Surface::UResolution` /
/// `VResolution`, approximated.
fn param_resolution(s: &dyn Surface, tol: f64, is_u: bool) -> f64 {
    let (r0, r1) = if is_u { s.u_range() } else { s.v_range() };
    let (o0, o1) = if is_u { s.v_range() } else { s.u_range() };
    let mid = |a: f64, b: f64| if a.is_finite() && b.is_finite() { 0.5 * (a + b) } else { 0.0 };
    let u = mid(r0, r1);
    let v = mid(o0, o1);
    let h = 1e-4;
    let p0 = s.d0(u, v);
    let p1 = if is_u { s.d0(u + h, v) } else { s.d0(u, v + h) };
    let dist = p0.distance(&p1);
    if dist > 1e-12 {
        tol * h / dist
    } else {
        tol
    }
}

/// Clamps a discrete range to a surface's geometric range, handling periodic
/// surfaces (whose discrete range may cross the seam). Source:
/// `BRepMesh_DefaultRangeSplitter::updateRange`.
fn update_range(
    geom_first: f64,
    geom_last: f64,
    is_periodic: bool,
    d_first: &mut f64,
    d_last: &mut f64,
) {
    if *d_first < geom_first || *d_last > geom_last {
        if is_periodic {
            if (*d_last - *d_first) > (geom_last - geom_first) {
                *d_last = *d_first + (geom_last - geom_first);
            }
        } else if (*d_first < geom_last) && (*d_last > geom_first) {
            // Protection against pcurves out of the surface domain (OCCT #23675).
            if geom_first > *d_first {
                *d_first = geom_first;
            }
            if geom_last < *d_last {
                *d_last = geom_last;
            }
        }
    }
}

/// Adds `param` to the set when it lies within `range`. Returns true on insert.
fn add_param(param: f64, range: (f64, f64), params: &mut ParamSet) -> bool {
    if param < range.0 || param > range.1 {
        return false;
    }
    params.insert(param);
    true
}

/// Initializes a parameter set from breakpoint intervals, optionally splitting
/// each interval at its midpoint. Source: `initParamsFromIntervals`.
fn init_params_from_intervals(
    intervals: &[f64],
    range: (f64, f64),
    is_split: bool,
    params: &mut ParamSet,
) -> bool {
    let mut added = false;
    for i in 0..intervals.len() {
        if add_param(intervals[i], range, params) {
            added = true;
        }
        if is_split && i + 1 < intervals.len() {
            let mid = 0.5 * (intervals[i] + intervals[i + 1]);
            if add_param(mid, range, params) {
                added = true;
            }
        }
    }
    added
}

/// Breaks the range into uniform breakpoints. The interval count comes from the
/// splitter's `get_undefined_interval_nb`; when it reports no interior structure
/// the two range endpoints are returned (OCCT falls back to `UIntervals`).
fn get_undefined_interval(
    s: &dyn RangeSplitter,
    is_u: bool,
    continuity: u8,
    range: (f64, f64),
) -> Vec<f64> {
    let intervals_nb = s.get_undefined_interval_nb(is_u, continuity);
    if intervals_nb > 1 {
        let diff = (range.1 - range.0) / intervals_nb as f64;
        let breakpoints: Vec<f64> =
            (1..intervals_nb).map(|i| range.0 + i as f64 * diff).collect();
        if !breakpoints.is_empty() {
            return breakpoints;
        }
    }
    vec![range.0, range.1]
}

/// Filters a parameter set so no two samples are closer than `min_dist`, then
/// keeps samples separated by at least `filter_dist` (except a trailing buffer).
/// Source: `BRepMesh_NURBSRangeSplitter::filterParameters`.
fn filter_parameters(params: &ParamSet, min_dist: f64, filter_dist: f64) -> Vec<f64> {
    let mut arr: Vec<f64> = params.iter().copied().collect();
    if arr.is_empty() {
        return Vec::new();
    }
    // Mandatory pre-filtering using the minimal distance.
    let mut len = 1;
    for j in 1..arr.len() {
        if arr[j] - arr[len - 1] > min_dist {
            if len < j {
                arr[len] = arr[j];
            }
            len += 1;
        }
    }
    // Series filtering.
    let mut result = Vec::new();
    let mut last_added = arr[0];
    let mut last_candidate = last_added;
    let mut candidate_defined = false;
    result.push(last_added);
    let mut j = 1;
    while j < len {
        let val = arr[j];
        if val - last_added > filter_dist {
            if candidate_defined {
                last_added = last_candidate;
                candidate_defined = false;
                // Re-examine the same sample against the promoted candidate.
            } else {
                last_added = val;
                j += 1;
            }
            result.push(last_added);
        } else {
            last_candidate = val;
            candidate_defined = true;
            j += 1;
        }
    }
    result.push(arr[len - 1]);
    result
}

/// Computes the grain (filter distance) and applies it to the source parameters.
/// Source: `BRepMesh_NURBSRangeSplitter::computeGrainAndFilterParameters`.
fn compute_grain_and_filter(
    s: &dyn RangeSplitter,
    source: &ParamSet,
    tol2d: f64,
    range_diff: f64,
    delta: f64,
    params: &MeshParameters,
) -> Vec<f64> {
    let surface = match s.surface() {
        Some(sf) => sf.as_ref(),
        None => return Vec::new(),
    };
    let mut min_diff = PCONFUSION;
    if delta < 1.0 {
        min_diff /= delta;
    }
    let min_size_2d = param_resolution(surface, params.min_size, true)
        .max(param_resolution(surface, params.min_size, false));
    min_diff = min_diff.max(min_size_2d);
    let diff_max_lim = 0.1 * range_diff;
    let diff_min_lim = (0.005 * range_diff).max(2.0 * tol2d);
    let diff = min_size_2d.max(diff_max_lim.min(diff_min_lim));
    filter_parameters(source, min_diff, diff)
}

/// Shared NURBS interior-node generation. Runs through `&dyn RangeSplitter` so
/// the `get_undefined_interval_nb` / `init_parameters` hooks dispatch to the
/// concrete splitter (NURBS, Undefined, BoundaryParams, Extrusion).
///
/// `ponytail:` the `AnalyticalFilter` step that removes too-dense iso-line
/// control parameters (OCCT) is omitted; the grid is built from the filtered
/// parameter sequences directly.
fn generate_nurbs_grid(s: &dyn RangeSplitter, params: &MeshParameters) -> Option<Vec<GpPnt2d>> {
    if !s.init_parameters() {
        return None;
    }
    let surface = s.surface()?;
    let range_u = s.range_u();
    let range_v = s.range_v();
    let delta = s.delta();
    let deflection = s.deflection();
    let continuity = surface.continuity();

    let mut u_params: ParamSet = s.parameters_u().cloned().unwrap_or_default();
    let mut v_params: ParamSet = s.parameters_v().cloned().unwrap_or_default();
    let u_intervals = get_undefined_interval(s, true, continuity, range_u);
    let v_intervals = get_undefined_interval(s, false, continuity, range_v);
    init_params_from_intervals(&u_intervals, range_u, false, &mut u_params);
    init_params_from_intervals(&v_intervals, range_v, false, &mut v_params);
    // ponytail: OCCT collects boundary parameters via grabParamsOfEdges (the
    // splitter has no model access here); add the endpoints so the grid reaches
    // the boundary.
    u_params.insert(range_u.0);
    u_params.insert(range_u.1);
    v_params.insert(range_v.0);
    v_params.insert(range_v.1);

    let tol_u = param_resolution(surface.as_ref(), deflection, true);
    let tol_v = param_resolution(surface.as_ref(), deflection, false);
    let mut u_seq =
        compute_grain_and_filter(s, &u_params, tol_u, range_u.1 - range_u.0, delta.0, params);
    let mut v_seq =
        compute_grain_and_filter(s, &v_params, tol_v, range_v.1 - range_v.0, delta.1, params);

    // `AnalyticalFilter`: refine U control params along V iso-lines, then V
    // control params along U iso-lines, so the grid converges onto curvature
    // extrema (OCCT NURBSRangeSplitter::GenerateSurfaceNodes).
    let angle_interior = if params.angle_interior >= 0.0 { params.angle_interior } else { 2.0 * params.angle };
    analytical_filter_insert(
        surface.as_ref(),
        false,
        &v_seq,
        &mut u_seq,
        deflection,
        angle_interior,
        params.min_size,
    );
    analytical_filter_insert(
        surface.as_ref(),
        true,
        &u_seq,
        &mut v_seq,
        deflection,
        angle_interior,
        params.min_size,
    );

    let mut nodes = Vec::with_capacity(u_seq.len() * v_seq.len());
    for &u in &u_seq {
        for &v in &v_seq {
            nodes.push(GpPnt2d::new(u, v));
        }
    }
    Some(nodes)
}

/// Squared distance from `mid` to the segment `p1..p2`
/// (`BRepMesh_GeomTool::SquareDeflectionOfSegment`).
fn sq_deflection_of_segment(p1: &GpPnt, p2: &GpPnt, mid: &GpPnt) -> f64 {
    let ab = p2.coord.subtracted(&p1.coord);
    let len2 = ab.square_modulus();
    if len2 <= f64::EPSILON {
        return mid.coord.subtracted(&p1.coord).square_modulus();
    }
    let am = mid.coord.subtracted(&p1.coord);
    let t = (am.dot(&ab) / len2).clamp(0.0, 1.0);
    let proj = p1.coord.added(&ab.multiplied(t));
    mid.coord.subtracted(&proj).square_modulus()
}

/// `BRepMesh_NURBSRangeSplitter::AnalyticalFilter` (insertion pass): refines the
/// control-parameter sequence of one direction along the iso-lines of the other
/// direction by inserting a midpoint wherever the segment's deflection or the
/// tangent angle exceeds the face deflection / interior angle. This is what makes
/// the interior node grid converge onto the surface's curvature extrema (the
/// bbox gate). The removal/thinning pass is omitted — it is a density optimisation
/// that does not move extrema.
fn analytical_filter_insert(
    surface: &dyn Surface,
    is_iso_u: bool,
    iso_params: &[f64],
    control_params: &mut Vec<f64>,
    deflection: f64,
    angle_interior: f64,
    min_size: f64,
) {
    if control_params.len() < 2 || iso_params.is_empty() {
        return;
    }
    let sq_max_deflection = deflection * deflection;
    let sq_min_size = min_size * min_size;

    // OCCT: IsoU scans every iso-line; IsoV skips the two outer ones.
    let (start, end) = if is_iso_u {
        (0, iso_params.len())
    } else {
        (1, iso_params.len().saturating_sub(1))
    };

    for &iso_param in &iso_params[start..end] {
        let iso_point = |t: f64| -> GpPnt {
            if is_iso_u {
                surface.d0(iso_param, t)
            } else {
                surface.d0(t, iso_param)
            }
        };
        let iso_tangent = |t: f64| -> GpVec {
            let (_, du, dv) = if is_iso_u {
                surface.d1(iso_param, t)
            } else {
                surface.d1(t, iso_param)
            };
            if is_iso_u {
                dv
            } else {
                du
            }
        };

        let mut prev_param = control_params[0];
        let mut prev_pnt = iso_point(prev_param);
        let mut prev_vec = iso_tangent(prev_param);

        let mut j = 1usize;
        while j < control_params.len() {
            let curr_param = control_params[j];
            let curr_pnt = iso_point(curr_param);
            let curr_vec = iso_tangent(curr_param);

            let mid_param = 0.5 * (prev_param + curr_param);
            let mid_pnt = iso_point(mid_param);
            let sq_dist = sq_deflection_of_segment(&prev_pnt, &curr_pnt, &mid_pnt);
            let angle = prev_vec.angle(&curr_vec);

            if (sq_dist > sq_max_deflection || angle > angle_interior) && sq_dist > sq_min_size {
                control_params.insert(j, mid_param);
                // Reprocess the inserted midpoint against `prev` (j stays).
                continue;
            }

            prev_param = curr_param;
            prev_pnt = curr_pnt;
            prev_vec = curr_vec;
            j += 1;
        }
    }
}

/// Range splitter — computes the discrete UV range of a face and, for the
/// analytical / NURBS splitters, generates interior (surface) nodes.
/// Source: `BRepMesh_DefaultRangeSplitter`.
pub trait RangeSplitter {
    /// Access to the embedded base splitter state.
    fn base(&self) -> &DefaultRangeSplitter;
    /// Mutable access to the embedded base splitter state.
    fn base_mut(&mut self) -> &mut DefaultRangeSplitter;

    // ---- lifecycle ----

    /// Resets the splitter for the given discrete face. Must be called before
    /// first use. Source: `Reset`.
    fn reset(&mut self, dface: &MeshFace, params: &MeshParameters) {
        self.reset_base(dface, params);
    }

    /// Base implementation of `Reset` (used by derived splitters).
    fn reset_base(&mut self, dface: &MeshFace, params: &MeshParameters) {
        self.base_mut().reset(dface, params);
    }

    /// Registers a boundary point. Source: `AddPoint`.
    fn add_point(&mut self, point: GpPnt2d) {
        self.add_point_base(point);
    }

    /// Base implementation of `AddPoint`.
    fn add_point_base(&mut self, point: GpPnt2d) {
        self.base_mut().add_point(point);
    }

    /// Computes the discrete range from the geometric range and the registered
    /// boundary points. Source: `AdjustRange`.
    fn adjust_range(&mut self) {
        self.adjust_range_base();
    }

    /// Base implementation of `AdjustRange` (used by derived splitters).
    fn adjust_range_base(&mut self) {
        let (mut du_first, mut du_second) = self.base().range_u();
        let (mut dv_first, mut dv_second) = self.base().range_v();
        let surface = self.surface().cloned();
        let Some(surf) = surface else {
            self.base_mut().set_valid(false);
            return;
        };

        let (gu0, gu1) = surf.u_range();
        update_range(gu0, gu1, surf.is_u_periodic(), &mut du_first, &mut du_second);
        if du_second < du_first {
            self.base_mut().set_valid(false);
            return;
        }

        let (gv0, gv1) = surf.v_range();
        update_range(gv0, gv1, surf.is_v_periodic(), &mut dv_first, &mut dv_second);
        if dv_second < dv_first {
            self.base_mut().set_valid(false);
            return;
        }

        self.base_mut().set_range_u((du_first, du_second));
        self.base_mut().set_range_v((dv_first, dv_second));

        let len_u = self.base().compute_length_u(surf.as_ref());
        let len_v = self.base().compute_length_v(surf.as_ref());
        let valid = len_u > PCONFUSION && len_v > PCONFUSION;
        self.base_mut().set_valid(valid);
        if valid {
            self.compute_tolerance(len_u, len_v);
            self.compute_delta(len_u, len_v);
        }
    }

    /// True when the computed range is valid. Source: `IsValid`.
    fn is_valid(&self) -> bool {
        self.base().is_valid()
    }

    /// Scales a point between real parametric space and the face basis.
    /// Source: `Scale`.
    fn scale(&self, point: GpPnt2d, to_face_basis: bool) -> GpPnt2d {
        self.base().scale(point, to_face_basis)
    }

    /// Returns the interior nodes generated from the surface data, or `None`
    /// when the splitter generates no nodes. Source: `GenerateSurfaceNodes`.
    fn generate_surface_nodes(&self, params: &MeshParameters) -> Option<Vec<GpPnt2d>> {
        self.base().generate_surface_nodes(params)
    }

    /// 3D point corresponding to the given UV parameter. Source: `Point`.
    fn point(&self, point2d: GpPnt2d) -> GpPnt {
        self.base().point(point2d)
    }

    // ---- accessors ----

    /// Discrete U range.
    fn range_u(&self) -> (f64, f64) {
        self.base().range_u()
    }
    /// Discrete V range.
    fn range_v(&self) -> (f64, f64) {
        self.base().range_v()
    }
    /// Scale factors per direction.
    fn delta(&self) -> (f64, f64) {
        self.base().delta()
    }
    /// Parametric tolerances per direction.
    fn tolerance_uv(&self) -> (f64, f64) {
        self.base().tolerance_uv()
    }
    /// The discrete face model.
    fn dface(&self) -> Option<&MeshFace> {
        self.base().dface()
    }
    /// The face surface.
    fn surface(&self) -> Option<&Arc<dyn Surface>> {
        self.base().surface()
    }
    /// The face deflection.
    fn deflection(&self) -> f64 {
        self.base().deflection()
    }

    // ---- UV-parameter maps (UVParamRangeSplitter and derivatives) ----

    /// Sorted U parameters collected from boundary points, when the splitter
    /// tracks them.
    fn parameters_u(&self) -> Option<&ParamSet> {
        None
    }
    /// Sorted V parameters collected from boundary points, when the splitter
    /// tracks them.
    fn parameters_v(&self) -> Option<&ParamSet> {
        None
    }

    // ---- virtual hooks overridden by derived splitters ----

    /// Computes the parametric tolerances. Source: `computeTolerance`.
    fn compute_tolerance(&mut self, len_u: f64, len_v: f64) {
        self.base_mut().compute_tolerance(len_u, len_v);
    }

    /// Computes the scale factors. Source: `computeDelta`.
    fn compute_delta(&mut self, len_u: f64, len_v: f64) {
        self.base_mut().compute_delta(len_u, len_v);
    }

    /// Number of intervals to subdivide an undefined range into. Source:
    /// `getUndefinedIntervalNb` — `(isU ? NbUPoles : NbVPoles) - 1` for NURBS,
    /// approximated here from the continuity degree.
    fn get_undefined_interval_nb(&self, _is_u: bool, continuity: u8) -> i32 {
        ((continuity as i32) + 1).clamp(2, 32)
    }

    /// Whether the surface parameter sets can be initialized. Source:
    /// `initParameters`.
    fn init_parameters(&self) -> bool {
        self.surface().is_some()
    }
}

/// Base range splitter. Source: `BRepMesh_DefaultRangeSplitter`.
pub struct DefaultRangeSplitter {
    dface: Option<MeshFace>,
    surface: Option<Arc<dyn Surface>>,
    deflection: f64,
    range_u: (f64, f64),
    range_v: (f64, f64),
    delta: (f64, f64),
    tolerance: (f64, f64),
    is_valid: bool,
}

impl DefaultRangeSplitter {
    /// Creates an empty splitter.
    pub fn new() -> Self {
        Self {
            dface: None,
            surface: None,
            deflection: 0.0,
            range_u: (1e100, -1e100),
            range_v: (1e100, -1e100),
            delta: (1.0, 1.0),
            tolerance: (CONFUSION, CONFUSION),
            is_valid: true,
        }
    }

    /// Resets the splitter. Source: `Reset`.
    pub fn reset(&mut self, dface: &MeshFace, _params: &MeshParameters) {
        self.dface = Some(dface.clone());
        self.surface = dface.surface();
        self.deflection = dface.deflection();
        self.range_u = (1e100, -1e100);
        self.range_v = (1e100, -1e100);
        self.delta = (1.0, 1.0);
        self.tolerance = (CONFUSION, CONFUSION);
        self.is_valid = true;
    }

    /// Registers a boundary point. Source: `AddPoint`.
    pub fn add_point(&mut self, point: GpPnt2d) {
        self.range_u.0 = self.range_u.0.min(point.x());
        self.range_u.1 = self.range_u.1.max(point.x());
        self.range_v.0 = self.range_v.0.min(point.y());
        self.range_v.1 = self.range_v.1.max(point.y());
    }

    /// True when the computed range is valid. Source: `IsValid`.
    pub fn is_valid(&self) -> bool {
        self.is_valid
    }

    /// Scales a point between real parametric space and the face basis.
    /// Source: `Scale`.
    pub fn scale(&self, point: GpPnt2d, to_face_basis: bool) -> GpPnt2d {
        if to_face_basis {
            GpPnt2d::new(
                (point.x() - self.range_u.0) / self.delta.0,
                (point.y() - self.range_v.0) / self.delta.1,
            )
        } else {
            GpPnt2d::new(
                point.x() * self.delta.0 + self.range_u.0,
                point.y() * self.delta.1 + self.range_v.0,
            )
        }
    }

    /// The base splitter generates no interior nodes (null list in OCCT).
    pub fn generate_surface_nodes(&self, _params: &MeshParameters) -> Option<Vec<GpPnt2d>> {
        None
    }

    /// 3D point at the given UV parameter. Source: `Point`.
    pub fn point(&self, point2d: GpPnt2d) -> GpPnt {
        self.surface
            .as_ref()
            .map(|s| s.d0(point2d.x(), point2d.y()))
            .unwrap_or_else(GpPnt::zero)
    }

    /// Discrete U range.
    pub fn range_u(&self) -> (f64, f64) {
        self.range_u
    }
    /// Discrete V range.
    pub fn range_v(&self) -> (f64, f64) {
        self.range_v
    }
    /// Scale factors per direction.
    pub fn delta(&self) -> (f64, f64) {
        self.delta
    }
    /// Parametric tolerances per direction.
    pub fn tolerance_uv(&self) -> (f64, f64) {
        self.tolerance
    }
    /// The discrete face model.
    pub fn dface(&self) -> Option<&MeshFace> {
        self.dface.as_ref()
    }
    /// The face surface.
    pub fn surface(&self) -> Option<&Arc<dyn Surface>> {
        self.surface.as_ref()
    }
    /// The face deflection.
    pub fn deflection(&self) -> f64 {
        self.deflection
    }

    /// Length along U of the discrete range, sampled on a 20-segment grid.
    /// Source: `computeLengthU`.
    pub fn compute_length_u(&self, s: &dyn Surface) -> f64 {
        let (u0, u1) = self.range_u;
        let (v0, v1) = self.range_v;
        let mut long = 0.0;
        let du = 0.05 * (u1 - u0);
        let v_ave = 0.5 * (v1 + v0);
        let mut p11 = s.d0(u0, v0);
        let mut p21 = s.d0(u0, v_ave);
        let mut p31 = s.d0(u0, v1);
        let mut u = u0 + du;
        for _ in 1..=20 {
            let p12 = s.d0(u, v0);
            let p22 = s.d0(u, v_ave);
            let p32 = s.d0(u, v1);
            long += p11.distance(&p12) + p21.distance(&p22) + p31.distance(&p32);
            p11 = p12;
            p21 = p22;
            p31 = p32;
            u += du;
        }
        long / 3.0
    }

    /// Length along V of the discrete range, sampled on a 20-segment grid.
    /// Source: `computeLengthV`.
    pub fn compute_length_v(&self, s: &dyn Surface) -> f64 {
        let (u0, u1) = self.range_u;
        let (v0, v1) = self.range_v;
        let mut long = 0.0;
        let dv = 0.05 * (v1 - v0);
        let u_ave = 0.5 * (u1 + u0);
        let mut p11 = s.d0(u0, v0);
        let mut p21 = s.d0(u_ave, v0);
        let mut p31 = s.d0(u1, v0);
        let mut v = v0 + dv;
        for _ in 1..=20 {
            let p12 = s.d0(u0, v);
            let p22 = s.d0(u_ave, v);
            let p32 = s.d0(u1, v);
            long += p11.distance(&p12) + p21.distance(&p22) + p31.distance(&p32);
            p11 = p12;
            p21 = p22;
            p31 = p32;
            v += dv;
        }
        long / 3.0
    }

    /// Computes the parametric tolerances. Source: `computeTolerance`.
    pub fn compute_tolerance(&mut self, _len_u: f64, _len_v: f64) {
        let diff_u = self.range_u.1 - self.range_u.0;
        let diff_v = self.range_v.1 - self.range_v.0;
        let face_tol = self
            .dface
            .as_ref()
            .map(|f| BRepTool::face_tolerance(f.face()))
            .unwrap_or(CONFUSION);
        let res_u = self
            .surface
            .as_ref()
            .map(|s| param_resolution(s.as_ref(), face_tol, true) * 1.1)
            .unwrap_or(face_tol);
        let res_v = self
            .surface
            .as_ref()
            .map(|s| param_resolution(s.as_ref(), face_tol, false) * 1.1)
            .unwrap_or(face_tol);
        const DEFLECTION_UV: f64 = 1e-5;
        self.tolerance.0 = (DEFLECTION_UV.min(res_u)).max(1e-7 * diff_u);
        self.tolerance.1 = (DEFLECTION_UV.min(res_v)).max(1e-7 * diff_v);
    }

    /// Computes the scale factors. Source: `computeDelta`.
    pub fn compute_delta(&mut self, len_u: f64, len_v: f64) {
        let diff_u = self.range_u.1 - self.range_u.0;
        let diff_v = self.range_v.1 - self.range_v.0;
        self.delta.0 = diff_u / (if len_u < self.tolerance.0 { 1.0 } else { len_u });
        self.delta.1 = diff_v / (if len_v < self.tolerance.1 { 1.0 } else { len_v });
    }

    pub(crate) fn set_range_u(&mut self, r: (f64, f64)) {
        self.range_u = r;
    }
    pub(crate) fn set_range_v(&mut self, r: (f64, f64)) {
        self.range_v = r;
    }
    pub(crate) fn set_delta(&mut self, d: (f64, f64)) {
        self.delta = d;
    }
    pub(crate) fn set_valid(&mut self, v: bool) {
        self.is_valid = v;
    }
}

impl RangeSplitter for DefaultRangeSplitter {
    fn base(&self) -> &DefaultRangeSplitter {
        self
    }
    fn base_mut(&mut self) -> &mut DefaultRangeSplitter {
        self
    }
}

/// Cylindrical surface splitter — U-periodic seam, interior nodes along the
/// parametric grid. Source: `BRepMesh_CylinderRangeSplitter`.
pub struct CylinderRangeSplitter {
    inner: DefaultRangeSplitter,
    du: f64,
}

impl CylinderRangeSplitter {
    /// Creates an empty splitter.
    pub fn new() -> Self {
        Self {
            inner: DefaultRangeSplitter::new(),
            du: 1.0,
        }
    }
}

impl RangeSplitter for CylinderRangeSplitter {
    fn base(&self) -> &DefaultRangeSplitter {
        &self.inner
    }
    fn base_mut(&mut self) -> &mut DefaultRangeSplitter {
        &mut self.inner
    }

    fn reset(&mut self, dface: &MeshFace, params: &MeshParameters) {
        self.reset_base(dface, params);
        let r = self.surface().map(|s| cylinder_radius(s.as_ref())).unwrap_or(1.0);
        let defl = self.deflection();
        self.du = arc_angular_step(r, defl, params.angle, params.min_size);
    }

    fn compute_delta(&mut self, _len_u: f64, len_v: f64) {
        let range_v = self.base().range_v();
        self.inner
            .set_delta((self.du / len_v.max(range_v.1 - range_v.0), 1.0));
    }

    fn generate_surface_nodes(&self, _params: &MeshParameters) -> Option<Vec<GpPnt2d>> {
        let range_u = self.base().range_u();
        let range_v = self.base().range_v();
        let radius = self.surface().map(|s| cylinder_radius(s.as_ref())).unwrap_or(0.0);
        let deflection = self.deflection();

        let su = range_u.1 - range_u.0;
        let sv = range_v.1 - range_v.0;
        let a_arc_len = su * radius;
        let mut nb_u = 0i32;
        let mut nb_v = 0i32;
        if a_arc_len > deflection {
            nb_u = (su / self.du) as i32;
            // ponytail: the OCCT V-step computation is commented out, so nbV stays 0
            // and no interior rows are produced.
        }
        let du = su / (nb_u + 1) as f64;
        let dv = sv / (nb_v + 1) as f64;

        let pas_max_v = range_v.1 - dv * 0.5;
        let pas_max_u = range_u.1 - du * 0.5;
        let mut nodes = Vec::new();
        let mut pas_v = range_v.0 + dv;
        while pas_v < pas_max_v {
            let mut pas_u = range_u.0 + du;
            while pas_u < pas_max_u {
                nodes.push(GpPnt2d::new(pas_u, pas_v));
                pas_u += du;
            }
            pas_v += dv;
        }
        Some(nodes)
    }
}

/// Conical surface splitter. Source: `BRepMesh_ConeRangeSplitter`.
pub struct ConeRangeSplitter {
    inner: DefaultRangeSplitter,
}

impl ConeRangeSplitter {
    /// Creates an empty splitter.
    pub fn new() -> Self {
        Self {
            inner: DefaultRangeSplitter::new(),
        }
    }

    /// Returns the split steps along U and V and the number of steps.
    /// Source: `GetSplitSteps`.
    pub fn get_split_steps(
        &self,
        params: &MeshParameters,
        steps_nb: &mut (i32, i32),
    ) -> (f64, f64) {
        let range_u = self.base().range_u();
        let range_v = self.base().range_v();
        let surface = self.surface().unwrap();
        let deflection = self.deflection();
        // ponytail: OCCT computes aRadius from gp_Cone RefRadius/SemiAngle; here it is
        // measured from the surface at the two V range ends.
        let a_radius = cone_radius_at(surface.as_ref(), range_v.0)
            .max(cone_radius_at(surface.as_ref(), range_v.1));

        let mut du = arc_angular_step(a_radius, deflection, params.angle, params.min_size);

        let a_diff_u = range_u.1 - range_u.0;
        let a_diff_v = range_v.1 - range_v.0;
        let a_scale = du * a_radius;
        let a_ratio = (a_diff_v / a_scale).ln().max(1.0);
        let nb_u = (a_diff_u / du) as i32;
        let nb_v = (a_diff_v / a_scale / a_ratio) as i32;

        du = a_diff_u / (nb_u + 1) as f64;
        let dv = a_diff_v / (nb_v + a_ratio as i32) as f64;

        steps_nb.0 = nb_u;
        steps_nb.1 = nb_v;
        (du, dv)
    }
}

impl RangeSplitter for ConeRangeSplitter {
    fn base(&self) -> &DefaultRangeSplitter {
        &self.inner
    }
    fn base_mut(&mut self) -> &mut DefaultRangeSplitter {
        &mut self.inner
    }

    fn generate_surface_nodes(&self, params: &MeshParameters) -> Option<Vec<GpPnt2d>> {
        let range_u = self.base().range_u();
        let range_v = self.base().range_v();
        let mut steps_nb = (0i32, 0i32);
        let (du, dv) = self.get_split_steps(params, &mut steps_nb);

        let pas_max_v = range_v.1 - dv * 0.5;
        let pas_max_u = range_u.1 - du * 0.5;
        let mut nodes = Vec::new();
        let mut pas_v = range_v.0 + dv;
        while pas_v < pas_max_v {
            let mut pas_u = range_u.0 + du;
            while pas_u < pas_max_u {
                nodes.push(GpPnt2d::new(pas_u, pas_v));
                pas_u += du;
            }
            pas_v += dv;
        }
        Some(nodes)
    }
}

/// Spherical surface splitter — staggered U/V grid. Source:
/// `BRepMesh_SphereRangeSplitter`.
pub struct SphereRangeSplitter {
    inner: DefaultRangeSplitter,
}

impl SphereRangeSplitter {
    /// Creates an empty splitter.
    pub fn new() -> Self {
        Self {
            inner: DefaultRangeSplitter::new(),
        }
    }

    /// Computes the step and upper bound for a range. Source: `computeStep`.
    fn compute_step(&self, range: (f64, f64), default_step: f64) -> (f64, f64) {
        let diff = range.1 - range.0;
        let step = diff / ((diff / default_step) as i32 + 1) as f64;
        (step, range.1 - PCONFUSION)
    }
}

impl RangeSplitter for SphereRangeSplitter {
    fn base(&self) -> &DefaultRangeSplitter {
        &self.inner
    }
    fn base_mut(&mut self) -> &mut DefaultRangeSplitter {
        &mut self.inner
    }

    fn generate_surface_nodes(&self, params: &MeshParameters) -> Option<Vec<GpPnt2d>> {
        let range_v = self.base().range_v();
        let range_u = self.base().range_u();
        let radius = self.surface().map(|s| sphere_radius(s.as_ref())).unwrap_or(0.0);
        let deflection = self.deflection();
        let a_step = 0.7 * arc_angular_step(radius, deflection, params.angle, params.min_size);

        let (step_v, max_v) = self.compute_step(range_v, a_step);
        let (step_u, max_u) = self.compute_step(range_u, a_step);
        let half_du = step_u * 0.5;

        let mut shift = false;
        let mut nodes = Vec::new();
        let mut pas_v = range_v.0 + step_v;
        while pas_v < max_v {
            shift = !shift;
            let d = if shift { half_du } else { 0.0 };
            let mut pas_u = range_u.0 + d;
            while pas_u < max_u {
                nodes.push(GpPnt2d::new(pas_u, pas_v));
                pas_u += step_u;
            }
            pas_v += step_v;
        }
        Some(nodes)
    }
}

/// UV range splitter — tracks the U/V parameters of boundary points. Source:
/// `BRepMesh_UVParamRangeSplitter`.
pub struct UVParamRangeSplitter {
    inner: DefaultRangeSplitter,
    u_params: ParamSet,
    v_params: ParamSet,
}

impl UVParamRangeSplitter {
    /// Creates an empty splitter.
    pub fn new() -> Self {
        Self {
            inner: DefaultRangeSplitter::new(),
            u_params: ParamSet::new(),
            v_params: ParamSet::new(),
        }
    }
}

impl RangeSplitter for UVParamRangeSplitter {
    fn base(&self) -> &DefaultRangeSplitter {
        &self.inner
    }
    fn base_mut(&mut self) -> &mut DefaultRangeSplitter {
        &mut self.inner
    }

    fn reset(&mut self, dface: &MeshFace, params: &MeshParameters) {
        self.reset_base(dface, params);
        self.u_params.clear();
        self.v_params.clear();
    }

    fn parameters_u(&self) -> Option<&ParamSet> {
        Some(&self.u_params)
    }
    fn parameters_v(&self) -> Option<&ParamSet> {
        Some(&self.v_params)
    }
}

/// Torus surface splitter — U/V periodic, boundary-parameter aware grid.
/// Source: `BRepMesh_TorusRangeSplitter`.
pub struct TorusRangeSplitter {
    inner: UVParamRangeSplitter,
}

impl TorusRangeSplitter {
    /// Creates an empty splitter.
    pub fn new() -> Self {
        Self {
            inner: UVParamRangeSplitter::new(),
        }
    }

    /// Fills a parameter sequence from the collected params, spaced at least
    /// `aStdStep` apart. Source: `fillParams`.
    fn fill_params(
        &self,
        params: &ParamSet,
        range: (f64, f64),
        steps_nb: i32,
        scale: f64,
    ) -> Vec<f64> {
        let mut arr: Vec<f64> = params.iter().copied().collect();
        let diff = (range.1 - range.0).abs();
        let mut step = calc_average_duv(&mut arr);
        step = step.max(diff / steps_nb as f64 / 2.0);

        let mut std_step = if arr.is_empty() { 0.0 } else { diff / arr.len() as f64 };
        if step > std_step {
            std_step = step;
        }
        std_step *= scale;

        let mut result = Vec::new();
        for &pp in &arr {
            let is_to_insert = result.iter().all(|&v: &f64| (v - pp).abs() > std_step);
            if is_to_insert {
                result.push(pp);
            }
        }
        result
    }
}

/// Average gap between consecutive sorted parameters. Source:
/// `FUN_CalcAverageDUV`.
fn calc_average_duv(p: &mut [f64]) -> f64 {
    p.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mut n = 0;
    let mut result = 0.0;
    for i in 1..p.len() {
        let d = (p[i] - p[i - 1]).abs();
        if d > 1e-7 {
            result += d;
            n += 1;
        }
    }
    if n > 0 {
        result / n as f64
    } else {
        -1.0
    }
}

impl RangeSplitter for TorusRangeSplitter {
    fn base(&self) -> &DefaultRangeSplitter {
        &self.inner.inner
    }
    fn base_mut(&mut self) -> &mut DefaultRangeSplitter {
        &mut self.inner.inner
    }

    fn reset(&mut self, dface: &MeshFace, params: &MeshParameters) {
        self.inner.reset(dface, params);
    }

    fn add_point(&mut self, point: GpPnt2d) {
        self.add_point_base(point);
        self.inner.u_params.insert(point.x());
        self.inner.v_params.insert(point.y());
    }

    fn parameters_u(&self) -> Option<&ParamSet> {
        Some(&self.inner.u_params)
    }
    fn parameters_v(&self) -> Option<&ParamSet> {
        Some(&self.inner.v_params)
    }

    fn generate_surface_nodes(&self, params: &MeshParameters) -> Option<Vec<GpPnt2d>> {
        let range_u = self.base().range_u();
        let range_v = self.base().range_v();
        let diff_u = range_u.1 - range_u.0;
        let diff_v = range_v.1 - range_v.0;

        let (r_major, r_minor) = self
            .surface()
            .map(|s| torus_radii(s.as_ref()))
            .unwrap_or((0.0, 0.0));
        let deflection = self.deflection();
        let r = r_minor;
        let R = r_major;

        let old_dv = arc_angular_step(r, deflection, params.angle, params.min_size);
        let dv = old_dv;

        let nb_v = (diff_v / dv) as i32;
        let nb_v = nb_v.max(2);
        let dv = diff_v / (nb_v + 1) as f64;

        let ru = R + r;
        let du = if ru > 1e-16 {
            let du0 = arc_angular_step(ru, deflection, params.angle, params.min_size);
            let aa = (du0 * du0 + old_dv * old_dv).sqrt();
            if aa < RESOLUTION {
                return None;
            }
            du0 * old_dv.min(du0) / aa
        } else {
            dv
        };

        let mut nb_u = (diff_u / du) as i32;
        nb_u = nb_u.max(2);
        let ratio_terms = if diff_v * r != 0.0 {
            nb_v as f64 * diff_u * R / (diff_v * r) / 5.0
        } else {
            0.0
        };
        nb_u = nb_u.max(ratio_terms as i32);
        let du = diff_u / (nb_u + 1) as f64;

        let param_u = if R < r {
            (0..=nb_u).map(|i| range_u.0 + i as f64 * du).collect()
        } else {
            self.fill_params(&self.inner.u_params, range_u, nb_u, 0.5)
        };
        let param_v = self.fill_params(&self.inner.v_params, range_v, nb_v, 2.0 / 3.0);

        let new_range_u = (range_u.0 + du * 0.1, range_u.1 - du * 0.1);
        let new_range_v = (range_v.0 + dv * 0.1, range_v.1 - dv * 0.1);

        let mut nodes = Vec::new();
        for &pas_u in &param_u {
            if pas_u >= new_range_u.0 && pas_u < new_range_u.1 {
                for &pas_v in &param_v {
                    if pas_v >= new_range_v.0 && pas_v < new_range_v.1 {
                        nodes.push(GpPnt2d::new(pas_u, pas_v));
                    }
                }
            }
        }
        Some(nodes)
    }
}

/// NURBS / Bezier / BSpline surface splitter — interval-based interior grid.
/// Source: `BRepMesh_NURBSRangeSplitter`.
pub struct NURBSRangeSplitter {
    inner: UVParamRangeSplitter,
    surface_type: SurfaceType,
}

impl NURBSRangeSplitter {
    /// Creates an empty splitter.
    pub fn new() -> Self {
        Self {
            inner: UVParamRangeSplitter::new(),
            surface_type: SurfaceType::OtherSurface,
        }
    }
}

impl RangeSplitter for NURBSRangeSplitter {
    fn base(&self) -> &DefaultRangeSplitter {
        &self.inner.inner
    }
    fn base_mut(&mut self) -> &mut DefaultRangeSplitter {
        &mut self.inner.inner
    }

    fn adjust_range(&mut self) {
        self.adjust_range_base();
        self.surface_type = self
            .surface()
            .map(|s| classify_surface(s.as_ref()))
            .unwrap_or(SurfaceType::OtherSurface);
        if self.surface_type == SurfaceType::BezierSurface {
            let (ru0, ru1) = self.base().range_u();
            let (rv0, rv1) = self.base().range_v();
            self.base_mut()
                .set_valid(ru0 >= -0.5 && ru1 <= 1.5 && rv0 >= -0.5 && rv1 <= 1.5);
        }
    }

    fn generate_surface_nodes(&self, params: &MeshParameters) -> Option<Vec<GpPnt2d>> {
        generate_nurbs_grid(self, params)
    }

    fn parameters_u(&self) -> Option<&ParamSet> {
        Some(&self.inner.u_params)
    }
    fn parameters_v(&self) -> Option<&ParamSet> {
        Some(&self.inner.v_params)
    }

    fn get_undefined_interval_nb(&self, _is_u: bool, continuity: u8) -> i32 {
        ((continuity as i32) + 1).clamp(2, 32)
    }
}

/// Splitter for surfaces that look like NURBS but expose no poles or other
/// interval characteristics — a single interval per direction. Source:
/// `BRepMesh_UndefinedRangeSplitter`.
pub struct UndefinedRangeSplitter {
    inner: NURBSRangeSplitter,
}

impl UndefinedRangeSplitter {
    /// Creates an empty splitter.
    pub fn new() -> Self {
        Self {
            inner: NURBSRangeSplitter::new(),
        }
    }
}

impl RangeSplitter for UndefinedRangeSplitter {
    fn base(&self) -> &DefaultRangeSplitter {
        &self.inner.inner.inner
    }
    fn base_mut(&mut self) -> &mut DefaultRangeSplitter {
        &mut self.inner.inner.inner
    }

    fn adjust_range(&mut self) {
        self.inner.adjust_range();
    }

    fn generate_surface_nodes(&self, params: &MeshParameters) -> Option<Vec<GpPnt2d>> {
        generate_nurbs_grid(self, params)
    }

    fn parameters_u(&self) -> Option<&ParamSet> {
        Some(&self.inner.inner.u_params)
    }
    fn parameters_v(&self) -> Option<&ParamSet> {
        Some(&self.inner.inner.v_params)
    }

    fn get_undefined_interval_nb(&self, _is_u: bool, _continuity: u8) -> i32 {
        1
    }
}

/// Splitter that seeds the U/V parameter sets from boundary points only. Source:
/// `BRepMesh_BoundaryParamsRangeSplitter`.
pub struct BoundaryParamsRangeSplitter {
    inner: NURBSRangeSplitter,
}

impl BoundaryParamsRangeSplitter {
    /// Creates an empty splitter.
    pub fn new() -> Self {
        Self {
            inner: NURBSRangeSplitter::new(),
        }
    }
}

impl RangeSplitter for BoundaryParamsRangeSplitter {
    fn base(&self) -> &DefaultRangeSplitter {
        &self.inner.inner.inner
    }
    fn base_mut(&mut self) -> &mut DefaultRangeSplitter {
        &mut self.inner.inner.inner
    }

    fn reset(&mut self, dface: &MeshFace, params: &MeshParameters) {
        self.reset_base(dface, params);
        // The derived boundary/torus splitters seed their own UV maps from
        // AddPoint; a fresh Reset must clear them (the base reset below only
        // resets the bottom DefaultRangeSplitter).
        self.inner.inner.u_params.clear();
        self.inner.inner.v_params.clear();
    }

    fn add_point(&mut self, point: GpPnt2d) {
        self.add_point_base(point);
        self.inner.inner.u_params.insert(point.x());
        self.inner.inner.v_params.insert(point.y());
    }

    fn adjust_range(&mut self) {
        self.inner.adjust_range();
    }

    fn generate_surface_nodes(&self, params: &MeshParameters) -> Option<Vec<GpPnt2d>> {
        generate_nurbs_grid(self, params)
    }

    fn parameters_u(&self) -> Option<&ParamSet> {
        Some(&self.inner.inner.u_params)
    }
    fn parameters_v(&self) -> Option<&ParamSet> {
        Some(&self.inner.inner.v_params)
    }

    /// Source: `BRepMesh_BoundaryParamsRangeSplitter::initParameters`.
    fn init_parameters(&self) -> bool {
        true
    }
}

/// Splitter for extrusion surfaces — the interval count follows the basis
/// curve, simplified here to a single interval. Source:
/// `BRepMesh_ExtrusionRangeSplitter`.
pub struct ExtrusionRangeSplitter {
    inner: NURBSRangeSplitter,
}

impl ExtrusionRangeSplitter {
    /// Creates an empty splitter.
    pub fn new() -> Self {
        Self {
            inner: NURBSRangeSplitter::new(),
        }
    }
}

impl RangeSplitter for ExtrusionRangeSplitter {
    fn base(&self) -> &DefaultRangeSplitter {
        &self.inner.inner.inner
    }
    fn base_mut(&mut self) -> &mut DefaultRangeSplitter {
        &mut self.inner.inner.inner
    }

    fn adjust_range(&mut self) {
        self.inner.adjust_range();
    }

    fn generate_surface_nodes(&self, params: &MeshParameters) -> Option<Vec<GpPnt2d>> {
        generate_nurbs_grid(self, params)
    }

    fn parameters_u(&self) -> Option<&ParamSet> {
        Some(&self.inner.inner.u_params)
    }
    fn parameters_v(&self) -> Option<&ParamSet> {
        Some(&self.inner.inner.v_params)
    }

    fn get_undefined_interval_nb(&self, _is_u: bool, _continuity: u8) -> i32 {
        1
    }
}

/// Creates the range splitter matching a surface's analytic type. Source:
/// `BRepMesh_MeshAlgoFactory` splitter selection (the `DeflectionControlMeshAlgo`
/// / `NodeInsertionMeshAlgo` wrapper choice is out of scope here).
pub fn create_range_splitter(surface: &dyn Surface) -> Box<dyn RangeSplitter> {
    match classify_surface(surface) {
        SurfaceType::Plane => Box::new(DefaultRangeSplitter::new()),
        SurfaceType::Sphere => Box::new(SphereRangeSplitter::new()),
        SurfaceType::Cylinder => Box::new(CylinderRangeSplitter::new()),
        SurfaceType::Cone => Box::new(ConeRangeSplitter::new()),
        SurfaceType::Torus => Box::new(TorusRangeSplitter::new()),
        SurfaceType::SurfaceOfRevolution => Box::new(BoundaryParamsRangeSplitter::new()),
        SurfaceType::SurfaceOfExtrusion => Box::new(ExtrusionRangeSplitter::new()),
        SurfaceType::BezierSurface | SurfaceType::BSplineSurface => Box::new(NURBSRangeSplitter::new()),
        SurfaceType::OffsetSurface | SurfaceType::OtherSurface => {
            Box::new(UndefinedRangeSplitter::new())
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;

    use occt_core::gp::{GpAx3, GpCylinder, GpPln, GpSphere, GpTorus};
    use occt_geom::{
        bspline_surface::{bspline_surface_uniform_knots, GeomBSplineSurface},
        GeomCylinder, GeomPlane, GeomSphere, GeomTorus,
    };

    use crate::builder::TopoBuilder;
    use crate::tgeometry::GeometryRegistry;

    /// Releases registry entries for a face so tests don't leave stale geometry.
    fn clear_face(mf: &MeshFace) {
        GeometryRegistry::global().clear_shape(&mf.face().0);
    }

    /// Builds a discrete face whose registered plane surface is replaced by the
    /// given analytic surface; deflection fixed to 0.001.
    fn make_dface(surface: Arc<dyn Surface>) -> MeshFace {
        let b = TopoBuilder::new();
        let face = b.make_face_plane(&GpPln::new(GpAx3::standard()));
        let mut mf = MeshFace::new(face);
        mf.set_surface(Some(surface));
        mf.set_deflection(0.001);
        mf
    }

    fn plane_surface() -> Arc<dyn Surface> {
        Arc::new(GeomPlane::new(GpPln::new(GpAx3::standard())))
    }

    fn cylinder_surface(r: f64) -> Arc<dyn Surface> {
        Arc::new(GeomCylinder::new(GpCylinder::new(GpAx3::standard(), r).unwrap()))
    }

    fn sphere_surface(r: f64) -> Arc<dyn Surface> {
        Arc::new(GeomSphere::new(GpSphere::new(GpAx3::standard(), r).unwrap()))
    }

    fn torus_surface(major: f64, minor: f64) -> Arc<dyn Surface> {
        Arc::new(GeomTorus::new(GpTorus::new(GpAx3::standard(), major, minor).unwrap()))
    }

    #[test]
    fn classify_analytic_surfaces() {
        assert_eq!(classify_surface(plane_surface().as_ref()), SurfaceType::Plane);
        assert_eq!(classify_surface(cylinder_surface(1.0).as_ref()), SurfaceType::Cylinder);
        assert_eq!(classify_surface(sphere_surface(1.0).as_ref()), SurfaceType::Sphere);
        assert_eq!(classify_surface(torus_surface(2.0, 1.0).as_ref()), SurfaceType::Torus);

        let (ku, kv) = bspline_surface_uniform_knots(4, 4, 3, 3);
        let poles = (0..4)
            .map(|i| (0..4).map(|j| GpPnt::new(i as f64, j as f64, 0.0)).collect())
            .collect();
        let bs: Arc<dyn Surface> = Arc::new(GeomBSplineSurface::new(poles, ku, kv, 3, 3).unwrap());
        assert_eq!(classify_surface(bs.as_ref()), SurfaceType::BSplineSurface);
    }

    #[test]
    fn factory_dispatches_by_surface_type() {
        let _ = create_range_splitter(plane_surface().as_ref());
        let _ = create_range_splitter(cylinder_surface(1.0).as_ref());
        let _ = create_range_splitter(sphere_surface(1.0).as_ref());
        let _ = create_range_splitter(torus_surface(2.0, 1.0).as_ref());
    }

    #[test]
    fn plane_range_splitter_computes_ranges_and_scale() {
        let dface = make_dface(plane_surface());
        let params = MeshParameters::default();
        let mut sp = DefaultRangeSplitter::new();
        sp.reset(&dface, &params);
        sp.add_point(GpPnt2d::new(0.0, 0.0));
        sp.add_point(GpPnt2d::new(1.0, 0.0));
        sp.add_point(GpPnt2d::new(0.0, 1.0));
        sp.adjust_range();
        assert!(sp.is_valid());

        let (u0, u1) = sp.range_u();
        let (v0, v1) = sp.range_v();
        assert!((u0 - 0.0).abs() < 1e-9 && (u1 - 1.0).abs() < 1e-9);
        assert!((v0 - 0.0).abs() < 1e-9 && (v1 - 1.0).abs() < 1e-9);

        // Face-basis scaling maps range origin to ~(0, 0).
        let s = sp.scale(GpPnt2d::new(u0, v0), true);
        assert!(s.x().abs() < 1e-9 && s.y().abs() < 1e-9);
        let back = sp.scale(GpPnt2d::new(0.0, 0.0), false);
        assert!((back.x() - u0).abs() < 1e-9 && (back.y() - v0).abs() < 1e-9);

        // The base splitter generates no interior nodes.
        assert!(sp.generate_surface_nodes(&params).is_none());

        // 3D point evaluation on the Z=0 plane.
        let p = sp.point(GpPnt2d::new(0.5, 0.5));
        assert!((p.z() - 0.0).abs() < 1e-9);
        clear_face(&dface);
    }

    #[test]
    fn cylinder_seam_range_clamps_to_period() {
        let dface = make_dface(cylinder_surface(1.0));
        let params = MeshParameters::default();
        let mut sp = CylinderRangeSplitter::new();
        sp.reset(&dface, &params);

        // Boundary points spanning more than one U-period cross the seam: the
        // discrete range must be clamped to a single period.
        sp.add_point(GpPnt2d::new(-0.2, 0.0));
        sp.add_point(GpPnt2d::new(6.5, 0.0));
        sp.add_point(GpPnt2d::new(0.0, -1.0));
        sp.add_point(GpPnt2d::new(1.0, 1.0));
        sp.adjust_range();
        assert!(sp.is_valid());

        let (u0, u1) = sp.range_u();
        let period = 2.0 * PI;
        assert!(u1 - u0 <= period + 1e-9, "seam range {u0}..{u1} exceeds one period");
        assert!((u1 - u0 - period).abs() < 1e-9);
        assert!((u0 - (-0.2)).abs() < 1e-9);

        // Cylinder delta: first component from the angular step, V delta = 1.
        let (du, dv) = sp.delta();
        assert!((dv - 1.0).abs() < 1e-12);
        assert!(du > 0.0);

        // Faithful to OCCT: the V-step code is commented out, so no interior rows.
        let nodes = sp.generate_surface_nodes(&params).unwrap_or_default();
        assert!(nodes.iter().all(|p| p.x() >= u0 - 1e-9 && p.x() <= u1 + 1e-9));
        assert!(nodes.iter().all(|p| p.y() >= -1.0 - 1e-9 && p.y() <= 1.0 + 1e-9));
        clear_face(&dface);
    }

    #[test]
    fn sphere_range_generates_staggered_nodes() {
        let dface = make_dface(sphere_surface(1.0));
        let params = MeshParameters::default();
        let mut sp = SphereRangeSplitter::new();
        sp.reset(&dface, &params);
        // Full sphere range.
        sp.add_point(GpPnt2d::new(0.0, -PI * 0.5));
        sp.add_point(GpPnt2d::new(2.0 * PI, PI * 0.5));
        sp.adjust_range();
        assert!(sp.is_valid());

        let nodes = sp.generate_surface_nodes(&params).expect("sphere nodes");
        assert!(nodes.len() > 100, "expected a dense staggered grid, got {}", nodes.len());

        let (u0, u1) = sp.range_u();
        let (v0, v1) = sp.range_v();
        for p in &nodes {
            assert!(p.x() >= u0 - 1e-6 && p.x() <= u1 + 1e-6);
            assert!(p.y() >= v0 - 1e-6 && p.y() <= v1 + 1e-6);
        }
        // Staggered rows: consecutive V rows start at the range origin and at
        // origin + half a U step, alternating — so at least two distinct
        // per-row minimum U values must exist.
        let mut map = std::collections::BTreeMap::<i64, f64>::new();
        for p in &nodes {
            let key = (p.y() * 1e6).round() as i64;
            let e = map.entry(key).or_insert(p.x());
            if p.x() < *e {
                *e = p.x();
            }
        }
        let mut row_mins: Vec<f64> = map.values().copied().collect();
        row_mins.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let staggered = row_mins.windows(2).any(|w| (w[1] - w[0]).abs() > 1e-3);
        assert!(staggered, "sphere rows should be staggered by half a step");
        clear_face(&dface);
    }

    #[test]
    fn torus_range_generates_nodes_within_range() {
        let dface = make_dface(torus_surface(2.0, 1.0));
        let params = MeshParameters::default();
        let mut sp = TorusRangeSplitter::new();
        sp.reset(&dface, &params);
        // A quarter of the torus in U, full V, with a dense boundary (as a
        // discretized face boundary would supply) so the density filter keeps
        // interior samples.
        for i in 0..=16 {
            let u = PI * 0.5 * i as f64 / 16.0;
            for j in 0..=32 {
                let v = 2.0 * PI * j as f64 / 32.0;
                sp.add_point(GpPnt2d::new(u, v));
            }
        }
        sp.adjust_range();
        assert!(sp.is_valid());

        let (u0, u1) = sp.range_u();
        let (v0, v1) = sp.range_v();
        let nodes = sp.generate_surface_nodes(&params).expect("torus nodes");
        assert!(nodes.len() > 100, "expected interior nodes, got {}", nodes.len());
        for p in &nodes {
            assert!(p.x() >= u0 - 1e-6 && p.x() <= u1 + 1e-6);
            assert!(p.y() >= v0 - 1e-6 && p.y() <= v1 + 1e-6);
        }
        clear_face(&dface);
    }

    #[test]
    fn nurbs_range_generates_grid_within_bounds() {
        let (ku, kv) = bspline_surface_uniform_knots(4, 4, 3, 3);
        let poles = (0..4)
            .map(|i| (0..4).map(|j| GpPnt::new(i as f64, j as f64, 0.0)).collect())
            .collect();
        let bs: Arc<dyn Surface> = Arc::new(GeomBSplineSurface::new(poles, ku, kv, 3, 3).unwrap());
        let dface = make_dface(bs);
        let params = MeshParameters::default();

        let mut sp = NURBSRangeSplitter::new();
        sp.reset(&dface, &params);
        sp.add_point(GpPnt2d::new(0.0, 0.0));
        sp.add_point(GpPnt2d::new(1.0, 1.0));
        sp.adjust_range();
        assert!(sp.is_valid());

        let (u0, u1) = sp.range_u();
        let (v0, v1) = sp.range_v();
        assert!((u0 - 0.0).abs() < 1e-9 && (u1 - 1.0).abs() < 1e-9);
        assert!((v0 - 0.0).abs() < 1e-9 && (v1 - 1.0).abs() < 1e-9);

        let nodes = sp.generate_surface_nodes(&params).expect("nurbs nodes");
        assert!(nodes.len() >= 4, "expected a grid, got {}", nodes.len());
        for p in &nodes {
            assert!(p.x() >= u0 - 1e-6 && p.x() <= u1 + 1e-6);
            assert!(p.y() >= v0 - 1e-6 && p.y() <= v1 + 1e-6);
        }
        clear_face(&dface);
    }

    #[test]
    fn boundary_params_splitter_collects_uv_params() {
        let dface = make_dface(plane_surface());
        let params = MeshParameters::default();
        // The plain UV splitter does not seed its maps from AddPoint (faithful
        // to OCCT); only the boundary/torus splitters do.
        let mut sp = UVParamRangeSplitter::new();
        sp.reset(&dface, &params);
        sp.add_point(GpPnt2d::new(0.5, 0.25));
        assert!(sp.parameters_u().unwrap().is_empty());

        let mut bp = BoundaryParamsRangeSplitter::new();
        bp.reset(&dface, &params);
        bp.add_point(GpPnt2d::new(0.5, 0.25));
        bp.add_point(GpPnt2d::new(0.5, 0.75));
        bp.add_point(GpPnt2d::new(1.5, 0.25));
        assert!(bp.parameters_u().unwrap().contains(&0.5));
        assert!(bp.parameters_u().unwrap().contains(&1.5));
        assert_eq!(bp.parameters_v().unwrap().len(), 2);
        // Reset clears the collected parameters.
        bp.reset(&dface, &params);
        assert!(bp.parameters_u().unwrap().is_empty());
        clear_face(&dface);
    }
}
