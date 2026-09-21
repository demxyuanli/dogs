use super::prelude::*;
use super::*;
use occt_geom::Curve;

/// Sorted set of real parameters used to track boundary U/V values. Mirrors
/// OCCT's `IMeshData::IMapOfReal`; `f64` does not implement `Ord` (NaN), so the
/// values are kept in a sorted `Vec` ordered by `partial_cmp`.
#[derive(Debug, Clone, Default)]

pub struct ParamSet {
    pub(super) values: Vec<f64>,
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
/// Exact transcription of `GeomAdaptor_Surface::Load`
/// (`GeomAdaptor_Surface.cxx:422-513`), which compares **`DynamicType`**
/// (exact class, not `IsKind`) in a fixed order: rectangular-trimmed (recurse
/// on the basis), Plane, Cylinder, Cone, Sphere, Torus, SurfaceOfRevolution,
/// SurfaceOfLinearExtrusion, **Bezier** (`cxx:480`, before BSpline),
/// BSpline, Offset, else `GeomAbs_OtherSurface`.
///
/// The previous body guessed Plane/Cylinder/Cone/Sphere/Torus from the
/// periodicity flags plus a sampled constant-radius test (`is_cylinder_like`),
/// which turned any V-bounded cylinder or cone into a `Sphere` and swapped the
/// whole internal node grid (audit A24).
pub fn classify_surface(s: &dyn Surface) -> SurfaceType {
    if let Some(basis) = s.rectangular_trimmed_basis() {
        return classify_surface(basis.as_ref());
    }
    if s.gp_pln().is_some() {
        return SurfaceType::Plane;
    }
    if s.gp_cylinder().is_some() {
        return SurfaceType::Cylinder;
    }
    if s.gp_cone().is_some() {
        return SurfaceType::Cone;
    }
    if s.gp_sphere().is_some() {
        return SurfaceType::Sphere;
    }
    if s.gp_torus().is_some() {
        return SurfaceType::Torus;
    }
    if s.is_surface_of_revolution() {
        return SurfaceType::SurfaceOfRevolution;
    }
    if s.is_surface_of_linear_extrusion() {
        return SurfaceType::SurfaceOfExtrusion;
    }
    if s.is_bezier_surface() {
        return SurfaceType::BezierSurface;
    }
    if s.is_bspline_surface() {
        return SurfaceType::BSplineSurface;
    }
    if s.is_offset_surface() {
        return SurfaceType::OffsetSurface;
    }
    SurfaceType::OtherSurface
}

/// Radius of a cylinder. `BRepAdaptor_Surface::Cylinder().Radius()` when the
/// face is a `Geom_CylindricalSurface`; otherwise the chord `(0,0)`--`(π,0)`.
pub(super) fn cylinder_radius(s: &dyn Surface) -> f64 {
    if let Some(basis) = s.rectangular_trimmed_basis() {
        return cylinder_radius(basis.as_ref());
    }
    if let Some(c) = s.gp_cylinder() {
        return c.radius();
    }
    s.d0(0.0, 0.0).distance(&s.d0(PI, 0.0)) * 0.5
}

/// Radius of a sphere. `BRepAdaptor_Surface::Sphere().Radius()` when the face
/// is a `Geom_SphericalSurface`; otherwise the chord `(0,0)`--`(π,0)`.
pub(super) fn sphere_radius(s: &dyn Surface) -> f64 {
    if let Some(basis) = s.rectangular_trimmed_basis() {
        return sphere_radius(basis.as_ref());
    }
    if let Some(sph) = s.gp_sphere() {
        return sph.radius();
    }
    s.d0(0.0, 0.0).distance(&s.d0(PI, 0.0)) * 0.5
}

/// `(major, minor)` radii. `BRepAdaptor_Surface::Torus()` when the face is a
/// `Geom_ToroidalSurface`; otherwise the same chord measure as before.
pub(super) fn torus_radii(s: &dyn Surface) -> (f64, f64) {
    if let Some(t) = s.gp_torus() {
        return (t.major_radius(), t.minor_radius());
    }
    let minor = s.d0(0.0, 0.0).distance(&s.d0(0.0, PI)) * 0.5;
    let major_plus_minor = s.d0(0.0, 0.0).distance(&s.d0(PI, 0.0)) * 0.5;
    (major_plus_minor - minor, minor)
}

/// Radius of the U-circle of a cone at parameter `v` measured from the surface.
pub(super) fn cone_radius_at(s: &dyn Surface, v: f64) -> f64 {
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

/// Parametric step along U or V that corresponds to a 3D displacement of `tol`.
/// BSpline uses `Geom_BSplineSurface::Resolution` (`cxx:2197-2221`,
/// `GeomAdaptor_Surface.cxx:1877-1880`). Other types keep the FD estimate.
pub(super) fn param_resolution(s: &dyn Surface, tol: f64, is_u: bool) -> f64 {
    if let Some((ru, rv)) = s.uv_resolution(tol) {
        return if is_u { ru } else { rv };
    }
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
pub(super) fn update_range(
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
pub(super) fn add_param(param: f64, range: (f64, f64), params: &mut ParamSet) -> bool {
    if param < range.0 || param > range.1 {
        return false;
    }
    params.insert(param);
    true
}

/// Initializes a parameter set from breakpoint intervals, optionally splitting
/// each interval at its midpoint. Source: `initParamsFromIntervals`.
pub(super) fn init_params_from_intervals(
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

/// `GeomAdaptor_Curve::NbIntervals` / `Intervals` on `[range]` (`cxx:371-413`).
/// `C0` or `S <= Continuity` is one span. Otherwise keep unique knots that
/// lie inside the trimmed range (CN uses `aCont = degree`, every knot breaks).
fn adaptor_curve_intervals(curve: &dyn Curve, continuity: u8, range: (f64, f64)) -> Vec<f64> {
    // `GeomAdaptor_Curve::NbIntervals` (`cxx:371-413`): BSpline uses
    // `LocalContinuity(First, Last)` then `BSplCLib::Intervals`.
    if let (Some(knots), Some(degree)) = (curve.bspline_knots(), curve.nurbs_degree()) {
        return occt_core::bspl::adaptor_intervals(
            knots,
            degree,
            curve.is_periodic(),
            continuity,
            range.0,
            range.1,
            curve.resolution(CONFUSION).min(PCONFUSION),
        );
    }
    if continuity == 0 || (!curve.is_periodic() && continuity <= curve.continuity()) {
        return vec![range.0, range.1];
    }
    let raw = curve.parameter_intervals(continuity);
    let mut iv = vec![range.0];
    for &k in &raw {
        if k > range.0 + 1e-14 && k < range.1 - 1e-14 && iv.last().map_or(true, |p| (k - p).abs() > 1e-14)
        {
            iv.push(k);
        }
    }
    if iv.last().map_or(true, |p| (range.1 - p).abs() > 1e-14) {
        iv.push(range.1);
    }
    iv
}

/// `GeomAdaptor_Surface` Offset interval shape (`cxx:664-682`, `722-740`).
/// C0->C1, C1->C2, C2->C3; C3/CN keep CN. Mesh `initParameters` passes CN.
fn offset_interval_continuity(s: u8) -> u8 {
    match s {
        0 => 2,
        2 => 4,
        4 => 5,
        _ => 6,
    }
}

/// `GeomAdaptor_Surface::NbUIntervals` / `UIntervals` (`cxx:643-818`).
/// BSpline uses `VIso(FirstVKnot)` trimmed to the face U range; V uses `UIso`.
/// Offset forwards to the basis adaptor (`cxx:664-684`), never Offset UIso.
/// Revolution U is `[UFirst,ULast]`; V uses a BSpline generatrix (`cxx:712-720`).
fn adaptor_surface_intervals(
    surf: &dyn Surface,
    is_u: bool,
    continuity: u8,
    range: (f64, f64),
) -> Vec<f64> {
    match classify_surface(surf) {
        SurfaceType::OffsetSurface => {
            let base_s = offset_interval_continuity(continuity);
            if let Some(basis) = surf.offset_basis_surface() {
                return adaptor_surface_intervals(basis.as_ref(), is_u, base_s, range);
            }
        }
        SurfaceType::SurfaceOfRevolution => {
            if !is_u {
                if let Some(c) = surf.revolution_basis_curve() {
                    if c.bspline_knots().is_some() {
                        return adaptor_curve_intervals(c.as_ref(), continuity, range);
                    }
                }
            }
            return vec![range.0, range.1];
        }
        SurfaceType::BSplineSurface => {}
        _ => {
            return if is_u {
                surf.u_intervals(continuity)
            } else {
                surf.v_intervals(continuity)
            };
        }
    }
    let iso = if is_u {
        let v0 = surf
            .v_intervals(0)
            .first()
            .copied()
            .unwrap_or(surf.v_range().0);
        surf.v_iso_curve(v0)
    } else {
        let u0 = surf
            .u_intervals(0)
            .first()
            .copied()
            .unwrap_or(surf.u_range().0);
        surf.u_iso_curve(u0)
    };
    match iso {
        Some(c) => adaptor_curve_intervals(c.as_ref(), continuity, range),
        None => {
            if is_u {
                surf.u_intervals(continuity)
            } else {
                surf.v_intervals(continuity)
            }
        }
    }
}

/// `BRepMesh_NURBSRangeSplitter::getUndefinedInterval` (`cxx:414-449`).
///
/// The interval count and break array come from the untrimmed surface adaptor
/// (`cxx:420-421`, `cxx:438-447`): `GetSurface()` is the `BRepAdaptor_Surface`
/// built with `R = false` by `IMeshData_Face` (`IMeshData_Face.hxx`), so
/// `NbUIntervals` / `UIntervals` see the full `Geom_Surface` bounds. The
/// adjusted range `GetRangeU()` (`cxx:460-461`) is only fed to the uniform
/// fallback grid (`cxx:430`).
pub(super) fn get_undefined_interval(
    s: &dyn RangeSplitter,
    is_u: bool,
    continuity: u8,
    range: (f64, f64),
) -> Vec<f64> {
    let Some(surf) = s.surface() else {
        return vec![range.0, range.1];
    };
    let geom_range = if is_u { surf.u_range() } else { surf.v_range() };
    let iv = adaptor_surface_intervals(surf.as_ref(), is_u, continuity, geom_range);
    let mut intervals_nb = iv.len().saturating_sub(1) as i32;
    if intervals_nb == 1 {
        intervals_nb = s.get_undefined_interval_nb(is_u, continuity);
        if intervals_nb > 1 {
            let diff = (range.1 - range.0) / intervals_nb as f64;
            return (1..intervals_nb)
                .map(|i| range.0 + i as f64 * diff)
                .collect();
        }
    }
    if iv.is_empty() {
        vec![range.0, range.1]
    } else {
        iv
    }
}

/// `toSplitIntervals`: split when `GeomLib::NormEstim` fails at a corner.
pub(super) fn to_split_intervals(surface: &dyn Surface, u_iv: &[f64], v_iv: &[f64]) -> bool {
    use occt_core::precision::CONFUSION;
    use super::super::geomlib_norm::norm_estim;
    for &u in u_iv {
        if !u.is_finite() {
            continue;
        }
        for &v in v_iv {
            if !v.is_finite() {
                continue;
            }
            if norm_estim(surface, u, v, CONFUSION).0 != 0 {
                return true;
            }
        }
    }
    false
}

/// Filters a parameter set so no two samples are closer than `min_dist`, then
/// keeps samples separated by at least `filter_dist` (except a trailing buffer).
/// Source: `BRepMesh_NURBSRangeSplitter::filterParameters`.
pub(super) fn filter_parameters(params: &ParamSet, min_dist: f64, filter_dist: f64) -> Vec<f64> {
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
    // `cxx:614`: `for (j = 2; j < aParamLength; j++)` -- the last
    // pre-filtered sample is not series-filtered; it is appended after
    // the loop (`cxx:637`). Including it here promotes the trailing
    // candidate and duplicates the end parameter.
    let mut j = 1;
    while j + 1 < len {
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
pub(super) fn compute_grain_and_filter(
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
/// `BRepMesh_NURBSRangeSplitter::GenerateSurfaceNodes` (`cxx:317-398`).
pub(super) fn generate_nurbs_grid(s: &dyn RangeSplitter, params: &MeshParameters) -> Option<Vec<GpPnt2d>> {
    if !s.init_parameters() {
        return None;
    }
    let surface = s.surface()?;
    let range_u = s.range_u();
    let range_v = s.range_v();
    let delta = s.delta();
    let deflection = s.deflection();
    // `initParameters` always requests `GeomAbs_CN` (`cxx:441`).
    const GEOM_ABS_CN: u8 = 6;

    let mut u_params: ParamSet = s.parameters_u().cloned().unwrap_or_default();
    let mut v_params: ParamSet = s.parameters_v().cloned().unwrap_or_default();
    // `BRepMesh_NURBSRangeSplitter::initParameters` (`cxx:438-477`) fills CN
    // intervals. `BoundaryParamsRangeSplitter::initParameters` returns true
    // without that fill (`cxx` override) — only AddPoint samples.
    if s.seeds_cn_intervals() {
        // `initParameters` (`cxx:465-487`): fail the whole grid when either
        // direction adds no interval sample. INTERNAL edge params are already
        // in `parameters_u/v` via `grabParamsOfEdges` (`cxx:487`); do not
        // invent range-end samples (not in cxx).
        let u_intervals = get_undefined_interval(s, true, GEOM_ABS_CN, range_u);
        let v_intervals = get_undefined_interval(s, false, GEOM_ABS_CN, range_v);
        let is_split = to_split_intervals(surface.as_ref(), &u_intervals, &v_intervals);
        let mut u_iv = ParamSet::new();
        let mut v_iv = ParamSet::new();
        if !init_params_from_intervals(&u_intervals, range_u, is_split, &mut u_iv) {
            return None;
        }
        if !init_params_from_intervals(&v_intervals, range_v, is_split, &mut v_iv) {
            return None;
        }
        for &p in u_iv.iter() {
            u_params.insert(p);
        }
        for &p in v_iv.iter() {
            v_params.insert(p);
        }
    }

    let tol_u = param_resolution(surface.as_ref(), deflection, true);
    let tol_v = param_resolution(surface.as_ref(), deflection, false);
    let mut u_seq =
        compute_grain_and_filter(s, &u_params, tol_u, range_u.1 - range_u.0, delta.0, params);
    let mut v_seq =
        compute_grain_and_filter(s, &v_params, tol_v, range_v.1 - range_v.0, delta.1, params);

    let angle_interior = if params.angle_interior >= 0.0 { params.angle_interior } else { 2.0 * params.angle };
    let mut u_remove: std::collections::HashSet<u64> = std::collections::HashSet::new();
    let mut v_remove: std::collections::HashSet<u64> = std::collections::HashSet::new();
    let mut u_fixed: std::collections::HashSet<u64> = std::collections::HashSet::new();
    let mut v_fixed: std::collections::HashSet<u64> = std::collections::HashSet::new();
    analytical_filter(
        surface.as_ref(),
        false,
        &v_seq,
        &mut u_seq,
        deflection,
        angle_interior,
        params.min_size,
        &mut u_remove,
        &mut v_fixed,
        &mut u_fixed,
    );
    analytical_filter(
        surface.as_ref(),
        true,
        &u_seq,
        &mut v_seq,
        deflection,
        angle_interior,
        params.min_size,
        &mut v_remove,
        &mut u_fixed,
        &mut v_fixed,
    );
    for k in &u_fixed {
        u_remove.remove(k);
    }
    for k in &v_fixed {
        v_remove.remove(k);
    }

    let mut nodes = Vec::with_capacity(u_seq.len() * v_seq.len());
    for &u in &u_seq {
        if u_remove.contains(&u.to_bits()) {
            continue;
        }
        for &v in &v_seq {
            if v_remove.contains(&v.to_bits()) {
                continue;
            }
            nodes.push(GpPnt2d::new(u, v));
        }
    }
    Some(nodes)
}

/// Squared distance from `mid` to the infinite line through `p1..p2`.
/// `BRepMesh_GeomTool::SquareDeflectionOfSegment` (`hxx:181-193`) uses
/// `gp_Lin::SquareDistance`, not a clamped chord.
pub(super) fn sq_deflection_of_segment(p1: &GpPnt, p2: &GpPnt, mid: &GpPnt) -> f64 {
    crate::meshing::geom_tool::GeomTool::square_deflection_of_segment(p1, p2, mid)
}

/// `BRepMesh_NURBSRangeSplitter::AnalyticalFilter` (`cxx:30-215`): insert
/// midpoints that miss deflection / angle, and mark too-dense control params
/// for removal after both iso passes.
#[allow(clippy::too_many_arguments)]
pub(super) fn analytical_filter(
    surface: &dyn Surface,
    is_iso_u: bool,
    iso_params: &[f64],
    control_params: &mut Vec<f64>,
    deflection: f64,
    angle_interior: f64,
    min_size: f64,
    control_remove: &mut std::collections::HashSet<u64>,
    iso_forbidden: &mut std::collections::HashSet<u64>,
    control_forbidden: &mut std::collections::HashSet<u64>,
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

    // `cxx:74-75` builds `GeomAdaptor_Curve(UIso/VIso)` then `D1`.
    // Offset UIso/VIso is AdvApprox (`Geom_OffsetSurface.cxx:601-687`).
    for &iso_param in &iso_params[start..end] {
        let iso_curve = if is_iso_u {
            surface.u_iso_curve(iso_param)
        } else {
            surface.v_iso_curve(iso_param)
        };
        let iso_point = |t: f64| -> GpPnt {
            if let Some(ref c) = iso_curve {
                c.d0(t)
            } else if is_iso_u {
                surface.d0(iso_param, t)
            } else {
                surface.d0(t, iso_param)
            }
        };
        let iso_tangent = |t: f64| -> GpVec {
            if let Some(ref c) = iso_curve {
                c.d1(t).1
            } else {
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
            let angle = if prev_vec.square_magnitude() > SQUARE_CONFUSION
                && curr_vec.square_magnitude() > SQUARE_CONFUSION
            {
                prev_vec.angle(&curr_vec).abs()
            } else {
                0.0
            };

            if (sq_dist > sq_max_deflection || angle > angle_interior) && sq_dist > sq_min_size {
                control_params.insert(j, mid_param);
                continue;
            }

            if (sq_dist < sq_max_deflection || angle < angle_interior)
                && control_params.len() > 3
                && j + 1 < control_params.len()
            {
                let next_param = control_params[j + 1];
                let next_pnt = iso_point(next_param);
                let next_vec = iso_tangent(next_param);
                let look_mid = iso_point(0.5 * (prev_param + next_param));
                let look_dist = sq_deflection_of_segment(&prev_pnt, &next_pnt, &look_mid);
                if look_dist < sq_max_deflection {
                    // `BRepMesh_NURBSRangeSplitter.cxx:165-167`: both
                    // `SquareMagnitude()` guards use `gp::Resolution()` =
                    // `RealSmall()` = `DBL_MIN` (`gp.hxx:60`).
                    let look_ok = prev_vec.square_magnitude() < REAL_SMALL
                        || next_vec.square_magnitude() < REAL_SMALL
                        || prev_vec.angle(&next_vec).abs() < angle_interior;
                    if look_ok {
                        control_remove.insert(curr_param.to_bits());
                        prev_param = next_param;
                        prev_pnt = next_pnt;
                        prev_vec = next_vec;
                        j += 2;
                        continue;
                    } else {
                        iso_forbidden.insert(iso_param.to_bits());
                        control_forbidden.insert(curr_param.to_bits());
                    }
                }
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

        // `BRepMesh_DefaultRangeSplitter::GetSurface()` (`hxx:81`) is the
        // discrete face's own adaptor, which `IMeshData_Face` builds
        // UNRESTRICTED (`IMeshData_Face.hxx:69`
        // `new BRepAdaptor_Surface(GetFace(), false)`). Its
        // `First/Last U|VParameter` therefore come from `Geom_Surface::Bounds`
        // (`BRepAdaptor_Surface.cxx:79` `Load(aSurface, trsf)`), NOT from
        // `BRepTools::UVBounds`. `AdjustRange` (`cxx:45-81`) must use that plain
        // surface range; clipping the discrete range to the face's pcurve box is
        // the `Restriction=true` adaptor's behaviour and collapses faces whose
        // stored pcurves lie off the surface domain (e.g. a plane face whose
        // pcurves sit in a shifted parameterization), which leaves
        // `computeLengthV` == 0 and `IsValid` == false.
        let (gu0, gu1) = surf.u_range();
        let (gv0, gv1) = surf.v_range();
        update_range(gu0, gu1, surf.is_u_periodic(), &mut du_first, &mut du_second);
        if du_second < du_first {
            self.base_mut().set_valid(false);
            return;
        }

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
    /// Mutable U parameter map. Used by `grabParamsOfEdges`.
    fn parameters_u_mut(&mut self) -> Option<&mut ParamSet> {
        None
    }
    /// Mutable V parameter map. Used by `grabParamsOfEdges`.
    fn parameters_v_mut(&mut self) -> Option<&mut ParamSet> {
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
    /// `getUndefinedIntervalNb` — `(isU ? NbUPoles : NbVPoles) - 1` for NURBS.
    fn get_undefined_interval_nb(&self, _is_u: bool, continuity: u8) -> i32 {
        ((continuity as i32) + 1).clamp(2, 32)
    }

    /// Whether the surface parameter sets can be initialized. Source:
    /// `initParameters`.
    fn init_parameters(&self) -> bool {
        self.surface().is_some()
    }

    /// True when `initParameters` seeds U/V from CN intervals.
    /// `BRepMesh_BoundaryParamsRangeSplitter` returns false.
    fn seeds_cn_intervals(&self) -> bool {
        true
    }
}

/// Base range splitter. Source: `BRepMesh_DefaultRangeSplitter`.
pub struct DefaultRangeSplitter {
    pub(super) dface: Option<MeshFace>,
    pub(super) surface: Option<Arc<dyn Surface>>,
    pub(super) deflection: f64,
    pub(super) range_u: (f64, f64),
    pub(super) range_v: (f64, f64),
    pub(super) delta: (f64, f64),
    pub(super) tolerance: (f64, f64),
    pub(super) is_valid: bool,
}
