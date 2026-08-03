//! Edge bean / face intersection ranges — the numerical heart of the precise
//! boolean.
//!
//! Port of `IntTools_BeanFaceIntersector` (TKBO). Given a curve *bean* (the
//! part of an edge over a parameter window) and a surface, compute the
//! parameter ranges of the curve that lie **on** the surface within the
//! combined edge/face tolerance. The real face boundary is *not* taken into
//! account (as in OCCT): the results are ranges of the curve whose points are
//! within `tol_e + tol_f` of the face surface, inside the surface parameter
//! window `[umin,umax]×[vmin,vmax]`.
//!
//! The `Perform()` flow mirrors the OCCT dispatch:
//!
//! 1. **Line/plane** → `ComputeLinePlane` (substitute the line into the plane
//!    equation → one root, or the whole range when the line lies in the plane).
//! 2. **`FastComputeAnalytic`** — analytic coincidence / no-intersection
//!    verdicts for conic curves against quadrics (plane/circle/ellipse,
//!    cylinder/line/circle, sphere/line). Returns `true` when it can conclude
//!    decisively; otherwise computation continues.
//! 3. **Range manager** (`IntTools_MarkedRangeSet` semantics) + coincidence
//!    scan (`TestComputeCoinside`).
//! 4. **Localized path** for high-degree NURBS surfaces: box culling of
//!    subdivided curve cells, then exact curve–surface intersection per
//!    surviving cell.
//! 5. **Exact path**: `IntCurveSurface_HInter` (`crate::intcurvesurface`) for
//!    discrete points / coincidence segments, then `Extrema_ExtCS`
//!    (`occt_geom::extrema_surf`) for near-surface spans, then boundary
//!    completion (`ComputeNearRangeBoundaries`).
//!
//! `ponytail:` `dyn Surface` cannot be downcast, so the analytic fast paths
//! classify curves/surfaces by sampling geometric invariants (the established
//! pattern in `crate::intcurvesurface`), and the local point-to-surface solve
//! of OCCT's `Extrema_GenLocateExtPS` is replaced by the grid+refine projector
//! `crate::brep_surface::surface_closest_params`.

use std::sync::Arc;

use occt_core::bnd::BndBox;
use occt_core::gp::{GpDir, GpLin, GpPnt, GpVec};
use occt_core::precision::{ANGULAR, CONFUSION, PCONFUSION};
use occt_geom::extrema_surf::curve_surface_extrema_all;
use occt_geom::{Curve, Surface};

use crate::brep_surface::{classify_surface, surface_closest_params, SurfaceKind};
use crate::intcurvesurface::perform_curve_surface;
use crate::inttools_data::IntRange;
use crate::inttools_range::IntContext;
use crate::meshing::range_splitter::{classify_surface as classify_surface_mesh, SurfaceType};

const PI: f64 = std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Curve / surface classification helpers
// ---------------------------------------------------------------------------

/// Coarse analytic kind of a curve (sampling classification).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum FastCurveKind {
    Line,
    Circle,
    Ellipse,
    Other,
}

/// Whether the curve is geometrically a straight line (unbounded, or all
/// samples collinear with the first–last chord).
fn is_line_like(c: &dyn Curve) -> bool {
    let (a, b) = (c.first_parameter(), c.last_parameter());
    if !a.is_finite() || !b.is_finite() {
        return true;
    }
    if (b - a).abs() <= 1e-15 {
        return false;
    }
    let p0 = c.d0(a);
    let pl = c.d0(b);
    let size = p0.distance(&pl);
    if size <= 1e-30 {
        return false;
    }
    let d0 = GpVec::from_pnts(&p0, &pl);
    let tol = 1e-6 * size;
    for i in 1..8 {
        let p = c.d0(a + (b - a) * i as f64 / 8.0);
        if GpVec::from_pnts(&p0, &p).crossed(&d0).magnitude() > tol * size {
            return false;
        }
    }
    true
}

/// First three non-collinear samples of a curve.
fn first_three_spanning(pts: &[GpPnt]) -> Option<(GpPnt, GpPnt, GpPnt)> {
    let p0 = pts[0];
    let mut i1 = None;
    for (i, p) in pts.iter().enumerate().skip(1) {
        if GpVec::from_pnts(&p0, p).magnitude() > 1e-9 {
            i1 = Some(i);
            break;
        }
    }
    let i1 = i1?;
    let p1 = pts[i1];
    let d0 = GpVec::from_pnts(&p0, &p1);
    for p in pts.iter().skip(i1 + 1) {
        if GpVec::from_pnts(&p0, p).crossed(&d0).magnitude() > 1e-9 * d0.magnitude().max(1e-9) {
            return Some((p0, p1, *p));
        }
    }
    None
}

fn det3(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

fn solve3(a: &[[f64; 3]; 3], rhs: &[f64; 3]) -> Option<[f64; 3]> {
    let d = det3(a);
    if d.abs() < 1e-20 {
        return None;
    }
    let mut x = [0.0; 3];
    for k in 0..3 {
        let mut m = *a;
        for i in 0..3 {
            m[i][k] = rhs[i];
        }
        x[k] = det3(&m) / d;
    }
    Some(x)
}

/// Circumcenter of three non-collinear points.
fn circumcenter(a: &GpPnt, b: &GpPnt, c: &GpPnt) -> Option<GpPnt> {
    let d1 = GpVec::from_pnts(a, b);
    let d2 = GpVec::from_pnts(a, c);
    let n = d1.crossed(&d2);
    if n.magnitude() < 1e-30 {
        return None;
    }
    let n2 = |p: &GpPnt| p.coord.dot(&p.coord);
    let mat = [
        [d1.xyz().x, d1.xyz().y, d1.xyz().z],
        [d2.xyz().x, d2.xyz().y, d2.xyz().z],
        [n.xyz().x, n.xyz().y, n.xyz().z],
    ];
    let rhs = [0.5 * (n2(b) - n2(a)), 0.5 * (n2(c) - n2(a)), a.coord.dot(&n.xyz())];
    let o = solve3(&mat, &rhs)?;
    Some(GpPnt::new(o[0], o[1], o[2]))
}

fn midpoint(a: &GpPnt, b: &GpPnt) -> GpPnt {
    GpPnt::new(0.5 * (a.x() + b.x()), 0.5 * (a.y() + b.y()), 0.5 * (a.z() + b.z()))
}

/// Whether the curve is a full circle-like closed curve: all samples coplanar
/// and equidistant from a common centre.
fn is_circle_like(c: &dyn Curve) -> bool {
    let (a, b) = (c.first_parameter(), c.last_parameter());
    if !a.is_finite() || !b.is_finite() || (b - a).abs() <= 1e-15 {
        return false;
    }
    let n = 8;
    let pts: Vec<GpPnt> = (0..=n).map(|i| c.d0(a + (b - a) * i as f64 / n as f64)).collect();
    let (p0, p1, p2) = match first_three_spanning(&pts) {
        Some(x) => x,
        None => return false,
    };
    let center = match circumcenter(&p0, &p1, &p2) {
        Some(c) => c,
        None => return false,
    };
    let radius = p0.distance(&center);
    if radius <= 1e-30 {
        return false;
    }
    let nrm = GpVec::from_pnts(&p0, &p1).crossed(&GpVec::from_pnts(&p0, &p2));
    let m = nrm.magnitude();
    if m <= 1e-30 {
        return false;
    }
    let nv = nrm.divided(m);
    let scale = radius.max(1.0);
    let tol = 1e-6 * scale;
    for p in pts {
        let v = GpVec::from_pnts(&center, &p);
        if v.dot(&nv).abs() > tol {
            return false;
        }
        if (v.magnitude() - radius).abs() > tol {
            return false;
        }
    }
    true
}

/// `(center, radius, unit plane normal)` of a circle-like curve.
fn circle_geometry(c: &dyn Curve) -> Option<(GpPnt, f64, GpDir)> {
    let (a, b) = (c.first_parameter(), c.last_parameter());
    if !a.is_finite() || !b.is_finite() {
        return None;
    }
    let n = 8;
    let pts: Vec<GpPnt> = (0..=n).map(|i| c.d0(a + (b - a) * i as f64 / n as f64)).collect();
    let (p0, p1, p2) = first_three_spanning(&pts)?;
    let center = circumcenter(&p0, &p1, &p2)?;
    let radius = p0.distance(&center);
    if radius <= 1e-30 {
        return None;
    }
    let nrm = GpVec::from_pnts(&p0, &p1).crossed(&GpVec::from_pnts(&p0, &p2));
    let m = nrm.magnitude();
    if m <= 1e-30 {
        return None;
    }
    let d = GpDir::from_vec(&nrm.divided(m)).ok()?;
    Some((center, radius, d))
}

/// Whether the curve is an ellipse-like closed curve (coplanar, periodic 2π,
/// not a circle).
fn is_ellipse_like(c: &dyn Curve) -> bool {
    if is_circle_like(c) {
        return false;
    }
    if !c.is_periodic() {
        return false;
    }
    let period = c.period();
    if (period - 2.0 * PI).abs() > 1e-6 {
        return false;
    }
    let (a, b) = (c.first_parameter(), c.last_parameter());
    if !a.is_finite() || !b.is_finite() {
        return false;
    }
    // All samples coplanar.
    let p0 = c.d0(a);
    let p1 = c.d0(a + (b - a) / 6.0);
    let p2 = c.d0(a + 2.0 * (b - a) / 6.0);
    let (p0, p1, p2) = match first_three_spanning(&[p0, p1, p2]) {
        Some(x) => x,
        None => return false,
    };
    let nrm = GpVec::from_pnts(&p0, &p1).crossed(&GpVec::from_pnts(&p0, &p2));
    let m = nrm.magnitude();
    if m <= 1e-30 {
        return false;
    }
    let nv = nrm.divided(m);
    for i in 0..24 {
        let u = a + (b - a) * i as f64 / 24.0;
        let v = GpVec::from_pnts(&p0, &c.d0(u));
        if v.dot(&nv).abs() > 1e-5 * v.magnitude().max(1.0) {
            return false;
        }
    }
    true
}

/// `(center, axis direction)` of an ellipse-like curve (natural 2π
/// parametrisation: opposite samples share the centre, the plane normal is the
/// axis).
fn ellipse_geometry(c: &dyn Curve) -> Option<(GpPnt, GpDir)> {
    let center = midpoint(&c.d0(0.0), &c.d0(PI));
    let p0 = c.d0(0.0);
    let pq = c.d0(PI / 2.0);
    let nrm = GpVec::from_pnts(&center, &p0).crossed(&GpVec::from_pnts(&center, &pq));
    let m = nrm.magnitude();
    if m <= 1e-30 {
        return None;
    }
    let d = GpDir::from_vec(&nrm.divided(m)).ok()?;
    // Validate the centre with the perpendicular pair.
    let c2 = midpoint(&c.d0(PI / 2.0), &c.d0(3.0 * PI / 2.0));
    if c2.distance(&center) > 1e-4 * center.distance(&p0).max(1.0) {
        return None;
    }
    Some((center, d))
}

/// `(location, unit direction)` of a line-like curve.
fn line_geometry(c: &dyn Curve) -> Option<(GpPnt, GpDir)> {
    let p0 = c.d0(0.0);
    let p1 = c.d0(1.0);
    let v = GpVec::from_pnts(&p0, &p1);
    let m = v.magnitude();
    if m <= 1e-30 {
        return None;
    }
    let d = GpDir::from_vec(&v.divided(m)).ok()?;
    Some((p0, d))
}

fn classify_fast_curve(c: &dyn Curve) -> FastCurveKind {
    if is_line_like(c) {
        return FastCurveKind::Line;
    }
    if is_circle_like(c) {
        return FastCurveKind::Circle;
    }
    if is_ellipse_like(c) {
        return FastCurveKind::Ellipse;
    }
    FastCurveKind::Other
}

/// `(location, x_dir, y_dir, normal)` of a planar surface (natural frame).
fn plane_geometry(s: &dyn Surface) -> Option<(GpPnt, GpDir, GpDir, GpDir)> {
    let (u0, _) = s.u_range();
    let (v0, _) = s.v_range();
    let u0 = if u0.is_finite() { u0 } else { 0.0 };
    let v0 = if v0.is_finite() { v0 } else { 0.0 };
    let o = s.d0(u0, v0);
    let xv = GpVec::from_pnts(&o, &s.d0(u0 + 1.0, v0));
    let yv = GpVec::from_pnts(&o, &s.d0(u0, v0 + 1.0));
    let nv = xv.crossed(&yv);
    let n = GpDir::from_vec(&nv).ok()?;
    let x = GpDir::from_vec(&xv).ok()?;
    let y = GpDir::from_vec(&yv).ok()?;
    Some((o, x, y, n))
}

/// `(center, radius)` of a spherical surface.
fn sphere_geometry(s: &dyn Surface) -> Option<(GpPnt, f64)> {
    let center = crate::brep_surface::sphere_center(s)?;
    let (u0, _) = s.u_range();
    let u0 = if u0.is_finite() { u0 } else { 0.0 };
    let r = s.d0(u0, 0.0).distance(&center);
    if r < 1e-12 {
        return None;
    }
    Some((center, r))
}

/// `(axis, radius)` of a cylindrical surface (axis direction from the V
/// advance, axis point = circle centre at V start, radius = half the opposite
/// sample distance).
fn cylinder_geometry(s: &dyn Surface) -> Option<(occt_core::gp::GpAx1, f64)> {
    let u0 = 0.0;
    let c0 = midpoint(&s.d0(u0, 0.0), &s.d0(u0 + PI, 0.0));
    let c1 = midpoint(&s.d0(u0, 1.0), &s.d0(u0 + PI, 1.0));
    let z = GpVec::from_pnts(&c0, &c1);
    let zm = z.magnitude();
    if zm < 1e-12 {
        return None;
    }
    let zd = GpDir::from_vec(&z.divided(zm)).ok()?;
    let r = s.d0(u0, 0.0).distance(&c0);
    if r < 1e-12 {
        return None;
    }
    Some((occt_core::gp::GpAx1::new(c0, zd), r))
}

/// Signed distance from `p` to a plane (`plane_geometry` frame).
fn plane_distance(ploc: &GpPnt, nrm: &GpDir, p: &GpPnt) -> f64 {
    GpVec::from_pnts(ploc, p).dot(&GpVec::from_xyz(nrm.xyz())).abs()
}

// ---------------------------------------------------------------------------
// MarkedRangeSet (port of IntTools_MarkedRangeSet)
// ---------------------------------------------------------------------------

/// A sorted set of parameter ranges, each carrying an integer flag. Ranges are
/// stored as consecutive `[boundaries[i], boundaries[i+1])` intervals; inserting
/// a range splits existing intervals and marks the overlap.
#[derive(Debug, Clone)]
struct MarkedRangeSet {
    boundaries: Vec<f64>,
    flags: Vec<i32>,
}

impl MarkedRangeSet {
    fn new() -> Self {
        Self { boundaries: Vec::new(), flags: Vec::new() }
    }

    /// `[first, last]` covered by one range carrying `init_flag`.
    fn set_boundaries(&mut self, first: f64, last: f64, init_flag: i32) {
        self.boundaries = vec![first, last];
        self.flags = vec![init_flag];
    }

    fn len(&self) -> usize {
        self.flags.len()
    }

    fn is_empty(&self) -> bool {
        self.flags.is_empty()
    }

    fn range(&self, i: usize) -> IntRange {
        IntRange::new_unchecked(self.boundaries[i], self.boundaries[i + 1])
    }

    fn flag(&self, i: usize) -> i32 {
        self.flags[i]
    }

    fn set_flag(&mut self, i: usize, flag: i32) {
        self.flags[i] = flag;
    }

    /// Index (0-based) of the range containing `value`, using the OCCT
    /// `GetIndex(value, UseLower)` semantics; `-1` when `value` is outside the
    /// set (or exactly at the last boundary with `UseLower`).
    fn get_index(&self, value: f64, use_lower: bool) -> isize {
        if self.boundaries.is_empty() {
            return -1;
        }
        if (use_lower && value < self.boundaries[0])
            || (!use_lower && value <= self.boundaries[0])
        {
            return -1;
        }
        for i in 1..self.boundaries.len() {
            if (use_lower && value < self.boundaries[i])
                || (!use_lower && value <= self.boundaries[i])
            {
                return (i - 1) as isize;
            }
        }
        -1
    }

    /// Indices (0-based) of every range containing `value` (a boundary value
    /// belongs to both adjacent ranges). Port of `GetIndices`.
    fn get_indices(&self, value: f64) -> Vec<usize> {
        let mut out = Vec::new();
        if self.boundaries.is_empty() || value < self.boundaries[0] {
            return out;
        }
        let mut found = false;
        for i in 1..self.boundaries.len() {
            if found {
                if value >= self.boundaries[i - 1] {
                    out.push(i - 1);
                } else {
                    break;
                }
            } else if value <= self.boundaries[i] {
                out.push(i - 1);
                found = true;
            }
        }
        out
    }

    /// Insert `[first, last]` with `flag`, splitting covered ranges. Returns
    /// `false` when the boundaries do not resolve (OCCT `InsertRange`).
    fn insert_range(&mut self, first: f64, last: f64, flag: i32) -> bool {
        if self.boundaries.is_empty() {
            return false;
        }
        let mut idx1 = self.get_index(first, true);
        if idx1 < 0 {
            return false;
        }
        let mut idx2 = self.get_index(last, false);
        if idx2 < 0 {
            return false;
        }
        if idx2 < idx1 {
            std::mem::swap(&mut idx1, &mut idx2);
            if last < first {
                return false;
            }
        }
        let idx1 = idx1 as usize;
        let mut idx2 = idx2 as usize;
        let are_equal = idx1 == idx2;
        let prev_flag = self.flags[idx1];

        self.boundaries.insert(idx1 + 1, first);
        self.flags.insert(idx1 + 1, flag);
        idx2 += 1;
        self.boundaries.insert(idx2 + 1, last);
        if are_equal {
            self.flags.insert(idx2 + 1, prev_flag);
        } else {
            self.flags.insert(idx2, flag);
        }
        if !are_equal {
            for i in (idx1 + 1)..idx2 {
                self.flags[i] = flag;
            }
        }
        true
    }
}

// ---------------------------------------------------------------------------
// BeanFaceIntersector
// ---------------------------------------------------------------------------

/// Computes the parameter ranges of a curve bean that lie on a surface within
/// the combined edge/face tolerance. Port of `IntTools_BeanFaceIntersector`.
pub struct BeanFaceIntersector {
    curve: Option<Arc<dyn Curve>>,
    surface: Option<Arc<dyn Surface>>,
    first_parameter: f64,
    last_parameter: f64,
    umin: f64,
    umax: f64,
    vmin: f64,
    vmax: f64,
    bean_tolerance: f64,
    face_tolerance: f64,
    curve_resolution: f64,
    criteria: f64,
    range_manager: MarkedRangeSet,
    context: Option<IntContext>,
    results: Vec<IntRange>,
    is_done: bool,
    min_sq_distance: f64,
}

impl BeanFaceIntersector {
    /// Create an empty intersector; call [`initialize`](Self::initialize)
    /// before [`perform`](Self::perform).
    pub fn new() -> Self {
        Self {
            curve: None,
            surface: None,
            first_parameter: 0.0,
            last_parameter: 0.0,
            umin: 0.0,
            umax: 0.0,
            vmin: 0.0,
            vmax: 0.0,
            bean_tolerance: 0.0,
            face_tolerance: 0.0,
            curve_resolution: PCONFUSION,
            criteria: CONFUSION,
            range_manager: MarkedRangeSet::new(),
            context: None,
            results: Vec::new(),
            is_done: false,
            min_sq_distance: f64::MAX,
        }
    }

    /// Initialise with a curve bean and a face surface plus the edge/face
    /// tolerances. The surface parameter window defaults to the surface's own
    /// ranges (override with [`set_surface_parameters`](Self::set_surface_parameters)).
    pub fn initialize(&mut self, curve: Arc<dyn Curve>, surface: Arc<dyn Surface>, tol_e: f64, tol_f: f64) {
        self.curve = Some(curve);
        self.surface = Some(surface);
        self.bean_tolerance = tol_e;
        self.face_tolerance = tol_f;
        self.criteria = tol_e + tol_f;
        let a = self.curve().first_parameter();
        let b = self.curve().last_parameter();
        let (first, last) = if a.is_finite() && b.is_finite() && b >= a { (a, b) } else { (0.0, 0.0) };
        self.first_parameter = first;
        self.last_parameter = last;
        self.curve_resolution = self.resolution(self.criteria);
        let (u0, u1) = self.surface().u_range();
        let (v0, v1) = self.surface().v_range();
        self.set_surface_parameters(u0, u1, v0, v1);
        self.results.clear();
        self.is_done = false;
        self.min_sq_distance = f64::MAX;
    }

    /// Restrict the curve bean to `[first, last]`.
    pub fn set_bean_parameters(&mut self, first: f64, last: f64) {
        self.first_parameter = first;
        self.last_parameter = last;
    }

    /// Restrict the surface parameter window to `[umin,umax]×[vmin,vmax]`.
    pub fn set_surface_parameters(&mut self, umin: f64, umax: f64, vmin: f64, vmax: f64) {
        self.umin = umin;
        self.umax = umax;
        self.vmin = vmin;
        self.vmax = vmax;
    }

    /// Attach an intersection context (cached, for interface compatibility; the
    /// projector used here operates directly on the supplied surface).
    pub fn set_context(&mut self, ctx: IntContext) {
        self.context = Some(ctx);
    }

    /// Whether the last [`perform`](Self::perform) completed.
    pub fn is_done(&self) -> bool {
        self.is_done
    }

    /// The parameter ranges of the curve lying on the surface.
    pub fn result(&self) -> Vec<IntRange> {
        self.results.clone()
    }

    /// The minimal (squared) distance found between the curve and the surface.
    pub fn minimal_square_distance(&self) -> f64 {
        self.min_sq_distance
    }

    // ---- accessors -------------------------------------------------------

    fn curve(&self) -> &dyn Curve {
        self.curve.as_ref().expect("BeanFaceIntersector: curve not set").as_ref()
    }

    fn surface(&self) -> &dyn Surface {
        self.surface.as_ref().expect("BeanFaceIntersector: surface not set").as_ref()
    }

    fn curve_d0(&self, u: f64) -> GpPnt {
        self.curve().d0(u)
    }

    fn surface_d0(&self, u: f64, v: f64) -> GpPnt {
        self.surface().d0(u, v)
    }

    /// Approximate `BRepAdaptor_Curve::Resolution(tol)`: the parameter
    /// increment over which the curve deviates from a straight chord by ~`tol`.
    /// Uses the current bean window (unbounded curves must have it set via
    /// [`set_bean_parameters`](Self::set_bean_parameters)).
    fn resolution(&self, tol: f64) -> f64 {
        let c = self.curve();
        let a = self.first_parameter;
        let b = self.last_parameter;
        let span = if a.is_finite() && b.is_finite() && b > a { b - a } else { 0.0 };
        if span.abs() <= 1e-30 || tol <= 0.0 {
            return PCONFUSION;
        }
        let n = 8;
        let h = span / n as f64 * 0.5;
        let mut max_speed = 0.0f64;
        for i in 0..=n {
            let u = a + span * i as f64 / n as f64;
            let speed = c.d0(u - h).distance(&c.d0(u + h)) / (2.0 * h);
            if speed > max_speed {
                max_speed = speed;
            }
        }
        if max_speed <= 1e-30 {
            return tol.max(PCONFUSION);
        }
        (tol / max_speed).max(PCONFUSION)
    }

    // ---- Distance ----------------------------------------------------------

    /// Closest surface parameters `(u, v)` of `p` and the surface distance.
    ///
    /// Planes and spheres are projected **exactly** (the reconstructed analytic
    /// frame), because the grid+refine projector `surface_closest_params` is
    /// only ~1e-3 accurate on analytic quadrics — far too coarse for the
    /// tolerance-driven range walk (`criteria` ~ 1e-7). All other surfaces use
    /// the grid+refine projector.
    fn closest_params_dist(&self, p: &GpPnt) -> (f64, f64, f64) {
        match classify_surface(self.surface()) {
            SurfaceKind::Plane => {
                if let Some((ploc, px, py, pn)) = plane_geometry(self.surface()) {
                    let d = GpVec::from_pnts(&ploc, p);
                    let u = d.dot(&GpVec::from_xyz(px.xyz()));
                    let v = d.dot(&GpVec::from_xyz(py.xyz()));
                    let dist = d.dot(&GpVec::from_xyz(pn.xyz())).abs();
                    return (u, v, dist);
                }
            }
            SurfaceKind::Sphere => {
                if let Some((sc, sr)) = sphere_geometry(self.surface()) {
                    let d = GpVec::from_pnts(&sc, p);
                    let m = d.magnitude();
                    if m > 1e-30 {
                        let u = d.y().atan2(d.x()).rem_euclid(2.0 * PI);
                        let v = (d.z() / m).clamp(-1.0, 1.0).asin();
                        return (u, v, (m - sr).abs());
                    }
                }
            }
            _ => {}
        }
        let (u, v) = surface_closest_params(self.surface(), p, 24, 24);
        let dist = self.surface_d0(u, v).distance(p);
        (u, v, dist)
    }

    /// Distance from the curve point at `arg` to the surface; the closest
    /// surface parameters (clamped to the window) are written to `u`/`v`.
    /// Port of `Distance(arg, u, v)`.
    fn distance_with_uv(&self, arg: f64, u: &mut f64, v: &mut f64) -> f64 {
        let p = self.curve_d0(arg);
        *u = self.umin;
        *v = self.vmin;
        let (su, sv, d) = self.closest_params_dist(&p);
        *u = su.clamp(self.umin, self.umax);
        *v = sv.clamp(self.vmin, self.vmax);
        d
    }

    /// Distance from the curve point at `arg` to the surface.
    fn distance(&self, arg: f64) -> f64 {
        let mut u = 0.0;
        let mut v = 0.0;
        self.distance_with_uv(arg, &mut u, &mut v)
    }

    // ---- Range expansion ----------------------------------------------------

    /// Expand a result range from a known on-surface point, increasing or
    /// decreasing the parameter. Port of `ComputeRangeFromStartPoint(bool, …)`.
    fn compute_range_from_start_point(&mut self, to_increase: bool, parameter: f64, u: f64, v: f64) {
        let found = self.range_manager.get_index(parameter, to_increase);
        if found < 0 {
            return;
        }
        self.compute_range_from_start_point_idx(to_increase, parameter, u, v, found as usize);
    }

    /// Expand a result range from a known on-surface point, walking the
    /// parameter by `curve_resolution` and bisecting on the distance until the
    /// curve leaves the surface. Port of
    /// `ComputeRangeFromStartPoint(bool, double, double, double, int)`.
    fn compute_range_from_start_point_idx(
        &mut self,
        to_increase: bool,
        parameter: f64,
        u_param: f64,
        v_param: f64,
        mut valid_index: usize,
    ) {
        if self.range_manager.flag(valid_index) > 0 {
            return;
        }
        let mut min_delta = self.curve_resolution * 0.5;
        let mut delta_restrictor = 0.1 * (self.last_parameter - self.first_parameter);
        if min_delta > delta_restrictor {
            min_delta = delta_restrictor * 0.5;
        }
        let ten_of_min_delta = min_delta * 10.0;
        let mut delta = self.curve_resolution;
        let mut cur_par = if to_increase { parameter + delta } else { parameter - delta };
        let mut prev_par = parameter;
        let mut current_range = self.range_manager.range(valid_index);
        let mut boundary_condition = if to_increase {
            cur_par > current_range.last
        } else {
            cur_par < current_range.first
        };
        if boundary_condition {
            cur_par = if to_increase { current_range.last } else { current_range.first };
            boundary_condition = false;
        }
        let mut loop_counter = 0;
        let _ = (u_param, v_param); // initial guess consumed by the local projector
        let mut another_solution_found = false;
        let mut is_boundary_index = false;
        let mut is_valid_index = true;

        while delta >= min_delta && loop_counter <= 10 {
            let a_point = self.curve_d0(cur_par);
            let (su, sv, sd) = self.closest_params_dist(&a_point);
            let point_found = if sd < self.criteria {
                true
            } else {
                self.distance(cur_par) < self.criteria
            };
            let _ = (su, sv);

            if point_found {
                prev_par = cur_par;
                another_solution_found = true;
                if boundary_condition && (is_boundary_index || !is_valid_index) {
                    break;
                }
            } else {
                delta_restrictor = delta;
            }

            delta = if point_found { delta * 2.0 } else { delta * 0.5 };
            delta = if delta < delta_restrictor { delta } else { delta_restrictor };
            cur_par = if to_increase { prev_par + delta } else { prev_par - delta };

            if cur_par == prev_par {
                break;
            }

            boundary_condition = if to_increase {
                cur_par > current_range.last
            } else {
                cur_par < current_range.first
            };
            is_boundary_index = false;
            is_valid_index = true;

            if boundary_condition {
                is_boundary_index = (!to_increase && valid_index == 0)
                    || (to_increase && valid_index + 1 == self.range_manager.len());
                if !is_boundary_index {
                    if point_found {
                        let adj_flag = if to_increase {
                            self.range_manager.flag(valid_index + 1)
                        } else {
                            self.range_manager.flag(valid_index - 1)
                        };
                        if adj_flag == 0 {
                            valid_index = if to_increase { valid_index + 1 } else { valid_index - 1 };
                            current_range = self.range_manager.range(valid_index);
                            if (to_increase && cur_par > current_range.last)
                                || (!to_increase && cur_par < current_range.first)
                            {
                                cur_par = 0.5 * (current_range.first + current_range.last);
                                delta *= 0.5;
                            }
                        } else {
                            is_valid_index = false;
                            cur_par = if to_increase { current_range.last } else { current_range.first };
                        }
                    }
                } else {
                    cur_par = if to_increase { current_range.last } else { current_range.first };
                }
                if delta < ten_of_min_delta {
                    loop_counter += 1;
                } else {
                    loop_counter = 0;
                }
            }
        }

        if another_solution_found {
            if to_increase {
                self.range_manager.insert_range(parameter, prev_par, 2);
            } else {
                self.range_manager.insert_range(prev_par, parameter, 2);
            }
        }
    }

    /// Insert a degenerate (point) result range when a crossing point added no
    /// interval — a tangency point. Port of the static `SetEmptyResultRange`.
    fn set_empty_result_range(&mut self, parameter: f64) {
        let indices = self.range_manager.get_indices(parameter);
        let mut add = !indices.is_empty();
        for &k in &indices {
            if self.range_manager.flag(k) == 2 {
                add = false;
                break;
            }
        }
        if add {
            self.range_manager.insert_range(parameter, parameter, 2);
        }
    }

    // ---- Main algorithm ------------------------------------------------------

    /// Launch the intersection. Sets `is_done` and fills `result()`.
    pub fn perform(&mut self) -> Result<(), String> {
        if self.curve.is_none() || self.surface.is_none() {
            return Err("BeanFaceIntersector::perform: curve and surface must be initialized".into());
        }
        self.is_done = false;
        self.results.clear();
        self.min_sq_distance = f64::MAX;

        // Fast computation of the Line/Plane case.
        if classify_fast_curve(self.curve()) == FastCurveKind::Line
            && classify_surface(self.surface()) == SurfaceKind::Plane
        {
            self.compute_line_plane();
            return Ok(());
        }

        // Fast analytic coincidence / no-intersection verdict.
        if self.fast_compute_analytic() {
            self.is_done = true;
            return Ok(());
        }

        // Range manager over the bean window.
        self.range_manager.set_boundaries(self.first_parameter, self.last_parameter, 0);

        // Coincidence scan over the whole bean.
        if self.test_compute_coinside() {
            self.results.push(IntRange::new_unchecked(self.first_parameter, self.last_parameter));
            self.is_done = true;
            return Ok(());
        }

        // Localized path for high-degree NURBS-like surfaces.
        let b_localize = self.surface_is_localizable();
        let is_localized = b_localize && self.compute_localized();

        if !is_localized {
            self.compute_around_exact_intersection();
            self.compute_using_extremum();
            self.compute_near_range_boundaries();
        }

        self.is_done = true;

        // Collect the flagged (== 2) ranges, merging adjacent ones.
        for i in 0..self.range_manager.len() {
            if self.range_manager.flag(i) != 2 {
                continue;
            }
            let r = self.range_manager.range(i);
            if let Some(last) = self.results.last_mut() {
                if (r.first - last.last).abs() > PCONFUSION {
                    self.results.push(r);
                } else {
                    last.last = last.last.max(r.last);
                }
            } else {
                self.results.push(r);
            }
        }
        Ok(())
    }

    /// Whether the localized path should be attempted for this surface type.
    fn surface_is_localizable(&self) -> bool {
        let finite = self.umin.is_finite()
            && self.umax.is_finite()
            && self.vmin.is_finite()
            && self.vmax.is_finite();
        if !finite {
            return false;
        }
        match classify_surface_mesh(self.surface()) {
            SurfaceType::BezierSurface | SurfaceType::OtherSurface => true,
            SurfaceType::BSplineSurface => true,
            _ => false,
        }
    }

    /// Line × plane intersection: substitute the line into the plane equation,
    /// then emit either a single root range (expanded by the tolerance-derived
    /// parameter width) or the whole range when the line lies in the plane.
    /// Port of `ComputeLinePlane`.
    fn compute_line_plane(&mut self) {
        let tol_ang = 1e-9;
        self.is_done = true;

        let (ploc, _px, _py, pn) = match plane_geometry(self.surface()) {
            Some(g) => g,
            None => return,
        };
        let (orig, ld) = match line_geometry(self.curve()) {
            Some(g) => g,
            None => return,
        };
        let nrm = GpVec::from_xyz(pn.xyz());
        let (a, b, c) = (nrm.x(), nrm.y(), nrm.z());
        let dcoef = -nrm.dot(&GpVec::from_pnts(&GpPnt::zero(), &ploc));
        let (al, bl, cl) = (ld.x(), ld.y(), ld.z());
        let direc = a * al + b * bl + c * cl;
        let dis = a * orig.x() + b * orig.y() + c * orig.z() + dcoef;

        let (mut parallel, mut inplane) = (false, false);
        if direc.abs() < tol_ang {
            parallel = true;
            inplane = dis.abs() < self.criteria;
        } else {
            let p1 = self.curve_d0(self.first_parameter);
            let p2 = self.curve_d0(self.last_parameter);
            let mut d1 = a * p1.x() + b * p1.y() + c * p1.z() + dcoef;
            if d1 < 0.0 {
                d1 = -d1;
            }
            let mut d2 = a * p2.x() + b * p2.y() + c * p2.z() + dcoef;
            if d2 < 0.0 {
                d2 = -d2;
            }
            if d1 <= self.criteria && d2 <= self.criteria {
                inplane = true;
            }
        }

        if inplane {
            self.results.push(IntRange::new_unchecked(self.first_parameter, self.last_parameter));
            return;
        }
        if parallel {
            return;
        }

        let t = -dis / direc;
        if t < self.first_parameter || t > self.last_parameter {
            return;
        }
        let pint = orig.translated_vec(&GpVec::from_xyz(ld.xyz()).multiplied_scalar(t));
        let (u, v) = plane_uv_of_point(&ploc, &ld, &pint, &_px, &_py);
        if self.umin > u || u > self.umax || self.vmin > v || v > self.vmax {
            return;
        }

        // Parameter half-width from the tolerances and the incidence angle.
        let angle = (PI * 0.5 - ld.angle(&pn)).abs();
        let a_dt = compute_int_range(self.bean_tolerance, self.face_tolerance, angle);
        let t1 = self.first_parameter.max(t - a_dt);
        let t2 = self.last_parameter.min(t + a_dt);
        self.results.push(IntRange::new_unchecked(t1, t2));
    }

    /// Fast analytic coincidence / no-intersection checks for conic curves
    /// against quadric surfaces. Returns `true` when a decisive verdict was
    /// reached; otherwise computation continues. Port of `FastComputeAnalytic`.
    fn fast_compute_analytic(&mut self) -> bool {
        let ck = classify_fast_curve(self.curve());
        if ck == FastCurveKind::Other {
            return false;
        }
        let sk = classify_surface(self.surface());
        let mut is_coincide = false;
        let mut has_intersection = true;

        match sk {
            SurfaceKind::Plane => {
                let (ploc, _px, _py, pn) = match plane_geometry(self.surface()) {
                    Some(g) => g,
                    None => return false,
                };
                let (adir, aloc) = match ck {
                    FastCurveKind::Circle => {
                        let (c, _r, n) = match circle_geometry(self.curve()) {
                            Some(g) => g,
                            None => return false,
                        };
                        (n, c)
                    }
                    FastCurveKind::Ellipse => {
                        let (c, n) = match ellipse_geometry(self.curve()) {
                            Some(g) => g,
                            None => return false,
                        };
                        (n, c)
                    }
                    _ => return false,
                };
                let angle = adir.angle(&pn);
                if angle > ANGULAR {
                    return false;
                }
                has_intersection = false;
                let dist = plane_distance(&ploc, &pn, &aloc);
                is_coincide = dist < self.criteria;
            }
            SurfaceKind::Sphere => {
                let (sc, sr) = match sphere_geometry(self.surface()) {
                    Some(g) => g,
                    None => return false,
                };
                if ck == FastCurveKind::Line {
                    let (lloc, ldir) = match line_geometry(self.curve()) {
                        Some(g) => g,
                        None => return false,
                    };
                    let lin = GpLin::from_pnt_dir(lloc, ldir);
                    let dist = lin.distance(&sc) - sr;
                    has_intersection = dist < self.criteria;
                } else {
                    return false;
                }
            }
            SurfaceKind::Cylinder | SurfaceKind::Other => {
                // `brep_surface::classify_surface` never returns `Cylinder`
                // directly; re-confirm via the mesh classifier.
                if classify_surface_mesh(self.surface()) != SurfaceType::Cylinder {
                    return false;
                }
                let (axis, radius) = match cylinder_geometry(self.surface()) {
                    Some(g) => g,
                    None => return false,
                };
                match ck {
                    FastCurveKind::Line => {
                        let (lloc, ldir) = match line_geometry(self.curve()) {
                            Some(g) => g,
                            None => return false,
                        };
                        if !ldir.is_parallel(axis.direction()) {
                            return false;
                        }
                        has_intersection = false;
                        let lin = GpLin::from_pnt_dir(lloc, ldir);
                        let dist = (lin.distance(axis.location()) - radius).abs();
                        is_coincide = dist < self.criteria;
                    }
                    FastCurveKind::Circle => {
                        let (cloc, cr, cn) = match circle_geometry(self.curve()) {
                            Some(g) => g,
                            None => return false,
                        };
                        let angle = axis.direction().angle(&cn);
                        if angle > ANGULAR {
                            return false;
                        }
                        let axis_lin = GpLin::from_pnt_dir(*axis.location(), *axis.direction());
                        let dist_loc = axis_lin.distance(&cloc);
                        let dist = dist_loc + (cr - radius).abs();
                        is_coincide = dist < self.criteria;
                        if !is_coincide {
                            has_intersection = (dist_loc - (cr + radius)) < self.criteria
                                && ((cr - radius).abs() - dist_loc) < self.criteria;
                        }
                    }
                    _ => return false,
                }
            }
            SurfaceKind::Cone | SurfaceKind::Torus => return false,
        }

        if is_coincide {
            self.results.push(IntRange::new_unchecked(self.first_parameter, self.last_parameter));
        }
        is_coincide || !has_intersection
    }

    /// Scan the whole bean for coincidence with the surface: sample 23 points,
    /// expand ranges from each. Port of `TestComputeCoinside`.
    fn test_compute_coinside(&mut self) -> bool {
        let cfp = self.first_parameter;
        let clp = self.last_parameter;
        let nb_seg = 23;
        let cdp = (clp - cfp) / nb_seg as f64;

        let mut u = 0.0;
        let mut v = 0.0;
        if self.distance_with_uv(cfp, &mut u, &mut v) > self.criteria {
            return false;
        }
        self.compute_range_from_start_point(true, cfp, u, v);

        let found = self.range_manager.get_index(clp, false);
        if found >= 0 && self.range_manager.flag(found as usize) == 2 {
            return true;
        }
        if self.distance_with_uv(clp, &mut u, &mut v) > self.criteria {
            return false;
        }
        self.compute_range_from_start_point(false, clp, u, v);

        for i in 1..nb_seg {
            let par = cfp + i as f64 * cdp;
            if self.distance_with_uv(par, &mut u, &mut v) > self.criteria {
                return false;
            }
            let n = self.range_manager.len();
            self.compute_range_from_start_point(false, par, u, v);
            self.compute_range_from_start_point(true, par, u, v);
            if n == self.range_manager.len() {
                self.set_empty_result_range(par);
            }
        }
        true
    }

    /// Exact curve–surface intersection via `IntCurveSurface_HInter`, then
    /// range expansion around every intersection point / segment. Port of
    /// `ComputeAroundExactIntersection`.
    fn compute_around_exact_intersection(&mut self) {
        let uv = self.finite_uv_bounds();
        let res = perform_curve_surface(
            self.curve(),
            self.surface(),
            (self.first_parameter, self.last_parameter),
            uv,
        );
        let hr = match res {
            Ok(h) => h,
            Err(_) => return,
        };

        // With more than one point, tighten the criteria to avoid merging
        // distinct crossings into one range.
        if hr.nb_points() > 1 {
            self.criteria = 3.0 * CONFUSION;
            self.curve_resolution = self.resolution(self.criteria);
        }

        for i in 0..hr.nb_points() {
            let p = hr.point(i);
            let mut w = p.param();
            if w < self.first_parameter || w > self.last_parameter {
                continue;
            }
            // Refine to the exact crossing (HInter's general path is only
            // accurate to ~1e-6, coarser than the result criteria).
            w = self.refine_crossing_param(w);
            let (su, sv) = {
                let pnt = self.curve_d0(w);
                let (a, b, _) = self.closest_params_dist(&pnt);
                (a, b)
            };
            let mut u = su;
            let mut v = sv;
            let u_not_valid = self.umin > u || u > self.umax;
            let v_not_valid = self.vmin > v || v > self.vmax;
            let mut solution_is_valid = !u_not_valid && !v_not_valid;

            if u_not_valid || v_not_valid {
                let mut b_u_corrected = true;
                if u_not_valid {
                    b_u_corrected = false;
                    solution_is_valid = false;
                    if self.surface().is_u_periodic() {
                        u = adjust_periodic(u, self.umin, self.umax, 2.0 * PI);
                        solution_is_valid = true;
                        b_u_corrected = true;
                    }
                }
                if b_u_corrected && v_not_valid {
                    solution_is_valid = false;
                    if self.surface().is_v_periodic() {
                        v = adjust_periodic(v, self.vmin, self.vmax, 2.0 * PI);
                        solution_is_valid = true;
                    }
                }
            }

            if !solution_is_valid {
                continue;
            }

            let n = self.range_manager.len();
            self.compute_range_from_start_point(false, w, u, v);
            self.compute_range_from_start_point(true, w, u, v);
            if n == self.range_manager.len() {
                self.set_empty_result_range(w);
            } else {
                self.min_sq_distance = 0.0;
            }
        }

        for i in 0..hr.nb_segments() {
            let seg = hr.segment(i);
            let p1 = seg.first_point();
            let p2 = seg.second_point();
            let first_param = if p1.param() < self.first_parameter {
                self.first_parameter
            } else {
                p1.param()
            };
            let last_param = if self.last_parameter < p2.param() {
                self.last_parameter
            } else {
                p2.param()
            };
            self.range_manager.insert_range(first_param, last_param, 2);
            self.compute_range_from_start_point(false, p1.param(), p1.u(), p1.v());
            self.compute_range_from_start_point(true, p2.param(), p2.u(), p2.v());
            self.min_sq_distance = 0.0;
        }
    }

    /// Complete result ranges whose start/end boundaries are near the surface
    /// but were missed by the discrete intersection points. Port of
    /// `ComputeNearRangeBoundaries`.
    fn compute_near_range_boundaries(&mut self) {
        let mut u = self.umin;
        let mut v = self.vmin;

        let n = self.range_manager.len();
        for i in 0..n {
            if self.range_manager.flag(i) > 0 {
                continue;
            }
            if i > 0 && self.range_manager.flag(i - 1) > 0 {
                continue;
            }
            let r = self.range_manager.range(i);
            if self.distance_with_uv(r.first, &mut u, &mut v) < self.criteria {
                let old_len = self.range_manager.len();
                if i > 0 {
                    self.compute_range_from_start_point_idx(false, r.first, u, v, i - 1);
                }
                let idx = i + (self.range_manager.len() - old_len);
                if idx < self.range_manager.len() {
                    self.compute_range_from_start_point_idx(true, r.first, u, v, idx);
                }
                if old_len == self.range_manager.len() {
                    self.set_empty_result_range(r.first);
                }
            }
        }

        if self.range_manager.is_empty() {
            return;
        }
        let last_idx = self.range_manager.len() - 1;
        if self.range_manager.flag(last_idx) == 0 {
            let r = self.range_manager.range(last_idx);
            if self.distance_with_uv(r.last, &mut u, &mut v) < self.criteria {
                let old_len = self.range_manager.len();
                self.compute_range_from_start_point_idx(false, r.last, u, v, last_idx);
                if old_len == self.range_manager.len() {
                    self.set_empty_result_range(r.last);
                }
            }
        }
    }

    /// Refine an HInter crossing parameter to the local minimum of the exact
    /// surface distance. The HInter general path converges with the coarse
    /// projector (`~1e-6` distance error), which can exceed the tightened
    /// result criteria (`3·Confusion`); refining makes the walk's starting
    /// point lie genuinely on the surface.
    fn refine_crossing_param(&self, w: f64) -> f64 {
        let window = (10.0 * self.curve_resolution).max(1e-5);
        let f = |u: f64| self.closest_params_dist(&self.curve_d0(u)).2;
        golden_1d(&f, w - window, w + window, 1e-10).0
    }

    /// Use curve–surface extrema to find near-surface spans that HInter's
    /// discrete points missed (tangencies / parallel spans). Port of
    /// `ComputeUsingExtremum`.
    fn compute_using_extremum(&mut self) {
        let tol = PCONFUSION;
        let mut i = 0usize;
        while i < self.range_manager.len() {
            if self.range_manager.flag(i) > 0 {
                i += 1;
                continue;
            }
            let r = self.range_manager.range(i);
            let anarg1 = r.first;
            let anarg2 = r.last;

            if anarg2 - anarg1 < PCONFUSION {
                if (i > 0 && self.range_manager.flag(i - 1) == 2)
                    || (i + 1 < self.range_manager.len() && self.range_manager.flag(i + 1) == 2)
                {
                    self.range_manager.set_flag(i, 1);
                    i += 1;
                    continue;
                }
            }

            let old_len = self.range_manager.len();
            let mut solution_found = false;

            // All curve–surface extrema, restricted to this range + surface window.
            let exts = curve_surface_extrema_all(self.curve(), self.surface());
            let mut candidates: Vec<(f64, f64, f64)> = Vec::new();
            for e in &exts {
                self.min_sq_distance = self.min_sq_distance.min(e.distance * e.distance);
                if e.distance * e.distance >= self.criteria * self.criteria {
                    continue;
                }
                if e.u1 < anarg1 - tol || e.u1 > anarg2 + tol {
                    continue;
                }
                let u = e.u2;
                let v = e.v2.unwrap_or(0.0);
                if u < self.umin || u > self.umax || v < self.vmin || v > self.vmax {
                    continue;
                }
                candidates.push((e.u1, u, v));
            }

            // Fall back to the single global extremum when the list is empty.
            if candidates.is_empty() {
                let e = occt_geom::extrema::curve_surface_extrema(self.curve(), self.surface(), 16);
                self.min_sq_distance = self.min_sq_distance.min(e.distance * e.distance);
                if e.distance * e.distance < self.criteria * self.criteria
                    && e.u1 >= anarg1 - tol
                    && e.u1 <= anarg2 + tol
                {
                    let u = e.u2;
                    let v = e.v2.unwrap_or(0.0);
                    if u >= self.umin && u <= self.umax && v >= self.vmin && v <= self.vmax {
                        let n = self.range_manager.len();
                        self.compute_range_from_start_point(false, e.u1, u, v);
                        self.compute_range_from_start_point(true, e.u1, u, v);
                        solution_found = true;
                        if n == self.range_manager.len() {
                            self.set_empty_result_range(e.u1);
                        }
                    }
                }
            } else {
                for (t, u, v) in candidates {
                    let n = self.range_manager.len();
                    self.compute_range_from_start_point(false, t, u, v);
                    self.compute_range_from_start_point(true, t, u, v);
                    solution_found = true;
                    if n == self.range_manager.len() {
                        self.set_empty_result_range(t);
                    }
                }
            }

            if !solution_found {
                self.range_manager.set_flag(i, 1);
            }
            let diff = self.range_manager.len() - old_len;
            if diff > 0 {
                i += diff;
            } else {
                i += 1;
            }
        }
    }

    /// Localized intersection for high-degree NURBS surfaces: box-cull the
    /// subdivided curve cells against the surface box, then run the exact
    /// curve–surface intersection per surviving cell. Port of `ComputeLocalized`
    /// + `LocalizeSolutions` (the subdivision machinery is collapsed into the
    /// per-cell box test + exact solve).
    fn compute_localized(&mut self) -> bool {
        let fbox = self.surface_box(self.umin, self.umax, self.vmin, self.vmax, self.criteria);
        let ebox = self.curve_box(self.first_parameter, self.last_parameter, self.bean_tolerance);

        // Whole-domain rejection: nothing to do.
        if ebox.is_out_box(&fbox) {
            for i in 0..self.range_manager.len() {
                self.range_manager.set_flag(i, 1);
            }
            return true;
        }

        // Subdivide the curve bean into cells; flag cells whose curve box is
        // disjoint from the surface box.
        let nb_c = 3usize;
        let first = self.first_parameter;
        let last = self.last_parameter;
        let span = last - first;
        let mut any_survive = false;
        let mut cells: Vec<(f64, f64)> = Vec::with_capacity(nb_c);
        for k in 0..nb_c {
            let a = first + k as f64 * span / nb_c as f64;
            let b = if k + 1 == nb_c { last } else { first + (k + 1) as f64 * span / nb_c as f64 };
            let cbox = self.curve_box(a, b, self.bean_tolerance);
            if cbox.is_out_box(&fbox) {
                self.range_manager.insert_range(a, b, 1);
            } else {
                cells.push((a, b));
                any_survive = true;
            }
        }
        if !any_survive {
            return true;
        }

        // Exact intersection per surviving cell.
        let uv = self.finite_uv_bounds();
        for (a, b) in cells {
            if let Ok(hr) = perform_curve_surface(self.curve(), self.surface(), (a, b), uv) {
                for i in 0..hr.nb_points() {
                    let p = hr.point(i);
                    let w = p.param();
                    if w < self.first_parameter || w > self.last_parameter {
                        continue;
                    }
                    let u = p.u().clamp(self.umin, self.umax);
                    let v = p.v().clamp(self.vmin, self.vmax);
                    let n = self.range_manager.len();
                    self.compute_range_from_start_point(false, w, u, v);
                    self.compute_range_from_start_point(true, w, u, v);
                    if n == self.range_manager.len() {
                        self.set_empty_result_range(w);
                    } else {
                        self.min_sq_distance = 0.0;
                    }
                }
            }
        }

        self.compute_near_range_boundaries();
        true
    }

    /// Finite `(u0, v0, u1, v1)` bounds for the HInter call; unbounded
    /// directions are clamped to a broad window (they only arise for analytic
    /// surfaces handled elsewhere).
    fn finite_uv_bounds(&self) -> (f64, f64, f64, f64) {
        let clamp = |a: f64, b: f64| {
            if a.is_finite() && b.is_finite() {
                (a, b)
            } else {
                (-1e4, 1e4)
            }
        };
        let (u0, u1) = clamp(self.umin, self.umax);
        let (v0, v1) = clamp(self.vmin, self.vmax);
        (u0, u1, v0, v1)
    }

    /// Axis-aligned bounding box of the surface patch over `[u0,u1]×[v0,v1]`,
    /// sampled on a grid and enlarged by `tol`.
    fn surface_box(&self, u0: f64, u1: f64, v0: f64, v1: f64, tol: f64) -> BndBox {
        let mut b = BndBox::new();
        let nu = 12usize;
        let nv = 12usize;
        for i in 0..=nu {
            for j in 0..=nv {
                let u = u0 + (u1 - u0) * i as f64 / nu as f64;
                let v = v0 + (v1 - v0) * j as f64 / nv as f64;
                b.add_point(&self.surface_d0(u, v));
            }
        }
        b.enlarge(tol);
        b
    }

    /// Axis-aligned bounding box of the curve over `[a,b]`, sampled and
    /// enlarged by `tol`.
    fn curve_box(&self, a: f64, b: f64, tol: f64) -> BndBox {
        let mut bb = BndBox::new();
        let n = 32usize;
        for i in 0..=n {
            let u = a + (b - a) * i as f64 / n as f64;
            bb.add_point(&self.curve_d0(u));
        }
        bb.enlarge(tol);
        bb
    }
}

impl Default for BeanFaceIntersector {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Static helpers
// ---------------------------------------------------------------------------

/// Wrap `value` into `[min, max]` by whole periods (positive modulo).
fn adjust_periodic(value: f64, min: f64, max: f64, period: f64) -> f64 {
    let p = period.abs();
    if p <= 1e-30 || !min.is_finite() || !max.is_finite() {
        return value;
    }
    let mut v = value;
    while v < min {
        v += p;
    }
    while v > max {
        v -= p;
    }
    v
}

/// `(u, v)` parameters of `p` in a plane's natural frame.
fn plane_uv_of_point(ploc: &GpPnt, _pn: &GpDir, p: &GpPnt, px: &GpDir, py: &GpDir) -> (f64, f64) {
    let d = GpVec::from_pnts(ploc, p);
    (d.dot(&GpVec::from_xyz(px.xyz())), d.dot(&GpVec::from_xyz(py.xyz())))
}

/// Port of `IntTools_Tools::ComputeIntRange`: the parameter half-width that
/// covers the tolerance band around a crossing at incidence `angle`.
fn compute_int_range(tol1: f64, tol2: f64, angle: f64) -> f64 {
    if (PI * 0.5 - angle).abs() < ANGULAR {
        return tol2;
    }
    let an_angle = if angle > PI * 0.5 { PI - angle } else { angle };
    let a1 = tol1 * (PI * 0.5 - an_angle).tan();
    let a2 = tol2 / an_angle.sin();
    a1 + a2
}

/// Golden-section minimization of `f` over `[lo, hi]`. Returns `(argmin, min)`.
fn golden_1d<F: Fn(f64) -> f64>(f: &F, lo: f64, hi: f64, eps: f64) -> (f64, f64) {
    const GOLD: f64 = 0.618_033_988_749_894_9;
    let mut a = lo;
    let mut b = hi;
    let mut c = b - GOLD * (b - a);
    let mut d = a + GOLD * (b - a);
    let mut fc = f(c);
    let mut fd = f(d);
    while (b - a) > eps {
        if fc < fd {
            b = d;
            d = c;
            fd = fc;
            c = b - GOLD * (b - a);
            fc = f(c);
        } else {
            a = c;
            c = d;
            fc = fd;
            d = a + GOLD * (b - a);
            fd = f(d);
        }
    }
    let x = 0.5 * (a + b);
    (x, f(x))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;

    use occt_core::gp::{GpAx2, GpAx3, GpCirc, GpDir, GpLin, GpPln, GpSphere as GpSphereT};
    use occt_geom::{GeomBSplineCurve, GeomLine, GeomPlane, GeomSphere};

    use crate::brep_extrema::test_box::unit_box;
    use crate::brep_tool::BRepTool;

    fn dir(x: f64, y: f64, z: f64) -> GpDir {
        GpDir::new(x, y, z).expect("dir")
    }

    fn sphere_at(origin: GpPnt, r: f64) -> GeomSphere {
        let ax3 = GpAx3::new(origin, dir(0.0, 0.0, 1.0), &dir(1.0, 0.0, 0.0)).unwrap();
        GeomSphere::new(GpSphereT::new(ax3, r).unwrap())
    }

    fn unit_sphere() -> Arc<dyn Surface> {
        Arc::new(sphere_at(GpPnt::zero(), 1.0))
    }

    /// A line edge through the box bottom face along +Z at (0.5, 0.5).
    fn z_line(px: f64, py: f64, z0: f64, z1: f64) -> Arc<dyn Curve> {
        let dir = if z1 >= z0 { dir(0.0, 0.0, 1.0) } else { dir(0.0, 0.0, -1.0) };
        Arc::new(GeomLine::new(GpLin::from_pnt_dir(GpPnt::new(px, py, z0), dir)))
    }

    /// The bottom face (z = 0) surface of the unit box.
    fn box_bottom_surface() -> Arc<dyn Surface> {
        let b = unit_box();
        let face = &b.faces[0];
        BRepTool::face_surface(face).expect("box bottom face surface")
    }

    fn run(curve: Arc<dyn Curve>, surface: Arc<dyn Surface>, u: (f64, f64, f64, f64)) -> Vec<IntRange> {
        let mut bfi = BeanFaceIntersector::new();
        let (cf, cl) = (curve.first_parameter(), curve.last_parameter());
        let (cf, cl) = if cf.is_finite() && cl.is_finite() { (cf, cl) } else { (-10.0, 10.0) };
        bfi.initialize(curve, surface, 1e-7, 1e-7);
        bfi.set_bean_parameters(cf, cl);
        bfi.set_surface_parameters(u.0, u.1, u.2, u.3);
        bfi.perform().expect("perform succeeds");
        assert!(bfi.is_done());
        bfi.result()
    }

    #[test]
    fn line_through_box_face_single_range() {
        let surf = box_bottom_surface();
        let curve = z_line(0.5, 0.5, -1.0, 1.0);
        let ranges = run(curve, surf, (0.0, 1.0, 0.0, 1.0));
        assert_eq!(ranges.len(), 1, "ranges: {ranges:?}");
        // Line from (0.5,0.5,-1) along +Z: the plane z = 0 is met at t = 1.
        let r = ranges[0];
        assert!(r.first <= 1.0 && r.last >= 1.0, "range {r:?} must bracket the crossing");
        assert!((r.first - 1.0).abs() < 1e-3, "first {r:?}");
        assert!((r.last - 1.0).abs() < 1e-3, "last {r:?}");
    }

    #[test]
    fn line_missing_box_face_empty() {
        let surf = box_bottom_surface();
        // Line outside the face's UV domain (x beyond the box).
        let curve = z_line(5.0, 5.0, -1.0, 1.0);
        let ranges = run(curve, surf, (0.0, 1.0, 0.0, 1.0));
        assert!(ranges.is_empty(), "ranges: {ranges:?}");
    }

    #[test]
    fn line_in_box_plane_covers_bean() {
        let surf = box_bottom_surface();
        // Line lying in z = 0, spanning the full bean window.
        let curve = Arc::new(GeomLine::new(GpLin::from_pnt_dir(
            GpPnt::new(0.0, 0.5, 0.0),
            dir(1.0, 0.0, 0.0),
        )));
        let mut bfi = BeanFaceIntersector::new();
        bfi.initialize(curve, surf, 1e-7, 1e-7);
        bfi.set_bean_parameters(-1.0, 2.0);
        bfi.set_surface_parameters(0.0, 1.0, 0.0, 1.0);
        bfi.perform().expect("perform");
        let ranges = bfi.result();
        assert_eq!(ranges.len(), 1, "ranges: {ranges:?}");
        let r = ranges[0];
        assert!((r.first - -1.0).abs() < 1e-9, "first {r:?}");
        assert!((r.last - 2.0).abs() < 1e-9, "last {r:?}");
    }

    #[test]
    fn bspline_crosses_sphere_two_ranges() {
        let surface = unit_sphere();
        // A B-spline through (-2,0,0) and (2,0,0) bowing toward z = 0.5; it
        // enters/exits the unit sphere near x = ±1.
        let poles = vec![
            GpPnt::new(-2.0, 0.0, 0.0),
            GpPnt::new(-0.5, 0.0, 0.6),
            GpPnt::new(0.5, 0.0, 0.6),
            GpPnt::new(2.0, 0.0, 0.0),
        ];
        let curve = Arc::new(
            GeomBSplineCurve::new(poles, vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0], 2).unwrap(),
        );
        let ranges = run(curve, surface, (0.0, 2.0 * PI, -PI / 2.0, PI / 2.0));
        assert!(!ranges.is_empty(), "ranges: {ranges:?}");
        // Every returned parameter must be within tolerance of the sphere.
        let mut bfi = BeanFaceIntersector::new();
        let surface2 = unit_sphere();
        let poles2 = vec![
            GpPnt::new(-2.0, 0.0, 0.0),
            GpPnt::new(-0.5, 0.0, 0.6),
            GpPnt::new(0.5, 0.0, 0.6),
            GpPnt::new(2.0, 0.0, 0.0),
        ];
        let curve2 = Arc::new(
            GeomBSplineCurve::new(poles2, vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0], 2).unwrap(),
        );
        bfi.initialize(curve2, surface2, 1e-7, 1e-7);
        bfi.set_bean_parameters(0.0, 1.0);
        bfi.set_surface_parameters(0.0, 2.0 * PI, -PI / 2.0, PI / 2.0);
        bfi.perform().expect("perform");
        for r in bfi.result() {
            let mut u = 0.0;
            let mut v = 0.0;
            let d = bfi.distance_with_uv(r.first, &mut u, &mut v);
            assert!(d < 1e-4, "range start {r:?} distance {d}");
            let d2 = bfi.distance_with_uv(r.last, &mut u, &mut v);
            assert!(d2 < 1e-4, "range end {r:?} distance {d2}");
        }
        // The minimal squared distance should be (essentially) zero.
        assert!(bfi.minimal_square_distance() < 1e-10, "min_sq {}", bfi.minimal_square_distance());
    }

    #[test]
    fn sphere_missed_by_bspline_is_empty() {
        let surface = unit_sphere();
        // A straight-ish B-spline passing well above the sphere.
        let poles = vec![
            GpPnt::new(-2.0, 0.0, 3.0),
            GpPnt::new(-0.5, 0.0, 3.0),
            GpPnt::new(0.5, 0.0, 3.0),
            GpPnt::new(2.0, 0.0, 3.0),
        ];
        let curve = Arc::new(
            GeomBSplineCurve::new(poles, vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0], 2).unwrap(),
        );
        let ranges = run(curve, surface, (0.0, 2.0 * PI, -PI / 2.0, PI / 2.0));
        assert!(ranges.is_empty(), "ranges: {ranges:?}");
    }

    #[test]
    fn minimal_square_distance_fast_reject_stays_large() {
        // A line at y = 3 against the unit sphere: `FastComputeAnalytic`
        // conclusively rejects the pair (no intersection) and returns early,
        // leaving `minimal_square_distance()` at its sentinel — exactly as OCCT
        // leaves `myMinSqDistance` at `RealLast()` on the fast path.
        let surface = unit_sphere();
        let curve = Arc::new(GeomLine::new(GpLin::from_pnt_dir(
            GpPnt::new(0.0, 3.0, 0.0),
            dir(1.0, 0.0, 0.0),
        )));
        let mut bfi = BeanFaceIntersector::new();
        bfi.initialize(curve, surface, 1e-7, 1e-7);
        bfi.set_bean_parameters(-3.0, 3.0);
        bfi.set_surface_parameters(0.0, 2.0 * PI, -PI / 2.0, PI / 2.0);
        bfi.perform().expect("perform");
        assert!(bfi.result().is_empty(), "line at y=3 misses the unit sphere");
        assert!(bfi.minimal_square_distance() > 1e100, "fast-reject leaves the sentinel");
    }

    #[test]
    fn circle_in_plane_full_range() {
        let ax3 = GpAx3::new(GpPnt::zero(), dir(0.0, 0.0, 1.0), &dir(1.0, 0.0, 0.0)).unwrap();
        let plane = Arc::new(GeomPlane::new(GpPln::new(ax3))) as Arc<dyn Surface>;
        let circle = Arc::new(occt_geom::GeomCircle::new(GpCirc::new(GpAx2::standard(), 1.0)))
            as Arc<dyn Curve>;
        let mut bfi = BeanFaceIntersector::new();
        bfi.initialize(circle, plane, 1e-7, 1e-7);
        bfi.set_bean_parameters(0.0, 2.0 * PI);
        bfi.set_surface_parameters(-5.0, 5.0, -5.0, 5.0);
        bfi.perform().expect("perform");
        let ranges = bfi.result();
        assert_eq!(ranges.len(), 1, "ranges: {ranges:?}");
        assert!((ranges[0].first - 0.0).abs() < 1e-9);
        assert!((ranges[0].last - 2.0 * PI).abs() < 1e-9);
    }

    #[test]
    fn line_crosses_bspline_surface_localized_path() {
        // A curved BSpline patch (poles at (u, v, u²) → z ≈ u²) over [0,1]²,
        // classified as a NURBS surface so the localized path runs. A vertical
        // line at (0.5, 0.5) crosses it somewhere over the patch; the result
        // must be one range whose midpoint lies on the surface.
        use occt_geom::bspline_surface::{bspline_surface_uniform_knots, GeomBSplineSurface};
        let (ku, kv) = bspline_surface_uniform_knots(3, 3, 2, 2);
        let poles: Vec<Vec<GpPnt>> = (0..3)
            .map(|i| {
                let u = i as f64 * 0.5;
                (0..3)
                    .map(|j| {
                        let v = j as f64 * 0.5;
                        GpPnt::new(u, v, u * u)
                    })
                    .collect()
            })
            .collect();
        let surface = Arc::new(GeomBSplineSurface::new(poles, ku, kv, 2, 2).unwrap()) as Arc<dyn Surface>;
        let curve = Arc::new(GeomLine::new(GpLin::from_pnt_dir(
            GpPnt::new(0.5, 0.5, -1.0),
            dir(0.0, 0.0, 1.0),
        )));
        let mut bfi = BeanFaceIntersector::new();
        bfi.initialize(curve, surface, 1e-7, 1e-7);
        bfi.set_bean_parameters(-1.0, 3.0);
        bfi.set_surface_parameters(0.0, 1.0, 0.0, 1.0);
        bfi.perform().expect("perform succeeds on BSpline surface");
        let ranges = bfi.result();
        assert_eq!(ranges.len(), 1, "ranges: {ranges:?}");
        let r = ranges[0];
        // The crossing point must genuinely lie on the surface.
        let mid = 0.5 * (r.first + r.last);
        let p = GeomLine::new(GpLin::from_pnt_dir(GpPnt::new(0.5, 0.5, -1.0), dir(0.0, 0.0, 1.0)))
            .d0(mid);
        let (_, _, d) = bfi.closest_params_dist(&p);
        assert!(d < 1e-4, "crossing at {mid} is {d} off the BSpline patch");
    }

    #[test]
    fn edge_face_cross_check_matches_approximate_solver() {
        // Line through the box bottom face: compare with the existing
        // approximate `inttools::edge_face_intersections`.
        let bx = unit_box();
        let face = &bx.faces[0];
        let b = crate::builder::TopoBuilder::new();
        let edge = b.make_edge_segment(&GpPnt::new(0.5, 0.5, -1.0), &GpPnt::new(0.5, 0.5, 1.0));
        let approx = crate::inttools::edge_face_intersections(&edge, face, 1e-9);
        assert_eq!(approx.len(), 1, "approx: {approx:?}");
        assert!(approx[0].1.distance(&GpPnt::new(0.5, 0.5, 0.0)) < 1e-6);

        let surf = BRepTool::face_surface(face).expect("face surface");
        let curve = Arc::new(GeomLine::new(GpLin::from_pnt_dir(
            GpPnt::new(0.5, 0.5, -1.0),
            dir(0.0, 0.0, 1.0),
        )));
        let mut bfi = BeanFaceIntersector::new();
        bfi.initialize(curve, surf, 1e-7, 1e-7);
        bfi.set_bean_parameters(-10.0, 10.0);
        bfi.set_surface_parameters(0.0, 1.0, 0.0, 1.0);
        bfi.perform().expect("perform");
        let ranges = bfi.result();
        assert_eq!(ranges.len(), 1, "ranges: {ranges:?}");
        // Approximate solver hits at edge parameter 1.0 (z = 0); the bean's
        // crossing is at curve parameter 1.0 too. Compare the 3D points.
        let mid = 0.5 * (ranges[0].first + ranges[0].last);
        let p = GeomLine::new(GpLin::from_pnt_dir(GpPnt::new(0.5, 0.5, -1.0), dir(0.0, 0.0, 1.0)))
            .d0(mid);
        assert!(p.distance(&approx[0].1) < 1e-3, "bean midpoint {p:?} vs approx {:?}", approx[0].1);
    }
}
