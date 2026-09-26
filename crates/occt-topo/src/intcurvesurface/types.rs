use super::prelude::*;


/// Tolerance for the analytic path (exact solves).

pub(super) const ANALYTIC_TOL: f64 = 1e-9;
/// Tolerance for the general sampling path.
pub(super) const GENERAL_TOL: f64 = 1e-6;
/// Parameter-merge window: points closer than this are merged as a tangency.
pub(super) const MERGE_TOL: f64 = 1e-7;
/// Angular tolerance for transition classification (mirrors `Precision::Angular`).
pub(super) const ANG_TOL: f64 = 1e-12;

// ---------------------------------------------------------------------------
// Data classes
// ---------------------------------------------------------------------------

/// Transition of the curve across the surface at an intersection point.
/// Mirrors `TopAbs_State`: `In` = entering the material, `Out` = leaving,
/// `On` = tangent (the curve only touches the surface), `Unknown` = undecided.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum State {
    In,
    Out,
    On,
    Unknown,
}

/// A point common to the curve and the surface.
///
/// Mirrors `IntCurveSurface_IntersectionPoint`: the parameter on the curve
/// (`param`, OCCT `W`), the surface parameters (`u`, `v`), the 3D point and
/// the curve transition at the crossing.
#[derive(Debug, Clone, Copy)]
pub struct IntersectionPoint {
    pub param: f64,
    pub u: f64,
    pub v: f64,
    pub pnt: GpPnt,
    pub state: State,
}

impl IntersectionPoint {
    pub fn new(param: f64, u: f64, v: f64, pnt: GpPnt, state: State) -> Self {
        Self { param, u, v, pnt, state }
    }
    /// Parameter on the curve.
    pub fn param(&self) -> f64 {
        self.param
    }
    /// U parameter on the surface.
    pub fn u(&self) -> f64 {
        self.u
    }
    /// V parameter on the surface.
    pub fn v(&self) -> f64 {
        self.v
    }
    /// The 3D intersection point.
    pub fn pnt(&self) -> GpPnt {
        self.pnt
    }
    /// Transition of the curve at the point.
    pub fn state(&self) -> State {
        self.state
    }
}

/// A curve span `[first, last]` lying on the surface (the coincidence case,
/// e.g. a line in a plane or a circle on a sphere).
///
/// Mirrors `IntCurveSurface_IntersectionSegment`.
#[derive(Debug, Clone, Copy)]
pub struct IntersectionSegment {
    pub(super) first: IntersectionPoint,
    pub(super) last: IntersectionPoint,
}

impl IntersectionSegment {
    pub fn new(first: IntersectionPoint, last: IntersectionPoint) -> Self {
        Self { first, last }
    }
    pub fn first_point(&self) -> &IntersectionPoint {
        &self.first
    }
    pub fn second_point(&self) -> &IntersectionPoint {
        &self.last
    }
}

/// Result of a curve–surface intersection: discrete points plus (degenerate)
/// coincident segments.
#[derive(Debug, Clone, Default)]
pub struct HInterResult {
    pub(super) points: Vec<IntersectionPoint>,
    pub(super) segments: Vec<IntersectionSegment>,
}

impl HInterResult {
    /// Number of isolated intersection points.
    pub fn nb_points(&self) -> usize {
        self.points.len()
    }
    /// The `i`-th intersection point (0-based).
    pub fn point(&self, i: usize) -> &IntersectionPoint {
        &self.points[i]
    }
    /// All intersection points.
    pub fn points(&self) -> &[IntersectionPoint] {
        &self.points
    }
    /// Number of coincidence segments.
    pub fn nb_segments(&self) -> usize {
        self.segments.len()
    }
    /// The `i`-th coincidence segment (0-based).
    pub fn segment(&self, i: usize) -> &IntersectionSegment {
        &self.segments[i]
    }
    /// All coincidence segments.
    pub fn segments(&self) -> &[IntersectionSegment] {
        &self.segments
    }
}

// ---------------------------------------------------------------------------
// Curve classification
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum CurveKind {
    Line,
    Circle,
    Ellipse,
    Other,
}

pub(super) fn classify_curve(c: &dyn Curve) -> CurveKind {
    if is_line_like(c) {
        return CurveKind::Line;
    }
    if c.is_periodic() && (c.period() - 2.0 * std::f64::consts::PI).abs() < 1e-9 {
        if is_circle_like(c) {
            CurveKind::Circle
        } else {
            CurveKind::Ellipse
        }
    } else {
        CurveKind::Other
    }
}

/// Whether the curve is geometrically a straight line: either unbounded (the
/// port's `GeomLine`) or with collinear samples.
pub(super) fn is_line_like(c: &dyn Curve) -> bool {
    let (a, b) = (c.first_parameter(), c.last_parameter());
    if a.is_finite() && b.is_finite() && (b - a).abs() <= 1e-15 {
        return false;
    }
    // Choose a sampling window: the natural range, or [-2, 2] for unbounded
    // curves (enough to expose a parabola/hyperbola's curvature).
    let (lo, hi) = if a.is_finite() && b.is_finite() {
        (a, b)
    } else {
        (-2.0, 2.0)
    };
    let p0 = c.d0(lo);
    let pl = c.d0(hi);
    let size = p0.distance(&pl);
    if size <= 1e-30 {
        return false;
    }
    let d0 = GpVec::from_pnts(&p0, &pl);
    let tol = 1e-7 * size;
    for i in 1..8 {
        let u = lo + (hi - lo) * i as f64 / 8.0;
        let p = c.d0(u);
        if GpVec::from_pnts(&p0, &p).crossed(&d0).magnitude() > tol * size {
            return false;
        }
    }
    true
}

/// Whether the curve is a circle: all samples are coplanar and equidistant
/// from one center.
pub(super) fn is_circle_like(c: &dyn Curve) -> bool {
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
    let nvec = nrm.divided(m);
    let scale = radius.max(1.0);
    let tol = 1e-6 * scale;
    for p in &pts {
        let v = GpVec::from_pnts(&center, p);
        if v.dot(&nvec).abs() > tol {
            return false;
        }
        if (v.magnitude() - radius).abs() > tol {
            return false;
        }
    }
    true
}

/// First three non-collinear points of a sample set.
pub(super) fn first_three_spanning(pts: &[GpPnt]) -> Option<(GpPnt, GpPnt, GpPnt)> {
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

pub(super) fn det3(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

pub(super) fn solve3(a: &[[f64; 3]; 3], rhs: &[f64; 3]) -> Option<[f64; 3]> {
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

/// Circumcenter of three non-collinear points (perpendicular-bisector system).
pub(super) fn circumcenter(a: &GpPnt, b: &GpPnt, c: &GpPnt) -> Option<GpPnt> {
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

pub(super) fn midpoint(a: &GpPnt, b: &GpPnt) -> GpPnt {
    GpPnt::new(0.5 * (a.x() + b.x()), 0.5 * (a.y() + b.y()), 0.5 * (a.z() + b.z()))
}

/// Reconstructed conic geometry. For a natural-parameter periodic 2π curve,
/// `d0(u) = o + x·r·cos(u) + y·r·sin(u)` (circle) or
/// `d0(u) = o + x·a·cos(u) − y·b·sin(u)` (ellipse).
pub(super) enum ConicGeom {
    Circle { o: GpPnt, r: f64, x: GpVec, y: GpVec },
    Ellipse { o: GpPnt, a: f64, b: f64, x: GpVec, y: GpVec },
}

/// Reconstruct a natural-parameter conic from the curve's `d0`, verifying the
/// parameterization so a non-natural curve (e.g. a trimmed reparametrization)
/// falls back to the general path.
pub(super) fn conic_geometry(c: &dyn Curve, kind: CurveKind) -> Option<ConicGeom> {
    match kind {
        CurveKind::Circle => {
            let p0 = c.d0(0.0);
            let p1 = c.d0(std::f64::consts::FRAC_PI_2);
            let p2 = c.d0(std::f64::consts::PI);
            let o = circumcenter(&p0, &p1, &p2)?;
            let r = p0.distance(&o);
            if r <= 1e-30 {
                return None;
            }
            let x = GpVec::from_pnts(&o, &p0).divided(r);
            let y = GpVec::from_pnts(&o, &p1).divided(r);
            // Verify natural parameterization at a few samples.
            for &u in &[0.3f64, 1.1f64, 5.0f64] {
                let expect = o.translated_vec(&x.multiplied_scalar(r * u.cos()))
                    .translated_vec(&y.multiplied_scalar(r * u.sin()));
                if c.d0(u).distance(&expect) > 1e-6 * r.max(1.0) {
                    return None;
                }
            }
            Some(ConicGeom::Circle { o, r, x, y })
        }
        CurveKind::Ellipse => {
            let p0 = c.d0(0.0);
            let pm = c.d0(std::f64::consts::PI);
            let pq = c.d0(std::f64::consts::FRAC_PI_2);
            let o = midpoint(&p0, &pm);
            let a = p0.distance(&o);
            let b = pq.distance(&o);
            if a <= 1e-30 || b <= 1e-30 {
                return None;
            }
            let x = GpVec::from_pnts(&o, &p0).divided(a);
            let y = GpVec::from_pnts(&o, &pq).divided(b); // P(pi/2) = o + y*b
            for &u in &[0.3f64, 1.1f64, 5.0f64] {
                let expect = o.translated_vec(&x.multiplied_scalar(a * u.cos()))
                    .translated_vec(&y.multiplied_scalar(b * u.sin()));
                if c.d0(u).distance(&expect) > 1e-6 * a.max(b).max(1.0) {
                    return None;
                }
            }
            Some(ConicGeom::Ellipse { o, a, b, x, y })
        }
        _ => None,
    }
}

/// Reconstruct a line as `d0(u) = loc + v·u` (v need not be unit: for a
/// `GeomLine` it is the unit direction; for a trimmed/other line-like curve it
/// is the chord direction scaled to the parameter range).
pub(super) fn line_geometry(c: &dyn Curve) -> Option<(GpPnt, GpVec)> {
    let p0 = c.d0(0.0);
    let p1 = c.d0(1.0);
    let v = GpVec::from_pnts(&p0, &p1);
    if v.magnitude() < 1e-30 {
        return None;
    }
    Some((p0, v))
}

// ---------------------------------------------------------------------------
// Surface classification and reconstruction
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum SurfaceKind {
    Plane,
    Sphere,
    Cylinder,
    Cone,
    Torus,
    Other,
}

/// Reconstructed quadric surface geometry. Each variant carries the data needed
/// for the implicit-equation solves and for exact `(u, v)` parameter recovery.
#[derive(Debug, Clone)]
pub(super) enum SurfaceGeom {
    Plane { o: GpPnt, x: GpVec, y: GpVec, n: GpVec },
    Sphere { o: GpPnt, r: f64, x: GpVec, y: GpVec, z: GpVec },
    Cylinder { a: GpPnt, z: GpVec, r: f64 },
    Cone { apex: GpPnt, z: GpVec, cosa: f64 },
    Torus { o: GpPnt, major: f64, minor: f64, x: GpVec, z: GpVec },
}

/// Finite, sane sampling bounds (unbounded directions clamp to ±1).
pub(super) fn sample_bounds(s: &dyn Surface) -> (f64, f64, f64, f64) {
    let (u0, u1) = s.u_range();
    let (v0, v1) = s.v_range();
    let clamp = |a: f64, b: f64| if a.is_finite() && b.is_finite() && b > a { (a, b) } else { (-1.0, 1.0) };
    let (u0, u1) = clamp(u0, u1);
    let (v0, v1) = clamp(v0, v1);
    (u0, u1, v0, v1)
}

/// The natural `(0, 0)`-anchored sampling point for a quadric's frame axes
/// (all ported quadrics start their periodic U at 0).
pub(super) fn quadric_origin_u(s: &dyn Surface) -> f64 {
    let (u0, _) = s.u_range();
    if u0.is_finite() { u0 } else { 0.0 }
}

pub(super) fn classify_surface_kind(s: &dyn Surface) -> SurfaceKind {
    if is_planar(s, 8, 8, 1e-6) {
        return SurfaceKind::Plane;
    }
    if let Some(geom) = build_sphere_geom(s) {
        if sphere_matches(s, &geom) {
            return SurfaceKind::Sphere;
        }
    }
    if let Some(geom) = build_cylinder_geom(s) {
        if cylinder_matches(s, &geom) {
            return SurfaceKind::Cylinder;
        }
    }
    if let Some(geom) = build_cone_geom(s) {
        if cone_matches(s, &geom) {
            return SurfaceKind::Cone;
        }
    }
    if let Some(geom) = build_torus_geom(s) {
        if torus_matches(s, &geom) {
            return SurfaceKind::Torus;
        }
    }
    SurfaceKind::Other
}

pub(super) fn build_surface_geom(s: &dyn Surface, kind: SurfaceKind) -> Option<SurfaceGeom> {
    match kind {
        SurfaceKind::Plane => {
            let u0 = match s.u_range().0 {
                a if a.is_finite() => a,
                _ => 0.0,
            };
            let v0 = match s.v_range().0 {
                a if a.is_finite() => a,
                _ => 0.0,
            };
            let o = s.d0(u0, v0);
            let x = GpVec::from_pnts(&o, &s.d0(u0 + 1.0, v0));
            let y = GpVec::from_pnts(&o, &s.d0(u0, v0 + 1.0));
            let n = x.crossed(&y);
            if n.magnitude() < 1e-12 {
                return None;
            }
            Some(SurfaceGeom::Plane { o, x, y, n })
        }
        SurfaceKind::Sphere => build_sphere_geom(s),
        SurfaceKind::Cylinder => build_cylinder_geom(s),
        SurfaceKind::Cone => build_cone_geom(s),
        SurfaceKind::Torus => build_torus_geom(s),
        SurfaceKind::Other => None,
    }
}

/// Sphere: center from `sphere_center`, radius from a sample, orthonormal
/// frame from the natural (0, 0), (π/2, 0), (0, π/2) samples.
pub(super) fn build_sphere_geom(s: &dyn Surface) -> Option<SurfaceGeom> {
    let o = sphere_center(s)?;
    let u0 = quadric_origin_u(s);
    let r = s.d0(u0, 0.0).distance(&o);
    if r < 1e-12 {
        return None;
    }
    let x = GpVec::from_pnts(&o, &s.d0(u0, 0.0)).divided(r);
    let y = GpVec::from_pnts(&o, &s.d0(u0 + std::f64::consts::FRAC_PI_2, 0.0)).divided(r);
    let z = GpVec::from_pnts(&o, &s.d0(u0, std::f64::consts::FRAC_PI_2)).divided(r);
    if !frame_ok(&x, &y, &z) {
        return None;
    }
    Some(SurfaceGeom::Sphere { o, r, x, y, z })
}

/// Verify a candidate sphere reconstruction against a sample grid.
pub(super) fn sphere_matches(s: &dyn Surface, geom: &SurfaceGeom) -> bool {
    let SurfaceGeom::Sphere { o, r, .. } = geom else { return false; };
    let (u0, u1, v0, v1) = sample_bounds(s);
    let tol = 1e-4 * r.abs().max(1.0);
    for i in 0..8 {
        for j in 0..8 {
            let u = u0 + (u1 - u0) * i as f64 / 7.0;
            let v = v0 + (v1 - v0) * j as f64 / 7.0;
            if (s.d0(u, v).distance(o) - r).abs() > tol {
                return false;
            }
        }
    }
    true
}

/// Cylinder: axis point = circle-center at the v start, axis direction from the
/// v advance, radius = half the opposite-sample distance.
pub(super) fn build_cylinder_geom(s: &dyn Surface) -> Option<SurfaceGeom> {
    let u0 = quadric_origin_u(s);
    let (_, _, v0, _) = sample_bounds(s);
    let c0 = midpoint(&s.d0(u0, v0), &s.d0(u0 + std::f64::consts::PI, v0));
    let c1 = midpoint(&s.d0(u0, v0 + 1.0), &s.d0(u0 + std::f64::consts::PI, v0 + 1.0));
    let zvec = GpVec::from_pnts(&c0, &c1);
    let zm = zvec.magnitude();
    if zm < 1e-12 {
        return None;
    }
    let z = zvec.divided(zm);
    let r = s.d0(u0, v0).distance(&c0);
    if r < 1e-12 {
        return None;
    }
    let x = GpVec::from_pnts(&c0, &s.d0(u0, v0)).divided(r);
    let y = GpVec::from_pnts(&c0, &s.d0(u0 + std::f64::consts::FRAC_PI_2, v0)).divided(r);
    // Reject a degenerate frame (x, y not perpendicular to z).
    if x.crossed(&y).dot(&z).abs() < 1e-6 {
        return None;
    }
    Some(SurfaceGeom::Cylinder { a: c0, z, r })
}

/// Verify a candidate cylinder reconstruction against a sample grid.
pub(super) fn cylinder_matches(s: &dyn Surface, geom: &SurfaceGeom) -> bool {
    let SurfaceGeom::Cylinder { a, z, r, .. } = geom else { return false; };
    let (u0, u1, v0, v1) = sample_bounds(s);
    let tol = 1e-4 * r.abs().max(1.0);
    for i in 0..8 {
        for j in 0..8 {
            let u = u0 + (u1 - u0) * i as f64 / 7.0;
            let v = v0 + (v1 - v0) * j as f64 / 7.0;
            let d = GpVec::from_pnts(a, &s.d0(u, v));
            let rho = d.coord.subtracted(&z.xyz().multiplied(d.dot(z))).modulus();
            if (rho - r).abs() > tol {
                return false;
            }
        }
    }
    true
}

/// Cone: axis from the circle-center line, apex at the extrapolated zero
/// radius, semi-angle from the radius growth rate.
pub(super) fn build_cone_geom(s: &dyn Surface) -> Option<SurfaceGeom> {
    let u0 = quadric_origin_u(s);
    let (_, _, v0, _) = sample_bounds(s);
    let c0 = midpoint(&s.d0(u0, v0), &s.d0(u0 + std::f64::consts::PI, v0));
    let c1 = midpoint(&s.d0(u0, v0 + 1.0), &s.d0(u0 + std::f64::consts::PI, v0 + 1.0));
    let zvec = GpVec::from_pnts(&c0, &c1);
    let zm = zvec.magnitude();
    if zm < 1e-12 {
        return None;
    }
    let z = zvec.divided(zm);
    let r0 = s.d0(u0, v0).distance(&c0);
    let r1 = s.d0(u0, v0 + 1.0).distance(&c1);
    let tan_alpha = (r1 - r0) / zm;
    if tan_alpha.abs() < 1e-12 {
        return None;
    }
    let cosa = 1.0 / (1.0 + tan_alpha * tan_alpha).sqrt();
    // Distance from apex to the v0 circle plane along the axis: r0 / tanα.
    let s0 = r0 / tan_alpha;
    let apex = c0.translated_vec(&z.multiplied_scalar(-s0));
    Some(SurfaceGeom::Cone { apex, z, cosa })
}

/// Verify a candidate cone reconstruction against a sample grid.
pub(super) fn cone_matches(s: &dyn Surface, geom: &SurfaceGeom) -> bool {
    let SurfaceGeom::Cone { apex, z, cosa, .. } = geom else { return false; };
    let (u0, u1, v0, v1) = sample_bounds(s);
    for i in 0..8 {
        for j in 0..8 {
            let u = u0 + (u1 - u0) * i as f64 / 7.0;
            let v = v0 + (v1 - v0) * j as f64 / 7.0;
            let d = GpVec::from_pnts(apex, &s.d0(u, v));
            let err = (d.dot(z) * d.dot(z) - d.square_magnitude() * cosa * cosa).abs();
            if err > 1e-4 * d.square_magnitude().max(1.0) {
                return false;
            }
        }
    }
    true
}

/// Torus: center from the v = 0 circle, major/minor radii from the v = 0 and
/// v = π circle radii, frame from the natural parameterization.
pub(super) fn build_torus_geom(s: &dyn Surface) -> Option<SurfaceGeom> {
    let u0 = quadric_origin_u(s);
    let pi = std::f64::consts::PI;
    let o = midpoint(&s.d0(u0, 0.0), &s.d0(u0 + pi, 0.0));
    let r_outer = s.d0(u0, 0.0).distance(&o);
    let o2 = midpoint(&s.d0(u0, pi), &s.d0(u0 + pi, pi));
    let r_inner = s.d0(u0, pi).distance(&o2);
    if o.distance(&o2) > 1e-6 * r_outer.max(1.0) {
        return None;
    }
    let major = 0.5 * (r_outer + r_inner);
    let minor = 0.5 * (r_outer - r_inner);
    if major <= 1e-12 || minor <= 1e-12 {
        return None;
    }
    let x = GpVec::from_pnts(&o, &s.d0(u0, 0.0)).divided(r_outer);
    let y = GpVec::from_pnts(&o, &s.d0(u0 + pi / 2.0, 0.0)).divided(r_outer);
    let z = GpVec::from_pnts(
        &o,
        &s.d0(u0, pi / 2.0),
    )
    .subtracted(&x.multiplied_scalar(major))
    .divided(minor);
    if !frame_ok(&x, &y, &z) {
        return None;
    }
    Some(SurfaceGeom::Torus { o, major, minor, x, z })
}

/// Verify a candidate torus reconstruction against a sample grid.
pub(super) fn torus_matches(s: &dyn Surface, geom: &SurfaceGeom) -> bool {
    let SurfaceGeom::Torus { o, major, minor, z, .. } = geom else { return false; };
    let (u0, u1, v0, v1) = sample_bounds(s);
    let tol = 1e-4 * (*major).max(*minor).max(1.0);
    for i in 0..8 {
        for j in 0..8 {
            let u = u0 + (u1 - u0) * i as f64 / 7.0;
            let v = v0 + (v1 - v0) * j as f64 / 7.0;
            let d = GpVec::from_pnts(o, &s.d0(u, v));
            let dz = d.dot(z);
            let rho = d.coord.subtracted(&z.xyz().multiplied(dz)).modulus();
            let err = ((rho - major) * (rho - major) + dz * dz - minor * minor).abs();
            if err > tol * (rho + major + minor).max(1.0) {
                return false;
            }
        }
    }
    true
}

/// Whether x, y, z form an orthonormal frame (up to scale/reflection noise).
pub(super) fn frame_ok(x: &GpVec, y: &GpVec, z: &GpVec) -> bool {
    (x.magnitude() - 1.0).abs() < 1e-6
        && (y.magnitude() - 1.0).abs() < 1e-6
        && (z.magnitude() - 1.0).abs() < 1e-6
        && x.dot(y).abs() < 1e-6
        && x.dot(z).abs() < 1e-6
        && y.dot(z).abs() < 1e-6
}

/// Wrap `x` into `[start, start + period)` (positive modulo).
pub(super) fn wrap_periodic(x: f64, start: f64, period: f64) -> f64 {
    let p = period.abs();
    if p <= 1e-30 {
        return x;
    }
    start + (x - start).rem_euclid(p)
}

/// Snap a cosine/sine value to ±1 when within `eps` of a pole, so `asin` does
/// not amplify a ~1e-16 frame error into a ~1e-8 angular error.
pub(super) fn snap_pole(s: f64) -> f64 {
    if s > 1.0 - 1e-9 {
        1.0
    } else if s < -1.0 + 1e-9 {
        -1.0
    } else {
        s
    }
}
