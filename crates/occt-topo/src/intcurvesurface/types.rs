use super::prelude::*;


/// Tolerance for the analytic path (exact solves).

pub(super) const ANALYTIC_TOL: f64 = 1e-9;
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

/// `Adaptor3d_Curve::GetType()` as reported by `GeomAdaptor_Curve` (a
/// `Geom_TrimmedCurve` reports the type of its basis curve).
///
/// `PerformBounds` (`IntCurveSurface_Inter.pxx:119-137`) switches on it: Line /
/// Circle / Ellipse / Parabola / Hyperbola go to `PerformConicSurf`, everything
/// else takes the polygon/polyhedron or curve-quadric arm.
///
/// The kind is read from the curve's own `gp_*` geometry (the same source as
/// `CurveTool::Line/Circle/Ellipse/Parabola/Hyperbola`), never from geometric
/// samples: a `Geom_BSplineCurve` that happens to be a straight line reports
/// `GeomAbs_BSplineCurve` and must take the general arm.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum CurveKind {
    Line,
    Circle,
    Ellipse,
    /// `PerformConicSurfParabola` (`IntCurveSurface_Inter.pxx:811-887`) routes
    /// Plane / Cylinder / Cone / Sphere through `IntAna_IntConicQuad`.
    Parabola,
    /// `PerformConicSurfHyperbola` (`IntCurveSurface_Inter.pxx:890-966`), same
    /// dispatch as [`CurveKind::Parabola`].
    Hyperbola,
    Other,
}

pub(super) fn classify_curve(c: &dyn Curve) -> CurveKind {
    if c.gp_line().is_some() {
        CurveKind::Line
    } else if c.gp_circ().is_some() {
        CurveKind::Circle
    } else if c.gp_ellipse().is_some() {
        CurveKind::Ellipse
    } else if c.gp_parabola().is_some() {
        CurveKind::Parabola
    } else if c.gp_hyperbola().is_some() {
        CurveKind::Hyperbola
    } else {
        CurveKind::Other
    }
}

/// Reconstructed conic geometry. For a natural-parameter periodic 2π curve,
/// `d0(u) = o + x·r·cos(u) + y·r·sin(u)` (circle) or
/// `d0(u) = o + x·a·cos(u) − y·b·sin(u)` (ellipse).
pub(super) enum ConicGeom {
    Circle { o: GpPnt, r: f64, x: GpVec, y: GpVec },
    Ellipse { o: GpPnt, a: f64, b: f64, x: GpVec, y: GpVec },
}

/// Read the natural-parameter conic geometry from the curve's own `gp_Circ` /
/// `gp_Elips` (the same source as `CurveTool::Circle` / `CurveTool::Ellipse`),
/// so `d0(u) = o + x·r·cos(u) + y·r·sin(u)` is `ElCLib::CircleValue` /
/// `ElCLib::EllipseValue`.
///
/// A reparametrized wrapper (whose `d0` is not the natural parameterization)
/// is still rejected by re-evaluating `d0` against it, so such a curve takes
/// the general arm.
pub(super) fn conic_geometry(c: &dyn Curve, kind: CurveKind) -> Option<ConicGeom> {
    match kind {
        CurveKind::Circle => {
            let circ = c.gp_circ()?;
            let pos = circ.position();
            let o = pos.location();
            let r = circ.radius();
            if r <= 1e-30 {
                return None;
            }
            let x = GpVec::from_xyz(pos.x_direction().xyz());
            let y = GpVec::from_xyz(pos.y_direction().xyz());
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
            let el = c.gp_ellipse()?;
            let pos = el.position();
            let o = pos.location();
            let a = el.major_radius();
            let b = el.minor_radius();
            if a <= 1e-30 || b <= 1e-30 {
                return None;
            }
            let x = GpVec::from_xyz(pos.x_direction().xyz());
            let y = GpVec::from_xyz(pos.y_direction().xyz());
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

/// The line's own `gp_Lin` (the same source as `CurveTool::Line`), evaluated as
/// `d0(u) = loc + u·dir` so a returned parameter is `ElCLib::LineParameter`.
pub(super) fn line_geometry(c: &dyn Curve) -> Option<(GpPnt, GpVec)> {
    let lin = c.gp_line()?;
    let v = GpVec::from_xyz(lin.direction().xyz());
    if v.magnitude() < 1e-30 {
        return None;
    }
    Some((lin.location(), v))
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
    Plane { o: GpPnt, n: GpVec },
    Sphere { o: GpPnt, r: f64 },
    Cylinder { a: GpPnt, z: GpVec, r: f64 },
    Cone { apex: GpPnt, z: GpVec, cosa: f64 },
    Torus { o: GpPnt, major: f64, minor: f64, x: GpVec, z: GpVec },
}

/// `Adaptor3d_Surface::GetType()` for the dispatcher: the concrete geometry
/// type, most specific first. `PerformBounds` (`IntCurveSurface_Inter.pxx:139-179`)
/// tests `SurfaceTool::GetType(theSurface)` against `GeomAbs_{Plane,Cylinder,
/// Cone,Sphere}`; every other type takes the polygon/polyhedron arm. A surface
/// that is only *geometrically* a quadric (a `Geom_BSplineSurface` patch) must
/// stay `Other`, exactly as `GetType` reports it — `build_surface_geom` reads
/// the analytic frame from the surface's own `gp_*` placement, which such a
/// patch does not have. Uses the
/// faithful `GeomBndLib_Surface::initFromSurface` classifier (it also recurses
/// into `Geom_RectangularTrimmedSurface` bases).
pub(super) fn classify_surface_kind(s: &dyn Surface) -> SurfaceKind {
    use crate::geom_bnd_lib_surface3d::SurfaceKind as TypeKind;
    match crate::geom_bnd_lib_surface3d::surface_kind(s) {
        TypeKind::Plane => SurfaceKind::Plane,
        TypeKind::Sphere => SurfaceKind::Sphere,
        TypeKind::Cylinder => SurfaceKind::Cylinder,
        TypeKind::Cone => SurfaceKind::Cone,
        TypeKind::Torus => SurfaceKind::Torus,
        _ => SurfaceKind::Other,
    }
}

/// Reconstruct the analytic quadric reference from the surface's own `gp_*`
/// placement — the same source as OCCT's
/// `SurfaceTool::Plane/Cylinder/Cone/Sphere/Torus(theSurface)`.
///
/// `classify_surface_kind` reports a kind only when the matching `gp_*` is
/// present (it recurses into a `Geom_RectangularTrimmedSurface` base), so every
/// arm reads the exact placement and no geometric sampling is involved.
pub(super) fn build_surface_geom(s: &dyn Surface, kind: SurfaceKind) -> Option<SurfaceGeom> {
    match kind {
        SurfaceKind::Plane => {
            let pln = s.gp_pln()?;
            let pos = pln.position();
            // `ElSLib::PlaneValue(U, V) = O + U·XDir + V·YDir`; normal = ZDir.
            Some(SurfaceGeom::Plane { o: pos.location(), n: GpVec::from_xyz(pos.direction().xyz()) })
        }
        SurfaceKind::Sphere => {
            let sph = s.gp_sphere()?;
            let pos = sph.position();
            // `ElSLib::SphereValue` is centered on the placement location.
            Some(SurfaceGeom::Sphere { o: pos.location(), r: sph.radius() })
        }
        SurfaceKind::Cylinder => {
            let cyl = s.gp_cylinder()?;
            let pos = cyl.position();
            // `ElSLib::CylinderValue` is axis-relative; ZDir is unit.
            Some(SurfaceGeom::Cylinder {
                a: pos.location(),
                z: GpVec::from_xyz(pos.direction().xyz()),
                r: cyl.radius(),
            })
        }
        SurfaceKind::Cone => {
            let cone = s.gp_cone()?;
            let pos = cone.position();
            let r = cone.radius();
            let alpha = cone.semi_angle();
            let tan_alpha = alpha.tan();
            if tan_alpha.abs() < 1e-12 {
                return None;
            }
            let z = GpVec::from_xyz(pos.direction().xyz());
            // `ElSLib::ConeValue` puts `v = 0` at `RefRadius` and grows the
            // radius on the +v side, so the apex is at `v = -R/sin(alpha)`,
            // i.e. `Location - Dir·(R/tan(alpha))`.
            let apex = pos.location().translated_vec(&z.multiplied_scalar(-r / tan_alpha));
            Some(SurfaceGeom::Cone { apex, z, cosa: alpha.cos() })
        }
        SurfaceKind::Torus => {
            let tor = s.gp_torus()?;
            let pos = tor.position();
            // `ElSLib::TorusValue(U, V) = O + (R + r·cos V)(XDir·cos U + YDir·sin U) + r·sin V·ZDir`.
            Some(SurfaceGeom::Torus {
                o: pos.location(),
                major: tor.major_radius(),
                minor: tor.minor_radius(),
                x: GpVec::from_xyz(pos.x_direction().xyz()),
                z: GpVec::from_xyz(pos.direction().xyz()),
            })
        }
        SurfaceKind::Other => None,
    }
}
