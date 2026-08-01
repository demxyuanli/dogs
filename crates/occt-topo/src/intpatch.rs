//! Surface–surface intersection curves.
//! Source: `IntPatch_Intersection`, `IntAna_IntConicQuad`, `IntWalk_PWalking` (TKMath / TKGeomAlgo).
//!
//! Ports OCCT's surface-intersection machinery to exact analytic closed forms
//! where the pair is amenable (plane∩sphere, sphere∩sphere, plane∩cylinder)
//! and a marching-squares grid tracer for the general case. The result of
//! intersecting two parametric surfaces is one or more 3D curves together with
//! the `(u, v)` parameters of each sampled curve point on both surfaces — the
//! data `BRepAlgoAPI_Section` needs to build p-curves on the trimmed faces.

use std::sync::Arc;

use occt_core::geom::polyline_simplify::rdp_simplify;
use occt_core::gp::{
    GpAx1, GpAx2, GpAx3, GpCirc, GpCone, GpCylinder, GpDir, GpElips, GpPln, GpPnt, GpSphere, GpTrsf,
    GpTorus, GpVec,
};
use occt_geom::{Curve, GeomCircle, GeomCylinder, GeomEllipse, GeomLine, GeomPlane, GeomSphere, GeomTorus, Surface};
use occt_math::spline_surface::{BSplineSurface, eval_bspline_surface, interpolate_grid};

use crate::brep_surface::{classify_surface, sphere_center, SurfaceKind};

/// A single 3D intersection curve between two surfaces, sampled into points
/// and per-surface parameters.
#[derive(Clone)]
pub struct IntersectionCurve {
    /// The 3D curve of intersection.
    pub curve: Arc<dyn Curve>,
    /// Sample points on the curve (curve evaluated on a uniform parameter grid).
    pub points: Vec<GpPnt>,
    /// The `(u, v)` parameters of each sample point on surface A.
    pub on_a: Vec<(f64, f64)>,
    /// The `(u, v)` parameters of each sample point on surface B.
    pub on_b: Vec<(f64, f64)>,
}

impl std::fmt::Debug for IntersectionCurve {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("IntersectionCurve")
            .field("points", &self.points)
            .field("on_a", &self.on_a)
            .field("on_b", &self.on_b)
            .finish()
    }
}

/// Result of intersecting two surfaces.
#[derive(Clone)]
pub enum SurfaceIntersection {
    /// One or more distinct intersection curves.
    Curves(Vec<IntersectionCurve>),
    /// The surfaces are geometrically coincident over a region.
    Coincident,
    /// No intersection found.
    None,
}

impl std::fmt::Debug for SurfaceIntersection {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SurfaceIntersection::Curves(c) => f.debug_tuple("Curves").field(c).finish(),
            SurfaceIntersection::Coincident => f.write_str("Coincident"),
            SurfaceIntersection::None => f.write_str("None"),
        }
    }
}

// ---------------------------------------------------------------------------
// Parameter helpers
// ---------------------------------------------------------------------------

/// Finite sampling bounds for a surface; unbounded ranges clamp to ±1
/// (matching `brep_surface::sample_bounds`).
pub fn sample_bounds(s: &dyn Surface) -> (f64, f64, f64, f64) {
    let (u0, u1) = s.u_range();
    let (v0, v1) = s.v_range();
    let clamp = |a: f64, b: f64| if a.is_finite() && b.is_finite() && b > a { (a, b) } else { (-1.0, 1.0) };
    let (u0, u1) = clamp(u0, u1);
    let (v0, v1) = clamp(v0, v1);
    (u0, u1, v0, v1)
}

/// Newton refinement of `(u, v)` minimizing `|s(u, v) − p|`.
///
/// Uses central finite differences of `d0` for the surface partials (many
/// ported analytic surfaces only implement `d0` and return zero `d1` vectors).
/// Returns `(u, v, s(u, v))`.
pub fn refine_point_on_surface(s: &dyn Surface, p: GpPnt, u0: f64, v0: f64, iters: usize) -> (f64, f64, GpPnt) {
    let (umin, umax) = s.u_range();
    let (vmin, vmax) = s.v_range();
    let clamp_u = |x: f64| if umin.is_finite() && umax.is_finite() { x.clamp(umin, umax) } else { x };
    let clamp_v = |x: f64| if vmin.is_finite() && vmax.is_finite() { x.clamp(vmin, vmax) } else { x };
    let (mut u, mut v) = (clamp_u(u0), clamp_v(v0));
    let h = 1e-6;
    for _ in 0..iters {
        let p0 = s.d0(u, v);
        let r = GpVec::from_pnts(&p0, &p);
        let pu = GpVec::from_pnts(&p0, &s.d0(u + h, v)).divided(h);
        let pv = GpVec::from_pnts(&p0, &s.d0(u, v + h)).divided(h);
        let (g11, g12, g22) = (pu.dot(&pu), pu.dot(&pv), pv.dot(&pv));
        let (b1, b2) = (pu.dot(&r), pv.dot(&r));
        let det = g11 * g22 - g12 * g12;
        if det.abs() < 1e-30 {
            break;
        }
        let du = (b1 * g22 - b2 * g12) / det;
        let dv = (b2 * g11 - b1 * g12) / det;
        u = clamp_u(u + du);
        v = clamp_v(v + dv);
        if du * du + dv * dv < 1e-24 {
            break;
        }
    }
    (u, v, s.d0(u, v))
}

/// Coarse grid search + refinement for the `(u, v)` parameters of the surface
/// point nearest `p`.
pub fn project_params(s: &dyn Surface, p: &GpPnt) -> (f64, f64) {
    let (u0, u1, v0, v1) = sample_bounds(s);
    let (nu, nv) = (24usize, 24usize);
    let mut bu = u0;
    let mut bv = v0;
    let mut bd = f64::INFINITY;
    for i in 0..=nu {
        for j in 0..=nv {
            let u = u0 + (u1 - u0) * i as f64 / nu as f64;
            let v = v0 + (v1 - v0) * j as f64 / nv as f64;
            let d = p.square_distance(&s.d0(u, v));
            if d < bd {
                bd = d;
                bu = u;
                bv = v;
            }
        }
    }
    let (u, v, _) = refine_point_on_surface(s, *p, bu, bv, 8);
    (u, v)
}

/// Approximate minimum distance from `p` to the surface `s`.
pub fn distance_to_surface(p: &GpPnt, s: &dyn Surface) -> f64 {
    let (u, v) = project_params(s, p);
    s.d0(u, v).distance(p)
}

/// The surface params of `p` when it lies within `tol` of the surface.
pub fn point_on_surface(s: &dyn Surface, p: &GpPnt, tol: f64) -> Option<(f64, f64)> {
    let (u, v) = project_params(s, p);
    if s.d0(u, v).distance(p) <= tol {
        Some((u, v))
    } else {
        None
    }
}

/// Extract a `GpPln` from a plane-like surface by sampling its geometry.
pub fn plane_from_surface(s: &dyn Surface) -> Option<GpPln> {
    if !crate::face_face::is_plane_like(s) {
        return None;
    }
    let (p, n) = crate::face_face::plane_geometry(s);
    let d = GpDir::from_vec(&n).ok()?;
    let z = GpDir::new(0.0, 0.0, 1.0).ok()?;
    let x = if d.is_normal(&z) { z } else { GpDir::new(1.0, 0.0, 0.0).ok()? };
    Some(GpPln::new(GpAx3::new(p, d, &x).ok()?))
}

// ---------------------------------------------------------------------------
// Curve construction helpers
// ---------------------------------------------------------------------------

/// A `GeomCircle` curve with the given center, plane normal and radius.
fn circle_curve(center: GpPnt, normal: GpDir, radius: f64) -> Result<Arc<dyn Curve>, String> {
    let z = GpDir::new(0.0, 0.0, 1.0).map_err(|_| "z dir")?;
    let x_dir = if normal.is_normal(&z) { z } else { GpDir::new(1.0, 0.0, 0.0).map_err(|_| "x dir")? };
    let ax2 = GpAx2::new(center, normal, x_dir).map_err(|e| e.to_string())?;
    Ok(Arc::new(GeomCircle::new(GpCirc::new(ax2, radius))))
}

/// An in-plane unit direction perpendicular to `normal`, for building frames.
fn perpendicular_dir(normal: &GpDir) -> Result<GpDir, String> {
    let z = GpDir::new(0.0, 0.0, 1.0).map_err(|_| "z dir")?;
    Ok(if normal.is_normal(&z) { z } else { GpDir::new(1.0, 0.0, 0.0).map_err(|_| "x dir")? })
}

/// Build an `IntersectionCurve` by sampling `curve` on a `n`-point grid and
/// projecting every sample onto both surfaces.
fn sample_curve_on(curve: Arc<dyn Curve>, a: &dyn Surface, b: &dyn Surface, n: usize) -> IntersectionCurve {
    let (f0, f1) = (curve.first_parameter(), curve.last_parameter());
    let n = n.max(2);
    let mut points = Vec::with_capacity(n);
    let mut on_a = Vec::with_capacity(n);
    let mut on_b = Vec::with_capacity(n);
    for i in 0..n {
        let t = if (f1 - f0).is_finite() {
            f0 + (f1 - f0) * i as f64 / (n - 1) as f64
        } else {
            -5.0 + 10.0 * i as f64 / (n - 1) as f64
        };
        let p = curve.d0(t);
        points.push(p);
        on_a.push(project_params(a, &p));
        on_b.push(project_params(b, &p));
    }
    IntersectionCurve { curve, points, on_a, on_b }
}

// ---------------------------------------------------------------------------
// B-spline surface support
// ---------------------------------------------------------------------------

/// A thin wrapper over `occt_math`'s tensor-product B-spline surface, exposing
/// it as an `occt_geom::Surface` on the unit square `[0, 1]²`.
/// Source: `Geom_BSplineSurface`.
#[derive(Debug, Clone)]
pub struct BsplinePatch {
    /// The underlying B-spline surface (poles, knots, degrees).
    pub bs: BSplineSurface,
}

impl BsplinePatch {
    /// Wrap an existing `BSplineSurface`.
    pub fn new(bs: BSplineSurface) -> Self {
        Self { bs }
    }

    /// Evaluate the patch at `(u, v)`.
    pub fn d0(&self, u: f64, v: f64) -> GpPnt {
        let a = eval_bspline_surface(&self.bs, u, v);
        GpPnt::new(a[0], a[1], a[2])
    }

    /// Convert into a boxed `Surface` handle.
    pub fn to_surface(self) -> Arc<dyn Surface> {
        Arc::new(BsplineSurfaceWrapper { patch: self })
    }
}

/// `occt_geom::Surface` implementation backed by a [`BsplinePatch`].
#[derive(Debug, Clone)]
pub struct BsplineSurfaceWrapper {
    pub patch: BsplinePatch,
}

impl BsplineSurfaceWrapper {
    pub fn new(bs: BSplineSurface) -> Self {
        Self { patch: BsplinePatch::new(bs) }
    }

    /// Evaluate the wrapped surface at `(u, v)`.
    pub fn d0(&self, u: f64, v: f64) -> GpPnt {
        self.patch.d0(u, v)
    }
}

impl Surface for BsplineSurfaceWrapper {
    fn d0(&self, u: f64, v: f64) -> GpPnt {
        let a = eval_bspline_surface(&self.patch.bs, u, v);
        GpPnt::new(a[0], a[1], a[2])
    }
    fn d1(&self, u: f64, v: f64) -> (GpPnt, GpVec, GpVec) {
        let p0 = self.patch.d0(u, v);
        let h = 1e-6;
        let pu = GpVec::from_pnts(&p0, &self.patch.d0(u + h, v)).divided(h);
        let pv = GpVec::from_pnts(&p0, &self.patch.d0(u, v + h)).divided(h);
        (p0, pu, pv)
    }
    fn u_range(&self) -> (f64, f64) { (0.0, 1.0) }
    fn v_range(&self) -> (f64, f64) { (0.0, 1.0) }
    fn continuity(&self) -> u8 {
        (self.patch.bs.deg_u.min(self.patch.bs.deg_v).saturating_sub(1)).min(3) as u8
    }
    fn transform(&mut self, t: &GpTrsf) {
        for row in &mut self.patch.bs.poles {
            for pole in row {
                let p = GpPnt::new(pole[0], pole[1], pole[2]).transformed(t);
                pole[0] = p.x();
                pole[1] = p.y();
                pole[2] = p.z();
            }
        }
    }
    fn clone_dyn(&self) -> Box<dyn Surface> { Box::new(self.clone()) }
}

/// Fit a B-spline patch through a `points[i][j]` grid (`i` = u rows, `j` = v
/// columns) with degrees `deg_u × deg_v`. Uses clamped uniform knots and
/// two-pass cubic collocation interpolation (from `occt_math`), so the surface
/// passes through every grid point at `u = i/(nu−1)`, `v = j/(nv−1)`.
pub fn fit_bspline_grid(points: &[Vec<GpPnt>], deg_u: usize, deg_v: usize) -> Result<BsplinePatch, String> {
    let nu = points.len();
    let nv = points.first().map(|r| r.len()).unwrap_or(0);
    if nu == 0 || nv == 0 {
        return Err("fit_bspline_grid: empty grid".into());
    }
    let rows: Vec<Vec<[f64; 3]>> = points
        .iter()
        .map(|row| row.iter().map(|p| [p.x(), p.y(), p.z()]).collect())
        .collect();
    let refs: Vec<&[[f64; 3]]> = rows.iter().map(|r| r.as_slice()).collect();
    let bs = interpolate_grid(&refs, nu, nv, deg_u, deg_v)?;
    Ok(BsplinePatch::new(bs))
}

/// Build a `Surface` from a grid of points, suitable as a general (B-spline)
/// face surface. See [`fit_bspline_grid`].
pub fn make_bspline_surface_from_grid(points: &[Vec<GpPnt>], deg_u: usize, deg_v: usize) -> Result<Arc<dyn Surface>, String> {
    Ok(fit_bspline_grid(points, deg_u, deg_v)?.to_surface())
}

// ---------------------------------------------------------------------------
// Analytic intersections
// ---------------------------------------------------------------------------

/// Plane ∩ sphere — a circle of intersection (when the plane cuts the sphere).
///
/// The plane is `n·x = d` with unit normal `n` and `d = n·location`. Distance
/// from the sphere center to the plane is `n·c − d`; when its magnitude is at
/// most `r` the intersection is a circle centered at `c − n·dist` with radius
/// `sqrt(r² − dist²)`, lying in the plane.
pub fn intersect_plane_sphere(pln: &GpPln, center: GpPnt, r: f64) -> Option<IntersectionCurve> {
    let n = GpVec::from_xyz(pln.axis().direction().xyz()).normalized();
    let d0 = n.dot(&GpVec::from_pnts(&GpPnt::zero(), &pln.location()));
    let dist = n.dot(&GpVec::from_pnts(&GpPnt::zero(), &center)) - d0;
    if dist.abs() > r + 1e-12 {
        return None;
    }
    let r_circle = (r * r - dist * dist).max(0.0).sqrt();
    let c = center.translated_vec(&n.multiplied_scalar(-dist));
    let normal = GpDir::from_vec(&n).ok()?;
    let curve = circle_curve(c, normal, r_circle).ok()?;
    let a: Arc<dyn Surface> = Arc::new(GeomPlane::new(pln.clone()));
    let ax3 = GpAx3::new(center, GpDir::new(0.0, 0.0, 1.0).ok()?, &GpDir::new(1.0, 0.0, 0.0).ok()?).ok()?;
    let b: Arc<dyn Surface> = Arc::new(GeomSphere::new(GpSphere::new(ax3, r).ok()?));
    Some(sample_curve_on(curve, a.as_ref(), b.as_ref(), 48))
}

/// Sphere ∩ sphere — a circle lying in the plane perpendicular to the line of
/// centers (analytic closed form).
///
/// Centers `c1, c2`, radii `r1, r2`. When `|r1 − r2| < d < r1 + r2` the
/// intersection is a circle centered on the axis at distance
/// `a = (r1² − r2² + d²) / (2d)` from `c1` with radius `sqrt(r1² − a²)`.
pub fn intersect_sphere_sphere(c1: GpPnt, r1: f64, c2: GpPnt, r2: f64) -> Option<IntersectionCurve> {
    let d = c1.distance(&c2);
    if d <= 1e-15 {
        return None; // concentric — coincident or disjoint, no single circle
    }
    // Externally/internal tangency collapses the circle to a single point
    // (radius → 0); treat within tolerance as a degenerate circle.
    if d > r1 + r2 + 1e-9 || d < (r1 - r2).abs() - 1e-9 {
        return None;
    }
    let u = GpVec::from_pnts(&c1, &c2).divided(d);
    let a = (r1 * r1 - r2 * r2 + d * d) / (2.0 * d);
    let r_circle = (r1 * r1 - a * a).max(0.0).sqrt();
    let c = c1.translated_vec(&u.multiplied_scalar(a));
    let normal = GpDir::from_vec(&u).ok()?;
    let curve = circle_curve(c, normal, r_circle).ok()?;
    // Build the two sphere surfaces for p-curve sampling.
    let z = GpDir::new(0.0, 0.0, 1.0).ok()?;
    let x = GpDir::new(1.0, 0.0, 0.0).ok()?;
    let ax3a = GpAx3::new(c1, z, &x).ok()?;
    let ax3b = GpAx3::new(c2, z, &x).ok()?;
    let a_surf: Arc<dyn Surface> = Arc::new(GeomSphere::new(GpSphere::new(ax3a, r1).ok()?));
    let b_surf: Arc<dyn Surface> = Arc::new(GeomSphere::new(GpSphere::new(ax3b, r2).ok()?));
    Some(sample_curve_on(curve, a_surf.as_ref(), b_surf.as_ref(), 48))
}

/// Plane ∩ cylinder — a circle (plane ⊥ axis), a pair of generatrix lines
/// (plane ∥ axis), or an ellipse (general orientation).
///
/// `ax` is the cylinder axis (a point on it + direction), `radius` the cylinder
/// radius. The result is `None` when the plane does not meet the cylinder.
pub fn intersect_plane_cylinder(pln: &GpPln, ax: &GpAx1, radius: f64) -> Option<IntersectionCurve> {
    let n = GpVec::from_xyz(pln.axis().direction().xyz()).normalized();
    let az = GpVec::from_xyz(ax.direction().xyz()).normalized();
    let n_par = n.dot(&az);
    let n_perp = n.subtracted(&az.multiplied_scalar(n_par));
    let n_perp_len = n_perp.magnitude();
    let axis_pt = *ax.location();
    // Plane equation: n·p = d0.
    let d0 = n.dot(&GpVec::from_pnts(&GpPnt::zero(), &pln.location()));
    // Offset from axis_pt to the plane along n:  n·(axis_pt + λ·n) = d0.
    let lambda = d0 - n.dot(&GpVec::from_pnts(&GpPnt::zero(), &axis_pt));

    // Case 1: plane perpendicular to the axis → circle.
    if n_perp_len <= 1e-9 {
        let t = lambda / n_par; // distance along the axis from axis_pt to the plane
        let c = axis_pt.translated_vec(&az.multiplied_scalar(t));
        let normal = GpDir::from_vec(&n).ok()?;
        let curve = circle_curve(c, normal, radius).ok()?;
        let a: Arc<dyn Surface> = Arc::new(GeomPlane::new(pln.clone()));
        let z = GpDir::from_vec(&az).ok()?;
        let x = perpendicular_dir(&z).ok()?;
        let ax3 = GpAx3::new(axis_pt, z, &x).ok()?;
        let b: Arc<dyn Surface> = Arc::new(GeomCylinder::new(GpCylinder::new(ax3, radius).ok()?));
        return Some(sample_curve_on(curve, a.as_ref(), b.as_ref(), 48));
    }

    // Case 2: plane parallel to the axis → two (or one) generatrix lines.
    if n_par.abs() <= 1e-9 {
        let dist = lambda.abs();
        if dist > radius + 1e-12 {
            return None;
        }
        let e = az.crossed(&n_perp).normalized();
        let half = (radius * radius - dist * dist).max(0.0).sqrt();
        let base = axis_pt.translated_vec(&n_perp.multiplied_scalar(-lambda));
        let a: Arc<dyn Surface> = Arc::new(GeomPlane::new(pln.clone()));
        let z = GpDir::from_vec(&az).ok()?;
        let x = perpendicular_dir(&z).ok()?;
        let ax3 = GpAx3::new(axis_pt, z, &x).ok()?;
        let b: Arc<dyn Surface> = Arc::new(GeomCylinder::new(GpCylinder::new(ax3, radius).ok()?));
        let dir = GpDir::from_vec(&az).ok()?;
        let p0 = base.translated_vec(&e.multiplied_scalar(half));
        let p1 = base.translated_vec(&e.multiplied_scalar(-half));
        let mut curves: Vec<Arc<dyn Curve>> = Vec::new();
        for p in [p0, p1] {
            let lin = occt_core::gp::GpLin::from_pnt_dir(p, dir);
            curves.push(Arc::new(GeomLine::new(lin)) as Arc<dyn Curve>);
        }
        // Return the first line with sampled points from both lines.
        let mut ic = sample_curve_on(curves[0].clone(), a.as_ref(), b.as_ref(), 16);
        let mut all_points = ic.points.clone();
        let mut on_a = ic.on_a.clone();
        let mut on_b = ic.on_b.clone();
        for c in curves.iter().skip(1) {
            for t in [-2.0, -1.0, 0.0, 1.0, 2.0] {
                let p = c.d0(t);
                all_points.push(p);
                on_a.push(project_params(a.as_ref(), &p));
                on_b.push(project_params(b.as_ref(), &p));
            }
        }
        ic.points = all_points;
        ic.on_a = on_a;
        ic.on_b = on_b;
        return Some(ic);
    }

    // Case 3: general — an ellipse.
    let t = lambda / n_par;
    let c = axis_pt.translated_vec(&az.multiplied_scalar(t));
    let e1 = n.crossed(&az).normalized();
    let e2 = n.crossed(&e1).normalized();
    let semi_a = radius / n_par.abs(); // along e2
    let semi_b = radius;               // along e1
    let (major, x_dir) = if semi_a >= semi_b { (semi_a, e2) } else { (semi_b, e1) };
    let minor = semi_a.min(semi_b);
    let normal = GpDir::from_vec(&n).ok()?;
    let xd = GpDir::from_vec(&x_dir).ok()?;
    let ax2 = GpAx2::new(c, normal, xd).ok()?;
    let elips = GpElips::new(ax2, major, minor);
    let curve: Arc<dyn Curve> = Arc::new(GeomEllipse::new(elips));
    let a: Arc<dyn Surface> = Arc::new(GeomPlane::new(pln.clone()));
    let z = GpDir::from_vec(&az).ok()?;
    let x = perpendicular_dir(&z).ok()?;
    let ax3 = GpAx3::new(axis_pt, z, &x).ok()?;
    let b: Arc<dyn Surface> = Arc::new(GeomCylinder::new(GpCylinder::new(ax3, radius).ok()?));
    Some(sample_curve_on(curve, a.as_ref(), b.as_ref(), 48))
}

/// Plane ∩ cone — a circle when the plane is perpendicular to the axis;
/// otherwise falls back to the general tracer (parabola/hyperbola/ellipse).
///
/// `apex` is the cone apex, `ax` its axis, `semi_angle` the half-angle between
/// the axis and a generator.
pub fn intersect_plane_cone(pln: &GpPln, apex: GpPnt, ax: &GpAx1, semi_angle: f64) -> Option<IntersectionCurve> {
    let n = GpVec::from_xyz(pln.axis().direction().xyz()).normalized();
    let az = GpVec::from_xyz(ax.direction().xyz()).normalized();
    let n_par = n.dot(&az);
    let n_perp = n.subtracted(&az.multiplied_scalar(n_par));
    if n_perp.magnitude() > 1e-9 {
        return None; // general conic — handled by trace_surface_curve
    }
    // Plane perpendicular to the axis: circle at t0 from the apex.
    let d0 = n.dot(&GpVec::from_pnts(&GpPnt::zero(), &pln.location()));
    let t0 = d0 - n.dot(&GpVec::from_pnts(&GpPnt::zero(), &apex));
    let t_dist = t0 / n_par;
    let r_circle = (t_dist * semi_angle.tan()).abs();
    let c = apex.translated_vec(&az.multiplied_scalar(t_dist));
    let normal = GpDir::from_vec(&n).ok()?;
    let curve = circle_curve(c, normal, r_circle).ok()?;
    let a: Arc<dyn Surface> = Arc::new(GeomPlane::new(pln.clone()));
    let z = GpDir::from_vec(&az).ok()?;
    let x = perpendicular_dir(&z).ok()?;
    let ax3 = GpAx3::new(apex, z, &x).ok()?;
    let cone = GpCone::new(ax3, r_circle.max(1e-9), semi_angle).ok()?;
    let b: Arc<dyn Surface> = Arc::new(occt_geom::GeomCone::new(cone));
    Some(sample_curve_on(curve, a.as_ref(), b.as_ref(), 48))
}

/// Plane ∩ torus — up to two circles when the plane contains the axis or is
/// perpendicular through the center; otherwise the general tracer.
pub fn intersect_plane_torus(pln: &GpPln, center: GpPnt, ax: &GpAx1, major: f64, minor: f64) -> Option<IntersectionCurve> {
    let n = GpVec::from_xyz(pln.axis().direction().xyz()).normalized();
    let az = GpVec::from_xyz(ax.direction().xyz()).normalized();
    let n_par = n.dot(&az);
    let n_perp = n.subtracted(&az.multiplied_scalar(n_par));
    let normal = GpDir::from_vec(&n).ok()?;

    // Case A: the plane contains the torus axis → two circles of radius `minor`
    // centered at ±major along the in-plane perpendicular to the axis.
    if n_par.abs() <= 1e-9 {
        let d0 = n.dot(&GpVec::from_pnts(&GpPnt::zero(), &pln.location()));
        let dc = n.dot(&GpVec::from_pnts(&GpPnt::zero(), &center));
        if (d0 - dc).abs() > 1e-9 {
            return None;
        }
        let e1 = n.crossed(&az).normalized();
        let c1 = center.translated_vec(&e1.multiplied_scalar(major));
        let c2 = center.translated_vec(&e1.multiplied_scalar(-major));
        let curve1 = circle_curve(c1, normal, minor).ok()?;
        let curve2 = circle_curve(c2, normal, minor).ok()?;
        let a: Arc<dyn Surface> = Arc::new(GeomPlane::new(pln.clone()));
        let z = GpDir::from_vec(&az).ok()?;
        let x = perpendicular_dir(&z).ok()?;
        let ax3 = GpAx3::new(center, z, &x).ok()?;
        let torus = GpTorus::new(ax3, major, minor).ok()?;
        let b: Arc<dyn Surface> = Arc::new(GeomTorus::new(torus));
        let ic1 = sample_curve_on(curve1, a.as_ref(), b.as_ref(), 32);
        let ic2 = sample_curve_on(curve2, a.as_ref(), b.as_ref(), 32);
        let mut merged = ic1.clone();
        merged.points.extend(ic2.points.iter().copied());
        merged.on_a.extend(ic2.on_a.iter().copied());
        merged.on_b.extend(ic2.on_b.iter().copied());
        return Some(merged);
    }

    // Case B: plane perpendicular to the axis and passing through the center →
    // two circles of radius major ± minor.
    if n_perp.magnitude() <= 1e-9 {
        let d0 = n.dot(&GpVec::from_pnts(&GpPnt::zero(), &pln.location()));
        let dc = n.dot(&GpVec::from_pnts(&GpPnt::zero(), &center));
        if (d0 - dc).abs() > 1e-9 {
            return None;
        }
        let mut curves = Vec::new();
        for r in [major + minor, (major - minor).abs()] {
            if r > 1e-12 {
                curves.push(circle_curve(center, normal, r).ok()?);
            }
        }
        let a: Arc<dyn Surface> = Arc::new(GeomPlane::new(pln.clone()));
        let z = GpDir::from_vec(&az).ok()?;
        let x = perpendicular_dir(&z).ok()?;
        let ax3 = GpAx3::new(center, z, &x).ok()?;
        let torus = GpTorus::new(ax3, major, minor).ok()?;
        let b: Arc<dyn Surface> = Arc::new(GeomTorus::new(torus));
        let ic1 = sample_curve_on(curves[0].clone(), a.as_ref(), b.as_ref(), 32);
        let mut merged = ic1.clone();
        if let Some(c2) = curves.get(1) {
            let ic2 = sample_curve_on(c2.clone(), a.as_ref(), b.as_ref(), 32);
            merged.points.extend(ic2.points.iter().copied());
            merged.on_a.extend(ic2.on_a.iter().copied());
            merged.on_b.extend(ic2.on_b.iter().copied());
        }
        return Some(merged);
    }

    None // general torus section — handled by trace_surface_curve
}

// ---------------------------------------------------------------------------
// General tracer
// ---------------------------------------------------------------------------

/// A piecewise-linear `Curve` through sampled points (parameter `[0, 1]`).
/// Used to carry grid-traced intersections.
#[derive(Debug, Clone)]
pub struct PolylineCurve {
    pub pts: Vec<GpPnt>,
}

impl Curve for PolylineCurve {
    fn d0(&self, u: f64) -> GpPnt {
        let n = self.pts.len();
        if n == 0 {
            return GpPnt::zero();
        }
        if n == 1 {
            return self.pts[0];
        }
        let t = u.clamp(0.0, 1.0) * (n - 1) as f64;
        let i = (t.floor() as usize).min(n - 2);
        let f = t - i as f64;
        let a = self.pts[i];
        let b = self.pts[i + 1];
        GpPnt::new(a.x() + f * (b.x() - a.x()), a.y() + f * (b.y() - a.y()), a.z() + f * (b.z() - a.z()))
    }
    fn d1(&self, u: f64) -> (GpPnt, GpVec) {
        let p0 = self.d0(u);
        let h = 1e-6;
        (p0, GpVec::from_pnts(&p0, &self.d0(u + h)).divided(h))
    }
    fn d2(&self, u: f64) -> (GpPnt, GpVec, GpVec) {
        let (p, d1) = self.d1(u);
        (p, d1, GpVec::zero())
    }
    fn first_parameter(&self) -> f64 { 0.0 }
    fn last_parameter(&self) -> f64 { 1.0 }
    fn continuity(&self) -> u8 { 0 }
    fn transform(&mut self, t: &GpTrsf) {
        for p in &mut self.pts {
            *p = p.transformed(t);
        }
    }
    fn reverse(&mut self) { self.pts.reverse(); }
    fn clone_dyn(&self) -> Box<dyn Curve> { Box::new(self.clone()) }
}

/// Crossing of a segment `a → b` (with signed field values `va, vb`) at the
/// zero level set.
fn edge_crossing(pa: &GpPnt, va: f64, pb: &GpPnt, vb: f64) -> Option<GpPnt> {
    if va * vb > 0.0 {
        return None;
    }
    if va.abs() < 1e-12 && vb.abs() < 1e-12 {
        return None;
    }
    if va.abs() < 1e-12 {
        return Some(*pa);
    }
    if vb.abs() < 1e-12 {
        return Some(*pb);
    }
    let t = va / (va - vb);
    Some(GpPnt::new(
        pa.x() + t * (pb.x() - pa.x()),
        pa.y() + t * (pb.y() - pa.y()),
        pa.z() + t * (pb.z() - pa.z()),
    ))
}

/// Chain a set of unordered 3D segments into polylines, merging endpoints that
/// are within `tol`, and return the flattened point stream.
fn chain_segments(segs: Vec<(GpPnt, GpPnt)>, tol: f64) -> Vec<GpPnt> {
    let mut remaining = segs;
    let mut out: Vec<GpPnt> = Vec::new();
    while !remaining.is_empty() {
        let (a, b) = remaining.remove(0);
        let mut poly = vec![a, b];
        let mut changed = true;
        while changed {
            changed = false;
            for i in (0..remaining.len()).rev() {
                let (c, d) = remaining[i];
                let first = poly[0];
                let last = *poly.last().unwrap();
                if c.distance(&last) <= tol {
                    poly.push(d);
                    remaining.swap_remove(i);
                    changed = true;
                } else if d.distance(&last) <= tol {
                    poly.push(c);
                    remaining.swap_remove(i);
                    changed = true;
                } else if c.distance(&first) <= tol {
                    poly.insert(0, d);
                    remaining.swap_remove(i);
                    changed = true;
                } else if d.distance(&first) <= tol {
                    poly.insert(0, c);
                    remaining.swap_remove(i);
                    changed = true;
                }
            }
        }
        out.extend(poly);
    }
    out
}

/// General fallback intersection: grid-march `a`'s parameter space and keep the
/// contour where the distance from `a(u, v)` to surface `b` equals `tol`.
///
/// Returns the traced 3D points (concatenated polylines). The level-set
/// contour is extracted with a marching-squares pass over the `(u, v)` grid.
pub fn trace_surface_curve(a: &dyn Surface, b: &dyn Surface, tol: f64) -> Vec<GpPnt> {
    trace_surface_curve_n(a, b, tol, 48, 48)
}

/// `trace_surface_curve` with an explicit marching-squares grid size.
fn trace_surface_curve_n(a: &dyn Surface, b: &dyn Surface, tol: f64, nu: usize, nv: usize) -> Vec<GpPnt> {
    let (u0, u1, v0, v1) = sample_bounds(a);
    let (nu, nv) = (nu.max(4), nv.max(4));
    let mut field = vec![vec![0.0f64; nv + 1]; nu + 1];
    for i in 0..=nu {
        for j in 0..=nv {
            let u = u0 + (u1 - u0) * i as f64 / nu as f64;
            let v = v0 + (v1 - v0) * j as f64 / nv as f64;
            let p = a.d0(u, v);
            field[i][j] = distance_to_surface(&p, b) - tol;
        }
    }

    // Marching squares: emit a segment per cell that the zero contour crosses.
    let mut segments: Vec<(GpPnt, GpPnt)> = Vec::new();
    for i in 0..nu {
        for j in 0..nv {
            // Corners CCW: 0=(i,j), 1=(i+1,j), 2=(i+1,j+1), 3=(i,j+1).
            let p = [
                a.d0(u0 + (u1 - u0) * i as f64 / nu as f64, v0 + (v1 - v0) * j as f64 / nv as f64),
                a.d0(u0 + (u1 - u0) * (i + 1) as f64 / nu as f64, v0 + (v1 - v0) * j as f64 / nv as f64),
                a.d0(u0 + (u1 - u0) * (i + 1) as f64 / nu as f64, v0 + (v1 - v0) * (j + 1) as f64 / nv as f64),
                a.d0(u0 + (u1 - u0) * i as f64 / nu as f64, v0 + (v1 - v0) * (j + 1) as f64 / nv as f64),
            ];
            let f = [field[i][j], field[i + 1][j], field[i + 1][j + 1], field[i][j + 1]];
            let e01 = edge_crossing(&p[0], f[0], &p[1], f[1]);
            let e12 = edge_crossing(&p[1], f[1], &p[2], f[2]);
            let e23 = edge_crossing(&p[2], f[2], &p[3], f[3]);
            let e30 = edge_crossing(&p[3], f[3], &p[0], f[0]);
            let mut pts = Vec::new();
            for e in [e01, e12, e23, e30] {
                if let Some(q) = e {
                    pts.push(q);
                }
            }
            if pts.len() == 2 {
                segments.push((pts[0], pts[1]));
            } else if pts.len() >= 4 {
                // Saddle — connect opposite crossings.
                segments.push((pts[0], pts[2]));
                segments.push((pts[1], pts[3]));
            }
        }
    }
    chain_segments(segments, tol * 2.0)
}

/// Chain an unordered set of 3D points into ordered polylines by greedy
/// nearest-neighbour walking from both ends of each chain. Points separated by
/// more than a (diagonal-relative) link distance start a new chain; junctions
/// are left as separate chains rather than forcing a single path through them.
///
/// The link tolerance is `max(tol, 2% of the bounding-box diagonal)`, so chains
/// whose overlapping endpoints lie within `tol` merge into a single polyline
/// while well-separated disjoint curves stay separate.
pub fn chain_intersection_points(pts: &[GpPnt], tol: f64) -> Vec<Vec<GpPnt>> {
    if pts.is_empty() {
        return Vec::new();
    }
    let mut lo = pts[0];
    let mut hi = pts[0];
    for p in pts {
        lo = GpPnt::new(lo.x().min(p.x()), lo.y().min(p.y()), lo.z().min(p.z()));
        hi = GpPnt::new(hi.x().max(p.x()), hi.y().max(p.y()), hi.z().max(p.z()));
    }
    let link_tol = tol.max(lo.distance(&hi) * 0.02);

    let mut used = vec![false; pts.len()];
    let mut chains: Vec<Vec<GpPnt>> = Vec::new();
    loop {
        let Some(start) = (0..pts.len()).find(|&i| !used[i]) else { break };
        used[start] = true;
        let mut chain: Vec<GpPnt> = vec![pts[start]];
        loop {
            let tail = *chain.last().unwrap();
            let head = chain[0];
            let mut best: Option<(f64, usize, bool)> = None; // (dist, idx, at_tail)
            for (i, p) in pts.iter().enumerate() {
                if used[i] {
                    continue;
                }
                let dt = p.distance(&tail);
                let dh = p.distance(&head);
                if best.map_or(true, |(bd, _, _)| dt < bd || dh < bd) {
                    if dt <= dh {
                        best = Some((dt, i, true));
                    } else {
                        best = Some((dh, i, false));
                    }
                }
            }
            match best {
                Some((d, i, at_tail)) if d <= link_tol => {
                    used[i] = true;
                    if at_tail {
                        chain.push(pts[i]);
                    } else {
                        chain.insert(0, pts[i]);
                    }
                }
                _ => break,
            }
        }
        chains.push(chain);
    }
    chains
}

/// Convert a chained intersection polyline into an [`IntersectionCurve`].
///
/// Long polylines (> 8 points) are first simplified with the
/// Ramer–Douglas–Peucker algorithm at `tol`; the resulting vertices become a
/// [`PolylineCurve`], and each vertex is projected onto both surfaces to fill
/// `on_a` / `on_b`.
pub fn polyline_to_curve(poly: &[GpPnt], a: &dyn Surface, b: &dyn Surface, tol: f64) -> Option<IntersectionCurve> {
    if poly.len() < 2 {
        return None;
    }
    let pts: Vec<GpPnt> = if poly.len() > 8 {
        let keep = rdp_simplify(poly, tol);
        keep.iter().map(|&i| poly[i]).collect()
    } else {
        poly.to_vec()
    };
    if pts.len() < 2 {
        return None;
    }
    let curve: Arc<dyn Curve> = Arc::new(PolylineCurve { pts: pts.clone() });
    let on_a: Vec<(f64, f64)> = pts.iter().map(|p| project_params(a, p)).collect();
    let on_b: Vec<(f64, f64)> = pts.iter().map(|p| project_params(b, p)).collect();
    Some(IntersectionCurve { curve, points: pts, on_a, on_b })
}

/// Trace the intersection curves of two general surfaces and return them as
/// ordered polylines (each `Vec<GpPnt>` is one chained curve). Exposed for the
/// curved boolean (`bop_curved`) which needs per-face intersection polylines.
pub fn intersection_curve_points(a: &dyn Surface, b: &dyn Surface, tol: f64, samples: usize) -> Vec<Vec<GpPnt>> {
    let grid = samples.clamp(24, 96);
    let pts = trace_surface_curve_n(a, b, tol, grid, grid);
    if pts.is_empty() {
        return Vec::new();
    }
    chain_intersection_points(&pts, tol)
}

/// Intersect two general (possibly non-analytic) surfaces using the grid
/// tracer, chaining the traced points into polylines and wrapping each in an
/// [`IntersectionCurve`].
pub fn intersect_general_surfaces(a: &dyn Surface, b: &dyn Surface, tol: f64) -> SurfaceIntersection {
    let pts = trace_surface_curve(a, b, tol);
    if pts.is_empty() {
        return SurfaceIntersection::None;
    }
    let chains = chain_intersection_points(&pts, tol);
    let mut curves = Vec::new();
    for chain in &chains {
        if let Some(ic) = polyline_to_curve(chain, a, b, tol) {
            curves.push(ic);
        }
    }
    if curves.is_empty() {
        SurfaceIntersection::None
    } else {
        SurfaceIntersection::Curves(curves)
    }
}

/// Quick overlap test: `true` when the grid tracer finds any intersection
/// point between the two surfaces.
pub fn surfaces_intersect_general(a: &dyn Surface, b: &dyn Surface, tol: f64) -> bool {
    !trace_surface_curve(a, b, tol).is_empty()
}

/// Sample points from the intersection curves of `a` and `b`, retaining only
/// those within `tol` of both surfaces (verification helper).
pub fn points_on_both(a: &dyn Surface, b: &dyn Surface, tol: f64, nsamples: usize) -> Vec<GpPnt> {
    let mut out = Vec::new();
    match surface_surface_intersection(a, b, tol) {
        SurfaceIntersection::Curves(curves) => {
            for ic in &curves {
                let (f0, f1) = (ic.curve.first_parameter(), ic.curve.last_parameter());
                let n = nsamples.max(4);
                for i in 0..n {
                    let t = if (f1 - f0).is_finite() {
                        f0 + (f1 - f0) * i as f64 / (n - 1) as f64
                    } else {
                        -2.0 + 4.0 * i as f64 / (n - 1) as f64
                    };
                    let p = ic.curve.d0(t);
                    if distance_to_surface(&p, a) <= tol && distance_to_surface(&p, b) <= tol {
                        out.push(p);
                    }
                }
            }
        }
        _ => {}
    }
    out
}

// ---------------------------------------------------------------------------
// Dispatcher
// ---------------------------------------------------------------------------

/// Intersect two surfaces, dispatching on their analytic types and falling back
/// to the grid tracer for general pairs.
pub fn surface_surface_intersection(a: &dyn Surface, b: &dyn Surface, tol: f64) -> SurfaceIntersection {
    let kind_a = classify_surface(a);
    let kind_b = classify_surface(b);

    match (kind_a, kind_b) {
        (SurfaceKind::Plane, SurfaceKind::Plane) => {
            let pa = match plane_from_surface(a) {
                Some(p) => p,
                None => return SurfaceIntersection::None,
            };
            let pb = match plane_from_surface(b) {
                Some(p) => p,
                None => return SurfaceIntersection::None,
            };
            let n1 = GpVec::from_xyz(pa.axis().direction().xyz());
            let n2 = GpVec::from_xyz(pb.axis().direction().xyz());
            if n1.cross_magnitude(&n2) <= tol {
                let dd = n1
                    .normalized()
                    .dot(&GpVec::from_pnts(&pa.location(), &pb.location()))
                    .abs();
                return if dd <= tol { SurfaceIntersection::Coincident } else { SurfaceIntersection::None };
            }
            return match crate::face_face::plane_plane_intersection(&pa, &pb) {
                Some((origin, dir)) => {
                    let lin = occt_core::gp::GpLin::from_pnt_dir(origin, GpDir::from_vec(&dir).unwrap_or_default());
                    let curve: Arc<dyn Curve> = Arc::new(GeomLine::new(lin));
                    SurfaceIntersection::Curves(vec![sample_curve_on(curve, a, b, 32)])
                }
                None => SurfaceIntersection::None,
            };
        }
        (SurfaceKind::Plane, SurfaceKind::Sphere) => {
            let pa = match plane_from_surface(a) {
                Some(p) => p,
                None => return SurfaceIntersection::None,
            };
            let (c, r) = match sphere_params(b) {
                Some(x) => x,
                None => return SurfaceIntersection::None,
            };
            return match intersect_plane_sphere(&pa, c, r) {
                Some(ic) => SurfaceIntersection::Curves(vec![ic]),
                None => SurfaceIntersection::None,
            };
        }
        (SurfaceKind::Sphere, SurfaceKind::Plane) => {
            let pb = match plane_from_surface(b) {
                Some(p) => p,
                None => return SurfaceIntersection::None,
            };
            let (c, r) = match sphere_params(a) {
                Some(x) => x,
                None => return SurfaceIntersection::None,
            };
            return match intersect_plane_sphere(&pb, c, r) {
                Some(ic) => SurfaceIntersection::Curves(vec![ic]),
                None => SurfaceIntersection::None,
            };
        }
        (SurfaceKind::Sphere, SurfaceKind::Sphere) => {
            let (c1, r1) = match sphere_params(a) {
                Some(x) => x,
                None => return SurfaceIntersection::None,
            };
            let (c2, r2) = match sphere_params(b) {
                Some(x) => x,
                None => return SurfaceIntersection::None,
            };
            return match intersect_sphere_sphere(c1, r1, c2, r2) {
                Some(ic) => SurfaceIntersection::Curves(vec![ic]),
                None => SurfaceIntersection::None,
            };
        }
        _ => {}
    }

    // General fallback: grid tracer chained into per-curve polylines.
    intersect_general_surfaces(a, b, tol)
}

/// Extract sphere center + radius from a surface by sampling.
pub fn sphere_params(s: &dyn Surface) -> Option<(GpPnt, f64)> {
    let c = sphere_center(s)?;
    let (u0, _, v0, _) = sample_bounds(s);
    let r = s.d0(u0, v0).distance(&c);
    Some((c, r))
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;
    use occt_core::gp::GpAx1;

    const TOL: f64 = 1e-6;

    fn unit_sphere(center: GpPnt) -> Arc<dyn Surface> {
        let ax3 = GpAx3::new(center, GpDir::new(0.0, 0.0, 1.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap()).unwrap();
        Arc::new(GeomSphere::new(GpSphere::new(ax3, 1.0).unwrap()))
    }

    fn plane_z(z: f64) -> GpPln {
        GpPln::new(GpAx3::new(GpPnt::new(0.0, 0.0, z), GpDir::new(0.0, 0.0, 1.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap()).unwrap())
    }

    #[test]
    fn plane_sphere_circle_radius_and_on_surface() {
        // Plane z = 0.5 cutting the unit sphere at origin: circle radius sqrt(1-0.25) ≈ 0.866.
        let pln = plane_z(0.5);
        let ic = intersect_plane_sphere(&pln, GpPnt::zero(), 1.0).expect("circle");
        assert_eq!(ic.points.len(), 48);
        for p in &ic.points {
            let r = p.distance(&GpPnt::zero());
            assert!((r - 1.0).abs() < 1e-9, "on sphere: r {r}");
            assert!((p.z() - 0.5).abs() < 1e-9, "in plane: z {}", p.z());
        }
        // Radius of the circle ≈ sqrt(3/4).
        let rad = ic.points[0].distance(&GpPnt::new(0.0, 0.0, 0.5));
        assert!((rad - (0.75f64.sqrt())).abs() < 1e-9, "circle radius {rad}");
        // Every sampled point is within tol of both surfaces.
        let s = unit_sphere(GpPnt::zero());
        let a: Arc<dyn Surface> = Arc::new(GeomPlane::new(pln.clone()));
        for p in &ic.points {
            assert!(distance_to_surface(p, a.as_ref()) < TOL);
            assert!(distance_to_surface(p, s.as_ref()) < TOL);
        }
    }

    #[test]
    fn plane_sphere_no_intersection() {
        let pln = plane_z(2.0); // above the unit sphere
        assert!(intersect_plane_sphere(&pln, GpPnt::zero(), 1.0).is_none());
    }

    #[test]
    fn plane_sphere_tangent() {
        let pln = plane_z(1.0);
        let ic = intersect_plane_sphere(&pln, GpPnt::zero(), 1.0).expect("tangent circle");
        // Tangent: radius 0, all points collapse to the tangent point.
        let p = ic.points[0];
        assert!(p.distance(&GpPnt::new(0.0, 0.0, 1.0)) < 1e-9, "tangent point {p:?}");
    }

    #[test]
    fn sphere_sphere_circle_on_both() {
        // Two unit spheres, centers 1.5 apart.
        let ic = intersect_sphere_sphere(GpPnt::zero(), 1.0, GpPnt::new(1.5, 0.0, 0.0), 1.0).expect("circle");
        let s1 = unit_sphere(GpPnt::zero());
        let s2 = unit_sphere(GpPnt::new(1.5, 0.0, 0.0));
        // Intersection plane at x = 0.75, radius sqrt(1 - 0.75²).
        let expected_r = (1.0 - 0.75f64 * 0.75).sqrt();
        for p in &ic.points {
            assert!((p.x() - 0.75).abs() < 1e-9, "plane x {}", p.x());
            let r = GpPnt::new(0.0, p.y(), p.z()).distance(&GpPnt::zero());
            assert!((r - expected_r).abs() < 1e-9, "radius {r} (expected {expected_r})");
            assert!(distance_to_surface(p, s1.as_ref()) < TOL);
            assert!(distance_to_surface(p, s2.as_ref()) < TOL);
        }
    }

    #[test]
    fn sphere_sphere_tangent_or_none() {
        // Tangent spheres (d = r1 + r2) → degenerate single point (radius 0).
        let ic = intersect_sphere_sphere(GpPnt::zero(), 1.0, GpPnt::new(2.0, 0.0, 0.0), 1.0).expect("tangent circle");
        let p = ic.points[0];
        assert!(p.distance(&GpPnt::new(1.0, 0.0, 0.0)) < 1e-6, "tangent point {p:?}");

        // Disjoint spheres → None.
        assert!(intersect_sphere_sphere(GpPnt::zero(), 1.0, GpPnt::new(5.0, 0.0, 0.0), 1.0).is_none());
    }

    #[test]
    fn plane_cylinder_circle() {
        // Plane z = 1 perpendicular to the Z-axis cylinder radius 2 → circle radius 2 at z=1.
        let pln = plane_z(1.0);
        let ax = GpAx1::new(GpPnt::zero(), GpDir::new(0.0, 0.0, 1.0).unwrap());
        let ic = intersect_plane_cylinder(&pln, &ax, 2.0).expect("circle");
        for p in &ic.points {
            assert!((p.z() - 1.0).abs() < 1e-9);
            let r = GpPnt::new(p.x(), p.y(), 0.0).distance(&GpPnt::zero());
            assert!((r - 2.0).abs() < 1e-9, "radius {r}");
        }
    }

    #[test]
    fn plane_cylinder_parallel_lines() {
        // Plane x = 0.5 parallel to the Z-axis cylinder radius 1 → two lines at y = ±sqrt(0.75).
        let pln = GpPln::new(GpAx3::new(
            GpPnt::new(0.5, 0.0, 0.0),
            GpDir::new(1.0, 0.0, 0.0).unwrap(),
            &GpDir::new(0.0, 1.0, 0.0).unwrap(),
        ).unwrap());
        let ax = GpAx1::new(GpPnt::zero(), GpDir::new(0.0, 0.0, 1.0).unwrap());
        let ic = intersect_plane_cylinder(&pln, &ax, 1.0).expect("lines");
        assert!(ic.points.len() >= 2);
        let mut ys: Vec<f64> = ic.points.iter().map(|p| p.y().abs()).collect();
        ys.sort_by(f64::total_cmp);
        let target = (0.75f64).sqrt();
        assert!((ys[0] - target).abs() < 1e-6, "line y {target}, got {}", ys[0]);
    }

    #[test]
    fn sphere_cylinder_via_trace() {
        // Unit sphere at origin ∩ cylinder radius 0.5 axis Z → circle at |z| = sqrt(0.75).
        let s = unit_sphere(GpPnt::zero());
        let ax3 = GpAx3::new(GpPnt::zero(), GpDir::new(0.0, 0.0, 1.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap()).unwrap();
        let c: Arc<dyn Surface> = Arc::new(GeomCylinder::new(GpCylinder::new(ax3, 0.5).unwrap()));
        let pts = trace_surface_curve(s.as_ref(), c.as_ref(), 0.03);
        assert!(!pts.is_empty(), "traced points");
        for p in pts.iter().take(20) {
            assert!(distance_to_surface(p, s.as_ref()) < 0.05);
            assert!(distance_to_surface(p, c.as_ref()) < 0.05);
            let r = GpPnt::new(p.x(), p.y(), 0.0).distance(&GpPnt::zero());
            assert!((r - 0.5).abs() < 0.05, "on cylinder r {r}");
        }
    }

    #[test]
    fn refine_point_on_sphere_converges() {
        let s = unit_sphere(GpPnt::zero());
        let target = GpPnt::new(0.0, 0.0, 1.0);
        let (u, v, p) = refine_point_on_surface(s.as_ref(), target, 0.5, 0.5, 20);
        assert!(p.distance(&target) < 1e-6, "refined {p:?} vs {target:?}");
        assert!((v - std::f64::consts::FRAC_PI_2).abs() < 1e-4, "v {v}");
        let _ = u;
    }

    #[test]
    fn points_on_both_filters() {
        let s1 = unit_sphere(GpPnt::zero());
        let s2 = unit_sphere(GpPnt::new(1.5, 0.0, 0.0));
        let pts = points_on_both(s1.as_ref(), s2.as_ref(), 1e-6, 32);
        assert!(!pts.is_empty());
        for p in &pts {
            assert!(distance_to_surface(p, s1.as_ref()) < 1e-6);
            assert!(distance_to_surface(p, s2.as_ref()) < 1e-6);
        }
    }

    #[test]
    fn surface_surface_intersection_sphere_sphere() {
        let s1 = unit_sphere(GpPnt::zero());
        let s2 = unit_sphere(GpPnt::new(1.5, 0.0, 0.0));
        match surface_surface_intersection(s1.as_ref(), s2.as_ref(), 1e-6) {
            SurfaceIntersection::Curves(curves) => {
                assert_eq!(curves.len(), 1);
                assert!(!curves[0].points.is_empty());
            }
            other => panic!("expected curves, got {other:?}"),
        }
    }

    #[test]
    fn surface_surface_intersection_disjoint_none() {
        let s1 = unit_sphere(GpPnt::zero());
        let s2 = unit_sphere(GpPnt::new(5.0, 0.0, 0.0));
        assert!(matches!(
            surface_surface_intersection(s1.as_ref(), s2.as_ref(), 1e-6),
            SurfaceIntersection::None
        ));
    }

    #[test]
    fn plane_cone_circle() {
        // Cone apex at origin, semi-angle 30°, plane z = 2 → circle radius 2·tan30 ≈ 1.155.
        let pln = plane_z(2.0);
        let ax = GpAx1::new(GpPnt::zero(), GpDir::new(0.0, 0.0, 1.0).unwrap());
        let semi = 30f64.to_radians();
        let ic = intersect_plane_cone(&pln, GpPnt::zero(), &ax, semi).expect("circle");
        let expected_r = 2.0 * semi.tan();
        for p in &ic.points {
            assert!((p.z() - 2.0).abs() < 1e-9);
            let r = GpPnt::new(p.x(), p.y(), 0.0).distance(&GpPnt::zero());
            assert!((r - expected_r).abs() < 1e-9, "radius {r} (expected {expected_r})");
        }
    }

    #[test]
    fn plane_torus_two_circles() {
        // Torus major 3 minor 1; plane y=0 (contains axis) → two circles radius 1 at x=±3.
        let pln = GpPln::new(GpAx3::new(
            GpPnt::zero(),
            GpDir::new(0.0, 1.0, 0.0).unwrap(),
            &GpDir::new(1.0, 0.0, 0.0).unwrap(),
        ).unwrap());
        let ax = GpAx1::new(GpPnt::zero(), GpDir::new(0.0, 0.0, 1.0).unwrap());
        let ic = intersect_plane_torus(&pln, GpPnt::zero(), &ax, 3.0, 1.0).expect("two circles");
        assert!(!ic.points.is_empty());
        let torus_ax3 = GpAx3::new(GpPnt::zero(), GpDir::new(0.0, 0.0, 1.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap()).unwrap();
        let t: Arc<dyn Surface> = Arc::new(GeomTorus::new(GpTorus::new(torus_ax3, 3.0, 1.0).unwrap()));
        for p in &ic.points {
            assert!((p.y()).abs() < 1e-9, "in plane y {}", p.y());
            assert!(distance_to_surface(p, t.as_ref()) < 1e-6);
            // Each point lies on one of the two circles (radius 1) centered at x = ±3.
            let d1 = p.distance(&GpPnt::new(3.0, 0.0, 0.0));
            let d2 = p.distance(&GpPnt::new(-3.0, 0.0, 0.0));
            assert!(
                (d1 - 1.0).abs() < 1e-9 || (d2 - 1.0).abs() < 1e-9,
                "on a minor circle: d1 {d1}, d2 {d2}"
            );
        }
    }

    #[test]
    fn trace_and_polyline_curve() {
        let c = PolylineCurve {
            pts: vec![GpPnt::zero(), GpPnt::new(1.0, 0.0, 0.0), GpPnt::new(2.0, 0.0, 0.0)],
        };
        assert!(c.d0(0.0).distance(&GpPnt::zero()) < 1e-12);
        assert!(c.d0(0.5).distance(&GpPnt::new(1.0, 0.0, 0.0)) < 1e-12);
        assert!(c.d0(1.0).distance(&GpPnt::new(2.0, 0.0, 0.0)) < 1e-12);
    }

    #[test]
    fn plane_plane_intersection_line() {
        let xy = plane_z(0.0);
        let yz = GpPln::new(GpAx3::new(
            GpPnt::zero(),
            GpDir::new(1.0, 0.0, 0.0).unwrap(),
            &GpDir::new(0.0, 1.0, 0.0).unwrap(),
        ).unwrap());
        let s1: Arc<dyn Surface> = Arc::new(GeomPlane::new(xy.clone()));
        let s2: Arc<dyn Surface> = Arc::new(GeomPlane::new(yz.clone()));
        match surface_surface_intersection(s1.as_ref(), s2.as_ref(), 1e-9) {
            SurfaceIntersection::Curves(curves) => {
                assert_eq!(curves.len(), 1);
                for p in &curves[0].points {
                    assert!(p.z().abs() < 1e-9);
                    assert!(p.x().abs() < 1e-9);
                }
            }
            other => panic!("expected line, got {other:?}"),
        }
    }

    // -- general (non-analytic) surface intersection extensions --

    fn sphere_r(center: GpPnt, r: f64) -> Arc<dyn Surface> {
        let ax3 = GpAx3::new(center, GpDir::new(0.0, 0.0, 1.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap()).unwrap();
        Arc::new(GeomSphere::new(GpSphere::new(ax3, r).unwrap()))
    }

    fn cylinder_z(r: f64) -> Arc<dyn Surface> {
        let ax3 = GpAx3::new(GpPnt::zero(), GpDir::new(0.0, 0.0, 1.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap()).unwrap();
        Arc::new(GeomCylinder::new(GpCylinder::new(ax3, r).unwrap()))
    }

    /// A degree-1 B-spline patch over a tilted plane `z = u + 0.5·v` on a
    /// 4×4 grid spanning `[0, 1]²`.
    fn tilted_plane_patch() -> Arc<dyn Surface> {
        let mut grid: Vec<Vec<GpPnt>> = Vec::new();
        for i in 0..=3usize {
            let mut row = Vec::new();
            for j in 0..=3usize {
                let u = i as f64 / 3.0;
                let v = j as f64 / 3.0;
                row.push(GpPnt::new(u, v, u + 0.5 * v));
            }
            grid.push(row);
        }
        make_bspline_surface_from_grid(&grid, 1, 1).expect("patch fit")
    }

    #[test]
    fn general_plane_plane_line() {
        let s1: Arc<dyn Surface> = Arc::new(GeomPlane::new(plane_z(0.0)));
        let yz = GpPln::new(GpAx3::new(
            GpPnt::zero(),
            GpDir::new(1.0, 0.0, 0.0).unwrap(),
            &GpDir::new(0.0, 1.0, 0.0).unwrap(),
        ).unwrap());
        let s2: Arc<dyn Surface> = Arc::new(GeomPlane::new(yz));
        // The analytic dispatcher still returns a line for two planes.
        match surface_surface_intersection(s1.as_ref(), s2.as_ref(), 1e-6) {
            SurfaceIntersection::Curves(curves) => {
                assert_eq!(curves.len(), 1);
                for p in &curves[0].points {
                    assert!(p.z().abs() < 1e-6);
                    assert!(p.x().abs() < 1e-6);
                }
            }
            other => panic!("expected line, got {other:?}"),
        }
    }

    #[test]
    fn general_sphere_sphere_via_trace() {
        // Two radius-2 spheres with centers 3 apart intersect in a circle;
        // the general tracer (bypassing the analytic sphere-sphere path) finds it.
        let a = sphere_r(GpPnt::zero(), 2.0);
        let b = sphere_r(GpPnt::new(3.0, 0.0, 0.0), 2.0);
        let tol = 0.05;
        match intersect_general_surfaces(a.as_ref(), b.as_ref(), tol) {
            SurfaceIntersection::Curves(curves) => {
                assert!(!curves.is_empty(), "expected at least one traced curve");
                for ic in &curves {
                    assert!(ic.points.len() >= 2, "chained points");
                    for p in &ic.points {
                        assert!(distance_to_surface(p, a.as_ref()) < 4.0 * tol, "on sphere a: {p:?}");
                        assert!(distance_to_surface(p, b.as_ref()) < 4.0 * tol, "on sphere b: {p:?}");
                    }
                }
            }
            other => panic!("expected curves, got {other:?}"),
        }
    }

    #[test]
    fn general_bspline_plane_intersection() {
        let patch = tilted_plane_patch();
        let pln: Arc<dyn Surface> = Arc::new(GeomPlane::new(plane_z(0.5)));
        let tol = 0.02;
        match intersect_general_surfaces(patch.as_ref(), pln.as_ref(), tol) {
            SurfaceIntersection::Curves(curves) => {
                assert!(!curves.is_empty(), "expected a traced line");
                let mut found = 0;
                for ic in &curves {
                    for p in &ic.points {
                        assert!(distance_to_surface(p, patch.as_ref()) < 6.0 * tol);
                        assert!(distance_to_surface(p, pln.as_ref()) < 6.0 * tol);
                        found += 1;
                    }
                }
                assert!(found >= 2, "at least two points on the intersection line");
            }
            other => panic!("expected curves, got {other:?}"),
        }
    }

    #[test]
    fn general_cylinder_plane_parallel_lines() {
        let cyl = cylinder_z(1.0);
        let pln_x = GpPln::new(GpAx3::new(
            GpPnt::new(0.5, 0.0, 0.0),
            GpDir::new(1.0, 0.0, 0.0).unwrap(),
            &GpDir::new(0.0, 1.0, 0.0).unwrap(),
        ).unwrap());
        let pln: Arc<dyn Surface> = Arc::new(GeomPlane::new(pln_x));
        let tol = 0.05;
        match intersect_general_surfaces(cyl.as_ref(), pln.as_ref(), tol) {
            SurfaceIntersection::Curves(curves) => {
                assert!(!curves.is_empty(), "expected generatrix lines");
                for ic in &curves {
                    for p in &ic.points {
                        assert!(distance_to_surface(p, cyl.as_ref()) < 4.0 * tol, "on cylinder: {p:?}");
                        assert!(distance_to_surface(p, pln.as_ref()) < 4.0 * tol, "on plane: {p:?}");
                    }
                }
            }
            other => panic!("expected curves, got {other:?}"),
        }
    }

    #[test]
    fn chain_intersection_merges() {
        // Two polylines whose endpoint regions overlap within tol merge to one chain.
        let chain_a: Vec<GpPnt> = (0..=10).map(|i| GpPnt::new(i as f64 * 0.1, 0.0, 0.0)).collect();
        let chain_b: Vec<GpPnt> = (0..=10).map(|i| GpPnt::new(0.95 + i as f64 * 0.1, 0.0, 0.0)).collect();
        let mut pts = chain_a.clone();
        pts.extend(chain_b);
        let chains = chain_intersection_points(&pts, 0.12);
        assert_eq!(chains.len(), 1, "overlapping polylines merge into one chain");
        assert_eq!(chains[0].len(), 22, "all 22 points chained in order");
        // Well-separated curves stay separate.
        let far: Vec<GpPnt> = (0..=5).map(|i| GpPnt::new(5.0 + i as f64 * 0.1, 0.0, 0.0)).collect();
        let mut pts2 = pts;
        pts2.extend(far);
        let chains2 = chain_intersection_points(&pts2, 0.12);
        assert_eq!(chains2.len(), 2, "a separated curve starts a new chain");
    }

    #[test]
    fn polyline_to_curve_on_params() {
        let s = unit_sphere(GpPnt::zero());
        let pln: Arc<dyn Surface> = Arc::new(GeomPlane::new(plane_z(0.5)));
        // Circle of intersection: radius sqrt(0.75) at z = 0.5.
        let r = (0.75f64).sqrt();
        let poly: Vec<GpPnt> = (0..16)
            .map(|i| {
                let a = 2.0 * std::f64::consts::PI * i as f64 / 16.0;
                GpPnt::new(r * a.cos(), r * a.sin(), 0.5)
            })
            .collect();
        let ic = polyline_to_curve(&poly, s.as_ref(), pln.as_ref(), 1e-4).expect("curve");
        assert_eq!(ic.points.len(), poly.len());
        assert_eq!(ic.on_a.len(), poly.len());
        assert_eq!(ic.on_b.len(), poly.len());
        for (p, (ua, va)) in ic.points.iter().zip(&ic.on_a) {
            let q = s.d0(*ua, *va);
            assert!(q.distance(p) < 1e-5, "on_a reconstructs {p:?} as {q:?}");
        }
        for (p, (ub, vb)) in ic.points.iter().zip(&ic.on_b) {
            let q = pln.d0(*ub, *vb);
            assert!(q.distance(p) < 1e-5, "on_b reconstructs {p:?} as {q:?}");
        }
    }

    #[test]
    fn general_nonintersecting_empty() {
        let s = unit_sphere(GpPnt::zero());
        let far_pln: Arc<dyn Surface> = Arc::new(GeomPlane::new(plane_z(10.0)));
        assert!(matches!(
            intersect_general_surfaces(s.as_ref(), far_pln.as_ref(), 0.05),
            SurfaceIntersection::None
        ));
        assert!(!surfaces_intersect_general(s.as_ref(), far_pln.as_ref(), 0.05));
        assert!(surfaces_intersect_general(s.as_ref(), unit_sphere(GpPnt::new(1.0, 0.0, 0.0)).as_ref(), 0.05));
    }

    #[test]
    fn fit_bspline_grid_corners() {
        let mut grid: Vec<Vec<GpPnt>> = Vec::new();
        for i in 0..=3usize {
            let mut row = Vec::new();
            for j in 0..=3usize {
                let u = i as f64 / 3.0;
                let v = j as f64 / 3.0;
                row.push(GpPnt::new(u * 2.0 - 1.0, v * 2.0 - 1.0, u * u + v));
            }
            grid.push(row);
        }
        let patch = fit_bspline_grid(&grid, 1, 1).expect("fit");
        let corners = [
            (0.0, 0.0, &grid[0][0]),
            (1.0, 0.0, &grid[3][0]),
            (0.0, 1.0, &grid[0][3]),
            (1.0, 1.0, &grid[3][3]),
        ];
        for (u, v, expected) in corners {
            let got = patch.d0(u, v);
            assert!(got.distance(expected) < 1e-9, "corner ({u},{v}): got {got:?} expected {expected:?}");
        }
    }
}
