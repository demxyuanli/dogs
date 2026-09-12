//! Closed-form quadric pair intersections (plane/sphere/cylinder/cone/torus).
//! Source: analytic helpers previously in `intpatch.rs`.

use std::sync::Arc;

use occt_core::gp::{
    GpAx1, GpAx2, GpAx3, GpCirc, GpCone, GpCylinder, GpDir, GpElips, GpPln, GpPnt, GpSphere, GpTorus,
    GpVec,
};
use occt_geom::{
    Curve, GeomCircle, GeomCylinder, GeomEllipse, GeomLine, GeomPlane, GeomSphere, GeomTorus,
    Surface,
};

use super::geom::project_params;
use super::IntersectionCurve;

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
pub(crate) fn sample_curve_on(
    curve: Arc<dyn Curve>,
    a: &dyn Surface,
    b: &dyn Surface,
    n: usize,
) -> IntersectionCurve {
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
