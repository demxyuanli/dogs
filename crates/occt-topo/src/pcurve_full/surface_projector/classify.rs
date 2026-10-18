//! Surface-type classification plus the sampled analytic frames (`projector_cylinder_params`/`cone_params`/`torus_params`) and the small vector helpers they share.

use super::*;



// ---------------------------------------------------------------------------
// Surface classification
// ---------------------------------------------------------------------------

/// Classify a surface by its concrete geometry type.
///
/// This mirrors `Adaptor3d_Surface::GetType()`: `ProjLib_ProjectedCurve.cxx`
/// (line 375) selects its analytic projection from `GeomAbs_SurfaceType`
/// (`case GeomAbs_Plane:` ... `case GeomAbs_Torus:`), and
/// `ShapeFix_Edge::FixAddPCurve` (`cxx:487`) tests `surf->IsKind(Geom_Plane)`.
/// Both are *type* tests, so a surface that is only geometrically a cylinder /
/// cone / torus (a `Geom_BSplineSurface` patch) must stay `Other` and be
/// projected by the general path; treating it as analytic would feed the
/// sampled analytic frame a parameterization it does not have.
pub fn classify_surface_kind(s: &dyn Surface) -> SurfaceKind {
    use crate::geom_bnd_lib_surface3d::GeomAbsSurfaceType as TypeKind;
    match crate::geom_bnd_lib_surface3d::surface_kind(s) {
        TypeKind::Plane => SurfaceKind::Plane,
        TypeKind::Sphere => SurfaceKind::Sphere,
        TypeKind::Cylinder => SurfaceKind::Cylinder,
        TypeKind::Cone => SurfaceKind::Cone,
        TypeKind::Torus => SurfaceKind::Torus,
        _ => SurfaceKind::Other,
    }
}

/// Cylinder parameters `(center, unit axis, radius)` from sampled invariants:
/// the axis is along `∂S/∂v` and two opposite points of a `u`-ring give the
/// radius (their midpoint lies on the axis).
pub(crate) fn projector_cylinder_params(s: &dyn Surface) -> Option<(GpPnt, GpVec, f64)> {
    let p0 = s.d0(0.0, 0.0);
    let p1 = s.d0(0.0, 1.0);
    let axis = GpVec::from_pnts(&p0, &p1);
    let m = axis.magnitude();
    if m < 1e-9 {
        return None;
    }
    let ax = axis.divided(m);
    let q0 = s.d0(0.0, 0.0);
    let q1 = s.d0(PI, 0.0);
    let r = q0.distance(&q1) / 2.0;
    if r <= 1e-9 {
        return None;
    }
    let center = sampled_mid(&q0, &q1);
    let (nu, nv) = (8, 6);
    for i in 0..=nu {
        for j in 0..=nv {
            let u = 2.0 * PI * i as f64 / nu as f64;
            let v = (j as f64 - nv as f64 / 2.0) / 2.0;
            let p = s.d0(u, v);
            if (dist_to_axis(&p, &center, &ax) - r).abs() > 1e-3 * r.max(1.0) {
                return None;
            }
        }
    }
    Some((center, ax, r))
}

/// Cone parameters `(apex, unit axis, semi-angle)`: two fixed-`u` generatrices
/// intersect at the apex; every sampled point makes the same angle with the axis.
pub(crate) fn cone_params(s: &dyn Surface) -> Option<(GpPnt, GpVec, f64)> {
    let g0a = s.d0(0.0, 0.0);
    let g0b = s.d0(0.0, 1.0);
    let d0 = GpVec::from_pnts(&g0a, &g0b);
    if d0.magnitude() < 1e-9 {
        return None;
    }
    let g1a = s.d0(FRAC_PI_2, 0.0);
    let g1b = s.d0(FRAC_PI_2, 1.0);
    let d1 = GpVec::from_pnts(&g1a, &g1b);
    if d1.magnitude() < 1e-9 {
        return None;
    }
    let apex = closest_point_lines(g0a, d0, g1a, d1)?;
    // The cone axis runs through the apex and the midpoint of two opposite
    // points on a u-ring. Sample the ring at v=1 (not v=0 — every u collapses
    // to the apex there, making the axis a zero vector).
    let m0 = sampled_mid(&s.d0(0.0, 1.0), &s.d0(PI, 1.0));
    let axis = GpVec::from_pnts(&apex, &m0);
    let am = axis.magnitude();
    if am < 1e-9 {
        return None;
    }
    let ax = axis.divided(am);
    // Use the v=1 point for the apex angle: `g0a` (v=0) may coincide with the
    // apex itself (a cone whose parameter v=0 sits at the vertex), making the
    // vector from apex to g0a a zero vector and the angle undefined.
    let alpha = GpVec::from_pnts(&apex, &g0b).angle(&ax);
    if alpha < 1e-6 || alpha > FRAC_PI_2 - 1e-6 {
        return None;
    }
    let (nu, nv) = (6, 5);
    for i in 0..=nu {
        for j in 1..=nv {
            let u = 2.0 * PI * i as f64 / nu as f64;
            // Sample the valid cone half only (v > 0): the surface's v-range is
            // unbounded but d0 is only well-defined on the cone sheet away from
            // the apex; negative v can land on the opposite (invalid) sheet.
            let v = j as f64 / nv as f64;
            let p = s.d0(u, v);
            if p.distance(&apex) < 1e-9 {
                continue;
            }
            let ang = GpVec::from_pnts(&apex, &p).angle(&ax);
            if (ang - alpha).abs() > 1e-3 * alpha.max(0.01) {
                return None;
            }
        }
    }
    Some((apex, ax, alpha))
}

/// Torus parameters `(center, axis, major radius, minor radius)`.
pub(crate) fn torus_params(s: &dyn Surface) -> Option<(GpPnt, GpVec, f64, f64)> {
    let top = s.d0(0.0, FRAC_PI_2);
    let bot = s.d0(0.0, -FRAC_PI_2);
    let axv = GpVec::from_pnts(&bot, &top);
    let m = axv.magnitude();
    if m < 1e-9 {
        return None;
    }
    let ax = axv.divided(m);
    let e0 = s.d0(0.0, 0.0);
    let e1 = s.d0(PI, 0.0);
    let center = sampled_mid(&e0, &e1);
    let rel = GpVec::from_pnts(&center, &top);
    let rel_axis = rel.dot(&ax);
    let in_plane = rel.subtracted(&ax.multiplied_scalar(rel_axis));
    let big_r = in_plane.magnitude();
    let small_r = rel_axis.abs();
    if big_r <= 1e-9 || small_r <= 1e-9 {
        return None;
    }
    let (nu, nv) = (6, 6);
    for i in 0..=nu {
        for j in 0..=nv {
            let u = 2.0 * PI * i as f64 / nu as f64;
            let v = 2.0 * PI * j as f64 / nv as f64;
            let p = s.d0(u, v);
            let rr = GpVec::from_pnts(&center, &p);
            let ra = rr.dot(&ax);
            // Signed radial distance `R + r·cos v` — not the unsigned magnitude.
            // A spindle torus (R < r) folds through itself, so its inner equator
            // (v ≈ π) has a negative rho; the unsigned |R + r·cos v| breaks the
            // tube equation there. Recover the sign from the torus identity
            // |rr|² = rho² + ra² ⇒ rho = (|rr|² + R² − r²) / (2R), which is exact
            // for both ring (R > r) and spindle (R < r) tori (ElSLib convention).
            let rho = (rr.square_magnitude() + big_r.powi(2) - small_r.powi(2)) / (2.0 * big_r);
            let err = ((rho - big_r).powi(2) + ra.powi(2)).sqrt();
            if (err - small_r).abs() > 1e-3 * small_r.max(1.0) {
                return None;
            }
        }
    }
    Some((center, ax, big_r, small_r))
}

// ---------------------------------------------------------------------------
// Small vector helpers
// ---------------------------------------------------------------------------

pub(crate) fn sampled_mid(a: &GpPnt, b: &GpPnt) -> GpPnt {
    GpPnt::new(0.5 * (a.x() + b.x()), 0.5 * (a.y() + b.y()), 0.5 * (a.z() + b.z()))
}

pub(crate) fn dist_to_axis(p: &GpPnt, c: &GpPnt, ax: &GpVec) -> f64 {
    let rel = GpVec::from_pnts(c, p);
    rel.subtracted(&ax.multiplied_scalar(rel.dot(ax))).magnitude()
}

pub(crate) fn pnt_add_vec(p: GpPnt, v: &GpVec) -> GpPnt {
    GpPnt::from_xyz(&p.coord.added(&v.coord))
}

/// Closest point between two (possibly skew) lines; `None` when parallel or
/// too far apart.
pub(crate) fn closest_point_lines(p1: GpPnt, d1: GpVec, p2: GpPnt, d2: GpVec) -> Option<GpPnt> {
    let r = GpVec::from_pnts(&p1, &p2);
    let a = d1.dot(&d1);
    let b = d1.dot(&d2);
    let c = d2.dot(&d2);
    let e = d1.dot(&r);
    let f = d2.dot(&r);
    let denom = a * c - b * b;
    if denom.abs() < 1e-12 {
        return None;
    }
    let t = (e * c - b * f) / denom;
    let s = (b * e - a * f) / denom;
    let q1 = pnt_add_vec(p1, &d1.multiplied_scalar(t));
    let q2 = pnt_add_vec(p2, &d2.multiplied_scalar(s));
    if q1.distance(&q2) < 1e-6 * (1.0 + p1.distance(&p2).max(1.0)) {
        Some(q1)
    } else {
        None
    }
}

// ---------------------------------------------------------------------------
// Analytic surface projector
// ---------------------------------------------------------------------------

