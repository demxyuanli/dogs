use super::prelude::*;
use super::*;

/// Surface `(u, v)` parameters of a point, exact for reconstructed quadrics
/// (so the reported parameters are precise and the UV-bounds validation is
/// meaningful), falling back to the sampling projector for other surfaces.
/// Periodic U/V directions are wrapped into the surface's natural range.
pub(super) fn surface_params(surface: &dyn Surface, geom: Option<&SurfaceGeom>, p: &GpPnt) -> (f64, f64) {
    let (u0, _, _, _) = sample_bounds(surface);
    if let Some(g) = geom {
        match g {
            // Plane and sphere frames are reconstructed from the surface's own
            // natural parameterization, so the recovered (u, v) are exact.
            SurfaceGeom::Plane { o, x, y, .. } => {
                let d = GpVec::from_pnts(o, p);
                return (d.dot(x), d.dot(y));
            }
            SurfaceGeom::Sphere { o, r, x, y, z } => {
                let d = GpVec::from_pnts(o, p);
                let u = wrap_periodic(d.dot(y).atan2(d.dot(x)), u0, 2.0 * std::f64::consts::PI);
                let s = snap_pole(d.dot(z) / r);
                return (u, s.asin());
            }
            // Cylinder/cone/torus frames use a reconstructed axis/reference that
            // can differ by a constant offset from the surface's natural
            // parameterization; the sampling projector keeps (u, v) consistent.
            _ => return surface_closest_params(surface, p, 24, 24),
        }
    }
    surface_closest_params(surface, p, 24, 24)
}

// ---------------------------------------------------------------------------
// Analytic solves
// ---------------------------------------------------------------------------

/// Result of an analytic curve-vs-quadric solve in the curve parameter.
pub(super) enum SolveResult {
    /// Discrete roots in the curve parameter.
    Roots(Vec<f64>),
    /// The curve lies on the surface (coincidence) — an `On` segment.
    OnSurface,
    /// The curve is disjoint (including parallel-and-offset).
    None,
}

/// Real roots of a·x² + b·x + c = 0 (stable), with the linear/parallel
/// degeneracy folded in.
pub(super) fn solve_quadratic(a: f64, b: f64, c: f64, surface_eps: f64) -> SolveResult {
    if a.abs() < 1e-300 {
        if b.abs() < 1e-300 {
            if c.abs() <= surface_eps {
                return SolveResult::OnSurface;
            }
            return SolveResult::None;
        }
        return SolveResult::Roots(vec![-c / b]);
    }
    let disc = b * b - 4.0 * a * c;
    if disc < 0.0 {
        return SolveResult::None;
    }
    if disc.abs() <= 1e-12 * a.abs().max(1.0) {
        return SolveResult::Roots(vec![-b / (2.0 * a)]);
    }
    let sq = disc.sqrt();
    let q = -0.5 * (b + b.signum() * sq);
    let r1 = q / a;
    let r2 = if q.abs() > 1e-300 { c / q } else { (-b - sq) / (2.0 * a) };
    SolveResult::Roots(vec![r1, r2])
}

pub(super) fn line_plane(loc: &GpPnt, v: &GpVec, o: &GpPnt, n: &GpVec) -> SolveResult {
    let denom = n.dot(v);
    if denom.abs() < ANG_TOL {
        // Parallel to the plane: intersects only if it lies in it.
        let d = GpVec::from_pnts(o, loc);
        if n.dot(&d).abs() <= ANALYTIC_TOL {
            SolveResult::OnSurface
        } else {
            SolveResult::None
        }
    } else {
        let t = n.dot(&GpVec::from_pnts(loc, o)) / denom;
        SolveResult::Roots(vec![t])
    }
}

pub(super) fn line_sphere(loc: &GpPnt, v: &GpVec, c: &GpPnt, r: f64) -> SolveResult {
    let d = GpVec::from_pnts(c, loc); // loc − c
    let a = v.square_magnitude();
    if a < 1e-30 {
        return SolveResult::None;
    }
    let b = 2.0 * d.dot(v);
    let cc = d.square_magnitude() - r * r;
    solve_quadratic(a, b, cc, 1e-9)
}

pub(super) fn line_cylinder(loc: &GpPnt, v: &GpVec, a: &GpPnt, z: &GpVec, r: f64) -> SolveResult {
    let d0 = GpVec::from_pnts(a, loc); // loc − a
    let d0z = d0.dot(z);
    let dz = v.dot(z);
    let a2 = v.square_magnitude() - dz * dz;
    let b = 2.0 * (d0.dot(v) - d0z * dz);
    let c = d0.square_magnitude() - d0z * d0z - r * r;
    solve_quadratic(a2, b, c, 1e-9)
}

pub(super) fn line_cone(loc: &GpPnt, v: &GpVec, apex: &GpPnt, z: &GpVec, cosa: f64) -> SolveResult {
    let d0 = GpVec::from_pnts(apex, loc); // loc − apex
    let d0z = d0.dot(z);
    let dz = v.dot(z);
    let v2 = v.square_magnitude();
    let a = dz * dz - v2 * cosa * cosa;
    let b = 2.0 * (d0z * dz - d0.dot(v) * cosa * cosa);
    let c = d0z * d0z - d0.square_magnitude() * cosa * cosa;
    solve_quadratic(a, b, c, 1e-9)
}

/// `a·cos(u) + b·sin(u) = c` roots in `[inf, sup]` via the trigonometric
/// reduction (`math_TrigonometricFunctionRoots`).
pub(super) fn trig_solve(a: f64, b: f64, c: f64, inf: f64, sup: f64, surface_eps: f64) -> SolveResult {
    if a.abs() < 1e-12 && b.abs() < 1e-12 {
        if c.abs() <= surface_eps {
            return SolveResult::OnSurface;
        }
        return SolveResult::None;
    }
    SolveResult::Roots(occt_math::trig_roots(a, b, -c, inf, sup))
}

// ---------------------------------------------------------------------------
// Main entry
// ---------------------------------------------------------------------------

/// Compute the intersection of `curve` with `surface`.
///
/// * `cu_range` — the curve parameter window to search (`(first, last)`); for
///   unbounded curves (lines) this must be finite.
/// * `uv_range` — the surface `(u0, v0, u1, v1)` window to search; unbounded
///   surface directions are clamped to it.
///
/// Returns discrete intersection points (sorted by curve parameter) and, for
/// the coincidence degeneracy (e.g. a line lying in a plane), `On` segments.
pub fn perform_curve_surface(
    curve: &dyn Curve,
    surface: &dyn Surface,
    cu_range: (f64, f64),
    uv_range: (f64, f64, f64, f64),
) -> Result<HInterResult, String> {
    let (t0, t1) = effective_curve_range(curve, cu_range)?;
    let uv = effective_uv_range(surface, uv_range)?;

    let curve_kind = classify_curve(curve);
    let surface_kind = classify_surface_kind(surface);
    let geom = build_surface_geom(surface, surface_kind);

    let mut points: Vec<IntersectionPoint> = Vec::new();
    let mut segments: Vec<IntersectionSegment> = Vec::new();

    // Analytic path for conic curves against supported quadrics.
    // Unbounded windows are valid here: `IntCurvesFace_Intersector.cxx:386-388`
    // calls `HICS.Perform(HLL, Hsurface)` with no polygon when the surface is
    // a plane/quadric. The sampling fallback needs a finite interval.
    let finite_window = t0.is_finite() && t1.is_finite();
    if curve_kind != CurveKind::Other {
        if let Some((mut pts, mut segs)) = geom
            .as_ref()
            .and_then(|g| perform_conic_surf(curve, curve_kind, surface, g, t0, t1, uv))
        {
            points.append(&mut pts);
            segments.append(&mut segs);
        } else if finite_window {
            let mut pts = general_curve_surface(curve, surface, geom.as_ref(), t0, t1, uv, GENERAL_TOL);
            points.append(&mut pts);
        }
    } else if finite_window {
        let mut pts = general_curve_surface(curve, surface, geom.as_ref(), t0, t1, uv, GENERAL_TOL);
        points.append(&mut pts);
    } else {
        return Err("perform_curve_surface: unbounded curve parameter range".into());
    }

    // Sort by curve parameter and merge near-duplicate (tangent) points.
    points.sort_by(|a, b| a.param.partial_cmp(&b.param).unwrap_or(Ordering::Equal));
    let mut deduped: Vec<IntersectionPoint> = Vec::with_capacity(points.len());
    for p in points {
        if let Some(last) = deduped.last_mut() {
            if (p.param - last.param).abs() < MERGE_TOL && last.pnt.distance(&p.pnt) <= 1e-6 {
                last.state = State::On;
                continue;
            }
        }
        deduped.push(p);
    }

    Ok(HInterResult { points: deduped, segments })
}

/// Clamp the curve parameter window to the curve's own range.
///
/// Unbounded curves (lines) may keep a non-finite search window: analytic
/// conic/quadric solves do not sample the curve
/// (`IntCurvesFace_Intersector.cxx:386-388`). The sampling fallback rejects
/// that window in `perform_curve_surface`.
pub(super) fn effective_curve_range(curve: &dyn Curve, cu_range: (f64, f64)) -> Result<(f64, f64), String> {
    let (ca, cb) = cu_range;
    let (lo, hi) = if ca <= cb { (ca, cb) } else { (cb, ca) };
    let fa = curve.first_parameter();
    let fb = curve.last_parameter();
    let a = if fa.is_finite() { fa.max(lo) } else { lo };
    let b = if fb.is_finite() { fb.min(hi) } else { hi };
    if a.is_finite() && b.is_finite() && b - a <= 1e-15 {
        return Err("perform_curve_surface: empty curve range".into());
    }
    Ok((a, b))
}

/// Clamp the surface `(u0, v0, u1, v1)` window to the surface's own ranges.
/// Unbounded directions take the caller's bounds directly.
pub(super) fn effective_uv_range(surface: &dyn Surface, uv: (f64, f64, f64, f64)) -> Result<(f64, f64, f64, f64), String> {
    let (u0, u1, v0, v1) = uv;
    let (su0, su1) = surface.u_range();
    let (sv0, sv1) = surface.v_range();
    let clamp_pair = |a: f64, b: f64, sa: f64, sb: f64| -> (f64, f64) {
        let (a, b) = if a <= b { (a, b) } else { (b, a) };
        let a = if sa.is_finite() { a.max(sa) } else { a };
        let b = if sb.is_finite() { b.min(sb) } else { b };
        (a, b)
    };
    let (u0, u1) = clamp_pair(u0, u1, su0, su1);
    let (v0, v1) = clamp_pair(v0, v1, sv0, sv1);
    if !u0.is_finite() || !u1.is_finite() || !v0.is_finite() || !v1.is_finite() {
        return Err("perform_curve_surface: unbounded surface UV range".into());
    }
    if u1 - u0 <= 1e-15 || v1 - v0 <= 1e-15 {
        return Err("perform_curve_surface: empty surface UV range".into());
    }
    Ok((u0, u1, v0, v1))
}

/// Analytic dispatch for conic curves against quadric surfaces. Returns
/// `None` when the combination has no analytic solve (caller falls back to the
/// general path).
pub(super) fn perform_conic_surf(
    curve: &dyn Curve,
    kind: CurveKind,
    surface: &dyn Surface,
    geom: &SurfaceGeom,
    t0: f64,
    t1: f64,
    uv: (f64, f64, f64, f64),
) -> Option<(Vec<IntersectionPoint>, Vec<IntersectionSegment>)> {
    let empty = (Vec::new(), Vec::new());
    let result = match (kind, geom) {
        (CurveKind::Line, SurfaceGeom::Plane { o, n, .. }) => {
            let (loc, v) = match line_geometry(curve) {
                Some(x) => x,
                None => return Some(empty),
            };
            line_plane(&loc, &v, o, n)
        }
        (CurveKind::Line, SurfaceGeom::Sphere { o, r, .. }) => {
            let (loc, v) = match line_geometry(curve) {
                Some(x) => x,
                None => return Some(empty),
            };
            line_sphere(&loc, &v, o, *r)
        }
        (CurveKind::Line, SurfaceGeom::Cylinder { a, z, r, .. }) => {
            let (loc, v) = match line_geometry(curve) {
                Some(x) => x,
                None => return Some(empty),
            };
            line_cylinder(&loc, &v, a, z, *r)
        }
        (CurveKind::Line, SurfaceGeom::Cone { apex, z, cosa, .. }) => {
            let (loc, v) = match line_geometry(curve) {
                Some(x) => x,
                None => return Some(empty),
            };
            line_cone(&loc, &v, apex, z, *cosa)
        }
        (CurveKind::Line, SurfaceGeom::Torus { o, major, minor, x, z, .. }) => {
            let (loc, v) = match line_geometry(curve) {
                Some(x) => x,
                None => return Some(empty),
            };
            let dir = match GpDir::from_vec(&v) {
                Ok(d) => d,
                Err(_) => return Some(empty),
            };
            let lin = GpLin::from_pnt_dir(loc, dir);
            let ax3 = match (
                GpDir::from_vec(z),
                GpDir::from_vec(x),
            ) {
                (Ok(zd), Ok(xd)) => match GpAx3::new(*o, zd, &xd) {
                    Ok(a) => a,
                    Err(_) => return Some(empty),
                },
                _ => return Some(empty),
            };
            let torus = match GpTorus::new(ax3, *major, *minor) {
                Ok(t) => t,
                Err(_) => return Some(empty),
            };
            let pts = line_torus_intersect(&lin, &torus);
            let v2 = v.square_magnitude();
            SolveResult::Roots(pts.into_iter().map(|p| GpVec::from_pnts(&loc, &p).dot(&v) / v2).collect())
        }
        (CurveKind::Circle, SurfaceGeom::Plane { o, n, .. }) => {
            let Some(ConicGeom::Circle { o: co, r, x, y }) = conic_geometry(curve, kind) else {
                return None;
            };
            let dq = GpVec::from_pnts(&co, o);
            trig_solve(r * n.dot(&x), r * n.dot(&y), n.dot(&dq), t0, t1, ANALYTIC_TOL)
        }
        (CurveKind::Circle, SurfaceGeom::Sphere { o, r: sph_r, .. }) => {
            let Some(ConicGeom::Circle { o: co, r, x, y }) = conic_geometry(curve, kind) else {
                return None;
            };
            let d = GpVec::from_pnts(o, &co); // circle center − sphere center
            let a = 2.0 * r * d.dot(&x);
            let b = 2.0 * r * d.dot(&y);
            let cc = sph_r * sph_r - r * r - d.square_magnitude();
            trig_solve(a, b, cc, t0, t1, ANALYTIC_TOL)
        }
        (CurveKind::Ellipse, SurfaceGeom::Plane { o, n, .. }) => {
            let Some(ConicGeom::Ellipse { o: eo, a, b, x, y }) = conic_geometry(curve, kind) else {
                return None;
            };
            let dq = GpVec::from_pnts(&eo, o);
            trig_solve(a * n.dot(&x), -b * n.dot(&y), n.dot(&dq), t0, t1, ANALYTIC_TOL)
        }
        _ => return None,
    };

    Some(match result {
        SolveResult::Roots(roots) => {
            let pts = roots
                .into_iter()
                .filter_map(|w| {
                    if w < t0 - ANALYTIC_TOL || w > t1 + ANALYTIC_TOL {
                        return None;
                    }
                    compute_append_point(curve, surface, Some(geom), w, uv, ANALYTIC_TOL)
                })
                .collect();
            (pts, Vec::new())
        }
        SolveResult::OnSurface => (Vec::new(), on_surface_segment(curve, surface, geom, t0, t1)),
        SolveResult::None => (Vec::new(), Vec::new()),
    })
}

/// Build an `On` coincidence segment spanning `[t0, t1]`.
pub(super) fn on_surface_segment(
    curve: &dyn Curve,
    surface: &dyn Surface,
    geom: &SurfaceGeom,
    t0: f64,
    t1: f64,
) -> Vec<IntersectionSegment> {
    let p0 = curve.d0(t0);
    let p1 = curve.d0(t1);
    let (u0, v0) = surface_params(surface, Some(geom), &p0);
    let (u1, v1) = surface_params(surface, Some(geom), &p1);
    let a = IntersectionPoint::new(t0, u0, v0, p0, State::On);
    let b = IntersectionPoint::new(t1, u1, v1, p1, State::On);
    vec![IntersectionSegment::new(a, b)]
}

/// Validate a candidate parameter and build an [`IntersectionPoint`]: the
/// curve point must lie on the surface (within `tol`) and its surface
/// parameters must fall inside `uv`.
pub(super) fn compute_append_point(
    curve: &dyn Curve,
    surface: &dyn Surface,
    geom: Option<&SurfaceGeom>,
    w: f64,
    uv: (f64, f64, f64, f64),
    tol: f64,
) -> Option<IntersectionPoint> {
    let pnt = curve.d0(w);
    let (su, sv) = surface_params(surface, geom, &pnt);
    let (u0, u1, v0, v1) = uv;
    if su < u0 - tol || su > u1 + tol || sv < v0 - tol || sv > v1 + tol {
        return None;
    }
    // Validate the point lies on the surface. Analytic points are checked
    // against the quadric's implicit equation (exact, and independent of the
    // reported (u, v)); the general path uses the normal (signed) distance.
    let on_surface = match geom {
        Some(g) => quadric_on_surface(g, &pnt, tol.max(1e-6)),
        None => signed_dist(surface, None, &pnt).abs() <= tol.max(1e-6),
    };
    if !on_surface {
        return None;
    }
    let state = compute_state(curve, surface, geom, w, su, sv);
    Some(IntersectionPoint::new(w, su, sv, pnt, state))
}

/// Exact outward normal of a reconstructed quadric at a point on it. Used for
/// the transition classification so tangencies are not corrupted by the
/// first-order error of the finite-difference normal (which tilts the normal by
/// a curvature term).
pub(super) fn quadric_normal(geom: &SurfaceGeom, p: &GpPnt) -> Option<GpVec> {
    match geom {
        SurfaceGeom::Plane { n, .. } => {
            let m = n.magnitude();
            if m < 1e-12 {
                None
            } else {
                Some(n.divided(m))
            }
        }
        SurfaceGeom::Sphere { o, .. } => {
            let d = GpVec::from_pnts(o, p);
            let m = d.magnitude();
            if m < 1e-12 {
                None
            } else {
                Some(d.divided(m))
            }
        }
        SurfaceGeom::Cylinder { a, z, .. } => {
            let d = GpVec::from_pnts(a, p);
            let radial = d.coord.subtracted(&z.xyz().multiplied(d.dot(z)));
            let m = radial.modulus();
            if m < 1e-12 {
                None
            } else {
                Some(GpVec::from_xyz(&radial.divided(m)))
            }
        }
        SurfaceGeom::Cone { apex, z, cosa, .. } => {
            let d = GpVec::from_pnts(apex, p);
            // f = |d|²·cos²α − (d·Z)²; ∇f ∝ cos²α·d − (d·Z)·Z.
            let grad = d.coord.multiplied(cosa * cosa).subtracted(&z.xyz().multiplied(d.dot(z)));
            let m = grad.modulus();
            if m < 1e-12 {
                None
            } else {
                Some(GpVec::from_xyz(&grad.divided(m)))
            }
        }
        SurfaceGeom::Torus { o, major, z, .. } => {
            let d = GpVec::from_pnts(o, p);
            let dz = d.dot(z);
            let radial = d.coord.subtracted(&z.xyz().multiplied(dz));
            let rho = radial.modulus();
            if rho < 1e-12 {
                return None;
            }
            // f = (ρ − R)² + z² − r²; ∇f ∝ (ρ − R)/ρ·d_perp + z·Z.
            let g = radial.multiplied((rho - major) / rho).added(&z.xyz().multiplied(dz));
            let m = g.modulus();
            if m < 1e-12 {
                None
            } else {
                Some(GpVec::from_xyz(&g.divided(m)))
            }
        }
    }
}

/// Whether a point lies on a reconstructed quadric, via the quadric's implicit
/// equation (relative tolerance). Independent of the surface parameterization,
/// so it stays exact for quadrics whose reconstructed reference differs from
/// the natural `(u, v)` origin.
pub(super) fn quadric_on_surface(geom: &SurfaceGeom, p: &GpPnt, tol: f64) -> bool {
    match geom {
        SurfaceGeom::Plane { o, n, .. } => {
            let d = GpVec::from_pnts(o, p);
            d.dot(n).abs() / n.magnitude() <= tol
        }
        SurfaceGeom::Sphere { o, r, .. } => (p.distance(o) - r).abs() <= tol,
        SurfaceGeom::Cylinder { a, z, r, .. } => {
            let d = GpVec::from_pnts(a, p);
            let rho = d.coord.subtracted(&z.xyz().multiplied(d.dot(z))).modulus();
            (rho - r).abs() <= tol
        }
        SurfaceGeom::Cone { apex, z, cosa, .. } => {
            let d = GpVec::from_pnts(apex, p);
            let lhs = d.dot(z) * d.dot(z);
            let rhs = d.square_magnitude() * cosa * cosa;
            (lhs - rhs).abs() <= tol * d.square_magnitude().max(1.0)
        }
        SurfaceGeom::Torus { o, major, minor, z, .. } => {
            let d = GpVec::from_pnts(o, p);
            let dz = d.dot(z);
            let rho = d.coord.subtracted(&z.xyz().multiplied(dz)).modulus();
            let lhs = (rho - major) * (rho - major) + dz * dz;
            (lhs - minor * minor).abs() <= tol * (rho + major + minor).max(1.0)
        }
    }
}

/// Transition of the curve at the crossing: aligned with the surface normal →
/// `Out`, opposed → `In`, perpendicular → `On`.
pub(super) fn compute_state(
    curve: &dyn Curve,
    surface: &dyn Surface,
    geom: Option<&SurfaceGeom>,
    w: f64,
    su: f64,
    sv: f64,
) -> State {
    let (_, d1) = curve.d1(w);
    let dm = d1.magnitude();
    if dm < 1e-12 {
        return State::Unknown;
    }
    let n = match geom.and_then(|g| quadric_normal(g, &curve.d0(w))) {
        Some(n) => n,
        None => surface_normal(surface, su, sv),
    };
    let nm = n.magnitude();
    if nm < 1e-12 {
        return State::Unknown;
    }
    let cos = d1.dot(&n) / (dm * nm);
    if -cos > ANG_TOL {
        State::In
    } else if cos > ANG_TOL {
        State::Out
    } else {
        State::On
    }
}

// ---------------------------------------------------------------------------
// General path (sampling + bisection)
// ---------------------------------------------------------------------------

/// Signed distance from a point to `surface`: positive when the point is on the
/// side of the surface normal, negative otherwise. When the reconstructed
/// quadric geometry is available its exact normal is used (the trait-level
/// `surface_normal` finite-difference fallback is ill-defined on unbounded
/// parameter directions); otherwise the sampling normal is used.
pub(super) fn signed_dist(surface: &dyn Surface, geom: Option<&SurfaceGeom>, p: &GpPnt) -> f64 {
    let (su, sv) = surface_closest_params(surface, p, 16, 16);
    let q = surface.d0(su, sv);
    let n = match geom.and_then(|g| quadric_normal(g, p)) {
        Some(n) => n,
        None => surface_normal(surface, su, sv),
    };
    let d = GpVec::from_pnts(&q, p);
    let nm = n.magnitude();
    if nm < 1e-12 {
        d.magnitude()
    } else {
        d.dot(&n) / nm
    }
}

/// General sampling path: uniform curve sampling, sign-change bisection for
/// transverse crossings, golden-section refinement of near-zero distance
/// minima for tangencies.
pub(super) fn general_curve_surface(
    curve: &dyn Curve,
    surface: &dyn Surface,
    geom: Option<&SurfaceGeom>,
    t0: f64,
    t1: f64,
    uv: (f64, f64, f64, f64),
    tol: f64,
) -> Vec<IntersectionPoint> {
    let n = 512usize;
    let mut points: Vec<IntersectionPoint> = Vec::new();
    let mut prev_u = t0;
    let mut prev_s = signed_dist(surface, geom, &curve.d0(t0));
    for i in 1..=n {
        let u = t0 + (t1 - t0) * i as f64 / n as f64;
        let s = signed_dist(surface, geom, &curve.d0(u));
        if prev_s * s < 0.0 {
            // Transverse crossing.
            if let Some(w) = bisect_root(curve, surface, geom, prev_u, u, prev_s, tol) {
                if let Some(pt) = compute_append_point(curve, surface, geom, w, uv, tol) {
                    points.push(pt);
                }
            }
        } else if s.abs() <= tol.max(1e-6) && prev_s.abs() > s.abs() {
            // Near-zero minimum → likely tangency.
            if let Some(w) = refine_min(curve, surface, geom, prev_u, u, tol) {
                if let Some(mut pt) = compute_append_point(curve, surface, geom, w, uv, tol) {
                    pt.state = State::On;
                    points.push(pt);
                }
            }
        }
        prev_u = u;
        prev_s = s;
    }
    points
}

/// Bisect a sign change of the signed distance over `[ua, ub]`.
pub(super) fn bisect_root(
    curve: &dyn Curve,
    surface: &dyn Surface,
    geom: Option<&SurfaceGeom>,
    ua: f64,
    ub: f64,
    sa: f64,
    tol: f64,
) -> Option<f64> {
    let mut lo = ua;
    let mut hi = ub;
    let mut flo = sa;
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        let sm = signed_dist(surface, geom, &curve.d0(mid));
        if flo * sm <= 0.0 {
            hi = mid;
        } else {
            lo = mid;
            flo = sm;
        }
        if (hi - lo).abs() < 1e-10 {
            break;
        }
        if sm.abs() < 1e-10 {
            break;
        }
    }
    let w = 0.5 * (lo + hi);
    let p = curve.d0(w);
    if signed_dist(surface, geom, &p).abs() <= tol.max(1e-6) {
        Some(w)
    } else {
        None
    }
}

/// Golden-section refinement of the unsigned distance minimum over `[lo, hi]`.
pub(super) fn refine_min(
    curve: &dyn Curve,
    surface: &dyn Surface,
    geom: Option<&SurfaceGeom>,
    lo: f64,
    hi: f64,
    tol: f64,
) -> Option<f64> {
    let f = |u: f64| signed_dist(surface, geom, &curve.d0(u)).abs();
    let (u, d) = golden_1d(&f, lo, hi, 1e-10);
    if d <= tol.max(1e-6) {
        Some(u)
    } else {
        None
    }
}

/// Golden-section minimization of `f` over `[lo, hi]`. Returns `(argmin, min)`.
pub(super) fn golden_1d<F: Fn(f64) -> f64>(f: &F, lo: f64, hi: f64, eps: f64) -> (f64, f64) {
    pub(super) const GOLD: f64 = 0.618_033_988_749_894_9;
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
