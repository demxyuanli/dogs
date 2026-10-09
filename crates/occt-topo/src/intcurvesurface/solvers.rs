use super::prelude::*;
use super::*;
use occt_math::MathFunctionSetRoot;

/// Surface `(u, v)` parameters of a point — `IntCurveSurface_InterUtils.pxx:899-925`
/// (`ComputeParamsOnQuadric`): switch on `Adaptor3d_Surface::GetType()` and take
/// `ElSLib::Parameters` on the surface's **own** placement
/// (`SurfaceTool::Plane/Cylinder/Cone/Sphere`). The port reads the placement off
/// the surface's own gp accessors, so the parameters share the surface's U/V
/// origin and are exact (a reconstructed frame can sit at a different origin).
///
/// OCCT's switch has no Torus case (`default: break`), but the line-torus arm
/// reaches this helper through `ProcessLinTorus`, whose source
/// `IntAna_IntLinTorus::Perform` stores the torus parameters from the same
/// `ElSLib::Parameters(gp_Torus, P)` call (`IntAna_IntLinTorus.cxx:98-113`;
/// consumed at `IntCurveSurface_InterUtils.pxx:1283-1316`), so Torus maps to
/// that source.
///
/// Returns `None` for a surface outside the five elementary types. OCCT never
/// selects this arm for such a surface (`Adaptor3d_Surface::GetType` routes it
/// to the polygon/polyhedron path), so the caller must not invent parameters.
pub(super) fn surface_params(surface: &dyn Surface, p: &GpPnt) -> Option<(f64, f64)> {
    if let Some(pln) = surface.gp_pln() {
        return Some(slib::plane_parameters(&pln.pos, p));
    }
    if let Some(cyl) = surface.gp_cylinder() {
        return Some(slib::cylinder_parameters(&cyl.pos, p));
    }
    if let Some(cone) = surface.gp_cone() {
        return Some(slib::cone_parameters(&cone.pos, cone.radius(), cone.semi_angle(), p));
    }
    if let Some(sph) = surface.gp_sphere() {
        return Some(slib::sphere_parameters(&sph.pos, p));
    }
    if let Some(tor) = surface.gp_torus() {
        return Some(slib::torus_parameters(&tor.pos, tor.major_radius(), tor.minor_radius(), p));
    }
    None
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
    // `IntCurveSurface_InterImpl::PerformBounds` (`Inter.pxx:105-182`) routes
    // every non-quadric surface (torus, B-spline, offset, ...) through the
    // polygon-polyhedron interference regardless of the curve type; only
    // plane / cylinder / cone / sphere have the analytic quadric arms.
    let polyhedron_path = matches!(surface_kind, SurfaceKind::Torus | SurfaceKind::Other);
    if curve_kind != CurveKind::Other {
        if let Some((mut pts, mut segs)) = geom
            .as_ref()
            .and_then(|g| perform_conic_surf(curve, curve_kind, surface, g, t0, t1, uv))
        {
            points.append(&mut pts);
            segments.append(&mut segs);
        } else if polyhedron_path && finite_window {
            let mut pts = polyhedron_curve_surface(curve, surface, uv);
            points.append(&mut pts);
        } else if finite_window {
            // `InternalPerformCurveQuadric` (`Inter.pxx:423-439` +
            // `IntCurveSurface_InterUtils.pxx:1243-1274`): non-conic curve
            // against a plane / cylinder / cone / sphere.
            let mut pts = quadric_curve_exact(curve, surface, geom.as_ref(), uv);
            points.append(&mut pts);
        }
    } else if polyhedron_path && finite_window {
        let mut pts = polyhedron_curve_surface(curve, surface, uv);
        points.append(&mut pts);
    } else if finite_window {
        // `InternalPerformCurveQuadric` (`Inter.pxx:423-439` +
        // `IntCurveSurface_InterUtils.pxx:1243-1274`).
        let mut pts = quadric_curve_exact(curve, surface, geom.as_ref(), uv);
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

/// Whether `(kind, geom)` is an analytic conic/quadric arm that OCCT routes
/// through `IntAna_IntConicQuad`.
///
/// `PerformConicSurfCircle` (`IntCurveSurface_Inter.pxx:711-759`),
/// `PerformConicSurfEllipse` (`:761-809`), `PerformConicSurfParabola`
/// (`:811-887`) and `PerformConicSurfHyperbola` (`:890-966`) each switch on
/// `SurfaceTool::GetType` and hand the conic plus the surface quadric to
/// `IntAna_IntConicQuad`. The arms already solved directly in
/// [`perform_conic_surf`] (`line` against Plane/Cylinder/Cone/Sphere/Torus,
/// `circle`/`ellipse` against Plane, `circle` against Sphere) are excluded:
/// they reach the same result with the same `IntAna` formulas, so they keep
/// their existing exact solves. Everything else in the four switches is
/// dispatched to `IntAna_IntConicQuad`.
fn uses_int_ana(kind: CurveKind, geom: &SurfaceGeom) -> bool {
    match (kind, geom) {
        (CurveKind::Circle, SurfaceGeom::Cylinder { .. })
        | (CurveKind::Circle, SurfaceGeom::Cone { .. })
        | (CurveKind::Ellipse, SurfaceGeom::Cylinder { .. })
        | (CurveKind::Ellipse, SurfaceGeom::Cone { .. })
        | (CurveKind::Ellipse, SurfaceGeom::Sphere { .. })
        | (CurveKind::Parabola, SurfaceGeom::Plane { .. })
        | (CurveKind::Parabola, SurfaceGeom::Cylinder { .. })
        | (CurveKind::Parabola, SurfaceGeom::Cone { .. })
        | (CurveKind::Parabola, SurfaceGeom::Sphere { .. })
        | (CurveKind::Hyperbola, SurfaceGeom::Plane { .. })
        | (CurveKind::Hyperbola, SurfaceGeom::Cylinder { .. })
        | (CurveKind::Hyperbola, SurfaceGeom::Cone { .. })
        | (CurveKind::Hyperbola, SurfaceGeom::Sphere { .. }) => true,
        _ => false,
    }
}

/// `SurfaceTool::Plane/Cylinder/Cone/Sphere(theSurface)` as an
/// `IntAna_Quadric` (`IntAna_IntConicQuad.hxx:41-48`): the implicit quadric
/// the conic arms pass to `IntAna_IntConicQuad`. The conic/plane overloads of
/// `IntAna_IntConicQuad` only forward to this same quadric solve
/// (`IntAna_IntConicQuad.cxx:562-575`), so the port calls the quadric
/// constructor directly.
fn surface_quadric(surface: &dyn Surface) -> Option<IntAnaQuadric> {
    if let Some(p) = surface.gp_pln() {
        return Some(IntAnaQuadric::from_plane(&p));
    }
    if let Some(c) = surface.gp_cylinder() {
        return Some(IntAnaQuadric::from_cylinder(&c));
    }
    if let Some(c) = surface.gp_cone() {
        return Some(IntAnaQuadric::from_cone(&c));
    }
    if let Some(s) = surface.gp_sphere() {
        return Some(IntAnaQuadric::from_sphere(&s));
    }
    None
}

/// `IntCurveSurface_InterUtils::ProcessIntAna` (`IntCurveSurface_InterUtils.pxx:1188-1227`)
/// plus `IntCurveSurface_HInter::AppendIntAna` (`IntCurveSurface_Inter.pxx:967-1002`):
/// every `IntAna_IntConicQuad` point becomes an intersection point through
/// `ComputeAppendPoint`, after its surface parameters are read by
/// `ComputeParamsOnQuadric`.
///
/// `ProcessIntAna` sets `theIsParallel` and appends no point when the conic is
/// in the quadric or parallel to it; `IntCurvesFace_Intersector` reads that
/// back as `IsParallel()` (`IntCurvesFace_Intersector.cxx:361`) and the port
/// represents it by the `On` segment (`int_curves_face.rs:171`).
fn process_int_ana(
    curve: &dyn Curve,
    surface: &dyn Surface,
    geom: &SurfaceGeom,
    ia: &IntAnaIntConicQuad,
    t0: f64,
    t1: f64,
    uv: (f64, f64, f64, f64),
) -> (Vec<IntersectionPoint>, Vec<IntersectionSegment>) {
    if !ia.is_done() {
        return (Vec::new(), Vec::new());
    }
    if ia.is_in_quadric() || ia.is_parallel() {
        return (Vec::new(), on_surface_segment(curve, surface, t0, t1));
    }
    let mut pts = Vec::new();
    for i in 1..=ia.nb_points() {
        let w = ia.param_on_conic(i);
        // The analytic solve searches the conic's own period, so the caller's
        // parameter window is the only filter left (same window the direct
        // arms apply to their roots).
        if w < t0 - ANALYTIC_TOL || w > t1 + ANALYTIC_TOL {
            continue;
        }
        if let Some(pt) = compute_append_point(curve, surface, Some(geom), w, uv, ANALYTIC_TOL) {
            pts.push(pt);
        }
    }
    (pts, Vec::new())
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

    // `IntAna_IntConicQuad` arms of the four `PerformConicSurf*` switches.
    if uses_int_ana(kind, geom) {
        let Some(quad) = surface_quadric(surface) else {
            return Some(empty);
        };
        let ia = match kind {
            CurveKind::Circle => match curve.gp_circ() {
                Some(c) => IntAnaIntConicQuad::circle_quadric(&c, &quad),
                None => return Some(empty),
            },
            CurveKind::Ellipse => match curve.gp_ellipse() {
                Some(e) => IntAnaIntConicQuad::ellipse_quadric(&e, &quad),
                None => return Some(empty),
            },
            CurveKind::Parabola => match curve.gp_parabola() {
                Some(p) => IntAnaIntConicQuad::parabola_quadric(&p, &quad),
                None => return Some(empty),
            },
            CurveKind::Hyperbola => match curve.gp_hyperbola() {
                Some(h) => IntAnaIntConicQuad::hyperbola_quadric(&h, &quad),
                None => return Some(empty),
            },
            _ => return Some(empty),
        };
        return Some(process_int_ana(curve, surface, geom, &ia, t0, t1, uv));
    }

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
            // `ElCLib::EllipseValue` uses `+Minor*sin(U)*YDir` (`cxx:176-189`),
            // so the Y coefficient carries the same sign as `b`.
            trig_solve(a * n.dot(&x), b * n.dot(&y), n.dot(&dq), t0, t1, ANALYTIC_TOL)
        }
        // `PerformConicSurfParabola` / `PerformConicSurfHyperbola`
        // (`IntCurveSurface_Inter.pxx:811-966`) default arm: the quadric
        // Plane / Cylinder / Cone / Sphere cases were taken above by
        // `uses_int_ana`, so any other surface falls through to the general
        // (polygon/polyhedron) arm, exactly as OCCT's `default` does.
        (CurveKind::Parabola, _) | (CurveKind::Hyperbola, _) => return None,
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
        SolveResult::OnSurface => (Vec::new(), on_surface_segment(curve, surface, t0, t1)),
        SolveResult::None => (Vec::new(), Vec::new()),
    })
}

/// Build an `On` coincidence segment spanning `[t0, t1]`.
pub(super) fn on_surface_segment(
    curve: &dyn Curve,
    surface: &dyn Surface,
    t0: f64,
    t1: f64,
) -> Vec<IntersectionSegment> {
    let p0 = curve.d0(t0);
    let p1 = curve.d0(t1);
    let Some((u0, v0)) = surface_params(surface, &p0) else {
        return Vec::new();
    };
    let Some((u1, v1)) = surface_params(surface, &p1) else {
        return Vec::new();
    };
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
    let Some((su, sv)) = surface_params(surface, &pnt) else {
        return None;
    };
    let (u0, u1, v0, v1) = uv;
    if su < u0 - tol || su > u1 + tol || sv < v0 - tol || sv > v1 + tol {
        return None;
    }
    // Validate the point lies on the surface. Analytic points are checked
    // against the quadric's implicit equation (exact, and independent of the
    // reported (u, v)). Without a reconstructed quadric there is nothing to
    // check against: OCCT's `ComputeAppendPoint` (`InterUtils.pxx:1116-1176`)
    // itself only validates the parameter windows and periodicity, never a
    // point-on-surface residual, and the `GetType`-based dispatch above never
    // reaches this with a non-elementary surface.
    let on_surface = match geom {
        Some(g) => quadric_on_surface(g, &pnt, tol.max(1e-6)),
        None => true,
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
// Quadric arm (exact roots of the implicit distance)
// ---------------------------------------------------------------------------

/// `InternalPerformCurveQuadric` (`Inter.pxx:423-439`) via
/// `IntCurveSurface_InterUtils::PerformCurveQuadric`
/// (`IntCurveSurface_InterUtils.pxx:1243-1274`): a non-conic curve against a
/// plane / cylinder / cone / sphere.
///
/// `TheQuadCurvExactHInter` runs `math_FunctionAllRoots` on the signed distance
/// `Q(w)` over each C1 interval of the curve; every returned root becomes an
/// intersection point through `ComputeAppendPoint`. OCCT searches the curve's
/// own range (no caller window) and lets `ComputeAppendPoint` filter by the
/// surface box, so no curve-parameter filter is applied here.
pub(super) fn quadric_curve_exact(
    curve: &dyn Curve,
    surface: &dyn Surface,
    geom: Option<&SurfaceGeom>,
    uv: (f64, f64, f64, f64),
) -> Vec<IntersectionPoint> {
    let mut out: Vec<IntersectionPoint> = Vec::new();
    let exact = TheQuadCurvExactHInter::new(surface, curve);
    if !exact.is_done() {
        return out;
    }
    let nb_roots = exact.nb_roots();
    for i in 1..=nb_roots {
        let w = exact.root(i);
        if let Some(pt) = compute_append_point(curve, surface, geom, w, uv, ANALYTIC_TOL) {
            out.push(pt);
        }
    }
    out
}

// ---------------------------------------------------------------------------
// General path (polyhedron interference + exact root refinement)
// ---------------------------------------------------------------------------

/// `IntCurveSurface_InterImpl::PerformBounds` default arm plus
/// `InternalPerformPolygonBounds` / `InternalPerform`
/// (`IntCurveSurface_Inter.pxx:105-182, 262-523`) for a non-quadric surface:
///
/// 1. `DecomposeSurfaceIntervals` splits the surface into C2 rectangles;
/// 2. `SamplePars` + `ThePolygonOfHInter` build the curve polygon per C2
///    interval of the curve;
/// 3. `ThePolyhedronOfHInter` builds the surface grid with the
///    `Adaptor3d_HSurfaceTool::NbSamplesU/V` counts;
/// 4. `TheInterferenceOfHInter` collects start points, which
///    `TheExactHInter` refines with `math_FunctionSetRoot`.
///
/// Unlike the sampling fallback this never projects a point onto the surface:
/// the `(u, v)` come straight out of the polyhedron
/// (`SectionPointToParameters`).
pub(super) fn polyhedron_curve_surface(
    curve: &dyn Curve,
    surface: &dyn Surface,
    uv: (f64, f64, f64, f64),
) -> Vec<IntersectionPoint> {
    /// `defl` / `NbMin` of `PerformBounds` (`Inter.pxx:155-157`).
    const DEFL: f64 = 0.1;
    const NB_MIN: usize = 10;

    let mut intervals = Vec::new();
    decompose_surface_intervals(surface, &mut intervals);

    let mut result: Vec<IntersectionPoint> = Vec::new();
    for iv in &intervals {
        let mut u1 = iv.u0;
        let mut u2 = iv.u1;
        let mut v1 = iv.v0;
        let mut v2 = iv.v1;
        clamp_uv_parameters(&mut u1, &mut u2, &mut v1, &mut v2);

        // The port intersects one face's UV box rather than the whole surface,
        // so the caller window is intersected with the C2 interval. Both the
        // polyhedron (inside `InternalPerformPolygonBounds`) and the Newton
        // domain (`ProcessSortedPoints`) receive this box.
        let (cu0, cv0, cu1, cv1) = uv;
        let (u1, u2) = (u1.max(cu0), u2.min(cu1));
        let (v1, v2) = (v1.max(cv0), v2.min(cv1));
        if u2 - u1 <= 1e-15 || v2 - v1 <= 1e-15 {
            continue;
        }

        // `InternalPerformPolygonBounds` (`Inter.pxx:461-477`).
        let mut nbsu = surface_nb_samples_u_range(surface, u1, u2);
        let mut nbsv = surface_nb_samples_v_range(surface, v1, v2);
        if nbsu > 40 {
            nbsu = 40;
        }
        if nbsv > 40 {
            nbsv = 40;
        }
        let polyhedron = ThePolyhedronOfHInter::new(surface, nbsu.max(1), nbsv.max(1), u1, v1, u2, v2);

        // `PerformBounds` (`Inter.pxx:139-178`): one polygon per C2 interval.
        let nb_intervals = curve.nb_intervals(2);
        let mut pars_list: Vec<Vec<f64>> = Vec::new();
        if nb_intervals > 1 {
            let tab_w = curve.parameter_intervals(2);
            for i in 0..nb_intervals as usize {
                pars_list.push(sample_pars(curve, tab_w[i], tab_w[i + 1], DEFL, NB_MIN));
            }
        } else {
            pars_list.push(sample_pars(
                curve,
                curve.first_parameter(),
                curve.last_parameter(),
                DEFL,
                NB_MIN,
            ));
        }

        for pars in &pars_list {
            if pars.len() < 2 {
                continue;
            }
            let polygon = ThePolygonOfHInter::with_params(curve, pars);
            let interference = TheInterferenceOfHInter::of_polygon_polyhedron(&polygon, &polyhedron);
            let mut start_points = SortedStartPoints::new();
            collect_interference_points(&interference, &polyhedron, &polygon, &mut start_points);
            sort_start_points(&mut start_points);

            // `InternalPerform` (`Inter.pxx:358-416`).
            let func = TheCSFunctionOfHInter::new(surface, curve);
            let mut exact = TheExactHInter::new(func, THE_TOLTANGENCY);
            let mut rsnld = MathFunctionSetRoot::with_iterations(exact.function(), 100);

            let mut pts: Vec<IntersectionPoint> = Vec::new();
            process_sorted_points(
                &mut exact,
                &mut rsnld,
                &start_points,
                u1,
                u2,
                v1,
                v2,
                polygon.inf_parameter(),
                polygon.sup_parameter(),
                curve,
                surface,
                &mut pts,
            );
            result.extend(pts);
        }
    }
    result
}
