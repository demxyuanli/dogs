use super::prelude::*;


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

/// Cylinder parameters `(center, unit axis, radius)` from sampled invariants:
/// the axis is along `∂S/∂v` and two opposite points of a `u`-ring give the
/// radius (their midpoint lies on the axis).
pub(super) fn cylinder_params(s: &dyn Surface) -> Option<(GpPnt, GpVec, f64)> {
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
    let center = mid(&q0, &q1);
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
    let m0 = mid(&s.d0(0.0, 1.0), &s.d0(PI, 1.0));
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
    let center = mid(&e0, &e1);
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

pub(super) fn mid(a: &GpPnt, b: &GpPnt) -> GpPnt {
    GpPnt::new(0.5 * (a.x() + b.x()), 0.5 * (a.y() + b.y()), 0.5 * (a.z() + b.z()))
}

pub(super) fn dist_to_axis(p: &GpPnt, c: &GpPnt, ax: &GpVec) -> f64 {
    let rel = GpVec::from_pnts(c, p);
    rel.subtracted(&ax.multiplied_scalar(rel.dot(ax))).magnitude()
}

pub(super) fn pnt_add_vec(p: GpPnt, v: &GpVec) -> GpPnt {
    GpPnt::from_xyz(&p.coord.added(&v.coord))
}

/// Closest point between two (possibly skew) lines; `None` when parallel or
/// too far apart.
pub(super) fn closest_point_lines(p1: GpPnt, d1: GpVec, p2: GpPnt, d2: GpVec) -> Option<GpPnt> {
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

/// Exact UV projection of a 3D point onto an analytic surface, recovered from
/// the surface's own parameterization so the results agree with
/// `edge_pcurve_on_face`.
pub(super) enum Projector {
    Plane { o: GpPnt, x: GpVec, y: GpVec },
    Cylinder { o: GpPnt, x: GpVec, y: GpVec, z: GpVec },
    Cone { o: GpPnt, x: GpVec, y: GpVec, z: GpVec, alpha: f64 },
    Sphere { o: GpPnt, x: GpVec, y: GpVec, z: GpVec, r: f64 },
    Torus { o: GpPnt, x: GpVec, y: GpVec, z: GpVec, r_big: f64, r_small: f64 },
}

impl Projector {
    pub(super) fn from_surface(s: &dyn Surface, kind: SurfaceKind) -> Option<Projector> {
        match kind {
            SurfaceKind::Plane => {
                let (o, x, y) = plane_axes(s)?;
                Some(Projector::Plane { o, x, y })
            }
            SurfaceKind::Cylinder => {
                let (o, x, y, z) = cylinder_axes(s)?;
                Some(Projector::Cylinder { o, x, y, z })
            }
            SurfaceKind::Cone => {
                let (apex, ax, alpha) = cone_params(s)?;
                let (x, y, z) = cone_frame(s, &apex, &ax)?;
                // The projection origin must be the surface's placement plane
                // point (v=0), not its apex: `cone_value` is parameterized from
                // the placement (ElSLib::ConeValue), so an apex-based origin
                // would offset the projected v by RefRadius/sin(α) and the UV
                // window would not be the exact inverse of `surface.d0`.
                let o = cone_placement(s, &apex, &z)?;
                Some(Projector::Cone { o, x, y, z, alpha })
            }
            SurfaceKind::Sphere => {
                let center = sphere_center(s)?;
                let r = s.d0(0.0, 0.0).distance(&center);
                if r < 1e-12 {
                    return None;
                }
                let x = GpVec::from_pnts(&center, &s.d0(0.0, 0.0)).normalized();
                let y = GpVec::from_pnts(&center, &s.d0(FRAC_PI_2, 0.0)).normalized();
                let z = x.crossed(&y).normalized();
                Some(Projector::Sphere { o: center, x, y, z, r })
            }
            SurfaceKind::Torus => {
                let (center, ax, r_big, r_small) = torus_params(s)?;
                let x = GpVec::from_pnts(&center, &s.d0(0.0, 0.0)).normalized();
                let y = GpVec::from_pnts(&center, &s.d0(FRAC_PI_2, 0.0)).normalized();
                // The axis must be the surface's *actual* revolution axis (the
                // `ElSLib::TorusParameters` `gp_Ax3` direction), recovered by
                // `torus_params` from `top − bot` — not `x × y`, which flips for
                // a left-handed STEP torus and would make `v = atan2(z·d, …)` land
                // on the wrong half of the tube (the screw's thread got v < 0 for
                // points above its centre).
                let z = ax;
                Some(Projector::Torus { o: center, x, y, z, r_big, r_small })
            }
            SurfaceKind::Other => None,
        }
    }

    /// Whether the surface's `u` parameter is periodic (revolution surfaces).
    pub(super) fn is_u_periodic(&self) -> bool {
        matches!(
            self,
            Projector::Cylinder { .. } | Projector::Cone { .. } | Projector::Sphere { .. } | Projector::Torus { .. }
        )
    }

    /// Whether the surface's `v` parameter is periodic (torus tube direction).
    pub(super) fn is_v_periodic(&self) -> bool {
        matches!(self, Projector::Torus { .. })
    }

    /// Project a single 3D point onto the face `(u, v)` domain. `u` is wrapped
    /// into the principal angle range `(-π, π]`; use [`Self::project_seq`] to
    /// unwrap a continuous sequence across the seam.
    pub(super) fn project(&self, p: &GpPnt) -> GpPnt2d {
        match self {
            Projector::Plane { o, x, y } => {
                let z = x.crossed(y);
                if let (Ok(zd), Ok(xd)) = (GpDir::from_vec(&z), GpDir::from_vec(x)) {
                    if let Ok(ax) = GpAx3::new(*o, zd, &xd) {
                        return occt_geom::projlib::project_pln_pnt(&GpPln::new(ax), p);
                    }
                }
                let d = p.coord.subtracted(&o.coord);
                let xx = x.xyz().dot(x.xyz());
                let yy = y.xyz().dot(y.xyz());
                let xy = x.xyz().dot(y.xyz());
                let dx = d.dot(x.xyz());
                let dy = d.dot(y.xyz());
                let det = xx * yy - xy * xy;
                if det.abs() < 1e-24 {
                    GpPnt2d::new(dx, dy)
                } else {
                    GpPnt2d::new((dx * yy - dy * xy) / det, (dy * xx - dx * xy) / det)
                }
            }
            Projector::Cylinder { o, x, y, z } => {
                let d = GpVec::from_pnts(o, p);
                GpPnt2d::new(d.dot(y).atan2(d.dot(x)), d.dot(z))
            }
            Projector::Cone { o, x, y, z, alpha } => {
                let d = GpVec::from_pnts(o, p);
                GpPnt2d::new(d.dot(y).atan2(d.dot(x)), d.dot(z) / alpha.cos())
            }
            Projector::Sphere { o, x, y, z, r } => {
                let d = GpVec::from_pnts(o, p);
                let v = (d.dot(z) / r).clamp(-1.0, 1.0).asin();
                GpPnt2d::new(d.dot(y).atan2(d.dot(x)), v)
            }
            Projector::Torus { o, x, y, z, r_big, r_small } => {
                let d = GpVec::from_pnts(o, p);
                let dz = d.dot(z);
                // Signed radial distance `R + r·cos v`. The unsigned |R + r·cos v|
                // folds a spindle torus (R < r) through its inner equator (v ≈ π),
                // where the radial distance is negative; `atan2(dz, |rho| − R)`
                // then returns v ≈ 0 instead of v ≈ π. Recover the sign from the
                // torus identity |d|² = rho² + dz² ⇒ rho = (|d|² + R² − r²)/(2R).
                let rho = (d.square_magnitude() + r_big.powi(2) - r_small.powi(2)) / (2.0 * r_big);
                GpPnt2d::new(d.dot(y).atan2(d.dot(x)), dz.atan2(rho - r_big))
            }
        }
    }

    /// Project a sequence of edge points, unwrapping the `u` (and, for a torus,
    /// the `v`) coordinate across the seam so the result stays continuous and on
    /// the same sheet of a self-intersecting surface.
    pub(super) fn project_seq(&self, curve: &dyn Curve, t_vals: &[f64]) -> Vec<GpPnt2d> {
        let mut out = Vec::with_capacity(t_vals.len());
        let mut prev_u = f64::NAN;
        let mut prev_v = f64::NAN;
        for &t in t_vals {
            let q = self.project(&curve.d0(t));
            let mut u = q.x();
            let mut v = q.y();
            if self.is_u_periodic() && prev_u.is_finite() {
                while u - prev_u > PI {
                    u -= 2.0 * PI;
                }
                while u - prev_u < -PI {
                    u += 2.0 * PI;
                }
            }
            if self.is_v_periodic() && prev_v.is_finite() {
                while v - prev_v > PI {
                    v -= 2.0 * PI;
                }
                while v - prev_v < -PI {
                    v += 2.0 * PI;
                }
            }
            prev_u = u;
            prev_v = v;
            out.push(GpPnt2d::new(u, v));
        }
        out
    }
}

/// Analytic projection of a 3D point onto an analytic surface
/// (plane/cylinder/cone/sphere/torus) — the OCCT `ProjLib` inverse. Returns
/// `None` for non-analytic surfaces (B-spline, …), where the caller must fall
/// back to a numeric inversion.
///
/// Unlike a Newton iterate seeded from the surface range centre, this is the
/// exact inverse of `d0` even on unbounded-v surfaces (a cylinder/cone point far
/// from v=0), where Newton from the centre cannot converge.
pub fn project_point_on_surface(s: &dyn Surface, p: &GpPnt) -> Option<GpPnt2d> {
    let kind = classify_surface_kind(s);
    let proj = Projector::from_surface(s, kind)?;
    Some(proj.project(p))
}

/// Continuous projection of a 3D curve onto an analytic surface — the OCCT
/// `ProjLib` `projected curve` route. Each sample is projected and then
/// unwrapped across the `u`/`v` seam so the polyline stays on the *same sheet*
/// of a self-intersecting surface (a spindle torus), which a point-wise
/// projection cannot guarantee. Restricted to the torus: for the other analytic
/// surfaces the point-wise grid search is already exact (their `d0` inverse has
/// no sheet ambiguity), so the caller keeps that path. Returns `None` otherwise.
/// `ShapeConstruct_ProjectCurveOnSurface.cxx:71`.
const THE_NCONTROL: usize = 23;

/// `generateCurvePoints` count (`cxx:531-556`).
fn generate_curve_point_count(curve: &dyn Curve, first: f64, last: f64) -> usize {
    let mut n = THE_NCONTROL;
    if curve.bspline_poles().is_some() {
        let knots = curve.parameter_intervals(0);
        let mut used = 0usize;
        for w in knots.windows(2) {
            if w[1] > first && w[0] < last {
                used += 1;
            }
        }
        let deg = curve.nurbs_degree().unwrap_or(1);
        let min_pnt = used * (deg + 1);
        while n < min_pnt {
            n += THE_NCONTROL - 1;
        }
    }
    n.max(2)
}

/// `ShapeConstruct_ProjectCurveOnSurface::projectAnalytic`
/// (`cxx:848-897`): the plane behind the surface itself, or behind one
/// `Geom_RectangularTrimmedSurface` / `Geom_OffsetSurface` wrapper. The wrapper
/// is resolved with `else if`, exactly as the C++ does, so only one level is
/// peeled.
///
/// The projected curve returned by `ProjLib_ProjectedCurve` on the wrapper is
/// the same `(u, v)` as the plane's because both wrappers keep the basis
/// parameterization: `Geom_RectangularTrimmedSurface::Value` clamps the
/// arguments but returns the basis point, and
/// `Geom_OffsetSurface::Value = basis + Offset*Normal`, so only the origin of
/// the frame moves (along the normal, which the `(u, v)` projection ignores).
fn plane_behind_wrapper(surf: &dyn Surface) -> Option<GpPln> {
    if let Some(pln) = surf.gp_pln() {
        return Some(pln);
    }
    if let Some(basis) = surf.rectangular_trimmed_basis() {
        if let Some(pln) = basis.gp_pln() {
            return Some(pln);
        }
    } else if let Some(basis) = surf.offset_basis_surface() {
        if let Some(pln) = basis.gp_pln() {
            return Some(pln);
        }
    }
    None
}

/// `ShapeConstruct_ProjectCurveOnSurface::projectAnalytic(aCurve3DTrim)`
/// (`cxx:848-897`): a `None` result means "no analytic plane projection", and
/// the caller falls back to the sampled approximation of `Perform`.
fn project_analytic(curve: &dyn Curve, surf: &dyn Surface) -> Option<Arc<dyn Curve2d>> {
    let pln = plane_behind_wrapper(surf)?;
    crate::pcurve::project_curve_on_plane(curve, &pln)
}

/// `GeomAdaptor_Surface::load` (`GeomAdaptor_Surface.cxx:423-430`) recurses
/// through every `Geom_RectangularTrimmedSurface` and stores the *basis*
/// geometry, so `SurfAdapt.Cylinder()` / `Cone()` / `Sphere()` / `Torus()`
/// answer for a trimmed face. The Rust `Surface` trait accessors do not peel
/// that wrapper, so the analytic arms of [`value_of_uv`] must look through it.
fn analytic_basis(s: &dyn Surface) -> Option<Arc<dyn Surface>> {
    let mut deepest: Option<Arc<dyn Surface>> = None;
    let mut cur = s.rectangular_trimmed_basis();
    while let Some(b) = cur {
        cur = b.rectangular_trimmed_basis();
        deepest = Some(b);
    }
    deepest
}

/// `ShapeAnalysis_Surface::ValueOfUV` (`ShapeAnalysis_Surface.cxx:1245-1515`):
/// the closest `(u, v)` of a 3D point on a surface.
///
/// `cxx:1262-1293` is a type switch on `Adaptor3d_Surface::GetType()`: for the
/// five analytic surfaces the answer is `ElSLib::Parameters` directly, with
/// `ShapeAnalysis::AdjustByPeriod` applied to each periodic parameter so the
/// value lands next to the middle of the surface range (`cxx:1271` cylinder,
/// `cxx:1276` cone, `cxx:1281` sphere, `cxx:1287-1290` torus U and V). The
/// switch is on the *type*, so a `Geom_BSplineSurface` that is only
/// geometrically a cylinder must not take these arms.
///
/// `cxx:1294-1476` is the general arm: `Extrema_ExtPS` seeded over the bounds
/// (`cxx:1346-1352`, extended by one resolution) and then polished with
/// `SurfaceNewton` (`cxx:1390-1401`).
/// [`occt_geom::extrema_surf`] ports both (`Extrema_GenExtPS` grid seeds +
/// Newton on `F = ((S-P)·Su, (S-P)·Sv)`), so it replaces the coarse
/// grid + axis-descent of `brep_surface::surface_closest_params`, whose
/// coordinate descent stalls ~1e-3 short along a curved valley.
///
/// This is the `mpt[i] == 3`-independent numeric inversion the whole
/// `approxPCurve` population loop (`cxx:1310-1464`) and `getLine`
/// (`cxx:902-1102`) run on. `ShapeAnalysis_Surface::NextValueOfUV`
/// (`cxx:1164-1243`) is the same solve seeded from the previous `(u, v)`; the
/// extrema path is used here because it is global (no seed continuity needed
/// for the sampled point sequence the caller builds).
///
/// PORTED (this function): the window extension
/// `du = min(myUDelt, UResolution(preci))` /
/// `dv = min(myVDelt, VResolution(preci))` (`cxx:1340-1348`, both zero for
/// offset surfaces, `myUDelt = myVDelt = 0.01` at `cxx:106-107`), fed to
/// `Extrema_ExtPS` over `[uf - du, ul + du] x [vf - dv, vl + dv]`
/// (`cxx:1352`), and `RestrictBounds` (`cxx:1323`, `cxx:64-91`) before it.
/// The closest point is therefore allowed to leave the face's parameter range,
/// which is what OCCT does when the minimum sits on the domain border.
///
/// PORTED (this function): the residue branch after the extrema minimum
/// (`cxx:1381-1437`) - the extrema distance is recomputed from the actual
/// surface value (`cxx:1374-1378`, the OCC486 workaround), the local-minimum
/// normal test `possLockal` (`cxx:1402-1419`), and [`uv_from_iso`]
/// (`cxx:1423` / `cxx:1459`, both in place of a general closest-point search),
/// which wins the result when it lands closer (`cxx:1428-1437`).
///
/// NOT PORTED (reported):
/// - `myExtOK` caching (`cxx:1327`): OCCT builds the `Extrema_ExtPS` object
///   once per surface; [`occt_geom::extrema_surf::point_surface_extrema_box`]
///   re-seeds its grid on every call. The result is the same, only the cost
///   differs. `value_of_uv` is a free function, so no cache state exists here.
/// - `ComputeBoxes` (`cxx:1898-1922`) and the `anIsoBox->Distance(aPBox) >
///   theMin` early-outs (`cxx:1607-1610`, `cxx:1676-1678`): they skip a
///   candidate whose bounding box is already farther than the best distance, so
///   they cannot change the result - only the cost.
///
/// PORTED (this function): `ShapeAnalysis_Surface::SurfaceNewton`
/// (`ShapeAnalysis_Surface.cxx:1065-1156`) at both of its `ValueOfUV` call
/// sites - the conic-extrusion attempt (`cxx:1309`) and the residue polish
/// (`cxx:1390-1401`).
pub(super) fn value_of_uv(s: &dyn Surface, p: &GpPnt, preci: f64) -> GpPnt2d {
    // `cxx:1254` `Bounds(uf, ul, vf, vl)` (`ShapeAnalysis_Surface.lxx:54-64`,
    // the surface's own parameter bounds buffered at construction).
    let (uf, ul) = s.u_range();
    let (vf, vl) = s.v_range();
    let two_pi = 2.0 * PI;
    // `cxx:1262-1293` reads `SurfAdapt.Cylinder()` / `Cone()` / `Sphere()` /
    // `Torus()`, and `GeomAdaptor_Surface::load` (`GeomAdaptor_Surface.cxx:423-430`)
    // answers those from the basis of a `Geom_RectangularTrimmedSurface`. The
    // period-adjustment target below stays the trimmed surface's own window
    // (`ShapeAnalysis_Surface.lxx:54-63`).
    let trimmed_basis = analytic_basis(s);
    let analytic: &dyn Surface = trimmed_basis.as_deref().unwrap_or(s);
    match classify_surface_kind(s) {
        // `cxx:1265-1268`: no periodic parameter, so no adjustment.
        SurfaceKind::Plane => {
            if let Some(pl) = analytic.gp_pln() {
                let (u, v) = occt_core::elib::slib::plane_parameters(&pl.pos, p);
                return GpPnt2d::new(u, v);
            }
        }
        // `cxx:1269-1275`.
        SurfaceKind::Cylinder => {
            if let Some(cy) = analytic.gp_cylinder() {
                let (mut u, v) = occt_core::elib::slib::cylinder_parameters(&cy.pos, p);
                u += crate::shhealing::adjust_by_period(u, 0.5 * (uf + ul), two_pi);
                return GpPnt2d::new(u, v);
            }
        }
        // `cxx:1276-1280`.
        SurfaceKind::Cone => {
            if let Some(co) = analytic.gp_cone() {
                let (mut u, v) =
                    occt_core::elib::slib::cone_parameters(&co.pos, co.radius, co.semi_angle, p);
                u += crate::shhealing::adjust_by_period(u, 0.5 * (uf + ul), two_pi);
                return GpPnt2d::new(u, v);
            }
        }
        // `cxx:1281-1285`.
        SurfaceKind::Sphere => {
            if let Some(sp) = analytic.gp_sphere() {
                let (mut u, v) = occt_core::elib::slib::sphere_parameters(&sp.pos, p);
                u += crate::shhealing::adjust_by_period(u, 0.5 * (uf + ul), two_pi);
                return GpPnt2d::new(u, v);
            }
        }
        // `cxx:1286-1293`, the only arm that adjusts V as well.
        SurfaceKind::Torus => {
            if let Some(to) = analytic.gp_torus() {
                let (mut u, mut v) = occt_core::elib::slib::torus_parameters(
                    &to.pos,
                    to.major_radius,
                    to.minor_radius,
                    p,
                );
                u += crate::shhealing::adjust_by_period(u, 0.5 * (uf + ul), two_pi);
                v += crate::shhealing::adjust_by_period(v, 0.5 * (vf + vl), two_pi);
                return GpPnt2d::new(u, v);
            }
        }
        SurfaceKind::Other => {}
    }
    // `cxx:1300-1301` `S = (uf + ul) / 2; T = (vf + vl) / 2;`: the general arm
    // starts from the middle of the window, which is also the answer it returns
    // when the extrema fail (`cxx:1454-1471`).
    let mut uf = uf;
    let mut ul = ul;
    let mut vf = vf;
    let mut vl = vl;
    let mut s_par = 0.5 * (uf + ul);
    let mut t_par = 0.5 * (vf + vl);
    // `cxx:1306-1321`: a `Geom_SurfaceOfExtrusion` with an infinite U range (a
    // conic extrusion - the basis curve is a line) is first solved by
    // `SurfaceNewton` from the middle of the window; on failure `uf, ul` become
    // `-500, 500` before `RestrictBounds`, which would otherwise take the
    // width-2000 window from the finite end.
    if s.is_surface_of_linear_extrusion()
        && occt_core::precision::Precision::is_infinite(uf)
        && occt_core::precision::Precision::is_infinite(ul)
    {
        // `cxx:1309-1313`: the Newton solution is returned outright.
        if let Some((solution, _)) =
            surface_newton(s, &GpPnt2d::new(s_par, t_par), p, preci)
        {
            return solution;
        }
        uf = -500.0;
        ul = 500.0;
    }
    // `cxx:1324` `RestrictBounds(uf, ul, vf, vl)`: an infinite range becomes
    // `[-1000, 1000]`, or a width-2000 window from the finite end.
    restrict_bounds(&mut uf, &mut ul);
    restrict_bounds(&mut vf, &mut vl);
    // `cxx:1332-1348`: the search window is the (restricted) bounds expanded by
    // `du = min(myUDelt, UResolution(preci))` / `dv = min(myVDelt, VResolution(preci))`;
    // the expansion is zero for an offset surface because evaluating one
    // outside its range can be undefined (id23943).
    //
    // `myUDelt` / `myVDelt` are NOT the constructor values
    // (`ShapeAnalysis_Surface.cxx:106-107`, `0.01`): the first general-arm
    // `ValueOfUV` forces `IsUClosed()` / `IsVClosed()` (`cxx:1332-1338`) which
    // overwrite them with `|ul - uf| / 20` (`cxx:672`), reduced by half the
    // resolution at the minimum boundary gap (`cxx:764`, `:786`, `:818`,
    // `:849`; `:959`, `:981`, `:1013`, `:1044`, `0` when the surface's own
    // `IsUClosed` already says closed).
    let mut du = 0.0;
    let mut dv = 0.0;
    if !s.is_offset_surface() {
        // `GeomAdaptor_Surface::UResolution` / `VResolution`
        // (`GeomAdaptor_Surface.cxx:1818`, `:1900`).
        du = sa_u_delt(s).min(occt_geom::approx_same_parameter::u_resolution(s, preci));
        dv = sa_v_delt(s).min(occt_geom::approx_same_parameter::v_resolution(s, preci));
    }
    let e = occt_geom::extrema_surf::point_surface_extrema_box(
        s,
        p,
        uf - du,
        ul + du,
        vf - dv,
        vl + dv,
    );
    let Some(v_min) = e.v2.filter(|v| e.u2.is_finite() && v.is_finite()) else {
        // `cxx:1448-1472`: no extrema solution, so `myGap = UVFromIso(P3D,
        // preci, UU, VV)` from the middle-of-window seed wins outright.
        let (_, iu, iv) = uv_from_iso(s, p, preci, s_par, t_par);
        return GpPnt2d::new(iu, iv);
    };
    s_par = e.u2;
    t_par = v_min;
    // `cxx:1376-1378`: the distance is recomputed from the actual surface value
    // (the OCC486 workaround - a surface of revolution can return a stale
    // distance from its extrema solution).
    let mut dis_surf = p.distance(&s.d0(s_par, t_par));
    if dis_surf > preci {
        // `cxx:1390-1401` PRO7226 #412920: the extrema foot is polished by
        // `SurfaceNewton` and the closer of the two keeps the answer.
        if let Some((pp, _)) = surface_newton(s, &GpPnt2d::new(s_par, t_par), p, preci) {
            let dist = p.distance(&s.value(pp.x(), pp.y()));
            if dist < dis_surf {
                dis_surf = dist;
                s_par = pp.x();
                t_par = pp.y();
            }
        }
        // `cxx:1402-1419`: a residue below `10 * preci` whose foot already lies
        // on the normal through `P3D` (`possLockal`) is a genuine local minimum,
        // and the iso search below is skipped for it.
        let mut poss_lockal = false;
        if dis_surf < 10.0 * preci && s.continuity() != 0 {
            // `GeomAbs_C0 == 0` (`GeomAbs_Shape`).
            let (pnt, d1u, d1v) = s.d1(s_par, t_par);
            let b = d1u.crossed(&d1v);
            let a = GpVec::from_pnts(&pnt, p);
            let ab = a.dot(&b);
            let nrm2 = b.square_magnitude();
            if nrm2 > 1e-10 {
                poss_lockal = a.square_magnitude() - (ab * ab) / nrm2
                    < occt_core::precision::CONFUSION * occt_core::precision::CONFUSION;
            }
        }
        if !poss_lockal {
            // `cxx:1423` (and `cxx:1459` on the failure path above).
            let (iso_dist, iu, iv) = uv_from_iso(s, p, preci, s_par, t_par);
            // `cxx:1428-1437`: the iso solution replaces the extrema one only
            // when it is closer.
            if dis_surf > iso_dist {
                s_par = iu;
                t_par = iv;
            }
        }
    }
    GpPnt2d::new(s_par, t_par)
}

/// [`value_of_uv`] plus the `myGap` that `ShapeAnalysis_Surface::ValueOfUV`
/// leaves in `Gap()`, which `ShapeConstruct_ProjectCurveOnSurface::getLine`
/// reads back at `cxx:973` / `:987` to grow the running `aTol2`.
///
/// `ValueOfUV` opens with `myGap = -1.` (`ShapeAnalysis_Surface.cxx:1249`) and
/// closes with `if (myGap <= 0) myGap = P3D.Distance(SurfAdapt.Value(S, T));`
/// (`cxx:1502-1513`). The five analytic arms (`cxx:1262-1293`) never assign it,
/// so the tail fills it with the distance between the sample and the surface
/// point of the returned parameters; the general arm stores the same quantity
/// as `min(disSurf, DistMinOnIso)` (`cxx:1432-1437`) or the `UVFromIso`
/// minimum (`cxx:1459`). The one path that returns with `myGap = -1.` is the
/// infinite-`U` `GeomAbs_SurfaceOfExtrusion` Newton arm (`cxx:1306-1313`),
/// which leaves before the tail - see [`conic_extrusion_newton_arm`].
pub(super) fn value_of_uv_with_gap(s: &dyn Surface, p: &GpPnt, preci: f64) -> (GpPnt2d, f64) {
    let q = value_of_uv(s, p, preci);
    let gap = if conic_extrusion_newton_arm(s, p, preci) {
        -1.0
    } else {
        s.d0(q.x(), q.y()).distance(p)
    };
    (q, gap)
}

/// `ShapeAnalysis_Surface::ValueOfUV` (`cxx:1306-1313`): a
/// `Geom_SurfaceOfExtrusion` with an infinite `U` range (the basis curve is a
/// line - the conic case) is solved by `SurfaceNewton` from the middle of the
/// window and `return solution;` leaves the function before the `myGap` tail
/// (`cxx:1502-1513`), so `Gap()` keeps the `-1.` of `cxx:1249`.
/// On Newton failure the arm continues into the general path and the tail runs.
fn conic_extrusion_newton_arm(s: &dyn Surface, p: &GpPnt, preci: f64) -> bool {
    if !s.is_surface_of_linear_extrusion() {
        return false;
    }
    let (uf, ul) = s.u_range();
    if !(occt_core::precision::Precision::is_infinite(uf)
        && occt_core::precision::Precision::is_infinite(ul))
    {
        return false;
    }
    let (vf, vl) = s.v_range();
    surface_newton(
        s,
        &GpPnt2d::new(0.5 * (uf + ul), 0.5 * (vf + vl)),
        p,
        preci,
    )
    .is_some()
}

/// `ShapeAnalysis_Surface::SurfaceNewton` (`ShapeAnalysis_Surface.cxx:1065-1156`).
///
/// The Newton step of the stationary-point system
/// `F = ((S(u, v) - P3D) . Su, (S(u, v) - P3D) . Sv)` solved from the second
/// derivatives `SurfAdapt.D2` (`cxx:1088`). PORTED here: the `D2` normal test
/// (`cxx:1091-1095`), the discriminant and step (`cxx:1101-1115`), the window
/// rejection `[uf - du, ul + du] x [vf - dv, vl + dv]` (`cxx:1117-1120`, the
/// window is the one `Bounds` + `UResolution`/`VResolution` give, `cxx:1073-1077`),
/// the convergence test `|du| + |dv| <= max(1e-12, (U + V) * 10e-16)`
/// (`cxx:1126-1130`) and the two solution rejections (`cxx:1134-1137`
/// `rs2 > rsfirst`, `cxx:1140-1143` `rs2 - rsn^2 / nrm2 > Tol^2`).
///
/// Returns `(sol, res)`; `res == 2` is OCCT's "strange attractor" guard
/// (`cxx:1151` `nrm2 < 0.01 * ru2 * rv2`) which `NextValueOfUV` treats as
/// "never trust this solution without the iso fallback" (`cxx:1220`). `None` is
/// the C++ `false` return (0).
pub(super) fn surface_newton(
    s: &dyn Surface,
    p2d_prev: &GpPnt2d,
    p3d: &GpPnt,
    preci: f64,
) -> Option<(GpPnt2d, i32)> {
    // `cxx:1072-1077`: `Bounds` widened by `SurfAdapt.UResolution(preci)` /
    // `VResolution(preci)` (`GeomAdaptor_Surface::UResolution`, the port in
    // `occt_geom::approx_same_parameter`). NOTE: the surface's own
    // `uv_resolution` is `None` for the analytic family (plane, cylinder, cone,
    // sphere, torus), which is why the resolution port is called directly here
    // instead of falling back to `Precision::Parametric(preci)`.
    let (uf, ul) = s.u_range();
    let (vf, vl) = s.v_range();
    let du_res = occt_geom::approx_same_parameter::u_resolution(s, preci);
    let dv_res = occt_geom::approx_same_parameter::v_resolution(s, preci);
    let u_f = uf - du_res;
    let u_l = ul + du_res;
    let v_f = vf - dv_res;
    let v_l = vl + dv_res;
    // `cxx:1079-1080`
    let tol = occt_core::precision::CONFUSION;
    let tol2 = tol * tol;
    let mut u = p2d_prev.x();
    let mut v = p2d_prev.y();
    // `cxx:1082` `rsfirst = P3D.XYZ() - Value(U, V).XYZ()`
    let rs_first = GpVec::from_pnts(&s.value(u, v), p3d);
    for _ in 0..25 {
        // `cxx:1088` `SurfAdapt.D2(U, V, pnt, ru, rv, ruu, rvv, ruv)`.
        let (pnt, ru, rv, ruu, rvv, ruv) = s.d2(u, v);
        // `cxx:1090-1096`.
        let ru2 = ru.dot(&ru);
        let rv2 = rv.dot(&rv);
        let n = ru.crossed(&rv);
        let nrm2 = n.square_magnitude();
        if nrm2 < 1e-10 || occt_core::precision::Precision::is_positive_infinite(nrm2) {
            break;
        }
        // `cxx:1098-1107`.
        let rs = GpVec::from_pnts(&pnt, p3d);
        let r_suu = rs.dot(&ruu);
        let r_svv = rs.dot(&rvv);
        let r_suv = rs.dot(&ruv);
        let d = -nrm2 + rv2 * r_suu + ru2 * r_svv - 2.0 * r_suv * ru.dot(&rv)
            + r_suv * r_suv
            - r_suu * r_svv;
        if d.abs() < 1e-10 {
            break;
        }
        // `cxx:1110-1115`.
        let fract = 1.0 / d;
        let step_u = rs.dot(
            &n.crossed(&rv)
                .added(&ru.multiplied_scalar(r_svv))
                .subtracted(&rv.multiplied_scalar(r_suv)),
        ) * fract;
        let step_v = rs.dot(
            &ru.crossed(&n)
                .added(&rv.multiplied_scalar(r_suu))
                .subtracted(&ru.multiplied_scalar(r_suv)),
        ) * fract;
        u += step_u;
        v += step_v;
        // `cxx:1117-1120`.
        if u < u_f || u > u_l || v < v_f || v > v_l {
            break;
        }
        // `cxx:1126-1130`.
        let a_resolution = 1e-12f64.max((u + v) * 10e-16);
        if step_u.abs() + step_v.abs() > a_resolution {
            continue;
        }
        // `cxx:1133-1137` PRO10109 4517.
        let rs2 = rs.square_magnitude();
        if rs2 > rs_first.square_magnitude() {
            break;
        }
        // `cxx:1140-1143`.
        let r_sn = rs.dot(&n);
        if rs2 - r_sn * r_sn / nrm2 > tol2 {
            break;
        }
        // `cxx:1147-1153`.
        let res = if nrm2 < 0.01 * ru2 * rv2 { 2 } else { 1 };
        return Some((GpPnt2d::new(u, v), res));
    }
    None
}

/// `GeomAdaptor_Surface::UContinuity` / `VContinuity` for a B-spline
/// (`GeomAdaptor_Surface.cxx:539-587` / `:591-639`) through
/// `LocalContinuity` (`cxx:110-165`): the direction is `GeomAbs_C0` when the
/// largest multiplicity of the unique knots strictly inside the window reaches
/// the degree (`cxx:139-153` wraps the end knots out of the range with the
/// `EpsKnot` nudges at `cxx:124-132`).
fn bspline_c0(direction_knots: &[f64], degree: usize, first: f64, last: f64) -> bool {
    let (knots, mults) = occt_geom::bspline_surface::GeomBSplineSurface::unique_knots_mults(
        direction_knots,
    );
    if knots.len() < 3 {
        return false;
    }
    let mut mult_max = 0i32;
    for (k, m) in knots.iter().zip(mults.iter()).skip(1).take(knots.len() - 2) {
        if *k > first && *k < last && *m > mult_max {
            mult_max = *m;
        }
    }
    (degree as i32) - mult_max <= 0
}

/// The `Geom_BSplineSurface` behind `s`, peeling `Geom_RectangularTrimmedSurface`
/// wrappers (`GeomAdaptor_Surface::load` dispatches a trimmed surface to its
/// basis type, so `SurfAdapt.BSpline()` is the basis, `Surface.cxx:423-430`).
fn bspline_behind(s: &dyn Surface, depth: usize) -> Option<occt_geom::bspline_surface::GeomBSplineSurface> {
    if let Some(basis) = s.rectangular_trimmed_basis() {
        if depth < 8 {
            return bspline_behind(basis.as_ref(), depth + 1);
        }
        return None;
    }
    s.osculating_bspline()
}

/// `ShapeAnalysis_Surface::NextValueOfUV` (`ShapeAnalysis_Surface.cxx:1164-1243`):
/// `SurfaceNewton` seeded from the previous `(u, v)` (the hint), used both by
/// `ShapeConstruct_ProjectCurveOnSurface::projectPoint` (`cxx:953`, `cxx:1410`)
/// and by the sample walk (`cxx:1358`, `cxx:1432`).
///
/// PORTED (this function): the type gate (`cxx:1174-1180`: `SurfaceNewton` is
/// only reachable for Bezier / B-spline / extrusion / revolution / offset
/// surfaces, every other type jumps straight to `ValueOfUV`, `cxx:1240`), the
/// C0-knot bail-out for a B-spline whose hint sits on a C0 knot
/// (`cxx:1182-1213`), and the `maxpreci` guard of the accepted solution
/// (`cxx:1216-1235`): `res == 2`, or a gap more than `Confusion` above
/// `maxpreci`, sends the result through `UVFromIso` and keeps it only when it is
/// at least as close as the iso minimum.
///
/// Returns `(uv, gap)` with the `myGap` OCCT leaves behind (`cxx:1226` / `:1232`
/// / through `ValueOfUV` at `cxx:1240`). The `ValueOfUV` fall-through returns
/// the distance recomputed at the returned `(u, v)`; see the note on
/// [`value_of_uv_with_gap`].
pub(super) fn next_value_of_uv(
    s: &dyn Surface,
    p2d_prev: &GpPnt2d,
    p3d: &GpPnt,
    preci: f64,
    maxpreci: f64,
) -> (GpPnt2d, f64) {
    use crate::geom_bnd_lib_surface3d::SurfaceKind as TypeKind;
    let takes_newton = matches!(
        crate::geom_bnd_lib_surface3d::surface_kind(s),
        TypeKind::BezierSurface
            | TypeKind::BSplineSurface
            | TypeKind::SurfaceOfExtrusion
            | TypeKind::SurfaceOfRevolution
            | TypeKind::OffsetSurface
    );
    if takes_newton {
        // `cxx:1182-1213`: near a C0 knot the Newton step is unreliable, so the
        // global solve is used instead.
        if let Some(bs) = bspline_behind(s, 0) {
            let (uf, ul) = s.u_range();
            let (vf, vl) = s.v_range();
            let on_u_knot = bspline_c0(&bs.knots_u, bs.deg_u, uf, ul)
                && bs
                    .knots_u
                    .iter()
                    .any(|k| (k - p2d_prev.x()).abs() < occt_core::precision::CONFUSION);
            let on_v_knot = bspline_c0(&bs.knots_v, bs.deg_v, vf, vl)
                && bs
                    .knots_v
                    .iter()
                    .any(|k| (k - p2d_prev.y()).abs() < occt_core::precision::CONFUSION);
            if on_u_knot || on_v_knot {
                return value_of_uv_with_gap(s, p3d, preci);
            }
        }
        // `cxx:1215-1237`.
        if let Some((sol, res)) = surface_newton(s, p2d_prev, p3d, preci) {
            let gap = p3d.distance(&s.value(sol.x(), sol.y()));
            if res == 2 || (maxpreci > 0.0 && gap - maxpreci > occt_core::precision::CONFUSION) {
                // `cxx:1226-1230`: `UVFromIso` writes back into `U` / `V`.
                let (my_gap, u, v) = uv_from_iso(s, p3d, preci, sol.x(), sol.y());
                if gap >= my_gap {
                    return (GpPnt2d::new(u, v), my_gap);
                }
            }
            return (sol, gap);
        }
    }
    // `cxx:1240`.
    value_of_uv_with_gap(s, p3d, preci)
}

/// `RestrictBounds` (`ShapeAnalysis_Surface.cxx:64-91`): clamp an infinite
/// parameter range to a finite window - `[-1000, 1000]` when both ends are
/// infinite, otherwise a width of 2000 from the finite end.
fn restrict_bounds(first: &mut f64, last: &mut f64) {
    let is_first_inf = occt_core::precision::Precision::is_negative_infinite(*first);
    let is_last_inf = occt_core::precision::Precision::is_positive_infinite(*last);
    if is_first_inf || is_last_inf {
        if is_first_inf && is_last_inf {
            *first = -1000.0;
            *last = 1000.0;
        } else if is_first_inf {
            *first = *last - 2000.0;
        } else {
            *last = *first + 2000.0;
        }
    }
}

// ---------------------------------------------------------------------------
// ShapeAnalysis_Surface::IsUClosed / IsVClosed
// ---------------------------------------------------------------------------

/// `Geom_Surface::IsUClosed` union with `IsUPeriodic`.
///
/// `ShapeAnalysis_Surface::IsUClosed` (`ShapeAnalysis_Surface.cxx:674`) and
/// `IsVClosed` (`cxx:881`) start from the concrete surface's own flag. For a
/// B-spline that flag is `Geom_BSplineSurface::IsUClosed`
/// (`Geom_BSplineSurface_1.cxx:1350-1369`): `IsUPeriodic()` first, then the
/// boundary `UIso` curves compared with `IsEqual(.., Confusion)`. The port's
/// [`Surface::is_u_closed`] implements a pole-row approximation of the second
/// term (`crates/occt-geom/src/bspline_surface.rs`) but omits the periodic
/// term, so the periodic flag is added here. For every other surface type
/// `is_u_closed` already defaults to `is_u_periodic`, so this is a no-op.
fn geom_u_closed(s: &dyn Surface) -> bool {
    s.is_u_closed() || s.is_u_periodic()
}

/// `Geom_Surface::IsVClosed` union with `IsVPeriodic`; see [`geom_u_closed`].
fn geom_v_closed(s: &dyn Surface) -> bool {
    s.is_v_closed() || s.is_v_periodic()
}

/// `Geom_BSplineSurface::UMultiplicity(1)` / `(NbUKnots)` / the same for V:
/// the first / last unique-knot multiplicities from the flat knot vector.
fn end_mults(flat: &[f64]) -> (i32, i32) {
    let (_, mults) = occt_geom::bspline_surface::GeomBSplineSurface::unique_knots_mults(flat);
    let first = mults.first().copied().unwrap_or(0);
    let last = mults.last().copied().unwrap_or(0);
    (first, last)
}

/// `ShapeAnalysis_Surface::IsUClosed` (`ShapeAnalysis_Surface.cxx:660-863`).
///
/// Returns `(myUCloseVal, anUmidVal, myUDelt)`. `myUCloseVal` is the maximum
/// 3D distance between the two U-boundary curves (`cxx:853-854`), `anUmidVal`
/// the distance from the first boundary sample to the mid-U sample of the
/// surface (`-1` = OCCT's "not computed", so the symmetry test
/// `myUCloseVal > sqrt(anUmidVal)` is skipped), and `myUDelt` the
/// `Extrema_ExtPS` window extension `ValueOfUV` reads as
/// `du = min(myUDelt, UResolution(preci))` (`cxx:1346`): a twentieth of the U
/// range (`cxx:672`), reduced by half the U resolution at the minimum boundary
/// gap for the sampled arms (`cxx:764`, `:786`, `:818`, `:849`).
///
/// The `AN_ASSERT(myUCloseVal >= 0)` guard of `cxx:853` is the `myUCloseVal < 0`
/// cache test of `cxx:667`; the cache is not reproduced (each call recomputes -
/// the value is preci-independent, so the boolean is unchanged; see the module
/// notes on the stateless port).
fn sa_u_close_state(s: &dyn Surface) -> (f64, f64, f64) {
    use crate::geom_bnd_lib_surface3d::SurfaceKind as TypeKind;
    // `cxx:668-671`: `Bounds` + `RestrictBounds`.
    let (mut uf, mut ul) = s.u_range();
    let (vf, vl) = s.v_range();
    restrict_bounds(&mut uf, &mut ul);
    // `cxx:672`.
    let mut delt = (ul - uf).abs() / 20.0;
    // `cxx:674-681`.
    if geom_u_closed(s) {
        return (0.0, -1.0, 0.0);
    }
    // `cxx:684-689`: a `Geom_RectangularTrimmedSurface` forces the `default`
    // (101-point) arm even though the adaptor already sees the basis type.
    let surftype = if s.rectangular_trimmed_basis().is_some() {
        TypeKind::OtherSurface
    } else {
        crate::geom_bnd_lib_surface3d::surface_kind(s)
    };
    let mut close = f64::MAX; // `RealLast()`
    let mut mid = -1.0;
    match surftype {
        // `cxx:692-694`.
        TypeKind::Plane => {}
        // `cxx:695-717`.
        TypeKind::SurfaceOfExtrusion => {
            if let Some(crv) = s.extrusion_basis_curve() {
                let f = crv.first_parameter();
                let l = crv.last_parameter();
                if !occt_core::precision::Precision::is_infinite(f)
                    && !occt_core::precision::Precision::is_infinite(l)
                {
                    let p1 = crv.d0(f);
                    let p2 = crv.d0(l);
                    close = p1.square_distance(&p2);
                    let pm = crv.d0(0.5 * (f + l));
                    mid = p1.square_distance(&pm);
                }
            }
        }
        // `cxx:719-788`.
        TypeKind::BSplineSurface => {
            if let Some(bs) = s.osculating_bspline() {
                let nbup = bs.poles.len();
                let mut distmin = f64::MAX;
                let (m1, m2) = end_mults(&bs.knots_u);
                if bs.u_periodic {
                    // `cxx:724-728`.
                    close = 0.0;
                    delt = 0.0;
                } else if nbup < 3 {
                    // `cxx:729-731`.
                } else if bs.is_rational()
                    || m1 != bs.deg_u as i32 + 1
                    || m2 != bs.deg_u as i32 + 1
                {
                    // `cxx:732-765`: sample along the V knots when the boundary
                    // is not a simple pole row (rational, or a boundary knot
                    // multiplicity that is not `degree + 1`).
                    let (vknots, _) = occt_geom::bspline_surface::GeomBSplineSurface::unique_knots_mults(&bs.knots_v);
                    let nbvk = vknots.len();
                    let mut v = vknots[0];
                    let p1 = s.value(uf, v);
                    let p2 = s.value(ul, v);
                    close = p1.square_distance(&p2);
                    let pm = s.value(0.5 * (uf + ul), v);
                    mid = p1.square_distance(&pm);
                    distmin = close;
                    for i in 1..nbvk {
                        // `cxx:747-749`: `v = 0.5 * (VKnot(i-1) + VKnot(i))`.
                        v = 0.5 * (vknots[i - 1] + vknots[i]);
                        let p1 = s.value(uf, v);
                        let p2 = s.value(ul, v);
                        let a_dist = p1.square_distance(&p2);
                        if a_dist > close {
                            close = a_dist;
                            let pm = s.value(0.5 * (uf + ul), v);
                            mid = p1.square_distance(&pm);
                        } else {
                            distmin = distmin.min(a_dist);
                        }
                    }
                    // `cxx:763-764`.
                    delt = delt.min(0.5 * occt_geom::approx_same_parameter::u_resolution(s, distmin.sqrt()));
                } else {
                    // `cxx:767-787`: the boundary pole rows.
                    let nbvp = bs.nb_poles_v();
                    close = bs.poles[0][0].square_distance(&bs.poles[nbup - 1][0]);
                    mid = bs.poles[0][0].square_distance(&bs.poles[nbup / 2][0]);
                    distmin = close;
                    for i in 1..nbvp {
                        let a_dist =
                            bs.poles[0][i].square_distance(&bs.poles[nbup - 1][i]);
                        if a_dist > close {
                            close = a_dist;
                            mid = bs.poles[0][i].square_distance(&bs.poles[nbup / 2][i]);
                        } else {
                            distmin = distmin.min(a_dist);
                        }
                    }
                    // `cxx:785-786`.
                    delt = delt.min(0.5 * occt_geom::approx_same_parameter::u_resolution(s, distmin.sqrt()));
                }
            }
        }
        // `cxx:790-820`.
        TypeKind::BezierSurface => {
            if let Some(bz) = s.osculating_bspline() {
                let nbup = bz.poles.len();
                if nbup >= 3 {
                    let nbvp = bz.nb_poles_v();
                    close = bz.poles[0][0].square_distance(&bz.poles[nbup - 1][0]);
                    mid = bz.poles[0][0].square_distance(&bz.poles[nbup / 2][0]);
                    let mut distmin = close;
                    for i in 0..nbvp {
                        let a_dist =
                            bz.poles[0][i].square_distance(&bz.poles[nbup - 1][i]);
                        if a_dist > close {
                            close = a_dist;
                            mid = bz.poles[0][i].square_distance(&bz.poles[nbup / 2][i]);
                        } else {
                            distmin = distmin.min(a_dist);
                        }
                    }
                    // `cxx:817-818`.
                    delt = delt.min(0.5 * occt_geom::approx_same_parameter::u_resolution(s, distmin.sqrt()));
                }
            }
        }
        // `cxx:822-851`: `Geom_RectangularTrimmedSurface` and
        // `Geom_OffsetSurface`, sampled at 101 V parameters.
        _ => {
            let nbpoints = 101i32;
            let p1 = s.value(uf, vf);
            let p2 = s.value(ul, vf);
            close = p1.square_distance(&p2);
            let pm = s.value(0.5 * (uf + ul), vf);
            mid = p1.square_distance(&pm);
            let mut distmin = close;
            for i in 1..nbpoints {
                let vparam = vf + (vl - vf) * i as f64 / (nbpoints - 1) as f64;
                let p1 = s.value(uf, vparam);
                let p2 = s.value(ul, vparam);
                let a_dist = p1.square_distance(&p2);
                if a_dist > close {
                    close = a_dist;
                    let pm = s.value(0.5 * (uf + ul), vparam);
                    mid = p1.square_distance(&pm);
                } else {
                    distmin = distmin.min(a_dist);
                }
            }
            // `cxx:848-849`.
            delt = delt.min(0.5 * occt_geom::approx_same_parameter::u_resolution(s, distmin.sqrt()));
        }
    }
    // `cxx:853-854`.
    close = close.sqrt();
    // `cxx:857-861`: the mid-line symmetry test discards a "closed" verdict
    // when the boundary gap exceeds the distance to the mid-U line.
    if mid > 0.0 && close > mid.sqrt() {
        return (f64::MAX, -1.0, delt);
    }
    (close, mid, delt)
}

/// `ShapeAnalysis_Surface::IsVClosed` (`ShapeAnalysis_Surface.cxx:865-1077`).
/// Same contract as [`sa_u_close_state`] with the V switch; `myVDelt` is the
/// `dv` counterpart of `myUDelt` (`cxx:879`).
fn sa_v_close_state(s: &dyn Surface) -> (f64, f64, f64) {
    use crate::geom_bnd_lib_surface3d::SurfaceKind as TypeKind;
    let (uf, ul) = s.u_range();
    let (mut vf, mut vl) = s.v_range();
    restrict_bounds(&mut vf, &mut vl);
    // `cxx:879`.
    let mut delt = (vl - vf).abs() / 20.0;
    // `cxx:881-888`.
    if geom_v_closed(s) {
        return (0.0, -1.0, 0.0);
    }
    // `cxx:891-896`.
    let surftype = if s.rectangular_trimmed_basis().is_some() {
        TypeKind::OtherSurface
    } else {
        crate::geom_bnd_lib_surface3d::surface_kind(s)
    };
    let mut close = f64::MAX;
    let mut mid = -1.0;
    match surftype {
        // `cxx:900-906`.
        TypeKind::Plane
        | TypeKind::Cone
        | TypeKind::Cylinder
        | TypeKind::Sphere
        | TypeKind::SurfaceOfExtrusion => {}
        // `cxx:908-916`.
        TypeKind::SurfaceOfRevolution => {
            if let Some(crv) = s.revolution_basis_curve() {
                let p1 = crv.d0(crv.first_parameter());
                let p2 = crv.d0(crv.last_parameter());
                close = p1.square_distance(&p2);
            }
        }
        // `cxx:917-984`.
        TypeKind::BSplineSurface => {
            if let Some(bs) = s.osculating_bspline() {
                let nbvp = bs.nb_poles_v();
                let mut distmin = f64::MAX;
                let (m1, m2) = end_mults(&bs.knots_v);
                if bs.v_periodic {
                    close = 0.0;
                    delt = 0.0;
                } else if nbvp < 3 {
                    // `cxx:928-931`.
                } else if bs.is_rational()
                    || m1 != bs.deg_v as i32 + 1
                    || m2 != bs.deg_v as i32 + 1
                {
                    // `cxx:932-960`.
                    let (uknots, _) = occt_geom::bspline_surface::GeomBSplineSurface::unique_knots_mults(&bs.knots_u);
                    let nbuk = uknots.len();
                    let mut u = uknots[0];
                    let p1 = s.value(u, vf);
                    let p2 = s.value(u, vl);
                    close = p1.square_distance(&p2);
                    let pm = s.value(u, 0.5 * (vf + vl));
                    mid = p1.square_distance(&pm);
                    distmin = close;
                    for i in 1..nbuk {
                        u = 0.5 * (uknots[i - 1] + uknots[i]);
                        let p1 = s.value(u, vf);
                        let p2 = s.value(u, vl);
                        let a_dist = p1.square_distance(&p2);
                        if a_dist > close {
                            close = a_dist;
                            let pm = s.value(u, 0.5 * (vf + vl));
                            mid = p1.square_distance(&pm);
                        } else {
                            distmin = distmin.min(a_dist);
                        }
                    }
                    // `cxx:958-959`.
                    delt = delt.min(0.5 * occt_geom::approx_same_parameter::v_resolution(s, distmin.sqrt()));
                } else {
                    // `cxx:962-983`: the boundary pole columns.
                    let nbup = bs.poles.len();
                    close = bs.poles[0][0].square_distance(&bs.poles[0][nbvp - 1]);
                    mid = bs.poles[0][0].square_distance(&bs.poles[0][nbvp / 2]);
                    distmin = close;
                    for i in 1..nbup {
                        let a_dist = bs.poles[i][0].square_distance(&bs.poles[i][nbvp - 1]);
                        if a_dist > close {
                            close = a_dist;
                            mid = bs.poles[i][0].square_distance(&bs.poles[i][nbvp / 2]);
                        } else {
                            distmin = distmin.min(a_dist);
                        }
                    }
                    // `cxx:980-981`.
                    delt = delt.min(0.5 * occt_geom::approx_same_parameter::v_resolution(s, distmin.sqrt()));
                }
            }
        }
        // `cxx:985-1015`.
        TypeKind::BezierSurface => {
            if let Some(bz) = s.osculating_bspline() {
                let nbvp = bz.nb_poles_v();
                if nbvp >= 3 {
                    let nbup = bz.nb_poles_u();
                    close = bz.poles[0][0].square_distance(&bz.poles[0][nbvp - 1]);
                    mid = bz.poles[0][0].square_distance(&bz.poles[0][nbvp / 2]);
                    let mut distmin = close;
                    for i in 1..nbup {
                        let a_dist = bz.poles[i][0].square_distance(&bz.poles[i][nbvp - 1]);
                        if a_dist > close {
                            close = a_dist;
                            mid = bz.poles[i][0].square_distance(&bz.poles[i][nbvp / 2]);
                        } else {
                            distmin = distmin.min(a_dist);
                        }
                    }
                    // `cxx:1012-1013`.
                    delt = delt.min(0.5 * occt_geom::approx_same_parameter::v_resolution(s, distmin.sqrt()));
                }
            }
        }
        // `cxx:1017-1046`.
        _ => {
            let nbpoints = 101i32;
            let p1 = s.value(uf, vf);
            let p2 = s.value(uf, vl);
            let pm = s.value(uf, 0.5 * (vf + vl));
            close = p1.square_distance(&p2);
            mid = p1.square_distance(&pm);
            let mut distmin = close;
            for i in 1..nbpoints {
                let uparam = uf + (ul - uf) * i as f64 / (nbpoints - 1) as f64;
                let p1 = s.value(uparam, vf);
                let p2 = s.value(uparam, vl);
                let a_dist = p1.square_distance(&p2);
                if a_dist > close {
                    close = a_dist;
                    let pm = s.value(uparam, 0.5 * (vf + vl));
                    mid = p1.square_distance(&pm);
                } else {
                    distmin = distmin.min(a_dist);
                }
            }
            // `cxx:1043-1044`.
            delt = delt.min(0.5 * occt_geom::approx_same_parameter::v_resolution(s, distmin.sqrt()));
        }
    }
    // `cxx:1048-1049`.
    close = close.sqrt();
    // `cxx:1052-1056`.
    if mid > 0.0 && close > mid.sqrt() {
        return (f64::MAX, -1.0, delt);
    }
    (close, mid, delt)
}

/// `ShapeAnalysis_Surface::IsUClosed(preci)` (`ShapeAnalysis_Surface.cxx:660-863`):
/// `myUCloseVal <= max(preci, Precision::Confusion())` (`cxx:658`, `:863`).
pub(crate) fn sa_is_u_closed(s: &dyn Surface, preci: f64) -> bool {
    sa_u_close_state(s).0 <= preci.max(occt_core::precision::CONFUSION)
}

/// `ShapeAnalysis_Surface::IsVClosed(preci)` (`ShapeAnalysis_Surface.cxx:865-1077`).
pub(crate) fn sa_is_v_closed(s: &dyn Surface, preci: f64) -> bool {
    sa_v_close_state(s).0 <= preci.max(occt_core::precision::CONFUSION)
}

/// `myUDelt` as `ShapeAnalysis_Surface::ValueOfUV` leaves it before reading
/// `du = min(myUDelt, SurfAdapt.UResolution(preci))` (`cxx:1332-1346`): the
/// first `ValueOfUV` of a general-arm surface forces `IsUClosed()` when
/// `myUCloseVal < 0`, so the field always carries the state computed here.
fn sa_u_delt(s: &dyn Surface) -> f64 {
    sa_u_close_state(s).2
}

/// `myVDelt`; see [`sa_u_delt`] and `cxx:1336-1347`.
fn sa_v_delt(s: &dyn Surface) -> f64 {
    sa_v_close_state(s).2
}

/// `ShapeAnalysis_Curve::Project(C3D, P3D, preci, proj, param, cf, cl, false)`
/// (`ShapeAnalysis_Curve.cxx:147-201`): an endpoint within
/// `Precision::Confusion()` is returned directly (`cxx:164-182`; the tolerance
/// is `Confusion` and not `preci` precisely because `AdjustToEnds` is false),
/// otherwise the range of a non-closed curve is widened by
/// `min(Resolution(preci), 0.1 * (umax - umin))` (`cxx:192-197`) and the minimum
/// over it is `ProjectAct` (`cxx:265-497`), the `Extrema_ExtPC` scan.
/// `ProjectAct(Adaptor3d_Curve, ..)` (`cxx:205-261`, taken by the offset arm at
/// `cxx:1806`) computes the same minimum, so one helper serves both arms and
/// `Adaptor3d_IsoCurve` is not needed.
///
/// Returns `(other, dist)` - OCCT's `param` / `distmin`.
///
/// UNPORTED (reported): `Geom_Curve::IsClosed()` is approximated by
/// `Curve::is_periodic()` plus coincident ends, because the trait exposes no
/// `IsClosed` (`Geom_Curve.hxx:155` is pure virtual and each curve implements
/// it, e.g. `Geom_BSplineCurve` through its first/last pole).
fn project_on_curve_range(
    iso: &dyn Curve,
    p: &GpPnt,
    preci: f64,
    cf: f64,
    cl: f64,
) -> Option<(f64, f64)> {
    let (u_min, u_max) = if cf < cl { (cf, cl) } else { (cl, cf) };
    let low = iso.d0(u_min);
    let d_low = low.distance(p);
    if d_low <= occt_core::precision::CONFUSION {
        return Some((u_min, d_low));
    }
    let high = iso.d0(u_max);
    let d_high = high.distance(p);
    if d_high <= occt_core::precision::CONFUSION {
        return Some((u_max, d_high));
    }
    let (mut a, mut b) = (u_min, u_max);
    let closed = iso.is_periodic() || low.distance(&high) <= occt_core::precision::CONFUSION;
    if !closed {
        let delta = iso.resolution(preci).min((u_max - u_min) * 0.1);
        a -= delta;
        b += delta;
    }
    let (t, _pt, d) = occt_geom::extrema_pc::extrema_ext_pc_min_in_range(iso, p, a, b)?;
    Some((t, d))
}

/// `ShapeAnalysis_Curve::NextProject` (`ShapeAnalysis_Curve.cxx:504-558`, the
/// `Geom_Curve` overload): the Newton walk `Extrema_LocateExtPC` from
/// `prev` (`cxx:572`), falling back to [`project_on_curve_range`] (OCCT
/// `Project(Adaptor3d_Curve,...)`, `cxx:205-261`) when the walk reports
/// `!IsDone`.
///
/// `AdjustToEnds` is false at the only call site (`cxx:2759`), so the bounded
/// endpoint snap uses `Precision::Confusion()` (`cxx:520-538`). As in OCCT,
/// a non-closed iso has its range widened by
/// `min(Resolution(preci), 0.1 * (umax - umin))` (`cxx:551-555`).
///
/// Returns `(param, dist)` - OCCT's `param` / the returned distance.
fn next_project_on_curve_range(
    iso: &dyn Curve,
    p: &GpPnt,
    preci: f64,
    prev: f64,
    cf: f64,
    cl: f64,
) -> Option<(f64, f64)> {
    let (mut u_min, mut u_max) = if cf < cl { (cf, cl) } else { (cl, cf) };
    if iso.first_parameter().is_finite() && iso.last_parameter().is_finite() {
        let low = iso.d0(u_min);
        let d_low = low.distance(p);
        if d_low <= occt_core::precision::CONFUSION {
            return Some((u_min, d_low));
        }
        let high = iso.d0(u_max);
        let d_high = high.distance(p);
        if d_high <= occt_core::precision::CONFUSION {
            return Some((u_max, d_high));
        }
    }
    let low = iso.d0(u_min);
    let high = iso.d0(u_max);
    let closed = iso.is_periodic() || low.distance(&high) <= occt_core::precision::CONFUSION;
    if !closed {
        let delta = iso.resolution(preci).min((u_max - u_min) * 0.1);
        u_min -= delta;
        u_max += delta;
    }
    // `Extrema_LocateExtPC(P3D, GAC, paramPrev, uMin, uMax, preci)` (`cxx:572`).
    if let Some((t, q)) =
        crate::int_tools_vertex_line::extrema_locate_ext_pc(iso, p, prev, u_min, u_max)
    {
        return Some((t, p.distance(&q)));
    }
    // `return Project(C3D, P3D, preci, proj, param, false)` (`cxx:578`), the
    // `Adaptor3d_Curve` overload on the already-widened `GAC` range.
    let d_low = iso.d0(u_min).distance(p);
    if d_low <= occt_core::precision::CONFUSION {
        return Some((u_min, d_low));
    }
    let d_high = iso.d0(u_max).distance(p);
    if d_high <= occt_core::precision::CONFUSION {
        return Some((u_max, d_high));
    }
    let (t, _pt, d) = occt_geom::extrema_pc::extrema_ext_pc_min_in_range(iso, p, u_min, u_max)?;
    if d < d_low + occt_core::precision::CONFUSION
        && d < d_high + occt_core::precision::CONFUSION
    {
        return Some((t, d));
    }
    if d_low < d_high {
        Some((u_min, d_low))
    } else {
        Some((u_max, d_high))
    }
}

/// `ShapeAnalysis_Surface::GetBoxUF/UL/VF/VL` (`ShapeAnalysis_Surface.cxx:1924-1946`)
/// reached through `ComputeBoxes` (`cxx:1897-1922`): the box of a bound iso
/// curve, `BndLib_Add3dCurve::Add(GeomAdaptor_Curve(iso), Precision::Confusion(),
/// B)` (`cxx:1904-1921`).
///
/// `BndLib_Add3dCurve::Add` (`BndLib_Add3dCurve.cxx:22-25`) uses the adaptor's
/// first/last parameter, i.e. the iso curve's own range; `box_curve` unwraps a
/// `Geom_TrimmedCurve` to its basis exactly like `GeomAdaptor_Curve::load`
/// (`GeomAdaptor_Curve.cxx:252-254`). The returned box carries `Gap = Confusion`
/// (`Bnd_Box.hxx:154` `Enlarge`), which is what `Bnd_Box::IsOut` tests
/// (`Bnd_Box.cxx`, `box3d.rs`).
fn iso_bound_box(iso: &dyn Curve) -> occt_core::bnd::BndBox {
    crate::geom_bnd_lib_curve3d::box_curve(
        iso,
        iso.first_parameter(),
        iso.last_parameter(),
        occt_core::precision::CONFUSION,
    )
}

/// Outputs of [`is_an_isoparametric`], the in-out reference parameters of
/// `ShapeConstruct_ProjectCurveOnSurface::isAnIsoparametric`
/// (`ShapeConstruct_ProjectCurveOnSurface.cxx:2461-2477`).
pub(super) struct IsoParamOut {
    /// `theIsTypeU`: the selected iso is a `UIso` (fixed `U`).
    pub is_type_u: bool,
    /// `theP1OnIso` / `theValueP1`: the first sample lies on some bound iso.
    pub p1_on_iso: bool,
    pub value_p1: GpPnt2d,
    /// `theP2OnIso` / `theValueP2`: same for the last sample.
    pub p2_on_iso: bool,
    pub value_p2: GpPnt2d,
    /// `theIsoPar2d3d`: the 3D sample parameters equal the iso parameters, so
    /// no 2D projection is needed.
    pub iso_par2d3d: bool,
    /// `theCIso`: the selected bound iso curve.
    pub c_iso: Option<Arc<dyn Curve>>,
    /// `theT1` / `theT2`: the iso curve range corresponding to the surface
    /// window (`V1,V2` for a `UIso`, `U1,U2` for a `VIso`).
    pub t1: f64,
    pub t2: f64,
    /// `theParamsOut`: the iso parameters of the interior samples found by the
    /// `NextProject` distance walk (`cxx:2757-2768`).
    pub params_out: Vec<f64>,
}

impl IsoParamOut {
    fn new(n: usize) -> Self {
        Self {
            is_type_u: false,
            p1_on_iso: false,
            value_p1: GpPnt2d::zero(),
            p2_on_iso: false,
            value_p2: GpPnt2d::zero(),
            iso_par2d3d: false,
            c_iso: None,
            t1: 0.0,
            t2: 0.0,
            params_out: vec![0.0; n],
        }
    }
}

/// `ShapeConstruct_ProjectCurveOnSurface::isAnIsoparametric`
/// (`ShapeConstruct_ProjectCurveOnSurface.cxx:2461-2796`): test whether the
/// sampled 3D curve runs along one of the four boundary isos of the surface
/// parameter window, and, if so, the exact UV values that make its pcurve a
/// straight `(u, v)` line.
///
/// The four candidates are `UIso(U1)` / `UIso(U2)` / `VIso(V1)` / `VIso(V2)`
/// (`cxx:2512-2556`, `mySurf->Bounds` at `cxx:2485`; a
/// `Geom_RectangularTrimmedSurface` re-reads its own bounds at
/// `cxx:2487-2492`). For each, the first and last samples are matched against
/// the iso ends (`cxx:2588-2610`, modes 1/2) or projected on the iso by
/// `ShapeAnalysis_Curve::Project` (mode 3, `cxx:2631-2641`), with the bound
/// box `Bnd_Box::IsOut` rejection at `cxx:2619` and the spherical-`VIso` skip
/// at `cxx:2614`. The best pair (smallest sum of squared distances,
/// `cxx:2682-2696`) is kept.
///
/// When both ends land on the iso the samples are validated directly against
/// the iso 3D points (`cxx:2714-2731`, `theIsoPar2d3d`) or, failing that, by
/// the `NextProject` distance walk (`cxx:2741-2769`).
///
/// `mySurf->Surface()->Bounds(uf, ul, vf, vl)` (`cxx:1160-1161`) and
/// `mySurf->Bounds(U1, U2, V1, V2)` (`cxx:2485`) are both the port's
/// `Surface::u_range` / `v_range` (`ShapeAnalysis_Surface::Bounds`,
/// `ShapeAnalysis_Surface.lxx:54-62`, returns `myUF/myUL/myVF/myVL`). Those
/// members are copied from `mySurf->Bounds` by both constructors
/// (`ShapeAnalysis_Surface.cxx:118`, `:142`) and are written nowhere else:
/// `ShapeAnalysis_Surface::SetDomain` (`cxx:1887-1895`) has no caller anywhere
/// under `src/`, and every `RestrictBounds` call (`cxx:671`, `:878`, `:1323`,
/// `:1616`, `:1711`, `:1741`, `:1779`, `:1806`) takes local copies of the
/// members. The two bounds are therefore equal on every face, `Shape-1` (60
/// faces, 16 spherical) included.
pub(super) fn is_an_isoparametric(
    surf: &dyn Surface,
    points: &[GpPnt],
    params: &[f64],
    out: &mut IsoParamOut,
) -> bool {
    // `cxx:2479`.
    let prec = occt_core::precision::CONFUSION;
    let nb_pnt = points.len();
    let mut iso_param = false;
    out.iso_par2d3d = false;
    let (u1, u2) = surf.u_range();
    let (v1, v2) = surf.v_range();
    // `cxx:2497-2502`.
    let mut mpt = [0i32; 2];
    let mut tpar = [0.0f64; 2];
    let mut iso_value = 0.0f64;
    let mut mindist2 = 4.0 * prec * prec;
    let mut mind2 = [mindist2; 2];
    out.p1_on_iso = false;
    out.p2_on_iso = false;
    // `cxx:2614`: `IsKind(STANDARD_TYPE(Geom_SphericalSurface))`, so a
    // `Geom_RectangularTrimmedSurface` wrapping a sphere does not qualify.
    let is_sphere = surf.gp_sphere().is_some();
    let mut chosen: Option<Arc<dyn Curve>> = None;

    for j in 0..4usize {
        // `cxx:2512-2556`.
        let (iso_u, iso_val, bound, tt1, tt2) = match j {
            0 => {
                if occt_core::precision::Precision::is_infinite(u1) {
                    continue;
                }
                (true, u1, u1, v1, v2)
            }
            1 => {
                if occt_core::precision::Precision::is_infinite(u2) {
                    continue;
                }
                (true, u2, u2, v1, v2)
            }
            2 => {
                if occt_core::precision::Precision::is_infinite(v1) {
                    continue;
                }
                (false, v1, v1, u1, u2)
            }
            _ => {
                if occt_core::precision::Precision::is_infinite(v2) {
                    continue;
                }
                (false, v2, v2, u1, u2)
            }
        };
        let ci = if iso_u {
            surf.u_iso_curve(bound)
        } else {
            surf.v_iso_curve(bound)
        };
        // `cxx:2558-2561`.
        let Some(ci) = ci else {
            continue;
        };
        let a_box = iso_bound_box(ci.as_ref());
        // `cxx:2578-2583`: a degenerate iso (all three probe points equal) is
        // skipped.
        let ext1 = ci.d0(tt1);
        let ext2 = ci.d0(tt2);
        let extmi = ci.d0(0.5 * (tt1 + tt2));
        if ext1.distance(&ext2) <= prec && ext1.distance(&extmi) <= prec {
            continue;
        }
        let mut pt_eq_ext1 = false;
        let mut pt_eq_ext2 = false;
        let mut currd2 = [0.0f64; 2];
        let mut tp = [0.0f64; 2];
        let mut mp = [0i32; 2];
        for i in 0..2usize {
            mp[i] = 0;
            // OCCT `k = (i == 0 ? 1 : theNbPnt)` (1-based) (`cxx:2594`).
            let k = if i == 0 { 0 } else { nb_pnt - 1 };
            currd2[i] = points[k].square_distance(&ext1);
            if currd2[i] <= prec * prec && !pt_eq_ext1 {
                mp[i] = 1;
                tp[i] = tt1;
                pt_eq_ext1 = true;
                continue;
            }
            currd2[i] = points[k].square_distance(&ext2);
            if currd2[i] <= prec * prec && !pt_eq_ext2 {
                mp[i] = 2;
                tp[i] = tt2;
                pt_eq_ext2 = true;
                continue;
            }
            // `cxx:2614-2617`: a spherical surface has no useful `VIso`
            // projection for its ends.
            if is_sphere && !iso_u {
                continue;
            }
            // `cxx:2619-2622`.
            if a_box.is_out(&points[k]) {
                continue;
            }
            let cf0 = ci.first_parameter();
            let cl0 = ci.last_parameter();
            let cf = if occt_core::precision::Precision::is_infinite(cf0) {
                -1000.0
            } else {
                cf0
            };
            let cl = if occt_core::precision::Precision::is_infinite(cl0) {
                1000.0
            } else {
                cl0
            };
            // `cxx:2635`: `sac.Project(cI, thePoints(k), prec, pt, t, Cf, Cl)`
            // with the default `AdjustToEnds = true`; since `prec` equals the
            // projection precision here the endpoint snap tolerance is the same
            // as the port's `AdjustToEnds = false` helper.
            if let Some((t, dist)) = project_on_curve_range(ci.as_ref(), &points[k], prec, cf, cl) {
                currd2[i] = dist * dist;
                if dist <= prec && t >= cf && t <= cl {
                    mp[i] = 3;
                    tp[i] = t;
                }
            }
        }
        // `cxx:2643-2647`.
        if mp[0] > 0 && mp[1] > 0 && (tp[0] - tp[1]).abs() < occt_core::precision::PCONFUSION {
            continue;
        }
        // `cxx:2649-2673`.
        if mp[0] > 0 && (!out.p1_on_iso || currd2[0] < mind2[0]) {
            out.p1_on_iso = true;
            mind2[0] = currd2[0];
            out.value_p1 = if iso_u {
                GpPnt2d::new(iso_val, tp[0])
            } else {
                GpPnt2d::new(tp[0], iso_val)
            };
        }
        if mp[1] > 0 && (!out.p2_on_iso || currd2[1] < mind2[1]) {
            out.p2_on_iso = true;
            mind2[1] = currd2[1];
            out.value_p2 = if iso_u {
                GpPnt2d::new(iso_val, tp[1])
            } else {
                GpPnt2d::new(tp[1], iso_val)
            };
        }
        // `cxx:2675-2696`.
        if mp[0] <= 0 || mp[1] <= 0 {
            continue;
        }
        let md2 = currd2[0] + currd2[1];
        if mindist2 <= md2 {
            continue;
        }
        mindist2 = md2;
        mpt[0] = mp[0];
        mpt[1] = mp[1];
        tpar[0] = tp[0];
        tpar[1] = tp[1];
        out.is_type_u = iso_u;
        iso_value = iso_val;
        out.t1 = tt1;
        out.t2 = tt2;
        chosen = Some(ci.clone());
    }

    // `cxx:2700-2771`.
    if mpt[0] > 0 && mpt[1] > 0 {
        let Some(ci) = chosen else {
            return false;
        };
        out.p1_on_iso = true;
        out.p2_on_iso = true;
        if out.is_type_u {
            out.value_p1 = GpPnt2d::new(iso_value, tpar[0]);
            out.value_p2 = GpPnt2d::new(iso_value, tpar[1]);
        } else {
            out.value_p1 = GpPnt2d::new(tpar[0], iso_value);
            out.value_p2 = GpPnt2d::new(tpar[1], iso_value);
        }
        if mpt[0] != 3 && mpt[1] != 3 {
            out.iso_par2d3d = true;
            for i in 1..nb_pnt.saturating_sub(1) {
                let t = if tpar[1] > tpar[0] {
                    params[i]
                } else {
                    out.t1 + out.t2 - params[i]
                };
                let pt = ci.d0(t);
                if points[i].distance(&pt) > prec {
                    out.iso_par2d3d = false;
                    break;
                }
            }
        }
        if out.iso_par2d3d {
            iso_param = true;
        } else {
            let mut prev_param = tpar[0];
            let cf0 = ci.first_parameter();
            let cl0 = ci.last_parameter();
            let cf = if occt_core::precision::Precision::is_infinite(cf0) {
                -1000.0
            } else {
                cf0
            };
            let cl = if occt_core::precision::Precision::is_infinite(cl0) {
                1000.0
            } else {
                cl0
            };
            let mut iso_by_distance = true;
            for i in 1..nb_pnt.saturating_sub(1) {
                match next_project_on_curve_range(ci.as_ref(), &points[i], prec, prev_param, cf, cl)
                {
                    Some((t, dist)) => {
                        prev_param = t;
                        out.params_out[i] = t;
                        if dist > prec || t < cf || t > cl {
                            iso_by_distance = false;
                            break;
                        }
                    }
                    None => {
                        iso_by_distance = false;
                        break;
                    }
                }
            }
            if iso_by_distance {
                iso_param = true;
            }
        }
        out.c_iso = Some(ci);
    }
    iso_param
}

/// `ShapeAnalysis_Surface::UVFromIso` (`ShapeAnalysis_Surface.cxx:1522-1842`):
/// the projection restricted to the two isos through the seed plus the four
/// window borders - candidates 0/1/2 are `UIso(myUF)` / `UIso(myUL)` / `UIso(U)`
/// and 3/4/5 `VIso(myVF)` / `VIso(myVL)` / `VIso(V)` (`cxx:1563-1603`) - then
/// refined by up to `MaxIters = 5` passes alternating the iso through the
/// running solution (`cxx:1688-1818`).
///
/// This is the "on the border" arm of `ValueOfUV`: when the extrema minimum
/// keeps a residue, the closest point of the face is often on one of the window
/// isos rather than at an interior stationary point.
///
/// Returns `(theMin, U, V)`; `u` / `v` are the seeds, matching OCCT where
/// `UU` / `VV` are in-out references and `theMin` is the returned residue.
fn uv_from_iso(s: &dyn Surface, p: &GpPnt, preci: f64, u: f64, v: f64) -> (f64, f64, f64) {
    // `cxx:1528-1537`: the seed itself is the first candidate.
    let mut uu = u;
    let mut vv = v;
    let mut the_min = s.d0(u, v).distance(p);
    if the_min < preci / 10.0 {
        // `cxx:1539-1542` "c etait deja OK"
        return (the_min, uu, vv);
    }
    // `cxx:1535-1538` reads `myUF` / `myUL` (the raw `Bounds`, not the
    // restricted window) and skips an infinite one (`cxx:1605`).
    let (uf, ul) = s.u_range();
    let (vf, vl) = s.v_range();
    let mut uv = true;
    for num in 0..6 {
        // `cxx:1563` `UV = (num < 3)`: an iso-U candidate fixes `par` to the U
        // of the iso and takes the projected parameter as V, and vice versa.
        uv = num < 3;
        let (par, iso) = match num {
            0 => (uf, s.u_iso_curve(uf)),
            1 => (ul, s.u_iso_curve(ul)),
            2 => (u, s.u_iso_curve(u)),
            3 => (vf, s.v_iso_curve(vf)),
            4 => (vl, s.v_iso_curve(vl)),
            _ => (v, s.v_iso_curve(v)),
        };
        if occt_core::precision::Precision::is_infinite(par) {
            continue;
        }
        // `cxx:1603` `iso.IsNull()`: no iso through that parameter.
        let Some(iso) = iso else {
            continue;
        };
        let (mut cf, mut cl) = (iso.first_parameter(), iso.last_parameter());
        restrict_bounds(&mut cf, &mut cl);
        if let Some((other, dist)) = project_on_curve_range(iso.as_ref(), p, preci, cf, cl) {
            if dist < the_min {
                the_min = dist;
                uu = if uv { par } else { other };
                vv = if uv { other } else { par };
            }
        }
    }
    // `cxx:1688-1818`: two projections per pass with opposite parities; OCCT
    // enters with `UV = false` (the last `UV = (num < 3)` has `num = 5`).
    let (mut prev_u, mut prev_v) = (u, v);
    let mut iters = 0;
    while (prev_u != uu || prev_v != vv) && iters < 5 && the_min > preci {
        prev_u = uu;
        prev_v = vv;
        for _ in 0..2 {
            let iso = if uv { s.u_iso_curve(uu) } else { s.v_iso_curve(vv) };
            if let Some(iso) = iso {
                let (mut cf, mut cl) = (iso.first_parameter(), iso.last_parameter());
                restrict_bounds(&mut cf, &mut cl);
                if let Some((other, dist)) = project_on_curve_range(iso.as_ref(), p, preci, cf, cl) {
                    if dist < the_min {
                        the_min = dist;
                        if uv {
                            vv = other;
                        } else {
                            uu = other;
                        }
                    }
                }
            }
            uv = !uv;
        }
        iters += 1;
    }
    (the_min, uu, vv)
}

/// `theIndCoord` accessor of `insertAdditionalPointOrAdjust` (`cxx:2039`):
/// 1 is U (the `x` coordinate), 2 is V (`thePoints2d.SetCoord(theIndCoord, ..)`).
fn iso_coord_of(p: &GpPnt2d, idx: usize) -> f64 {
    if idx == 1 { p.x() } else { p.y() }
}

fn set_iso_coord_of(p: &mut GpPnt2d, idx: usize, v: f64) {
    if idx == 1 { p.set_x(v) } else { p.set_y(v); }
}

/// `ShapeConstruct_ProjectCurveOnSurface::insertAdditionalPointOrAdjust`
/// (`ShapeConstruct_ProjectCurveOnSurface.cxx:2039-2147`): a period jump
/// between `thePrevCoord` and `theCurCoord` is either bridged by the projected
/// 3D midpoint (when the curve really crosses the seam) or removed by shifting
/// the current coordinate one period back onto the previous one.
/// `index` is 0-based here (OCCT's `theIndex` is 1-based); it is advanced past
/// the inserted midpoint. `params` / `pts2d` / `pts3d` grow together, as in
/// OCCT (`cxx:2111-2131` inserts the 3D midpoint in `thePoints` as well).
#[allow(clippy::too_many_arguments)]
fn insert_additional_point_or_adjust(
    to_adjust: &mut bool,
    idx_coord: usize,
    period: f64,
    tol_on_period: f64,
    cur_coord_raw: f64,
    cur_coord: &mut f64,
    prev_coord: f64,
    preci: f64,
    curve: &dyn Curve,
    index: &mut usize,
    params: &mut Vec<f64>,
    pts3d: &mut Vec<GpPnt>,
    pts2d: &mut Vec<GpPnt2d>,
    surf: &dyn Surface,
) {
    // `cxx:2053-2054`
    let corrected = occt_core::bspl::locate::in_period(
        cur_coord_raw,
        prev_coord - 0.5 * period,
        prev_coord + 0.5 * period,
    );

    if !*to_adjust {
        // `cxx:2057-2064`
        let cur_par = params[*index];
        let prev_par = params[*index - 1];
        let mut mid_par = 0.5 * (prev_par + cur_par);
        let mut mid_p3d = curve.d0(mid_par);
        let mut mid_p2d = value_of_uv(surf, &mid_p3d, preci);
        let mut mid_coord = occt_core::bspl::locate::in_period(
            iso_coord_of(&mid_p2d, idx_coord),
            prev_coord - 0.5 * period,
            prev_coord + 0.5 * period,
        );
        let mut first_coord = prev_coord;
        let mut last_coord = corrected;
        if last_coord < first_coord {
            std::mem::swap(&mut first_coord, &mut last_coord);
        }
        if last_coord - first_coord <= tol_on_period {
            // `cxx:2074-2077`
            *to_adjust = true;
        } else if first_coord <= mid_coord && mid_coord <= last_coord {
            // `cxx:2078-2081`
            *to_adjust = true;
        } else {
            // `cxx:2082-2141`: bisect in the 3D parameter until the midpoint
            // projection leaves the half-period window around `thePrevCoord`.
            let mut success = true;
            let mut first_t = prev_par;
            let mut last_t = cur_par;
            mid_coord = iso_coord_of(&mid_p2d, idx_coord);
            while (mid_coord - prev_coord).abs() >= 0.5 * period - tol_on_period
                || (cur_coord_raw - mid_coord).abs() >= 0.5 * period - tol_on_period
            {
                if mid_par - first_t <= occt_core::precision::PCONFUSION
                    || last_t - mid_par <= occt_core::precision::PCONFUSION
                {
                    success = false;
                    break;
                }
                if (mid_coord - prev_coord).abs() >= 0.5 * period - tol_on_period {
                    last_t = 0.5 * (first_t + last_t);
                } else {
                    first_t = 0.5 * (first_t + last_t);
                }
                mid_par = 0.5 * (first_t + last_t);
                mid_p3d = curve.d0(mid_par);
                mid_p2d = value_of_uv(surf, &mid_p3d, preci);
                mid_coord = iso_coord_of(&mid_p2d, idx_coord);
            }
            if success {
                // `cxx:2110-2135`: insert the extra sample at `theIndex` in all
                // three arrays.
                params.insert(*index, mid_par);
                pts3d.insert(*index, mid_p3d);
                pts2d.insert(*index, mid_p2d);
                // `cxx:2135` `theIndex++`: in OCCT the `for` loop then increments
                // `aPntIter` again (`cxx:1542`), so the next iteration lands back
                // on the shifted original sample. The Rust caller increments `i`
                // itself, so no advance is made here; that reproduces the same
                // revisit - and, with it, the same `minX` / `maxX` accumulation
                // of the original sample's coordinate.
            } else {
                *to_adjust = true;
            }
        }
    }
    if *to_adjust {
        // `cxx:2142-2146`
        *cur_coord = corrected;
        set_iso_coord_of(&mut pts2d[*index], idx_coord, corrected);
    }
}

/// `approxPCurve` (`ShapeConstruct_ProjectCurveOnSurface.cxx:1499-1690`):
/// "Handle U-closed surfaces" / "Handle V-closed surfaces". A closed (not
/// necessarily periodic) surface's point-wise projection returns every sample
/// reduced into the parameter window, so a curve that runs once around the
/// closed direction carries a full-period jump in the middle of its `(u, v)`
/// sequence. Interpolating through that sequence produces a curve that sweeps
/// the whole window instead of running along the surface (measured on
/// ATU01038: 0.77 to 16.0 against a 0.01 projection tolerance), and the
/// mesher's UV frontier then degenerates.
///
/// `theIsFromCache` / `aSavedPoint` are produced by the sample loop
/// (`cxx:1295-1465`, ported in [`project_curve_on_surface_perform`]): in the
/// `!isRecompute` endpoint arm they are `getLine`'s `theIsFromCache` / the
/// `thePoints2d(1)` projection (`cxx:1390-1399`), and in the `isRecompute`
/// endpoint arm the `myCache` lookup sets them from the cache anchor
/// (`cxx:1401-1424`). `is_from_cache` switches off the "pull the first sample
/// into the window" and the final recentring steps (`cxx:1508`,
/// `cxx:1521-1536`, `cxx:1571`, `cxx:1599`, `cxx:1613-1628`, `cxx:1664`) so
/// that the cached parameter branch is preserved.
///
/// NOT PORTED HERE (reported): the `isoParam` arms of the sampling loop
/// (`cxx:1317-1375`), so this function's inputs are always the `else` arm of
/// `cxx:1376-1385`.
///
/// The `needResolveUJump` / `needResolveVJump` flags (`cxx:1436-1451`, latched
/// in the sampling loop of the caller) arrive as `need_resolve_u_jump` /
/// `need_resolve_v_jump` and widen the two triggers below exactly as in OCCT
/// (`cxx:1504`, `cxx:1596-1598`).
///
/// The `Handle AdjustOverDegen` block (`cxx:1746-1890`) runs at the end of this
/// function, after the V-closed seam adjustment, exactly as in OCCT.
fn resolve_closed_surface_period_jump(
    surf: &dyn Surface,
    curve: &dyn Curve,
    singularities: &[super::singularities::Singularity],
    preci: f64,
    need_resolve_u_jump: bool,
    need_resolve_v_jump: bool,
    is_from_cache: bool,
    saved_point: GpPnt2d,
    params: &mut Vec<f64>,
    pts3d: &mut Vec<GpPnt>,
    pts2d: &mut Vec<GpPnt2d>,
) {
    let (uf, ul) = surf.u_range();
    let (vf, vl) = surf.v_range();
    let up = ul - uf;
    let vp = vl - vf;
    // `cxx:1501-1502`
    let tol_on_u_period = occt_core::precision::CONFUSION * up;
    let tol_on_v_period = occt_core::precision::CONFUSION * vp;

    if sa_is_u_closed(surf, preci) || need_resolve_u_jump {
        // `cxx:1505-1519`: pull the first sample into `[uf, ul]` unless it came
        // from the cache, in which case that representative is kept.
        let mut first_x = pts2d[0].x();
        if !is_from_cache {
            while first_x < uf {
                first_x += up;
                pts2d[0].set_x(first_x);
            }
            while first_x > ul {
                first_x -= up;
                pts2d[0].set_x(first_x);
            }
        }
        // `cxx:1521-1536`: with a cached anchor the period window is placed
        // around it and the first sample is folded into that window.
        if surf.is_u_periodic() && is_from_cache {
            let mut a_min_param = uf;
            let mut a_max_param = ul;
            while a_min_param > saved_point.x() {
                a_min_param -= up;
                a_max_param -= up;
            }
            while a_max_param < saved_point.x() {
                a_min_param += up;
                a_max_param += up;
            }
            first_x += super::projection_cache::adjust_to_period(first_x, a_min_param, a_max_param);
            pts2d[0].set_x(first_x);
        }
        let mut prev_x = first_x;
        let mut min_x = first_x;
        let mut max_x = first_x;
        let mut to_adjust = false;
        // `cxx:1544-1569`
        let mut i = 1usize;
        while i < pts2d.len() {
            let raw_x = pts2d[i].x();
            let mut cur_x = raw_x;
            if (cur_x - prev_x).abs() > 0.5 * up {
                insert_additional_point_or_adjust(
                    &mut to_adjust,
                    1,
                    up,
                    tol_on_u_period,
                    raw_x,
                    &mut cur_x,
                    prev_x,
                    preci,
                    curve,
                    &mut i,
                    params,
                    pts3d,
                    pts2d,
                    surf,
                );
            }
            prev_x = cur_x;
            if min_x > cur_x {
                min_x = cur_x;
            } else if max_x < cur_x {
                max_x = cur_x;
            }
            i += 1;
        }
        // `cxx:1571-1591`: recentre the sequence on the surface window. Skipped
        // when the anchor came from the cache (`cxx:1571`).
        if !is_from_cache {
            let mid_x = 0.5 * (min_x + max_x);
            let shift_x = if mid_x > ul {
                -up
            } else if mid_x < uf {
                up
            } else {
                0.0
            };
            if shift_x != 0.0 {
                for p in pts2d.iter_mut() {
                    p.set_x(p.x() + shift_x);
                }
            }
        }
    }

    // `cxx:1596-1598`; the `Geom_SphericalSurface` arm is the extra trigger
    // OCCT adds because a sphere is closed in V without being V-periodic.
    if sa_is_v_closed(surf, preci) || need_resolve_v_jump || surf.gp_sphere().is_some() {
        // `cxx:1599-1611`
        let mut first_y = pts2d[0].y();
        if !is_from_cache {
            while first_y < vf {
                first_y += vp;
                pts2d[0].set_y(first_y);
            }
            while first_y > vl {
                first_y -= vp;
                pts2d[0].set_y(first_y);
            }
        }
        // `cxx:1613-1628`
        if surf.is_v_periodic() && is_from_cache {
            let mut a_min_param = vf;
            let mut a_max_param = vl;
            while a_min_param > saved_point.y() {
                a_min_param -= vp;
                a_max_param -= vp;
            }
            while a_max_param < saved_point.y() {
                a_min_param += vp;
                a_max_param += vp;
            }
            first_y += super::projection_cache::adjust_to_period(first_y, a_min_param, a_max_param);
            pts2d[0].set_y(first_y);
        }
        let mut prev_y = first_y;
        let mut min_y = first_y;
        let mut max_y = first_y;
        let mut to_adjust = false;
        // `cxx:1636-1662`
        let mut i = 1usize;
        while i < pts2d.len() {
            let raw_y = pts2d[i].y();
            let mut cur_y = raw_y;
            if (cur_y - prev_y).abs() > 0.5 * vp {
                insert_additional_point_or_adjust(
                    &mut to_adjust,
                    2,
                    vp,
                    tol_on_v_period,
                    raw_y,
                    &mut cur_y,
                    prev_y,
                    preci,
                    curve,
                    &mut i,
                    params,
                    pts3d,
                    pts2d,
                    surf,
                );
            }
            prev_y = cur_y;
            if min_y > cur_y {
                min_y = cur_y;
            } else if max_y < cur_y {
                max_y = cur_y;
            }
            i += 1;
        }
        // `cxx:1664-1685`. Skipped when the anchor came from the cache.
        if !is_from_cache {
            let mid_y = 0.5 * (min_y + max_y);
            let shift_y = if mid_y > vl {
                -vp
            } else if mid_y < vf {
                vp
            } else {
                0.0
            };
            if shift_y != 0.0 {
                for p in pts2d.iter_mut() {
                    p.set_y(p.y() + shift_y);
                }
            }
        }

        // `cxx:1687-1746`: "Handle V-closed seam adjustment". A pair that still
        // straddles the seam by more than half a period really crosses it, so
        // whichever of the two samples is nearest to a seam boundary is pinned
        // exactly onto it (`vf` or `vl`); leaving it a tolerance away from the
        // seam makes the mesher's UV frontier cut the surface instead of
        // closing along it.
        for i in 1..pts2d.len() {
            let prev_y = pts2d[i - 1].y();
            let cur_y = pts2d[i].y();
            if (cur_y - prev_y).abs() > 0.5 * vp {
                let dist_prev_vf = (prev_y - vf).abs();
                let dist_prev_vl = (prev_y - vl).abs();
                let dist_curr_vf = (cur_y - vf).abs();
                let dist_curr_vl = (cur_y - vl).abs();

                let mut prev_on_first = true;
                let mut prev_on_last = false;
                let mut curr_on_first = false;
                let mut curr_on_last = false;
                let mut the_min = dist_prev_vf;
                if dist_prev_vl < the_min {
                    the_min = dist_prev_vl;
                    prev_on_first = false;
                    prev_on_last = true;
                }
                if dist_curr_vf < the_min {
                    the_min = dist_curr_vf;
                    prev_on_first = false;
                    prev_on_last = false;
                    curr_on_first = true;
                }
                if dist_curr_vl < the_min {
                    prev_on_first = false;
                    prev_on_last = false;
                    curr_on_first = false;
                    curr_on_last = true;
                }

                if prev_on_first {
                    pts2d[i - 1].set_y(vf);
                } else if prev_on_last {
                    pts2d[i - 1].set_y(vl);
                } else if curr_on_first {
                    pts2d[i].set_y(vf);
                } else if curr_on_last {
                    pts2d[i].set_y(vl);
                }
            }
        }
    }

    // `cxx:1746-1890`: "Handle AdjustOverDegen".
    super::singularities::adjust_over_degenerated(surf, singularities, preci, pts3d, pts2d);
}

/// `ShapeConstruct_ProjectCurveOnSurface::Perform` (`cxx:707-774`):
/// `projectAnalytic` first (`cxx:731-736`), then `generateCurvePoints`
/// (`cxx:750-752`) into the `approxPCurve` chain (`cxx:760`). `preci` is
/// `myPreci` (`cxx:533` `myProjector->Init(sas, preci)`), which `approxPCurve`
/// forwards to `getLine` (`cxx:1119`). `tol_first` / `tol_last` are `Perform`'s
/// `theTolFirst` / `theTolLast`, the tolerances of the edge's end vertices
/// (`ShapeFix_Edge.cxx:521-531`); a negative value means `Precision::Confusion()`
/// (`cxx:1475-1476`).
///
/// The `myCache` endpoint cache and the B-spline corner cache are PORTED
/// (`pcurve_full/projection_cache.rs`): `getLine` seeds its four probes from them, anchors
/// `fixPeriodicityTroubles` on the cached point (`cxx:957-961`) and refills
/// `myCache` from the pcurve it returns (`cxx:1122-1142`); the fallback path
/// refills it at `cxx:1892-1903`. `cache` is `myCache`, cleared by `SetSurface`
/// (`cxx:676-683`) - the caller keeps one per reprojected edge chain.
///
/// The sample loop (`cxx:1311-1465`) is PORTED: interior samples walk from the
/// previous parameter point with `NextValueOfUV(p2d, p3d, myPreci,
/// Confusion + 1000 * gap)` (`cxx:1432`), the endpoints are re-projected from
/// `myCache` only when `getLine` fixed a period jump (`cxx:1401-1424`), and the
/// next seed is the extrapolation `2 * p2d - thePoints2d(prev)` (`cxx:1455-1463`).
///
/// PARK (unported, reported):
/// - `isBSplineCurveInvalid` / `PerformByProjLib` (`cxx:428`, `cxx:739-748`).
/// - the `isoParam` arms of the sample loop (`cxx:1317-1375`) and the
///   `p1OnIso` / `p2OnIso` endpoint overrides (`cxx:1379-1385`): both require
///   `isAnIsoparametric` (`cxx:1170`, impl `cxx:2461-2796`), an edge that lies
///   on a boundary iso of the surface parameter window and becomes an exact iso
///   pcurve (a straight line in `(u, v)`), with `theParamsOut` reparameterized
///   through `ShapeAnalysis_Curve::Project` / `NextProject`. Four things block
///   it:
///   1. `mySurf->UIso(U1)` / `UIso(U2)` / `VIso(V1)` / `VIso(V2)`
///      (`cxx:2519`, `:2530`, `:2540`, `:2549`). The trait's
///      `u_iso_curve` / `v_iso_curve` return `None` for the whole analytic
///      family (plane, cylinder, cone, sphere, torus), so those four arms would
///      silently never fire; implementing them means editing
///      `crates/occt-geom/src/**`, outside this port's allowed files.
///   2. `ShapeAnalysis_Surface::GetBoxUF/GetBoxUL/GetBoxVF/GetBoxVL`
///      (`cxx:2522`, `:2532`, `:2543`, `:2552`), used as the
///      `aBox->IsOut(thePoints(k))` rejection (`cxx:2619`).
///   3. `ShapeAnalysis_Curve::Project` (`ShapeAnalysis_Curve.cxx:126-202`,
///      `:205-261`, `ProjectAct` `:265-...`, extrema + Newton on a
///      `GeomAdaptor_Curve`) and `NextProject` (`:504-...`, `:561-...`); there
///      is no `ShapeAnalysis_Curve` module in this port.
///   4. `theCIso->D0` / `FirstParameter` / `LastParameter` are used on the
///      returned handle (`cxx:2578-2584`, `:2727`, `:2744`), which is available
///      through [`Curve`], but point 1 means the handle never exists.
///   Note that `theP1OnIso` / `theP2OnIso` can be set by `isAnIsoparametric`
///   even when it returns `isoParam == false` (`cxx:2647-2672`), so this port's
///   endpoints always take the `myCache` / `thePoints2d` arm. Measured on
///   ATU01038: the `getLine` arm above already collapses the iso-edge tolerance
///   band (top edge tolerance now equals OCCT's 0.007844513808 to 1e-13).
/// - `needResolveUJump` / `needResolveVJump` (`cxx:1436-1451`, consumed at
///   `cxx:1504` / `cxx:1596-1598`) are PORTED: latched in the sample loop above
///   and passed to [`resolve_closed_surface_period_jump`]. Their
///   `IsUClosed(myPreci)` / `IsVClosed(myPreci)` gates are the ported
///   [`sa_is_u_closed`] / [`sa_is_v_closed`]
///   (`ShapeAnalysis_Surface.cxx:660-863`, `:865-1077`).
/// - `theC2D = interpolatePCurve(...)` (`cxx:771`) is ported; `approximatePCurve`
///   (`cxx:2220`) has no caller in this version.
pub fn project_curve_on_surface_perform(
    curve: &dyn Curve,
    surf: &dyn Surface,
    first: f64,
    last: f64,
    preci: f64,
    tol_first: f64,
    tol_last: f64,
    cache: &mut super::projection_cache::ProjectorCache,
) -> Option<Arc<dyn Curve2d>> {
    if !first.is_finite() || !last.is_finite() || last - first < 1e-15 {
        return None;
    }
    if let Some(c2d) = project_analytic(curve, surf) {
        return Some(c2d);
    }
    let n = generate_curve_point_count(curve, first, last);
    let mut t_vals: Vec<f64> = (0..n)
        .map(|i| {
            if i == 0 {
                first
            } else if i + 1 == n {
                last
            } else {
                first + (last - first) * (i as f64) / ((n - 1) as f64)
            }
        })
        .collect();
    let mut pts3d: Vec<GpPnt> = t_vals.iter().map(|&t| curve.d0(t)).collect();

    // `approxPCurve` (`cxx:1106-1143`) tries the straight-pcurve shortcut
    // `getLine` before everything else (`cxx:1119`); a hit gives the exact
    // `Geom2d_Line` / 2-pole pcurve OCCT stores for an isoparametric edge and
    // refills `myCache` from its two endpoints (`cxx:1122-1142`).
    let mut gl = super::projection_cache::GetLineOut::default();
    if let Some(c2d) = super::projection_cache::get_line(surf, &pts3d, &t_vals, preci, cache, &mut gl) {
        let change_cycle = cache.change_cycle(&pts3d[0], &pts3d[pts3d.len() - 1]);
        cache.store(&pts3d, &[gl.first2d, gl.last2d], change_cycle);
        return Some(c2d);
    }

    // `approxPCurve` (`cxx:1147-1287`): the isoparametric test and its pre-block.
    // `cxx:1153-1157`: `isAnalytic` is false only for the exact dynamic types
    // `Geom_BezierSurface` / `Geom_BSplineSurface`. `DynamicType()`, not
    // `IsKind`, is what OCCT tests, so a `Geom_RectangularTrimmedSurface`
    // wrapping either patch stays analytic.
    let is_analytic = !(surf.is_bspline_surface() || surf.osculating_bspline().is_some());
    let (uf, ul) = surf.u_range();
    let (vf, vl) = surf.v_range();
    let mut iso_out = IsoParamOut::new(n);
    let iso_param = is_an_isoparametric(surf, &pts3d, &t_vals, &mut iso_out);
    let mut iso_value = 0.0f64;
    let mut iso_par1 = 0.0f64;
    let mut iso_par2 = 0.0f64;
    let mut value_p1 = iso_out.value_p1;
    let mut value_p2 = iso_out.value_p2;
    let mut iso_closed = false;
    if iso_param {
        // `cxx:1190-1209`.
        let (parf, parl) = if iso_out.is_type_u {
            iso_value = value_p1.x();
            iso_par1 = value_p1.y();
            iso_par2 = value_p2.y();
            iso_closed = sa_is_v_closed(surf, preci);
            (vf, vl)
        } else {
            iso_value = value_p1.y();
            iso_par1 = value_p1.x();
            iso_par2 = value_p2.x();
            iso_closed = sa_is_u_closed(surf, preci);
            (uf, ul)
        };
        if !iso_out.iso_par2d3d && !is_analytic {
            // `cxx:1210-1244`: only for a non-analytic surface whose samples
            // are not already the iso points. OCCT computes `Cf` / `Cl` from
            // `cIso` here but never uses them (`cxx:1211-1222` is dead), so the
            // port omits that read.
            let tdeb = iso_out.params_out[1];
            if iso_closed && (iso_par1 == parf || iso_par1 == parl) {
                // `cxx:1224-1244`.
                if (tdeb - parf).abs() < (tdeb - parl).abs() {
                    iso_par1 = parf;
                } else {
                    iso_par1 = parl;
                }
                if iso_out.is_type_u {
                    value_p1.set_y(iso_par1);
                } else {
                    value_p1.set_x(iso_par1);
                }
            }
            if iso_closed && (iso_par2 == parf || iso_par2 == parl) {
                // `cxx:1245-1265`.
                let tfin = iso_out.params_out[n - 2];
                if (tfin - parf).abs() < (tfin - parl).abs() {
                    iso_par2 = parf;
                } else {
                    iso_par2 = parl;
                }
                if iso_out.is_type_u {
                    value_p2.set_y(iso_par2);
                } else {
                    value_p2.set_x(iso_par2);
                }
            }
            if !iso_closed {
                // `cxx:1267-1285`: swap the ends when the outer end is farther
                // from its sample.
                let tfin = iso_out.params_out[n - 2];
                if (tdeb - iso_par1).abs() > (tdeb - iso_par2).abs()
                    && (tfin - iso_par2).abs() > (tfin - iso_par1).abs()
                {
                    std::mem::swap(&mut value_p1, &mut value_p2);
                    if iso_out.is_type_u {
                        iso_value = value_p1.x();
                        iso_par1 = value_p1.y();
                        iso_par2 = value_p2.y();
                    } else {
                        iso_value = value_p1.y();
                        iso_par1 = value_p1.x();
                        iso_par2 = value_p2.x();
                    }
                }
            }
        }
    }

    // `approxPCurve` (`ShapeConstruct_ProjectCurveOnSurface.cxx:1290-1465`): the
    // point-sampling loop. Every interior sample is projected by
    // `NextValueOfUV` walking from the previous parameter point (`cxx:1432`)
    // instead of by an independent global projection. The two endpoints keep the
    // projections `getLine` left in `thePoints2d` (`cxx:1388-1399`) unless
    // `getLine` reported a fixed period jump (`cxx:1401-1424`), in which case
    // they are re-projected from `myCache` or globally. The `isoParam` arms
    // (`cxx:1317-1375`) and the `p1OnIso` / `p2OnIso` endpoint overrides
    // (`cxx:1379-1385`) come from `isAnIsoparametric`, ported above.
    //
    // `needResolveUJump` / `needResolveVJump` (`cxx:1436-1451`) ARE ported here
    // and consumed by `resolve_closed_surface_period_jump` (`cxx:1504`,
    // `cxx:1596-1598`).
    //
    // `cxx:1291`: the running seed tolerance, `myPreci` before the first sample.
    let mut gap = preci;
    // `cxx:1290-1291`: `Up` / `Vp` are the unrestricted surface bounds'
    // widths (`cxx:1160-1161`).
    let up = ul - uf;
    let vp = vl - vf;
    // `cxx:1296-1304`: the sample order reverses when the incoming cache is
    // anchored on the LAST point of this curve.
    let a_change_cycle = cache.change_cycle(&pts3d[0], &pts3d[n - 1]);
    let mut is_from_cache = false;
    let mut a_saved_point = GpPnt2d::zero();
    let mut pts: Vec<GpPnt2d> = vec![GpPnt2d::zero(); n];
    // `cxx:1020-1021`: `getLine` wrote `thePoints2d(1)` / `thePoints2d(aNb)`
    // before it failed; the `!isRecompute` endpoint arm reads them back.
    pts[0] = gl.first2d;
    pts[n - 1] = gl.last2d;
    // `cxx:1466-1467`: the singularities of the surface. OCCT triggers
    // `ComputeSingularities` lazily from `NbSingularities(myPreci)` inside the
    // loop below (`ShapeAnalysis_Surface.cxx:305-309`) and again at `cxx:1469`;
    // it is a pure function of the surface and its bounds, so it is computed
    // once here and shared by the latch, `correct_degenerated_points` and
    // `resolve_closed_surface_period_jump`.
    let singularities = super::singularities::compute_singularities(surf);
    // `cxx:1308-1310` `p2d`: loop-carried. OCCT leaves it uninitialized, but
    // every first-iteration arm assigns it before use.
    let mut p2d = GpPnt2d::zero();
    // `cxx:1308-1311`: the lattice of the jump test below.
    let mut need_resolve_u_jump = false;
    let mut need_resolve_v_jump = false;
    let mut prev_p3d = GpPnt::zero();
    let mut prev_p2d = GpPnt2d::zero();
    for ii in 1..=n {
        // `cxx:1313`: the sample order may be reversed.
        let a_pnt_index = if a_change_cycle { n - ii + 1 } else { ii };
        let p3d = pts3d[a_pnt_index - 1];
        if iso_param {
            // `cxx:1317-1375`: the isoparametric arms.
            let mut t_par = 0.0f64;
            if iso_out.iso_par2d3d {
                // `cxx:1319-1328`.
                t_par = if iso_par2 > iso_par1 {
                    t_vals[a_pnt_index - 1]
                } else {
                    iso_out.t1 + iso_out.t2 - t_vals[a_pnt_index - 1]
                };
            } else if !is_analytic {
                // `cxx:1330-1346`: the projected iso parameters `pout`.
                t_par = if a_pnt_index == 1 {
                    iso_par1
                } else if a_pnt_index == n {
                    iso_par2
                } else {
                    iso_out.params_out[a_pnt_index - 1]
                };
            }
            if !iso_out.iso_par2d3d && is_analytic {
                // `cxx:1347-1362`: an analytic surface whose samples are not
                // already the iso points still projects by `NextValueOfUV`;
                // only the two ends take the iso values.
                if a_pnt_index == 1 {
                    p2d = value_p1;
                } else if a_pnt_index == n {
                    p2d = value_p2;
                } else {
                    let (q, g) = next_value_of_uv(
                        surf,
                        &p2d,
                        &p3d,
                        preci,
                        occt_core::precision::CONFUSION + 1000.0 * gap,
                    );
                    p2d = q;
                    gap = g;
                }
            } else if iso_out.is_type_u {
                // `cxx:1364-1376`.
                p2d = GpPnt2d::new(iso_value, t_par);
            } else {
                p2d = GpPnt2d::new(t_par, iso_value);
            }
        } else if a_pnt_index == 1 && iso_out.p1_on_iso {
            // `cxx:1379-1381`: endpoint override from `isAnIsoparametric`,
            // which can be set even when `isoParam` is false.
            p2d = value_p1;
        } else if a_pnt_index == n && iso_out.p2_on_iso {
            // `cxx:1382-1384`.
            p2d = value_p2;
        } else if a_pnt_index == 1 || a_pnt_index == n {
            if !gl.is_recompute {
                // `cxx:1388-1399`: no period jump among the `getLine` probes, so
                // the endpoint keeps the seeded-Newton projection there.
                p2d = pts[a_pnt_index - 1];
                gap = gl.last_gap;
                if a_pnt_index == 1 {
                    is_from_cache = gl.is_from_cache;
                    a_saved_point = p2d;
                }
                continue;
            }
            // `cxx:1401-1424`: a period jump was fixed, so re-project the
            // endpoint from `myCache` (walking Newton) or, without a cache hit,
            // globally.
            match cache.find(&p3d, preci * preci) {
                Some(a_cache_pnt) => {
                    let (q, g) = next_value_of_uv(
                        surf,
                        &a_cache_pnt.second,
                        &p3d,
                        preci,
                        occt_core::precision::CONFUSION + gap,
                    );
                    p2d = q;
                    gap = g;
                    if a_pnt_index == 1 {
                        is_from_cache = true;
                        a_saved_point = a_cache_pnt.second;
                    }
                }
                None => {
                    let (q, g) = value_of_uv_with_gap(surf, &p3d, preci);
                    p2d = q;
                    gap = g;
                }
            }
        } else {
            // `cxx:1432`: the interior sample walks from the previous one.
            let (q, g) = next_value_of_uv(
                surf,
                &p2d,
                &p3d,
                preci,
                occt_core::precision::CONFUSION + 1000.0 * gap,
            );
            p2d = q;
            gap = g;
        }
        pts[a_pnt_index - 1] = p2d;
        // `cxx:1436-1451`: `needResolveUJump` / `needResolveVJump`. A sample
        // whose U (V) parameter moved a full period while its 3D point stayed
        // within `myPreci` of the previous one is a seam crossing the point-wise
        // projection resolved the wrong way; the flags widen the consumers below
        // (`cxx:1504`, `cxx:1596-1598`) from "the surface is closed" to "this
        // pcurve jumped a period". Only an interior sample of a longer sequence
        // (`theNbPnt > 23`, `ii > 2`, `ii < theNbPnt`) that lands on a
        // `Geom_BSplineSurface` with a singularity and is not already
        // U(V)-closed qualifies.
        if n > 23 && ii > 2 && ii < n {
            let is_bspline = crate::geom_bnd_lib_surface3d::surface_kind(surf)
                == crate::geom_bnd_lib_surface3d::SurfaceKind::BSplineSurface;
            if is_bspline
                && prev_p3d.distance(&p3d) < preci
                && super::singularities::nb_singularities(&singularities, preci) > 0
            {
                if (p2d.x() - prev_p2d.x()).abs() > 0.95 * up
                    && !sa_is_u_closed(surf, preci)
                {
                    need_resolve_u_jump = true;
                }
                if (p2d.y() - prev_p2d.y()).abs() > 0.95 * vp
                    && !sa_is_v_closed(surf, preci)
                {
                    need_resolve_v_jump = true;
                }
            }
        }
        // `cxx:1452-1453`.
        prev_p3d = p3d;
        prev_p2d = p2d;
        // `cxx:1455-1463`: the seed for the next sample is the linear
        // extrapolation of this sample against the previously stored one.
        if ii > 1 {
            let base = if a_change_cycle {
                pts[a_pnt_index]
            } else {
                pts[a_pnt_index - 2]
            };
            p2d = GpPnt2d::new(2.0 * p2d.x() - base.x(), 2.0 * p2d.y() - base.y());
        }
    }
    // `approxPCurve` (`cxx:1466-1470`): the singularities of the surface
    // (`ShapeAnalysis_Surface::ComputeSingularities`) feed both the
    // degenerate-point projection of the sampled pcurve and the singularity
    // loop that corrects its ends. The degenerate-point projection is guarded
    // by `if (!isoPar2d3d)` at
    // `ShapeConstruct_ProjectCurveOnSurface.cxx:1467`.
    if !iso_out.iso_par2d3d {
        super::singularities::correct_degenerated_points(
            surf,
            curve,
            &singularities,
            preci,
            tol_first,
            tol_last,
            &pts3d,
            &t_vals,
            &mut pts,
        );
    }
    // `approxPCurve` (`cxx:1499-1890`): before any curve is built, remove the
    // full-period jumps that the point-wise projection left in the sampled UV
    // sequence of a closed surface, then run "Handle AdjustOverDegen".
    resolve_closed_surface_period_jump(
        surf,
        curve,
        &singularities,
        preci,
        need_resolve_u_jump,
        need_resolve_v_jump,
        is_from_cache,
        a_saved_point,
        &mut t_vals,
        &mut pts3d,
        &mut pts,
    );
    // `Perform:cxx:771`: `theC2D = interpolatePCurve(aNbPini, aPoints2d, aParams)`,
    // the single other source of `theC2D` besides `getLine`. A `None` result is
    // OCCT's null `theC2D` (`cxx:773`, `ShapeExtend_FAIL1`).
    let c2d = super::projection_cache::interpolate_pcurve(&pts, &t_vals, preci)?;
    // `cxx:1892-1903`: `myCache` is refilled from the pcurve just built, so the
    // next pcurve on this surface is projected from these endpoints.
    let change_cycle = cache.change_cycle(&pts3d[0], &pts3d[pts3d.len() - 1]);
    cache.store(&pts3d, &pts, change_cycle);
    Some(c2d)
}

pub fn project_curve_on_surface(
    s: &dyn Surface,
    curve: &dyn Curve,
    t_vals: &[f64],
) -> Option<Vec<GpPnt2d>> {
    let kind = classify_surface_kind(s);
    if kind != SurfaceKind::Torus {
        return None;
    }
    let proj = Projector::from_surface(s, kind)?;
    Some(proj.project_seq(curve, t_vals))
}

/// Recover the plane's location and in-plane axes from its surface `d0`.
pub(super) fn plane_axes(surf: &dyn Surface) -> Option<(GpPnt, GpVec, GpVec)> {
    let o = surf.d0(0.0, 0.0);
    let x = GpVec::from_pnts(&o, &surf.d0(1.0, 0.0));
    let y = GpVec::from_pnts(&o, &surf.d0(0.0, 1.0));
    if x.square_magnitude() < 1e-24 || y.square_magnitude() < 1e-24 {
        return None;
    }
    Some((o, x, y))
}

/// Recover a cylinder's axis point, unit axes and radius from its surface `d0`.
pub(super) fn cylinder_axes(surf: &dyn Surface) -> Option<(GpPnt, GpVec, GpVec, GpVec)> {
    let a = surf.d0(0.0, 0.0);
    let b = surf.d0(PI, 0.0);
    let c = surf.d0(FRAC_PI_2, 0.0);
    let o = mid(&a, &b);
    let r = 0.5 * a.distance(&b);
    if r < 1e-12 {
        return None;
    }
    let x = GpVec::new((a.x() - o.x()) / r, (a.y() - o.y()) / r, (a.z() - o.z()) / r);
    let y = GpVec::new((c.x() - o.x()) / r, (c.y() - o.y()) / r, (c.z() - o.z()) / r);
    let z = GpVec::from_pnts(&a, &surf.d0(0.0, 1.0));
    if z.square_magnitude() < 1e-24 {
        return None;
    }
    Some((o, x, y, z))
}

/// Orthonormal radial frame of a cone: `x` is the `u = 0` radial direction
/// (the generatrix direction projected perpendicular to the axis).
pub(super) fn cone_frame(s: &dyn Surface, _apex: &GpPnt, z: &GpVec) -> Option<(GpVec, GpVec, GpVec)> {
    let dz = GpVec::from_pnts(&s.d0(0.0, 0.0), &s.d0(0.0, 1.0));
    let x = dz.subtracted(&z.multiplied_scalar(dz.dot(z))).normalized();
    if x.square_magnitude() < 1e-20 {
        return None;
    }
    let y = z.crossed(&x).normalized();
    Some((x, y, *z))
}

/// The cone's placement point (the v=0 ring plane, `ElSLib::ConeValue` base).
///
/// `s.d0(0,0)` lies on the v=0 ring: `placement + RefRadius·x`. Subtracting
/// its radial component about the axis (from the apex) recovers the placement
/// point on the axis.
pub(super) fn cone_placement(s: &dyn Surface, apex: &GpPnt, z: &GpVec) -> Option<GpPnt> {
    let p0 = s.d0(0.0, 0.0);
    let d = GpVec::from_pnts(apex, &p0);
    let radial = d.subtracted(&z.multiplied_scalar(d.dot(z)));
    Some(GpPnt::from_xyz(&p0.coord.subtracted(&radial.coord)))
}

// ---------------------------------------------------------------------------
// Analytic isoparametric pcurves
// ---------------------------------------------------------------------------

/// `Adaptor3d_Curve::GetType() == GeomAbs_Line`: a `Geom_Line` (or a trimmed
/// one). This is the *concrete class* test of
/// `GeomAdaptor_Curve::load` (`GeomAdaptor_Curve.cxx:252-311`), which unwraps a
/// `Geom_TrimmedCurve` to its basis (the port's trimmed wrappers forward the
/// query) and then tests `gp_line()`. It replaces the earlier
/// "unbounded parameter range" heuristic, which also matched an offset curve.
pub(super) fn is_line(curve: &dyn Curve) -> bool {
    curve.gp_line().is_some()
}

/// `Adaptor3d_Curve::GetType() == GeomAbs_Circle`: a `Geom_Circle` (or a
/// trimmed one). Deliberately NOT a geometric test - a `Geom_BSplineCurve`
/// whose image is a circle stays `GeomAbs_BSplineCurve` in OCCT
/// (`GeomAdaptor_Curve::load`, `GeomAdaptor_Curve.cxx:252-311`; the same order
/// as `edge_edge/find_solutions.rs:78-98`), and
/// `ProjLib_ProjectedCurve::Project` (`ProjLib_ProjectedCurve.cxx:242-270`)
/// only has analytic overloads for `GeomAbs_Line` (`:247`), `GeomAbs_Circle`
/// (`:250`), `GeomAbs_Ellipse` (`:253`), `GeomAbs_Hyperbola` (`:256`) and
/// `GeomAbs_Parabola` (`:259`): `GeomAbs_BSplineCurve` / `BezierCurve` /
/// `OffsetCurve` / `OtherCurve` break out (`:262-266`) into the general
/// approximation (`ProjLib_ComputeApprox`, `ProjLib_ProjectedCurve.cxx:720`).
pub(super) fn is_circle(curve: &dyn Curve) -> bool {
    curve.gp_circ().is_some()
}

/// Fit a clamped B-spline through `pts` with the exact parameter range [a, b].
/// Degree 1 is a polyline through the sample points (endpoints interpolated).
pub(super) fn bspline_from_samples(pts: &[GpPnt2d], a: f64, b: f64, degree: usize) -> Result<Geom2dBSplineCurve, String> {
    if pts.len() < 2 {
        return Err("pcurve_full: need at least 2 points".into());
    }
    if degree >= pts.len() {
        return Err("pcurve_full: degree must be below pole count".into());
    }
    let mut knots = occt_core::bspl::knots::build_uniform_knots(pts.len(), degree);
    for k in knots.iter_mut() {
        *k = a + (b - a) * *k;
    }
    let xs: Vec<f64> = pts.iter().map(|p| p.x()).collect();
    let ys: Vec<f64> = pts.iter().map(|p| p.y()).collect();
    Geom2dBSplineCurve::new(xs, ys, knots, degree).map_err(|e| e.to_string())
}

/// Build a `Geom2dLine` (unit-speed) or a degree-1 B-spline (otherwise) through
/// the first and last projected points over the parameter range `[a, b]`.
pub(super) fn line_from_pts(pts: &[GpPnt2d], a: f64, b: f64) -> Option<Arc<dyn Curve2d>> {
    let q0 = pts[0];
    let q1 = pts[pts.len() - 1];
    let d = GpVec2d::new(q1.x() - q0.x(), q1.y() - q0.y());
    let mag = d.magnitude();
    if mag < 1e-15 {
        return None;
    }
    let dir = GpDir2d::from_vec2d(&d).ok()?;
    let len = b - a;
    if (mag - len).abs() <= 1e-6 * len.max(1.0) {
        let loc = GpPnt2d::new(q0.x() - a * dir.x, q0.y() - a * dir.y);
        Some(Arc::new(Geom2dLine::new(GpAx2d::new(loc, dir))))
    } else {
        let bs = bspline_from_samples(pts, a, b, 1).ok()?;
        Some(Arc::new(bs))
    }
}

/// Analytic pcurve of a line edge on an analytic revolution surface: a
/// generatrix projects to a `u`-isoline (varying `v`) or `v`-isoline.
/// `None` when the projection is not an isoparametric line (the caller then
/// samples instead).
pub(super) fn iso_line_pcurve(curve: &dyn Curve, proj: &Projector, a: f64, b: f64) -> Option<Arc<dyn Curve2d>> {
    if b - a < 1e-15 {
        return None;
    }
    let t_vals = [a, 0.5 * (a + b), b];
    let pts = proj.project_seq(curve, &t_vals);
    let tol = 1e-6 * (1.0 + pts[0].x().abs() + pts[0].y().abs());
    let u_const = pts.iter().all(|q| (q.x() - pts[0].x()).abs() < tol);
    let v_const = pts.iter().all(|q| (q.y() - pts[0].y()).abs() < tol);
    if !u_const && !v_const {
        return None;
    }
    line_from_pts(&pts, a, b)
}

/// Analytic pcurve of a circle edge on an analytic revolution surface: a
/// latitude (constant `v`) or meridian (constant `u`) projects to a straight
/// UV line; a U-periodic full circle unwraps across the seam. `None` when the
/// projection is not an isoparametric line.
pub(super) fn iso_circle_pcurve(curve: &dyn Curve, proj: &Projector, a: f64, b: f64) -> Option<Arc<dyn Curve2d>> {
    if b - a < 1e-15 {
        return None;
    }
    let n = 9;
    let t_vals: Vec<f64> = (0..n).map(|i| a + (b - a) * i as f64 / (n - 1) as f64).collect();
    let pts = proj.project_seq(curve, &t_vals);
    let tol = 1e-6 * (1.0 + pts[0].x().abs() + pts[0].y().abs());
    let u_const = pts.iter().all(|q| (q.x() - pts[0].x()).abs() < tol);
    let v_const = pts.iter().all(|q| (q.y() - pts[0].y()).abs() < tol);
    if u_const == v_const {
        // Degenerate (both) or non-isoparametric (neither): let the sampling
        // path handle it.
        return None;
    }
    line_from_pts(&pts, a, b)
}

// ---------------------------------------------------------------------------
// Main entry points
// ---------------------------------------------------------------------------

/// A 2D curve linearly reparameterized from `[p0, p1]` onto `[a, b]`: a parameter
/// `t ∈ [a, b]` maps to `inner.d0(p0 + (t − a)·(p1 − p0)/(b − a))`. Aligns a
/// STEP-imported pcurve (parameterized over its natural range) with the edge's
/// 3D range.
#[derive(Clone)]
pub struct ReparamCurve2d {
    pub(super) inner: Arc<dyn Curve2d>,
    pub(super) a: f64,
    pub(super) b: f64,
    pub(super) p0: f64,
    pub(super) p1: f64,
}

impl ReparamCurve2d {
    pub fn new(inner: Arc<dyn Curve2d>, a: f64, b: f64, p0: f64, p1: f64) -> Self {
        Self { inner, a, b, p0, p1 }
    }

    pub(super) fn scale(&self) -> f64 {
        if (self.b - self.a).abs() < 1e-30 {
            1.0
        } else {
            (self.p1 - self.p0) / (self.b - self.a)
        }
    }

    pub(super) fn map(&self, t: f64) -> f64 {
        self.p0 + (t - self.a) * self.scale()
    }
}

impl Curve2d for ReparamCurve2d {
    fn d0(&self, u: f64) -> GpPnt2d {
        self.inner.d0(self.map(u))
    }
    fn d1(&self, u: f64) -> (GpPnt2d, GpVec2d) {
        let (p, d) = self.inner.d1(self.map(u));
        let s = self.scale();
        (p, GpVec2d::new(d.x() * s, d.y() * s))
    }
    fn d2(&self, u: f64) -> (GpPnt2d, GpVec2d, GpVec2d) {
        let (p, d1, d2) = self.inner.d2(self.map(u));
        let s = self.scale();
        (
            p,
            GpVec2d::new(d1.x() * s, d1.y() * s),
            GpVec2d::new(d2.x() * s * s, d2.y() * s * s),
        )
    }
    fn first_parameter(&self) -> f64 {
        self.a
    }
    fn last_parameter(&self) -> f64 {
        self.b
    }
    fn continuity(&self) -> u8 {
        self.inner.continuity()
    }
    fn transform(&mut self, t: &GpTrsf2d) {
        let mut inner = self.inner.clone_dyn();
        inner.transform(t);
        self.inner = Arc::from(inner);
    }
    fn reverse(&mut self) {
        let mut inner = self.inner.clone_dyn();
        inner.reverse();
        self.inner = Arc::from(inner);
    }
    fn clone_dyn(&self) -> Box<dyn Curve2d> {
        Box::new(self.clone())
    }
    // Preserve Geom2d_Line kind through the wrapper so
    // `Adaptor3d_CurveOnSurface::EvalKPart` (`cxx:1552-1732`) and
    // `GeomLib::SameRange` Line translate (`cxx:862-870`) keep iso/line
    // detection. A bare default `is_line=false` densifies Offset seams.
    fn is_line(&self) -> bool {
        self.inner.is_line()
    }
    fn gp_lin2d(&self) -> Option<occt_core::gp::GpLin2d> {
        self.inner.gp_lin2d()
    }
    // Kind delegation. `Geom2dAdaptor_Curve::load` unwraps the stored curve
    // before classifying it (`Geom2dAdaptor_Curve.cxx:226-230`), so the
    // wrapper's sample count (`Geom2dAdaptor_Curve::NbSamples`,
    // `Geom2dAdaptor_Curve.cxx:1351-1394`) and bounding-box sample budget
    // (`Geom2dAdaptor_Curve::GetType`, mirrored by
    // `geom_bnd_lib_sample2d::sample_kind_of`) must come from the inner curve.
    // A bare trait default makes a reparameterized B-spline pcurve report
    // `GeomAbs_OtherCurve` and fall back to 20 samples.
    fn gp_circ2d(&self) -> Option<occt_core::gp::GpCirc2d> {
        self.inner.gp_circ2d()
    }
    fn bezier_nb_poles(&self) -> Option<usize> {
        self.inner.bezier_nb_poles()
    }
    fn bspline_nb_knots(&self) -> Option<usize> {
        self.inner.bspline_nb_knots()
    }
    fn bspline_degree(&self) -> Option<usize> {
        self.inner.bspline_degree()
    }
    fn offset_basis(&self) -> Option<&dyn Curve2d> {
        self.inner.offset_basis()
    }
}

/// Wrap a 2D curve in a linear reparameterization from `[p0, p1]` onto `[a, b]`.
pub fn reparam_curve2d(
    inner: Arc<dyn Curve2d>,
    a: f64,
    b: f64,
    p0: f64,
    p1: f64,
) -> Arc<dyn Curve2d> {
    Arc::new(ReparamCurve2d::new(inner, a, b, p0, p1))
}

/// Parameter of the 2D curve closest to `uv` (the 1-D projection of a surface
/// point onto the pcurve). A bounded curve is searched on `[first, last]`; an
/// unbounded line is projected in closed form.
pub fn project_uv_on_curve2d(curve: &dyn Curve2d, uv: GpPnt2d) -> f64 {
    let (f, l) = (curve.first_parameter(), curve.last_parameter());
    if f.is_finite() && l.is_finite() {
        let mut best = f;
        let mut best_d = f64::INFINITY;
        for i in 0..=64 {
            let u = f + (l - f) * i as f64 / 64.0;
            let d = curve.d0(u).distance(&uv);
            if d < best_d {
                best_d = d;
                best = u;
            }
        }
        let mut step = (l - f) / 64.0;
        for _ in 0..32 {
            step *= 0.5;
            let mut u = best;
            let mut dd = best_d;
            for &du in &[-step, step] {
                let nu = (u + du).clamp(f, l);
                let d = curve.d0(nu).distance(&uv);
                if d < dd {
                    dd = d;
                    u = nu;
                }
            }
            best = u;
            best_d = dd;
        }
        best
    } else {
        let p0 = curve.d0(0.0);
        let dir = curve.d0(1.0).coord.subtracted(&p0.coord);
        let len2 = dir.square_modulus();
        if len2 < 1e-30 {
            return 0.0;
        }
        uv.coord.subtracted(&p0.coord).dot(&dir) / len2
    }
}
