use super::prelude::*;
use super::*;

/// Parameter of the pcurve whose 3D image on `surface` is closest to `p`
/// (`Extrema_LocateExtPC` on a pcurve). The distance `F(u)=|C(u)−p|²` has an
/// extremum when `g(u)=dF/du=(p−C(u))·C'(u)=0`; a Newton step from `u_guess`
/// (the linear-scaled parameter) finds that root locally, as
/// `math_FunctionRoot` does. No global scan — the pcurve may have several
/// stationary points and OCCT keeps the one nearest the guess.
pub fn project_point_on_pcurve(
    curve: &dyn Curve2d,
    surface: &dyn Surface,
    p: &GpPnt,
    u_guess: f64,
) -> f64 {
    let dist = |u: f64| {
        let uv = curve.d0(u);
        surface.d0(uv.x(), uv.y()).distance(p)
    };
    let (f, l) = (curve.first_parameter(), curve.last_parameter());
    let (lo, hi) = if f.is_finite() && l.is_finite() {
        (f, l)
    } else {
        (u_guess - 1e6, u_guess + 1e6)
    };
    let mut u = u_guess.clamp(lo, hi);
    let mut best = dist(u);
    let mut step = 1.0;
    for _ in 0..64 {
        let mut improved = false;
        for &du in &[-step, step] {
            let nu = (u + du).clamp(lo, hi);
            let d = dist(nu);
            if d < best {
                best = d;
                u = nu;
                improved = true;
            }
        }
        if !improved {
            step *= 0.5;
            if step < 1e-12 {
                break;
            }
        }
    }
    u
}

/// Construct the pcurve of `edge` on `face` — the edge's 3D curve projected
/// into the face surface's `(u, v)` parameter domain, with analytic cases for
/// every analytic surface and a sampling fallback for B-spline faces.
///
/// The returned pcurve is parameterized over the edge's range `[a, b]` so
/// `d0(a)` / `d0(b)` land on the projected endpoints.
/// Select which of two seam pcurves is the forward one (its 2D direction matches
/// the edge's 3D direction), port of `ShapeAnalysis_Curve::SelectForwardSeam`.
/// Returns 1 (first pcurve) or 2 (second pcurve).
pub fn select_forward_seam(c1: &dyn Curve2d, c2: &dyn Curve2d) -> usize {
    // `ShapeAnalysis_Curve.cxx:852-969`. A `Geom2d_Line` uses Location/Direction;
    // a bounded curve uses the Start→End chord (same signs for an isoline).
    let line_of = |c: &dyn Curve2d| -> Option<(GpPnt2d, GpVec2d)> {
        let (a, b) = (c.first_parameter(), c.last_parameter());
        let (a, b) = if a.is_finite() && b.is_finite() { (a, b) } else { (0.0, 1.0) };
        let p0 = c.d0(a);
        let p1 = c.d0(b);
        let v = GpVec2d::new(p1.x() - p0.x(), p1.y() - p0.y());
        if v.square_magnitude() < occt_core::precision::RESOLUTION {
            return None;
        }
        Some((p0, v))
    };
    let Some((loc1, d1)) = line_of(c1) else {
        return 0;
    };
    let Some((loc2, _)) = line_of(c2) else {
        return 0;
    };
    let ux = d1.x();
    let uy = d1.y();
    if ux > 0.0 {
        if loc1.y() < loc2.y() { 1 } else { 2 }
    } else if ux < 0.0 {
        if loc1.y() < loc2.y() { 2 } else { 1 }
    } else if uy > 0.0 {
        if loc1.x() > loc2.x() { 1 } else { 2 }
    } else if uy < 0.0 {
        if loc1.x() > loc2.x() { 2 } else { 1 }
    } else {
        0
    }
}

pub fn make_pcurve_full(edge: &Edge, face: &Face) -> Result<Arc<dyn Curve2d>, String> {
    // A STEP-imported pcurve (BRep_TEdge's stored `(face -> Geom2d_Curve)`)
    // wins over projection: it is the exact trimming curve of the face, so a
    // B-spline face whose 3D trimming curve does not lie on the surface still
    // gets the correct 2D boundary. A seam edge stores two (forward/reversed
    // sides); pick the one matching the edge's traversal orientation.
    let face_key = GeometryRegistry::shape_key(&face.0);
    let pcs = GeometryRegistry::global().edge_pcurves(&edge.0, face_key);
    if pcs.len() >= 2 {
        let idx = if edge.0.orientation().is_reversed() { 1 } else { 0 };
        return Ok(pcs[idx].clone());
    }
    if let Some(pc) = pcs.first() {
        return Ok(pc.clone());
    }
    let Some(curve) = GeometryRegistry::global().edge_curve(&edge.0) else {
        return Err("make_pcurve_full: edge has no 3D curve".into());
    };
    let (a, b) = GeometryRegistry::global().edge_parameters(&edge.0);
    if !a.is_finite() || !b.is_finite() || b - a < 1e-15 {
        return Err("make_pcurve_full: edge range is empty or unbounded".into());
    }
    let Some(surf) = GeometryRegistry::global().face_surface(&face.0) else {
        return Err("make_pcurve_full: face has no surface".into());
    };
    let kind = classify_surface_kind(surf.as_ref());

    // Plane faces: reuse the Phase-16a implementation verbatim, so
    // make_pcurve_full and pcurve::make_pcurve_on_face agree exactly.
    if kind == SurfaceKind::Plane {
        return crate::pcurve::make_pcurve_on_face(edge, face);
    }

    // Analytic revolution surfaces: exact isoparametric pcurves.
    if let Some(proj) = Projector::from_surface(surf.as_ref(), kind) {
        if is_line(curve.as_ref()) {
            if let Some(pc) = iso_line_pcurve(curve.as_ref(), &proj, a, b) {
                return Ok(pc);
            }
        }
        if is_circle(curve.as_ref()) {
            if let Some(pc) = iso_circle_pcurve(curve.as_ref(), &proj, a, b) {
                return Ok(pc);
            }
        }
    }

    // General B-spline surface (or a non-isoparametric edge): sample the
    // projection and fit a clamped degree-1 B-spline over the edge range.
    let pts = edge_pcurve_on_face(edge, face, 48);
    if pts.len() < 2 {
        return Err("make_pcurve_full: too few projected samples".into());
    }
    let bs = bspline_from_samples(&pts, a, b, 1)?;
    Ok(Arc::new(bs))
}

/// Direction correction of a pcurve: `+1` when the pcurve parameter follows the
/// edge's 3D direction, `-1` when it must be reversed.
///
/// The surface-mapped pcurve tangent `∂S/∂u·u' + ∂S/∂v·v'` is compared against
/// the edge's 3D tangent (flipped for a REVERSED edge). Surfaces whose `d1`
/// returns zero vectors fall back to comparing endpoint displacements.
pub fn pc_curve_orientation(edge: &Edge, face: &Face, curve: &dyn Curve2d) -> f64 {
    let (a, b) = GeometryRegistry::global().edge_parameters(&edge.0);
    if !a.is_finite() || !b.is_finite() || b - a < 1e-15 {
        return 1.0;
    }
    let tm = 0.5 * (a + b);
    let (Some(c3d), Some(surf)) = (
        GeometryRegistry::global().edge_curve(&edge.0),
        GeometryRegistry::global().face_surface(&face.0),
    ) else {
        return 1.0;
    };

    let (_, mut d3d) = c3d.d1(tm);
    if edge.orientation().is_reversed() {
        d3d = d3d.reversed();
    }
    let q = curve.d0(tm);
    let (_, du, dv) = surf.d1(q.x(), q.y());
    let d2d = curve.d1(tm).1;
    let mapped = du.multiplied_scalar(d2d.x()).add(&dv.multiplied_scalar(d2d.y()));
    let dot = d3d.dot(&mapped);
    if dot.abs() < 1e-12 {
        // Degenerate surface tangent (analytic `d1` may be zero): compare the
        // endpoint displacement of the edge curve against the UV-mapped one.
        let p0 = c3d.d0(a);
        let p1 = c3d.d0(b);
        let mut disp3 = GpVec::from_pnts(&p0, &p1);
        if edge.orientation().is_reversed() {
            disp3 = disp3.reversed();
        }
        let q0 = curve.d0(a);
        let q1 = curve.d0(b);
        let suv0 = surf.d0(q0.x(), q0.y());
        let suv1 = surf.d0(q1.x(), q1.y());
        let disp_uv = GpVec::from_pnts(&suv0, &suv1);
        return if disp3.dot(&disp_uv) >= 0.0 { 1.0 } else { -1.0 };
    }
    if dot > 0.0 {
        1.0
    } else {
        -1.0
    }
}

/// Wrap `x` into the half-open interval `[lo, hi)` by adding/subtracting the
/// period `hi − lo` (a no-op when the interval is degenerate).
pub(super) fn wrap_into(x: f64, lo: f64, hi: f64) -> f64 {
    let p = hi - lo;
    if p <= 0.0 || !p.is_finite() {
        return x;
    }
    let mut r = (x - lo) % p;
    if r < 0.0 {
        r += p;
    }
    lo + r
}

/// The `u`-period of a periodic surface (its `u`-range span).
pub(super) fn u_period(s: &dyn Surface) -> f64 {
    let (u0, u1) = s.u_range();
    if u0.is_finite() && u1.is_finite() {
        u1 - u0
    } else {
        2.0 * PI
    }
}

/// The `v`-period of a periodic surface (its `v`-range span).
pub(super) fn v_period(s: &dyn Surface) -> f64 {
    let (v0, v1) = s.v_range();
    if v0.is_finite() && v1.is_finite() {
        v1 - v0
    } else {
        2.0 * PI
    }
}

/// Bring a pcurve inside the face's UV bounds.
///
/// For periodic surface dimensions (the `u` of cylinder/cone/sphere/torus and
/// the `v` of a torus) the curve is translated by whole periods so its midpoint
/// lies in range. If part of the curve still leaves a non-periodic bounded
/// **UNPORTED (A15/T-51 余项)**: OCCT has no trimming step on this position.
/// `BOPTools_AlgoTools2D::AdjustPCurveOnSurf` (`BOPTools_AlgoTools2D.cxx:247-400`)
/// only translates the pcurve by whole surface periods and leaves its range
/// alone; the faithful port is
/// [`crate::algo_tools::AlgoTools2D::adjust_pcurve_on_surf`], which
/// `pave_blocks/p01.rs` now calls. This trim helper is kept only as a utility for
/// callers that explicitly want the pcurve clipped to the face rectangle; do not
/// use it in an `AdjustPCurveOnFace` position.
/// dimension, it is trimmed to the in-bounds parameter subrange.
pub fn trim_pcurve_to_face(curve: &Arc<dyn Curve2d>, face: &Face, _tol: f64) -> Result<Arc<dyn Curve2d>, String> {
    let (umin, umax, vmin, vmax) = face_uv_bounds(face);
    let surf = GeometryRegistry::global()
        .face_surface(&face.0)
        .ok_or("trim_pcurve_to_face: face has no surface")?;
    let u_periodic = surf.is_u_periodic();
    let v_periodic = surf.is_v_periodic();
    let up = if u_periodic { u_period(surf.as_ref()) } else { 0.0 };
    let vp = if v_periodic { v_period(surf.as_ref()) } else { 0.0 };

    let (mut c0, mut c1) = (curve.first_parameter(), curve.last_parameter());
    if !c0.is_finite() || !c1.is_finite() || c1 - c0 < 1e-15 {
        // Unbounded 2D curve (a Geom2dLine): evaluate a unit window around 0.
        c0 = 0.0;
        c1 = 1.0;
    }

    // Sample the curve and compute the whole-period shift that brings the
    // midpoint inside the periodic bounds.
    let n = 33;
    let t_vals: Vec<f64> = (0..n).map(|i| c0 + (c1 - c0) * i as f64 / (n - 1) as f64).collect();
    let qm = curve.d0(0.5 * (c0 + c1));
    let mut du = 0.0;
    let mut dv = 0.0;
    if u_periodic && up > 0.0 && umin.is_finite() {
        du = wrap_into(qm.x(), umin, umin + up) - qm.x();
    }
    if v_periodic && vp > 0.0 && vmin.is_finite() {
        dv = wrap_into(qm.y(), vmin, vmin + vp) - qm.y();
    }

    let mut tr = GpTrsf2d::identity();
    tr.set_translation_vec(&GpVec2d::new(du, dv));
    let shifted: Arc<dyn Curve2d> = Arc::from(curve.transformed(&tr));

    let inside = |q: &GpPnt2d| {
        (if umin.is_finite() { q.x() >= umin - 1e-7 } else { true })
            && (if umax.is_finite() { q.x() <= umax + 1e-7 } else { true })
            && (if vmin.is_finite() { q.y() >= vmin - 1e-7 } else { true })
            && (if vmax.is_finite() { q.y() <= vmax + 1e-7 } else { true })
    };

    let shifted_pts: Vec<GpPnt2d> = t_vals.iter().map(|&t| shifted.d0(t)).collect();
    if shifted_pts.iter().all(inside) {
        return Ok(shifted);
    }

    // Some samples still violate a non-periodic bounded dimension: trim to the
    // in-bounds parameter subrange around the midpoint.
    let in_bounds: Vec<bool> = shifted_pts.iter().map(inside).collect();
    let Some(first) = in_bounds.iter().position(|&b| b) else {
        return Err("trim_pcurve_to_face: pcurve lies entirely outside the face UV bounds".into());
    };
    let last = in_bounds.iter().rposition(|&b| b).unwrap();
    let t0 = t_vals[first];
    let t1 = t_vals[last];
    let kept: Vec<GpPnt2d> = shifted_pts[first..=last].to_vec();
    let bs = bspline_from_samples(&kept, t0, t1, 1)?;
    Ok(Arc::new(bs))
}

/// Evaluate a pcurve at `param`, validating that the parameter is inside the
/// curve's range (a `Geom2dLine` accepts any parameter).
pub fn pcurve_point_on_face(curve: &dyn Curve2d, param: f64) -> Result<GpPnt2d, String> {
    let (a, b) = (curve.first_parameter(), curve.last_parameter());
    if a.is_finite() && b.is_finite() && (param < a - 1e-9 || param > b + 1e-9) {
        return Err(format!(
            "pcurve_point_on_face: parameter {param} outside range [{a}, {b}]"
        ));
    }
    Ok(curve.d0(param))
}
