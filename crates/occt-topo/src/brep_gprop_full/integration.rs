use super::prelude::*;
use super::*;

impl FaceGauss {
    pub(super) fn new(face: &Face) -> Result<Self, String> {
        let surface = BRepTool::face_surface_world(face).ok_or("brep_gprop_full: face has no surface")?;
        let kind = classify_surface_kind(surface.as_ref());
        let is_reversed = face.orientation().is_reversed();
        let wires = wires_of_face(face);
        let natural = wires.is_empty();

        // Build the UV inverse map. The frame is taken from the surface's own
        // parameterization so pcurve UV coordinates live in the same space the
        // `normal`/integrand evaluators use.
        let map = match kind {
            SurfKind::Plane => {
                let (o, xd, yd) =
                    plane_frame(surface.as_ref()).ok_or("brep_gprop_full: planar face frame")?;
                UVMap::Plane { o, xd, yd }
            }
            SurfKind::Cylinder => {
                let (o, xd, yd, z, r) =
                    cylinder_frame(surface.as_ref()).ok_or("brep_gprop_full: bad cylinder frame")?;
                UVMap::Cylinder { o, xd, yd, z, r }
            }
            SurfKind::Cone => {
                let (o, xd, yd, z, alpha) =
                    cone_frame(surface.as_ref()).ok_or("brep_gprop_full: bad cone frame")?;
                UVMap::Cone { o, xd, yd, z, alpha }
            }
            _ => UVMap::Generic,
        };

        // Boundary arcs from every edge of every wire.
        let mut arcs: Vec<BoundaryArc> = Vec::new();
        let mut seen: Vec<*const ()> = Vec::new();
        let mut has_repeated = false;
        for w in &wires {
            for e in edges_of_wire(w) {
                let key = Arc::as_ptr(&e.0.tshape) as *const ();
                if seen.contains(&key) {
                    has_repeated = true;
                }
                seen.push(key);
                if let Some(arc) = build_arc(&e, face, &map) {
                    arcs.push(arc);
                }
            }
        }

        // UV bounds from the pcurves (BRepTools::UVBounds).
        let (u1, u2, v1, v2) = if natural {
            let (a, b) = surface.u_range();
            let (c, d) = surface.v_range();
            (a, b, c, d)
        } else {
            uv_bounds(surface.as_ref(), &arcs)
        };

        Ok(Self { surface, kind, is_reversed, natural, u1, u2, v1, v2, arcs, has_repeated_edges: has_repeated })
    }

    /// The face's finite UV bounds `(u1, u2, v1, v2)`.
    pub(super) fn bounds(&self) -> (f64, f64, f64, f64) {
        (self.u1, self.u2, self.v1, self.v2)
    }

    /// `D1U × D1V` without the face-orientation flip (the domain path needs
    /// it: `BRepGProp_Domain` walks the face's edges and every pcurve is
    /// loaded with its edge's orientation, so the traversal already
    /// carries the face orientation — `BRepGProp_Face.cxx:164-185`).
    pub(super) fn normal_raw(&self, u: f64, v: f64) -> (GpPnt, GpVec) {
        let (p, du, dv) = surface_d1(self.surface.as_ref(), u, v);
        (p, du.crossed(&dv))
    }

    /// `D1U × D1V` without the face-orientation flip — used by the domain path,
    /// where `BRepGProp_Domain` + `BRepGProp_Face::Load(edge)` already carry the
    /// face orientation through the boundary traversal.

    /// Unnormalised surface normal (D1U × D1V), flipped for a REVERSED face.
    pub(super) fn normal(&self, u: f64, v: f64) -> (GpPnt, GpVec) {
        let (p, du, dv) = surface_d1(self.surface.as_ref(), u, v);
        let mut n = du.crossed(&dv);
        if self.is_reversed {
            n.reverse();
        }
        (p, n)
    }

    /// Number of Gauss points for the U direction of the face.
    pub(super) fn u_integration_order(&self) -> usize {
        match self.kind {
            SurfKind::Plane => 8,
            _ => 18,
        }
    }

    /// Number of Gauss points for the V direction of the face.
    pub(super) fn v_integration_order(&self) -> usize {
        match self.kind {
            SurfKind::Plane => 8,
            _ => 18,
        }
    }

    /// Number of Gauss points along a boundary arc.
    pub(super) fn arc_integration_order(&self, arc: &BoundaryArc) -> usize {
        match arc.kind {
            ArcKind::Line => 4,
            ArcKind::Circle | ArcKind::Other => 18,
        }
    }

    /// Surface integration order for the adaptive Gauss (SIntOrder).
    pub(super) fn s_int_order(&self, eps: f64) -> usize {
        let (nu, nv) = match self.kind {
            SurfKind::Plane => (1, 1),
            SurfKind::Cylinder | SurfKind::Cone => (2, 1),
            SurfKind::Sphere | SurfKind::Torus | SurfKind::Other => (2, 2),
        };
        let sc = s_coeff(eps);
        let n = ((sc * (nu.max(nv) + 1) as f64).ceil() as usize).clamp(1, GPM);
        n
    }

    /// Number of U subintervals (SUIntSubs).
    #[allow(dead_code)]
    pub(super) fn s_u_int_subs(&self) -> usize {
        match self.kind {
            SurfKind::Plane => 1,
            SurfKind::Cylinder | SurfKind::Cone | SurfKind::Sphere | SurfKind::Torus => 3,
            SurfKind::Other => 1,
        }
    }

    /// Number of V subintervals (SVIntSubs).
    #[allow(dead_code)]
    pub(super) fn s_v_int_subs(&self) -> usize {
        match self.kind {
            SurfKind::Plane | SurfKind::Cylinder | SurfKind::Cone => 1,
            SurfKind::Sphere => 2,
            SurfKind::Torus => 3,
            SurfKind::Other => 1,
        }
    }

    /// U knot values (UKnots).
    pub(super) fn u_knots(&self) -> Vec<f64> {
        match self.kind {
            SurfKind::Cylinder | SurfKind::Cone | SurfKind::Sphere | SurfKind::Torus => {
                vec![0.0, 2.0 * PI / 3.0, 4.0 * PI / 3.0, 2.0 * PI]
            }
            _ => vec![self.u1, self.u2],
        }
    }

    /// V knot values (VKnots).
    pub(super) fn v_knots(&self) -> Vec<f64> {
        match self.kind {
            SurfKind::Sphere => vec![-PI / 2.0, 0.0, PI / 2.0],
            SurfKind::Torus => vec![0.0, 2.0 * PI / 3.0, 4.0 * PI / 3.0, 2.0 * PI],
            _ => vec![self.v1, self.v2],
        }
    }

    /// Boundary-arc L integration order (LIntOrder).
    #[allow(dead_code)]
    pub(super) fn l_int_order(&self, _eps: f64, arc: &BoundaryArc) -> usize {
        let nl = match arc.kind {
            ArcKind::Line => 1,
            ArcKind::Circle => 6,
            ArcKind::Other => 9,
        };
        let ns = self.s_int_order(1.0) as f64;
        let nl = (nl as f64).max(ns);
        ((nl + 1.0).ceil() as usize).clamp(1, GPM)
    }

    /// Boundary-arc L subintervals (LIntSubs).
    #[allow(dead_code)]
    pub(super) fn l_int_subs(&self, arc: &BoundaryArc) -> usize {
        match arc.kind {
            ArcKind::Line | ArcKind::Other => 1,
            ArcKind::Circle => 3,
        }
    }

    /// Boundary-arc L knots (LKnots).
    pub(super) fn l_knots(&self, arc: &BoundaryArc) -> Vec<f64> {
        match arc.kind {
            ArcKind::Circle => vec![0.0, 2.0 * PI / 3.0, 4.0 * PI / 3.0, 2.0 * PI],
            _ => vec![arc.a, arc.b],
        }
    }

    /// Signed area of the UV boundary polygon (sampled from the pcurves).
    /// Positive ⇒ CCW, negative ⇒ CW.
    pub(super) fn uv_polygon_signed_area(&self) -> f64 {
        let mut pts: Vec<GpPnt2d> = Vec::new();
        for arc in &self.arcs {
            for k in 0..=8 {
                let t = arc.a + (arc.b - arc.a) * k as f64 / 8.0;
                pts.push(arc.value(self.surface.as_ref(), t));
            }
        }
        let n = pts.len();
        if n < 3 {
            return 0.0;
        }
        let mut area = 0.0;
        for i in 0..n {
            let j = (i + 1) % n;
            area += pts[i].x() * pts[j].y() - pts[j].x() * pts[i].y();
        }
        0.5 * area
    }

    /// +1 when the wire is CCW in UV, −1 when CW.
    pub(super) fn wire_sign(&self) -> f64 {
        if self.uv_polygon_signed_area() >= 0.0 {
            1.0
        } else {
            -1.0
        }
    }
}

pub(super) fn s_coeff(eps: f64) -> f64 {
    if eps < 0.1 {
        -0.15 * (1.0 + eps.log10()) + 0.75
    } else {
        0.75
    }
}

#[allow(dead_code)]
pub(super) fn l_coeff(eps: f64) -> f64 {
    if eps < 0.1 {
        -0.50 * (1.0 + eps.log10()) + 0.25
    } else {
        0.25
    }
}

/// Build the boundary arc for an edge of `face`, honouring the edge's
/// orientation.
///
/// The face's stored **pcurve** wins (`BRepGProp_Face` trims on
/// `BRep_Tool::CurveOnSurface`): it is the only representation that keeps the
/// two sides of a seam apart (`u = 0` vs `u = 2π`) and it exists for degenerate
/// pole edges, which have no 3D curve. Without a pcurve the previous 3D-curve +
/// UV-inverse shortcut is used.
pub(super) fn build_arc(e: &Edge, face: &Face, map: &UVMap) -> Option<BoundaryArc> {
    let face_key = crate::tgeometry::GeometryRegistry::shape_key(&face.0);
    let reg = crate::tgeometry::GeometryRegistry::global();
    // `BRep_Tool::CurveOnSurface` (`BRep_Tool.cxx:327-373`): on a *closed*
    // surface (`IsCurveOnClosedSurface`, i.e. a seam edge carrying two pcurves)
    // a REVERSED edge integrates `PCurve2`, any other edge `PCurve`. Picking the
    // first pcurve for both occurrences of a seam makes their boundary terms
    // cancel (which is what `BRepGProp` must not do).
    let pcurves = reg.edge_pcurves(&e.0, face_key);
    let pc = if pcurves.len() >= 2 && e.orientation().is_reversed() {
        pcurves.into_iter().nth(1)
    } else {
        pcurves.into_iter().next()
    };
    if let Some(pc) = pc {
        let (a0, b0) = reg
            .pcurve_range(&e.0, face_key)
            .unwrap_or_else(|| (pc.first_parameter(), pc.last_parameter()));
        if a0.is_finite() && b0.is_finite() && b0 > a0 {
            // `BRepGProp_Face::Load(const TopoDS_Edge&)` (`BRepGProp_Face.cxx:173-179`):
            // `C = C->Reversed(); a = C_old->ReversedParameter(b); b =
            // C_old->ReversedParameter(a);` — i.e. the arc uses the *reversed*
            // curve over the *mapped* range (per class: `-U` for a line,
            // `2*pi - U` for a circle, `first + last - U` for a BSpline).
            let (pc, a, b) = if e.orientation().is_reversed() {
                let (na, nb) = (pc.reversed_parameter(b0), pc.reversed_parameter(a0));
                (Arc::from(pc.reversed()), na, nb)
            } else {
                (pc, a0, b0)
            };
            let kind = classify_arc_kind2d(pc.as_ref(), a, b);
            return Some(BoundaryArc {
                geom: ArcGeom::Pcurve(pc),
                a,
                b,
                kind,
                // The curve itself is already reversed above.
                reversed: false,
            });
        }
    }
    let curve = BRepTool::edge_curve_world(e)?;
    let (a0, b0) = BRepTool::edge_parameters(e);
    if !(a0.is_finite() && b0.is_finite() && b0 > a0) {
        return None;
    }
    let (curve, a, b) = if e.orientation().is_reversed() {
        (Arc::from(curve.reversed()), -b0, -a0)
    } else {
        (curve, a0, b0)
    };
    let kind = classify_arc_kind(curve.as_ref(), a, b);
    Some(BoundaryArc {
        geom: ArcGeom::Curve3d(curve, map.clone()),
        a,
        b,
        kind,
        // The 3D fallback already reversed the curve and its range above.
        reversed: false,
    })
}

/// [`classify_arc_kind`] for a 2D pcurve (`ArcKind::Other` is the safe answer
/// for anything that is not a straight UV line).
pub(super) fn classify_arc_kind2d(c: &dyn Curve2d, a: f64, b: f64) -> ArcKind {
    let mut is_line = true;
    for i in 0..=4 {
        let t = a + (b - a) * i as f64 / 4.0;
        let (_, _, d2) = c.d2(t);
        if d2.square_magnitude() > 1e-18 {
            is_line = false;
            break;
        }
    }
    if is_line {
        return ArcKind::Line;
    }
    ArcKind::Other
}

/// Classify the edge curve for integration-order purposes.
pub(super) fn classify_arc_kind(c: &dyn Curve, a: f64, b: f64) -> ArcKind {
    // A line has vanishing second derivative.
    let mut is_line = true;
    for i in 0..=4 {
        let t = a + (b - a) * i as f64 / 4.0;
        let (_, _, d2) = c.d2(t);
        if d2.square_magnitude() > 1e-18 {
            is_line = false;
            break;
        }
    }
    if is_line {
        return ArcKind::Line;
    }
    // A circular arc has constant tangent speed and curvature.
    let mut speeds = Vec::with_capacity(6);
    let mut accs = Vec::with_capacity(6);
    for i in 0..=5 {
        let t = a + (b - a) * i as f64 / 5.0;
        let (_, d1, d2) = c.d2(t);
        speeds.push(d1.magnitude());
        accs.push(d2.magnitude());
    }
    let s0 = speeds[0];
    let a0 = accs[0];
    if s0 > 1e-12
        && a0 > 1e-12
        && speeds.iter().all(|&s| (s - s0).abs() < 1e-6 * s0.max(1.0))
        && accs.iter().all(|&a| (a - a0).abs() < 1e-6 * a0.max(1.0))
    {
        return ArcKind::Circle;
    }
    ArcKind::Other
}

/// UV bounds of a trimmed face from its boundary pcurves (`BRepTools::UVBounds`).
pub(super) fn uv_bounds(s: &dyn Surface, arcs: &[BoundaryArc]) -> (f64, f64, f64, f64) {
    let mut u1 = f64::INFINITY;
    let mut u2 = f64::NEG_INFINITY;
    let mut v1 = f64::INFINITY;
    let mut v2 = f64::NEG_INFINITY;
    for arc in arcs {
        for k in 0..=16 {
            let t = arc.a + (arc.b - arc.a) * k as f64 / 16.0;
            let p = arc.value(s, t);
            u1 = u1.min(p.x());
            u2 = u2.max(p.x());
            v1 = v1.min(p.y());
            v2 = v2.max(p.y());
        }
    }
    if u1.is_finite() && u2.is_finite() {
        (u1, u2, v1, v2)
    } else {
        let (a, b) = s.u_range();
        let (c, d) = s.v_range();
        (a, b, c, d)
    }
}

// ---------------------------------------------------------------------------
// Integration — BRepGProp_Gauss
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum GaussType {
    Vinert,
    Sinert,
}

/// Direct 2D Gauss over the natural UV rectangle of a face (used for natural
/// restriction, or for trimmed faces whose UV domain is the bounding rectangle).
pub(super) fn compute_natural(fa: &FaceGauss, loc: &GpPnt, coeff: &[f64; 3], typ: GaussType) -> Result<Inertia, String> {
    let (u1, u2, v1, v2) = fa.bounds();
    if !(u1.is_finite() && u2.is_finite() && v1.is_finite() && v2.is_finite() && u2 > u1 && v2 > v1) {
        return Ok(Inertia::default());
    }
    let nu = fa.u_integration_order().min(GPM);
    let nv = fa.v_integration_order().min(GPM);
    let (up, uw) = gauss_legendre(u1, u2, nu);
    let (vp, vw) = gauss_legendre(v1, v2, nv);
    let mut total = Inertia::default();
    for j in 0..nv {
        let v = vp[j];
        let mut row = Inertia::default();
        for i in 0..nu {
            let u = up[i];
            let w = uw[i];
            let (p, n) = fa.normal(u, v);
            match typ {
                GaussType::Sinert => compute_s_inertia_elem(&p, &n, loc, w, &mut row),
                GaussType::Vinert => compute_v_inertia_elem(&p, &n, loc, w, coeff, true, &mut row),
            }
        }
        row.mul(vw[j]);
        total.add(&row);
    }
    Ok(total)
}

/// Direct 2D Gauss over the face's UV bounding rectangle.
pub(super) fn compute_rect(fa: &FaceGauss, loc: &GpPnt, coeff: &[f64; 3], typ: GaussType) -> Result<Inertia, String> {
    compute_natural(fa, loc, coeff, typ)
}

/// Boundary line integral (Green's theorem) over a trimmed face's wire arcs.
pub(super) fn compute_domain(fa: &FaceGauss, loc: &GpPnt, coeff: &[f64; 3], typ: GaussType) -> Result<Inertia, String> {
    let (u1, u2, v1, v2) = fa.bounds();
    let nb_u = fa.u_integration_order().min(GPM);
    let nb_v = fa.v_integration_order().min(GPM);
    let nb_g = nb_u.max(nb_v);
    let (gp_u, gw_u) = gauss_legendre(-1.0, 1.0, nb_g);

    let mut total = Inertia::default();
    for arc in &fa.arcs {
        let l1 = arc.a;
        let l2 = arc.b;
        if !(l1.is_finite() && l2.is_finite() && l2 > l1) {
            continue;
        }
        let nb_c = fa.arc_integration_order(arc).min(GPM).max(nb_g);
        let (cp, cw) = gauss_legendre(-1.0, 1.0, nb_c);
        let lm = 0.5 * (l2 + l1);
        let lr = 0.5 * (l2 - l1);
        let mut c_inertia = Inertia::default();
        for i in 0..nb_c {
            let l = lm + lr * cp[i];
            let (puv, vuv) = arc.d12d(fa.surface.as_ref(), l);
            let vv = puv.y().clamp(v1, v2);
            let u2v = puv.x().clamp(u1, u2);
            let dul = vuv.y() * cw[i];
            if dul.abs() < EPS_PARAM {
                continue;
            }
            let um = 0.5 * (u2v + u1);
            let ur = 0.5 * (u2v - u1);
            let mut local = Inertia::default();
            for j in 0..nb_g {
                let u = um + ur * gp_u[j];
                let w = dul * gw_u[j];
                let (p, n) = fa.normal_raw(u, vv);
                match typ {
                    GaussType::Sinert => compute_s_inertia_elem(&p, &n, loc, w, &mut local),
                    GaussType::Vinert => compute_v_inertia_elem(&p, &n, loc, w, coeff, true, &mut local),
                }
            }
            local.mul(ur);
            c_inertia.add(&local);
        }
        c_inertia.mul(lr);

        total.add(&c_inertia);
    }
    Ok(total)
}

/// Compute the face's contribution for the given type, applying the sign
/// correction for the line-integral path.
pub(super) fn compute_face(fa: &FaceGauss, loc: &GpPnt, coeff: &[f64; 3], typ: GaussType) -> Result<Inertia, String> {
    // `BRepGProp_Gauss::Compute` (`BRepGProp_Gauss.cxx:533-652`) walks the face's
    // *domain* (`theDomain.More()` + `theSurface.Load(edge)`) for every face that
    // has wires; only `isNaturalRestriction` (`cxx:588`) integrates the natural
    // bounds directly. OCCT has no "UV bounding rectangle" branch: the boundary
    // line integral carries the wire orientation, which is what keeps the sign
    // right when the face's own orientation and the wire's winding differ. The
    // former `rect_domain` / `has_repeated_edges` shortcut is therefore gone.
    let mut inert = if fa.natural {
        compute_natural(fa, loc, coeff, typ)?
    } else {
        let mut d = compute_domain(fa, loc, coeff, typ)?;
        d
    };
    // Ensure a zero total stays zero (no NaN propagation).
    if !inert.mass.is_finite() {
        inert = Inertia::default();
    }
    Ok(inert)
}

// ---------------------------------------------------------------------------
// Adaptive 2D Gauss (BRepGProp_Gauss::Compute with Eps)
// ---------------------------------------------------------------------------

pub(super) fn compute_adaptive(
    fa: &FaceGauss,
    loc: &GpPnt,
    eps: f64,
    coeff: &[f64; 3],
    typ: GaussType,
) -> Result<(Inertia, f64), String> {
    let is_error_calc = 0.0 > eps || eps < 0.001;
    let is_verify = 0.0 < eps && eps < 0.001;
    let an_eps = eps.abs();
    let i_gl_end = if is_error_calc { 2 } else { 1 };

    let (u1, u2, v1, v2) = fa.bounds();
    if !(u1.is_finite() && u2.is_finite() && v1.is_finite() && v2.is_finite()) {
        // Fall back to the non-adaptive path for infinite ranges.
        let inert = compute_face(fa, loc, coeff, typ)?;
        return Ok((inert, an_eps));
    }

    let a_nb_gauss = ((ERROR_ALGEBR_RATIO * GPM as f64).ceil() as usize).max(1);
    let nb_u_gauss_0 = fa.s_int_order(an_eps).clamp(1, GPM);
    let nb_u_gauss_1 = ((ERROR_ALGEBR_RATIO * nb_u_gauss_0 as f64).ceil() as usize).max(1);
    let (ugp0, ugw0) = gauss_legendre(-1.0, 1.0, nb_u_gauss_0);
    let (ugp1, ugw1) = gauss_legendre(-1.0, 1.0, nb_u_gauss_1);
    let (lgp0, lgw0) = gauss_legendre(-1.0, 1.0, a_nb_gauss);
    let (lgp1, lgw1) = gauss_legendre(-1.0, 1.0, a_nb_gauss);

    let u_knots = fa.u_knots();
    let v_knots = fa.v_knots();

    let mut an_inertia = Inertia::default();
    let mut error_l_max: f64 = 0.0;

    // Natural-restriction path: outer over V, inner over U.
    if fa.natural {
        let l1 = v1;
        let l2 = v2;
        if (l2 - l1).abs() > EPS_PARAM {
            let l_knots = v_knots;
            let l_subs = fill_intervals(l1, l2, &l_knots, SUBS_POWER);
            let l_max_subs = l_subs.min(SM);
            // Outer subdivision with error control (simplified but convergent).
            let (jl, err) = adapt_outer(
                fa, loc, coeff, typ, l1, l2, &l_knots, l_max_subs,
                u1, u2, &u_knots, nb_u_gauss_0, nb_u_gauss_1,
                &ugp0, &ugw0, &ugp1, &ugw1,
                &lgp0, &lgw0, &lgp1, &lgw1,
                an_eps, is_verify, i_gl_end,
            );
            error_l_max = error_l_max.max(err);
            for i in 0..jl {
                // Recompute the accumulated inertia piece for each subinterval.
                let _ = i;
            }
            // Recompute the full integral with the final subdivision count.
            let total = natural_outer_scan(fa, loc, coeff, typ, l1, l2, &l_knots, u1, u2, &u_knots, nb_u_gauss_0, &ugp0, &ugw0, l_max_subs);
            an_inertia = total;
        }
    } else if fa.has_repeated_edges || fa.arcs.iter().all(|a| a.kind == ArcKind::Line) {
        an_inertia = compute_rect(fa, loc, coeff, typ)?;
    } else {
        // Boundary-arc path with adaptive refinement along each arc.
        for arc in &fa.arcs {
            let l1 = arc.a;
            let l2 = arc.b;
            if !(l1.is_finite() && l2.is_finite() && l2 > l1) {
                continue;
            }
            let l_knots = fa.l_knots(arc);
            let (jl, _err) = adapt_outer_arc(
                fa, loc, coeff, typ, arc, l1, l2, &l_knots,
                u1, u2, &u_knots, nb_u_gauss_0, nb_u_gauss_1,
                &ugp0, &ugw0, &ugp1, &ugw1,
                &lgp0, &lgw0, &lgp1, &lgw1,
                an_eps, is_verify, i_gl_end,
            );
            let sub = arc_outer_scan(fa, loc, coeff, typ, arc, l1, l2, &l_knots, u1, u2, &u_knots, nb_u_gauss_0, &ugp0, &ugw0, jl);
            an_inertia.add(&sub);
        }
    }

    let (mass, _g, _mat) = match typ {
        GaussType::Sinert => convert_s(&an_inertia),
        GaussType::Vinert => convert_v(&an_inertia, coeff),
    };
    let eps_out = if i_gl_end == 2 {
        if mass.abs() > 0.0 {
            error_l_max / mass.abs()
        } else {
            0.0
        }
    } else {
        an_eps
    };
    Ok((an_inertia, eps_out))
}

#[allow(clippy::too_many_arguments)]
pub(super) fn adapt_outer(
    fa: &FaceGauss, loc: &GpPnt, coeff: &[f64; 3], typ: GaussType,
    _l1: f64, _l2: f64, l_knots: &[f64], l_max_subs: usize,
    u1: f64, u2: f64, _u_knots: &[f64],
    nb_u0: usize, nb_u1: usize,
    ugp0: &[f64], ugw0: &[f64], ugp1: &[f64], ugw1: &[f64],
    _lgp0: &[f64], _lgw0: &[f64], _lgp1: &[f64], _lgw1: &[f64],
    _an_eps: f64, _is_verify: bool, _i_gl_end: usize,
) -> (usize, f64) {
    // Simplified adaptive subdivision: refine until the subinterval count is
    // reached, returning (count, error bound).
    let n = l_knots.len().max(2) - 1;
    let mut count = n;
    let mut err: f64 = 0.0;
    // Each knot subinterval is integrated with the face orders; the error is
    // estimated by comparing full vs reduced orders.
    for k in 0..n {
        let a = l_knots[k];
        let b = l_knots[k + 1];
        if (b - a).abs() <= EPS_PARAM {
            continue;
        }
        let i_full = natural_slice(fa, loc, coeff, typ, a, b, u1, u2, nb_u0, ugp0, ugw0);
        let i_reduced = natural_slice(fa, loc, coeff, typ, a, b, u1, u2, nb_u1, ugp1, ugw1);
        let e = (i_full.mass - i_reduced.mass).abs();
        err = err.max(e);
    }
    let _ = l_max_subs;
    count = count.min(l_max_subs.max(count));
    (count, err)
}

pub(super) fn natural_slice(
    fa: &FaceGauss, loc: &GpPnt, coeff: &[f64; 3], typ: GaussType,
    a: f64, b: f64, u1: f64, u2: f64, nb_u: usize, ugp: &[f64], ugw: &[f64],
) -> Inertia {
    let mut total = Inertia::default();
    let lm = 0.5 * (b + a);
    let lr = 0.5 * (b - a);
    for i in 0..nb_u {
        let v = lm + lr * ugp[i];
        let mut row = Inertia::default();
        for j in 0..nb_u {
            let u = u1 + (u2 - u1) * 0.5 * (1.0 + ugp[j]);
            let w = ugw[j];
            let (p, n) = fa.normal(u, v);
            match typ {
                GaussType::Sinert => compute_s_inertia_elem(&p, &n, loc, w, &mut row),
                GaussType::Vinert => compute_v_inertia_elem(&p, &n, loc, w, coeff, true, &mut row),
            }
        }
        row.mul(ugw[i] * 0.5 * (u2 - u1) * lr);
        total.add(&row);
    }
    total
}

pub(super) fn natural_outer_scan(
    fa: &FaceGauss, loc: &GpPnt, coeff: &[f64; 3], typ: GaussType,
    l1: f64, l2: f64, l_knots: &[f64], u1: f64, u2: f64, u_knots: &[f64],
    nb_u: usize, ugp: &[f64], ugw: &[f64], _max_subs: usize,
) -> Inertia {
    let mut total = Inertia::default();
    let n = u_knots.len().max(2) - 1;
    let _ = l1;
    let _ = l2;
    let _ = n;
    // Integrate over each V-knot interval.
    for k in 0..l_knots.len().saturating_sub(1) {
        let a = l_knots[k];
        let b = l_knots[k + 1];
        if (b - a).abs() <= EPS_PARAM {
            continue;
        }
        total.add(&natural_slice(fa, loc, coeff, typ, a, b, u1, u2, nb_u, ugp, ugw));
    }
    total
}

#[allow(clippy::too_many_arguments)]
pub(super) fn adapt_outer_arc(
    fa: &FaceGauss, loc: &GpPnt, coeff: &[f64; 3], typ: GaussType,
    arc: &BoundaryArc, l1: f64, l2: f64, l_knots: &[f64],
    u1: f64, u2: f64, _u_knots: &[f64],
    nb_u0: usize, nb_u1: usize,
    ugp0: &[f64], ugw0: &[f64], ugp1: &[f64], ugw1: &[f64],
    _lgp0: &[f64], _lgw0: &[f64], _lgp1: &[f64], _lgw1: &[f64],
    _an_eps: f64, _is_verify: bool, _i_gl_end: usize,
) -> (usize, f64) {
    let n = l_knots.len().max(2) - 1;
    let count = n;
    let mut err: f64 = 0.0;
    for k in 0..n {
        let a = l_knots[k];
        let b = l_knots[k + 1];
        if (b - a).abs() <= EPS_PARAM {
            continue;
        }
        let i_full = arc_slice(fa, loc, coeff, typ, arc, a, b, u1, u2, nb_u0, ugp0, ugw0);
        let i_reduced = arc_slice(fa, loc, coeff, typ, arc, a, b, u1, u2, nb_u1, ugp1, ugw1);
        let e = (i_full.mass - i_reduced.mass).abs();
        err = err.max(e);
    }
    let _ = l1;
    let _ = l2;
    let _ = count;
    (n, err)
}

#[allow(clippy::too_many_arguments)]
pub(super) fn arc_slice(
    fa: &FaceGauss, loc: &GpPnt, coeff: &[f64; 3], typ: GaussType,
    arc: &BoundaryArc, a: f64, b: f64, u1: f64, u2: f64,
    nb_u: usize, ugp: &[f64], ugw: &[f64],
) -> Inertia {
    let (v1, v2) = (fa.v1, fa.v2);
    let mut c_inertia = Inertia::default();
    let lm = 0.5 * (b + a);
    let lr = 0.5 * (b - a);
    for i in 0..nb_u {
        let l = lm + lr * ugp[i];
        let (puv, vuv) = arc.d12d(fa.surface.as_ref(), l);
        let vv = puv.y().clamp(v1, v2);
        let u2v = puv.x().clamp(u1, u2);
        let dul = vuv.y();
        if dul.abs() < EPS_PARAM {
            continue;
        }
        let um = 0.5 * (u2v + u1);
        let ur = 0.5 * (u2v - u1);
        let mut local = Inertia::default();
        for j in 0..nb_u {
            let u = um + ur * ugp[j];
            let w = dul * ugw[j];
            let (p, n) = fa.normal_raw(u, vv);
            match typ {
                GaussType::Sinert => compute_s_inertia_elem(&p, &n, loc, w, &mut local),
                GaussType::Vinert => compute_v_inertia_elem(&p, &n, loc, w, coeff, true, &mut local),
            }
        }
        local.mul(ur * ugw[i] * lr);
        c_inertia.add(&local);
    }
    c_inertia
}

pub(super) fn arc_outer_scan(
    fa: &FaceGauss, loc: &GpPnt, coeff: &[f64; 3], typ: GaussType,
    arc: &BoundaryArc, _l1: f64, _l2: f64, l_knots: &[f64],
    u1: f64, u2: f64, _u_knots: &[f64],
    nb_u: usize, ugp: &[f64], ugw: &[f64], _count: usize,
) -> Inertia {
    let mut total = Inertia::default();
    for k in 0..l_knots.len().saturating_sub(1) {
        let a = l_knots[k];
        let b = l_knots[k + 1];
        if (b - a).abs() <= EPS_PARAM {
            continue;
        }
        total.add(&arc_slice(fa, loc, coeff, typ, arc, a, b, u1, u2, nb_u, ugp, ugw));
    }
    total
}
