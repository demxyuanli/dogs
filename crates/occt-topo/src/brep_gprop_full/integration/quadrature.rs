use super::*;

// ---------------------------------------------------------------------------
// Integration — BRepGProp_Gauss
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(in crate::brep_gprop_full) enum GaussType {
    Vinert,
    Sinert,
}

/// Direct 2D Gauss over the natural UV rectangle of a face (used for natural
/// restriction, or for trimmed faces whose UV domain is the bounding rectangle).
pub(in crate::brep_gprop_full) fn compute_natural(fa: &FaceGauss, loc: &GpPnt, coeff: &[f64; 3], typ: GaussType) -> Result<Inertia, String> {
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
pub(in crate::brep_gprop_full) fn compute_rect(fa: &FaceGauss, loc: &GpPnt, coeff: &[f64; 3], typ: GaussType) -> Result<Inertia, String> {
    compute_natural(fa, loc, coeff, typ)
}

/// Boundary line integral (Green's theorem) over a trimmed face's wire arcs.
/// `BRepGProp_Face::UKnots` (`BRepGProp_Face.cxx:343-390`): the surface's U iso
/// knots `BRepGProp_Gauss::Compute` splits the U integration at. A plane (and any
/// type not listed) reports its own U bounds, giving a single sub-interval;
/// cylinders, cones, spheres and tori report `{0, 2*pi/3, 4*pi/3, 2*pi}`
/// (`SUIntSubs()` = 4, `BRepGProp_Face.cxx:273-330`). Trailing `NaN`s never
/// compare inside `(u1, u2)`.
fn u_iso_knots(fa: &FaceGauss) -> [f64; 4] {
    use crate::brep_surface::SurfaceKind;
    let (u1, u2) = fa.surface.u_range();
    match crate::brep_surface::classify_surface(fa.surface.as_ref()) {
        SurfaceKind::Cylinder | SurfaceKind::Cone | SurfaceKind::Sphere | SurfaceKind::Torus => [
            0.0,
            2.0 * std::f64::consts::PI / 3.0,
            4.0 * std::f64::consts::PI / 3.0,
            2.0 * std::f64::consts::PI,
        ],
        _ => [u1, u2, f64::NAN, f64::NAN],
    }
}

pub(in crate::brep_gprop_full) fn compute_domain(fa: &FaceGauss, loc: &GpPnt, coeff: &[f64; 3], typ: GaussType) -> Result<Inertia, String> {
    // `BRepGProp_Face::Bounds` (`BRepGProp_Face.cxx:150-156`) gives `u1 = BU1` and
    // the `u2`/`v` clamps from the *surface* parameters (natural range for a
    // cylinder/cone/sphere/torus, face UV bounds only for a plane), not from the
    // face's own UV bounds.
    let (u1, u2, v1, v2) = fa.surface_ranges();
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
            // `BRepGProp_Gauss::Compute` splits the U direction at the surface's
            // iso knots — `FillIntervalBounds(u1, u2, UKnots, NumSubs, anInertiaU,
            // U1, U2, ErrU, aDummy)` (`BRepGProp_Gauss.cxx:772-773`) with `UKnots`
            // from `BRepGProp_Face::UKnots` (`BRepGProp_Face.cxx:343-390`) — and
            // integrates each sub-interval with its own Gauss rule (`:829-871`
            // accumulates `anUI.Mass = mult(aLocal[0].Mass, ur)` per sub-interval).
            let mut pieces: Vec<(f64, f64)> = Vec::new();
            let mut run = u1;
            for k in u_iso_knots(fa) {
                if k > u1 + EPS_PARAM && k < u2v - EPS_PARAM {
                    pieces.push((run, k));
                    run = k;
                }
            }
            pieces.push((run, u2v));
            for (su1, su2) in pieces {
                let um = 0.5 * (su2 + su1);
                let ur = 0.5 * (su2 - su1);
                if ur.abs() < EPS_PARAM {
                    continue;
                }
                let mut local = Inertia::default();
                for j in 0..nb_g {
                    let u = um + ur * gp_u[j];
                    let w = dul * gw_u[j];
                    let (p, n) = fa.normal(u, vv);
                    match typ {
                        GaussType::Sinert => compute_s_inertia_elem(&p, &n, loc, w, &mut local),
                        GaussType::Vinert => compute_v_inertia_elem(&p, &n, loc, w, coeff, true, &mut local),
                    }
                }
                local.mul(ur);
                c_inertia.add(&local);
            }
        }
        c_inertia.mul(lr);

        total.add(&c_inertia);
    }
    Ok(total)
}

/// Compute the face's contribution for the given type, applying the sign
/// correction for the line-integral path.
pub(in crate::brep_gprop_full) fn compute_face(fa: &FaceGauss, loc: &GpPnt, coeff: &[f64; 3], typ: GaussType) -> Result<Inertia, String> {
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
        let d = compute_domain(fa, loc, coeff, typ)?;
        d
    };
    // Ensure a zero total stays zero (no NaN propagation).
    if !inert.mass.is_finite() {
        inert = Inertia::default();
    }
    Ok(inert)
}
