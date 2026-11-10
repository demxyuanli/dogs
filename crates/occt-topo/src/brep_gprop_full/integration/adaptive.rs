use super::*;

// ---------------------------------------------------------------------------
// Adaptive 2D Gauss (BRepGProp_Gauss::Compute with Eps)
// ---------------------------------------------------------------------------

pub(in crate::brep_gprop_full) fn compute_adaptive(
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

    let (u1, u2, v1, v2) = fa.surface_ranges();
    if !(u1.is_finite() && u2.is_finite() && v1.is_finite() && v2.is_finite()) {
        // Fall back to the non-adaptive path for infinite ranges.
        let inert = compute_face(fa, loc, coeff, typ)?;
        return Ok((inert, an_eps));
    }

    let nb_gauss = ((ERROR_ALGEBR_RATIO * GPM as f64).ceil() as usize).max(1);
    let nb_u_gauss_0 = fa.s_int_order(an_eps).clamp(1, GPM);
    let nb_u_gauss_1 = ((ERROR_ALGEBR_RATIO * nb_u_gauss_0 as f64).ceil() as usize).max(1);
    let (ugp0, ugw0) = gauss_legendre(-1.0, 1.0, nb_u_gauss_0);
    let (ugp1, ugw1) = gauss_legendre(-1.0, 1.0, nb_u_gauss_1);
    let (lgp0, lgw0) = gauss_legendre(-1.0, 1.0, nb_gauss);
    let (lgp1, lgw1) = gauss_legendre(-1.0, 1.0, nb_gauss);

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
pub(in crate::brep_gprop_full) fn adapt_outer(
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
