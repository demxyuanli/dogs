use super::*;

#[allow(clippy::too_many_arguments)]
pub(in crate::brep_gprop_full) fn adapt_outer_arc(
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
pub(in crate::brep_gprop_full) fn arc_slice(
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
            let (p, n) = fa.normal(u, vv);
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

pub(in crate::brep_gprop_full) fn arc_outer_scan(
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
