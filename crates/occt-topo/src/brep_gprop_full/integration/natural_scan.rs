use super::*;

pub(in crate::brep_gprop_full) fn natural_slice(
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

pub(in crate::brep_gprop_full) fn natural_outer_scan(
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
