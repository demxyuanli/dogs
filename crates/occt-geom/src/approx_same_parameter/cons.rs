//! `ConsValue` / `ConsD1` and `ComputeTolReached`.

use super::*;

pub(crate) fn cons_value(c2d: &dyn Curve2d, s: &dyn Surface, t: f64) -> GpPnt {
    let uv = c2d.d0(t);
    s.d0(uv.x(), uv.y())
}

pub(crate) fn cons_d1(c2d: &dyn Curve2d, s: &dyn Surface, t: f64) -> (GpPnt, GpVec) {
    let (uv, duv) = c2d.d1(t);
    let (p, du, dv) = s.d1(uv.x(), uv.y());
    let d = du.multiplied_scalar(duv.x()).added(&dv.multiplied_scalar(duv.y()));
    (p, d)
}

pub(crate) fn compute_tol_reached(c3d: &dyn Curve, c2d: &dyn Curve2d, s: &dyn Surface, first: f64, last: f64) -> f64 {
    let nbp = 2 * NB_SAMPLES;
    let mut d2: f64 = 0.0;
    for i in 0..=nbp {
        let t = i as f64 / nbp as f64;
        let u = first * (1.0 - t) + last * t;
        let pc3d = c3d.d0(u);
        let pcons = cons_value(c2d, s, u);
        if !pcons.x().is_finite() || !pcons.y().is_finite() || !pcons.z().is_finite() {
            return INFINITE;
        }
        d2 = d2.max(pc3d.square_distance(&pcons));
    }
    (1.05 * d2.sqrt()).max(CONFUSION)
}

