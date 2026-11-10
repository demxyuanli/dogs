//! `DistData` and `BuildInitialDistribution` (`cxx:547-648`).

use super::*;

pub(crate) struct DistData {
    pub(crate) pc3d: [f64; MAX_ARRAY],
    pub(crate) pc2d: [f64; MAX_ARRAY],
    pub(crate) nb: usize,
    pub(crate) c3d_f: f64,
    pub(crate) c3d_l: f64,
    pub(crate) c2d_f: f64,
    pub(crate) c2d_l: f64,
    pub(crate) tol: f64,
}

pub(crate) fn project_point_on_curve(
    init: f64,
    point: &GpPnt,
    tolerance: f64,
    curve: &dyn Curve,
    first: f64,
    last: f64,
) -> Option<f64> {
    let mut param = init;
    for _ in 0..30 {
        let (point_val, d1, d2) = curve.d2(param);
        let vector = GpVec::from_pnts(&point_val, point);
        let func = vector.dot(&d1);
        if func.abs() < tolerance * d1.magnitude() {
            return Some(param);
        }
        let func_derivative = vector.dot(&d2) - d1.dot(&d1);
        if func_derivative.abs() > 1.0e-12 {
            param -= func / func_derivative;
        }
        param = param.max(first).min(last);
    }
    None
}

pub(crate) fn increase_initial_nb_samples(data: &mut DistData, c2d: &dyn Curve2d) -> bool {
    let intervals = c2d.parameter_intervals(2);
    let length = intervals.len();
    if length == 0 {
        return true;
    }
    let at = |i1: usize| intervals[i1 - 1];
    let mut inter = 1usize;
    let mut nb_int = length;
    while inter <= nb_int && at(inter) <= data.c3d_f + PCONFUSION {
        inter += 1;
    }
    while nb_int > 0 && at(nb_int) >= data.c3d_l - PCONFUSION {
        nb_int -= 1;
    }
    let mut new_par: Vec<f64> = Vec::new();
    new_par.push(data.c3d_f);
    let mut ii = 1usize;
    while inter <= nb_int || (ii < NB_SAMPLES && inter <= length) {
        if at(inter) < data.pc2d[ii] {
            new_par.push(at(inter));
            if data.pc2d[ii] - at(inter) <= PCONFUSION {
                ii += 1;
                if ii > NB_SAMPLES {
                    ii = NB_SAMPLES;
                }
            }
            inter += 1;
        } else {
            if (at(inter) - data.pc2d[ii]) > PCONFUSION {
                new_par.push(data.pc2d[ii]);
            }
            ii += 1;
        }
    }
    data.nb = new_par.len();
    if data.nb > MAX_ARRAY - 1 {
        return false;
    }
    for i in 1..data.nb {
        let v = new_par[i];
        data.pc2d[i] = v;
        data.pc3d[i] = v;
    }
    data.pc3d[data.nb] = data.c3d_l;
    data.pc2d[data.nb] = data.c2d_l;
    true
}

pub(crate) fn build_initial(data: &mut DistData, c2d: &dyn Curve2d) -> bool {
    let deltacons = (data.c2d_l - data.c2d_f) / NB_SAMPLES as f64;
    let deltac3d = (data.c3d_l - data.c3d_f) / NB_SAMPLES as f64;
    let mut wcons = data.c2d_f;
    let mut wc3d = data.c3d_f;
    for ii in 0..NB_SAMPLES {
        data.pc2d[ii] = wcons;
        data.pc3d[ii] = wc3d;
        wcons += deltacons;
        wc3d += deltac3d;
    }
    data.nb = NB_SAMPLES;
    data.pc2d[data.nb] = data.c2d_l;
    data.pc3d[data.nb] = data.c3d_l;
    if c2d.continuity() < 2 {
        return increase_initial_nb_samples(data, c2d);
    }
    true
}

