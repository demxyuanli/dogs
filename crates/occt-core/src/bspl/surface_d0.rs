//! `BSplSLib::D0` / `HomogeneousD0` (`BSplSLib.cxx:737-851`).

use crate::bspl::bohm::eval;
use crate::bspl::prepare_eval::prepare_eval;
use crate::gp::GpPnt;

/// `BSplSLib::HomogeneousD0` (`cxx:782-851`). Returns `(W, N)` in homogeneous
/// form; for non-rational surfaces `W` is 1 and `N` is the cartesian point.
pub fn homogeneous_d0(
    u: f64,
    v: f64,
    u_index: i32,
    v_index: i32,
    poles: &[Vec<GpPnt>],
    weights: Option<&[Vec<f64>]>,
    u_knots: &[f64],
    v_knots: &[f64],
    u_mults: Option<&[i32]>,
    v_mults: Option<&[i32]>,
    u_degree: i32,
    v_degree: i32,
    u_rat: bool,
    v_rat: bool,
    u_per: bool,
    v_per: bool,
) -> (f64, GpPnt) {
    let mut prep = prepare_eval(
        u, v, u_index, v_index, u_degree, v_degree, u_rat, v_rat, u_per, v_per, poles, weights,
        u_knots, v_knots, u_mults, v_mults,
    );
    let dim = if prep.rational { 4 } else { 3 };
    let dim2 = dim * (prep.d2 + 1);
    eval(prep.u1, prep.d1, &prep.knots1, dim2, &mut prep.poles);
    eval(prep.u2, prep.d2, &prep.knots2, dim, &mut prep.poles);
    if prep.rational {
        let w = prep.poles.get(3).copied().unwrap_or(1.0);
        (
            w,
            GpPnt::new(
                prep.poles.first().copied().unwrap_or(0.0),
                prep.poles.get(1).copied().unwrap_or(0.0),
                prep.poles.get(2).copied().unwrap_or(0.0),
            ),
        )
    } else {
        (
            1.0,
            GpPnt::new(
                prep.poles.first().copied().unwrap_or(0.0),
                prep.poles.get(1).copied().unwrap_or(0.0),
                prep.poles.get(2).copied().unwrap_or(0.0),
            ),
        )
    }
}

/// `BSplSLib::D0` (`cxx:737-778`): `HomogeneousD0` then divide by `W`.
pub fn d0(
    u: f64,
    v: f64,
    u_index: i32,
    v_index: i32,
    poles: &[Vec<GpPnt>],
    weights: Option<&[Vec<f64>]>,
    u_knots: &[f64],
    v_knots: &[f64],
    u_mults: Option<&[i32]>,
    v_mults: Option<&[i32]>,
    u_degree: i32,
    v_degree: i32,
    u_rat: bool,
    v_rat: bool,
    u_per: bool,
    v_per: bool,
) -> GpPnt {
    let (w, p) = homogeneous_d0(
        u, v, u_index, v_index, poles, weights, u_knots, v_knots, u_mults, v_mults, u_degree,
        v_degree, u_rat, v_rat, u_per, v_per,
    );
    if w.abs() < 1e-30 {
        return p;
    }
    GpPnt::new(p.x() / w, p.y() / w, p.z() / w)
}
