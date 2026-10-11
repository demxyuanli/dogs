//! `BSplSLib::D2` (`BSplSLib.cxx:1063-1247`).

use crate::bspl::bohm::{bohm, eval};
use crate::bspl::prepare_eval::prepare_eval;
use crate::bspl::rational_derivative::rational_derivative;
use crate::gp::{GpPnt, GpVec};

fn xyz_at(buf: &[f64], off: usize) -> GpVec {
    if off + 2 >= buf.len() {
        return GpVec::zero();
    }
    GpVec::new(buf[off], buf[off + 1], buf[off + 2])
}

fn pnt_at(buf: &[f64], off: usize) -> GpPnt {
    if off + 2 >= buf.len() {
        return GpPnt::zero();
    }
    GpPnt::new(buf[off], buf[off + 1], buf[off + 2])
}

fn zero_or(buf: &[f64], off: Option<usize>) -> GpVec {
    match off {
        Some(o) => xyz_at(buf, o),
        None => GpVec::zero(),
    }
}

/// Second-order surface jet from `BSplSLib::D2`.
pub type SurfD2 = (GpPnt, GpVec, GpVec, GpVec, GpVec, GpVec);

/// `BSplSLib::D2` (`cxx:1063-1247`).
pub fn d2(
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
) -> SurfD2 {
    let mut prep = prepare_eval(
        u, v, u_index, v_index, u_degree, v_degree, u_rat, v_rat, u_per, v_per, poles, weights,
        u_knots, v_knots, u_mults, v_mults,
    );
    let dim = if prep.rational { 4 } else { 3 };
    let dim2 = dim * (prep.d2 + 1);
    bohm(
        prep.u1,
        prep.d1,
        2,
        &prep.knots1,
        dim2,
        &mut prep.poles,
    );
    let stride = dim2 as usize;
    bohm(prep.u2, prep.d2, 2, &prep.knots2, dim, &mut prep.poles);
    bohm(
        prep.u2,
        prep.d2,
        1,
        &prep.knots2,
        dim,
        &mut prep.poles[stride..],
    );
    if prep.d1 > 1 {
        eval(
            prep.u2,
            prep.d2,
            &prep.knots2,
            dim,
            &mut prep.poles[2 * stride..],
        );
    }

    if prep.rational {
        // `RationalDerivative(..., N=2, M=2, All=true)` — ders layout stride
        // `3 * (M+1) = 9` (`cxx:1129-1135` / `:1184-1190`).
        let mut ders = [0.0f64; 27];
        rational_derivative(prep.d1, prep.d2, 2, 2, &prep.poles, &mut ders, true);
        return if prep.ufirst {
            (
                pnt_at(&ders, 0),
                xyz_at(&ders, 9),
                xyz_at(&ders, 3),
                xyz_at(&ders, 18),
                xyz_at(&ders, 6),
                xyz_at(&ders, 12),
            )
        } else {
            (
                pnt_at(&ders, 0),
                xyz_at(&ders, 3),
                xyz_at(&ders, 9),
                xyz_at(&ders, 6),
                xyz_at(&ders, 18),
                xyz_at(&ders, 12),
            )
        };
    }

    let result = &prep.poles;
    let dim2u = stride;
    let vuv = ((prep.d2 << 1) + prep.d2 + 6) as usize;
    if prep.ufirst {
        // `cxx:1149-1168`
        let vuu = if u_degree <= 1 {
            None
        } else {
            Some(2 * dim2u)
        };
        let vvv = if v_degree <= 1 { None } else { Some(6) };
        (
            pnt_at(result, 0),
            xyz_at(result, dim2u),
            xyz_at(result, 3),
            zero_or(result, vuu),
            zero_or(result, vvv),
            xyz_at(result, vuv),
        )
    } else {
        // `cxx:1204-1223`
        let vuu = if u_degree <= 1 { None } else { Some(6) };
        let vvv = if v_degree <= 1 {
            None
        } else {
            Some(2 * dim2u)
        };
        (
            pnt_at(result, 0),
            xyz_at(result, 3),
            xyz_at(result, dim2u),
            zero_or(result, vuu),
            zero_or(result, vvv),
            xyz_at(result, vuv),
        )
    }
}
