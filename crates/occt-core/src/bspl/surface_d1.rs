//! `BSplSLib::D1` / `HomogeneousD1` (`BSplSLib.cxx:855-1059`).

use crate::bspl::bohm::{bohm, eval};
use crate::bspl::prepare_eval::prepare_eval;
use crate::bspl::rational_derivative::rational_derivative;
use crate::gp::{GpPnt, GpVec};

/// Homogeneous first jet: `(N, Nu, Nv, D, Du, Dv)` from `HomogeneousD1`.
pub type HomogSurfD1 = (GpPnt, GpVec, GpVec, f64, f64, f64);

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

/// First-order surface jet from `BSplSLib::D1`.
pub type SurfD1 = (GpPnt, GpVec, GpVec);

/// `BSplSLib::D1` (`cxx:855-970`).
pub fn d1(
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
) -> SurfD1 {
    let mut prep = prepare_eval(
        u, v, u_index, v_index, u_degree, v_degree, u_rat, v_rat, u_per, v_per, poles, weights,
        u_knots, v_knots, u_mults, v_mults,
    );
    let dim = if prep.rational { 4 } else { 3 };
    let dim2 = dim * (prep.d2 + 1);
    bohm(
        prep.u1,
        prep.d1,
        1,
        &prep.knots1,
        dim2,
        &mut prep.poles,
    );
    let stride = dim2 as usize;
    bohm(prep.u2, prep.d2, 1, &prep.knots2, dim, &mut prep.poles);
    eval(
        prep.u2,
        prep.d2,
        &prep.knots2,
        dim,
        &mut prep.poles[stride..],
    );

    if prep.rational {
        // `RationalDerivative(..., N=1, M=1, All=true)` — stride `3*(M+1)=6`.
        let mut ders = [0.0f64; 12];
        rational_derivative(prep.d1, prep.d2, 1, 1, &prep.poles, &mut ders, true);
        return if prep.ufirst {
            // `cxx:914-916`
            (pnt_at(&ders, 0), xyz_at(&ders, 6), xyz_at(&ders, 3))
        } else {
            // `cxx:941-943`
            (pnt_at(&ders, 0), xyz_at(&ders, 3), xyz_at(&ders, 6))
        };
    }

    let result = &prep.poles;
    if prep.ufirst {
        // `cxx:926-928`
        (
            pnt_at(result, 0),
            xyz_at(result, stride),
            xyz_at(result, 3),
        )
    } else {
        // `cxx:953-955`
        (
            pnt_at(result, 0),
            xyz_at(result, 3),
            xyz_at(result, stride),
        )
    }
}

/// `BSplSLib::HomogeneousD1` (`cxx:974-1059`). Returns homogeneous
/// `(N, Nu, Nv, D, Du, Dv)`; for non-rational surfaces `D=1`, `Du=Dv=0`.
pub fn homogeneous_d1(
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
) -> HomogSurfD1 {
    let mut prep = prepare_eval(
        u, v, u_index, v_index, u_degree, v_degree, u_rat, v_rat, u_per, v_per, poles, weights,
        u_knots, v_knots, u_mults, v_mults,
    );
    let dim = if prep.rational { 4 } else { 3 };
    let dim2 = dim * (prep.d2 + 1);
    bohm(
        prep.u1,
        prep.d1,
        1,
        &prep.knots1,
        dim2,
        &mut prep.poles,
    );
    let stride = dim2 as usize;
    bohm(prep.u2, prep.d2, 1, &prep.knots2, dim, &mut prep.poles);
    eval(
        prep.u2,
        prep.d2,
        &prep.knots2,
        dim,
        &mut prep.poles[stride..],
    );
    let result = &prep.poles;
    let (off_u, off_v) = if prep.ufirst {
        (stride, dim as usize)
    } else {
        (dim as usize, stride)
    };
    let n = pnt_at(result, 0);
    let nu = xyz_at(result, off_u);
    let nv = xyz_at(result, off_v);
    if prep.rational {
        let d = result.get(3).copied().unwrap_or(1.0);
        let du = result.get(off_u + 3).copied().unwrap_or(0.0);
        let dv = result.get(off_v + 3).copied().unwrap_or(0.0);
        (n, nu, nv, d, du, dv)
    } else {
        (n, nu, nv, 1.0, 0.0, 0.0)
    }
}
