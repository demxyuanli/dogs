//! `BSplSLib::D3` (`BSplSLib.cxx:1251-1515`).

use crate::bspl::bohm::{bohm, eval};
use crate::bspl::prepare_eval::prepare_eval;
use crate::bspl::rational_derivative::rational_derivative;
use crate::elib::surface_eval_d3::SurfD3;
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

/// `BSplSLib::D3` (`cxx:1251-1515`).
pub fn surface_d3(
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
) -> SurfD3 {
    let mut prep = prepare_eval(
        u, v, u_index, v_index, u_degree, v_degree, u_rat, v_rat, u_per, v_per, poles, weights,
        u_knots, v_knots, u_mults, v_mults,
    );
    let dim = if prep.rational { 4 } else { 3 };
    let dim2 = dim * (prep.d2 + 1);
    bohm(
        prep.u1,
        prep.d1,
        3,
        &prep.knots1,
        dim2,
        &mut prep.poles,
    );
    let stride = dim2 as usize;
    bohm(prep.u2, prep.d2, 3, &prep.knots2, dim, &mut prep.poles);
    bohm(
        prep.u2,
        prep.d2,
        2,
        &prep.knots2,
        dim,
        &mut prep.poles[stride..],
    );
    if prep.d1 > 1 {
        bohm(
            prep.u2,
            prep.d2,
            1,
            &prep.knots2,
            dim,
            &mut prep.poles[2 * stride..],
        );
    }
    if prep.d1 > 2 {
        eval(
            prep.u2,
            prep.d2,
            &prep.knots2,
            dim,
            &mut prep.poles[3 * stride..],
        );
    }

    if prep.rational {
        // `RationalDerivative(..., All=true)` default (`BSplSLib.hxx:177`).
        let mut ders = [0.0f64; 48];
        rational_derivative(prep.d1, prep.d2, 3, 3, &prep.poles, &mut ders, true);
        return if prep.ufirst {
            // `cxx:1326-1335`
            SurfD3 {
                point: pnt_at(&ders, 0),
                d1u: xyz_at(&ders, 12),
                d1v: xyz_at(&ders, 3),
                d2u: xyz_at(&ders, 24),
                d2v: xyz_at(&ders, 6),
                d2uv: xyz_at(&ders, 15),
                d3u: xyz_at(&ders, 36),
                d3v: xyz_at(&ders, 9),
                d3_uuv: xyz_at(&ders, 27),
                d3_uvv: xyz_at(&ders, 18),
            }
        } else {
            // `cxx:1413-1422`
            SurfD3 {
                point: pnt_at(&ders, 0),
                d1u: xyz_at(&ders, 3),
                d1v: xyz_at(&ders, 12),
                d2u: xyz_at(&ders, 6),
                d2v: xyz_at(&ders, 24),
                d2uv: xyz_at(&ders, 15),
                d3u: xyz_at(&ders, 9),
                d3v: xyz_at(&ders, 36),
                d3_uuv: xyz_at(&ders, 18),
                d3_uvv: xyz_at(&ders, 27),
            }
        };
    }

    let result = &prep.poles;
    let dim2u = stride;
    if prep.ufirst {
        // `cxx:1353-1392`
        let vuu = if u_degree <= 1 {
            None
        } else {
            Some(2 * dim2u)
        };
        let vuuv = if u_degree <= 1 {
            None
        } else {
            Some(2 * dim2u + 3)
        };
        let vvv = if v_degree <= 1 { None } else { Some(6) };
        let vuvv = if v_degree <= 1 {
            None
        } else {
            Some(dim2u + 6)
        };
        let vuuu = if u_degree <= 2 {
            None
        } else {
            Some(3 * dim2u)
        };
        let vvvv = if v_degree <= 2 { None } else { Some(9) };
        let vuv = ((prep.d2 << 1) + prep.d2 + 6) as usize;
        SurfD3 {
            point: pnt_at(result, 0),
            d1u: xyz_at(result, dim2u),
            d1v: xyz_at(result, 3),
            d2u: zero_or(result, vuu),
            d2v: zero_or(result, vvv),
            d2uv: xyz_at(result, vuv),
            d3u: zero_or(result, vuuu),
            d3v: zero_or(result, vvvv),
            d3_uuv: zero_or(result, vuuv),
            d3_uvv: zero_or(result, vuvv),
        }
    } else {
        // `cxx:1440-1479`
        let vuu = if u_degree <= 1 { None } else { Some(6) };
        let vuuv = if u_degree <= 1 {
            None
        } else {
            Some(dim2u + 6)
        };
        let vvv = if v_degree <= 1 {
            None
        } else {
            Some(2 * dim2u)
        };
        let vuvv = if v_degree <= 1 {
            None
        } else {
            Some(2 * dim2u + 3)
        };
        let vuuu = if u_degree <= 2 { None } else { Some(9) };
        let vvvv = if v_degree <= 2 {
            None
        } else {
            Some(3 * dim2u)
        };
        let vuv = ((prep.d2 << 1) + prep.d2 + 6) as usize;
        SurfD3 {
            point: pnt_at(result, 0),
            d1u: xyz_at(result, 3),
            d1v: xyz_at(result, dim2u),
            d2u: zero_or(result, vuu),
            d2v: zero_or(result, vvv),
            d2uv: xyz_at(result, vuv),
            d3u: zero_or(result, vuuu),
            d3v: zero_or(result, vvvv),
            d3_uuv: zero_or(result, vuuv),
            d3_uvv: zero_or(result, vuvv),
        }
    }
}
