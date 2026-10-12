//! `BSplCLib` curve evaluation: `PrepareEval_T` + `BSplCLib_D0/D1/D2/D3`.
//! Source: `BSplCLib_CurveComputation.pxx:719-1010`.
//!
//! `Geom_BSplineCurve::EvalD0/EvalD1/EvalD2/EvalD3` (`Geom_BSplineCurve_1.cxx`)
//! dispatch here, so the derivative vectors of a B-spline must come from the
//! same `BSplCLib::Eval` / `BSplCLib::Bohm` arithmetic OCCT uses: evaluating a
//! degree-(d-1)/(d-2) control-point difference curve is mathematically equal
//! but rounds differently, which is enough to drift an edge tessellation.

use crate::bspl::bohm::{bohm, eval};
use crate::bspl::build_knots::build_knots;
use crate::bspl::knots::pole_index;
use crate::bspl::locate::locate_parameter;
use crate::bspl::plib_rational;
use crate::gp::{GpPnt, GpVec};

/// `BSplCLib::IsRational` (`BSplCLib.cxx`).
pub(crate) fn is_rational_span(weights: &[f64], i1: i32, i2: i32) -> bool {
    let f = 1i32;
    let l = weights.len() as i32;
    if l <= 0 {
        return false;
    }
    let i3 = i2 - f;
    let mut i = i1 - f;
    while i < i3 {
        let a = weights[(f + i.rem_euclid(l) - 1) as usize];
        let b = weights[(f + (i + 1).rem_euclid(l) - 1) as usize];
        if a != b {
            return true;
        }
        i += 1;
    }
    false
}

/// `BSplCLib_BuildEval` (`pxx:719-763`): the local pole window (homogeneous
/// `(x, y, z, w)` when the span is rational).
fn build_eval(
    degree: i32,
    index: i32,
    poles: &[GpPnt],
    weights: Option<&[f64]>,
) -> (i32, Vec<f64>) {
    let p_lower = 1i32;
    let p_upper = poles.len() as i32;
    let mut ip = p_lower + index - 1;
    if weights.is_none() {
        let mut lp = vec![0.0; ((degree + 1) * 3).max(0) as usize];
        let mut pole = 0usize;
        for _ in 0..=degree {
            ip += 1;
            if ip > p_upper {
                ip = p_lower;
            }
            let p = poles[(ip - 1) as usize];
            lp[pole] = p.x();
            lp[pole + 1] = p.y();
            lp[pole + 2] = p.z();
            pole += 3;
        }
        (3, lp)
    } else {
        let wts = weights.unwrap();
        let mut lp = vec![0.0; ((degree + 1) * 4).max(0) as usize];
        let mut pole = 0usize;
        for _ in 0..=degree {
            ip += 1;
            if ip > p_upper {
                ip = p_lower;
            }
            let p = poles[(ip - 1) as usize];
            let w = wts[(ip - 1) as usize];
            lp[pole + 3] = w;
            lp[pole] = p.x() * w;
            lp[pole + 1] = p.y() * w;
            lp[pole + 2] = p.z() * w;
            pole += 4;
        }
        (4, lp)
    }
}

/// Local knot window + homogeneous pole window of one span (`PrepareEval_T`).
pub struct CurvePrepared {
    /// `u` after `LocateParameter` (period-wrapped for a periodic curve).
    pub u: f64,
    /// 3, or 4 for a rational span.
    pub dim: i32,
    pub rational: bool,
    pub knots: Vec<f64>,
    pub poles: Vec<f64>,
}

/// `PrepareEval_T` (`pxx:777-830`).
///
/// `index` is the incoming knot guess (see `Geom_BSplineCurve::EvalDx`, which
/// passes the span it located in `myKnots`). `mults == None` is the flat-knot
/// convention (`BSplCLib::NoMults`).
pub fn prepare_eval(
    u: f64,
    index: i32,
    degree: i32,
    periodic: bool,
    poles: &[GpPnt],
    weights: Option<&[f64]>,
    knots: &[f64],
    mults: Option<&[i32]>,
) -> CurvePrepared {
    let (mut knot_index, uu) = locate_parameter(degree, knots, mults, u, periodic, index);
    let local_knots = build_knots(degree, knot_index, periodic, knots, mults);
    if mults.is_none() {
        knot_index -= 1 + degree;
    } else if let Some(m) = mults {
        knot_index = pole_index(degree, knot_index, periodic, m);
    }
    let mut rational = weights.is_some();
    if rational {
        let w_lower = 1 + knot_index;
        rational = is_rational_span(weights.unwrap(), w_lower, w_lower + degree);
    }
    let (dim, poles) = if rational {
        build_eval(degree, knot_index, poles, weights)
    } else {
        build_eval(degree, knot_index, poles, None)
    };
    CurvePrepared {
        u: uu,
        dim,
        rational,
        knots: local_knots,
        poles,
    }
}

/// `BSplCLib::D0` (`pxx:843-882`).
pub fn d0(
    u: f64,
    index: i32,
    degree: i32,
    periodic: bool,
    poles: &[GpPnt],
    weights: Option<&[f64]>,
    knots: &[f64],
    mults: Option<&[i32]>,
) -> GpPnt {
    let mut p = prepare_eval(u, index, degree, periodic, poles, weights, knots, mults);
    eval(p.u, degree, &p.knots, p.dim, &mut p.poles);
    if p.rational {
        let w = p.poles[3];
        GpPnt::new(p.poles[0] / w, p.poles[1] / w, p.poles[2] / w)
    } else {
        GpPnt::new(p.poles[0], p.poles[1], p.poles[2])
    }
}

/// `BSplCLib::D1` (`pxx:896-935`).
pub fn d1(
    u: f64,
    index: i32,
    degree: i32,
    periodic: bool,
    poles: &[GpPnt],
    weights: Option<&[f64]>,
    knots: &[f64],
    mults: Option<&[i32]>,
) -> (GpPnt, GpVec) {
    let mut p = prepare_eval(u, index, degree, periodic, poles, weights, knots, mults);
    bohm(p.u, degree, 1, &p.knots, p.dim, &mut p.poles);
    if p.rational {
        let ders = plib_rational::rational_derivatives(degree, 1, &p.poles);
        (
            GpPnt::new(ders[0], ders[1], ders[2]),
            GpVec::new(ders[3], ders[4], ders[5]),
        )
    } else {
        let dim = p.dim as usize;
        (
            GpPnt::new(p.poles[0], p.poles[1], p.poles[2]),
            GpVec::new(p.poles[dim], p.poles[dim + 1], p.poles[dim + 2]),
        )
    }
}

/// `BSplCLib::D2` (`pxx:948-998`).
pub fn d2(
    u: f64,
    index: i32,
    degree: i32,
    periodic: bool,
    poles: &[GpPnt],
    weights: Option<&[f64]>,
    knots: &[f64],
    mults: Option<&[i32]>,
) -> (GpPnt, GpVec, GpVec) {
    let mut p = prepare_eval(u, index, degree, periodic, poles, weights, knots, mults);
    bohm(p.u, degree, 2, &p.knots, p.dim, &mut p.poles);
    if p.rational {
        let ders = plib_rational::rational_derivatives(degree, 2, &p.poles);
        (
            GpPnt::new(ders[0], ders[1], ders[2]),
            GpVec::new(ders[3], ders[4], ders[5]),
            GpVec::new(ders[6], ders[7], ders[8]),
        )
    } else {
        let dim = p.dim as usize;
        let pnt = GpPnt::new(p.poles[0], p.poles[1], p.poles[2]);
        let d1 = GpVec::new(p.poles[dim], p.poles[dim + 1], p.poles[dim + 2]);
        if degree < 2 {
            return (pnt, d1, GpVec::zero());
        }
        (
            pnt,
            d1,
            GpVec::new(
                p.poles[2 * dim],
                p.poles[2 * dim + 1],
                p.poles[2 * dim + 2],
            ),
        )
    }
}

/// `BSplCLib::D3` (`pxx:1013-1055`).
pub fn d3(
    u: f64,
    index: i32,
    degree: i32,
    periodic: bool,
    poles: &[GpPnt],
    weights: Option<&[f64]>,
    knots: &[f64],
    mults: Option<&[i32]>,
) -> (GpPnt, GpVec, GpVec, GpVec) {
    let mut p = prepare_eval(u, index, degree, periodic, poles, weights, knots, mults);
    bohm(p.u, degree, 3, &p.knots, p.dim, &mut p.poles);
    if p.rational {
        let ders = plib_rational::rational_derivatives(degree, 3, &p.poles);
        (
            GpPnt::new(ders[0], ders[1], ders[2]),
            GpVec::new(ders[3], ders[4], ders[5]),
            GpVec::new(ders[6], ders[7], ders[8]),
            GpVec::new(ders[9], ders[10], ders[11]),
        )
    } else {
        let dim = p.dim as usize;
        let pnt = GpPnt::new(p.poles[0], p.poles[1], p.poles[2]);
        let d1 = GpVec::new(p.poles[dim], p.poles[dim + 1], p.poles[dim + 2]);
        if degree < 2 {
            return (pnt, d1, GpVec::zero(), GpVec::zero());
        }
        let d2 = GpVec::new(
            p.poles[2 * dim],
            p.poles[2 * dim + 1],
            p.poles[2 * dim + 2],
        );
        if degree < 3 {
            return (pnt, d1, d2, GpVec::zero());
        }
        (
            pnt,
            d1,
            d2,
            GpVec::new(
                p.poles[3 * dim],
                p.poles[3 * dim + 1],
                p.poles[3 * dim + 2],
            ),
        )
    }
}
