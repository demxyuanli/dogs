//! `BSplCLib::PrepareEval` + `BSplCLib::DN` for 3D poles.
//! Source: `BSplCLib_CurveComputation.pxx:720-830` and `pxx:1085-1134`.

use crate::bspl::bohm::bohm;
use crate::bspl::build_knots::build_knots;
use crate::bspl::knots::pole_index;
use crate::bspl::locate::locate_parameter;
use crate::bspl::plib_rational;
use crate::gp::{GpPnt, GpVec};

fn is_rational_span(weights: &[f64], i1: i32, i2: i32) -> bool {
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

fn build_eval(degree: i32, index: i32, poles: &[GpPnt], weights: Option<&[f64]>) -> (i32, Vec<f64>) {
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

/// `PrepareEval_T` (`pxx:777-830`) then `BSplCLib_DN` (`pxx:1085-1134`).
///
/// `index` is the incoming knot guess (`0` from `Geom_BSplineCurve::EvalDN`).
/// `mults == None` is `BSplCLib::NoMults` (flat knots).
///
/// `n == 0` reproduces `BSplCLib::D0` (`BSplCLib_1.cxx:248-267`): `Bohm` with
/// derivative order 0 leaves the value in slot 0 and `RationalDerivative(…, 0)`
/// divides by the weight, so the point comes back as a vector. This is the arm
/// the periodic evaluation in `GeomBSplineCurve::{d0,d1,d2}` uses.
pub fn dn(
    u: f64,
    n: i32,
    index: i32,
    degree: i32,
    periodic: bool,
    poles: &[GpPnt],
    weights: Option<&[f64]>,
    knots: &[f64],
    mults: Option<&[i32]>,
) -> GpVec {
    if degree < 0 || n < 0 {
        return GpVec::zero();
    }
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
    let (dim, mut local_poles) = if rational {
        build_eval(degree, knot_index, poles, weights)
    } else {
        build_eval(degree, knot_index, poles, None)
    };
    bohm(uu, degree, n, &local_knots, dim, &mut local_poles);
    if rational {
        let v = plib_rational::rational_derivative(degree, n, &local_poles);
        GpVec::new(v[0], v[1], v[2])
    } else if n > degree {
        GpVec::zero()
    } else {
        let off = (n * 3) as usize;
        if off + 2 < local_poles.len() {
            GpVec::new(local_poles[off], local_poles[off + 1], local_poles[off + 2])
        } else {
            GpVec::zero()
        }
    }
}
