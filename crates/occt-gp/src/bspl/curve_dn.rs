//! `BSplCLib::PrepareEval` + `BSplCLib::DN` for 3D poles.
//! Source: `BSplCLib_CurveComputation.pxx:720-830` and `pxx:1085-1134`.

use crate::bspl::bohm::bohm;
use crate::bspl::curve_eval::prepare_curve_eval;
use crate::bspl::plib_rational;
use crate::gp::{GpPnt, GpVec};

/// `PrepareEval_T` (`pxx:777-830`) then `BSplCLib_DN` (`pxx:1085-1134`).
///
/// `index` is the incoming knot guess (`0` from `Geom_BSplineCurve::EvalDN`).
/// `mults == None` is `BSplCLib::NoMults` (flat knots).
///
/// `n == 0` reproduces `BSplCLib::D0` (`BSplCLib_1.cxx:248-267`): `Bohm` with
/// derivative order 0 leaves the value in slot 0 and `RationalDerivative(…, 0)`
/// divides by the weight, so the point comes back as a vector. This is the arm
/// the periodic evaluation in `GeomBSplineCurve::{d0,d1,d2}` uses.
pub fn curve_dn(
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
    let mut p = prepare_curve_eval(u, index, degree, periodic, poles, weights, knots, mults);
    bohm(p.u, degree, n, &p.knots, p.dim, &mut p.poles);
    if p.rational {
        let v = plib_rational::rational_derivative_xyz(degree, n, &p.poles);
        GpVec::new(v[0], v[1], v[2])
    } else if n > degree {
        GpVec::zero()
    } else {
        let off = (n * p.dim) as usize;
        if off + 2 < p.poles.len() {
            GpVec::new(p.poles[off], p.poles[off + 1], p.poles[off + 2])
        } else {
            GpVec::zero()
        }
    }
}
