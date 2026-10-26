//! `PLib::EvalLength` — arc length of a polynomial in `n` dimensions.
//!
//! Source: `FoundationClasses/TKMath/PLib/PLib.cxx:2244-2339` (OCCT 8.0.0).
//!
//! Translation notes:
//! - OCCT takes `double& PolynomialCoeff` (the first coefficient) and indexes it
//!   as a raw array; the port takes the coefficient slice of one element, laid
//!   out degree-major / dimension-minor exactly like the OCCT
//!   `NCollection_Array1<double>` block (`FEmTool_Curve.cxx:400-465`).
//! - The Gauss rule is fetched on `[-1, 1]` and the `C1` / `C2` change of
//!   variable stays where OCCT has it (`cxx:2262-2264`), so the accumulation
//!   order is unchanged.
//! - `math::GaussPoints` / `math::GaussWeights` (`cxx:2258-2261`) are
//!   `crate::gauss::gauss_legendre`.
//! - The tolerance overload (`cxx:2308-2339`, `Tol` + `Error`) is UNPORTED: no
//!   caller in the `FEmTool` / `AppDef` chain uses it.

use crate::gauss::gauss_legendre;

/// `PLib::EvalLength(Degree, Dimension, PolynomialCoeff, U1, U2, Length)`
/// (`PLib.cxx:2244-2305`).
pub fn eval_length(degree: i32, dimension: i32, poly: &[f64], u1: f64, u2: f64) -> f64 {
    let dim = dimension.max(0) as usize;
    let deg = degree.max(0);
    if dim == 0 || poly.len() < (deg as usize + 1) * dim {
        return 0.0;
    }

    // `int NbGaussPoints = 4 * std::min((Degree / 4) + 1, 10);` (`cxx:2256`).
    let nb_gauss_points = 4 * ((deg / 4) + 1).min(10);
    let (gauss_points, gauss_weights) = gauss_legendre(-1.0, 1.0, nb_gauss_points as usize);

    let c1 = (u2 + u1) / 2.0;
    let c2 = (u2 - u1) / 2.0;

    let degdim = (deg * dimension) as usize;
    let mut sum = 0.0;

    // `cxx:2272-2300`: symmetric accumulation over `NbGaussPoints / 2` pairs.
    for j in 0..(nb_gauss_points / 2) as usize {
        let tran = c2 * gauss_points[j];
        let x1 = c1 + tran;
        let x2 = c1 - tran;

        let mut der1 = 0.0;
        let mut der2 = 0.0;
        for idim in 0..dim {
            let mut d1 = deg as f64 * poly[idim + degdim];
            let mut d2 = d1;
            for i in (1..deg).rev() {
                let dd = i as f64 * poly[idim + i as usize * dim];
                d1 = d1 * x1 + dd;
                d2 = d2 * x2 + dd;
            }
            der1 += d1 * d1;
            der2 += d2 * d2;
        }

        sum += gauss_weights[j] * c2 * (der1.sqrt() + der2.sqrt());
    }

    sum
}
