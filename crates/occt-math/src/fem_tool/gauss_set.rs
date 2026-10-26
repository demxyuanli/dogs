//! `math_GaussSetIntegration` — Gauss integration of a vector function.
//!
//! Source: `FoundationClasses/TKMath/math/math_GaussSetIntegration.cxx:27-100`
//! (OCCT 8.0.0).
//!
//! Translation notes:
//! - The OCCT class validates `NbVar == 1` and `Order <= math::GaussPointsMax()`
//!   (`cxx:41-42`). The port reports those as `None` instead of raising
//!   `Standard_NotImplemented`.
//! - `F.Value(Tval, Val)` / `F.Value(Tval, FVal1)` are replaced by a closure that
//!   writes `nb_equations` values; a `false` return aborts the integration with
//!   `None`, as OCCT leaves `Done` false.
//! - The odd-order centre point and the symmetric pairs keep the OCCT order of
//!   accumulation (`cxx:64-97`).

use crate::gauss::gauss_legendre;

/// `math::GaussPointsMax()` (`math.cxx:24-27`).
const GAUSS_POINTS_MAX: usize = 61;

/// `math_GaussSetIntegration(F, Lower, Upper, Order)` (`cxx:27-100`): returns
/// the `nb_equations` values of `F` integrated over `[lower, upper]`.
pub fn gauss_set_integration<F>(
    nb_equations: usize,
    order: usize,
    lower: f64,
    upper: f64,
    f: F,
) -> Option<Vec<f64>>
where
    F: Fn(f64, &mut [f64]) -> bool,
{
    // `cxx:41-42`: `NbVar != 1 || Order > math::GaussPointsMax()`.
    if order == 0 || order > GAUSS_POINTS_MAX {
        return None;
    }

    let (gauss_p, gauss_w) = gauss_legendre(-1.0, 1.0, order);

    // `cxx:60-61`: change of variable to [Lower, Upper].
    let xm = 0.5 * (lower + upper);
    let xr = 0.5 * (upper - lower);

    let mut val = vec![0.0f64; nb_equations];
    let mut f_val1 = vec![0.0f64; nb_equations];
    let mut f_val2 = vec![0.0f64; nb_equations];

    // `ind = Ordre / 2`, `ind1 = (Ordre + 1) / 2` in 1-based OCCT terms.
    let ind = order / 2;
    let ind1 = (order + 1) / 2;

    if ind1 > ind {
        // Odd order: the centre point is evaluated once, without `Xr`.
        if !f(xm, &mut val) {
            return None;
        }
        let w = gauss_w[ind1 - 1];
        for v in val.iter_mut() {
            *v *= w;
        }
    }

    for i in 1..=ind {
        if !f(xm + xr * gauss_p[i - 1], &mut f_val1) {
            return None;
        }
        if !f(xm - xr * gauss_p[i - 1], &mut f_val2) {
            return None;
        }
        let w = gauss_w[i - 1];
        for k in 0..nb_equations {
            val[k] += (f_val1[k] + f_val2[k]) * w;
        }
    }

    for v in val.iter_mut() {
        *v *= xr;
    }
    Some(val)
}
