//! `PLib::HermiteCoefficients` — the Hermite polynomial basis coefficients
//! `GeomFill_PolynomialConvertor` multiplies its constraint vectors with.
//!
//! Source: `FoundationClasses/TKMath/PLib/PLib.cxx:1399-1477` (OCCT 8.0.0).
//!
//! Translation notes:
//! - `math_Matrix MatrixCoefs` becomes one `Vec<f64>` per constraint row, in
//!   the same row order (`MatrixCoefs.SetRow(ii, Coeff)`, `cxx:1473`).
//! - `math_Gauss ResolCoeff(MAT, 1.0e-10)` (`cxx:1464`) is [`Gauss`] with the
//!   same `MinPivot`, so a singular `MAT` yields `false` exactly as
//!   `!ResolCoeff.IsDone()` does.
//! - `Coeff`, `B` and `MAT` are 1-based in OCCT; the port keeps the same
//!   indices in [`MathVector`] / [`MathMatrix`] handles.

use crate::newton_minimum::Gauss;
use crate::{MathMatrix, MathVector};

/// `PLib::HermiteCoefficients(FirstParameter, LastParameter, FirstOrder,
/// LastOrder, MatrixCoefs)` (`PLib.cxx:1404-1477`).
///
/// Returns the `FirstOrder + LastOrder + 2` coefficient rows, or `None` when
/// OCCT returns `false` (invalid orders, out-of-range parameters, degenerate
/// interval, or a singular constraint matrix).
pub fn hermite_coefficients(
    first_parameter: f64,
    last_parameter: f64,
    first_order: i32,
    last_order: i32,
) -> Option<Vec<Vec<f64>>> {
    let nb_coeff = (first_order + last_order + 2) as usize;
    if first_order < 0 || last_order < 0 {
        return None;
    }
    let d1 = first_parameter.abs();
    let mut d2 = last_parameter.abs();
    if d1 > 100.0 || d2 > 100.0 {
        return None;
    }
    d2 += d1;
    if d2 < 0.01 {
        return None;
    }
    if (last_parameter - first_parameter).abs() / d2 < 0.01 {
        return None;
    }

    // Calcul de la matrice a inverser (MAT).
    let ordre = [first_order + 1, last_order + 1];
    let mut t_borne = first_parameter;
    let mut iof = 0usize;
    let mut coeff = MathVector::with_init(1, nb_coeff, 1.0);
    let mut mat = MathMatrix::with_init(1, nb_coeff, 1, nb_coeff, 0.0);

    for cote in 0..=1 {
        coeff.init(1.0);
        for pp in 1..=ordre[cote] as usize {
            let ii = pp + iof;
            let mut prod = 1.0;
            for jj in pp..=nb_coeff {
                // All the work is in the next three lines (`cxx:1453-1456`).
                mat.set_value(ii, jj, coeff.value(jj) * prod);
                coeff.set_value(jj, coeff.value(jj) * (jj - pp) as f64);
                prod *= t_borne;
            }
        }
        t_borne = last_parameter;
        iof = ordre[0] as usize;
    }

    // Resolution du systeme.
    let resol = Gauss::new(&mat, 1.0e-10);
    if !resol.is_done() {
        return None;
    }

    let mut b = MathVector::with_init(1, nb_coeff, 0.0);
    let mut rows = Vec::with_capacity(nb_coeff);
    for ii in 1..=nb_coeff {
        b.set_value(ii, 1.0);
        let row = resol.solve(&b)?;
        rows.push((1..=nb_coeff).map(|j| row.value(j)).collect());
        b.set_value(ii, 0.0);
    }
    Some(rows)
}
