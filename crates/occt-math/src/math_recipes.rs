//! `DACTCL_Decompose` / `DACTCL_Solve` — the symmetric banded solver that
//! `AppParCurves_LeastSquare` drives its normal equations with.
//!
//! Source: `FoundationClasses/TKMath/math/math_Recipes.cxx:750-901`
//! (declarations `math_Recipes.hxx:88-118`, OCCT 8.0.0). Only these two
//! entries of `math_Recipes` are ported here; the rest of that file (SVD,
//! Jacobi, `Math_*` helpers) is still UNPORTED.
//!
//! Translation notes:
//! - `a` stores the upper triangle of the banded symmetric matrix in the
//!   packed order produced by `AppParCurves_LeastSquare::MakeTAA`, and `indx`
//!   holds the position of each diagonal element in that packing. Both keep
//!   OCCT's 1-based indexing; the integer locals are `isize` so that the
//!   address arithmetic is the same signed arithmetic OCCT performs, and a
//!   computed position below 1 panics like `Standard_OutOfRange`.
//! - OCCT returns `math_Status_SingularMatrix` when a pivot falls below
//!   `MinPivot` (`cxx:818-821`, `cxx:875-878`); [`MathStatus`](crate::MathStatus)
//!   has no such code, so [`MathStatus::FunctionError`](crate::MathStatus::FunctionError)
//!   is returned and the caller's `if (!Error)` test keeps OCCT's meaning.
//! - The default `MinPivot` is `1.e-20` (`math_Recipes.hxx:100`).

use crate::{MathIntVector, MathStatus, MathVector};

/// 1-based position from signed address arithmetic (`Standard_OutOfRange`
/// when OCCT would index out of the array).
fn pos(v: isize) -> usize {
    if v < 1 {
        panic!("Standard_OutOfRange: DACTCL");
    }
    v as usize
}

/// `DACTCL_Decompose(a, indx, MinPivot)` (`math_Recipes.cxx:750-830`).
///
/// Decomposes the packed symmetric banded matrix `a` in place, in the order
/// given by `indx`.
pub fn dactcl_decompose(a: &mut MathVector, indx: &MathIntVector, min_pivot: f64) -> MathStatus {
    let neq = indx.len() as isize;
    let mut jr: isize = 0;

    for j in 1..=neq {
        let mut diag = false;
        let jd = indx.value(j as usize) as isize;
        let jh = jd - jr;
        let is = j - jh + 2;
        if jh - 2 == 0 {
            diag = true;
        }
        if jh - 2 > 0 {
            let ie = j - 1;
            let mut k = jr + 2;
            let mut id = indx.value(pos(is - 1)) as isize;
            // Reduction des coefficients non diagonaux:
            for i in is..=ie {
                let ir = id;
                id = indx.value(i as usize) as isize;
                let mut ih = id - ir - 1;
                let mh = i - is + 1;
                if ih > mh {
                    ih = mh;
                }
                if ih > 0 {
                    let mut dot = 0.0;
                    let idot1 = k - ih - 1;
                    let idot2 = id - ih - 1;
                    for idot in 1..=ih {
                        dot += a.value(pos(idot1 + idot)) * a.value(pos(idot2 + idot));
                    }
                    let ak = a.value(pos(k));
                    a.set_value(pos(k), ak - dot);
                }
                k += 1;
            }
            diag = true;
        }

        if diag {
            // Reduction des coefficients diagonaux:
            let ir = jr + 1;
            let ie = jd - 1;
            let k = j - jd;
            for i in ir..=ie {
                let id = indx.value(pos(k + i)) as usize;
                let mut aa = a.value(id);
                if aa < 0.0 {
                    aa = -aa;
                }
                if aa <= min_pivot {
                    return MathStatus::FunctionError;
                }
                let d = a.value(i as usize);
                a.set_value(i as usize, d / a.value(id));
                let jdv = a.value(jd as usize) - d * a.value(i as usize);
                a.set_value(jd as usize, jdv);
            }
        }
        jr = jd;
    }
    MathStatus::Ok
}

/// `DACTCL_Solve(a, b, indx, MinPivot)` (`math_Recipes.cxx:832-901`).
///
/// Solves `a * x = b` with the decomposition left by [`dactcl_decompose`];
/// `b` is replaced by the solution.
pub fn dactcl_solve(
    a: &MathVector,
    b: &mut MathVector,
    indx: &MathIntVector,
    min_pivot: f64,
) -> MathStatus {
    let neq = indx.len() as isize;
    let mut jr: isize = 0;

    for j in 1..=neq {
        let jd = indx.value(j as usize) as isize;
        let jh = jd - jr;
        let is = j - jh + 2;

        // Reduction du second membre:
        let mut dot = 0.0;
        let idot1 = jr;
        let idot2 = is - 2;
        let jh1 = jh - 1;
        for idot in 1..=jh1 {
            dot += a.value(pos(idot1 + idot)) * b.value(pos(idot2 + idot));
        }
        b.set_value(j as usize, b.value(j as usize) - dot);

        jr = jd;
    }

    // Division par les pivots diagonaux:
    for i in 1..=neq {
        let id = indx.value(i as usize) as usize;
        let mut aa = a.value(id);
        if aa < 0.0 {
            aa = -aa;
        }
        if aa <= min_pivot {
            return MathStatus::FunctionError;
        }
        b.set_value(i as usize, b.value(i as usize) / a.value(id));
    }

    // Substitution arriere:
    let mut jd = indx.value(neq as usize) as isize;
    for j in (1..neq).rev() {
        let d = b.value((j + 1) as usize);
        jr = indx.value(j as usize) as isize;
        if jd - jr > 1 {
            let is = j - jd + jr + 2;
            let k = jr - is + 1;
            for i in is..=j {
                let v = b.value(i as usize) - a.value(pos(i + k)) * d;
                b.set_value(i as usize, v);
            }
        }
        jd = jr;
    }
    MathStatus::Ok
}
