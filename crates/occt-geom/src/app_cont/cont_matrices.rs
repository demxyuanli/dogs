//! `AppCont_ContMatrices` - the Bernstein mass matrix family used by
//! `AppCont_LeastSquare`.
//!
//! Source: `ModelingData/TKGeomBase/AppCont/AppCont_ContMatrices.pxx:24-180`
//! plus the generated tables `AppCont_ContMatrices_{BB,InvM,IBP,IBT,VB}.pxx`
//! (OCCT 8.0.0).
//!
//! Translation notes:
//! - OCCT stores the five matrices as flat `static const double[]` tables
//!   addressed by the offset helpers `MMatrixOffset` / `IBPMatrixOffset` /
//!   `IBTMatrixOffset` / `Som`. Those tables are pure precomputations of the
//!   Bernstein basis on a Gauss-Legendre rule, so the port evaluates the same
//!   mathematics at run time instead of shipping 120k table lines:
//!   - `MMatrix` is the Bernstein mass matrix
//!     `M(i,j) = C(d,i-1) C(d,j-1) / ((2d+1) C(2d, i+j-2))`, `d = classe-1`
//!     (verified against `AppCont_ContMatrices_BB.pxx` for classe 2..4);
//!   - `InvMMatrix` is `M^-1`;
//!   - `IBPMatrix` is `(MP^T MP)^-1`, `MP` the `M` columns `2..classe-1`;
//!   - `IBTMatrix` is `(MT^T MT)^-1`, `MT` the `M` columns `3..classe-2`;
//!   - `VBernstein` is `B(i,j) = B_{i-1}^{classe-1}(t_j)` with
//!     `t_j = 0.5 * (1 + g_j)`, `g_j` the `j`-th Gauss-Legendre point on
//!     `[-1, 1]` in ascending order (verified against
//!     `AppCont_ContMatrices_VB.pxx` for classe 2).
//! - The legatee-inverse branches at `pxx:480-499` of `AppCont_LeastSquare.cxx`
//!   (`MP.Transposed() * MP` inverted in place) and the table branches are the
//!   same matrices, so the port always uses the closed form and keeps OCCT's
//!   `classe <= 26` / `classe <= 24` range guards.
//! - `Standard_DimensionError` becomes a panic, the port's counterpart of the
//!   OCCT exception (same convention as `AppParCurves_MultiPoint`).

use occt_math::{MathMatrix, gauss_legendre};

/// `THE_MAX_CLASSE` (`pxx:29`).
pub const THE_MAX_CLASSE: i32 = 26;
/// `THE_MAX_MMATRIX_CLASSE` (`pxx:30`).
pub const THE_MAX_MMATRIX_CLASSE: i32 = 24;
/// `THE_MAX_NBPOINTS` (`pxx:31`).
pub const THE_MAX_NBPOINTS: i32 = 24;

/// `C(n, k)` as a `f64`, `0` outside `0 <= k <= n`.
fn binomial(n: usize, k: usize) -> f64 {
    if k > n {
        return 0.0;
    }
    let k = k.min(n - k);
    let mut acc = 1.0f64;
    for i in 0..k {
        acc *= (n - i) as f64;
        acc /= (i + 1) as f64;
    }
    acc
}

/// `B_i^d(t)` - the `i`-th Bernstein basis polynomial of degree `d`.
fn bernstein(d: usize, i: usize, t: f64) -> f64 {
    if i > d {
        return 0.0;
    }
    binomial(d, i) * t.powi(i as i32) * (1.0 - t).powi((d - i) as i32)
}

/// `MMatrix` (`pxx:71-87`): the Bernstein mass matrix, `classe x classe`.
pub fn m_matrix(classe: i32) -> MathMatrix {
    if classe > THE_MAX_MMATRIX_CLASSE {
        panic!("Standard_DimensionError: MMatrix: classe > 24");
    }
    let c = classe as usize;
    let d = classe - 1;
    let den = (2 * d + 1) as f64;
    let mut m = MathMatrix::with_init(1, c, 1, c, 0.0);
    for i in 1..=c {
        for j in 1..=c {
            let num = binomial(d as usize, i - 1) * binomial(d as usize, j - 1);
            let cbin = binomial((2 * d) as usize, i + j - 2);
            m.set_value(i, j, num / (den * cbin));
        }
    }
    m
}

/// `InvMMatrix` (`pxx:90-105`): `M^-1`.
pub fn inv_m_matrix(classe: i32) -> MathMatrix {
    if classe > THE_MAX_MMATRIX_CLASSE {
        panic!("Standard_DimensionError: InvMMatrix: classe > 24");
    }
    m_matrix(classe)
        .inverted()
        .expect("Standard_DimensionError: InvMMatrix singular")
}

/// `(P^T P)^-1` for the `classe x n` slice `P` of `M` spanning columns
/// `lo..=hi` - the `MP` / `MT` of `AppCont_LeastSquare.cxx:487-493`.
///
/// The gram matrix is accumulated directly (`math_Matrix::Transposed` is
/// square-only in this port, and `P` is not), which is the same product.
///
/// The result carries the `lo..=hi` bounds, not `1..=n`. OCCT declares the
/// destination as `math_Matrix IBP(bdeb, bfin, bdeb, bfin)`
/// (`AppCont_LeastSquare.cxx:473`) and `IBPMatrix` / `IBTMatrix` write it at
/// `i + LowerRow() - 1` (`AppCont_ContMatrices.pxx:123-126`, `:144-147`), so
/// the caller then reads `IBP(i, j)` with the same `bdeb..=bfin` indices
/// (`cxx:506`). Re-basing to `1..=n` made that read index past the end for
/// every `classe != 5`.
fn gram_inverse(classe: i32, lo: i32, hi: i32) -> MathMatrix {
    let m = m_matrix(classe);
    let mut gram = MathMatrix::with_init(lo as usize, hi as usize, lo as usize, hi as usize, 0.0);
    for i in lo..=hi {
        for j in lo..=hi {
            let mut sum = 0.0;
            for k in 1..=classe {
                sum += m.value(k as usize, i as usize) * m.value(k as usize, j as usize);
            }
            gram.set_value(i as usize, j as usize, sum);
        }
    }
    gram.invert()
        .expect("Standard_DimensionError: IBPMatrix/IBTMatrix singular");
    gram
}

/// `IBPMatrix` (`pxx:112-127`): `(MP^T MP)^-1` for the pass-point constraints,
/// size `(classe-2) x (classe-2)`.
pub fn ibp_matrix(classe: i32) -> MathMatrix {
    if classe > THE_MAX_CLASSE {
        panic!("Standard_DimensionError: IBPMatrix: classe > 26");
    }
    gram_inverse(classe, 2, classe - 1)
}

/// `IBTMatrix` (`pxx:133-150`): `(MT^T MT)^-1` for the tangency constraints,
/// size `(classe-4) x (classe-4)`.
pub fn ibt_matrix(classe: i32) -> MathMatrix {
    if classe > THE_MAX_CLASSE {
        panic!("Standard_DimensionError: IBTMatrix: classe > 26");
    }
    gram_inverse(classe, 3, classe - 2)
}

/// `VBernstein` (`pxx:152-178`): `B_{i-1}^{classe-1}` evaluated at the
/// `nbpoints` Gauss-Legendre points mapped to `[0, 1]`, `classe x nbpoints`.
pub fn v_bernstein(classe: i32, nbpoints: i32) -> MathMatrix {
    if classe > THE_MAX_CLASSE {
        panic!("Standard_DimensionError: VBernstein: classe > 26");
    }
    if nbpoints > THE_MAX_NBPOINTS {
        panic!("Standard_DimensionError: VBernstein: nbpoints > 24");
    }
    let (pts, _) = gauss_legendre(-1.0, 1.0, nbpoints as usize);
    let c = classe as usize;
    let n = nbpoints as usize;
    let d = classe as usize - 1;
    let mut m = MathMatrix::with_init(1, c, 1, n, 0.0);
    for i in 1..=c {
        for j in 1..=n {
            let t = 0.5 * (1.0 + pts[j - 1]);
            m.set_value(i, j, bernstein(d, i - 1, t));
        }
    }
    m
}
