//! `PLib_HermitJacobi` — the Hermite + W(t)*Jacobi basis used by `FEmTool`
//! and `AppDef_LinearCriteria`.
//!
//! Source: `FoundationClasses/TKMath/PLib/PLib_HermitJacobi.hxx:25-141`,
//! `PLib_HermitJacobi.cxx:24-372` (OCCT 8.0.0).
//!
//! P(t) = H(t) + W(t) * Q(t), `W(t) = (1 - t^2)^(2*niv_constr + 2)`; the
//! Hermite part carries `2*niv_constr + 2` coefficients, the Jacobi part the
//! remaining `work_degree - (2*niv_constr + 1)` ones.
//!
//! Translation notes:
//! - The `NCollection_Array1<double>&` outputs of `ToCoefficients` /
//!   `D0..D3` become return values / slices whose index 0 is the OCCT lower
//!   index, so `ibegHJC` / `ibegC` drop out.
//! - `MaxError` / `ReduceDegree` / `AverageError` take `double& HermJacCoeff`
//!   in OCCT; the port reads the array as a slice (Jacobson's
//!   `ReduceDegree` never writes the coefficients, it only inspects them).
//! - `GetHermiteMatrix_C{0,1,2}` (`cxx:53-84`) are `PLib::HermiteCoefficients
//!   (-1, 1, niv, niv)` computed once; [`hermite_matrix`] caches them the same
//!   way.
//! - `NCollection_LocalArray<double> jac0(4 * 20)` (`cxx:203-207`) only bounds
//!   the scratch size (degree <= 30); the port allocates exactly
//!   `jac_degree + 1` values.

use std::sync::OnceLock;

use occt_core::bspl::coefs::{eval_polynomial, no_derivative_eval_polynomial};
use occt_core::bspl::plib_jacobi::JacobiPolynomial;
use occt_core::kernel::geomabs::Shape;

use crate::plib_hermite::hermite_coefficients;

/// `WCoeff_C0` (`cxx:27`): `W(t) = 1 - t^2`.
const W_COEFF_C0: [f64; 3] = [1.0, 0.0, -1.0];
/// `WCoeff_C1` (`cxx:31`): `W(t) = 1 - 2t^2 + t^4`.
const W_COEFF_C1: [f64; 5] = [1.0, 0.0, -2.0, 0.0, 1.0];
/// `WCoeff_C2` (`cxx:35`): `W(t) = 1 - 3t^2 + 3t^4 - t^6`.
const W_COEFF_C2: [f64; 7] = [1.0, 0.0, -3.0, 0.0, 3.0, 0.0, -1.0];

/// `GetWCoefficients(niv_constr)` (`cxx:39-51`).
fn w_coefficients(niv_constr: i32) -> &'static [f64] {
    match niv_constr {
        1 => &W_COEFF_C1,
        2 => &W_COEFF_C2,
        _ => &W_COEFF_C0,
    }
}

/// `GetHermiteMatrix(niv_constr)` (`cxx:53-96`): the `(degree_h + 1)` rows of
/// Hermite polynomial coefficients on `[-1, 1]`.
fn hermite_matrix(niv_constr: i32) -> &'static Vec<Vec<f64>> {
    static C0: OnceLock<Vec<Vec<f64>>> = OnceLock::new();
    static C1: OnceLock<Vec<Vec<f64>>> = OnceLock::new();
    static C2: OnceLock<Vec<Vec<f64>>> = OnceLock::new();
    let cell = match niv_constr {
        1 => &C1,
        2 => &C2,
        _ => &C0,
    };
    cell.get_or_init(|| {
        // `PLib::HermiteCoefficients(-1., 1., niv, niv, aResult)` cannot fail
        // for these arguments; the fallback keeps the signature panic-free.
        hermite_coefficients(-1.0, 1.0, niv_constr, niv_constr).unwrap_or_default()
    })
}

/// `class PLib_HermitJacobi` (`hxx:58-141`).
#[derive(Clone)]
pub struct HermitJacobi {
    /// `myJacobi`
    jacobi: JacobiPolynomial,
}

impl HermitJacobi {
    /// `PLib_HermitJacobi(WorkDegree, ConstraintOrder)` (`cxx:102-104`).
    pub fn new(work_degree: i32, constraint_order: Shape) -> Result<Self, ()> {
        Ok(Self {
            jacobi: JacobiPolynomial::new(work_degree, constraint_order)?,
        })
    }

    /// `WorkDegree()` (`hxx:124`).
    pub fn work_degree(&self) -> i32 {
        self.jacobi.work_degree()
    }

    /// `NivConstr()` (`hxx:127`).
    pub fn niv_constr(&self) -> i32 {
        self.jacobi.niv_constr()
    }

    /// `MaxError(Dimension, HermJacCoeff, NewDegree)` (`cxx:108-112`).
    pub fn max_error(&self, dimension: i32, hermit_jac_coeff: &[f64], new_degree: i32) -> f64 {
        self.jacobi.max_error(dimension, hermit_jac_coeff, new_degree)
    }

    /// `ReduceDegree(Dimension, MaxDegree, Tol, HermJacCoeff, NewDegree,
    /// MaxError)` (`cxx:116-124`): returns `(NewDegree, MaxError)`.
    pub fn reduce_degree(
        &self,
        dimension: i32,
        max_degree: i32,
        tol: f64,
        hermit_jac_coeff: &[f64],
    ) -> (i32, f64) {
        self.jacobi
            .reduce_degree(dimension, max_degree, tol, hermit_jac_coeff)
    }

    /// `AverageError(Dimension, HermJacCoeff, NewDegree)` (`cxx:128-132`).
    pub fn average_error(&self, dimension: i32, hermit_jac_coeff: &[f64], new_degree: i32) -> f64 {
        self.jacobi
            .average_error(dimension, hermit_jac_coeff, new_degree)
    }

    /// `ToCoefficients(Dimension, Degree, HermJacCoeff, Coefficients)`
    /// (`cxx:138-188`).
    ///
    /// `hermit_jac_coeff` is laid out degree-major / dimension-minor, `(Degree +
    /// 1) * Dimension` values; the returned vector has the same layout and is
    /// the polynomial in the canonical basis.
    pub fn to_coefficients(
        &self,
        dimension: i32,
        degree: i32,
        hermit_jac_coeff: &[f64],
    ) -> Vec<f64> {
        let dim = dimension.max(0) as usize;
        let niv = self.niv_constr();
        let degree_h = 2 * niv + 1;
        let matrix = hermite_matrix(niv);

        let total = ((degree + 1).max(0) as usize) * dim;
        let mut aux_coeff = vec![0.0f64; total];

        for k in 0..=degree_h {
            let kdim = (k as usize) * dim;
            for i in 0..=niv {
                let h1 = matrix
                    .get(i as usize)
                    .and_then(|row| row.get(k as usize))
                    .copied()
                    .unwrap_or(0.0);
                let h2 = matrix
                    .get((i + niv + 1) as usize)
                    .and_then(|row| row.get(k as usize))
                    .copied()
                    .unwrap_or(0.0);
                let i1 = (i as usize) * dim;
                let i2 = ((i + niv + 1) as usize) * dim;
                for d in 0..dim {
                    let a = hermit_jac_coeff.get(i1 + d).copied().unwrap_or(0.0);
                    let b = hermit_jac_coeff.get(i2 + d).copied().unwrap_or(0.0);
                    if kdim + d < aux_coeff.len() {
                        aux_coeff[kdim + d] += a * h1 + b * h2;
                    }
                }
            }
        }

        let kdim = total;
        for k in (((degree_h + 1) as usize) * dim)..kdim {
            aux_coeff[k] = hermit_jac_coeff.get(k).copied().unwrap_or(0.0);
        }

        if degree > degree_h {
            self.jacobi.to_coefficients(dimension, degree, &aux_coeff)
        } else {
            aux_coeff
        }
    }

    /// `D0123(NDeriv, U, BasisValue, BasisD1, BasisD2, BasisD3)`
    /// (`cxx:196-334`, FORTRAN `MPOBAS`).
    ///
    /// Each output slice must hold `work_degree + 1` values; only the first
    /// `NDeriv + 1` derivative levels are written, as in OCCT.
    pub fn d0123(
        &self,
        n_deriv: i32,
        u: f64,
        basis_value: &mut [f64],
        basis_d1: &mut [f64],
        basis_d2: &mut [f64],
        basis_d3: &mut [f64],
    ) {
        let niv = self.niv_constr();
        let work_degree = self.work_degree();
        let degree_h = 2 * niv + 1;
        let jac_degree = work_degree - degree_h - 1;
        let matrix = hermite_matrix(niv);

        let mut jac0 = vec![0.0f64; (jac_degree.max(0) + 1) as usize];
        let mut jac1 = vec![0.0f64; (jac_degree.max(0) + 1) as usize];
        let mut jac2 = vec![0.0f64; (jac_degree.max(0) + 1) as usize];
        let mut jac3 = vec![0.0f64; (jac_degree.max(0) + 1) as usize];
        let mut w_values = vec![0.0f64; (n_deriv + 1).max(0) as usize];

        // Evaluation des polynomes d'hermite (`cxx:222-245`).
        let n_rows = (degree_h + 1) as usize;
        let n_cols = (n_deriv + 1).max(0) as usize;
        let mut hermit_values = vec![vec![0.0f64; n_cols]; n_rows];
        for i in 0..n_rows {
            let coeffs = matrix.get(i).cloned().unwrap_or_default();
            if n_deriv == 0 {
                let r = no_derivative_eval_polynomial(u, degree_h, 1, &coeffs);
                hermit_values[i][0] = r[0];
            } else {
                let r = eval_polynomial(u, n_deriv, degree_h, 1, &coeffs);
                for k in 0..=n_deriv as usize {
                    hermit_values[i][k] = r.get(k).copied().unwrap_or(0.0);
                }
            }
        }

        // Evaluation des polynomes de Jacobi + de W(t) (`cxx:247-284`).
        if jac_degree >= 0 {
            match n_deriv {
                0 => self.jacobi.d0(u, &mut jac0),
                1 => self.jacobi.d1(u, &mut jac0, &mut jac1),
                2 => self.jacobi.d2(u, &mut jac0, &mut jac1, &mut jac2),
                3 => self.jacobi.d3(u, &mut jac0, &mut jac1, &mut jac2, &mut jac3),
                _ => {}
            }

            let w_coeff = w_coefficients(niv);
            if n_deriv == 0 {
                let r = no_derivative_eval_polynomial(u, degree_h + 1, 1, w_coeff);
                w_values[0] = r[0];
            } else {
                let r = eval_polynomial(u, n_deriv, degree_h + 1, 1, w_coeff);
                for k in 0..=n_deriv as usize {
                    w_values[k] = r.get(k).copied().unwrap_or(0.0);
                }
            }
        }

        // Evaluation a l'ordre 0 (`cxx:286-294`).
        for i in 0..n_rows {
            if i < basis_value.len() {
                basis_value[i] = hermit_values[i][0];
            }
        }
        let w0 = w_values[0];
        let mut j = 0usize;
        let mut i = (degree_h + 1).max(0);
        while i <= work_degree {
            if let Some(slot) = basis_value.get_mut(i as usize) {
                *slot = w0 * jac0.get(j).copied().unwrap_or(0.0);
            }
            j += 1;
            i += 1;
        }

        // Evaluation a l'ordre 1 (`cxx:296-334`).
        if n_deriv >= 1 {
            let w1 = w_values[1];
            for i in 0..n_rows {
                if i < basis_d1.len() {
                    basis_d1[i] = hermit_values[i][1];
                }
            }
            let mut j = 0usize;
            let mut i = (degree_h + 1).max(0);
            while i <= work_degree {
                if let Some(slot) = basis_d1.get_mut(i as usize) {
                    *slot = w0 * jac1.get(j).copied().unwrap_or(0.0)
                        + w1 * jac0.get(j).copied().unwrap_or(0.0);
                }
                j += 1;
                i += 1;
            }

            if n_deriv >= 2 {
                let w2 = w_values[2];
                for i in 0..n_rows {
                    if i < basis_d2.len() {
                        basis_d2[i] = hermit_values[i][2];
                    }
                }
                let mut j = 0usize;
                let mut i = (degree_h + 1).max(0);
                while i <= work_degree {
                    if let Some(slot) = basis_d2.get_mut(i as usize) {
                        *slot = w0 * jac2.get(j).copied().unwrap_or(0.0)
                            + 2.0 * w1 * jac1.get(j).copied().unwrap_or(0.0)
                            + w2 * jac0.get(j).copied().unwrap_or(0.0);
                    }
                    j += 1;
                    i += 1;
                }

                if n_deriv == 3 {
                    let w3 = w_values[3];
                    for i in 0..n_rows {
                        if i < basis_d3.len() {
                            basis_d3[i] = hermit_values[i][3];
                        }
                    }
                    let mut j = 0usize;
                    let mut i = (degree_h + 1).max(0);
                    while i <= work_degree {
                        if let Some(slot) = basis_d3.get_mut(i as usize) {
                            let a = jac3.get(j).copied().unwrap_or(0.0);
                            let b = jac0.get(j).copied().unwrap_or(0.0);
                            let c = jac2.get(j).copied().unwrap_or(0.0);
                            let d = jac1.get(j).copied().unwrap_or(0.0);
                            *slot = w0 * a + w3 * b + 3.0 * (w1 * c + w2 * d);
                        }
                        j += 1;
                        i += 1;
                    }
                }
            }
        }
    }

    /// `D0(U, BasisValue)` (`cxx:338-341`).
    pub fn d0(&self, u: f64, basis_value: &mut [f64]) {
        let mut empty1: [f64; 0] = [];
        let mut empty2: [f64; 0] = [];
        let mut empty3: [f64; 0] = [];
        self.d0123(0, u, basis_value, &mut empty1, &mut empty2, &mut empty3);
    }

    /// `D1(U, BasisValue, BasisD1)` (`cxx:345-351`).
    pub fn d1(&self, u: f64, basis_value: &mut [f64], basis_d1: &mut [f64]) {
        let mut empty2: [f64; 0] = [];
        let mut empty3: [f64; 0] = [];
        self.d0123(1, u, basis_value, basis_d1, &mut empty2, &mut empty3);
    }

    /// `D2(U, BasisValue, BasisD1, BasisD2)` (`cxx:355-360`).
    pub fn d2(
        &self,
        u: f64,
        basis_value: &mut [f64],
        basis_d1: &mut [f64],
        basis_d2: &mut [f64],
    ) {
        let mut empty3: [f64; 0] = [];
        self.d0123(2, u, basis_value, basis_d1, basis_d2, &mut empty3);
    }

    /// `D3(U, BasisValue, BasisD1, BasisD2, BasisD3)` (`cxx:364-371`).
    pub fn d3(
        &self,
        u: f64,
        basis_value: &mut [f64],
        basis_d1: &mut [f64],
        basis_d2: &mut [f64],
        basis_d3: &mut [f64],
    ) {
        self.d0123(3, u, basis_value, basis_d1, basis_d2, basis_d3);
    }
}
