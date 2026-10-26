//! `FEmTool_ElementsOfRefMatrix` — the vector of upper-triangle products of
//! the Hermite/Jacobi basis at one parameter.
//!
//! Source: `ModelingData/TKGeomBase/FEmTool/FEmTool_ElementsOfRefMatrix.hxx:28-72`
//! and `FEmTool_ElementsOfRefMatrix.cxx:22-83` (OCCT 8.0.0).
//!
//! Translation notes:
//! - `Value(X, F)` is `math_FunctionSet::Value`; the port takes the single
//!   parameter directly and writes `nb_equations()` entries of `f`.
//! - `Basis` / `Aux` are `NCollection_Array1<double>(0, WorkDegree)`; the port
//!   uses `vec![0.0; work_degree + 1]`.
//! - `D1` / `D2` / `D3` write the derivative into the **second** argument
//!   (`cxx:61-66`: `myBase.D1(u, Aux, Basis)`), so `Basis` holds the derivative
//!   of order `myDerOrder`; the port keeps that argument order.

use crate::plib_hermit_jacobi::HermitJacobi;

/// `DerOrder < 0 || DerOrder > 3` (`cxx:26-29`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ConstructionError;

/// `FEmTool_ElementsOfRefMatrix` (`hxx:28-72`).
#[derive(Clone)]
pub struct ElementsOfRefMatrix {
    base: HermitJacobi,
    der_order: i32,
    nb_equations: usize,
}

impl ElementsOfRefMatrix {
    /// `FEmTool_ElementsOfRefMatrix(TheBase, DerOrder)` (`cxx:22-34`).
    pub fn new(base: &HermitJacobi, der_order: i32) -> Result<Self, ConstructionError> {
        if !(0..=3).contains(&der_order) {
            return Err(ConstructionError);
        }
        let work_degree = base.work_degree();
        Ok(Self {
            base: base.clone(),
            der_order,
            nb_equations: ((work_degree + 2) * (work_degree + 1) / 2).max(0) as usize,
        })
    }

    /// `NbVariables()` (`cxx:36-38`): always 1.
    pub fn nb_variables(&self) -> usize {
        1
    }

    /// `NbEquations()` (`cxx:40-43`).
    pub fn nb_equations(&self) -> usize {
        self.nb_equations
    }

    /// `WorkDegree()` of the underlying base.
    pub fn work_degree(&self) -> i32 {
        self.base.work_degree()
    }

    /// `Value(X, F)` (`cxx:45-82`).
    pub fn value(&self, x: f64, f: &mut [f64]) -> bool {
        if f.len() < self.nb_equations {
            return false; // `cxx:47-50` raises `Standard_OutOfRange`.
        }

        let work_degree = self.base.work_degree();
        let n = (work_degree + 1).max(0) as usize;
        let mut basis = vec![0.0f64; n];
        let mut aux = vec![0.0f64; n];
        let mut aux1 = vec![0.0f64; n];
        let mut aux2 = vec![0.0f64; n];

        match self.der_order {
            0 => self.base.d0(x, &mut basis),
            1 => self.base.d1(x, &mut aux, &mut basis),
            // OCCT (`cxx:62-63`) passes the same `Aux` array for every
            // derivative slot; only the last argument (`Basis`) is read
            // afterwards, so distinct scratch buffers are equivalent.
            2 => self.base.d2(x, &mut aux, &mut aux1, &mut basis),
            _ => self.base.d3(x, &mut aux, &mut aux1, &mut aux2, &mut basis),
        }

        let mut ii = 0usize;
        for i in 0..n {
            for j in i..n {
                f[ii] = basis[i] * basis[j];
                ii += 1;
            }
        }
        true
    }
}
