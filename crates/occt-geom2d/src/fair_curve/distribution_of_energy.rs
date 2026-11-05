//! `FairCurve_DistributionOfEnergy` (`FairCurve_DistributionOfEnergy.hxx:30-59`,
//! `.cxx:19-57`): the common `math_FunctionSet` state shared by the local
//! energy distributions integrated by `math_GaussSetIntegration`.
//!
//! Rust note: `MyPoles` is an `NCollection_HArray1<gp_Pnt2d>` handle
//! (`FairCurve_DistributionOfEnergy.hxx:54`) that the owning `FairCurve_Energy`
//! edits through `ComputePoles`; it is modelled as a shared `Rc<RefCell<..>>`
//! holding the 1-based OCCT array (element `i` is `Vec[i - 1]`).

use std::cell::RefCell;
use std::rc::Rc;

use occt_core::gp::pnt2d::GpPnt2d;

/// Shared `NCollection_HArray1<gp_Pnt2d>` handle (1-based, `Vec[i - 1]`).
pub type SharedPoles = Rc<RefCell<Vec<GpPnt2d>>>;

/// `MyPoles->Value(i)` for the 1-based index `i`.
#[inline]
pub(crate) fn pole(poles: &[GpPnt2d], i: i32) -> GpPnt2d {
    poles[(i - 1) as usize]
}

/// `FairCurve_DistributionOfEnergy` (`FairCurve_DistributionOfEnergy.hxx:30-59`).
#[derive(Debug, Clone)]
pub struct DistributionOfEnergy {
    /// `MyBSplOrder`.
    pub(crate) bspl_order: i32,
    /// `MyFlatKnots->Array1()`.
    pub(crate) flat_knots: Vec<f64>,
    /// `MyPoles`.
    pub(crate) poles: SharedPoles,
    /// `MyDerivativeOrder`.
    pub(crate) derivative_order: i32,
    /// `MyNbVar`.
    pub(crate) nb_var: i32,
    /// `MyNbEqua`.
    pub(crate) nb_equa: i32,
    /// `MyNbValAux`.
    pub(crate) nb_val_aux: i32,
}

impl DistributionOfEnergy {
    /// `FairCurve_DistributionOfEnergy(BSplOrder, FlatKnots, Poles,
    /// DerivativeOrder, NbValAux = 0)` (`.cxx:27-42`).
    pub(crate) fn new(
        bspl_order: i32,
        flat_knots: Vec<f64>,
        poles: SharedPoles,
        derivative_order: i32,
        nb_val_aux: i32,
    ) -> Self {
        let mut value = Self {
            bspl_order,
            flat_knots,
            poles,
            derivative_order: 0,
            nb_var: 1,
            nb_equa: 0,
            nb_val_aux,
        };
        // `MyNbVar = 1; SetDerivativeOrder(DerivativeOrder);` (`cxx:40-41`).
        value.set_derivative_order(derivative_order);
        value
    }

    /// `SetDerivativeOrder(DerivativeOrder)` (`.cxx:44-57`).
    pub fn set_derivative_order(&mut self, derivative_order: i32) {
        self.nb_equa = 1;
        if derivative_order >= 1 {
            self.nb_equa += 2 * self.poles.borrow().len() as i32 + self.nb_val_aux;
        }
        if derivative_order >= 2 {
            self.nb_equa = self.nb_equa + (self.nb_equa - 1) * self.nb_equa / 2;
        }
        self.derivative_order = derivative_order;
    }

    /// `NbVariables()` (`.cxx:19-22`).
    pub fn nb_variables(&self) -> usize {
        self.nb_var as usize
    }

    /// `NbEquations()` (`.cxx:24-27`).
    pub fn nb_equations(&self) -> usize {
        self.nb_equa as usize
    }
}
