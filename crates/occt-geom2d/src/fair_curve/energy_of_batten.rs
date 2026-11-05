//! `FairCurve_EnergyOfBatten` (`FairCurve_EnergyOfBatten.hxx:30-58`,
//! `.cxx:30-134`): the batten objective, the sum of the sagging and tension
//! distributions integrated over `[0, 1]` by `math_GaussSetIntegration`.
//!
//! Rust note: OCCT's `Compute` / `ComputePoles` / `Variable` are virtual
//! overrides called from the base class. The port keeps the base
//! [`Energy`] state inside a `RefCell` (the `math_MultipleVarFunction`
//! interface of this crate takes `&self`) and calls the base's
//! `*_computed` entry points after its own `ComputePoles`, so exactly one
//! `ComputePoles` runs per evaluation, as in `.cxx:60-73`.

use std::cell::{Cell, RefCell};

use occt_math::{
    GaussSetIntegration, MathMatrix, MathMultipleVarFunction,
    MathMultipleVarFunctionWithGradient, MathMultipleVarFunctionWithHessian, MathVector,
};

use super::analysis_code::FairCurveAnalysisCode;
use super::batten_law::BattenLaw;
use super::distribution_of_energy::SharedPoles;
use super::distribution_of_sagging::DistributionOfSagging;
use super::distribution_of_tension::DistributionOfTension;
use super::energy::Energy;

/// `FairCurve_EnergyOfBatten` (`FairCurve_EnergyOfBatten.hxx:30-58`).
pub struct EnergyOfBatten {
    /// The `FairCurve_Energy` base sub-object.
    energy: RefCell<Energy>,
    /// `MyLengthSliding`.
    my_length_sliding: Cell<f64>,
    /// `OriginalSliding`.
    original_sliding: f64,
    /// `MyBattenLaw`.
    my_batten_law: BattenLaw,
    /// `MyTension`.
    my_tension: RefCell<DistributionOfTension>,
    /// `MySagging`.
    my_sagging: RefCell<DistributionOfSagging>,
    /// `MyStatus`.
    my_status: Cell<FairCurveAnalysisCode>,
    /// The shared pole array (also held by the base and the distributions).
    poles: SharedPoles,
}

impl EnergyOfBatten {
    /// `FairCurve_EnergyOfBatten(BSplOrder, FlatKnots, Poles, ContrOrder1,
    /// ContrOrder2, Law, LengthSliding, FreeSliding, Angle1, Angle2)`
    /// (`.cxx:29-51`).
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        bspl_order: i32,
        flat_knots: Vec<f64>,
        poles: SharedPoles,
        contr_order1: i32,
        contr_order2: i32,
        law: &BattenLaw,
        length_sliding: f64,
        free_sliding: bool,
        angle1: f64,
        angle2: f64,
    ) -> Self {
        let nb_val_aux = if free_sliding { 1 } else { 0 };
        Self {
            energy: RefCell::new(Energy::new(
                poles.clone(),
                contr_order1,
                contr_order2,
                free_sliding,
                angle1,
                angle2,
                2,
                0.0,
                0.0,
            )),
            my_length_sliding: Cell::new(length_sliding),
            original_sliding: length_sliding,
            my_batten_law: law.clone(),
            my_tension: RefCell::new(DistributionOfTension::new(
                bspl_order,
                flat_knots.clone(),
                poles.clone(),
                1,
                length_sliding,
                law,
                nb_val_aux,
                false,
            )),
            my_sagging: RefCell::new(DistributionOfSagging::new(
                bspl_order,
                flat_knots,
                poles.clone(),
                1,
                law,
                nb_val_aux,
            )),
            my_status: Cell::new(FairCurveAnalysisCode::Ok),
            poles,
        }
    }

    /// `Variable(X)` (`.cxx:53-63`).
    pub fn variable(&self, x: &mut MathVector) -> bool {
        let ok = self.energy.borrow().variable(x);
        if self.energy.borrow().with_aux_value() {
            let upper = x.upper();
            x.set_value(upper, self.my_length_sliding.get());
        }
        ok
    }

    /// `ComputePoles(X)` (`.cxx:66-73`).
    pub fn compute_poles(&self, x: &MathVector) {
        self.energy.borrow_mut().compute_poles(x);
        if self.energy.borrow().with_aux_value() {
            self.my_length_sliding.set(x.value(x.upper()));
        }
    }

    /// `Status()` (`.hxx:52-55`).
    pub fn status(&self) -> FairCurveAnalysisCode {
        self.my_status.get()
    }

    /// `Compute(DerivativeOrder, Result)` (`.cxx:76-133`).
    fn compute(&self, derivative_order: i32, result: &mut MathVector) -> bool {
        // `math_Vector Debut(1, 1, 0.), Fin(1, 1, 1.);
        //  math_IntegerVector MyOrder(1, 1, 24);` (`.cxx:80-81`).
        let mut debut = MathVector::with_init(1, 1, 0.0);
        let mut fin = MathVector::with_init(1, 1, 1.0);
        let my_order = 24usize;
        let mut ok = false;

        // Blindage contre les longueur de glissement trop exotique (`.cxx:85-97`).
        self.my_status.set(FairCurveAnalysisCode::Ok);
        if self.my_length_sliding.get() > 10.0 * self.original_sliding {
            self.my_status.set(FairCurveAnalysisCode::InfiniteSliding);
            return false;
        }
        if self.my_length_sliding.get() < self.original_sliding / 100.0 {
            self.my_length_sliding.set(self.original_sliding / 100.0);
        }

        // Mise a jour des objets sous-fonction (`.cxx:99-102`).
        let mut tension = self.my_tension.borrow_mut();
        let mut sagging = self.my_sagging.borrow_mut();
        tension.set_derivative_order(derivative_order);
        tension.set_length_sliding(self.my_length_sliding.get());
        sagging.set_derivative_order(derivative_order);
        self.my_batten_law.set_sliding(self.my_length_sliding.get());

        // on decoupe afin d'avoir au moins 2 points d'integration par poles
        // 24 points de Gauss => 12 poles maximum (`.cxx:106-110`).
        let n_poles = self.poles.borrow().len() as i32;
        let nb_interv = (n_poles - 1) / 12 + 1;
        let delta = 1.0 / nb_interv as f64;
        result.init(0.0);

        for ii in 1..=nb_interv {
            debut.set_value(1, ((ii - 1) as f64) * delta);
            fin.set_value(1, (ii as f64) * delta);

            let sum_tension = GaussSetIntegration::new(&mut *tension, &debut, &fin, my_order);
            ok = sum_tension.is_done();
            if !ok {
                return ok;
            }

            let sum_sagging = GaussSetIntegration::new(&mut *sagging, &debut, &fin, my_order);
            ok = sum_sagging.is_done();

            if ok {
                let total = sum_sagging.value().added(sum_tension.value());
                result.add(&total);
            }
        }

        ok
    }
}

impl MathMultipleVarFunction for EnergyOfBatten {
    fn nb_variables(&self) -> usize {
        self.energy.borrow().nb_variables() as usize
    }

    /// `Value(X, E)` (`.cxx:60-68`).
    fn value(&self, x: &MathVector, f: &mut f64) -> bool {
        self.compute_poles(x);
        let mut energy = self.energy.borrow_mut();
        let mut compute = |order: i32, result: &mut MathVector| self.compute(order, result);
        energy.value_computed(f, &mut compute)
    }

    fn as_hessian(&self) -> Option<&dyn MathMultipleVarFunctionWithHessian> {
        Some(self)
    }
}

impl MathMultipleVarFunctionWithGradient for EnergyOfBatten {
    /// `Gradient(X, G)` (`.cxx:71-79`).
    fn gradient(&self, x: &MathVector, g: &mut MathVector) -> bool {
        let mut e = 0.0;
        self.values(x, &mut e, g)
    }

    /// `Values(X, E, G)` (`.cxx:135-149`).
    fn values(&self, x: &MathVector, f: &mut f64, g: &mut MathVector) -> bool {
        self.compute_poles(x);
        let mut energy = self.energy.borrow_mut();
        let mut compute = |order: i32, result: &mut MathVector| self.compute(order, result);
        energy.values_computed(f, g, &mut compute)
    }
}

impl MathMultipleVarFunctionWithHessian for EnergyOfBatten {
    /// `Values(X, E, G, H)` (`.cxx:152-165`).
    fn values_hessian(
        &self,
        x: &MathVector,
        f: &mut f64,
        g: &mut MathVector,
        h: &mut MathMatrix,
    ) -> bool {
        self.compute_poles(x);
        let mut energy = self.energy.borrow_mut();
        let mut compute = |order: i32, result: &mut MathVector| self.compute(order, result);
        energy.values_hessian_computed(f, g, h, &mut compute)
    }
}
