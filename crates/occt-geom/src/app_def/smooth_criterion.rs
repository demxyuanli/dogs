//! `AppDef_SmoothCriterion` — the functional the variational optimizer
//! minimizes.
//!
//! Source: `ModelingData/TKGeomBase/AppDef/AppDef_SmoothCriterion.hxx:34-101`
//! (OCCT 8.0.0).
//!
//! Translation notes:
//! - `occ::handle<FEmTool_Curve>&` in/out parameters are shared
//!   `Rc<RefCell<Curve>>` here, so the criterion and the driver see one curve.
//! - `AssemblyTable()` returns the `NCollection_HArray2` content as
//!   `[dimension][element]` lists of global variable numbers, and
//!   `DependenceTable()` as the `[dimension][dimension]` coupling flags; both
//!   are what [`occt_math::fem_tool::Assembly::new`] consumes.
//! - `EstLength()` returns `double&` in OCCT; the port returns `&mut f64`.

use std::cell::RefCell;
use std::rc::Rc;

use occt_math::fem_tool::Curve;
use occt_math::matrix::MathMatrix;
use occt_math::vector::MathVector;

/// `AppDef_SmoothCriterion` (`hxx:34-101`).
pub trait SmoothCriterion {
    /// `SetParameters(Parameters)`.
    fn set_parameters(&mut self, parameters: &MathVector);

    /// `SetCurve(C)`.
    fn set_curve(&mut self, c: &Rc<RefCell<Curve>>);

    /// `GetCurve(C)`.
    fn get_curve(&self) -> Option<Rc<RefCell<Curve>>>;

    /// `SetEstimation(E1, E2, E3)`.
    fn set_estimation(&mut self, e1: f64, e2: f64, e3: f64);

    /// `EstLength()`.
    fn est_length(&mut self) -> &mut f64;

    /// `GetEstimation(E1, E2, E3)`.
    fn estimation(&self) -> (f64, f64, f64);

    /// `AssemblyTable()`.
    fn assembly_table(&self) -> Vec<Vec<Vec<i32>>>;

    /// `DependenceTable()`.
    fn dependence_table(&self) -> Vec<Vec<i32>>;

    /// `QualityValues(J1min, J2min, J3min, J1, J2, J3)`.
    fn quality_values(
        &mut self,
        j1min: f64,
        j2min: f64,
        j3min: f64,
        j1: &mut f64,
        j2: &mut f64,
        j3: &mut f64,
    ) -> i32;

    /// `ErrorValues(MaxError, QuadraticError, AverageError)`.
    fn error_values(&mut self, max_error: &mut f64, quadratic_error: &mut f64, average_error: &mut f64);

    /// `Hessian(Element, Dimension1, Dimension2, H)`.
    fn hessian(&mut self, element: usize, dim1: i32, dim2: i32, h: &mut MathMatrix);

    /// `Gradient(Element, Dimension, G)`.
    fn gradient(&mut self, element: usize, dim: i32, g: &mut MathVector);

    /// `InputVector(X, AssTable)` — converts the assembly vector back into the
    /// curve.
    fn input_vector(&mut self, x: &MathVector, ass_table: &[Vec<Vec<i32>>]);

    /// `SetWeight(QuadraticWeight, QualityWeight, percentJ1, percentJ2, percentJ3)`.
    fn set_weight(&mut self, quadratic_weight: f64, quality_weight: f64, percent: [f64; 3]);

    /// `GetWeight(QuadraticWeight, QualityWeight)`.
    fn weight(&self) -> (f64, f64);

    /// `SetWeight(Weight)` — the per-point weights.
    fn set_pnt_weight(&mut self, weight: &MathVector);
}
