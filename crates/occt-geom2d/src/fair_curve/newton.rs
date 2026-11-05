//! `FairCurve_Newton` (`FairCurve_Newton.hxx:27-56`, `.cxx:22-50`): the
//! `math_NewtonMinimum` used by the fair-curve objectives with the spatial
//! convergence test of `.cxx:45-50`.

use occt_math::{MathMultipleVarFunctionWithHessian, MathVector, NewtonConvergence, NewtonMinimum};

/// `FairCurve_Newton` (`FairCurve_Newton.hxx:27-56`).
pub struct FairCurveNewton {
    /// The `math_NewtonMinimum` base sub-object.
    nm: NewtonMinimum,
    /// `mySpTol`.
    my_sp_tol: f64,
}

impl FairCurveNewton {
    /// `FairCurve_Newton(theFunction, theSpatialTolerance = 1.0e-7,
    /// theCriteriumTolerance = 1.0e-7, theNbIterations = 40,
    /// theConvexity = 1.0e-6, theWithSingularity = true)`
    /// (`.cxx:22-35`).
    pub fn new(
        nb_variables: usize,
        spatial_tolerance: f64,
        criterium_tolerance: f64,
        nb_iterations: usize,
        convexity: f64,
        with_singularity: bool,
    ) -> Self {
        Self {
            nm: NewtonMinimum::new(
                nb_variables,
                criterium_tolerance,
                nb_iterations,
                convexity,
                with_singularity,
            ),
            my_sp_tol: spatial_tolerance,
        }
    }

    /// `IsConverged()` (`.cxx:45-50`) applied to a Newton state snapshot.
    pub fn is_converged_state(&self, state: &NewtonConvergence) -> bool {
        let n = state.step_norm;
        (n <= 0.01 * self.my_sp_tol)
            || (n <= self.my_sp_tol
                && (state.minimum - state.previous_minimum).abs()
                    <= state.x_tol * state.previous_minimum.abs())
    }

    /// `IsConverged()` (`.cxx:45-50`) on the current state.
    pub fn is_converged(&self) -> bool {
        self.is_converged_state(&self.nm.convergence_state())
    }

    /// `Perform(theFunction, theStartingPoint)`: `math_NewtonMinimum::Perform`
    /// (`math_NewtonMinimum.cxx:77-263`) driven by the redefined
    /// `IsConverged()`.
    pub fn perform(
        &mut self,
        f: &dyn MathMultipleVarFunctionWithHessian,
        starting_point: &MathVector,
    ) {
        let sp_tol = self.my_sp_tol;
        self.nm.perform_with(f, starting_point, |state| {
            let n = state.step_norm;
            (n <= 0.01 * sp_tol)
                || (n <= sp_tol
                    && (state.minimum - state.previous_minimum).abs()
                        <= state.x_tol * state.previous_minimum.abs())
        });
    }

    /// `IsDone()` (`math_NewtonMinimum.lxx:25-28`).
    pub fn is_done(&self) -> bool {
        self.nm.is_done()
    }

    /// `Location()` (`math_NewtonMinimum.lxx:30-34`).
    pub fn location(&self) -> &MathVector {
        self.nm.location()
    }
}
