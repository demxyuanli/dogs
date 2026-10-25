//! `AppParCurves_Gradient` — the parametrisation refinement driver:
//! a Rogers & Fog projection pass followed by `math_BFGS` on
//! [`ParFunction`](super::function::ParFunction).
//!
//! Source: `ModelingData/TKGeomBase/AppParCurves/AppParCurves_Gradient.gxx:43-217`
//! (instantiated for `AppDef` as `AppDef_MyGradientOfCompute`,
//! `AppDef_MyGradientOfCompute_0.cxx`).
//!
//! Translation notes:
//! - `BSplCLib::PolesCoefficients` / `CoefsD1` are ported in
//!   `occt_core::bspl::coefs`, including the `1 / (k-1)!` scaling of the cache
//!   coefficients (`BSplCLib_CurveComputation.pxx:1584-1589`).
//! - `AppParCurves_Gradient_BFGS` (`AppDef_Gradient_BFGSOfTheGradient_0.cxx:26-54`)
//!   is `math_BFGS(F.NbVariables(), Eps, NbIterations, Eps)` with an overriding
//!   `IsSolutionReached`; [`run_bfgs`] reproduces both through
//!   `BFGS::minimize_with`.
//! - `Parameters` is an in-out argument in OCCT; the port takes `&mut` for the
//!   same reason.

use std::cell::RefCell;

use occt_core::bspl::coefs::{coefs_d1, coefs_d1_2d, poles_coefficients, poles_coefficients2d};
use occt_core::gp::{GpPnt, GpPnt2d, GpVec, GpVec2d};
use occt_math::{BFGS, MathVector};

use crate::approx::line_tool::LineTool;
use crate::app_par_curves::constraint::ConstraintCouple;
use crate::app_par_curves::function::ParFunction;
use crate::app_par_curves::multi_curve::MultiCurve;

/// `RealEpsilon()` (`Standard_Real.hxx`): `DBL_EPSILON`.
const REAL_EPSILON: f64 = f64::EPSILON;

/// `AppParCurves_Gradient` (`gxx:43-217`).
#[derive(Debug, Clone)]
pub struct Gradient {
    /// `ParError` (`gxx:52`): the worst error at each point.
    par_error: MathVector,
    /// `AvError` (`gxx:53`).
    av_error: f64,
    /// `MError3d` (`gxx:54`).
    m_error3d: f64,
    /// `MError2d` (`gxx:55`).
    m_error2d: f64,
    /// `SCU` (`gxx:99`, `gxx:183`, `gxx:189`, `gxx:198`).
    scu: MultiCurve,
    /// `Done`.
    done: bool,
}

impl Gradient {
    /// `AppParCurves_Gradient(SSP, FirstPoint, LastPoint, TheConstraints,
    /// Parameters, Deg, Tol3d, Tol2d, NbIterations)` (`gxx:43-217`).
    ///
    /// `parameters` is `math_Vector&` in OCCT and is left at the refined
    /// parametrisation, hence the `&mut`.
    pub fn new<L: LineTool>(
        ssp: &L,
        first_point: i32,
        last_point: i32,
        the_constraints: &[ConstraintCouple],
        parameters: &mut MathVector,
        deg: i32,
        tol3d: f64,
        tol2d: f64,
        nb_iterations: i32,
    ) -> Self {
        let mut this = Self {
            par_error: MathVector::with_init(
                first_point.max(0) as usize,
                last_point.max(0) as usize,
                0.0,
            ),
            av_error: 0.0,
            m_error3d: 0.0,
            m_error2d: 0.0,
            scu: MultiCurve::new(),
            done: false,
        };
        this.perform(
            ssp,
            first_point,
            last_point,
            the_constraints,
            parameters,
            deg,
            tol3d,
            tol2d,
            nb_iterations,
        );
        this
    }

    /// `Perform` — the constructor body (`gxx:57-217`).
    #[allow(clippy::too_many_arguments)]
    fn perform<L: LineTool>(
        &mut self,
        ssp: &L,
        first_point: i32,
        last_point: i32,
        the_constraints: &[ConstraintCouple],
        parameters: &mut MathVector,
        deg_param: i32,
        tol3d: f64,
        tol2d: f64,
        nb_iterations: i32,
    ) {
        let nb_p3d = ssp.nb_p3d();
        let nb_p2d = ssp.nb_p2d();
        let mut my_nb_p3d = nb_p3d;
        let mut my_nb_p2d = nb_p2d;
        let nb_p = nb_p3d + nb_p2d;
        self.done = false;

        if nb_p3d == 0 {
            my_nb_p3d = 1;
        }
        if nb_p2d == 0 {
            my_nb_p2d = 1;
        }

        // `AppParCurves_ParFunction MyF(SSP, FirstPoint, LastPoint,
        // TheConstraints, Parameters, Deg)` (`gxx:90`).
        let mut my_f = ParFunction::new(
            ssp,
            first_point,
            last_point,
            the_constraints,
            parameters,
            deg_param,
        );

        let mut fval = 0.0;
        if !my_f.value(parameters, &mut fval) {
            self.done = false;
            return;
        }

        self.scu = my_f.curve_value().clone();
        let deg = self.scu.nb_poles() as i32 - 1;
        let npol = (deg + 1).max(1) as usize;

        let mut the_coef: Vec<GpPnt> = vec![GpPnt::new(0.0, 0.0, 0.0); npol * my_nb_p3d as usize];
        let mut the_coef2d: Vec<GpPnt2d> =
            vec![GpPnt2d::new(0.0, 0.0); npol * my_nb_p2d as usize];

        // `gxx:107-126`: store the cache coefficients of every curve.
        let mut i2 = 0usize;
        for k in 1..=nb_p3d {
            let tab_pole = self.scu.curve3d(k);
            let tab_coef = poles_coefficients(&tab_pole);
            for j in 0..npol {
                the_coef[j + i2] = tab_coef[j];
            }
            i2 += npol;
        }
        i2 = 0;
        for k in 1..=nb_p2d {
            let tab_pole = self.scu.curve2d(nb_p3d + k);
            let tab_coef2d = poles_coefficients2d(&tab_pole);
            for j in 0..npol {
                the_coef2d[j + i2] = tab_coef2d[j];
            }
            i2 += npol;
        }

        let mut tab_coef: Vec<GpPnt> = vec![GpPnt::new(0.0, 0.0, 0.0); npol];
        let mut tab_coef2d: Vec<GpPnt2d> = vec![GpPnt2d::new(0.0, 0.0); npol];

        // Iteration de Projection (`gxx:132-175`).
        for j in first_point + 1..=last_point - 1 {
            let mut uf = parameters.value(j as usize);
            let (tab_p, tab_p2d) = ssp.value(j);

            let mut fu = 0.0;
            let mut dfu = 0.0;
            i2 = 0;
            for k in 1..=nb_p3d {
                for l in 0..npol {
                    tab_coef[l] = the_coef[l + i2];
                }
                i2 += npol;
                let (pt, v1) = coefs_d1(uf, &tab_coef);
                let target = tab_p[(k - 1) as usize];
                let my_v = GpVec::new(target.x() - pt.x(), target.y() - pt.y(), target.z() - pt.z());
                fu += my_v.dot(&v1);
                dfu += v1.square_magnitude();
            }
            i2 = 0;
            for k in 1..=nb_p2d {
                for l in 0..npol {
                    tab_coef2d[l] = the_coef2d[l + i2];
                }
                i2 += npol;
                let (pt2d, v12d) = coefs_d1_2d(uf, &tab_coef2d);
                let target = tab_p2d[(k - 1) as usize];
                let my_v2d = GpVec2d::new(target.x() - pt2d.x(), target.y() - pt2d.y());
                fu += my_v2d.dot(&v12d);
                dfu += v12d.square_magnitude();
            }

            if dfu >= REAL_EPSILON {
                let mut du = fu / dfu;
                du = du.signum() * du.abs().min(5.0e-02);
                uf += du;
                parameters.set_value(j as usize, uf);
            }
        }

        if !my_f.value(parameters, &mut fval) {
            self.scu = MultiCurve::new();
            self.done = false;
            return;
        }
        self.m_error3d = my_f.max_error3d();
        self.m_error2d = my_f.max_error2d();

        if self.m_error3d <= tol3d && self.m_error2d <= tol2d {
            self.done = true;
            self.scu = my_f.curve_value().clone();
        } else if nb_iterations != 0 {
            // NbIterations de gradient conjugue (`gxx:191-199`).
            let eps = 1.0e-07;
            let function = RefCell::new(my_f);
            run_bfgs(&function, parameters, tol3d, tol2d, eps, nb_iterations);
            *parameters = function.borrow().new_parameters().clone();
            my_f = function.into_inner();
            self.scu = my_f.curve_value().clone();
        }

        // `gxx:201-217`.
        self.av_error = 0.0;
        for j in first_point..=last_point {
            for k in 1..=nb_p {
                let e = my_f.error(j, k);
                let cur = self.par_error.value(j as usize);
                self.par_error.set_value(j as usize, cur.max(e));
            }
            self.av_error += self.par_error.value(j as usize);
        }
        self.av_error /= (last_point - first_point + 1) as f64;

        self.m_error3d = my_f.max_error3d();
        self.m_error2d = my_f.max_error2d();
        if self.m_error3d <= tol3d && self.m_error2d <= tol2d {
            self.done = true;
        }
    }

    /// `Value()` (`gxx:219-222`).
    pub fn value(&self) -> MultiCurve {
        self.scu.clone()
    }

    /// `IsDone()` (`gxx:224-227`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `Error(Index)` (`gxx:229-232`).
    pub fn error(&self, index: i32) -> f64 {
        self.par_error.value(index as usize)
    }

    /// `AverageError()` (`gxx:234-237`).
    pub fn average_error(&self) -> f64 {
        self.av_error
    }

    /// `MaxError3d()` (`gxx:239-242`).
    pub fn max_error3d(&self) -> f64 {
        self.m_error3d
    }

    /// `MaxError2d()` (`gxx:244-247`).
    pub fn max_error2d(&self) -> f64 {
        self.m_error2d
    }
}

/// `AppParCurves_Gradient_BFGS` (`AppDef_Gradient_BFGSOfTheGradient_0.cxx:26-54`).
///
/// The constructor is `math_BFGS(F.NbVariables(), Eps, NbIterations, Eps)`
/// followed by `Perform(F, StartingPoint)`; `IsSolutionReached` replaces the
/// base test with `2|TheMinimum - PreviousMinimum| <= 1e-10 (|TheMinimum| +
/// |PreviousMinimum|) + 1e-12 || (MaxError3d <= Tol3d && MaxError2d <= Tol2d)`.
/// `function` is shared by the value, gradient and stop-test closures, which is
/// the port's stand-in for OCCT's single mutable `math_MultipleVarFunction`.
fn run_bfgs<L: LineTool>(
    function: &RefCell<ParFunction<L>>,
    parameters: &MathVector,
    tol3d: f64,
    tol2d: f64,
    eps: f64,
    nb_iterations: i32,
) {
    let mut bfgs = BFGS::new();
    // `math_BFGS(NbVariables, Tolerance = Eps, NbIterations, ZEPS = Eps)`
    // (`math_BFGS.cxx:458-476`).
    bfgs.tolerance = eps;
    bfgs.zeps = eps;
    bfgs.max_iter = nb_iterations.max(0) as usize;

    let f = |x: &MathVector| -> f64 {
        let mut val = 0.0;
        let _ = function.borrow_mut().value(x, &mut val);
        val
    };
    let g = |x: &MathVector| -> MathVector {
        let mut grad = MathVector::new(x.lower(), x.upper());
        let _ = function.borrow_mut().gradient(x, &mut grad);
        grad
    };
    let stop = |the_minimum: f64, previous_minimum: f64| -> bool {
        let r1 = 2.0 * (the_minimum - previous_minimum).abs()
            <= 1.0e-10 * (the_minimum.abs() + previous_minimum.abs()) + 1.0e-12;
        let f_ref = function.borrow();
        let r2 = f_ref.max_error3d() <= tol3d && f_ref.max_error2d() <= tol2d;
        r1 || r2
    };

    let _ = bfgs.minimize_with(f, g, parameters, Some(&stop));
}






