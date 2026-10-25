//! `AppParCurves_BSpGradient` — the B-spline form of the parametrisation
//! refinement driver: a Rogers & Fog projection pass followed by `math_BFGS`
//! on [`BspFunction`](super::bsp_function::BspFunction).
//!
//! Source: `ModelingData/TKGeomBase/AppParCurves/AppParCurves_BSpGradient.gxx:82-410`,
//! declaration `AppDef_MyBSplGradientOfBSplineCompute.hxx:40-134`; instantiated
//! for `AppDef` as `AppDef_MyBSplGradientOfBSplineCompute`
//! (`AppDef_MyBSplGradientOfBSplineCompute_0.cxx:24-39`).
//!
//! Translation notes:
//! - `AppParCurves_BSpGradient_BFGS`
//!   (`AppDef_BSpGradient_BFGSOfMyBSplGradientOfBSplineCompute_0.cxx:26-56`) is
//!   `math_BFGS(F.NbVariables(), Eps, NbIterations, Eps)` with an overriding
//!   `IsSolutionReached`; [`run_bfgs`] reproduces both through
//!   `BFGS::minimize_with`.
//! - `Parameters` is an in-out argument in OCCT (`gxx:88`, `gxx:363`), so the
//!   port takes `&mut` for the same reason.
//! - The two `FirstConstraint` / `LastConstraint` statics (`gxx:36-77`) are
//!   file-local in the `.gxx`, so they stay file-local here.
//! - `TabV` / `TabV2d` (`gxx:179-180`) are declared and never read; the port
//!   drops them.

use std::cell::RefCell;

use occt_core::gp::{GpPnt, GpPnt2d, GpVec, GpVec2d};
use occt_math::{BFGS, MathVector};

use crate::approx::line_tool::LineTool;
use crate::app_par_curves::bsp_function::BspFunction;
use crate::app_par_curves::constraint::{Constraint, ConstraintCouple};
use crate::app_par_curves::least_square::LeastSquare;
use crate::app_par_curves::multi_bspline_curve::MultiBSpCurve;

/// `RealEpsilon()` (`Standard_Real.hxx`): `DBL_EPSILON`.
const REAL_EPSILON: f64 = f64::EPSILON;

/// `Cons >= AppParCurves_TangencyPoint`, using the enum's declaration order
/// (`AppParCurves_Constraint.hxx:22-30`).
fn at_least_tangency(c: Constraint) -> bool {
    matches!(c, Constraint::TangencyPoint | Constraint::CurvaturePoint)
}

/// `FirstConstraint(TheConstraints, FirstPoint)` (`gxx:36-55`).
fn first_constraint(the_constraints: &[ConstraintCouple], first_point: i32) -> Constraint {
    let mut cons = Constraint::NoConstraint;
    for couple in the_constraints {
        cons = couple.constraint();
        if couple.index() == first_point {
            break;
        }
    }
    cons
}

/// `LastConstraint(TheConstraints, LastPoint)` (`gxx:57-77`).
fn last_constraint(the_constraints: &[ConstraintCouple], last_point: i32) -> Constraint {
    let mut cons = Constraint::NoConstraint;
    for couple in the_constraints {
        cons = couple.constraint();
        if couple.index() == last_point {
            break;
        }
    }
    cons
}

/// `AppParCurves_BSpGradient` (`hxx:124-134`).
#[derive(Debug, Clone)]
pub struct BspGradient {
    /// `SCU` (`hxx:125`).
    scu: MultiBSpCurve,
    /// `ParError` (`hxx:126`).
    par_error: MathVector,
    /// `AvError` (`hxx:127`).
    av_error: f64,
    /// `MError3d` (`hxx:128`).
    m_error3d: f64,
    /// `MError2d` (`hxx:129`).
    m_error2d: f64,
    /// `mylambda1` (`hxx:130`).
    my_lambda1: f64,
    /// `mylambda2` (`hxx:131`).
    my_lambda2: f64,
    /// `myIsLambdaDefined` (`hxx:132`).
    my_is_lambda_defined: bool,
    /// `Done` (`hxx:133`).
    done: bool,
}

impl BspGradient {
    /// `AppParCurves_BSpGradient(SSP, FirstPoint, LastPoint, TheConstraints,
    /// Parameters, Knots, Mults, Deg, Tol3d, Tol2d, NbIterations)` (`gxx:82-108`).
    #[allow(clippy::too_many_arguments)]
    pub fn new<L: LineTool>(
        ssp: &L,
        first_point: i32,
        last_point: i32,
        the_constraints: &[ConstraintCouple],
        parameters: &mut MathVector,
        knots: &[f64],
        mults: &[i32],
        deg: i32,
        tol3d: f64,
        tol2d: f64,
        nb_iterations: i32,
    ) -> Self {
        let mut this = Self {
            scu: MultiBSpCurve::new(),
            par_error: MathVector::with_init(
                first_point.max(0) as usize,
                last_point.max(0) as usize,
                0.0,
            ),
            av_error: 0.0,
            m_error3d: 0.0,
            m_error2d: 0.0,
            my_lambda1: 0.0,
            my_lambda2: 0.0,
            my_is_lambda_defined: false,
            done: false,
        };
        this.perform(
            ssp,
            first_point,
            last_point,
            the_constraints,
            parameters,
            knots,
            mults,
            deg,
            tol3d,
            tol2d,
            nb_iterations,
        );
        this
    }

    /// `AppParCurves_BSpGradient(..., NbIterations, lambda1, lambda2)`
    /// (`gxx:110-140`): the tangency multipliers come from the caller instead
    /// of being searched by the fitting least squares.
    #[allow(clippy::too_many_arguments)]
    pub fn with_lambdas<L: LineTool>(
        ssp: &L,
        first_point: i32,
        last_point: i32,
        the_constraints: &[ConstraintCouple],
        parameters: &mut MathVector,
        knots: &[f64],
        mults: &[i32],
        deg: i32,
        tol3d: f64,
        tol2d: f64,
        nb_iterations: i32,
        lambda1: f64,
        lambda2: f64,
    ) -> Self {
        let mut this = Self {
            scu: MultiBSpCurve::new(),
            par_error: MathVector::with_init(
                first_point.max(0) as usize,
                last_point.max(0) as usize,
                0.0,
            ),
            av_error: 0.0,
            m_error3d: 0.0,
            m_error2d: 0.0,
            my_lambda1: lambda1,
            my_lambda2: lambda2,
            my_is_lambda_defined: true,
            done: false,
        };
        this.perform(
            ssp,
            first_point,
            last_point,
            the_constraints,
            parameters,
            knots,
            mults,
            deg,
            tol3d,
            tol2d,
            nb_iterations,
        );
        this
    }

    /// `Perform` (`gxx:142-381`).
    #[allow(clippy::too_many_arguments)]
    fn perform<L: LineTool>(
        &mut self,
        ssp: &L,
        first_point: i32,
        last_point: i32,
        the_constraints: &[ConstraintCouple],
        parameters: &mut MathVector,
        knots: &[f64],
        mults: &[i32],
        deg: i32,
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

        // `gxx:187-191`: `nbpoles = -Deg - 1; nbpoles += sum(Mults)`.
        let mut nbpoles = -deg - 1;
        for m in mults {
            nbpoles += *m;
        }
        let npol = nbpoles.max(0) as usize;

        let first_cons = first_constraint(the_constraints, first_point);
        let last_cons = last_constraint(the_constraints, last_point);

        // `AppParCurves_BSpParFunction MyF(SSP, FirstPoint, LastPoint,
        // TheConstraints, Parameters, Knots, Mults, nbpoles)` (`gxx:201-202`).
        let mut my_f = BspFunction::new(
            ssp,
            first_point,
            last_point,
            the_constraints,
            parameters,
            knots,
            mults,
            nbpoles,
        );

        if at_least_tangency(first_cons) || at_least_tangency(last_cons) {
            if !self.my_is_lambda_defined {
                // `gxx:208-222`.
                let thefitt = LeastSquare::with_knots_parameters(
                    ssp,
                    knots,
                    mults,
                    first_point,
                    last_point,
                    first_cons,
                    last_cons,
                    parameters,
                    nbpoles,
                );
                if at_least_tangency(first_cons) {
                    self.my_lambda1 = thefitt.first_lambda();
                    my_f.set_first_lambda(self.my_lambda1);
                }
                if at_least_tangency(last_cons) {
                    self.my_lambda2 = thefitt.last_lambda();
                    my_f.set_last_lambda(self.my_lambda2);
                }
            } else {
                // `gxx:224-228`.
                my_f.set_first_lambda(self.my_lambda1);
                my_f.set_last_lambda(self.my_lambda2);
            }
        }

        // `gxx:230-233`.
        let mut fval = 0.0;
        my_f.value(parameters, &mut fval);
        self.m_error3d = my_f.max_error3d();
        self.m_error2d = my_f.max_error2d();
        self.scu = my_f.curve_value().clone();

        if self.m_error3d > tol3d || self.m_error2d > tol2d {
            // Storage of curve poles for projection (`gxx:237-255`).
            let mut the_poles: Vec<GpPnt> =
                vec![GpPnt::new(0.0, 0.0, 0.0); npol * my_nb_p3d as usize];
            let mut the_poles2d: Vec<GpPnt2d> =
                vec![GpPnt2d::new(0.0, 0.0); npol * my_nb_p2d as usize];

            let mut i2 = 0usize;
            for k in 1..=nb_p3d {
                let tab_pole = self.scu.base().curve3d(k);
                for j in 0..npol {
                    the_poles[j + i2] = tab_pole[j];
                }
                i2 += npol;
            }
            i2 = 0;
            for k in 1..=nb_p2d {
                let tab_pole = self.scu.base().curve2d(nb_p3d + k);
                for j in 0..npol {
                    the_poles2d[j + i2] = tab_pole[j];
                }
                i2 += npol;
            }

            // `gxx:260-330`: one projection step per interior point.
            for j in first_point + 1..=last_point - 1 {
                let mut uf = parameters.value(j as usize);
                let (tab_p, tab_p2d) = ssp.value(j);

                let mut fu = 0.0;
                let mut dfu = 0.0;
                let mut i2 = 0usize;

                let index_deb = my_f.k_index().value(j as usize) + 1;
                let index_fin = index_deb + deg;

                for k in 1..=nb_p3d {
                    let mut a = 0.0;
                    let mut b = 0.0;
                    let mut c = 0.0;
                    let mut d = 0.0;
                    let mut e = 0.0;
                    let mut f = 0.0;
                    for l in index_deb..=index_fin {
                        let pt = the_poles[(l + i2 as i32) as usize - 1];
                        let px = pt.x();
                        let py = pt.y();
                        let pz = pt.z();
                        let aa = my_f.function_matrix().value(j as usize, l as usize);
                        let da = my_f.derivative_function_matrix().value(j as usize, l as usize);
                        a += aa * px;
                        d += da * px;
                        b += aa * py;
                        e += da * py;
                        c += aa * pz;
                        f += da * pz;
                    }
                    let pt = GpPnt::new(a, b, c);
                    let v1 = GpVec::new(d, e, f);
                    i2 += npol;

                    let target = tab_p[(k - 1) as usize];
                    let my_v = GpVec::new(target.x() - pt.x(), target.y() - pt.y(), target.z() - pt.z());
                    fu += my_v.dot(&v1);
                    dfu += v1.square_magnitude();
                }
                i2 = 0;
                for k in 1..=nb_p2d {
                    let mut a = 0.0;
                    let mut b = 0.0;
                    let mut d = 0.0;
                    let mut e = 0.0;
                    for l in index_deb..=index_fin {
                        let pt2d = the_poles2d[(l + i2 as i32) as usize - 1];
                        let px = pt2d.x();
                        let py = pt2d.y();
                        let aa = my_f.function_matrix().value(j as usize, l as usize);
                        let da = my_f.derivative_function_matrix().value(j as usize, l as usize);
                        a += aa * px;
                        d += da * px;
                        b += aa * py;
                        e += da * py;
                    }
                    let pt2d = GpPnt2d::new(a, b);
                    let v12d = GpVec2d::new(d, e);
                    i2 += npol;

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

            // `gxx:343-347`.
            my_f.value(parameters, &mut fval);
            self.m_error3d = my_f.max_error3d();
            self.m_error2d = my_f.max_error2d();
        }

        // `gxx:349-358`.
        if self.m_error3d <= tol3d && self.m_error2d <= tol2d {
            self.done = true;
        } else if nb_iterations != 0 {
            let eps = 1.0e-07;
            let function = RefCell::new(my_f);
            run_bfgs(&function, parameters, tol3d, tol2d, eps, nb_iterations);
            my_f = function.into_inner();
        }

        // `gxx:361`.
        self.scu = my_f.curve_value().clone();

        // `gxx:363-375`.
        self.av_error = 0.0;
        for j in first_point..=last_point {
            parameters.set_value(j as usize, my_f.new_parameters().value(j as usize));
            for k in 1..=nb_p {
                let e = my_f.error(j, k);
                let cur = self.par_error.value(j as usize);
                self.par_error.set_value(j as usize, cur.max(e));
            }
            self.av_error += self.par_error.value(j as usize);
        }
        self.av_error /= (last_point - first_point + 1) as f64;

        // `gxx:377-381`.
        self.m_error3d = my_f.max_error3d();
        self.m_error2d = my_f.max_error2d();
        if self.m_error3d <= tol3d && self.m_error2d <= tol2d {
            self.done = true;
        }
    }

    /// `Value()` (`gxx:383-386`).
    pub fn value(&self) -> MultiBSpCurve {
        self.scu.clone()
    }

    /// `IsDone()` (`gxx:388-391`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `Error(Index)` (`gxx:393-396`).
    pub fn error(&self, index: i32) -> f64 {
        self.par_error.value(index as usize)
    }

    /// `AverageError()` (`gxx:398-401`).
    pub fn average_error(&self) -> f64 {
        self.av_error
    }

    /// `MaxError3d()` (`gxx:403-406`).
    pub fn max_error3d(&self) -> f64 {
        self.m_error3d
    }

    /// `MaxError2d()` (`gxx:408-411`).
    pub fn max_error2d(&self) -> f64 {
        self.m_error2d
    }
}

/// `AppParCurves_BSpGradient_BFGS`
/// (`AppDef_BSpGradient_BFGSOfMyBSplGradientOfBSplineCompute_0.cxx:26-56`).
///
/// The constructor is `math_BFGS(F.NbVariables(), Eps, NbIterations, Eps)`
/// followed by `Perform(F, StartingPoint)`; `IsSolutionReached` replaces the
/// base test with `2|TheMinimum - PreviousMinimum| <= 1e-10 (|TheMinimum| +
/// |PreviousMinimum|) + 1e-12 || (MaxError3d <= Tol3d && MaxError2d <= Tol2d)`.
/// `function` is shared by the value, gradient and stop-test closures, which is
/// the port's stand-in for OCCT's single mutable `math_MultipleVarFunction`.
fn run_bfgs<L: LineTool>(
    function: &RefCell<BspFunction<L>>,
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






