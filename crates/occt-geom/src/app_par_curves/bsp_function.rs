//! `AppParCurves_BSpFunction` — the parametrisation objective and gradient for
//! the B-spline least-squares driver.
//!
//! Source: `ModelingData/TKGeomBase/AppParCurves/AppParCurves_BSpFunction.gxx:27-338`
//! and `AppDef_BSplineFunctionOfMyGradientOfCompute.hxx` (OCCT 8.0.0).
//!
//! Translation notes:
//! - `Value`'s and `Perform`'s constrained arms are **empty in OCCT**
//!   (`gxx:205-209`, `gxx:238-240`): the B-spline driver never enables the
//!   Uzawa resolution, it only re-fits with the `lambda1` / `lambda2` penalty
//!   pair (`gxx:193`, `gxx:221`). The port keeps the `Contraintes` scan and the
//!   empty arms instead of dropping them.
//! - Because those arms are empty, `MyF` is never filled, so `Error` returns
//!   `MyLeastSquare.Distance()` in the unconstrained case (`gxx:301-304`) and
//!   `sqrt(MyF(IPoint, CurveIndex))` — all zeroes — otherwise.
//! - `tabdim` (`gxx:105-137`) is only read inside those dead arms; the port
//!   keeps it for the record.

use occt_math::{MathIntVector, MathMatrix, MathVector};

use crate::app_def::MultiLine;
use crate::approx::line_tool::LineTool;
use crate::app_par_curves::constraint::{Constraint, ConstraintCouple};
use crate::app_par_curves::least_square::LeastSquare;
use crate::app_par_curves::multi_bspline_curve::MultiBSpCurve;

use super::function::at_least_pass;

/// `AppParCurves_BSpFunction` (`gxx:27-338`).
///
/// `A`, `DA`, `NbCu` and `myConstraints` are initialised by the constructor
/// and never read again in the `.gxx` either: they belong to the empty
/// constrained arms. Kept as the OCCT record.
#[allow(dead_code)]
pub struct BspFunction<L: LineTool = MultiLine> {
    /// `Done` (`gxx:196`, `gxx:224`): written on the failure paths only.
    #[allow(dead_code)]
    done: bool,
    /// `MyMultiLine`.
    my_multi_line: L,
    /// `MyMultiBSpCurve`.
    my_multi_bsp_curve: MultiBSpCurve,
    /// `myParameters`.
    my_parameters: MathVector,
    /// `FVal`.
    f_val: f64,
    /// `ValGrad_F`.
    val_grad_f: MathVector,
    /// `MyF` (`gxx:40`): never filled, see the module notes.
    my_f: MathMatrix,
    /// `PTLX`.
    ptlx: MathMatrix,
    /// `PTLY`.
    ptly: MathMatrix,
    /// `PTLZ`.
    ptlz: MathMatrix,
    /// `A`.
    a: MathMatrix,
    /// `DA`.
    da: MathMatrix,
    /// `MyLeastSquare`.
    my_least_square: LeastSquare,
    /// `Contraintes`.
    contraintes: bool,
    /// `NbP`.
    nbp: i32,
    /// `NbCu`.
    nb_cu: i32,
    /// `Adeb`.
    adeb: i32,
    /// `Afin`.
    afin: i32,
    /// `nbpoles`.
    nbpoles: i32,
    /// `tabdim` (`gxx:105`).
    tabdim: Vec<i32>,
    /// `ERR3d`.
    err3d: f64,
    /// `ERR2d`.
    err2d: f64,
    /// `FirstP`.
    first_p: i32,
    /// `LastP`.
    last_p: i32,
    /// `myConstraints`.
    my_constraints: Vec<ConstraintCouple>,
    /// `mylambda1`.
    mylambda1: f64,
    /// `mylambda2`.
    mylambda2: f64,
}

impl<L: LineTool> BspFunction<L> {
    /// `AppParCurves_BSpFunction(SSP, FirstPoint, LastPoint, TheConstraints,
    /// Parameters, Knots, Mults, NbPol)` (`gxx:27-140`).
    pub fn new(
        ssp: &L,
        first_point: i32,
        last_point: i32,
        the_constraints: &[ConstraintCouple],
        parameters: &MathVector,
        knots: &[f64],
        mults: &[i32],
        nb_pol: i32,
    ) -> Self {
        let nb_cu = ssp.nb_p3d() + ssp.nb_p2d();
        let ncu = nb_cu.max(0) as usize;
        let rl = first_point.max(0) as usize;
        let ru = last_point.max(0) as usize;
        let ncol = nb_pol.max(0) as usize;

        let mut this = Self {
            done: false,
            my_multi_line: ssp.clone(),
            my_multi_bsp_curve: MultiBSpCurve::with_nb_poles(ncol),
            my_parameters: MathVector::new(parameters.lower(), parameters.upper()),
            f_val: 0.0,
            val_grad_f: MathVector::new(rl, ru),
            my_f: MathMatrix::with_init(rl, ru, 1, ncu, 0.0),
            ptlx: MathMatrix::with_init(rl, ru, 1, ncu, 0.0),
            ptly: MathMatrix::with_init(rl, ru, 1, ncu, 0.0),
            ptlz: MathMatrix::with_init(rl, ru, 1, ncu, 0.0),
            a: MathMatrix::new(rl, ru, 1, ncol),
            da: MathMatrix::new(rl, ru, 1, ncol),
            my_least_square: LeastSquare::with_knots(
                ssp,
                knots,
                mults,
                first_point,
                last_point,
                Self::first_constraint(the_constraints, first_point),
                Self::last_constraint(the_constraints, last_point),
                nb_pol,
            ),
            contraintes: false,
            nbp: 0,
            nb_cu,
            adeb: 0,
            afin: 0,
            nbpoles: nb_pol,
            tabdim: vec![0; ncu],
            err3d: 0.0,
            err2d: 0.0,
            first_p: first_point,
            last_p: last_point,
            my_constraints: the_constraints.to_vec(),
            mylambda1: 0.0,
            mylambda2: 0.0,
        };

        // `gxx:55-57`.
        for i in parameters.lower()..=parameters.upper() {
            this.my_parameters.set_value(i, parameters.value(i));
        }

        this.first_p = first_point;
        this.last_p = last_point;
        this.nbp = this.last_p - this.first_p + 1;
        this.adeb = this.first_p;
        this.afin = this.last_p;
        this.nbpoles = nb_pol;
        this.my_multi_bsp_curve.set_knots(knots);
        this.my_multi_bsp_curve.set_multiplicities(mults);
        this.contraintes = false;

        // `gxx:74-92`.
        for couple in the_constraints {
            let cons = couple.constraint();
            let myindex = couple.index();
            if myindex == this.first_p {
                if at_least_pass(cons) {
                    this.adeb += 1;
                }
            } else if myindex == this.last_p {
                if at_least_pass(cons) {
                    this.afin -= 1;
                }
            } else if at_least_pass(cons) {
                this.contraintes = true;
            }
        }

        let nb3d = ssp.nb_p3d();
        let _nb2d = ssp.nb_p2d();

        if this.contraintes {
            for i in 1..=nb_cu {
                if i <= nb3d {
                    this.tabdim[(i - 1) as usize] = 3;
                } else {
                    this.tabdim[(i - 1) as usize] = 2;
                }
            }

            // `gxx:119-138`. `TabP2d(j)` is read with the full 1..NbCu index in
            // OCCT while the array only holds `nb2d` slots; the intended slot
            // is `j - nb3d`, the same spelling the right-hand sides use.
            for i in this.first_p..=this.last_p {
                let (tab_p, tab_p2d) = ssp.value(i);
                for j in 1..=nb_cu {
                    let idx = j as usize;
                    if this.tabdim[(j - 1) as usize] == 3 {
                        let p = tab_p[idx - 1];
                        this.ptlx.set_value(i as usize, idx, p.x());
                        this.ptly.set_value(i as usize, idx, p.y());
                        this.ptlz.set_value(i as usize, idx, p.z());
                    } else {
                        let slot = (j - nb3d).max(1) as usize;
                        let p = tab_p2d[slot - 1];
                        this.ptlx.set_value(i as usize, idx, p.x());
                        this.ptly.set_value(i as usize, idx, p.y());
                    }
                }
            }
        }

        this
    }

    /// `FirstConstraint(TheConstraints, FirstPoint)` (`gxx:142-162`).
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

    /// `LastConstraint(TheConstraints, LastPoint)` (`gxx:164-184`).
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

    /// `SetFirstLambda(l1)` (`gxx:242-245`).
    pub fn set_first_lambda(&mut self, l1: f64) {
        self.mylambda1 = l1;
    }

    /// `SetLastLambda(l2)` (`gxx:247-250`).
    pub fn set_last_lambda(&mut self, l2: f64) {
        self.mylambda2 = l2;
    }

    /// `Value(X, F)` (`gxx:186-211`).
    pub fn value(&mut self, x: &MathVector, f: &mut f64) -> bool {
        self.my_parameters = x.clone();

        // `MyLeastSquare.Perform(myParameters, mylambda1, mylambda2)`
        // (`gxx:193`).
        self.my_least_square
            .perform_parameters_lambda(&self.my_parameters, self.mylambda1, self.mylambda2);
        if !self.my_least_square.is_done() {
            self.done = false;
            return false;
        }
        if !self.contraintes {
            let (fval, err3d, err2d) = self.my_least_square.error();
            self.f_val = fval;
            self.err3d = err3d;
            self.err2d = err2d;
            *f = fval;
        }
        // else: empty arm (`gxx:205-209`).
        true
    }

    /// `Perform(X)` (`gxx:213-240`).
    fn perform(&mut self, x: &MathVector) {
        self.my_parameters = x.clone();
        self.my_least_square
            .perform_parameters_lambda(&self.my_parameters, self.mylambda1, self.mylambda2);

        if !self.my_least_square.is_done() {
            self.done = false;
            return;
        }

        for j in self.my_parameters.lower()..=self.my_parameters.upper() {
            self.val_grad_f.set_value(j, 0.0);
        }

        if !self.contraintes {
            let (f, err3d, err2d) = self.my_least_square.error_gradient(&mut self.val_grad_f);
            self.f_val = f;
            self.err3d = err3d;
            self.err2d = err2d;
        }
        // else: empty arm (`gxx:238-240`).
    }

    /// `NbVariables()` (`gxx:252-255`).
    pub fn nb_variables(&self) -> usize {
        self.nbp.max(0) as usize
    }

    /// `Gradient(X, G)` (`gxx:257-265`).
    pub fn gradient(&mut self, x: &MathVector, g: &mut MathVector) -> bool {
        self.perform(x);
        *g = self.val_grad_f.clone();
        true
    }

    /// `Values(X, F, G)` (`gxx:267-290`).
    pub fn values(&mut self, x: &MathVector, f: &mut f64, g: &mut MathVector) -> bool {
        self.perform(x);
        *f = self.f_val;
        *g = self.val_grad_f.clone();
        true
    }

    /// `CurveValue()` (`gxx:292-297`).
    pub fn curve_value(&mut self) -> &MultiBSpCurve {
        if !self.contraintes {
            self.my_multi_bsp_curve = self.my_least_square.bspline_value().clone();
        }
        &self.my_multi_bsp_curve
    }

    /// `Error(IPoint, CurveIndex)` (`gxx:299-307`).
    pub fn error(&mut self, ipoint: i32, curve_index: i32) -> f64 {
        if !self.contraintes {
            return self
                .my_least_square
                .distance()
                .value(ipoint as usize, curve_index as usize);
        }
        self.my_f
            .value(ipoint as usize, curve_index as usize)
            .sqrt()
    }

    /// `MaxError3d()` (`gxx:309-312`).
    pub fn max_error3d(&self) -> f64 {
        self.err3d
    }

    /// `MaxError2d()` (`gxx:314-317`).
    pub fn max_error2d(&self) -> f64 {
        self.err2d
    }

    /// `NewParameters()` (`gxx:319-322`).
    pub fn new_parameters(&self) -> &MathVector {
        &self.my_parameters
    }

    /// `FunctionMatrix()` (`gxx:324-327`).
    pub fn function_matrix(&self) -> &MathMatrix {
        self.my_least_square.function_matrix()
    }

    /// `DerivativeFunctionMatrix()` (`gxx:329-332`).
    pub fn derivative_function_matrix(&self) -> &MathMatrix {
        self.my_least_square.derivative_function_matrix()
    }

    /// `Index()` (`gxx:334-337`).
    pub fn k_index(&self) -> &MathIntVector {
        self.my_least_square.k_index()
    }

    /// The stored `MultiLine` (OCCT: `MyMultiLine`).
    pub fn multi_line(&self) -> &L {
        &self.my_multi_line
    }

    /// The stored knot vector (OCCT: `MyMultiBSpCurve.Knots()`).
    pub fn nb_poles(&self) -> i32 {
        self.nbpoles
    }
}








