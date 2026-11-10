//! `AppParCurves_Function` — the objective `F = sum ||C(ui) - Pt(ui)||^2` and
//! its gradient, as a function of the parametrisation `ui`.
//!
//! Source: `ModelingData/TKGeomBase/AppParCurves/AppParCurves_Function.gxx:35-695`
//! and the generated declaration
//! `AppDef_ParFunctionOfMyGradientOfCompute.hxx:37-131` (OCCT 8.0.0). The
//! `.gxx` is instantiated for each `AppDef` / `GeomInt` / `BRepApprox` driver;
//! they differ only in the substituted `MultiLine` / `LineTool` names, so the
//! port keeps the single `AppDef` instantiation and collapses the generics
//! onto [`crate::app_def::MultiLine`] / [`crate::app_def::MyLineTool`].
//!
//! Layout:
//! - [`mod@value`] — `Value` (`gxx:189-293`),
//! - [`mod@gradient`] — `Perform` (`gxx:295-646`), the shared body of
//!   `Gradient` / `Values`,
//! - [`mod@matrix_ops`] — the bound-correct `math_Matrix` products.
//!
//! Translation notes:
//! - `MyLeastSquare` is `AppDef_ParLeastSquareOfMyGradientOfCompute`, i.e.
//!   [`LeastSquare`] with the Bezier constructors; `BezierValue()` therefore
//!   never raises `Standard_NoSuchObject` here.
//! - `ResolCons` is `AppDef_ResConstraintOfMyGradientOfCompute`, ported as
//!   [`ResolConstraint`].
//! - The C++ `Done` member (`hxx:105`) is written on every failure path and
//!   never read or exported, so the port keeps it as a write-only mirror.
//! - `MyMultiCurve`'s size is `Deg + 1` (`gxx:44`); `Perform` overwrites it
//!   from `BezierValue()` (`gxx:337`, `gxx:403`).
//! - `PTLX` / `PTLY` / `PTLZ` are only filled when `Contraintes` is true
//!   (`gxx:109-142`), so the unconstrained paths never read them.

use occt_math::{MathMatrix, MathVector};

use crate::app_def::MultiLine;
use crate::approx::line_tool::LineTool;
use crate::app_par_curves::constraint::{Constraint, ConstraintCouple};
use crate::app_par_curves::least_square::LeastSquare;
use crate::app_par_curves::multi_curve::MultiCurve;

mod gradient;
mod matrix_ops;
mod value;

/// C++ `Cons >= AppParCurves_PassPoint` (`Constraint` declaration order).
pub(super) fn at_least_pass(c: Constraint) -> bool {
    matches!(
        c,
        Constraint::PassPoint | Constraint::TangencyPoint | Constraint::CurvaturePoint
    )
}

/// `AppParCurves_Function` (`AppDef_ParFunctionOfMyGradientOfCompute.hxx:37-131`).
#[derive(Debug, Clone)]
pub struct ParFunction<L: LineTool = MultiLine> {
    /// `Done` (`hxx:105`): written on the failure paths only.
    #[allow(dead_code)]
    done: bool,
    /// `MyMultiLine` (`hxx:106`).
    my_multi_line: L,
    /// `MyMultiCurve` (`hxx:107`).
    my_multi_curve: MultiCurve,
    /// `Degre` (`hxx:108`).
    degre: i32,
    /// `myParameters` (`hxx:109`).
    my_parameters: MathVector,
    /// `FVal` (`hxx:110`): uninitialised in OCCT until `Value` / `Perform`.
    f_val: f64,
    /// `ValGrad_F` (`hxx:111`).
    val_grad_f: MathVector,
    /// `MyF` (`hxx:112`): the per-point, per-curve squared distance.
    my_f: MathMatrix,
    /// `PTLX` (`hxx:113`).
    ptlx: MathMatrix,
    /// `PTLY` (`hxx:114`).
    ptly: MathMatrix,
    /// `PTLZ` (`hxx:115`).
    ptlz: MathMatrix,
    /// `A` (`hxx:116`): the Bernstein matrix.
    a: MathMatrix,
    /// `DA` (`hxx:117`): its derivative.
    da: MathMatrix,
    /// `MyLeastSquare` (`hxx:118`).
    my_least_square: LeastSquare,
    /// `Contraintes` (`hxx:119`): true when an interior point is constrained.
    contraintes: bool,
    /// `NbP` (`hxx:120`).
    nbp: i32,
    /// `NbCu` (`hxx:121`).
    nb_cu: i32,
    /// `Adeb` (`hxx:122`).
    adeb: i32,
    /// `Afin` (`hxx:123`).
    afin: i32,
    /// `tabdim` (`hxx:124`): 3 for the 3d curves, 2 for the 2d ones, indexed
    /// from 0 (`gxx:107`, `gxx:115-117`).
    tabdim: Vec<i32>,
    /// `ERR3d` (`hxx:126`).
    err3d: f64,
    /// `ERR2d` (`hxx:127`).
    err2d: f64,
    /// `FirstP` (`hxx:128`).
    first_p: i32,
    /// `LastP` (`hxx:129`).
    last_p: i32,
    /// `myConstraints` (`hxx:130`).
    my_constraints: Vec<ConstraintCouple>,
}

impl<L: LineTool> ParFunction<L> {
    /// `AppDef_ParFunctionOfMyGradientOfCompute(SSP, FirstPoint, LastPoint,
    /// TheConstraints, Parameters, Deg)` (`gxx:35-144`).
    pub fn new(
        ssp: &L,
        first_point: i32,
        last_point: i32,
        the_constraints: &[ConstraintCouple],
        parameters: &MathVector,
        deg: i32,
    ) -> Self {
        let nb_cu = ssp.nb_p3d() + ssp.nb_p2d();
        let first_cons = Self::first_constraint(the_constraints, first_point);
        let last_cons = Self::last_constraint(the_constraints, last_point);
        let rl = first_point.max(0) as usize;
        let ru = last_point.max(0) as usize;
        let ncu = nb_cu.max(0) as usize;
        let ncol = (deg + 1).max(0) as usize;

        let mut this = Self {
            done: false,
            my_multi_line: ssp.clone(),
            my_multi_curve: MultiCurve::with_nb_poles(ncol),
            degre: deg,
            my_parameters: MathVector::new(parameters.lower(), parameters.upper()),
            f_val: 0.0,            val_grad_f: MathVector::new(rl, ru),
            my_f: MathMatrix::with_init(rl, ru, 1, ncu, 0.0),
            ptlx: MathMatrix::with_init(rl, ru, 1, ncu, 0.0),
            ptly: MathMatrix::with_init(rl, ru, 1, ncu, 0.0),
            ptlz: MathMatrix::with_init(rl, ru, 1, ncu, 0.0),
            a: MathMatrix::new(rl, ru, 1, ncol),
            da: MathMatrix::new(rl, ru, 1, ncol),
            my_least_square: LeastSquare::new(
                ssp,
                first_point,
                last_point,
                first_cons,
                last_cons,
                deg + 1,
            ),
            contraintes: false,
            nbp: 0,
            nb_cu,
            adeb: 0,
            afin: 0,
            tabdim: vec![0; ncu],
            err3d: 0.0,
            err2d: 0.0,
            first_p: first_point,
            last_p: last_point,
            my_constraints: the_constraints.to_vec(),
        };

        // `for (i = Parameters.Lower(); ...) myParameters(i) = Parameters(i);`
        // (`gxx:60-62`); the member was sized from `Parameters`' own bounds.
        for i in parameters.lower()..=parameters.upper() {
            this.my_parameters.set_value(i, parameters.value(i));
        }

        this.first_p = first_point;
        this.last_p = last_point;
        this.nbp = this.last_p - this.first_p + 1;
        this.adeb = this.first_p;
        this.afin = this.last_p;
        this.degre = deg;
        this.contraintes = false;

        // `gxx:71-94`: the points carrying a constraint move the fitted
        // window and/or raise `Contraintes`.
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

        // `gxx:97-99`: `nb3d` / `nb2d` only size the least-squares arrays,
        // which the `LeastSquare` constructor already did; `nb2d` itself is
        // not read again here.
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

            // `gxx:119-141`: cache the constrained points. `TabP2d(j)` is read
            // with the full 1..NbCu index in OCCT while the array only holds
            // `nb2d` slots; the intended slot is `j - nb3d` (the same
            // spelling the right-hand sides use at `gxx:270`), which is what
            // the port reads.
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

    /// `FirstConstraint(TheConstraints, FirstPoint)` (`gxx:144-165`).
    pub fn first_constraint(the_constraints: &[ConstraintCouple], first_point: i32) -> Constraint {
        let mut cons = Constraint::NoConstraint;
        for couple in the_constraints {
            cons = couple.constraint();
            if couple.index() == first_point {
                break;
            }
        }
        cons
    }

    /// `LastConstraint(TheConstraints, LastPoint)` (`gxx:167-185`).
    pub fn last_constraint(the_constraints: &[ConstraintCouple], last_point: i32) -> Constraint {
        let mut cons = Constraint::NoConstraint;
        for couple in the_constraints {
            cons = couple.constraint();
            if couple.index() == last_point {
                break;
            }
        }
        cons
    }

    /// `Value(X, F)` (`gxx:189-293`).
    pub fn value(&mut self, x: &MathVector, f: &mut f64) -> bool {
        value::value(self, x, f)
    }

    /// `Gradient(X, G)` (`gxx:651-661`).
    pub fn gradient(&mut self, x: &MathVector, g: &mut MathVector) -> bool {
        self.perform(x);
        *g = self.val_grad_f.clone();
        true
    }

    /// `Values(X, F, G)` (`gxx:663-668`).
    pub fn values(&mut self, x: &MathVector, f: &mut f64, g: &mut MathVector) -> bool {
        self.perform(x);
        *f = self.f_val;
        *g = self.val_grad_f.clone();
        true
    }

    /// `NbVariables()` (`gxx:648-651`): the number of MultiPoints.
    pub fn nb_variables(&self) -> usize {
        self.nbp.max(0) as usize
    }

    /// `NewParameters()` (`gxx:692-695`).
    pub fn new_parameters(&self) -> &MathVector {
        &self.my_parameters
    }

    /// `CurveValue()` (`gxx:670-677`): recomputes the Bezier poles when the
    /// objective is unconstrained.
    pub fn curve_value(&mut self) -> &MultiCurve {
        if !self.contraintes {
            self.my_multi_curve = self.my_least_square.bezier_value();
        }
        &self.my_multi_curve
    }

    /// `Error(IPoint, CurveIndex)` (`gxx:679-682`).
    pub fn error(&self, ipoint: i32, curve_index: i32) -> f64 {
        self.my_f
            .value(ipoint as usize, curve_index as usize)
            .sqrt()
    }

    /// `MaxError3d()` (`gxx:684-687`).
    pub fn max_error3d(&self) -> f64 {
        self.err3d
    }

    /// `MaxError2d()` (`gxx:689-692`).
    pub fn max_error2d(&self) -> f64 {
        self.err2d
    }

    /// `Perform(X)` (`gxx:295-646`): shared by `Gradient` and `Values`.
    fn perform(&mut self, x: &MathVector) {
        gradient::perform(self, x);
    }
}








