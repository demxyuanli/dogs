//! `Approx_BSplComputeLine` (`AppDef_BSplineCompute`) — the B-spline
//! degree-escalating approximation driver that `AppBlend_AppSurf` drives per
//! surface row.
//!
//! Source: `ModelingData/TKGeomBase/Approx/Approx_BSplComputeLine.gxx:1-1458`
//! (instantiated for `AppDef` as `AppDef_BSplineCompute`), declaration
//! `AppDef_BSplineCompute.hxx:47-257`.
//!
//! Layout:
//! - [`mod@perform`] — `Perform` (`gxx:636-788`) + `FindRealConstraints`
//!   (`gxx:1206-1281`),
//! - [`mod@compute`] — `Compute` (`gxx:896-1080`),
//! - [`mod@parameters`] — `Parameters` (`gxx:811-894`),
//! - [`mod@interpol`] — `Interpol` (`gxx:1284-1422`),
//! - [`mod@tangency`] — `FirstTangencyVector` / `LastTangencyVector` /
//!   `TangencyVector` / `SearchFirstLambda` / `SearchLastLambda`
//!   (`gxx:197-490`, `gxx:1427-1458`).
//!
//! Translation notes:
//! - `NCollection_HArray1<double>` (`myParameters`, `myfirstParam`, `myknots`)
//!   and `HArray1<int>` (`mymults`) become `(lower, Vec)` pairs so OCCT's
//!   `Lower()`-based indexing stays visible, as in
//!   [`ComputeLine`](super::compute_line::ComputeLine).
//! - The generic `LineTool` parameter is fixed to
//!   [`MyLineTool`](crate::app_def::MyLineTool), matching the `AppDef`
//!   instantiation; `MultiLine` is `AppDef_MultiLine`.
//! - `Approx_BSpParLeastSquareOfMyBSplGradient` is
//!   [`LeastSquare`](crate::app_par_curves::LeastSquare),
//!   `Approx_MyBSplGradient` is [`BspGradient`](crate::app_par_curves::BspGradient)
//!   and `Approx_MyGradientbis` is [`Gradient`](crate::app_par_curves::Gradient).
//! - `ComputeCurve` (`hxx:192-194`) is declared but has no definition in the
//!   `.gxx` and no call site; it stays UNPORTED by omission.
//! - The `mydebug` / `DUMP` blocks (`gxx:33-191`, `gxx:638-641`, `gxx:783-786`)
//!   are `OCCT_DEBUG && DRAW && !WNT` only and are not ported.

mod compute;
mod interpol;
mod parameters;
mod perform;
mod tangency;

use crate::approx::line_tool::LineTool;
use crate::app_par_curves::constraint::{Constraint, ConstraintCouple};
use crate::app_par_curves::multi_bspline_curve::MultiBSpCurve;
use occt_core::approx::ApproxParametrizationType;
use occt_math::MathVector;

/// `RealLast()` (`Standard_Real.hxx:99`).
const REAL_LAST: f64 = f64::MAX;

/// `Approx_BSplComputeLine` (`hxx:47-257`).
#[derive(Debug, Clone)]
pub struct BsplComputeLine {
    /// The TheMultiLineTool template parameter (LineTool); the line itself
    /// only enters through Perform / Interpol.
    /// `TheMultiBSpCurve` (`hxx:229`).
    the_multi_bsp_curve: MultiBSpCurve,
    /// `alldone` (`hxx:230`).
    alldone: bool,
    /// `tolreached` (`hxx:231`).
    tolreached: bool,
    /// `Par` (`hxx:232`).
    par: ApproxParametrizationType,
    /// `myParameters` (`hxx:233`).
    my_parameters: Option<(i32, Vec<f64>)>,
    /// `myfirstParam` (`hxx:234`).
    my_first_param: Option<(i32, Vec<f64>)>,
    /// `myknots` (`hxx:235`).
    my_knots: Option<(i32, Vec<f64>)>,
    /// `mymults` (`hxx:236`).
    my_mults: Option<(i32, Vec<i32>)>,
    /// `myhasknots` (`hxx:237`).
    my_has_knots: bool,
    /// `myhasmults` (`hxx:238`).
    my_has_mults: bool,
    /// `myConstraints` (`hxx:239`): always the two couples
    /// `(FirstPoint, realfirstC)` and `(LastPoint, reallastC)` inserted by
    /// `Perform` (`gxx:656-658`).
    my_constraints: Vec<ConstraintCouple>,
    /// `mydegremin` (`hxx:240`).
    my_degremin: i32,
    /// `mydegremax` (`hxx:241`).
    my_degremax: i32,
    /// `mytol3d` (`hxx:242`).
    my_tol3d: f64,
    /// `mytol2d` (`hxx:243`).
    my_tol2d: f64,
    /// `currenttol3d` (`hxx:244`).
    current_tol3d: f64,
    /// `currenttol2d` (`hxx:245`).
    current_tol2d: f64,
    /// `mycut` (`hxx:246`).
    my_cut: bool,
    /// `mysquares` (`hxx:247`).
    my_squares: bool,
    /// `myitermax` (`hxx:248`).
    my_iter_max: i32,
    /// `myfirstC` (`hxx:249`).
    my_first_c: Constraint,
    /// `mylastC` (`hxx:250`).
    my_last_c: Constraint,
    /// `realfirstC` (`hxx:251`).
    real_first_c: Constraint,
    /// `reallastC` (`hxx:252`).
    real_last_c: Constraint,
    /// `mycont` (`hxx:253`).
    my_cont: i32,
    /// `mylambda1` (`hxx:254`).
    my_lambda1: f64,
    /// `mylambda2` (`hxx:255`).
    my_lambda2: f64,
    /// `myPeriodic` (`hxx:256`).
    my_periodic: bool,
}

/// The shared member initialisation of the four constructors
/// (`gxx:505-531`, `:544-568`, `:585-598`, `:615-630`).
#[allow(clippy::too_many_arguments)]
fn empty(
    degreemin: i32,
    degreemax: i32,
    tolerance3d: f64,
    tolerance2d: f64,
    nb_iterations: i32,
    cutting: bool,
    parametrization: ApproxParametrizationType,
    squares: bool,
) -> BsplComputeLine {
    BsplComputeLine {
        the_multi_bsp_curve: MultiBSpCurve::new(),
        alldone: false,
        tolreached: false,
        par: parametrization,
        my_parameters: None,
        my_first_param: None,
        my_knots: None,
        my_mults: None,
        my_has_knots: false,
        my_has_mults: false,
        my_constraints: vec![ConstraintCouple::new(), ConstraintCouple::new()],
        my_degremin: degreemin,
        my_degremax: degreemax,
        my_tol3d: tolerance3d,
        my_tol2d: tolerance2d,
        current_tol3d: REAL_LAST,
        current_tol2d: REAL_LAST,
        my_cut: cutting,
        my_squares: squares,
        my_iter_max: nb_iterations,
        my_first_c: Constraint::TangencyPoint,
        my_last_c: Constraint::TangencyPoint,
        real_first_c: Constraint::NoConstraint,
        real_last_c: Constraint::NoConstraint,
        my_cont: -1,
        my_lambda1: 0.0,
        my_lambda2: 0.0,
        my_periodic: false,
    }
}

impl BsplComputeLine {
    /// `Approx_BSplComputeLine(Line, Parameters, degreemin, degreemax,
    /// Tolerance3d, Tolerance2d, NbIterations, cutting, Squares)`
    /// (`gxx:496-531`): the first parametrisation is caller-supplied and
    /// `Perform` runs immediately.
    #[allow(clippy::too_many_arguments)]
    pub fn with_line_and_parameters<L: LineTool>(
        line: &L,
        parameters: &MathVector,
        degreemin: i32,
        degreemax: i32,
        tolerance3d: f64,
        tolerance2d: f64,
        nb_iterations: i32,
        cutting: bool,
        squares: bool,
    ) -> Self {
        let mut this = empty(
            degreemin,
            degreemax,
            tolerance3d,
            tolerance2d,
            nb_iterations,
            cutting,
            ApproxParametrizationType::IsoParametric,
            squares,
        );
        this.my_first_param = Some((
            parameters.lower() as i32,
            (parameters.lower()..=parameters.upper())
                .map(|i| parameters.value(i))
                .collect(),
        ));
        this.perform(line);
        this
    }

    /// `Approx_BSplComputeLine(Parameters, degreemin, degreemax, Tolerance3d,
    /// Tolerance2d, NbIterations, cutting, Squares)` (`gxx:534-568`): no
    /// `Perform` — the line arrives later.
    #[allow(clippy::too_many_arguments)]
    pub fn from_parameters(
        parameters: &MathVector,
        degreemin: i32,
        degreemax: i32,
        tolerance3d: f64,
        tolerance2d: f64,
        nb_iterations: i32,
        cutting: bool,
        squares: bool,
    ) -> Self {
        let mut this = empty(
            degreemin,
            degreemax,
            tolerance3d,
            tolerance2d,
            nb_iterations,
            cutting,
            ApproxParametrizationType::IsoParametric,
            squares,
        );
        this.my_first_param = Some((
            parameters.lower() as i32,
            (parameters.lower()..=parameters.upper())
                .map(|i| parameters.value(i))
                .collect(),
        ));
        this
    }

    /// `Approx_BSplComputeLine(degreemin, degreemax, Tolerance3d, Tolerance2d,
    /// NbIterations, cutting, parametrization, Squares)` (`gxx:572-599`): no
    /// `Perform`.
    #[allow(clippy::too_many_arguments)]
    pub fn with_settings(
        degreemin: i32,
        degreemax: i32,
        tolerance3d: f64,
        tolerance2d: f64,
        nb_iterations: i32,
        cutting: bool,
        parametrization: ApproxParametrizationType,
        squares: bool,
    ) -> Self {
        empty(
            degreemin,
            degreemax,
            tolerance3d,
            tolerance2d,
            nb_iterations,
            cutting,
            parametrization,
            squares,
        )
    }

    /// `Approx_BSplComputeLine(Line, degreemin, degreemax, Tolerance3d,
    /// Tolerance2d, NbIterations, cutting, parametrization, Squares)`
    /// (`gxx:602-632`).
    #[allow(clippy::too_many_arguments)]
    pub fn from_line<L: LineTool>(
        line: &L,
        degreemin: i32,
        degreemax: i32,
        tolerance3d: f64,
        tolerance2d: f64,
        nb_iterations: i32,
        cutting: bool,
        parametrization: ApproxParametrizationType,
        squares: bool,
    ) -> Self {
        let mut this = empty(
            degreemin,
            degreemax,
            tolerance3d,
            tolerance2d,
            nb_iterations,
            cutting,
            parametrization,
            squares,
        );
        this.perform(line);
        this
    }

    /// `Perform(Line)` (`gxx:636-788`).
    pub fn perform<L: LineTool>(&mut self, line: &L) {
        perform::perform(self, line);
    }

    /// `Interpol(Line)` (`gxx:1284-1422`).
    pub fn interpol<L: LineTool>(&mut self, line: &L) {
        interpol::interpol(self, line);
    }

    /// `Init(degreemin, degreemax, Tolerance3d, Tolerance2d, NbIterations,
    /// cutting, parametrization, Squares)` (`gxx:1126-1144`).
    #[allow(clippy::too_many_arguments)]
    pub fn init(
        &mut self,
        degreemin: i32,
        degreemax: i32,
        tolerance3d: f64,
        tolerance2d: f64,
        nb_iterations: i32,
        cutting: bool,
        parametrization: ApproxParametrizationType,
        squares: bool,
    ) {
        self.my_degremin = degreemin;
        self.my_degremax = degreemax;
        self.my_tol3d = tolerance3d;
        self.my_tol2d = tolerance2d;
        self.par = parametrization;
        self.my_squares = squares;
        self.my_cut = cutting;
        self.my_iter_max = nb_iterations;
    }

    /// `SetParameters(ThePar)` (`gxx:1084-1091`).
    pub fn set_parameters(&mut self, the_par: &MathVector) {
        self.my_first_param = Some((
            the_par.lower() as i32,
            (the_par.lower()..=the_par.upper())
                .map(|i| the_par.value(i))
                .collect(),
        ));
    }

    /// `SetKnots(Knots)` (`gxx:1095-1103`).
    pub fn set_knots(&mut self, knots: &[f64]) {
        self.my_has_knots = true;
        self.my_knots = Some((1, knots.to_vec()));
    }

    /// `SetKnotsAndMultiplicities(Knots, Mults)` (`gxx:1107-1123`).
    pub fn set_knots_and_multiplicities(&mut self, knots: &[f64], mults: &[i32]) {
        self.my_has_knots = true;
        self.my_has_mults = true;
        self.my_knots = Some((1, knots.to_vec()));
        self.my_mults = Some((1, mults.to_vec()));
    }

    /// `SetDegrees(degreemin, degreemax)` (`gxx:1148-1152`).
    pub fn set_degrees(&mut self, degreemin: i32, degreemax: i32) {
        self.my_degremin = degreemin;
        self.my_degremax = degreemax;
    }

    /// `SetTolerances(Tolerance3d, Tolerance2d)` (`gxx:1156-1160`).
    pub fn set_tolerances(&mut self, tolerance3d: f64, tolerance2d: f64) {
        self.my_tol3d = tolerance3d;
        self.my_tol2d = tolerance2d;
    }

    /// `SetConstraints(FirstC, LastC)` (`gxx:1164-1169`).
    pub fn set_constraints(&mut self, first_c: Constraint, last_c: Constraint) {
        self.my_first_c = first_c;
        self.my_last_c = last_c;
    }

    /// `SetPeriodic(thePeriodic)` (`gxx:1173-1176`).
    pub fn set_periodic(&mut self, the_periodic: bool) {
        self.my_periodic = the_periodic;
    }

    /// `SetContinuity(C)` (`gxx:1203-1206`).
    pub fn set_continuity(&mut self, c: i32) {
        self.my_cont = c;
    }

    /// `IsAllApproximated()` (`gxx:1180-1183`).
    pub fn is_all_approximated(&self) -> bool {
        self.alldone
    }

    /// `IsToleranceReached()` (`gxx:1187-1190`).
    pub fn is_tolerance_reached(&self) -> bool {
        self.tolreached
    }

    /// `Error(tol3d, tol2d)` (`gxx:1194-1198`).
    pub fn error(&self) -> (f64, f64) {
        (self.current_tol3d, self.current_tol2d)
    }

    /// `Parameters()` (`gxx:792-795`): the parametrisation of the result.
    pub fn parameters(&self) -> &[f64] {
        let Some((_, v)) = &self.my_parameters else {
            panic!("Standard_NullObject: Approx_BSplComputeLine::Parameters");
        };
        v
    }

    /// `Value()` (`gxx:799-802`).
    pub fn value(&self) -> &MultiBSpCurve {
        &self.the_multi_bsp_curve
    }

    /// `ChangeValue()` (`gxx:806-809`).
    pub fn change_value(&mut self) -> &mut MultiBSpCurve {
        &mut self.the_multi_bsp_curve
    }
}









