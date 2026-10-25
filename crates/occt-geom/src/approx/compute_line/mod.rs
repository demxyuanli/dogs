//! `Approx_ComputeLine` (`AppDef_Compute`) — the Bezier degree-escalating
//! approximation driver that `AppBlend_AppSurf` drives per section.
//!
//! Source: `ModelingData/TKGeomBase/Approx/Approx_ComputeLine.gxx:1-1749`
//! (instantiated for `AppDef` as `AppDef_Compute`, `AppDef_Compute_0.cxx`).
//!
//! Translation notes:
//! - `NCollection_HArray1<double>` members (`myfirstParam`, `myParameters`, the
//!   `myPar` entries) become `(lower, Vec<f64>)` pairs so OCCT's `Lower()`-
//!   based indexing stays visible.
//! - The `TheMultiLineTool` template parameter of the `.gxx` is
//!   [`LineTool`](crate::approx::line_tool::LineTool): its associated `Line`
//!   is OCCT's `TheMultiLine`. The default instantiation is
//!   `LineTool = MyLineTool` / `Line = AppDef_MultiLine`, i.e. `AppDef_Compute`
//!   (`AppDef_Compute_0.cxx`); `BRepApprox_TheComputeLineOfApprox` names the
//!   same struct with the `ApproxInt` tool instead.
//! - `Approx_ParLeastSquareOfMyGradient` is
//!   [`LeastSquare`](crate::app_par_curves::LeastSquare) and
//!   `Approx_MyGradient` is [`Gradient`](crate::app_par_curves::Gradient).

pub mod check_multi_curve;
pub mod compute;
pub mod lambda;
pub mod parameters;
pub mod perform;

use occt_core::approx::ApproxParametrizationType;

use crate::app_par_curves::constraint::{Constraint, ConstraintCouple};
use crate::app_par_curves::multi_bspline_curve::MultiBSpCurve;
use crate::app_par_curves::multi_curve::MultiCurve;
use crate::approx::line_tool::LineTool;
use occt_math::MathVector;

/// `Approx_ComputeLine` (`Approx_ComputeLine.hxx:29-232`).
///
/// `LT` is the `TheMultiLineTool` template parameter (`LineTool`); the default
/// names the `AppDef_Compute` instantiation.
#[derive(Debug, Clone)]
pub struct ComputeLine {
    /// `TheMultiLine` never enters the member list (the line arrives per
    /// `myConstraints`: always two couples, `(FirstPoint, myfirstC)` and
    /// `(LastPoint, mylastC)` (`gxx:860-864`).
    my_constraints: Vec<ConstraintCouple>,
    /// `myfirstParam`: the caller-supplied parametrisation, consumed by the
    /// first cut (`gxx:1189-1193`).
    my_first_param: Option<(i32, Vec<f64>)>,
    /// `myfirstC`.
    my_first_c: Constraint,
    /// `mylastC`.
    my_last_c: Constraint,
    /// `Par`.
    par: ApproxParametrizationType,
    /// `mydegremin`.
    my_degremin: i32,
    /// `mydegremax`.
    my_degremax: i32,
    /// `mytol3d`.
    my_tol3d: f64,
    /// `mytol2d`.
    my_tol2d: f64,
    /// `mysquares`: selects the least-square path over the BFGS one.
    my_squares: bool,
    /// `mycut`.
    my_cut: bool,
    /// `myitermax`.
    my_iter_max: i32,
    /// `alldone`.
    alldone: bool,
    /// `tolreached`.
    tolreached: bool,
    /// `TheMultiCurve` (`gxx:880`): the current best approximation.
    the_multi_curve: MultiCurve,
    /// `currenttol3d`.
    current_tol3d: f64,
    /// `currenttol2d`.
    current_tol2d: f64,
    /// `myParameters`.
    my_parameters: Option<(i32, Vec<f64>)>,
    /// `myMultiCurves`.
    my_multi_curves: Vec<MultiCurve>,
    /// `myPar`: the parametrisation of each stored curve.
    my_par: Vec<(i32, Vec<f64>)>,
    /// `Tolers3d`.
    tol_ers3d: Vec<f64>,
    /// `Tolers2d`.
    tol_ers2d: Vec<f64>,
    /// `myMultiLineNb`.
    my_multi_line_nb: i32,
    /// `myIsClear`.
    my_is_clear: bool,
    /// `myspline` (`SplineValue`).
    my_spline: MultiBSpCurve,
}

impl Default for ComputeLine {
    fn default() -> Self {
        Self::with_settings(
            3,
            8,
            1.0e-3,
            1.0e-3,
            0,
            false,
            ApproxParametrizationType::ChordLength,
            true,
        )
    }
}

impl ComputeLine {
    /// The shared member initialisation of the four constructors
    /// (`gxx:735-745`, `:765-778`, `:790-801`, `:815-829`).
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
    ) -> Self {
        Self {
            my_constraints: vec![ConstraintCouple::new(), ConstraintCouple::new()],
            my_first_param: None,
            my_first_c: Constraint::TangencyPoint,
            my_last_c: Constraint::TangencyPoint,
            par: parametrization,
            my_degremin: degreemin,
            my_degremax: degreemax,
            my_tol3d: tolerance3d,
            my_tol2d: tolerance2d,
            my_squares: squares,
            my_cut: cutting,
            my_iter_max: nb_iterations,
            alldone: false,
            tolreached: false,
            the_multi_curve: MultiCurve::new(),
            current_tol3d: 0.0,
            current_tol2d: 0.0,
            my_parameters: None,
            my_multi_curves: Vec::new(),
            my_par: Vec::new(),
            tol_ers3d: Vec::new(),
            tol_ers2d: Vec::new(),
            my_multi_line_nb: 0,
            my_is_clear: false,
            my_spline: MultiBSpCurve::new(),
        }
    }

    /// `Approx_ComputeLine(Line, Parameters, degreemin, degreemax, Tolerance3d,
    /// Tolerance2d, NbIterations, cutting, Squares)` (`gxx:717-747`).
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
        let mut this = Self::empty(
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

    /// `Approx_ComputeLine(Parameters, degreemin, degreemax, Tolerance3d,
    /// Tolerance2d, NbIterations, cutting, Squares)` (`gxx:749-777`): no
    /// `Perform` — the line arrives later through `Perform`.
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
        let mut this = Self::empty(
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

    /// `Approx_ComputeLine(degreemin, degreemax, Tolerance3d, Tolerance2d,
    /// NbIterations, cutting, parametrization, Squares)` (`gxx:779-802`):
    /// no `Perform`.
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
        Self::empty(
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

    /// `Approx_ComputeLine(Line, degreemin, degreemax, Tolerance3d,
    /// Tolerance2d, NbIterations, cutting, parametrization, Squares)`
    /// (`gxx:804-830`).
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
        let mut this = Self::empty(
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

    /// `Parameters(Index)` (`gxx:1221-1224`).
    pub fn parameters(&self, index: i32) -> &[f64] {
        &self.my_par[(index - 1) as usize].1
    }

    /// `NbMultiCurves()` (`gxx:1226-1229`).
    pub fn nb_multi_curves(&self) -> usize {
        self.my_multi_curves.len()
    }

    /// `ChangeValue(Index)` (`gxx:1231-1234`).
    pub fn change_value(&mut self, index: i32) -> &mut MultiCurve {
        &mut self.my_multi_curves[(index - 1) as usize]
    }

    /// `Value(Index)` (`gxx:1236-1239`).
    pub fn value(&self, index: i32) -> &MultiCurve {
        &self.my_multi_curves[(index - 1) as usize]
    }

    /// `SplineValue()` (`gxx:1241-1247`).
    pub fn spline_value(&mut self) -> &MultiBSpCurve {
        let mut trans = crate::approx::m_curves_to_bsp_curve::MCurvesToBSpCurve::new();
        trans.perform_seq(&self.my_multi_curves);
        self.my_spline = trans.value().clone();
        &self.my_spline
    }

    /// `Init(degreemin, degreemax, Tolerance3d, Tolerance2d, NbIterations,
    /// cutting, parametrization, Squares)` (`gxx:1692-1709`).
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

    /// `SetDegrees` (`gxx:1711-1715`).
    pub fn set_degrees(&mut self, degreemin: i32, degreemax: i32) {
        self.my_degremin = degreemin;
        self.my_degremax = degreemax;
    }

    /// `SetTolerances` (`gxx:1717-1721`).
    pub fn set_tolerances(&mut self, tolerance3d: f64, tolerance2d: f64) {
        self.my_tol3d = tolerance3d;
        self.my_tol2d = tolerance2d;
    }

    /// `SetConstraints` (`gxx:1723-1728`).
    pub fn set_constraints(&mut self, first_c: Constraint, last_c: Constraint) {
        self.my_first_c = first_c;
        self.my_last_c = last_c;
    }

    /// `IsAllApproximated` (`gxx:1730-1733`).
    pub fn is_all_approximated(&self) -> bool {
        self.alldone
    }

    /// `IsToleranceReached` (`gxx:1735-1738`).
    pub fn is_tolerance_reached(&self) -> bool {
        self.tolreached
    }

    /// `Error(Index, tol3d, tol2d)` (`gxx:1740-1744`).
    pub fn error(&self, index: i32) -> (f64, f64) {
        (
            self.tol_ers3d[(index - 1) as usize],
            self.tol_ers2d[(index - 1) as usize],
        )
    }

    /// `Parametrization()` (`gxx:1746-1749`).
    pub fn parametrization(&self) -> ApproxParametrizationType {
        self.par
    }
}







