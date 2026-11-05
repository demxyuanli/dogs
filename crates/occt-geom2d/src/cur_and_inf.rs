//! `LProp_CurAndInf` (`LProp_CurAndInf.cxx:19-131`, `.hxx:24-54`) and
//! `GeomLProp_CurAndInf2d` (`GeomLProp_CurAndInf2d.cxx:24-197`, `.hxx:...`)
//! with the numeric helpers `GeomLProp_FuncCurExt2d`
//! (`GeomLProp_FuncCurExt2d.pxx:25-151`), `GeomLProp_FuncCurNul2d`
//! (`GeomLProp_FuncCurNul2d.pxx:24-103`) and `GeomLProp_NumericCurInf2d`
//! (`GeomLProp_NumericCurInf2d.pxx:33-134`).
//!
//! `Perform` collects the curvature extrema (`dKC/dU = 0`) and the inflection
//! points (`F = (V1 ^ V2) / (|V1| |V2|) = 0`) of a bounded 2D curve, sorted by
//! parameter. Analytic curves answer in closed form, everything else through
//! `math_FunctionRoots` + `math_BracketedRoot`.
//!
//! `MAT2d_CutCurve` (`MAT2d_CutCurve.cxx:39-85`) is the current consumer: it
//! splits a curve at those parameters.

use std::f64::consts::PI;
use std::sync::Arc;

use occt_core::elib::clib;
use occt_core::gp::GpVec2d;
use occt_core::kernel::geomabs::CurveType;
use occt_core::math_bracketed_root::{
    bracketed_root, BRACKETED_ROOT_NB_ITERATIONS, BRACKETED_ROOT_ZEPS,
};
use occt_core::math_fn::{MathFunction, MathFunctionWithDerivative};
use occt_core::math_function_roots::FunctionRoots;
use occt_core::precision::{PCONFUSION, RESOLUTION};

use crate::adaptor_curve;
use crate::curve::Curve2d;

/// `LProp_CIType` (`LProp_CIType.hxx:...`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LPropCIType {
    /// `LProp_Inflection`.
    Inflection,
    /// `LProp_MinCur`: extremum of minimum radius of curvature.
    MinCur,
    /// `LProp_MaxCur`.
    MaxCur,
}

/// `LProp_CurAndInf` (`LProp_CurAndInf.hxx:24-54`): the sorted list of the
/// characteristic parameters of a curve.
#[derive(Debug, Clone, Default)]
pub struct LPropCurAndInf {
    /// `theParams` (`.hxx:50`).
    the_params: Vec<f64>,
    /// `theTypes` (`.hxx:52`).
    the_types: Vec<LPropCIType>,
}

impl LPropCurAndInf {
    /// `LProp_CurAndInf()` (`cxx:22`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `AddInflection(Param)` (`cxx:26-49`).
    pub fn add_inflection(&mut self, param: f64) {
        insert_sorted(&mut self.the_params, &mut self.the_types, param, LPropCIType::Inflection);
    }

    /// `AddExtCur(Param, IsMin)` (`cxx:55-87`).
    pub fn add_ext_cur(&mut self, param: f64, is_min: bool) {
        let type_point = if is_min {
            LPropCIType::MinCur
        } else {
            LPropCIType::MaxCur
        };
        insert_sorted(&mut self.the_params, &mut self.the_types, param, type_point);
    }

    /// `Clear()` (`cxx:92-96`).
    pub fn clear(&mut self) {
        self.the_params.clear();
        self.the_types.clear();
    }

    /// `IsEmpty()` (`cxx:101-104`).
    pub fn is_empty(&self) -> bool {
        self.the_params.is_empty()
    }

    /// `NbPoints()` (`cxx:109-112`).
    pub fn nb_points(&self) -> usize {
        self.the_params.len()
    }

    /// `Parameter(N)` (`cxx:116-123`); panics where OCCT raises
    /// `Standard_OutOfRange`.
    pub fn parameter(&self, n: i32) -> f64 {
        if n < 1 || (n as usize) > self.nb_points() {
            panic!("LProp_CurAndInf::Parameter: Standard_OutOfRange");
        }
        self.the_params[(n as usize) - 1]
    }

    /// `Type(N)` (`cxx:127-134`); panics where OCCT raises
    /// `Standard_OutOfRange`.
    pub fn type_of(&self, n: i32) -> LPropCIType {
        if n < 1 || (n as usize) > self.nb_points() {
            panic!("LProp_CurAndInf::Type: Standard_OutOfRange");
        }
        self.the_types[(n as usize) - 1]
    }
}

/// The shared tail of `AddInflection` (`cxx:27-49`) and `AddExtCur`
/// (`cxx:64-86`): append when the list is empty or the parameter is the new
/// last one, otherwise insert before the first larger parameter.
fn insert_sorted(
    params: &mut Vec<f64>,
    types: &mut Vec<LPropCIType>,
    param: f64,
    type_point: LPropCIType,
) {
    if params.is_empty() {
        params.push(param);
        types.push(type_point);
        return;
    }
    if param > *params.last().unwrap() {
        params.push(param);
        types.push(type_point);
        return;
    }
    for i in 0..params.len() {
        if params[i] > param {
            params.insert(i, param);
            types.insert(i, type_point);
            break;
        }
    }
}

/// `GeomLProp_FuncCurExt2d` (`FuncCurExt2d.pxx:25-151`): `F = dKC/dU`, whose
/// roots are the curvature extrema.
#[derive(Clone)]
pub struct GeomLPropFuncCurExt2d {
    /// `myCurve` (`.pxx:147`).
    my_curve: Arc<dyn Curve2d>,
    /// `myEpsX` (`.pxx:148`).
    my_eps_x: f64,
}

impl GeomLPropFuncCurExt2d {
    /// `GeomLProp_FuncCurExt2d(theCurve, theTol)` (`.pxx:36-39`).
    pub fn new(the_curve: Arc<dyn Curve2d>, the_tol: f64) -> Self {
        Self { my_curve: the_curve, my_eps_x: the_tol }
    }

    /// `Value(theU, theF)` (`.pxx:44-63`).
    fn eval_value(&self, the_u: f64) -> Option<f64> {
        let (_, v1, v2, v3) = self.my_curve.d3(the_u);
        let cross12 = v1.crossed(&v2);
        let cross13 = v1.crossed(&v3);
        let dot12 = v1.dot(&v2);
        let sq_norm1 = v1.square_magnitude();
        let norm1 = sq_norm1.sqrt();
        let norm1_cb = sq_norm1 * norm1;
        let norm1_pow5 = norm1_cb * sq_norm1;

        if norm1_pow5 < RESOLUTION {
            return None;
        }
        Some(cross13 / norm1_cb - 3.0 * cross12 * dot12 / norm1_pow5)
    }

    /// `Values(theU, theF, theD)` (`.pxx:86-100`): forward (or backward at the
    /// upper end) finite difference with step `myEpsX / 100`.
    fn eval_values(&self, the_u: f64) -> Option<(f64, f64)> {
        let mut delta_u = self.my_eps_x / 100.0;
        if the_u + delta_u > self.my_curve.last_parameter() {
            delta_u = -delta_u;
        }
        let f = self.eval_value(the_u)?;
        let f_next = self.eval_value(the_u + delta_u)?;
        Some((f, (f_next - f) / delta_u))
    }

    /// `IsMinKC(theU)` (`.pxx:107-145`): true when `|KC(U)| > |KC(U + eps)|`,
    /// i.e. `U` is a local maximum of the curvature (minimum radius).
    pub fn is_min_kc(&self, the_u: f64) -> bool {
        let mut delta_u = self.my_eps_x;

        // `.pxx:112-124`: curvature at `theU`.
        let (_, v1, v2, _) = self.my_curve.d3(the_u);
        let a_norm1_cb = curvature_denominator(&v1);
        if a_norm1_cb < RESOLUTION {
            return false;
        }
        let a_curv_curr = v1.crossed(&v2) / a_norm1_cb;

        if the_u + delta_u > self.my_curve.last_parameter() {
            delta_u = -delta_u;
        }

        // `.pxx:131-142`: curvature at the probe point.
        let (_, v1, v2, _) = self.my_curve.d3(the_u + delta_u);
        let a_norm1_cb = curvature_denominator(&v1);
        if a_norm1_cb < RESOLUTION {
            return false;
        }
        let a_curv_probe = v1.crossed(&v2) / a_norm1_cb;

        a_curv_curr.abs() > a_curv_probe.abs()
    }
}

/// `||V1||^3` (`.pxx:55`).
fn curvature_denominator(v1: &GpVec2d) -> f64 {
    let sq = v1.square_magnitude();
    sq * sq.sqrt()
}

impl MathFunction for GeomLPropFuncCurExt2d {
    fn value(&mut self, x: f64, f: &mut f64) -> bool {
        match self.eval_value(x) {
            Some(v) => {
                *f = v;
                true
            }
            None => false,
        }
    }
}

impl MathFunctionWithDerivative for GeomLPropFuncCurExt2d {
    fn derivative(&mut self, x: f64, d: &mut f64) -> bool {
        let mut f = 0.0;
        self.values(x, &mut f, d)
    }

    fn values(&mut self, x: f64, f: &mut f64, d: &mut f64) -> bool {
        match self.eval_values(x) {
            Some((fv, dv)) => {
                *f = fv;
                *d = dv;
                true
            }
            None => false,
        }
    }
}

/// `GeomLProp_FuncCurNul2d` (`FuncCurNul2d.pxx:24-103`): the normalized
/// curvature `F = (V1 ^ V2) / (|V1| |V2|)`, whose roots are the inflections.
#[derive(Clone)]
pub struct GeomLPropFuncCurNul2d {
    /// `myCurve` (`.pxx:100`).
    my_curve: Arc<dyn Curve2d>,
}

impl GeomLPropFuncCurNul2d {
    /// `GeomLProp_FuncCurNul2d(theCurve)` (`.pxx:32-35`).
    pub fn new(the_curve: Arc<dyn Curve2d>) -> Self {
        Self { my_curve: the_curve }
    }

    /// `Values(theU, theF, theD)` (`.pxx:66-98`).
    fn eval_values(&self, the_u: f64) -> Option<(f64, f64)> {
        let (_, v1, v2, v3) = self.my_curve.d3(the_u);
        let cross12 = v1.crossed(&v2);
        let cross13 = v1.crossed(&v3);
        let dot12 = v1.dot(&v2);
        let dot23 = v2.dot(&v3);
        let norm1 = v1.magnitude();
        let norm2 = v2.magnitude();

        // `.pxx:80-82`: a vanishing second derivative reports `F = D = 0`.
        if norm2 < 1.0e-4 {
            return Some((0.0, 0.0));
        }
        if norm1 * norm2 < RESOLUTION {
            return None;
        }
        let norm_prod = norm1 * norm2;
        let f = cross12 / norm_prod;
        let d = (cross13 - cross12 * dot12 / (norm1 * norm1) - cross12 * dot23 / (norm2 * norm2))
            / norm_prod;
        Some((f, d))
    }
}

impl MathFunction for GeomLPropFuncCurNul2d {
    fn value(&mut self, x: f64, f: &mut f64) -> bool {
        let mut d = 0.0;
        self.values(x, f, &mut d)
    }
}

impl MathFunctionWithDerivative for GeomLPropFuncCurNul2d {
    fn derivative(&mut self, x: f64, d: &mut f64) -> bool {
        let mut f = 0.0;
        self.values(x, &mut f, d)
    }

    fn values(&mut self, x: f64, f: &mut f64, d: &mut f64) -> bool {
        match self.eval_values(x) {
            Some((fv, dv)) => {
                *f = fv;
                *d = dv;
                true
            }
            None => false,
        }
    }
}

/// `GeomLProp_NumericCurInf2d` (`NumericCurInf2d.pxx:33-134`).
#[derive(Debug, Clone, Default)]
pub struct GeomLPropNumericCurInf2d {
    /// `myIsDone` (`.pxx:131`).
    my_is_done: bool,
}

impl GeomLPropNumericCurInf2d {
    /// `GeomLProp_NumericCurInf2d()` (`.pxx:36-39`).
    pub fn new() -> Self {
        Self { my_is_done: false }
    }

    /// `PerformCurExt(theCurve, theResult)` (`.pxx:43-47`).
    pub fn perform_cur_ext(&mut self, the_curve: &dyn Curve2d, the_result: &mut LPropCurAndInf) {
        let (first, last) = (the_curve.first_parameter(), the_curve.last_parameter());
        self.perform_cur_ext_range(the_curve, first, last, the_result);
    }

    /// `PerformCurExt(theCurve, theUMin, theUMax, theResult)` (`.pxx:58-86`).
    pub fn perform_cur_ext_range(
        &mut self,
        the_curve: &dyn Curve2d,
        the_u_min: f64,
        the_u_max: f64,
        the_result: &mut LPropCurAndInf,
    ) {
        self.my_is_done = true;

        let an_eps_h = 1.0e-4 * (the_u_max - the_u_min);
        let a_tol = PCONFUSION;
        let a_nb_samples = 100;

        let mut a_func = GeomLPropFuncCurExt2d::new(Arc::from(the_curve.clone_dyn()), an_eps_h);
        let a_root_finder = FunctionRoots::new(
            &mut a_func,
            the_u_min,
            the_u_max,
            a_nb_samples,
            an_eps_h,
            an_eps_h,
            an_eps_h,
            0.0,
        );

        if a_root_finder.is_done() {
            for j in 1..=a_root_finder.nb_solutions() {
                let mut a_param = a_root_finder.value(j);
                // Refine the root with bracketed bisection (`.pxx:72-76`).
                let mut refined = a_func.clone();
                let a_refinement = bracketed_root(
                    &mut refined,
                    a_param - an_eps_h,
                    a_param + an_eps_h,
                    a_tol,
                    BRACKETED_ROOT_NB_ITERATIONS,
                    BRACKETED_ROOT_ZEPS,
                );
                if a_refinement.done {
                    a_param = a_refinement.root;
                }
                let is_min_radius = a_func.is_min_kc(a_param);
                the_result.add_ext_cur(a_param, is_min_radius);
            }
        } else {
            self.my_is_done = false;
        }
    }

    /// `PerformInf(theCurve, theResult)` (`.pxx:92-95`).
    pub fn perform_inf(&mut self, the_curve: &dyn Curve2d, the_result: &mut LPropCurAndInf) {
        let (first, last) = (the_curve.first_parameter(), the_curve.last_parameter());
        self.perform_inf_range(the_curve, first, last, the_result);
    }

    /// `PerformInf(theCurve, theUMin, theUMax, theResult)` (`.pxx:104-126`).
    pub fn perform_inf_range(
        &mut self,
        the_curve: &dyn Curve2d,
        the_u_min: f64,
        the_u_max: f64,
        the_result: &mut LPropCurAndInf,
    ) {
        self.my_is_done = true;
        let mut a_func = GeomLPropFuncCurNul2d::new(Arc::from(the_curve.clone_dyn()));
        let an_eps_x = 1.0e-6;
        let an_eps_f = 1.0e-6;
        let a_nb_samples = 30;

        let a_root_finder = FunctionRoots::new(
            &mut a_func,
            the_u_min,
            the_u_max,
            a_nb_samples,
            an_eps_x,
            an_eps_f,
            an_eps_x,
            0.0,
        );

        if a_root_finder.is_done() {
            for j in 1..=a_root_finder.nb_solutions() {
                the_result.add_inflection(a_root_finder.value(j));
            }
        } else {
            self.my_is_done = false;
        }
    }

    /// `IsDone()` (`.pxx:129`).
    pub fn is_done(&self) -> bool {
        self.my_is_done
    }
}

/// `GeomLProp_CurAndInf2d` (`GeomLProp_CurAndInf2d.hxx:...`), a
/// `LProp_CurAndInf` with the analytic / numeric `Perform` drivers.
#[derive(Debug, Clone, Default)]
pub struct GeomLPropCurAndInf2d {
    /// The `LProp_CurAndInf` base fields.
    base: LPropCurAndInf,
    /// `myIsDone` (`cxx:26`).
    my_is_done: bool,
}

impl GeomLPropCurAndInf2d {
    /// `GeomLProp_CurAndInf2d()` (`cxx:24-27`).
    pub fn new() -> Self {
        Self { base: LPropCurAndInf::new(), my_is_done: false }
    }

    /// The inherited `LProp_CurAndInf` accessors.
    pub fn points(&self) -> &LPropCurAndInf {
        &self.base
    }

    /// `Perform(C)` (`cxx:32-37`).
    pub fn perform(&mut self, c: &dyn Curve2d) {
        self.base.clear();
        self.my_is_done = true;
        self.perform_cur_ext(c);
        self.perform_inf(c);
    }

    /// `PerformCurExt(C)` (`cxx:42-47`).
    pub fn perform_cur_ext(&mut self, c: &dyn Curve2d) {
        self.base.clear();
        self.my_is_done = true;
        self.do_perform_cur_ext(c);
    }

    /// `PerformInf(C)` (`cxx:131-136`).
    pub fn perform_inf(&mut self, c: &dyn Curve2d) {
        self.base.clear();
        self.my_is_done = true;
        self.do_perform_inf(c);
    }

    /// `IsDone()` (`cxx:191-194`).
    pub fn is_done(&self) -> bool {
        self.my_is_done
    }

    /// `performCurExt(theCurve)` (`cxx:52-126`).
    fn do_perform_cur_ext(&mut self, the_curve: &dyn Curve2d) {
        let a_curve_type = adaptor_curve::type2d(the_curve);
        let (a_u_first, a_u_last) = (the_curve.first_parameter(), the_curve.last_parameter());
        let mut a_numeric_inf = GeomLPropNumericCurInf2d::new();

        match a_curve_type {
            // `cxx:59-62`: a line has no curvature extremum.
            CurveType::Line => {}
            CurveType::Circle => {}
            CurveType::Ellipse => {
                // `cxx:63-78`: the four quarter points of the parameter circle,
                // mapped into the curve's period and clipped to its range.
                let a_u_f_plus2pi = a_u_first + 2.0 * PI;
                let a_quarter_pts = [0.0, PI / 2.0, PI, 3.0 * PI / 2.0];
                let a_is_min = [true, false, true, false];
                for i in 0..4 {
                    let a_u = clib::in_period(a_quarter_pts[i], a_u_first, a_u_f_plus2pi);
                    if a_u_first <= a_u && a_u <= a_u_last {
                        self.base.add_ext_cur(a_u, a_is_min[i]);
                    }
                }
            }
            CurveType::Hyperbola => {
                // `cxx:79-87`.
                if a_u_first <= 0.0 && a_u_last >= 0.0 {
                    self.base.add_ext_cur(0.0, true);
                }
            }
            CurveType::Parabola => {
                // `cxx:88-96`.
                if a_u_first <= 0.0 && a_u_last >= 0.0 {
                    self.base.add_ext_cur(0.0, true);
                }
            }
            CurveType::BSplineCurve => {
                // `cxx:97-118`.
                if the_curve.continuity() >= 5 {
                    a_numeric_inf.perform_cur_ext(the_curve, &mut self.base);
                    self.my_is_done = a_numeric_inf.is_done();
                } else {
                    self.my_is_done = true;
                    let intervals = the_curve.parameter_intervals(5);
                    for pair in intervals.windows(2) {
                        a_numeric_inf.perform_cur_ext_range(
                            the_curve,
                            pair[0],
                            pair[1],
                            &mut self.base,
                        );
                        if !a_numeric_inf.is_done() {
                            self.my_is_done = false;
                        }
                    }
                }
            }
            _ => {
                // `cxx:120-125`.
                a_numeric_inf.perform_cur_ext(the_curve, &mut self.base);
                self.my_is_done = a_numeric_inf.is_done();
            }
        }
    }

    /// `performInf(theCurve)` (`cxx:138-187`).
    fn do_perform_inf(&mut self, the_curve: &dyn Curve2d) {
        let a_curve_type = adaptor_curve::type2d(the_curve);
        let mut a_numeric_inf = GeomLPropNumericCurInf2d::new();

        match a_curve_type {
            // `cxx:147-156`: conics have no inflection point.
            CurveType::Line
            | CurveType::Circle
            | CurveType::Ellipse
            | CurveType::Hyperbola
            | CurveType::Parabola => {}
            CurveType::BSplineCurve => {
                // `cxx:157-179`.
                if the_curve.continuity() >= 5 {
                    a_numeric_inf.perform_inf(the_curve, &mut self.base);
                    self.my_is_done = a_numeric_inf.is_done();
                } else {
                    self.my_is_done = true;
                    let intervals = the_curve.parameter_intervals(5);
                    for pair in intervals.windows(2) {
                        a_numeric_inf.perform_inf_range(the_curve, pair[0], pair[1], &mut self.base);
                        if !a_numeric_inf.is_done() {
                            self.my_is_done = false;
                        }
                    }
                }
            }
            _ => {
                // `cxx:181-186`.
                a_numeric_inf.perform_inf(the_curve, &mut self.base);
                self.my_is_done = a_numeric_inf.is_done();
            }
        }
    }
}
