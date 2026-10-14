//! `GeomLib_CheckCurveOnSurface`
//! (`GeomLib_CheckCurveOnSurface.cxx:1-772`,
//! `GeomLib_CheckCurveOnSurface.hxx:23-78`): the maximum distance between a 3D
//! curve and a 2D curve lifted onto a surface.
//!
//! The class maximises `|C3d(t) - S(C2d(t))|` by minimising
//! `TargetFunc(t) = -|C3d(t) - S(C2d(t))|^2` over each sub-interval of
//! `[First, Last]`: `math_PSO` seeds the optimum (`cxx:645-694`) and
//! `math_NewtonMinimum` refines it (`cxx:695-772`). `TargetFunc` (`cxx:52-205`)
//! supplies the exact value, 1st and 2nd derivatives from `D1`/`D2` of the two
//! curves (`cxx:104-147`); `FillSubIntervals` (`cxx:430-640`) sub-divides at
//! the knots of the 3D and 2D B-spline curves.
//!
//! Rust notes:
//! * `Adaptor3d_CurveOnSurface` is passed as `&dyn Curve` (the port's
//!   `meshing::CurveOnSurface`) plus its `GetCurve()` 2D handle (`cxx:363`,
//!   `:381`), the same convention as [`crate::geom_lib::build_curve3d`].
//! * OCCT's `myCurve` is a `GeomAdaptor_Curve(theCurve, theFirst, theLast)`
//!   (`IntTools_Tools.cxx:749-750`). The port's `myCurve` is the curve itself,
//!   whose `first_parameter()` / `last_parameter()` carry the same values at
//!   every OCCT call site (`IntTools_FaceFace.cxx:636-637`,
//!   `BOPTools_AlgoTools2D_1.cxx:90-93`).
//! * `OSD_ThreadPool` (`cxx:386-410`) is not used: `perform` always runs the
//!   per-interval loop serially. `myIsParallel` is kept so `SetParallel`
//!   (`hxx:51-53`) round-trips; the serial result is identical.
//! * `Standard_Failure` guards (`cxx:78-90`, `:110-143`, `:422-425`,
//!   `:628-635`, `:701-771`) have no Rust counterpart and are not needed: the
//!   evaluations below are total.

use occt_core::gp::GpVec;
use occt_core::precision::PCONFUSION;
use occt_geom2d::curve::Curve2d;
use occt_math::{
    MathMatrix, MathMultipleVarFunction, MathMultipleVarFunctionWithGradient,
    MathMultipleVarFunctionWithHessian, MathStatus, MathVector, NewtonMinimum, Pso,
    PsoParticlesPool,
};

use crate::curve::Curve;

/// `RealLast()`.
const REAL_LAST: f64 = f64::MAX;

/// `aMaxKnots` (`cxx:436`): the largest knot array fed to `FillSubIntervals`.
const A_MAX_KNOTS: usize = 101;

/// `GeomLib_CheckCurveOnSurface_TargetFunc` (`cxx:52-205`).
struct TargetFunc<'a> {
    /// `myCurve1` (`cxx:58`): the 3D curve (`Adaptor3d_Curve`).
    curve1: &'a dyn Curve,
    /// `myCurve2` (`cxx:59`): the curve on surface (`Adaptor3d_CurveOnSurface`).
    curve2: &'a dyn Curve,
    /// `myFirst` (`cxx:60`).
    first: f64,
    /// `myLast` (`cxx:61`).
    last: f64,
}

impl<'a> TargetFunc<'a> {
    /// `CheckParameter` (`cxx:196-201`).
    fn check_parameter(&self, param: f64) -> bool {
        self.first <= param && param <= self.last
    }

    /// `FirstParameter()` (`cxx:185-187`).
    fn first_parameter(&self) -> f64 {
        self.first
    }

    /// `LastParameter()` (`cxx:189-191`).
    fn last_parameter(&self) -> f64 {
        self.last
    }

    /// `Value(const double theX, double& theFVal)` (`cxx:74-93`):
    /// `F(X) = -|C3d(X) - C2d(X)|^2`.
    fn value_scalar(&self, x: f64, fval: &mut f64) -> bool {
        if !self.check_parameter(x) {
            return false;
        }
        let p1 = self.curve1.d0(x);
        let p2 = self.curve2.d0(x);
        *fval = -1.0 * p1.square_distance(&p2);
        true
    }

    /// `Derive(const double theX, double& theDeriv1, double* const theDeriv2)`
    /// (`cxx:107-147`).
    fn derive(&self, x: f64, deriv1: &mut f64, deriv2: Option<&mut f64>) -> bool {
        if !self.check_parameter(x) {
            return false;
        }
        //
        // `if (!theDeriv2) { D1 } else { D2 }` (`cxx:121-130`).
        let (p1, dc1, dcc1, p2, dc2, dcc2) = match deriv2 {
            None => {
                let (p1, dc1) = self.curve1.d1(x);
                let (p2, dc2) = self.curve2.d1(x);
                (p1, dc1, GpVec::zero(), p2, dc2, GpVec::zero())
            }
            Some(_) => {
                let (p1, dc1, dcc1) = self.curve1.d2(x);
                let (p2, dc2, dcc2) = self.curve2.d2(x);
                (p1, dc1, dcc1, p2, dc2, dcc2)
            }
        };

        // `gp_Vec aVec1(aP1, aP2)` is `aP2 - aP1`; `aVec2 = aDC2 - aDC1`.
        let vec1 = GpVec::new(p2.x() - p1.x(), p2.y() - p1.y(), p2.z() - p1.z());
        let vec2 = dc2.subtracted(&dc1);
        //
        *deriv1 = -2.0 * vec1.dot(&vec2);

        if let Some(out) = deriv2 {
            let vec3 = dcc2.subtracted(&dcc1);
            *out = -2.0 * (vec2.square_magnitude() + vec1.dot(&vec3));
        }
        true
    }
}

impl MathMultipleVarFunction for TargetFunc<'_> {
    /// `NbVariables()` (`cxx:68`).
    fn nb_variables(&self) -> usize {
        1
    }

    /// `Value(const math_Vector& theX, double& theFVal)` (`cxx:72`).
    fn value(&self, x: &MathVector, f: &mut f64) -> bool {
        self.value_scalar(x.value(1), f)
    }
}

impl MathMultipleVarFunctionWithGradient for TargetFunc<'_> {
    /// `Gradient` (`cxx:100-104`).
    fn gradient(&self, x: &MathVector, g: &mut MathVector) -> bool {
        let mut d1 = 0.0;
        if !self.derive(x.value(1), &mut d1, None) {
            return false;
        }
        g.set_value(1, d1);
        true
    }

    /// `Values(theX, theVal, theGrad)` (`cxx:150-164`).
    fn values(&self, x: &MathVector, f: &mut f64, g: &mut MathVector) -> bool {
        if !self.value(x, f) {
            return false;
        }
        self.gradient(x, g)
    }
}

impl MathMultipleVarFunctionWithHessian for TargetFunc<'_> {
    /// `Values(theX, theVal, theGrad, theHessian)` (`cxx:166-183`).
    fn values_hessian(
        &self,
        x: &MathVector,
        f: &mut f64,
        g: &mut MathVector,
        h: &mut MathMatrix,
    ) -> bool {
        if !self.value(x, f) {
            return false;
        }
        //
        let mut d1 = 0.0;
        let mut d2 = 0.0;
        if !self.derive(x.value(1), &mut d1, Some(&mut d2)) {
            return false;
        }
        g.set_value(1, d1);
        h.set_value(1, 1, d2);
        true
    }
}

/// Run lengths of equal consecutive knots, i.e. `Geom_BSplineCurve::Multiplicities()`
/// derived from the flat `Knots()` sequence
/// (`Geom_BSplineCurve_1.cxx:148-151`, `:380-383`).
fn knot_multiplicities(flat: &[f64]) -> Vec<usize> {
    let mut mults = Vec::new();
    let mut i = 0usize;
    while i < flat.len() {
        let mut j = i + 1;
        while j < flat.len() && flat[j] == flat[i] {
            j += 1;
        }
        mults.push(j - i);
        i = j;
    }
    mults
}

/// `BSplCLib::FirstUKnotIndex` (`BSplCLib.cxx:111-122`): the 1-based index in
/// `Mults` of the first knot whose cumulative multiplicity exceeds `Degree`.
///
/// The walk terminates in range: `sum(Mults) = NbKnots()` is at least
/// `Degree + 2` for any valid B-spline, so `SigmaMult` must exceed `Degree`.
fn first_uknot_index(degree: usize, mults: &[usize]) -> usize {
    let mut index = 0usize;
    let mut sigma_mult = mults[index];

    while sigma_mult <= degree {
        index += 1;
        sigma_mult += mults[index];
    }
    index + 1
}

/// `BSplCLib::LastUKnotIndex` (`BSplCLib.cxx:126-137`).
fn last_uknot_index(degree: usize, mults: &[usize]) -> usize {
    let mut index = mults.len() - 1;
    let mut sigma_mult = mults[index];

    while sigma_mult <= degree {
        index -= 1;
        sigma_mult += mults[index];
    }
    index + 1
}

/// The uniform `aMaxKnots`-point split used when a curve has more knots than
/// the guard allows (`cxx:491-501`, `cxx:541-551`).
fn uniform_knots(first: f64, last: f64) -> Vec<f64> {
    let mut arr = vec![0.0; A_MAX_KNOTS];
    arr[0] = first;
    arr[A_MAX_KNOTS - 1] = last;
    let dt = (last - first) / (A_MAX_KNOTS as f64 - 1.0);
    let mut t = first + dt;
    for slot in arr.iter_mut().take(A_MAX_KNOTS - 1).skip(1) {
        *slot = t;
        t += dt;
    }
    arr
}

/// Knot array of one curve for `FillSubIntervals` (`cxx:448-556`).
fn curve_interval_knots(knots: Option<&[f64]>, degree: Option<usize>, first: f64, last: f64) -> Vec<f64> {
    // `if (!aBS3DCurv.IsNull()) { ... } else { anArrTemp }` (`cxx:456-506`).
    let Some(knots) = knots else {
        return vec![first, last];
    };

    // `if (aBS3DCurv->NbKnots() <= aMaxKnots)` (`cxx:460`).
    if knots.len() <= A_MAX_KNOTS {
        return knots.to_vec();
    }

    // OCCT's `isTrimmed3D` / `isTrimmed2D` (`cxx:442`) are initialised to
    // `false` and never assigned anywhere in the file, so the trimmed arms
    // (`cxx:466-481`, `cxx:517-531`) are unreachable; only the `else` arm
    // (`cxx:483-486`, `cxx:534-537`) is ported. `Degree()` is `Some` whenever
    // `Knots()` is `Some` (`crates/occt-geom/src/curve.rs:75-83`).
    let mults = knot_multiplicities(knots);
    let degree = degree.unwrap_or(0);
    let knot_count = last_uknot_index(degree, &mults) - first_uknot_index(degree, &mults) + 1;

    if knot_count <= A_MAX_KNOTS {
        knots.to_vec()
    } else {
        uniform_knots(first, last)
    }
}

/// `FillSubIntervals` (`cxx:430-640`).
///
/// Returns the number of sub-intervals `[aNbSubIntervals]`; when
/// `sub_intervals` is `Some` the `aNbSubIntervals + 1` break points are pushed
/// into it in OCCT's 1-based order (`cxx:570`, `:587`, `:605`, `:617`).
/// `nb_particles` is reset to 3 and then raised to the largest curve degree
/// (`cxx:441`, `:619-627`).
///
/// UNPORTED: the `try`/`catch` (`cxx:559-637`) that turns an array-range
/// failure into `aNbSubIntervals = 0` — the accesses below are in range by
/// construction.
fn fill_sub_intervals(
    curve3d: &dyn Curve,
    curve2d: &dyn Curve2d,
    first: f64,
    last: f64,
    nb_particles: &mut usize,
    mut sub_intervals: Option<&mut Vec<f64>>,
) -> usize {
    let an_arr_knots3d = curve_interval_knots(
        curve3d.bspline_knots(),
        curve3d.nurbs_degree(),
        first,
        last,
    );
    let an_arr_knots2d = curve_interval_knots(
        curve2d.bspline_knots2d(),
        curve2d.bspline_degree(),
        first,
        last,
    );

    //
    // `theNbParticles = 3;` (`cxx:441`).
    *nb_particles = 3;

    let mut a_nb_sub_intervals = 1usize;

    let an_ind_max3d = an_arr_knots3d.len();
    let an_ind_max2d = an_arr_knots2d.len();

    let mut an_index3d = 0usize;
    let mut an_index2d = 0usize;

    if let Some(ref mut out) = sub_intervals {
        out.push(first); // `theSubIntervals->ChangeValue(1)` (`cxx:570`).
    }

    while an_index3d < an_ind_max3d && an_index2d < an_ind_max2d {
        let a_val3d = an_arr_knots3d[an_index3d];
        let a_val2d = an_arr_knots2d[an_index2d];
        let a_delta = a_val3d - a_val2d;

        if a_delta < PCONFUSION {
            // `aVal3D <= aVal2D` (`cxx:579-595`).
            if a_val3d > first && a_val3d < last {
                a_nb_sub_intervals += 1;
                if let Some(ref mut out) = sub_intervals {
                    out.push(a_val3d);
                }
            }
            an_index3d += 1;

            if -a_delta < PCONFUSION {
                // `aVal3D == aVal2D` (`cxx:592-595`).
                an_index2d += 1;
            }
        } else {
            // `aVal2D < aVal3D` (`cxx:597-610`).
            if a_val2d > first && a_val2d < last {
                a_nb_sub_intervals += 1;
                if let Some(ref mut out) = sub_intervals {
                    out.push(a_val2d);
                }
            }
            an_index2d += 1;
        }
    }

    if let Some(ref mut out) = sub_intervals {
        out.push(last); // `ChangeValue(aNbSubIntervals + 1)` (`cxx:617`).
    }

    //
    // `theNbParticles = std::max(theNbParticles, Degree())` (`cxx:619-627`).
    if let Some(degree) = curve3d.nurbs_degree() {
        *nb_particles = (*nb_particles).max(degree);
    }
    if let Some(degree) = curve2d.bspline_degree() {
        *nb_particles = (*nb_particles).max(degree);
    }

    a_nb_sub_intervals
}

/// `PSO_Perform` (`cxx:645-694`).
fn pso_perform(
    func: &TargetFunc<'_>,
    par_inf: &MathVector,
    par_sup: &MathVector,
    epsilon: f64,
    nb_particles: usize,
    best_value: &mut f64,
    out_param: &mut MathVector,
) -> bool {
    let a_delta_param = par_sup.value(1) - par_inf.value(1);
    if a_delta_param < PCONFUSION {
        return false;
    }

    let mut a_step_par = MathVector::new(1, 1);
    a_step_par.set_value(1, epsilon * a_delta_param);

    let mut a_particles = PsoParticlesPool::new(nb_particles, 1);

    // The pool is seeded from `aNbControlPoints = 3 * theNbParticles` samples of
    // a regular grid whose last sample is `theParSup(1)` (`cxx:664-689`).
    let a_nb_control_points = 3 * nb_particles;
    let a_step = a_delta_param / (a_nb_control_points as f64 - 1.0);
    let mut a_count = 1usize;
    let mut a_prm = par_inf.value(1);
    while a_count <= a_nb_control_points {
        let mut a_val = REAL_LAST;
        if func.value_scalar(a_prm, &mut a_val) {
            // `GetWorstParticle()` (`math_PSOParticlesPool.cxx:63-66`).
            let a_particle = a_particles.worst_particle_index();
            if a_val <= a_particles.particle(a_particle).best_distance {
                let particle = a_particles.particle_mut(a_particle);
                particle.position[0] = a_prm;
                particle.best_position[0] = a_prm;
                particle.distance = a_val;
                particle.best_distance = a_val;
            }
        }
        a_count += 1;
        a_prm = if a_count == a_nb_control_points {
            par_sup.value(1)
        } else {
            a_prm + a_step
        };
    }

    // `math_PSO aPSO(&theFunction, theParInf, theParSup, aStepPar)` uses the
    // `math_PSO.hxx:74-75` defaults `theNbParticles = 32` and `theNbIter = 100`;
    // `Perform(pool, theNbParticles, ...)` overrides the particle count only.
    let mut a_pso = Pso::new(par_inf.clone(), par_sup.clone(), a_step_par, 32, 100);
    a_pso.perform_with_pool(func, &mut a_particles, nb_particles, best_value, out_param);

    true
}

/// `MinComputing` (`cxx:695-772`).
///
/// Returns the optimal value (the inverse of the square of the maximal
/// distance) and its parameter in `best_value` / `best_parameter`.
fn min_computing(
    func: &TargetFunc<'_>,
    epsilon: f64,
    nb_particles: usize,
    best_value: &mut f64,
    best_parameter: &mut f64,
) -> bool {
    let mut par_inf = MathVector::new(1, 1);
    let mut par_sup = MathVector::new(1, 1);
    let mut out_param = MathVector::new(1, 1);
    par_inf.set_value(1, func.first_parameter());
    par_sup.set_value(1, func.last_parameter());
    *best_parameter = par_inf.value(1);
    *best_value = REAL_LAST;

    if !pso_perform(
        func,
        &par_inf,
        &par_sup,
        epsilon,
        nb_particles,
        best_value,
        &mut out_param,
    ) {
        return false;
    }

    *best_parameter = out_param.value(1);

    // `math_NewtonMinimum aMinSol(theFunction)` (`cxx:734`): the default
    // constructor is `Precision::Confusion()`, 40 iterations, convexity `1.0e-6`
    // and singularity treatment on (`math_NewtonMinimum.hxx:41-45`).
    let mut a_min_sol = NewtonMinimum::with_defaults(1);
    a_min_sol.perform(func, &out_param);

    if a_min_sol.is_done() && a_min_sol.get_status() == MathStatus::Ok {
        // math_NewtonMinimum has precised the value. We take it.
        out_param.set_value(1, a_min_sol.location().value(1));
        *best_parameter = out_param.value(1);
        *best_value = a_min_sol.minimum();
    } else {
        // Use math_PSO again but on the smaller range (`cxx:742-762`).
        let a_step = epsilon * (par_sup.value(1) - par_inf.value(1));
        par_inf.set_value(1, *best_parameter - 0.5 * a_step);
        par_sup.set_value(1, *best_parameter + 0.5 * a_step);

        let mut a_value = REAL_LAST;
        if pso_perform(
            func,
            &par_inf,
            &par_sup,
            epsilon,
            nb_particles,
            &mut a_value,
            &mut out_param,
        ) {
            if a_value < *best_value {
                *best_value = a_value;
                *best_parameter = out_param.value(1);
            }
        }
    }

    true
}

/// `GeomLib_CheckCurveOnSurface` (`GeomLib_CheckCurveOnSurface.hxx:23-78`).
///
/// OCCT's `myCurve` is a `Handle(Adaptor3d_Curve)`; here it is a borrow of the
/// caller's `Adaptor3d_Curve` (`crates/occt-geom/src/curve.rs`), so the class
/// does not clone the curve it is asked to check.
pub struct CheckCurveOnSurface<'a> {
    /// `myCurve` (`cxx:292`): `Adaptor3d_Curve`; `None` is OCCT's null handle.
    curve: Option<&'a dyn Curve>,
    /// `myErrorStatus` (`cxx:289`).
    error_status: i32,
    /// `myMaxDistance` (`cxx:290`).
    max_distance: f64,
    /// `myMaxParameter` (`cxx:291`).
    max_parameter: f64,
    /// `myTolRange` (`cxx:292`).
    tol_range: f64,
    /// `myIsParallel` (`cxx:293`).
    is_parallel: bool,
}

impl Default for CheckCurveOnSurface<'_> {
    fn default() -> Self {
        Self::new()
    }
}

impl<'a> CheckCurveOnSurface<'a> {
    /// `GeomLib_CheckCurveOnSurface()` (`cxx:286-293`).
    pub fn new() -> Self {
        Self {
            curve: None,
            error_status: 0,
            max_distance: REAL_LAST,
            max_parameter: 0.0,
            tol_range: PCONFUSION,
            is_parallel: false,
        }
    }

    /// `GeomLib_CheckCurveOnSurface(theCurve, theTolRange)` (`cxx:297-308`).
    ///
    /// `theCurve` is OCCT's `GeomAdaptor_Curve(theCurve, theFirst, theLast)`
    /// (`IntTools_Tools.cxx:749-750`); the port passes the curve itself, whose
    /// `first_parameter()` / `last_parameter()` are those bounds at every OCCT
    /// call site.
    pub fn with_curve(curve: &'a dyn Curve, tol_range: f64) -> Self {
        Self {
            curve: Some(curve),
            error_status: 0,
            max_distance: REAL_LAST,
            max_parameter: 0.0,
            tol_range,
            is_parallel: false,
        }
    }

    /// `Init()` (`cxx:312-319`): initialises all members by default values.
    pub fn init(&mut self) {
        self.curve = None;
        self.error_status = 0;
        self.max_distance = REAL_LAST;
        self.max_parameter = 0.0;
        self.tol_range = PCONFUSION;
    }

    /// `Init(theCurve, theTolRange)` (`cxx:323-331`).
    pub fn init_with(&mut self, curve: &'a dyn Curve, tol_range: f64) {
        self.curve = Some(curve);
        self.error_status = 0;
        self.max_distance = REAL_LAST;
        self.max_parameter = 0.0;
        self.tol_range = tol_range;
    }

    /// `SetParallel` (`hxx:51`).
    pub fn set_parallel(&mut self, is_parallel: bool) {
        self.is_parallel = is_parallel;
    }

    /// `IsParallel` (`hxx:56`).
    pub fn is_parallel(&self) -> bool {
        self.is_parallel
    }

    /// `IsDone` (`hxx:58`): true if the max distance has been found.
    pub fn is_done(&self) -> bool {
        self.error_status == 0
    }

    /// `ErrorStatus` (`hxx:65`): 0 OK, 1 null curve, 2 invalid parametric
    /// range, 3 error in calculations.
    pub fn error_status(&self) -> i32 {
        self.error_status
    }

    /// `MaxDistance` (`hxx:68`).
    pub fn max_distance(&self) -> f64 {
        self.max_distance
    }

    /// `MaxParameter` (`hxx:71`): the parameter at which the distance is
    /// maximal.
    pub fn max_parameter(&self) -> f64 {
        self.max_parameter
    }

    /// `Perform(theCurveOnSurface)` (`cxx:335-426`).
    ///
    /// `curve_on_surface` is the `Adaptor3d_CurveOnSurface`; `curve2d` is its
    /// `GetCurve()` handle (`cxx:363`, `:381`).
    pub fn perform(&mut self, curve_on_surface: &dyn Curve, curve2d: &dyn Curve2d) {
        let Some(my_curve) = self.curve else {
            self.error_status = 1;
            return;
        };

        // `if (myCurve.IsNull() || theCurveOnSurface.IsNull())` (`cxx:338-342`).
        if (my_curve.first_parameter() - curve_on_surface.first_parameter() > self.tol_range)
            || (my_curve.last_parameter() - curve_on_surface.last_parameter() < -self.tol_range)
        {
            self.error_status = 2;
            return;
        }

        // `const double anEpsilonRange = 1.e-3;` (`cxx:351`).
        let an_epsilon_range = 1.0e-3;

        let mut a_nb_particles = 3usize;

        // Polynomial function with degree n has not more than n-1 maxima and
        // minima, so one particle per monotonicity interval means n particles.
        let a_nb_sub_intervals = fill_sub_intervals(
            my_curve,
            curve2d,
            my_curve.first_parameter(),
            my_curve.last_parameter(),
            &mut a_nb_particles,
            None,
        );

        if a_nb_sub_intervals == 0 {
            self.error_status = 3;
            return;
        }

        let mut an_intervals: Vec<f64> = Vec::with_capacity(a_nb_sub_intervals + 2);
        fill_sub_intervals(
            my_curve,
            curve2d,
            my_curve.first_parameter(),
            my_curve.last_parameter(),
            &mut a_nb_particles,
            Some(&mut an_intervals),
        );

        // `GeomLib_CheckCurveOnSurface_Local::operator()` (`cxx:228-253`).
        // `anIntervals(1) == first`, `anIntervals(aNbSubIntervals + 1) == last`.
        let mut arr_of_dist = vec![REAL_LAST; a_nb_sub_intervals];
        let mut arr_of_param = vec![0.0f64; a_nb_sub_intervals];

        for an_i in 0..a_nb_sub_intervals {
            let a_func = TargetFunc {
                curve1: my_curve,
                curve2: curve_on_surface,
                first: an_intervals[an_i],
                last: an_intervals[an_i + 1],
            };

            let mut a_min_dist = REAL_LAST;
            let mut a_par = 0.0;
            if !min_computing(
                &a_func,
                an_epsilon_range,
                a_nb_particles,
                &mut a_min_dist,
                &mut a_par,
            ) {
                arr_of_dist[an_i] = REAL_LAST;
                arr_of_param[an_i] = a_func.first_parameter();
                continue;
            }

            arr_of_dist[an_i] = a_min_dist;
            arr_of_param[an_i] = a_par;
        }

        // `OptimalValues` (`cxx:256-269`): the smallest `-distance^2`.
        let mut a_minimal_value = arr_of_dist[0];
        let mut a_parameter = arr_of_param[0];
        for i in 1..a_nb_sub_intervals {
            if arr_of_dist[i] < a_minimal_value {
                a_minimal_value = arr_of_dist[i];
                a_parameter = arr_of_param[i];
            }
        }

        self.max_distance = a_minimal_value.abs().sqrt();
        self.max_parameter = a_parameter;
    }
}
