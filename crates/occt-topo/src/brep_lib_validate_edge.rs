//! `BRepLib_ValidateEdge` (`BRepLib_ValidateEdge.cxx:1-244`,
//! `BRepLib_ValidateEdge.hxx:26-101`): the maximum distance between a 3D curve
//! and a curve on surface.
//!
//! Ported: the constructor (`cxx:24-38`), `CheckTolerance` (`cxx:42-45`),
//! `GetMaxDistance` (`cxx:49-53`), `UpdateTolerance` (`cxx:57-64`),
//! `correctTolerance` (`cxx:68-77`), `SetExitIfToleranceExceeded` (`cxx:81-85`),
//! `Process` (`cxx:89-99`), `processApprox` (`cxx:103-230`) and `processExact`
//! (`cxx:234-244`).
//!
//! UNPORTED: none.
//!
//! OCCT holds the two curves as `Handle(Adaptor3d_Curve)` /
//! `Handle(Adaptor3d_CurveOnSurface)`; here they are borrows of the port's
//! `Adaptor3d_Curve` (`occt_geom::curve::Curve`) and of
//! `Adaptor3d_CurveOnSurface` (`occt_geom::curve_on_surface::CurveOnSurface`,
//! which implements `Curve`). `myOtherCurve->GetSurface()` (`cxx:70`) is the
//! surface the curve-on-surface was loaded with; the exact arm additionally
//! needs `GetCurve()` (`cxx:238`), so both are carried separately.

use occt_core::precision::{PCONFUSION, CONFUSION};
use occt_geom::check_curve_on_surface::CheckCurveOnSurface;
use occt_geom::curve::Curve;
use occt_geom::extrema_pc::LocateExtPC;
use occt_geom::surface::Surface;
use occt_geom2d::curve::Curve2d;

use crate::brep_check::{prec_curve, prec_surface};

/// `BRepLib_ValidateEdge` (`BRepLib_ValidateEdge.hxx:26-101`).
pub struct ValidateEdge<'a> {
    /// `myReferenceCurve` (`hxx:91`).
    reference_curve: &'a dyn Curve,
    /// `myOtherCurve` (`hxx:92`), the `Adaptor3d_CurveOnSurface`.
    other_curve: &'a dyn Curve,
    /// `myOtherCurve->GetSurface()` (`cxx:70`); `None` is OCCT's null handle.
    other_surface: Option<&'a dyn Surface>,
    /// `myOtherCurve->GetCurve()` (`cxx:238`).
    other_curve2d: Option<&'a dyn Curve2d>,
    /// `mySameParameter` (`hxx:93`).
    same_parameter: bool,
    /// `myControlPointsNumber` (`hxx:94`), `22` when constructed.
    control_points_number: i32,
    /// `myToleranceForChecking` (`hxx:95`).
    tolerance_for_checking: f64,
    /// `myCalculatedDistance` (`hxx:96`).
    calculated_distance: f64,
    /// `myExitIfToleranceExceeded` (`hxx:97`).
    exit_if_tolerance_exceeded: bool,
    /// `myIsDone` (`hxx:98`).
    is_done: bool,
    /// `myIsExactMethod` (`hxx:99`).
    is_exact_method: bool,
    /// `myIsMultiThread` (`hxx:100`).
    is_multi_thread: bool,
}

impl<'a> ValidateEdge<'a> {
    /// `BRepLib_ValidateEdge(theReferenceCurve, theOtherCurve, theSameParameter)`
    /// (`cxx:24-38`). The `GetSurface()` / `GetCurve()` handles of
    /// `theOtherCurve` are passed beside it.
    pub fn new(
        reference_curve: &'a dyn Curve,
        other_curve: &'a dyn Curve,
        other_surface: Option<&'a dyn Surface>,
        other_curve2d: Option<&'a dyn Curve2d>,
        same_parameter: bool,
    ) -> Self {
        Self {
            reference_curve,
            other_curve,
            other_surface,
            other_curve2d,
            same_parameter,
            control_points_number: 22,
            tolerance_for_checking: 0.0,
            calculated_distance: 0.0,
            exit_if_tolerance_exceeded: false,
            is_done: false,
            is_exact_method: false,
            is_multi_thread: false,
        }
    }

    /// `SetExactMethod` (`hxx:39`).
    pub fn set_exact_method(&mut self, is_exact: bool) {
        self.is_exact_method = is_exact;
    }

    /// `IsExactMethod` (`hxx:42`).
    pub fn is_exact_method(&self) -> bool {
        self.is_exact_method
    }

    /// `SetParallel` (`hxx:45`).
    pub fn set_parallel(&mut self, is_multi_thread: bool) {
        self.is_multi_thread = is_multi_thread;
    }

    /// `IsParallel` (`hxx:48`).
    pub fn is_parallel(&self) -> bool {
        self.is_multi_thread
    }

    /// `SetControlPointsNumber` (`hxx:51-54`).
    pub fn set_control_points_number(&mut self, control_points_number: i32) {
        self.control_points_number = control_points_number;
    }

    /// `CheckTolerance` (`cxx:42-45`).
    pub fn check_tolerance(&self, tolerance_to_check: f64) -> bool {
        self.correct_tolerance(tolerance_to_check) > self.calculated_distance
    }

    /// `GetMaxDistance` (`cxx:49-53`).
    pub fn get_max_distance(&self) -> f64 {
        self.calculated_distance * 1.00001
    }

    /// `UpdateTolerance` (`cxx:57-64`): `theToleranceToUpdate` is grown to the
    /// corrected distance when the latter is larger.
    pub fn update_tolerance(&self, tolerance_to_update: &mut f64) {
        let a_corrected_tolerance = self.calculated_distance * 1.00001;
        if a_corrected_tolerance > *tolerance_to_update {
            *tolerance_to_update = a_corrected_tolerance;
        }
    }

    /// `correctTolerance(theTolerance)` (`cxx:68-77`).
    fn correct_tolerance(&self, tolerance: f64) -> f64 {
        // `const Handle(Adaptor3d_Surface)& aSurface = myOtherCurve->GetSurface();`
        // (`cxx:70`).
        // `double aCurvePrecision = BRepCheck::PrecCurve(*myReferenceCurve);` (`cxx:71`).
        let a_curve_precision = prec_curve(self.reference_curve);
        // `double aSurfacePrecision = BRepCheck::PrecSurface(aSurface);` (`cxx:72`).
        let a_surface_precision = match self.other_surface {
            Some(s) => prec_surface(s),
            // A null `Adaptor3d_Surface` handle is what OCCT dereferences; the
            // port has no null surface to inspect, so `PrecSurface` sees no
            // characteristic quantities and answers `RealEpsilon()`.
            None => occt_core::precision::REAL_EPSILON,
        };
        // `double aToleranceDelta = (aCurvePrecision > aSurfacePrecision)
        //  ? aCurvePrecision : aSurfacePrecision;` (`cxx:73-74`).
        let a_tolerance_delta = if a_curve_precision > a_surface_precision {
            a_curve_precision
        } else {
            a_surface_precision
        };
        // `return theTolerance + aToleranceDelta;` (`cxx:75-76`).
        tolerance + a_tolerance_delta
    }

    /// `SetExitIfToleranceExceeded(theToleranceForChecking)` (`cxx:81-85`).
    pub fn set_exit_if_tolerance_exceeded(&mut self, tolerance_for_checking: f64) {
        self.exit_if_tolerance_exceeded = true;
        self.tolerance_for_checking = self.correct_tolerance(tolerance_for_checking);
    }

    /// `Process()` (`cxx:89-99`).
    pub fn process(&mut self) {
        if self.is_exact_method && self.same_parameter {
            self.process_exact();
        } else {
            self.process_approx();
        }
    }

    /// `IsDone()` (`hxx:69`).
    pub fn is_done(&self) -> bool {
        self.is_done
    }

    /// `processApprox()` (`cxx:103-230`).
    fn process_approx(&mut self) {
        // `cxx:105-111`.
        self.is_done = true;
        let a_square_tolerance_for_checking =
            self.tolerance_for_checking * self.tolerance_for_checking;
        let a_reference_first_param = self.reference_curve.first_parameter();
        let a_reference_last_param = self.reference_curve.last_parameter();
        let an_other_first_param = self.other_curve.first_parameter();
        let an_other_last_param = self.other_curve.last_parameter();
        let mut a_max_square_distance = 0.0;

        // `cxx:113-117`.
        let a_control_points_number = if self.control_points_number < 1 {
            1
        } else {
            self.control_points_number
        };
        let an_is_projection = !self.same_parameter
            || (an_other_first_param - a_reference_first_param).abs() > PCONFUSION
            || (an_other_last_param - a_reference_last_param).abs() > PCONFUSION;

        if !an_is_projection {
            // `cxx:121-139`.
            for i in 0..=a_control_points_number {
                let a_control_point_param = ((a_control_points_number - i) as f64
                    * a_reference_first_param
                    + i as f64 * a_reference_last_param)
                    / a_control_points_number as f64;
                let a_reference_point = self.reference_curve.value(a_control_point_param);
                let an_other_point = self.other_curve.value(a_control_point_param);
                let a_square_distance = a_reference_point.square_distance(&an_other_point);
                if a_square_distance > a_max_square_distance {
                    a_max_square_distance = a_square_distance;
                }
                // `cxx:134-138`: stop for best performance.
                if self.exit_if_tolerance_exceeded
                    && a_max_square_distance > a_square_tolerance_for_checking
                {
                    self.calculated_distance = a_max_square_distance.sqrt();
                    return;
                }
            }
        } else {
            // `cxx:143-154`: the first end.
            let mut a_reference_point = self.reference_curve.value(a_reference_first_param);
            let mut an_other_point = self.other_curve.value(an_other_first_param);
            let mut a_square_distance = a_reference_point.square_distance(&an_other_point);
            if a_square_distance > a_max_square_distance {
                a_max_square_distance = a_square_distance;
            }
            if self.exit_if_tolerance_exceeded
                && a_max_square_distance > a_square_tolerance_for_checking
            {
                self.calculated_distance = a_max_square_distance.sqrt();
                return;
            }

            // `cxx:156-167`: the last end.
            a_reference_point = self.reference_curve.value(a_reference_last_param);
            an_other_point = self.other_curve.value(an_other_last_param);
            a_square_distance = a_reference_point.square_distance(&an_other_point);
            if a_square_distance > a_max_square_distance {
                a_max_square_distance = a_square_distance;
            }
            if self.exit_if_tolerance_exceeded
                && a_max_square_distance > a_square_tolerance_for_checking
            {
                self.calculated_distance = a_max_square_distance.sqrt();
                return;
            }

            // `cxx:169-177`: the two projectors, initialized over the curves'
            // own ranges with `Resolution(Precision::Confusion())`.
            let mut a_reference_extrema = LocateExtPC::default();
            a_reference_extrema.initialize(
                self.reference_curve,
                a_reference_first_param,
                a_reference_last_param,
                self.reference_curve.resolution(CONFUSION),
            );
            let mut an_other_extrema = LocateExtPC::default();
            an_other_extrema.initialize(
                self.other_curve,
                an_other_first_param,
                an_other_last_param,
                self.other_curve.resolution(CONFUSION),
            );

            // `cxx:178-227`.
            for i in 1..a_control_points_number {
                let a_reference_param = ((a_control_points_number - i) as f64
                    * a_reference_first_param
                    + i as f64 * a_reference_last_param)
                    / a_control_points_number as f64;
                let a_reference_extrema_point = self.reference_curve.value(a_reference_param);
                let an_other_param = ((a_control_points_number - i) as f64
                    * an_other_first_param
                    + i as f64 * an_other_last_param)
                    / a_control_points_number as f64;
                let an_other_extrema_point = self.other_curve.value(an_other_param);

                // `cxx:188-206`.
                a_reference_extrema.perform(&an_other_extrema_point, a_reference_param);
                if a_reference_extrema.is_done() {
                    if a_reference_extrema.square_distance() > a_max_square_distance {
                        a_max_square_distance = a_reference_extrema.square_distance();
                    }
                    if self.exit_if_tolerance_exceeded
                        && a_max_square_distance > a_square_tolerance_for_checking
                    {
                        self.calculated_distance = a_max_square_distance.sqrt();
                        return;
                    }
                } else {
                    self.is_done = false;
                    return;
                }

                // `cxx:208-226`.
                an_other_extrema.perform(&a_reference_extrema_point, an_other_param);
                if an_other_extrema.is_done() {
                    if an_other_extrema.square_distance() > a_max_square_distance {
                        a_max_square_distance = an_other_extrema.square_distance();
                    }
                    if self.exit_if_tolerance_exceeded
                        && a_max_square_distance > a_square_tolerance_for_checking
                    {
                        self.calculated_distance = a_max_square_distance.sqrt();
                        return;
                    }
                } else {
                    self.is_done = false;
                    return;
                }
            }
        }

        // `cxx:229`.
        self.calculated_distance = a_max_square_distance.sqrt();
    }

    /// `processExact()` (`cxx:234-244`).
    fn process_exact(&mut self) {
        // `GeomLib_CheckCurveOnSurface aCheckCurveOnSurface(myReferenceCurve);`
        // (`cxx:236`): the constructor's tolerance is `Precision::PConfusion()`.
        let mut a_check_curve_on_surface =
            CheckCurveOnSurface::with_curve(self.reference_curve, PCONFUSION);
        // `aCheckCurveOnSurface.SetParallel(myIsMultiThread);` (`cxx:237`).
        a_check_curve_on_surface.set_parallel(self.is_multi_thread);
        // `aCheckCurveOnSurface.Perform(myOtherCurve);` (`cxx:238`): the port's
        // `Perform` takes the curve-on-surface and its `GetCurve()` apart.
        let curve2d = match self.other_curve2d {
            Some(c) => c,
            None => return,
        };
        a_check_curve_on_surface.perform(self.other_curve, curve2d);
        // `myIsDone = aCheckCurveOnSurface.IsDone();` (`cxx:239`).
        self.is_done = a_check_curve_on_surface.is_done();
        if self.is_done {
            // `myCalculatedDistance = aCheckCurveOnSurface.MaxDistance();` (`cxx:242`).
            self.calculated_distance = a_check_curve_on_surface.max_distance();
        }
    }
}
