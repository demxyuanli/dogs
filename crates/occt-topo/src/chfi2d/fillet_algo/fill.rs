//! `ChFi2d_FilletAlgo`'s iterative core: `FillPoint`, `FillDiff`,
//! `PerformNewton` and `ProcessPoint`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi2d/ChFi2d_FilletAlgo.cxx:257-563`
//! (OCCT 8.0.0).
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | `FillPoint` (`cxx:257-333`) | 77 | [`FilletAlgo::fill_point`] |
//! | `FillDiff` (`cxx:335-348`) | 14 | [`FilletAlgo::fill_diff`] |
//! | `ProcessPoint` (`cxx:415-461`) | 47 | [`FilletAlgo::process_point`] |
//! | `PerformNewton` (`cxx:463-563`) | 101 | [`FilletAlgo::perform_newton`] |
//!
//! `IsRadiusIntersected` (`cxx:162-254`) lives in [`super::radius_intersected`].
//!
//! Translation notes:
//! - `FillPoint` reads `aPoint` / `aCenter` after its loop even when the loop
//!   broke on the limit test before ever assigning them. In C++ that is an
//!   uninitialised read; the port starts both at `(0, 0)` and the divergence is
//!   recorded in `specs/_tasks.md` (T-354),
//! - `aPerp.Normalize()` raises `gp_VectorWithNullMagnitude` in OCCT when
//!   `|aVec| <= gp::Resolution()`; the preceding
//!   `aVec.SquareMagnitude() < Precision::Confusion()` guard makes that
//!   unreachable for any real curve, so the port keeps the unnormalised vector,
//! - the recursion guard `myDegreeOfRecursion` is a member, shared with every
//!   nested `ProcessPoint` call, exactly as in OCCT.

use occt_core::gp::{GpPnt2d, GpVec2d};
use occt_core::precision::CONFUSION;
use occt_geom2d::geom2d_api_project_point_on_curve::ProjectPointOnCurve2d;

use super::fillet_point::FilletPoint;
use super::radius_intersected::is_radius_intersected;
use super::FilletAlgo;

impl FilletAlgo {
    /// `FillPoint(thePoint, theLimit)` (`cxx:257-333`): walks the first curve
    /// away from `thePoint`'s parameter in doubling steps until the radius
    /// vector leaves the second curve, then fills the point's function values.
    pub(super) fn fill_point(&self, the_point: &mut FilletPoint, the_limit: f64) {
        let mut a_step = CONFUSION;
        // OCCT leaves both uninitialised when the loop breaks on the first
        // limit test; see the module note.
        let mut a_center = GpPnt2d::zero();
        let mut a_point = GpPnt2d::zero();
        let mut a_param = the_point.param();
        if the_limit < a_param {
            a_step = -a_step;
        }
        let curve1 = self.curve1.as_ref().expect("FillPoint: curve1");
        let curve2 = self.curve2.as_ref().expect("FillPoint: curve2");
        let mut a_valid = false;
        while !a_valid {
            if (a_param - a_step - the_limit) * (a_param - the_limit) <= 0.0 {
                break; // limit was exceeded
            }
            a_step *= 2.0;
            let (p, a_vec) = curve1.d1(a_param);
            a_point = p;
            if a_vec.square_magnitude() < CONFUSION {
                a_param += a_step;
                continue;
            }

            let perp = GpVec2d::new(
                if self.start_side { -1.0 } else { 1.0 } * a_vec.y(),
                if self.start_side { 1.0 } else { -1.0 } * a_vec.x(),
            );
            // `aPerp.Normalize(); aPerp.Multiply(myRadius);`
            let perp = perp.normalized().unwrap_or(perp);
            let perp = perp.multiplied_scalar(self.radius);
            a_center = a_point.translated_vec(&perp);

            // `Geom2dAPI_ProjectPointOnCurve aProjInt(aPoint, myCurve2,
            // myStart2, myEnd2);` (`cxx:288-293`).
            let proj = ProjectPointOnCurve2d::with_curve_range(
                &a_point,
                curve2,
                self.start2,
                self.end2,
            );
            match proj.nearest_point() {
                None => {
                    a_valid = true;
                    break;
                }
                Some(q) if a_point.distance(&q) > CONFUSION => {
                    a_valid = true;
                    break;
                }
                _ => {}
            }
            a_param += a_step;
        }
        if a_valid {
            the_point.set_param(a_param);
            the_point.set_center(a_center);
            a_valid = !is_radius_intersected(
                curve2,
                self.start2,
                self.end2,
                &a_point,
                &a_center,
                true,
            );
        }

        // `Geom2dAPI_ProjectPointOnCurve aProj(aCenter, myCurve2);`
        // (`cxx:305`): unbounded, i.e. the curve's own range.
        let a_proj = ProjectPointOnCurve2d::with_curve(&a_center, curve2);
        let a_nb = a_proj.nb_points();
        let mut a = a_nb;
        while a > 0 {
            let Some(a_proj_point) = a_proj.point(a) else {
                a -= 1;
                continue;
            };
            if a_point.square_distance(&a_proj_point) < CONFUSION {
                a -= 1;
                continue;
            }

            let mut a_valid2 = a_valid;
            if a_valid2 {
                a_valid2 = !is_radius_intersected(
                    curve1,
                    self.start1,
                    self.end1,
                    &a_center,
                    &a_proj_point,
                    false,
                );
            }

            // checking the right parameter (`cxx:319-324`).
            let Some(mut a_param_proj) = a_proj.parameter(a) else {
                a -= 1;
                continue;
            };
            while curve2.is_periodic() && a_param_proj < self.start2 {
                a_param_proj += curve2.period();
            }

            let Some(d) = a_proj.distance(a) else {
                a -= 1;
                continue;
            };
            the_point.append_value(
                d * d - self.radius * self.radius,
                a_param_proj >= self.start2 && a_param_proj <= self.end2 && a_valid2,
            );
            if (d - self.radius).abs() < CONFUSION {
                the_point.set_param2(a_param_proj);
            }
            a -= 1;
        }
    }

    /// `FillDiff(thePoint, theDiffStep, theFront)` (`cxx:335-348`): the
    /// derivative by a symmetric difference with a step in the front or
    /// reversed direction, retried on the opposite side when the two value
    /// sequences do not match up.
    pub(super) fn fill_diff(&self, the_point: &mut FilletPoint, the_diff_step: f64, front: bool) {
        let a_delta = if front { the_diff_step } else { -the_diff_step };
        let mut a_diff = FilletPoint::new(the_point.param() + a_delta);
        self.fill_point(&mut a_diff, a_delta * 999.0);
        if !the_point.calculate_diff(&a_diff) {
            // `aDiff->setParam(...)` only moves the parameter; the values the
            // first `FillPoint` appended stay in place.
            a_diff.set_param(the_point.param() - a_delta);
            self.fill_point(&mut a_diff, -a_delta * 999.0);
            the_point.calculate_diff(&a_diff);
        }
    }

    /// `ProcessPoint(theLeft, theRight, theParameter)` (`cxx:415-461`): inserts
    /// one more point between the left and the right one and re-runs Newton on
    /// both halves. `false` is `theParameter` outside `[left, right)`.
    pub(super) fn process_point(
        &mut self,
        the_left: &FilletPoint,
        the_right: &FilletPoint,
        the_parameter: f64,
    ) -> bool {
        if !(the_parameter >= the_left.param() && the_parameter < the_right.param()) {
            return false;
        }
        let mut the_parameter = the_parameter;
        let a_dx = the_right.param() - the_left.param();
        if the_parameter - the_left.param() < a_dx / 100.0 {
            the_parameter = the_left.param() + a_dx / 100.0;
        }
        if the_right.param() - the_parameter < a_dx / 100.0 {
            the_parameter = the_right.param() - a_dx / 100.0;
        }

        // Protection on infinite loops (`cxx:431-445`).
        self.degree_of_recursion += 1;
        let mut diffx = 0.001 * a_dx;
        if self.degree_of_recursion > 100 {
            diffx *= 10.0;
            if self.degree_of_recursion > 1000 {
                diffx *= 10.0;
                if self.degree_of_recursion > 3000 {
                    return true;
                }
            }
        }

        let mut a_point1 = the_left.copy();
        let mut a_point2 = FilletPoint::new(the_parameter);
        self.fill_point(&mut a_point2, a_point1.param());
        self.fill_diff(&mut a_point2, diffx, true);

        a_point1.filter_points(&a_point2);
        self.perform_newton(&a_point1, &a_point2);
        a_point2.filter_points(the_right);
        self.perform_newton(&a_point2, the_right);
        true
    }

    /// `PerformNewton(theLeft, theRight)` (`cxx:463-563`): stores the left point
    /// when it is already a root of the fillet function, otherwise divides the
    /// interval at the root of the local quadratic model.
    pub(super) fn perform_newton(&mut self, the_left: &FilletPoint, the_right: &FilletPoint) {
        // Check the left: if this is a solution, store it and remove it from the
        // list of researching points of the left.
        if let Some(a) = the_left.has_solution(self.radius) {
            if the_left.is_valid(a) {
                self.result_params.push(the_left.param());
                self.result_orientation.push(self.start_side);
            }
            return;
        }

        let a_dx = the_right.param() - the_left.param();
        if a_dx < 1.0e-6 * CONFUSION {
            if let Some(a) = the_right.has_solution(self.radius) {
                if the_right.is_valid(a) {
                    self.result_params.push(the_right.param());
                    self.result_orientation.push(self.start_side);
                }
            }
            return;
        }

        for a in 1..=the_left.nb_values() {
            let a_near = the_left.near(a);
            let a_a = (the_right.diff(a_near) - the_left.diff(a)) / a_dx;
            let a_b = the_left.diff(a) - a_a * the_left.param();
            let a_c = the_left.value(a) - the_left.diff(a) * the_left.param()
                + a_a * the_left.param() * the_left.param() / 2.0;
            let mut a_det = a_b * a_b - 2.0 * a_a * a_c;

            if a_a.abs() < CONFUSION {
                // linear case (`cxx:502-519`).
                if a_b.abs() > 10e-20 {
                    let a_x0 = -a_c / a_b; // use extremum
                    if a_x0 > the_left.param() && a_x0 < the_right.param() {
                        self.process_point(the_left, the_right, a_x0);
                    }
                } else {
                    // linear division otherwise.
                    let mid = the_left.param() + a_dx / 2.0;
                    self.process_point(the_left, the_right, mid);
                }
            } else if a_b.abs() > (a_det * 1000000.).abs() {
                // possible floating point operations accuracy errors
                // (`cxx:524-529`).
                let mid = the_left.param() + a_dx / 2.0;
                self.process_point(the_left, the_right, mid);
            } else if a_det > 0.0 {
                // two solutions (`cxx:531-547`).
                a_det = a_det.sqrt();
                let mut a_res = self.process_point(the_left, the_right, (-a_b + a_det) / a_a);
                if !a_res {
                    a_res = self.process_point(the_left, the_right, (-a_b - a_det) / a_a);
                }
                if !a_res {
                    let mid = the_left.param() + a_dx / 2.0;
                    self.process_point(the_left, the_right, mid);
                }
            } else {
                // use extremum (`cxx:548-561`).
                let a_x0 = -a_b / a_a;
                if a_x0 > the_left.param() && a_x0 < the_right.param() {
                    self.process_point(the_left, the_right, a_x0);
                } else {
                    let mid = the_left.param() + a_dx / 2.0;
                    self.process_point(the_left, the_right, mid);
                }
            }
        }
    }
}
