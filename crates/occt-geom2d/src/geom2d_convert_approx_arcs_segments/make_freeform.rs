//! `Geom2dConvert_ApproxArcsSegments::makeFreeform` and the file-scope helpers
//! it uses.
//!
//! Source: `Geom2dConvert_ApproxArcsSegments.cxx:290-542` (`makeFreeform`,
//! `getLinearParts`, `findInflection`), `:545-740` (`makeApproximation`,
//! `calculateBiArcs`, `calculateLines`) and `:826-919` (`checkContinuity`,
//! `getParameter`, both `isInflectionPoint` overloads).

use std::sync::Arc;

use occt_core::gp::{GpVec2d, GpXY};
use occt_core::precision::CONFUSION;

use crate::curve::Curve2d;
use crate::gcpnts_quasi_uniform_deflection::QuasiUniformDeflection2d;
use crate::geom2d_convert_ppoint::PPoint;
use crate::trimmed::Geom2dTrimmedCurve;

use super::{ApproxArcsSegments, Status, MAX_POINTS, MY_CURVATURE_TOLERANCE};

/// `RealLast()` (`Standard_Real.hxx:179-182`).
const REAL_LAST: f64 = f64::MAX;

/// `checkContinuity` (`cxx:826-838`): the tangents at the join are parallel
/// within `theAngleTol`.
fn check_continuity(
    curve1: Option<&Arc<dyn Curve2d>>,
    curve2: Option<&Arc<dyn Curve2d>>,
    angle_tol: f64,
) -> bool {
    let (Some(c1), Some(c2)) = (curve1, curve2) else {
        return false;
    };
    let (_, v11) = c1.d1(c1.last_parameter());
    let (_, v21) = c2.d1(c2.first_parameter());
    v11.is_parallel(&v21, angle_tol)
}

/// `getParameter` (`cxx:843-890`): the parameter of the point of the adaptor
/// curve nearest to `xy1`, by half division over `[first, last]`.
fn get_parameter(xy1: &GpXY, first_param: f64, last_param: f64, c: &dyn Curve2d) -> PPoint {
    let mut prev_param = first_param;
    let mut af1 = first_param;
    let mut af2 = last_param;
    let mut min_dist = REAL_LAST;

    for _ in 1..=MAX_POINTS {
        let at_af1 = PPoint::on_curve(af1, c);
        let dist1 = xy1.subtracted(&at_af1.point().xy()).modulus();
        if dist1 < CONFUSION {
            return at_af1;
        }

        let at_af2 = PPoint::on_curve(af2, c);
        let dist2 = xy1.subtracted(&at_af2.point().xy()).modulus();
        if dist2 < CONFUSION {
            return at_af2;
        }

        if min_dist <= dist2 - CONFUSION && min_dist <= dist1 - CONFUSION {
            break;
        }
        if dist1 < dist2 - CONFUSION {
            prev_param = af1;
            min_dist = dist1;
            af2 = (af1 + af2) * 0.5;
        } else {
            prev_param = af2;
            min_dist = dist2;
            af1 = (af1 + af2) * 0.5;
        }
    }
    PPoint::on_curve(prev_param, c)
}

/// `isInflectionPoint(theParam, theCurve)` (`cxx:894-903`): the curvature at
/// `theParam` is below `MyCurvatureTolerance`.
fn is_inflection_point(param: f64, c: &dyn Curve2d) -> bool {
    let (_, d1, d2) = c.d2(param);
    let sq_mod = d1.xy().square_modulus();
    let curvature = d1.xy().crossed(&d2.xy()).abs() / (sq_mod * sq_mod.sqrt());
    curvature < MY_CURVATURE_TOLERANCE
}

/// `isInflectionPoint(theParam, theFirstInfl, theCurve, theAnglTol)`
/// (`cxx:906-919`): the curvature test plus the angle with the chord from the
/// first inflection point.
fn is_inflection_point_rel(
    param: f64,
    first_infl: &PPoint,
    c: &dyn Curve2d,
    angle_tol: f64,
) -> bool {
    let (p1, d1, d2) = c.d2(param);
    let sq_mod = d1.xy().square_modulus();
    let curvature = d1.xy().crossed(&d2.xy()).abs() / (sq_mod * sq_mod.sqrt());
    let chord = GpVec2d::new(
        p1.x() - first_infl.point().x(),
        p1.y() - first_infl.point().y(),
    );
    let mut cont_angle = chord.angle(&d1).abs();
    cont_angle = cont_angle.min((std::f64::consts::PI - cont_angle).abs());
    curvature < MY_CURVATURE_TOLERANCE && cont_angle < angle_tol
}

impl<'a> ApproxArcsSegments<'a> {
    /// `makeFreeform` (`cxx:290-388`).
    pub(super) fn make_freeform(&mut self) -> bool {
        let mut seq_param_points: Vec<PPoint> = Vec::new();
        let mut prev_param = self.ext[0];

        self.get_linear_parts(&mut seq_param_points);
        let is_no_inf_points = seq_param_points.is_empty();

        let mut seq_linear_parts: Vec<Option<Arc<dyn Curve2d>>> = Vec::new();
        let mut is_done = true;

        // `cxx:301-315`: the linear parts between odd and even indices.
        let mut i = 1usize;
        while i < seq_param_points.len() {
            let mut param0 = seq_param_points[i - 1];
            let mut param1 = seq_param_points[i];
            let line_curve = if param0 != param1 {
                self.make_line(&mut param0, &mut param1, false)
            } else {
                None
            };
            seq_linear_parts.push(line_curve);
            i += 2;
        }

        let mut i = 1usize;
        while i < seq_param_points.len() {
            // `cxx:317-336`: the non-linear part preceding the linear part.
            if seq_param_points[i - 1] != prev_param {
                let last_ind = self.seq_curves.len();
                let mut first = prev_param;
                let mut last = seq_param_points[i - 1];
                is_done = self.make_approximation(&mut first, &mut last);
                if is_done && last_ind != 0 && self.seq_curves.len() > last_ind {
                    is_done = check_continuity(
                        self.seq_curves[last_ind - 1].as_ref(),
                        self.seq_curves[last_ind].as_ref(),
                        self.angle_tolerance,
                    );
                }
                if !is_done {
                    self.status = Status::Error;
                    break;
                }
            }

            let a_curve = seq_linear_parts[(i + 1) / 2 - 1].clone();
            if a_curve.is_some() {
                // `cxx:337-340`.
                self.seq_curves.push(a_curve);
            } else {
                // `cxx:341-360`.
                let mut param0 = seq_param_points[i - 1];
                let mut param1 = seq_param_points[i];
                let last_ind = self.seq_curves.len();
                is_done = self.make_approximation(&mut param0, &mut param1);
                if is_done && last_ind != 0 && self.seq_curves.len() > last_ind {
                    is_done = check_continuity(
                        self.seq_curves[last_ind - 1].as_ref(),
                        self.seq_curves[last_ind].as_ref(),
                        self.angle_tolerance,
                    );
                }
                if !is_done {
                    self.status = Status::Error;
                    break;
                }
            }
            prev_param = seq_param_points[i];
            i += 2;
        }

        // `cxx:366-387`: the non-linear part following the last linear part.
        if is_done && prev_param != self.ext[1] {
            let closing = self.ext[0].point().xy().subtracted(&self.ext[1].point().xy());
            if is_no_inf_points && closing.modulus() < self.tolerance {
                // A closed edge such as an ellipse is approximated in halves.
                let mid_param = 0.5 * (self.ext[0].parameter() + self.ext[1].parameter());
                let mid = PPoint::on_curve(mid_param, self.curve);
                let mut first = self.ext[0];
                let mut last = mid;
                is_done = self.make_approximation(&mut first, &mut last);
                // `myExt[0]` is passed by reference in OCCT, so its D1 follows.
                self.ext[0] = first;
                if is_done {
                    let mut first2 = mid;
                    let mut last2 = self.ext[1];
                    is_done = self.make_approximation(&mut first2, &mut last2);
                }
            } else {
                let mut first = prev_param;
                let mut last = self.ext[1];
                is_done = self.make_approximation(&mut first, &mut last);
            }
            if !is_done {
                self.status = Status::Error;
            }
        }

        !self.seq_curves.is_empty() && self.status != Status::Error
    }

    /// `getLinearParts` (`cxx:393-518`).
    fn get_linear_parts(&mut self, the_seq_par: &mut Vec<PPoint>) {
        self.seq_params.clear();

        // `cxx:399-423`: the sample points from the deflection sampler, or the
        // uniform fallback when it fails.
        let mut sampler = QuasiUniformDeflection2d::new();
        sampler.initialize_whole(self.curve, self.tolerance * 0.5, 1);
        if sampler.is_done() {
            for i in 1..=sampler.nb_points() {
                if let Some(param) = sampler.parameter(i) {
                    let p = PPoint::on_curve(param, self.curve);
                    self.seq_params.push(p);
                }
            }
        } else {
            let step = (self.ext[1].parameter() - self.ext[0].parameter()) / MAX_POINTS as f64;
            for i in 1..=MAX_POINTS {
                let param = self.ext[0].parameter() + step * i as f64;
                let p = PPoint::on_curve(param, self.curve);
                self.seq_params.push(p);
            }
        }

        // `cxx:426-448`: may the whole curve be linearised?
        let mut dir = self.ext[1].point().xy().subtracted(&self.ext[0].point().xy());
        let mod2 = dir.square_modulus();
        if mod2 > CONFUSION {
            let mut is_linear = true;
            dir = dir.divided(mod2.sqrt());
            for p in &self.seq_params {
                let d = dir.crossed(&p.point().xy().subtracted(&self.ext[0].point().xy()));
                if d > self.tolerance * 0.5 || d < -self.tolerance * 0.5 {
                    is_linear = false;
                    break;
                }
            }
            if is_linear {
                the_seq_par.push(self.ext[0]);
                the_seq_par.push(self.ext[1]);
                return;
            }
        }

        // `cxx:450-517`: the inflection-point intervals.
        let mut ind_start_linear = 0i32;
        let mut last_infl_param = self.ext[0];
        let mut first_infl_param = self.ext[0];

        for i in 1..=self.seq_params.len() {
            let cur_param = self.seq_params[i - 1];
            if ind_start_linear != 0 {
                let mut is_still_inflection_first = is_inflection_point_rel(
                    first_infl_param.parameter(),
                    &cur_param,
                    self.curve,
                    self.angle_tolerance,
                );
                if is_inflection_point_rel(
                    cur_param.parameter(),
                    &first_infl_param,
                    self.curve,
                    self.angle_tolerance,
                ) {
                    last_infl_param = self.seq_params[i - 1];
                    while !is_still_inflection_first {
                        ind_start_linear += 1;
                        if ind_start_linear >= i as i32 {
                            ind_start_linear = 0;
                            break;
                        }
                        first_infl_param = self.seq_params[(ind_start_linear - 1) as usize];
                        is_still_inflection_first = is_inflection_point_rel(
                            first_infl_param.parameter(),
                            &cur_param,
                            self.curve,
                            self.angle_tolerance,
                        );
                    }
                } else {
                    // `cxx:480-497`: the interval is only kept when it is more
                    // than ten times the tolerance.
                    last_infl_param = self.find_inflection(&last_infl_param, &cur_param);
                    if !is_inflection_point_rel(
                        first_infl_param.parameter(),
                        &last_infl_param,
                        self.curve,
                        self.angle_tolerance,
                    ) {
                        first_infl_param = self.find_inflection(&last_infl_param, &first_infl_param);
                    }
                    let dist = first_infl_param
                        .point()
                        .xy()
                        .subtracted(&last_infl_param.point().xy())
                        .modulus();
                    if first_infl_param.parameter() < last_infl_param.parameter()
                        && dist > 10.0 * self.tolerance
                    {
                        the_seq_par.push(first_infl_param);
                        the_seq_par.push(last_infl_param);
                    }
                    ind_start_linear = 0;
                }
            } else if is_inflection_point(cur_param.parameter(), self.curve) {
                last_infl_param = cur_param;
                if i > 1 {
                    first_infl_param = self.find_inflection(&cur_param, &self.seq_params[i - 2]);
                }
                ind_start_linear = i as i32;
            }
        }

        let dist = first_infl_param
            .point()
            .xy()
            .subtracted(&self.ext[1].point().xy())
            .modulus();
        if ind_start_linear != 0 && dist > 10.0 * self.tolerance {
            the_seq_par.push(first_infl_param);
            the_seq_par.push(self.ext[1]);
        }
    }

    /// `findInflection` (`cxx:521-542`): three bisections towards the
    /// inflection point.
    fn find_inflection(&self, is_infl: &PPoint, no_infl: &PPoint) -> PPoint {
        let mut lower = is_infl.parameter();
        let mut upper = no_infl.parameter();
        let mut test = 0.0;
        for _ in 0..3 {
            test = (lower + upper) * 0.5;
            if is_inflection_point_rel(test, is_infl, self.curve, self.angle_tolerance) {
                lower = test;
            } else {
                upper = test;
            }
        }
        PPoint::on_curve(test, self.curve)
    }

    /// `makeApproximation` (`cxx:545-585`).
    fn make_approximation(&mut self, first_param: &mut PPoint, last_param: &mut PPoint) -> bool {
        let mut is_done = false;
        if *last_param != *first_param {
            let distance = first_param
                .point()
                .xy()
                .subtracted(&last_param.point().xy())
                .modulus();
            if distance < self.tolerance {
                let a_curve = self.make_line(first_param, last_param, true);
                is_done = a_curve.is_some();
                if is_done && !self.seq_curves.is_empty() {
                    is_done = check_continuity(
                        self.seq_curves.last().and_then(|c| c.as_ref()),
                        a_curve.as_ref(),
                        self.angle_tolerance,
                    );
                }
                if is_done || distance < CONFUSION {
                    self.seq_curves.push(a_curve);
                    return is_done;
                }
            }

            // `cxx:569-585`: the biarc, retried on the split interval.
            is_done = self.calculate_bi_arcs(first_param, last_param);
            if !is_done {
                let mut mid = PPoint::on_curve(
                    first_param.parameter()
                        + (last_param.parameter() - first_param.parameter()) * 0.55,
                    self.curve,
                );
                is_done = self.make_approximation(first_param, &mut mid);
                if is_done {
                    let mut last = *last_param;
                    is_done = self.make_approximation(&mut mid, &mut last);
                }
            }
        }
        is_done
    }

    /// `calculateBiArcs` (`cxx:589-702`).
    fn calculate_bi_arcs(&mut self, first_param: &mut PPoint, last_param: &mut PPoint) -> bool {
        let a_resolution = occt_core::precision::REAL_SMALL;
        if first_param.d1().square_magnitude() < a_resolution
            || last_param.d1().square_magnitude() < a_resolution
        {
            return false;
        }

        let a_pnt = [first_param.point(), last_param.point()];
        let mut a_dir = [first_param.d1(), last_param.d1()];

        // `cxx:604-621`: a single arc when the two tangent angles match.
        let a_delta = a_pnt[1].xy().subtracted(&a_pnt[0].xy());
        let a_delta_vec = GpVec2d::from_xy(a_delta);
        let mut an_angle1 = a_delta_vec.angle(&a_dir[0]);
        if an_angle1 < 0.0 {
            an_angle1 = -an_angle1;
        }
        let mut an_angle2 = a_delta_vec.angle(&a_dir[1]);
        if an_angle2 < 0.0 {
            an_angle2 = -an_angle2;
        }

        let mut an_angle_threshold = occt_core::precision::ANGULAR * 10.0;
        if *first_param != self.ext[0] || *last_param != self.ext[1] {
            an_angle_threshold = self.angle_tolerance * 0.1;
        }
        if (an_angle1 - an_angle2).abs() < an_angle_threshold {
            // The last parameter is protected from the D1 update of `makeArc`
            // when the arc is rejected (`cxx:626-629`).
            let mut a_last_param = *last_param;
            let mut a_curve: Option<Geom2dTrimmedCurve> = None;
            if !self.make_arc(first_param, &mut a_last_param, true, &mut a_curve) {
                return false;
            }
            let a_curve = a_curve.unwrap();
            if self.check_curve(&a_curve, first_param.parameter(), a_last_param.parameter()) {
                *last_param = a_last_param;
                self.seq_curves.push(Some(Arc::new(a_curve)));
                return true;
            }
        }

        // `cxx:649-677`: the intersection of the two tangent lines.
        a_dir[1] = GpVec2d::from_xy(a_dir[1].xy().reversed());
        let a_prod = [
            a_pnt[0].xy().crossed(&a_dir[0].xy()),
            a_pnt[1].xy().crossed(&a_dir[1].xy()),
            a_dir[1].xy().crossed(&a_dir[0].xy()),
        ];
        let a_int_point = GpXY::new(
            (a_prod[0] * a_dir[1].x() - a_prod[1] * a_dir[0].x()) / a_prod[2],
            (a_prod[0] * a_dir[1].y() - a_prod[1] * a_dir[0].y()) / a_prod[2],
        );
        let a_diff = [
            a_int_point.subtracted(&a_pnt[0].xy()),
            a_int_point.subtracted(&a_pnt[1].xy()),
        ];
        if a_diff[0].dot(&a_dir[0].xy()) < 0.0 || a_diff[1].dot(&a_dir[1].xy()) < 0.0 {
            return false;
        }

        let ad1 = a_diff[0].modulus();
        let ad2 = a_diff[1].modulus();
        let ad12 = a_delta.modulus();
        if ad1 + ad2 == 0.0 {
            return false;
        }
        let a_b1 = ad1 / (ad1 + ad2);
        if (a_b1 - 0.5).abs() < 0.0001 {
            return false;
        }

        let a_xy = [
            a_pnt[0]
                .xy()
                .added(&a_dir[0].xy().multiplied_scalar(ad12 * ad1 / (ad12 + ad1 + ad2))),
            a_pnt[1]
                .xy()
                .added(&a_dir[1].xy().multiplied_scalar(ad12 * ad2 / (ad12 + ad1 + ad2))),
        ];
        let a_xy_mid_arc = a_xy[0].added(&a_xy[1].subtracted(&a_xy[0]).multiplied_scalar(a_b1));
        let mut a_param_mid_arc = get_parameter(
            &a_xy_mid_arc,
            first_param.parameter(),
            last_param.parameter(),
            self.curve,
        );

        // `cxx:680-691`: the first arc of the biarc.
        let mut a_curve1: Option<Geom2dTrimmedCurve> = None;
        if !self.make_arc(first_param, &mut a_param_mid_arc, true, &mut a_curve1) {
            return false;
        }
        let a_curve1 = a_curve1.unwrap();
        if !self.check_curve(&a_curve1, first_param.parameter(), a_param_mid_arc.parameter()) {
            return false;
        }

        // `cxx:692-700`: the second arc of the biarc.
        let mut a_curve2: Option<Geom2dTrimmedCurve> = None;
        if self.make_arc(last_param, &mut a_param_mid_arc, false, &mut a_curve2) {
            let a_curve2 = a_curve2.unwrap();
            if self.check_curve(&a_curve2, a_param_mid_arc.parameter(), last_param.parameter()) {
                self.seq_curves.push(Some(Arc::new(a_curve1)));
                self.seq_curves.push(Some(Arc::new(a_curve2)));
                return true;
            }
        }
        false
    }

    /// `calculateLines` (`cxx:706-740`). Dead in OCCT (no caller, the two
    /// `makeFreeform` loops build their lines directly), kept for parity with
    /// the class.
    #[allow(dead_code)]
    fn calculate_lines(&mut self, first_param: &mut PPoint, last_param: &mut PPoint) -> bool {
        let mut prev_param = *first_param;
        for i in 1..=self.seq_params.len() {
            let mut cur_param = self.seq_params[i - 1];
            if cur_param.parameter() < prev_param.parameter() {
                continue;
            }
            if cur_param.parameter() > last_param.parameter() {
                break;
            }

            if cur_param != prev_param {
                let distance = cur_param
                    .point()
                    .xy()
                    .subtracted(&prev_param.point().xy())
                    .modulus();
                if distance > self.tolerance {
                    let a_curve = self.make_line(&mut prev_param, &mut cur_param, false);
                    if a_curve.is_none() {
                        return false;
                    }
                    self.seq_curves.push(a_curve);
                    prev_param = self.seq_params[i - 1];
                }
            }
        }
        true
    }
}
