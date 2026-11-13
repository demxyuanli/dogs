//! Reparametrization helpers of `GeomFill_GordonUtilities.pxx` (`pxx:137-185`,
//! `:1069-1548`).

use occt_core::bspl::banded_interp::{interpolate_contact, knot_sequence};
use occt_core::bspl::function_ops::function_reparameterise;
use occt_core::bspl::knots::{nb_poles, reparametrize};
use occt_core::gp::GpPnt;
use occt_core::precision::{CONFUSION, PCONFUSION, REAL_SMALL};

use crate::bspline_curve::GeomBSplineCurve;
use crate::curve::Curve;
use crate::geom_fill::gordon::utilities::params::{
    column_means, interpolation_degree, parameter_resolution, push_constraint, row_means,
    sampling_parameters,
};

/// `ScalarLaw` (`pxx:1069-1161`): monotone piecewise-linear target -> source map.
#[derive(Clone, Debug, Default)]
pub struct ScalarLaw {
    source: Vec<f64>,
    target: Vec<f64>,
}

impl ScalarLaw {
    /// `ScalarLaw::Init` (`pxx:1074-1095`).
    pub fn init(&mut self, source: &[f64], target: &[f64]) -> bool {
        if source.len() < 2 || source.len() != target.len() {
            return false;
        }
        self.source.clear();
        self.target.clear();
        for idx in 0..source.len() {
            if !self.source.is_empty()
                && (source[idx] <= *self.source.last().unwrap() + PCONFUSION
                    || target[idx] <= *self.target.last().unwrap() + PCONFUSION)
            {
                return false;
            }
            self.source.push(source[idx]);
            self.target.push(target[idx]);
        }
        true
    }

    /// `ScalarLaw::Value` (`pxx:1097-1125`).
    pub fn value(&self, target_parameter: f64) -> f64 {
        if self.source.is_empty() {
            return 0.0;
        }
        if target_parameter <= self.target[0] {
            return self.source[0];
        }
        if target_parameter >= *self.target.last().unwrap() {
            return *self.source.last().unwrap();
        }
        for param_idx in 1..self.target.len() {
            if target_parameter <= self.target[param_idx] {
                let t0 = self.target[param_idx - 1];
                let t1 = self.target[param_idx];
                let weight = (target_parameter - t0) / (t1 - t0);
                return self.source[param_idx - 1]
                    + weight * (self.source[param_idx] - self.source[param_idx - 1]);
            }
        }
        *self.source.last().unwrap()
    }

    /// `ScalarLaw::Evaluate` derivative arm (`pxx:1127-1156`).
    pub fn evaluate_derivative(&self, target_parameter: f64) -> f64 {
        for param_idx in 1..self.target.len() {
            if target_parameter <= self.target[param_idx] {
                return (self.source[param_idx] - self.source[param_idx - 1])
                    / (self.target[param_idx] - self.target[param_idx - 1]);
            }
        }
        let n = self.source.len();
        if n < 2 {
            return 0.0;
        }
        (self.source[n - 1] - self.source[n - 2]) / (self.target[n - 1] - self.target[n - 2])
    }
}

/// `exposePeriod` (`pxx:137-167`).
pub fn expose_period(curve: &mut GeomBSplineCurve) -> bool {
    if !curve.periodic() {
        return true;
    }
    let first = curve.first_parameter();
    let last = curve.last_parameter();
    if last <= first + PCONFUSION {
        return false;
    }
    if curve.segment(first, last, PCONFUSION).is_err() {
        return false;
    }
    if curve.periodic() {
        curve.set_not_periodic();
    }
    !curve.periodic()
}

/// `makeUnitCurves` (`pxx:169-185`).
pub fn make_unit_curves(curves: &mut [GeomBSplineCurve]) -> bool {
    for curve in curves.iter_mut() {
        if !expose_period(curve) {
            return false;
        }
        let mut knots = curve.distinct_knots().to_vec();
        reparametrize(0.0, 1.0, &mut knots);
        if curve.set_knots(&knots).is_err() {
            return false;
        }
    }
    true
}

/// Local port of `GeomLib_Interpolate` (`GeomLib_Interpolate.cxx`) used by
/// `rebuildCurve` (`pxx:1448-1454`). Builds clamped flat knots and calls
/// `BSplCLib::Interpolate`.
fn geom_lib_interpolate(
    degree: i32,
    points: &[GpPnt],
    parameters: &[f64],
) -> Option<GeomBSplineCurve> {
    let num_points = points.len();
    if num_points < 2 || parameters.len() != num_points || degree < 1 {
        return None;
    }
    let num_flat = num_points + degree as usize + 1;
    let mut flat_knots = vec![0.0; num_flat];
    for ii in 0..=degree as usize {
        flat_knots[ii] = parameters[0];
    }
    for ii in 1..num_points - 1 {
        flat_knots[ii + degree as usize] = parameters[ii];
    }
    for ii in num_points..num_flat {
        flat_knots[ii] = parameters[num_points - 1];
    }
    let mut poles_buf = vec![0.0; num_points * 3];
    for (i, p) in points.iter().enumerate() {
        poles_buf[i * 3] = p.x();
        poles_buf[i * 3 + 1] = p.y();
        poles_buf[i * 3 + 2] = p.z();
    }
    let contact = vec![0i32; num_points];
    if interpolate_contact(
        degree as usize,
        &flat_knots,
        parameters,
        &contact,
        &mut poles_buf,
        3,
    )
    .is_err()
    {
        return None;
    }
    let poles: Vec<GpPnt> = (0..num_points)
        .map(|i| {
            GpPnt::new(
                poles_buf[i * 3],
                poles_buf[i * 3 + 1],
                poles_buf[i * 3 + 2],
            )
        })
        .collect();
    GeomBSplineCurve::new(poles, flat_knots, degree as usize).ok()
}

/// `reparametrizeExactly` (`pxx:1163-1294`).
pub fn reparametrize_exactly(
    curve: &GeomBSplineCurve,
    source_params: &[f64],
    target_params: &[f64],
    law: &ScalarLaw,
) -> Option<GeomBSplineCurve> {
    let nb_constraints = source_params.len();
    if nb_constraints < 2 || target_params.len() != nb_constraints {
        return None;
    }
    let degree = curve.degree() as i32;
    let mut new_knot_values = Vec::new();
    let mut new_knot_mults = Vec::new();
    let mut constraint_idx = 0usize;
    let distinct = curve.distinct_knots();
    let mults = curve.multiplicities();
    for knot_idx in 0..distinct.len() {
        let source_knot = distinct[knot_idx];
        while constraint_idx < nb_constraints
            && source_params[constraint_idx] < source_knot - PCONFUSION
        {
            new_knot_values.push(target_params[constraint_idx]);
            new_knot_mults.push(degree);
            constraint_idx += 1;
        }
        let mut target_knot = target_params[nb_constraints - 1];
        for param_idx in 0..nb_constraints - 1 {
            if source_knot <= source_params[param_idx + 1] + PCONFUSION {
                let source_first = source_params[param_idx];
                let source_last = source_params[param_idx + 1];
                let target_first = target_params[param_idx];
                let target_last = target_params[param_idx + 1];
                target_knot = target_first
                    + (source_knot - source_first) * (target_last - target_first)
                        / (source_last - source_first);
                break;
            }
        }
        let mut multiplicity = mults[knot_idx];
        if constraint_idx < nb_constraints
            && (source_params[constraint_idx] - source_knot).abs() <= PCONFUSION
        {
            if constraint_idx > 0 && constraint_idx + 1 < nb_constraints {
                multiplicity = multiplicity.max(degree);
            }
            constraint_idx += 1;
        }
        new_knot_values.push(target_knot);
        new_knot_mults.push(multiplicity);
    }
    let nb_new_poles = nb_poles(degree, false, &new_knot_mults);
    if nb_new_poles <= 0 {
        return None;
    }
    let source_flat = knot_sequence(distinct, mults, degree);
    let target_flat = knot_sequence(&new_knot_values, &new_knot_mults, degree);
    let eval = |der: i32, _se: [f64; 2], p: f64| -> Result<f64, ()> {
        if der == 0 {
            Ok(law.value(p))
        } else if der == 1 {
            Ok(law.evaluate_derivative(p))
        } else {
            Ok(0.0)
        }
    };
    if !curve.is_rational() {
        let mut poles_buf = vec![0.0; curve.nb_poles() * 3];
        for (i, p) in curve.poles.iter().enumerate() {
            poles_buf[i * 3] = p.x();
            poles_buf[i * 3 + 1] = p.y();
            poles_buf[i * 3 + 2] = p.z();
        }
        let mut new_poles_buf = vec![0.0; nb_new_poles as usize * 3];
        function_reparameterise(
            &eval,
            degree,
            &source_flat,
            &poles_buf,
            3,
            &target_flat,
            degree,
            &mut new_poles_buf,
        )
        .ok()?;
        let new_poles: Vec<GpPnt> = (0..nb_new_poles as usize)
            .map(|i| {
                GpPnt::new(
                    new_poles_buf[i * 3],
                    new_poles_buf[i * 3 + 1],
                    new_poles_buf[i * 3 + 2],
                )
            })
            .collect();
        return GeomBSplineCurve::from_poles_knots_mults(
            new_poles,
            new_knot_values,
            new_knot_mults,
            degree as usize,
        )
        .ok();
    }
    let weights = curve.weights.as_ref()?;
    let mut homogeneous = vec![0.0; 4 * curve.nb_poles()];
    for (pole_idx, pole) in curve.poles.iter().enumerate() {
        let w = weights[pole_idx];
        let base = 4 * pole_idx;
        homogeneous[base] = pole.x() * w;
        homogeneous[base + 1] = pole.y() * w;
        homogeneous[base + 2] = pole.z() * w;
        homogeneous[base + 3] = w;
    }
    let mut new_homogeneous = vec![0.0; 4 * nb_new_poles as usize];
    function_reparameterise(
        &eval,
        degree,
        &source_flat,
        &homogeneous,
        4,
        &target_flat,
        degree,
        &mut new_homogeneous,
    )
    .ok()?;
    let mut new_poles = Vec::with_capacity(nb_new_poles as usize);
    let mut new_weights = Vec::with_capacity(nb_new_poles as usize);
    for pole_idx in 0..nb_new_poles as usize {
        let base = 4 * pole_idx;
        let w = new_homogeneous[base + 3];
        if w <= REAL_SMALL {
            return None;
        }
        new_poles.push(GpPnt::new(
            new_homogeneous[base] / w,
            new_homogeneous[base + 1] / w,
            new_homogeneous[base + 2] / w,
        ));
        new_weights.push(w);
    }
    let flat = knot_sequence(&new_knot_values, &new_knot_mults, degree);
    GeomBSplineCurve::rational(new_poles, new_weights, flat, degree as usize).ok()
}

/// `rebuildCurve` (`pxx:1339-1469`).
pub fn rebuild_curve(
    curve: &mut GeomBSplineCurve,
    source_params: &[f64],
    target_params: &[f64],
    allow_approximate_rational: bool,
    did_approximate: &mut bool,
    max_deviation: &mut f64,
) -> bool {
    let mut source = Vec::new();
    let mut target = Vec::new();
    let source_tol = parameter_resolution(curve.first_parameter(), curve.last_parameter());
    let target_tol = parameter_resolution(0.0, 1.0);
    push_constraint(
        &mut source,
        &mut target,
        curve.first_parameter(),
        0.0,
        source_tol,
        target_tol,
    );
    for idx in 0..source_params.len() {
        let source_param = source_params[idx];
        let target_param = target_params[idx];
        if (source_param - curve.first_parameter()).abs() <= source_tol
            || (source_param - curve.last_parameter()).abs() <= source_tol
            || target_param.abs() <= target_tol
            || (target_param - 1.0).abs() <= target_tol
        {
            continue;
        }
        push_constraint(
            &mut source,
            &mut target,
            source_param,
            target_param,
            source_tol,
            target_tol,
        );
    }
    push_constraint(
        &mut source,
        &mut target,
        curve.last_parameter(),
        1.0,
        source_tol,
        target_tol,
    );
    let mut law = ScalarLaw::default();
    if !law.init(&source, &target) {
        return false;
    }
    if let Some(exact) = reparametrize_exactly(curve, &source, &target, &law) {
        *curve = exact;
        return true;
    }
    if curve.is_rational() {
        if !allow_approximate_rational {
            return false;
        }
        *did_approximate = true;
    }
    let samples = sampling_parameters(target_params, curve.nb_poles());
    let mut points = Vec::new();
    let mut parameters = Vec::new();
    for &target_param in &samples {
        let point = curve.d0(law.value(target_param));
        if let Some(last) = points.last() {
            if point.distance(last) <= CONFUSION {
                continue;
            }
        }
        points.push(point);
        parameters.push(target_param);
    }
    let nb_points = points.len() as i32;
    if nb_points < interpolation_degree() {
        return false;
    }
    let Some(new_curve) = geom_lib_interpolate(interpolation_degree(), &points, &parameters)
    else {
        return false;
    };
    for &target_param in &samples {
        let source_point = curve.d0(law.value(target_param));
        let new_point = new_curve.d0(target_param);
        *max_deviation = max_deviation.max(source_point.distance(&new_point));
    }
    *curve = new_curve;
    true
}

/// `equalizeParameters` (`pxx:1471-1548`).
pub fn equalize_parameters(
    profiles: &mut [GeomBSplineCurve],
    guides: &mut [GeomBSplineCurve],
    profile_params: &mut [Vec<f64>],
    guide_params: &mut [Vec<f64>],
    allow_approximate_rational: bool,
    did_approximate: &mut bool,
    max_deviation: &mut f64,
) -> bool {
    let profile_targets = column_means(profile_params);
    let guide_targets = row_means(guide_params);
    for profile_idx in 0..profiles.len() {
        let mut old_profile_params = vec![0.0; guides.len()];
        let mut is_already_target = true;
        for col in 0..guides.len() {
            old_profile_params[col] = profile_params[profile_idx][col];
            is_already_target = is_already_target
                && (old_profile_params[col] - profile_targets[col]).abs() <= PCONFUSION;
        }
        if is_already_target {
            continue;
        }
        if !rebuild_curve(
            &mut profiles[profile_idx],
            &old_profile_params,
            &profile_targets,
            allow_approximate_rational,
            did_approximate,
            max_deviation,
        ) {
            return false;
        }
        for col in 0..profile_targets.len() {
            profile_params[profile_idx][col] = profile_targets[col];
        }
    }
    for guide_idx in 0..guides.len() {
        let mut old_guide_params = vec![0.0; profiles.len()];
        let mut is_already_target = true;
        for row in 0..profiles.len() {
            old_guide_params[row] = guide_params[row][guide_idx];
            is_already_target = is_already_target
                && (old_guide_params[row] - guide_targets[row]).abs() <= PCONFUSION;
        }
        if is_already_target {
            continue;
        }
        if !rebuild_curve(
            &mut guides[guide_idx],
            &old_guide_params,
            &guide_targets,
            allow_approximate_rational,
            did_approximate,
            max_deviation,
        ) {
            return false;
        }
        for row in 0..guide_targets.len() {
            guide_params[row][guide_idx] = guide_targets[row];
        }
    }
    true
}
