//! `Law_Interpolate` and `Law_Interpol` �?interpolation of a 1D evolution law.
//!
//! Source: `Law_Interpolate.hxx`/`Law_Interpolate.cxx:30-605` and
//! `Law_Interpol.hxx`/`Law_Interpol.cxx:30-181` (OCCT 8.0.0).
//!
//! Translation notes:
//! - OCCT passes `NCollection_Array1` around by pointer and fills the arrays
//!   in place (`BuildPeriodicTangent` writes `myTangents` /
//!   `myTangentFlags`). The port keeps the same scalars in a struct and takes
//!   `&mut self` where OCCT mutates the members.
//! - The work arrays inside `PerformNonPeriodic` / `PerformPeriodic` are
//!   indexed 1-based in OCCT (`parameters`, `flatknots`, `mults`,
//!   `contact_order_array`, `poles`). The port allocates them with a leading
//!   sentinel so the OCCT index `I` is the vector index `I`; the slices handed
//!   to `BSplCLib` are `[1..]`.
//! - `BSplCLib::Interpolate` answers `inversion_problem`; its failure leaves
//!   `myIsDone` false (`Law_Interpolate.cxx:390`, `:477`, `:573`). The port's
//!   `interpolate_contact` `Err` is treated the same way.
//! - `Law_Interpolate::ClearTangents` is declared but never defined in OCCT
//!   (`Law_Interpolate.hxx:70` has no `.cxx` counterpart and no caller), so it
//!   is not ported.
//! - `PLib::EvalLagrange` is called with `ArrayDimension == 1` throughout; the
//!   port routes it through `occt_core`'s 2D `eval_lagrange` with a zero second
//!   coordinate, which is the same arithmetic.

use std::sync::Arc;

use occt_core::bspl::banded_interp::interpolate_contact;
use occt_core::bspl::plib_eval::eval_lagrange;
use occt_core::kernel::geomabs::Shape;
use occt_core::precision::REAL_SMALL;

use crate::law::bspline::{LawBSpline, LawBSpFunc};
use crate::law::function::LawFunction;

/// `PLib::EvalLagrange` with `ArrayDimension == 1` (`PLib.cxx:1122-1250`):
/// `values` / `params` hold `degree + 1` entries, the answer holds
/// `min(derivative_request, degree) + 1` entries.
fn eval_lagrange_1d(
    parameter: f64,
    derivative_request: usize,
    degree: usize,
    values: &[f64],
    params: &[f64],
) -> Result<Vec<f64>, i32> {
    let pairs: Vec<[f64; 2]> = values.iter().map(|&v| [v, 0.0]).collect();
    let result = eval_lagrange(parameter, derivative_request, degree, &pairs, params)?;
    Ok(result.iter().map(|r| r[0]).collect())
}

/// `CheckParameters` (`Law_Interpolate.cxx:32-41`).
fn check_parameters(parameters: &[f64]) -> bool {
    let mut result = true;
    for ii in 0..parameters.len().saturating_sub(1) {
        if !result {
            break;
        }
        let distance = parameters[ii + 1] - parameters[ii];
        result = distance >= REAL_SMALL;
    }
    result
}

/// `BuildParameters` (`Law_Interpolate.cxx:45-70`).
fn build_parameters(periodic_flag: bool, points: &[f64]) -> Vec<f64> {
    let num_parameters = if periodic_flag {
        points.len() + 1
    } else {
        points.len()
    };
    let mut parameters = vec![0.0; num_parameters];
    // `cxx:55`: `Parameters(1) = 0.`
    if num_parameters > 0 {
        parameters[0] = 0.0;
    }
    // `cxx:58-63`: `ii` runs over `[Lower, Upper)`, `index` starts at 2.
    if !points.is_empty() {
        let mut index = 1usize;
        for ii in 0..points.len().saturating_sub(1) {
            let distance = (points[ii] - points[ii + 1]).abs();
            parameters[index] = parameters[ii] + distance;
            index += 1;
        }
        // `cxx:65-69`.
        if periodic_flag && index < num_parameters {
            let ii = points.len() - 1;
            let distance = (points[points.len() - 1] - points[0]).abs();
            parameters[index] = parameters[ii] + distance;
        }
    }
    parameters
}

/// `BuildPeriodicTangent` (`Law_Interpolate.cxx:74-110`).
fn build_periodic_tangent(
    points: &[f64],
    tangents: &mut [f64],
    tangent_flags: &mut [bool],
    parameters: &[f64],
) -> Result<(), String> {
    if points.len() < 2 {
        // `cxx:82-85`.
        tangent_flags[0] = true;
        tangents[0] = 0.0;
    } else if !tangent_flags[0] {
        // `cxx:86-108`: the closing tangent comes from a degree-2 Lagrange
        // interpolation through the last point, the closing point and the
        // second point.
        let degree = 2;
        let period = parameters[parameters.len() - 1] - parameters[0];
        let point_array = [
            points[points.len() - 1],
            points[0],
            points[1],
        ];
        let parameter_array = [
            parameters[parameters.len() - 2] - period,
            parameters[0],
            parameters[1],
        ];
        tangent_flags[0] = true;
        let eval_result = eval_lagrange_1d(
            parameter_array[1],
            1,
            degree,
            &point_array,
            &parameter_array,
        )
        .map_err(|_| "Standard_ConstructionError: Law_Interpolate".to_string())?;
        tangents[0] = eval_result[1];
    }
    Ok(())
}

/// `BuildTangents` (`Law_Interpolate.cxx:114-158`).
fn build_tangents(
    points: &[f64],
    tangents: &mut [f64],
    tangent_flags: &mut [bool],
    parameters: &[f64],
) -> Result<(), String> {
    let mut degree = 3;
    if points.len() < 3 {
        // `cxx:122-125`.
        return Err("Standard_ConstructionError: Law_Interpolate".into());
    }
    if points.len() == 3 {
        degree = 2;
    }
    if !tangent_flags[0] {
        // `cxx:130-144`: the first `degree + 1` points and parameters.
        let point_array = &points[..=degree];
        let parameter_array = &parameters[..=degree];
        tangent_flags[0] = true;
        let eval_result = eval_lagrange_1d(parameters[0], 1, degree, point_array, parameter_array)
            .map_err(|_| "Standard_ConstructionError: Law_Interpolate".to_string())?;
        tangents[0] = eval_result[1];
    }
    let last = tangent_flags.len() - 1;
    if !tangent_flags[last] {
        // `cxx:145-156`: the last `degree + 1` points and parameters.
        let start = points.len() - 1 - degree;
        let point_array = &points[start..];
        let parameter_array = &parameters[start..];
        tangent_flags[last] = true;
        let eval_result = eval_lagrange_1d(
            parameters[parameters.len() - 1],
            1,
            degree,
            point_array,
            parameter_array,
        )
        .map_err(|_| "Standard_ConstructionError: Law_Interpolate".to_string())?;
        tangents[last] = eval_result[1];
    }
    Ok(())
}

/// `Law_Interpolate` (`Law_Interpolate.hxx:34-100`).
pub struct LawInterpolate {
    tolerance: f64,
    points: Vec<f64>,
    is_done: bool,
    curve: Option<LawBSpline>,
    tangents: Vec<f64>,
    tangent_flags: Vec<bool>,
    parameters: Vec<f64>,
    periodic: bool,
    tangent_request: bool,
}

impl LawInterpolate {
    /// `Law_Interpolate(Points, PeriodicFlag, Tolerance)`
    /// (`Law_Interpolate.cxx:163-179`).
    pub fn new(points: Vec<f64>, periodic_flag: bool, tolerance: f64) -> Self {
        let num_points = points.len();
        let parameters = build_parameters(periodic_flag, &points);
        Self {
            tolerance,
            points,
            is_done: false,
            curve: None,
            tangents: vec![0.0; num_points],
            tangent_flags: vec![false; num_points],
            parameters,
            periodic: periodic_flag,
            tangent_request: false,
        }
    }

    /// `Law_Interpolate(Points, Parameters, PeriodicFlag, Tolerance)`
    /// (`Law_Interpolate.cxx:184-210`).
    pub fn new_with_parameters(
        points: Vec<f64>,
        parameters: Vec<f64>,
        periodic_flag: bool,
        tolerance: f64,
    ) -> Result<Self, String> {
        // `cxx:193-199`.
        if periodic_flag && points.len() + 1 != parameters.len() {
            return Err("Standard_ConstructionError: Law_Interpolate".into());
        }
        // `cxx:204-208`.
        if !check_parameters(&parameters) {
            return Err("Standard_ConstructionError: Law_Interpolate".into());
        }
        let num_points = points.len();
        Ok(Self {
            tolerance,
            points,
            is_done: false,
            curve: None,
            tangents: vec![0.0; num_points],
            tangent_flags: vec![false; num_points],
            parameters,
            periodic: periodic_flag,
            tangent_request: false,
        })
    }

    /// `Law_Interpolate::Load(Tangents, TangentFlags)`
    /// (`Law_Interpolate.cxx:213-230`).
    pub fn load(&mut self, tangents: &[f64], tangent_flags: &[bool]) -> Result<(), String> {
        // `cxx:219-221`: the request flag is set before the length check.
        self.tangent_request = true;
        if tangents.len() != self.points.len() || tangent_flags.len() != self.points.len() {
            return Err("Standard_ConstructionError: Law_Interpolate".into());
        }
        self.tangents = tangents.to_vec();
        self.tangent_flags = tangent_flags.to_vec();
        Ok(())
    }

    /// `Law_Interpolate::Load(InitialTangent, FinalTangent)`
    /// (`Law_Interpolate.cxx:234-242`).
    pub fn load_end_tangents(&mut self, initial_tangent: f64, final_tangent: f64) {
        self.tangent_request = true;
        if self.points.is_empty() {
            return;
        }
        let last = self.points.len() - 1;
        self.tangent_flags[0] = true;
        self.tangents[0] = initial_tangent;
        self.tangent_flags[last] = true;
        self.tangents[last] = final_tangent;
    }

    /// `Law_Interpolate::Perform` (`Law_Interpolate.cxx:246-256`).
    pub fn perform(&mut self) -> Result<(), String> {
        if self.periodic {
            self.perform_periodic()
        } else {
            self.perform_non_periodic()
        }
    }

    /// `Law_Interpolate::Curve` (`Law_Interpolate.cxx:588-595`).
    pub fn curve(&self) -> Result<&LawBSpline, String> {
        if !self.is_done {
            return Err("StdFail_NotDone: Law_Interpolate".into());
        }
        self.curve
            .as_ref()
            .ok_or_else(|| "StdFail_NotDone: Law_Interpolate".to_string())
    }

    /// `Law_Interpolate::IsDone` (`Law_Interpolate.cxx:597-601`).
    pub fn is_done(&self) -> bool {
        self.is_done
    }

    /// The tolerance the interpolation was built with.
    pub fn tolerance(&self) -> f64 {
        self.tolerance
    }

    /// `Law_Interpolate::PerformNonPeriodic`
    /// (`Law_Interpolate.cxx:399-580`).
    fn perform_non_periodic(&mut self) -> Result<(), String> {
        let num_points = self.points.len();
        let mut num_distinct_knots = num_points;
        let mut num_poles = num_points;
        let degree;
        if num_poles == 2 && !self.tangent_request {
            // `cxx:408-411`.
            degree = 1;
        } else if num_poles == 3 && !self.tangent_request {
            // `cxx:412-416`.
            degree = 2;
            num_distinct_knots = 2;
        } else {
            // `cxx:418-430`.
            degree = 3;
            num_poles += 2;
            if self.tangent_request {
                // `cxx:422-429`: `ii` runs over `[Lower + 1, Upper)`.
                for ii in 1..num_points.saturating_sub(1) {
                    if self.tangent_flags[ii] {
                        num_poles += 1;
                    }
                }
            }
        }
        let mut parameters = vec![0.0; num_poles + 1];
        let mut flatknots = vec![0.0; num_poles + degree + 2];
        let mut mults = vec![0i32; num_distinct_knots + 1];
        let mut knots = vec![0.0; num_distinct_knots + 1];
        let mut contact_order_array = vec![0i32; num_poles + 1];
        let mut poles = vec![0.0; num_poles + 1];
        let my_parameters = |ii: usize| self.parameters[ii - 1];
        let my_points = |ii: usize| self.points[ii - 1];

        // `cxx:438-441`.
        for ii in 1..=degree + 1 {
            flatknots[ii] = my_parameters(1);
            flatknots[ii + num_poles] = my_parameters(num_points);
        }
        // `cxx:442-445`: `contact_order_array` is already zero.
        // `cxx:446-450`.
        for ii in 2..num_distinct_knots {
            mults[ii] = 1;
        }
        // `cxx:451-452`.
        mults[1] = (degree + 1) as i32;
        mults[num_distinct_knots] = (degree + 1) as i32;

        match degree {
            1 => {
                // `cxx:455-465`.
                for ii in 1..=num_poles {
                    poles[ii] = my_points(ii);
                }
                let curve = LawBSpline::new(
                    poles[1..=num_poles].to_vec(),
                    self.parameters.clone(),
                    mults[1..=num_distinct_knots].to_vec(),
                    degree as i32,
                    false,
                )?;
                self.curve = Some(curve);
                self.is_done = true;
            }
            2 => {
                // `cxx:466-483`.
                knots[1] = my_parameters(1);
                knots[2] = my_parameters(num_poles);
                for ii in 1..=num_poles {
                    poles[ii] = my_points(ii);
                }
                let inversion_problem = interpolate_contact(
                    degree as usize,
                    &flatknots[1..],
                    &self.parameters,
                    &contact_order_array[1..],
                    &mut poles[1..],
                    1,
                );
                if inversion_problem.is_ok() {
                    let curve = LawBSpline::new(
                        poles[1..=num_poles].to_vec(),
                        knots[1..=num_distinct_knots].to_vec(),
                        mults[1..=num_distinct_knots].to_vec(),
                        degree as i32,
                        false,
                    )?;
                    self.curve = Some(curve);
                    self.is_done = true;
                }
            }
            _ => {
                // `cxx:484-577`.
                if num_points >= 3 {
                    build_tangents(
                        &self.points,
                        &mut self.tangents,
                        &mut self.tangent_flags,
                        &self.parameters,
                    )?;
                }
                contact_order_array[2] = 1;
                parameters[1] = my_parameters(1);
                parameters[2] = my_parameters(1);
                poles[1] = my_points(1);
                poles[2] = self.tangents[0];
                let mut mult_index = 2usize;
                let mut index = 3usize;
                let mut index3 = degree + 2;
                if self.tangent_request {
                    // `cxx:509-537`: `index1` / `index2` both start at
                    // `Lower + 1`.
                    let mut index1 = 2usize;
                    let mut index2 = 2usize;
                    for ii in 2..num_points {
                        parameters[index] = my_parameters(ii);
                        poles[index] = my_points(index2);
                        flatknots[index3] = my_parameters(ii);
                        index += 1;
                        index3 += 1;
                        if self.tangent_flags[index1 - 1] {
                            mults[mult_index] += 1;
                            contact_order_array[index] = 1;
                            flatknots[index3] = my_parameters(ii);
                            parameters[index] = my_parameters(ii);
                            poles[index] = self.tangents[ii - 1];
                            index += 1;
                            index3 += 1;
                        }
                        mult_index += 1;
                        index1 += 1;
                        index2 += 1;
                    }
                } else {
                    // `cxx:538-557`.
                    let mut index1 = 2usize;
                    for ii in 1..=num_points {
                        parameters[index1] = my_parameters(ii);
                        index1 += 1;
                    }
                    let mut idx = 3usize;
                    for ii in 2..num_points {
                        poles[idx] = my_points(ii);
                        idx += 1;
                    }
                    let mut idx = degree + 1;
                    for ii in 1..=num_points {
                        flatknots[idx] = my_parameters(ii);
                        idx += 1;
                    }
                }
                poles[num_poles - 1] = self.tangents[num_points - 1];
                contact_order_array[num_poles - 1] = 1;
                parameters[num_poles] = my_parameters(num_points);
                parameters[num_poles - 1] = my_parameters(num_points);
                poles[num_poles] = my_points(num_points);

                let inversion_problem = interpolate_contact(
                    degree as usize,
                    &flatknots[1..],
                    &parameters[1..=num_poles],
                    &contact_order_array[1..],
                    &mut poles[1..],
                    1,
                );
                if inversion_problem.is_ok() {
                    let curve = LawBSpline::new(
                        poles[1..=num_poles].to_vec(),
                        self.parameters.clone(),
                        mults[1..=num_distinct_knots].to_vec(),
                        degree as i32,
                        false,
                    )?;
                    self.curve = Some(curve);
                    self.is_done = true;
                }
            }
        }
        Ok(())
    }

    /// `Law_Interpolate::PerformPeriodic` (`Law_Interpolate.cxx:260-395`).
    fn perform_periodic(&mut self) -> Result<(), String> {
        let num_points = self.points.len();
        let period = self.parameters[self.parameters.len() - 1] - self.parameters[0];
        let num_distinct_knots = num_points + 1;
        let mut num_poles = num_points + 1;
        let half_order = 2usize;
        let degree = 3;
        num_poles += 2;
        if self.tangent_request {
            // `cxx:281-287`: `ii` runs over `[Lower + 1, Upper]`.
            for ii in 2..=num_points {
                if self.tangent_flags[ii - 1] {
                    num_poles += 1;
                }
            }
        }
        let mut parameters = vec![0.0; num_poles + 1];
        let mut flatknots = vec![0.0; num_poles + degree + 2];
        let mut mults = vec![0i32; num_distinct_knots + 1];
        let mut contact_order_array = vec![0i32; num_poles + 1];
        let mut poles = vec![0.0; num_poles + 1];
        let my_parameters = |ii: usize| self.parameters[ii - 1];
        let my_points = |ii: usize| self.points[ii - 1];
        let upper = num_points + 1;

        // `cxx:295-301`.
        for ii in 1..=half_order {
            flatknots[ii] = my_parameters(upper - 1) - period;
            flatknots[ii + half_order] = my_parameters(1);
            flatknots[num_poles + ii] = my_parameters(upper);
            flatknots[num_poles + half_order + ii] = my_parameters(half_order) + period;
        }
        // `cxx:302-306`: `contact_order_array` is already zero.
        // `cxx:307-311`.
        for ii in 2..num_distinct_knots {
            mults[ii] = 1;
        }
        mults[1] = half_order as i32;
        mults[num_distinct_knots] = half_order as i32;

        build_periodic_tangent(
            &self.points,
            &mut self.tangents,
            &mut self.tangent_flags,
            &self.parameters,
        )?;

        // `cxx:318-322`.
        contact_order_array[2] = 1;
        parameters[1] = my_parameters(1);
        parameters[2] = my_parameters(1);
        poles[1] = my_points(1);
        poles[2] = self.tangents[0];
        let mut mult_index = 2usize;
        let mut index = 3usize;
        let mut index1 = degree + 2;
        if self.tangent_request {
            // `cxx:327-348`.
            for ii in 2..=num_points {
                parameters[index] = my_parameters(ii);
                flatknots[index1] = my_parameters(ii);
                poles[index] = my_points(ii);
                index += 1;
                index1 += 1;
                if self.tangent_flags[ii - 1] {
                    mults[mult_index] += 1;
                    contact_order_array[index] = 1;
                    parameters[index] = my_parameters(ii);
                    flatknots[index1] = my_parameters(ii);
                    poles[index] = self.tangents[ii - 1];
                    index += 1;
                    index1 += 1;
                }
                mult_index += 1;
            }
        } else {
            // `cxx:349-370`.
            let mut idx = degree + 1;
            let mut idx1 = 2usize;
            for ii in 1..=upper {
                parameters[idx1] = my_parameters(ii);
                flatknots[idx] = my_parameters(ii);
                idx += 1;
                idx1 += 1;
            }
            let mut idx = 3usize;
            for ii in 2..=num_points {
                poles[idx] = my_points(ii);
                idx += 1;
            }
        }
        // `cxx:371-381`: for a periodic curve only the tangent of the first
        // point is used, because the curve closes on the first point.
        contact_order_array[num_poles - 1] = 1;
        parameters[num_poles - 1] = my_parameters(upper);
        poles[num_poles - 1] = self.tangents[0];
        parameters[num_poles] = my_parameters(upper);
        poles[num_poles] = my_points(1);

        let inversion_problem = interpolate_contact(
            degree as usize,
            &flatknots[1..],
            &parameters[1..=num_poles],
            &contact_order_array[1..],
            &mut poles[1..],
            1,
        );
        if inversion_problem.is_ok() {
            // `cxx:391-395`: the closing pole is dropped.
            let curve = LawBSpline::new(
                poles[1..=num_poles - 2].to_vec(),
                self.parameters.clone(),
                mults[1..=num_distinct_knots].to_vec(),
                degree as i32,
                self.periodic,
            )?;
            self.curve = Some(curve);
            self.is_done = true;
        }
        Ok(())
    }
}

/// `Law_Interpol` (`Law_Interpol.hxx:24-100`): a `Law_BSpFunc` built by
/// interpolation of a `(parameter, value)` table.
#[derive(Debug, Clone, Default)]
pub struct LawInterpol {
    func: LawBSpFunc,
}

impl LawInterpol {
    /// `Law_Interpol::Law_Interpol` (`Law_Interpol.cxx:32`).
    pub fn new() -> Self {
        Self::default()
    }

    /// The `(parameter, value)` table split into the `par` / `rad` arrays that
    /// `Law_Interpol.cxx:39-62` builds. For a periodic table the last value is
    /// dropped, so `rad` holds `nbp - 1` entries as OCCT does.
    fn split(par_and_rad: &[(f64, f64)], periodic: bool) -> (Vec<f64>, Vec<f64>) {
        let nbp = par_and_rad.len();
        let mut par = Vec::with_capacity(nbp);
        let mut rad = Vec::with_capacity(if periodic { nbp.saturating_sub(1) } else { nbp });
        for (i, &(x, y)) in par_and_rad.iter().enumerate() {
            par.push(x);
            if !periodic || i + 1 != nbp {
                rad.push(y);
            }
        }
        (par, rad)
    }

    /// `Law_Interpol::Set(ParAndRad, Periodic)` (`Law_Interpol.cxx:36-67`).
    pub fn set(&mut self, par_and_rad: &[(f64, f64)], periodic: bool) -> Result<(), String> {
        let (par, rad) = Self::split(par_and_rad, periodic);
        self.build(rad, par, periodic, None)
    }

    /// `Law_Interpol::Set(ParAndRad, Dd, Df, Periodic)`
    /// (`Law_Interpol.cxx:108-142`).
    pub fn set_with_derivatives(
        &mut self,
        par_and_rad: &[(f64, f64)],
        dd: f64,
        df: f64,
        periodic: bool,
    ) -> Result<(), String> {
        let (par, rad) = Self::split(par_and_rad, periodic);
        self.build(rad, par, periodic, Some((dd, df)))
    }

    /// `Law_Interpol::SetInRelative(ParAndRad, Ud, Uf, Periodic)`
    /// (`Law_Interpol.cxx:71-104`): the parameters are remapped from their own
    /// span `[wd, wf]` onto `[Ud, Uf]`.
    pub fn set_in_relative(
        &mut self,
        par_and_rad: &[(f64, f64)],
        ud: f64,
        uf: f64,
        periodic: bool,
    ) -> Result<(), String> {
        let (par, rad) = Self::split_relative(par_and_rad, ud, uf, periodic);
        self.build(rad, par, periodic, None)
    }

    /// `Law_Interpol::SetInRelative(ParAndRad, Ud, Uf, Dd, Df, Periodic)`
    /// (`Law_Interpol.cxx:146-181`).
    pub fn set_in_relative_with_derivatives(
        &mut self,
        par_and_rad: &[(f64, f64)],
        ud: f64,
        uf: f64,
        dd: f64,
        df: f64,
        periodic: bool,
    ) -> Result<(), String> {
        let (par, rad) = Self::split_relative(par_and_rad, ud, uf, periodic);
        self.build(rad, par, periodic, Some((dd, df)))
    }

    /// `Law_Interpol.cxx:75-99`: `wd = ParAndRad(Lower).X()`,
    /// `wf = ParAndRad(Upper).X()`, and every parameter becomes
    /// `(Uf * (x - wd) + Ud * (wf - x)) / (wf - wd)`.
    fn split_relative(
        par_and_rad: &[(f64, f64)],
        ud: f64,
        uf: f64,
        periodic: bool,
    ) -> (Vec<f64>, Vec<f64>) {
        let nbp = par_and_rad.len();
        let wd = par_and_rad[0].0;
        let wf = par_and_rad[nbp - 1].0;
        let mut par = Vec::with_capacity(nbp);
        let mut rad = Vec::with_capacity(if periodic { nbp.saturating_sub(1) } else { nbp });
        for (i, &(x, y)) in par_and_rad.iter().enumerate() {
            par.push((uf * (x - wd) + ud * (wf - x)) / (wf - wd));
            if !periodic || i + 1 != nbp {
                rad.push(y);
            }
        }
        (par, rad)
    }

    fn build(
        &mut self,
        rad: Vec<f64>,
        par: Vec<f64>,
        periodic: bool,
        tangents: Option<(f64, f64)>,
    ) -> Result<(), String> {
        let mut inter = LawInterpolate::new_with_parameters(
            rad,
            par,
            periodic,
            occt_core::precision::CONFUSION,
        )?;
        if let Some((dd, df)) = tangents {
            inter.load_end_tangents(dd, df);
        }
        inter.perform()?;
        let curve = inter.curve()?.clone();
        self.func.set_curve(curve);
        Ok(())
    }

    /// `Law_Interpol::Curve` (`Law_Interpol.cxx:105` via `Law_BSpFunc`).
    pub fn curve(&self) -> Option<&LawBSpline> {
        self.func.curve()
    }

    /// `Law_BSpFunc::SetCurve` as reached from `Law_Interpol`.
    pub fn set_curve(&mut self, curve: LawBSpline) {
        self.func.set_curve(curve);
    }

    /// The wrapped `Law_BSpFunc`.
    pub fn func(&self) -> &LawBSpFunc {
        &self.func
    }
}

impl LawFunction for LawInterpol {
    fn continuity(&self) -> Shape {
        self.func.continuity()
    }

    fn nb_intervals(&self, s: Shape) -> Result<i32, String> {
        self.func.nb_intervals(s)
    }

    fn intervals(&self, s: Shape) -> Result<Vec<f64>, String> {
        self.func.intervals(s)
    }

    fn value(&self, x: f64) -> f64 {
        self.func.value(x)
    }

    fn d1(&self, x: f64) -> (f64, f64) {
        self.func.d1(x)
    }

    fn d2(&self, x: f64) -> (f64, f64, f64) {
        self.func.d2(x)
    }

    fn trim(&self, p_first: f64, p_last: f64, tol: f64) -> Result<Arc<dyn LawFunction>, String> {
        self.func.trim(p_first, p_last, tol)
    }

    fn bounds(&self) -> (f64, f64) {
        self.func.bounds()
    }
}
