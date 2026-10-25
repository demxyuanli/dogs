//! `Geom2dAPI_Interpolate` — 2D B-spline interpolation through points.
//!
//! Source: `ModelingAlgorithms/TKGeomAlgo/Geom2dAPI/Geom2dAPI_Interpolate.cxx`
//! (OCCT 8.0.0) and `Geom2dAPI_Interpolate.hxx`.
//!
//! Ported: `CheckPoints` (`cxx:34-45`), `CheckTangents` (`cxx:48-65`),
//! `CheckParameters` (`cxx:70-80`), `BuildParameters` (`cxx:85-112`),
//! `BuildTangents` (`cxx:158-216`), `ScaleTangents` (`cxx:225-279`), both
//! constructors (`cxx:284-353`), both `Load` overloads (`cxx:357-418`) and
//! `PerformNonPeriodic` (`cxx:616-810`).
//!
//! UNPORTED: `PerformPeriodic` (`cxx:436-612`) and `BuildPeriodicTangent`
//! (`cxx:114-154`). [`Interpolate::perform`] on a periodic instance returns
//! `Err` carrying that marker and leaves `is_done()` false. The only in-tree
//! caller, `GeomFill_Gordon::reparamCurve`, always builds the interpolator with
//! `PeriodicFlag = false` (`GeomFill_Gordon.cxx:306`), so the missing branch is
//! not on the reachable path there.
//!
//! `Standard_ConstructionError` (`cxx:300`, `:339`, `:346`, `:368`, `:388`,
//! `:408`, `:127`, `:171`) becomes `Err`.

use occt_core::bspl::banded_interp::{interpolate_contact, knot_sequence};
use occt_core::bspl::plib_eval::eval_lagrange;
use occt_core::gp::{GpPnt2d, GpVec2d};
use occt_core::precision::REAL_SMALL;

use crate::bspline_curve::Geom2dBSplineCurve;

const CONSTRUCTION_ERROR: &str = "Geom2dAPI_Interpolate: Standard_ConstructionError";

/// `Geom2dAPI_Interpolate` (`Geom2dAPI_Interpolate.hxx`).
pub struct Interpolate {
    /// `myTolerance`.
    tolerance: f64,
    /// `myPoints`.
    points: Vec<GpPnt2d>,
    /// `myParameters`.
    parameters: Vec<f64>,
    /// `myTangents`.
    tangents: Vec<GpVec2d>,
    /// `myTangentFlags`.
    tangent_flags: Vec<bool>,
    /// `myPeriodic`.
    periodic: bool,
    /// `myTangentRequest`.
    tangent_request: bool,
    /// `myCurve`.
    curve: Option<Geom2dBSplineCurve>,
    /// `myIsDone`.
    is_done: bool,
}

impl Interpolate {
    /// `Geom2dAPI_Interpolate(thePoints, thePeriodicFlag, theTolerance)`
    /// (`cxx:284-308`).
    pub fn new(points: Vec<GpPnt2d>, periodic: bool, tolerance: f64) -> Result<Self, String> {
        if !check_points(&points, tolerance) {
            return Err(CONSTRUCTION_ERROR.into());
        }
        let parameters = build_parameters(periodic, &points);
        Ok(Self::with_state(points, parameters, periodic, tolerance))
    }

    /// `Geom2dAPI_Interpolate(thePoints, theParameters, thePeriodicFlag,
    /// theTolerance)` (`cxx:312-353`).
    pub fn with_parameters(
        points: Vec<GpPnt2d>,
        parameters: Vec<f64>,
        periodic: bool,
        tolerance: f64,
    ) -> Result<Self, String> {
        if periodic && points.len() + 1 != parameters.len() {
            return Err(CONSTRUCTION_ERROR.into());
        }
        if !check_points(&points, tolerance) {
            return Err(CONSTRUCTION_ERROR.into());
        }
        if !check_parameters(&parameters) {
            return Err(CONSTRUCTION_ERROR.into());
        }
        Ok(Self::with_state(points, parameters, periodic, tolerance))
    }

    fn with_state(
        points: Vec<GpPnt2d>,
        parameters: Vec<f64>,
        periodic: bool,
        tolerance: f64,
    ) -> Self {
        let tangents = vec![GpVec2d::zero(); points.len()];
        let tangent_flags = vec![false; points.len()];
        Self {
            tolerance,
            points,
            parameters,
            tangents,
            tangent_flags,
            periodic,
            tangent_request: false,
            curve: None,
            is_done: false,
        }
    }

    /// `Load(theTangents, theTangentFlags, theScale)` (`cxx:357-393`).
    pub fn load_tangents(
        &mut self,
        tangents: &[GpVec2d],
        flags: &[bool],
        scale: bool,
    ) -> Result<(), String> {
        self.tangent_request = true;
        if tangents.len() != self.points.len() || flags.len() != self.points.len() {
            return Err(CONSTRUCTION_ERROR.into());
        }
        self.tangent_flags = flags.to_vec();
        if !check_tangents(tangents, &self.tangent_flags, self.tolerance) {
            return Err(CONSTRUCTION_ERROR.into());
        }
        self.tangents = tangents.to_vec();
        if scale {
            self.scale_tangents()?;
        }
        Ok(())
    }

    /// `Load(theInitialTangent, theFinalTangent, theScale)` (`cxx:397-418`).
    pub fn load_end_tangents(
        &mut self,
        initial: GpVec2d,
        final_tangent: GpVec2d,
        scale: bool,
    ) -> Result<(), String> {
        let n = self.points.len();
        self.tangent_request = true;
        self.tangent_flags[0] = true;
        self.tangent_flags[n - 1] = true;
        self.tangents[0] = initial;
        self.tangents[n - 1] = final_tangent;
        if !check_tangents(&self.tangents, &self.tangent_flags, self.tolerance) {
            return Err(CONSTRUCTION_ERROR.into());
        }
        if scale {
            self.scale_tangents()?;
        }
        Ok(())
    }

    /// `Perform()` (`cxx:423-432`).
    pub fn perform(&mut self) {
        if self.periodic {
            // UNPORTED: `PerformPeriodic` (`cxx:436-612`).
            self.is_done = false;
            self.curve = None;
            return;
        }
        self.perform_non_periodic();
    }

    /// `Curve()` (`cxx:816-826`); `None` where OCCT throws `StdFail_NotDone`.
    pub fn curve(&self) -> Option<&Geom2dBSplineCurve> {
        if self.is_done {
            self.curve.as_ref()
        } else {
            None
        }
    }

    /// `IsDone()` (`cxx:834-837`).
    pub fn is_done(&self) -> bool {
        self.is_done
    }

    /// `ScaleTangents` (`cxx:225-279`) over the stored tangent arrays.
    fn scale_tangents(&mut self) -> Result<(), String> {
        let n = self.points.len();
        let degree = if n == 2 {
            1
        } else if n >= 3 {
            2
        } else {
            0
        };
        let mut index = 0usize;
        for ii in 0..self.tangent_flags.len() {
            if !self.tangent_flags[ii] {
                continue;
            }
            let last = (index + degree).min(n - 1);
            let values: Vec<[f64; 2]> = self.points[index..=last]
                .iter()
                .map(|p| [p.x(), p.y()])
                .collect();
            let params = &self.parameters[index..=last];
            let eval = eval_lagrange(self.parameters[ii], 1, degree, &values, params)
                .map_err(|_| CONSTRUCTION_ERROR.to_string())?;
            let tangent = self.tangents[ii];
            let mut value = [0.0f64; 2];
            for jj in 0..2 {
                value[0] += if jj == 0 { tangent.x().abs() } else { tangent.y().abs() };
                value[1] += eval[1][jj].abs();
            }
            let ratio = value[1] / value[0];
            self.tangents[ii] = GpVec2d::new(ratio * tangent.x(), ratio * tangent.y());
            if ii != 0 {
                index += 1;
            }
            if index > n - 1 - degree {
                index = n - 1 - degree;
            }
        }
        Ok(())
    }

    /// `PerformNonPeriodic` (`cxx:616-810`). Arrays are 1-based as in OCCT.
    fn perform_non_periodic(&mut self) {
        self.is_done = false;
        self.curve = None;

        let num_points = self.points.len();
        let mut num_distinct_knots = num_points;
        let mut num_poles = num_points;
        let degree: usize;
        if num_poles == 2 && !self.tangent_request {
            degree = 1;
        } else if num_poles == 3 && !self.tangent_request {
            degree = 2;
            num_distinct_knots = 2;
        } else {
            degree = 3;
            num_poles += 2;
            if self.tangent_request {
                // `cxx:640-647`: interior flags only.
                for ii in 1..num_points.saturating_sub(1) {
                    if self.tangent_flags[ii] {
                        num_poles += 1;
                    }
                }
            }
        }

        let mut parameters = vec![0.0f64; num_poles + 1];
        let mut flatknots = vec![0.0f64; num_poles + degree + 2];
        let mut mults = vec![0i32; num_distinct_knots + 1];
        let mut knots = vec![0.0f64; num_distinct_knots + 1];
        let mut contact = vec![0i32; num_poles + 1];
        let mut poles = vec![GpPnt2d::zero(); num_poles + 1];

        for ii in 1..=degree + 1 {
            flatknots[ii] = self.parameters[0];
            flatknots[ii + num_poles] = self.parameters[num_points - 1];
        }
        // `cxx:668-671`: interior multiplicities 1, both ends `degree + 1`.
        for ii in 2..num_distinct_knots {
            mults[ii] = 1;
        }
        mults[1] = degree as i32 + 1;
        mults[num_distinct_knots] = degree as i32 + 1;

        match degree {
            1 => {
                for ii in 1..=num_poles {
                    poles[ii] = self.points[ii - 1];
                }
                self.finish_curve(&poles[1..], &self.parameters.clone(), &mults[1..], 1);
            }
            2 => {
                knots[1] = self.parameters[0];
                knots[2] = self.parameters[2];
                for ii in 1..=num_poles {
                    poles[ii] = self.points[ii - 1];
                }
                let ok = solve_interpolation(
                    degree,
                    &flatknots[1..],
                    &self.parameters,
                    &contact[1..],
                    &mut poles,
                );
                if ok.is_ok() {
                    self.finish_curve(&poles[1..], &knots[1..], &mults[1..], 2);
                }
            }
            _ => {
                if num_points >= 3 {
                    // `cxx:706-714`: endpoint tangents when not supplied.
                    if self.build_tangents().is_err() {
                        return;
                    }
                }
                contact[2] = 1;
                parameters[1] = self.parameters[0];
                parameters[2] = self.parameters[0];
                poles[1] = self.points[0];
                poles[2] = GpPnt2d::new(self.tangents[0].x(), self.tangents[0].y());
                if self.tangent_request {
                    let mut mult_index = 2usize;
                    let mut index = 3usize;
                    let mut index1 = 2usize;
                    let mut index2 = 2usize;
                    let mut index3 = degree + 2;
                    for ii in 1..num_points.saturating_sub(1) {
                        parameters[index] = self.parameters[ii];
                        poles[index] = self.points[index2 - 1];
                        flatknots[index3] = self.parameters[ii];
                        index += 1;
                        index3 += 1;
                        if self.tangent_flags[index1 - 1] {
                            mults[mult_index] += 1;
                            contact[index] = 1;
                            flatknots[index3] = self.parameters[ii];
                            parameters[index] = self.parameters[ii];
                            poles[index] =
                                GpPnt2d::new(self.tangents[ii].x(), self.tangents[ii].y());
                            index += 1;
                            index3 += 1;
                        }
                        mult_index += 1;
                        index1 += 1;
                        index2 += 1;
                    }
                } else {
                    let mut index1 = 2usize;
                    for ii in 0..num_points {
                        parameters[index1] = self.parameters[ii];
                        index1 += 1;
                    }
                    let mut index = 3usize;
                    for ii in 1..num_points.saturating_sub(1) {
                        poles[index] = self.points[ii];
                        index += 1;
                    }
                    let mut index = degree + 1;
                    for ii in 0..num_points {
                        flatknots[index] = self.parameters[ii];
                        index += 1;
                    }
                }
                poles[num_poles - 1] =
                    GpPnt2d::new(self.tangents[num_points - 1].x(), self.tangents[num_points - 1].y());
                contact[num_poles - 1] = 1;
                parameters[num_poles] = self.parameters[num_points - 1];
                parameters[num_poles - 1] = self.parameters[num_points - 1];
                poles[num_poles] = self.points[num_points - 1];
                let ok = solve_interpolation(
                    degree,
                    &flatknots[1..],
                    &parameters[1..],
                    &contact[1..],
                    &mut poles,
                );
                if ok.is_ok() {
                    let curve_knots: Vec<f64> = self.parameters.clone();
                    self.finish_curve(&poles[1..], &curve_knots, &mults[1..], 3);
                }
            }
        }
    }

    /// `BuildTangents` (`cxx:158-216`) on the stored arrays.
    fn build_tangents(&mut self) -> Result<(), String> {
        let n = self.points.len();
        if n < 3 {
            return Err(CONSTRUCTION_ERROR.into());
        }
        let degree = if n == 3 { 2 } else { 3 };
        if !self.tangent_flags[0] {
            self.tangent_flags[0] = true;
            let values: Vec<[f64; 2]> = self.points[..=degree]
                .iter()
                .map(|p| [p.x(), p.y()])
                .collect();
            let eval = eval_lagrange(
                self.parameters[0],
                1,
                degree,
                &values,
                &self.parameters[..=degree],
            )
            .map_err(|_| CONSTRUCTION_ERROR.to_string())?;
            self.tangents[0] = GpVec2d::new(eval[1][0], eval[1][1]);
        }
        if !self.tangent_flags[n - 1] {
            self.tangent_flags[n - 1] = true;
            let base = n - 1 - degree;
            let values: Vec<[f64; 2]> = self.points[base..=n - 1]
                .iter()
                .map(|p| [p.x(), p.y()])
                .collect();
            let eval = eval_lagrange(
                self.parameters[n - 1],
                1,
                degree,
                &values,
                &self.parameters[base..=n - 1],
            )
            .map_err(|_| CONSTRUCTION_ERROR.to_string())?;
            self.tangents[n - 1] = GpVec2d::new(eval[1][0], eval[1][1]);
        }
        Ok(())
    }

    /// `new Geom2d_BSplineCurve(poles, knots, mults, degree)`
    /// (`Geom2d_BSplineCurve.cxx`, `UpdateKnots`).
    fn finish_curve(
        &mut self,
        poles: &[GpPnt2d],
        knots: &[f64],
        mults: &[i32],
        degree: usize,
    ) {
        let xs: Vec<f64> = poles.iter().map(|p| p.x()).collect();
        let ys: Vec<f64> = poles.iter().map(|p| p.y()).collect();
        let flat = knot_sequence(knots, mults, degree as i32);
        match Geom2dBSplineCurve::from_flat(xs, ys, None, flat, degree, false) {
            Ok(curve) => {
                self.curve = Some(curve);
                self.is_done = true;
            }
            Err(_) => {
                self.curve = None;
                self.is_done = false;
            }
        }
    }
}

/// `BSplCLib::Interpolate(degree, flatknots, parameters, contact,
/// poles, inversion_problem)` (`Geom2dAPI_Interpolate.cxx:698`, `:797`).
fn solve_interpolation(
    degree: usize,
    flat: &[f64],
    parameters: &[f64],
    contact: &[i32],
    poles: &mut [GpPnt2d],
) -> Result<(), i32> {
    let mut flat_poles = Vec::with_capacity(poles.len() * 2);
    for p in poles.iter() {
        flat_poles.push(p.x());
        flat_poles.push(p.y());
    }
    interpolate_contact(degree, flat, parameters, contact, &mut flat_poles, 2)?;
    for (i, p) in poles.iter_mut().enumerate() {
        *p = GpPnt2d::new(flat_poles[2 * i], flat_poles[2 * i + 1]);
    }
    Ok(())
}

/// `CheckPoints` (`cxx:34-45`).
fn check_points(points: &[GpPnt2d], tolerance: f64) -> bool {
    let tol2 = tolerance * tolerance;
    for i in 0..points.len().saturating_sub(1) {
        if points[i].square_distance(&points[i + 1]) < tol2 {
            return false;
        }
    }
    true
}

/// `CheckTangents` (`cxx:48-65`).
fn check_tangents(tangents: &[GpVec2d], flags: &[bool], tolerance: f64) -> bool {
    let tol2 = tolerance * tolerance;
    for (i, flag) in flags.iter().enumerate() {
        if *flag {
            let t = tangents[i];
            if t.x() * t.x() + t.y() * t.y() < tol2 {
                return false;
            }
        }
    }
    true
}

/// `CheckParameters` (`cxx:70-80`).
fn check_parameters(parameters: &[f64]) -> bool {
    for i in 0..parameters.len().saturating_sub(1) {
        if parameters[i + 1] - parameters[i] < REAL_SMALL {
            return false;
        }
    }
    true
}

/// `BuildParameters` (`cxx:85-112`): chord-length parameters, plus the
/// wrap-around parameter on a periodic curve.
fn build_parameters(periodic: bool, points: &[GpPnt2d]) -> Vec<f64> {
    let n = points.len();
    let mut parameters = vec![0.0f64; n];
    for i in 0..n.saturating_sub(1) {
        parameters[i + 1] = parameters[i] + points[i].distance(&points[i + 1]);
    }
    if periodic {
        let wrap = points[n - 1].distance(&points[0]);
        parameters.push(parameters[n - 1] + wrap);
    }
    parameters
}
