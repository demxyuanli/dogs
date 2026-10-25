//! `GeomAPI_Interpolate` — the B-spline interpolator through points with
//! explicit parameters and optional tangents.
//!
//! Source: `ModelingAlgorithms/TKGeomAlgo/GeomAPI/GeomAPI_Interpolate.cxx`
//! (OCCT 8.0.0).
//!
//! The class is the driver `GeomFill_GordonBuilder` uses to skin the profile
//! and guide sections (`GeomFill_GordonBuilder.cxx:241-249`, `:378-381`): its
//! `Perform` builds the flat knot vector for the requested degree
//! (1/2/3 depending on the point count), fills the contact-order array, and
//! hands the whole system to `BSplCLib::Interpolate`
//! ([`occt_core::bspl::banded_interp::interpolate_contact`]).
//!
//! Translation notes:
//! - OCCT indexes every array from 1; the port keeps `Vec` storage and reads
//!   through 1-based helpers so the index arithmetic of the original stays
//!   readable. The helpers are local to [`perform`].
//! - `Standard_ConstructionError` raised by the two constructors and the two
//!   `Load` overloads (`cxx:297`, `:330-335`, `:343-346`, `:365-387`,
//!   `:402-407`) becomes `Err`, `StdFail_NotDone` from `Curve()` (`cxx:817-821`)
//!   becomes `None`.

use crate::bspline_curve::GeomBSplineCurve;
use occt_core::gp::{GpPnt, GpVec};

mod perform;
mod tangents;

/// `GeomAPI_Interpolate` (`GeomAPI_Interpolate.hxx:36-160`).
pub struct GeomAPIInterpolate {
    /// `myTolerance`.
    tolerance: f64,
    /// `myPoints` — the interpolation points, 0-based.
    points: Vec<GpPnt>,
    /// `myParameters` — one parameter per point (plus the wrap-around value on
    /// a periodic curve), 0-based.
    parameters: Vec<f64>,
    /// `myPeriodic`.
    periodic: bool,
    /// `myTangents`.
    tangents: Vec<GpVec>,
    /// `myTangentFlags`.
    tangent_flags: Vec<bool>,
    /// `myTangentRequest`.
    tangent_request: bool,
    /// `myCurve`.
    curve: Option<GeomBSplineCurve>,
    /// `myIsDone`.
    is_done: bool,
}

/// `CheckPoints` (`GeomAPI_Interpolate.cxx:30-43`): consecutive points must be
/// at least `tolerance` apart.
fn check_points(points: &[GpPnt], tolerance: f64) -> bool {
    let tolerance_squared = tolerance * tolerance;
    for i in 0..points.len().saturating_sub(1) {
        if points[i].square_distance(&points[i + 1]) < tolerance_squared {
            return false;
        }
    }
    true
}

/// `CheckTangents` (`GeomAPI_Interpolate.cxx:47-64`): a flagged tangent must
/// not be shorter than `tolerance`.
fn check_tangents(tangents: &[GpVec], flags: &[bool], tolerance: f64) -> bool {
    let tolerance_squared = tolerance * tolerance;
    for (i, flag) in flags.iter().enumerate() {
        if *flag && tangents[i].square_magnitude() < tolerance_squared {
            return false;
        }
    }
    true
}

/// `CheckParameters` (`GeomAPI_Interpolate.cxx:68-80`): the parameters must be
/// strictly increasing by at least `RealSmall()`.
fn check_parameters(parameters: &[f64]) -> bool {
    for i in 0..parameters.len().saturating_sub(1) {
        if parameters[i + 1] - parameters[i] < occt_core::precision::REAL_SMALL {
            return false;
        }
    }
    true
}

/// `BuildParameters` (`GeomAPI_Interpolate.cxx:84-110`): chord-length
/// parameters on `[0, total]`, with one extra wrap-around value for a periodic
/// curve.
fn build_parameters(periodic: bool, points: &[GpPnt]) -> Vec<f64> {
    let n = points.len();
    let num_parameters = n + usize::from(periodic);
    let mut parameters = vec![0.0f64; num_parameters];
    parameters[0] = 0.0;
    // `cxx:97-103`: `ParametersPtr->Value(ii)` is the previous parameter, i.e.
    // the slot `index - 2` in 0-based storage.
    for ii in 0..n.saturating_sub(1) {
        let distance = points[ii].distance(&points[ii + 1]);
        parameters[ii + 1] = parameters[ii] + distance;
    }
    if periodic && n > 0 {
        let distance = points[n - 1].distance(&points[0]);
        parameters[n] = parameters[n - 1] + distance;
    }
    parameters
}

impl GeomAPIInterpolate {
    /// `GeomAPI_Interpolate(Points, PeriodicFlag, Tolerance)`
    /// (`GeomAPI_Interpolate.cxx:282-307`): the parameters are the chord
    /// lengths of the point sequence.
    pub fn new(points: Vec<GpPnt>, periodic: bool, tolerance: f64) -> Result<Self, String> {
        if !check_points(&points, tolerance) {
            return Err("Standard_ConstructionError: GeomAPI_Interpolate".into());
        }
        let parameters = build_parameters(periodic, &points);
        let n = points.len();
        Ok(Self {
            tolerance,
            points,
            parameters,
            periodic,
            tangents: vec![GpVec::zero(); n],
            tangent_flags: vec![false; n],
            tangent_request: false,
            curve: None,
            is_done: false,
        })
    }

    /// `GeomAPI_Interpolate(Points, Parameters, PeriodicFlag, Tolerance)`
    /// (`GeomAPI_Interpolate.cxx:310-350`). A periodic curve needs
    /// `points + 1` parameters (`cxx:326-331`).
    pub fn with_parameters(
        points: Vec<GpPnt>,
        parameters: Vec<f64>,
        periodic: bool,
        tolerance: f64,
    ) -> Result<Self, String> {
        if !check_points(&points, tolerance) {
            return Err("Standard_ConstructionError: GeomAPI_Interpolate".into());
        }
        if periodic && points.len() + 1 != parameters.len() {
            return Err("Standard_ConstructionError: GeomAPI_Interpolate".into());
        }
        if !check_parameters(&parameters) {
            return Err("Standard_ConstructionError: GeomAPI_Interpolate".into());
        }
        let n = points.len();
        Ok(Self {
            tolerance,
            points,
            parameters,
            periodic,
            tangents: vec![GpVec::zero(); n],
            tangent_flags: vec![false; n],
            tangent_request: false,
            curve: None,
            is_done: false,
        })
    }

    /// `Load(Tangents, TangentFlags, Scale)` (`GeomAPI_Interpolate.cxx:353-388`).
    pub fn load_tangents(
        &mut self,
        tangents: Vec<GpVec>,
        flags: Vec<bool>,
        scale: bool,
    ) -> Result<(), String> {
        self.tangent_request = true;
        if tangents.len() != self.points.len() || flags.len() != self.points.len() {
            return Err("Standard_ConstructionError: GeomAPI_Interpolate::Load".into());
        }
        if !check_tangents(&tangents, &flags, self.tolerance) {
            return Err("Standard_ConstructionError: GeomAPI_Interpolate::Load".into());
        }
        self.tangents = tangents;
        self.tangent_flags = flags;
        if scale {
            self.scale_tangents();
        }
        Ok(())
    }

    /// `Load(InitialTangent, FinalTangent, Scale)`
    /// (`GeomAPI_Interpolate.cxx:394-417`).
    pub fn load_end_tangents(
        &mut self,
        initial: GpVec,
        fin: GpVec,
        scale: bool,
    ) -> Result<(), String> {
        self.tangent_request = true;
        let n = self.points.len();
        if n < 2 {
            return Err("Standard_ConstructionError: GeomAPI_Interpolate::Load".into());
        }
        self.tangent_flags[0] = true;
        self.tangent_flags[n - 1] = true;
        self.tangents[0] = initial;
        self.tangents[n - 1] = fin;
        if !check_tangents(&self.tangents, &self.tangent_flags, self.tolerance) {
            return Err("Standard_ConstructionError: GeomAPI_Interpolate::Load".into());
        }
        if scale {
            self.scale_tangents();
        }
        Ok(())
    }

    /// `Perform()` (`GeomAPI_Interpolate.cxx:420-430`).
    pub fn perform(&mut self) {
        if self.periodic {
            self.perform_periodic();
        } else {
            self.perform_non_periodic();
        }
    }

    /// `Curve()` (`GeomAPI_Interpolate.cxx:814-823`); `None` where OCCT throws
    /// `StdFail_NotDone`.
    pub fn curve(&self) -> Option<&GeomBSplineCurve> {
        if self.is_done {
            self.curve.as_ref()
        } else {
            None
        }
    }

    /// `IsDone()` (`GeomAPI_Interpolate.cxx:834-837`).
    pub fn is_done(&self) -> bool {
        self.is_done
    }

    /// `myPoints` (`GeomAPI_Interpolate.hxx:74`).
    pub fn points(&self) -> &[GpPnt] {
        &self.points
    }

    /// `myParameters` (`GeomAPI_Interpolate.hxx:78`).
    pub fn parameters(&self) -> &[f64] {
        &self.parameters
    }
}
