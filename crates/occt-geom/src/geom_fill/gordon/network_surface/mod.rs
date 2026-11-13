//! `GeomFill_NetworkSurface` — pole-based Gordon network surface builder.
//!
//! Source: `ModelingAlgorithms/TKGeomAlgo/GeomFill/GeomFill_NetworkSurface.cxx`
//! and `.hxx` (OCCT 8.0.1).
//!
//! Builds the final surface for an already prepared Gordon network of
//! non-periodic B-spline profiles and guides. Curve families are made
//! compatible by `GeomFill_Profiler` when polynomial; rational families must
//! already share degree, poles, and knots.
//!
//! `Surface()` in OCCT throws `StdFail_NotDone` when `!IsDone()`; the port
//! returns [`None`].

mod helpers;
mod perform;

use occt_core::gp::GpPnt;

use crate::bspline_curve::GeomBSplineCurve;
use crate::bspline_surface::GeomBSplineSurface;

use helpers::{
    check_parameters, check_weights, has_periodic_curve, prepare_curve_family,
    reparametrize_surface,
};
use perform::{apply_periodicity, make_network_surface};

/// `GeomFill_NetworkSurface::ResultStatus` (`GeomFill_NetworkSurface.hxx:57-70`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ResultStatus {
    /// `NotStarted` — Perform has not been called since initialization.
    NotStarted,
    /// `Done` — surface has been constructed.
    Done,
    /// `InvalidInput` — prepared network does not satisfy builder requirements.
    InvalidInput,
    /// `CurveCompatibilityFailed` — families could not share a compatible basis.
    CurveCompatibilityFailed,
    /// `SkinningFailed` — profile or guide skin interpolation failed.
    SkinningFailed,
    /// `ReferenceSurfaceFailed` — intersection-grid reference surface failed.
    ReferenceSurfaceFailed,
    /// `KnotAlignmentFailed` — intermediate surfaces could not share knots.
    KnotAlignmentFailed,
    /// `RationalDegreeOverflow` — exact rational product degree exceeds MaxDegree.
    RationalDegreeOverflow,
    /// `RationalConstructionFailed` — exact rational numerator/denominator failed.
    RationalConstructionFailed,
    /// `ConstructionFailed` — internal B-spline construction failed.
    ConstructionFailed,
    /// `PeriodicityFailed` — closed seam could not become periodic.
    PeriodicityFailed,
}

/// `GeomFill_NetworkSurface` (`GeomFill_NetworkSurface.hxx:51-120`).
pub struct NetworkSurface {
    profiles: Vec<GeomBSplineCurve>,
    guides: Vec<GeomBSplineCurve>,
    profile_parameters: Vec<f64>,
    guide_parameters: Vec<f64>,
    /// `[guide][profile]` contact grid (`myIntersectionPoints`).
    intersection_points: Vec<Vec<GpPnt>>,
    /// `[guide][profile]` contact weights (`myIntersectionWeights`).
    intersection_weights: Vec<Vec<f64>>,
    surface: Option<GeomBSplineSurface>,
    tolerance: f64,
    is_u_closed: bool,
    is_v_closed: bool,
    status: ResultStatus,
}

impl Default for NetworkSurface {
    fn default() -> Self {
        Self::new()
    }
}

impl NetworkSurface {
    /// `GeomFill_NetworkSurface::GeomFill_NetworkSurface()` (`cxx:1171`).
    pub fn new() -> Self {
        Self {
            profiles: Vec::new(),
            guides: Vec::new(),
            profile_parameters: Vec::new(),
            guide_parameters: Vec::new(),
            intersection_points: Vec::new(),
            intersection_weights: Vec::new(),
            surface: None,
            tolerance: 0.0,
            is_u_closed: false,
            is_v_closed: false,
            status: ResultStatus::NotStarted,
        }
    }

    /// `GeomFill_NetworkSurface::Init` (`cxx:1175-1197`).
    ///
    /// `intersection_points` / `intersection_weights` are indexed
    /// `[guide_index][profile_index]` matching OCCT
    /// `Array2` row = guide, column = profile.
    #[allow(clippy::too_many_arguments)]
    pub fn init(
        &mut self,
        profiles: Vec<GeomBSplineCurve>,
        guides: Vec<GeomBSplineCurve>,
        profile_parameters: Vec<f64>,
        guide_parameters: Vec<f64>,
        intersection_points: Vec<Vec<GpPnt>>,
        intersection_weights: Vec<Vec<f64>>,
        tolerance: f64,
        is_u_closed: bool,
        is_v_closed: bool,
    ) {
        self.profiles = profiles;
        self.guides = guides;
        self.profile_parameters = profile_parameters;
        self.guide_parameters = guide_parameters;
        self.intersection_points = intersection_points;
        self.intersection_weights = intersection_weights;
        self.tolerance = tolerance;
        self.is_u_closed = is_u_closed;
        self.is_v_closed = is_v_closed;
        self.status = ResultStatus::NotStarted;
        self.surface = None;
    }

    /// `GeomFill_NetworkSurface::Perform` (`cxx:1201-1258`).
    pub fn perform(&mut self) {
        self.status = ResultStatus::NotStarted;
        self.surface = None;

        if !is_ready_to_build(
            &self.profiles,
            &self.guides,
            &self.profile_parameters,
            &self.guide_parameters,
            &self.intersection_points,
            &self.intersection_weights,
        ) {
            self.status = ResultStatus::InvalidInput;
            return;
        }

        // OCCT wraps the body in `try { OCC_CATCH_SIGNALS ... } catch
        // (Standard_Failure)` (`cxx:1217-1255`). Fallible Rust ports map to the
        // same status codes; unexpected failure is `ConstructionFailed`.
        if !prepare_curve_family(&mut self.profiles) || !prepare_curve_family(&mut self.guides) {
            self.status = ResultStatus::CurveCompatibilityFailed;
            return;
        }

        let mut surface_status = ResultStatus::ConstructionFailed;
        let surface = make_network_surface(
            &self.profiles,
            &self.guides,
            &self.profile_parameters,
            &self.guide_parameters,
            &self.intersection_points,
            &self.intersection_weights,
            self.is_u_closed,
            &mut surface_status,
        );
        let Some(mut surface) = surface else {
            self.status = surface_status;
            return;
        };

        if !apply_periodicity(
            &mut surface,
            self.is_u_closed,
            self.is_v_closed,
            self.tolerance,
        ) {
            self.surface = None;
            self.status = ResultStatus::PeriodicityFailed;
            return;
        }
        reparametrize_surface(&mut surface);
        self.surface = Some(surface);
        self.status = ResultStatus::Done;
    }

    /// `IsDone()` (`hxx:99`).
    pub fn is_done(&self) -> bool {
        self.status == ResultStatus::Done
    }

    /// `Status()` (`hxx:102`).
    pub fn status(&self) -> ResultStatus {
        self.status
    }

    /// `Surface()` (`cxx:1262-1269`); `None` where OCCT throws `StdFail_NotDone`.
    pub fn surface(&self) -> Option<&GeomBSplineSurface> {
        if self.is_done() {
            self.surface.as_ref()
        } else {
            None
        }
    }
}

/// `isReadyToBuild` (`cxx:857-873`).
fn is_ready_to_build(
    profiles: &[GeomBSplineCurve],
    guides: &[GeomBSplineCurve],
    profile_parameters: &[f64],
    guide_parameters: &[f64],
    intersection_points: &[Vec<GpPnt>],
    intersection_weights: &[Vec<f64>],
) -> bool {
    let nb_profiles = profiles.len();
    let nb_guides = guides.len();
    if nb_profiles < 2 || nb_guides < 2 {
        return false;
    }
    if has_periodic_curve(profiles) || has_periodic_curve(guides) {
        return false;
    }
    if profiles.len() != profile_parameters.len() || guides.len() != guide_parameters.len() {
        return false;
    }
    // `ColLength == guide params`, `RowLength == profile params` (`cxx:867-870`).
    if intersection_points.len() != guide_parameters.len()
        || intersection_weights.len() != guide_parameters.len()
    {
        return false;
    }
    if intersection_points
        .iter()
        .any(|row| row.len() != profile_parameters.len())
        || intersection_weights
            .iter()
            .any(|row| row.len() != profile_parameters.len())
    {
        return false;
    }
    check_parameters(profile_parameters)
        && check_parameters(guide_parameters)
        && check_weights(intersection_weights)
}
