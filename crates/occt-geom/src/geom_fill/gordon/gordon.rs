//! `GeomFill_Gordon` — high-level Gordon surface driver (OCCT 8.0.1).
//!
//! Source: `ModelingAlgorithms/TKGeomAlgo/GeomFill/GeomFill_Gordon.cxx`
//! and `GeomFill_Gordon.hxx`.
//!
//! `Perform` prepares the network through [`NetworkPreparation`], then builds
//! the surface with [`NetworkSurface`]. Approximate fallback uses
//! `GeomAPI_PointsToBSplineSurface` in OCCT (`cxx:148-158`); that API is not
//! ported here, so the Approx path returns `None` / `ApproximationFailed`.

use occt_core::convert::ParameterisationType;
use occt_core::gp::GpPnt;
use occt_core::precision::{CONFUSION, PCONFUSION};

use crate::bspline_curve::{GeomBSplineCurve, MAX_DEGREE};
use crate::bspline_surface::GeomBSplineSurface;
use crate::convert_bspl::curve_to_bspline_curve;
use crate::curve::Curve;
use crate::geom_fill::gordon::network_surface::{
    NetworkSurface, ResultStatus as NetworkResultStatus,
};
use crate::geom_fill::gordon::utilities::{column_means, row_means, NetworkPreparation};
use crate::surface::Surface;

/// `GeomFill_Gordon::ResultStatus` (`GeomFill_Gordon.hxx:67-88`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ResultStatus {
    NotStarted,
    Done,
    InvalidInput,
    ConversionFailed,
    IntersectionFailed,
    OrderingFailed,
    ReparametrizationFailed,
    CompatibilityFailed,
    CurveCompatibilityFailed,
    RationalReparametrizationFailed,
    SkinningFailed,
    ReferenceSurfaceFailed,
    KnotAlignmentFailed,
    RationalDegreeOverflow,
    RationalConstructionFailed,
    PeriodicityFailed,
    ApproximationFailed,
    ConstructionFailed,
}

/// `GeomFill_Gordon::ApproximationMode` (`GeomFill_Gordon.hxx:91-95`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ApproximationMode {
    ExactOnly,
    AllowApproximateFallback,
}

/// `GeomFill_Gordon::BuildStage` (`GeomFill_Gordon.hxx:98-108`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BuildStage {
    NotStarted,
    InputConversion,
    ContactDiscovery,
    NetworkOrdering,
    Reparametrization,
    ExactConstruction,
    Validation,
    Approximation,
}

/// `GeomFill_Gordon::BuildReport` (`GeomFill_Gordon.hxx:111-121`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BuildReport {
    pub status: ResultStatus,
    pub failed_stage: BuildStage,
    pub is_approximate: bool,
    pub max_contact_gap: f64,
    pub max_reparametrization_deviation: f64,
    pub max_profile_deviation: f64,
    pub max_guide_deviation: f64,
    pub max_approximation_deviation: f64,
}

impl Default for BuildReport {
    fn default() -> Self {
        Self {
            status: ResultStatus::NotStarted,
            failed_stage: BuildStage::NotStarted,
            is_approximate: false,
            max_contact_gap: 0.0,
            max_reparametrization_deviation: 0.0,
            max_profile_deviation: 0.0,
            max_guide_deviation: 0.0,
            max_approximation_deviation: 0.0,
        }
    }
}

/// `GeomFill_Gordon` (`GeomFill_Gordon.hxx:61-181`).
pub struct Gordon {
    /// `myInputProfiles`.
    input_profiles: Vec<GeomBSplineCurve>,
    /// `myInputGuides`.
    input_guides: Vec<GeomBSplineCurve>,
    /// `myProfiles`.
    profiles: Vec<GeomBSplineCurve>,
    /// `myGuides`.
    guides: Vec<GeomBSplineCurve>,
    /// `myProfileParams`, indexed `[profile][guide]`.
    profile_params: Vec<Vec<f64>>,
    /// `myGuideParams`, indexed `[profile][guide]`.
    guide_params: Vec<Vec<f64>>,
    /// `mySurface`.
    surface: Option<GeomBSplineSurface>,
    /// `myTolerance`.
    tolerance: f64,
    /// `myIsUClosed`.
    is_u_closed: bool,
    /// `myIsVClosed`.
    is_v_closed: bool,
    /// `myToUseParallel`.
    to_use_parallel: bool,
    /// `myApproximationMode`.
    approximation_mode: ApproximationMode,
    /// `myReport`.
    report: BuildReport,
}

impl Default for Gordon {
    fn default() -> Self {
        Self::new()
    }
}

impl Gordon {
    /// `GeomFill_Gordon::GeomFill_Gordon()` (`cxx:577`).
    pub fn new() -> Self {
        Self {
            input_profiles: Vec::new(),
            input_guides: Vec::new(),
            profiles: Vec::new(),
            guides: Vec::new(),
            profile_params: Vec::new(),
            guide_params: Vec::new(),
            surface: None,
            tolerance: 0.0,
            is_u_closed: false,
            is_v_closed: false,
            to_use_parallel: false,
            approximation_mode: ApproximationMode::ExactOnly,
            report: BuildReport::default(),
        }
    }

    /// `GeomFill_Gordon::Init` (`cxx:581-625`).
    ///
    /// Accepts `&dyn Curve` slices (same public shape as the prior port) and
    /// converts each curve with `GeomConvert::CurveToBSplineCurve`
    /// (`curve_to_bspline_curve`).
    pub fn init(&mut self, profiles: &[&dyn Curve], guides: &[&dyn Curve], tolerance: f64) {
        self.tolerance = tolerance;
        self.is_u_closed = false;
        self.is_v_closed = false;
        reset_report(&mut self.report);
        self.surface = None;

        let nb_prof = profiles.len();
        let nb_guid = guides.len();

        if nb_prof < 2 || nb_guid < 2 {
            // `cxx:594-601`: store counts but skip conversion.
            self.input_profiles = Vec::new();
            self.input_guides = Vec::new();
            self.profiles = Vec::new();
            self.guides = Vec::new();
            self.profile_params = Vec::new();
            self.guide_params = Vec::new();
            return;
        }

        // `cxx:606-618`.
        let mut new_profiles = Vec::with_capacity(nb_prof);
        for curve in profiles {
            match curve_to_bspline_curve(*curve, ParameterisationType::TgtThetaOver2) {
                Ok(bspline) => new_profiles.push(bspline),
                Err(_) => {
                    self.input_profiles = Vec::new();
                    self.input_guides = Vec::new();
                    self.profiles = Vec::new();
                    self.guides = Vec::new();
                    self.profile_params = Vec::new();
                    self.guide_params = Vec::new();
                    return;
                }
            }
        }
        let mut new_guides = Vec::with_capacity(nb_guid);
        for curve in guides {
            match curve_to_bspline_curve(*curve, ParameterisationType::TgtThetaOver2) {
                Ok(bspline) => new_guides.push(bspline),
                Err(_) => {
                    self.input_profiles = Vec::new();
                    self.input_guides = Vec::new();
                    self.profiles = Vec::new();
                    self.guides = Vec::new();
                    self.profile_params = Vec::new();
                    self.guide_params = Vec::new();
                    return;
                }
            }
        }
        self.input_profiles = new_profiles;
        self.input_guides = new_guides;
        // `cxx:620-624`.
        self.profile_params = vec![vec![0.0; nb_guid]; nb_prof];
        self.guide_params = vec![vec![0.0; nb_guid]; nb_prof];
        self.profiles = Vec::new();
        self.guides = Vec::new();
    }

    /// `SetParallelMode` (`hxx:139`).
    pub fn set_parallel_mode(&mut self, to_use_parallel: bool) {
        self.to_use_parallel = to_use_parallel;
    }

    /// `SetApproximationMode` (`hxx:143`).
    pub fn set_approximation_mode(&mut self, mode: ApproximationMode) {
        self.approximation_mode = mode;
    }

    /// `GetApproximationMode` (`hxx:146`).
    pub fn get_approximation_mode(&self) -> ApproximationMode {
        self.approximation_mode
    }

    /// `IsParallelMode` (`hxx:149`).
    pub fn is_parallel_mode(&self) -> bool {
        self.to_use_parallel
    }

    /// `IsDone` (`hxx:152`).
    pub fn is_done(&self) -> bool {
        self.report.status == ResultStatus::Done
    }

    /// `IsApproximate` (`hxx:156`).
    pub fn is_approximate(&self) -> bool {
        self.report.is_approximate
    }

    /// `Status` (`hxx:159`).
    pub fn status(&self) -> ResultStatus {
        self.report.status
    }

    /// `Report` (`hxx:162`).
    pub fn report(&self) -> &BuildReport {
        &self.report
    }

    /// `Surface()` (`cxx:629-636`); `None` where OCCT throws `StdFail_NotDone`.
    pub fn surface(&self) -> Option<&GeomBSplineSurface> {
        if self.is_done() {
            self.surface.as_ref()
        } else {
            None
        }
    }

    /// `GeomFill_Gordon::Perform` (`cxx:443-573`).
    pub fn perform(&mut self) {
        reset_report(&mut self.report);
        self.surface = None;
        self.is_u_closed = false;
        self.is_v_closed = false;

        if self.input_profiles.len() < 2 || self.input_guides.len() < 2 {
            set_status(
                &mut self.report,
                ResultStatus::InvalidInput,
                BuildStage::InputConversion,
            );
            return;
        }

        // `cxx:456-463`.
        let mut profiles = match copy_working_curves(&self.input_profiles) {
            Some(c) => c,
            None => {
                set_status(
                    &mut self.report,
                    ResultStatus::ConversionFailed,
                    BuildStage::InputConversion,
                );
                return;
            }
        };
        let mut guides = match copy_working_curves(&self.input_guides) {
            Some(c) => c,
            None => {
                set_status(
                    &mut self.report,
                    ResultStatus::ConversionFailed,
                    BuildStage::InputConversion,
                );
                return;
            }
        };

        let nb_prof = self.input_profiles.len();
        let nb_guid = self.input_guides.len();
        let mut profile_params = vec![vec![0.0; nb_guid]; nb_prof];
        let mut guide_params = vec![vec![0.0; nb_guid]; nb_prof];

        // NetworkPreparation borrows the working arrays; recreate between
        // stages that need immutable access to those arrays.
        {
            let mut network = NetworkPreparation::new(
                &mut profiles,
                &mut guides,
                &mut profile_params,
                &mut guide_params,
                self.tolerance,
                self.to_use_parallel,
            );

            if !network.convert_curves_to_unit_interval() {
                set_status(
                    &mut self.report,
                    ResultStatus::ConversionFailed,
                    BuildStage::InputConversion,
                );
                return;
            }

            if !network.fill_intersection_parameters() {
                set_status(
                    &mut self.report,
                    ResultStatus::IntersectionFailed,
                    BuildStage::ContactDiscovery,
                );
                return;
            }
            self.report.max_contact_gap = network.max_contact_gap();

            if !network.reorder_network() {
                set_status(
                    &mut self.report,
                    ResultStatus::OrderingFailed,
                    BuildStage::NetworkOrdering,
                );
                return;
            }
            if !network
                .normalize_intersection_parameters(&mut self.is_u_closed, &mut self.is_v_closed)
            {
                set_status(
                    &mut self.report,
                    ResultStatus::CompatibilityFailed,
                    BuildStage::NetworkOrdering,
                );
                return;
            }
        }
        self.is_u_closed = self.is_u_closed || are_closed_curves(&profiles, self.tolerance);

        let mut has_approximate_reparametrization = false;
        let allow_approx =
            self.approximation_mode == ApproximationMode::AllowApproximateFallback;
        let rational_network = has_rational_curve_pair(&profiles, &guides);
        let (intersection_points, intersection_weights) = {
            let mut network = NetworkPreparation::new(
                &mut profiles,
                &mut guides,
                &mut profile_params,
                &mut guide_params,
                self.tolerance,
                self.to_use_parallel,
            );
            if !network.equalize_intersection_parameters(
                allow_approx,
                &mut has_approximate_reparametrization,
            ) {
                let status = if rational_network {
                    ResultStatus::RationalReparametrizationFailed
                } else {
                    ResultStatus::ReparametrizationFailed
                };
                set_status(&mut self.report, status, BuildStage::Reparametrization);
                return;
            }
            self.report.max_reparametrization_deviation =
                network.max_reparametrization_deviation();

            if !network
                .normalize_intersection_parameters(&mut self.is_u_closed, &mut self.is_v_closed)
            {
                set_status(
                    &mut self.report,
                    ResultStatus::CompatibilityFailed,
                    BuildStage::Reparametrization,
                );
                return;
            }
            (network.intersection_points(), network.intersection_weights())
        };
        self.is_u_closed = self.is_u_closed || are_closed_curves(&profiles, self.tolerance);

        // `cxx:521-525`.
        let guide_param_values = column_means(&profile_params);
        let profile_param_values = row_means(&guide_params);

        let mut network_surface = NetworkSurface::new();
        network_surface.init(
            profiles.clone(),
            guides.clone(),
            profile_param_values.clone(),
            guide_param_values.clone(),
            intersection_points.clone(),
            intersection_weights,
            self.tolerance,
            self.is_u_closed,
            self.is_v_closed,
        );
        network_surface.perform();

        if !network_surface.is_done() {
            let exact_status = map_network_status(network_surface.status());
            if self.approximation_mode == ApproximationMode::AllowApproximateFallback
                && can_approximate_after_exact_failure(exact_status)
            {
                // `cxx:545-556`.
                match approximate_prepared_network(
                    &profiles,
                    &guides,
                    &intersection_points,
                    self.tolerance,
                ) {
                    Some(surface) => {
                        self.report.is_approximate = true;
                        self.report.max_approximation_deviation = max_profile_deviation(
                            &surface,
                            &profiles,
                            &profile_param_values,
                        )
                        .max(max_guide_deviation(
                            &surface,
                            &guides,
                            &guide_param_values,
                        ));
                        self.surface = Some(surface);
                        self.profiles = profiles;
                        self.guides = guides;
                        self.profile_params = profile_params;
                        self.guide_params = guide_params;
                        set_status(
                            &mut self.report,
                            ResultStatus::Done,
                            BuildStage::Approximation,
                        );
                        return;
                    }
                    None => {
                        set_status(
                            &mut self.report,
                            ResultStatus::ApproximationFailed,
                            BuildStage::Approximation,
                        );
                        return;
                    }
                }
            }
            set_status(
                &mut self.report,
                exact_status,
                BuildStage::ExactConstruction,
            );
            return;
        }

        let surface = match network_surface.surface().cloned() {
            Some(s) => s,
            None => {
                set_status(
                    &mut self.report,
                    ResultStatus::ConstructionFailed,
                    BuildStage::ExactConstruction,
                );
                return;
            }
        };
        self.report.max_profile_deviation =
            max_profile_deviation(&surface, &profiles, &profile_param_values);
        self.report.max_guide_deviation =
            max_guide_deviation(&surface, &guides, &guide_param_values);
        self.report.is_approximate = has_approximate_reparametrization;
        self.report.max_approximation_deviation = if self.report.is_approximate {
            self.report
                .max_profile_deviation
                .max(self.report.max_guide_deviation)
        } else {
            0.0
        };
        self.surface = Some(surface);
        self.profiles = profiles;
        self.guides = guides;
        self.profile_params = profile_params;
        self.guide_params = guide_params;
        set_status(&mut self.report, ResultStatus::Done, BuildStage::Validation);
    }
}

//=================================================================================================
// Anonymous helpers (`GeomFill_Gordon.cxx:49-437`).

/// `hasRationalCurve` over one family (`cxx:49-59`).
fn has_rational_curve(curves: &[GeomBSplineCurve]) -> bool {
    curves.iter().any(|c| c.is_rational())
}

/// `hasRationalCurve` over both families (`cxx:64-68`).
fn has_rational_curve_pair(
    profiles: &[GeomBSplineCurve],
    guides: &[GeomBSplineCurve],
) -> bool {
    has_rational_curve(profiles) || has_rational_curve(guides)
}

/// `mapNetworkStatus` (`cxx:73-101`).
fn map_network_status(status: NetworkResultStatus) -> ResultStatus {
    match status {
        NetworkResultStatus::Done => ResultStatus::Done,
        NetworkResultStatus::InvalidInput => ResultStatus::CompatibilityFailed,
        NetworkResultStatus::CurveCompatibilityFailed => ResultStatus::CurveCompatibilityFailed,
        NetworkResultStatus::SkinningFailed => ResultStatus::SkinningFailed,
        NetworkResultStatus::ReferenceSurfaceFailed => ResultStatus::ReferenceSurfaceFailed,
        NetworkResultStatus::KnotAlignmentFailed => ResultStatus::KnotAlignmentFailed,
        NetworkResultStatus::RationalDegreeOverflow => ResultStatus::RationalDegreeOverflow,
        NetworkResultStatus::RationalConstructionFailed => {
            ResultStatus::RationalConstructionFailed
        }
        NetworkResultStatus::PeriodicityFailed => ResultStatus::PeriodicityFailed,
        NetworkResultStatus::NotStarted | NetworkResultStatus::ConstructionFailed => {
            ResultStatus::ConstructionFailed
        }
    }
}

/// `canApproximateAfterExactFailure` (`cxx:106-115`).
fn can_approximate_after_exact_failure(status: ResultStatus) -> bool {
    matches!(
        status,
        ResultStatus::CurveCompatibilityFailed
            | ResultStatus::SkinningFailed
            | ResultStatus::ReferenceSurfaceFailed
            | ResultStatus::KnotAlignmentFailed
            | ResultStatus::RationalDegreeOverflow
            | ResultStatus::RationalConstructionFailed
            | ResultStatus::ConstructionFailed
    )
}

/// `fallbackSampleCount` (`cxx:118-129`).
fn fallback_sample_count(curves: &[GeomBSplineCurve]) -> usize {
    let mut nb_samples = curves.len().max(2);
    for curve in curves {
        nb_samples = nb_samples.max(curve.nb_poles());
    }
    nb_samples.max(MAX_DEGREE / 2)
}

/// `unitParameters` (`cxx:134-143`).
fn unit_parameters(nb_samples: usize) -> Vec<f64> {
    if nb_samples == 0 {
        return Vec::new();
    }
    if nb_samples == 1 {
        return vec![0.0];
    }
    (0..nb_samples)
        .map(|i| i as f64 / (nb_samples - 1) as f64)
        .collect()
}

/// `approximateGrid` (`cxx:148-158`).
///
/// UNPORTED: `GeomAPI_PointsToBSplineSurface::Init` with `Approx_IsoParametric`,
/// degree range `[1, MaxDegree]`, `GeomAbs_C1`, and
/// `max(tolerance, Precision::Confusion())` (`cxx:151-157`). No
/// `PointsToBSplineSurface` driver exists in this crate yet (`surface_fit.rs`
/// only provides plane/sphere/grid helpers).
fn approximate_grid(_points: &[Vec<GpPnt>], _tolerance: f64) -> Option<GeomBSplineSurface> {
    // UNPORTED: GeomAPI_PointsToBSplineSurface (GeomFill_Gordon.cxx:148-158).
    None
}

/// `profileSampleGrid` (`cxx:164-183`).
fn profile_sample_grid(profiles: &[GeomBSplineCurve], nb_u_samples: usize) -> Vec<Vec<GpPnt>> {
    let u_params = unit_parameters(nb_u_samples);
    let mut points = vec![vec![GpPnt::zero(); profiles.len()]; nb_u_samples];
    for (profile_idx, profile) in profiles.iter().enumerate() {
        for (u_idx, &u) in u_params.iter().enumerate() {
            points[u_idx][profile_idx] = profile.d0(u);
        }
    }
    points
}

/// `guideSampleGrid` (`cxx:188-207`).
fn guide_sample_grid(guides: &[GeomBSplineCurve], nb_v_samples: usize) -> Vec<Vec<GpPnt>> {
    let v_params = unit_parameters(nb_v_samples);
    let mut points = vec![vec![GpPnt::zero(); nb_v_samples]; guides.len()];
    for (guide_idx, guide) in guides.iter().enumerate() {
        for (v_idx, &v) in v_params.iter().enumerate() {
            points[guide_idx][v_idx] = guide.d0(v);
        }
    }
    points
}

/// `approximatePreparedNetwork` (`cxx:212-256`).
fn approximate_prepared_network(
    profiles: &[GeomBSplineCurve],
    guides: &[GeomBSplineCurve],
    intersections: &[Vec<GpPnt>],
    tolerance: f64,
) -> Option<GeomBSplineSurface> {
    let nb_u_samples = fallback_sample_count(profiles);
    let nb_v_samples = fallback_sample_count(guides);

    let profile_surface = approximate_grid(&profile_sample_grid(profiles, nb_u_samples), tolerance)?;
    let guide_surface = approximate_grid(&guide_sample_grid(guides, nb_v_samples), tolerance)?;
    let reference_surface = approximate_grid(intersections, tolerance)?;

    let u_params = unit_parameters(nb_u_samples);
    let v_params = unit_parameters(nb_v_samples);
    let mut points = vec![vec![GpPnt::zero(); nb_v_samples]; nb_u_samples];
    for (u_idx, &u) in u_params.iter().enumerate() {
        for (v_idx, &v) in v_params.iter().enumerate() {
            let p = profile_surface.d0(u, v);
            let g = guide_surface.d0(u, v);
            let r = reference_surface.d0(u, v);
            // `cxx:249-252`: P + G - R.
            points[u_idx][v_idx] = GpPnt::new(
                p.x() + g.x() - r.x(),
                p.y() + g.y() - r.y(),
                p.z() + g.z() - r.z(),
            );
        }
    }
    approximate_grid(&points, tolerance)
}

/// `resetReport` (`cxx:261-264`).
fn reset_report(report: &mut BuildReport) {
    *report = BuildReport::default();
}

/// `setStatus` (`cxx:269-277`).
fn set_status(report: &mut BuildReport, new_status: ResultStatus, stage: BuildStage) {
    report.status = new_status;
    report.failed_stage = if new_status == ResultStatus::Done {
        BuildStage::NotStarted
    } else {
        stage
    };
}

/// `areClosedCurves` (`cxx:282-301`).
fn are_closed_curves(curves: &[GeomBSplineCurve], tolerance: f64) -> bool {
    if curves.is_empty() {
        return false;
    }
    let tol = tolerance.max(CONFUSION);
    for curve in curves {
        let start = curve.d0(curve.first_parameter());
        let end = curve.d0(curve.last_parameter());
        if start.distance(&end) > tol {
            return false;
        }
    }
    true
}

/// `copyWorkingCurves` (`cxx:306-332`).
fn copy_working_curves(input: &[GeomBSplineCurve]) -> Option<Vec<GeomBSplineCurve>> {
    let mut working = Vec::with_capacity(input.len());
    for curve in input {
        working.push(curve.clone());
    }
    Some(working)
}

/// `surfaceStartParameter` (`cxx:337-340`).
fn surface_start_parameter(surface: &GeomBSplineSurface, is_u: bool) -> f64 {
    if is_u {
        *surface.uknots().first().unwrap_or(&0.0)
    } else {
        *surface.vknots().first().unwrap_or(&0.0)
    }
}

/// `surfaceEndParameter` (`cxx:345-349`).
fn surface_end_parameter(surface: &GeomBSplineSurface, is_u: bool) -> f64 {
    if is_u {
        *surface.uknots().last().unwrap_or(&1.0)
    } else {
        *surface.vknots().last().unwrap_or(&1.0)
    }
}

/// `surfaceParameter` (`cxx:354-360`).
fn surface_parameter(surface: &GeomBSplineSurface, is_u: bool, unit_parameter: f64) -> f64 {
    let first = surface_start_parameter(surface, is_u);
    first + unit_parameter * (surface_end_parameter(surface, is_u) - first)
}

/// `maxCurveDeviation` (`cxx:365-401`).
fn max_curve_deviation(
    surface: &GeomBSplineSurface,
    curve: &GeomBSplineCurve,
    fixed_parameter: f64,
    is_profile: bool,
) -> f64 {
    let mut max_deviation: f64 = 0.0;
    let samples_per_span = curve.degree()
        + if is_profile {
            surface.deg_u()
        } else {
            surface.deg_v()
        };
    let knots = curve.distinct_knots();
    for knot_idx in 0..knots.len().saturating_sub(1) {
        let first = knots[knot_idx];
        let last = knots[knot_idx + 1];
        if last - first <= PCONFUSION {
            continue;
        }
        for sample_idx in 0..=samples_per_span {
            let ratio = sample_idx as f64 / samples_per_span.max(1) as f64;
            let parameter = first + ratio * (last - first);
            let u = if is_profile {
                surface_parameter(surface, true, parameter)
            } else {
                surface_parameter(surface, true, fixed_parameter)
            };
            let v = if is_profile {
                surface_parameter(surface, false, fixed_parameter)
            } else {
                surface_parameter(surface, false, parameter)
            };
            max_deviation = max_deviation.max(curve.d0(parameter).distance(&surface.d0(u, v)));
        }
    }
    max_deviation
}

/// `maxProfileDeviation` (`cxx:406-419`).
fn max_profile_deviation(
    surface: &GeomBSplineSurface,
    profiles: &[GeomBSplineCurve],
    profile_params: &[f64],
) -> f64 {
    let mut max_deviation: f64 = 0.0;
    for (profile_idx, profile) in profiles.iter().enumerate() {
        let v = surface_parameter(surface, false, profile_params[profile_idx]);
        max_deviation = max_deviation.max(max_curve_deviation(surface, profile, v, true));
    }
    max_deviation
}

/// `maxGuideDeviation` (`cxx:424-437`).
fn max_guide_deviation(
    surface: &GeomBSplineSurface,
    guides: &[GeomBSplineCurve],
    guide_params: &[f64],
) -> f64 {
    let mut max_deviation: f64 = 0.0;
    for (guide_idx, guide) in guides.iter().enumerate() {
        let u = surface_parameter(surface, true, guide_params[guide_idx]);
        max_deviation = max_deviation.max(max_curve_deviation(surface, guide, u, false));
    }
    max_deviation
}
