//! `GeomFill_Gordon` — high-level Gordon surface from a curve network.
//!
//! Source: `ModelingAlgorithms/TKGeomAlgo/GeomFill/GeomFill_Gordon.cxx`
//! (OCCT 8.0.0), `GeomFill_Gordon.hxx`.
//!
//! The class converts the input curves to B-splines, finds every
//! profile/guide intersection, sorts the network, snaps and reparametrizes the
//! intersection parameters, unifies the curves with `GeomFill_Profiler`, and
//! delegates the boolean-sum construction to [`super::GordonBuilder`].
//!
//! Translation notes:
//! - `Surface()` throws `StdFail_NotDone` before `Perform` (`cxx:610-618`);
//!   the port returns `None`.
//! - `OSD_Parallel::For` (`cxx:725-738`) only selects the executor; the port
//!   runs the intersection loop serially. Every intersection is independent, so
//!   the results do not depend on that policy.
//! - `Init` with fewer than two profiles or guides stores the counts but no
//!   curves, and `Perform` returns before touching them (`cxx:452-459`,
//!   `cxx:488-493`); the port stores empty vectors, which `perform` rejects the
//!   same way.
//! - A curve that `GeomConvert::CurveToBSplineCurve` cannot convert leaves
//!   `IsDone()` false (`cxx:466-477`), like a failed conversion.

use occt_core::convert::ParameterisationType;
use occt_core::gp::GpPnt;
use occt_core::precision::PCONFUSION;

use crate::bspline_curve::GeomBSplineCurve;
use crate::bspline_surface::GeomBSplineSurface;
use crate::convert_bspl::curve_to_bspline_curve;
use crate::curve::Curve;
use crate::extrema_cc::curve_curve_extrema_all;
use crate::geom_fill::gordon::builder::GordonBuilder;
use crate::geom_fill::gordon::reparam_curve::reparam_curve;
use crate::geom_fill::profiler::Profiler;

/// `THE_SNAP_TOL` (`GeomFill_Gordon.cxx:46-47`): snap tolerance for boundary
/// parameters.
const THE_SNAP_TOL: f64 = 1.0e-5;
/// `THE_REL_TOL_CLOSED` (`cxx:48-49`): relative closedness tolerance.
const THE_REL_TOL_CLOSED: f64 = 1.0e-8;
/// `THE_NETWORK_CHECK_TOL_FACTOR` (`cxx:85-86`): relative tolerance for the
/// post-reparametrization compatibility check.
const THE_NETWORK_CHECK_TOL_FACTOR: f64 = 1.0e-5;

/// `GeomFill_Gordon` (`GeomFill_Gordon.hxx:57-136`).
pub struct Gordon {
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
    /// `myScale`.
    scale: f64,
    /// `myIsUClosed`.
    is_u_closed: bool,
    /// `myIsVClosed`.
    is_v_closed: bool,
    /// `myToUseParallel`.
    to_use_parallel: bool,
    /// `myIsDone`.
    is_done: bool,
}

impl Default for Gordon {
    fn default() -> Self {
        Self::new()
    }
}

impl Gordon {
    /// `GeomFill_Gordon::GeomFill_Gordon()` (`cxx:434`).
    pub fn new() -> Self {
        Self {
            profiles: Vec::new(),
            guides: Vec::new(),
            profile_params: Vec::new(),
            guide_params: Vec::new(),
            surface: None,
            tolerance: 0.0,
            scale: 1.0,
            is_u_closed: false,
            is_v_closed: false,
            to_use_parallel: false,
            is_done: false,
        }
    }

    /// `GeomFill_Gordon::Init` (`cxx:438-483`).
    pub fn init(&mut self, profiles: &[&dyn Curve], guides: &[&dyn Curve], tolerance: f64) {
        self.tolerance = tolerance;
        self.is_done = false;
        self.scale = 1.0;
        self.is_u_closed = false;
        self.is_v_closed = false;
        self.surface = None;

        let nb_prof = profiles.len();
        let nb_guid = guides.len();
        if nb_prof < 2 || nb_guid < 2 {
            self.profiles = Vec::new();
            self.guides = Vec::new();
            return;
        }

        // `cxx:462-478`: `GeomConvert::CurveToBSplineCurve` on every input.
        let mut new_profiles = Vec::with_capacity(nb_prof);
        for curve in profiles {
            match curve_to_bspline_curve(*curve, ParameterisationType::TgtThetaOver2) {
                Ok(bspline) => new_profiles.push(bspline),
                Err(_) => {
                    self.profiles = Vec::new();
                    self.guides = Vec::new();
                    return;
                }
            }
        }
        let mut new_guides = Vec::with_capacity(nb_guid);
        for curve in guides {
            match curve_to_bspline_curve(*curve, ParameterisationType::TgtThetaOver2) {
                Ok(bspline) => new_guides.push(bspline),
                Err(_) => {
                    self.profiles = Vec::new();
                    self.guides = Vec::new();
                    return;
                }
            }
        }
        self.profiles = new_profiles;
        self.guides = new_guides;
        // `cxx:479-481`.
        self.profile_params = vec![vec![0.0; nb_guid]; nb_prof];
        self.guide_params = vec![vec![0.0; nb_guid]; nb_prof];
    }

    /// `SetParallelMode` (`GeomFill_Gordon.hxx:77`).
    pub fn set_parallel_mode(&mut self, to_use_parallel: bool) {
        self.to_use_parallel = to_use_parallel;
    }

    /// `IsParallelMode` (`GeomFill_Gordon.hxx:79`).
    pub fn is_parallel_mode(&self) -> bool {
        self.to_use_parallel
    }

    /// `IsDone` (`GeomFill_Gordon.hxx:82`).
    pub fn is_done(&self) -> bool {
        self.is_done
    }

    /// `Surface()` (`cxx:610-618`); `None` where OCCT throws `StdFail_NotDone`.
    pub fn surface(&self) -> Option<&GeomBSplineSurface> {
        if self.is_done {
            self.surface.as_ref()
        } else {
            None
        }
    }

    /// `GeomFill_Gordon::Perform` (`cxx:487-606`).
    pub fn perform(&mut self) {
        self.is_done = false;
        if self.profiles.len() < 2 || self.guides.len() < 2 {
            return;
        }
        if !self.convert_to_bspline() {
            return;
        }
        if !self.compute_intersections() {
            return;
        }
        if !self.sort_network() {
            return;
        }
        self.compute_scale();
        self.eliminate_inaccuracies();
        self.detect_closedness();
        if !self.reparametrize() {
            return;
        }
        if !self.check_network_compatibility() {
            return;
        }

        // `cxx:534-558`: unify the profiles, then the guides.
        let mut profile_profiler = Profiler::default();
        for curve in &self.profiles {
            if profile_profiler.add_curve(curve).is_err() {
                return;
            }
        }
        if profile_profiler.perform(PCONFUSION).is_err() {
            return;
        }
        self.profiles = profile_profiler.curves().to_vec();

        let mut guide_profiler = Profiler::default();
        for curve in &self.guides {
            if guide_profiler.add_curve(curve).is_err() {
                return;
            }
        }
        if guide_profiler.perform(PCONFUSION).is_err() {
            return;
        }
        self.guides = guide_profiler.curves().to_vec();

        // `cxx:561-587`: averaged intersection parameters for the builder.
        let nb_prof = self.profiles.len();
        let nb_guid = self.guides.len();
        let mut guide_param_values = vec![0.0f64; nb_guid];
        for guide_idx in 0..nb_guid {
            let mut sum = 0.0;
            for profile_idx in 0..nb_prof {
                sum += self.profile_params[profile_idx][guide_idx];
            }
            guide_param_values[guide_idx] = sum / nb_prof as f64;
        }
        let mut profile_param_values = vec![0.0f64; nb_prof];
        for profile_idx in 0..nb_prof {
            let mut sum = 0.0;
            for guide_idx in 0..nb_guid {
                sum += self.guide_params[profile_idx][guide_idx];
            }
            profile_param_values[profile_idx] = sum / nb_guid as f64;
        }

        // `cxx:589-601`.
        let mut builder = GordonBuilder::new();
        builder.set_parallel_mode(self.to_use_parallel);
        builder.init(
            self.profiles.clone(),
            self.guides.clone(),
            profile_param_values,
            guide_param_values,
            self.tolerance,
            self.is_u_closed,
            self.is_v_closed,
        );
        builder.perform();
        if !builder.is_done() {
            return;
        }
        self.surface = builder.surface().cloned();
        self.is_done = true;
    }

    /// `convertToBSpline` (`cxx:622-658`): reject rational/periodic curves and
    /// reparametrize every curve to `[0, 1]`.
    fn convert_to_bspline(&mut self) -> bool {
        for curve in self.profiles.iter_mut() {
            if curve.is_rational() || curve.is_periodic() {
                return false;
            }
            let mut knots = curve.distinct_knots().to_vec();
            occt_core::bspl::knots::reparametrize(0.0, 1.0, &mut knots);
            if curve.set_knots(&knots).is_err() {
                return false;
            }
        }
        for curve in self.guides.iter_mut() {
            if curve.is_rational() || curve.is_periodic() {
                return false;
            }
            let mut knots = curve.distinct_knots().to_vec();
            occt_core::bspl::knots::reparametrize(0.0, 1.0, &mut knots);
            if curve.set_knots(&knots).is_err() {
                return false;
            }
        }
        true
    }

    /// `computeIntersections` (`cxx:663-741`).
    fn compute_intersections(&mut self) -> bool {
        let nb_prof = self.profiles.len();
        let nb_guid = self.guides.len();
        for profile_idx in 0..nb_prof {
            for guide_idx in 0..nb_guid {
                let extrema =
                    curve_curve_extrema_all(&self.profiles[profile_idx], &self.guides[guide_idx]);
                // `cxx:687-691`: the closest extremum below the tolerance.
                let mut best: Option<(f64, f64, f64)> = None;
                for ext in &extrema {
                    if ext.distance < self.tolerance
                        && best.map_or(true, |b| ext.distance < b.2)
                    {
                        best = Some((ext.u1, ext.u2, ext.distance));
                    }
                }
                // `cxx:695-710` retries through `TotalNearestPoints`, which
                // `GeomAPI_ExtremaCurveCurve::TotalPerform`
                // (`GeomAPI_ExtremaCurveCurve.cxx:281-320`) resolves to the same
                // closest extremum; it cannot succeed once the loop above
                // rejected every extremum (no extremum at all, or all of them
                // farther than the tolerance).
                match best {
                    Some((prof_param, guid_param, _)) => {
                        self.profile_params[profile_idx][guide_idx] = prof_param;
                        self.guide_params[profile_idx][guide_idx] = guid_param;
                    }
                    None => return false,
                }
            }
        }
        true
    }

    /// `computeScale` (`cxx:745-798`): bounding-box diagonal of the curve
    /// endpoints, clamped below by 1.
    fn compute_scale(&mut self) {
        if self.profiles.is_empty() && self.guides.is_empty() {
            self.scale = 1.0;
            return;
        }
        let seed_pnt = if !self.profiles.is_empty() {
            self.profiles[0].d0(self.profiles[0].first_parameter())
        } else {
            self.guides[0].d0(self.guides[0].first_parameter())
        };
        let mut min = seed_pnt;
        let mut max = seed_pnt;
        let mut update_bbox = |pnt: &GpPnt| {
            min.set_x(min.x().min(pnt.x()));
            min.set_y(min.y().min(pnt.y()));
            min.set_z(min.z().min(pnt.z()));
            max.set_x(max.x().max(pnt.x()));
            max.set_y(max.y().max(pnt.y()));
            max.set_z(max.z().max(pnt.z()));
        };
        for curve in self.profiles.iter().chain(self.guides.iter()) {
            update_bbox(&curve.d0(curve.first_parameter()));
            update_bbox(&curve.d0(curve.last_parameter()));
        }
        let diag = max.distance(&min);
        self.scale = diag.max(1.0);
    }

    /// `detectClosedness` (`cxx:802-858`).
    fn detect_closedness(&mut self) {
        let nb_prof = self.profiles.len();
        let nb_guid = self.guides.len();
        let close_tol = THE_REL_TOL_CLOSED * self.scale;

        // U direction: first and last guides coincide and the first/last grid
        // columns match (`cxx:814-831`).
        {
            let first_guide = &self.guides[0];
            let last_guide = &self.guides[nb_guid - 1];
            let guid_first_start = first_guide.d0(first_guide.first_parameter());
            let guid_first_end = first_guide.d0(first_guide.last_parameter());
            let guid_last_start = last_guide.d0(last_guide.first_parameter());
            let guid_last_end = last_guide.d0(last_guide.last_parameter());
            let is_guides_equal = guid_first_start.distance(&guid_last_start) <= close_tol
                && guid_first_end.distance(&guid_last_end) <= close_tol;
            let mut is_grid_closed = true;
            for profile_idx in 0..nb_prof {
                let p1 = self.profiles[profile_idx].d0(self.profile_params[profile_idx][0]);
                let p2 =
                    self.profiles[profile_idx].d0(self.profile_params[profile_idx][nb_guid - 1]);
                if p1.distance(&p2) > close_tol {
                    is_grid_closed = false;
                    break;
                }
            }
            self.is_u_closed = is_guides_equal && is_grid_closed;
        }

        // V direction: first and last profiles coincide and the first/last grid
        // rows match (`cxx:835-856`).
        {
            let first_profile = &self.profiles[0];
            let last_profile = &self.profiles[nb_prof - 1];
            let prof_first_start = first_profile.d0(first_profile.first_parameter());
            let prof_first_end = first_profile.d0(first_profile.last_parameter());
            let prof_last_start = last_profile.d0(last_profile.first_parameter());
            let prof_last_end = last_profile.d0(last_profile.last_parameter());
            let is_profiles_equal = prof_first_start.distance(&prof_last_start) <= close_tol
                && prof_first_end.distance(&prof_last_end) <= close_tol;
            let mut is_grid_closed = true;
            for guide_idx in 0..nb_guid {
                let p1 = self.guides[guide_idx].d0(self.guide_params[0][guide_idx]);
                let p2 =
                    self.guides[guide_idx].d0(self.guide_params[nb_prof - 1][guide_idx]);
                if p1.distance(&p2) > close_tol {
                    is_grid_closed = false;
                    break;
                }
            }
            self.is_v_closed = is_profiles_equal && is_grid_closed;
        }
    }

    /// `eliminateInaccuracies` (`cxx:862-906`): snap parameters near curve
    /// boundaries to the exact boundary values.
    fn eliminate_inaccuracies(&mut self) {
        let nb_prof = self.profiles.len();
        let nb_guid = self.guides.len();
        for profile_idx in 0..nb_prof {
            let first = self.profiles[profile_idx].first_parameter();
            let last = self.profiles[profile_idx].last_parameter();
            let range = last - first;
            for guide_idx in 0..nb_guid {
                let param = &mut self.profile_params[profile_idx][guide_idx];
                if (*param - first).abs() < THE_SNAP_TOL * range {
                    *param = first;
                } else if (*param - last).abs() < THE_SNAP_TOL * range {
                    *param = last;
                }
            }
        }
        for guide_idx in 0..nb_guid {
            let first = self.guides[guide_idx].first_parameter();
            let last = self.guides[guide_idx].last_parameter();
            let range = last - first;
            for profile_idx in 0..nb_prof {
                let param = &mut self.guide_params[profile_idx][guide_idx];
                if (*param - first).abs() < THE_SNAP_TOL * range {
                    *param = first;
                } else if (*param - last).abs() < THE_SNAP_TOL * range {
                    *param = last;
                }
            }
        }
    }

    /// `checkNetworkCompatibility` (`cxx:912-937`).
    fn check_network_compatibility(&self) -> bool {
        let check_tol = self.scale * THE_NETWORK_CHECK_TOL_FACTOR;
        for profile_idx in 0..self.profiles.len() {
            for guide_idx in 0..self.guides.len() {
                let prof_pt = self.profiles[profile_idx]
                    .d0(self.profile_params[profile_idx][guide_idx]);
                let guid_pt = self.guides[guide_idx]
                    .d0(self.guide_params[profile_idx][guide_idx]);
                if prof_pt.distance(&guid_pt) > check_tol {
                    return false;
                }
            }
        }
        true
    }

    /// `sortNetwork` (`cxx:942-1136`).
    fn sort_network(&mut self) -> bool {
        let nb_prof = self.profiles.len();
        let nb_guid = self.guides.len();
        if nb_prof < 2 || nb_guid < 2 {
            return false;
        }

        // `cxx:950-990`: the profile/guide pair whose intersection is the
        // minimum parameter on both curves.
        let mut start_prof = 0usize;
        let mut start_guid = 0usize;
        let mut guide_must_reverse = false;
        let mut is_start_found = false;
        for profile_idx in 0..nb_prof {
            let mut min_guide_idx = 0usize;
            let mut min_val = self.profile_params[profile_idx][0];
            for guide_idx in 1..nb_guid {
                if self.profile_params[profile_idx][guide_idx] < min_val {
                    min_val = self.profile_params[profile_idx][guide_idx];
                    min_guide_idx = guide_idx;
                }
            }
            let mut min_profile_idx = 0usize;
            let mut min_guide_val = self.guide_params[0][min_guide_idx];
            for profile_scan in 1..nb_prof {
                if self.guide_params[profile_scan][min_guide_idx] < min_guide_val {
                    min_guide_val = self.guide_params[profile_scan][min_guide_idx];
                    min_profile_idx = profile_scan;
                }
            }
            if min_profile_idx == profile_idx {
                start_prof = profile_idx;
                start_guid = min_guide_idx;
                guide_must_reverse = false;
                is_start_found = true;
                break;
            }
        }

        // `cxx:993-1027`: fallback for an oppositely oriented guide.
        if !is_start_found {
            for profile_idx in 0..nb_prof {
                let mut min_guide_idx = 0usize;
                let mut min_val = self.profile_params[profile_idx][0];
                for guide_idx in 1..nb_guid {
                    if self.profile_params[profile_idx][guide_idx] < min_val {
                        min_val = self.profile_params[profile_idx][guide_idx];
                        min_guide_idx = guide_idx;
                    }
                }
                let mut max_profile_idx = 0usize;
                let mut max_guide_val = self.guide_params[0][min_guide_idx];
                for profile_scan in 1..nb_prof {
                    if self.guide_params[profile_scan][min_guide_idx] > max_guide_val {
                        max_guide_val = self.guide_params[profile_scan][min_guide_idx];
                        max_profile_idx = profile_scan;
                    }
                }
                if max_profile_idx == profile_idx {
                    start_prof = profile_idx;
                    start_guid = min_guide_idx;
                    guide_must_reverse = true;
                    is_start_found = true;
                    break;
                }
            }
        }
        if !is_start_found {
            return false;
        }

        // `cxx:1087-1091`: move the starting curves to position 0.
        self.swap_profiles(0, start_prof);
        self.swap_guides(0, start_guid);
        if guide_must_reverse {
            self.reverse_guide(0);
        }

        // `cxx:1093-1106`: bubble-sort guides by their parameter on the first
        // profile.
        let mut sort_end = nb_guid;
        while sort_end > 1 {
            for guide_idx in 1..sort_end {
                if self.profile_params[0][guide_idx - 1]
                    > self.profile_params[0][guide_idx]
                {
                    self.swap_guides(guide_idx - 1, guide_idx);
                }
            }
            sort_end -= 1;
        }

        // `cxx:1108-1114`: bubble-sort profiles by their parameter on the first
        // guide.
        let mut sort_end = nb_prof;
        while sort_end > 1 {
            for profile_idx in 1..sort_end {
                if self.guide_params[profile_idx - 1][0]
                    > self.guide_params[profile_idx][0]
                {
                    self.swap_profiles(profile_idx - 1, profile_idx);
                }
            }
            sort_end -= 1;
        }

        // `cxx:1117-1134`: reverse the individual curves whose parameters run
        // backwards.
        for profile_idx in 1..nb_prof {
            if self.profile_params[profile_idx][0] > self.profile_params[profile_idx][nb_guid - 1]
            {
                self.reverse_profile(profile_idx);
            }
        }
        for guide_idx in 1..nb_guid {
            if self.guide_params[0][guide_idx] > self.guide_params[nb_prof - 1][guide_idx] {
                self.reverse_guide(guide_idx);
            }
        }
        true
    }

    /// `swapProfiles` (`cxx:1035-1048`).
    fn swap_profiles(&mut self, i1: usize, i2: usize) {
        if i1 == i2 {
            return;
        }
        self.profiles.swap(i1, i2);
        self.profile_params.swap(i1, i2);
        self.guide_params.swap(i1, i2);
    }

    /// `swapGuides` (`cxx:1050-1061`).
    fn swap_guides(&mut self, g1: usize, g2: usize) {
        if g1 == g2 {
            return;
        }
        self.guides.swap(g1, g2);
        for profile_idx in 0..self.profile_params.len() {
            self.profile_params[profile_idx].swap(g1, g2);
            self.guide_params[profile_idx].swap(g1, g2);
        }
    }

    /// `reverseProfile` (`cxx:1063-1072`): mirror the parameters and reverse
    /// the curve.
    fn reverse_profile(&mut self, profile_idx: usize) {
        let first = self.profiles[profile_idx].first_parameter();
        let last = self.profiles[profile_idx].last_parameter();
        for guide_idx in 0..self.guides.len() {
            let param = &mut self.profile_params[profile_idx][guide_idx];
            *param = -*param + first + last;
        }
        self.profiles[profile_idx].reverse();
    }

    /// `reverseGuide` (`cxx:1074-1082`).
    fn reverse_guide(&mut self, guide_idx: usize) {
        let first = self.guides[guide_idx].first_parameter();
        let last = self.guides[guide_idx].last_parameter();
        for profile_idx in 0..self.profiles.len() {
            let param = &mut self.guide_params[profile_idx][guide_idx];
            *param = -*param + first + last;
        }
        self.guides[guide_idx].reverse();
    }

    /// `reparametrize` (`cxx:1141-1213`).
    fn reparametrize(&mut self) -> bool {
        let nb_prof = self.profiles.len();
        let nb_guid = self.guides.len();

        let mut avg_profile_params = vec![0.0f64; nb_guid];
        for guide_idx in 0..nb_guid {
            let mut sum = 0.0;
            for profile_idx in 0..nb_prof {
                sum += self.profile_params[profile_idx][guide_idx];
            }
            avg_profile_params[guide_idx] = sum / nb_prof as f64;
        }
        let mut avg_guide_params = vec![0.0f64; nb_prof];
        for profile_idx in 0..nb_prof {
            let mut sum = 0.0;
            for guide_idx in 0..nb_guid {
                sum += self.guide_params[profile_idx][guide_idx];
            }
            avg_guide_params[profile_idx] = sum / nb_guid as f64;
        }

        for profile_idx in 0..nb_prof {
            let old_params = self.profile_params[profile_idx].clone();
            let new_params = avg_profile_params.clone();
            if reparam_curve(&mut self.profiles[profile_idx], &old_params, &new_params).is_err() {
                return false;
            }
            self.profile_params[profile_idx] = avg_profile_params.clone();
        }
        for guide_idx in 0..nb_guid {
            let mut old_params = Vec::with_capacity(nb_prof);
            for profile_idx in 0..nb_prof {
                old_params.push(self.guide_params[profile_idx][guide_idx]);
            }
            let new_params = avg_guide_params.clone();
            if reparam_curve(&mut self.guides[guide_idx], &old_params, &new_params).is_err() {
                return false;
            }
            for profile_idx in 0..nb_prof {
                self.guide_params[profile_idx][guide_idx] = avg_guide_params[profile_idx];
            }
        }
        true
    }
}
