//! Contact discovery of `GeomFill_GordonUtilities.pxx` (`pxx:212-616`,
//! `:1550-1615`).

use occt_core::gp::{GpPnt, GpXyz};
use occt_core::precision::{CONFUSION, PCONFUSION, REAL_SMALL};

use crate::bspline_curve::GeomBSplineCurve;
use crate::curve::Curve;
use crate::geom_api::{project_point_on_curve, ExtremaCurveCurve};
use crate::geom_fill::gordon::utilities::params::{
    boundary_parameter_resolution, parameter_resolution,
};
use crate::geom_fill::gordon::utilities::{Contact, ContactSet};

/// `GeomLib_Tool::Parameter` via `GeomAPI_ProjectPointOnCurve`
/// (`pxx:265-266`, `:285-288`, `:296-299`).
fn curve_parameter(curve: &GeomBSplineCurve, point: &GpPnt, tolerance: f64) -> Option<f64> {
    let poc = project_point_on_curve(curve, point, tolerance)?;
    if poc.distance <= tolerance {
        Some(poc.parameter)
    } else {
        None
    }
}

/// `appendContact` (`pxx:212-245`).
pub fn append_contact(
    profile: &GeomBSplineCurve,
    guide: &GeomBSplineCurve,
    profile_param: f64,
    guide_param: f64,
    tolerance: f64,
    contacts: &mut ContactSet,
) {
    let gap = profile.d0(profile_param).distance(&guide.d0(guide_param));
    if gap > tolerance {
        return;
    }
    let profile_tol = parameter_resolution(profile.first_parameter(), profile.last_parameter());
    let guide_tol = parameter_resolution(guide.first_parameter(), guide.last_parameter());
    for contact in &contacts.candidates {
        if (contact.profile_param - profile_param).abs() <= profile_tol
            && (contact.guide_param - guide_param).abs() <= guide_tol
        {
            return;
        }
    }
    contacts.candidates.push(Contact {
        is_valid: true,
        profile_param,
        guide_param,
        gap,
    });
}

/// `appendExtremaContact` (`pxx:247-270`).
pub fn append_extrema_contact(
    profile: &GeomBSplineCurve,
    guide: &GeomBSplineCurve,
    profile_point: &GpPnt,
    guide_point: &GpPnt,
    tolerance: f64,
    contacts: &mut ContactSet,
) {
    if profile_point.distance(guide_point) > tolerance {
        return;
    }
    let mid = GpXyz::new(
        (profile_point.x() + guide_point.x()) * 0.5,
        (profile_point.y() + guide_point.y()) * 0.5,
        (profile_point.z() + guide_point.z()) * 0.5,
    );
    let shared = GpPnt::from_xyz(&mid);
    if let (Some(profile_param), Some(guide_param)) = (
        curve_parameter(profile, &shared, tolerance),
        curve_parameter(guide, &shared, tolerance),
    ) {
        append_contact(
            profile,
            guide,
            profile_param,
            guide_param,
            tolerance,
            contacts,
        );
    }
}

/// `appendBoundaryContacts` (`pxx:272-304`).
pub fn append_boundary_contacts(
    profile: &GeomBSplineCurve,
    guide: &GeomBSplineCurve,
    tolerance: f64,
    contacts: &mut ContactSet,
) {
    let profile_bounds = [profile.first_parameter(), profile.last_parameter()];
    let guide_bounds = [guide.first_parameter(), guide.last_parameter()];
    for &profile_param in &profile_bounds {
        if let Some(guide_param) =
            curve_parameter(guide, &profile.d0(profile_param), tolerance)
        {
            append_contact(
                profile,
                guide,
                profile_param,
                guide_param,
                tolerance,
                contacts,
            );
        }
    }
    for &guide_param in &guide_bounds {
        if let Some(profile_param) =
            curve_parameter(profile, &guide.d0(guide_param), tolerance)
        {
            append_contact(
                profile,
                guide,
                profile_param,
                guide_param,
                tolerance,
                contacts,
            );
        }
    }
}

/// `appendLinearSegmentContact` (`pxx:306-375`).
pub fn append_linear_segment_contact(
    profile: &GeomBSplineCurve,
    guide: &GeomBSplineCurve,
    tolerance: f64,
    contacts: &mut ContactSet,
) -> bool {
    if profile.degree() != 1
        || guide.degree() != 1
        || profile.nb_poles() != 2
        || guide.nb_poles() != 2
    {
        return false;
    }
    let profile_first = profile.first_parameter();
    let profile_last = profile.last_parameter();
    let guide_first = guide.first_parameter();
    let guide_last = guide.last_parameter();
    if profile_last <= profile_first + PCONFUSION || guide_last <= guide_first + PCONFUSION {
        return true;
    }
    let p0 = profile.d0(profile_first);
    let p1 = profile.d0(profile_last);
    let q0 = guide.d0(guide_first);
    let q1 = guide.d0(guide_last);
    let u = GpXyz::new(p1.x() - p0.x(), p1.y() - p0.y(), p1.z() - p0.z());
    let v = GpXyz::new(q1.x() - q0.x(), q1.y() - q0.y(), q1.z() - q0.z());
    let w = GpXyz::new(p0.x() - q0.x(), p0.y() - q0.y(), p0.z() - q0.z());
    let a = u.square_modulus();
    let b = u.dot(&v);
    let c = v.square_modulus();
    let d = u.dot(&w);
    let e = v.dot(&w);
    if a <= REAL_SMALL || c <= REAL_SMALL {
        return true;
    }
    let denom = a * c - b * b;
    if denom.abs() <= REAL_SMALL {
        return false;
    }
    let mut s = ((b * e - c * d) / denom).clamp(0.0, 1.0);
    let mut t = (b * s + e) / c;
    if t < 0.0 {
        t = 0.0;
        s = (-d / a).clamp(0.0, 1.0);
    } else if t > 1.0 {
        t = 1.0;
        s = ((b - d) / a).clamp(0.0, 1.0);
    }
    append_contact(
        profile,
        guide,
        profile_first + s * (profile_last - profile_first),
        guide_first + t * (guide_last - guide_first),
        tolerance,
        contacts,
    );
    true
}

/// Closest extrema points — stands in for UNPORTED
/// `GeomAPI_ExtremaCurveCurve::TotalNearestPoints` (`pxx:399`,
/// `GeomAPI_ExtremaCurveCurve.cxx:281-320`): same closest index as `myIndex`.
fn total_nearest_points(extrema: &ExtremaCurveCurve) -> Option<(GpPnt, GpPnt)> {
    let idx = extrema.index();
    if idx < 1 || idx > extrema.nb_extrema() {
        return None;
    }
    Some(extrema.points(idx))
}

/// `collectContacts` (`pxx:377-419`).
pub fn collect_contacts(
    profile: &GeomBSplineCurve,
    guide: &GeomBSplineCurve,
    tolerance: f64,
) -> ContactSet {
    let mut contacts = ContactSet::default();
    if append_linear_segment_contact(profile, guide, tolerance, &mut contacts) {
        return contacts;
    }
    let extrema = ExtremaCurveCurve::new(profile, guide);
    for ext_idx in 1..=extrema.nb_extrema() {
        let (profile_point, guide_point) = extrema.points(ext_idx);
        append_extrema_contact(
            profile,
            guide,
            &profile_point,
            &guide_point,
            tolerance,
            &mut contacts,
        );
    }
    if let Some((profile_point, guide_point)) = total_nearest_points(&extrema) {
        append_extrema_contact(
            profile,
            guide,
            &profile_point,
            &guide_point,
            tolerance,
            &mut contacts,
        );
    }
    append_boundary_contacts(profile, guide, tolerance, &mut contacts);
    contacts.candidates.sort_by(|left, right| {
        left.gap
            .partial_cmp(&right.gap)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| {
                left.profile_param
                    .partial_cmp(&right.profile_param)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .then_with(|| {
                left.guide_param
                    .partial_cmp(&right.guide_param)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
    });
    contacts
}

/// `contactIndex` (`pxx:421-425`).
pub fn contact_index(profile_idx: usize, guide_idx: usize, nb_guides: usize) -> usize {
    profile_idx * nb_guides + guide_idx
}

/// `chooseNearestContacts` (`pxx:427-449`).
pub fn choose_nearest_contacts(
    contacts: &[ContactSet],
    nb_profiles: usize,
    nb_guides: usize,
    profile_params: &mut [Vec<f64>],
    guide_params: &mut [Vec<f64>],
) -> bool {
    for profile_idx in 0..nb_profiles {
        for guide_idx in 0..nb_guides {
            let set = &contacts[contact_index(profile_idx, guide_idx, nb_guides)];
            let Some(contact) = set.candidates.first() else {
                return false;
            };
            profile_params[profile_idx][guide_idx] = contact.profile_param;
            guide_params[profile_idx][guide_idx] = contact.guide_param;
        }
    }
    true
}

/// `maxNearestContactGap` (`pxx:451-464`).
pub fn max_nearest_contact_gap(contacts: &[ContactSet]) -> f64 {
    let mut max_gap: f64 = 0.0;
    for set in contacts {
        if let Some(first) = set.candidates.first() {
            max_gap = max_gap.max(first.gap);
        }
    }
    max_gap
}

/// `fillContactTable` (`pxx:466-514`).
///
/// `OSD_Parallel::For` with `!theToUseParallel` (`pxx:486-506`) is sequential
/// when parallel mode is off; this port always runs the pair loop sequentially.
pub fn fill_contact_table(
    profiles: &[GeomBSplineCurve],
    guides: &[GeomBSplineCurve],
    profile_params: &mut [Vec<f64>],
    guide_params: &mut [Vec<f64>],
    contact_sets: &mut Vec<ContactSet>,
    tolerance: f64,
    _to_use_parallel: bool,
) -> bool {
    let nb_profiles = profiles.len();
    let nb_guides = guides.len();
    let nb_pairs = nb_profiles * nb_guides;
    contact_sets.clear();
    contact_sets.resize_with(nb_pairs, ContactSet::default);
    let mut is_done = true;
    for pair_idx in 0..nb_pairs {
        if !is_done {
            break;
        }
        let profile_idx = pair_idx / nb_guides;
        let guide_idx = pair_idx % nb_guides;
        let set = collect_contacts(&profiles[profile_idx], &guides[guide_idx], tolerance);
        if set.candidates.is_empty() {
            is_done = false;
        }
        contact_sets[contact_index(profile_idx, guide_idx, nb_guides)] = set;
    }
    is_done
        && choose_nearest_contacts(
            contact_sets,
            nb_profiles,
            nb_guides,
            profile_params,
            guide_params,
        )
}

/// `snappedParameter` (`pxx:516-533`).
pub fn snapped_parameter(curve: &GeomBSplineCurve, parameter: f64, tolerance: f64) -> f64 {
    let first = curve.first_parameter();
    let last = curve.last_parameter();
    let tol = boundary_parameter_resolution(first, last, tolerance);
    if (parameter - first).abs() < tol {
        return first;
    }
    if (parameter - last).abs() < tol {
        return last;
    }
    parameter
}

/// `normalizeContactTable` (`pxx:535-616`).
pub fn normalize_contact_table(
    profiles: &[GeomBSplineCurve],
    guides: &[GeomBSplineCurve],
    profile_params: &mut [Vec<f64>],
    guide_params: &mut [Vec<f64>],
    tolerance: f64,
    is_u_closed: &mut bool,
    is_v_closed: &mut bool,
) -> bool {
    let tol = tolerance.max(CONFUSION);
    for profile_idx in 0..profiles.len() {
        for guide_idx in 0..guides.len() {
            profile_params[profile_idx][guide_idx] = snapped_parameter(
                &profiles[profile_idx],
                profile_params[profile_idx][guide_idx],
                tolerance,
            );
            guide_params[profile_idx][guide_idx] = snapped_parameter(
                &guides[guide_idx],
                guide_params[profile_idx][guide_idx],
                tolerance,
            );
        }
    }
    let mut u_gap: f64 = 0.0;
    let mut v_gap: f64 = 0.0;
    let last_profile = profiles.len() - 1;
    let last_guide = guides.len() - 1;
    for profile_idx in 0..profiles.len() {
        for guide_idx in 0..guides.len() {
            let profile_point =
                profiles[profile_idx].d0(profile_params[profile_idx][guide_idx]);
            let guide_point = guides[guide_idx].d0(guide_params[profile_idx][guide_idx]);
            if profile_point.distance(&guide_point) > tol {
                return false;
            }
            if guide_idx == 0 || guide_idx == last_guide {
                let opposite = profiles[profile_idx]
                    .d0(profile_params[profile_idx][last_guide - guide_idx]);
                u_gap = f64::max(u_gap, profile_point.distance(&opposite));
            }
            if profile_idx == 0 || profile_idx == last_profile {
                let opposite = guides[guide_idx]
                    .d0(guide_params[last_profile - profile_idx][guide_idx]);
                v_gap = f64::max(v_gap, guide_point.distance(&opposite));
            }
        }
    }
    *is_u_closed = u_gap <= tol;
    *is_v_closed = v_gap <= tol;
    if *is_u_closed {
        for profile_idx in 0..profiles.len() {
            profile_params[profile_idx][0] = profiles[profile_idx].first_parameter();
            profile_params[profile_idx][last_guide] = profiles[profile_idx].last_parameter();
        }
    }
    if *is_v_closed {
        for guide_idx in 0..guides.len() {
            guide_params[0][guide_idx] = guides[guide_idx].first_parameter();
            guide_params[last_profile][guide_idx] = guides[guide_idx].last_parameter();
        }
    }
    true
}

/// `contactPointGrid` (`pxx:1550-1570`): guide/profile order.
pub fn contact_point_grid(
    profiles: &[GeomBSplineCurve],
    profile_params: &[Vec<f64>],
) -> Vec<Vec<GpPnt>> {
    let nb_guides = profile_params.first().map_or(0, |r| r.len());
    let nb_profiles = profiles.len();
    let mut points = vec![vec![GpPnt::new(0.0, 0.0, 0.0); nb_profiles]; nb_guides];
    for profile_idx in 0..nb_profiles {
        for guide_idx in 0..nb_guides {
            points[guide_idx][profile_idx] =
                profiles[profile_idx].d0(profile_params[profile_idx][guide_idx]);
        }
    }
    points
}

/// `curveWeightAt` (`pxx:1572-1594`).
pub fn curve_weight_at(curve: &GeomBSplineCurve, parameter: f64) -> f64 {
    if !curve.is_rational() {
        return 1.0;
    }
    let weights = match curve.weights.as_ref() {
        Some(w) => w,
        None => return 1.0,
    };
    let mut prepared = occt_core::bspl::curve_eval::prepare_curve_eval(
        parameter,
        0,
        curve.degree() as i32,
        curve.periodic(),
        &curve.poles,
        Some(weights),
        &curve.knots,
        None,
    );
    occt_core::bspl::bohm::eval(
        prepared.u,
        curve.degree() as i32,
        &prepared.knots,
        prepared.dim,
        &mut prepared.poles,
    );
    if prepared.rational {
        prepared.poles[3]
    } else {
        1.0
    }
}

/// `contactWeightGrid` (`pxx:1596-1615`): guide/profile order.
pub fn contact_weight_grid(
    profiles: &[GeomBSplineCurve],
    profile_params: &[Vec<f64>],
) -> Vec<Vec<f64>> {
    let nb_guides = profile_params.first().map_or(0, |r| r.len());
    let nb_profiles = profiles.len();
    let mut weights = vec![vec![1.0; nb_profiles]; nb_guides];
    for profile_idx in 0..nb_profiles {
        for guide_idx in 0..nb_guides {
            weights[guide_idx][profile_idx] =
                curve_weight_at(&profiles[profile_idx], profile_params[profile_idx][guide_idx]);
        }
    }
    weights
}
