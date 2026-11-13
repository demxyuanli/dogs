//! Network reorder helpers of `GeomFill_GordonUtilities.pxx` (`pxx:618-1067`).

use occt_core::precision::{ANGULAR, CONFUSION, REAL_SMALL};

use crate::bspline_curve::GeomBSplineCurve;
use crate::curve::Curve;
use crate::geom_fill::gordon::utilities::contacts::contact_index;
use crate::geom_fill::gordon::utilities::params::{
    column_means, parameter_resolution, row_means, OrderedCurve,
};
use crate::geom_fill::gordon::utilities::{Contact, ContactSet};

/// `orderedCurves` (`pxx:618-657`).
pub fn ordered_curves(ranks: &[f64]) -> Vec<OrderedCurve> {
    let mut order: Vec<OrderedCurve> = ranks
        .iter()
        .enumerate()
        .map(|(index, &rank)| OrderedCurve { index, rank })
        .collect();
    order.sort_by(|a, b| a.rank.partial_cmp(&b.rank).unwrap_or(std::cmp::Ordering::Equal));
    let tol = parameter_resolution(0.0, 1.0);
    if order.len() > 2 && (order[0].rank - order[1].rank).abs() <= tol {
        let mut opened = Vec::with_capacity(order.len());
        opened.push(order[0]);
        opened.extend_from_slice(&order[2..]);
        opened.push(order[1]);
        return opened;
    }
    if order.len() > 2
        && (order[order.len() - 1].rank - order[order.len() - 2].rank).abs() <= tol
    {
        let mut opened = Vec::with_capacity(order.len());
        opened.push(order[order.len() - 2]);
        opened.extend_from_slice(&order[..order.len() - 2]);
        opened.push(order[order.len() - 1]);
        return opened;
    }
    order
}

/// `isMonotoneLayout` (`pxx:659-687`).
pub fn is_monotone_layout(profile_params: &[Vec<f64>], guide_params: &[Vec<f64>]) -> bool {
    let tol = parameter_resolution(0.0, 1.0);
    for row in profile_params {
        for col in 1..row.len() {
            if row[col] <= row[col - 1] + tol {
                return false;
            }
        }
    }
    if guide_params.is_empty() {
        return true;
    }
    let nb_cols = guide_params[0].len();
    for col in 0..nb_cols {
        for row in 1..guide_params.len() {
            if guide_params[row][col] <= guide_params[row - 1][col] + tol {
                return false;
            }
        }
    }
    true
}

/// `openClosedBoundaryCandidates` (`pxx:689-742`).
pub fn open_closed_boundary_candidates(
    profiles: &[GeomBSplineCurve],
    guides: &[GeomBSplineCurve],
    profile_params: &mut [Vec<f64>],
    guide_params: &mut [Vec<f64>],
    tolerance: f64,
) {
    let tol = tolerance.max(CONFUSION);
    let last_guide = guides.len() - 1;
    let last_prof = profiles.len() - 1;
    let mut is_u_repeated = true;
    for profile_idx in 0..profiles.len() {
        let first = profiles[profile_idx].d0(profile_params[profile_idx][0]);
        let last = profiles[profile_idx].d0(profile_params[profile_idx][last_guide]);
        if first.distance(&last) > tol {
            is_u_repeated = false;
            break;
        }
    }
    if is_u_repeated {
        for profile_idx in 0..profiles.len() {
            profile_params[profile_idx][0] = profiles[profile_idx].first_parameter();
            profile_params[profile_idx][last_guide] = profiles[profile_idx].last_parameter();
        }
    }
    let mut is_v_repeated = true;
    for guide_idx in 0..guides.len() {
        let first = guides[guide_idx].d0(guide_params[0][guide_idx]);
        let last = guides[guide_idx].d0(guide_params[last_prof][guide_idx]);
        if first.distance(&last) > tol {
            is_v_repeated = false;
            break;
        }
    }
    if is_v_repeated {
        for guide_idx in 0..guides.len() {
            guide_params[0][guide_idx] = guides[guide_idx].first_parameter();
            guide_params[last_prof][guide_idx] = guides[guide_idx].last_parameter();
        }
    }
}

/// `orientedCurve` (`pxx:744-754`).
pub fn oriented_curve(curve: &GeomBSplineCurve, is_forward: bool) -> GeomBSplineCurve {
    let mut out = curve.clone();
    if !is_forward {
        out.reverse();
    }
    out
}

/// `orientedParameter` (`pxx:756-763`).
pub fn oriented_parameter(curve: &GeomBSplineCurve, parameter: f64, is_forward: bool) -> f64 {
    if is_forward {
        parameter
    } else {
        curve.first_parameter() + curve.last_parameter() - parameter
    }
}

/// `orderedContactScore` (`pxx:765-795`).
#[allow(clippy::too_many_arguments)]
pub fn ordered_contact_score(
    contact: &Contact,
    profile_param: f64,
    guide_param: f64,
    profile_target: f64,
    guide_target: f64,
    previous_profile_param: f64,
    previous_guide_param: f64,
    previous_profile_target: f64,
    previous_guide_target: f64,
    has_previous_profile: bool,
    has_previous_guide: bool,
    tolerance: f64,
) -> f64 {
    let mut score = contact.gap / tolerance.max(CONFUSION)
        + (profile_param - profile_target).abs()
        + (guide_param - guide_target).abs();
    if has_previous_profile {
        score += ((profile_param - previous_profile_param)
            - (profile_target - previous_profile_target))
            .abs();
    }
    if has_previous_guide {
        score +=
            ((guide_param - previous_guide_param) - (guide_target - previous_guide_target)).abs();
    }
    score
}

/// `seamTangentMismatch` (`pxx:797-836`).
pub fn seam_tangent_mismatch(
    curve: &GeomBSplineCurve,
    first_parameter: f64,
    last_parameter: f64,
    is_forward: bool,
    tolerance: f64,
    mismatch: &mut f64,
) -> bool {
    let first_raw = if is_forward {
        first_parameter
    } else {
        curve.first_parameter() + curve.last_parameter() - first_parameter
    };
    let last_raw = if is_forward {
        last_parameter
    } else {
        curve.first_parameter() + curve.last_parameter() - last_parameter
    };
    let (first_pnt, mut first_tangent) = curve.d1(first_raw);
    let (last_pnt, mut last_tangent) = curve.d1(last_raw);
    if first_pnt.distance(&last_pnt) > tolerance.max(CONFUSION) {
        return false;
    }
    if !is_forward {
        first_tangent.reverse();
        last_tangent.reverse();
    }
    let first_mag = first_tangent.magnitude();
    let last_mag = last_tangent.magnitude();
    if first_mag <= REAL_SMALL || last_mag <= REAL_SMALL {
        return false;
    }
    let cos = (first_tangent.dot(&last_tangent) / (first_mag * last_mag)).clamp(-1.0, 1.0);
    *mismatch = 1.0 - cos;
    true
}

/// `chooseOrderedContact` (`pxx:838-931`).
#[allow(clippy::too_many_arguments)]
pub fn choose_ordered_contact(
    contacts: &ContactSet,
    profile: &GeomBSplineCurve,
    guide: &GeomBSplineCurve,
    profile_forward: bool,
    guide_forward: bool,
    profile_target: f64,
    guide_target: f64,
    previous_profile_param: f64,
    previous_guide_param: f64,
    has_previous_profile: bool,
    has_previous_guide: bool,
    previous_profile_target: f64,
    previous_guide_target: f64,
    first_profile_param: f64,
    first_guide_param: f64,
    is_profile_seam: bool,
    is_guide_seam: bool,
    tolerance: f64,
    selected: &mut Contact,
) -> bool {
    let param_tol = parameter_resolution(0.0, 1.0);
    let mut best_tangent_mismatch = f64::MAX;
    let mut best = f64::MAX;
    for candidate in &contacts.candidates {
        let profile_param =
            oriented_parameter(profile, candidate.profile_param, profile_forward);
        let guide_param = oriented_parameter(guide, candidate.guide_param, guide_forward);
        if (has_previous_profile && profile_param <= previous_profile_param + param_tol)
            || (has_previous_guide && guide_param <= previous_guide_param + param_tol)
        {
            continue;
        }
        let score = ordered_contact_score(
            candidate,
            profile_param,
            guide_param,
            profile_target,
            guide_target,
            previous_profile_param,
            previous_guide_param,
            previous_profile_target,
            previous_guide_target,
            has_previous_profile,
            has_previous_guide,
            tolerance,
        );
        let mut tangent_mismatch = 0.0;
        let mut has_tangent_match = false;
        if is_profile_seam {
            let mut profile_mismatch = 0.0;
            if seam_tangent_mismatch(
                profile,
                first_profile_param,
                profile_param,
                profile_forward,
                tolerance,
                &mut profile_mismatch,
            ) {
                tangent_mismatch += profile_mismatch;
                has_tangent_match = true;
            }
        }
        if is_guide_seam {
            let mut guide_mismatch = 0.0;
            if seam_tangent_mismatch(
                guide,
                first_guide_param,
                guide_param,
                guide_forward,
                tolerance,
                &mut guide_mismatch,
            ) {
                tangent_mismatch += guide_mismatch;
                has_tangent_match = true;
            }
        }
        let take = if has_tangent_match {
            tangent_mismatch < best_tangent_mismatch - ANGULAR
                || ((tangent_mismatch - best_tangent_mismatch).abs() <= ANGULAR && score < best)
        } else {
            best_tangent_mismatch == f64::MAX && score < best
        };
        if take {
            best_tangent_mismatch = if has_tangent_match {
                tangent_mismatch
            } else {
                f64::MAX
            };
            best = score;
            *selected = *candidate;
            selected.profile_param = profile_param;
            selected.guide_param = guide_param;
        }
    }
    selected.is_valid
}

/// `reorderNetwork` (`pxx:933-1067`).
pub fn reorder_network(
    profiles: &mut Vec<GeomBSplineCurve>,
    guides: &mut Vec<GeomBSplineCurve>,
    profile_params: &mut Vec<Vec<f64>>,
    guide_params: &mut Vec<Vec<f64>>,
    contact_sets: &[ContactSet],
    tolerance: f64,
) -> bool {
    if profiles.len() < 2 || guides.len() < 2 {
        return false;
    }
    let profile_order = ordered_curves(&row_means(guide_params));
    let guide_order = ordered_curves(&column_means(profile_params));

    let mut new_profiles = Vec::with_capacity(profiles.len());
    let mut profile_forward = Vec::with_capacity(profiles.len());
    for order in &profile_order {
        let src = order.index;
        let first_guide = guide_order[0].index;
        let last_guide = guide_order[guide_order.len() - 1].index;
        let is_forward =
            profile_params[src][first_guide] <= profile_params[src][last_guide];
        new_profiles.push(oriented_curve(&profiles[src], is_forward));
        profile_forward.push(is_forward);
    }

    let mut new_guides = Vec::with_capacity(guides.len());
    let mut guide_forward = Vec::with_capacity(guides.len());
    for order in &guide_order {
        let src = order.index;
        let first_profile = profile_order[0].index;
        let last_profile = profile_order[profile_order.len() - 1].index;
        let is_forward =
            guide_params[first_profile][src] <= guide_params[last_profile][src];
        new_guides.push(oriented_curve(&guides[src], is_forward));
        guide_forward.push(is_forward);
    }

    let nb_prof = profile_order.len();
    let nb_guid = guide_order.len();
    let mut new_profile_params = vec![vec![0.0; nb_guid]; nb_prof];
    let mut new_guide_params = vec![vec![0.0; nb_guid]; nb_prof];

    for profile_idx in 0..nb_prof {
        let src_profile = profile_order[profile_idx].index;
        for guide_idx in 0..nb_guid {
            let src_guide = guide_order[guide_idx].index;
            let profile_target = if nb_guid == 1 {
                0.0
            } else {
                guide_idx as f64 / (nb_guid - 1) as f64
            };
            let guide_target = if nb_prof == 1 {
                0.0
            } else {
                profile_idx as f64 / (nb_prof - 1) as f64
            };
            let previous_profile_param = if guide_idx == 0 {
                0.0
            } else {
                new_profile_params[profile_idx][guide_idx - 1]
            };
            let previous_guide_param = if profile_idx == 0 {
                0.0
            } else {
                new_guide_params[profile_idx - 1][guide_idx]
            };
            let previous_profile_target = if guide_idx < 2 {
                0.0
            } else {
                (guide_idx - 1) as f64 / (nb_guid - 1) as f64
            };
            let previous_guide_target = if profile_idx < 2 {
                0.0
            } else {
                (profile_idx - 1) as f64 / (nb_prof - 1) as f64
            };
            let is_profile_seam = guide_idx + 1 == nb_guid;
            let is_guide_seam = profile_idx + 1 == nb_prof;
            let first_profile_param = if is_profile_seam {
                new_profile_params[profile_idx][0]
            } else {
                0.0
            };
            let first_guide_param = if is_guide_seam {
                new_guide_params[0][guide_idx]
            } else {
                0.0
            };
            let mut contact = Contact::default();
            if !choose_ordered_contact(
                &contact_sets[contact_index(src_profile, src_guide, guides.len())],
                &profiles[src_profile],
                &guides[src_guide],
                profile_forward[profile_idx],
                guide_forward[guide_idx],
                profile_target,
                guide_target,
                previous_profile_param,
                previous_guide_param,
                guide_idx != 0,
                profile_idx != 0,
                previous_profile_target,
                previous_guide_target,
                first_profile_param,
                first_guide_param,
                is_profile_seam,
                is_guide_seam,
                tolerance,
                &mut contact,
            ) {
                return false;
            }
            new_profile_params[profile_idx][guide_idx] = contact.profile_param;
            new_guide_params[profile_idx][guide_idx] = contact.guide_param;
        }
    }

    open_closed_boundary_candidates(
        &new_profiles,
        &new_guides,
        &mut new_profile_params,
        &mut new_guide_params,
        tolerance,
    );
    if !is_monotone_layout(&new_profile_params, &new_guide_params) {
        return false;
    }
    *profiles = new_profiles;
    *guides = new_guides;
    *profile_params = new_profile_params;
    *guide_params = new_guide_params;
    true
}
