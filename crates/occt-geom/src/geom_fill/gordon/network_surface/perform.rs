//! Construction helpers of `GeomFill_NetworkSurface.cxx` (OCCT 8.0.1):
//! skins, seam completion, alignment, corrected profile, periodicity.

use occt_core::gp::GpPnt;
use occt_core::precision::{CONFUSION, PCONFUSION};

use crate::bspline_curve::GeomBSplineCurve;
use crate::bspline_surface::GeomBSplineSurface;
use crate::curve::Curve;
use crate::geom_api_interpolate::GeomAPIInterpolate;

use super::helpers::{
    gp_resolution, has_rational_curve, interpolate_poles, interpolate_surface_weighted,
    make_surface, make_triple_product_basis, multiply_by_denominators, pole_weight,
    prepare_curve_family, unite_knots, SkinningBasis,
};
use super::ResultStatus;

/// `makeProfileSkin` (`cxx:394-441`).
pub(super) fn make_profile_skin(
    profiles: &[GeomBSplineCurve],
    profile_parameters: &[f64],
    v_degree: i32,
    v_knots: &[f64],
    v_mults: &[i32],
    v_flat_knots: &[f64],
) -> Option<GeomBSplineSurface> {
    let base = profiles.first()?;
    let nb_u = base.nb_poles();
    let nb_v = profiles.len();
    let mut poles = vec![vec![GpPnt::zero(); nb_v]; nb_u];
    let mut weights = vec![vec![1.0; nb_v]; nb_u];

    for u_pole_idx in 1..=nb_u as i32 {
        let mut column = Vec::with_capacity(nb_v);
        let mut column_weights = Vec::with_capacity(nb_v);
        for profile in profiles {
            column.push(*profile.pole((u_pole_idx - 1) as usize));
            column_weights.push(pole_weight(profile, u_pole_idx));
        }
        if !interpolate_poles(
            v_degree,
            v_flat_knots,
            profile_parameters,
            &mut column,
            &mut column_weights,
        ) {
            return None;
        }
        for v_idx in 0..nb_v {
            poles[(u_pole_idx - 1) as usize][v_idx] = column[v_idx];
            weights[(u_pole_idx - 1) as usize][v_idx] = column_weights[v_idx];
        }
    }

    make_surface(
        poles,
        Some(weights),
        base.distinct_knots().to_vec(),
        v_knots.to_vec(),
        base.multiplicities().to_vec(),
        v_mults.to_vec(),
        base.degree(),
        v_degree as usize,
    )
}

/// `makeGuideSkin` (`cxx:443-490`).
pub(super) fn make_guide_skin(
    guides: &[GeomBSplineCurve],
    guide_parameters: &[f64],
    u_degree: i32,
    u_knots: &[f64],
    u_mults: &[i32],
    u_flat_knots: &[f64],
) -> Option<GeomBSplineSurface> {
    let base = guides.first()?;
    let nb_u = guides.len();
    let nb_v = base.nb_poles();
    let mut poles = vec![vec![GpPnt::zero(); nb_v]; nb_u];
    let mut weights = vec![vec![1.0; nb_v]; nb_u];

    for v_pole_idx in 1..=nb_v as i32 {
        let mut row = Vec::with_capacity(nb_u);
        let mut row_weights = Vec::with_capacity(nb_u);
        for guide in guides {
            row.push(*guide.pole((v_pole_idx - 1) as usize));
            row_weights.push(pole_weight(guide, v_pole_idx));
        }
        if !interpolate_poles(
            u_degree,
            u_flat_knots,
            guide_parameters,
            &mut row,
            &mut row_weights,
        ) {
            return None;
        }
        for u_idx in 0..nb_u {
            poles[u_idx][(v_pole_idx - 1) as usize] = row[u_idx];
            weights[u_idx][(v_pole_idx - 1) as usize] = row_weights[u_idx];
        }
    }

    make_surface(
        poles,
        Some(weights),
        u_knots.to_vec(),
        base.distinct_knots().to_vec(),
        u_mults.to_vec(),
        base.multiplicities().to_vec(),
        u_degree as usize,
        base.degree(),
    )
}

/// `alignSurfaces` (`cxx:594-689`).
pub(super) fn align_surfaces(
    profile_surface: &mut GeomBSplineSurface,
    guide_surface: &mut GeomBSplineSurface,
    reference_surface: &mut GeomBSplineSurface,
) -> bool {
    let u_degree = profile_surface
        .deg_u()
        .max(guide_surface.deg_u())
        .max(reference_surface.deg_u());
    let v_degree = profile_surface
        .deg_v()
        .max(guide_surface.deg_v())
        .max(reference_surface.deg_v());

    if profile_surface.increase_degree(u_degree, v_degree).is_err()
        || guide_surface.increase_degree(u_degree, v_degree).is_err()
        || reference_surface.increase_degree(u_degree, v_degree).is_err()
    {
        return false;
    }

    let Some((u_knots12, u_mults12)) = unite_knots(
        u_degree as i32,
        profile_surface.uknots(),
        profile_surface.umults(),
        guide_surface.uknots(),
        guide_surface.umults(),
    ) else {
        return false;
    };
    let Some((u_knots, u_mults)) = unite_knots(
        u_degree as i32,
        &u_knots12,
        &u_mults12,
        reference_surface.uknots(),
        reference_surface.umults(),
    ) else {
        return false;
    };

    let Some((v_knots12, v_mults12)) = unite_knots(
        v_degree as i32,
        profile_surface.vknots(),
        profile_surface.vmults(),
        guide_surface.vknots(),
        guide_surface.vmults(),
    ) else {
        return false;
    };
    let Some((v_knots, v_mults)) = unite_knots(
        v_degree as i32,
        &v_knots12,
        &v_mults12,
        reference_surface.vknots(),
        reference_surface.vmults(),
    ) else {
        return false;
    };

    if profile_surface
        .insert_u_knots(&u_knots, &u_mults, PCONFUSION, false)
        .is_err()
        || guide_surface
            .insert_u_knots(&u_knots, &u_mults, PCONFUSION, false)
            .is_err()
        || reference_surface
            .insert_u_knots(&u_knots, &u_mults, PCONFUSION, false)
            .is_err()
        || profile_surface
            .insert_v_knots(&v_knots, &v_mults, PCONFUSION, false)
            .is_err()
        || guide_surface
            .insert_v_knots(&v_knots, &v_mults, PCONFUSION, false)
            .is_err()
        || reference_surface
            .insert_v_knots(&v_knots, &v_mults, PCONFUSION, false)
            .is_err()
    {
        return false;
    }
    true
}

/// `makeCorrectedProfileSkin` (`cxx:692-836`).
pub(super) fn make_corrected_profile_skin(
    profile_surface: &GeomBSplineSurface,
    guide_surface: &GeomBSplineSurface,
    reference_surface: &GeomBSplineSurface,
    status: &mut ResultStatus,
) -> Option<GeomBSplineSurface> {
    if profile_surface.nb_poles_u() != guide_surface.nb_poles_u()
        || profile_surface.nb_poles_u() != reference_surface.nb_poles_u()
        || profile_surface.nb_poles_v() != guide_surface.nb_poles_v()
        || profile_surface.nb_poles_v() != reference_surface.nb_poles_v()
    {
        *status = ResultStatus::RationalConstructionFailed;
        return None;
    }

    if !profile_surface.is_u_rational()
        && !profile_surface.is_v_rational()
        && !guide_surface.is_u_rational()
        && !guide_surface.is_v_rational()
        && !reference_surface.is_u_rational()
        && !reference_surface.is_v_rational()
    {
        let nu = profile_surface.nb_poles_u();
        let nv = profile_surface.nb_poles_v();
        let mut result_poles = vec![vec![GpPnt::zero(); nv]; nu];
        for i in 0..nu {
            for j in 0..nv {
                let xyz = profile_surface.pole(i, j).xyz().added(guide_surface.pole(i, j).xyz())
                    .subtracted(reference_surface.pole(i, j).xyz());
                result_poles[i][j] = GpPnt::from_xyz(&xyz);
            }
        }
        return make_surface(
            result_poles,
            None,
            profile_surface.uknots().to_vec(),
            profile_surface.vknots().to_vec(),
            profile_surface.umults().to_vec(),
            profile_surface.vmults().to_vec(),
            profile_surface.deg_u(),
            profile_surface.deg_v(),
        );
    }

    let Some((u_degree, u_knots, u_mults, u_flat_knots)) = make_triple_product_basis(
        profile_surface.uknots(),
        profile_surface.umults(),
        profile_surface.deg_u() as i32,
    ) else {
        *status = ResultStatus::RationalDegreeOverflow;
        return None;
    };
    let Some((v_degree, v_knots, v_mults, v_flat_knots)) = make_triple_product_basis(
        profile_surface.vknots(),
        profile_surface.vmults(),
        profile_surface.deg_v() as i32,
    ) else {
        *status = ResultStatus::RationalDegreeOverflow;
        return None;
    };

    let nb_u_poles = u_flat_knots.len() as i32 - u_degree - 1;
    let nb_v_poles = v_flat_knots.len() as i32 - v_degree - 1;
    if nb_u_poles < 2 || nb_v_poles < 2 {
        *status = ResultStatus::RationalConstructionFailed;
        return None;
    }

    let Some((profile_num, profile_den)) = multiply_by_denominators(
        profile_surface,
        guide_surface,
        reference_surface,
        u_degree,
        v_degree,
        &u_flat_knots,
        &v_flat_knots,
    ) else {
        *status = ResultStatus::RationalConstructionFailed;
        return None;
    };
    let Some((guide_num, guide_den)) = multiply_by_denominators(
        guide_surface,
        profile_surface,
        reference_surface,
        u_degree,
        v_degree,
        &u_flat_knots,
        &v_flat_knots,
    ) else {
        *status = ResultStatus::RationalConstructionFailed;
        return None;
    };
    let Some((reference_num, reference_den)) = multiply_by_denominators(
        reference_surface,
        profile_surface,
        guide_surface,
        u_degree,
        v_degree,
        &u_flat_knots,
        &v_flat_knots,
    ) else {
        *status = ResultStatus::RationalConstructionFailed;
        return None;
    };

    let nu = profile_num.len();
    let nv = profile_num.first().map_or(0, |r| r.len());
    let mut result_poles = vec![vec![GpPnt::zero(); nv]; nu];
    let mut result_weights = vec![vec![0.0; nv]; nu];
    for i in 0..nu {
        for j in 0..nv {
            let weight = profile_den[i][j];
            if weight <= CONFUSION
                || (guide_den[i][j] - weight).abs() > CONFUSION
                || (reference_den[i][j] - weight).abs() > CONFUSION
            {
                *status = ResultStatus::RationalConstructionFailed;
                return None;
            }
            let numerator = profile_num[i][j]
                .xyz()
                .added(guide_num[i][j].xyz())
                .subtracted(reference_num[i][j].xyz());
            result_weights[i][j] = weight;
            result_poles[i][j] = GpPnt::from_xyz(&numerator.divided(weight));
        }
    }

    make_surface(
        result_poles,
        Some(result_weights),
        u_knots,
        v_knots,
        u_mults,
        v_mults,
        u_degree as usize,
        v_degree as usize,
    )
}

/// `canSetPeriodic` (`cxx:838-855`).
pub(super) fn can_set_periodic(
    surface: &GeomBSplineSurface,
    along_u: bool,
    tolerance: f64,
) -> bool {
    let nb_checks = if along_u {
        surface.nb_poles_v()
    } else {
        surface.nb_poles_u()
    };
    for check_idx in 0..nb_checks {
        let (first, last) = if along_u {
            (
                surface.pole(0, check_idx),
                surface.pole(surface.nb_poles_u() - 1, check_idx),
            )
        } else {
            (
                surface.pole(check_idx, 0),
                surface.pole(check_idx, surface.nb_poles_v() - 1),
            )
        };
        if first.distance(&last) > tolerance {
            return false;
        }
    }
    true
}

/// `applyPeriodicity` (`cxx:1141-1165`).
pub(super) fn apply_periodicity(
    surface: &mut GeomBSplineSurface,
    is_u_closed: bool,
    is_v_closed: bool,
    tolerance: f64,
) -> bool {
    let tolerance = tolerance.max(CONFUSION);
    if is_u_closed {
        if !can_set_periodic(surface, true, tolerance) {
            return false;
        }
        surface.set_u_periodic();
    }
    if is_v_closed {
        if !can_set_periodic(surface, false, tolerance) {
            return false;
        }
        surface.set_v_periodic();
    }
    true
}

/// `makeProfileSeamGuide` (`cxx:876-895`).
pub(super) fn make_profile_seam_guide(
    profiles: &[GeomBSplineCurve],
    profile_parameters: &[f64],
) -> Option<GeomBSplineCurve> {
    let mut points = Vec::with_capacity(profiles.len());
    for profile in profiles {
        points.push(profile.d0(profile.first_parameter()));
    }
    let mut interpolator =
        GeomAPIInterpolate::with_parameters(points, profile_parameters.to_vec(), false, CONFUSION)
            .ok()?;
    interpolator.perform();
    interpolator.curve().cloned()
}

/// `completeClosedGuideSeam` (`cxx:898-995`).
#[allow(clippy::too_many_arguments)]
pub(super) fn complete_closed_guide_seam(
    profiles: &[GeomBSplineCurve],
    profile_parameters: &[f64],
    guides: &[GeomBSplineCurve],
    guide_parameters: &[f64],
    intersection_points: &[Vec<GpPnt>],
    intersection_weights: &[Vec<f64>],
) -> Option<(
    Vec<GeomBSplineCurve>,
    Vec<f64>,
    Vec<Vec<GpPnt>>,
    Vec<Vec<f64>>,
)> {
    let profile_start = profiles[0].first_parameter();
    let profile_end = profiles[0].last_parameter();
    let needs_start_seam = guide_parameters[0] > profile_start + PCONFUSION;
    let needs_end_seam = *guide_parameters.last().unwrap() < profile_end - PCONFUSION;
    if !needs_start_seam && !needs_end_seam {
        return Some((
            guides.to_vec(),
            guide_parameters.to_vec(),
            intersection_points.to_vec(),
            intersection_weights.to_vec(),
        ));
    }

    let seam_guide = if needs_start_seam {
        Some(make_profile_seam_guide(profiles, profile_parameters)?)
    } else {
        None
    };

    let guide_offset = usize::from(needs_start_seam);
    let nb_guides = guides.len() + guide_offset + usize::from(needs_end_seam);
    let mut closed_guides = Vec::with_capacity(nb_guides);
    let mut closed_guide_parameters = Vec::with_capacity(nb_guides);
    if needs_start_seam {
        closed_guides.push(seam_guide.clone().unwrap());
        closed_guide_parameters.push(profile_start);
    }
    for (idx, guide) in guides.iter().enumerate() {
        closed_guides.push(guide.clone());
        closed_guide_parameters.push(guide_parameters[idx]);
    }
    if needs_end_seam {
        if needs_start_seam {
            closed_guides.push(seam_guide.unwrap());
        } else {
            closed_guides.push(guides[0].clone());
        }
        closed_guide_parameters.push(profile_end);
    }

    let nb_profiles = intersection_points.first().map_or(0, |r| r.len());
    let mut closed_points = vec![vec![GpPnt::zero(); nb_profiles]; nb_guides];
    let mut closed_weights = vec![vec![1.0; nb_profiles]; nb_guides];
    if needs_start_seam {
        for profile_idx in 0..nb_profiles {
            closed_points[0][profile_idx] =
                profiles[profile_idx].d0(profiles[profile_idx].first_parameter());
            closed_weights[0][profile_idx] = pole_weight(&profiles[profile_idx], 1);
        }
    }
    for guide_idx in 0..intersection_points.len() {
        for profile_idx in 0..nb_profiles {
            closed_points[guide_idx + guide_offset][profile_idx] =
                intersection_points[guide_idx][profile_idx];
            closed_weights[guide_idx + guide_offset][profile_idx] =
                intersection_weights[guide_idx][profile_idx];
        }
    }
    if needs_end_seam {
        for profile_idx in 0..nb_profiles {
            closed_points[nb_guides - 1][profile_idx] = closed_points[0][profile_idx];
            closed_weights[nb_guides - 1][profile_idx] = closed_weights[0][profile_idx];
        }
    }
    Some((
        closed_guides,
        closed_guide_parameters,
        closed_points,
        closed_weights,
    ))
}

/// `makeNetworkSurface` (`cxx:997-1139`).
pub(super) fn make_network_surface(
    profiles: &[GeomBSplineCurve],
    guides: &[GeomBSplineCurve],
    profile_parameters: &[f64],
    guide_parameters: &[f64],
    intersection_points: &[Vec<GpPnt>],
    intersection_weights: &[Vec<f64>],
    is_u_closed: bool,
    status: &mut ResultStatus,
) -> Option<GeomBSplineSurface> {
    *status = ResultStatus::ConstructionFailed;

    let should_complete_u_seam = is_u_closed
        && (guide_parameters[0] > profiles[0].first_parameter() + PCONFUSION
            || *guide_parameters.last().unwrap() < profiles[0].last_parameter() - PCONFUSION);

    let (
        closed_guides,
        closed_guide_parameters,
        closed_intersection_points,
        closed_intersection_weights,
    ) = if should_complete_u_seam {
        match complete_closed_guide_seam(
            profiles,
            profile_parameters,
            guides,
            guide_parameters,
            intersection_points,
            intersection_weights,
        ) {
            Some(v) => v,
            None => {
                *status = ResultStatus::SkinningFailed;
                return None;
            }
        }
    } else {
        (
            guides.to_vec(),
            guide_parameters.to_vec(),
            intersection_points.to_vec(),
            intersection_weights.to_vec(),
        )
    };

    let mut guides_work = closed_guides;
    if should_complete_u_seam && !prepare_curve_family(&mut guides_work) {
        *status = ResultStatus::CurveCompatibilityFailed;
        return None;
    }

    let guide_parameters_work = if should_complete_u_seam {
        &closed_guide_parameters
    } else {
        guide_parameters
    };
    let intersection_points_work = if should_complete_u_seam {
        &closed_intersection_points
    } else {
        intersection_points
    };
    let intersection_weights_work = if should_complete_u_seam {
        &closed_intersection_weights
    } else {
        intersection_weights
    };

    let u_basis = SkinningBasis::init(guide_parameters_work);
    let v_basis = SkinningBasis::init(profile_parameters);

    let Some(mut profile_surface) = make_profile_skin(
        profiles,
        profile_parameters,
        v_basis.degree,
        &v_basis.knots,
        &v_basis.mults,
        &v_basis.flat_knots,
    ) else {
        *status = ResultStatus::SkinningFailed;
        return None;
    };

    let Some(mut guide_surface) = make_guide_skin(
        &guides_work,
        guide_parameters_work,
        u_basis.degree,
        &u_basis.knots,
        &u_basis.mults,
        &u_basis.flat_knots,
    ) else {
        *status = ResultStatus::SkinningFailed;
        return None;
    };

    // Reference surface from the contact grid (`cxx:1077-1132`).
    let mut reference_poles = intersection_points_work.to_vec();
    let mut reference_weights = intersection_weights_work.to_vec();
    let is_reference_polynomial =
        has_rational_curve(profiles) != has_rational_curve(&guides_work);
    if is_reference_polynomial {
        for row in &mut reference_weights {
            for w in row.iter_mut() {
                *w = 1.0;
            }
        }
    } else {
        for i in 0..reference_poles.len() {
            for j in 0..reference_poles[i].len() {
                let w = reference_weights[i][j];
                let xyz = reference_poles[i][j].xyz().multiplied(w);
                reference_poles[i][j] = GpPnt::from_xyz(&xyz);
            }
        }
    }

    if !interpolate_surface_weighted(
        u_basis.degree,
        v_basis.degree,
        &u_basis.flat_knots,
        &v_basis.flat_knots,
        guide_parameters_work,
        profile_parameters,
        &mut reference_poles,
        &mut reference_weights,
    ) {
        *status = ResultStatus::ReferenceSurfaceFailed;
        return None;
    }
    if !is_reference_polynomial {
        for i in 0..reference_poles.len() {
            for j in 0..reference_poles[i].len() {
                let w = reference_weights[i][j];
                if w <= gp_resolution() {
                    *status = ResultStatus::ReferenceSurfaceFailed;
                    return None;
                }
                let xyz = reference_poles[i][j].xyz().divided(w);
                reference_poles[i][j] = GpPnt::from_xyz(&xyz);
            }
        }
    }

    let Some(mut reference_surface) = make_surface(
        reference_poles,
        Some(reference_weights),
        u_basis.knots,
        v_basis.knots,
        u_basis.mults,
        v_basis.mults,
        u_basis.degree as usize,
        v_basis.degree as usize,
    ) else {
        *status = ResultStatus::ReferenceSurfaceFailed;
        return None;
    };

    if !align_surfaces(
        &mut profile_surface,
        &mut guide_surface,
        &mut reference_surface,
    ) {
        *status = ResultStatus::KnotAlignmentFailed;
        return None;
    }
    make_corrected_profile_skin(
        &profile_surface,
        &guide_surface,
        &reference_surface,
        status,
    )
}
