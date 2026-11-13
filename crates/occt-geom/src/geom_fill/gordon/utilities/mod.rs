//! `GordonUtilities` from `GeomFill_GordonUtilities.pxx` (OCCT 8.0.1).
//!
//! Included by `GeomFill_Gordon.cxx` via `#include "GeomFill_GordonUtilities.pxx"`.
//! Split by theme: [`params`], [`contacts`], [`reorder`], [`reparam`].

pub mod contacts;
pub mod params;
pub mod reparam;
pub mod reorder;

use crate::bspline_curve::GeomBSplineCurve;
use occt_core::gp::GpPnt;

pub use contacts::{
    append_boundary_contacts, append_contact, append_extrema_contact,
    append_linear_segment_contact, choose_nearest_contacts, collect_contacts, contact_index,
    contact_point_grid, contact_weight_grid, curve_weight_at, fill_contact_table,
    max_nearest_contact_gap, normalize_contact_table, snapped_parameter,
};
pub use params::{
    add_interval_midpoints, append_unit_parameter, boundary_parameter_resolution, column_means,
    compact_parameters, interpolation_degree, parameter_resolution, push_constraint, row_means,
    sampling_parameters, OrderedCurve,
};
pub use reparam::{
    equalize_parameters, expose_period, make_unit_curves, rebuild_curve, reparametrize_exactly,
    ScalarLaw,
};
pub use reorder::{
    choose_ordered_contact, is_monotone_layout, open_closed_boundary_candidates,
    ordered_contact_score, ordered_curves, oriented_curve, oriented_parameter, reorder_network,
    seam_tangent_mismatch,
};

/// `GordonUtilities::Contact` (`pxx:20-27`).
#[derive(Clone, Copy, Debug)]
pub struct Contact {
    pub is_valid: bool,
    pub profile_param: f64,
    pub guide_param: f64,
    pub gap: f64,
}

impl Default for Contact {
    fn default() -> Self {
        Self {
            is_valid: false,
            profile_param: 0.0,
            guide_param: 0.0,
            gap: f64::MAX,
        }
    }
}

/// `GordonUtilities::ContactSet` (`pxx:29-33`).
#[derive(Clone, Debug, Default)]
pub struct ContactSet {
    pub candidates: Vec<Contact>,
}

/// `GordonUtilities::NetworkPreparation` (`pxx:1619-1725`).
pub struct NetworkPreparation<'a> {
    profiles: &'a mut Vec<GeomBSplineCurve>,
    guides: &'a mut Vec<GeomBSplineCurve>,
    profile_params: &'a mut Vec<Vec<f64>>,
    guide_params: &'a mut Vec<Vec<f64>>,
    contact_sets: Vec<ContactSet>,
    tolerance: f64,
    max_contact_gap: f64,
    max_reparametrization_deviation: f64,
    to_use_parallel: bool,
}

impl<'a> NetworkPreparation<'a> {
    /// `NetworkPreparation::NetworkPreparation` (`pxx:1624-1637`).
    pub fn new(
        profiles: &'a mut Vec<GeomBSplineCurve>,
        guides: &'a mut Vec<GeomBSplineCurve>,
        profile_params: &'a mut Vec<Vec<f64>>,
        guide_params: &'a mut Vec<Vec<f64>>,
        tolerance: f64,
        to_use_parallel: bool,
    ) -> Self {
        Self {
            profiles,
            guides,
            profile_params,
            guide_params,
            contact_sets: Vec::new(),
            tolerance,
            max_contact_gap: 0.0,
            max_reparametrization_deviation: 0.0,
            to_use_parallel,
        }
    }

    /// `ConvertCurvesToUnitInterval` (`pxx:1639-1643`).
    pub fn convert_curves_to_unit_interval(&mut self) -> bool {
        make_unit_curves(self.profiles) && make_unit_curves(self.guides)
    }

    /// `FillIntersectionParameters` (`pxx:1645-1657`).
    pub fn fill_intersection_parameters(&mut self) -> bool {
        let is_done = fill_contact_table(
            self.profiles,
            self.guides,
            self.profile_params,
            self.guide_params,
            &mut self.contact_sets,
            self.tolerance,
            self.to_use_parallel,
        );
        self.max_contact_gap = if is_done {
            max_nearest_contact_gap(&self.contact_sets)
        } else {
            0.0
        };
        is_done
    }

    /// `ReorderNetwork` (`pxx:1659-1668`).
    pub fn reorder_network(&mut self) -> bool {
        reorder_network(
            self.profiles,
            self.guides,
            self.profile_params,
            self.guide_params,
            &self.contact_sets,
            self.tolerance,
        )
    }

    /// `NormalizeIntersectionParameters` (`pxx:1670-1680`).
    pub fn normalize_intersection_parameters(
        &mut self,
        is_u_closed: &mut bool,
        is_v_closed: &mut bool,
    ) -> bool {
        normalize_contact_table(
            self.profiles,
            self.guides,
            self.profile_params,
            self.guide_params,
            self.tolerance,
            is_u_closed,
            is_v_closed,
        )
    }

    /// `EqualizeIntersectionParameters` (`pxx:1682-1695`).
    pub fn equalize_intersection_parameters(
        &mut self,
        allow_approximate_rational: bool,
        did_approximate: &mut bool,
    ) -> bool {
        self.max_reparametrization_deviation = 0.0;
        equalize_parameters(
            self.profiles,
            self.guides,
            self.profile_params,
            self.guide_params,
            allow_approximate_rational,
            did_approximate,
            &mut self.max_reparametrization_deviation,
        )
    }

    /// `MaxContactGap` (`pxx:1697-1698`).
    pub fn max_contact_gap(&self) -> f64 {
        self.max_contact_gap
    }

    /// `MaxReparametrizationDeviation` (`pxx:1700-1701`).
    pub fn max_reparametrization_deviation(&self) -> f64 {
        self.max_reparametrization_deviation
    }

    /// `IntersectionPoints` (`pxx:1703-1707`).
    pub fn intersection_points(&self) -> Vec<Vec<GpPnt>> {
        contact_point_grid(self.profiles, self.profile_params)
    }

    /// `IntersectionWeights` (`pxx:1709-1713`).
    pub fn intersection_weights(&self) -> Vec<Vec<f64>> {
        contact_weight_grid(self.profiles, self.profile_params)
    }
}
