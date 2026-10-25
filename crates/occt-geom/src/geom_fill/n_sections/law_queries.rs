//! The read-only `SectionLaw` queries of `GeomFill_NSections`.
//!
//! Source: `GeomFill_NSections.cxx:640-872` (OCCT 8.0.0): `SectionShape`,
//! `Knots`, `Mults`, `IsRational`, `IsUPeriodic`, `IsVPeriodic`, `NbIntervals`,
//! `Intervals`, `GetTolerance`, `BarycentreOfSurf`, `MaximalSection` and
//! `GetMinimalWeight`.
//!
//! Translation notes:
//! - `NbIntervals` / `Intervals` build a `GeomAdaptor_Surface` on `mySurface`
//!   and ask for `NbVIntervals` / `VIntervals`; the port reaches the same
//!   methods through the [`Surface`] trait `GeomBSplineSurface` implements.
//! - OCCT leaves the caller's `NCollection_Array1` untouched when `mySurface`
//!   is null (`cxx:643-646`, `:658-661`, `:669-672`); the port answers the
//!   neutral value (`(0, 0, 0)`, empty vector, 0 intervals).
//! - `BarycentreOfSurf` samples 21 x 21 points but does **not** reset `U`
//!   between the `V` rows (`cxx:806-812`); the port keeps that control flow.

use occt_core::gp::GpPnt;

use crate::geom_fill::n_sections::NSections;
use crate::gcpnts::curve_length;
use crate::surface::Surface;

impl NSections {
    /// `GeomFill_NSections::SectionShape` (`cxx:640-650`).
    pub(super) fn shape_values(&self) -> (i32, i32, i32) {
        match &self.surface {
            None => (0, 0, 0),
            Some(surface) => (
                surface.nb_poles_u() as i32,
                surface.distinct_knots_and_mults_u().0.len() as i32,
                surface.deg_u() as i32,
            ),
        }
    }

    /// `GeomFill_NSections::Knots` (`cxx:654-661`).
    pub(super) fn knot_values(&self) -> Vec<f64> {
        match &self.surface {
            None => Vec::new(),
            Some(surface) => surface.distinct_knots_and_mults_u().0,
        }
    }

    /// `GeomFill_NSections::Mults` (`cxx:665-672`).
    pub(super) fn mult_values(&self) -> Vec<i32> {
        match &self.surface {
            None => Vec::new(),
            Some(surface) => surface.distinct_knots_and_mults_u().1,
        }
    }

    /// `GeomFill_NSections::IsRational` (`cxx:676-684`): `IsURational()`.
    pub(super) fn is_rational_flag(&self) -> bool {
        match &self.surface {
            None => false,
            Some(surface) => surface.is_u_rational(),
        }
    }

    /// `GeomFill_NSections::IsUPeriodic` (`cxx:690-697`).
    pub(super) fn is_u_periodic_flag(&self) -> bool {
        match &self.surface {
            None => false,
            Some(surface) => surface.u_periodic(),
        }
    }

    /// `GeomFill_NSections::IsVPeriodic` (`cxx:703-711`).
    pub(super) fn is_v_periodic_flag(&self) -> bool {
        match &self.surface {
            None => false,
            Some(surface) => surface.v_periodic(),
        }
    }

    /// `GeomFill_NSections::NbIntervals` (`cxx:715-723`).
    pub(super) fn v_interval_count(&self, shape: u8) -> i32 {
        match &self.surface {
            None => 0,
            Some(surface) => surface.nb_v_intervals(shape),
        }
    }

    /// `GeomFill_NSections::Intervals` (`cxx:730-737`).
    pub(super) fn v_interval_bounds(&self, shape: u8) -> Vec<f64> {
        match &self.surface {
            None => Vec::new(),
            Some(surface) => surface.v_intervals(shape),
        }
    }

    /// `GeomFill_NSections::GetTolerance` (`cxx:773-785`).
    pub(super) fn tolerance_values(&self, bound_tol: f64, surf_tol: f64) -> Vec<f64> {
        let n = match &self.surface {
            None => 0,
            Some(surface) => surface.nb_poles_u(),
        };
        let mut tol3d = vec![surf_tol; n];
        if bound_tol < surf_tol && n > 0 {
            tol3d[0] = bound_tol;
            tol3d[n - 1] = bound_tol;
        }
        tol3d
    }

    /// `GeomFill_NSections::BarycentreOfSurf` (`cxx:791-815`).
    pub(super) fn barycentre(&self) -> Result<GpPnt, String> {
        let Some(surface) = self.surface.as_ref() else {
            return Ok(GpPnt::zero());
        };
        let (u0, u1) = surface.u_range();
        let (v0, v1) = surface.v_range();
        let delta_v = (v1 - v0) / 20.0;
        let delta_u = (u1 - u0) / 20.0;
        let mut sum = occt_core::gp::GpXyz::zero();
        let mut v = v0;
        let mut u = u0;
        for _j in 0..=20 {
            for _i in 0..=20 {
                let p = surface.d0(u, v);
                sum = sum.added(p.xyz());
                u += delta_u;
            }
            v += delta_v;
        }
        let bary = sum.divided((21 * 21) as f64);
        Ok(GpPnt::from_xyz(&bary))
    }

    /// `GeomFill_NSections::MaximalSection` (`cxx:819-832`).
    pub(super) fn maximal_section_length(&self) -> f64 {
        let mut max = 0.0;
        for section in &self.sections {
            let length = curve_length(section.as_ref());
            if length > max {
                max = length;
            }
        }
        max
    }

    /// `GeomFill_NSections::GetMinimalWeight` (`cxx:837-872`).
    pub(super) fn minimal_weights(&self) -> Result<Vec<f64>, String> {
        let Some(surface) = self.surface.as_ref() else {
            return Ok(Vec::new());
        };
        let n_u = surface.nb_poles_u();
        // `cxx:849`: `mySurface->IsURational()`.
        if !surface.is_u_rational() {
            return Ok(vec![1.0; n_u]);
        }
        let n_v = surface.nb_poles_v();
        let mut weights = Vec::with_capacity(n_u);
        for i in 0..n_u {
            let mut min = surface.weight(i, 0);
            for j in 1..n_v {
                let w = surface.weight(i, j);
                if min > w {
                    min = w;
                }
            }
            weights.push(min);
        }
        Ok(weights)
    }
}
