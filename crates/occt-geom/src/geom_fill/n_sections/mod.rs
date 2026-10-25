//! `GeomFill_NSections` -- a section law defined by N sections.
//!
//! Source: `TKGeomAlgo/GeomFill/GeomFill_NSections.hxx:38-219` and
//! `GeomFill_NSections.cxx:201-1047` (OCCT 8.0.0).
//!
//! Translation notes:
//! - `NCollection_Sequence<Handle(Geom_Curve)>` is `Vec<Arc<dyn Curve>>`,
//!   `NCollection_Sequence<gp_Trsf>` is `Vec<GpTrsf>`,
//!   `NCollection_Sequence<double>` is `Vec<f64>` and
//!   `Handle(Geom_BSplineSurface)` is `Option<Arc<GeomBSplineSurface>>`.
//! - The five constructors (`cxx:201-286`) all end in `ComputeSurface()`. OCCT
//!   leaves a half-built law behind when the profiler or the segmentation
//!   raises; the port returns `Result` so the same failure is reported.
//! - `NCollection_Array1` out-parameters (the `Poles` / `Weights` families)
//!   become `Vec`s, as in [`SectionLaw`].
//! - `D0` / `D1` / `D2` return OCCT's `bool`, mapped to `Ok` / `Err`; the
//!   "null surface" (`cxx:293-296`), "V degree < 2" (`cxx:410-413`) and
//!   "null weight" (`cxx:378-380`, `cxx:466-468`) answers are errors here.
//! - `ComputeSurface` without a reference surface ends in
//!   `GeomFill_AppSurf::Perform(line, section, true)` (`cxx:578-579`), i.e. the
//!   `SpApprox` branch of `AppBlend_AppSurf` (`AppBlend_AppSurf.gxx:420-471`,
//!   `AppDef_BSplineCompute` / `Approx_BSplComputeLine.gxx:1-1458`), which is
//!   ported; the `myRefSurf` branch (`cxx:592-623`) is ported in full.

mod derivatives;
mod law_queries;
mod section_curves;

use std::sync::Arc;

use occt_core::convert::ParameterisationType;
use occt_core::gp::{GpPnt, GpTrsf, GpVec};
use occt_core::kernel::geomabs::Shape;
use occt_core::precision::PCONFUSION;

use crate::bspline_surface::GeomBSplineSurface;
use crate::convert_bspl::curve_to_bspline_curve;
use crate::curve::Curve;
use crate::geom_fill::app_surf::AppSurf;
use crate::geom_fill::line::Line;
use crate::geom_fill::section_generator::SectionGenerator;
use crate::geom_fill::section_law::SectionLaw;

/// `GeomFill_NSections` (`GeomFill_NSections.hxx:36-219`).
pub struct NSections {
    /// `UFirst` (`hxx:208`).
    u_first: f64,
    /// `ULast` (`hxx:209`).
    u_last: f64,
    /// `VFirst` (`hxx:210`).
    v_first: f64,
    /// `VLast` (`hxx:211`).
    v_last: f64,
    /// `mySections` (`hxx:212`).
    sections: Vec<Arc<dyn Curve>>,
    /// `myTrsfs` (`hxx:213`).
    trsfs: Vec<GpTrsf>,
    /// `myParams` (`hxx:214`).
    params: Vec<f64>,
    /// `mySurface` (`hxx:215`).
    surface: Option<Arc<GeomBSplineSurface>>,
    /// `myRefSurf` (`hxx:216`).
    ref_surf: Option<Arc<GeomBSplineSurface>>,
}

impl NSections {
    /// `GeomFill_NSections(NC)` (`cxx:203-210`).
    pub fn new(sections: Vec<Arc<dyn Curve>>) -> Result<Self, String> {
        Self::with_bounds(sections, Vec::new(), 0.0, 1.0, 0.0, 1.0)
    }

    /// `GeomFill_NSections(NC, NP)` (`cxx:216-226`).
    pub fn with_params(
        sections: Vec<Arc<dyn Curve>>,
        params: Vec<f64>,
    ) -> Result<Self, String> {
        Self::with_bounds(sections, params, 0.0, 1.0, 0.0, 1.0)
    }

    /// `GeomFill_NSections(NC, NP, UF, UL)` (`cxx:232-243`): `VFirst = 0`,
    /// `VLast = 1`.
    pub fn with_parametric_bounds(
        sections: Vec<Arc<dyn Curve>>,
        params: Vec<f64>,
        u_first: f64,
        u_last: f64,
    ) -> Result<Self, String> {
        Self::with_bounds(sections, params, u_first, u_last, 0.0, 1.0)
    }

    /// `GeomFill_NSections(NC, NP, UF, UL, VF, VL)` (`cxx:247-263`).
    pub fn with_bounds(
        sections: Vec<Arc<dyn Curve>>,
        params: Vec<f64>,
        u_first: f64,
        u_last: f64,
        v_first: f64,
        v_last: f64,
    ) -> Result<Self, String> {
        let mut law = Self {
            u_first,
            u_last,
            v_first,
            v_last,
            sections,
            trsfs: Vec::new(),
            params,
            surface: None,
            ref_surf: None,
        };
        law.compute_surface()?;
        Ok(law)
    }

    /// `GeomFill_NSections(NC, Trsfs, NP, UF, UL, VF, VL, Surf)`
    /// (`cxx:267-286`).
    #[allow(clippy::too_many_arguments)]
    pub fn with_reference_surface(
        sections: Vec<Arc<dyn Curve>>,
        trsfs: Vec<GpTrsf>,
        params: Vec<f64>,
        u_first: f64,
        u_last: f64,
        v_first: f64,
        v_last: f64,
        reference_surface: Arc<GeomBSplineSurface>,
    ) -> Result<Self, String> {
        let mut law = Self {
            u_first,
            u_last,
            v_first,
            v_last,
            sections,
            trsfs,
            params,
            surface: None,
            ref_surf: Some(reference_surface),
        };
        law.compute_surface()?;
        Ok(law)
    }

    /// `GeomFill_NSections::SetSurface(RefSurf)` (`cxx:492-496`): only stores
    /// the reference surface, `mySurface` is not recomputed.
    pub fn set_surface(&mut self, reference_surface: Arc<GeomBSplineSurface>) {
        self.ref_surf = Some(reference_surface);
    }

    /// `mySection` accessor: the stored input sections (`hxx:212`).
    pub fn sections(&self) -> &[Arc<dyn Curve>] {
        &self.sections
    }

    /// `myParams` accessor (`hxx:214`).
    pub fn params(&self) -> &[f64] {
        &self.params
    }

    /// `GeomFill_NSections::ComputeSurface` (`cxx:500-636`).
    pub fn compute_surface(&mut self) -> Result<(), String> {
        let computed = match self.ref_surf.clone() {
            None => self.approximate_surface()?,
            Some(reference) => {
                Some(Self::segment_reference_surface(&reference, self.u_first, self.u_last)?)
            }
        };
        let Some(mut surface) = computed else {
            self.surface = None;
            return Ok(());
        };
        // `cxx:625-628`: "On augmente le degre pour que le positionnement D2
        // soit correct".
        if surface.deg_v() < 2 {
            surface.increase_degree(surface.deg_u(), 2)?;
        }
        self.surface = Some(Arc::new(surface));
        Ok(())
    }

    /// The `myRefSurf.IsNull()` arm of `ComputeSurface` (`cxx:505-590`).
    fn approximate_surface(&self) -> Result<Option<GeomBSplineSurface>, String> {
        let nb_curves = self.sections.len();
        // `cxx:511-514`: a single curve cannot build a surface.
        if nb_curves <= 1 {
            return Ok(None);
        }

        let mut section = SectionGenerator::new();
        for curv in &self.sections {
            // `cxx:526-531`: keep a B-spline section as it is, otherwise
            // `GeomConvert::CurveToBSplineCurve(curv, Convert_QuasiAngular)`.
            // `curve_to_bspline_curve` copies a B-spline input unchanged and
            // converts every other type, so it covers both arms.
            let mut curv_bs = curve_to_bspline_curve(curv.as_ref(), ParameterisationType::QuasiAngular)
                .map_err(|_| "GeomFill_NSections: GeomConvert::CurveToBSplineCurve".to_string())?;
            // `cxx:533-535`: `BSplCLib::Reparametrize(UFirst, ULast, BSK)`
            // followed by `curvBS->SetKnots(BSK)`.
            let (mut knots, _) = curv_bs.distinct_knots_and_mults();
            occt_core::bspl::knots::reparametrize(self.u_first, self.u_last, &mut knots);
            curv_bs
                .set_knots(&knots)
                .map_err(|e| format!("Geom_BSplineCurve::SetKnots: {e}"))?;
            section.add_curve(&curv_bs)?;
        }

        // `cxx:561-571`: the parameters, when the law carries them.
        if !self.params.is_empty() {
            section.set_param(&self.params);
        }
        section.perform(PCONFUSION)?;

        let line = Line::with_nb_points(nb_curves as i32);
        // `cxx:573-579`: `GeomFill_AppSurf(2, 6, 1e-6, 1e-6, 0, knownP)` with
        // `Continuity(C1)` and `Perform(line, section, true)`.
        let mut approx = AppSurf::with_settings(2, 6, 1.0e-6, 1.0e-6, 0, !self.params.is_empty());
        approx.set_continuity(Shape::C1);
        approx.perform(&line, &section, true)?;

        if !approx.is_done() {
            // OCCT reads `SurfPoles()` on a not-done `AppBlend_AppSurf` and
            // raises `StdFail_NotDone` (`AppBlend_AppSurf.gxx:1073-1090`).
            return Ok(None);
        }
        Ok(Self::surface_from_approximation(&approx, section.is_periodic()))
    }

    /// The `new Geom_BSplineSurface(anApprox.SurfPoles(), ...)` of
    /// `cxx:581-589`.
    ///
    /// The 10-argument `Geom_BSplineSurface` constructor used there
    /// (`Geom_BSplineSurface.hxx:213-222`) takes `UPeriodic` and leaves
    /// `VPeriodic` at its `false` default.
    fn surface_from_approximation(
        approx: &AppSurf,
        u_periodic: bool,
    ) -> Option<GeomBSplineSurface> {
        let deg_u = approx.u_degree() as usize;
        let deg_v = approx.v_degree() as usize;
        let poles = approx.surf_poles().to_vec();
        let weights = approx.surf_weights().to_vec();
        let u_knots = approx.surf_u_knots().to_vec();
        let v_knots = approx.surf_v_knots().to_vec();
        let u_mults = approx.surf_u_mults().to_vec();
        let v_mults = approx.surf_v_mults().to_vec();

        let flat_u = if u_periodic {
            occt_core::bspl::knots::knot_sequence_periodic(&u_knots, &u_mults, deg_u as i32)
        } else {
            occt_core::bspl::banded_interp::knot_sequence(&u_knots, &u_mults, deg_u as i32)
        };
        let flat_v = occt_core::bspl::banded_interp::knot_sequence(&v_knots, &v_mults, deg_v as i32);
        // `Geom_BSplineSurface::Rational(Weights, myURational, myVRational)`
        // (`Geom_BSplineSurface.cxx:110-138`): identical weights leave the
        // surface non-rational.
        let weights = if rational_weights(&weights) {
            Some(weights)
        } else {
            None
        };
        Some(GeomBSplineSurface::from_parts(
            poles, weights, flat_u, flat_v, deg_u, deg_v, u_periodic, false,
        ))
    }

    /// The `myRefSurf` arm of `ComputeSurface` (`cxx:592-623`): snap the
    /// requested U bounds onto reference knots and cut the reference surface.
    fn segment_reference_surface(
        reference: &GeomBSplineSurface,
        u_first: f64,
        u_last: f64,
    ) -> Result<GeomBSplineSurface, String> {
        let mut surface = reference.clone();
        let u_knots = reference.uknots();
        let mut ui1 = u_first;
        let mut ui2 = u_last;

        let (i1, i2) = reference.locate_u(ui1, PCONFUSION);
        if let Some(&k) = knot_at(u_knots, i1) {
            if (ui1 - k).abs() <= PCONFUSION {
                ui1 = k;
            }
        }
        if let Some(&k) = knot_at(u_knots, i2) {
            if (ui1 - k).abs() <= PCONFUSION {
                ui1 = k;
            }
        }
        let (i1, i2) = reference.locate_u(ui2, PCONFUSION);
        if let Some(&k) = knot_at(u_knots, i1) {
            if (ui2 - k).abs() <= PCONFUSION {
                ui2 = k;
            }
        }
        if let Some(&k) = knot_at(u_knots, i2) {
            if (ui2 - k).abs() <= PCONFUSION {
                ui2 = k;
            }
        }

        // `cxx:619-620`: `VKnot(FirstVKnotIndex())` / `VKnot(LastVKnotIndex())`.
        let v_knots = reference.vknots();
        let v0 = *v_knots
            .first()
            .ok_or_else(|| "GeomFill_NSections: empty V knot sequence".to_string())?;
        let v1 = *v_knots
            .last()
            .ok_or_else(|| "GeomFill_NSections: empty V knot sequence".to_string())?;
        surface.check_and_segment(ui1, ui2, v0, v1)?;
        Ok(surface)
    }
}

/// `myUKnots.Value(index)` for the 1-based `index` returned by `LocateU`
/// (`Geom_BSplineSurface_1.cxx:1464-1514`).
fn knot_at(knots: &[f64], index: i32) -> Option<&f64> {
    if index < 1 {
        return None;
    }
    knots.get((index - 1) as usize)
}

/// `Geom_BSplineSurface::Rational(Weights, URational, VRational)`
/// (`Geom_BSplineSurface.cxx:110-138`).
fn rational_weights(weights: &[Vec<f64>]) -> bool {
    let n_v = weights.first().map_or(0, |row| row.len());
    // Vrational: consecutive U indices inside one column.
    for j in 0..n_v {
        for i in 0..weights.len().saturating_sub(1) {
            let a = weights[i][j];
            if (a - weights[i + 1][j]).abs() > occt_core::precision::epsilon(a.abs()) {
                return true;
            }
        }
    }
    // Urational: consecutive V indices inside one row.
    for row in weights {
        for j in 0..row.len().saturating_sub(1) {
            let a = row[j];
            if (a - row[j + 1]).abs() > occt_core::precision::epsilon(a.abs()) {
                return true;
            }
        }
    }
    false
}

impl SectionLaw for NSections {
    /// `GeomFill_NSections::D0` (`cxx:288-311`).
    fn d0(&self, param: f64) -> Result<(Vec<GpPnt>, Vec<f64>), String> {
        self.poles_at(param)
    }

    /// `GeomFill_NSections::D1` (`cxx:315-390`).
    fn d1(&self, param: f64) -> Result<(Vec<GpPnt>, Vec<GpVec>, Vec<f64>, Vec<f64>), String> {
        self.poles_and_derivative_at(param)
    }

    /// `GeomFill_NSections::D2` (`cxx:394-480`).
    fn d2(
        &self,
        param: f64,
    ) -> Result<
        (
            Vec<GpPnt>,
            Vec<GpVec>,
            Vec<GpVec>,
            Vec<f64>,
            Vec<f64>,
            Vec<f64>,
        ),
        String,
    > {
        self.poles_and_second_derivative_at(param)
    }

    /// `GeomFill_NSections::BSplineSurface` (`cxx:484-488`).
    fn bspline_surface(&self) -> Option<Arc<GeomBSplineSurface>> {
        self.surface.clone()
    }

    /// `GeomFill_NSections::SectionShape` (`cxx:640-650`).
    fn section_shape(&self) -> (i32, i32, i32) {
        self.shape_values()
    }

    /// `GeomFill_NSections::Knots` (`cxx:654-661`).
    fn knots(&self) -> Vec<f64> {
        self.knot_values()
    }

    /// `GeomFill_NSections::Mults` (`cxx:665-672`).
    fn mults(&self) -> Vec<i32> {
        self.mult_values()
    }

    /// `GeomFill_NSections::IsRational` (`cxx:676-684`).
    fn is_rational(&self) -> bool {
        self.is_rational_flag()
    }

    /// `GeomFill_NSections::IsUPeriodic` (`cxx:690-697`).
    fn is_u_periodic(&self) -> bool {
        self.is_u_periodic_flag()
    }

    /// `GeomFill_NSections::IsVPeriodic` (`cxx:703-711`).
    fn is_v_periodic(&self) -> bool {
        self.is_v_periodic_flag()
    }

    /// `GeomFill_NSections::NbIntervals` (`cxx:715-723`).
    fn nb_intervals(&self, s: u8) -> i32 {
        self.v_interval_count(s)
    }

    /// `GeomFill_NSections::Intervals` (`cxx:730-737`).
    fn intervals(&self, s: u8) -> Vec<f64> {
        self.v_interval_bounds(s)
    }

    /// `GeomFill_NSections::SetInterval` (`cxx:746-750`): "nothing to do,
    /// mySurface is assumed Cn in V".
    fn set_interval(&self, _first: f64, _last: f64) {}

    /// `GeomFill_NSections::GetInterval` (`cxx:756-759`).
    fn get_interval(&self) -> (f64, f64) {
        (self.v_first, self.v_last)
    }

    /// `GeomFill_NSections::GetDomain` (`cxx:764-767`).
    fn get_domain(&self) -> (f64, f64) {
        (self.v_first, self.v_last)
    }

    /// `GeomFill_NSections::GetTolerance` (`cxx:773-785`).
    fn get_tolerance(&self, bound_tol: f64, surf_tol: f64, _angle_tol: f64) -> Vec<f64> {
        self.tolerance_values(bound_tol, surf_tol)
    }

    /// `GeomFill_NSections::BarycentreOfSurf` (`cxx:791-815`).
    fn barycentre_of_surf(&self) -> Result<GpPnt, String> {
        self.barycentre()
    }

    /// `GeomFill_NSections::MaximalSection` (`cxx:819-832`).
    fn maximal_section(&self) -> f64 {
        self.maximal_section_length()
    }

    /// `GeomFill_NSections::GetMinimalWeight` (`cxx:837-872`).
    fn get_minimal_weight(&self) -> Result<Vec<f64>, String> {
        self.minimal_weights()
    }

    /// `GeomFill_NSections::IsConstant` (`cxx:876-938`).
    fn is_constant(&self, error: &mut f64) -> bool {
        *error = 0.0;
        self.constant_section_test()
    }

    /// `GeomFill_NSections::ConstantSection` (`cxx:942-950`).
    fn constant_section(&self) -> Result<Arc<dyn Curve>, String> {
        self.constant_section_copy()
    }

    /// `GeomFill_NSections::IsConicalLaw` (`cxx:954-1011`).
    fn is_conical_law(&self, error: &mut f64) -> bool {
        *error = 0.0;
        self.conical_law_test()
    }

    /// `GeomFill_NSections::CirclSection` (`cxx:1016-1047`).
    fn circl_section(&self, param: f64) -> Result<Arc<dyn Curve>, String> {
        self.circle_section_at(param)
    }
}
