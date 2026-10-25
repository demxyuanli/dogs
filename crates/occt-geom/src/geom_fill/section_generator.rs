//! `GeomFill_SectionGenerator` — the `GeomFill_Profiler` view `AppBlend_AppSurf`
//! evaluates section by section.
//!
//! Source: `TKGeomAlgo/GeomFill/GeomFill_SectionGenerator.hxx:38-74` and
//! `GeomFill_SectionGenerator.cxx:22-109` (OCCT 8.0.0).
//!
//! Translation notes:
//! - The OCCT class derives from `GeomFill_Profiler`; the port holds a
//!   [`Profiler`] and forwards [`SectionGenerator::add_curve`] /
//!   [`SectionGenerator::perform`] to it.
//! - `myParams` (`NCollection_HArray1<double>`) is a `Vec<f64>`. OCCT's
//!   constructor only fills it when the profile already has more than one
//!   curve, which cannot happen at construction time (the sequence is empty),
//!   so [`SectionGenerator::new`] starts with no parameters; they arrive
//!   through [`SectionGenerator::set_param`] (`cxx:38-47`).
//! - `WeightsArray()` (`Geom_BSplineCurve.hxx:805`) is the allocated
//!   weight array, i.e. all ones when the curve is not rational.

use crate::bspline_curve::GeomBSplineCurve;
use crate::curve::Curve;
use crate::geom_fill::profiler::Profiler;

/// `GeomFill_SectionGenerator` (`GeomFill_SectionGenerator.hxx:38-74`).
#[derive(Default)]
pub struct SectionGenerator {
    /// The `GeomFill_Profiler` base (`hxx:33`).
    profiler: Profiler,
    /// `myParams` (`hxx:71`).
    params: Vec<f64>,
}

impl SectionGenerator {
    /// `GeomFill_SectionGenerator::GeomFill_SectionGenerator` (`cxx:22-35`).
    pub fn new() -> Self {
        Self::default()
    }

    /// The `GeomFill_Profiler::AddCurve` base arm
    /// (`GeomFill_Profiler.cxx:125-162`).
    pub fn add_curve(&mut self, curve: &dyn Curve) -> Result<(), String> {
        self.profiler.add_curve(curve)
    }

    /// The `GeomFill_Profiler::Perform` base arm
    /// (`GeomFill_Profiler.cxx:166-246`).
    pub fn perform(&mut self, ptol: f64) -> Result<(), String> {
        self.profiler.perform(ptol)
    }

    /// The `GeomFill_Profiler` base object.
    pub fn profiler(&self) -> &Profiler {
        &self.profiler
    }

    /// `GeomFill_Profiler::IsPeriodic` inherited by
    /// `GeomFill_SectionGenerator` (`GeomFill_Profiler.hxx:56`); used as the
    /// `UPeriodic` flag of the approximated surface
    /// (`GeomFill_NSections.cxx:589`).
    pub fn is_periodic(&self) -> bool {
        self.profiler.is_periodic()
    }

    /// `GeomFill_SectionGenerator::SetParam` (`cxx:38-47`).
    pub fn set_param(&mut self, params: &[f64]) {
        self.params = params.to_vec();
    }

    /// `GeomFill_SectionGenerator::GetShape` (`cxx:50-60`):
    /// `(NbPoles, NbKnots, Degree, NbPoles2d)` of the first section.
    pub fn get_shape(&self) -> Result<(i32, i32, i32, i32), String> {
        let c = self
            .profiler
            .curves()
            .first()
            .ok_or_else(|| "Standard_OutOfRange: GeomFill_SectionGenerator::GetShape".to_string())?;
        Ok((c.nb_poles() as i32, c.nb_knots() as i32, c.degree() as i32, 0))
    }

    /// `GeomFill_SectionGenerator::Knots` (`cxx:62-67`).
    pub fn knots(&self) -> Result<Vec<f64>, String> {
        let c = self
            .profiler
            .curves()
            .first()
            .ok_or_else(|| "Standard_OutOfRange: GeomFill_SectionGenerator::Knots".to_string())?;
        Ok(c.distinct_knots_and_mults().0)
    }

    /// `GeomFill_SectionGenerator::Mults` (`cxx:70-75`).
    pub fn mults(&self) -> Result<Vec<i32>, String> {
        let c = self
            .profiler
            .curves()
            .first()
            .ok_or_else(|| "Standard_OutOfRange: GeomFill_SectionGenerator::Mults".to_string())?;
        Ok(c.distinct_knots_and_mults().1)
    }

    /// `GeomFill_SectionGenerator::Section(P, Poles, Poles2d, Weigths)`
    /// (`cxx:92-102`). `P` is 1-based, as in OCCT.
    pub fn section(&self, p: i32) -> Result<(Vec<occt_core::gp::GpPnt>, Vec<f64>), String> {
        let curve = self.section_curve(p)?;
        Ok((curve.poles.to_vec(), weights_array(curve)))
    }

    /// `GeomFill_SectionGenerator::Section(P, Poles, DPoles, Poles2d,
    /// DPoles2d, Weigths, DWeigths)` (`cxx:77-90`): the base section plus
    /// `false`, because the derivatives are not computed.
    pub fn section_with_derivatives(
        &self,
        p: i32,
    ) -> Result<(Vec<occt_core::gp::GpPnt>, Vec<f64>), String> {
        self.section(p)
    }

    /// `GeomFill_SectionGenerator::Parameter` (`cxx:104-107`).
    pub fn parameter(&self, p: i32) -> Result<f64, String> {
        self.params
            .get((p - 1).max(0) as usize)
            .copied()
            .ok_or_else(|| "Standard_OutOfRange: GeomFill_SectionGenerator::Parameter".to_string())
    }

    fn section_curve(&self, p: i32) -> Result<&GeomBSplineCurve, String> {
        if p < 1 {
            return Err("Standard_OutOfRange: GeomFill_SectionGenerator::Section".into());
        }
        self.profiler
            .curves()
            .get((p - 1) as usize)
            .ok_or_else(|| "Standard_OutOfRange: GeomFill_SectionGenerator::Section".to_string())
    }
}

/// `Geom_BSplineCurve::WeightsArray` (`Geom_BSplineCurve.hxx:805`): the
/// stored weight array, all ones when the curve is not rational.
fn weights_array(curve: &GeomBSplineCurve) -> Vec<f64> {
    match &curve.weights {
        Some(w) => w.clone(),
        None => vec![1.0; curve.poles.len()],
    }
}
