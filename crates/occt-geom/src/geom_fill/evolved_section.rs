//! `GeomFill_EvolvedSection` — a section law that scales one section curve by
//! a `Law_Function` along the sweep.
//!
//! Source: `TKGeomAlgo/GeomFill/GeomFill_EvolvedSection.hxx:35-166` and
//! `GeomFill_EvolvedSection.cxx:33-347` (OCCT 8.0.0): constructor
//! `cxx:33-51`, `D0` `cxx:55-72`, `D1` `cxx:74-97`, `D2` `cxx:99-129`,
//! `BSplineSurface` `cxx:131-162`, `SectionShape` `cxx:164-171`, `Knots`
//! `cxx:173-177`, `Mults` `cxx:179-186`, `IsRational` `cxx:188-191`,
//! `IsUPeriodic` `cxx:193-197`, `IsVPeriodic` `cxx:199-208`, `NbIntervals`
//! `cxx:210-215`, `Intervals` `cxx:217-224`, `SetInterval` `cxx:226-231`,
//! `GetInterval` `cxx:233-241`, `GetDomain` `cxx:243-249`, `GetTolerance`
//! `cxx:251-267`, `BarycentreOfSurf` `cxx:269-292`, `MaximalSection`
//! `cxx:294-312`, `GetMinimalWeight` `cxx:314-318`, `IsConstant` `cxx:320-326`,
//! `ConstantSection` `cxx:328-347`.
//!
//! Translation notes:
//! - `mySection = C->Copy()` (`cxx:38`) is the caller's `Arc`, exactly as the
//!   other section laws of the port keep it.
//! - `TLaw` is rewritten by `SetInterval` (`cxx:230`), so it lives behind a
//!   `Mutex` and the law takes `&self`.
//! - The `Geom_BSplineCurve` down_cast (`cxx:41-50`) is
//!   [`is_bspline`] / [`bspline_copy`]; the not-a-B-spline arm converts with
//!   `Convert_QuasiAngular` and removes the wrap-around knot of a periodic
//!   result, like `GeomFill_UniformSection`.
//! - `BSplineSurface` returns a null handle (`cxx:157-161`; the construction is
//!   commented out in OCCT) and hence `None`.
//! - `ConstantSection` throws `StdFail_NotDone` (`cxx:332-334`) and every
//!   `Law_Function` failure is an OCCT throw: both become `Err` / panics at the
//!   same call sites as the OCCT exceptions.

use std::sync::{Arc, Mutex};

use occt_core::gp::{GpPnt, GpVec};
use occt_core::kernel::geomabs::Shape;
use occt_core::precision::{CONFUSION, PCONFUSION};

use crate::bspline_curve::GeomBSplineCurve;
use crate::bspline_surface::GeomBSplineSurface;
use crate::convert_bspl::{bspline_copy, curve_to_bspline_curve, is_bspline};
use crate::curve::Curve;
use crate::gcpnts::curve_length;
use crate::geom_fill::section_law::SectionLaw;
use crate::law::LawFunction;
use occt_core::convert::ParameterisationType;

/// `GeomFill_EvolvedSection` (`GeomFill_EvolvedSection.hxx:35-166`).
pub struct EvolvedSection {
    /// `First` (`hxx:158`): the law bound.
    first: f64,
    /// `Last` (`hxx:159`): the law bound.
    last: f64,
    /// `mySection` (`hxx:160`): the copied section curve.
    section: Arc<dyn Curve>,
    /// `myLaw` (`hxx:161`): the trimmed input law.
    my_law: Arc<dyn LawFunction>,
    /// `TLaw` (`hxx:162`): `myLaw` trimmed by `SetInterval`.
    t_law: Mutex<Arc<dyn LawFunction>>,
    /// `myCurve` (`hxx:163`): the B-spline representation of the section.
    curve: GeomBSplineCurve,
}

impl EvolvedSection {
    /// `GeomFill_EvolvedSection(C, L)` (`cxx:33-51`).
    pub fn new(c: Arc<dyn Curve>, l: Arc<dyn LawFunction>) -> Result<Self, String> {
        let (first, last) = l.bounds();
        let law = l
            .trim(first, last, 1.0e-20)
            .map_err(|_| "Standard_DomainError: GeomFill_EvolvedSection".to_string())?;
        let curve = if is_bspline(c.as_ref()) {
            bspline_copy(c.as_ref())
                .map_err(|_| "Standard_DomainError: GeomFill_EvolvedSection".to_string())?
        } else {
            let mut curve = curve_to_bspline_curve(c.as_ref(), ParameterisationType::QuasiAngular)
                .map_err(|_| "Standard_DomainError: GeomFill_EvolvedSection".to_string())?;
            if curve.is_periodic() {
                let m = curve.degree() as i32 / 2 + 1;
                curve.remove_knot(1, m, CONFUSION).map_err(|_| {
                    "Standard_ConstructionError: GeomFill_EvolvedSection".to_string()
                })?;
            }
            curve
        };
        Ok(Self {
            first,
            last,
            section: c,
            t_law: Mutex::new(law.clone()),
            my_law: law,
            curve,
        })
    }

    fn trimmed_law(&self) -> Arc<dyn LawFunction> {
        self.t_law
            .lock()
            .unwrap_or_else(|e| e.into_inner())
            .clone()
    }

    /// `Geom_BSplineCurve::WeightsArray()` (`Geom_BSplineCurve.hxx:805`).
    fn weights_array(&self) -> Vec<f64> {
        match &self.curve.weights {
            Some(w) => w.clone(),
            None => vec![1.0; self.curve.poles.len()],
        }
    }
}

impl SectionLaw for EvolvedSection {
    /// `GeomFill_EvolvedSection::D0` (`cxx:55-72`): every pole scaled by
    /// `TLaw->Value(U)`.
    fn d0(&self, u: f64) -> Result<(Vec<GpPnt>, Vec<f64>), String> {
        let val = self.trimmed_law().value(u);
        let poles = self
            .curve
            .poles
            .iter()
            .map(|p| GpPnt::new(p.x() * val, p.y() * val, p.z() * val))
            .collect();
        Ok((poles, self.weights_array()))
    }

    /// `GeomFill_EvolvedSection::D1` (`cxx:74-97`).
    fn d1(&self, u: f64) -> Result<(Vec<GpPnt>, Vec<GpVec>, Vec<f64>, Vec<f64>), String> {
        let (val, dval) = self.trimmed_law().d1(u);
        let dpoles = self
            .curve
            .poles
            .iter()
            .map(|p| GpVec::new(p.x() * dval, p.y() * dval, p.z() * dval))
            .collect();
        let poles = self
            .curve
            .poles
            .iter()
            .map(|p| GpPnt::new(p.x() * val, p.y() * val, p.z() * val))
            .collect();
        let n = self.curve.poles.len();
        Ok((poles, dpoles, self.weights_array(), vec![0.0; n]))
    }

    /// `GeomFill_EvolvedSection::D2` (`cxx:99-129`).
    fn d2(
        &self,
        u: f64,
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
        let (val, dval, d2val) = self.trimmed_law().d2(u);
        let dpoles = self
            .curve
            .poles
            .iter()
            .map(|p| GpVec::new(p.x() * dval, p.y() * dval, p.z() * dval))
            .collect();
        let d2poles = self
            .curve
            .poles
            .iter()
            .map(|p| GpVec::new(p.x() * d2val, p.y() * d2val, p.z() * d2val))
            .collect();
        let poles = self
            .curve
            .poles
            .iter()
            .map(|p| GpPnt::new(p.x() * val, p.y() * val, p.z() * val))
            .collect();
        let n = self.curve.poles.len();
        Ok((
            poles,
            dpoles,
            d2poles,
            self.weights_array(),
            vec![0.0; n],
            vec![0.0; n],
        ))
    }

    /// `GeomFill_EvolvedSection::BSplineSurface` (`cxx:131-162`): null.
    fn bspline_surface(&self) -> Option<Arc<GeomBSplineSurface>> {
        None
    }

    /// `GeomFill_EvolvedSection::SectionShape` (`cxx:164-171`).
    fn section_shape(&self) -> (i32, i32, i32) {
        (
            self.curve.nb_poles() as i32,
            self.curve.distinct_knots_and_mults().0.len() as i32,
            self.curve.degree() as i32,
        )
    }

    /// `GeomFill_EvolvedSection::Knots` (`cxx:173-177`).
    fn knots(&self) -> Vec<f64> {
        self.curve.distinct_knots_and_mults().0
    }

    /// `GeomFill_EvolvedSection::Mults` (`cxx:179-186`).
    fn mults(&self) -> Vec<i32> {
        self.curve.distinct_knots_and_mults().1
    }

    /// `GeomFill_EvolvedSection::IsRational` (`cxx:188-191`).
    fn is_rational(&self) -> bool {
        self.curve.is_rational()
    }

    /// `GeomFill_EvolvedSection::IsUPeriodic` (`cxx:193-197`).
    fn is_u_periodic(&self) -> bool {
        self.curve.is_periodic()
    }

    /// `GeomFill_EvolvedSection::IsVPeriodic` (`cxx:199-208`):
    /// `|myLaw->Value(First) - myLaw->Value(Last)| < Precision::Confusion()`.
    fn is_v_periodic(&self) -> bool {
        (self.my_law.value(self.first) - self.my_law.value(self.last)).abs() < CONFUSION
    }

    /// `GeomFill_EvolvedSection::NbIntervals` (`cxx:210-215`).
    fn nb_intervals(&self, s: u8) -> i32 {
        self.my_law
            .nb_intervals(shape_of(s))
            .expect("Law_Function::NbIntervals")
    }

    /// `GeomFill_EvolvedSection::Intervals` (`cxx:217-224`).
    fn intervals(&self, s: u8) -> Vec<f64> {
        self.my_law
            .intervals(shape_of(s))
            .expect("Law_Function::Intervals")
    }

    /// `GeomFill_EvolvedSection::SetInterval` (`cxx:226-231`):
    /// `TLaw = myLaw->Trim(F, L, Precision::PConfusion())`.
    fn set_interval(&self, first: f64, last: f64) {
        let trimmed = self
            .my_law
            .trim(first, last, PCONFUSION)
            .expect("Law_Function::Trim");
        *self.t_law.lock().unwrap_or_else(|e| e.into_inner()) = trimmed;
    }

    /// `GeomFill_EvolvedSection::GetInterval` (`cxx:233-241`): `TLaw->Bounds`.
    fn get_interval(&self) -> (f64, f64) {
        self.trimmed_law().bounds()
    }

    /// `GeomFill_EvolvedSection::GetDomain` (`cxx:243-249`): `First`, `Last`.
    fn get_domain(&self) -> (f64, f64) {
        (self.first, self.last)
    }

    /// `GeomFill_EvolvedSection::GetTolerance` (`cxx:251-267`): `SurfTol`
    /// everywhere, `BoundTol` at the two ends when it is the smaller one.
    fn get_tolerance(&self, bound_tol: f64, surf_tol: f64, _angle_tol: f64) -> Vec<f64> {
        let n = self.curve.poles.len();
        let mut tol3d = vec![surf_tol; n];
        if bound_tol < surf_tol && n > 0 {
            tol3d[0] = bound_tol;
            tol3d[n - 1] = bound_tol;
        }
        tol3d
    }

    /// `GeomFill_EvolvedSection::BarycentreOfSurf` (`cxx:269-292`): the mean of
    /// 21 section points times the mean of 21 law values, over 21 * 21.
    fn barycentre_of_surf(&self) -> Result<GpPnt, String> {
        let mut u = self.section.first_parameter();
        let delta = (self.curve.last_parameter() - u) / 20.0;
        let mut bary = occt_core::gp::GpXyz::zero();
        for _ in 0..=20 {
            let p = self.curve.d0(u);
            bary = bary.added(p.xyz());
            u += delta;
        }

        u = self.first;
        let delta = (self.last - self.first) / 20.0;
        let mut b = 0.0;
        for _ in 0..=20 {
            b += self.my_law.value(u);
            u += delta;
        }
        let c = bary.multiplied(b / (21.0 * 21.0));
        Ok(GpPnt::new(c.x(), c.y(), c.z()))
    }

    /// `GeomFill_EvolvedSection::MaximalSection` (`cxx:294-312`): the section
    /// length times the largest law value on 21 samples.
    fn maximal_section(&self) -> f64 {
        let l = curve_length(self.section.as_ref());

        let delta = (self.last - self.first) / 20.0;
        let mut max = 0.0;
        let mut u = self.first;
        for _ in 0..=20 {
            let val = self.my_law.value(u);
            if val > max {
                max = val;
            }
            u += delta;
        }
        l * max
    }

    /// `GeomFill_EvolvedSection::GetMinimalWeight` (`cxx:314-318`).
    fn get_minimal_weight(&self) -> Result<Vec<f64>, String> {
        Ok(self.weights_array())
    }

    /// `GeomFill_EvolvedSection::IsConstant` (`cxx:320-326`): false, `Error = 0`.
    fn is_constant(&self, error: &mut f64) -> bool {
        *error = 0.0;
        false
    }

    /// `GeomFill_EvolvedSection::ConstantSection` (`cxx:328-347`): only valid
    /// for a constant law, which `IsConstant` never reports.
    fn constant_section(&self) -> Result<Arc<dyn Curve>, String> {
        Err("StdFail_NotDone: The Law is not Constant!".into())
    }
}

/// The `GeomAbs_Shape` of the `u8` code the trait passes around.
fn shape_of(s: u8) -> Shape {
    match s {
        0 => Shape::C0,
        1 => Shape::G1,
        2 => Shape::C1,
        3 => Shape::G2,
        4 => Shape::C2,
        5 => Shape::C3,
        _ => Shape::CN,
    }
}
