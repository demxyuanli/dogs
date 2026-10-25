//! `GeomFill_UniformSection` — a section law that repeats one curve.
//!
//! Source: `TKGeomAlgo/GeomFill/GeomFill_UniformSection.hxx` and
//! `GeomFill_UniformSection.cxx:34-295` (OCCT 8.0.0).
//!
//! Translation notes:
//! - `mySection = C->Copy()` and `myCurve = down_cast<Geom_BSplineCurve>(C)`
//!   (`cxx:39-53`). When `C` is not a B-spline, OCCT converts it with
//!   `GeomConvert::CurveToBSplineCurve(C, Convert_QuasiAngular)` and removes the
//!   wrap-around knot of a periodic result (`cxx:44-52`); the port calls
//!   `convert_bspl::curve_to_bspline_curve` and `GeomBSplineCurve::remove_knot`.
//! - `D0`/`D1`/`D2` fill the caller's arrays from `myCurve->Poles()` and
//!   `myCurve->WeightsArray()`; the port returns them (see
//!   [`SectionLaw`]).
//! - `BSplineSurface` (`cxx:110-140`) builds the degree-1-in-V surface through
//!   `GeomBSplineSurface` directly; the port returns the same poles and knot
//!   arrays as a [`GeomBSplineSurface`].

use std::sync::Arc;

use occt_core::convert::ParameterisationType;
use occt_core::gp::{GpPnt, GpVec};
use occt_core::precision::CONFUSION;

use crate::bspline_curve::GeomBSplineCurve;
use crate::bspline_surface::GeomBSplineSurface;
use crate::convert_bspl::curve_to_bspline_curve;
use crate::curve::Curve;
use crate::gcpnts::curve_length;
use crate::geom_fill::section_law::SectionLaw;

/// `GeomFill_UniformSection` (`GeomFill_UniformSection.hxx:35-100`).
pub struct UniformSection {
    /// `mySection` (`hxx:80`): the copied input curve.
    section: Arc<dyn Curve>,
    /// `myCurve` (`hxx:81`): the B-spline representation.
    curve: GeomBSplineCurve,
    /// `First` (`hxx:82`).
    first: f64,
    /// `Last` (`hxx:83`).
    last: f64,
}

impl UniformSection {
    /// `GeomFill_UniformSection::GeomFill_UniformSection` (`cxx:34-54`).
    pub fn new(
        c: Arc<dyn Curve>,
        first_parameter: f64,
        last_parameter: f64,
    ) -> Result<Self, String> {
        let mut curve = match curve_to_bspline_curve(c.as_ref(), ParameterisationType::QuasiAngular)
        {
            Ok(b) => b,
            Err(_) => return Err("Standard_DomainError: GeomFill_UniformSection".into()),
        };
        if curve.is_periodic() {
            let m = curve.degree() as i32 / 2 + 1;
            curve
                .remove_knot(1, m, CONFUSION)
                .map_err(|_| "Standard_ConstructionError: GeomFill_UniformSection".to_string())?;
        }
        Ok(Self {
            section: c,
            curve,
            first: first_parameter,
            last: last_parameter,
        })
    }

    /// `myCurve` (`hxx:81`).
    pub fn curve(&self) -> &GeomBSplineCurve {
        &self.curve
    }

    /// The stored weights array (`Geom_BSplineCurve.hxx:805`).
    fn weights_array(&self) -> Vec<f64> {
        match &self.curve.weights {
            Some(w) => w.clone(),
            None => vec![1.0; self.curve.poles.len()],
        }
    }
}

impl SectionLaw for UniformSection {
    /// `GeomFill_UniformSection::D0` (`cxx:57-65`).
    fn d0(&self, _param: f64) -> Result<(Vec<GpPnt>, Vec<f64>), String> {
        Ok((self.curve.poles.clone(), self.weights_array()))
    }

    /// `GeomFill_UniformSection::D1` (`cxx:70-84`): the derivatives are zero.
    fn d1(&self, _param: f64) -> Result<(Vec<GpPnt>, Vec<GpVec>, Vec<f64>, Vec<f64>), String> {
        let zeros = vec![GpVec::zero(); self.curve.poles.len()];
        let dweights = vec![0.0; self.curve.poles.len()];
        Ok((
            self.curve.poles.clone(),
            zeros,
            self.weights_array(),
            dweights,
        ))
    }

    /// `GeomFill_UniformSection::D2` (`cxx:86-105`): both derivatives zero.
    fn d2(
        &self,
        _param: f64,
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
        let n = self.curve.poles.len();
        Ok((
            self.curve.poles.clone(),
            vec![GpVec::zero(); n],
            vec![GpVec::zero(); n],
            self.weights_array(),
            vec![0.0; n],
            vec![0.0; n],
        ))
    }

    /// `GeomFill_UniformSection::BSplineSurface` (`cxx:110-140`).
    fn bspline_surface(&self) -> Option<Arc<GeomBSplineSurface>> {
        let (u_knots, u_mults) = self.curve.distinct_knots_and_mults();
        let poles = self
            .curve
            .poles
            .iter()
            .map(|p| vec![*p, *p])
            .collect::<Vec<_>>();
        let v_knots = vec![self.first, self.last];
        let v_mults = vec![2, 2];
        let mut surface = GeomBSplineSurface::from_poles_knots_mults(
            poles,
            u_knots,
            v_knots,
            u_mults,
            v_mults,
            self.curve.degree(),
            1,
        )
        .ok()?;
        // `cxx:136-138`: the constructor receives `myCurve->IsPeriodic()` as
        // `UPeriodic`; the port sets the same flag after construction.
        surface.set_u_periodic_flag(self.curve.is_periodic());
        Some(Arc::new(surface))
    }

    /// `GeomFill_UniformSection::SectionShape` (`cxx:147-153`).
    fn section_shape(&self) -> (i32, i32, i32) {
        (
            self.curve.nb_poles() as i32,
            self.curve.distinct_knots_and_mults().0.len() as i32,
            self.curve.degree() as i32,
        )
    }

    /// `GeomFill_UniformSection::Knots` (`cxx:155-158`).
    fn knots(&self) -> Vec<f64> {
        self.curve.distinct_knots_and_mults().0
    }

    /// `GeomFill_UniformSection::Mults` (`cxx:160-164`).
    fn mults(&self) -> Vec<i32> {
        self.curve.distinct_knots_and_mults().1
    }

    /// `GeomFill_UniformSection::IsRational` (`cxx:166-170`).
    fn is_rational(&self) -> bool {
        self.curve.is_rational()
    }

    /// `GeomFill_UniformSection::IsUPeriodic` (`cxx:172-176`).
    fn is_u_periodic(&self) -> bool {
        self.curve.is_periodic()
    }

    /// `GeomFill_UniformSection::IsVPeriodic` (`cxx:178-182`): true.
    fn is_v_periodic(&self) -> bool {
        true
    }

    /// `GeomFill_UniformSection::NbIntervals` (`cxx:192-196`): 1.
    fn nb_intervals(&self, _s: u8) -> i32 {
        1
    }

    /// `GeomFill_UniformSection::Intervals` (`cxx:201-209`).
    fn intervals(&self, _s: u8) -> Vec<f64> {
        vec![self.first, self.last]
    }

    /// `GeomFill_UniformSection::SetInterval` (`cxx:211-216`): does nothing.
    fn set_interval(&self, _first: f64, _last: f64) {}

    /// `GeomFill_UniformSection::GetInterval` (`cxx:218-224`).
    fn get_interval(&self) -> (f64, f64) {
        (self.first, self.last)
    }

    /// `GeomFill_UniformSection::GetDomain` (`cxx:226-231`).
    fn get_domain(&self) -> (f64, f64) {
        (self.first, self.last)
    }

    /// `GeomFill_UniformSection::GetTolerance` (`cxx:236-251`).
    fn get_tolerance(&self, bound_tol: f64, surf_tol: f64, _angle_tol: f64) -> Vec<f64> {
        let n = self.curve.poles.len();
        let mut tol3d = vec![surf_tol; n];
        if bound_tol < surf_tol && n > 0 {
            tol3d[0] = bound_tol;
            tol3d[n - 1] = bound_tol;
        }
        tol3d
    }

    /// `GeomFill_UniformSection::BarycentreOfSurf` (`cxx:255-270`).
    fn barycentre_of_surf(&self) -> Result<GpPnt, String> {
        let u0 = self.section.first_parameter();
        let delta = (self.curve.last_parameter() - u0) / 20.0;
        let mut sum = occt_core::gp::GpXyz::zero();
        let mut u = u0;
        for _ in 0..=20 {
            let p = self.curve.d0(u);
            sum = sum.added(&occt_core::gp::GpXyz::new(p.x(), p.y(), p.z()));
            u += delta;
        }
        let c = sum.divided(21.0);
        Ok(GpPnt::new(c.x(), c.y(), c.z()))
    }

    /// `GeomFill_UniformSection::MaximalSection` (`cxx:272-277`):
    /// `GCPnts_AbscissaPoint::Length(GeomAdaptor_Curve(mySection))`.
    fn maximal_section(&self) -> f64 {
        curve_length(self.section.as_ref())
    }

    /// `GeomFill_UniformSection::GetMinimalWeight` (`cxx:279-282`).
    fn get_minimal_weight(&self) -> Result<Vec<f64>, String> {
        Ok(self.weights_array())
    }

    /// `GeomFill_UniformSection::IsConstant` (`cxx:284-288`).
    fn is_constant(&self, error: &mut f64) -> bool {
        *error = 0.0;
        true
    }

    /// `GeomFill_UniformSection::ConstantSection` (`cxx:290-295`).
    fn constant_section(&self) -> Result<Arc<dyn Curve>, String> {
        Ok(Arc::clone(&self.section))
    }
}
