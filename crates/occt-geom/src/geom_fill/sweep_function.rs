//! `GeomFill_SweepFunction` — the `Approx_SweepFunction` a general sweep is
//! approximated from.
//!
//! Source: `TKGeomAlgo/GeomFill/GeomFill_SweepFunction.hxx:39-171` and
//! `GeomFill_SweepFunction.cxx:36-372` (OCCT 8.0.0).
//!
//! Translation notes:
//! - The `NCollection_Array1` out-parameters become the returned vectors of
//!   [`ApproxSweepFunction`]; the loop bounds `L = Poles.Length()` follow from
//!   `SectionShape`.
//! - `myf` / `myfOnS` / `myRatio` are the OCCT members (`hxx:162-164`); every
//!   section parameter is `T = myfOnS + (Param - myf) * myRatio`
//!   (`cxx:60`, `cxx:98`, `cxx:147`).
//! - `GeomLib::FuseIntervals` (`cxx:252`, `cxx:297`) is
//!   [`crate::geom_lib::fuse_intervals`] with `Precision::PConfusion() * 0.99`.
//! - `M`, `V`, `DM`, `DV`, `D2M`, `D2V` are OCCT members (`hxx:165-169`) but
//!   are pure scratch values recomputed in every call; the port keeps them as
//!   locals.

use std::sync::Arc;

use occt_core::approx::{ApproxSweepFunction, SweepD1, SweepD2};
use occt_core::gp::{GpPnt, GpPnt2d, GpVec};
use occt_core::precision::PCONFUSION;

use crate::geom_fill::location_law::LocationLaw;
use crate::geom_fill::section_law::SectionLaw;
use crate::geom_lib::fuse_intervals;

/// `GeomFill_SweepFunction` (`GeomFill_SweepFunction.hxx:39-171`).
pub struct SweepFunction {
    /// `myLoc` (`hxx:159`).
    loc: Arc<dyn LocationLaw>,
    /// `mySec` (`hxx:160`).
    sec: Arc<dyn SectionLaw>,
    /// `myf` (`hxx:161`).
    f: f64,
    /// `myfOnS` (`hxx:162`).
    f_on_s: f64,
    /// `myRatio` (`hxx:163`).
    ratio: f64,
}

impl SweepFunction {
    /// `GeomFill_SweepFunction::GeomFill_SweepFunction` (`cxx:36-47`).
    pub fn new(
        section: Arc<dyn SectionLaw>,
        location: Arc<dyn LocationLaw>,
        first_parameter: f64,
        first_parameter_on_s: f64,
        ratio_parameter_on_s: f64,
    ) -> Self {
        Self {
            loc: location,
            sec: section,
            f: first_parameter,
            f_on_s: first_parameter_on_s,
            ratio: ratio_parameter_on_s,
        }
    }

    /// `mySec` (`hxx:160`).
    pub fn section(&self) -> &Arc<dyn SectionLaw> {
        &self.sec
    }

    /// `myLoc` (`hxx:159`).
    pub fn location(&self) -> &Arc<dyn LocationLaw> {
        &self.loc
    }
}

impl ApproxSweepFunction for SweepFunction {
    /// `GeomFill_SweepFunction::D0` (`cxx:52-81`).
    fn d0(
        &self,
        param: f64,
        _first: f64,
        _last: f64,
    ) -> Result<(Vec<GpPnt>, Vec<GpPnt2d>, Vec<f64>), String> {
        let t = self.f_on_s + (param - self.f) * self.ratio;
        let (m, v, poles2d) = self.loc.d0_with_poles(param)?;
        let (mut poles, weights) = self.sec.d0(t)?;
        for pole in poles.iter_mut() {
            let p = occt_core::gp::GpXyz::new(pole.x(), pole.y(), pole.z()).multiplied_mat(&m);
            *pole = GpPnt::new(p.x() + v.x(), p.y() + v.y(), p.z() + v.z());
        }
        Ok((poles, poles2d, weights))
    }

    /// `GeomFill_SweepFunction::D1` (`cxx:85-126`).
    fn d1(&self, param: f64, _first: f64, _last: f64) -> Result<SweepD1, String> {
        let t = self.f_on_s + (param - self.f) * self.ratio;
        let loc = self.loc.d1(param)?;
        let sec = self.sec.d1(t)?;
        let (mut poles, mut dpoles, mut weights, mut dweights) = sec;
        let l = poles.len().min(dpoles.len());
        for ii in 0..l {
            let p = occt_core::gp::GpXyz::new(poles[ii].x(), poles[ii].y(), poles[ii].z());
            let mut pprim = occt_core::gp::GpXyz::new(
                dpoles[ii].x(),
                dpoles[ii].y(),
                dpoles[ii].z(),
            );
            pprim = pprim.multiply_scalar(self.ratio);
            dweights[ii] *= self.ratio;
            pprim = pprim.multiplied_mat(&loc.m);
            pprim = pprim.added(&loc.dm.multiplied(&p));
            pprim = pprim.added(&occt_core::gp::GpXyz::new(loc.dv.x(), loc.dv.y(), loc.dv.z()));
            dpoles[ii] = GpVec::from_xyz(&pprim);
            let moved = p.multiplied_mat(&loc.m);
            poles[ii] = GpPnt::new(
                moved.x() + loc.v.x(),
                moved.y() + loc.v.y(),
                moved.z() + loc.v.z(),
            );
        }
        let _ = &mut weights;
        Ok(SweepD1 {
            poles,
            dpoles,
            poles2d: loc.poles2d,
            dpoles2d: loc.dpoles2d,
            weights,
            dweights,
        })
    }

    /// `GeomFill_SweepFunction::D2` (`cxx:131-186`).
    fn d2(&self, param: f64, _first: f64, _last: f64) -> Result<SweepD2, String> {
        let t = self.f_on_s + (param - self.f) * self.ratio;
        let squareratio = self.ratio * self.ratio;
        let loc = self.loc.d2(param)?;
        let sec = self.sec.d2(t)?;
        let (mut poles, mut dpoles, mut d2poles, mut weights, mut dweights, mut d2weights) = sec;
        let l = poles.len().min(dpoles.len());
        for ii in 0..l {
            let p = occt_core::gp::GpXyz::new(poles[ii].x(), poles[ii].y(), poles[ii].z());
            let mut pprim =
                occt_core::gp::GpXyz::new(dpoles[ii].x(), dpoles[ii].y(), dpoles[ii].z());
            pprim = pprim.multiply_scalar(self.ratio);
            dweights[ii] *= self.ratio;
            let mut psecn =
                occt_core::gp::GpXyz::new(d2poles[ii].x(), d2poles[ii].y(), d2poles[ii].z());
            psecn = psecn.multiply_scalar(squareratio);
            d2weights[ii] *= squareratio;

            psecn = psecn.multiplied_mat(&loc.m);
            psecn = psecn.added(&loc.dm.multiplied(&pprim).multiply_scalar(2.0));
            psecn = psecn.added(&loc.d2m.multiplied(&p));
            psecn = psecn.added(&occt_core::gp::GpXyz::new(
                loc.d2v.x(),
                loc.d2v.y(),
                loc.d2v.z(),
            ));
            d2poles[ii] = GpVec::from_xyz(&psecn);

            pprim = pprim.multiplied_mat(&loc.m);
            pprim = pprim.added(&loc.dm.multiplied(&p));
            pprim = pprim.added(&occt_core::gp::GpXyz::new(loc.dv.x(), loc.dv.y(), loc.dv.z()));
            dpoles[ii] = GpVec::from_xyz(&pprim);

            let moved = p.multiplied_mat(&loc.m);
            poles[ii] = GpPnt::new(
                moved.x() + loc.v.x(),
                moved.y() + loc.v.y(),
                moved.z() + loc.v.z(),
            );
        }
        let _ = &mut weights;
        Ok(SweepD2 {
            poles,
            dpoles,
            d2poles,
            poles2d: loc.poles2d,
            dpoles2d: loc.dpoles2d,
            d2poles2d: loc.d2poles2d,
            weights,
            dweights,
            d2weights,
        })
    }

    /// `GeomFill_SweepFunction::Nb2dCurves` (`cxx:190-193`).
    fn nb_2d_curves(&self) -> i32 {
        self.loc.nb_2d_curves()
    }

    /// `GeomFill_SweepFunction::SectionShape` (`cxx:197-200`).
    fn section_shape(&self) -> (i32, i32, i32) {
        self.sec.section_shape()
    }

    /// `GeomFill_SweepFunction::Knots` (`cxx:204-207`).
    fn knots(&self) -> Vec<f64> {
        self.sec.knots()
    }

    /// `GeomFill_SweepFunction::Mults` (`cxx:211-214`).
    fn mults(&self) -> Vec<i32> {
        self.sec.mults()
    }

    /// `GeomFill_SweepFunction::IsRational` (`cxx:218-221`).
    fn is_rational(&self) -> bool {
        self.sec.is_rational()
    }

    /// `GeomFill_SweepFunction::NbIntervals` (`cxx:225-254`).
    fn nb_intervals(&self, s: u8) -> i32 {
        let nb_sec = self.sec.nb_intervals(s);
        let nb_loc = self.loc.nb_intervals(s);
        if nb_sec == 1 {
            return nb_loc;
        }
        if nb_loc == 1 {
            return nb_sec;
        }
        let int_s: Vec<f64> = self
            .sec
            .intervals(s)
            .into_iter()
            .map(|t| (t - self.f_on_s) / self.ratio + self.f)
            .collect();
        let int_l = self.loc.intervals(s);
        let inter = fuse_intervals(&int_s, &int_l, PCONFUSION * 0.99, false);
        inter.len() as i32 - 1
    }

    /// `GeomFill_SweepFunction::Intervals` (`cxx:258-300`).
    fn intervals(&self, s: u8) -> Vec<f64> {
        let nb_sec = self.sec.nb_intervals(s);
        let nb_loc = self.loc.nb_intervals(s);
        if nb_sec == 1 {
            return self.loc.intervals(s);
        }
        if nb_loc == 1 {
            return self
                .sec
                .intervals(s)
                .into_iter()
                .map(|t| (t - self.f_on_s) / self.ratio + self.f)
                .collect();
        }
        let int_s: Vec<f64> = self
            .sec
            .intervals(s)
            .into_iter()
            .map(|t| (t - self.f_on_s) / self.ratio + self.f)
            .collect();
        let int_l = self.loc.intervals(s);
        fuse_intervals(&int_s, &int_l, PCONFUSION * 0.99, false)
    }

    /// `GeomFill_SweepFunction::SetInterval` (`cxx:304-311`).
    fn set_interval(&self, first: f64, last: f64) {
        self.loc.set_interval(first, last);
        let uf = self.f + (first - self.f) * self.ratio;
        let ul = self.f + (last - self.f) * self.ratio;
        self.sec.set_interval(uf, ul);
    }

    /// `GeomFill_SweepFunction::GetTolerance` (`cxx:315-321`).
    fn get_tolerance(&self, bound_tol: f64, surf_tol: f64, angle_tol: f64) -> Vec<f64> {
        self.sec.get_tolerance(bound_tol, surf_tol, angle_tol)
    }

    /// `GeomFill_SweepFunction::Resolution` (`cxx:325-332`).
    fn resolution(&self, index: i32, tol: f64) -> Result<(f64, f64), String> {
        self.loc.resolution(index, tol)
    }

    /// `GeomFill_SweepFunction::SetTolerance` (`cxx:336-340`).
    fn set_tolerance(&self, tol3d: f64, tol2d: f64) {
        self.sec.set_tolerance(tol3d, tol2d);
        self.loc.set_tolerance(tol3d, tol2d);
    }

    /// `GeomFill_SweepFunction::BarycentreOfSurf` (`cxx:344-355`).
    fn barycentre_of_surf(&self) -> Result<GpPnt, String> {
        let bary = self.sec.barycentre_of_surf()?;
        let (am, translate) = self.loc.get_average_law();
        let p = occt_core::gp::GpXyz::new(bary.x(), bary.y(), bary.z()).multiplied_mat(&am);
        Ok(GpPnt::new(
            p.x() + translate.x(),
            p.y() + translate.y(),
            p.z() + translate.z(),
        ))
    }

    /// `GeomFill_SweepFunction::MaximalSection` (`cxx:359-364`).
    fn maximal_section(&self) -> Result<f64, String> {
        Ok(self.sec.maximal_section() * self.loc.get_maximal_norm())
    }

    /// `GeomFill_SweepFunction::GetMinimalWeight` (`cxx:368-372`).
    fn get_minimal_weight(&self) -> Result<Vec<f64>, String> {
        self.sec.get_minimal_weight()
    }
}
