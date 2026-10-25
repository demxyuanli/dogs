//! `GeomFill_LocationGuide` — the location law that keeps the section in
//! contact with a guide curve through the revolution surface of the section.
//!
//! Source: `TKGeomAlgo/GeomFill/GeomFill_LocationGuide.cxx:1-1500` and
//! `GeomFill_LocationGuide.hxx:39-219` (OCCT 8.0.0).
//!
//! Layout:
//! - `mod.rs` (this file): the state, the constructor, `Set` / `SetCurve` /
//!   `Copy` / `SetTrsf`, the interval and norm queries, and
//!   `ComputeAutomaticLaw`,
//! - `rotation.rs`: `SetRotation` (`cxx:184-467`), `InitX` (`cxx:1352-1425`)
//!   and the `math_FunctionSetRoot` solve shared with `D0`,
//! - `eval.rs`: `D0` / `D1` / `D2` (`cxx:584-1144`).
//!
//! Translation notes:
//! - `occ::handle<Adaptor3d_Curve>` is `Arc<dyn Curve>`; `myCurve` and
//!   `myTrimmed` follow the [`TrihedronLawCurve`] convention (the trim keeps
//!   the parameters of `myCurve`).
//! - `myGuide` is `myLaw->Guide()` (`cxx:162`); when the guide is not periodic
//!   the constructor extends it by one percent on each side (`cxx:164-173`,
//!   `Trim(f, l, delta * 1.e-7)`) with [`GeomTrimmedCurveBasis`], the
//!   parameter-preserving `Geom_TrimmedCurve` view.
//! - `Extrema_ExtCS` (`cxx:324`) is `curve_surface_extrema_all` over the guide
//!   and the revolution surface; the minimum square distance scan of
//!   `cxx:375-389` is reproduced (first minimum wins, as `<` does in OCCT).
//! - `GeomFill_FunctionGuide` + `math_FunctionSetRoot` (`cxx:339-360`,
//!   `cxx:617-645`) are [`FunctionGuide`] and [`MathFunctionSetRoot`].
//! - `D2`'s `if (WithTrans)` block (`cxx:925-930`) multiplies `M`, `DM` and
//!   `D2M` *before* they are assigned anywhere in that function, so it has no
//!   observable effect and the non-rotation frame it returns is untransformed.
//!   That OCCT inconsistency is kept: see [`crate::geom_fill::location_guide::eval`].
//! - `D1` / `D2` return `false` as soon as `rotation` is set (`cxx:768-770`,
//!   `cxx:933-935`); the C1 / C2 approximation of a guided sweep therefore
//!   fails exactly as it does in OCCT.

pub mod eval;
pub mod rotation;

use std::sync::{Arc, Mutex, MutexGuard};

use occt_core::gp::{GpMat, GpPnt2d, GpVec};
use occt_core::precision::{CONFUSION, PCONFUSION};

use crate::curve::Curve;
use crate::geom_fill::location_law::LocationLaw;
use crate::geom_fill::pipe_error::GeomFillPipeError;
use crate::geom_fill::section_law::SectionLaw;
use crate::geom_fill::trihedron_with_guide::TrihedronWithGuide;
use crate::geom_lib::fuse_intervals;
use crate::trimmed::GeomTrimmedCurveBasis;

/// `myNbPts` (`cxx:162`): the number of guide samples of `SetRotation`.
pub const NB_PTS: usize = 21;

/// The mutable part of `GeomFill_LocationGuide` (`hxx:195-218`). OCCT hands one
/// law object to several holders and mutates it through any of them, so the
/// port keeps the whole state behind a `Mutex` and takes `&self`.
pub struct LgState {
    /// `myLaw` (`hxx:196`).
    pub law: Arc<dyn TrihedronWithGuide>,
    /// `mySec` (`hxx:197`).
    pub sec: Option<Arc<dyn SectionLaw>>,
    /// `myCurve` (`hxx:198`).
    pub curve: Option<Arc<dyn Curve>>,
    /// `myGuide` (`hxx:200`), possibly extended by the constructor.
    pub guide: Option<Arc<dyn Curve>>,
    /// `myTrimmed` (`hxx:201`).
    pub trimmed: Option<Arc<dyn Curve>>,
    /// `myNbPts` (`hxx:202`).
    pub nb_pts: usize,
    /// `rotation` (`hxx:203`): "contact or not".
    pub rotation: bool,
    /// `OrigParam1` / `OrigParam2` (`hxx:204-205`).
    pub orig_param1: f64,
    pub orig_param2: f64,
    /// `Uf` / `Ul` (`hxx:206-207`): the section curve parameter range.
    pub uf: f64,
    pub ul: f64,
    /// `myFirstS` / `myLastS` (`hxx:208-209`).
    pub first_s: f64,
    pub last_s: f64,
    /// `ratio` (`hxx:210`).
    pub ratio: f64,
    /// `WithTrans` / `Trans` (`hxx:211-212`).
    pub with_trans: bool,
    pub trans: GpMat,
    /// `TolRes` (`hxx:213`), 1-based in OCCT, stored 0-based here.
    pub tol_res: [f64; 3],
    /// `Inf` / `Sup` (`hxx:214-215`).
    pub inf: [f64; 3],
    pub sup: [f64; 3],
    /// `X` / `R` (`hxx:216-217`).
    pub x: [f64; 3],
    pub r: [f64; 3],
    /// `myPoles2d` (`hxx:190`): `(1, 2, 1, myNbPts)`; entry `k` is the OCCT
    /// column `k + 1`, `[0]` the `(t, w)` pair and `[1]` the `(angle, v)` pair.
    pub poles2d: Vec<[GpPnt2d; 2]>,
    /// `myStatus` (`hxx:218`).
    pub status: GeomFillPipeError,
}

/// `GeomFill_LocationGuide` (`GeomFill_LocationGuide.hxx:40-219`).
pub struct LocationGuide {
    st: Mutex<LgState>,
}

impl LocationGuide {
    /// `GeomFill_LocationGuide(Triedre)` (`cxx:144-181`).
    ///
    /// OCCT dereferences `myLaw->Guide()` immediately; a law without a guide
    /// has no faithful continuation, so this returns `Err` instead of
    /// dereferencing a null handle.
    pub fn new(law: Arc<dyn TrihedronWithGuide>) -> Result<Self, String> {
        let guide = law
            .guide()
            .ok_or_else(|| "Standard_NullObject: GeomFill_LocationGuide::GeomFill_LocationGuide".to_string())?;
        // `cxx:164-173`: extend a non-periodic guide by one percent per side.
        let guide = if guide.is_periodic() {
            guide
        } else {
            let f = guide.first_parameter();
            let l = guide.last_parameter();
            let delta = (l - f) / 100.0;
            Arc::new(GeomTrimmedCurveBasis::new(guide, f - delta, l + delta))
        };
        let mut tol_res = [0.0; 3];
        tol_res[0] = 1.0e-6;
        Ok(Self {
            st: Mutex::new(LgState {
                law,
                sec: None,
                curve: None,
                guide: Some(guide),
                trimmed: None,
                nb_pts: NB_PTS,
                rotation: false,
                orig_param1: 0.0,
                orig_param2: 1.0,
                uf: 0.0,
                ul: 0.0,
                first_s: -505.0e77,
                last_s: -505.0e77,
                ratio: 0.0,
                with_trans: false,
                trans: GpMat::identity(),
                tol_res,
                inf: [0.0; 3],
                sup: [0.0; 3],
                x: [0.0; 3],
                r: [0.0; 3],
                poles2d: vec![[GpPnt2d::zero(); 2]; NB_PTS],
                status: GeomFillPipeError::PipeOk,
            }),
        })
    }

    /// The shared state (`hxx:190-218`).
    pub(crate) fn lock(&self) -> MutexGuard<'_, LgState> {
        self.st.lock().unwrap_or_else(|e| e.into_inner())
    }

    /// `GeomFill_LocationGuide::Set` (`cxx:472-503`). Returns `LastAngle`;
    /// `myStatus` reports an impossible contact exactly as OCCT does, without
    /// turning it into an error.
    pub fn set(
        &self,
        section: Option<Arc<dyn SectionLaw>>,
        rotat: bool,
        s_first: f64,
        s_last: f64,
        prec_angle: f64,
    ) -> Result<f64, String> {
        {
            let mut st = self.lock();
            st.status = GeomFillPipeError::PipeOk;
            st.first_s = s_first;
            st.last_s = s_last;
            if let Some(curve) = &st.curve {
                st.ratio = (s_last - s_first) / (curve.last_parameter() - curve.first_parameter());
            } else {
                st.ratio = 0.0;
            }
            st.sec = section;
        }
        if rotat {
            self.set_rotation(prec_angle)
        } else {
            self.lock().rotation = false;
            Ok(prec_angle)
        }
    }

    /// `GeomFill_LocationGuide::EraseRotation` (`cxx:507-512`).
    pub fn erase_rotation(&self) {
        let mut st = self.lock();
        st.rotation = false;
        if st.status == GeomFillPipeError::ImpossibleContact {
            st.status = GeomFillPipeError::PipeOk;
        }
    }

    /// `GeomFill_LocationGuide::SetOrigine` (`cxx:1432-1436`).
    pub fn set_origine(&self, param1: f64, param2: f64) {
        let mut st = self.lock();
        st.orig_param1 = param1;
        st.orig_param2 = param2;
    }

    /// `GeomFill_LocationGuide::Section` (`cxx:1312-1315`).
    pub fn section(&self) -> Result<Arc<dyn Curve>, String> {
        let st = self.lock();
        let sec = st
            .sec
            .clone()
            .ok_or_else(|| "Standard_NullObject: GeomFill_LocationGuide::Section".to_string())?;
        drop(st);
        sec.constant_section()
    }

    /// `GeomFill_LocationGuide::Guide` (`cxx:1318-1321`).
    pub fn guide(&self) -> Option<Arc<dyn Curve>> {
        self.lock().guide.clone()
    }

    /// `GeomFill_LocationGuide::ComputeAutomaticLaw` (`cxx:1438-1472`): the
    /// `(parameter, width)` pairs for the automatic-radius law, plus the first
    /// law status that stopped the scan.
    pub fn compute_automatic_law(&self) -> (GeomFillPipeError, Vec<GpPnt2d>) {
        let st = self.lock();
        let nb_pts = st.nb_pts;
        let (Some(curve), Some(law)) = (st.curve.clone(), Some(st.law.clone())) else {
            return (st.status, Vec::new());
        };
        let f = curve.first_parameter();
        let l = curve.last_parameter();
        let mut par_and_rad = Vec::with_capacity(nb_pts);
        for ii in 1..=nb_pts {
            let t = ((nb_pts - ii) as f64 * f + (ii - 1) as f64 * l) / (nb_pts - 1) as f64;
            let p = curve.d0(t);
            if law.d0(t).is_err() {
                return (law.error_status(), par_and_rad);
            }
            let point_on_guide = law.current_point_on_guide();
            par_and_rad.push(GpPnt2d::new(t, p.distance(&point_on_guide)));
        }
        (st.status, par_and_rad)
    }

    /// `myStatus` (`hxx:218`).
    pub fn error_status(&self) -> GeomFillPipeError {
        self.lock().status
    }

    /// `rotation` (`hxx:203`).
    pub fn is_rotating(&self) -> bool {
        self.lock().rotation
    }

    /// `TolRes` / `Inf` / `Sup` / `X` / `R` (`hxx:213-217`), 1-based accessors.
    pub fn tol_res(&self) -> [f64; 3] {
        self.lock().tol_res
    }

    pub fn bounds(&self) -> ([f64; 3], [f64; 3]) {
        let st = self.lock();
        (st.inf, st.sup)
    }
}

impl LocationLaw for LocationGuide {
    /// `GeomFill_LocationGuide::SetCurve` (`cxx:534-552`).
    fn set_curve(&self, curve: Arc<dyn Curve>) -> bool {
        let prec_angle;
        {
            let mut st = self.lock();
            st.curve = Some(curve.clone());
            st.trimmed = Some(curve.clone());
            st.law.set_curve(curve);
            st.law.origine(st.orig_param1, st.orig_param2);
            st.status = st.law.error_status();
            if !st.rotation {
                return st.status == GeomFillPipeError::PipeOk;
            }
            prec_angle = st.poles2d[0][0].x();
        }
        match self.set_rotation(prec_angle) {
            Ok(_) => self.error_status() == GeomFillPipeError::PipeOk,
            Err(_) => false,
        }
    }

    /// `GeomFill_LocationGuide::GetCurve` (`cxx:556-559`).
    fn get_curve(&self) -> Option<Arc<dyn Curve>> {
        self.lock().curve.clone()
    }

    /// `GeomFill_LocationGuide::SetTrsf` (`cxx:563-580`).
    fn set_trsf(&self, transfo: &GpMat) {
        let mut st = self.lock();
        st.trans = *transfo;
        let aux = GpMat::identity().subtract(transfo);
        st.with_trans = false;
        for ii in 1..=3 {
            for jj in 1..=3 {
                if aux.value(ii, jj).abs() > 1.0e-14 {
                    st.with_trans = true;
                }
            }
        }
    }

    /// `GeomFill_LocationGuide::Copy` (`cxx:516-527`).
    fn copy(&self) -> Arc<dyn LocationLaw> {
        let (law, sec, rotation, first_s, last_s, prec_angle, trans) = {
            let st = self.lock();
            (
                st.law.copy_with_guide(),
                st.sec.clone(),
                st.rotation,
                st.first_s,
                st.last_s,
                st.poles2d[0][0].x(),
                st.trans,
            )
        };
        let copy = LocationGuide::new(law).expect("GeomFill_LocationGuide::Copy: law without guide");
        {
            let st = self.lock();
            copy.set_origine(st.orig_param1, st.orig_param2);
        }
        let _ = copy.set(sec, rotation, first_s, last_s, prec_angle);
        copy.set_trsf(&trans);
        Arc::new(copy)
    }

    /// `GeomFill_LocationGuide::D0` — see `eval.rs`.
    fn d0(&self, param: f64) -> Result<(GpMat, GpVec), String> {
        self.d0_frame(param)
    }

    /// `GeomFill_LocationGuide::D0(Param, M, V, Poles2d)` (`cxx:657-727`): the
    /// `Poles2d` array is ignored by OCCT, so it comes back empty.
    fn d0_with_poles(&self, param: f64) -> Result<(GpMat, GpVec, Vec<GpPnt2d>), String> {
        let (m, v) = self.d0_frame(param)?;
        Ok((m, v, Vec::new()))
    }

    /// `GeomFill_LocationGuide::D1` (`cxx:733-891`).
    fn d1(&self, param: f64) -> Result<crate::geom_fill::location_law::LocationD1, String> {
        self.d1_frame(param)
    }

    /// `GeomFill_LocationGuide::D2` (`cxx:896-1145`).
    fn d2(&self, param: f64) -> Result<crate::geom_fill::location_law::LocationD2, String> {
        self.d2_frame(param)
    }

    /// `GeomFill_LocationGuide::HasFirstRestriction` (`cxx:1149-1152`).
    fn has_first_restriction(&self) -> bool {
        false
    }

    /// `GeomFill_LocationGuide::HasLastRestriction` (`cxx:1155-1158`).
    fn has_last_restriction(&self) -> bool {
        false
    }

    /// `GeomFill_LocationGuide::TraceNumber` (`cxx:1162-1165`).
    fn trace_number(&self) -> i32 {
        0
    }

    /// `GeomFill_LocationGuide::ErrorStatus` (`cxx:1170-1173`).
    fn error_status(&self) -> GeomFillPipeError {
        LocationGuide::error_status(self)
    }

    /// `GeomFill_LocationGuide::NbIntervals` (`cxx:1177-1200`).
    fn nb_intervals(&self, s: u8) -> i32 {
        let st = self.lock();
        let nb_sec = st.trimmed.as_ref().map_or(1, |c| c.nb_intervals(s));
        let nb_law = st.law.nb_intervals(s);
        if nb_sec == 1 {
            return nb_law;
        }
        if nb_law == 1 {
            return nb_sec;
        }
        let int_c = st
            .trimmed
            .as_ref()
            .map(|c| c.parameter_intervals(s))
            .unwrap_or_default();
        let int_l = st.law.intervals(s);
        let inter = fuse_intervals(&int_c, &int_l, PCONFUSION * 0.99, false);
        inter.len() as i32 - 1
    }

    /// `GeomFill_LocationGuide::Intervals` (`cxx:1204-1232`).
    fn intervals(&self, s: u8) -> Vec<f64> {
        let st = self.lock();
        let nb_sec = st.trimmed.as_ref().map_or(1, |c| c.nb_intervals(s));
        let nb_law = st.law.nb_intervals(s);
        if nb_sec == 1 {
            return st.law.intervals(s);
        }
        if nb_law == 1 {
            return st
                .trimmed
                .as_ref()
                .map(|c| c.parameter_intervals(s))
                .unwrap_or_default();
        }
        let int_c = st
            .trimmed
            .as_ref()
            .map(|c| c.parameter_intervals(s))
            .unwrap_or_default();
        let int_l = st.law.intervals(s);
        fuse_intervals(&int_c, &int_l, PCONFUSION * 0.99, false)
    }

    /// `GeomFill_LocationGuide::SetInterval` (`cxx:1237-1240`):
    /// `myTrimmed = myCurve->Trim(First, Last, 0)`.
    fn set_interval(&self, first: f64, last: f64) {
        let mut st = self.lock();
        st.law.set_interval(first, last).ok();
        if let Some(curve) = st.curve.clone() {
            st.trimmed = Some(Arc::new(GeomTrimmedCurveBasis::new(curve, first, last)));
        }
    }

    /// `GeomFill_LocationGuide::GetInterval` (`cxx:1244-1248`).
    fn get_interval(&self) -> (f64, f64) {
        let st = self.lock();
        match &st.trimmed {
            Some(c) => (c.first_parameter(), c.last_parameter()),
            None => (0.0, 0.0),
        }
    }

    /// `GeomFill_LocationGuide::GetDomain` (`cxx:1252-1257`).
    fn get_domain(&self) -> (f64, f64) {
        let st = self.lock();
        match &st.curve {
            Some(c) => (c.first_parameter(), c.last_parameter()),
            None => (0.0, 0.0),
        }
    }

    /// `GeomFill_LocationGuide::SetTolerance` (`cxx:1260-1264`):
    /// `TolRes(1) = myGuide->Resolution(Tol3d)` and
    /// `Resolution(1, Tol3d, TolRes(2), TolRes(3))`.
    fn set_tolerance(&self, tol3d: f64, _tol2d: f64) {
        let mut st = self.lock();
        if let Some(guide) = &st.guide {
            st.tol_res[0] = guide.resolution(tol3d);
        }
        st.tol_res[1] = tol3d / 100.0;
        st.tol_res[2] = tol3d / 100.0;
    }

    /// `GeomFill_LocationGuide::Resolution` (`cxx:1269-1275`): `Tol / 100`
    /// on both 2d axes; `Index` is ignored.
    fn resolution(&self, _index: i32, tol: f64) -> Result<(f64, f64), String> {
        Ok((tol / 100.0, tol / 100.0))
    }

    /// `GeomFill_LocationGuide::GetMaximalNorm` (`cxx:1282-1285`).
    fn get_maximal_norm(&self) -> f64 {
        1.0
    }

    /// `GeomFill_LocationGuide::GetAverageLaw` (`cxx:1289-1307`).
    fn get_average_law(&self) -> (GpMat, GpVec) {
        let st = self.lock();
        let (v1, v2, v3) = st.law.get_average_law();
        let am = GpMat::from_cols(v1.xyz(), v2.xyz(), v3.xyz());
        let Some(trimmed) = &st.trimmed else {
            return (am, GpVec::zero());
        };
        let mut av = GpVec::zero();
        let delta = (trimmed.last_parameter() - trimmed.first_parameter()) / 10.0;
        let mut u = trimmed.first_parameter();
        for _ in 0..=st.nb_pts {
            let p = trimmed.d0(u);
            av = av.add(&GpVec::new(p.x(), p.y(), p.z()));
            u += delta;
        }
        (am, av.multiply_scalar(1.0 / (st.nb_pts + 1) as f64))
    }

    /// `GeomFill_LocationGuide::IsTranslation` (`cxx:1343-1346`): false.
    fn is_translation(&self, _error: &mut f64) -> bool {
        false
    }

    /// `GeomFill_LocationGuide::IsRotation` (`cxx:1326-1329`): false.
    fn is_rotation(&self, _error: &mut f64) -> bool {
        false
    }

    /// `GeomFill_LocationGuide::Rotation` (`cxx:1334-1337`):
    /// `Standard_NotImplemented`.
    fn rotation(&self) -> Result<occt_core::gp::GpPnt, String> {
        Err("Standard_NotImplemented: GeomFill_LocationGuide::Rotation".into())
    }
}

/// `Precision::Confusion()` — the tolerance `SetRotation` passes to
/// `Extrema_ExtCS` (`cxx:324`).
pub const EXTREMA_TOL: f64 = CONFUSION;
