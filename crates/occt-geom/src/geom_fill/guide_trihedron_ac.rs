//! `GeomFill_GuideTrihedronAC` — the trihedron driven by a guide curve
//! through the arc-length ratio between the path and the guide.
//!
//! Source: `TKGeomAlgo/GeomFill/GeomFill_GuideTrihedronAC.hxx` and
//! `GeomFill_GuideTrihedronAC.cxx:36-395` (OCCT 8.0.0): constructor
//! `cxx:36-45`, `D0` `cxx:58-83`, `D1` `cxx:90-151`, `D2` `cxx:158-262`,
//! `Copy` `cxx:266-272`, `SetCurve` `cxx:276-288`, `NbIntervals`
//! `cxx:292-307`, `Intervals` `cxx:310-329`, `SetInterval` `cxx:334-348`,
//! `GetAverageLaw` `cxx:352-372`, `IsConstant` `cxx:376-379`,
//! `IsOnlyBy3dCurve` `cxx:384-387`, `Origine` `cxx:392-395`.
//!
//! Translation notes:
//! - `myGuideAC` / `myCurveAC` are `Approx_CurvlinFunc` case 1
//!   ([`CurvlinFunc::from_curve`]); `GetSParameter` / `GetUParameter` are
//!   [`CurvlinFunc::s_parameter`] / [`CurvlinFunc::u_parameter`].
//! - `myGuide` / `myTrimG` / `myCurPointOnGuide` come from
//!   [`TrihedronWithGuide`]; `myTrimG = myGuide->Trim(...)` is
//!   [`GeomTrimmedCurveBasis`], the parameter-preserving trim.
//! - `GeomLib::FuseIntervals` defaults are `Epspar = 1.0e-9`,
//!   `IsAdjustToFirstInterval = false` (`GeomLib.hxx:206-210`), i.e.
//!   [`crate::geom_lib::fuse_intervals`] with those arguments.
//! - The two `TODO: finding #8` comments (`cxx:72-75`, `cxx:123-125`) are kept
//!   verbatim: OCCT normalizes without a zero guard, and the port must not
//!   invent one.

use std::sync::{Arc, Mutex};

use occt_core::gp::{GpPnt, GpVec};
use occt_core::kernel::geomabs::Shape;
use occt_core::precision::PCONFUSION;

use crate::approx_curvlin_func::CurvlinFunc;
use crate::curve::Curve;
use crate::geom_fill::trihedron_law::{TrihedronD1, TrihedronD2, TrihedronLaw, TrihedronLawCurve};
use crate::geom_fill::trihedron_with_guide::{TrihedronWithGuide, TrihedronWithGuideState};
use crate::geom_lib::fuse_intervals;
use crate::trimmed::GeomTrimmedCurveBasis;

/// `Lguide` / `L` tolerance and the `1.e-7` of the two `Approx_CurvlinFunc`
/// constructors (`cxx:41`, `cxx:283`).
const AC_TOL: f64 = 1.0e-7;
/// Default `Epspar` of `GeomLib::FuseIntervals` (`GeomLib.hxx:209`).
const FUSE_EPS: f64 = 1.0e-9;

/// `GeomFill_GuideTrihedronAC` (`GeomFill_GuideTrihedronAC.hxx`).
pub struct GuideTrihedronAC {
    /// `myCurve` / `myTrimmed` (`GeomFill_TrihedronLaw.hxx:111-112`).
    data: TrihedronLawCurve,
    /// `myGuide` (`hxx`): non-null, the constructor argument.
    guide: Arc<dyn Curve>,
    /// `myTrimG` / `myCurPointOnGuide` of the base class.
    gstate: Mutex<TrihedronWithGuideState>,
    /// `myCurveAC` (`hxx`).
    curve_ac: Mutex<Option<CurvlinFunc>>,
    /// `myGuideAC` (`hxx`).
    guide_ac: Mutex<Option<CurvlinFunc>>,
    /// `Lguide` (`hxx`): the guide length, set by the constructor.
    lguide: f64,
    /// `L` (`hxx`): the path length, set by `SetCurve` (`cxx:285`).
    l: Mutex<f64>,
    /// `UTol` (`hxx`): the `myCurve->Trim` tolerance of `SetInterval`
    /// (`cxx:337`); the port's `TrihedronLawCurve::set_interval` stores the
    /// parameter range only, so the value is kept for fidelity.
    #[allow(dead_code)]
    utol: f64,
    /// `STol` (`hxx`): set by the constructor, never read again.
    #[allow(dead_code)]
    stol: f64,
    /// `Orig1` / `Orig2` (`hxx`): the multi-edge path origins (`cxx:43-44`).
    orig: Mutex<(f64, f64)>,
}

impl GuideTrihedronAC {
    /// `GeomFill_GuideTrihedronAC(guide)` (`cxx:36-45`): `myCurve` null,
    /// `myGuide = myTrimG = guide`, `myGuideAC = new Approx_CurvlinFunc(guide,
    /// 1.e-7)`, `Lguide = GetLength()`, `UTol = STol = PConfusion`,
    /// `Orig1 = 0`, `Orig2 = 1`.
    pub fn new(guide: Arc<dyn Curve>) -> Self {
        let guide_ac = CurvlinFunc::from_curve(guide.clone(), AC_TOL);
        let lguide = guide_ac.length();
        Self {
            data: TrihedronLawCurve::new(),
            guide: guide.clone(),
            gstate: Mutex::new(TrihedronWithGuideState {
                guide: Some(guide.clone()),
                trim_g: Some(guide),
                cur_point_on_guide: GpPnt::new(0.0, 0.0, 0.0),
            }),
            curve_ac: Mutex::new(None),
            guide_ac: Mutex::new(Some(guide_ac)),
            lguide,
            l: Mutex::new(0.0),
            utol: PCONFUSION,
            stol: PCONFUSION,
            orig: Mutex::new((0.0, 1.0)),
        }
    }

    /// `Lguide` (`hxx`).
    pub fn guide_length(&self) -> f64 {
        self.lguide
    }

    fn origin(&self) -> (f64, f64) {
        *self.orig.lock().unwrap_or_else(|e| e.into_inner())
    }

    fn length(&self) -> f64 {
        *self.l.lock().unwrap_or_else(|e| e.into_inner())
    }

    /// `GeomFill_GuideTrihedronAC::Copy` body (`cxx:268-272`): a new law on the
    /// same guide, the same path and the same `Orig1` / `Orig2`.
    fn duplicate(&self) -> GuideTrihedronAC {
        let opy = GuideTrihedronAC::new(self.guide.clone());
        if let Some(curve) = self.data.curve() {
            TrihedronLaw::set_curve(&opy, curve);
        }
        let (o1, o2) = self.origin();
        opy.origine(o1, o2);
        opy
    }

    /// `myCurveAC->GetSParameter(Param)` (`cxx:61`).
    fn s_parameter(&self, param: f64) -> Result<f64, String> {
        let mut guard = self.curve_ac.lock().unwrap_or_else(|e| e.into_inner());
        match guard.as_mut() {
            Some(ac) => Ok(ac.s_parameter(param)),
            None => Err("GeomFill_GuideTrihedronAC: myCurveAC is null".into()),
        }
    }

    /// `myGuideAC->GetUParameter(*myGuide, S, 1)` (`cxx:63`).
    fn guide_parameter(&self, s: f64) -> Result<f64, String> {
        let mut guard = self.guide_ac.lock().unwrap_or_else(|e| e.into_inner());
        match guard.as_mut() {
            Some(ac) => Ok(ac.u_parameter(s)),
            None => Err("GeomFill_GuideTrihedronAC: myGuideAC is null".into()),
        }
    }

    /// `OrigG = Orig1 + s * (Orig2 - Orig1)` (`cxx:62`).
    fn orig_g(&self, s: f64) -> f64 {
        let (o1, o2) = self.origin();
        o1 + s * (o2 - o1)
    }

    fn set_cur_point(&self, pg: GpPnt) {
        self.gstate
            .lock()
            .unwrap_or_else(|e| e.into_inner())
            .cur_point_on_guide = pg;
    }

    fn trimmed_guide(&self) -> Result<Arc<dyn Curve>, String> {
        self.trim_g()
            .ok_or_else(|| "GeomFill_GuideTrihedronAC: myTrimG is null".to_string())
    }
}

impl TrihedronWithGuide for GuideTrihedronAC {
    fn guide_state(&self) -> &Mutex<TrihedronWithGuideState> {
        &self.gstate
    }

    /// `GeomFill_GuideTrihedronAC::Origine` (`cxx:392-395`).
    fn origine(&self, param1: f64, param2: f64) {
        *self.orig.lock().unwrap_or_else(|e| e.into_inner()) = (param1, param2);
    }

    /// `GeomFill_GuideTrihedronAC::Copy` (`cxx:266-272`), typed as a guide law
    /// for `GeomFill_LocationGuide::Copy` (`cxx:519-521`).
    fn copy_with_guide(&self) -> Arc<dyn TrihedronWithGuide> {
        Arc::new(self.duplicate())
    }
}

impl TrihedronLaw for GuideTrihedronAC {
    /// `GeomFill_GuideTrihedronAC::SetCurve` (`cxx:276-288`): the path and its
    /// `Approx_CurvlinFunc`, hence `L`.
    fn set_curve(&self, curve: Arc<dyn Curve>) -> bool {
        self.data.set_curve(curve.clone());
        let ac = CurvlinFunc::from_curve(curve, AC_TOL);
        *self.l.lock().unwrap_or_else(|e| e.into_inner()) = ac.length();
        *self.curve_ac.lock().unwrap_or_else(|e| e.into_inner()) = Some(ac);
        true
    }

    /// `GeomFill_GuideTrihedronAC::Copy` (`cxx:266-272`).
    fn copy(&self) -> Arc<dyn TrihedronLaw> {
        Arc::new(self.duplicate())
    }

    /// `GeomFill_GuideTrihedronAC::D0` (`cxx:58-83`).
    fn d0(&self, param: f64) -> Result<(GpVec, GpVec, GpVec), String> {
        let s = self.s_parameter(param)?;
        let orig_g = self.orig_g(s);
        let t_g = self.guide_parameter(orig_g)?;

        let curve = self.data.curve_or_err("GeomFill_GuideTrihedronAC")?;
        let (_p, to) = curve.d1(param);
        let pg = self.trimmed_guide()?.d0(t_g);
        self.set_cur_point(pg);

        let n = GpVec::from_pnts(&_p, &pg);
        // TODO: finding #8 - no zero-magnitude guard before
        // Normalized()/Magnitude() (`cxx:72-75`).
        let normal = n.normalized();
        let b = to.crossed(&normal);
        let bi_normal = b.divided(b.magnitude());
        let tangent = normal.crossed(&bi_normal).normalized();
        Ok((tangent, normal, bi_normal))
    }

    /// `GeomFill_GuideTrihedronAC::D1` (`cxx:90-151`).
    fn d1(&self, param: f64) -> Result<TrihedronD1, String> {
        let s = self.s_parameter(param)?;
        let orig_g = self.orig_g(s);
        let t_g = self.guide_parameter(orig_g)?;

        let curve = self.data.curve_or_err("GeomFill_GuideTrihedronAC")?;
        let (p, to, d_to) = curve.d2(param);
        let (pg, tg) = {
            let guide = self.trimmed_guide()?;
            let (pg, tg, _) = guide.d2(t_g);
            (pg, tg)
        };
        self.set_cur_point(pg);

        let mut n = GpVec::from_pnts(&p, &pg);
        let mut norm = n.magnitude();
        if norm < 1.0e-12 {
            norm = 1.0;
        }
        n = n.divided(norm);

        // TODO: finding #8 - no zero-magnitude guard before TG.Magnitude()/L
        // division (`cxx:123-125`).
        let (o1, o2) = self.origin();
        let l = self.length();
        let dtg = (o2 - o1) * (to.magnitude() / tg.magnitude()) * (self.lguide / l);
        let dn = tg
            .multiplied_scalar(dtg)
            .subtracted(&to)
            .divided(norm);

        let normal = n;
        let b = to.crossed(&normal);
        let norm_b = b.magnitude();
        let b = b.divided(norm_b);
        let bi_normal = b;
        let tangent = normal.crossed(&bi_normal).normalized();

        let d_normal = n
            .multiplied_scalar(-n.dot(&dn))
            .added(&dn);
        let b_prim = d_to
            .crossed(&normal)
            .added(&to.crossed(&d_normal));
        let d_bi_normal = b
            .multiplied_scalar(-b.dot(&b_prim))
            .added(&b_prim)
            .divided(norm_b);
        let d_tangent = normal
            .crossed(&d_bi_normal)
            .added(&d_normal.crossed(&bi_normal));

        Ok(TrihedronD1 {
            tangent,
            dtangent: d_tangent,
            normal,
            dnormal: d_normal,
            binormal: bi_normal,
            dbinormal: d_bi_normal,
        })
    }

    /// `GeomFill_GuideTrihedronAC::D2` (`cxx:158-262`).
    fn d2(&self, param: f64) -> Result<TrihedronD2, String> {
        let s = self.s_parameter(param)?;
        let orig_g = self.orig_g(s);
        let t_g = self.guide_parameter(orig_g)?;

        let curve = self.data.curve_or_err("GeomFill_GuideTrihedronAC")?;
        let (p, to, d_to, d2_to) = curve.d3(param);
        let (pg, tg, dtg, _d2tg) = self.trimmed_guide()?.d3(t_g);
        self.set_cur_point(pg);

        let n_to = to.magnitude();
        let n2_to = to.square_magnitude();
        let n_tg = tg.magnitude();
        let n2_tp = tg.square_magnitude();
        let (o1, o2) = self.origin();
        let l = self.length();
        let ratio = self.lguide / l;
        let dtg_dt = (o2 - o1) * (n_to / n_tg) * ratio;

        let n = GpVec::from_pnts(&p, &pg);
        let norm = n.magnitude();
        let dn = tg.multiplied_scalar(dtg_dt).subtracted(&to);

        // `cxx:195-200`.
        let d2tp_dt2 = (o2 - o1)
            * ratio
            * (d_to.dot(&to) / (n_to * n_tg) - n2_to * tg.dot(&dtg) * ratio / (n2_tp * n2_tp));
        let d2n = dtg
            .multiplied_scalar(dtg_dt * dtg_dt)
            .added(&tg.multiplied_scalar(d2tp_dt2))
            .subtracted(&d_to);

        let (mut n, mut dn, mut d2n) = (n, dn, d2n);
        if norm > 1.0e-9 {
            n = n.divided(norm);
            dn = dn.divided(norm);
            d2n = d2n.divided(norm);
        }
        let normal = n;

        let mut tn = to.crossed(&normal);
        let norma = tn.magnitude();
        if norma > 1.0e-9 {
            tn = tn.divided(norma);
        }
        let bi_normal = tn;
        let tangent = normal.crossed(&bi_normal);

        let ndn = n.dot(&dn);
        let d_normal = n.multiplied_scalar(-ndn).added(&dn);
        let mut d_tn = d_to
            .crossed(&normal)
            .added(&to.crossed(&d_normal))
            .divided(norma);
        let tndtn = tn.dot(&d_tn);
        let d_bi_normal = tn.multiplied_scalar(-tndtn).added(&d_tn);
        let d_tangent = normal
            .crossed(&d_bi_normal)
            .added(&d_normal.crossed(&bi_normal));

        let tn2 = tn.square_magnitude();
        let d2_normal = dn
            .multiplied_scalar(-2.0 * ndn)
            .added(&n.multiplied_scalar(3.0 * ndn * ndn - (dn.square_magnitude() + n.dot(&d2n))))
            .added(&d2n);

        d_tn = d2_to
            .crossed(&normal)
            .added(&d_to.crossed(&d_normal).multiplied_scalar(2.0))
            .added(&to.crossed(&d2_normal))
            .divided(norma);
        let d2_bi_normal = d_tn
            .multiplied_scalar(-2.0 * tndtn)
            .added(&tn.multiplied_scalar(3.0 * tndtn * tndtn - (tn2 + tn.dot(&d_tn))))
            .added(&d_tn);

        let d2_tangent = d2_normal
            .crossed(&bi_normal)
            .added(&d_normal.crossed(&d_bi_normal).multiplied_scalar(2.0))
            .added(&normal.crossed(&d2_bi_normal));

        Ok(TrihedronD2 {
            tangent,
            dtangent: d_tangent,
            d2tangent: d2_tangent,
            normal,
            dnormal: d_normal,
            d2normal: d2_normal,
            binormal: bi_normal,
            dbinormal: d_bi_normal,
            d2binormal: d2_bi_normal,
        })
    }

    /// `GeomFill_GuideTrihedronAC::NbIntervals` (`cxx:292-307`).
    fn nb_intervals(&self, s: u8) -> i32 {
        let shape = s;
        let disc_c = match self.curve_ac.lock().unwrap_or_else(|e| e.into_inner()).as_ref() {
            Some(ac) => ac.intervals(shape_of(shape)),
            None => vec![0.0, 1.0],
        };
        let disc_g = match self.guide_ac.lock().unwrap_or_else(|e| e.into_inner()).as_ref() {
            Some(ac) => ac.intervals(shape_of(shape)),
            None => vec![0.0, 1.0],
        };
        let seq = fuse_intervals(&disc_c, &disc_g, FUSE_EPS, false);
        seq.len().saturating_sub(1) as i32
    }

    /// `GeomFill_GuideTrihedronAC::Intervals` (`cxx:310-329`): the fused knots
    /// remapped through `myCurveAC->GetUParameter(*myCurve, knot, 1)`.
    fn intervals(&self, s: u8) -> Vec<f64> {
        let shape = s;
        let disc_c = match self.curve_ac.lock().unwrap_or_else(|e| e.into_inner()).as_ref() {
            Some(ac) => ac.intervals(shape_of(shape)),
            None => vec![0.0, 1.0],
        };
        let disc_g = match self.guide_ac.lock().unwrap_or_else(|e| e.into_inner()).as_ref() {
            Some(ac) => ac.intervals(shape_of(shape)),
            None => vec![0.0, 1.0],
        };
        let seq = fuse_intervals(&disc_c, &disc_g, FUSE_EPS, false);
        let mut guard = self.curve_ac.lock().unwrap_or_else(|e| e.into_inner());
        let Some(ac) = guard.as_mut() else {
            return seq;
        };
        seq.into_iter().map(|knot| ac.u_parameter(knot)).collect()
    }

    /// `GeomFill_GuideTrihedronAC::SetInterval` (`cxx:334-348`).
    fn set_interval(&self, first: f64, last: f64) -> Result<(), String> {
        self.data.set_interval(first, last)?;
        let s_f = self.s_parameter(first)?;
        let s_l = self.s_parameter(last)?;
        let u = self.orig_g(s_f);
        let s_f = self.guide_parameter(u)?;
        let u = self.orig_g(s_l);
        let s_l = self.guide_parameter(u)?;
        let guide = self.guide.clone();
        let trim_g: Arc<dyn Curve> = Arc::new(GeomTrimmedCurveBasis::new(guide, s_f, s_l));
        self.gstate
            .lock()
            .unwrap_or_else(|e| e.into_inner())
            .trim_g = Some(trim_g);
        Ok(())
    }

    /// `GeomFill_TrihedronLaw::GetInterval`
    /// (`GeomFill_TrihedronLaw.cxx:62-67`) through the `myTrimmed` trim.
    fn get_interval(&self) -> (f64, f64) {
        self.data.interval()
    }

    /// `GeomFill_GuideTrihedronAC::GetAverageLaw` (`cxx:352-372`): 20 samples
    /// starting at `myCurve->FirstParameter()`, step `Delta = (Last -
    /// First) / 20.001`.
    fn get_average_law(&self) -> (GpVec, GpVec, GpVec) {
        let mut a_tangent = GpVec::zero();
        let mut a_normal = GpVec::zero();
        let mut a_bi_normal = GpVec::zero();
        let Some(curve) = self.data.curve() else {
            return (a_tangent, a_normal, a_bi_normal);
        };
        let delta = (curve.last_parameter() - curve.first_parameter()) / 20.001;
        for ii in 1..=20 {
            let t = curve.first_parameter() + (ii - 1) as f64 * delta;
            if let Ok((tv, n, b)) = self.d0(t) {
                a_tangent = a_tangent.added(&tv);
                a_normal = a_normal.added(&n);
                a_bi_normal = a_bi_normal.added(&b);
            }
        }
        a_tangent = a_tangent.divided(20.0);
        a_normal = a_normal.divided(20.0);
        a_bi_normal = a_bi_normal.divided(20.0);
        (a_tangent, a_normal, a_bi_normal)
    }

    /// `GeomFill_GuideTrihedronAC::IsConstant` (`cxx:376-379`): false.
    fn is_constant(&self) -> bool {
        false
    }

    /// `GeomFill_GuideTrihedronAC::IsOnlyBy3dCurve` (`cxx:384-387`): false.
    fn is_only_by_3d_curve(&self) -> bool {
        false
    }
}

/// The `GeomAbs_Shape` of the `u8` shape code the trait passes around.
fn shape_of(s: u8) -> Shape {
    match s {
        0 => Shape::C0,
        1 => Shape::G1,
        2 => Shape::C1,
        3 => Shape::G2,
        4 => Shape::C2,
        5 => Shape::C3,
        6 => Shape::CN,
        _ => Shape::CN,
    }
}
