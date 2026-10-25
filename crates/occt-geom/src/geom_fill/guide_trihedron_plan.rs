//! `GeomFill_GuideTrihedronPlan` — the trihedron of a sweep along a guide
//! curve met by the planes normal to the spine.
//!
//! Source: `TKGeomAlgo/GeomFill/GeomFill_GuideTrihedronPlan.hxx:38-131` and
//! `GeomFill_GuideTrihedronPlan.cxx:51-582` (OCCT 8.0.0): `InGoodPeriod`
//! `cxx:51-65`, constructor `cxx:69-86`, `Init` `cxx:92-191`, `SetCurve`
//! `cxx:195-205`, `D0` `cxx:214-262`, `D1` `cxx:268-352`, `D2` `cxx:358-389`,
//! `Copy` `cxx:393-398`, `ErrorStatus` `cxx:402-405`, `NbIntervals`
//! `cxx:411-432`, `Intervals` `cxx:436-455`, `SetInterval` `cxx:459-462`,
//! `GetAverageLaw` `cxx:466-489`, `IsConstant` `cxx:493-506`,
//! `IsOnlyBy3dCurve` `cxx:510-513`, `Origine` `cxx:519`, `InitX` `cxx:525-582`.
//!
//! Translation notes:
//! - `myGuide` / `myTrimG` / `myCurPointOnGuide` belong to
//!   [`TrihedronWithGuide`]; `myTrimG = myGuide->Trim(...)` is
//!   [`GeomTrimmedCurveBasis`], the parameter-preserving trim.
//! - `X` / `Inf` / `Sup` / `Pole` / `myStatus` are written by `Init` and
//!   `InitX`, which OCCT calls through `D0` / `D1` on a shared handle: they
//!   live behind one `Mutex` as [`PlanState`].
//! - `frenet->D0` failing leaves `Tangent` untouched (`cxx:98` declares the
//!   three vectors *before* the loop, and `cxx:131` ignores the result), so
//!   `Init` keeps the previous tangent on failure, exactly like OCCT. In `D0`
//!   / `D1` the same failure has no out-parameters to preserve; the port
//!   reports it as the law's error.
//! - `D2` (`cxx:358-389`) evaluates `myCurve->D2` and `frenet->D2` into its
//!   out-parameters and then returns `false` unconditionally (`cxx:388`).
//!   Nothing of that frame reaches a caller, so the port reports the failure
//!   directly.
//! - `IntCurveSurface_HInter` (`cxx:135`) lives in `occt-topo`
//!   (`crates/occt-topo/src/intcurvesurface`) and cannot be reached from
//!   `occt-geom`. For the `GeomAbs_Plane` case that call resolves to
//!   `IntCurveSurface_InterImpl::InternalPerformCurveQuadric`
//!   (`IntCurveSurface_Inter.pxx:178`, `:423-439`), i.e.
//!   `IntCurveSurface_QuadricCurveExactInterUtils::PerformIntersection`
//!   (`IntCurveSurface_QuadricCurveExactInterUtils.pxx:50-136`): for every
//!   C1 interval of the guide a `math_FunctionSample` of
//!   `IntCurveSurface_TheHCurveTool::NbSamples` points feeds
//!   `math_FunctionAllRoots` on
//!   `IntCurveSurface_TheQuadCurvFuncOfTheQuadCurvExactHInter`, whose value is
//!   `IntSurf_Quadric::Distance` of the curve point (`IntSurf_Quadric.cxx:173-179`
//!   gives the plane form `n.P + d`). [`plane_intersections`] reproduces that
//!   branch with the same function, sample count and tolerances. The conic
//!   short-cuts of `IntCurveSurface_InterImpl::PerformBounds`
//!   (`IntCurveSurface_Inter.pxx:119-137`, `PerformConicSurf*`) solve the same
//!   signed-distance equation analytically and are not ported separately.

use std::sync::{Arc, Mutex};

use occt_core::elib::clib::in_period;
use occt_core::gp::{GpPnt, GpPnt2d, GpVec};
use occt_core::kernel::geomabs::Shape;
use occt_math::{FunctionWithDerivative, MathFunctionRoot};

use crate::curve::Curve;
use crate::geom_fill::frenet::Frenet;
use crate::geom_fill::pipe_error::GeomFillPipeError;
use crate::geom_fill::plane_curve_inter::plane_intersections;
use crate::geom_fill::plan_func::PlanFunc;
use crate::geom_fill::trihedron_law::{TrihedronD1, TrihedronD2, TrihedronLaw, TrihedronLawCurve};
use crate::geom_fill::trihedron_with_guide::{TrihedronWithGuide, TrihedronWithGuideState};
use crate::trimmed::GeomTrimmedCurveBasis;

/// `myNbPts` (`cxx:80`): the number of poles computed by `Init`.
const NB_PTS: usize = 20;
/// `Iter` of the two `math_FunctionRoot` calls (`cxx:229`, `cxx:285`).
const ROOT_ITERATIONS: i32 = 50;

/// `X(1)` / `Inf(1)` / `Sup(1)` / `Pole` / `myStatus` (`hxx:121-129`).
struct PlanState {
    /// `Pole` (`hxx:121`): `Pole(1, ii)` stores `(path param, guide param)`.
    /// Stored 0-based, index `ii - 1` is OCCT's `Pole(1, ii)`.
    pole: Vec<GpPnt2d>,
    /// `X(1)` (`hxx:122`): the initial guide parameter of `InitX`.
    x: f64,
    /// `Inf(1)` (`hxx:124`).
    inf: f64,
    /// `Sup(1)` (`hxx:125`).
    sup: f64,
    /// `myStatus` (`hxx:129`).
    status: GeomFillPipeError,
}

/// `GeomFill_GuideTrihedronPlan` (`GeomFill_GuideTrihedronPlan.hxx:38-131`).
pub struct GuideTrihedronPlan {
    /// `myCurve` / `myTrimmed` (`GeomFill_TrihedronLaw.hxx:111-112`).
    data: TrihedronLawCurve,
    /// `myGuide` of the base class: non-null, the constructor argument.
    guide: Arc<dyn Curve>,
    /// `myTrimG` / `myCurPointOnGuide` of the base class.
    gstate: Mutex<TrihedronWithGuideState>,
    /// `X` / `Inf` / `Sup` / `Pole` / `myStatus`.
    state: Mutex<PlanState>,
    /// `XTol(1)` (`hxx:123`): `myGuide->Resolution(1.e-6)` (`cxx:85`).
    xtol: f64,
    /// `frenet` (`hxx:127`): the Frenet frame of the path.
    frenet: Arc<Frenet>,
}

impl GuideTrihedronPlan {
    /// `GeomFill_GuideTrihedronPlan(guide)` (`cxx:69-86`): `myCurve` null,
    /// `myGuide = myTrimG = guide`, `myNbPts = 20`, `XTol.Init(1.e-6)` then
    /// `XTol(1) = myGuide->Resolution(1.e-6)`, `myStatus = PipeOk`.
    pub fn new(guide: Arc<dyn Curve>) -> Self {
        let xtol = guide.resolution(1.0e-6);
        Self {
            data: TrihedronLawCurve::new(),
            gstate: Mutex::new(TrihedronWithGuideState {
                guide: Some(guide.clone()),
                trim_g: Some(guide.clone()),
                cur_point_on_guide: GpPnt::new(0.0, 0.0, 0.0),
            }),
            guide,
            state: Mutex::new(PlanState {
                pole: vec![GpPnt2d::zero(); NB_PTS],
                x: 0.0,
                inf: 0.0,
                sup: 0.0,
                status: GeomFillPipeError::PipeOk,
            }),
            xtol,
            frenet: Arc::new(Frenet::new()),
        }
    }

    fn lock_state(&self) -> std::sync::MutexGuard<'_, PlanState> {
        self.state.lock().unwrap_or_else(|e| e.into_inner())
    }

    fn set_status(&self, status: GeomFillPipeError) {
        self.lock_state().status = status;
    }

    fn root_range(&self) -> (f64, f64) {
        let state = self.lock_state();
        (state.inf, state.sup)
    }

    fn guess(&self) -> f64 {
        self.lock_state().x
    }

    fn pole_point(&self, index: usize) -> GpPnt2d {
        self.lock_state().pole[index]
    }

    fn set_pole_point(&self, index: usize, value: GpPnt2d) {
        self.lock_state().pole[index] = value;
    }

    fn set_range(&self, inf: f64, sup: f64) {
        let mut state = self.lock_state();
        state.inf = inf;
        state.sup = sup;
    }

    fn trimmed_guide(&self) -> Result<Arc<dyn Curve>, String> {
        self.trim_g()
            .ok_or_else(|| "GeomFill_GuideTrihedronPlan: myTrimG is null".to_string())
    }

    /// `GeomFill_GuideTrihedronPlan::Copy` body (`cxx:393-398`): a new law on
    /// the same guide carrying the same path.
    fn duplicate(&self) -> GuideTrihedronPlan {
        let copy = GuideTrihedronPlan::new(self.guide.clone());
        if let Some(curve) = self.data.curve() {
            TrihedronLaw::set_curve(&copy, curve);
        }
        copy
    }

    /// `GeomFill_GuideTrihedronPlan::Init` (`cxx:92-191`): the `myNbPts`
    /// guide parameters of the intersection of the normal planes with the
    /// guide.
    fn init(&self) {
        self.set_status(GeomFillPipeError::PipeOk);
        let Some(curve) = self.data.curve() else {
            return;
        };
        let guide = self.guide.clone();
        let f = curve.first_parameter();
        let l = curve.last_parameter();

        self.frenet.set_curve(curve.clone());
        let mut delta_g = (guide.last_parameter() - guide.first_parameter()) / 2.0;
        let g_first = guide.first_parameter();
        let g_last = guide.last_parameter();

        self.set_range(g_first - delta_g, g_last + delta_g);

        let trim_g: Arc<dyn Curve> = if !guide.is_periodic() {
            // `cxx:116-118`: the third argument (DeltaG * 1.e-7) is the trim
            // tolerance of the adaptor trim, which does not move the bounds.
            Arc::new(GeomTrimmedCurveBasis::new(
                guide.clone(),
                g_first - delta_g / 100.0,
                g_last + delta_g / 100.0,
            ))
        } else {
            guide.clone()
        };
        self.gstate
            .lock()
            .unwrap_or_else(|e| e.into_inner())
            .trim_g = Some(trim_g.clone());

        delta_g /= 3.0;
        let mut w = 0.0;
        let mut tangent = GpVec::zero();
        for ii in 1..=NB_PTS {
            let t = ((NB_PTS - ii) as f64 * f + (ii - 1) as f64 * l) / (NB_PTS - 1) as f64;
            let p = curve.d0(t);
            // `cxx:131`: the result is ignored; a failing `D0` leaves the
            // vectors at their previous value (`cxx:98`).
            if let Ok((tng, _n, _b)) = self.frenet.d0(t) {
                tangent = tng;
            }

            let points = plane_intersections(trim_g.as_ref(), &p, &tangent);
            if points.is_empty() {
                w = if (g_last - w).abs() > (g_first - w).abs() {
                    g_first
                } else {
                    g_last
                };
                self.set_status(GeomFillPipeError::PlaneNotIntersectGuide);
            } else {
                let mut best = points[0];
                let mut dmin = p.distance(&best.1);
                for &(pw, pp) in points.iter().skip(1) {
                    let d = p.distance(&pp);
                    if d < dmin {
                        best = (pw, pp);
                        dmin = d;
                    }
                }
                w = best.0;
            }

            if ii > 1 {
                let prev = self.pole_point(ii - 2).y();
                if (w - prev).abs() > delta_g && guide.is_periodic() {
                    in_good_period(prev, guide.period(), &mut w);
                }
            }

            self.set_pole_point(ii - 1, GpPnt2d::new(t, w));
        }
    }

    /// `GeomFill_GuideTrihedronPlan::InitX` (`cxx:525-582`): the initial guide
    /// parameter, interpolated between the two neighbouring poles.
    fn init_x(&self, param: f64) {
        let guide = &self.guide;
        let mut state = self.lock_state();
        let pole = &state.pole;

        let mut ideb = 0usize;
        let mut ifin = NB_PTS - 1;

        let valeur = pole[ideb].x();
        if param == valeur {
            ifin = ideb + 1;
        }

        let valeur = pole[ifin].x();
        if param == valeur {
            ideb = ifin - 1;
        }

        while ideb + 1 != ifin {
            let idemi = (ideb + ifin) / 2;
            let valeur = pole[idemi].x();
            if valeur < param {
                ideb = idemi;
            } else if valeur > param {
                ifin = idemi;
            } else {
                ideb = idemi;
                ifin = ideb + 1;
            }
        }

        let t1 = pole[ideb].x();
        let t2 = pole[ifin].x();
        let diff = t2 - t1;
        let mut x = if diff > 1.0e-7 {
            let b = (param - t1) / diff;
            let a = (t2 - param) / diff;
            pole[ideb].y() * a + pole[ifin].y() * b
        } else {
            (pole[ideb].y() + pole[ifin].y()) / 2.0
        };
        if guide.is_periodic() {
            x = in_period(x, guide.first_parameter(), guide.last_parameter());
        }
        state.x = x;
    }
}

impl TrihedronWithGuide for GuideTrihedronPlan {
    fn guide_state(&self) -> &Mutex<TrihedronWithGuideState> {
        &self.gstate
    }

    /// `GeomFill_GuideTrihedronPlan::Origine` (`cxx:519`): nothing.
    fn origine(&self, _param1: f64, _param2: f64) {}

    /// `GeomFill_GuideTrihedronPlan::Copy` (`cxx:393-398`), typed as a guide
    /// law for `GeomFill_LocationGuide::Copy` (`cxx:519-521`).
    fn copy_with_guide(&self) -> Arc<dyn TrihedronWithGuide> {
        Arc::new(self.duplicate())
    }
}

impl TrihedronLaw for GuideTrihedronPlan {
    /// `GeomFill_GuideTrihedronPlan::SetCurve` (`cxx:195-205`): stores the
    /// path and runs `Init` when it is not null.
    fn set_curve(&self, curve: Arc<dyn Curve>) -> bool {
        self.data.set_curve(curve);
        self.init();
        true
    }

    /// `GeomFill_GuideTrihedronPlan::Copy` (`cxx:393-398`).
    fn copy(&self) -> Arc<dyn TrihedronLaw> {
        Arc::new(self.duplicate())
    }

    /// `GeomFill_GuideTrihedronPlan::ErrorStatus` (`cxx:402-405`).
    fn error_status(&self) -> GeomFillPipeError {
        self.lock_state().status
    }

    /// `GeomFill_GuideTrihedronPlan::D0` (`cxx:214-262`).
    fn d0(&self, param: f64) -> Result<(GpVec, GpVec, GpVec), String> {
        let curve = self.data.curve_or_err("GeomFill_GuideTrihedronPlan")?;
        let p = curve.d0(param);
        let (tangent, _normal, _binormal) = self.frenet.d0(param)?;

        self.init_x(param);
        let mut e = PlanFunc::new(p, tangent, self.guide.clone());
        let (inf, sup) = self.root_range();
        let result =
            MathFunctionRoot::new_with_bounds(&mut e, self.guess(), self.xtol, inf, sup, ROOT_ITERATIONS);

        if result.is_done() {
            let res = result.root();
            let pprime = self.trimmed_guide()?.d0(res);
            let n = GpVec::from_pnts(&p, &pprime);
            let normal = n.normalized();
            let bi_normal = tangent.crossed(&normal).normalized();
            Ok((tangent, normal, bi_normal))
        } else {
            self.set_status(GeomFillPipeError::PlaneNotIntersectGuide);
            Err("GeomFill_GuideTrihedronPlan::D0: plane does not intersect the guide".into())
        }
    }

    /// `GeomFill_GuideTrihedronPlan::D1` (`cxx:268-352`).
    fn d1(&self, param: f64) -> Result<TrihedronD1, String> {
        let curve = self.data.curve_or_err("GeomFill_GuideTrihedronPlan")?;
        let (p, to) = curve.d1(param);
        let frame = self.frenet.d1(param)?;

        self.init_x(param);
        let mut e = PlanFunc::new(p, frame.tangent, self.guide.clone());
        let (inf, sup) = self.root_range();
        let result =
            MathFunctionRoot::new_with_bounds(&mut e, self.guess(), self.xtol, inf, sup, ROOT_ITERATIONS);

        if !result.is_done() {
            self.set_status(GeomFillPipeError::PlaneNotIntersectGuide);
            return Err("GeomFill_GuideTrihedronPlan::D1: plane does not intersect the guide".into());
        }

        let res = result.root();
        let (pg, tg) = self.trimmed_guide()?.d1(res);
        let mut n = GpVec::from_pnts(&p, &pg);
        let mut norm = n.magnitude();
        if norm < 1.0e-12 {
            norm = 1.0;
        }
        n = n.divided(norm);

        let normal = n;
        let bi_normal = frame.tangent.crossed(&normal);

        let mut dedx = 0.0;
        FunctionWithDerivative::derivative(&mut e, res, &mut dedx);
        let mut dedt = 0.0;
        e.dedt(res, &to, &frame.dtangent, &mut dedt);
        let dtg_dt = -dedt / dedx;

        let dn = tg.multiplied_scalar(dtg_dt).subtracted(&to);
        let d_normal = dn
            .subtracted(&n.multiplied_scalar(n.dot(&dn)))
            .divided(norm);
        let d_bi_normal = frame
            .tangent
            .crossed(&d_normal)
            .added(&frame.dtangent.crossed(&normal));

        Ok(TrihedronD1 {
            tangent: frame.tangent,
            dtangent: frame.dtangent,
            normal,
            dnormal: d_normal,
            binormal: bi_normal,
            dbinormal: d_bi_normal,
        })
    }

    /// `GeomFill_GuideTrihedronPlan::D2` (`cxx:358-389`): OCCT fills the frame
    /// from `myCurve->D2` / `frenet->D2` and returns `false` regardless.
    fn d2(&self, _param: f64) -> Result<TrihedronD2, String> {
        Err("GeomFill_GuideTrihedronPlan::D2: returned false by OCCT (cxx:388)".into())
    }

    /// `GeomFill_GuideTrihedronPlan::NbIntervals` (`cxx:411-432`): the shape is
    /// promoted by one order before asking `myCurve`.
    fn nb_intervals(&self, s: u8) -> i32 {
        match self.data.curve() {
            Some(curve) => curve.nb_intervals(promoted_shape(s) as u8),
            None => 1,
        }
    }

    /// `GeomFill_GuideTrihedronPlan::Intervals` (`cxx:436-455`).
    fn intervals(&self, s: u8) -> Vec<f64> {
        match self.data.curve() {
            Some(curve) => curve.parameter_intervals(promoted_shape(s) as u8),
            None => vec![0.0, 0.0],
        }
    }

    /// `GeomFill_GuideTrihedronPlan::SetInterval` (`cxx:459-462`):
    /// `myTrimmed = myCurve->Trim(First, Last, Precision::Confusion())`.
    fn set_interval(&self, first: f64, last: f64) -> Result<(), String> {
        self.data.set_interval(first, last)
    }

    /// `GeomFill_TrihedronLaw::GetInterval`
    /// (`GeomFill_TrihedronLaw.cxx:62-67`).
    fn get_interval(&self) -> (f64, f64) {
        self.data.interval()
    }

    /// `GeomFill_GuideTrihedronPlan::GetAverageLaw` (`cxx:466-489`).
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

    /// `GeomFill_GuideTrihedronPlan::IsConstant` (`cxx:493-506`).
    fn is_constant(&self) -> bool {
        let Some(curve) = self.data.curve() else {
            return false;
        };
        if !curve.is_line() || !self.guide.is_line() {
            return false;
        }
        let (Some(curve_line), Some(guide_line)) = (curve.gp_line(), self.guide.gp_line()) else {
            return false;
        };
        let angle = curve_line.angle(&guide_line);
        angle < 1.0e-12 || (2.0 * std::f64::consts::PI - angle) < 1.0e-12
    }

    /// `GeomFill_GuideTrihedronPlan::IsOnlyBy3dCurve` (`cxx:510-513`): false.
    fn is_only_by_3d_curve(&self) -> bool {
        false
    }
}

/// `InGoodPeriod(Prec, Period, Current)` (`cxx:51-65`): moves `Current` into
/// the period centred on `Prec`.
fn in_good_period(prec: f64, period: f64, current: &mut f64) {
    let mut diff = *current - prec;
    let nb = (diff / period).trunc() as i32 as f64;
    *current -= nb * period;
    diff = *current - prec;
    if diff > period / 2.0 {
        *current -= period;
    } else if diff < -period / 2.0 {
        *current += period;
    }
}

/// The `GeomAbs_Shape` promotion of the `NbIntervals` / `Intervals` switches
/// (`cxx:415-428`, `cxx:440-453`).
fn promoted_shape(s: u8) -> Shape {
    match s {
        0 => Shape::C1,
        2 => Shape::C2,
        4 => Shape::C3,
        _ => Shape::CN,
    }
}
