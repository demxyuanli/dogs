//! `GeomFill_CorrectedFrenet` — the angle-law machinery: `Init`,
//! `InitInterval`, `FindPlane`, `smoothlaw`, `CalcAngleAT` and `GetAngleAT`.
//!
//! Source: `TKGeomAlgo/GeomFill/GeomFill_CorrectedFrenet.cxx:53-1870`
//! (OCCT 8.0.0).
//!
//! Translation notes:
//! - `occ::handle<Law_Function>` members (`EvolAroundT`, `TLaw`) and the OCC78
//!   sample arrays (`HArrPoles` / `HArrAngle` / `HArrTangent` / `HArrNormal`)
//!   are written from the non-const `SetCurve`, which the port calls through a
//!   shared handle, so they live in [`CorrectedState`] behind a `Mutex`.
//! - `BndLib_Add3dCurve::Add(CS, First, Last, 1.e-2, bbox)` (`cxx:460-462`)
//!   is `GeomBndLib_OtherCurve`'s 33-sample envelope, the same construction as
//!   `occt_topo::geom_bnd_lib_curve3d::box_other`; only `GetGap()` is read.
//! - `GeomLib::AxeOfInertia` is [`crate::geom_lib::geom_lib_axe_of_inertia`].

use std::sync::Arc;

use occt_core::gp::{GpAx3, GpPln, GpPnt, GpVec};
use occt_core::kernel::geomabs::Shape;
use occt_core::precision::{CONFUSION, PCONFUSION, RESOLUTION, SQUARE_CONFUSION};

use crate::curve::Curve;
use crate::extrema_pc::{extrema_curve_type, GeomAbsCurveType};
use crate::geom_fill::corrected_frenet::CorrectedFrenet;
use crate::geom_fill::trihedron_law::TrihedronLaw;
use crate::law::composite::LawComposite;
use crate::law::function::{LawConstant, LawFunction};
use crate::law::interpolate::LawInterpolate;
use crate::law::bspline::{LawBSpFunc, LawBSpline};

/// The mutable OCCT members of `GeomFill_CorrectedFrenet`
/// (`GeomFill_CorrectedFrenet.hxx:145-155`).
pub(crate) struct CorrectedState {
    /// `isFrenet` (`hxx:146`).
    pub is_frenet: bool,
    /// `EvolAroundT` (`hxx:147`).
    pub evol_around_t: Option<Arc<LawComposite>>,
    /// `TLaw` (`hxx:148`): `EvolAroundT`, or a trim of it after
    /// `SetInterval`.
    pub t_law: Option<Arc<dyn LawFunction>>,
    /// `AT` (`hxx:149`).
    pub at: GpVec,
    /// `AN` (`hxx:150`).
    pub an: GpVec,
    /// `HArrPoles` (`hxx:152`).
    pub poles: Vec<f64>,
    /// `HArrAngle` (`hxx:153`).
    pub angle: Vec<f64>,
    /// `HArrTangent` (`hxx:154`).
    pub tangent: Vec<GpVec>,
    /// `HArrNormal` (`hxx:155`).
    pub normal: Vec<GpVec>,
}

impl Default for CorrectedState {
    fn default() -> Self {
        Self {
            is_frenet: true,
            evol_around_t: None,
            t_law: None,
            at: GpVec::zero(),
            an: GpVec::zero(),
            poles: Vec::new(),
            angle: Vec::new(),
            tangent: Vec::new(),
            normal: Vec::new(),
        }
    }
}

/// `ComputeTorsion` (`cxx:53-78`).
pub(crate) fn compute_torsion(param: f64, curve: &dyn Curve) -> f64 {
    let (_, dc1, dc2, dc3) = curve.d3(param);
    let cross = dc1.crossed(&dc2);
    let norm = cross.magnitude();
    let mixed = cross.dot(&dc3);
    let tol = RESOLUTION;
    let square_norm = norm * norm;
    if square_norm <= tol {
        0.0
    } else {
        mixed / square_norm
    }
}

/// `smoothlaw` (`cxx:86-187`): reduce the number of knots of an angle law
/// while keeping the sampled values within `tol`.
pub(crate) fn smoothlaw(law: &mut LawBSpline, points: &[f64], params: &[f64], tol_in: f64) {
    let mut bs = law.copy();
    let mut tol = tol_in / 10.0;

    // First pass: tight tolerance (`cxx:95-106`).
    let mut removed = false;
    for ii in (2..bs.nb_knots() as i32).rev() {
        if bs.remove_knot(ii, 0, tol).unwrap_or(false) {
            removed = true;
        }
    }

    let mut use_bs;
    if removed {
        // Control (`cxx:108-129`).
        let mut max_dev: f64 = 0.0;
        let mut ok = true;
        for (point, param) in points.iter().zip(params.iter()) {
            let d = (bs.value(*param) - *point).abs();
            if d > max_dev {
                max_dev = d;
            }
            ok = max_dev <= tol_in;
            if !ok {
                break;
            }
        }
        if !ok {
            return;
        }
        use_bs = true;
        tol = (tol_in - max_dev) / 2.0;
    } else {
        use_bs = false;
        tol = tol_in / 2.0;
    }

    // Second pass: loose tolerance (`cxx:140-150`).
    let mut removed2 = false;
    for ii in (2..bs.nb_knots() as i32).rev() {
        if bs.remove_knot(ii, 0, tol).unwrap_or(false) {
            removed2 = true;
        }
    }

    if removed2 {
        // Control (`cxx:152-172`).
        let mut max_dev: f64 = 0.0;
        let mut ok = true;
        for (point, param) in points.iter().zip(params.iter()) {
            let d = (bs.value(*param) - *point).abs();
            if d > max_dev {
                max_dev = d;
            }
            ok = max_dev <= tol_in;
            if !ok {
                break;
            }
        }
        if ok {
            use_bs = true;
        }
    }

    if use_bs {
        *law = bs;
    }
}

/// `FindPlane` (`cxx:189-302`): the plane the curve lies on, when one is
/// found.
pub(crate) fn find_plane(curve: &dyn Curve) -> Option<GpPln> {
    let mut found = true;
    let mut tab_p: Option<Vec<GpPnt>> = None;
    let mut plane: Option<GpPln> = None;

    match extrema_curve_type(curve) {
        GeomAbsCurveType::Line => {
            found = false;
        }
        GeomAbsCurveType::Circle => {
            if let Some(c) = curve.gp_circ() {
                plane = Some(GpPln::new(GpAx3::from_ax2(&c.position())));
            } else {
                found = false;
            }
        }
        GeomAbsCurveType::Ellipse => {
            if let Some(e) = curve.gp_ellipse() {
                plane = Some(GpPln::new(GpAx3::from_ax2(e.position())));
            } else {
                found = false;
            }
        }
        GeomAbsCurveType::Hyperbola => {
            if let Some(h) = curve.gp_hyperbola() {
                plane = Some(GpPln::new(GpAx3::from_ax2(h.position())));
            } else {
                found = false;
            }
        }
        GeomAbsCurveType::Parabola => {
            if let Some(p) = curve.gp_parabola() {
                plane = Some(GpPln::new(GpAx3::from_ax2(p.position())));
            } else {
                found = false;
            }
        }
        GeomAbsCurveType::BezierCurve => match curve.bezier_poles() {
            Some(poles) if poles.len() > 2 => tab_p = Some(poles.to_vec()),
            _ => found = false,
        },
        GeomAbsCurveType::BSplineCurve => match curve.bspline_poles() {
            Some(poles) if poles.len() > 2 => tab_p = Some(poles.to_vec()),
            _ => found = false,
        },
        _ => {
            // Sampling (`cxx:256-268`).
            let nbp = 15 + curve.nb_intervals(Shape::C3 as u8);
            let f = curve.first_parameter();
            let l = curve.last_parameter();
            let inv = 1.0 / (nbp - 1) as f64;
            let mut samples = Vec::with_capacity(nbp.max(0) as usize);
            for ii in 1..=nbp {
                let t = (f * (nbp - ii) as f64 + l * (ii - 1) as f64) * inv;
                samples.push(curve.d0(t));
            }
            tab_p = Some(samples);
        }
    }

    if let Some(points) = tab_p {
        // Mean plane and control (`cxx:271-298`).
        let (axe, is_singular) = crate::geom_lib::geom_lib_axe_of_inertia(&points, CONFUSION);
        if is_singular {
            found = false;
        } else {
            plane = Some(GpPln::new(GpAx3::from_ax2(&axe)));
        }
        if found {
            match &plane {
                Some(pln) => {
                    let (a, b, c, d) = pln.coefficients();
                    for p in &points {
                        let xyz = p.xyz();
                        let dist = a * xyz.x + b * xyz.y + c * xyz.z + d;
                        if dist.abs() > CONFUSION {
                            found = false;
                            break;
                        }
                    }
                }
                None => found = false,
            }
        }
    }

    if found {
        plane
    } else {
        None
    }
}

/// `corr2PI_PI` (`cxx:637-640`).
fn corr_2pi_pi(ang: f64) -> f64 {
    if ang < std::f64::consts::PI {
        ang
    } else {
        ang - 2.0 * std::f64::consts::PI
    }
}

/// `diffAng` (`cxx:642-647`).
fn diff_ang(a: f64, ao: f64) -> f64 {
    let da = (a - ao) - ((a - ao) / 2.0 / std::f64::consts::PI).floor() * 2.0 * std::f64::consts::PI;
    if da >= 0.0 {
        corr_2pi_pi(da)
    } else {
        -corr_2pi_pi(-da)
    }
}

/// The out-parameters of `InitInterval` (`cxx:435-448`) beyond the law.
pub(crate) struct InitIntervalResult {
    /// `FuncInt`.
    pub func: Arc<dyn LawFunction>,
    /// The `isZero` return value.
    pub is_zero: bool,
}

impl CorrectedFrenet {
    /// `GeomFill_CorrectedFrenet::Init` (`cxx:364-430`).
    pub(crate) fn init(&self) {
        let curve = match self.curve_arc() {
            Some(c) => c,
            None => return,
        };
        let mut composite = LawComposite::new();

        let nb_i = self.frenet.nb_intervals(Shape::C0 as u8);
        let t = self.frenet.intervals(Shape::C0 as u8);
        if (nb_i as usize) + 1 > t.len() {
            return;
        }

        let mut seq_poles: Vec<f64> = Vec::new();
        let mut seq_angle: Vec<f64> = Vec::new();
        let mut seq_tangent: Vec<GpVec> = Vec::new();
        let mut seq_normal: Vec<GpVec> = Vec::new();

        let (first, last) = self.data.interval();
        let (mut tangent, mut normal, _bn) = match self.frenet.d0(first) {
            Ok(f) => f,
            Err(_) => return,
        };
        let mut start_ang = 0.0;
        let mut at = self.state.lock().unwrap_or_else(|e| e.into_inner()).at;
        let mut an = self.state.lock().unwrap_or_else(|e| e.into_inner()).an;

        const NB_STEP: f64 = 10.0;
        let av_step = (last - first) / NB_STEP;
        let mut is_frenet = self.is_frenet();

        for i in 1..=nb_i as usize {
            let span = t[i] - t[i - 1];
            let nb_step = ((span / av_step) as i32).max(3);
            let step = span / nb_step as f64;
            let result = self.init_interval(
                t[i - 1],
                t[i],
                step,
                &mut start_ang,
                &mut tangent,
                &mut normal,
                &mut at,
                &mut an,
                &curve,
                &mut seq_poles,
                &mut seq_angle,
                &mut seq_tangent,
                &mut seq_normal,
            );
            let func = match result {
                Ok(r) => {
                    if !r.is_zero {
                        is_frenet = false;
                    }
                    r.func
                }
                Err(_) => {
                    is_frenet = false;
                    Arc::new(LawConstant::new()) as Arc<dyn LawFunction>
                }
            };
            composite.add_law(func);
        }
        if curve.is_periodic() {
            composite.set_periodic();
        }

        let composite = Arc::new(composite);
        let mut state = self.state.lock().unwrap_or_else(|e| e.into_inner());
        state.evol_around_t = Some(composite.clone());
        state.t_law = Some(composite as Arc<dyn LawFunction>);
        state.at = at;
        state.an = an;
        let i_end = seq_poles.len();
        if i_end != 0 {
            state.poles = seq_poles;
            state.angle = seq_angle;
            state.tangent = seq_tangent;
            state.normal = seq_normal;
        } else {
            state.poles.clear();
            state.angle.clear();
            state.tangent.clear();
            state.normal.clear();
        }
        drop(state);
        self.set_is_frenet(is_frenet);
    }

    /// `GeomFill_CorrectedFrenet::InitInterval` (`cxx:435-597`): the angle
    /// law of one interval of continuity of the Frenet law.
    #[allow(clippy::too_many_arguments)]
    fn init_interval(
        &self,
        first: f64,
        last: f64,
        step: f64,
        start_ang: &mut f64,
        prev_tangent: &mut GpVec,
        prev_normal: &mut GpVec,
        a_t: &mut GpVec,
        a_n: &mut GpVec,
        curve: &Arc<dyn Curve>,
        seq_poles: &mut Vec<f64>,
        seq_angle: &mut Vec<f64>,
        seq_tangent: &mut Vec<GpVec>,
        seq_normal: &mut Vec<GpVec>,
    ) -> Result<InitIntervalResult, String> {
        let mut parameters: Vec<f64> = Vec::new();
        let mut evol_at: Vec<f64> = Vec::new();
        let mut param = first;
        let mut is_zero = true;
        let mut is_const = true;

        // `frenet->SetInterval(First, Last)` (`cxx:458`).
        self.frenet.set_interval(first, last)?;

        let snglr = crate::geom_fill::snglr_func::SnglrFunc::new(curve.clone());
        let length_min = snglr_gap(&snglr, first, last) * 1.0e-4;

        *a_t = GpVec::zero();
        *a_n = GpVec::zero();

        let mut angle_at: f64 = 0.0;
        let mut curr_step = step;

        let is_planar = if self.for_evaluation() {
            false
        } else {
            find_plane(curve.as_ref()).is_some()
        };

        let mut i = 1;
        let mut curr_param = param;
        let d_last = last - PCONFUSION;

        while param < last {
            if curr_param > d_last {
                if (d_last - param).abs() < SQUARE_CONFUSION {
                    param = curr_param;
                }
                curr_step = d_last - param;
                curr_param = last;
            }
            if is_planar {
                curr_param = last;
            }

            let (tangent, normal, _bn) = self.frenet.d0(curr_param)?;
            if prev_tangent.angle(&tangent) < std::f64::consts::PI / 3.0 || i == 1 {
                parameters.push(curr_param);
                // OCC78 (`cxx:503-508`).
                seq_poles.push(param);
                seq_angle.push(if i > 1 { evol_at[i - 2] } else { *start_ang });
                seq_tangent.push(*prev_tangent);
                seq_normal.push(*prev_normal);
                angle_at = self.calc_angle_at(&tangent, &normal, prev_tangent, prev_normal);

                if is_const && i > 1 && angle_at.abs() > PCONFUSION {
                    is_const = false;
                }

                angle_at += if i > 1 { evol_at[i - 2] } else { *start_ang };
                evol_at.push(angle_at);
                *prev_normal = normal;

                if is_zero && angle_at.abs() > PCONFUSION {
                    is_zero = false;
                }

                *a_t = a_t.added(&tangent);

                let cross = tangent.crossed(&normal);
                let tcross = tangent.crossed(&cross);
                *a_n = cross
                    .multiplied_scalar(angle_at.sin())
                    .added(&tcross.multiplied_scalar(1.0 - angle_at.cos()))
                    .added(&normal)
                    .added(&*a_n);
                *prev_tangent = tangent;
                param = curr_param;
                i += 1;

                // Evaluate the next step (`cxx:542-552`).
                let (pon_c, d1) = snglr.d1(param);
                let l = (pon_c.xyz().modulus() / 2.0).max(length_min);
                let mut norm = d1.magnitude();
                if norm < CONFUSION {
                    norm = CONFUSION;
                }
                curr_step = l / norm;
                if curr_step > step {
                    curr_step = step;
                }
            } else {
                curr_step /= 2.0;
            }

            curr_param = param + curr_step;
        }

        if !is_planar {
            let n = (parameters.len() - 1) as f64;
            *a_t = a_t.divided(n);
            *a_n = a_n.divided(n);
        }
        *start_ang = angle_at;

        // Interpolation (`cxx:569-596`).
        if is_const || is_planar {
            let mut law = LawConstant::new();
            law.set(angle_at, first, last);
            Ok(InitIntervalResult {
                func: Arc::new(law) as Arc<dyn LawFunction>,
                is_zero,
            })
        } else {
            let mut law_at = LawInterpolate::new_with_parameters(
                evol_at.clone(),
                parameters.clone(),
                false,
                PCONFUSION,
            )?;
            law_at.perform()?;
            let mut bs = law_at.curve()?.copy();
            smoothlaw(&mut bs, &evol_at, &parameters, 0.1);
            Ok(InitIntervalResult {
                func: Arc::new(LawBSpFunc::new_with_curve(bs, first, last)) as Arc<dyn LawFunction>,
                is_zero,
            })
        }
    }

    /// `GeomFill_CorrectedFrenet::CalcAngleAT` (`cxx:607-632`): the angle
    /// between the normal of `Tangent` and the normal of `prevTangent`.
    pub(crate) fn calc_angle_at(
        &self,
        tangent: &GpVec,
        normal: &GpVec,
        prev_tangent: &GpVec,
        prev_normal: &GpVec,
    ) -> f64 {
        let angle = tangent.angle(prev_tangent);
        let normal_rot = if angle.abs() > occt_core::precision::ANGULAR
            && angle.abs() < std::f64::consts::PI - occt_core::precision::ANGULAR
        {
            let cross = tangent.crossed(prev_tangent).normalized();
            normal
                .added(&cross.crossed(normal).multiplied_scalar(angle.sin()))
                .added(&cross.crossed(&cross.crossed(normal)).multiplied_scalar(1.0 - angle.cos()))
        } else {
            *normal
        };
        let mut angle_at = normal_rot.angle(prev_normal);
        if angle_at > occt_core::precision::ANGULAR
            && std::f64::consts::PI - angle_at > occt_core::precision::ANGULAR
        {
            if normal_rot
                .crossed(prev_normal)
                .is_opposite_ang(prev_tangent, occt_core::precision::ANGULAR)
            {
                angle_at = -angle_at;
            }
        }
        angle_at
    }

    /// `GeomFill_CorrectedFrenet::GetAngleAT` (`cxx:654-696`): the corrected
    /// angle at `param`, fixing the OCC78 discontinuity.
    pub(crate) fn get_angle_at(&self, param: f64) -> f64 {
        let state = self.state.lock().unwrap_or_else(|e| e.into_inner());
        let Some(t_law) = state.t_law.clone() else {
            return 0.0;
        };
        let poles = &state.poles;
        if poles.is_empty() {
            return t_law.value(param);
        }
        let len = poles.len() as i32;
        let mut i_b = 1i32;
        let mut i_e = len;
        let mut i_c = (i_e + i_b) / 2;
        if param == poles[(i_b - 1) as usize] {
            return t_law.value(param);
        }
        if param > poles[(i_e - 1) as usize] {
            i_c = i_e;
        }
        if i_c < i_e {
            while i_c >= 1
                && (poles[(i_c - 1) as usize] > param
                    || (i_c < len && param > poles[i_c as usize]))
            {
                if poles[(i_c - 1) as usize] < param {
                    i_b = i_c;
                } else {
                    i_e = i_c;
                }
                i_c = (i_e + i_b) / 2;
            }
            if poles[(i_c - 1) as usize] == param
                || (i_c < len && param == poles[i_c as usize])
            {
                return t_law.value(param);
            }
        }
        let ang_p = t_law.value(param);
        let ang_po = state.angle[(i_c - 1) as usize];
        let d_ang = ang_p - ang_po;
        let (tangent, normal, _bn) = match self.frenet.d0(param) {
            Ok(f) => f,
            Err(_) => return ang_p,
        };
        let d_ang_t = self.calc_angle_at(
            &tangent,
            &normal,
            &state.tangent[(i_c - 1) as usize],
            &state.normal[(i_c - 1) as usize],
        );
        let da = diff_ang(d_ang_t, d_ang);
        if da.abs() > std::f64::consts::PI / 2.0 {
            ang_po + d_ang_t
        } else {
            ang_p
        }
    }
}

/// `BndLib_Add3dCurve::Add(CS, First, Last, 1.e-2, bbox)` reduced to
/// `GetGap()`: `GeomBndLib_OtherCurve`'s 33-sample envelope
/// (`GeomBndLib_OtherCurve.cxx:72-87`) enlarged by the tolerance.
fn snglr_gap(curve: &dyn Curve, first: f64, last: f64) -> f64 {
    let mut points: Vec<GpPnt> = Vec::with_capacity(67);
    points.push(curve.d0(first));
    let dp = last - first;
    if dp.abs() > PCONFUSION {
        let n = 33;
        let step = dp / (2 * n) as f64;
        let mut p = first;
        for _ in 1..=n {
            p += step;
            points.push(curve.d0(p));
            p += step;
            points.push(curve.d0(p));
        }
    } else {
        points.push(curve.d0(last));
    }
    let mut bb = occt_core::bnd::BndBox::new();
    for p in &points {
        bb.add_point(p);
    }
    bb.enlarge(1.0e-2);
    bb.gap()
}
