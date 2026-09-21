//! `Extrema_ExtPC` / `Extrema_GGExtPC::Perform` + `IntervalPerform`.
//!
//! Source: `Extrema_GGExtPC.hxx` (Perform ~132-528, IntervalPerform ~590-614,
//! AddSol ~617-632), `Extrema_GGenExtPC.hxx` (Perform ~160-173),
//! `Extrema_GFuncExtPC.hxx` (Value/Values/GetStateNumber),
//! `Extrema_CurveTool::DeflCurvIntervals` (`Extrema_CurveTool.cxx:41-91`).
//!
//! Used by `ShapeAnalysis_Curve::ProjectAct` (`ShapeAnalysis_Curve.cxx:277`).
//! Singular DN / MaxDerivOrder: `GFuncExtPC::Value` (`hxx:167-241`) and
//! `Values` three-point DF (`hxx:297-338`).

use occt_core::bspl::knots::unique_knots_mults;
use occt_core::bspl::locate::{first_u_knot_index, last_u_knot_index};
use occt_core::gcpnts::{perform_tangential_curve, CurveSecondDeriv};
use occt_core::gp::{GpPnt, GpVec};
use occt_core::math_fn::{MathFunction, MathFunctionWithDerivative};
use occt_core::math_function_roots::FunctionRoots;
use occt_core::precision::{CONFUSION, PCONFUSION, Precision};

use crate::curve::Curve;

const EPS_DL_DT: f64 = 1.0e-3;
const MAX_DEFL: f64 = 1.0e3;
const MIN_DEFL: f64 = 1.0e-3;
const DEFL_LENGTH_SAMPLES: i32 = 23;

const TOL_F: f64 = 1.0e-10;
const MIN_TOL: f64 = 1.0e-20;
const TOL_FACTOR: f64 = 1.0e-12;
const MIN_STEP: f64 = 1.0e-7;
const MAX_ORDER: i32 = 3;
const MAX_SAMPLE: i32 = 17;

/// One `Extrema_ExtPC` solution (`SquareDistance` / `IsMin` / `Point`).
#[derive(Debug, Clone, Copy)]
pub struct ExtPcSolution {
    pub u: f64,
    pub point: GpPnt,
    pub sq_dist: f64,
    pub is_min: bool,
}

/// `Extrema_GFuncExtPC` for 3D point–curve (`Extrema_PCFOfEPCOfExtPC`).
struct ExtPcf<'a> {
    curve: &'a dyn Curve,
    point: GpPnt,
    u: f64,
    pc: GpPnt,
    d1f: f64,
    u_inf: f64,
    u_sup: f64,
    /// `myMaxDerivOrder` (`hxx:67` / Initialize GetType switch).
    max_deriv_order: i32,
    /// `myTol` — SearchOfTolerance or MinTol (`hxx:92-98`).
    tol: f64,
    sols: Vec<ExtPcSolution>,
}

/// `GFuncExtPC` Initialize GetType arm (`hxx:87-99`, `115-128`).
fn gfunc_max_deriv_order(curve: &dyn Curve) -> i32 {
    if curve.bspline_poles().is_some()
        || curve.bspline_knots().is_some()
        || curve.bezier_poles().is_some()
    {
        return MAX_ORDER;
    }
    // GeomAbs_OffsetCurve / GeomAbs_OtherCurve (no elementary line/circle).
    if curve.is_line() || curve.gp_circ().is_some() || curve.circle_radius().is_some() {
        0
    } else {
        MAX_ORDER
    }
}

/// `GFuncExtPC::SearchOfTolerance` (`hxx:428-457`) on curve First/Last.
fn gfunc_search_of_tolerance(curve: &dyn Curve) -> f64 {
    const N_POINT: i32 = 10;
    let u_inf = curve.first_parameter();
    let u_sup = curve.last_parameter();
    let a_step = (u_sup - u_inf) / (N_POINT as f64);
    let mut a_num = 0;
    let mut a_max = -Precision::INFINITE;
    loop {
        let mut u = u_inf + (a_num as f64) * a_step;
        if u > u_sup {
            u = u_sup;
        }
        let (_p, v_der) = curve.d1(u);
        if !(Precision::is_infinite(v_der.x()) || Precision::is_infinite(v_der.y())) {
            let vm = v_der.magnitude();
            if vm > a_max {
                a_max = vm;
            }
        }
        a_num += 1;
        if a_num >= N_POINT + 1 {
            break;
        }
    }
    (a_max * TOL_FACTOR).max(MIN_TOL)
}

impl<'a> ExtPcf<'a> {
    fn new(curve: &'a dyn Curve, point: &GpPnt, u_inf: f64, u_sup: f64) -> Self {
        let max_deriv_order = gfunc_max_deriv_order(curve);
        let tol = if max_deriv_order != 0 {
            gfunc_search_of_tolerance(curve)
        } else {
            MIN_TOL
        };
        Self {
            curve,
            point: *point,
            u: 0.0,
            pc: GpPnt::zero(),
            d1f: 0.0,
            u_inf,
            u_sup,
            max_deriv_order,
            tol,
            sols: Vec::new(),
        }
    }

    fn clear_sols(&mut self) {
        self.sols.clear();
    }

    fn sub_interval_initialize(&mut self, u_first: f64, u_last: f64) {
        self.u_inf = u_first;
        self.u_sup = u_last;
    }

    /// Subinterval span for singular delta (`hxx:170-176`, `301-307`).
    fn singular_du(&self) -> f64 {
        if !self.u_sup.is_finite() || !self.u_inf.is_finite() {
            0.0
        } else {
            self.u_sup - self.u_inf
        }
    }

    /// `GFuncExtPC::Value` (`hxx:146-253`) including MaxDerivOrder singular D1.
    fn value_at(&mut self, the_u: f64, the_f: &mut f64) -> bool {
        self.u = the_u;
        let (pc, mut d1c) = self.curve.d1(the_u);
        self.pc = pc;
        if Precision::is_infinite(d1c.x()) || Precision::is_infinite(d1c.y()) {
            *the_f = f64::INFINITY;
            return false;
        }
        let mut ndu = d1c.magnitude();

        if self.max_deriv_order != 0 && ndu <= self.tol {
            // Singular case (`hxx:167-241`).
            const DIVISION_FACTOR: f64 = 1.0e-3;
            let a_delta = (self.singular_du() * DIVISION_FACTOR).max(MIN_STEP);

            let mut n = 1;
            let mut v = GpVec::zero();
            let mut is_derive_found = false;
            while !is_derive_found && n < self.max_deriv_order {
                n += 1;
                v = self.curve.eval_dn(the_u, n);
                ndu = v.magnitude();
                is_derive_found = ndu > self.tol;
            }

            if is_derive_found {
                let u = if self.u - self.u_inf < a_delta {
                    self.u + a_delta
                } else {
                    self.u - a_delta
                };
                let p1 = self.curve.d0(self.u.min(u));
                let p2 = self.curve.d0(self.u.max(u));
                let v1 = GpVec::from_pnts(&p1, &p2);
                if v.dot(&v1) < 0.0 {
                    d1c = v.reversed();
                } else {
                    d1c = v;
                }
            } else {
                // Three-point Taylor (`hxx:213-238`); Ptemp default origin.
                let (p1, p2, p3, grown) = if self.u - self.u_inf < 2.0 * a_delta {
                    (
                        self.curve.d0(self.u),
                        self.curve.d0(self.u + a_delta),
                        self.curve.d0(self.u + 2.0 * a_delta),
                        true,
                    )
                } else {
                    (
                        self.curve.d0(self.u - 2.0 * a_delta),
                        self.curve.d0(self.u - a_delta),
                        self.curve.d0(self.u),
                        false,
                    )
                };
                let v1 = GpVec::new(p1.x(), p1.y(), p1.z());
                let v2 = GpVec::new(p2.x(), p2.y(), p2.z());
                let v3 = GpVec::new(p3.x(), p3.y(), p3.z());
                d1c = if grown {
                    v1.multiplied_scalar(-3.0)
                        .added(&v2.multiplied_scalar(4.0))
                        .subtracted(&v3)
                } else {
                    v1.subtracted(&v2.multiplied_scalar(4.0))
                        .added(&v3.multiplied_scalar(3.0))
                };
            }
            ndu = d1c.magnitude();
        }

        if ndu <= MIN_TOL {
            return false;
        }
        let ppc = GpVec::from_pnts(&self.point, &self.pc);
        *the_f = ppc.dot(&d1c) / ndu;
        true
    }

    /// `GFuncExtPC::Values` (`hxx:274-352`) including singular DF (`hxx:297-338`).
    fn values_at(&mut self, the_u: f64, the_f: &mut f64, the_df: &mut f64) -> bool {
        let pc_old = self.pc;
        let p_old = self.point;
        if !self.value_at(the_u, the_f) {
            return false;
        }
        self.u = the_u;
        self.pc = pc_old;
        self.point = p_old;
        let (pc, d1c, d2c) = self.curve.d2(the_u);
        self.pc = pc;
        let ndu = d1c.magnitude();
        if ndu <= self.tol {
            // Singular DF via three F samples (`hxx:297-338`).
            const DIVISION_FACTOR: f64 = 0.01;
            let a_delta = (self.singular_du() * DIVISION_FACTOR).max(MIN_STEP);
            let f1;
            let f2;
            let f3;
            if self.u - self.u_inf < 2.0 * a_delta {
                f1 = *the_f;
                let u2 = self.u + a_delta;
                let u3 = self.u + a_delta * 2.0;
                let mut f2v = 0.0;
                let mut f3v = 0.0;
                if !(self.value_at(u2, &mut f2v) && self.value_at(u3, &mut f3v)) {
                    return false;
                }
                f2 = f2v;
                f3 = f3v;
                *the_df = (-3.0 * f1 + 4.0 * f2 - f3) / (2.0 * a_delta);
            } else {
                f3 = *the_f;
                let u1 = self.u - a_delta * 2.0;
                let u2 = self.u - a_delta;
                let mut f1v = 0.0;
                let mut f2v = 0.0;
                if !(self.value_at(u2, &mut f2v) && self.value_at(u1, &mut f1v)) {
                    return false;
                }
                f1 = f1v;
                f2 = f2v;
                *the_df = (f1 - 4.0 * f2 + 3.0 * f3) / (2.0 * a_delta);
            }
            self.u = the_u;
            self.pc = pc_old;
            self.point = p_old;
        } else {
            let ppc = GpVec::from_pnts(&self.point, &self.pc);
            *the_df = ndu + (ppc.dot(&d2c) / ndu) - *the_f * (d1c.dot(&d2c)) / (ndu * ndu);
        }
        self.d1f = *the_df;
        true
    }
}

impl MathFunction for ExtPcf<'_> {
    fn value(&mut self, x: f64, f: &mut f64) -> bool {
        self.value_at(x, f)
    }

    fn get_state_number(&mut self) -> i32 {
        // `GFuncExtPC::GetStateNumber` (`hxx:357-378`).
        let sq = self.pc.square_distance(&self.point);
        let mut ff = 0.0;
        let mut dd = 0.0;
        let _ = self.values_at(self.u, &mut ff, &mut dd);
        let is_min = self.d1f > 0.0;
        self.sols.push(ExtPcSolution {
            u: self.u,
            point: self.pc,
            sq_dist: sq,
            is_min,
        });
        0
    }
}

impl MathFunctionWithDerivative for ExtPcf<'_> {
    fn derivative(&mut self, x: f64, d: &mut f64) -> bool {
        let mut f = 0.0;
        self.values_at(x, &mut f, d)
    }

    fn values(&mut self, x: f64, f: &mut f64, d: &mut f64) -> bool {
        self.values_at(x, f, d)
    }
}

struct ExtPcState {
    sols: Vec<ExtPcSolution>,
    done: bool,
    dist1: f64,
    dist2: f64,
    pf: GpPnt,
    pl: GpPnt,
    tol_u: f64,
    u_inf: f64,
    u_sup: f64,
    sample: i32,
    int_u_inf: f64,
    int_u_sup: f64,
}

impl ExtPcState {
    fn new(curve: &dyn Curve, u_inf: f64, u_sup: f64) -> Self {
        Self {
            sols: Vec::new(),
            done: false,
            dist1: f64::MAX,
            dist2: f64::MAX,
            pf: GpPnt::zero(),
            pl: GpPnt::zero(),
            tol_u: curve.resolution(CONFUSION),
            u_inf,
            u_sup,
            sample: 17,
            int_u_inf: u_inf,
            int_u_sup: u_sup,
        }
    }

    /// `Extrema_GGExtPC::AddSol` (`hxx:617-632`).
    fn add_sol(&mut self, u: f64, p: GpPnt, sq_dist: f64, is_min: bool) {
        for s in &self.sols {
            if (s.u - u).abs() <= self.tol_u {
                return;
            }
        }
        self.sols.push(ExtPcSolution {
            u,
            point: p,
            sq_dist,
            is_min,
        });
    }

    /// `IntervalPerform` (`hxx:590-614`) via `math_FunctionRoots` on `GFuncExtPC`.
    fn interval_perform(&mut self, curve: &dyn Curve, point: &GpPnt) {
        let mut f = ExtPcf::new(curve, point, self.int_u_inf, self.int_u_sup);
        f.sub_interval_initialize(self.int_u_inf, self.int_u_sup);
        f.clear_sols();
        let roots = FunctionRoots::new(
            &mut f,
            self.int_u_inf,
            self.int_u_sup,
            self.sample,
            self.tol_u,
            TOL_F,
            TOL_F,
            0.0,
        );
        self.done = roots.is_done() && !roots.is_all_null();
        if !self.done {
            return;
        }
        let period = if curve.is_periodic() {
            curve.period()
        } else {
            0.0
        };
        for s in f.sols {
            let mut u = s.u;
            if period.abs() > 0.0 {
                u = elclib_in_period(u, self.u_inf, self.u_inf + period);
            }
            if u >= self.u_inf - self.tol_u && u <= self.u_sup + self.tol_u {
                self.add_sol(u, s.point, s.sq_dist, s.is_min);
            }
        }
    }
}

/// `ElCLib::InPeriod` (`ElCLib.cxx`).
fn elclib_in_period(u: f64, u_first: f64, u_last: f64) -> f64 {
    let mut period = u_last - u_first;
    if period <= 0.0 {
        return u;
    }
    let mut x = u;
    while x < u_first {
        x += period;
    }
    while x > u_last || (x - u_last).abs() <= PCONFUSION {
        if (x - u_last).abs() <= PCONFUSION {
            return u_last;
        }
        x -= period;
        if x < u_first {
            break;
        }
    }
    if x < u_first {
        x += period;
    }
    x
}

fn perform_bspline(state: &mut ExtPcState, curve: &dyn Curve, point: &GpPnt) {
    let flat = match curve.bspline_knots() {
        Some(k) if !k.is_empty() => k,
        _ => {
            perform_general(state, curve, point);
            return;
        }
    };
    let degree = curve.nurbs_degree().unwrap_or(1) as i32;
    let (knots, mults) = unique_knots_mults(flat);
    if knots.len() < 2 {
        perform_general(state, curve, point);
        return;
    }
    let first_idx = first_u_knot_index(degree, &mults);
    let last_idx = last_u_knot_index(degree, &mults);
    let a_tol_coeff = (state.u_sup - state.u_inf) * PCONFUSION;
    let mut period_jump = 0.0;
    if curve.is_periodic() {
        let period = curve.period();
        if period.abs() > 0.0 {
            let mut shift = ((state.u_inf - knots[(first_idx - 1) as usize]) / period) as i32;
            if state.u_inf < knots[(first_idx - 1) as usize] - a_tol_coeff {
                shift -= 1;
            }
            period_jump = period * shift as f64;
        }
    }

    let mut first_used = first_idx;
    let mut last_used = last_idx;
    for an_idx in first_idx..=last_idx {
        let a_knot = knots[(an_idx - 1) as usize] + period_jump;
        if state.u_inf >= a_knot - a_tol_coeff {
            first_used = an_idx;
        } else {
            break;
        }
    }
    for an_idx in (first_idx..=last_idx).rev() {
        let a_knot = knots[(an_idx - 1) as usize] + period_jump;
        if state.u_sup <= a_knot + a_tol_coeff {
            last_used = an_idx;
        } else {
            break;
        }
    }
    if first_used == last_used {
        first_used = first_idx;
        last_used = first_idx + 1;
    }

    state.sample = degree + 1;

    if state.sample == 2 {
        // Degree-1 knot-span arm (`GGExtPC.hxx:233-297`).
        let mut tmin = 0.0;
        let mut distmin = f64::MAX;
        let mut pmin = GpPnt::zero();
        let mut a_min2 = 0.0;
        for an_idx in first_used..last_used {
            let mut a_f = knots[(an_idx - 1) as usize] + period_jump;
            let mut a_l = knots[an_idx as usize] + period_jump;
            if an_idx == first_used {
                a_f = state.u_inf;
            } else if an_idx == last_used - 1 {
                a_l = state.u_sup;
            }
            let a_p1 = curve.d0(a_f);
            let a_p2 = curve.d0(a_l);
            let a_base1 = GpVec::from_pnts(point, &a_p1);
            let a_base2 = GpVec::from_pnts(point, &a_p2);
            let a_v = GpVec::from_pnts(&a_p2, &a_p1);
            let a_val1 = a_v.dot(&a_base1);
            let a_val2 = a_v.dot(&a_base2);
            let a_min1 = if an_idx == first_used {
                point.square_distance(&a_p1)
            } else {
                let prev = a_min2;
                if distmin > prev {
                    distmin = prev;
                    tmin = a_f;
                    pmin = a_p1;
                }
                prev
            };
            a_min2 = point.square_distance(&a_p2);
            let a_min_sq = a_min1.min(a_min2);
            let a_min_der = a_val1.abs().min(a_val2.abs());
            if !(Precision::is_infinite(a_val1) || Precision::is_infinite(a_val2))
                && (a_val1 * a_val2 <= 0.0
                    || a_min_sq < 100.0 * CONFUSION * CONFUSION
                    || 2.0 * a_min_der < CONFUSION)
            {
                state.int_u_inf = a_f;
                state.int_u_sup = a_l;
                state.interval_perform(curve, point);
            }
        }
        if distmin.is_finite() {
            let mut is_to_add = true;
            for s in &state.sols {
                if !(distmin < s.sq_dist && (s.u - tmin).abs() > state.tol_u) {
                    is_to_add = false;
                    break;
                }
            }
            if is_to_add {
                state.add_sol(tmin, pmin, distmin, true);
            }
        }
    } else {
        // Degree>=2 sample + IntervalPerform (`GGExtPC.hxx:300-385`).
        let n_spans = (last_used - first_used).max(1) as usize;
        let cap = (state.sample as usize) * n_spans + 1;
        let mut vals = Vec::with_capacity(cap);
        let mut params = Vec::with_capacity(cap);
        for an_idx in first_used..last_used {
            let mut a_f = knots[(an_idx - 1) as usize] + period_jump;
            let mut a_l = knots[an_idx as usize] + period_jump;
            if an_idx == first_used {
                a_f = state.u_inf;
            }
            if an_idx == last_used - 1 {
                a_l = state.u_sup;
            }
            let a_step = (a_l - a_f) / state.sample as f64;
            for a_pnt_idx in 0..state.sample {
                let a_cur = a_f + a_step * a_pnt_idx as f64;
                vals.push(curve.d0(a_cur).square_distance(point));
                params.push(a_cur);
            }
        }
        vals.push(curve.d0(state.u_sup).square_distance(point));
        params.push(state.u_sup);

        for an_idx in 1..vals.len().saturating_sub(1) {
            if vals[an_idx] <= CONFUSION * CONFUSION {
                state.add_sol(
                    params[an_idx],
                    curve.d0(params[an_idx]),
                    vals[an_idx],
                    true,
                );
            }
            if (vals[an_idx] >= vals[an_idx + 1] && vals[an_idx] >= vals[an_idx - 1])
                || (vals[an_idx] <= vals[an_idx + 1] && vals[an_idx] <= vals[an_idx - 1])
            {
                state.int_u_inf = params[an_idx - 1];
                state.int_u_sup = params[an_idx + 1];
                state.interval_perform(curve, point);
            }
        }

        if state.dist1 > CONFUSION * CONFUSION && state.dist1.is_finite() && params.len() >= 2 {
            let (p1, v1) = curve.d1(params[0]);
            let (p2, v2) = curve.d1(params[1]);
            let a_base1 = GpVec::from_pnts(point, &p1);
            let a_base2 = GpVec::from_pnts(point, &p2);
            let a_val1 = v1.dot(&a_base1);
            let a_val2 = v2.dot(&a_base2);
            if !(Precision::is_infinite(a_val1) || Precision::is_infinite(a_val2))
                && (a_val1 * a_val2 <= 0.0
                    || a_base1.dot(&a_base2) <= 0.0
                    || 2.0 * a_val1.abs() < CONFUSION)
            {
                state.int_u_inf = params[0];
                state.int_u_sup = params[1];
                state.interval_perform(curve, point);
            }
        }
        if state.dist2 > CONFUSION * CONFUSION && state.dist2.is_finite() && params.len() >= 2 {
            let n = params.len();
            let (p1, v1) = curve.d1(params[n - 2]);
            let (p2, v2) = curve.d1(params[n - 1]);
            let a_base1 = GpVec::from_pnts(point, &p1);
            let a_base2 = GpVec::from_pnts(point, &p2);
            let a_val1 = v1.dot(&a_base1);
            let a_val2 = v2.dot(&a_base2);
            if !(Precision::is_infinite(a_val1) || Precision::is_infinite(a_val2))
                && (a_val1 * a_val2 <= 0.0
                    || a_base1.dot(&a_base2) <= 0.0
                    || 2.0 * a_val2.abs() < CONFUSION)
            {
                state.int_u_inf = params[n - 2];
                state.int_u_sup = params[n - 1];
                state.interval_perform(curve, point);
            }
        }
    }
    state.done = true;
}

/// Adapter for `GCPnts_TangentialDeflection::PerformCurve`.
struct CurveTd<'a>(&'a dyn Curve);

impl CurveSecondDeriv for CurveTd<'_> {
    fn point(&self, u: f64) -> GpPnt {
        self.0.d0(u)
    }
    fn d2(&self, u: f64) -> (GpPnt, GpVec, GpVec) {
        self.0.d2(u)
    }
}

/// `Extrema_CurveTool::DeflCurvIntervals` (`Extrema_CurveTool.cxx:41-91`).
///
/// Builds curvature-deflection sample breakpoints via
/// `GCPnts_TangentialDeflection` when a single C2 span would otherwise leave
/// `[First, Last]` as the only ExtPC interval.
fn defl_curv_intervals(curve: &dyn Curve) -> Vec<f64> {
    let tf = curve.first_parameter();
    let tl = curve.last_parameter();
    if !(tf.is_finite() && tl.is_finite()) || (tl - tf).abs() <= 0.0 {
        return vec![tf, tl];
    }

    // Chord-sum from First to each of 22 interior/end samples (`cxx:51-57`).
    // OCCT does not update the running point — distances are from `Value(tf)`.
    let a_p = curve.d0(tf);
    let mut length = 0.0;
    let nbpnts = DEFL_LENGTH_SAMPLES;
    for i in 2..=nbpnts {
        let t = (tf * (nbpnts - i) as f64 + (i - 1) as f64 * tl) / (nbpnts - 1) as f64;
        length += a_p.distance(&curve.d0(t));
    }

    let span = tl - tf;
    let dldt = length / span;
    if length <= CONFUSION || dldt < EPS_DL_DT || span > 10000.0 {
        return vec![tf, tl];
    }

    let a_defl = (0.01 * length / (2.0 * std::f64::consts::PI)).max(MIN_DEFL);
    if a_defl > MAX_DEFL {
        return vec![tf, tl];
    }
    let a_min_len = (0.00001 * length).max(CONFUSION);
    let a_tol = (0.00001 * span).max(PCONFUSION);

    // `GCPnts_TangentialDeflection(C, PI/6, aDefl, 2, aTol, aMinLen)` then
    // `PerformCurve` CN intervals (`TangentialDeflection.cxx:554-556`).
    let cn = {
        let mut iv = curve.parameter_intervals(6);
        if iv.len() < 2 {
            iv = vec![tf, tl];
        }
        iv
    };
    // Offset / OtherCurve stay at MinimumOfPoints=2; BSpline/Bezier bump is
    // inside PerformCurve when GetType matches (`cxx:568-578`).
    let degree_min_nb = curve
        .nurbs_degree()
        .map(|d| (d as usize + 1).max(2))
        .unwrap_or(2);

    let (params, _) = perform_tangential_curve(
        &CurveTd(curve),
        tf,
        tl,
        std::f64::consts::PI / 6.0,
        a_defl,
        2,
        a_tol,
        a_min_len,
        &cn,
        degree_min_nb,
    );
    if params.len() < 2 {
        vec![tf, tl]
    } else {
        params
    }
}

/// Default / OtherCurve / Offset path (`GGExtPC.hxx:390-469`).
/// Bezier early-return is `hxx:182-188`.
fn perform_general(state: &mut ExtPcState, curve: &dyn Curve, point: &GpPnt) {
    if let Some(poles) = curve.bezier_poles() {
        state.int_u_inf = state.u_inf;
        state.int_u_sup = state.u_sup;
        state.sample = (poles.len() as i32) * 2;
        state.interval_perform(curve, point);
        return;
    }

    // `GGExtPC.hxx:395-405`: C2 intervals when NbIntervals>1, else DeflCurvIntervals.
    let c2 = curve.parameter_intervals(2);
    let n_c2 = c2.len().saturating_sub(1);
    let intervals = if n_c2 > 1 {
        c2
    } else {
        defl_curv_intervals(curve)
    };
    let n = intervals.len().saturating_sub(1).max(1);
    let mut maxint = 0.0;
    for i in 0..n {
        let dt = intervals[i + 1] - intervals[i];
        if maxint < dt {
            maxint = dt;
        }
    }
    if maxint <= 0.0 {
        maxint = (state.u_sup - state.u_inf).abs().max(CONFUSION);
    }
    let is_periodic = curve.is_periodic();
    let mut int_ext_done = false;
    let mut s2 = 0.0;
    for i in 0..n {
        state.int_u_inf = intervals[i];
        state.int_u_sup = intervals[i + 1];
        // `RealToInt` truncates toward zero (`GGExtPC.hxx:424`).
        state.sample = ((MAX_SAMPLE as f64) * (state.int_u_sup - state.int_u_inf) / maxint) as i32;
        state.sample = state.sample.max(3);

        let mut an_inf = state.int_u_inf;
        let mut a_sup = state.int_u_sup;
        if is_periodic {
            let period = curve.period();
            an_inf = elclib_in_period(state.int_u_inf, state.u_inf, state.u_inf + period);
            a_sup = state.int_u_sup + (an_inf - state.int_u_inf);
        }
        if state.u_inf > a_sup || state.u_sup < an_inf {
            continue;
        }
        if state.u_inf >= an_inf {
            an_inf = state.u_inf;
        }
        if state.u_sup <= a_sup {
            a_sup = state.u_sup;
        }
        if (a_sup - an_inf) <= state.tol_u {
            continue;
        }
        state.int_u_inf = an_inf;
        state.int_u_sup = a_sup;

        if i != 0 {
            let (pp, v1) = curve.d1(state.int_u_inf);
            let s1 = GpVec::from_pnts(point, &pp).dot(&v1);
            if s1 * s2 < 0.0 {
                state.add_sol(state.int_u_inf, pp, pp.square_distance(point), s1 < 0.0);
            }
        }
        if i + 1 != n {
            let (pp, v1) = curve.d1(state.int_u_sup);
            s2 = GpVec::from_pnts(point, &pp).dot(&v1);
        }
        state.interval_perform(curve, point);
        int_ext_done = int_ext_done || state.done;
    }
    state.done = int_ext_done;
}

fn postprocess_ends(state: &mut ExtPcState) {
    // `GGExtPC.hxx:474-502` for BSpline / Offset / OtherCurve.
    if state.dist1 >= CONFUSION * CONFUSION && state.dist2 >= CONFUSION * CONFUSION {
        return;
    }
    let mut is_first = false;
    let mut is_last = false;
    for s in &state.sols {
        if (s.u - state.u_inf).abs() < state.tol_u {
            is_first = true;
        } else if (state.u_sup - s.u).abs() < state.tol_u {
            is_last = true;
        }
    }
    if !is_first && state.dist1 < CONFUSION * CONFUSION {
        state.sols.insert(
            0,
            ExtPcSolution {
                u: state.u_inf,
                point: state.pf,
                sq_dist: state.dist1,
                is_min: true,
            },
        );
    }
    if !is_last && state.dist2 < CONFUSION * CONFUSION {
        state.add_sol(state.u_sup, state.pl, state.dist2, true);
    }
    state.done = true;
}

/// `Extrema_ExtPC(P, C, Uinf, Usup)` — range-aware Perform.
///
/// Analytic Line/Circle stay in `point_curve_extrema_all`; this path covers
/// BSpline / Bezier / OtherCurve matching `GGExtPC::Perform`.
pub fn extrema_ext_pc_range(
    curve: &dyn Curve,
    point: &GpPnt,
    u_inf: f64,
    u_sup: f64,
) -> Vec<ExtPcSolution> {
    let (mut u_inf, mut u_sup) = (u_inf, u_sup);
    if u_inf > u_sup {
        std::mem::swap(&mut u_inf, &mut u_sup);
    }
    let mut state = ExtPcState::new(curve, u_inf, u_sup);
    if Precision::is_infinite(u_inf) {
        state.dist1 = f64::MAX;
    } else {
        state.pf = curve.d0(u_inf);
        state.dist1 = point.square_distance(&state.pf);
    }
    if Precision::is_infinite(u_sup) {
        state.dist2 = f64::MAX;
    } else {
        state.pl = curve.d0(u_sup);
        state.dist2 = point.square_distance(&state.pl);
    }

    // `Extrema_GGExtPC` curve-type switch: elementary curves go to
    // `Extrema_ExtPElC`, BSpline/Bezier/OtherCurve to the `default:` arm
    // (`Extrema_GGExtPC.hxx:390-502`).
    if let Some(pairs) = super::point_curve::ext_pelc_all(curve, point, u_inf, u_sup) {
        // `ExtPElC` filters to `[Uinf, Usup]` itself (`Extrema_ExtPElC.cxx:180`)
        // and sets `myIsMin` per solution (`cxx:77`, `:184`); `myDone = true`
        // regardless of the number of solutions (`cxx:189`).
        for (e, is_min) in pairs {
            state.add_sol(e.u1, e.p2, e.distance * e.distance, is_min);
        }
        state.done = true;
    } else if curve.bspline_knots().is_some() {
        perform_bspline(&mut state, curve, point);
        postprocess_ends(&mut state);
    } else {
        perform_general(&mut state, curve, point);
        postprocess_ends(&mut state);
    }
    state.sols
}

/// Closest minimum among `IsMin` solutions (`ShapeAnalysis_Curve.cxx:280-301`).
pub fn extrema_ext_pc_min_in_range(
    curve: &dyn Curve,
    point: &GpPnt,
    u_inf: f64,
    u_sup: f64,
) -> Option<(f64, GpPnt, f64)> {
    let sols = extrema_ext_pc_range(curve, point, u_inf, u_sup);
    let mut best: Option<ExtPcSolution> = None;
    for s in sols {
        if !s.is_min {
            continue;
        }
        if best.map_or(true, |b| s.sq_dist < b.sq_dist) {
            best = Some(s);
        }
    }
    best.map(|s| (s.u, s.point, s.sq_dist.sqrt()))
}
