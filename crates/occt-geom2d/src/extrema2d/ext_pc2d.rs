//! `Extrema_ExtPC2d` -- all extremal distances between a 2D point and a 2D curve.
//!
//! Source: `Extrema_ExtPC2d.hxx:31-39` (the alias
//! `Extrema_GGExtPC<Adaptor2d_Curve2d, Extrema_Curve2dTool, Extrema_ExtPElC2d,
//! gp_Pnt2d, gp_Vec2d, Extrema_POnCurv2d, NCollection_Sequence<...>,
//! Extrema_EPCOfExtPC2d>`) and `Extrema_GGExtPC.hxx`: `Initialize` `:114-127`,
//! `Perform` `:131-528`, `IntervalPerform` `:590-614`, `AddSol` `:617-632`,
//! postprocessing `:473-527`.
//!
//! The engine behind `IntervalPerform` is `Extrema_EPCOfExtPC2d` =
//! `Extrema_GGenExtPC` (`Extrema_GGenExtPC.hxx:36-232`), whose `Perform`
//! (`:160-173`) runs `math_FunctionRoots` over the sub-interval on the
//! `Extrema_PCFOfEPCOfExtPC2d` function ported in `super::gfunc_ext_pc2d`.
//!
//! Curve-type dispatch (`hxx:160-181`): Circle / Ellipse / Parabola / Hyperbola
//! / Line go to `Extrema_ExtPElC2d` (the analytic arms already ported in
//! `super::analytic_solvers`); Bezier runs one `IntervalPerform` with
//! `mysample = NbPoles()*2` and returns; BSpline runs the knot-span arms
//! (`hxx:190-388`); every other type runs the `default:` arm (`hxx:390-469`).
//!
//! `myIsMin` of `Extrema_ExtPElC2d` is `(NoSol == 0)` for the
//! Ellipse/Hyperbola/Parabola arms (`Extrema_ExtPElC2d.cxx:208`, `:277`,
//! `:343`), i.e. always `false`; this port keeps that value rather than
//! reporting a minimum the OCCT run would not report.

use occt_core::bspl::knots::unique_knots_mults;
use occt_core::bspl::locate::{first_u_knot_index, last_u_knot_index};
use occt_core::elib::clib2d::in_period;
use occt_core::gcpnts::{perform_tangential_curve, CurveSecondDeriv};
use occt_core::gp::{GpPnt, GpPnt2d, GpVec, GpVec2d};
use occt_core::kernel::geomabs::CurveType;
use occt_core::math_function_roots::FunctionRoots;
use occt_core::precision::Precision;

use crate::adaptor_curve::basis2d;

use super::analytic_solvers::{
    circle_all2d, ellipse_all2d, hyperbola_all2d, line_all2d, parabola_all2d,
};
use super::curve2d_tool::{
    adaptor_intervals2d, adaptor_is_periodic2d, adaptor_period2d, adaptor_resolution2d,
    adaptor_type2d,
};
use super::gfunc_ext_pc2d::FuncExtPc2d;
use super::prelude::*;

/// `GGExtPC::Perform` `mysample` reset (`Extrema_GGExtPC.hxx:139`).
const SAMPLE_DEFAULT: i32 = 17;
/// `default:` arm `aMaxSample` (`Extrema_GGExtPC.hxx:391`).
const MAX_SAMPLE: i32 = 17;
/// `t3d` handed to `Extrema_ExtPElC2d` (`Extrema_GGExtPC.hxx:140`).
const T3D: f64 = CONFUSION;
/// `myTolF` default (`Extrema_GGExtPC.hxx:87`, `:99`, `:113`).
const TOL_F: f64 = 1.0e-10;
/// `Extrema_Curve2dTool::DeflCurvIntervals` `epsd` (`Extrema_Curve2dTool.cxx:29`).
const EPS_DL_DT: f64 = 1.0e-3;
/// `maxdefl` (`Extrema_Curve2dTool.cxx:30`).
const MAX_DEFL: f64 = 1.0e3;
/// `mindefl` (`Extrema_Curve2dTool.cxx:31`).
const MIN_DEFL: f64 = 1.0e-3;
/// `nbpnts` (`Extrema_Curve2dTool.cxx:33`).
const DEFL_LENGTH_SAMPLES: i32 = 23;
/// `GeomAbs_Shape` numerics used by `Extrema_Curve2dTool::{NbIntervals,Intervals}`.
const GEOM_ABS_C2: u8 = 4;
const GEOM_ABS_CN: u8 = 6;

/// One `Extrema_ExtPC2d` solution: `mypoint(i)` / `mySqDist(i)` / `myismin(i)`.
#[derive(Debug, Clone, Copy)]
pub struct ExtPc2dSolution {
    /// `Point(i).Parameter()`.
    pub u: f64,
    /// `Point(i).Value()`.
    pub point: GpPnt2d,
    /// `SquareDistance(i)`.
    pub sq_dist: f64,
    /// `IsMin(i)`.
    pub is_min: bool,
}

/// `Extrema_GGExtPC::Initialize` (`Extrema_GGExtPC.hxx:114-127`) plus the
/// `Extrema_ExtPC2d(P, C, Uinf, Usup)` construction (`:88-99`).
struct Pc2dState<'a> {
    c: &'a dyn Curve2d,
    point: GpPnt2d,
    curve_type: CurveType,
    sols: Vec<ExtPc2dSolution>,
    done: bool,
    /// `mydist1` / `mydist2`.
    dist1: f64,
    dist2: f64,
    /// `Pf` / `Pl`.
    pf: GpPnt2d,
    pl: GpPnt2d,
    /// `mytolu` = `Extrema_Curve2dTool::Resolution(C, Precision::Confusion())`.
    tol_u: f64,
    /// `myuinf` / `myusup`.
    u_inf: f64,
    u_sup: f64,
    /// `mysample`.
    sample: i32,
    /// `myintuinf` / `myintusup`.
    int_u_inf: f64,
    int_u_sup: f64,
}

impl<'a> Pc2dState<'a> {
    fn new(c: &'a dyn Curve2d, p: &GpPnt2d, u_inf: f64, u_sup: f64) -> Self {
        Self {
            c,
            point: *p,
            curve_type: adaptor_type2d(c),
            sols: Vec::new(),
            done: false,
            dist1: f64::MAX,
            dist2: f64::MAX,
            pf: GpPnt2d::zero(),
            pl: GpPnt2d::zero(),
            // `hxx:123`.
            tol_u: adaptor_resolution2d(c, CONFUSION),
            u_inf,
            u_sup,
            sample: SAMPLE_DEFAULT,
            int_u_inf: u_inf,
            int_u_sup: u_sup,
        }
    }

    /// `GGExtPC::AddSol` (`Extrema_GGExtPC.hxx:617-632`).
    fn add_sol(&mut self, u: f64, p: GpPnt2d, sq_dist: f64, is_min: bool) {
        for s in &self.sols {
            if (s.u - u).abs() <= self.tol_u {
                return;
            }
        }
        self.sols.push(ExtPc2dSolution { u, point: p, sq_dist, is_min });
    }

    /// `GGExtPC::IntervalPerform` (`Extrema_GGExtPC.hxx:590-614`) on top of
    /// `Extrema_GGenExtPC::Perform` (`Extrema_GGenExtPC.hxx:160-173`).
    fn interval_perform(&mut self) {
        let c = self.c;
        let mut f = FuncExtPc2d::new(c, &self.point, self.int_u_inf, self.int_u_sup);
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
        // `hxx:598-609`.
        let period = if adaptor_is_periodic2d(c) { adaptor_period2d(c) } else { 0.0 };
        for s in f.sols() {
            let mut u = s.u;
            if period.abs() > 0.0 {
                u = in_period(u, self.u_inf, self.u_inf + period);
            }
            if u >= self.u_inf - self.tol_u && u <= self.u_sup + self.tol_u {
                self.add_sol(u, s.point, s.sq_dist, s.is_min);
            }
        }
    }
}

/// `Extrema_ExtPC2d` (`Extrema_ExtPC2d.hxx:31-39`) result set.
pub struct ExtremaExtPC2d<'a> {
    state: Pc2dState<'a>,
}

impl<'a> ExtremaExtPC2d<'a> {
    /// `Extrema_ExtPC2d(theP, theC, theTolF)` (`Extrema_GGExtPC.hxx:96-105`):
    /// the curve's own range.
    pub fn new(c: &'a dyn Curve2d, p: &GpPnt2d) -> Self {
        Self::new_range(c, p, c.first_parameter(), c.last_parameter())
    }

    /// `Extrema_ExtPC2d(theP, theC, theUinf, theUsup, theTolF)`
    /// (`Extrema_GGExtPC.hxx:84-93`).
    pub fn new_range(c: &'a dyn Curve2d, p: &GpPnt2d, u_inf: f64, u_sup: f64) -> Self {
        let mut out = Self { state: Pc2dState::new(c, p, u_inf, u_sup) };
        out.perform();
        out
    }

    /// `IsDone()` (`Extrema_GGExtPC.hxx:530-531`).
    pub fn is_done(&self) -> bool {
        self.state.done
    }

    /// `NbExt()` (`Extrema_GGExtPC.hxx:544-549`); OCCT throws when not done,
    /// the port returns 0.
    pub fn nb_ext(&self) -> usize {
        if !self.state.done {
            return 0;
        }
        self.state.sols.len()
    }

    /// `SquareDistance(theN)` (`Extrema_GGExtPC.hxx:535-542`), `theN` 1-based.
    pub fn square_distance(&self, n: usize) -> f64 {
        self.state.sols[n - 1].sq_dist
    }

    /// `IsMin(theN)` (`Extrema_GGExtPC.hxx:554-560`), `theN` 1-based.
    pub fn is_min(&self, n: usize) -> bool {
        self.state.sols[n - 1].is_min
    }

    /// `Point(theN)` (`Extrema_GGExtPC.hxx:563-570`): `(parameter, point)`.
    pub fn point(&self, n: usize) -> (f64, GpPnt2d) {
        let s = self.state.sols[n - 1];
        (s.u, s.point)
    }

    /// The `mypoint` / `mySqDist` / `myismin` sequences, in solver order.
    pub fn solutions(&self) -> &[ExtPc2dSolution] {
        &self.state.sols
    }

    /// `Extrema_GGExtPC::Perform` (`Extrema_GGExtPC.hxx:131-528`).
    fn perform(&mut self) {
        let st = &mut self.state;
        st.sols.clear();
        st.sample = SAMPLE_DEFAULT;
        let p = st.point;
        let c = st.c;

        // `hxx:142-157`.
        if Precision::is_infinite(st.u_inf) {
            st.dist1 = f64::MAX;
        } else {
            st.pf = c.d0(st.u_inf);
            st.dist1 = p.square_distance(&st.pf);
        }
        if Precision::is_infinite(st.u_sup) {
            st.dist2 = f64::MAX;
        } else {
            st.pl = c.d0(st.u_sup);
            st.dist2 = p.square_distance(&st.pl);
        }

        match st.curve_type {
            CurveType::Circle
            | CurveType::Ellipse
            | CurveType::Parabola
            | CurveType::Hyperbola
            | CurveType::Line => {
                // `hxx:162-181` + the elementary postprocessing (`hxx:504-527`).
                let (done, sols) = elementary2d(c, &p, st.curve_type, st.u_inf, st.u_sup);
                st.done = done;
                if done {
                    let period = if adaptor_is_periodic2d(c) { adaptor_period2d(c) } else { 0.0 };
                    for (u, pt, sq_dist, is_min) in sols {
                        // `PC.SetValues(U, myExtPElC.Point(i).Value())`: the
                        // parameter is remapped, the point is not re-evaluated.
                        let uu = if period.abs() > 0.0 {
                            in_period(u, st.u_inf, st.u_inf + period)
                        } else {
                            u
                        };
                        if uu >= st.u_inf - st.tol_u && uu <= st.u_sup + st.tol_u {
                            st.sols.push(ExtPc2dSolution { u: uu, point: pt, sq_dist, is_min });
                        }
                    }
                }
            }
            CurveType::BezierCurve => {
                // `hxx:182-188`: `IntervalPerform` then `return`.
                st.int_u_inf = st.u_inf;
                st.int_u_sup = st.u_sup;
                let np = basis2d(c).bezier_nb_poles().unwrap_or(0) as i32;
                st.sample = np * 2;
                st.interval_perform();
                return;
            }
            CurveType::BSplineCurve => {
                perform_bspline(st);
                // `hxx:387`.
                st.done = true;
            }
            _ => {
                // `hxx:390-469`.
                perform_general(st);
            }
        }

        postprocess_ends(st);
    }
}

/// `myExtPElC.Perform(theP, C, t3d, myuinf, myusup)` for each elementary type
/// (`Extrema_GGExtPC.hxx:162-181` -> `Extrema_ExtPElC2d::Perform`):
/// Line `Extrema_ExtPElC2d.cxx:56-79`, Circle `:94-152`, Ellipse `:168-212`,
/// Hyperbola `:228-284`, Parabola `:298-350`.
///
/// Returns `(IsDone(), (U, Point, SquareDistance, IsMin))`. `IsMin` is `true`
/// for the Line (`:78`) and Circle (`:140`, `:150`) arms and `(NoSol == 0)` —
/// i.e. always `false` — for Ellipse (`:208`), Hyperbola (`:277`) and Parabola
/// (`:343`).
#[allow(clippy::type_complexity)]
fn elementary2d(
    c: &dyn Curve2d,
    p: &GpPnt2d,
    curve_type: CurveType,
    u_inf: f64,
    u_sup: f64,
) -> (bool, Vec<(f64, GpPnt2d, f64, bool)>) {
    let collect = |v: Vec<super::analytic_solvers::Extrema2d>, is_min: bool| {
        v.into_iter()
            .map(|e| (e.u1, e.p2, e.p2.square_distance(p), is_min))
            .collect()
    };
    match curve_type {
        CurveType::Line => match c.gp_lin2d() {
            Some(l) => (true, collect(line_all2d(&l, p, u_inf, u_sup), true)),
            None => (true, Vec::new()),
        },
        CurveType::Circle => match c.gp_circ2d() {
            Some(gc) => {
                // `Extrema_ExtPElC2d.cxx:100-103`: not done when P is the centre.
                if gc.location().is_equal(p, T3D) {
                    return (false, Vec::new());
                }
                (true, collect(circle_all2d(&gc, p, u_inf, u_sup), true))
            }
            None => (true, Vec::new()),
        },
        CurveType::Ellipse => match c.gp_elips2d() {
            Some(e) => {
                // `:186-189`.
                if e.pos.point.is_equal(p, T3D)
                    && (e.major_radius - e.minor_radius).abs() <= T3D
                {
                    return (false, Vec::new());
                }
                (true, collect(ellipse_all2d(&e, p, u_inf, u_sup), false))
            }
            None => (true, Vec::new()),
        },
        CurveType::Hyperbola => match c.gp_hypr2d() {
            Some(h) => (true, collect(hyperbola_all2d(&h, p, u_inf, u_sup), false)),
            None => (true, Vec::new()),
        },
        CurveType::Parabola => match c.gp_parab2d() {
            Some(pa) => {
                // `:310-317`: `math_DirectPolynomialRoots(1/(4F), ...)` is not
                // done for `F == 0`, which `parabola_all2d` reports as empty.
                if pa.focal.abs() < RESOLUTION {
                    return (false, Vec::new());
                }
                (true, collect(parabola_all2d(&pa, p, u_inf, u_sup), false))
            }
            None => (true, Vec::new()),
        },
        _ => (true, Vec::new()),
    }
}

/// `GeomAbs_BSplineCurve` arm of `GGExtPC::Perform` (`Extrema_GGExtPC.hxx:190-388`).
fn perform_bspline(st: &mut Pc2dState<'_>) {
    let c = st.c;
    let basis = basis2d(c);
    let flat = match basis.bspline_knots2d() {
        Some(k) if !k.is_empty() => k,
        // `TheCurveTool::BSpline(aCurve)` always answers for an adaptor whose
        // `GetType()` is `GeomAbs_BSplineCurve`; fall back like the 3D port.
        _ => {
            perform_general(st);
            return;
        }
    };
    let degree = basis.bspline_degree().unwrap_or(1) as i32;
    let (knots, mults) = unique_knots_mults(flat);
    if knots.len() < 2 {
        perform_general(st);
        return;
    }
    let first_idx = first_u_knot_index(degree, &mults);
    let last_idx = last_u_knot_index(degree, &mults);
    let tol_coeff = (st.u_sup - st.u_inf) * PCONFUSION;

    // `hxx:196-204`.
    let mut period_jump = 0.0;
    if adaptor_is_periodic2d(c) {
        let period = adaptor_period2d(c);
        if period.abs() > 0.0 {
            let mut shift = ((st.u_inf - knots[(first_idx - 1) as usize]) / period) as i32;
            if st.u_inf < knots[(first_idx - 1) as usize] - tol_coeff {
                shift -= 1;
            }
            period_jump = period * shift as f64;
        }
    }

    // `hxx:206-227`.
    let mut first_used = first_idx;
    let mut last_used = last_idx;
    for an_idx in first_idx..=last_idx {
        let knot = knots[(an_idx - 1) as usize] + period_jump;
        if st.u_inf >= knot - tol_coeff {
            first_used = an_idx;
        } else {
            break;
        }
    }
    for an_idx in (first_idx..=last_idx).rev() {
        let knot = knots[(an_idx - 1) as usize] + period_jump;
        if st.u_sup <= knot + tol_coeff {
            last_used = an_idx;
        } else {
            break;
        }
    }
    if first_used == last_used {
        first_used = first_idx;
        last_used = first_idx + 1;
    }

    // `hxx:231`.
    st.sample = degree + 1;

    if st.sample == 2 {
        // Degree-1 knot-span arm (`hxx:233-297`).
        let mut tmin = 0.0;
        let mut distmin = f64::MAX;
        let mut pmin = GpPnt2d::zero();
        let mut min2 = 0.0;
        for an_idx in first_used..last_used {
            let mut f = knots[(an_idx - 1) as usize] + period_jump;
            let mut l = knots[an_idx as usize] + period_jump;
            if an_idx == first_used {
                f = st.u_inf;
            } else if an_idx == last_used - 1 {
                l = st.u_sup;
            }
            let p1_val = c.d0(f);
            let p2_val = c.d0(l);
            let p = st.point;
            let base1 = vec_from_pnts2d(&p, &p1_val);
            let base2 = vec_from_pnts2d(&p, &p2_val);
            let v = vec_from_pnts2d(&p2_val, &p1_val);
            let val1 = v.dot(&base1);
            let val2 = v.dot(&base2);
            let min1 = if an_idx == first_used {
                p.square_distance(&p1_val)
            } else {
                let prev = min2;
                if distmin > prev {
                    distmin = prev;
                    tmin = f;
                    pmin = p1_val;
                }
                prev
            };
            min2 = p.square_distance(&p2_val);
            let min_sq = min1.min(min2);
            let min_der = val1.abs().min(val2.abs());
            if !(Precision::is_infinite(val1) || Precision::is_infinite(val2))
                && (val1 * val2 <= 0.0
                    || min_sq < 100.0 * SQUARE_CONFUSION
                    || 2.0 * min_der < CONFUSION)
            {
                st.int_u_inf = f;
                st.int_u_sup = l;
                st.interval_perform();
            }
        }
        if !Precision::is_infinite(distmin) {
            let mut is_to_add = true;
            for s in &st.sols {
                if !(distmin < s.sq_dist && (s.u - tmin).abs() > st.tol_u) {
                    is_to_add = false;
                    break;
                }
            }
            if is_to_add {
                st.add_sol(tmin, pmin, distmin, true);
            }
        }
    } else {
        // Degree >= 2 sample + IntervalPerform (`hxx:300-385`).
        let n_spans = (last_used - first_used).max(1) as usize;
        let cap = (st.sample as usize) * n_spans + 1;
        let mut vals = Vec::with_capacity(cap);
        let mut params = Vec::with_capacity(cap);
        for an_idx in first_used..last_used {
            let mut f = knots[(an_idx - 1) as usize] + period_jump;
            let mut l = knots[an_idx as usize] + period_jump;
            if an_idx == first_used {
                f = st.u_inf;
            }
            if an_idx == last_used - 1 {
                l = st.u_sup;
            }
            let step = (l - f) / st.sample as f64;
            for pnt_idx in 0..st.sample {
                let cur = f + step * pnt_idx as f64;
                vals.push(c.d0(cur).square_distance(&st.point));
                params.push(cur);
            }
        }
        vals.push(c.d0(st.u_sup).square_distance(&st.point));
        params.push(st.u_sup);

        for an_idx in 1..vals.len().saturating_sub(1) {
            if vals[an_idx] <= SQUARE_CONFUSION {
                st.sols.push(ExtPc2dSolution {
                    u: params[an_idx],
                    point: c.d0(params[an_idx]),
                    sq_dist: vals[an_idx],
                    is_min: true,
                });
            }
            if (vals[an_idx] >= vals[an_idx + 1] && vals[an_idx] >= vals[an_idx - 1])
                || (vals[an_idx] <= vals[an_idx + 1] && vals[an_idx] <= vals[an_idx - 1])
            {
                st.int_u_inf = params[an_idx - 1];
                st.int_u_sup = params[an_idx + 1];
                st.interval_perform();
            }
        }

        if st.dist1 > SQUARE_CONFUSION
            && !Precision::is_positive_infinite(st.dist1)
            && params.len() >= 2
        {
            let (p1, v1) = c.d1(params[0]);
            let (p2, v2) = c.d1(params[1]);
            let base1 = vec_from_pnts2d(&st.point, &p1);
            let base2 = vec_from_pnts2d(&st.point, &p2);
            let val1 = v1.dot(&base1);
            let val2 = v2.dot(&base2);
            if !(Precision::is_infinite(val1) || Precision::is_infinite(val2))
                && (val1 * val2 <= 0.0
                    || base1.dot(&base2) <= 0.0
                    || 2.0 * val1.abs() < CONFUSION)
            {
                st.int_u_inf = params[0];
                st.int_u_sup = params[1];
                st.interval_perform();
            }
        }
        if st.dist2 > SQUARE_CONFUSION
            && !Precision::is_positive_infinite(st.dist2)
            && params.len() >= 2
        {
            let n = params.len();
            let (p1, v1) = c.d1(params[n - 2]);
            let (p2, v2) = c.d1(params[n - 1]);
            let base1 = vec_from_pnts2d(&st.point, &p1);
            let base2 = vec_from_pnts2d(&st.point, &p2);
            let val1 = v1.dot(&base1);
            let val2 = v2.dot(&base2);
            if !(Precision::is_infinite(val1) || Precision::is_infinite(val2))
                && (val1 * val2 <= 0.0
                    || base1.dot(&base2) <= 0.0
                    || 2.0 * val2.abs() < CONFUSION)
            {
                st.int_u_inf = params[n - 2];
                st.int_u_sup = params[n - 1];
                st.interval_perform();
            }
        }
    }
}

/// Adapter over a 2D curve for `GCPnts_TangentialDeflection`, mirroring the
/// lift `gp_Pnt2d`/`gp_Vec2d` -> `gp_Pnt`/`gp_Vec` with `Z = 0`
/// (`GCPnts_TangentialDeflection.cxx:46-67`).
struct CurveTd2d<'a>(&'a dyn Curve2d);

impl CurveSecondDeriv for CurveTd2d<'_> {
    fn point(&self, u: f64) -> GpPnt {
        let p = self.0.d0(u);
        GpPnt::new(p.x(), p.y(), 0.0)
    }
    fn d2(&self, u: f64) -> (GpPnt, GpVec, GpVec) {
        let (p, d1, d2) = self.0.d2(u);
        (
            GpPnt::new(p.x(), p.y(), 0.0),
            GpVec::new(d1.x(), d1.y(), 0.0),
            GpVec::new(d2.x(), d2.y(), 0.0),
        )
    }
}

/// `Extrema_Curve2dTool::DeflCurvIntervals` (`Extrema_Curve2dTool.cxx:26-73`).
///
/// The trailing `GCPnts_TangentialDeflection` call goes through the shared
/// `PerformCurve` port. `GCPnts_TangentialDeflection::initialize`
/// (`GCPnts_TangentialDeflection.cxx:415-452`) routes `GeomAbs_Line` /
/// `GeomAbs_Circle` to `PerformLinear` / `PerformCircular` and a 2-pole
/// BSpline/Bezier to `PerformLinear`; those arms are unreachable here because
/// `GGExtPC::Perform` only calls this from its `default:` arm, i.e. for
/// Parabola / Hyperbola / OffsetCurve / OtherCurve.
fn defl_curv_intervals2d(c: &dyn Curve2d) -> Vec<f64> {
    let tf = c.first_parameter();
    let tl = c.last_parameter();
    if !(tf.is_finite() && tl.is_finite()) || (tl - tf).abs() <= 0.0 {
        return vec![tf, tl];
    }

    // Chord sum from `Value(tf)`, OCCT does not advance the running point
    // (`Extrema_Curve2dTool.cxx:37-42`).
    let p = c.d0(tf);
    let mut length = 0.0;
    let nbpnts = DEFL_LENGTH_SAMPLES;
    for i in 2..=nbpnts {
        let t = (tf * (nbpnts - i) as f64 + (i - 1) as f64 * tl) / (nbpnts - 1) as f64;
        length += p.distance(&c.d0(t));
    }

    let span = tl - tf;
    let dldt = length / span;
    if length <= CONFUSION || dldt < EPS_DL_DT || span > 10000.0 {
        return vec![tf, tl];
    }

    let defl = (0.01 * length / (2.0 * std::f64::consts::PI)).max(MIN_DEFL);
    if defl > MAX_DEFL {
        return vec![tf, tl];
    }
    let min_len = (0.00001 * length).max(CONFUSION);
    let tol = (0.00001 * span).max(PCONFUSION);

    let mut cn = adaptor_intervals2d(c, GEOM_ABS_CN);
    if cn.len() < 2 {
        cn = vec![tf, tl];
    }
    let degree_min_nb = basis2d(c)
        .bspline_degree()
        .map(|d| (d + 1).max(2))
        .unwrap_or(2);
    let (params, _) = perform_tangential_curve(
        &CurveTd2d(c),
        tf,
        tl,
        std::f64::consts::PI / 6.0,
        defl,
        2,
        tol,
        min_len,
        &cn,
        degree_min_nb,
    );
    if params.len() < 2 {
        vec![tf, tl]
    } else {
        params
    }
}

/// `default:` arm of `GGExtPC::Perform` (`Extrema_GGExtPC.hxx:390-469`).
fn perform_general(st: &mut Pc2dState<'_>) {
    let c = st.c;

    // `hxx:395-405`: C2 intervals when `NbIntervals(C2) > 1`, else
    // `DeflCurvIntervals`.
    let c2 = adaptor_intervals2d(c, GEOM_ABS_C2);
    let n_c2 = c2.len().saturating_sub(1);
    let intervals = if n_c2 > 1 { c2 } else { defl_curv_intervals2d(c) };
    let n = intervals.len().saturating_sub(1).max(1);
    let mut maxint = 0.0;
    for i in 0..n {
        let dt = intervals[i + 1] - intervals[i];
        if maxint < dt {
            maxint = dt;
        }
    }
    if maxint <= 0.0 {
        maxint = (st.u_sup - st.u_inf).abs().max(CONFUSION);
    }
    let is_periodic = adaptor_is_periodic2d(c);
    let mut int_ext_done = false;
    let mut s2 = 0.0;
    for i in 0..n {
        // `hxx:422-424`: the raw interval bounds are what `IntervalPerform`
        // searches; the clamped copies below are only the validity test.
        let int_inf = intervals[i];
        let int_sup = intervals[i + 1];
        st.int_u_inf = int_inf;
        st.int_u_sup = int_sup;
        // `RealToInt` truncates toward zero (`hxx:424`).
        st.sample = ((MAX_SAMPLE as f64) * (int_sup - int_inf) / maxint) as i32;
        st.sample = st.sample.max(3);

        let mut an_inf = int_inf;
        let mut sup = int_sup;
        if is_periodic {
            let period = adaptor_period2d(c);
            an_inf = in_period(int_inf, st.u_inf, st.u_inf + period);
            sup = int_sup + (an_inf - int_inf);
        }
        if st.u_inf > sup || st.u_sup < an_inf {
            continue;
        }
        if st.u_inf >= an_inf {
            an_inf = st.u_inf;
        }
        if st.u_sup <= sup {
            sup = st.u_sup;
        }
        if (sup - an_inf) <= st.tol_u {
            continue;
        }

        if i != 0 {
            let (pp, v1) = c.d1(int_inf);
            let s1 = vec_from_pnts2d(&st.point, &pp).dot(&v1);
            if s1 * s2 < 0.0 {
                st.sols.push(ExtPc2dSolution {
                    u: int_inf,
                    point: pp,
                    sq_dist: pp.square_distance(&st.point),
                    is_min: s1 < 0.0,
                });
            }
        }
        if i + 1 != n {
            let (pp, v1) = c.d1(int_sup);
            s2 = vec_from_pnts2d(&st.point, &pp).dot(&v1);
        }
        st.interval_perform();
        int_ext_done = int_ext_done || st.done;
    }
    st.done = int_ext_done;
}

/// Postprocessing (`Extrema_GGExtPC.hxx:473-527`).
fn postprocess_ends(st: &mut Pc2dState<'_>) {
    if matches!(
        st.curve_type,
        CurveType::BSplineCurve | CurveType::OffsetCurve | CurveType::OtherCurve
    ) {
        // `hxx:474-502`.
        if st.dist1 >= SQUARE_CONFUSION && st.dist2 >= SQUARE_CONFUSION {
            return;
        }
        let mut is_first = false;
        let mut is_last = false;
        for s in &st.sols {
            if (s.u - st.u_inf).abs() < st.tol_u {
                is_first = true;
            } else if (st.u_sup - s.u).abs() < st.tol_u {
                is_last = true;
            }
        }
        if !is_first && st.dist1 < SQUARE_CONFUSION {
            st.sols.insert(
                0,
                ExtPc2dSolution {
                    u: st.u_inf,
                    point: st.pf,
                    sq_dist: st.dist1,
                    is_min: true,
                },
            );
        }
        if !is_last && st.dist2 < SQUARE_CONFUSION {
            st.sols.push(ExtPc2dSolution {
                u: st.u_sup,
                point: st.pl,
                sq_dist: st.dist2,
                is_min: true,
            });
        }
        st.done = true;
    }
}

/// `TheVector(P1, P2)` for `gp_Vec2d`.
fn vec_from_pnts2d(p1: &GpPnt2d, p2: &GpPnt2d) -> GpVec2d {
    GpVec2d::new(p2.x() - p1.x(), p2.y() - p1.y())
}
