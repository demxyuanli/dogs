//! `Extrema_GFuncExtPC` 2D instantiation (`Extrema_PCFOfEPCOfExtPC2d`).
//!
//! Source: `Extrema_PCFOfEPCOfExtPC2d.hxx:29-36` (the template alias) and
//! `Extrema_GFuncExtPC.hxx`: `Initialize` `:104-129`, `SetPoint` `:132-141`,
//! `Value` `:146-253`, `Values` `:274-352`, `GetStateNumber` `:355-379`,
//! `SubIntervalInitialize` `:419-424`, `SearchOfTolerance` `:428-456`,
//! constants `:459-463`.
//!
//! Driven by `Extrema_GGExtPC::IntervalPerform` (`Extrema_GGExtPC.hxx:590-614`)
//! through `Extrema_EPCOfExtPC2d` = `Extrema_GGenExtPC`
//! (`Extrema_GGenExtPC.hxx:36-232`): `Perform` (`:160-173`) runs
//! `math_FunctionRoots` over the sub-interval and the roots are collected by
//! `GetStateNumber`.
//!
//! The `myIsMin`/`mySqDist`/`mypoint` sequences of `Extrema_GFuncExtPC` are the
//! collected root list here (`sols`), which is what `Extrema_GGenExtPC::NbExt`
//! / `Point` / `SquareDistance` / `IsMin` read back.

use occt_core::gp::{GpPnt2d, GpVec2d};
use occt_core::kernel::geomabs::CurveType;
use occt_core::math_fn::{MathFunction, MathFunctionWithDerivative};
use occt_core::precision::Precision;

use super::curve2d_tool::adaptor_type2d;
use super::prelude::*;

/// `GFuncExtPC::TolFactor` (`Extrema_GFuncExtPC.hxx:460`).
const TOL_FACTOR: f64 = 1.0e-12;
/// `GFuncExtPC::MinTol` (`Extrema_GFuncExtPC.hxx:461`).
const MIN_TOL: f64 = 1.0e-20;
/// `GFuncExtPC::MinStep` (`Extrema_GFuncExtPC.hxx:462`).
const MIN_STEP: f64 = 1.0e-7;
/// `GFuncExtPC::MaxOrder` (`Extrema_GFuncExtPC.hxx:463`).
const MAX_ORDER: i32 = 3;
/// `SearchOfTolerance` sample count (`Extrema_GFuncExtPC.hxx:430`).
const N_POINT: i32 = 10;
/// `Value` singular-D1 division factor (`Extrema_GFuncExtPC.hxx:169`).
const DIVISION_FACTOR_D1: f64 = 1.0e-3;
/// `Values` singular-DF division factor (`Extrema_GFuncExtPC.hxx:300`).
const DIVISION_FACTOR_DF: f64 = 0.01;
/// Taylor three-point coefficients (`Extrema_GFuncExtPC.hxx:236-238`, `:326-336`).
const TAYLOR_A: f64 = 3.0;
const TAYLOR_B: f64 = 4.0;

/// One extremum accepted by `GetStateNumber` (`Extrema_GFuncExtPC.hxx:355-379`).
#[derive(Debug, Clone, Copy)]
pub(super) struct GFuncSolution2d {
    /// `myU`, the parameter of the root.
    pub u: f64,
    /// `myPc`, the point on the curve.
    pub point: GpPnt2d,
    /// `mySqDist`, the square distance to `myP`.
    pub sq_dist: f64,
    /// `myIsMin == 1` i.e. `myD1f > 0` (`:370-373`).
    pub is_min: bool,
}

/// `Extrema_GFuncExtPC<Adaptor2d_Curve2d, Extrema_Curve2dTool, ...>`.
pub(super) struct FuncExtPc2d<'a> {
    curve: &'a dyn Curve2d,
    point: GpPnt2d,
    u: f64,
    pc: GpPnt2d,
    d1f: f64,
    u_inf: f64,
    u_sup: f64,
    /// `myMaxDerivOrder` (`hxx:68`); the `Initialize` GetType switch (`:115-128`).
    max_deriv_order: i32,
    /// `myTol` — `SearchOfTolerance()` or `MinTol` (`:92-98`).
    tol: f64,
    sols: Vec<GFuncSolution2d>,
}

impl<'a> FuncExtPc2d<'a> {
    /// `Extrema_GFuncExtPC(theP, theC)` (`hxx:76-101`) plus
    /// `Extrema_GGenExtPC::Initialize` (`Extrema_GGenExtPC.hxx:105-113`,
    /// `:128-132`): the tolerance is computed on the curve's own range and only
    /// then is the sub-interval narrowed to `[the_u_inf, the_u_sup]`.
    pub(super) fn new(
        curve: &'a dyn Curve2d,
        point: &GpPnt2d,
        u_inf_par: f64,
        u_sup_par: f64,
    ) -> Self {
        let mut out = Self {
            curve,
            point: *point,
            u: 0.0,
            pc: GpPnt2d::zero(),
            d1f: 0.0,
            u_inf: 0.0,
            u_sup: 0.0,
            max_deriv_order: 0,
            tol: MIN_TOL,
            sols: Vec::new(),
        };
        // `hxx:87-99` / `:115-128`.
        out.max_deriv_order = match adaptor_type2d(curve) {
            CurveType::BezierCurve
            | CurveType::BSplineCurve
            | CurveType::OffsetCurve
            | CurveType::OtherCurve => MAX_ORDER,
            _ => 0,
        };
        out.sub_interval_initialize(curve.first_parameter(), curve.last_parameter());
        out.tol = if out.max_deriv_order != 0 {
            out.search_of_tolerance()
        } else {
            MIN_TOL
        };
        out.sub_interval_initialize(u_inf_par, u_sup_par);
        out
    }

    /// `SubIntervalInitialize` (`hxx:419-424`).
    pub(super) fn sub_interval_initialize(&mut self, u_first: f64, u_last: f64) {
        self.u_inf = u_first;
        self.u_sup = u_last;
    }

    /// `SetPoint` (`hxx:132-141`): also clears the three result sequences.
    pub(super) fn set_point(&mut self, p: &GpPnt2d) {
        self.point = *p;
        self.sols.clear();
    }

    /// `Point(i).Parameter()` / `SquareDistance(i)` / `IsMin(i)` read-back of
    /// the collected `myPoint` / `mySqDist` / `myIsMin` (`:381-417`), 1-based.
    pub(super) fn sols(&self) -> &[GFuncSolution2d] {
        &self.sols
    }

    /// Sub-interval span for the singular delta (`hxx:170-176`, `:301-307`).
    fn singular_du(&self) -> f64 {
        if !self.u_sup.is_finite() || !self.u_inf.is_finite() {
            0.0
        } else {
            self.u_sup - self.u_inf
        }
    }

    /// `GFuncExtPC::Value` (`hxx:146-253`), including the `MaxOrder` singular D1
    /// reconstruction (`hxx:167-241`).
    fn value_at(&mut self, u_par: f64, f: &mut f64) -> bool {
        self.u = u_par;
        let (pc, mut d1c) = self.curve.d1(u_par);
        self.pc = pc;
        if Precision::is_infinite(d1c.x()) || Precision::is_infinite(d1c.y()) {
            *f = Precision::INFINITE;
            return false;
        }
        let mut ndu = d1c.magnitude();

        if self.max_deriv_order != 0 && ndu <= self.tol {
            // Singular case (`hxx:167-241`).
            let delta = (self.singular_du() * DIVISION_FACTOR_D1).max(MIN_STEP);

            // `do { V = DN(++n); ... } while (!found && n < MaxOrder)`.
            let mut n = 1;
            let mut v = GpVec2d::zero();
            let mut is_derive_found = false;
            while !is_derive_found && n < self.max_deriv_order {
                n += 1;
                v = self.curve.eval_dn(u_par, n);
                ndu = v.magnitude();
                is_derive_found = ndu > self.tol;
            }

            if is_derive_found {
                let u = if self.u - self.u_inf < delta {
                    self.u + delta
                } else {
                    self.u - delta
                };
                let p1 = self.curve.d0(self.u.min(u));
                let p2 = self.curve.d0(self.u.max(u));
                let v1 = vec_from_pnts2d(&p1, &p2);
                d1c = if v.dot(&v1) < 0.0 { v.reversed() } else { v };
            } else {
                // Three-point Taylor (`hxx:213-238`); `Ptemp` stays at the origin.
                let (p1, p2, p3, grown) = if self.u - self.u_inf < 2.0 * delta {
                    (
                        self.curve.d0(self.u),
                        self.curve.d0(self.u + delta),
                        self.curve.d0(self.u + 2.0 * delta),
                        true,
                    )
                } else {
                    (
                        self.curve.d0(self.u - 2.0 * delta),
                        self.curve.d0(self.u - delta),
                        self.curve.d0(self.u),
                        false,
                    )
                };
                let v1 = vec_from_pnts2d(&GpPnt2d::zero(), &p1);
                let v2 = vec_from_pnts2d(&GpPnt2d::zero(), &p2);
                let v3 = vec_from_pnts2d(&GpPnt2d::zero(), &p3);
                d1c = if grown {
                    v1.multiplied_scalar(-TAYLOR_A)
                        .added(&v2.multiplied_scalar(TAYLOR_B))
                        .subtracted(&v3)
                } else {
                    v1.subtracted(&v2.multiplied_scalar(TAYLOR_B))
                        .added(&v3.multiplied_scalar(TAYLOR_A))
                };
            }
            ndu = d1c.magnitude();
        }

        if ndu <= MIN_TOL {
            // Warning: 1st derivative is equal to zero!
            return false;
        }

        let ppc = vec_from_pnts2d(&self.point, &self.pc);
        *f = ppc.dot(&d1c) / ndu;
        true
    }

    /// `GFuncExtPC::Values` (`hxx:274-352`), including the singular DF
    /// three-point form (`hxx:297-338`).
    fn values_at(&mut self, u_par: f64, f: &mut f64, df: &mut f64) -> bool {
        let pc_old = self.pc;
        let p_old = self.point;
        if !self.value_at(u_par, f) {
            return false;
        }
        self.u = u_par;
        self.pc = pc_old;
        self.point = p_old;

        let (pc, d1c, d2c) = self.curve.d2(u_par);
        self.pc = pc;
        let ndu = d1c.magnitude();
        if ndu <= self.tol {
            // Singular DF (`hxx:297-338`).
            let delta = (self.singular_du() * DIVISION_FACTOR_DF).max(MIN_STEP);
            let (f1, f2, f3);
            if self.u - self.u_inf < 2.0 * delta {
                f1 = *f;
                let u2 = self.u + delta;
                let u3 = self.u + delta * 2.0;
                let mut f2v = 0.0;
                let mut f3v = 0.0;
                if !(self.value_at(u2, &mut f2v) && self.value_at(u3, &mut f3v)) {
                    return false;
                }
                f2 = f2v;
                f3 = f3v;
                *df = (-TAYLOR_A * f1 + TAYLOR_B * f2 - f3) / (2.0 * delta);
            } else {
                f3 = *f;
                let u1 = self.u - delta * 2.0;
                let u2 = self.u - delta;
                let mut f1v = 0.0;
                let mut f2v = 0.0;
                if !(self.value_at(u2, &mut f2v) && self.value_at(u1, &mut f1v)) {
                    return false;
                }
                f1 = f1v;
                f2 = f2v;
                *df = (f1 - TAYLOR_B * f2 + TAYLOR_A * f3) / (2.0 * delta);
            }
            self.u = u_par;
            self.pc = pc_old;
            self.point = p_old;
        } else {
            let ppc = vec_from_pnts2d(&self.point, &self.pc);
            *df = ndu + (ppc.dot(&d2c) / ndu) - *f * (d1c.dot(&d2c)) / (ndu * ndu);
        }

        self.d1f = *df;
        true
    }

    /// `GFuncExtPC::SearchOfTolerance` (`hxx:428-456`) over the sub-interval
    /// currently set (i.e. the curve's own `[First, Last]` at construction).
    fn search_of_tolerance(&mut self) -> f64 {
        let step = (self.u_sup - self.u_inf) / (N_POINT as f64);
        let mut num = 0;
        let mut max_val = -Precision::INFINITE;
        loop {
            let mut u = self.u_inf + (num as f64) * step;
            if u > self.u_sup {
                u = self.u_sup;
            }
            let (_p, v_der) = self.curve.d1(u);
            // `continue` inside the do-while still runs `++aNum` (`hxx:447`).
            if !(Precision::is_infinite(v_der.x()) || Precision::is_infinite(v_der.y())) {
                let vm = v_der.magnitude();
                if vm > max_val {
                    max_val = vm;
                }
            }
            num += 1;
            if num >= N_POINT + 1 {
                break;
            }
        }
        (max_val * TOL_FACTOR).max(MIN_TOL)
    }
}

/// `TheVector(P1, P2)` for `gp_Vec2d`.
fn vec_from_pnts2d(p1: &GpPnt2d, p2: &GpPnt2d) -> GpVec2d {
    GpVec2d::new(p2.x() - p1.x(), p2.y() - p1.y())
}

impl MathFunction for FuncExtPc2d<'_> {
    fn value(&mut self, x: f64, f: &mut f64) -> bool {
        self.value_at(x, f)
    }

    /// `GFuncExtPC::GetStateNumber` (`hxx:355-379`).
    fn get_state_number(&mut self) -> i32 {
        let sq = self.pc.square_distance(&self.point);
        // It is necessary to always compute myD1f (`hxx:365-368`).
        let mut ff = 0.0;
        let mut dd = 0.0;
        let _ = self.values_at(self.u, &mut ff, &mut dd);
        let is_min = self.d1f > 0.0;
        self.sols.push(GFuncSolution2d {
            u: self.u,
            point: self.pc,
            sq_dist: sq,
            is_min,
        });
        0
    }
}

impl MathFunctionWithDerivative for FuncExtPc2d<'_> {
    fn derivative(&mut self, x: f64, d: &mut f64) -> bool {
        let mut f = 0.0;
        self.values_at(x, &mut f, d)
    }

    fn values(&mut self, x: f64, f: &mut f64, d: &mut f64) -> bool {
        self.values_at(x, f, d)
    }
}
