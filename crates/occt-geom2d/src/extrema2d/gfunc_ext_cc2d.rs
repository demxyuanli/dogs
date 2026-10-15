//! `Extrema_GFuncExtCC` 2D instantiation (`Extrema_CCLocFOfLocECC2d`).
//!
//! Source: `Extrema_CCLocFOfLocECC2d.hxx:29-36` (the template alias) and
//! `Extrema_GFuncExtCC.hxx`: `SearchOfTolerance` `:150-197`,
//! `Extrema_GFuncExtCC(theTol)` `:204-230`,
//! `Extrema_GFuncExtCC(theC1, theC2, theTol)` `:235-291`, `SetCurve` `:299-351`,
//! `Value` `:359-550`, `Derivatives` `:551-565`, `Values` `:568-827`,
//! `GetStateNumber` `:832-870`, `Points` `:885-895`,
//! `SubIntervalInitialize` `:899-925`, constants `:109-112`.
//!
//! This is the `math_FunctionSetWithDerivatives` whose zero is searched by
//! `Extrema_GenLocateExtCC` (`Extrema_LocECC2d`), i.e. the engine behind
//! `Extrema_LocateExtCC2d` (`Extrema_LocateExtCC2d.cxx:39`).
//!
//! Two notes on the instantiation, neither a control-flow difference:
//!
//! * `TheCurveTool1` and `TheCurveTool2` are both `Extrema_Curve2dTool` here, so
//!   the methods OCCT spells `TheCurveTool1::...` on the *second* curve
//!   (`SearchOfTolerance` `:189`, which the template always routes through
//!   `TheCurveTool1`) are evaluated on that curve's tool.
//! * The `myC1`/`myC2` `void*` identity tests of `SearchOfTolerance`
//!   (`:162`, `:167`) select the first or the second sub-interval; the port
//!   selects it by rank (1 or 2), the same decision the pointer comparison makes.

use super::prelude::*;

use occt_core::kernel::geomabs::CurveType;
use occt_math::{MathFunctionSetWithDerivatives, MathMatrix, MathVector};

use super::curve2d_tool::adaptor_type2d;
use crate::geom2d_int::POnCurv2d;

/// `GFuncExtCC::THE_MIN_TOL` (`Extrema_GFuncExtCC.hxx:109`).
const THE_MIN_TOL: f64 = 1.0e-20;
/// `GFuncExtCC::THE_TOL_FACTOR` (`Extrema_GFuncExtCC.hxx:110`).
const THE_TOL_FACTOR: f64 = 1.0e-12;
/// `GFuncExtCC::THE_MIN_STEP` (`Extrema_GFuncExtCC.hxx:111`).
const THE_MIN_STEP: f64 = 1.0e-7;
/// `GFuncExtCC::THE_MAX_ORDER` (`Extrema_GFuncExtCC.hxx:112`).
const THE_MAX_ORDER: i32 = 3;
/// `RealLast()` (`Standard_Real.hxx:179-183`).
const REAL_LAST: f64 = f64::MAX;
/// `RealFirst()` (`Standard_Real.hxx:167-171`).
const REAL_FIRST: f64 = -f64::MAX;

/// `SearchOfTolerance`'s `NPoint` (`Extrema_GFuncExtCC.hxx:158`).
const N_POINT: i32 = 10;

/// `Adaptor2d_Curve2d` reports an unbounded curve through
/// `Precision::Infinite()` (= `2e100`, `Precision.hxx:371`) -- e.g.
/// `Geom2d_Line::FirstParameter()` (`Geom2d_Line.cxx:142-145`) -- while the
/// port's unbounded curves report IEEE `+/-infinity`, which would make
/// `myUsupremum - myUinfium` infinite and the `RealLast()`/`RealFirst()` tests
/// below (`Extrema_GFuncExtCC.hxx:388`, `:472`) meaningless. Substitute the
/// OCCT value; the same substitution is documented in `super::general_extrema`.
fn border(v: f64) -> f64 {
    if v.is_finite() {
        v
    } else if v < 0.0 {
        -INFINITE
    } else {
        INFINITE
    }
}

/// `Extrema_GFuncExtCC<Adaptor2d_Curve2d, Extrema_Curve2dTool,
/// Adaptor2d_Curve2d, Extrema_Curve2dTool, Extrema_POnCurv2d, gp_Pnt2d,
/// gp_Vec2d, NCollection_Sequence<Extrema_POnCurv2d>>`
/// (`Extrema_CCLocFOfLocECC2d.hxx:29-36`).
pub(super) struct GFuncExtCC2d<'a> {
    /// `myC1` (`hxx:114`).
    my_c1: &'a dyn Curve2d,
    /// `myC2` (`hxx:115`).
    my_c2: &'a dyn Curve2d,
    /// `myTol` (`hxx:116`).
    my_tol: f64,
    /// `myU` (`hxx:117`).
    my_u: f64,
    /// `myV` (`hxx:118`).
    my_v: f64,
    /// `myP1` (`hxx:119`).
    my_p1: GpPnt2d,
    /// `myP2` (`hxx:120`).
    my_p2: GpPnt2d,
    /// `myDu` (`hxx:121`).
    my_du: GpVec2d,
    /// `myDv` (`hxx:122`).
    my_dv: GpVec2d,
    /// `mySqDist` (`hxx:123`).
    my_sq_dist: Vec<f64>,
    /// `myPoints` (`hxx:124`); two consecutive entries per extremum.
    my_points: Vec<POnCurv2d>,
    /// `myTolC1` (`hxx:125`).
    my_tol_c1: f64,
    /// `myTolC2` (`hxx:126`).
    my_tol_c2: f64,
    /// `myMaxDerivOrderC1` (`hxx:127`).
    my_max_deriv_order_c1: i32,
    /// `myMaxDerivOrderC2` (`hxx:128`).
    my_max_deriv_order_c2: i32,
    /// `myUinfium` (`hxx:129`).
    my_uinfium: f64,
    /// `myUsupremum` (`hxx:130`).
    my_usupremum: f64,
    /// `myVinfium` (`hxx:131`).
    my_vinfium: f64,
    /// `myVsupremum` (`hxx:132`).
    my_vsupremum: f64,
}

impl<'a> GFuncExtCC2d<'a> {
    /// `Extrema_GFuncExtCC(theC1, theC2, theTol = 1.0e-10)`
    /// (`Extrema_GFuncExtCC.hxx:235-291`).
    pub(super) fn new(c1: &'a dyn Curve2d, c2: &'a dyn Curve2d, the_tol: f64) -> Self {
        let mut out = Self {
            my_c1: c1,
            my_c2: c2,
            my_tol: the_tol,
            my_u: 0.0,
            my_v: 0.0,
            my_p1: GpPnt2d::zero(),
            my_p2: GpPnt2d::zero(),
            my_du: GpVec2d::zero(),
            my_dv: GpVec2d::zero(),
            my_sq_dist: Vec::new(),
            my_points: Vec::new(),
            my_tol_c1: THE_MIN_TOL,
            my_tol_c2: THE_MIN_TOL,
            my_max_deriv_order_c1: 0,
            my_max_deriv_order_c2: 0,
            my_uinfium: 0.0,
            my_usupremum: 0.0,
            my_vinfium: 0.0,
            my_vsupremum: 0.0,
        };
        // `hxx:256-261`: SubIntervalInitialize over both curves' own ranges.
        out.sub_interval_initialize(
            c1.first_parameter(),
            c1.last_parameter(),
            c2.first_parameter(),
            c2.last_parameter(),
        );

        // `hxx:263-276` / `hxx:278-291`: the GetType switches.
        out.my_max_deriv_order_c1 = match adaptor_type2d(c1) {
            CurveType::BezierCurve
            | CurveType::BSplineCurve
            | CurveType::OffsetCurve
            | CurveType::OtherCurve => THE_MAX_ORDER,
            _ => 0,
        };
        out.my_tol_c1 = if out.my_max_deriv_order_c1 != 0 {
            out.search_of_tolerance(1)
        } else {
            THE_MIN_TOL
        };
        out.my_max_deriv_order_c2 = match adaptor_type2d(c2) {
            CurveType::BezierCurve
            | CurveType::BSplineCurve
            | CurveType::OffsetCurve
            | CurveType::OtherCurve => THE_MAX_ORDER,
            _ => 0,
        };
        out.my_tol_c2 = if out.my_max_deriv_order_c2 != 0 {
            out.search_of_tolerance(2)
        } else {
            THE_MIN_TOL
        };
        out
    }

    /// `SubIntervalInitialize(theUfirst, theUlast)`
    /// (`Extrema_GFuncExtCC.hxx:899-925`).
    fn sub_interval_initialize(&mut self, u_inf: f64, u_sup: f64, v_inf: f64, v_sup: f64) {
        self.my_uinfium = border(u_inf);
        self.my_usupremum = border(u_sup);
        self.my_vinfium = border(v_inf);
        self.my_vsupremum = border(v_sup);
    }

    /// `SearchOfTolerance(theC)` (`Extrema_GFuncExtCC.hxx:150-197`); `the_rank`
    /// selects the sub-interval and the curve the `D1` samples come from, which
    /// is what the `theC == myC1` / `theC == myC2` pointer tests decide.
    ///
    /// Unlike `Extrema_GFuncExtPC::SearchOfTolerance`
    /// (`Extrema_GFuncExtPC.hxx:434-444`) this one has no infinity guard: a
    /// sample whose `D1` is not finite contributes its (non-finite) magnitude.
    fn search_of_tolerance(&self, the_rank: i32) -> f64 {
        let (a_start_param, an_end_param) = if the_rank == 1 {
            (self.my_uinfium, self.my_usupremum)
        } else {
            (self.my_vinfium, self.my_vsupremum)
        };

        let a_step = (an_end_param - a_start_param) / (N_POINT as f64);

        let mut a_num = 0;
        let mut a_max = -INFINITE;

        loop {
            let mut u = a_start_param + (a_num as f64) * a_step;
            if u > an_end_param {
                u = an_end_param;
            }

            // `TheCurveTool1::D1` (`hxx:189`) -- both tools are
            // `Extrema_Curve2dTool` in this instantiation.
            let (_, v_der) = if the_rank == 1 {
                self.my_c1.d1(u)
            } else {
                self.my_c2.d1(u)
            };
            let vm = v_der.magnitude();
            if vm > a_max {
                a_max = vm;
            }

            a_num += 1;
            if a_num >= N_POINT + 1 {
                break;
            }
        }

        (a_max * THE_TOL_FACTOR).max(THE_MIN_TOL)
    }

    /// The `myMaxDerivOrderCx != 0 && |D1| <= myTolCx` arm of `Value`
    /// (`Extrema_GFuncExtCC.hxx:381-457` for the first curve, `:465-540` for the
    /// second): rebuild `myDu` / `myDv`. Returns the reconstructed derivative
    /// and its magnitude.
    fn reconstruct_d1(&self, the_rank: i32) -> (GpVec2d, f64) {
        let (u, u_inf, u_sup, tol, max_order) = if the_rank == 1 {
            (
                self.my_u,
                self.my_uinfium,
                self.my_usupremum,
                self.my_tol_c1,
                self.my_max_deriv_order_c1,
            )
        } else {
            (
                self.my_v,
                self.my_vinfium,
                self.my_vsupremum,
                self.my_tol_c2,
                self.my_max_deriv_order_c2,
            )
        };
        let c: &dyn Curve2d = if the_rank == 1 { self.my_c1 } else { self.my_c2 };

        const DIVISION_FACTOR: f64 = 1.0e-3; // `hxx:386` / `:469`
        let d = if u_sup >= REAL_LAST || u_inf <= REAL_FIRST {
            0.0
        } else {
            u_sup - u_inf
        };
        let a_delta = (d * DIVISION_FACTOR).max(THE_MIN_STEP);

        // `do { V = DN(++n); ... } while (!IsDeriveFound && n < myMaxDerivOrderCx)`.
        let mut n = 1;
        let mut v;
        let mut is_derive_found;
        loop {
            n += 1;
            v = c.eval_dn(u, n);
            is_derive_found = v.magnitude() > tol;
            if is_derive_found || n >= max_order {
                break;
            }
        }

        if is_derive_found {
            let uu = if u - u_inf < a_delta { u + a_delta } else { u - a_delta };
            let p1 = c.d0(u.min(uu));
            let p2 = c.d0(u.max(uu));
            let v1 = GpVec2d::new(p2.x() - p1.x(), p2.y() - p1.y());
            let out = if v.dot(&v1) < 0.0 { v.reversed() } else { v };
            let mag = out.magnitude();
            (out, mag)
        } else {
            // Three-point Taylor; `Ptemp` stays at the default `gp_Pnt2d(0, 0)`.
            let (grown, p1, p2, p3) = if u - u_inf < 2.0 * a_delta {
                (true, c.d0(u), c.d0(u + a_delta), c.d0(u + 2.0 * a_delta))
            } else {
                (
                    false,
                    c.d0(u - 2.0 * a_delta),
                    c.d0(u - a_delta),
                    c.d0(u),
                )
            };
            let v1 = GpVec2d::new(p1.x(), p1.y());
            let v2 = GpVec2d::new(p2.x(), p2.y());
            let v3 = GpVec2d::new(p3.x(), p3.y());
            let out = if grown {
                v1.multiplied_scalar(-3.0)
                    .added(&v2.multiplied_scalar(4.0))
                    .subtracted(&v3)
            } else {
                v1.subtracted(&v2.multiplied_scalar(4.0))
                    .added(&v3.multiplied_scalar(3.0))
            };
            let mag = out.magnitude();
            (out, mag)
        }
    }

    /// `Extrema_GFuncExtCC::Value(theUV, theF)`
    /// (`Extrema_GFuncExtCC.hxx:359-550`).
    fn value_impl(&mut self, the_u: f64, the_v: f64, the_f: &mut MathVector) -> bool {
        self.my_u = the_u;
        self.my_v = the_v;

        let (p1, du) = self.my_c1.d1(the_u);
        let (p2, dv) = self.my_c2.d1(the_v);
        self.my_p1 = p1;
        self.my_p2 = p2;
        self.my_du = du;
        self.my_dv = dv;

        let p1p2 = GpVec2d::new(p2.x() - p1.x(), p2.y() - p1.y());

        let mut ndu = du.magnitude();
        if self.my_max_deriv_order_c1 != 0 && ndu <= self.my_tol_c1 {
            let (v, mag) = self.reconstruct_d1(1);
            self.my_du = v;
            ndu = mag;
        }

        if ndu <= THE_MIN_TOL {
            // Warning: 1st derivative of the first curve is equal to zero!
            return false;
        }

        let mut ndv = dv.magnitude();
        if self.my_max_deriv_order_c2 != 0 && ndv <= self.my_tol_c2 {
            let (v, mag) = self.reconstruct_d1(2);
            self.my_dv = v;
            ndv = mag;
        }

        if ndv <= THE_MIN_TOL {
            return false;
        }

        the_f.set_value(1, p1p2.dot(&self.my_du) / ndu);
        the_f.set_value(2, p1p2.dot(&self.my_dv) / ndv);
        true
    }

    /// `Extrema_GFuncExtCC::Values(theUV, theF, theDF)`
    /// (`Extrema_GFuncExtCC.hxx:568-827`).
    fn values_impl(&mut self, the_uv: &MathVector, the_f: &mut MathVector, the_df: &mut MathMatrix) -> bool {
        if !self.value_impl(the_uv.value(1), the_uv.value(2), the_f) {
            return false;
        }
        let the_f1 = the_f.value(1);
        let the_f2 = the_f.value(2);

        let (p1, du, duu) = self.my_c1.d2(self.my_u);
        let (p2, dv, dvv) = self.my_c2.d2(self.my_v);
        self.my_p1 = p1;
        self.my_p2 = p2;

        let my_u_old = self.my_u;
        let my_v_old = self.my_v;
        let my_p1_old = self.my_p1;
        let my_p2_old = self.my_p2;
        let my_du_old = self.my_du;
        let my_dv_old = self.my_dv;

        const DIVISION_FACTOR: f64 = 0.01; // `hxx:614`

        let du_span = if self.my_usupremum >= REAL_LAST || self.my_uinfium <= REAL_FIRST {
            0.0
        } else {
            self.my_usupremum - self.my_uinfium
        };
        let a_delta_u = (du_span * DIVISION_FACTOR).max(THE_MIN_STEP);

        let dv_span = if self.my_vsupremum >= REAL_LAST || self.my_vinfium <= REAL_FIRST {
            0.0
        } else {
            self.my_vsupremum - self.my_vinfium
        };
        let a_delta_v = (dv_span * DIVISION_FACTOR).max(THE_MIN_STEP);

        let p1p2 = GpVec2d::new(self.my_p2.x() - self.my_p1.x(), self.my_p2.y() - self.my_p1.y());

        if self.my_max_deriv_order_c1 != 0 && du.magnitude() <= self.my_tol_c1 {
            // `hxx:636-722`.
            let mut ff1 = MathVector::new(1, 2);
            let mut ff2 = MathVector::new(1, 2);
            let mut ff3 = MathVector::new(1, 2);

            let (v1a, v2a, v3a);
            if my_u_old - self.my_uinfium < 2.0 * a_delta_u {
                v1a = the_f1;
                let uv2_1 = my_u_old + a_delta_u;
                let uv3_1 = my_u_old + 2.0 * a_delta_u;
                if !(self.value_impl(uv2_1, my_v_old, &mut ff2)
                    && self.value_impl(uv3_1, my_v_old, &mut ff3))
                {
                    return false;
                }
                v2a = ff2.value(1);
                v3a = ff3.value(1);
            } else {
                v3a = the_f1;
                let uv2_1 = my_u_old - a_delta_u;
                let uv1_1 = my_u_old - 2.0 * a_delta_u;
                if !(self.value_impl(uv2_1, my_v_old, &mut ff2)
                    && self.value_impl(uv1_1, my_v_old, &mut ff1))
                {
                    return false;
                }
                v1a = ff1.value(1);
                v2a = ff2.value(1);
            }
            let d11 = if my_u_old - self.my_uinfium < 2.0 * a_delta_u {
                (-3.0 * v1a + 4.0 * v2a - v3a) / (2.0 * a_delta_u)
            } else {
                (v1a - 4.0 * v2a + 3.0 * v3a) / (2.0 * a_delta_u)
            };
            the_df.set_value(1, 1, d11);

            self.my_u = my_u_old;
            self.my_v = my_v_old;

            let (w1, w2, w3);
            if my_v_old - self.my_vinfium < 2.0 * a_delta_v {
                w1 = the_f1;
                let uv2_2 = my_v_old + a_delta_v;
                let uv3_2 = my_v_old + 2.0 * a_delta_v;
                if !(self.value_impl(my_u_old, uv2_2, &mut ff2)
                    && self.value_impl(my_u_old, uv3_2, &mut ff3))
                {
                    return false;
                }
                w2 = ff2.value(1);
                w3 = ff3.value(1);
            } else {
                w3 = the_f1;
                let uv2_2 = my_v_old - a_delta_v;
                let uv1_2 = my_v_old - 2.0 * a_delta_v;
                if !(self.value_impl(my_u_old, uv2_2, &mut ff2)
                    && self.value_impl(my_u_old, uv1_2, &mut ff1))
                {
                    return false;
                }
                w1 = ff1.value(1);
                w2 = ff2.value(1);
            }
            let d12 = if my_v_old - self.my_vinfium < 2.0 * a_delta_v {
                (-3.0 * w1 + 4.0 * w2 - w3) / (2.0 * a_delta_v)
            } else {
                (w1 - 4.0 * w2 + 3.0 * w3) / (2.0 * a_delta_v)
            };
            the_df.set_value(1, 2, d12);

            self.my_u = my_u_old;
            self.my_v = my_v_old;
            self.my_p1 = my_p1_old;
            self.my_p2 = my_p2_old;
            self.my_du = my_du_old;
            self.my_dv = my_dv_old;
        } else {
            // `hxx:724-728`.
            let ndu = self.my_du.magnitude();
            the_df.set_value(
                1,
                1,
                -ndu + (p1p2.dot(&duu) / ndu)
                    - the_f1 * (self.my_du.dot(&duu) / (ndu * ndu)),
            );
            the_df.set_value(1, 2, self.my_dv.dot(&self.my_du) / ndu);
        }

        if self.my_max_deriv_order_c2 != 0 && dv.magnitude() <= self.my_tol_c2 {
            // `hxx:731-820`.
            let mut ff1 = MathVector::new(1, 2);
            let mut ff2 = MathVector::new(1, 2);
            let mut ff3 = MathVector::new(1, 2);

            let (v1a, v2a, v3a);
            if my_v_old - self.my_vinfium < 2.0 * a_delta_v {
                v1a = the_f2;
                let uv2_2 = my_v_old + a_delta_v;
                let uv3_2 = my_v_old + 2.0 * a_delta_v;
                if !(self.value_impl(my_u_old, uv2_2, &mut ff2)
                    && self.value_impl(my_u_old, uv3_2, &mut ff3))
                {
                    return false;
                }
                v2a = ff2.value(2);
                v3a = ff3.value(2);
            } else {
                v3a = the_f2;
                let uv2_2 = my_v_old - a_delta_v;
                let uv1_2 = my_v_old - 2.0 * a_delta_v;
                if !(self.value_impl(my_u_old, uv2_2, &mut ff2)
                    && self.value_impl(my_u_old, uv1_2, &mut ff1))
                {
                    return false;
                }
                v1a = ff1.value(2);
                v2a = ff2.value(2);
            }
            let d22 = if my_v_old - self.my_vinfium < 2.0 * a_delta_v {
                (-3.0 * v1a + 4.0 * v2a - v3a) / (2.0 * a_delta_v)
            } else {
                (v1a - 4.0 * v2a + 3.0 * v3a) / (2.0 * a_delta_v)
            };
            the_df.set_value(2, 2, d22);

            self.my_u = my_u_old;
            self.my_v = my_v_old;

            let (w1, w2, w3);
            if my_u_old - self.my_uinfium < 2.0 * a_delta_u {
                w1 = the_f2;
                let uv2_1 = my_u_old + a_delta_u;
                let uv3_1 = my_u_old + 2.0 * a_delta_u;
                if !(self.value_impl(uv2_1, my_v_old, &mut ff2)
                    && self.value_impl(uv3_1, my_v_old, &mut ff3))
                {
                    return false;
                }
                w2 = ff2.value(2);
                w3 = ff3.value(2);
            } else {
                w3 = the_f2;
                let uv2_1 = my_u_old - a_delta_u;
                let uv1_1 = my_u_old - 2.0 * a_delta_u;
                if !(self.value_impl(uv2_1, my_v_old, &mut ff2)
                    && self.value_impl(uv1_1, my_v_old, &mut ff1))
                {
                    return false;
                }
                w1 = ff1.value(2);
                w2 = ff2.value(2);
            }
            let d21 = if my_u_old - self.my_uinfium < 2.0 * a_delta_u {
                (-3.0 * w1 + 4.0 * w2 - w3) / (2.0 * a_delta_u)
            } else {
                (w1 - 4.0 * w2 + 3.0 * w3) / (2.0 * a_delta_u)
            };
            the_df.set_value(2, 1, d21);

            self.my_u = my_u_old;
            self.my_v = my_v_old;
            self.my_p1 = my_p1_old;
            self.my_p2 = my_p2_old;
            self.my_du = my_du_old;
            self.my_dv = my_dv_old;
        } else {
            // `hxx:822-827`.
            let ndv = self.my_dv.magnitude();
            the_df.set_value(
                2,
                2,
                ndv + (p1p2.dot(&dvv) / ndv)
                    - the_f2 * (self.my_dv.dot(&dvv) / (ndv * ndv)),
            );
            the_df.set_value(2, 1, -self.my_du.dot(&self.my_dv) / ndv);
        }

        true
    }

    /// `NbExt()` (`Extrema_GFuncExtCC.hxx:84-86`).
    pub(super) fn nb_ext(&self) -> usize {
        self.my_sq_dist.len()
    }

    /// `SquareDistance(theN)` (`Extrema_GFuncExtCC.hxx:88-89`), `the_n` 1-based.
    pub(super) fn square_distance(&self, the_n: usize) -> f64 {
        self.my_sq_dist[the_n - 1]
    }

    /// `Points(theN, P1, P2)` (`Extrema_GFuncExtCC.hxx:885-895`).
    pub(super) fn points(&self, the_n: usize) -> (POnCurv2d, POnCurv2d) {
        (self.my_points[2 * the_n - 2], self.my_points[2 * the_n - 1])
    }
}

impl<'a> MathFunctionSetWithDerivatives for GFuncExtCC2d<'a> {
    /// `NbVariables()` (`Extrema_GFuncExtCC.hxx:66`).
    fn nb_variables(&self) -> usize {
        2
    }

    /// `NbEquations()` (`Extrema_GFuncExtCC.hxx:69`).
    fn nb_equations(&self) -> usize {
        2
    }

    /// `Value(theUV, theF)` (`Extrema_GFuncExtCC.hxx:359-550`).
    fn value(&mut self, the_uv: &MathVector, the_f: &mut MathVector) -> bool {
        self.value_impl(the_uv.value(1), the_uv.value(2), the_f)
    }

    /// `Derivatives(theUV, theDF)` (`Extrema_GFuncExtCC.hxx:551-565`).
    fn derivatives(&mut self, the_uv: &MathVector, the_df: &mut MathMatrix) -> bool {
        let mut f = MathVector::new(1, 2);
        self.values_impl(the_uv, &mut f, the_df)
    }

    /// `Values(theUV, theF, theDF)` (`Extrema_GFuncExtCC.hxx:568-827`).
    fn values(&mut self, the_uv: &MathVector, the_f: &mut MathVector, the_df: &mut MathMatrix) -> bool {
        self.values_impl(the_uv, the_f, the_df)
    }

    /// `GetStateNumber()` (`Extrema_GFuncExtCC.hxx:832-870`).
    fn get_state_number(&mut self) -> i32 {
        let mut du = self.my_du;
        let mut dv = self.my_dv;
        let p1p2 = GpVec2d::new(self.my_p2.x() - self.my_p1.x(), self.my_p2.y() - self.my_p1.y());

        let m = du.magnitude();
        if m > self.my_tol_c1 {
            du = du.multiplied_scalar(1.0 / m);
        }
        let m = dv.magnitude();
        if m > self.my_tol_c2 {
            dv = dv.multiplied_scalar(1.0 / m);
        }

        if p1p2.dot(&du).abs() <= self.my_tol && p1p2.dot(&dv).abs() <= self.my_tol {
            self.my_sq_dist.push(self.my_p1.square_distance(&self.my_p2));
            self.my_points.push(POnCurv2d::new(self.my_u, self.my_p1));
            self.my_points.push(POnCurv2d::new(self.my_v, self.my_p2));
        }
        0
    }
}
