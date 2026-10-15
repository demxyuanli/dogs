//! `Extrema_GenLocateExtCC` 3D instantiation (`Extrema_LocECC`) and its only
//! consumer `Extrema_LocateExtCC`.
//!
//! Source: `Extrema_LocECC.hxx:27-30` (the template alias),
//! `Extrema_GenLocateExtCC.hxx:49-110` (the constructor that performs the
//! search) and `Extrema_LocateExtCC.cxx:24-74` (the class the port's
//! `curve_curve::locate_extcc` wraps).
//!
//! The search itself is `math_FunctionSetRoot` over the
//! [`GFuncExtCC`](super::gfunc_ext_cc::GFuncExtCC) function set, exactly
//! as `Extrema_GenLocateExtCC` does (`Extrema_GenLocateExtCC.hxx:95-96`).
//! The 2D twin of this file is `occt_geom2d::extrema2d::gen_locate_ext_cc2d`.

use occt_core::precision::CONFUSION;
use occt_math::{MathFunctionSetRoot, MathFunctionSetWithDerivatives, MathVector};

use super::prelude::*;

use super::gfunc_ext_cc::GFuncExtCC;
use super::p_on_curv::POnCurv;

/// `RealLast()` (`Standard_Real.hxx:179-183`), the `mySqDist` initial value
/// (`Extrema_GenLocateExtCC.hxx:57`).
const REAL_LAST: f64 = f64::MAX;

/// `Extrema_GenLocateExtCC<Adaptor3d_Curve, Extrema_CurveTool,
/// Extrema_POnCurv, Extrema_CCLocFOfLocECC>` = `Extrema_LocECC`
/// (`Extrema_LocECC.hxx:27-30`), whose whole body is the constructor
/// (`Extrema_GenLocateExtCC.hxx:49-110`). Like OCCT, the function set
/// (`TheCCLocF F`, `hxx:78`) is a local of that constructor: its accumulated
/// results are what the three members below hold.
pub(super) struct GenLocateExtCC {
    /// `myDone` (`hxx:139`).
    my_done: bool,
    /// `mySqDist` (`hxx:140`).
    my_sq_dist: f64,
    /// `myPoint1` (`hxx:141`).
    my_point1: POnCurv,
    /// `myPoint2` (`hxx:142`).
    my_point2: POnCurv,
}

impl GenLocateExtCC {
    /// `Extrema_GenLocateExtCC(theC1, theC2, theU0, theV0, theTolU, theTolV)`
    /// (`Extrema_GenLocateExtCC.hxx:49-110`).
    pub(super) fn new(
        c1: &dyn Curve,
        c2: &dyn Curve,
        u0: f64,
        v0: f64,
        tol_u: f64,
        tol_v: f64,
    ) -> Self {
        let mut my_done = false;
        let mut my_sq_dist = REAL_LAST;

        // `hxx:59-77`: the axis-aligned bounds, each swapped when reversed.
        let mut u_inf = c1.first_parameter();
        let mut u_sup = c1.last_parameter();
        if u_inf > u_sup {
            std::mem::swap(&mut u_inf, &mut u_sup);
        }
        let mut v_inf = c2.first_parameter();
        let mut v_sup = c2.last_parameter();
        if v_inf > v_sup {
            std::mem::swap(&mut v_inf, &mut v_sup);
        }

        let mut f = GFuncExtCC::new(c1, c2, 1.0e-10); // `hxx:78`

        let tol = MathVector::from_slice(&[tol_u, tol_v]); // `hxx:79-81`
        let tol_f = 1.0e-10; // `hxx:82`

        let mut start = MathVector::from_slice(&[u0, v0]); // `hxx:84-89`
        let uu_inf = MathVector::from_slice(&[u_inf, v_inf]); // `hxx:90-91`
        let uu_sup = MathVector::from_slice(&[u_sup, v_sup]); // `hxx:92-93`

        // `math_FunctionSetRoot S(F, Tol); S.Perform(F, Start, Uuinf, Uusup);`
        // (`hxx:95-96`).
        let mut s = MathFunctionSetRoot::new(&f, &tol, 100);
        s.perform_with_bounds(&mut f, &start, &uu_inf, &uu_sup, false);

        let mut my_point1 = POnCurv::new(0.0, GpPnt::zero());
        let mut my_point2 = POnCurv::new(0.0, GpPnt::zero());

        if s.is_done() && f.nb_ext() > 0 {
            // `hxx:99-109`.
            my_sq_dist = f.square_distance(1);
            let (p1, p2) = f.points(1);
            my_point1 = p1;
            my_point2 = p2;

            start.set_value(1, my_point1.parameter());
            start.set_value(2, my_point2.parameter());
            let mut ff = MathVector::new(1, 2);
            f.value(&start, &mut ff);
            if ff.value(1) < tol_f && ff.value(2) < tol_f {
                my_done = true;
            }
        }

        Self {
            my_done,
            my_sq_dist,
            my_point1,
            my_point2,
        }
    }

    /// `IsDone()` (`hxx:113`).
    pub(super) fn is_done(&self) -> bool {
        self.my_done
    }

    /// `SquareDistance()` (`hxx:116-123`); OCCT throws `StdFail_NotDone` when
    /// not done, the port reports the raw value.
    pub(super) fn square_distance(&self) -> f64 {
        self.my_sq_dist
    }

    /// `Point(theP1, theP2)` (`hxx:128-136`); same not-done convention.
    pub(super) fn point(&self) -> (POnCurv, POnCurv) {
        (self.my_point1, self.my_point2)
    }
}

/// `Extrema_LocateExtCC` (`Extrema_LocateExtCC.hxx:29-59`,
/// `Extrema_LocateExtCC.cxx:24-74`): the two tolerances are the curves'
/// `Extrema_CurveTool::Resolution(.., Precision::Confusion())`.
pub struct ExtremaLocateExtCC {
    /// `myDone` (`hxx:56`).
    my_done: bool,
    /// `mySqDist` (`hxx:57`); `RealLast()` before the search.
    my_sq_dist: f64,
    /// `myPoint1` (`hxx:58`).
    my_point1: POnCurv,
    /// `myPoint2` (`hxx:59`).
    my_point2: POnCurv,
}

impl ExtremaLocateExtCC {
    /// `Extrema_LocateExtCC(C1, C2, U0, V0)`
    /// (`Extrema_LocateExtCC.cxx:24-47`).
    pub fn new(c1: &dyn Curve, c2: &dyn Curve, u0: f64, v0: f64) -> Self {
        // `cxx:30-31`.
        let tol_u = c1.resolution(CONFUSION);
        let tol_v = c2.resolution(CONFUSION);

        // `cxx:36`: `Extrema_LocECC Xtrem(C1, C2, U0, V0, TolU, TolV)`.
        let xtrem = GenLocateExtCC::new(c1, c2, u0, v0, tol_u, tol_v);

        // `cxx:39-46`.
        let my_done = xtrem.is_done();
        if my_done {
            let (p1, p2) = xtrem.point();
            Self {
                my_done,
                my_sq_dist: xtrem.square_distance(),
                my_point1: p1,
                my_point2: p2,
            }
        } else {
            Self {
                my_done,
                my_sq_dist: REAL_LAST,
                my_point1: POnCurv::new(0.0, GpPnt::zero()),
                my_point2: POnCurv::new(0.0, GpPnt::zero()),
            }
        }
    }

    /// `IsDone()` (`cxx:49-53`).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// `SquareDistance()` (`cxx:55-63`); OCCT throws `StdFail_NotDone`.
    pub fn square_distance(&self) -> f64 {
        self.my_sq_dist
    }

    /// `Point(P1, P2)` (`cxx:65-74`); same not-done convention.
    pub fn point(&self) -> (POnCurv, POnCurv) {
        (self.my_point1, self.my_point2)
    }
}
