//! `Extrema_GenLocateExtCC` 2D instantiation (`Extrema_LocECC2d`) and its only
//! consumer `Extrema_LocateExtCC2d`.
//!
//! Source: `Extrema_LocECC2d.hxx:27-31` (the template alias),
//! `Extrema_GenLocateExtCC.hxx:47-110` (the constructor that performs the
//! search) and `Extrema_LocateExtCC2d.cxx:26-82` (the class the port's
//! `curve_curve::locate_extcc2d` wraps).
//!
//! The search itself is `math_FunctionSetRoot` over the
//! [`GFuncExtCC2d`](super::gfunc_ext_cc2d::GFuncExtCC2d) function set, exactly
//! as `Extrema_GenLocateExtCC` does (`Extrema_GenLocateExtCC.hxx:95-108`).

use occt_core::precision::CONFUSION;
use occt_math::{MathFunctionSetRoot, MathFunctionSetWithDerivatives, MathVector};

use super::prelude::*;

use super::curve2d_tool::adaptor_resolution2d;
use super::gfunc_ext_cc2d::GFuncExtCC2d;
use crate::geom2d_int::POnCurv2d;

/// `RealLast()` (`Standard_Real.hxx:179-183`), the `mySqDist` initial value
/// (`Extrema_GenLocateExtCC.hxx:57`).
const REAL_LAST: f64 = f64::MAX;

/// `Extrema_GenLocateExtCC<Adaptor2d_Curve2d, Extrema_Curve2dTool,
/// Extrema_POnCurv2d, Extrema_CCLocFOfLocECC2d>` = `Extrema_LocECC2d`
/// (`Extrema_LocECC2d.hxx:27-31`), whose whole body is the constructor
/// (`Extrema_GenLocateExtCC.hxx:47-110`). Like OCCT, the function set
/// (`TheCCLocF F`, `hxx:78`) is a local of that constructor: its accumulated
/// results are what the three members below hold.
pub(super) struct GenLocateExtCC2d {
    /// `myDone` (`hxx:138`).
    done: bool,
    /// `mySqDist` (`hxx:139`).
    sq_dist: f64,
    /// `myPoint1` (`hxx:140`).
    point1: POnCurv2d,
    /// `myPoint2` (`hxx:141`).
    point2: POnCurv2d,
}

impl GenLocateExtCC2d {
    /// `Extrema_GenLocateExtCC(theC1, theC2, theU0, theV0, theTolU, theTolV)`
    /// (`Extrema_GenLocateExtCC.hxx:47-110`).
    pub(super) fn new(
        c1: &dyn Curve2d,
        c2: &dyn Curve2d,
        u0: f64,
        v0: f64,
        tol_u: f64,
        tol_v: f64,
    ) -> Self {
        let mut done = false;
        let mut sq_dist = REAL_LAST;

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

        let mut f = GFuncExtCC2d::new(c1, c2, 1.0e-10); // `hxx:78`

        let tol = MathVector::from_slice(&[tol_u, tol_v]); // `hxx:79-81`
        let tol_f = 1.0e-10; // `hxx:82`

        let mut start = MathVector::from_slice(&[u0, v0]); // `hxx:84-88`
        let uu_inf = MathVector::from_slice(&[u_inf, v_inf]); // `hxx:90-91`
        let uu_sup = MathVector::from_slice(&[u_sup, v_sup]); // `hxx:92-93`

        // `math_FunctionSetRoot S(F, Tol); S.Perform(F, Start, Uuinf, Uusup);`
        // (`hxx:95-96`).
        let mut s = MathFunctionSetRoot::new(&f, &tol, 100);
        s.perform_with_bounds(&mut f, &start, &uu_inf, &uu_sup, false);

        let mut point1 = POnCurv2d::new(0.0, GpPnt2d::zero());
        let mut point2 = POnCurv2d::new(0.0, GpPnt2d::zero());

        if s.is_done() && f.nb_ext() > 0 {
            // `hxx:99-108`.
            sq_dist = f.square_distance(1);
            let (p1, p2) = f.points(1);
            point1 = p1;
            point2 = p2;

            start.set_value(1, point1.parameter());
            start.set_value(2, point2.parameter());
            let mut ff = MathVector::new(1, 2);
            f.value(&start, &mut ff);
            if ff.value(1) < tol_f && ff.value(2) < tol_f {
                done = true;
            }
        }

        Self {
            done,
            sq_dist,
            point1,
            point2,
        }
    }

    /// `IsDone()` (`hxx:112`).
    pub(super) fn is_done(&self) -> bool {
        self.done
    }

    /// `SquareDistance()` (`hxx:115-123`); OCCT throws `StdFail_NotDone` when
    /// not done, the port reports the raw value.
    pub(super) fn square_distance(&self) -> f64 {
        self.sq_dist
    }

    /// `Point(theP1, theP2)` (`hxx:127-134`); same not-done convention.
    pub(super) fn point(&self) -> (POnCurv2d, POnCurv2d) {
        (self.point1, self.point2)
    }
}

/// `Extrema_LocateExtCC2d` (`Extrema_LocateExtCC2d.hxx:31-62`,
/// `Extrema_LocateExtCC2d.cxx:26-82`): the two tolerances are the curves'
/// `Extrema_Curve2dTool::Resolution(.., Precision::Confusion())`.
pub struct ExtremaLocateExtCC2d {
    /// `myDone` (`hxx:58`).
    done: bool,
    /// `mySqDist` (`hxx:59`); `RealLast()` before the search.
    sq_dist: f64,
    /// `myPoint1` (`hxx:60`).
    point1: POnCurv2d,
    /// `myPoint2` (`hxx:61`).
    point2: POnCurv2d,
}

impl ExtremaLocateExtCC2d {
    /// `Extrema_LocateExtCC2d(C1, C2, U0, V0)`
    /// (`Extrema_LocateExtCC2d.cxx:26-50`).
    pub fn new(c1: &dyn Curve2d, c2: &dyn Curve2d, u0: f64, v0: f64) -> Self {
        // `cxx:32-33`.
        let tol_u = adaptor_resolution2d(c1, CONFUSION);
        let tol_v = adaptor_resolution2d(c2, CONFUSION);

        // `cxx:39`: `Extrema_LocECC2d Xtrem(C1, C2, U0, V0, TolU, TolV)`.
        let xtrem = GenLocateExtCC2d::new(c1, c2, u0, v0, tol_u, tol_v);

        // `cxx:42-49`.
        let done = xtrem.is_done();
        if done {
            let (p1, p2) = xtrem.point();
            Self {
                done,
                sq_dist: xtrem.square_distance(),
                point1: p1,
                point2: p2,
            }
        } else {
            Self {
                done,
                sq_dist: REAL_LAST,
                point1: POnCurv2d::new(0.0, GpPnt2d::zero()),
                point2: POnCurv2d::new(0.0, GpPnt2d::zero()),
            }
        }
    }

    /// `IsDone()` (`cxx:56-59`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `SquareDistance()` (`cxx:64-70`); OCCT throws `StdFail_NotDone`.
    pub fn square_distance(&self) -> f64 {
        self.sq_dist
    }

    /// `Point(P1, P2)` (`cxx:76-82`); same not-done convention.
    pub fn point(&self) -> (POnCurv2d, POnCurv2d) {
        (self.point1, self.point2)
    }
}
