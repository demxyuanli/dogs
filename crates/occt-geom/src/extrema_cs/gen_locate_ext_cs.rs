//! `Extrema_GenLocateExtCS` — local curve/surface extrema from an initial
//! `(T, U, V)` triple. Source: `Extrema_GenLocateExtCS.cxx:28-139`,
//! `.hxx:31-77`.
//!
//! ```text
//! Perform -> math_FunctionSetRoot over Extrema_FuncExtCS bounded by the
//!            curve and surface ranges (cxx:63-89), then keep solution #1
//!            (cxx:96-98)
//! ```
//!
//! There is no global search here: the caller supplies the seed, exactly as
//! `Extrema_GenExtCS::GlobMinConicS` does when it re-polishes a "bad" solution
//! (`Extrema_GenExtCS.cxx:688-706`).

use occt_math::{MathFunctionSetRoot, MathVector};

use super::func_ext_cs::FuncExtCS;
use super::prelude::*;

/// `Extrema_GenLocateExtCS` (`Extrema_GenLocateExtCS.hxx:31-77`).
pub struct GenLocateExtCS<'a> {
    /// `myDone`.
    done: bool,
    /// `mySqDist`.
    sq_dist: f64,
    /// `myPoint1` / `myPoint2`.
    point1: ExtremaPOnCurv,
    point2: ExtremaPOnSurf,
    /// Keeps the lifetime of the curve/surface handles used by the last
    /// `Perform` alive (OCCT stores them inside `myF`).
    _marker: std::marker::PhantomData<&'a ()>,
}

impl<'a> GenLocateExtCS<'a> {
    /// `Extrema_GenLocateExtCS()` (`cxx:28-33`): `myDone = false`,
    /// `mySqDist = RealLast()`.
    pub fn new() -> Self {
        GenLocateExtCS {
            done: false,
            sq_dist: f64::MAX,
            point1: ExtremaPOnCurv::new(
                0.0,
                GpPnt::new(0.0, 0.0, 0.0),
            ),
            point2: ExtremaPOnSurf::new(0.0, 0.0, GpPnt::new(0.0, 0.0, 0.0)),
            _marker: std::marker::PhantomData,
        }
    }

    /// `Extrema_GenLocateExtCS(C, S, T, U, V, Tol1, Tol2)` (`cxx:36-47`).
    pub fn with_curve_surface(
        c: &'a dyn Curve,
        s: &'a dyn Surface,
        t: f64,
        u: f64,
        v: f64,
        tol1: f64,
        tol2: f64,
    ) -> Self {
        let mut me = GenLocateExtCS::new();
        me.perform(c, s, t, u, v, tol1, tol2);
        me
    }

    /// `Perform(C, S, T, U, V, Tol1, Tol2)` (`cxx:50-99`).
    pub fn perform(
        &mut self,
        c: &'a dyn Curve,
        s: &'a dyn Surface,
        t: f64,
        u: f64,
        v: f64,
        tol1: f64,
        tol2: f64,
    ) {
        self.done = false;

        let mut f = FuncExtCS::with_curve_surface(c, s);

        let mut tol = MathVector::new(1, 3);
        tol.set_value(1, tol1);
        tol.set_value(2, tol2);
        tol.set_value(3, tol2);

        let mut start = MathVector::new(1, 3);
        start.set_value(1, t);
        start.set_value(2, u);
        start.set_value(3, v);

        let mut b_inf = MathVector::new(1, 3);
        b_inf.set_value(1, c.first_parameter());
        b_inf.set_value(2, s.u_range().0);
        b_inf.set_value(3, s.v_range().0);

        let mut b_sup = MathVector::new(1, 3);
        b_sup.set_value(1, c.last_parameter());
        b_sup.set_value(2, s.u_range().1);
        b_sup.set_value(3, s.v_range().1);

        let mut sr = MathFunctionSetRoot::new(&f, &tol, 100);
        sr.perform_with_bounds(&mut f, &mut start, &b_inf, &b_sup, false);
        if !sr.is_done() {
            return;
        }

        // `cxx:96-98`: keep solution #1.
        self.sq_dist = f.square_distance(1);
        self.point1 = f.point_on_curve(1);
        self.point2 = f.point_on_surface(1);
        self.done = true;
    }

    /// `IsDone()` (`cxx:102-105`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `SquareDistance()` (`cxx:110-117`): OCCT throws `StdFail_NotDone`.
    pub fn square_distance(&self) -> f64 {
        assert!(
            self.is_done(),
            "Extrema_GenLocateExtCS::SquareDistance: StdFail_NotDone"
        );
        self.sq_dist
    }

    /// `PointOnCurve()` (`cxx:120-127`).
    pub fn point_on_curve(&self) -> ExtremaPOnCurv {
        assert!(
            self.is_done(),
            "Extrema_GenLocateExtCS::PointOnCurve: StdFail_NotDone"
        );
        self.point1
    }

    /// `PointOnSurface()` (`cxx:132-139`).
    pub fn point_on_surface(&self) -> ExtremaPOnSurf {
        assert!(
            self.is_done(),
            "Extrema_GenLocateExtCS::PointOnSurface: StdFail_NotDone"
        );
        self.point2
    }
}

impl Default for GenLocateExtCS<'_> {
    fn default() -> Self {
        GenLocateExtCS::new()
    }
}
