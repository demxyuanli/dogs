//! `Extrema_EPCOfExtPC2d` — all extremal distances between a point and a 2D
//! curve, sampled into `NbSample` sub-intervals.
//!
//! Source: `ModelingData/TKGeomBase/Extrema/Extrema_GGenExtPC.hxx:36-232`
//! (OCCT 8.0.0); `Extrema_EPCOfExtPC2d` is the template instantiation
//! `Extrema_EPCOfExtPC2d.hxx:28-32`. Ported: `Initialize` `:98-110` (full
//! range), the explicit-range `Initialize` `:119-133`, `Perform` `:160-173` and
//! the read-backs `:176-225`.
//!
//! This is the engine `BRepBlend_BlendTool::Project`
//! (`BRepBlend_BlendTool.cxx:52-55`) drives directly, *not* `Extrema_ExtPC2d`:
//! `Extrema_ExtPC2d` dispatches on the curve type and takes the analytic
//! conic arms (`Extrema_GGExtPC.hxx:160-181`), while this class always runs the
//! sampled root search over the whole curve. The two report different extrema
//! counts on conics, so [`super::ext_pc2d::ExtremaExtPC2d`] is not a substitute.
//!
//! `math_FunctionRoots` collects the roots through
//! `Extrema_GFuncExtPC::GetStateNumber`; the port's
//! [`super::gfunc_ext_pc2d::FuncExtPc2d`] implements that by pushing onto
//! `sols`, and `NbExt` / `SquareDistance` / `IsMin` / `Point` read it back.

use occt_core::gp::GpPnt2d;
use occt_core::math_function_roots::FunctionRoots;

use super::gfunc_ext_pc2d::FuncExtPc2d;
use super::prelude::*;

/// `math_FunctionRoots(F, A, B, NbSample, EpsX, EpsF, EpsNull)` — the
/// 7-argument form `Perform` uses (`Extrema_GGenExtPC.hxx:166`), so `K` takes
/// its `math_FunctionRoots` default (`math_FunctionRoots.hxx:105`).
const K_DEFAULT: f64 = 0.0;

/// `Extrema_EPCOfExtPC2d` = `Extrema_GGenExtPC<...>` over the 2D point-curve
/// function of `super::gfunc_ext_pc2d`.
pub struct EPcOfExtPc2d<'a> {
    /// `myF` (`hxx:234`).
    f: FuncExtPc2d<'a>,
    /// `myDone` (`hxx:227`).
    done: bool,
    /// `mynbsample` (`hxx:229`).
    nb_sample: i32,
    /// `myumin` (`hxx:230`).
    u_min: f64,
    /// `myusup` (`hxx:231`).
    u_sup: f64,
    /// `mytolu` (`hxx:232`).
    tol_u: f64,
    /// `mytolF` (`hxx:233`).
    tol_f: f64,
}

impl<'a> EPcOfExtPc2d<'a> {
    /// `Extrema_GGenExtPC(P, C, NbSample, TolU, TolF)` (`hxx:62-71`): sampled
    /// over the curve's own `[First, Last]`.
    pub fn new(p: &GpPnt2d, c: &'a dyn Curve2d, nb_sample: i32, tol_u: f64, tol_f: f64) -> Self {
        let first = c.first_parameter();
        let last = c.last_parameter();
        let mut out = Self::with_function(p, c, first, last);
        // `Initialize(theC, theNbSample, theTolU, theTolF)` (`hxx:68`).
        out.nb_sample = nb_sample;
        out.tol_u = tol_u;
        out.tol_f = tol_f;
        out.u_min = first;
        out.u_sup = last;
        // `Perform(theP)` (`hxx:69`).
        out.perform(p);
        out
    }

    /// `Extrema_GGenExtPC(P, C, NbSample, Umin, Usup, TolU, TolF)`
    /// (`hxx:80-91`).
    pub fn new_range(
        p: &GpPnt2d,
        c: &'a dyn Curve2d,
        nb_sample: i32,
        u_min: f64,
        u_sup: f64,
        tol_u: f64,
        tol_f: f64,
    ) -> Self {
        let first = c.first_parameter();
        let last = c.last_parameter();
        let mut out = Self::with_function(p, c, first, last);
        // `Initialize(theC, theNbSample, theUmin, theUsup, theTolU, theTolF)`
        // (`hxx:89`).
        out.nb_sample = nb_sample;
        out.tol_u = tol_u;
        out.tol_f = tol_f;
        out.u_min = u_min;
        out.u_sup = u_sup;
        out.perform(p);
        out
    }

    /// `myF(theP, theC)` (`hxx:67`, `:86`): `Extrema_GFuncExtPC(theP, theC)`
    /// (`Extrema_GFuncExtPC.hxx:76-101`) computes the singular tolerance over
    /// the curve's own range, and the sampled range starts there too.
    fn with_function(p: &GpPnt2d, c: &'a dyn Curve2d, first: f64, last: f64) -> Self {
        Self {
            f: FuncExtPc2d::new(c, p, first, last),
            done: false,
            nb_sample: 0,
            u_min: first,
            u_sup: last,
            tol_u: 0.0,
            tol_f: 0.0,
        }
    }

    /// `Perform(theP)` (`hxx:160-173`).
    pub fn perform(&mut self, p: &GpPnt2d) {
        self.f.set_point(p);
        self.f.sub_interval_initialize(self.u_min, self.u_sup);
        self.done = false;

        let roots = FunctionRoots::new(
            &mut self.f,
            self.u_min,
            self.u_sup,
            self.nb_sample,
            self.tol_u,
            self.tol_f,
            self.tol_f,
            K_DEFAULT,
        );
        // `if (!S.IsDone() || S.IsAllNull()) return;` (`hxx:168-171`).
        if !roots.is_done() || roots.is_all_null() {
            return;
        }
        self.done = true;
    }

    /// `IsDone()` (`hxx:176`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `NbExt()` (`hxx:180-187`). OCCT raises `StdFail_NotDone` when the search
    /// did not finish; every OCCT caller tests `IsDone()` first
    /// (`BRepBlend_BlendTool.cxx:64-66`), so the port returns the collected
    /// count.
    pub fn nb_ext(&self) -> usize {
        self.f.sols().len()
    }

    /// `SquareDistance(theN)` (`hxx:192-198`), 1-based.
    ///
    /// OCCT raises `Standard_OutOfRange` outside `1..=NbExt()`; the port returns
    /// `0.0` for an index the search never produced, mirroring
    /// [`super::ext_pc2d::ExtremaExtPC2d::square_distance`].
    pub fn square_distance(&self, n: usize) -> f64 {
        match n.checked_sub(1).and_then(|i| self.f.sols().get(i)) {
            Some(s) => s.sq_dist,
            None => 0.0,
        }
    }

    /// `IsMin(theN)` (`hxx:203-211`), 1-based, with the same out-of-range
    /// convention as [`Self::square_distance`].
    pub fn is_min(&self, n: usize) -> bool {
        match n.checked_sub(1).and_then(|i| self.f.sols().get(i)) {
            Some(s) => s.is_min,
            None => false,
        }
    }

    /// `Point(theN)` (`hxx:215-225`), 1-based: the parameter on the curve and
    /// the point on it, with the same out-of-range convention as
    /// [`Self::square_distance`].
    pub fn point(&self, n: usize) -> (f64, GpPnt2d) {
        match n.checked_sub(1).and_then(|i| self.f.sols().get(i)) {
            Some(s) => (s.u, s.point),
            None => (0.0, GpPnt2d::zero()),
        }
    }
}
