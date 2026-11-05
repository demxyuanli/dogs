//! `Extrema_GenLocateExtPS` + `Extrema_FuncPSDist` — local point / surface
//! extremum around a starting `(U, V)`.
//!
//! Source: `ModelingData/TKGeomBase/Extrema/Extrema_GenLocateExtPS.cxx`
//! (OCCT 8.0.0, 253 lines), `Extrema_GenLocateExtPS.hxx:30-84`,
//! `Extrema_FuncPSDist.cxx` (91 lines) and `Extrema_FuncPSDist.hxx:30-70`.
//! Line references below are to the `.cxx` unless another file is named.
//!
//! Two branches, selected by `isDistanceCriteria` (`hxx:47-49`, default
//! `false`):
//!
//! - `true`, the distance criteria (`cxx:146-178`): `math_BFGS` on
//!   [`ExtremaFuncPsDist`], falling back to `math_FRPR` when the BFGS status is
//!   not `math_OK`. The objective is the square distance from the point to the
//!   surface; no normal equation is solved.
//! - `false`, the normal projection criteria (`cxx:180-223`): the
//!   `Extrema_FuncPSNorm` two-equation system (`(S-P).Su = 0`, `(S-P).Sv = 0`)
//!   through `math_FunctionSetRoot`, retried with
//!   `math_NewtonFunctionSetRoot` when the first solver did not converge or
//!   when [`correct_tol`] widened the tolerances. The extremum kept is the one
//!   of smallest square distance among the stored ones (`cxx:211-222`).
//!
//! Divergences that are not heuristics:
//! - `math_BFGS::Perform` / `math_FRPR::Perform` report `math_FunctionError`
//!   through `math_Status`; the port's solvers take infallible closures, so
//!   [`ExtremaFuncPsDist`] answers `f64::INFINITY` outside the parametric box
//!   instead of refusing the value (`Extrema_FuncPSDist.cxx:41-44`, `:53-56`,
//!   `:67-70`), and [`GenLocateExtPs::perform`] then checks the returned
//!   location against that same box. An optimum outside the box is the
//!   "not done" answer OCCT gives when the objective refuses the value.
//! - `Extrema_POnSurf::Parameter(U, V)` (`Extrema_POnSurf.hxx:61-65`) is the
//!   `(f64, f64)` of [`ExtremaPOnSurf::parameter`].
//! - `StdFail_NotDone` of `SquareDistance` / `Point` (`cxx:236-253`) is the
//!   `Option` of the two accessors; OCCT raises, the callers of this class in
//!   the port only read them after `is_done`.

use crate::approx_same_parameter;
use crate::extrema_surf::{ExtremaFuncPsNorm, ExtremaPOnSurf};
use crate::surface::Surface;

use occt_core::gp::{GpPnt, GpVec};
use occt_core::kernel::string_tools::real_to_int;
use occt_core::precision::{self, CONFUSION, PCONFUSION};
use occt_math::frpr::Frpr;
use occt_math::{MathFunctionSetRoot, MathNewtonFunctionSetRoot, MathVector, BFGS};

// ---------------------------------------------------------------------------
// Extrema_FuncPSDist (Extrema_FuncPSDist.hxx:30-70, .cxx:29-91)
// ---------------------------------------------------------------------------

/// `Extrema_FuncPSDist` (`hxx:30-70`): the square distance from one point to a
/// point of the surface, and its gradient in `(U, V)`.
///
/// The class is a `math_MultipleVarFunctionWithGradient` in OCCT; the port's
/// `math_*` minimization entry points take plain closures, so the three bodies
/// stay inherent methods and the call sites wrap them.
pub struct ExtremaFuncPsDist<'a> {
    /// `mySurf` (`hxx:66`).
    surf: &'a dyn Surface,
    /// `myP` (`hxx:67`).
    p: GpPnt,
}

impl<'a> ExtremaFuncPsDist<'a> {
    /// `Extrema_FuncPSDist(S, P)` (`cxx:29-34`).
    pub fn new(surf: &'a dyn Surface, p: GpPnt) -> Self {
        ExtremaFuncPsDist { surf, p }
    }

    /// `NbVariables()` (`cxx:38-41`): `2`.
    pub fn nb_variables(&self) -> usize {
        2
    }

    /// `IsInside(X)` (`cxx:85-91`): both parameters inside the surface box.
    pub fn is_inside(&self, x: &MathVector) -> bool {
        let (umin, umax) = self.surf.u_range();
        let (vmin, vmax) = self.surf.v_range();
        x.value(1) >= umin && x.value(1) <= umax && x.value(2) >= vmin && x.value(2) <= vmax
    }

    /// `Value(X, F)` (`cxx:44-51`): `F = |S(U, V) - P|^2`.
    ///
    /// `false` of the `.cxx` (outside the box) is `f64::INFINITY` — see the
    /// module note.
    pub fn value(&self, x: &MathVector) -> f64 {
        if !self.is_inside(x) {
            return f64::INFINITY;
        }
        self.surf
            .value(x.value(1), x.value(2))
            .square_distance(&self.p)
    }

    /// `Gradient(X, G)` (`cxx:55-69`): `G(1) = P1P2.Du`, `G(2) = P1P2.Dv`
    /// with `P1P2 = gp_Vec(S(U, V), P)`.
    ///
    /// Outside the box the `.cxx` answers `false`; the port answers a zero
    /// gradient so the solver's own progress test sees no descent.
    pub fn gradient(&self, x: &MathVector) -> MathVector {
        let mut g = MathVector::new(1, 2);
        if !self.is_inside(x) {
            g.init(0.0);
            return g;
        }
        let (s, du, dv) = self.surf.d1(x.value(1), x.value(2));
        let p1p2 = GpVec::from_pnts(&s, &self.p);
        g.set_value(1, p1p2.dot(&du));
        g.set_value(2, p1p2.dot(&dv));
        g
    }

    /// `Values(X, F, G)` (`cxx:73-82`): the two expressions above, computed
    /// from the one `D1` evaluation of the `.cxx` body.
    pub fn values(&self, x: &MathVector) -> (f64, MathVector) {
        let mut g = MathVector::new(1, 2);
        if !self.is_inside(x) {
            g.init(0.0);
            return (f64::INFINITY, g);
        }
        let (s, du, dv) = self.surf.d1(x.value(1), x.value(2));
        let p1p2 = GpVec::from_pnts(&s, &self.p);
        g.set_value(1, p1p2.dot(&du));
        g.set_value(2, p1p2.dot(&dv));
        (s.square_distance(&self.p), g)
    }
}

// ---------------------------------------------------------------------------
// Extrema_GenLocateExtPS (Extrema_GenLocateExtPS.hxx:30-84, .cxx:29-253)
// ---------------------------------------------------------------------------

/// `static void CorrectTol(U0, V0, Tol)` (`cxx:29-59`): widen the parametric
/// tolerances for large `U0` / `V0`, one decade per order of magnitude of
/// `Epsilon(U0)` above `Epsilon(1)`.
fn correct_tol(u0: f64, v0: f64, tol: &mut MathVector) {
    let tol_ref = PCONFUSION;
    let eps_ref = precision::epsilon(1.0);
    let tolog10 = 0.43429;

    let epsu = precision::epsilon(u0);
    if epsu > eps_ref {
        let n = real_to_int(tolog10 * ((epsu / eps_ref).ln()) + 1.0) + 1;
        let mut t = tol_ref;
        for _ in 1..=n {
            t *= 10.0;
        }
        tol.set_value(1, tol.value(1).max(t));
    }

    let epsv = precision::epsilon(v0);
    if epsv > eps_ref {
        let n = real_to_int(tolog10 * ((epsv / eps_ref).ln()) + 1.0) + 1;
        let mut t = tol_ref;
        for _ in 1..=n {
            t *= 10.0;
        }
        tol.set_value(2, tol.value(2).max(t));
    }
}

/// `Extrema_GenLocateExtPS` (`hxx:30-84`): the local extremum of the distance
/// from a point to a surface, searched around a starting `(U, V)`.
pub struct GenLocateExtPs<'a> {
    /// `mySurf` (`hxx:72`).
    surf: &'a dyn Surface,
    /// `myTolU` (`hxx:73`).
    tol_u: f64,
    /// `myTolV` (`hxx:73`).
    tol_v: f64,
    /// `myDone` (`hxx:76`).
    done: bool,
    /// `mySqDist` (`hxx:79`).
    sq_dist: f64,
    /// `myPoint` (`hxx:80`).
    point: ExtremaPOnSurf,
}

impl<'a> GenLocateExtPs<'a> {
    /// `Extrema_GenLocateExtPS(S, TolU = PConfusion, TolV = PConfusion)`
    /// (`cxx:107-117`).
    pub fn new(surf: &'a dyn Surface, tol_u: f64, tol_v: f64) -> Self {
        GenLocateExtPs {
            surf,
            tol_u,
            tol_v,
            done: false,
            sq_dist: -1.0,
            point: ExtremaPOnSurf::default(),
        }
    }

    /// `Extrema_GenLocateExtPS(S)` with both tolerances at
    /// `Precision::PConfusion()` (`hxx:37-39`).
    pub fn with_default_tolerance(surf: &'a dyn Surface) -> Self {
        Self::new(surf, PCONFUSION, PCONFUSION)
    }

    /// `IsMinDist(P, S, U0, V0)` (`cxx:63-105`): true when none of the eight
    /// neighbours of `(U0, V0)` at one `10 * Confusion` resolution step is
    /// closer to `P`.
    pub fn is_min_dist(p: &GpPnt, surf: &dyn Surface, u0: f64, v0: f64) -> bool {
        let du = approx_same_parameter::u_resolution(surf, 10.0 * CONFUSION).max(10.0 * PCONFUSION);
        let dv = approx_same_parameter::v_resolution(surf, 10.0 * CONFUSION).max(10.0 * PCONFUSION);
        let p0 = surf.value(u0, v0);
        let d0 = p.square_distance(&p0);
        let (umin, umax) = surf.u_range();
        let (vmin, vmax) = surf.v_range();
        for iu in -1..=1 {
            let mut u = u0 + f64::from(iu) * du;
            if !surf.is_u_periodic() {
                u = u.max(umin).min(umax);
            }
            for iv in -1..=1 {
                if iu == 0 && iv == 0 {
                    continue;
                }
                let mut v = v0 + f64::from(iv) * dv;
                if !surf.is_v_periodic() {
                    v = v.max(vmin).min(vmax);
                }
                let d = p.square_distance(&surf.value(u, v));
                if d < d0 {
                    return false;
                }
            }
        }
        true
    }

    /// `Perform(P, U0, V0, isDistanceCriteria = false)` (`cxx:121-225`).
    ///
    /// `isDistanceCriteria == true` runs the `math_BFGS` / `math_FRPR` branch
    /// over [`ExtremaFuncPsDist`] (`cxx:146-178`); `false` runs the
    /// `math_FunctionSetRoot` / `math_NewtonFunctionSetRoot` branch over
    /// [`ExtremaFuncPsNorm`] (`cxx:180-223`).
    pub fn perform(&mut self, p: &GpPnt, u0: f64, v0: f64, is_distance_criteria: bool) {
        self.done = false;

        // `math_Vector aTol(1, 2), aStart(1, 2), aBoundInf(1, 2),
        // aBoundSup(1, 2);` (`cxx:130-138`).
        let mut tol = MathVector::new(1, 2);
        tol.set_value(1, self.tol_u);
        tol.set_value(2, self.tol_v);
        let mut start = MathVector::new(1, 2);
        start.set_value(1, u0);
        start.set_value(2, v0);

        // `aBoundInf(1) = mySurf.FirstUParameter(); ... aBoundSup(2) =
        // mySurf.LastVParameter();` (`cxx:140-144`).
        let (umin, umax) = self.surf.u_range();
        let (vmin, vmax) = self.surf.v_range();
        let mut inf_bound = MathVector::new(1, 2);
        inf_bound.set_value(1, umin);
        inf_bound.set_value(2, vmin);
        let mut sup_bound = MathVector::new(1, 2);
        sup_bound.set_value(1, umax);
        sup_bound.set_value(2, vmax);

        if is_distance_criteria {
            // `Extrema_FuncPSDist F(mySurf, theP); math_BFGS aSolver(2,
            // aRelTol); aSolver.Perform(F, aStart);` (`cxx:148-156`).
            let a_rel_tol = 1.0e-8;
            let f = ExtremaFuncPsDist::new(self.surf, *p);

            let mut bfgs = BFGS::new();
            bfgs.tolerance = a_rel_tol;
            let bfgs_result = bfgs.minimize(|x| f.value(x), |x| f.gradient(x), &start);

            let location = match bfgs_result {
                Ok(location) => {
                    // `aSolver.Location(aResPnt); mySqDist =
                    // aSolver.Minimum();` (`cxx:169-172`).
                    self.sq_dist = bfgs.minimum();
                    location
                }
                Err(_) => {
                    // `math_FRPR aSolver1(F, aRelTol); aSolver1.Perform(F,
                    // aStart); if (!aSolver1.IsDone()) return;
                    // aSolver1.Location(aResPnt); mySqDist =
                    // aSolver1.Minimum();` (`cxx:158-168`).
                    let mut frpr = Frpr::new();
                    frpr.tolerance = a_rel_tol;
                    match frpr.perform(|x| f.value(x), |x| f.gradient(x), &start) {
                        Ok(result) => {
                            self.sq_dist = result.minimum;
                            result.location
                        }
                        Err(_) => return,
                    }
                }
            };

            // The two solvers are unconstrained in the `.cxx`; the objective
            // refuses values outside the box (see the module note), which the
            // port reproduces by rejecting the returned location.
            if !f.is_inside(&location) {
                return;
            }

            // `myPoint.SetParameters(aResPnt(1), aResPnt(2),
            // mySurf.Value(aResPnt(1), aResPnt(2))); myDone = true;`
            // (`cxx:175-177`).
            self.point = ExtremaPOnSurf::new(
                location.value(1),
                location.value(2),
                self.surf.value(location.value(1), location.value(2)),
            );
            self.done = true;
        } else {
            // `Extrema_FuncPSNorm F(theP, mySurf);` (`cxx:182`).
            let mut f = ExtremaFuncPsNorm::new(*p, self.surf);

            // `if (mySurf.GetType() == GeomAbs_BSplineSurface) { aTol(1) =
            // myTolU; aTol(2) = myTolV; CorrectTol(theU0, theV0, aTol); }`
            // (`cxx:184-189`). `GetType()` is the
            // `Surface::is_bspline_surface` query of the port.
            if self.surf.is_bspline_surface() {
                tol.set_value(1, self.tol_u);
                tol.set_value(2, self.tol_v);
                correct_tol(u0, v0, &mut tol);
            }

            // `bool isCorrectTol = (abs(aTol(1) - myTolU) > PConfusion ||
            // abs(aTol(2) - myTolV) > PConfusion);` (`cxx:191-192`).
            let is_correct_tol = (tol.value(1) - self.tol_u).abs() > PCONFUSION
                || (tol.value(2) - self.tol_v).abs() > PCONFUSION;

            // `math_FunctionSetRoot aSR(F, aTol); aSR.Perform(F, aStart,
            // aBoundInf, aBoundSup);` (`cxx:194-195`).
            let mut sr = MathFunctionSetRoot::new(&f, &tol, 100);
            sr.perform_with_bounds(&mut f, &start, &inf_bound, &sup_bound, false);

            if !sr.is_done() || is_correct_tol {
                // `if (isCorrectTol) { aTol(1) = myTolU; aTol(2) = myTolV; }
                // math_NewtonFunctionSetRoot aNSR(F, aTol,
                // Precision::Confusion()); aNSR.Perform(F, aStart,
                // aBoundInf, aBoundSup); if (!aSR.IsDone() &&
                // !aNSR.IsDone()) return;` (`cxx:197-209`).
                if is_correct_tol {
                    tol.set_value(1, self.tol_u);
                    tol.set_value(2, self.tol_v);
                }
                let mut nsr =
                    MathNewtonFunctionSetRoot::new(&mut f, &tol, CONFUSION, 100);
                nsr.perform_bounded(&mut f, &start, &inf_bound, &sup_bound);
                if !sr.is_done() && !nsr.is_done() {
                    return;
                }
            }

            // `double aNbExt = F.NbExt(); mySqDist = F.SquareDistance(1);
            // myPoint = F.Point(1); for (i = 2; i <= aNbExt; ++i) ...`
            // (`cxx:211-223`).
            let nb_ext = f.nb_ext();
            if nb_ext == 0 {
                // `F.SquareDistance(1)` raises Standard_OutOfRange in OCCT
                // when no station was stored; the port answers "not done".
                return;
            }
            self.sq_dist = f.square_distance(1);
            self.point = f.point(1);
            for i in 2..=nb_ext {
                if f.square_distance(i) < self.sq_dist {
                    self.sq_dist = f.square_distance(i);
                    self.point = f.point(i);
                }
            }
            self.done = true;
        }
    }

    /// `Perform(P, U0, V0)` with the `hxx:47-49` default,
    /// `isDistanceCriteria = false` — the normal projection criteria.
    pub fn perform_default(&mut self, p: &GpPnt, u0: f64, v0: f64) {
        self.perform(p, u0, v0, false);
    }

    /// `IsDone()` (`cxx:229-232`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `SquareDistance()` (`cxx:236-243`): `None` where OCCT raises
    /// `StdFail_NotDone`.
    pub fn square_distance(&self) -> Option<f64> {
        if self.done {
            Some(self.sq_dist)
        } else {
            None
        }
    }

    /// `Point()` (`cxx:247-253`): `None` where OCCT raises
    /// `StdFail_NotDone`.
    pub fn point(&self) -> Option<ExtremaPOnSurf> {
        if self.done {
            Some(self.point)
        } else {
            None
        }
    }
}
