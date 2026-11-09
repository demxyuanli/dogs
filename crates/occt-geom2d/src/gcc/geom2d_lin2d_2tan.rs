//! `Geom2dGcc_Lin2d2Tan` (`Geom2dGcc_Lin2d2Tan.{hxx:48-121,cxx:31-465}`): lines
//! tangent to two circles/curves or to a circle/curve and a point.
//!
//! The class is a front end: a circle/circle or circle/point pair goes to the
//! analytic [`GccAnaLin2d2Tan`]; every other pair is solved by sampling the two
//! parameter ranges and calling [`Geom2dGccLin2d2TanIter`] at each grid point,
//! keeping the results that [`Geom2dGccLin2d2Tan::add`] does not reject.
//!
//! | OCCT | arguments | `cxx` |
//! |---|---|---|
//! | `(QualifiedCurve, QualifiedCurve, Tolang)` | `:31-120` |
//! | `(QualifiedCurve, gp_Pnt2d, Tolang)` | `:122-191` |
//! | `(QualifiedCurve, QualifiedCurve, Tolang, Param1, Param2)` | `:195-262` |
//! | `(QualifiedCurve, gp_Pnt2d, Tolang, Param1)` | `:264-326` |
//! | `IsDone`..`Add` | `:328-465` |
//!
//! Translation notes:
//! - `Geom2dAdaptor_Curve C1 = Qualified1.Qualified()` (`cxx:44-45`) is rebuilt
//!   from the stored `Arc<dyn Curve2d>` with
//!   [`Geom2dAdaptorCurve::from_curve`], so `FirstParameter` / `LastParameter`
//!   / `NbSamples` / `EpsX` see the same curve OCCT wraps,
//! - the sampling loops are `for (i = 0; i <= aNbSamples1 && NbrSol < 4; i++)`
//!   / `for (j = 0; j <= aNbSamples2 && NbrSol < 4; j++)` (`cxx:98-115`) and the
//!   point variants stop at `< 2` (`cxx:175`); `Param1` advances once per outer
//!   iteration, `Param2` restarts at `a2FPar` and advances once per inner
//!   iteration; kept,
//! - `aStep = (aLastPar - aFirstPar) / aNbSamples` divides by
//!   `CurveTool::NbSamples` = 20 (`cxx:87`, `:170`); a zero would be a division
//!   by zero in both OCCT and the port,
//! - `Add` (`cxx:407-465`) rejects a candidate when `|aPar1arg - pararg1(i)| <=
//!   theTol && |aPar2arg - pararg2(i)| <= theTol` for any earlier solution, and
//!   when `|aLinDir x gp_Dir2d(aVTan)| > theTol` at either tangency point; the
//!   second check is skipped for the null second adaptor of the point overload
//!   (`cxx:443`),
//! - the `Modified by Sergey KHROMOV` guards (`cxx:69-70`, `:233-234`,
//!   `:295-296`) only add the `if (WellDone)` wrapper around the copy loops;
//!   kept,
//! - the three-argument constructors only enter the sampling branch when the
//!   operand is not a circle, so a curve/circle pair also samples the circle's
//!   range: the circle arrives as a full `Geom2d_Circle` adaptor, not as a
//!   `GpCirc2d`,
//! - `StdFail_NotDone` / `Standard_OutOfRange` are panics; `ThisSolution`
//!   checks only the index (`cxx:337-343`), `WhichQualifier` / `Tangency1` /
//!   `Tangency2` also check `WellDone`.
//!
//! UNPORTED: none - the four constructors, the six accessors and `Add`.

use occt_core::gp::{GpCirc2d, GpDir2d, GpLin2d, GpPnt2d};

use crate::adaptor_curve::{basis2d, type2d, Geom2dAdaptorCurve};
use crate::gcc::curve_tool;
use crate::gcc::gcc_ent::{GccEntPosition, GccEntQualifiedCirc};
use crate::gcc::lin2d_2tan::{GccAnaLin2d2Tan, MAX_SOL};
use crate::gcc::lin2d_2tan_iter::Geom2dGccLin2d2TanIter;
use crate::gcc::qualified_curve::{Geom2dGccQCurve, Geom2dGccQualifiedCurve};
use occt_core::kernel::geomabs::CurveType;

/// `Geom2dGcc_Lin2d2Tan` (`hxx:48-121`).
#[derive(Clone)]
pub struct Geom2dGccLin2d2Tan {
    /// `WellDone` (`hxx:109`).
    well_done: bool,
    /// `NbrSol` (`hxx:110`).
    nbr_sol: usize,
    /// `linsol` (`hxx:111`).
    linsol: [GpLin2d; MAX_SOL],
    /// `qualifier1` (`hxx:112`).
    qualifier1: [GccEntPosition; MAX_SOL],
    /// `qualifier2` (`hxx:113`).
    qualifier2: [GccEntPosition; MAX_SOL],
    /// `pnttg1sol` (`hxx:114`).
    pnttg1sol: [GpPnt2d; MAX_SOL],
    /// `pnttg2sol` (`hxx:115`).
    pnttg2sol: [GpPnt2d; MAX_SOL],
    /// `par1sol` (`hxx:116`).
    par1sol: [f64; MAX_SOL],
    /// `par2sol` (`hxx:117`).
    par2sol: [f64; MAX_SOL],
    /// `pararg1` (`hxx:118`).
    pararg1: [f64; MAX_SOL],
    /// `pararg2` (`hxx:119`).
    pararg2: [f64; MAX_SOL],
}

/// `Geom2d_Circle` / `Geom2d_Circle` casts of the analytic arm (`cxx:59-62`).
fn as_circle(q: &Geom2dGccQualifiedCurve) -> GpCirc2d {
    basis2d(q.qualified().as_ref())
        .gp_circ2d()
        .expect("Geom2d_Circle down_cast")
}

/// `GccEnt_QualifiedCirc(c1, Qualified1.Qualifier())`.
fn as_qualified_circ(q: &Geom2dGccQualifiedCurve) -> GccEntQualifiedCirc {
    GccEntQualifiedCirc::new(as_circle(q), q.qualifier())
}

/// `Geom2dAdaptor_Curve C1 = Qualified1.Qualified()` (`cxx:44`).
fn as_adaptor(q: &Geom2dGccQualifiedCurve) -> Geom2dAdaptorCurve {
    Geom2dAdaptorCurve::from_curve(q.qualified().clone())
}

/// `gp_Dir2d(Vect)` (`gp_Dir2d.cxx:47-54`).
fn dir_of(x: f64, y: f64) -> GpDir2d {
    GpDir2d::new(x, y).expect("gp_VectorWithNullMagnitude in Geom2dGcc_Lin2d2Tan")
}

impl Geom2dGccLin2d2Tan {
    /// The `: linsol(1, 4), qualifier1(1, 4), ... :` initializer list
    /// (`cxx:36-43`, `:126-133`, `:201-208`, `:268-275`).
    fn empty() -> Self {
        Self {
            well_done: false,
            nbr_sol: 0,
            linsol: [GpLin2d::default(); MAX_SOL],
            qualifier1: [GccEntPosition::NoQualifier; MAX_SOL],
            qualifier2: [GccEntPosition::NoQualifier; MAX_SOL],
            pnttg1sol: [GpPnt2d::zero(); MAX_SOL],
            pnttg2sol: [GpPnt2d::zero(); MAX_SOL],
            par1sol: [0.0; MAX_SOL],
            par2sol: [0.0; MAX_SOL],
            pararg1: [0.0; MAX_SOL],
            pararg2: [0.0; MAX_SOL],
        }
    }

    /// `(QualifiedCurve, QualifiedCurve, Tolang)` (`cxx:31-120`).
    pub fn new_two_curves(
        qualified1: &Geom2dGccQualifiedCurve,
        qualified2: &Geom2dGccQualifiedCurve,
        tolang: f64,
    ) -> Self {
        let mut s = Self::empty();
        let type1 = type2d(qualified1.qualified().as_ref());
        let type2 = type2d(qualified2.qualified().as_ref());
        s.nbr_sol = 0;
        if type1 == CurveType::Circle && type2 == CurveType::Circle {
            let lin = GccAnaLin2d2Tan::from_two_circles(
                &as_qualified_circ(qualified1),
                &as_qualified_circ(qualified2),
                tolang,
            );
            s.well_done = lin.is_done();
            if s.well_done {
                s.nbr_sol = lin.nb_solutions();
                for i in 1..=s.nbr_sol {
                    let k = i - 1;
                    s.linsol[k] = lin.this_solution(i);
                    let (par_sol, par_arg, pnt) = lin.tangency1(i);
                    s.par1sol[k] = par_sol;
                    s.pararg1[k] = par_arg;
                    s.pnttg1sol[k] = pnt;
                    let (par_sol, par_arg, pnt) = lin.tangency2(i);
                    s.par2sol[k] = par_sol;
                    s.pararg2[k] = par_arg;
                    s.pnttg2sol[k] = pnt;
                    let (q1, q2) = lin.which_qualifier(i);
                    s.qualifier1[k] = q1;
                    s.qualifier2[k] = q2;
                }
            }
        } else {
            let c1 = as_adaptor(qualified1);
            let c2 = as_adaptor(qualified2);
            let qc1 = Geom2dGccQCurve::new(c1.clone(), qualified1.qualifier());
            let a1f_par = curve_tool::first_parameter(&c1);
            let a1l_par = curve_tool::last_parameter(&c1);
            let a_nb_samples1 = curve_tool::nb_samples(&c1);
            let a_step1 = (a1l_par - a1f_par) / f64::from(a_nb_samples1);
            let mut param1 = a1f_par;
            let qc2 = Geom2dGccQCurve::new(c2.clone(), qualified2.qualifier());
            let a2f_par = curve_tool::first_parameter(&c2);
            let a2l_par = curve_tool::last_parameter(&c2);
            let a_nb_samples2 = curve_tool::nb_samples(&c2);
            let a_step2 = (a2l_par - a2f_par) / f64::from(a_nb_samples2);
            let mut i = 0i32;
            while i <= a_nb_samples1 && s.nbr_sol < MAX_SOL {
                let mut param2 = a2f_par;
                let mut j = 0i32;
                while j <= a_nb_samples2 && s.nbr_sol < MAX_SOL {
                    let lin =
                        Geom2dGccLin2d2TanIter::new_curve_curve(&qc1, &qc2, param1, param2, tolang);
                    if lin.is_done()
                        && s.add(s.nbr_sol + 1, &lin, tolang, &c1, Some(&c2))
                    {
                        s.nbr_sol += 1;
                    }
                    param2 += a_step2;
                    j += 1;
                }
                param1 += a_step1;
                i += 1;
            }
            s.well_done = s.nbr_sol > 0;
        }
        s
    }

    /// `(QualifiedCurve, ThePoint, Tolang)` (`cxx:122-191`).
    pub fn new_curve_point(
        qualified1: &Geom2dGccQualifiedCurve,
        the_point: &GpPnt2d,
        tolang: f64,
    ) -> Self {
        let mut s = Self::empty();
        let type1 = type2d(qualified1.qualified().as_ref());
        s.nbr_sol = 0;
        if type1 == CurveType::Circle {
            let lin = GccAnaLin2d2Tan::from_circle_and_point(
                &as_qualified_circ(qualified1),
                the_point,
                tolang,
            );
            s.well_done = lin.is_done();
            if s.well_done {
                s.nbr_sol = lin.nb_solutions();
                for i in 1..=s.nbr_sol {
                    let k = i - 1;
                    s.linsol[k] = lin.this_solution(i);
                    let (par_sol, par_arg, pnt) = lin.tangency1(i);
                    s.par1sol[k] = par_sol;
                    s.pararg1[k] = par_arg;
                    s.pnttg1sol[k] = pnt;
                    let (par_sol, par_arg, pnt) = lin.tangency2(i);
                    s.par2sol[k] = par_sol;
                    s.pararg2[k] = par_arg;
                    s.pnttg2sol[k] = pnt;
                    let (q1, q2) = lin.which_qualifier(i);
                    s.qualifier1[k] = q1;
                    s.qualifier2[k] = q2;
                }
            }
        } else {
            let c1 = as_adaptor(qualified1);
            let qc1 = Geom2dGccQCurve::new(c1.clone(), qualified1.qualifier());
            let a_first_par = curve_tool::first_parameter(&c1);
            let a_last_par = curve_tool::last_parameter(&c1);
            let a_nb_samples = curve_tool::nb_samples(&c1);
            let a_step = (a_last_par - a_first_par) / f64::from(a_nb_samples);
            let mut param1 = a_first_par;
            let mut i = 0i32;
            while i <= a_nb_samples && s.nbr_sol < 2 {
                let lin = Geom2dGccLin2d2TanIter::new_curve_point(&qc1, the_point, param1, tolang);
                if lin.is_done() && s.add(s.nbr_sol + 1, &lin, tolang, &c1, None) {
                    s.nbr_sol += 1;
                }
                param1 += a_step;
                i += 1;
            }
            s.well_done = s.nbr_sol > 0;
        }
        s
    }

    /// `(QualifiedCurve, QualifiedCurve, Tolang, Param1, Param2)`
    /// (`cxx:195-262`).
    pub fn new_two_curves_from(
        qualified1: &Geom2dGccQualifiedCurve,
        qualified2: &Geom2dGccQualifiedCurve,
        tolang: f64,
        param1: f64,
        param2: f64,
    ) -> Self {
        let mut s = Self::empty();
        let type1 = type2d(qualified1.qualified().as_ref());
        let type2 = type2d(qualified2.qualified().as_ref());
        s.nbr_sol = 0;
        if type1 == CurveType::Circle && type2 == CurveType::Circle {
            let lin = GccAnaLin2d2Tan::from_two_circles(
                &as_qualified_circ(qualified1),
                &as_qualified_circ(qualified2),
                tolang,
            );
            s.well_done = lin.is_done();
            if s.well_done {
                s.nbr_sol = lin.nb_solutions();
                for i in 1..=s.nbr_sol {
                    let k = i - 1;
                    s.linsol[k] = lin.this_solution(i);
                    let (par_sol, par_arg, pnt) = lin.tangency1(i);
                    s.par1sol[k] = par_sol;
                    s.pararg1[k] = par_arg;
                    s.pnttg1sol[k] = pnt;
                    let (par_sol, par_arg, pnt) = lin.tangency2(i);
                    s.par2sol[k] = par_sol;
                    s.pararg2[k] = par_arg;
                    s.pnttg2sol[k] = pnt;
                    let (q1, q2) = lin.which_qualifier(i);
                    s.qualifier1[k] = q1;
                    s.qualifier2[k] = q2;
                }
            }
        } else {
            let qc1 = Geom2dGccQCurve::new(as_adaptor(qualified1), qualified1.qualifier());
            let qc2 = Geom2dGccQCurve::new(as_adaptor(qualified2), qualified2.qualifier());
            let lin =
                Geom2dGccLin2d2TanIter::new_curve_curve(&qc1, &qc2, param1, param2, tolang);
            s.well_done = lin.is_done();
            if s.well_done {
                s.nbr_sol = 1;
                s.linsol[0] = lin.this_solution();
                let (par_sol, par_arg, pnt) = lin.tangency1();
                s.par1sol[0] = par_sol;
                s.pararg1[0] = par_arg;
                s.pnttg1sol[0] = pnt;
                let (par_sol, par_arg, pnt) = lin.tangency2();
                s.par2sol[0] = par_sol;
                s.pararg2[0] = par_arg;
                s.pnttg2sol[0] = pnt;
                let (q1, q2) = lin.which_qualifier();
                s.qualifier1[0] = q1;
                s.qualifier2[0] = q2;
            }
        }
        s
    }

    /// `(QualifiedCurve, ThePoint, Tolang, Param1)` (`cxx:264-326`).
    pub fn new_curve_point_from(
        qualified1: &Geom2dGccQualifiedCurve,
        the_point: &GpPnt2d,
        tolang: f64,
        param1: f64,
    ) -> Self {
        let mut s = Self::empty();
        let type1 = type2d(qualified1.qualified().as_ref());
        s.nbr_sol = 0;
        if type1 == CurveType::Circle {
            let lin = GccAnaLin2d2Tan::from_circle_and_point(
                &as_qualified_circ(qualified1),
                the_point,
                tolang,
            );
            s.well_done = lin.is_done();
            if s.well_done {
                s.nbr_sol = lin.nb_solutions();
                for i in 1..=s.nbr_sol {
                    let k = i - 1;
                    s.linsol[k] = lin.this_solution(i);
                    let (par_sol, par_arg, pnt) = lin.tangency1(i);
                    s.par1sol[k] = par_sol;
                    s.pararg1[k] = par_arg;
                    s.pnttg1sol[k] = pnt;
                    let (par_sol, par_arg, pnt) = lin.tangency2(i);
                    s.par2sol[k] = par_sol;
                    s.pararg2[k] = par_arg;
                    s.pnttg2sol[k] = pnt;
                    let (q1, q2) = lin.which_qualifier(i);
                    s.qualifier1[k] = q1;
                    s.qualifier2[k] = q2;
                }
            }
        } else {
            let qc1 = Geom2dGccQCurve::new(as_adaptor(qualified1), qualified1.qualifier());
            let lin = Geom2dGccLin2d2TanIter::new_curve_point(&qc1, the_point, param1, tolang);
            s.well_done = lin.is_done();
            if s.well_done {
                s.nbr_sol = 1;
                s.linsol[0] = lin.this_solution();
                let (par_sol, par_arg, pnt) = lin.tangency1();
                s.par1sol[0] = par_sol;
                s.pararg1[0] = par_arg;
                s.pnttg1sol[0] = pnt;
                let (par_sol, par_arg, pnt) = lin.tangency2();
                s.par2sol[0] = par_sol;
                s.pararg2[0] = par_arg;
                s.pnttg2sol[0] = pnt;
                let (q1, q2) = lin.which_qualifier();
                s.qualifier1[0] = q1;
                s.qualifier2[0] = q2;
            }
        }
        s
    }

    /// `Add(theIndex, theLin, theTol, theC1, theC2)` (`cxx:407-465`).
    fn add(
        &mut self,
        the_index: usize,
        the_lin: &Geom2dGccLin2d2TanIter,
        the_tol: f64,
        the_c1: &Geom2dAdaptorCurve,
        the_c2: Option<&Geom2dAdaptorCurve>,
    ) -> bool {
        let a_lin = the_lin.this_solution();
        let (a_par1sol, a_par1arg, a_pnt1sol) = the_lin.tangency1();
        let (a_par2sol, a_par2arg, a_pnt2sol) = the_lin.tangency2();

        for i in 1..the_index {
            if (a_par1arg - self.pararg1[i - 1]).abs() <= the_tol
                && (a_par2arg - self.pararg2[i - 1]).abs() <= the_tol
            {
                return false;
            }
        }

        let a_lin_dir = *a_lin.direction();
        let (_, a_v_tan) = curve_tool::d1(the_c1, a_par1arg)
            .expect("Geom2dAdaptor_Curve::D1 in Geom2dGcc_Lin2d2Tan::Add");
        if a_lin_dir.crossed(&dir_of(a_v_tan.x(), a_v_tan.y())).abs() > the_tol {
            return false;
        }

        if let Some(c2) = the_c2 {
            if c2.curve().is_some() {
                let (_, a_v_tan) = curve_tool::d1(c2, a_par2arg)
                    .expect("Geom2dAdaptor_Curve::D1 in Geom2dGcc_Lin2d2Tan::Add");
                if a_lin_dir.crossed(&dir_of(a_v_tan.x(), a_v_tan.y())).abs() > the_tol {
                    return false;
                }
            }
        }

        let k = the_index - 1;
        self.linsol[k] = a_lin;
        self.par1sol[k] = a_par1sol;
        self.pararg1[k] = a_par1arg;
        self.pnttg1sol[k] = a_pnt1sol;
        self.par2sol[k] = a_par2sol;
        self.pararg2[k] = a_par2arg;
        self.pnttg2sol[k] = a_pnt2sol;

        let (q1, q2) = the_lin.which_qualifier();
        self.qualifier1[k] = q1;
        self.qualifier2[k] = q2;

        true
    }

    /// `IsDone()` (`cxx:328-331`).
    pub fn is_done(&self) -> bool {
        self.well_done
    }

    /// `NbSolutions()` (`cxx:333-336`).
    pub fn nb_solutions(&self) -> usize {
        self.nbr_sol
    }

    /// `ThisSolution(Index)` (`cxx:338-344`): only the index is checked.
    pub fn this_solution(&self, index: usize) -> GpLin2d {
        assert!(
            index >= 1 && index <= self.nbr_sol,
            "Standard_OutOfRange in Geom2dGcc_Lin2d2Tan"
        );
        self.linsol[index - 1]
    }

    /// `WhichQualifier(Index, Qualif1, Qualif2)` (`cxx:346-364`).
    pub fn which_qualifier(&self, index: usize) -> (GccEntPosition, GccEntPosition) {
        self.check(index);
        (self.qualifier1[index - 1], self.qualifier2[index - 1])
    }

    /// `Tangency1(Index, ParSol, ParArg, PntSol)` (`cxx:366-384`).
    pub fn tangency1(&self, index: usize) -> (f64, f64, GpPnt2d) {
        self.check(index);
        (
            self.par1sol[index - 1],
            self.pararg1[index - 1],
            self.pnttg1sol[index - 1],
        )
    }

    /// `Tangency2(Index, ParSol, ParArg, PntSol)` (`cxx:386-404`).
    pub fn tangency2(&self, index: usize) -> (f64, f64, GpPnt2d) {
        self.check(index);
        (
            self.par2sol[index - 1],
            self.pararg2[index - 1],
            self.pnttg2sol[index - 1],
        )
    }

    /// `if (!WellDone) throw StdFail_NotDone()` and
    /// `if (Index <= 0 || Index > NbrSol) throw Standard_OutOfRange()`.
    fn check(&self, index: usize) {
        assert!(self.well_done, "StdFail_NotDone in Geom2dGcc_Lin2d2Tan");
        assert!(
            index >= 1 && index <= self.nbr_sol,
            "Standard_OutOfRange in Geom2dGcc_Lin2d2Tan"
        );
    }
}
