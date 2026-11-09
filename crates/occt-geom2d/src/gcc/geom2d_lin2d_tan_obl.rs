//! `Geom2dGcc_Lin2dTanObl` (`Geom2dGcc_Lin2dTanObl.{hxx:39-166,cxx:31-287}`):
//! lines tangent to a circle/curve and making a given angle with a reference
//! line.
//!
//! The class is a front end: a circle operand goes to the analytic
//! [`GccAnaLin2dTanObl`]; every other operand is solved by sampling the
//! parameter range and calling [`Geom2dGccLin2dTanOblIter`] at each grid point,
//! keeping the results that [`Geom2dGccLin2dTanObl::add`] does not reject.
//!
//! | OCCT | arguments | `cxx` |
//! |---|---|---|
//! | `(QualifiedCurve, gp_Lin2d, TolAng, Angle)` | `:31-101` |
//! | `(QualifiedCurve, gp_Lin2d, TolAng, Param1, Angle)` | `:103-160` |
//! | `IsDone`..`Add` | `:162-287` |
//!
//! Translation notes:
//! - the circle arm (`cxx:54-72`, `:129-146`) reads `WellDone = Lin.IsDone()`
//!   and, only inside `if (WellDone)`, copies `NbrSol = Lin.NbSolutions()`
//!   solutions through `ThisSolution` / `Tangency1` / `Intersection2` /
//!   `WhichQualifier`,
//! - the sampling loop is `for (i = 0; i <= aNbSamples && NbrSol < 2; i++)`
//!   with `aStep = (aLastPar - aFirstPar) / aNbSamples` (`cxx:80-97`), i.e. 21
//!   `Lin2dTanOblIter` calls at most; `WellDone` is then `NbrSol > 0`,
//! - **OCCT defect kept**: the five-argument constructor calls
//!   `Geom2dGcc_Lin2dTanOblIter Lin(Qc1, TheLine, TolAng, Param1, Angle)`
//!   (`cxx:150`) while the iter signature is
//!   `(Qualified1, TheLin, Param1, TolAng, Angle)` (`Lin2dTanOblIter.cxx:38-41`),
//!   so the angular tolerance and the start parameter are swapped; the port
//!   passes `tol_ang` into the parameter slot and `param1` into the tolerance
//!   slot,
//! - **OCCT defect kept**: the same constructor never assigns `NbrSol`
//!   (`cxx:147-159` leaves it at the `:127` initialisation), so a successful
//!   curve arm reports `WellDone = true` with `NbrSol = 0`; `NbSolutions()` is
//!   `0` and `ThisSolution(1)` throws `Standard_OutOfRange`,
//! - `Add` (`cxx:240-287`) re-reads the candidate solution and its two
//!   tangency/intersection data sets, rejects a candidate whose `pararg1` /
//!   `pararg2` repeat within `theTol` (`cxx:258-263`) and one whose direction
//!   is not parallel to the curve tangent at `pararg1` within `theTol`
//!   (`cxx:270-273`); `gp_Dir2d(aVTan)` normalises the tangent, so a null
//!   tangent panics,
//! - `Add` calls `theLin.Intersection2` (`cxx:256`) unconditionally, so a
//!   parallel candidate propagates `Geom2dGcc_IsParallel` out of the
//!   constructor; kept,
//! - `ThisSolution` (`cxx:171-180`) and `NbSolutions` (`cxx:168-171`) check
//!   neither `WellDone` nor the solution count beyond `1 <= Index <= NbrSol`;
//!   `WhichQualifier` / `Tangency1` / `Intersection2` check both,
//! - `Paral2` (`hxx:149`) is initialised to `false` (`cxx:36`, `:108`) and never
//!   read again; `IsParallel2` is declared (`hxx:143`) but has no definition in
//!   `Geom2dGcc_Lin2dTanObl.cxx`, so there is no behaviour to port and the dead
//!   member is omitted rather than carried as an unread field.
//!
//! UNPORTED: none - both constructors, the six accessors and `Add`.

use occt_core::gp::{GpCirc2d, GpDir2d, GpLin2d, GpPnt2d};
use occt_core::kernel::geomabs::CurveType;

use crate::adaptor_curve::{basis2d, type2d, Geom2dAdaptorCurve};
use crate::gcc::curve_tool;
use crate::gcc::gcc_ent::{GccEntPosition, GccEntQualifiedCirc};
use crate::gcc::lin2d_tan_obl::{GccAnaLin2dTanObl, MAX_SOL};
use crate::gcc::lin2d_tan_obl_iter::Geom2dGccLin2dTanOblIter;
use crate::gcc::qualified_curve::{Geom2dGccQCurve, Geom2dGccQualifiedCurve};

/// `Geom2d_Circle` cast of the analytic arm (`cxx:57-58`, `:131-132`).
fn as_circle(q: &Geom2dGccQualifiedCurve) -> GpCirc2d {
    basis2d(q.qualified().as_ref())
        .gp_circ2d()
        .expect("Geom2d_Circle down_cast")
}

/// `GccEnt_QualifiedCirc(c1, Qualified1.Qualifier())` (`cxx:59`, `:133`).
fn as_qualified_circ(q: &Geom2dGccQualifiedCurve) -> GccEntQualifiedCirc {
    GccEntQualifiedCirc::new(as_circle(q), q.qualifier())
}

/// `Geom2dAdaptor_Curve C1 = Qualified1.Qualified()` (`cxx:46`, `:118`).
fn as_adaptor(q: &Geom2dGccQualifiedCurve) -> Geom2dAdaptorCurve {
    Geom2dAdaptorCurve::from_curve(q.qualified().clone())
}

/// `gp_Dir2d(Vect)` (`gp_Dir2d.cxx:47-54`).
fn dir_of(x: f64, y: f64) -> GpDir2d {
    GpDir2d::new(x, y).expect("gp_VectorWithNullMagnitude in Geom2dGcc_Lin2dTanObl")
}

/// `Geom2dGcc_Lin2dTanObl` (`hxx:39-166`).
#[derive(Clone)]
pub struct Geom2dGccLin2dTanObl {
    /// `WellDone` (`hxx:151`).
    well_done: bool,
    /// `NbrSol` (`hxx:152`).
    nbr_sol: usize,
    /// `linsol` (`hxx:153`).
    linsol: [GpLin2d; MAX_SOL],
    /// `qualifier1` (`hxx:154`).
    qualifier1: [GccEntPosition; MAX_SOL],
    /// `pnttg1sol` (`hxx:155`).
    pnttg1sol: [GpPnt2d; MAX_SOL],
    /// `pntint2sol` (`hxx:156`).
    pntint2sol: [GpPnt2d; MAX_SOL],
    /// `par1sol` (`hxx:157`).
    par1sol: [f64; MAX_SOL],
    /// `par2sol` (`hxx:158`).
    par2sol: [f64; MAX_SOL],
    /// `pararg1` (`hxx:159`).
    pararg1: [f64; MAX_SOL],
    /// `pararg2` (`hxx:160`).
    pararg2: [f64; MAX_SOL],
}

impl Geom2dGccLin2dTanObl {
    /// The `: linsol(1, 2), qualifier1(1, 2), ... :` initializer list
    /// (`cxx:35-44`, `:107-116`).
    fn empty() -> Self {
        Self {
            well_done: false,
            nbr_sol: 0,
            linsol: [GpLin2d::default(); MAX_SOL],
            qualifier1: [GccEntPosition::NoQualifier; MAX_SOL],
            pnttg1sol: [GpPnt2d::zero(); MAX_SOL],
            pntint2sol: [GpPnt2d::zero(); MAX_SOL],
            par1sol: [0.0; MAX_SOL],
            par2sol: [0.0; MAX_SOL],
            pararg1: [0.0; MAX_SOL],
            pararg2: [0.0; MAX_SOL],
        }
    }

    /// The circle arm shared by both constructors (`cxx:54-72`,
    /// `cxx:129-146`).
    fn from_analytic(&mut self, qualified1: &Geom2dGccQualifiedCurve, the_line: &GpLin2d, angle: f64) {
        let lin = GccAnaLin2dTanObl::tangent_to_circle(&as_qualified_circ(qualified1), the_line, angle);
        self.well_done = lin.is_done();
        if self.well_done {
            self.nbr_sol = lin.nb_solutions();
            for i in 1..=self.nbr_sol {
                let k = i - 1;
                self.linsol[k] = lin.this_solution(i);
                let (par_sol, par_arg, pnt) = lin.tangency1(i);
                self.par1sol[k] = par_sol;
                self.pararg1[k] = par_arg;
                self.pnttg1sol[k] = pnt;
                let (par_sol, par_arg, pnt) = lin.intersection2(i);
                self.par2sol[k] = par_sol;
                self.pararg2[k] = par_arg;
                self.pntint2sol[k] = pnt;
                self.qualifier1[k] = lin.which_qualifier(i);
            }
        }
    }

    /// `(QualifiedCurve, gp_Lin2d, TolAng, Angle)` (`cxx:31-101`).
    pub fn new(
        qualified1: &Geom2dGccQualifiedCurve,
        the_line: &GpLin2d,
        tol_ang: f64,
        angle: f64,
    ) -> Self {
        let mut s = Self::empty();
        let type1 = type2d(qualified1.qualified().as_ref());
        // `cxx:53-54`.
        s.well_done = false;
        s.nbr_sol = 0;
        if type1 == CurveType::Circle {
            s.from_analytic(qualified1, the_line, angle);
        } else {
            // `cxx:74-100`.
            let c1 = as_adaptor(qualified1);
            let qc1 = Geom2dGccQCurve::new(c1.clone(), qualified1.qualifier());
            let a_first_par = curve_tool::first_parameter(&c1);
            let a_last_par = curve_tool::last_parameter(&c1);
            let a_nb_samples = curve_tool::nb_samples(&c1);
            let a_step = (a_last_par - a_first_par) / f64::from(a_nb_samples);
            let mut param1 = a_first_par;
            let mut i = 0;
            while i <= a_nb_samples && s.nbr_sol < 2 {
                let lin = Geom2dGccLin2dTanOblIter::new(&qc1, the_line, param1, tol_ang, angle);
                if lin.is_done() && s.add(s.nbr_sol + 1, &lin, tol_ang, &c1) {
                    s.nbr_sol += 1;
                }
                param1 += a_step;
                i += 1;
            }
            s.well_done = s.nbr_sol > 0;
        }
        s
    }

    /// `(QualifiedCurve, gp_Lin2d, TolAng, Param1, Angle)` (`cxx:103-160`).
    pub fn new_from(
        qualified1: &Geom2dGccQualifiedCurve,
        the_line: &GpLin2d,
        tol_ang: f64,
        param1: f64,
        angle: f64,
    ) -> Self {
        let mut s = Self::empty();
        let type1 = type2d(qualified1.qualified().as_ref());
        // `cxx:126-127`.
        s.well_done = false;
        s.nbr_sol = 0;
        if type1 == CurveType::Circle {
            s.from_analytic(qualified1, the_line, angle);
        } else {
            // `cxx:148-159`.
            let c1 = as_adaptor(qualified1);
            let qc1 = Geom2dGccQCurve::new(c1, qualified1.qualifier());
            // `Geom2dGcc_Lin2dTanOblIter Lin(Qc1, TheLine, TolAng, Param1,
            // Angle)` (`cxx:150`): the OCCT argument order against the iter
            // signature `(Qualified1, TheLin, Param1, TolAng, Angle)`.
            let lin = Geom2dGccLin2dTanOblIter::new(&qc1, the_line, tol_ang, param1, angle);
            s.well_done = lin.is_done();
            if s.well_done {
                s.linsol[0] = lin.this_solution();
                let (par_sol, par_arg, pnt) = lin.tangency1();
                s.par1sol[0] = par_sol;
                s.pararg1[0] = par_arg;
                s.pnttg1sol[0] = pnt;
                let (par_sol, par_arg, pnt) = lin.intersection2();
                s.par2sol[0] = par_sol;
                s.pararg2[0] = par_arg;
                s.pntint2sol[0] = pnt;
                s.qualifier1[0] = lin.which_qualifier();
            }
        }
        s
    }

    /// `Add(theIndex, theLin, theTol, theC1)` (`cxx:240-287`).
    fn add(
        &mut self,
        the_index: usize,
        the_lin: &Geom2dGccLin2dTanOblIter,
        the_tol: f64,
        the_c1: &Geom2dAdaptorCurve,
    ) -> bool {
        let a_lin = the_lin.this_solution();
        let (a_par1sol, a_par1arg, a_pnt1sol) = the_lin.tangency1();
        let (a_par2sol, a_par2arg, a_pnt2sol) = the_lin.intersection2();

        // `cxx:258-263`.
        for i in 1..the_index {
            if (a_par1arg - self.pararg1[i - 1]).abs() <= the_tol
                && (a_par2arg - self.pararg2[i - 1]).abs() <= the_tol
            {
                return false;
            }
        }

        // `cxx:265-274`.
        let a_lin_dir = *a_lin.direction();
        let (_, a_v_tan) =
            curve_tool::d1(the_c1, a_par1arg).expect("Geom2dAdaptor_Curve::D1 in Add");
        if a_lin_dir.crossed(&dir_of(a_v_tan.x(), a_v_tan.y())).abs() > the_tol {
            return false;
        }

        let k = the_index - 1;
        self.linsol[k] = a_lin;
        self.par1sol[k] = a_par1sol;
        self.pararg1[k] = a_par1arg;
        self.pnttg1sol[k] = a_pnt1sol;
        self.par2sol[k] = a_par2sol;
        self.pararg2[k] = a_par2arg;
        self.pntint2sol[k] = a_pnt2sol;
        self.qualifier1[k] = the_lin.which_qualifier();
        true
    }

    /// `IsDone()` (`cxx:162-165`).
    pub fn is_done(&self) -> bool {
        self.well_done
    }

    /// `NbSolutions()` (`cxx:167-171`): no `WellDone` check.
    pub fn nb_solutions(&self) -> usize {
        self.nbr_sol
    }

    /// `ThisSolution(Index)` (`cxx:173-180`): only the index is checked.
    pub fn this_solution(&self, index: usize) -> GpLin2d {
        assert!(
            index >= 1 && index <= self.nbr_sol,
            "Standard_OutOfRange in Geom2dGcc_Lin2dTanObl"
        );
        self.linsol[index - 1]
    }

    /// `WhichQualifier(Index, Qualif1)` (`cxx:182-195`).
    pub fn which_qualifier(&self, index: usize) -> GccEntPosition {
        self.check(index);
        self.qualifier1[index - 1]
    }

    /// `Tangency1(Index, ParSol, ParArg, PntSol)` (`cxx:197-215`).
    pub fn tangency1(&self, index: usize) -> (f64, f64, GpPnt2d) {
        self.check(index);
        (
            self.par1sol[index - 1],
            self.pararg1[index - 1],
            self.pnttg1sol[index - 1],
        )
    }

    /// `Intersection2(Index, ParSol, ParArg, PntSol)` (`cxx:217-235`).
    pub fn intersection2(&self, index: usize) -> (f64, f64, GpPnt2d) {
        self.check(index);
        (
            self.par2sol[index - 1],
            self.pararg2[index - 1],
            self.pntint2sol[index - 1],
        )
    }

    /// `if (!WellDone) throw StdFail_NotDone()` and
    /// `if (Index <= 0 || Index > NbrSol) throw Standard_OutOfRange()`.
    fn check(&self, index: usize) {
        assert!(self.well_done, "StdFail_NotDone in Geom2dGcc_Lin2dTanObl");
        assert!(
            index >= 1 && index <= self.nbr_sol,
            "Standard_OutOfRange in Geom2dGcc_Lin2dTanObl"
        );
    }
}
