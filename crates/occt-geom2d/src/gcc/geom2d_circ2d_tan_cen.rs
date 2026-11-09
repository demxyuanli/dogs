//! `Geom2dGcc_Circ2dTanCen` (`Geom2dGcc_Circ2dTanCen.{hxx:51-149,cxx:32-203}`):
//! circles centred on a point and tangent to a line, a circle or a general
//! curve.
//!
//! The class is a front end: a line or a circle operand goes to the analytic
//! [`GccAnaCirc2dTanCen`], anything else to the geometric
//! [`Geom2dGccCirc2dTanCenGeo`]; every solution is copied into the `cirsol` /
//! `qualifier1` / ... arrays of `hxx:140-149`.
//!
//! | OCCT | `cxx` |
//! |---|---|
//! | `(QualifiedCurve, Geom2d_Point, Tolerance)` | `:32-119` |
//! | `IsDone`..`IsTheSame1` | `:121-203` |
//!
//! Translation notes:
//! - the circle arm passes `Tolerance` (`cxx:59`), the line arm does not
//!   (`cxx:81`, `GccAna_Circ2dTanCen(Linetan, Pcenter)`), and the geometric arm
//!   passes it again (`cxx:108`); kept,
//! - both analytic arms call `Tangency1` unconditionally (`cxx:73`, `:95`),
//!   which throws `StdFail_NotDone` when `IsTheSame1(j)` is set (the centre
//!   coincides with the circle centre, `GccAna_Circ2dTanCen.cxx:180-188`); the
//!   port keeps the panic through `expect`,
//! - the geometric arm writes `TheSame1(j) = 0` and calls `Tangency1` **before**
//!   `WhichQualifier` (`cxx:112-116`), the analytic arms the other way round
//!   (`cxx:64-73`, `:86-95`); kept,
//! - `Geom2dGcc_QualifiedCurve` only stores the curve in this port, so
//!   `Geom2dAdaptor_Curve C1 = Qualified1.Qualified()` (`cxx:42`) becomes
//!   [`Geom2dAdaptorCurve::from_curve`] on that curve, exactly like the
//!   `Geom2dAdaptor_Curve(const Handle(Geom2d_Curve)&)` constructor OCCT runs
//!   when the front end takes a `QualifiedCurve`,
//! - `StdFail_NotDone` / `Standard_OutOfRange` are panics, matching the rest of
//!   the `gcc` module.
//!
//! UNPORTED: none - the constructor and all five accessors.

use occt_core::gp::{GpCirc2d, GpLin2d, GpPnt2d};
use occt_core::kernel::geomabs::CurveType;

use crate::adaptor_curve::{basis2d, type2d, Geom2dAdaptorCurve};
use crate::gcc::circ2d_tan_cen::{GccAnaCirc2dTanCen, MAX_SOL};
use crate::gcc::circ2d_tan_cen_geo::Geom2dGccCirc2dTanCenGeo;
use crate::gcc::gcc_ent::{GccEntPosition, GccEntQualifiedCirc};
use crate::gcc::qualified_curve::{Geom2dGccQCurve, Geom2dGccQualifiedCurve};

/// `Geom2dGcc_Circ2dTanCen` (`hxx:51-149`).
#[derive(Clone)]
pub struct Geom2dGccCirc2dTanCen {
    /// `WellDone` (`hxx:141`).
    well_done: bool,
    /// `NbrSol` (`hxx:142`).
    nbr_sol: usize,
    /// `cirsol` (`hxx:143`).
    cirsol: [GpCirc2d; MAX_SOL],
    /// `qualifier1` (`hxx:144`).
    qualifier1: [GccEntPosition; MAX_SOL],
    /// `TheSame1` (`hxx:145`).
    the_same1: [i32; MAX_SOL],
    /// `pnttg1sol` (`hxx:146`).
    pnttg1sol: [GpPnt2d; MAX_SOL],
    /// `par1sol` (`hxx:147`).
    par1sol: [f64; MAX_SOL],
    /// `pararg1` (`hxx:148`).
    pararg1: [f64; MAX_SOL],
}

/// `Geom2d_Circle` / `Geom2d_Line` casts of the front end (`cxx:57`, `:79`).
fn as_circle(q: &Geom2dGccQualifiedCurve) -> GpCirc2d {
    basis2d(q.qualified().as_ref())
        .gp_circ2d()
        .expect("Geom2d_Circle down_cast")
}

fn as_line(q: &Geom2dGccQualifiedCurve) -> GpLin2d {
    basis2d(q.qualified().as_ref())
        .gp_lin2d()
        .expect("Geom2d_Line down_cast")
}

fn as_qualified_circ(q: &Geom2dGccQualifiedCurve) -> GccEntQualifiedCirc {
    GccEntQualifiedCirc::new(as_circle(q), q.qualifier())
}

/// `Geom2dGcc_QCurve Qc1(C1, Qualified1.Qualifier())` (`cxx:107`).
fn as_qcurve(q: &Geom2dGccQualifiedCurve) -> Geom2dGccQCurve {
    Geom2dGccQCurve::new(
        Geom2dAdaptorCurve::from_curve(q.qualified().clone()),
        q.qualifier(),
    )
}

/// `C1.GetType()` (`cxx:44`).
fn curve_type(q: &Geom2dGccQualifiedCurve) -> CurveType {
    type2d(q.qualified().as_ref())
}

/// `TheSame1(j) = 1` / `TheSame1(j) = 0` (`cxx:66-73`) on an `int` array.
fn flag(value: bool) -> i32 {
    if value {
        1
    } else {
        0
    }
}

impl Geom2dGccCirc2dTanCen {
    /// The `: cirsol(1, 2), qualifier1(1, 2), ... :` initializer list
    /// (`cxx:36-40`); `WellDone` and `NbrSol` are read before being written by
    /// the constructor.
    fn empty() -> Self {
        Self {
            well_done: false,
            nbr_sol: 0,
            cirsol: [GpCirc2d::default(); MAX_SOL],
            qualifier1: [GccEntPosition::NoQualifier; MAX_SOL],
            the_same1: [0; MAX_SOL],
            pnttg1sol: [GpPnt2d::zero(); MAX_SOL],
            par1sol: [0.0; MAX_SOL],
            pararg1: [0.0; MAX_SOL],
        }
    }

    /// `Geom2dGcc_Circ2dTanCen(Qualified1, PCenter, Tolerance)` (`cxx:32-119`).
    pub fn new(
        qualified1: &Geom2dGccQualifiedCurve,
        pcenter: &GpPnt2d,
        tolerance: f64,
    ) -> Self {
        let mut s = Self::empty();
        let type1 = curve_type(qualified1);
        s.nbr_sol = 0;
        if type1 == CurveType::Line || type1 == CurveType::Circle {
            if type1 == CurveType::Circle {
                // `cxx:56-76`: the circle arm.
                let circ = GccAnaCirc2dTanCen::tangent_to_circle(
                    &as_qualified_circ(qualified1),
                    pcenter,
                    tolerance,
                );
                s.well_done = circ.is_done();
                s.nbr_sol = circ.nb_solutions();
                for j in 1..=s.nbr_sol {
                    let k = j - 1;
                    s.cirsol[k] = circ.this_solution(j);
                    s.qualifier1[k] = circ.which_qualifier(j);
                    s.the_same1[k] = flag(circ.is_the_same1(j));
                    let (par_sol, par_arg, pnt) = circ
                        .tangency1(j)
                        .expect("StdFail_NotDone in Geom2dGcc_Circ2dTanCen");
                    s.par1sol[k] = par_sol;
                    s.pararg1[k] = par_arg;
                    s.pnttg1sol[k] = pnt;
                }
            } else {
                // `cxx:78-97`: the line arm.
                let circ = GccAnaCirc2dTanCen::tangent_to_line(&as_line(qualified1), pcenter);
                s.well_done = circ.is_done();
                s.nbr_sol = circ.nb_solutions();
                for j in 1..=s.nbr_sol {
                    let k = j - 1;
                    s.cirsol[k] = circ.this_solution(j);
                    s.qualifier1[k] = circ.which_qualifier(j);
                    s.the_same1[k] = flag(circ.is_the_same1(j));
                    let (par_sol, par_arg, pnt) = circ
                        .tangency1(j)
                        .expect("StdFail_NotDone in Geom2dGcc_Circ2dTanCen");
                    s.par1sol[k] = par_sol;
                    s.pararg1[k] = par_arg;
                    s.pnttg1sol[k] = pnt;
                }
            }
        } else {
            // `cxx:105-117`: the geometric arm.
            let geo = Geom2dGccCirc2dTanCenGeo::new(&as_qcurve(qualified1), pcenter, tolerance);
            s.well_done = geo.is_done();
            s.nbr_sol = geo.nb_solutions();
            for j in 1..=s.nbr_sol {
                let k = j - 1;
                s.cirsol[k] = geo.this_solution(j);
                s.the_same1[k] = 0;
                let (par_sol, par_arg, pnt) = geo.tangency1(j);
                s.par1sol[k] = par_sol;
                s.pararg1[k] = par_arg;
                s.pnttg1sol[k] = pnt;
                s.qualifier1[k] = geo.which_qualifier(j);
            }
        }
        s
    }

    /// `IsDone()` (`cxx:121-124`).
    pub fn is_done(&self) -> bool {
        self.well_done
    }

    /// `NbSolutions()` (`cxx:126-129`).
    pub fn nb_solutions(&self) -> usize {
        self.nbr_sol
    }

    /// `ThisSolution(Index)` (`cxx:131-142`).
    pub fn this_solution(&self, index: usize) -> GpCirc2d {
        self.check(index);
        self.cirsol[index - 1]
    }

    /// `WhichQualifier(Index, Qualif1)` (`cxx:144-158`).
    pub fn which_qualifier(&self, index: usize) -> GccEntPosition {
        self.check(index);
        self.qualifier1[index - 1]
    }

    /// `Tangency1(Index, ParSol, ParArg, PntSol)` (`cxx:160-186`).
    ///
    /// `None` where OCCT raises `StdFail_NotDone` for `TheSame1(Index) != 0`.
    pub fn tangency1(&self, index: usize) -> Option<(f64, f64, GpPnt2d)> {
        self.check(index);
        if self.the_same1[index - 1] != 0 {
            return None;
        }
        Some((
            self.par1sol[index - 1],
            self.pararg1[index - 1],
            self.pnttg1sol[index - 1],
        ))
    }

    /// `IsTheSame1(Index)` (`cxx:188-203`).
    pub fn is_the_same1(&self, index: usize) -> bool {
        self.check(index);
        self.the_same1[index - 1] != 0
    }

    /// `if (!WellDone) throw StdFail_NotDone()` and
    /// `if (Index <= 0 || Index > NbrSol) throw Standard_OutOfRange()`.
    fn check(&self, index: usize) {
        assert!(
            self.well_done,
            "StdFail_NotDone in Geom2dGcc_Circ2dTanCen"
        );
        assert!(
            index >= 1 && index <= self.nbr_sol,
            "Standard_OutOfRange in Geom2dGcc_Circ2dTanCen"
        );
    }
}
