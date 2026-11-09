//! `Geom2dGcc_Circ2d2TanOn` (`Geom2dGcc_Circ2d2TanOn.{hxx:46-205,cxx:32-1021}`):
//! circles tangent to two operands with the centre on a third one.
//!
//! The class is a front end: it dispatches on the three operand types between
//! the analytic [`GccAnaCirc2d2TanOn`], the geometric
//! [`Geom2dGccCirc2d2TanOnGeo`] and the iterative
//! [`Geom2dGccCirc2d2TanOnIter`], then hands every solution back through the
//! `cirsol` / `qualifier1` / ... arrays of `hxx:189-205`.
//!
//! | OCCT | arguments | `cxx` |
//! |---|---|---|
//! | `(QualifiedCurve, QualifiedCurve, Adaptor, Tol, Param1, Param2, ParamOn)` | `:32-440` |
//! | `(QualifiedCurve, Geom2d_Point, Adaptor, Tol, Param1, ParamOn)` | `:442-686` |
//! | `(Geom2d_Point, Geom2d_Point, Adaptor, Tol)` | `:688-768` |
//! | `Results(GccAna_Circ2d2TanOn)` | `:771-797` |
//! | `Results(Geom2dGcc_Circ2d2TanOnGeo)` | `:799-826` |
//! | `IsDone`..`IsTheSame2` | `:828-1021` |
//!
//! Translation notes:
//! - the `cxx:345-438` arm (at least one tangency operand is a general curve)
//!   runs the `Type3 == Line || Type3 == Circle` block *and then* the
//!   `(Qc1, Qc2, OnCurve)` iterative solver unconditionally, overwriting
//!   `WellDone`, `NbrSol` and every array entry; the port keeps both calls in
//!   that order so the second one wins exactly as in OCCT,
//! - `Invert` is set only inside the `if (WellDone)` blocks (`cxx:189`, `:207`,
//!   `:322`), so a failed solver leaves it false; the port keeps that,
//! - the `cxx:350-379` arm (carrier is a circle) does not call
//!   `Circ.WhichQualifier(...)`, unlike the line arm (`cxx:406`) and the final
//!   block (`cxx:434`); kept as written,
//! - `Geom2dGcc_QualifiedCurve` only stores the curve in this port, so
//!   `Geom2dAdaptor_Curve C1 = Qualified1.Qualified()` (`cxx:54`) becomes
//!   [`Geom2dAdaptorCurve::from_curve`] on that curve: the adaptor then carries
//!   the curve's own `FirstParameter` / `LastParameter`, exactly like the
//!   `Geom2dAdaptor_Curve(const Handle(Geom2d_Curve)&)` constructor OCCT runs
//!   when the facade takes a `QualifiedCurve`,
//! - `Results` calls `Tangency1` / `Tangency2` unconditionally (`cxx:793-794`,
//!   `:822-823`), and both throw `StdFail_NotDone` when `TheSame` is set, so
//!   such a solution aborts the constructor; the port panics accordingly,
//! - `GccAna_Circ2d2TanOn::CenterOn3` answers `pnttg1sol` (`GccAna_Circ2d2TanOn.cxx:484-485`),
//!   so the `pntcen` of the analytic arm receives the first tangency point, not
//!   the centre; the port keeps that,
//! - `StdFail_NotDone` / `Standard_OutOfRange` / `GccEnt_BadQualifier` are
//!   panics, matching the rest of the `gcc` module.
//!
//! UNPORTED: none - all three constructors, both `Results` and all nine
//! accessors.

use occt_core::gp::{GpCirc2d, GpLin2d, GpPnt2d};
use occt_core::kernel::geomabs::CurveType;

use crate::adaptor_curve::{basis2d, type2d, Geom2dAdaptorCurve};
use crate::gcc::circ2d_2tan_on::{GccAnaCirc2d2TanOn, MAX_SOL};
use crate::gcc::circ2d_2tan_on_geo::Geom2dGccCirc2d2TanOnGeo;
use crate::gcc::circ2d_2tan_on_iter::Geom2dGccCirc2d2TanOnIter;
use crate::gcc::gcc_ent::{GccEntPosition, GccEntQualifiedCirc, GccEntQualifiedLin};
use crate::gcc::qualified_curve::{Geom2dGccQCurve, Geom2dGccQualifiedCurve};

/// `Geom2dGcc_Circ2d2TanOn` (`Geom2dGcc_Circ2d2TanOn.hxx:46-205`).
#[derive(Clone)]
pub struct Geom2dGccCirc2d2TanOn {
    /// `WellDone` (`hxx:189`).
    well_done: bool,
    /// `cirsol` (`hxx:190`).
    cirsol: [GpCirc2d; MAX_SOL],
    /// `NbrSol` (`hxx:191`).
    nbr_sol: usize,
    /// `qualifier1` (`hxx:192`).
    qualifier1: [GccEntPosition; MAX_SOL],
    /// `qualifier2` (`hxx:193`).
    qualifier2: [GccEntPosition; MAX_SOL],
    /// `TheSame1` (`hxx:194`).
    the_same1: [i32; MAX_SOL],
    /// `TheSame2` (`hxx:195`).
    the_same2: [i32; MAX_SOL],
    /// `pnttg1sol` (`hxx:196`).
    pnttg1sol: [GpPnt2d; MAX_SOL],
    /// `pnttg2sol` (`hxx:197`).
    pnttg2sol: [GpPnt2d; MAX_SOL],
    /// `pntcen` (`hxx:198`).
    pntcen: [GpPnt2d; MAX_SOL],
    /// `par1sol` (`hxx:199`).
    par1sol: [f64; MAX_SOL],
    /// `par2sol` (`hxx:200`).
    par2sol: [f64; MAX_SOL],
    /// `pararg1` (`hxx:201`).
    pararg1: [f64; MAX_SOL],
    /// `pararg2` (`hxx:202`).
    pararg2: [f64; MAX_SOL],
    /// `parcen3` (`hxx:203`).
    parcen3: [f64; MAX_SOL],
    /// `Invert` (`hxx:204`).
    invert: bool,
}

/// `Geom2d_Circle` / `Geom2d_Line` casts of the facade (`cxx:78`, `:123`,
/// `:166`, `:269`, `:304`, `:516`, `:578`, `:693`, ...).
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

fn as_qualified_lin(q: &Geom2dGccQualifiedCurve) -> GccEntQualifiedLin {
    GccEntQualifiedLin::new(as_line(q), q.qualifier())
}

/// `Geom2dGcc_QCurve Qc1(C1, Qualified1.Qualifier())` (`cxx:347-348`,
/// `:595`, `:658`).
fn as_qcurve(q: &Geom2dGccQualifiedCurve) -> Geom2dGccQCurve {
    Geom2dGccQCurve::new(
        Geom2dAdaptorCurve::from_curve(q.qualified().clone()),
        q.qualifier(),
    )
}

/// `C1.GetType()` (`cxx:58-60`).
fn curve_type(q: &Geom2dGccQualifiedCurve) -> CurveType {
    type2d(q.qualified().as_ref())
}

/// `Type3 == GeomAbs_Line || Type3 == GeomAbs_Circle`.
fn line_or_circle(t: CurveType) -> bool {
    t == CurveType::Line || t == CurveType::Circle
}

/// `TheSame1(j) = 1` / `TheSame1(j) = 0` (`cxx:778-785`) on an `int` array.
fn flag(value: bool) -> i32 {
    if value {
        1
    } else {
        0
    }
}

/// `CCon->Circ2d()` on the `OnCurve` adaptor.
fn on_circle(on_curve: &Geom2dAdaptorCurve) -> GpCirc2d {
    on_curve.circle().expect("Geom2d_Circle down_cast on OnCurve")
}

/// `CCon->Lin2d()` on the `OnCurve` adaptor.
fn on_line(on_curve: &Geom2dAdaptorCurve) -> GpLin2d {
    on_curve.line().expect("Geom2d_Line down_cast on OnCurve")
}

impl Geom2dGccCirc2d2TanOn {
    /// The `: WellDone(false), cirsol(1, 8), ... :` initializer list every
    /// constructor shares (`cxx:40-54`, `:449-463`, `:691-705`).
    fn empty() -> Self {
        Self {
            well_done: false,
            cirsol: [GpCirc2d::default(); MAX_SOL],
            nbr_sol: 0,
            qualifier1: [GccEntPosition::NoQualifier; MAX_SOL],
            qualifier2: [GccEntPosition::NoQualifier; MAX_SOL],
            the_same1: [0; MAX_SOL],
            the_same2: [0; MAX_SOL],
            pnttg1sol: [GpPnt2d::zero(); MAX_SOL],
            pnttg2sol: [GpPnt2d::zero(); MAX_SOL],
            pntcen: [GpPnt2d::zero(); MAX_SOL],
            par1sol: [0.0; MAX_SOL],
            par2sol: [0.0; MAX_SOL],
            pararg1: [0.0; MAX_SOL],
            pararg2: [0.0; MAX_SOL],
            parcen3: [0.0; MAX_SOL],
            invert: false,
        }
    }

    /// `Geom2dGcc_Circ2d2TanOn(Qualified1, Qualified2, OnCurve, Tolerance,
    /// Param1, Param2, ParamOn)` (`cxx:32-440`).
    pub fn new(
        qualified1: &Geom2dGccQualifiedCurve,
        qualified2: &Geom2dGccQualifiedCurve,
        on_curve: &Geom2dAdaptorCurve,
        tolerance: f64,
        param1: f64,
        param2: f64,
        param_on: f64,
    ) -> Self {
        let mut s = Self::empty();
        let type1 = curve_type(qualified1);
        let type2 = curve_type(qualified2);
        let type3 = on_curve
            .get_type()
            .unwrap_or(CurveType::OtherCurve);
        if line_or_circle(type1) && line_or_circle(type2) {
            if line_or_circle(type3) {
                // `cxx:70-254`: the analytic solver.
                let (solver, invert) = match (type1, type2, type3) {
                    (CurveType::Circle, CurveType::Circle, CurveType::Circle) => (
                        GccAnaCirc2d2TanOn::new_circ_circ_on_circ(
                            &as_qualified_circ(qualified1),
                            &as_qualified_circ(qualified2),
                            &on_circle(on_curve),
                            tolerance,
                        ),
                        false,
                    ),
                    (CurveType::Circle, CurveType::Circle, CurveType::Line) => (
                        GccAnaCirc2d2TanOn::new_circ_circ_on_line(
                            &as_qualified_circ(qualified1),
                            &as_qualified_circ(qualified2),
                            &on_line(on_curve),
                            tolerance,
                        ),
                        false,
                    ),
                    (CurveType::Circle, CurveType::Line, CurveType::Circle) => (
                        GccAnaCirc2d2TanOn::new_circ_lin_on_circ(
                            &as_qualified_circ(qualified1),
                            &as_qualified_lin(qualified2),
                            &on_circle(on_curve),
                            tolerance,
                        ),
                        false,
                    ),
                    (CurveType::Circle, CurveType::Line, CurveType::Line) => (
                        GccAnaCirc2d2TanOn::new_circ_lin_on_line(
                            &as_qualified_circ(qualified1),
                            &as_qualified_lin(qualified2),
                            &on_line(on_curve),
                            tolerance,
                        ),
                        false,
                    ),
                    (CurveType::Line, CurveType::Circle, CurveType::Circle) => (
                        GccAnaCirc2d2TanOn::new_circ_lin_on_circ(
                            &as_qualified_circ(qualified2),
                            &as_qualified_lin(qualified1),
                            &on_circle(on_curve),
                            tolerance,
                        ),
                        true,
                    ),
                    (CurveType::Line, CurveType::Circle, CurveType::Line) => (
                        GccAnaCirc2d2TanOn::new_circ_lin_on_line(
                            &as_qualified_circ(qualified2),
                            &as_qualified_lin(qualified1),
                            &on_line(on_curve),
                            tolerance,
                        ),
                        true,
                    ),
                    (CurveType::Line, CurveType::Line, CurveType::Circle) => (
                        GccAnaCirc2d2TanOn::new_lin_lin_on_circ(
                            &as_qualified_lin(qualified1),
                            &as_qualified_lin(qualified2),
                            &on_circle(on_curve),
                        ),
                        false,
                    ),
                    (CurveType::Line, CurveType::Line, CurveType::Line) => (
                        GccAnaCirc2d2TanOn::new_lin_lin_on_line(
                            &as_qualified_lin(qualified1),
                            &as_qualified_lin(qualified2),
                            &on_line(on_curve),
                            tolerance,
                        ),
                        false,
                    ),
                    _ => unreachable!("the three types are line or circle"),
                };
                s.well_done = solver.is_done();
                if s.well_done {
                    s.nbr_sol = solver.nb_solutions();
                    for i in 1..=s.nbr_sol {
                        let (q1, q2) = solver.which_qualifier(i);
                        s.qualifier1[i - 1] = q1;
                        s.qualifier2[i - 1] = q2;
                    }
                    s.results_ana(&solver);
                    s.invert = invert;
                }
            } else {
                // `cxx:258-343`: the geometric solver.
                let (solver, invert) = match (type1, type2) {
                    (CurveType::Circle, CurveType::Circle) => (
                        Geom2dGccCirc2d2TanOnGeo::new_circ_circ(
                            &as_qualified_circ(qualified1),
                            &as_qualified_circ(qualified2),
                            on_curve,
                            tolerance,
                        ),
                        false,
                    ),
                    (CurveType::Circle, CurveType::Line) => (
                        Geom2dGccCirc2d2TanOnGeo::new_circ_lin(
                            &as_qualified_circ(qualified1),
                            &as_qualified_lin(qualified2),
                            on_curve,
                            tolerance,
                        ),
                        false,
                    ),
                    (CurveType::Line, CurveType::Circle) => (
                        Geom2dGccCirc2d2TanOnGeo::new_circ_lin(
                            &as_qualified_circ(qualified2),
                            &as_qualified_lin(qualified1),
                            on_curve,
                            tolerance,
                        ),
                        true,
                    ),
                    (CurveType::Line, CurveType::Line) => (
                        Geom2dGccCirc2d2TanOnGeo::new_lin_lin(
                            &as_qualified_lin(qualified1),
                            &as_qualified_lin(qualified2),
                            on_curve,
                            tolerance,
                        ),
                        false,
                    ),
                    _ => unreachable!("the two types are line or circle"),
                };
                s.well_done = solver.is_done();
                if s.well_done {
                    s.nbr_sol = solver.nb_solutions();
                    for i in 1..=s.nbr_sol {
                        let (q1, q2) = solver.which_qualifier(i);
                        s.qualifier1[i - 1] = q1;
                        s.qualifier2[i - 1] = q2;
                    }
                    s.results_geo(&solver);
                    s.invert = invert;
                }
            }
        } else {
            // `cxx:345-438`.
            let qc1 = as_qcurve(qualified1);
            let qc2 = as_qcurve(qualified2);
            if line_or_circle(type3) {
                if type3 == CurveType::Circle {
                    let circ = on_circle(on_curve);
                    let iter = Geom2dGccCirc2d2TanOnIter::new_cu_cu_on_ci(
                        &qc1, &qc2, &circ, param1, param2, param_on, tolerance,
                    );
                    s.copy_iter_cu_cu_on_ci(&iter);
                } else {
                    let line = on_line(on_curve);
                    let iter = Geom2dGccCirc2d2TanOnIter::new_cu_cu_on_li(
                        &qc1, &qc2, &line, param1, param2, param_on, tolerance,
                    );
                    s.copy_iter_with_qualifier(&iter);
                }
            }
            // `cxx:413-438`: unconditional, and therefore the winning call.
            let iter = Geom2dGccCirc2d2TanOnIter::new_cu_cu_on_cu(
                &qc1, &qc2, on_curve, param1, param2, param_on, tolerance,
            );
            s.copy_iter_with_qualifier(&iter);
        }
        s
    }

    /// `Geom2dGcc_Circ2d2TanOn(Qualified1, Point, OnCurve, Tolerance, Param1,
    /// ParamOn)` (`cxx:442-686`).
    pub fn new_with_point(
        qualified1: &Geom2dGccQualifiedCurve,
        point: &GpPnt2d,
        on_curve: &Geom2dAdaptorCurve,
        tolerance: f64,
        param1: f64,
        param_on: f64,
    ) -> Self {
        let mut s = Self::empty();
        let type1 = curve_type(qualified1);
        let type3 = on_curve
            .get_type()
            .unwrap_or(CurveType::OtherCurve);
        if line_or_circle(type1) {
            if line_or_circle(type3) {
                // `cxx:472-551`: the analytic solver.
                let solver = match (type1, type3) {
                    (CurveType::Circle, CurveType::Circle) => {
                        GccAnaCirc2d2TanOn::new_circ_pnt_on_circ(
                            &as_qualified_circ(qualified1),
                            point,
                            &on_circle(on_curve),
                            tolerance,
                        )
                    }
                    (CurveType::Circle, CurveType::Line) => {
                        GccAnaCirc2d2TanOn::new_circ_pnt_on_line(
                            &as_qualified_circ(qualified1),
                            point,
                            &on_line(on_curve),
                            tolerance,
                        )
                    }
                    (CurveType::Line, CurveType::Circle) => {
                        GccAnaCirc2d2TanOn::new_lin_pnt_on_circ(
                            &as_qualified_lin(qualified1),
                            point,
                            &on_circle(on_curve),
                            tolerance,
                        )
                    }
                    (CurveType::Line, CurveType::Line) => {
                        GccAnaCirc2d2TanOn::new_lin_pnt_on_line(
                            &as_qualified_lin(qualified1),
                            point,
                            &on_line(on_curve),
                            tolerance,
                        )
                    }
                    _ => unreachable!("type1 and type3 are line or circle"),
                };
                s.well_done = solver.is_done();
                if s.well_done {
                    s.nbr_sol = solver.nb_solutions();
                    for i in 1..=s.nbr_sol {
                        let (q1, q2) = solver.which_qualifier(i);
                        s.qualifier1[i - 1] = q1;
                        s.qualifier2[i - 1] = q2;
                    }
                    s.results_ana(&solver);
                }
            } else {
                // `cxx:556-593`: the geometric solver.
                let solver = if type1 == CurveType::Circle {
                    Geom2dGccCirc2d2TanOnGeo::new_circ_pnt(
                        &as_qualified_circ(qualified1),
                        point,
                        on_curve,
                        tolerance,
                    )
                } else {
                    Geom2dGccCirc2d2TanOnGeo::new_lin_pnt(
                        &as_qualified_lin(qualified1),
                        point,
                        on_curve,
                        tolerance,
                    )
                };
                s.well_done = solver.is_done();
                if s.well_done {
                    s.nbr_sol = solver.nb_solutions();
                    for i in 1..=s.nbr_sol {
                        let (q1, q2) = solver.which_qualifier(i);
                        s.qualifier1[i - 1] = q1;
                        s.qualifier2[i - 1] = q2;
                    }
                    s.results_geo(&solver);
                }
            }
        } else {
            // `cxx:595-685`.
            let qc1 = as_qcurve(qualified1);
            if line_or_circle(type3) {
                if type3 == CurveType::Circle {
                    let circ = on_circle(on_curve);
                    let iter = Geom2dGccCirc2d2TanOnIter::new_cu_pt_on_ci(
                        &qc1, point, &circ, param1, param_on, tolerance,
                    );
                    s.copy_iter_cu_pt_on_ci(&iter);
                } else {
                    let line = on_line(on_curve);
                    let iter = Geom2dGccCirc2d2TanOnIter::new_cu_pt_on_li(
                        &qc1, point, &line, param1, param_on, tolerance,
                    );
                    s.copy_iter_cu_pt_on_ci(&iter);
                }
            } else {
                let iter = Geom2dGccCirc2d2TanOnIter::new_cu_pt_on_cu(
                    &qc1, point, on_curve, param1, param_on, tolerance,
                );
                s.copy_iter_full(&iter);
            }
        }
        s
    }

    /// `Geom2dGcc_Circ2d2TanOn(Point1, Point2, OnCurve, Tolerance)`
    /// (`cxx:688-768`).
    pub fn new_with_two_points(
        point1: &GpPnt2d,
        point2: &GpPnt2d,
        on_curve: &Geom2dAdaptorCurve,
        tolerance: f64,
    ) -> Self {
        let mut s = Self::empty();
        let type3 = on_curve
            .get_type()
            .unwrap_or(CurveType::OtherCurve);
        if line_or_circle(type3) {
            let solver = if type3 == CurveType::Circle {
                GccAnaCirc2d2TanOn::new_pnt_pnt_on_circ(
                    point1,
                    point2,
                    &on_circle(on_curve),
                    tolerance,
                )
            } else {
                GccAnaCirc2d2TanOn::new_pnt_pnt_on_line(
                    point1,
                    point2,
                    &on_line(on_curve),
                    tolerance,
                )
            };
            s.well_done = solver.is_done();
            if s.well_done {
                s.nbr_sol = solver.nb_solutions();
                for i in 1..=s.nbr_sol {
                    let (q1, q2) = solver.which_qualifier(i);
                    s.qualifier1[i - 1] = q1;
                    s.qualifier2[i - 1] = q2;
                }
                s.results_ana(&solver);
            }
        } else {
            let solver =
                Geom2dGccCirc2d2TanOnGeo::new_pnt_pnt(point1, point2, on_curve, tolerance);
            s.well_done = solver.is_done();
            if s.well_done {
                s.nbr_sol = solver.nb_solutions();
                for i in 1..=s.nbr_sol {
                    let (q1, q2) = solver.which_qualifier(i);
                    s.qualifier1[i - 1] = q1;
                    s.qualifier2[i - 1] = q2;
                }
                s.results_geo(&solver);
            }
        }
        s
    }

    /// `Results(const GccAna_Circ2d2TanOn&)` (`cxx:771-797`).
    fn results_ana(&mut self, circ: &GccAnaCirc2d2TanOn) {
        for j in 1..=self.nbr_sol {
            self.cirsol[j - 1] = circ.this_solution(j);
            self.the_same1[j - 1] = flag(circ.is_the_same1(j));
            self.the_same2[j - 1] = flag(circ.is_the_same2(j));
            let (q1, q2) = circ.which_qualifier(j);
            self.qualifier1[j - 1] = q1;
            self.qualifier2[j - 1] = q2;
            let (ps, pa, p) = circ
                .tangency1(j)
                .expect("StdFail_NotDone in GccAna_Circ2d2TanOn::Tangency1");
            self.par1sol[j - 1] = ps;
            self.pararg1[j - 1] = pa;
            self.pnttg1sol[j - 1] = p;
            let (ps, pa, p) = circ
                .tangency2(j)
                .expect("StdFail_NotDone in GccAna_Circ2d2TanOn::Tangency2");
            self.par2sol[j - 1] = ps;
            self.pararg2[j - 1] = pa;
            self.pnttg2sol[j - 1] = p;
            let (pa, p) = circ.center_on3(j);
            self.parcen3[j - 1] = pa;
            self.pntcen[j - 1] = p;
        }
    }

    /// `Results(const Geom2dGcc_Circ2d2TanOnGeo&)` (`cxx:799-826`).
    fn results_geo(&mut self, circ: &Geom2dGccCirc2d2TanOnGeo) {
        for j in 1..=self.nbr_sol {
            self.cirsol[j - 1] = circ.this_solution(j);
            self.the_same1[j - 1] = flag(circ.is_the_same1(j));
            self.the_same2[j - 1] = flag(circ.is_the_same2(j));
            let (q1, q2) = circ.which_qualifier(j);
            self.qualifier1[j - 1] = q1;
            self.qualifier2[j - 1] = q2;
            let (ps, pa, p) = circ
                .tangency1(j)
                .expect("StdFail_NotDone in Geom2dGcc_Circ2d2TanOnGeo::Tangency1");
            self.par1sol[j - 1] = ps;
            self.pararg1[j - 1] = pa;
            self.pnttg1sol[j - 1] = p;
            let (ps, pa, p) = circ
                .tangency2(j)
                .expect("StdFail_NotDone in Geom2dGcc_Circ2d2TanOnGeo::Tangency2");
            self.par2sol[j - 1] = ps;
            self.pararg2[j - 1] = pa;
            self.pnttg2sol[j - 1] = p;
            let (pa, p) = circ.center_on3(j);
            self.parcen3[j - 1] = pa;
            self.pntcen[j - 1] = p;
        }
    }

    /// The `cxx:350-379` arm: a single iterative solution over a circle
    /// carrier. `IsTheSame1` / `IsTheSame2` and both tangencies, but no
    /// `WhichQualifier`.
    fn copy_iter_cu_cu_on_ci(&mut self, iter: &Geom2dGccCirc2d2TanOnIter) {
        self.well_done = iter.is_done();
        if self.well_done {
            self.nbr_sol = 1;
            self.cirsol[0] = iter.this_solution();
            self.the_same1[0] = flag(iter.is_the_same1());
            self.the_same2[0] = flag(iter.is_the_same2());
            let (ps, pa, p) = iter.tangency1();
            self.par1sol[0] = ps;
            self.pararg1[0] = pa;
            self.pnttg1sol[0] = p;
            let (ps, pa, p) = iter.tangency2();
            self.par2sol[0] = ps;
            self.pararg2[0] = pa;
            self.pnttg2sol[0] = p;
        }
    }

    /// The `cxx:381-410` and `cxx:413-438` arms: one solution, `IsTheSame1` /
    /// `IsTheSame2`, `WhichQualifier` and both tangencies.
    fn copy_iter_with_qualifier(&mut self, iter: &Geom2dGccCirc2d2TanOnIter) {
        self.well_done = iter.is_done();
        if self.well_done {
            self.nbr_sol = 1;
            self.cirsol[0] = iter.this_solution();
            self.the_same1[0] = flag(iter.is_the_same1());
            self.the_same2[0] = flag(iter.is_the_same2());
            let (q1, q2) = iter.which_qualifier();
            self.qualifier1[0] = q1;
            self.qualifier2[0] = q2;
            let (ps, pa, p) = iter.tangency1();
            self.par1sol[0] = ps;
            self.pararg1[0] = pa;
            self.pnttg1sol[0] = p;
            let (ps, pa, p) = iter.tangency2();
            self.par2sol[0] = ps;
            self.pararg2[0] = pa;
            self.pnttg2sol[0] = p;
        }
    }

    /// The `cxx:600-655` arms of the point constructor: one solution,
    /// `IsTheSame1`, `WhichQualifier` and both tangencies. `TheSame2` is not
    /// touched.
    fn copy_iter_cu_pt_on_ci(&mut self, iter: &Geom2dGccCirc2d2TanOnIter) {
        self.well_done = iter.is_done();
        if self.well_done {
            self.nbr_sol = 1;
            self.cirsol[0] = iter.this_solution();
            self.the_same1[0] = flag(iter.is_the_same1());
            let (q1, q2) = iter.which_qualifier();
            self.qualifier1[0] = q1;
            self.qualifier2[0] = q2;
            let (ps, pa, p) = iter.tangency1();
            self.par1sol[0] = ps;
            self.pararg1[0] = pa;
            self.pnttg1sol[0] = p;
            let (ps, pa, p) = iter.tangency2();
            self.par2sol[0] = ps;
            self.pararg2[0] = pa;
            self.pnttg2sol[0] = p;
        }
    }

    /// The `cxx:657-683` arm of the point constructor: one solution with both
    /// `IsTheSame` queries, `WhichQualifier` and both tangencies.
    fn copy_iter_full(&mut self, iter: &Geom2dGccCirc2d2TanOnIter) {
        self.well_done = iter.is_done();
        if self.well_done {
            self.nbr_sol = 1;
            self.cirsol[0] = iter.this_solution();
            self.the_same1[0] = flag(iter.is_the_same1());
            self.the_same2[0] = flag(iter.is_the_same2());
            let (q1, q2) = iter.which_qualifier();
            self.qualifier1[0] = q1;
            self.qualifier2[0] = q2;
            let (ps, pa, p) = iter.tangency1();
            self.par1sol[0] = ps;
            self.pararg1[0] = pa;
            self.pnttg1sol[0] = p;
            let (ps, pa, p) = iter.tangency2();
            self.par2sol[0] = ps;
            self.pararg2[0] = pa;
            self.pnttg2sol[0] = p;
        }
    }

    /// `IsDone()` (`cxx:828-831`).
    pub fn is_done(&self) -> bool {
        self.well_done
    }

    /// `NbSolutions()` (`cxx:833-836`).
    pub fn nb_solutions(&self) -> usize {
        self.nbr_sol
    }

    /// `ThisSolution(Index)` (`cxx:838-848`).
    pub fn this_solution(&self, index: usize) -> GpCirc2d {
        self.check(index);
        self.cirsol[index - 1]
    }

    /// `WhichQualifier(Index, Qualif1, Qualif2)` (`cxx:850-877`).
    pub fn which_qualifier(&self, index: usize) -> (GccEntPosition, GccEntPosition) {
        self.check(index);
        if self.invert {
            (self.qualifier2[index - 1], self.qualifier1[index - 1])
        } else {
            (self.qualifier1[index - 1], self.qualifier2[index - 1])
        }
    }

    /// `Tangency1(Index, ParSol, ParArg, PntSol)` (`cxx:879-918`).
    ///
    /// `None` where OCCT raises `StdFail_NotDone` for the selected `TheSame`.
    pub fn tangency1(&self, index: usize) -> Option<(f64, f64, GpPnt2d)> {
        self.check(index);
        if self.invert {
            if self.the_same2[index - 1] == 0 {
                Some((
                    self.par2sol[index - 1],
                    self.pararg2[index - 1],
                    self.pnttg2sol[index - 1],
                ))
            } else {
                None
            }
        } else if self.the_same1[index - 1] == 0 {
            Some((
                self.par1sol[index - 1],
                self.pararg1[index - 1],
                self.pnttg1sol[index - 1],
            ))
        } else {
            None
        }
    }

    /// `Tangency2(Index, ParSol, ParArg, PntSol)` (`cxx:920-962`).
    pub fn tangency2(&self, index: usize) -> Option<(f64, f64, GpPnt2d)> {
        self.check(index);
        if !self.invert {
            if self.the_same2[index - 1] == 0 {
                Some((
                    self.par2sol[index - 1],
                    self.pararg2[index - 1],
                    self.pnttg2sol[index - 1],
                ))
            } else {
                None
            }
        } else if self.the_same1[index - 1] == 0 {
            Some((
                self.par1sol[index - 1],
                self.pararg1[index - 1],
                self.pnttg1sol[index - 1],
            ))
        } else {
            None
        }
    }

    /// `CenterOn3(Index, ParArg, PntSol)` (`cxx:964-980`).
    pub fn center_on3(&self, index: usize) -> (f64, GpPnt2d) {
        self.check(index);
        (self.parcen3[index - 1], self.pntcen[index - 1])
    }

    /// `IsTheSame1(Index)` (`cxx:982-1000`).
    pub fn is_the_same1(&self, index: usize) -> bool {
        self.check(index);
        if self.invert {
            self.the_same2[index - 1] != 0
        } else {
            self.the_same1[index - 1] != 0
        }
    }

    /// `IsTheSame2(Index)` (`cxx:1002-1021`).
    pub fn is_the_same2(&self, index: usize) -> bool {
        self.check(index);
        if !self.invert {
            self.the_same2[index - 1] != 0
        } else {
            self.the_same1[index - 1] != 0
        }
    }

    /// `if (!WellDone) throw StdFail_NotDone()` and
    /// `if (Index <= 0 || Index > NbrSol) throw Standard_OutOfRange()`, the
    /// guard of every accessor above (`cxx:839-847`, `:855-862`, ...).
    fn check(&self, index: usize) {
        assert!(
            self.well_done,
            "StdFail_NotDone in Geom2dGcc_Circ2d2TanOn"
        );
        assert!(
            index >= 1 && index <= self.nbr_sol,
            "Standard_OutOfRange in Geom2dGcc_Circ2d2TanOn"
        );
    }
}

