//! `Geom2dGcc_Circ2d2TanRad`: the `Geom2dAdaptor_Curve` front end of
//! [`super::circ2d_2tan_rad::GccAnaCirc2d2TanRad`] and of
//! [`super::circ2d_2tan_rad_geo::Geom2dGccCirc2d2TanRadGeo`].
//!
//! | OCCT | lines |
//! |---|---|
//! | `Geom2dGcc_Circ2d2TanRad.cxx` (curve, curve) | `cxx:52-289` |
//! | `Geom2dGcc_Circ2d2TanRad.cxx` (curve, point) | `cxx:291-370` |
//! | `Geom2dGcc_Circ2d2TanRad.cxx` (point, point) | `cxx:372-410` |
//! | `Results(GccAna_Circ2d2TanRad)` | `cxx:412-437` |
//! | `Results(Geom2dGcc_Circ2d2TanRadGeo)` | `cxx:439-462` |
//! | `IsDone`..`IsTheSame2` | `cxx:463-641` |
//!
//! Translation notes:
//! - the constructor dispatches on the arithmetic types of the two operands
//!   (`cxx:84-190`): when both are a line or a circle the analytic solver
//!   answers, otherwise the `GccGeo` arm (`cxx:192-288`) picks the two-operand
//!   general-curve constructor and records `Invert` when it had to swap so the
//!   line or circle came first,
//! - `Radius < 0.0` is tested before any dispatch (`cxx:68-70`, `:306-310`,
//!   `:388-392`), so a negative radius raises `Standard_NegativeValue` even
//!   when a qualifier is also bad,
//! - the `GccGeo` arm rejects an `Enclosing` line operand with
//!   `GccEnt_BadQualifier` (`cxx:200-205`, `:239-243`),
//! - the curve/point constructor has no bad-qualifier guard at all
//!   (`cxx:311-369`): a bad qualifier surfaces from the solver it delegates to,
//! - the point/point constructor sizes its arrays `(1, 2)` (`cxx:377-390`)
//!   rather than `aNbSolMAX`, so at most two solutions are stored,
//! - `Invert` swaps the operand order in `WhichQualifier`, `Tangency1`,
//!   `Tangency2`, `IsTheSame1` and `IsTheSame2` (`cxx:500-641`).
//!
//! UNPORTED: none - the three constructors, both `Results` helpers and the six
//! accessors.

use crate::adaptor_curve::{basis2d, type2d};
use crate::gcc::circ2d_2tan_rad::GccAnaCirc2d2TanRad;
use crate::gcc::circ2d_2tan_rad_geo::Geom2dGccCirc2d2TanRadGeo;
use crate::gcc::gcc_ent::{GccEntPosition, GccEntQualifiedCirc, GccEntQualifiedLin};
use crate::gcc::qualified_curve::{Geom2dGccQCurve, Geom2dGccQualifiedCurve};
use crate::adaptor_curve::Geom2dAdaptorCurve;
use occt_core::gp::circ2d::GpCirc2d;
use occt_core::gp::lin2d::GpLin2d;
use occt_core::gp::pnt2d::GpPnt2d;
use occt_core::kernel::geomabs::CurveType;

/// The solver a `Geom2dGcc_Circ2d2TanRad` constructor selected.
///
/// OCCT keeps both members (`GccAna_Circ2d2TanRad` used by the analytic arms and
/// `Geom2dGcc_Circ2d2TanRadGeo` by the `GccGeo` arm) and fills one of them per
/// constructor; the port stores the one that answered, because every accessor
/// reads through `Results` in both cases.
enum Solver {    Ana(GccAnaCirc2d2TanRad),
    Geo(Geom2dGccCirc2d2TanRadGeo),
}

impl Solver {
    fn is_done(&self) -> bool {
        match self {
            Solver::Ana(s) => s.is_done(),
            Solver::Geo(s) => s.is_done(),
        }
    }

    fn nb_solutions(&self) -> usize {
        match self {
            Solver::Ana(s) => s.nb_solutions(),
            Solver::Geo(s) => s.nb_solutions(),
        }
    }

    fn this_solution(&self, index: usize) -> GpCirc2d {
        match self {
            Solver::Ana(s) => s.this_solution(index),
            Solver::Geo(s) => s.this_solution(index),
        }
    }

    fn which_qualifier(&self, index: usize) -> (GccEntPosition, GccEntPosition) {
        match self {
            Solver::Ana(s) => s.which_qualifier(index),
            Solver::Geo(s) => s.which_qualifier(index),
        }
    }

    /// `Circ.Tangency1(j, par1sol(j), pararg1(j), pnttg1sol(j))`
    /// (`cxx:434`, `:459`): the `Geom2dGcc_Circ2d2TanRadGeo` arm never sets a
    /// non-zero `TheSame1`, so it cannot raise `StdFail_NotDone` there.
    fn tangency1(&self, index: usize) -> Option<(f64, f64, GpPnt2d)> {
        match self {
            Solver::Ana(s) => s.tangency1(index),
            Solver::Geo(s) => Some(s.tangency1(index)),
        }
    }

    fn tangency2(&self, index: usize) -> Option<(f64, f64, GpPnt2d)> {
        match self {
            Solver::Ana(s) => s.tangency2(index),
            Solver::Geo(s) => Some(s.tangency2(index)),
        }
    }

    fn is_the_same1(&self, index: usize) -> bool {
        match self {
            Solver::Ana(s) => s.is_the_same1(index),
            Solver::Geo(s) => s.is_the_same1(index),
        }
    }

    fn is_the_same2(&self, index: usize) -> bool {
        match self {
            Solver::Ana(s) => s.is_the_same2(index),
            Solver::Geo(s) => s.is_the_same2(index),
        }
    }
}

/// `Geom2dGcc_Circ2d2TanRad` (`Geom2dGcc_Circ2d2TanRad.hxx:38-208`).
pub struct Geom2dGccCirc2d2TanRad {
    /// `Invert`: set when the constructor had to swap its two arguments so the
    /// circle one came first (`cxx:152`, `:236`, `:259`, `:264`).
    invert: bool,
    /// The solver that produced the results.
    solver: Solver,
}

/// `gp_Circ2d c(CCC->Circ2d())` after the `Geom2d_Circle` down cast
/// (`cxx:91-92`, `:221-222`, `:263-264`).
fn as_circle(q: &Geom2dGccQualifiedCurve) -> GpCirc2d {
    basis2d(q.qualified().as_ref())
        .gp_circ2d()
        .expect("Geom2d_Circle down_cast")
}

/// `gp_Lin2d l(LLL->Lin2d())` after the `Geom2d_Line` down cast
/// (`cxx:136-137`, `:198-199`, `:236-237`).
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

/// `Geom2dGcc_QCurve Qc(C, Qualified.Qualifier())` (`cxx:206-207`, `:225-226`,
/// `:245-246`, `:266-267`, `:275-276`, `:359`).
fn as_qcurve(q: &Geom2dGccQualifiedCurve) -> Geom2dGccQCurve {
    Geom2dGccQCurve::new(
        Geom2dAdaptorCurve::from_curve(q.qualified().clone()),
        q.qualifier(),
    )
}

/// `C.GetType()` (`cxx:78-79`).
fn curve_type(q: &Geom2dGccQualifiedCurve) -> CurveType {
    type2d(q.qualified().as_ref())
}

fn line_or_circle(t: CurveType) -> bool {
    t == CurveType::Line || t == CurveType::Circle
}

fn check_radius(radius: f64) {
    // `if (Radius < 0.) throw Standard_NegativeValue()` (`cxx:68-70`,
    // `:306-310`, `:388-392`): the raw sign test, so NaN passes.
    assert!(
        !(radius < 0.0),
        "Standard_NegativeValue in Geom2dGcc_Circ2d2TanRad"
    );
}

fn bad_qualifier() -> ! {
    panic!("GccEnt_BadQualifier in Geom2dGcc_Circ2d2TanRad");
}

impl Geom2dGccCirc2d2TanRad {
    /// `Geom2dGcc_Circ2d2TanRad(Qualified1, Qualified2, Radius, Tolerance)`
    /// (`cxx:52-289`).
    pub fn new(
        qualified1: &Geom2dGccQualifiedCurve,
        qualified2: &Geom2dGccQualifiedCurve,
        radius: f64,
        tolerance: f64,
    ) -> Self {
        check_radius(radius);
        let type1 = curve_type(qualified1);
        let type2 = curve_type(qualified2);
        // `cxx:86`: both operands analytic -> `GccAna`.
        if line_or_circle(type1) && line_or_circle(type2) {
            if type1 == CurveType::Circle {
                let qc1 = as_qualified_circ(qualified1);
                if type2 == CurveType::Circle {
                    let qc2 = as_qualified_circ(qualified2);
                    return Self {
                        invert: false,
                        solver: Solver::Ana(GccAnaCirc2d2TanRad::new_circ_circ(
                            &qc1, &qc2, radius, tolerance,
                        )),
                    };
                }
                // `if (!Qualified2.IsEnclosing()) ... else throw`
                // (`cxx:116-133`).
                if qualified2.is_enclosing() {
                    bad_qualifier();
                }
                let ql2 = as_qualified_lin(qualified2);
                return Self {
                    invert: false,
                    solver: Solver::Ana(GccAnaCirc2d2TanRad::new_circ_lin(
                        &qc1, &ql2, radius, tolerance,
                    )),
                };
            }
            // `cxx:139-144`: a line operand cannot be `Enclosing`.
            if qualified1.is_enclosing() {
                bad_qualifier();
            }
            let ql1 = as_qualified_lin(qualified1);
            if type2 == CurveType::Circle {
                // `Invert = true` (`cxx:152`): the circle argument goes first.
                let qc2 = as_qualified_circ(qualified2);
                return Self {
                    invert: true,
                    solver: Solver::Ana(GccAnaCirc2d2TanRad::new_circ_lin(
                        &qc2, &ql1, radius, tolerance,
                    )),
                };
            }
            if qualified2.is_enclosing() {
                bad_qualifier();
            }
            let ql2 = as_qualified_lin(qualified2);
            return Self {
                invert: false,
                solver: Solver::Ana(GccAnaCirc2d2TanRad::new_lin_lin(
                    &ql1, &ql2, radius, tolerance,
                )),
            };
        }
        // `cxx:192-288`: `GccGeo`.
        if type1 == CurveType::Line {
            if qualified1.is_enclosing() {
                bad_qualifier();
            }
            let ql1 = as_qualified_lin(qualified1);
            let qc2 = as_qcurve(qualified2);
            return Self {
                invert: false,
                solver: Solver::Geo(Geom2dGccCirc2d2TanRadGeo::new_lin_qcurve(
                    &ql1, &qc2, radius, tolerance,
                )),
            };
        }
        if type1 == CurveType::Circle {
            let qc1 = as_qualified_circ(qualified1);
            let qc2 = as_qcurve(qualified2);
            return Self {
                invert: false,
                solver: Solver::Geo(Geom2dGccCirc2d2TanRadGeo::new_circ_qcurve(
                    &qc1, &qc2, radius, tolerance,
                )),
            };
        }
        if type2 == CurveType::Line {
            // `Invert = true` (`cxx:236`).
            if qualified2.is_enclosing() {
                bad_qualifier();
            }
            let ql2 = as_qualified_lin(qualified2);
            let qc1 = as_qcurve(qualified1);
            return Self {
                invert: true,
                solver: Solver::Geo(Geom2dGccCirc2d2TanRadGeo::new_lin_qcurve(
                    &ql2, &qc1, radius, tolerance,
                )),
            };
        }
        if type2 == CurveType::Circle {
            // `Invert = true` (`cxx:259`).
            let qc2 = as_qualified_circ(qualified2);
            let qc1 = as_qcurve(qualified1);
            return Self {
                invert: true,
                solver: Solver::Geo(Geom2dGccCirc2d2TanRadGeo::new_circ_qcurve(
                    &qc2, &qc1, radius, tolerance,
                )),
            };
        }
        let qc1 = as_qcurve(qualified1);
        let qc2 = as_qcurve(qualified2);
        Self {
            invert: false,
            solver: Solver::Geo(Geom2dGccCirc2d2TanRadGeo::new_two_qcurves(
                &qc1, &qc2, radius, tolerance,
            )),
        }
    }

    /// `Geom2dGcc_Circ2d2TanRad(Qualified1, Point, Radius, Tolerance)`
    /// (`cxx:291-370`).
    pub fn new_curve_point(
        qualified1: &Geom2dGccQualifiedCurve,
        point: &GpPnt2d,
        radius: f64,
        tolerance: f64,
    ) -> Self {
        check_radius(radius);
        let type1 = curve_type(qualified1);
        if type1 == CurveType::Circle {
            let qc1 = as_qualified_circ(qualified1);
            return Self {
                invert: false,
                solver: Solver::Ana(GccAnaCirc2d2TanRad::new_circ_pnt(
                    &qc1, point, radius, tolerance,
                )),
            };
        }
        if type1 == CurveType::Line {
            let ql1 = as_qualified_lin(qualified1);
            return Self {
                invert: false,
                solver: Solver::Ana(GccAnaCirc2d2TanRad::new_lin_pnt(
                    &ql1, point, radius, tolerance,
                )),
            };
        }
        // `cxx:356-367`: `GccGeo`.
        let qc1 = as_qcurve(qualified1);
        Self {
            invert: false,
            solver: Solver::Geo(Geom2dGccCirc2d2TanRadGeo::new_qcurve_pnt(
                &qc1, point, radius, tolerance,
            )),
        }
    }

    /// `Geom2dGcc_Circ2d2TanRad(Point1, Point2, Radius, Tolerance)`
    /// (`cxx:372-410`).
    pub fn new_point_point(
        point1: &GpPnt2d,
        point2: &GpPnt2d,
        radius: f64,
        tolerance: f64,
    ) -> Self {
        check_radius(radius);
        Self {
            invert: false,
            solver: Solver::Ana(GccAnaCirc2d2TanRad::new_pnt_pnt(
                point1, point2, radius, tolerance,
            )),
        }
    }

    /// `IsDone()` (`cxx:463-466`).
    pub fn is_done(&self) -> bool {
        self.solver.is_done()
    }

    /// `NbSolutions()` (`cxx:468-471`).
    pub fn nb_solutions(&self) -> usize {
        self.solver.nb_solutions()
    }

    /// `ThisSolution(Index)` (`cxx:473-486`).
    pub fn this_solution(&self, index: usize) -> GpCirc2d {
        self.solver.this_solution(index)
    }

    /// `WhichQualifier(Index, Qualif1, Qualif2)` (`cxx:488-512`).
    pub fn which_qualifier(&self, index: usize) -> (GccEntPosition, GccEntPosition) {
        let (q1, q2) = self.solver.which_qualifier(index);
        if self.invert {
            (q2, q1)
        } else {
            (q1, q2)
        }
    }

    /// `Tangency1(Index, ParSol, ParArg, PntSol)` (`cxx:514-556`).
    ///
    /// `None` where OCCT raises `StdFail_NotDone` (`cxx:539`, `:549`).
    pub fn tangency1(&self, index: usize) -> Option<(f64, f64, GpPnt2d)> {
        if self.invert {
            self.solver.tangency2(index)
        } else {
            self.solver.tangency1(index)
        }
    }

    /// `Tangency2(Index, ParSol, ParArg, PntSol)` (`cxx:558-600`).
    pub fn tangency2(&self, index: usize) -> Option<(f64, f64, GpPnt2d)> {
        if self.invert {
            self.solver.tangency1(index)
        } else {
            self.solver.tangency2(index)
        }
    }

    /// `IsTheSame1(Index)` (`cxx:602-621`).
    pub fn is_the_same1(&self, index: usize) -> bool {
        if self.invert {
            self.solver.is_the_same2(index)
        } else {
            self.solver.is_the_same1(index)
        }
    }

    /// `IsTheSame2(Index)` (`cxx:623-641`).
    pub fn is_the_same2(&self, index: usize) -> bool {
        if self.invert {
            self.solver.is_the_same1(index)
        } else {
            self.solver.is_the_same2(index)
        }
    }
}
