//! `Geom2dGcc_Circ2d2TanRad`: the `Geom2dAdaptor_Curve` front end of
//! [`super::circ2d_2tan_rad::GccAnaCirc2d2TanRad`].
//!
//! | OCCT | lines |
//! |---|---|
//! | `Geom2dGcc_Circ2d2TanRad.cxx` (curve, curve) | `cxx:52-289` |
//! | `Results(GccAna_Circ2d2TanRad)` | `cxx:412-437` |
//! | `IsDone`..`IsTheSame2` | `cxx:463-641` |
//!
//! Translation notes:
//! - `cxx:78-190` (both arguments a line or a circle) is the only arm the
//!   `ChFi2d` pipeline reaches: `ChFi2d_Builder::BuildFilletEdge` rejects any
//!   other p-curve with `ChFi2d_NotAuthorized` through `IsLineOrCircle`
//!   (`ChFi2d_Builder.cxx:1269-1288`),
//! - the `GccGeo` arm (`cxx:192-288`), i.e. `Geom2dGcc_Circ2d2TanRadGeo`
//!   (`Geom2dGcc_Circ2d2TanRadGeo.cxx:1-1141`) and everything it pulls in
//!   (`Geom2dGcc_QCurve`, `Geom2dGcc_FunctionTanCuCu`, `Geom2dInt_GInter`), is
//!   **UNPORTED**; reaching it panics with the OCCT file and line range rather
//!   than silently answering `!IsDone()`,
//! - the point overloads (`cxx:291-409`) are **UNPORTED** together with
//!   `GccAna_Circ2d2TanRad_2.cxx` / `_5.cxx`, which they delegate to.

use crate::adaptor_curve::{basis2d, type2d};
use crate::gcc::circ2d_2tan_rad::GccAnaCirc2d2TanRad;
use crate::gcc::gcc_ent::{GccEntPosition, GccEntQualifiedCirc, GccEntQualifiedLin};
use crate::gcc::qualified_curve::Geom2dGccQualifiedCurve;
use occt_core::gp::circ2d::GpCirc2d;
use occt_core::gp::pnt2d::GpPnt2d;
use occt_core::kernel::geomabs::CurveType;

/// `Geom2dGcc_Circ2d2TanRad` (`Geom2dGcc_Circ2d2TanRad.hxx:38-208`).
pub struct Geom2dGccCirc2d2TanRad {
    /// `Invert`: set when the constructor had to swap its two arguments so the
    /// circle one came first (`cxx:152`, `:236`, `:259`).
    invert: bool,
    /// The analytic solver that produced the results. OCCT additionally keeps
    /// a `Geom2dGcc_Circ2d2TanRadGeo` member for the `GccGeo` arm, which is
    /// UNPORTED (see the module notes).
    solver: GccAnaCirc2d2TanRad,
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
        // `Geom2dAdaptor_Curve C1 = Qualified1.Qualified(); ... GetType()`
        // (`cxx:74-78`).
        let type1 = type2d(qualified1.qualified().as_ref());
        let type2 = type2d(qualified2.qualified().as_ref());
        let line_or_circle = |t: CurveType| t == CurveType::Line || t == CurveType::Circle;
        if !line_or_circle(type1) || !line_or_circle(type2) {
            panic!(
                "UNPORTED Geom2dGcc_Circ2d2TanRadGeo arm of Geom2dGcc_Circ2d2TanRad.cxx:192-288"
            );
        }

        // `cxx:84-190`: dispatch on the argument types, swapping so the circle
        // comes first and recording `Invert`.
        let circle_of = |q: &Geom2dGccQualifiedCurve| -> GccEntQualifiedCirc {
            let c = basis2d(q.qualified().as_ref())
                .gp_circ2d()
                .expect("Geom2d_Circle basis");
            GccEntQualifiedCirc::new(c, q.qualifier())
        };
        let line_of = |q: &Geom2dGccQualifiedCurve| -> GccEntQualifiedLin {
            let l = basis2d(q.qualified().as_ref())
                .gp_lin2d()
                .expect("Geom2d_Line basis");
            GccEntQualifiedLin::new(l, q.qualifier())
        };

        if type1 == CurveType::Circle {
            let qc1 = circle_of(qualified1);
            if type2 == CurveType::Circle {
                let qc2 = circle_of(qualified2);
                return Self {
                    invert: false,
                    solver: GccAnaCirc2d2TanRad::new_circ_circ(&qc1, &qc2, radius, tolerance),
                };
            }
            // `if (!Qualified2.IsEnclosing()) ... else throw GccEnt_BadQualifier()`
            // (`cxx:116-133`).
            if qualified2.is_enclosing() {
                panic!("GccEnt_BadQualifier in Geom2dGcc_Circ2d2TanRad");
            }
            let ql2 = line_of(qualified2);
            return Self {
                invert: false,
                solver: GccAnaCirc2d2TanRad::new_circ_lin(&qc1, &ql2, radius, tolerance),
            };
        }

        let ql1 = line_of(qualified1);
        if qualified1.is_enclosing() {
            panic!("GccEnt_BadQualifier in Geom2dGcc_Circ2d2TanRad");
        }
        if type2 == CurveType::Circle {
            // `cxx:146-163`: `Invert = true`, the circle argument goes first.
            let qc2 = circle_of(qualified2);
            return Self {
                invert: true,
                solver: GccAnaCirc2d2TanRad::new_circ_lin(&qc2, &ql1, radius, tolerance),
            };
        }
        if qualified2.is_enclosing() {
            panic!("GccEnt_BadQualifier in Geom2dGcc_Circ2d2TanRad");
        }
        let ql2 = line_of(qualified2);
        Self {
            invert: false,
            solver: GccAnaCirc2d2TanRad::new_lin_lin(&ql1, &ql2, radius, tolerance),
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
