//! `GccAna_Circ2d2TanOn`: circles tangent to two arguments with the centre on
//! a third one.
//!
//! | OCCT | argument list | lines | port |
//! |---|---|---|---|
//! | `GccAna_Circ2d2TanOn.cxx` | circ, circ, on line | `cxx:42-378` | `circ_circ_on_line.rs` |
//! | `GccAna_Circ2d2TanOn_1.cxx` | circ, lin, on line | `cxx:47-345` | `circ_lin_on_line.rs` |
//! | `GccAna_Circ2d2TanOn_2.cxx` | lin, lin, on line | `cxx:30-189` | `lin_lin_on_line.rs` |
//! | `GccAna_Circ2d2TanOn_3.cxx` | circ, pnt, on line | `cxx:48-276` | `circ_pnt_on_line.rs` |
//! | `GccAna_Circ2d2TanOn_4.cxx` | lin, pnt, on line | `cxx:51-241` | `lin_pnt_on_line.rs` |
//! | `GccAna_Circ2d2TanOn_5.cxx` | pnt, pnt, on line | `cxx:38-99` | `pnt_pnt_on_line.rs` |
//! | `GccAna_Circ2d2TanOn_6.cxx` | circ, circ, on circ | `cxx:45-401` | `circ_circ_on_circ.rs` |
//! | `GccAna_Circ2d2TanOn_7.cxx` | circ, lin, on circ | `cxx:34-327` | `circ_lin_on_circ.rs` |
//! | `GccAna_Circ2d2TanOn_8.cxx` | circ, pnt, on circ | `cxx:42-315` | `circ_pnt_on_circ.rs` |
//! | `GccAna_Circ2d2TanOn_9.cxx` | lin, lin, on circ | `cxx:31-379` | `lin_lin_on_circ.rs` |
//! | `GccAna_Circ2d2TanOn_10.cxx` | lin, pnt, on circ | `cxx:47-235` | `lin_pnt_on_circ.rs` |
//! | `GccAna_Circ2d2TanOn_11.cxx` | pnt, pnt, on circ | `cxx:47-123` | `pnt_pnt_on_circ.rs` |
//! | `IsDone`..`IsTheSame2` | | `cxx:381-523` | this file |
//!
//! Translation notes:
//! - OCCT sizes the fifteen `NCollection_Array1` fields per constructor:
//!   `(1, 8)` in circ/circ/on circ (`_6.cxx:49-62`), `(1, 4)` in the eight
//!   other `(circ|lin|pnt)` / `(circ|lin|pnt)` over line or circle
//!   constructors, `(1, 2)` in lin/lin/on line, pnt/pnt/on line and
//!   pnt/pnt/on circ. The port sizes every field at [`MAX_SOL`] (= 8, the
//!   largest) so the accessors do not depend on which constructor ran; each
//!   constructor keeps its own OCCT upper bound where the source clamps
//!   (`_8.cxx:263-266` `NbrSol == 4`, `_2.cxx:96` `NbrSol < 2`),
//! - the local `NCollection_Array1<double> Radius(1, 2)` / `Rbid(1, 2)` /
//!   `RBid(1, 2)` / `Rradius(1, 2)` of the constructors become
//!   `[f64; MAX_SOL]`; OCCT writes `Radius(3)`/`Radius(4)` past the end of the
//!   `(1, 2)` arrays in the general case of `GccAna_Circ2d2TanOn.cxx:299` and
//!   `_6.cxx:320` (`nnsol` reaches 4 with two `Unqualified` arguments), and
//!   reads an uninitialised `Rbid(2)`/`Rbid(3)` there; the port's zero-filled
//!   array is the only defined reading of that uninitialised slot,
//! - `NCollection_Array1` is 1-based; the port keeps OCCT's 1-based loop
//!   variables and offsets only when indexing the Rust arrays,
//! - the `GccEnt_BadQualifier`, `StdFail_NotDone` and `Standard_OutOfRange`
//!   raises have no Rust equivalent at the call sites of this package and panic
//!   with the OCCT exception name, matching [`super::circ2d_3tan`],
//! - `CenterOn3` returns `pnttg1sol(Index)` (`cxx:484-485`), not `pntcen`; the
//!   port keeps that and the accessor is documented where it is declared.
//!
//! UNPORTED: none - all twelve constructors and all nine accessors.

mod circ_circ_on_circ;
mod circ_circ_on_line;
mod circ_lin_on_circ;
mod circ_lin_on_line;
mod circ_pnt_on_circ;
mod circ_pnt_on_line;
mod lin_lin_on_circ;
mod lin_lin_on_line;
mod lin_pnt_on_circ;
mod lin_pnt_on_line;
mod pnt_pnt_on_circ;
mod pnt_pnt_on_line;

use occt_core::gp::ax22d::GpAx22d;
use occt_core::gp::circ2d::GpCirc2d;
use occt_core::gp::dir2d::GpDir2d;
use occt_core::gp::lin2d::GpLin2d;
use occt_core::gp::pnt2d::GpPnt2d;

use super::gcc_ent::{GccEntPosition, GccEntQualifiedCirc, GccEntQualifiedLin};

/// The largest `NCollection_Array1` upper bound of the twelve constructors
/// (`cirsol(1, 8)` of `_6.cxx:49`).
pub const MAX_SOL: usize = 8;

/// `GccAna_Circ2d2TanOn` (`GccAna_Circ2d2TanOn.hxx:45-260`).
#[derive(Debug, Clone, Copy)]
pub struct GccAnaCirc2d2TanOn {
    pub(crate) well_done: bool,
    pub(crate) nbr_sol: usize,
    pub(crate) cirsol: [GpCirc2d; MAX_SOL],
    pub(crate) qualifier1: [GccEntPosition; MAX_SOL],
    pub(crate) qualifier2: [GccEntPosition; MAX_SOL],
    pub(crate) the_same1: [i32; MAX_SOL],
    pub(crate) the_same2: [i32; MAX_SOL],
    pub(crate) pnttg1sol: [GpPnt2d; MAX_SOL],
    pub(crate) pnttg2sol: [GpPnt2d; MAX_SOL],
    pub(crate) pntcen: [GpPnt2d; MAX_SOL],
    pub(crate) par1sol: [f64; MAX_SOL],
    pub(crate) par2sol: [f64; MAX_SOL],
    pub(crate) pararg1: [f64; MAX_SOL],
    pub(crate) pararg2: [f64; MAX_SOL],
    pub(crate) parcen3: [f64; MAX_SOL],
}

impl GccAnaCirc2d2TanOn {
    /// The `: cirsol(1, ...), ... :` initializer list every constructor shares,
    /// with `TheSame1`/`TheSame2` at zero, `WellDone = false` and `NbrSol = 0`.
    pub(crate) fn empty() -> Self {
        Self {
            well_done: false,
            nbr_sol: 0,
            cirsol: [dummy_circle(); MAX_SOL],
            qualifier1: [GccEntPosition::Unqualified; MAX_SOL],
            qualifier2: [GccEntPosition::Unqualified; MAX_SOL],
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
        }
    }

    /// `IsDone()` (`cxx:383-386`).
    pub fn is_done(&self) -> bool {
        self.well_done
    }

    /// `NbSolutions()` (`cxx:388-391`).
    pub fn nb_solutions(&self) -> usize {
        self.nbr_sol
    }

    /// `ThisSolution(Index)` (`cxx:393-396`).
    ///
    /// OCCT indexes `cirsol` without the `WellDone`/range guards of the other
    /// accessors (`cxx:393-396` has neither throw); the port keeps the access
    /// unguarded to preserve the `Standard_OutOfRange`-free behaviour.
    pub fn this_solution(&self, index: usize) -> GpCirc2d {
        self.cirsol[index - 1]
    }

    /// `WhichQualifier(Index, Qualif1, Qualif2)` (`cxx:398-415`).
    pub fn which_qualifier(&self, index: usize) -> (GccEntPosition, GccEntPosition) {
        self.check_done();
        self.check_index(index);
        (self.qualifier1[index - 1], self.qualifier2[index - 1])
    }

    /// `Tangency1(Index, ParSol, ParArg, PntSol)` (`cxx:417-445`).
    ///
    /// `None` where OCCT raises `StdFail_NotDone` for `TheSame1(Index) != 0`.
    pub fn tangency1(&self, index: usize) -> Option<(f64, f64, GpPnt2d)> {
        self.check_done();
        self.check_index(index);
        if self.the_same1[index - 1] != 0 {
            return None;
        }
        Some((
            self.par1sol[index - 1],
            self.pararg1[index - 1],
            self.pnttg1sol[index - 1],
        ))
    }

    /// `Tangency2(Index, ParSol, ParArg, PntSol)` (`cxx:447-472`).
    ///
    /// `None` where OCCT raises `StdFail_NotDone` for `TheSame2(Index) != 0`.
    pub fn tangency2(&self, index: usize) -> Option<(f64, f64, GpPnt2d)> {
        self.check_done();
        self.check_index(index);
        if self.the_same2[index - 1] != 0 {
            return None;
        }
        Some((
            self.par2sol[index - 1],
            self.pararg2[index - 1],
            self.pnttg2sol[index - 1],
        ))
    }

    /// `CenterOn3(Index, ParArg, PntArg)` (`cxx:474-488`).
    ///
    /// The point OCCT hands back is `pnttg1sol(Index)` (`cxx:484-485`), not the
    /// stored centre `pntcen(Index)`; the port returns the same field.
    pub fn center_on3(&self, index: usize) -> (f64, GpPnt2d) {
        self.check_done();
        self.check_index(index);
        (self.parcen3[index - 1], self.pnttg1sol[index - 1])
    }

    /// `IsTheSame1(Index)` (`cxx:490-507`).
    pub fn is_the_same1(&self, index: usize) -> bool {
        self.check_done();
        self.check_index(index);
        self.the_same1[index - 1] != 0
    }

    /// `IsTheSame2(Index)` (`cxx:509-523`).
    pub fn is_the_same2(&self, index: usize) -> bool {
        self.check_done();
        self.check_index(index);
        self.the_same2[index - 1] != 0
    }

    /// `if (!WellDone) throw StdFail_NotDone()`.
    fn check_done(&self) {
        assert!(self.well_done, "StdFail_NotDone in GccAna_Circ2d2TanOn");
    }

    /// `if (Index <= 0 || Index > NbrSol) throw Standard_OutOfRange()`.
    fn check_index(&self, index: usize) {
        assert!(
            index >= 1 && index <= self.nbr_sol,
            "Standard_OutOfRange in GccAna_Circ2d2TanOn"
        );
    }

    /// The twelve OCCT constructors (`GccAna_Circ2d2TanOn.hxx:51-144`).
    pub fn new_circ_circ_on_line(
        qualified1: &GccEntQualifiedCirc,
        qualified2: &GccEntQualifiedCirc,
        on_line: &GpLin2d,
        tolerance: f64,
    ) -> Self {
        circ_circ_on_line::new(qualified1, qualified2, on_line, tolerance)
    }

    /// `GccAna_Circ2d2TanOn(QualifiedCirc, QualifiedLin, gp_Lin2d, double)`
    /// (`_1.cxx:47-345`).
    pub fn new_circ_lin_on_line(
        qualified1: &GccEntQualifiedCirc,
        qualified2: &GccEntQualifiedLin,
        on_line: &GpLin2d,
        tolerance: f64,
    ) -> Self {
        circ_lin_on_line::new(qualified1, qualified2, on_line, tolerance)
    }

    /// `GccAna_Circ2d2TanOn(QualifiedLin, QualifiedLin, gp_Lin2d, double)`
    /// (`_2.cxx:30-189`).
    pub fn new_lin_lin_on_line(
        qualified1: &GccEntQualifiedLin,
        qualified2: &GccEntQualifiedLin,
        on_line: &GpLin2d,
        tolerance: f64,
    ) -> Self {
        lin_lin_on_line::new(qualified1, qualified2, on_line, tolerance)
    }

    /// `GccAna_Circ2d2TanOn(QualifiedCirc, gp_Pnt2d, gp_Lin2d, double)`
    /// (`_3.cxx:48-276`).
    pub fn new_circ_pnt_on_line(
        qualified1: &GccEntQualifiedCirc,
        point2: &GpPnt2d,
        on_line: &GpLin2d,
        tolerance: f64,
    ) -> Self {
        circ_pnt_on_line::new(qualified1, point2, on_line, tolerance)
    }

    /// `GccAna_Circ2d2TanOn(QualifiedLin, gp_Pnt2d, gp_Lin2d, double)`
    /// (`_4.cxx:51-241`).
    pub fn new_lin_pnt_on_line(
        qualified1: &GccEntQualifiedLin,
        point2: &GpPnt2d,
        on_line: &GpLin2d,
        tolerance: f64,
    ) -> Self {
        lin_pnt_on_line::new(qualified1, point2, on_line, tolerance)
    }

    /// `GccAna_Circ2d2TanOn(gp_Pnt2d, gp_Pnt2d, gp_Lin2d, double)`
    /// (`_5.cxx:38-99`).
    pub fn new_pnt_pnt_on_line(
        point1: &GpPnt2d,
        point2: &GpPnt2d,
        on_line: &GpLin2d,
        tolerance: f64,
    ) -> Self {
        pnt_pnt_on_line::new(point1, point2, on_line, tolerance)
    }

    /// `GccAna_Circ2d2TanOn(QualifiedCirc, QualifiedCirc, gp_Circ2d, double)`
    /// (`_6.cxx:45-401`).
    pub fn new_circ_circ_on_circ(
        qualified1: &GccEntQualifiedCirc,
        qualified2: &GccEntQualifiedCirc,
        on_circ: &GpCirc2d,
        tolerance: f64,
    ) -> Self {
        circ_circ_on_circ::new(qualified1, qualified2, on_circ, tolerance)
    }

    /// `GccAna_Circ2d2TanOn(QualifiedCirc, QualifiedLin, gp_Circ2d, double)`
    /// (`_7.cxx:34-327`).
    pub fn new_circ_lin_on_circ(
        qualified1: &GccEntQualifiedCirc,
        qualified2: &GccEntQualifiedLin,
        on_circ: &GpCirc2d,
        tolerance: f64,
    ) -> Self {
        circ_lin_on_circ::new(qualified1, qualified2, on_circ, tolerance)
    }

    /// `GccAna_Circ2d2TanOn(QualifiedCirc, gp_Pnt2d, gp_Circ2d, double)`
    /// (`_8.cxx:42-315`).
    pub fn new_circ_pnt_on_circ(
        qualified1: &GccEntQualifiedCirc,
        point2: &GpPnt2d,
        on_circ: &GpCirc2d,
        tolerance: f64,
    ) -> Self {
        circ_pnt_on_circ::new(qualified1, point2, on_circ, tolerance)
    }

    /// `GccAna_Circ2d2TanOn(QualifiedLin, QualifiedLin, gp_Circ2d, double)`
    /// (`_9.cxx:31-379`; the fourth argument is unnamed and unused).
    pub fn new_lin_lin_on_circ(
        qualified1: &GccEntQualifiedLin,
        qualified2: &GccEntQualifiedLin,
        on_circ: &GpCirc2d,
    ) -> Self {
        lin_lin_on_circ::new(qualified1, qualified2, on_circ)
    }

    /// `GccAna_Circ2d2TanOn(QualifiedLin, gp_Pnt2d, gp_Circ2d, double)`
    /// (`_10.cxx:47-235`).
    pub fn new_lin_pnt_on_circ(
        qualified1: &GccEntQualifiedLin,
        point2: &GpPnt2d,
        on_circ: &GpCirc2d,
        tolerance: f64,
    ) -> Self {
        lin_pnt_on_circ::new(qualified1, point2, on_circ, tolerance)
    }

    /// `GccAna_Circ2d2TanOn(gp_Pnt2d, gp_Pnt2d, gp_Circ2d, double)`
    /// (`_11.cxx:47-123`).
    pub fn new_pnt_pnt_on_circ(
        point1: &GpPnt2d,
        point2: &GpPnt2d,
        on_circ: &GpCirc2d,
        tolerance: f64,
    ) -> Self {
        pnt_pnt_on_circ::new(point1, point2, on_circ, tolerance)
    }
}

/// `gp_Circ2d(gp_Ax2d(gp_Pnt2d(0, 0), gp::DX2d()), 0.0)`: only ever read back
/// after a slot has been written by a constructor.
fn dummy_circle() -> GpCirc2d {
    GpCirc2d::new(axis_at(GpPnt2d::zero()), 0.0)
}

/// `gp_Dir2d dirx(gp_Dir2d::D::X)`: the `dirx` of every constructor.
pub(crate) fn dirx() -> GpDir2d {
    GpDir2d::new(1.0, 0.0).expect("gp::DX2d()")
}

/// `gp_Ax2d(Center, dirx)`.
pub(crate) fn axis_at(center: GpPnt2d) -> GpAx22d {
    GpAx22d::from_xdir(center, dirx())
}

/// `gp_Dir2d(p2.XY() - p1.XY())`, panicking where OCCT would raise
/// `Standard_ConstructionError` on coincident points.
pub(crate) fn dir_between(from: &GpPnt2d, to: &GpPnt2d) -> GpDir2d {
    GpDir2d::new(to.x() - from.x(), to.y() - from.y())
        .expect("gp_Dir2d from two distinct points")
}

/// `gp_Dir2d(-d.Y(), d.X())`: the left normal every constructor builds by hand.
pub(crate) fn left_normal(d: &GpDir2d) -> GpDir2d {
    GpDir2d::new(-d.y, d.x).expect("gp_Dir2d of the left normal")
}

/// `gp_Pnt2d(p.XY() + k * dir.XY())`.
pub(crate) fn offset_xy(p: &GpPnt2d, k: f64, dir: &GpDir2d) -> GpPnt2d {
    GpPnt2d::new(p.x() + k * dir.x, p.y() + k * dir.y)
}

/// The side test every constructor spells out by hand:
/// `((o.X() - c.X()) * n.X()) + ((o.Y() - c.Y()) * n.Y())`.
pub(crate) fn side(origin: &GpPnt2d, center: &GpPnt2d, normal: &GpDir2d) -> f64 {
    ((origin.x() - center.x()) * normal.x) + ((origin.y() - center.y()) * normal.y)
}

/// The three-branch circle qualifier every constructor repeats after it stores
/// a solution:
///
/// ```text
/// if (!Qualified.IsUnqualified())     qualifier = Qualified.Qualifier();
/// else if (|dist + r - R| < Tol)      qualifier = GccEnt_enclosed;
/// else if (|dist - R - r| < Tol)      qualifier = GccEnt_outside;
/// else                                qualifier = GccEnt_enclosing;
/// ```
pub(crate) fn qualifier_of_circle(
    qualified: &GccEntQualifiedCirc,
    dist: f64,
    radius: f64,
    r_arg: f64,
    tol: f64,
) -> GccEntPosition {
    if !qualified.is_unqualified() {
        qualified.qualifier()
    } else if (dist + radius - r_arg).abs() < tol {
        GccEntPosition::Enclosed
    } else if (dist - r_arg - radius).abs() < tol {
        GccEntPosition::Outside
    } else {
        GccEntPosition::Enclosing
    }
}

/// `if (!Qualified.IsUnqualified()) qualifier = Qualified.Qualifier();
/// else if (dc.Dot(normL) > 0.0) qualifier = GccEnt_outside;
/// else qualifier = GccEnt_enclosed;`
pub(crate) fn qualifier_of_line(
    qualified: &GccEntQualifiedLin,
    dc_dot_norm: f64,
) -> GccEntPosition {
    if !qualified.is_unqualified() {
        qualified.qualifier()
    } else if dc_dot_norm > 0.0 {
        GccEntPosition::Outside
    } else {
        GccEntPosition::Enclosed
    }
}
