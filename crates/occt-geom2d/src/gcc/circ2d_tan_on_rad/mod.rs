//! `GccAna_Circ2dTanOnRad`: circles of a given radius tangent to one argument
//! and centred on a second one.
//!
//! | OCCT | argument list | lines | port |
//! |---|---|---|---|
//! | `GccAna_Circ2dTanOnRad.cxx` | circ, on line | `cxx:55-322` | `circ_on_line.rs` |
//! | `GccAna_Circ2dTanOnRad_1.cxx` | lin, on line | `cxx:42-153` | `lin_on_line.rs` |
//! | `GccAna_Circ2dTanOnRad_2.cxx` | pnt, on line | `cxx:43-179` | `pnt_on_line.rs` |
//! | `GccAna_Circ2dTanOnRad_3.cxx` | circ, on circ | `cxx:40-317` | `circ_on_circ.rs` |
//! | `GccAna_Circ2dTanOnRad_4.cxx` | lin, on circ | `cxx:38-175` | `lin_on_circ.rs` |
//! | `GccAna_Circ2dTanOnRad_5.cxx` | pnt, on circ | `cxx:37-124` | `pnt_on_circ.rs` |
//! | `IsDone`..`IsTheSame1` | | `cxx:320-412` | this file |
//!
//! Translation notes:
//! - OCCT sizes the eight `NCollection_Array1` fields per constructor: `(1, 4)`
//!   in circ/on line and circ/on circ, `(1, 2)` in the four line/point
//!   constructors. The port sizes every field at [`MAX_SOL`] (= 4, the largest)
//!   so the accessors do not depend on which constructor ran; each constructor
//!   still keeps its own OCCT upper bound where the source clamps,
//! - the local 1-based arrays of the constructors (`Coef(1, 2)`, `eps(1, 2)`,
//!   `Center(1, 4)`, `dir1on(1, 4)`) become fixed-size Rust arrays indexed by
//!   `index - 1`; the `gp_Dir2d` slots OCCT leaves uninitialised are seeded
//!   with `gp::DX2d()` and only read after a write, exactly like the OCCT
//!   slots they mirror,
//! - `Tol = std::abs(Tolerance)` (`cxx:73`) is the tolerance of every test in
//!   this class; the raw `Tolerance` argument is never used,
//! - `GccEnt_BadQualifier`, `Standard_NegativeValue`, `StdFail_NotDone` and
//!   `Standard_OutOfRange` panic with the OCCT exception name, matching
//!   [`super::circ2d_2tan_on`],
//! - `CenterOn3` hands back `pnttg1sol(Index)` (`cxx:391`), not the stored
//!   centre `pntcen3(Index)`; the port keeps that,
//! - `Tangency1` (`cxx:358-377`) has **no** `TheSame1` guard in this class,
//!   unlike `GccAna_Circ2d2TanOn::Tangency1`; the port returns the stored
//!   parameters unconditionally.
//!
//! UNPORTED: none - all six constructors and all seven accessors.

mod circ_on_circ;
mod circ_on_line;
mod lin_on_circ;
mod lin_on_line;
mod pnt_on_circ;
mod pnt_on_line;

use occt_core::gp::ax22d::GpAx22d;
use occt_core::gp::circ2d::GpCirc2d;
use occt_core::gp::pnt2d::GpPnt2d;

use super::gcc_ent::{GccEntPosition, GccEntQualifiedCirc, GccEntQualifiedLin};

/// The largest `NCollection_Array1` upper bound of the six constructors
/// (`cirsol(1, 4)` of `GccAna_Circ2dTanOnRad.cxx:60`).
pub const MAX_SOL: usize = 4;

/// `GccAna_Circ2dTanOnRad` (`GccAna_Circ2dTanOnRad.hxx:57-206`).
#[derive(Debug, Clone, Copy)]
pub struct GccAnaCirc2dTanOnRad {
    pub(crate) well_done: bool,
    pub(crate) nbr_sol: usize,
    pub(crate) cirsol: [GpCirc2d; MAX_SOL],
    pub(crate) qualifier1: [GccEntPosition; MAX_SOL],
    pub(crate) the_same1: [i32; MAX_SOL],
    pub(crate) pnttg1sol: [GpPnt2d; MAX_SOL],
    pub(crate) pntcen3: [GpPnt2d; MAX_SOL],
    pub(crate) par1sol: [f64; MAX_SOL],
    pub(crate) pararg1: [f64; MAX_SOL],
    pub(crate) parcen3: [f64; MAX_SOL],
}

impl GccAnaCirc2dTanOnRad {
    /// The `: cirsol(1, ...), ... :` initializer list every constructor shares,
    /// with `WellDone = false` and `NbrSol = 0` (`cxx:55-68`).
    pub(crate) fn empty() -> Self {
        Self {
            well_done: false,
            nbr_sol: 0,
            cirsol: [dummy_circle(); MAX_SOL],
            qualifier1: [GccEntPosition::Unqualified; MAX_SOL],
            the_same1: [0; MAX_SOL],
            pnttg1sol: [GpPnt2d::zero(); MAX_SOL],
            pntcen3: [GpPnt2d::zero(); MAX_SOL],
            par1sol: [0.0; MAX_SOL],
            pararg1: [0.0; MAX_SOL],
            parcen3: [0.0; MAX_SOL],
        }
    }

    /// `IsDone()` (`cxx:320-323`).
    pub fn is_done(&self) -> bool {
        self.well_done
    }

    /// `NbSolutions()` (`cxx:325-328`).
    pub fn nb_solutions(&self) -> usize {
        self.nbr_sol
    }

    /// `ThisSolution(Index)` (`cxx:330-341`).
    ///
    /// OCCT guards the index but not `WellDone` here; the port keeps that.
    pub fn this_solution(&self, index: usize) -> GpCirc2d {
        assert!(
            index >= 1 && index <= self.nbr_sol,
            "Standard_OutOfRange in GccAna_Circ2dTanOnRad"
        );
        self.cirsol[index - 1]
    }

    /// `WhichQualifier(Index, Qualif1)` (`cxx:343-356`).
    pub fn which_qualifier(&self, index: usize) -> GccEntPosition {
        self.check_done();
        self.check_index(index);
        self.qualifier1[index - 1]
    }

    /// `Tangency1(Index, ParSol, ParArg, PntSol)` (`cxx:358-377`).
    ///
    /// `(ParSol, ParArg, PntSol)`; OCCT has no `TheSame1` guard here, so a
    /// solution stored as "the same as the argument" still reports the fields.
    pub fn tangency1(&self, index: usize) -> (f64, f64, GpPnt2d) {
        self.check_done();
        self.check_index(index);
        (
            self.par1sol[index - 1],
            self.pararg1[index - 1],
            self.pnttg1sol[index - 1],
        )
    }

    /// `CenterOn3(Index, ParArg, PntSol)` (`cxx:380-393`).
    ///
    /// The point OCCT hands back is `pnttg1sol(Index)` (`cxx:391`), not the
    /// stored centre `pntcen3(Index)`; the port returns the same field.
    pub fn center_on3(&self, index: usize) -> (f64, GpPnt2d) {
        self.check_done();
        self.check_index(index);
        (self.parcen3[index - 1], self.pnttg1sol[index - 1])
    }

    /// `IsTheSame1(Index)` (`cxx:395-412`).
    pub fn is_the_same1(&self, index: usize) -> bool {
        self.check_done();
        self.check_index(index);
        self.the_same1[index - 1] != 0
    }

    /// `if (!WellDone) throw StdFail_NotDone()`.
    fn check_done(&self) {
        assert!(self.well_done, "StdFail_NotDone in GccAna_Circ2dTanOnRad");
    }

    /// `if (Index <= 0 || Index > NbrSol) throw Standard_OutOfRange()`.
    fn check_index(&self, index: usize) {
        assert!(
            index >= 1 && index <= self.nbr_sol,
            "Standard_OutOfRange in GccAna_Circ2dTanOnRad"
        );
    }

    /// The six OCCT constructors (`GccAna_Circ2dTanOnRad.hxx:69-121`).
    pub fn new_circ_on_line(
        qualified1: &GccEntQualifiedCirc,
        on_line: &occt_core::gp::lin2d::GpLin2d,
        radius: f64,
        tolerance: f64,
    ) -> Self {
        circ_on_line::new(qualified1, on_line, radius, tolerance)
    }

    /// `GccAna_Circ2dTanOnRad(QualifiedLin, gp_Lin2d, double, double)`
    /// (`_1.cxx:42-153`).
    pub fn new_lin_on_line(
        qualified1: &GccEntQualifiedLin,
        on_line: &occt_core::gp::lin2d::GpLin2d,
        radius: f64,
        tolerance: f64,
    ) -> Self {
        lin_on_line::new(qualified1, on_line, radius, tolerance)
    }

    /// `GccAna_Circ2dTanOnRad(gp_Pnt2d, gp_Lin2d, double, double)`
    /// (`_2.cxx:43-179`).
    pub fn new_pnt_on_line(
        point1: &GpPnt2d,
        on_line: &occt_core::gp::lin2d::GpLin2d,
        radius: f64,
        tolerance: f64,
    ) -> Self {
        pnt_on_line::new(point1, on_line, radius, tolerance)
    }

    /// `GccAna_Circ2dTanOnRad(QualifiedCirc, gp_Circ2d, double, double)`
    /// (`_3.cxx:40-317`).
    pub fn new_circ_on_circ(
        qualified1: &GccEntQualifiedCirc,
        on_circ: &GpCirc2d,
        radius: f64,
        tolerance: f64,
    ) -> Self {
        circ_on_circ::new(qualified1, on_circ, radius, tolerance)
    }

    /// `GccAna_Circ2dTanOnRad(QualifiedLin, gp_Circ2d, double, double)`
    /// (`_4.cxx:38-175`).
    pub fn new_lin_on_circ(
        qualified1: &GccEntQualifiedLin,
        on_circ: &GpCirc2d,
        radius: f64,
        tolerance: f64,
    ) -> Self {
        lin_on_circ::new(qualified1, on_circ, radius, tolerance)
    }

    /// `GccAna_Circ2dTanOnRad(gp_Pnt2d, gp_Circ2d, double, double)`
    /// (`_5.cxx:37-124`).
    pub fn new_pnt_on_circ(
        point1: &GpPnt2d,
        on_circ: &GpCirc2d,
        radius: f64,
        tolerance: f64,
    ) -> Self {
        pnt_on_circ::new(point1, on_circ, radius, tolerance)
    }
}

/// `gp_Circ2d(gp_Ax2d(gp_Pnt2d(0, 0), gp::DX2d()), 0.0)`: only ever read back
/// after a slot has been written by a constructor.
fn dummy_circle() -> GpCirc2d {
    GpCirc2d::new(GpAx22d::from_xdir(GpPnt2d::zero(), super::circ2d_2tan_on::dirx()), 0.0)
}
