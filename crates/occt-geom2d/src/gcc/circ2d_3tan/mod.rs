//! `GccAna_Circ2d3Tan`: circles tangent to three points / lines / circles.
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | `GccAna_Circ2d3Tan.cxx` (circ, circ, circ) | `cxx:31-1017` | `three_circles.rs` |
//! | `GccAna_Circ2d3Tan_1.cxx` (circ, circ, lin) | `cxx:38-415` | `two_circles_line.rs` |
//! | `GccAna_Circ2d3Tan_2.cxx` (circ, lin, lin) | `cxx:38-323` | `circle_two_lines.rs` |
//! | `GccAna_Circ2d3Tan_3.cxx` (lin, lin, lin) | `cxx:38-547` | `three_lines.rs` |
//! | `GccAna_Circ2d3Tan_4.cxx` (circ, circ, pnt) | `cxx:41-465` | `two_circles_point.rs` |
//! | `GccAna_Circ2d3Tan_5.cxx` (circ, lin, pnt) | `cxx:38-309` | `circle_line_point.rs` |
//! | `GccAna_Circ2d3Tan_6.cxx` (lin, lin, pnt) | `cxx:34-251` | `two_lines_point.rs` |
//! | `GccAna_Circ2d3Tan_7.cxx` (circ, pnt, pnt) | `cxx:36-241` | `circle_two_points.rs` |
//! | `GccAna_Circ2d3Tan_8.cxx` (lin, pnt, pnt) | `cxx:34-177` | `line_two_points.rs` |
//! | `GccAna_Circ2d3Tan_9.cxx` (pnt, pnt, pnt) | `cxx:40-160` | `three_points.rs` |
//! | `IsDone`..`IsTheSame3` | `cxx:1020-1210` | this file |
//!
//! Translation notes:
//! - OCCT sizes the fourteen `NCollection_Array1` fields per constructor:
//!   `(1, 16)` in circ/circ/circ and circ/circ/lin, `(1, 8)` in circ/lin/lin,
//!   `(1, 4)` in lin/lin/lin and circ/lin/pnt, `(1, 20)` in circ/circ/pnt
//!   (`static int MaxSol = 20`, `_4.cxx:35`), `(1, 2)` in lin/lin/pnt,
//!   circ/pnt/pnt and lin/pnt/pnt, `(1, 1)` in pnt/pnt/pnt. The port sizes every
//!   field at [`MAX_SOL`] (= 20, the largest) so the accessors do not depend on
//!   which constructor ran; each constructor keeps its own OCCT upper bound as
//!   the `upper` local used by its `NbrSol == MaxSol` clamps,
//! - `NCollection_Array1` is 1-based and so is `FirstSol`/`FirstSol1` in the
//!   circ/circ/circ solver; the port keeps OCCT's 1-based loop variables and
//!   offsets only when indexing the Rust arrays,
//! - the `GccEnt_BadQualifier`, `StdFail_NotDone` and `Standard_OutOfRange`
//!   raises have no Rust equivalent at the call sites of this package and panic
//!   with the OCCT exception name, matching [`super::circ2d_2tan_rad`],
//! - `_3.cxx:519-521` appends filtered copies of already stored solutions while
//!   its arrays are `(1, 4)`; OCCT writes past the upper bound there. The port
//!   stops at [`MAX_SOL`] (a memory-safety guard only, no test-visible rule was
//!   added).
//!
//! UNPORTED: none - all ten constructors and all nine accessors.

mod circle_line_point;
mod circle_two_lines;
mod circle_two_points;
mod line_two_points;
mod three_circles;
mod three_lines;
mod three_points;
mod two_circles_line;
mod two_circles_point;
mod two_lines_point;

use occt_core::gp::ax22d::GpAx22d;
use occt_core::gp::circ2d::GpCirc2d;
use occt_core::gp::dir2d::GpDir2d;
use occt_core::gp::pnt2d::GpPnt2d;

use super::gcc_ent::{GccEntPosition, GccEntQualifiedCirc, GccEntQualifiedLin};

/// The largest `NCollection_Array1` upper bound of the ten constructors
/// (`MaxSol = 20` of `_4.cxx:35`).
pub const MAX_SOL: usize = 20;

/// `GccAna_Circ2d3Tan` (`GccAna_Circ2d3Tan.hxx:53-259`).
#[derive(Debug, Clone, Copy)]
pub struct GccAnaCirc2d3Tan {
    pub(crate) well_done: bool,
    pub(crate) nbr_sol: usize,
    pub(crate) cirsol: [GpCirc2d; MAX_SOL],
    pub(crate) qualifier1: [GccEntPosition; MAX_SOL],
    pub(crate) qualifier2: [GccEntPosition; MAX_SOL],
    pub(crate) qualifier3: [GccEntPosition; MAX_SOL],
    pub(crate) the_same1: [i32; MAX_SOL],
    pub(crate) the_same2: [i32; MAX_SOL],
    pub(crate) the_same3: [i32; MAX_SOL],
    pub(crate) pnttg1sol: [GpPnt2d; MAX_SOL],
    pub(crate) pnttg2sol: [GpPnt2d; MAX_SOL],
    pub(crate) pnttg3sol: [GpPnt2d; MAX_SOL],
    pub(crate) par1sol: [f64; MAX_SOL],
    pub(crate) par2sol: [f64; MAX_SOL],
    pub(crate) par3sol: [f64; MAX_SOL],
    pub(crate) pararg1: [f64; MAX_SOL],
    pub(crate) pararg2: [f64; MAX_SOL],
    pub(crate) pararg3: [f64; MAX_SOL],
}

impl GccAnaCirc2d3Tan {
    /// The `: cirsol(1, ...), ... :` initializer list every constructor shares,
    /// with `WellDone = false` and `NbrSol = 0`.
    pub(crate) fn empty() -> Self {
        Self {
            well_done: false,
            nbr_sol: 0,
            cirsol: [dummy_circle(); MAX_SOL],
            qualifier1: [GccEntPosition::Unqualified; MAX_SOL],
            qualifier2: [GccEntPosition::Unqualified; MAX_SOL],
            qualifier3: [GccEntPosition::Unqualified; MAX_SOL],
            the_same1: [0; MAX_SOL],
            the_same2: [0; MAX_SOL],
            the_same3: [0; MAX_SOL],
            pnttg1sol: [GpPnt2d::zero(); MAX_SOL],
            pnttg2sol: [GpPnt2d::zero(); MAX_SOL],
            pnttg3sol: [GpPnt2d::zero(); MAX_SOL],
            par1sol: [0.0; MAX_SOL],
            par2sol: [0.0; MAX_SOL],
            par3sol: [0.0; MAX_SOL],
            pararg1: [0.0; MAX_SOL],
            pararg2: [0.0; MAX_SOL],
            pararg3: [0.0; MAX_SOL],
        }
    }

    /// `IsDone()` (`cxx:1022-1025`).
    pub fn is_done(&self) -> bool {
        self.well_done
    }

    /// `NbSolutions()` (`cxx:1027-1030`).
    pub fn nb_solutions(&self) -> usize {
        self.nbr_sol
    }

    /// `ThisSolution(Index)` (`cxx:1032-1045`).
    pub fn this_solution(&self, index: usize) -> GpCirc2d {
        self.check_done();
        self.check_index(index);
        self.cirsol[index - 1]
    }

    /// `WhichQualifier(Index, Qualif1, Qualif2, Qualif3)` (`cxx:1047-1066`).
    pub fn which_qualifier(
        &self,
        index: usize,
    ) -> (GccEntPosition, GccEntPosition, GccEntPosition) {
        self.check_done();
        self.check_index(index);
        (
            self.qualifier1[index - 1],
            self.qualifier2[index - 1],
            self.qualifier3[index - 1],
        )
    }

    /// `Tangency1(Index, ParSol, ParArg, PntSol)` (`cxx:1068-1094`).
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

    /// `Tangency2(Index, ParSol, ParArg, PntSol)` (`cxx:1096-1122`).
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

    /// `Tangency3(Index, ParSol, ParArg, PntSol)` (`cxx:1124-1150`).
    pub fn tangency3(&self, index: usize) -> Option<(f64, f64, GpPnt2d)> {
        self.check_done();
        self.check_index(index);
        if self.the_same3[index - 1] != 0 {
            return None;
        }
        Some((
            self.par3sol[index - 1],
            self.pararg3[index - 1],
            self.pnttg3sol[index - 1],
        ))
    }

    /// `IsTheSame1(Index)` (`cxx:1152-1172`).
    pub fn is_the_same1(&self, index: usize) -> bool {
        self.check_done();
        self.check_index(index);
        self.the_same1[index - 1] != 0
    }

    /// `IsTheSame2(Index)` (`cxx:1174-1194`).
    pub fn is_the_same2(&self, index: usize) -> bool {
        self.check_done();
        self.check_index(index);
        self.the_same2[index - 1] != 0
    }

    /// `IsTheSame3(Index)` (`cxx:1196-1210`).
    pub fn is_the_same3(&self, index: usize) -> bool {
        self.check_done();
        self.check_index(index);
        self.the_same3[index - 1] != 0
    }

    /// `if (!WellDone) throw StdFail_NotDone()`.
    fn check_done(&self) {
        assert!(self.well_done, "StdFail_NotDone in GccAna_Circ2d3Tan");
    }

    /// `if (Index <= 0 || Index > NbrSol) throw Standard_OutOfRange()`.
    fn check_index(&self, index: usize) {
        assert!(
            index >= 1 && index <= self.nbr_sol,
            "Standard_OutOfRange in GccAna_Circ2d3Tan"
        );
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

/// The three-branch qualifier resolution every constructor repeats after it
/// stores a solution:
///
/// ```text
/// if (!Qualified.IsUnqualified())     qualifier = Qualified.Qualifier();
/// else if (|dist + r - R| <= Tol)     qualifier = GccEnt_enclosed;
/// else if (|dist - R - r| <= Tol)     qualifier = GccEnt_outside;
/// else                                qualifier = GccEnt_enclosing;
/// ```
///
/// (`GccAna_Circ2d3Tan.cxx:898-913` and the same block in every other file,
/// always with `dist = Center.Distance(that circle's centre)`).
pub(crate) fn qualifier_of_circle(
    qualified: &GccEntQualifiedCirc,
    dist: f64,
    radius: f64,
    r_arg: f64,
    tol: f64,
) -> GccEntPosition {
    if !qualified.is_unqualified() {
        qualified.qualifier()
    } else if (dist + radius - r_arg).abs() <= tol {
        GccEntPosition::Enclosed
    } else if (dist - r_arg - radius).abs() <= tol {
        GccEntPosition::Outside
    } else {
        GccEntPosition::Enclosing
    }
}

/// The ten OCCT constructors (`GccAna_Circ2d3Tan.hxx:64-215`), each named after
/// its argument list.
impl GccAnaCirc2d3Tan {
    /// `GccAna_Circ2d3Tan(Qualified1, Qualified2, Qualified3, Tolerance)`
    /// (`GccAna_Circ2d3Tan.cxx:31-1017`): three circles.
    pub fn new_circ_circ_circ(
        qualified1: &GccEntQualifiedCirc,
        qualified2: &GccEntQualifiedCirc,
        qualified3: &GccEntQualifiedCirc,
        tolerance: f64,
    ) -> Self {
        three_circles::new(qualified1, qualified2, qualified3, tolerance)
    }

    /// `GccAna_Circ2d3Tan(Qualified1, Qualified2, Qualified3, Tolerance)`
    /// (`GccAna_Circ2d3Tan_1.cxx:38-415`): two circles and a line.
    pub fn new_circ_circ_lin(
        qualified1: &GccEntQualifiedCirc,
        qualified2: &GccEntQualifiedCirc,
        qualified3: &GccEntQualifiedLin,
        tolerance: f64,
    ) -> Self {
        two_circles_line::new(qualified1, qualified2, qualified3, tolerance)
    }

    /// `GccAna_Circ2d3Tan(Qualified1, Qualified2, Qualified3, Tolerance)`
    /// (`GccAna_Circ2d3Tan_2.cxx:38-323`): one circle and two lines.
    pub fn new_circ_lin_lin(
        qualified1: &GccEntQualifiedCirc,
        qualified2: &GccEntQualifiedLin,
        qualified3: &GccEntQualifiedLin,
        tolerance: f64,
    ) -> Self {
        circle_two_lines::new(qualified1, qualified2, qualified3, tolerance)
    }

    /// `GccAna_Circ2d3Tan(Qualified1, Qualified2, Qualified3)`
    /// (`GccAna_Circ2d3Tan_3.cxx:38-547`): three lines.
    pub fn new_lin_lin_lin(
        qualified1: &GccEntQualifiedLin,
        qualified2: &GccEntQualifiedLin,
        qualified3: &GccEntQualifiedLin,
    ) -> Self {
        three_lines::new(qualified1, qualified2, qualified3)
    }

    /// `GccAna_Circ2d3Tan(Qualified1, Qualified2, Point3, Tolerance)`
    /// (`GccAna_Circ2d3Tan_4.cxx:41-465`): two circles and a point.
    pub fn new_circ_circ_pnt(
        qualified1: &GccEntQualifiedCirc,
        qualified2: &GccEntQualifiedCirc,
        point3: &GpPnt2d,
        tolerance: f64,
    ) -> Self {
        two_circles_point::new(qualified1, qualified2, point3, tolerance)
    }

    /// `GccAna_Circ2d3Tan(Qualified1, Qualified2, Point3, Tolerance)`
    /// (`GccAna_Circ2d3Tan_5.cxx:38-309`): one circle, one line and a point.
    pub fn new_circ_lin_pnt(
        qualified1: &GccEntQualifiedCirc,
        qualified2: &GccEntQualifiedLin,
        point3: &GpPnt2d,
        tolerance: f64,
    ) -> Self {
        circle_line_point::new(qualified1, qualified2, point3, tolerance)
    }

    /// `GccAna_Circ2d3Tan(Qualified1, Qualified2, Point3, Tolerance)`
    /// (`GccAna_Circ2d3Tan_6.cxx:34-251`): two lines and a point.
    pub fn new_lin_lin_pnt(
        qualified1: &GccEntQualifiedLin,
        qualified2: &GccEntQualifiedLin,
        point3: &GpPnt2d,
        tolerance: f64,
    ) -> Self {
        two_lines_point::new(qualified1, qualified2, point3, tolerance)
    }

    /// `GccAna_Circ2d3Tan(Qualified1, Point2, Point3, Tolerance)`
    /// (`GccAna_Circ2d3Tan_7.cxx:36-241`): one circle and two points.
    pub fn new_circ_pnt_pnt(
        qualified1: &GccEntQualifiedCirc,
        point2: &GpPnt2d,
        point3: &GpPnt2d,
        tolerance: f64,
    ) -> Self {
        circle_two_points::new(qualified1, point2, point3, tolerance)
    }

    /// `GccAna_Circ2d3Tan(Qualified1, Point2, Point3, Tolerance)`
    /// (`GccAna_Circ2d3Tan_8.cxx:34-177`): one line and two points.
    pub fn new_lin_pnt_pnt(
        qualified1: &GccEntQualifiedLin,
        point2: &GpPnt2d,
        point3: &GpPnt2d,
        tolerance: f64,
    ) -> Self {
        line_two_points::new(qualified1, point2, point3, tolerance)
    }

    /// `GccAna_Circ2d3Tan(Point1, Point2, Point3, Tolerance)`
    /// (`GccAna_Circ2d3Tan_9.cxx:40-160`): three points.
    pub fn new_pnt_pnt_pnt(
        point1: &GpPnt2d,
        point2: &GpPnt2d,
        point3: &GpPnt2d,
        tolerance: f64,
    ) -> Self {
        three_points::new(point1, point2, point3, tolerance)
    }
}

