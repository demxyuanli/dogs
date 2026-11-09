//! `GccAna_Circ2d2TanRad`: circles tangent to two arguments with a given
//! radius.
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | `GccAna_Circ2d2TanRad.cxx` (circle, circle) | `cxx:43-802` | `circ_circ.rs` |
//! | `GccAna_Circ2d2TanRad_1.cxx` (circle, line) | `cxx:42-626` | `circ_lin.rs` |
//! | `GccAna_Circ2d2TanRad_2.cxx` (circle, point) | `cxx:43-291` | `circ_pnt.rs` |
//! | `GccAna_Circ2d2TanRad_3.cxx` (line, point) | `cxx:46-271` | `lin_pnt.rs` |
//! | `GccAna_Circ2d2TanRad_4.cxx` (line, line) | `cxx:42-244` | `lin_lin.rs` |
//! | `GccAna_Circ2d2TanRad_5.cxx` (point, point) | `cxx:35-115` | `pnt_pnt.rs` |
//! | `IsDone`..`IsTheSame2` | `cxx:804-934` | this file |
//!
//! Translation notes:
//! - `NCollection_Array1<...>(1, 8)` / `(1, 4)` / `(1, 2)` become fixed arrays
//!   of [`MAX_SOL`] entries; OCCT indexes them from 1, the port from 0, and the
//!   public queries keep OCCT's 1-based `Index` and translate at the boundary,
//! - the `StdFail_NotDone`, `Standard_OutOfRange`, `Standard_NegativeValue` and
//!   `GccEnt_BadQualifier` raises have no Rust equivalent at the call sites of
//!   this package, so they panic with the OCCT exception name, which is how the
//!   port already models `StdFail_NotDone` in `IntAna2d_AnaIntersection`,
//! - the three `gp_Pnt2d` constructors (`_2.cxx`, `_3.cxx`, `_5.cxx`) declare
//!   `(1, 4)`, `(1, 2)` and `(1, 2)` upper bounds; the port sizes them at
//!   [`MAX_SOL`] like the rest and each constructor keeps its own OCCT clamp
//!   where the source has one (`_3.cxx:217` `NbrSol < 2`).

mod circ_circ;
mod circ_lin;
mod circ_pnt;
mod lin_lin;
mod lin_pnt;
mod pnt_pnt;

use occt_core::gp::ax22d::GpAx22d;
use occt_core::gp::circ2d::GpCirc2d;
use occt_core::gp::dir2d::GpDir2d;
use occt_core::gp::pnt2d::GpPnt2d;

use crate::gcc::gcc_ent::GccEntPosition;

/// `NCollection_Array1<...>(1, 8)`: the largest solution count any of the
/// ported constructors declares.
pub const MAX_SOL: usize = 8;

/// `GccAna_Circ2d2TanRad` (`GccAna_Circ2d2TanRad.hxx:39-190`).
///
/// The `qualifier1`, `qualifier2`, `TheSame1`, `TheSame2`, `cirsol`,
/// `pnttg1sol`, `pnttg2sol`, `par1sol`, `par2sol`, `pararg1` and `pararg2`
/// arrays are all `(1, 8)` in the circle/circle constructor, `(1, 8)` in the
/// circle/line one and `(1, 4)` in the line/line one; the port sizes them at
/// `MAX_SOL` for every constructor so the accessors do not depend on which
/// constructor ran.
pub struct GccAnaCirc2d2TanRad {
    pub(crate) qualifier1: [GccEntPosition; MAX_SOL],
    pub(crate) qualifier2: [GccEntPosition; MAX_SOL],
    pub(crate) the_same1: [i32; MAX_SOL],
    pub(crate) the_same2: [i32; MAX_SOL],
    pub(crate) cirsol: [GpCirc2d; MAX_SOL],
    pub(crate) pnttg1sol: [GpPnt2d; MAX_SOL],
    pub(crate) pnttg2sol: [GpPnt2d; MAX_SOL],
    pub(crate) par1sol: [f64; MAX_SOL],
    pub(crate) par2sol: [f64; MAX_SOL],
    pub(crate) pararg1: [f64; MAX_SOL],
    pub(crate) pararg2: [f64; MAX_SOL],
    pub(crate) well_done: bool,
    pub(crate) nbr_sol: usize,
}

impl GccAnaCirc2d2TanRad {
    /// The `: qualifier1(1, 8), ...` initializer list every constructor shares
    /// (`GccAna_Circ2d2TanRad.cxx:46-57`, `_1.cxx:46-57`, `_4.cxx:45-56`).
    pub(crate) fn empty() -> Self {
        let dummy_circ = GpCirc2d::new(GpAx22d::default(), 0.0);
        Self {
            qualifier1: [GccEntPosition::Unqualified; MAX_SOL],
            qualifier2: [GccEntPosition::Unqualified; MAX_SOL],
            the_same1: [0; MAX_SOL],
            the_same2: [0; MAX_SOL],
            cirsol: [dummy_circ; MAX_SOL],
            pnttg1sol: [GpPnt2d::zero(); MAX_SOL],
            pnttg2sol: [GpPnt2d::zero(); MAX_SOL],
            par1sol: [0.0; MAX_SOL],
            par2sol: [0.0; MAX_SOL],
            pararg1: [0.0; MAX_SOL],
            pararg2: [0.0; MAX_SOL],
            well_done: false,
            nbr_sol: 0,
        }
    }

    /// `IsDone()` (`cxx:804-807`).
    pub fn is_done(&self) -> bool {
        self.well_done
    }

    /// `NbSolutions()` (`cxx:809-812`).
    pub fn nb_solutions(&self) -> usize {
        self.nbr_sol
    }

    /// `ThisSolution(Index)` (`cxx:814-826`).
    pub fn this_solution(&self, index: usize) -> GpCirc2d {
        self.check_done();
        self.check_index(index);
        self.cirsol[index - 1]
    }

    /// `WhichQualifier(Index, Qualif1, Qualif2)` (`cxx:828-841`).
    pub fn which_qualifier(&self, index: usize) -> (GccEntPosition, GccEntPosition) {
        self.check_done();
        self.check_index(index);
        (self.qualifier1[index - 1], self.qualifier2[index - 1])
    }

    /// `Tangency1(Index, ParSol, ParArg, PntSol)` (`cxx:843-870`).
    ///
    /// Returns `None` where OCCT raises `StdFail_NotDone` for
    /// `TheSame1(Index) != 0`.
    pub fn tangency1(&self, index: usize) -> Option<(f64, f64, GpPnt2d)> {
        self.check_done();
        self.check_index(index);
        if self.the_same1[index - 1] != 0 {
            return None;
        }
        Some((self.par1sol[index - 1], self.pararg1[index - 1], self.pnttg1sol[index - 1]))
    }

    /// `Tangency2(Index, ParSol, ParArg, PntSol)` (`cxx:872-898`).
    ///
    /// Returns `None` where OCCT raises `StdFail_NotDone` for
    /// `TheSame2(Index) != 0`.
    pub fn tangency2(&self, index: usize) -> Option<(f64, f64, GpPnt2d)> {
        self.check_done();
        self.check_index(index);
        if self.the_same2[index - 1] != 0 {
            return None;
        }
        Some((self.par2sol[index - 1], self.pararg2[index - 1], self.pnttg2sol[index - 1]))
    }

    /// `IsTheSame1(Index)` (`cxx:900-916`).
    pub fn is_the_same1(&self, index: usize) -> bool {
        self.check_done();
        self.check_index(index);
        self.the_same1[index - 1] != 0
    }

    /// `IsTheSame2(Index)` (`cxx:918-934`).
    pub fn is_the_same2(&self, index: usize) -> bool {
        self.check_done();
        self.check_index(index);
        self.the_same2[index - 1] != 0
    }

    /// `if (!WellDone) throw StdFail_NotDone()`.
    fn check_done(&self) {
        assert!(self.well_done, "StdFail_NotDone in GccAna_Circ2d2TanRad");
    }

    /// `if (Index <= 0 || Index > NbrSol) throw Standard_OutOfRange()`.
    fn check_index(&self, index: usize) {
        assert!(
            index >= 1 && index <= self.nbr_sol,
            "Standard_OutOfRange in GccAna_Circ2d2TanRad"
        );
    }

    /// `ElCLib::Parameter(cirsol(i), pnttg1sol(i))` / `...2...` of the closing
    /// loop (`GccAna_Circ2d2TanRad.cxx:790-800`, `_1.cxx:616-625`,
    /// `_4.cxx:237-243`).
    ///
    /// The line variant of `ElCLib::Parameter` is
    /// [`occt_core::elib::clib2d::parameter_lin2d`] and the circle variant is
    /// [`occt_core::elib::clib2d::parameter_circ2d`]; the solver always builds
    /// `cirsol` on the `gp::DX2d()` axis, so `par1sol`/`par2sol` both use the
    /// circle overload.
    pub(crate) fn fill_solution_parameters(&mut self) {
        use occt_core::elib::clib2d::parameter_circ2d;
        for i in 0..self.nbr_sol {
            self.par1sol[i] = parameter_circ2d(&self.cirsol[i], &self.pnttg1sol[i]);
            self.par2sol[i] = parameter_circ2d(&self.cirsol[i], &self.pnttg2sol[i]);
        }
    }
}

/// `gp::DX2d()`: the `gp_Dir2d(1, 0)` every constructor stores in `dirx`.
pub(crate) fn dx2d() -> GpDir2d {
    GpDir2d::new(1.0, 0.0).expect("gp::DX2d()")
}

/// `gp_Ax2d(center, dirx)`: the axis every solution circle is built on.
pub(crate) fn axis_at(center: GpPnt2d) -> GpAx22d {
    GpAx22d::from_xdir(center, dx2d())
}

/// `gp_Circ2d(C1.XAxis(), r)` (`GccAna_Circ2d2TanRad_2.cxx:113`, `:136`, `:155`,
/// `:194`, `:200`, `:203`): the auxiliary circle keeps `C1`'s centre and X
/// direction and only replaces the radius. Distinct from [`axis_at`], which
/// always uses `gp::DX2d()`.
pub(crate) fn circ_on_x_axis(c1: &GpCirc2d, radius: f64) -> GpCirc2d {
    let ax = c1.x_axis();
    GpCirc2d::new(GpAx22d::from_xdir(ax.loc, ax.vdir), radius)
}

/// `gp_Pnt2d(center.XY() + k * dir.XY())`.
pub(crate) fn offset_by(center: &GpPnt2d, k: f64, dir: &GpDir2d) -> GpPnt2d {
    GpPnt2d::new(center.x() + k * dir.x, center.y() + k * dir.y)
}

/// `gp_Dir2d(p2.XY() - p1.XY())`.
pub(crate) fn direction_between(from: &GpPnt2d, to: &GpPnt2d) -> GpDir2d {
    GpDir2d::new(to.x() - from.x(), to.y() - from.y()).expect("gp_Dir2d from two distinct points")
}
