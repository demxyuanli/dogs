//! `GccAna_Lin2dTanPar` (`GccAna_Lin2dTanPar.cxx:37-193`): lines parallel to
//! a reference line, either through a point or tangent to a circle.
//!
//! | OCCT | lines |
//! |---|---|
//! | `GccAna_Lin2dTanPar(gp_Pnt2d, gp_Lin2d)` | `cxx:37-53` |
//! | `GccAna_Lin2dTanPar(GccEnt_QualifiedCirc, gp_Lin2d)` | `cxx:66-129` |
//! | `IsDone`..`Tangency1` | `cxx:131-193` |
//!
//! Translation notes:
//! - the `NCollection_Array1<...>(1, 1)` / `(1, 2)` fields become fixed arrays
//!   of [`MAX_SOL`]; OCCT indexes them from 1, the port from 0, and the public
//!   queries keep OCCT's 1-based `Index` and translate at the boundary,
//! - `GccEnt_BadQualifier`, `StdFail_NotDone` and `Standard_OutOfRange` have no
//!   Rust equivalent at the call sites of this package, so they panic with the
//!   OCCT exception name, matching [`super::circ2d_2tan_rad`].

use occt_core::elib::clib2d::parameter_circ2d;
use occt_core::gp::lin2d::GpLin2d;
use occt_core::gp::pnt2d::GpPnt2d;
use occt_core::gp::xy::GpXY;

use super::gcc_ent::{GccEntPosition, GccEntQualifiedCirc};

/// `NCollection_Array1<...>(1, 2)`: the largest solution count of the two
/// constructors.
pub const MAX_SOL: usize = 2;

/// `GccAna_Lin2dTanPar` (`GccAna_Lin2dTanPar.hxx`).
#[derive(Debug, Clone, Copy)]
pub struct GccAnaLin2dTanPar {
    pub(crate) linsol: [GpLin2d; MAX_SOL],
    pub(crate) qualifier1: [GccEntPosition; MAX_SOL],
    pub(crate) pnttg1sol: [GpPnt2d; MAX_SOL],
    pub(crate) par1sol: [f64; MAX_SOL],
    pub(crate) pararg1: [f64; MAX_SOL],
    pub(crate) well_done: bool,
    pub(crate) nbr_sol: usize,
}

impl GccAnaLin2dTanPar {
    /// The `: linsol(1, 1), qualifier1(1, 1), ...` initializer list the two
    /// constructors share (`cxx:38-42`, `cxx:67-71`).
    pub(crate) fn empty() -> Self {
        Self {
            linsol: [GpLin2d::default(); MAX_SOL],
            qualifier1: [GccEntPosition::Unqualified; MAX_SOL],
            pnttg1sol: [GpPnt2d::zero(); MAX_SOL],
            par1sol: [0.0; MAX_SOL],
            pararg1: [0.0; MAX_SOL],
            well_done: false,
            nbr_sol: 0,
        }
    }

    /// `GccAna_Lin2dTanPar(ThePoint, Lin1)` (`cxx:37-53`): the line through the
    /// point with the direction of `Lin1`.
    pub fn through_point(the_point: &GpPnt2d, lin1: &GpLin2d) -> Self {
        let mut s = Self::empty();
        s.linsol[0] = GpLin2d::from_pnt_dir(*the_point, *lin1.direction());
        // `qualifier1(1) = GccEnt_noqualifier` (`cxx:47`).
        s.qualifier1[0] = GccEntPosition::NoQualifier;
        s.pnttg1sol[0] = *the_point;
        s.par1sol[0] = 0.0;
        s.pararg1[0] = 0.0;
        s.nbr_sol = 1;
        s.well_done = true;
        s
    }

    /// `GccAna_Lin2dTanPar(Qualified1, Lin1)` (`cxx:66-129`): the lines parallel
    /// to `Lin1` and tangent to `Qualified1`.
    ///
    /// `xy = (-R1 * ydir, R1 * xdir)` is the radius vector turned by 90 degrees
    /// (`cxx:114`); the loop flips `signe` before every use (`cxx:117`), so the
    /// unqualified case emits the `outside` line first and the `enclosing` one
    /// second, exactly as the qualifier table at `cxx:107-113` declares.
    pub fn tangent_to_circle(qualified1: &GccEntQualifiedCirc, lin1: &GpLin2d) -> Self {
        let mut s = Self::empty();
        let nbsol: usize;
        let mut signe: i32;
        // `cxx:78-83`.
        if !(qualified1.is_enclosed()
            || qualified1.is_enclosing()
            || qualified1.is_outside()
            || qualified1.is_unqualified())
        {
            panic!("GccEnt_BadQualifier in GccAna_Lin2dTanPar");
        }
        let c1 = qualified1.qualified();
        let xdir = lin1.direction().x;
        let ydir = lin1.direction().y;

        // `cxx:88-113`.
        if qualified1.is_enclosed() {
            panic!("GccEnt_BadQualifier in GccAna_Lin2dTanPar");
        } else if qualified1.is_enclosing() {
            nbsol = 1;
            signe = 1;
            s.qualifier1[0] = GccEntPosition::Enclosing;
        } else if qualified1.is_outside() {
            nbsol = 1;
            signe = -1;
            s.qualifier1[0] = GccEntPosition::Outside;
        } else {
            nbsol = 2;
            signe = -1;
            s.qualifier1[0] = GccEntPosition::Outside;
            s.qualifier1[1] = GccEntPosition::Enclosing;
        }

        // `gp_XY xy(-C1.Radius() * ydir, C1.Radius() * xdir)` (`cxx:114`).
        let xy = GpXY::new(-c1.radius() * ydir, c1.radius() * xdir);
        for _ in 1..=nbsol {
            signe = -signe;
            s.nbr_sol += 1;
            let i = s.nbr_sol - 1;
            let pnt = GpPnt2d::from_xy(
                c1.location().xy().added(&xy.multiplied_scalar(signe as f64)),
            );
            s.linsol[i] = GpLin2d::from_pnt_dir(pnt, *lin1.direction());
            s.pnttg1sol[i] = pnt;
            s.par1sol[i] = 0.0;
            s.pararg1[i] = parameter_circ2d(&c1, &pnt);
            s.well_done = true;
        }
        s
    }

    /// `IsDone()` (`cxx:131-134`).
    pub fn is_done(&self) -> bool {
        self.well_done
    }

    /// `NbSolutions()` (`cxx:136-143`).
    pub fn nb_solutions(&self) -> usize {
        self.check_done();
        self.nbr_sol
    }

    /// `ThisSolution(Index)` (`cxx:145-156`).
    pub fn this_solution(&self, index: usize) -> GpLin2d {
        self.check_done();
        self.check_index(index);
        self.linsol[index - 1]
    }

    /// `WhichQualifier(Index, Qualif1)` (`cxx:158-172`).
    pub fn which_qualifier(&self, index: usize) -> GccEntPosition {
        self.check_done();
        self.check_index(index);
        self.qualifier1[index - 1]
    }

    /// `Tangency1(Index, ParSol, ParArg, Pnt)` (`cxx:174-193`).
    pub fn tangency1(&self, index: usize) -> (f64, f64, GpPnt2d) {
        self.check_done();
        self.check_index(index);
        (self.par1sol[index - 1], self.pararg1[index - 1], self.pnttg1sol[index - 1])
    }

    /// `if (!WellDone) throw StdFail_NotDone()`.
    fn check_done(&self) {
        assert!(self.well_done, "StdFail_NotDone in GccAna_Lin2dTanPar");
    }

    /// `if (Index <= 0 || Index > NbrSol) throw Standard_OutOfRange()`.
    fn check_index(&self, index: usize) {
        assert!(
            index >= 1 && index <= self.nbr_sol,
            "Standard_OutOfRange in GccAna_Lin2dTanPar"
        );
    }
}
