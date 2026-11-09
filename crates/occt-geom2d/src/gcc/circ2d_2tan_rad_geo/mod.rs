//! `Geom2dGcc_Circ2d2TanRadGeo`
//! (`Geom2dGcc_Circ2d2TanRadGeo.{hxx:56-181,cxx:51-1141}`): circles of a
//! prescribed radius tangent to two operands where at least one is a general
//! curve.
//!
//! The user gave the radius, so no bisector is built: every operand is offset
//! by `-+Radius` (`Adaptor2d_OffsetCurve`) and the centre is the intersection of
//! the two offsets, or of one offset with the carrier of the other operand.
//! Four constructors cover the operand pairs the front end can build.
//!
//! | OCCT constructor | arguments | `cxx` | file |
//! |---|---|---|---|
//! | 1 | `GccEnt_QualifiedLin, Geom2dGcc_QCurve` | `:51-249` | [`lin_qcurve`] |
//! | 2 | `GccEnt_QualifiedCirc, Geom2dGcc_QCurve` | `:265-474` | [`circ_qcurve`] |
//! | 3 | `Geom2dGcc_QCurve, gp_Pnt2d` | `:490-608` | [`qcurve_pnt`] |
//! | 4 | `Geom2dGcc_QCurve, Geom2dGcc_QCurve` | `:785-1006` | [`two_qcurves`] |
//!
//! Translation notes:
//! - `aNbSolMAX = 16` (`cxx:35`), twice the `Circ2dTanOnRadGeo` bound,
//! - the `cote1` / `cote2` cascade (`cxx:110-189`, `:321-397`, `:840-916`) is
//!   copied branch for branch in [`cotes`], including the fact that it is *not*
//!   one `if / else if` chain: four standalone tests follow the first chain,
//!   each guarded by conjunctions that keep them mutually exclusive, so the
//!   first matching branch wins as written,
//! - the inner sweep of the curve side is `for jcote2 = 1; jcote2 <= nbrcote2;
//!   jcote2++` with no `NbrSol` bound in constructors 2 and 4 (`cxx:411`,
//!   `:935`), while the point loop that follows is bounded (`cxx:429`,
//!   `:953`); constructor 1 bounds the outer sweep only (`cxx:191-197`),
//! - `par1sol` / `par2sol` are `ElCLib::Parameter` of the solution circle at the
//!   tangency points, evaluated after both tangency points are written
//!   (`cxx:238-242`), unlike the `Circ2dTanOnRadGeo` family,
//! - OCCT defect kept: the curve/curve constructor writes
//!   `qualifier1(NbrSol) = Qualified1.Qualifier();` twice (`cxx:987-988`) and
//!   never writes `qualifier2`, so [`Geom2dGccCirc2d2TanRadGeo::which_qualifier`]
//!   answers the untouched initial slot for the second operand,
//! - OCCT defect kept: `PrecRoot` is asked to refine in place
//!   (`PrecRoot(C1, C2, aU0, aV0, aU0, aV0)`, `cxx:980`), and the refined `aU0`
//!   is used for the centre (`cxx:984`) while `pararg1` / `pnttg1sol` still use
//!   the unrefined `Intp.Point(i)` parameters (`cxx:991-994`),
//! - the `#ifdef OCCT_DEBUG` blocks (`cxx:436-440`, `:812-817`, `:925-931`,
//!   `:939-948`) are not compiled in a release build and are not ported,
//! - the accessors guard `WellDone` first and the index second
//!   (`cxx:1018-1141`).
//!
//! UNPORTED: none - all four constructors, the eight accessors, the `PrecRoot`
//! static (`cxx:618-770`) and the field set.

use occt_core::gp::{GpCirc2d, GpPnt2d};

use crate::gcc::gcc_ent::{GccEntPosition, GccEntQualifiedCirc, GccEntQualifiedLin};
use crate::gcc::qualified_curve::Geom2dGccQCurve;

mod circ_qcurve;
mod lin_qcurve;
mod qcurve_pnt;
mod two_qcurves;

/// `aNbSolMAX = 16` (`cxx:35`).
pub(super) const MAX_SOL: usize = 16;

/// The qualifier categories the `cote` cascade and the bad-qualifier guard test.
#[derive(Clone, Copy, Debug, Default)]
pub(super) struct Categories {
    pub enclosed: bool,
    pub enclosing: bool,
    pub outside: bool,
    pub unqualified: bool,
}

impl Categories {
    /// The predicates of a line operand; a line has no `Enclosing`.
    pub(super) fn of_lin(q: &GccEntQualifiedLin) -> Self {
        Self {
            enclosed: q.is_enclosed(),
            enclosing: false,
            outside: q.is_outside(),
            unqualified: q.is_unqualified(),
        }
    }

    /// The predicates of a circle operand.
    pub(super) fn of_circ(q: &GccEntQualifiedCirc) -> Self {
        Self {
            enclosed: q.is_enclosed(),
            enclosing: q.is_enclosing(),
            outside: q.is_outside(),
            unqualified: q.is_unqualified(),
        }
    }

    /// The predicates of a general-curve operand.
    pub(super) fn of_qcurve(q: &Geom2dGccQCurve) -> Self {
        Self {
            enclosed: q.is_enclosed(),
            enclosing: q.is_enclosing(),
            outside: q.is_outside(),
            unqualified: q.is_unqualified(),
        }
    }

    /// `IsEnclosed() || IsEnclosing() || IsOutside() || IsUnqualified()`.
    pub(super) fn all_of(&self) -> bool {
        self.enclosed || self.enclosing || self.outside || self.unqualified
    }

    /// The same test without `IsEnclosing()`, used for line operands
    /// (`cxx:90-92`).
    pub(super) fn without_enclosing(&self) -> bool {
        self.enclosed || self.outside || self.unqualified
    }
}

/// The `cote1` / `cote2` cascade (`cxx:110-189`, `:321-397`, `:840-916`).
///
/// Returns `(nbrcote1, cote1, nbrcote2, cote2)`. The OCCT source writes four
/// `if / else if` arms and then four more tests that are *not* chained onto
/// them; every arm assigns all four values, so the result is the first matching
/// arm of the written order. A pair with no arm (`Enclosing`) leaves every count
/// at zero and the sweeps empty.
pub(super) fn cotes(
    c1: Categories,
    c2: Categories,
    radius: f64,
) -> (usize, [f64; 2], usize, [f64; 2]) {
    let mut nbrcote1 = 0;
    let mut nbrcote2 = 0;
    let mut cote1 = [0.0f64; 2];
    let mut cote2 = [0.0f64; 2];
    if c1.enclosed && c2.enclosed {
        nbrcote1 = 1;
        nbrcote2 = 1;
        cote1[0] = radius;
        cote2[0] = radius;
    } else if c1.enclosed && c2.outside {
        nbrcote1 = 1;
        nbrcote2 = 1;
        cote1[0] = radius;
        cote2[0] = -radius;
    } else if c1.outside && c2.enclosed {
        nbrcote1 = 1;
        nbrcote2 = 1;
        cote1[0] = -radius;
        cote2[0] = radius;
    } else if c1.outside && c2.outside {
        nbrcote1 = 1;
        nbrcote2 = 1;
        cote1[0] = -radius;
        cote2[0] = -radius;
    }
    if c1.enclosed && c2.unqualified {
        nbrcote1 = 1;
        nbrcote2 = 2;
        cote1[0] = radius;
        cote2[0] = radius;
        cote2[1] = -radius;
    }
    if c1.unqualified && c2.enclosed {
        nbrcote1 = 2;
        nbrcote2 = 1;
        cote1[0] = radius;
        cote1[1] = -radius;
        cote2[0] = radius;
    } else if c1.outside && c2.unqualified {
        nbrcote1 = 1;
        nbrcote2 = 2;
        cote1[0] = -radius;
        cote2[0] = radius;
        cote2[1] = -radius;
    }
    if c1.unqualified && c2.outside {
        nbrcote1 = 2;
        nbrcote2 = 1;
        cote1[0] = radius;
        cote1[1] = -radius;
        cote2[0] = -radius;
    } else if c1.unqualified && c2.unqualified {
        nbrcote1 = 2;
        nbrcote2 = 2;
        cote1[0] = radius;
        cote1[1] = -radius;
        cote2[0] = radius;
        cote2[1] = -radius;
    }
    (nbrcote1, cote1, nbrcote2, cote2)
}

/// `Geom2dGcc_Circ2d2TanRadGeo` (`Geom2dGcc_Circ2d2TanRadGeo.hxx:56-181`).
#[derive(Debug, Clone, Copy)]
pub struct Geom2dGccCirc2d2TanRadGeo {
    /// `WellDone` (`hxx:167`).
    well_done: bool,
    /// `NbrSol` (`hxx:168`).
    nbr_sol: usize,
    /// `cirsol` (`hxx:169`).
    cirsol: [GpCirc2d; MAX_SOL],
    /// `qualifier1` (`hxx:170`).
    qualifier1: [GccEntPosition; MAX_SOL],
    /// `qualifier2` (`hxx:171`).
    qualifier2: [GccEntPosition; MAX_SOL],
    /// `TheSame1` (`hxx:172`).
    the_same1: [i32; MAX_SOL],
    /// `TheSame2` (`hxx:173`).
    the_same2: [i32; MAX_SOL],
    /// `pnttg1sol` (`hxx:174`).
    pnttg1sol: [GpPnt2d; MAX_SOL],
    /// `pnttg2sol` (`hxx:175`).
    pnttg2sol: [GpPnt2d; MAX_SOL],
    /// `par1sol` (`hxx:176`).
    par1sol: [f64; MAX_SOL],
    /// `par2sol` (`hxx:177`).
    par2sol: [f64; MAX_SOL],
    /// `pararg1` (`hxx:178`).
    pararg1: [f64; MAX_SOL],
    /// `pararg2` (`hxx:179`).
    pararg2: [f64; MAX_SOL],
}

impl Default for Geom2dGccCirc2d2TanRadGeo {
    fn default() -> Self {
        Self::empty()
    }
}

impl Geom2dGccCirc2d2TanRadGeo {
    /// The `: cirsol(1, aNbSolMAX), ... :` initializer list every constructor
    /// shares (`cxx:60-71`, `:274-285`, `:500-510`, `:794-805`).
    pub(super) fn empty() -> Self {
        Self {
            well_done: false,
            nbr_sol: 0,
            cirsol: [GpCirc2d::default(); MAX_SOL],
            qualifier1: [GccEntPosition::NoQualifier; MAX_SOL],
            qualifier2: [GccEntPosition::NoQualifier; MAX_SOL],
            the_same1: [0; MAX_SOL],
            the_same2: [0; MAX_SOL],
            pnttg1sol: [GpPnt2d::zero(); MAX_SOL],
            pnttg2sol: [GpPnt2d::zero(); MAX_SOL],
            par1sol: [0.0; MAX_SOL],
            par2sol: [0.0; MAX_SOL],
            pararg1: [0.0; MAX_SOL],
            pararg2: [0.0; MAX_SOL],
        }
    }

    /// `IsDone()` (`cxx:1008-1011`).
    pub fn is_done(&self) -> bool {
        self.well_done
    }

    /// `NbSolutions()` (`cxx:1013-1016`).
    pub fn nb_solutions(&self) -> usize {
        self.nbr_sol
    }

    /// `ThisSolution(Index)` (`cxx:1018-1031`).
    pub fn this_solution(&self, index: usize) -> GpCirc2d {
        self.check(index);
        self.cirsol[index - 1]
    }

    /// `WhichQualifier(Index, Qualif1, Qualif2)` (`cxx:1033-1049`).
    pub fn which_qualifier(&self, index: usize) -> (GccEntPosition, GccEntPosition) {
        self.check(index);
        (self.qualifier1[index - 1], self.qualifier2[index - 1])
    }

    /// `Tangency1(Index, ParSol, ParArg, PntSol)` (`cxx:1051-1078`): a non-zero
    /// `TheSame1` throws `StdFail_NotDone` in the `else` arm.
    pub fn tangency1(&self, index: usize) -> (f64, f64, GpPnt2d) {
        self.check(index);
        assert!(
            self.the_same1[index - 1] == 0,
            "StdFail_NotDone in Geom2dGcc_Circ2d2TanRadGeo::Tangency1"
        );
        (
            self.par1sol[index - 1],
            self.pararg1[index - 1],
            self.pnttg1sol[index - 1],
        )
    }

    /// `Tangency2(Index, ParSol, ParArg, PntSol)` (`cxx:1080-1106`).
    pub fn tangency2(&self, index: usize) -> (f64, f64, GpPnt2d) {
        self.check(index);
        assert!(
            self.the_same2[index - 1] == 0,
            "StdFail_NotDone in Geom2dGcc_Circ2d2TanRadGeo::Tangency2"
        );
        (
            self.par2sol[index - 1],
            self.pararg2[index - 1],
            self.pnttg2sol[index - 1],
        )
    }

    /// `IsTheSame1(Index)` (`cxx:1108-1123`).
    pub fn is_the_same1(&self, index: usize) -> bool {
        self.check(index);
        self.the_same1[index - 1] != 0
    }

    /// `IsTheSame2(Index)` (`cxx:1125-1141`).
    pub fn is_the_same2(&self, index: usize) -> bool {
        self.check(index);
        self.the_same2[index - 1] != 0
    }

    /// `if (!WellDone) throw StdFail_NotDone()` followed by
    /// `if (Index <= 0 || Index > NbrSol) throw Standard_OutOfRange()`.
    fn check(&self, index: usize) {
        assert!(
            self.well_done,
            "StdFail_NotDone in Geom2dGcc_Circ2d2TanRadGeo"
        );
        assert!(
            index >= 1 && index <= self.nbr_sol,
            "Standard_OutOfRange in Geom2dGcc_Circ2d2TanRadGeo"
        );
    }

    /// `Geom2dGcc_Circ2d2TanRadGeo(GccEnt_QualifiedLin, Geom2dGcc_QCurve,
    /// Radius, Tolerance)` (`Geom2dGcc_Circ2d2TanRadGeo.hxx:73-77`).
    pub fn new_lin_qcurve(
        qualified1: &GccEntQualifiedLin,
        qualified2: &Geom2dGccQCurve,
        radius: f64,
        tolerance: f64,
    ) -> Self {
        lin_qcurve::new(qualified1, qualified2, radius, tolerance)
    }

    /// `Geom2dGcc_Circ2d2TanRadGeo(GccEnt_QualifiedCirc, Geom2dGcc_QCurve,
    /// Radius, Tolerance)` (`Geom2dGcc_Circ2d2TanRadGeo.hxx:65-68`).
    pub fn new_circ_qcurve(
        qualified1: &GccEntQualifiedCirc,
        qualified2: &Geom2dGccQCurve,
        radius: f64,
        tolerance: f64,
    ) -> Self {
        circ_qcurve::new(qualified1, qualified2, radius, tolerance)
    }

    /// `Geom2dGcc_Circ2d2TanRadGeo(Geom2dGcc_QCurve, const gp_Pnt2d&, Radius,
    /// Tolerance)` (`Geom2dGcc_Circ2d2TanRadGeo.hxx:91-94`).
    pub fn new_qcurve_pnt(
        qualified1: &Geom2dGccQCurve,
        point2: &GpPnt2d,
        radius: f64,
        tolerance: f64,
    ) -> Self {
        qcurve_pnt::new(qualified1, point2, radius, tolerance)
    }

    /// `Geom2dGcc_Circ2d2TanRadGeo(Geom2dGcc_QCurve, Geom2dGcc_QCurve, Radius,
    /// Tolerance)` (`Geom2dGcc_Circ2d2TanRadGeo.hxx:82-85`).
    pub fn new_two_qcurves(
        qualified1: &Geom2dGccQCurve,
        qualified2: &Geom2dGccQCurve,
        radius: f64,
        tolerance: f64,
    ) -> Self {
        two_qcurves::new(qualified1, qualified2, radius, tolerance)
    }
}

/// `!(Qualified1 ok) || !(Qualified2 ok)` -> `GccEnt_BadQualifier`
/// (`cxx:90-98`, `:304-311`, `:527-532`, `:825-831`).
pub(super) fn check_qualifier(q1_ok: bool, q2_ok: bool, who: &str) {
    assert!(q1_ok && q2_ok, "GccEnt_BadQualifier in {who}");
}

/// `if (Radius < 0.0) throw Standard_NegativeValue()` (`cxx:106-109`,
/// `:314-317`, `:535-539`, `:834-837`): the test is on the raw sign, so a NaN
/// radius passes it.
pub(super) fn check_radius(radius: f64, who: &str) {
    assert!(!(radius < 0.0), "Standard_NegativeValue in {who}");
}
