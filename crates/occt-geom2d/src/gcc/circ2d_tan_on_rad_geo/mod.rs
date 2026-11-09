//! `Geom2dGcc_Circ2dTanOnRadGeo`
//! (`Geom2dGcc_Circ2dTanOnRadGeo.{hxx:44-203,cxx:61-912}`): circles of a
//! prescribed radius tangent to one operand with the centre on a second one.
//!
//! Unlike `Geom2dGcc_Circ2d2TanOnGeo` the radius is an input, so no bisector is
//! built: the centre locus is the operand offset by `-+Radius` (a
//! `Adaptor2d_OffsetCurve`) intersected with the carrier. Six constructors
//! cover the carrier being a line, a circle or a general curve, and the
//! operand being a circle, a line, a point or a general curve.
//!
//! | OCCT constructor | arguments | `cxx` | file |
//! |---|---|---|---|
//! | 1 | `Geom2dGcc_QCurve, gp_Lin2d` | `:61-167` | [`on_line`] |
//! | 2 | `Geom2dGcc_QCurve, gp_Circ2d` | `:185-298` | [`on_circ`] |
//! | 3 | `GccEnt_QualifiedCirc, Adaptor` | `:317-445` | [`circ_on_curve`] |
//! | 4 | `GccEnt_QualifiedLin, Adaptor` | `:466-583` | [`lin_on_curve`] |
//! | 5 | `Geom2dGcc_QCurve, Adaptor` | `:604-717` | [`qcurve_on_curve`] |
//! | 6 | `gp_Pnt2d, Adaptor` | `:736-816` | [`pnt_on_curve`] |
//!
//! Translation notes:
//! - `Adaptor2d_OffsetCurve(HCu1, -Coef(j))` and `Adaptor2d_OffsetCurve(HCu1,
//!   -cote1(j))` (`cxx:136`, `:267`, `:675`) evaluate through
//!   `Geom2d_OffsetCurveUtils::CalculateD0` (`Adaptor2d_OffsetCurve.cxx:295-309`),
//!   which is what [`crate::offset::Geom2dOffsetCurve`] implements, so the port
//!   wraps the operand in that curve and hands it to the intersectors as a
//!   `dyn Curve2d`; `FirstParameter` / `LastParameter` delegate to the basis in
//!   both,
//! - `Adaptor2d_OffsetCurve C2(HOnCurv)` (`cxx:683`) has the default offset
//!   `0.0`, i.e. `Value` is the basis value and both parameters are the basis
//!   parameters, so the port intersects the carrier as is,
//! - `par1sol` is written **before** `pnttg1sol` in all six constructors
//!   (`cxx:159-160`, `:290-291`, `:436-437`, `:575-576`, `:709-710`,
//!   `:808-809`), so it is `ElCLib::Parameter` of the **previous** iteration's
//!   tangency point, and of the zero-initialised slot on the first one; the
//!   port keeps that order,
//! - `CenterOn3` reports `pnttg1sol(Index)` instead of `pntcen3(Index)`
//!   (`cxx:890`), so the front end's `pntcen3` receives the tangency point; the
//!   port keeps that,
//! - `ThisSolution` has no `WellDone` guard (`cxx:830-840`) while the other
//!   accessors do; the port keeps that,
//! - `Standard_NegativeValue` tests the raw `Radius < 0.0` (`cxx:103`, `:229`,
//!   `:362`, `:509`, `:646`, `:770`), and `GccEnt_BadQualifier` is thrown
//!   before it,
//! - the `IsEnclosing` qualifier passes the guard but no `cote1` branch
//!   (`cxx:110-127`, `:236-253`, `:368-385`, `:514-533`, `:651-668`), so
//!   `nbrcote1` stays `0`, the loop never runs and `WellDone` stays false,
//! - `WellDone = true` sits inside `if (Intp.IsDone())` but outside
//!   `if (!Intp.IsEmpty())` (`cxx:161-166`), so an empty intersection still
//!   marks the solver done.
//!
//! UNPORTED: none - all six constructors, the eight accessors and the field
//! set.

use occt_core::gp::{GpCirc2d, GpLin2d, GpPnt2d};

use crate::adaptor_curve::Geom2dAdaptorCurve;
use crate::gcc::gcc_ent::{GccEntPosition, GccEntQualifiedCirc, GccEntQualifiedLin};
use crate::gcc::qualified_curve::Geom2dGccQCurve;

mod circ_on_curve;
mod lin_on_curve;
mod on_circ;
mod on_line;
mod pnt_on_curve;
mod qcurve_on_curve;

/// `aNbSolMAX = 8` (`cxx:44`).
pub const MAX_SOL: usize = 8;

/// `thefirst = -100000.` / `thelast = 100000.` (`cxx:88-89`, `:213-214`,
/// `:345-346`, `:490-491`, `:629-630`, `:763-764`).

/// `Geom2dGcc_Circ2dTanOnRadGeo` (`Geom2dGcc_Circ2dTanOnRadGeo.hxx:44-203`).
#[derive(Clone)]
pub struct Geom2dGccCirc2dTanOnRadGeo {
    /// `WellDone` (`hxx:186`).
    well_done: bool,
    /// `NbrSol` (`hxx:187`).
    nbr_sol: usize,
    /// `cirsol` (`hxx:188`).
    cirsol: [GpCirc2d; MAX_SOL],
    /// `qualifier1` (`hxx:189`).
    qualifier1: [GccEntPosition; MAX_SOL],
    /// `TheSame1` (`hxx:190`).
    the_same1: [i32; MAX_SOL],
    /// `pnttg1sol` (`hxx:191`).
    pnttg1sol: [GpPnt2d; MAX_SOL],
    /// `pntcen3` (`hxx:192`).
    pntcen3: [GpPnt2d; MAX_SOL],
    /// `par1sol` (`hxx:193`).
    par1sol: [f64; MAX_SOL],
    /// `pararg1` (`hxx:194`).
    pararg1: [f64; MAX_SOL],
    /// `parcen3` (`hxx:195`).
    parcen3: [f64; MAX_SOL],
}

impl Default for Geom2dGccCirc2dTanOnRadGeo {
    fn default() -> Self {
        Self::empty()
    }
}

impl Geom2dGccCirc2dTanOnRadGeo {
    /// The `: cirsol(1, aNbSolMAX), ... :` initializer list every constructor
    /// shares (`cxx:73-79`, `:197-205`, `:328-336`, `:477-482`, `:614-621`,
    /// `:747-754`).
    pub(super) fn empty() -> Self {
        Self {
            well_done: false,
            nbr_sol: 0,
            cirsol: [GpCirc2d::default(); MAX_SOL],
            qualifier1: [GccEntPosition::NoQualifier; MAX_SOL],
            the_same1: [0; MAX_SOL],
            pnttg1sol: [GpPnt2d::zero(); MAX_SOL],
            pntcen3: [GpPnt2d::zero(); MAX_SOL],
            par1sol: [0.0; MAX_SOL],
            pararg1: [0.0; MAX_SOL],
            parcen3: [0.0; MAX_SOL],
        }
    }

    /// `IsDone()` (`cxx:820-823`).
    pub fn is_done(&self) -> bool {
        self.well_done
    }

    /// `NbSolutions()` (`cxx:825-828`).
    pub fn nb_solutions(&self) -> usize {
        self.nbr_sol
    }

    /// `ThisSolution(Index)` (`cxx:830-840`): OCCT guards the index only.
    pub fn this_solution(&self, index: usize) -> GpCirc2d {
        assert!(
            index >= 1 && index <= self.nbr_sol,
            "Standard_OutOfRange in Geom2dGcc_Circ2dTanOnRadGeo"
        );
        self.cirsol[index - 1]
    }

    /// `WhichQualifier(Index, Qualif1)` (`cxx:842-855`).
    pub fn which_qualifier(&self, index: usize) -> GccEntPosition {
        self.check(index);
        self.qualifier1[index - 1]
    }

    /// `Tangency1(Index, ParSol, ParArg, PntSol)` (`cxx:857-875`): no
    /// `TheSame1` guard here, unlike the front end's.
    pub fn tangency1(&self, index: usize) -> (f64, f64, GpPnt2d) {
        self.check(index);
        (
            self.par1sol[index - 1],
            self.pararg1[index - 1],
            self.pnttg1sol[index - 1],
        )
    }

    /// `CenterOn3(Index, ParArg, PntSol)` (`cxx:877-892`).
    ///
    /// The point is `pnttg1sol(Index)` (`cxx:890`), not the stored centre.
    pub fn center_on3(&self, index: usize) -> (f64, GpPnt2d) {
        self.check(index);
        (self.parcen3[index - 1], self.pnttg1sol[index - 1])
    }

    /// `IsTheSame1(Index)` (`cxx:894-912`).
    pub fn is_the_same1(&self, index: usize) -> bool {
        self.check(index);
        self.the_same1[index - 1] != 0
    }

    /// `if (!WellDone) throw StdFail_NotDone()` followed by
    /// `if (Index <= 0 || Index > NbrSol) throw Standard_OutOfRange()`.
    fn check(&self, index: usize) {
        assert!(
            self.well_done,
            "StdFail_NotDone in Geom2dGcc_Circ2dTanOnRadGeo"
        );
        assert!(
            index >= 1 && index <= self.nbr_sol,
            "Standard_OutOfRange in Geom2dGcc_Circ2dTanOnRadGeo"
        );
    }

    /// The six public constructors (`Geom2dGcc_Circ2dTanOnRadGeo.hxx:57-139`).
    pub fn new_qcurve_on_line(
        qualified1: &Geom2dGccQCurve,
        on_line: &GpLin2d,
        radius: f64,
        tolerance: f64,
    ) -> Self {
        on_line::new(qualified1, on_line, radius, tolerance)
    }

    pub fn new_qcurve_on_circ(
        qualified1: &Geom2dGccQCurve,
        on_circ: &GpCirc2d,
        radius: f64,
        tolerance: f64,
    ) -> Self {
        on_circ::new(qualified1, on_circ, radius, tolerance)
    }

    pub fn new_circ_on_curve(
        qualified1: &GccEntQualifiedCirc,
        on_curve: &Geom2dAdaptorCurve,
        radius: f64,
        tolerance: f64,
    ) -> Self {
        circ_on_curve::new(qualified1, on_curve, radius, tolerance)
    }

    pub fn new_lin_on_curve(
        qualified1: &GccEntQualifiedLin,
        on_curve: &Geom2dAdaptorCurve,
        radius: f64,
        tolerance: f64,
    ) -> Self {
        lin_on_curve::new(qualified1, on_curve, radius, tolerance)
    }

    pub fn new_qcurve_on_curve(
        qualified1: &Geom2dGccQCurve,
        on_curve: &Geom2dAdaptorCurve,
        radius: f64,
        tolerance: f64,
    ) -> Self {
        qcurve_on_curve::new(qualified1, on_curve, radius, tolerance)
    }

    pub fn new_pnt_on_curve(
        point1: &GpPnt2d,
        on_curve: &Geom2dAdaptorCurve,
        radius: f64,
        tolerance: f64,
    ) -> Self {
        pnt_on_curve::new(point1, on_curve, radius, tolerance)
    }
}

/// `gp_Ax2d(Center, gp_Dir2d::D::X)` / `gp_Circ2d(gp_Ax2d(Center, dirx),
/// Radius)` (`cxx:155`, `:286`, `:417`, `:558`, `:706`, `:803`).
pub(super) use crate::gcc::tan_geo_common::{
    circle_at, circle_domain, curve_domain, offset_of, value_on_circ, value_on_line, THEFIRST,
    THELAST,
};

/// `if (Qualified1.IsEnclosed()) { nbrcote1 = 1; cote1(1) = Radius; } else if
/// (Qualified1.IsOutside()) { nbrcote1 = 1; cote1(1) = -Radius; } else if
/// (Qualified1.IsUnqualified()) { nbrcote1 = 2; cote1(1) = Radius; cote1(2) =
/// -Radius; }` (`cxx:110-127`, `:236-253`, `:368-385`, `:514-533`,
/// `:651-668`).
///
/// `IsEnclosing` matches none of the branches, so the count stays zero and the
/// caller's loop never runs.
pub(super) fn cotes(
    is_enclosed: bool,
    is_outside: bool,
    is_unqualified: bool,
    radius: f64,
) -> (usize, [f64; 2]) {
    if is_enclosed {
        (1, [radius, 0.0])
    } else if is_outside {
        (1, [-radius, 0.0])
    } else if is_unqualified {
        (2, [radius, -radius])
    } else {
        (0, [0.0, 0.0])
    }
}

