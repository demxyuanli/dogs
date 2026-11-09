//! `Geom2dGcc_Circ2d2TanOnGeo`
//! (`Geom2dGcc_Circ2d2TanOnGeo.{hxx:33-204,cxx:50-1297}`): circles tangent to
//! two operands with the centre on a third one, solved geometrically by
//! intersecting the bisector of the two operands with the carrier curve.
//!
//! | OCCT | arguments | `cxx` | here |
//! |---|---|---|---|
//! | `(QualifiedCirc, QualifiedCirc, Adaptor)` | `:50-309` | `circ_circ.rs` |
//! | `(QualifiedCirc, QualifiedLin, Adaptor)` | `:325-522` | `circ_lin.rs` |
//! | `(QualifiedLin, QualifiedLin, Adaptor)` | `:537-705` | `lin_lin.rs` |
//! | `(QualifiedCirc, gp_Pnt2d, Adaptor)` | `:716-899` | `circ_pnt.rs` |
//! | `(QualifiedLin, gp_Pnt2d, Adaptor)` | `:912-1055` | `lin_pnt.rs` |
//! | `(gp_Pnt2d, gp_Pnt2d, Adaptor)` | `:1064-1143` | `pnt_pnt.rs` |
//! | `IsDone`..`IsTheSame2` | `:1145-1297` | this file |
//!
//! Translation notes:
//! - `Adaptor2d_OffsetCurve Cu2(HCu2, 0.)` (`cxx:99`, `:374`, `:587`, `:763`,
//!   `:955`, `:1102`) is an offset adaptor built with a zero offset, so
//!   `Value(U)`, `FirstParameter()` and `LastParameter()` are the ones of the
//!   basis `Geom2dAdaptor_Curve`; the port uses the copied adaptor directly
//!   instead of porting `Adaptor2d_OffsetCurve`,
//! - the six `NCollection_Array1(1, aNbSolMAX)` families become `[_; MAX_SOL]`;
//!   `cirsol(1, 8)` is the largest bound of the class (`cxx:47`),
//! - the three local `NCollection_Array1<double> Rbid/RBid/Radius(1, 2)` of the
//!   circle/circle constructor (`cxx:74-76`) become `[f64; MAX_SOL]`: OCCT
//!   writes `Radius(3)` / `Radius(4)` when both qualifiers are unqualified
//!   (`nnsol` reaches 4 at `cxx:230-233`) and reads an uninitialised
//!   `Rbid(2)` there, so the zero-filled slot is the only defined reading,
//! - `Rbid` / `RBid` / `Radius` survive between the points of the `j` loop
//!   (`cxx:74-76` declares them outside every loop), so a point that fails the
//!   `Qualified2` guard keeps the value the previous point stored; the port
//!   keeps the same arrays across iterations,
//! - `throw GccEnt_BadQualifier()` (`cxx:82`, `:354`, `:564`, `:745`, `:941`)
//!   and `throw Standard_ConstructionError()` (`cxx:159`, `:402`, `:816`,
//!   `:981`) are panics, matching the rest of the `gcc` module,
//! - `switch (Sol->ArcType())` with a `GccInt_Pnt` solution falls into
//!   `default:` and throws `Standard_ConstructionError`; the port matches the
//!   [`GccIntBisec`](crate::gcc_int::GccIntBisec) variant directly and keeps
//!   that panic,
//! - `cxx:247-275` computes `qualifier1` / `qualifier2` from `Radius(i)` --
//!   the *bisector* index, not the solution index `k` -- while `cirsol` uses
//!   `Radius(k)`; the port keeps `radius[i - 1]` there,
//! - `cxx:199-200` writes `Rbid(1)` twice for an unqualified first circle, so
//!   `Rbid(1)` ends up as `|dist1 - R1|` and `Rbid(2)` is never written; the
//!   port keeps both assignments and the zero-filled second slot,
//! - `cxx:886-887` writes `par1sol` twice (`0.` then `ElCLib::Parameter`);
//!   only the second survives, so the port writes once.
//!
//! UNPORTED: none - all six constructors and all nine accessors.

mod circ_circ;
mod circ_lin;
mod circ_pnt;
mod lin_lin;
mod lin_pnt;
mod pnt_pnt;

use std::f64::consts::PI;
use std::sync::Arc;

use occt_core::elib::clib2d::{
    circle_value_ax22d, ellipse_value_ax22d, hyperbola_value_ax22d, parameter_circ2d,
    parameter_lin2d, parabola_value_ax22d,
};
use occt_core::gp::{GpAx22d, GpCirc2d, GpDir2d, GpLin2d, GpPnt2d};
use occt_core::intres2d::IntRes2dDomain;

use crate::adaptor_curve::Geom2dAdaptorCurve;
use crate::curve::Curve2d;
use crate::gcc::gcc_ent::{GccEntPosition, GccEntQualifiedCirc, GccEntQualifiedLin};

/// `aNbSolMAX = 8` (`cxx:48`).
pub const MAX_SOL: usize = 8;

/// `thefirst = -100000.` (`cxx:69` and the five twins).
const THEFIRST: f64 = -100000.0;

/// `thelast = 100000.` (`cxx:70` and the five twins).
const THELAST: f64 = 100000.0;

/// `Geom2dGcc_Circ2d2TanOnGeo` (`Geom2dGcc_Circ2d2TanOnGeo.hxx:43-204`).
#[derive(Clone)]
pub struct Geom2dGccCirc2d2TanOnGeo {
    /// `WellDone` (`hxx:190`).
    well_done: bool,
    /// `NbrSol` (`hxx:191`).
    nbr_sol: usize,
    /// `cirsol` (`hxx:192`).
    cirsol: [GpCirc2d; MAX_SOL],
    /// `qualifier1` (`hxx:193`).
    qualifier1: [GccEntPosition; MAX_SOL],
    /// `qualifier2` (`hxx:194`).
    qualifier2: [GccEntPosition; MAX_SOL],
    /// `TheSame1` (`hxx:195`).
    the_same1: [i32; MAX_SOL],
    /// `TheSame2` (`hxx:196`).
    the_same2: [i32; MAX_SOL],
    /// `pnttg1sol` (`hxx:197`).
    pnttg1sol: [GpPnt2d; MAX_SOL],
    /// `pnttg2sol` (`hxx:198`).
    pnttg2sol: [GpPnt2d; MAX_SOL],
    /// `pntcen` (`hxx:199`).
    pntcen: [GpPnt2d; MAX_SOL],
    /// `par1sol` (`hxx:200`).
    par1sol: [f64; MAX_SOL],
    /// `par2sol` (`hxx:201`).
    par2sol: [f64; MAX_SOL],
    /// `pararg1` (`hxx:202`).
    pararg1: [f64; MAX_SOL],
    /// `pararg2` (`hxx:203`).
    pararg2: [f64; MAX_SOL],
    /// `parcen3` (`hxx:204`).
    parcen3: [f64; MAX_SOL],
}

impl Geom2dGccCirc2d2TanOnGeo {
    /// The `: cirsol(1, aNbSolMAX), ... :` initializer list every constructor
    /// shares (`cxx:54-67` and the five twins), with `WellDone = false` and
    /// `NbrSol = 0`.
    fn empty() -> Self {
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
            pntcen: [GpPnt2d::zero(); MAX_SOL],
            par1sol: [0.0; MAX_SOL],
            par2sol: [0.0; MAX_SOL],
            pararg1: [0.0; MAX_SOL],
            pararg2: [0.0; MAX_SOL],
            parcen3: [0.0; MAX_SOL],
        }
    }

    /// `IsDone()` (`cxx:1145-1148`).
    pub fn is_done(&self) -> bool {
        self.well_done
    }

    /// `NbSolutions()` (`cxx:1150-1153`).
    pub fn nb_solutions(&self) -> usize {
        self.nbr_sol
    }

    /// `ThisSolution(Index)` (`cxx:1155-1167`).
    pub fn this_solution(&self, index: usize) -> GpCirc2d {
        self.check(index);
        self.cirsol[index - 1]
    }

    /// `WhichQualifier(Index, Qualif1, Qualif2)` (`cxx:1169-1182`).
    pub fn which_qualifier(&self, index: usize) -> (GccEntPosition, GccEntPosition) {
        self.check(index);
        (self.qualifier1[index - 1], self.qualifier2[index - 1])
    }

    /// `Tangency1(Index, ParSol, ParArg, PntSol)` (`cxx:1184-1215`).
    ///
    /// `None` where OCCT raises `StdFail_NotDone` for `TheSame1(Index) != 0`.
    pub fn tangency1(&self, index: usize) -> Option<(f64, f64, GpPnt2d)> {
        self.check(index);
        if self.the_same1[index - 1] != 0 {
            return None;
        }
        Some((
            self.par1sol[index - 1],
            self.pararg1[index - 1],
            self.pnttg1sol[index - 1],
        ))
    }

    /// `Tangency2(Index, ParSol, ParArg, PntSol)` (`cxx:1217-1243`).
    pub fn tangency2(&self, index: usize) -> Option<(f64, f64, GpPnt2d)> {
        self.check(index);
        if self.the_same2[index - 1] != 0 {
            return None;
        }
        Some((
            self.par2sol[index - 1],
            self.pararg2[index - 1],
            self.pnttg2sol[index - 1],
        ))
    }

    /// `CenterOn3(Index, ParArg, PntSol)` (`cxx:1245-1258`).
    pub fn center_on3(&self, index: usize) -> (f64, GpPnt2d) {
        self.check(index);
        (self.parcen3[index - 1], self.pntcen[index - 1])
    }

    /// `IsTheSame1(Index)` (`cxx:1260-1278`).
    pub fn is_the_same1(&self, index: usize) -> bool {
        self.check(index);
        self.the_same1[index - 1] != 0
    }

    /// `IsTheSame2(Index)` (`cxx:1280-1297`).
    pub fn is_the_same2(&self, index: usize) -> bool {
        self.check(index);
        self.the_same2[index - 1] != 0
    }

    /// `if (!WellDone) throw StdFail_NotDone()` followed by
    /// `if (Index <= 0 || Index > NbrSol) throw Standard_OutOfRange()`, the
    /// guard every accessor above shares (`cxx:1157-1163`, `:1174-1181`, ...).
    fn check(&self, index: usize) {
        assert!(self.well_done, "StdFail_NotDone in Geom2dGcc_Circ2d2TanOnGeo");
        assert!(
            index >= 1 && index <= self.nbr_sol,
            "Standard_OutOfRange in Geom2dGcc_Circ2d2TanOnGeo"
        );
    }

    /// The six public constructors (`Geom2dGcc_Circ2d2TanOnGeo.hxx:52-97`).
    pub fn new_circ_circ(
        qualified1: &GccEntQualifiedCirc,
        qualified2: &GccEntQualifiedCirc,
        on_curve: &Geom2dAdaptorCurve,
        tolerance: f64,
    ) -> Self {
        circ_circ::new(qualified1, qualified2, on_curve, tolerance)
    }

    pub fn new_circ_lin(
        qualified1: &GccEntQualifiedCirc,
        qualified2: &GccEntQualifiedLin,
        on_curve: &Geom2dAdaptorCurve,
        tolerance: f64,
    ) -> Self {
        circ_lin::new(qualified1, qualified2, on_curve, tolerance)
    }

    pub fn new_lin_lin(
        qualified1: &GccEntQualifiedLin,
        qualified2: &GccEntQualifiedLin,
        on_curve: &Geom2dAdaptorCurve,
        tolerance: f64,
    ) -> Self {
        lin_lin::new(qualified1, qualified2, on_curve, tolerance)
    }

    pub fn new_circ_pnt(
        qualified1: &GccEntQualifiedCirc,
        point2: &GpPnt2d,
        on_curve: &Geom2dAdaptorCurve,
        tolerance: f64,
    ) -> Self {
        circ_pnt::new(qualified1, point2, on_curve, tolerance)
    }

    pub fn new_lin_pnt(
        qualified1: &GccEntQualifiedLin,
        point2: &GpPnt2d,
        on_curve: &Geom2dAdaptorCurve,
        tolerance: f64,
    ) -> Self {
        lin_pnt::new(qualified1, point2, on_curve, tolerance)
    }

    pub fn new_pnt_pnt(
        point1: &GpPnt2d,
        point2: &GpPnt2d,
        on_curve: &Geom2dAdaptorCurve,
        tolerance: f64,
    ) -> Self {
        pnt_pnt::new(point1, point2, on_curve, tolerance)
    }
}

/// `gp_Ax2d(Center, dirx)` with `dirx = gp_Dir2d::D::X` (`cxx:90`, `:357`,
/// `:569`, `:751`, `:943`, `:1094`).
fn axis_at(center: GpPnt2d) -> GpAx22d {
    GpAx22d::from_xdir(center, GpDir2d::new(1.0, 0.0).expect("gp::DX2d()"))
}

/// `gp_Circ2d(gp_Ax2d(Center, dirx), Radius)`.
fn circle_at(center: GpPnt2d, radius: f64) -> GpCirc2d {
    GpCirc2d::new(axis_at(center), radius)
}

/// `gp_Dir2d(from.XY() - to.XY())`: panics where OCCT raises
/// `Standard_ConstructionError` on coincident points.
fn dir_from(from: &GpPnt2d, to: &GpPnt2d) -> GpDir2d {
    GpDir2d::new(from.x() - to.x(), from.y() - to.y()).expect("gp_Dir2d from two distinct points")
}

/// `gp_Dir2d(sign * gp_XY(-dir.Y(), dir.X()))` (`cxx:510`, `:693`, `:1040`):
/// a zero `sign` selects the zero vector and panics as OCCT does.
fn signed_normal(sign: f64, dir: &GpDir2d) -> GpDir2d {
    GpDir2d::new(sign * -dir.y(), sign * dir.x()).expect("gp_Dir2d of the signed normal")
}

/// `ElCLib::Parameter(cirsol, P)` (`ElCLib.cxx:1096-1100`).
fn par_on_circle(circle: &GpCirc2d, p: &GpPnt2d) -> f64 {
    parameter_circ2d(circle, p)
}

/// `ElCLib::Parameter(line, P)` (`ElCLib.cxx:1081-1085`).
fn par_on_line(line: &GpLin2d, p: &GpPnt2d) -> f64 {
    parameter_lin2d(line, p)
}

/// `IntRes2d_Domain D2(C2.Value(firstparam), firstparam, Tol, C2.Value(lastparam),
/// lastparam, Tol)` (`cxx:102-106` and the five twins) together with the curve
/// the six `Perform` calls hand to the intersector.
///
/// `C2` is `Adaptor2d_OffsetCurve(Cu2, 0.)`, so its `Value` / `FirstParameter`
/// / `LastParameter` are those of the basis adaptor.
fn carrier_domain(
    on_curve: &Geom2dAdaptorCurve,
    tol: f64,
) -> (IntRes2dDomain, Arc<dyn Curve2d>) {
    let first = on_curve.first_parameter().max(THEFIRST);
    let last = on_curve.last_parameter().min(THELAST);
    let p_first = on_curve
        .value(first)
        .expect("Geom2dAdaptor_Curve::Value on the carrier");
    let p_last = on_curve
        .value(last)
        .expect("Geom2dAdaptor_Curve::Value on the carrier");
    let domain = IntRes2dDomain::bounded(&p_first, first, tol, &p_last, last, tol);
    let curve = on_curve
        .curve()
        .expect("Geom2dAdaptor_Curve::Curve on the carrier")
        .clone();
    (domain, curve)
}

/// `IntRes2d_Domain D1(ElCLib::Value(0., Circ), 0., Tol1,
/// ElCLib::Value(2. * M_PI, Circ), 2. * M_PI, Tol2); D1.SetEquivalentParameters(0.,
/// 2. * M_PI);` (`cxx:119-123`, `:779-783`, `:792-796`).
fn circle_domain(circle: &GpCirc2d, tol1: f64, tol2: f64) -> IntRes2dDomain {
    let first = circle_value_ax22d(0.0, circle.position(), circle.radius());
    let last = circle_value_ax22d(2.0 * PI, circle.position(), circle.radius());
    let mut domain = IntRes2dDomain::bounded(&first, 0.0, tol1, &last, 2.0 * PI, tol2);
    domain.set_equivalent_parameters(0.0, 2.0 * PI);
    domain
}

/// The same construction for an `gp_Elips2d` (`cxx:130-134`).
fn ellipse_domain(
    ellipse: &occt_core::gp::GpElips2d,
    tol1: f64,
    tol2: f64,
) -> IntRes2dDomain {
    let first = ellipse_value_ax22d(
        0.0,
        ellipse.axis(),
        ellipse.major_radius,
        ellipse.minor_radius,
    );
    let last = ellipse_value_ax22d(
        2.0 * PI,
        ellipse.axis(),
        ellipse.major_radius,
        ellipse.minor_radius,
    );
    let mut domain = IntRes2dDomain::bounded(&first, 0.0, tol1, &last, 2.0 * PI, tol2);
    domain.set_equivalent_parameters(0.0, 2.0 * PI);
    domain
}

/// `IntRes2d_Domain D1(ElCLib::Value(-4., Hypr), -4., Tol1, ElCLib::Value(4.,
/// Hypr), 4., Tol2);` (`cxx:141-146`).
fn hyperbola_domain(
    hyperbola: &occt_core::gp::GpHypr2d,
    tol1: f64,
    tol2: f64,
) -> IntRes2dDomain {
    let first = hyperbola_value_ax22d(
        -4.0,
        hyperbola.axis(),
        hyperbola.major_radius,
        hyperbola.minor_radius,
    );
    let last = hyperbola_value_ax22d(
        4.0,
        hyperbola.axis(),
        hyperbola.major_radius,
        hyperbola.minor_radius,
    );
    IntRes2dDomain::bounded(&first, -4.0, tol1, &last, 4.0, tol2)
}

/// `IntRes2d_Domain D1(ElCLib::Value(-40, Parab), -40, Tol1,
/// ElCLib::Value(40, Parab), 40, Tol1);` (`cxx:392-397`, `:971-976`) -- note
/// the last tolerance is `Tol1`, not `Tol2`.
fn parabola_domain(parabola: &occt_core::gp::GpParab2d, tol1: f64) -> IntRes2dDomain {
    let first = parabola_value_ax22d(-40.0, parabola.axis(), parabola.focal);
    let last = parabola_value_ax22d(40.0, parabola.axis(), parabola.focal);
    IntRes2dDomain::bounded(&first, -40.0, tol1, &last, 40.0, tol1)
}

/// `IntRes2d_Domain D1;` -- the infinite domain of the line arms
/// (`cxx:153`, `:387`, `:590`, `:955`, `:1110`).
fn empty_domain() -> IntRes2dDomain {
    IntRes2dDomain::new()
}
