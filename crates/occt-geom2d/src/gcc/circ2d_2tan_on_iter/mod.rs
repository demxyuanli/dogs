//! `Geom2dGcc_Circ2d2TanOnIter`
//! (`Geom2dGcc_Circ2d2TanOnIter.{hxx:34-270,cxx:45-1438}`): the iterative solver
//! behind `Geom2dGcc_Circ2d2TanOn` when at least one operand is a general curve.
//!
//! A single solution is computed: a circle tangent to two operands whose center
//! lies on a third. The two tangent operands can be a line, a circle, a point or
//! a `Geom2dGcc_QCurve`; the carrier can be a line, a circle or a
//! `Geom2dAdaptor_Curve`. The twelve OCCT constructors are named
//! `new_<arg1>_<arg2>_on_<carrier>` with `cu` (curve), `ci` (circle), `li`
//! (line) and `pt` (point):
//!
//! | OCCT | `cxx` | here |
//! |---|---|---|
//! | `(QualifiedLin, QCurve, gp_Lin2d)` | `:45-151` | [`new_li_cu_on_li`](Geom2dGccCirc2d2TanOnIter::new_li_cu_on_li) |
//! | `(QCurve, QCurve, gp_Lin2d)` | `:153-272` | [`new_cu_cu_on_li`](Geom2dGccCirc2d2TanOnIter::new_cu_cu_on_li) |
//! | `(QCurve, gp_Pnt2d, gp_Lin2d)` | `:274-364` | [`new_cu_pt_on_li`](Geom2dGccCirc2d2TanOnIter::new_cu_pt_on_li) |
//! | `(QualifiedCirc, QCurve, gp_Lin2d)` | `:366-480` | [`new_ci_cu_on_li`](Geom2dGccCirc2d2TanOnIter::new_ci_cu_on_li) |
//! | `(QualifiedCirc, QCurve, gp_Circ2d)` | `:482-595` | [`new_ci_cu_on_ci`](Geom2dGccCirc2d2TanOnIter::new_ci_cu_on_ci) |
//! | `(QualifiedLin, QCurve, gp_Circ2d)` | `:597-704` | [`new_li_cu_on_ci`](Geom2dGccCirc2d2TanOnIter::new_li_cu_on_ci) |
//! | `(QCurve, QCurve, gp_Circ2d)` | `:706-827` | [`new_cu_cu_on_ci`](Geom2dGccCirc2d2TanOnIter::new_cu_cu_on_ci) |
//! | `(QCurve, gp_Pnt2d, gp_Circ2d)` | `:829-918` | [`new_cu_pt_on_ci`](Geom2dGccCirc2d2TanOnIter::new_cu_pt_on_ci) |
//! | `(QCurve, QCurve, Geom2dAdaptor_Curve)` | `:920-1035` | [`new_cu_cu_on_cu`](Geom2dGccCirc2d2TanOnIter::new_cu_cu_on_cu) |
//! | `(QualifiedCirc, QCurve, Geom2dAdaptor_Curve)` | `:1037-1146` | [`new_ci_cu_on_cu`](Geom2dGccCirc2d2TanOnIter::new_ci_cu_on_cu) |
//! | `(QualifiedLin, QCurve, Geom2dAdaptor_Curve)` | `:1148-1253` | [`new_li_cu_on_cu`](Geom2dGccCirc2d2TanOnIter::new_li_cu_on_cu) |
//! | `(QCurve, gp_Pnt2d, Geom2dAdaptor_Curve)` | `:1255-1344` | [`new_cu_pt_on_cu`](Geom2dGccCirc2d2TanOnIter::new_cu_pt_on_cu) |
//!
//! Translation notes:
//! - `GccEnt_BadQualifier` (`cxx:69`, `:177`, `:292`, `:390`, `:505`, `:619`,
//!   `:730`, `:850`, `:944`, `:1059`, `:1169`, `:1274`) and `StdFail_NotDone`
//!   (`cxx:1361`, `:1375`, `:1393`, `:1407`, `:1421`, `:1431`) are panics here,
//!   matching the rest of the `gcc` module;
//! - `Func.Value(Ufirst, Umin)` right after `Root.Perform` (`cxx:100` and the
//!   eleven twins) only refills the output vector OCCT hands in as `Umin`; the
//!   port's system is stateless, so the call is dropped;
//! - `new_ci_cu_on_li` (`cxx:436`) evaluates `OnLine` at `Ufirst(1)` -- the
//!   circle parameter, not the line parameter `Ufirst(3)`. The port keeps
//!   `x.value(1)`;
//! - `new_ci_cu_on_ci` (`cxx:538`) calls the
//!   `(gp_Circ2d, Geom2dAdaptor_Curve, gp_Circ2d)` constructor, whose `TheType`
//!   is `Geom2dGcc_CuCuOnCi` (see `function_tan_cu_cu_on_cu.rs`), so the central
//!   curve is never used; the port keeps that dispatch;
//! - `new_cu_cu_on_li` returns early when `(dist1 + dist2) / 2 < Tol`
//!   (`cxx:220-223`), before the `IsDone` bookkeeping;
//! - `new_li_cu_on_li` sets all four tolerances from `1.e-15` (`cxx:90-93`),
//!   where every other constructor uses `Tol / 10.` for the fourth and `Tol`
//!   is `|Tolerance|`;
//! - `IsUnqualified()` seeds `qualifier1` / `qualifier2` with
//!   `GccEnt_noqualifier`, and `WhichQualifier` only throws when `!IsDone()`.
//!
//! UNPORTED: none. The `.cxx` has no other member.

use std::f64::consts::PI;

use occt_core::elib::clib2d::{
    circle_d1_ax22d, circle_d2_ax22d, circle_parameter_ax22d, circle_value_ax22d, line_d1_ax2d,
    line_value_ax2d,
};
use occt_core::gp::{GpAx22d, GpCirc2d, GpDir2d, GpLin2d, GpPnt2d, GpVec2d};
use occt_core::precision::REAL_SMALL;
use occt_math::{MathFunctionSetRoot, MathVector};

use crate::adaptor_curve::Geom2dAdaptorCurve;
use crate::gcc::curve_tool;
use crate::gcc::function_tan_cu_cu_on_cu::Geom2dGccFunctionTanCuCuOnCu;
use crate::gcc::gcc_ent::{GccEntPosition, GccEntQualifiedCirc, GccEntQualifiedLin};
use crate::gcc::qualified_curve::Geom2dGccQCurve;

/// `RealFirst()` (`Standard_Real.hxx:88-91`) = `-DBL_MAX`.
const REAL_FIRST: f64 = -f64::MAX;

/// `RealLast()` (`Standard_Real.hxx:94-97`) = `DBL_MAX`.
const REAL_LAST: f64 = f64::MAX;

/// `2.e-15 * M_PI`, the angular tolerance of `cxx:412`, `:530`, `:644`, `:754`,
/// `:868`, `:1083`, `:1289`.
const ANGULAR_TOL: f64 = 2.0e-15 * PI;

/// `1.e-15`, the linear tolerance of `cxx:90`, `:416`, `:642`, `:1193`, `:1293`.
const LINEAR_TOL: f64 = 1.0e-15;

/// `gp_Ax2d(P, gp_Dir2d(gp_Dir2d::D::X))`, the axis of every `cirsol`.
fn axis_at(center: GpPnt2d) -> GpAx22d {
    GpAx22d::from_xdir(center, GpDir2d::new(1.0, 0.0).expect("gp::DX2d()"))
}

/// `gp_Vec2d(P1, P2)` (`gp_Vec2d.cxx:83-90`).
fn vec_between(from: &GpPnt2d, to: &GpPnt2d) -> GpVec2d {
    GpVec2d::new(to.x() - from.x(), to.y() - from.y())
}

/// `Vec.Angle(Tan)` guarded by `norme >= gp::Resolution()`
/// (`cxx:447-455` and its eleven twins).
fn angle_or_zero(vec: &GpVec2d, tan: &GpVec2d) -> f64 {
    if vec.magnitude() >= REAL_SMALL && tan.magnitude() >= REAL_SMALL {
        vec.angle(tan)
    } else {
        0.0
    }
}

/// `point3.XY().Dot(gp_XY(-Dir.Y(), Dir.X()))` of `cxx:683`, `:1229`.
fn signed_side(point: &GpPnt2d, line: &GpLin2d) -> f64 {
    let d = line.direction();
    point.x() * (-d.y()) + point.y() * d.x()
}

/// `ElCLib::Parameter(cirsol, P)` (`ElCLib.cxx:1096-1100`).
fn par_on_solution(sol: &GpCirc2d, p: &GpPnt2d) -> f64 {
    circle_parameter_ax22d(sol.position(), p)
}

/// The `Qualified1.IsUnqualified() || (IsEnclosing() && Rsol >= R1 && dist <=
/// Rsol) || (IsOutside() && dist >= Rsol) || (IsEnclosed() && Rsol <= R1 &&
/// dist <= Rsol)` test of `cxx:457-460` and its twins.
fn accept_circ(q: GccEntPosition, r1: f64, rsol: f64, dist: f64) -> bool {
    q.is_unqualified()
        || (q.is_enclosing() && rsol >= r1 && dist <= rsol)
        || (q.is_outside() && dist >= rsol)
        || (q.is_enclosed() && rsol <= r1 && dist <= rsol)
}

/// The `Qualified1.IsUnqualified() || (IsOutside() && pscal <= 0.) ||
/// (IsEnclosed() && pscal >= 0.)` test of `cxx:683-686` and its twins.
fn accept_lin(q: GccEntPosition, pscal: f64) -> bool {
    q.is_unqualified() || (q.is_outside() && pscal <= 0.) || (q.is_enclosed() && pscal >= 0.)
}

/// The `IsUnqualified() || (IsEnclosing() && angle <= 0.) || (IsOutside() &&
/// angle >= 0.) || (IsEnclosed() && angle <= 0.)` test of `cxx:457-465` and its
/// twins, used for a curve operand.
fn accept_curve(q: GccEntPosition, angle: f64) -> bool {
    q.is_unqualified()
        || (q.is_enclosing() && angle <= 0.)
        || (q.is_outside() && angle >= 0.)
        || (q.is_enclosed() && angle <= 0.)
}

/// `Geom2dGcc_Circ2d2TanOnIter` (`Geom2dGcc_Circ2d2TanOnIter.hxx:34-270`).
#[derive(Clone)]
pub struct Geom2dGccCirc2d2TanOnIter {
    /// `WellDone` (`hxx:254`).
    well_done: bool,
    /// `cirsol` (`hxx:255`).
    cirsol: GpCirc2d,
    /// `qualifier1` (`hxx:256`).
    qualifier1: GccEntPosition,
    /// `qualifier2` (`hxx:257`).
    qualifier2: GccEntPosition,
    /// `TheSame1` (`hxx:258`).
    the_same1: bool,
    /// `TheSame2` (`hxx:259`).
    the_same2: bool,
    /// `pnttg1sol` (`hxx:260`).
    pnttg1sol: GpPnt2d,
    /// `pnttg2sol` (`hxx:261`).
    pnttg2sol: GpPnt2d,
    /// `pntcen` (`hxx:262`).
    pntcen: GpPnt2d,
    /// `par1sol` (`hxx:263`).
    par1sol: f64,
    /// `par2sol` (`hxx:264`).
    par2sol: f64,
    /// `pararg1` (`hxx:265`).
    pararg1: f64,
    /// `pararg2` (`hxx:266`).
    pararg2: f64,
    /// `parcen3` (`hxx:267`).
    parcen3: f64,
}

impl Geom2dGccCirc2d2TanOnIter {
    /// The member state every constructor starts from (`cxx:52-65` and twins).
    fn not_done() -> Self {
        Self {
            well_done: false,
            cirsol: GpCirc2d::default(),
            qualifier1: GccEntPosition::NoQualifier,
            qualifier2: GccEntPosition::NoQualifier,
            the_same1: false,
            the_same2: false,
            pnttg1sol: GpPnt2d::zero(),
            pnttg2sol: GpPnt2d::zero(),
            pntcen: GpPnt2d::zero(),
            par1sol: 0.0,
            par2sol: 0.0,
            pararg1: 0.0,
            pararg2: 0.0,
            parcen3: 0.0,
        }
    }

    /// `Root.Perform(Func, Ufirst, Umin, Umax)` followed by `Root.Root(Ufirst)`
    /// (`cxx:101-105` and twins), answering `None` when `!Root.IsDone()`.
    fn solve(
        func: &mut Geom2dGccFunctionTanCuCuOnCu,
        tol: &MathVector,
        ufirst: &MathVector,
        umin: &MathVector,
        umax: &MathVector,
    ) -> Option<MathVector> {
        let mut root = MathFunctionSetRoot::new(func, tol, 100);
        root.perform_with_bounds(func, ufirst, umin, umax, false);
        if root.is_done() {
            Some(root.root())
        } else {
            None
        }
    }

    /// The twelve `qualifier1 = ...; qualifier2 = ...; ...; WellDone = true;`
    /// tails (`cxx:466-478` and twins).
    #[allow(clippy::too_many_arguments)]
    fn commit(
        &mut self,
        sol: GpCirc2d,
        q1: GccEntPosition,
        q2: GccEntPosition,
        pnt1: GpPnt2d,
        pararg1: f64,
        par1sol: f64,
        pnt2: GpPnt2d,
        pararg2: f64,
        par2sol: f64,
        center: GpPnt2d,
        parcen3: f64,
    ) {
        self.cirsol = sol;
        self.qualifier1 = q1;
        self.qualifier2 = q2;
        self.pnttg1sol = pnt1;
        self.pararg1 = pararg1;
        self.par1sol = par1sol;
        self.pnttg2sol = pnt2;
        self.pararg2 = pararg2;
        self.par2sol = par2sol;
        self.pntcen = center;
        self.parcen3 = parcen3;
        self.well_done = true;
    }

    /// `IsDone()` (`cxx:1346-1349`).
    pub fn is_done(&self) -> bool {
        self.well_done
    }

    /// `ThisSolution()` (`cxx:1351-1354`), panicking where OCCT raises
    /// `StdFail_NotDone`.
    pub fn this_solution(&self) -> GpCirc2d {
        if !self.well_done {
            panic!("StdFail_NotDone in Geom2dGcc_Circ2d2TanOnIter");
        }
        self.cirsol
    }

    /// `WhichQualifier(Qualif1, Qualif2)` (`cxx:1356-1368`).
    pub fn which_qualifier(&self) -> (GccEntPosition, GccEntPosition) {
        if !self.well_done {
            panic!("StdFail_NotDone in Geom2dGcc_Circ2d2TanOnIter");
        }
        (self.qualifier1, self.qualifier2)
    }

    /// `Tangency1(ParSol, ParArg, PntSol)` (`cxx:1370-1388`).
    pub fn tangency1(&self) -> (f64, f64, GpPnt2d) {
        if !self.well_done {
            panic!("StdFail_NotDone in Geom2dGcc_Circ2d2TanOnIter");
        }
        if !self.the_same1 {
            (0.0, 0.0, self.pnttg1sol)
        } else {
            panic!("StdFail_NotDone in Geom2dGcc_Circ2d2TanOnIter");
        }
    }

    /// `Tangency2(ParSol, ParArg, PntSol)` (`cxx:1390-1404`).
    pub fn tangency2(&self) -> (f64, f64, GpPnt2d) {
        if !self.well_done {
            panic!("StdFail_NotDone in Geom2dGcc_Circ2d2TanOnIter");
        }
        (0.0, 0.0, self.pnttg2sol)
    }

    /// `CenterOn3(ParArg, PntSol)` (`cxx:1406-1417`).
    pub fn center_on3(&self) -> (f64, GpPnt2d) {
        if !self.well_done {
            panic!("StdFail_NotDone in Geom2dGcc_Circ2d2TanOnIter");
        }
        (0.0, self.pntcen)
    }

    /// `IsTheSame1()` (`cxx:1419-1426`).
    pub fn is_the_same1(&self) -> bool {
        if !self.well_done {
            panic!("StdFail_NotDone in Geom2dGcc_Circ2d2TanOnIter");
        }
        self.the_same1
    }

    /// `IsTheSame2()` (`cxx:1428-1438`).
    pub fn is_the_same2(&self) -> bool {
        if !self.well_done {
            panic!("StdFail_NotDone in Geom2dGcc_Circ2d2TanOnIter");
        }
        self.the_same2
    }
}

mod on_circle;
mod on_curve;
mod on_line;
