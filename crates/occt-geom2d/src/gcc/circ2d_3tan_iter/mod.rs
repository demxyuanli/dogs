//! `Geom2dGcc_Circ2d3TanIter`
//! (`Geom2dGcc_Circ2d3TanIter.{hxx:32-216,cxx:44-1568}`): the iterative solver
//! behind `Geom2dGcc_Circ2d3Tan` when the three operands are not all of the
//! analytic kind.
//!
//! One solution is computed. The three operands are qualified circles, lines,
//! `Geom2dGcc_QCurve`s or plain points; each curve operand contributes one
//! unknown, each circle one angle in `[0, 2*pi)`, each line one unbounded
//! parameter, and each point a radius-0 circle whose angle starts at `pi`. Once
//! `math_FunctionSetRoot` converges the three tangency points are handed to the
//! analytic `GccAna_Circ2d3Tan(point, point, point)` and the resulting circle is
//! accepted only if every operand's qualifier agrees and every tangency dot is
//! `<= 1.e-12`.
//!
//! | OCCT constructor | `cxx` | here |
//! |---|---|---|
//! | `(QualifiedCirc, QualifiedCirc, QCurve)` | `:335-484` | [`new_ci_ci_cu`](Geom2dGccCirc2d3TanIter::new_ci_ci_cu) |
//! | `(QualifiedCirc, QCurve, QCurve)` | `:187-332` | [`new_ci_cu_cu`](Geom2dGccCirc2d3TanIter::new_ci_cu_cu) |
//! | `(QualifiedCirc, QualifiedLin, QCurve)` | `:1171-1315` | [`new_ci_li_cu`](Geom2dGccCirc2d3TanIter::new_ci_li_cu) |
//! | `(QualifiedCirc, QCurve, gp_Pnt2d)` | `:1318-1453` | [`new_ci_cu_pt`](Geom2dGccCirc2d3TanIter::new_ci_cu_pt) |
//! | `(QualifiedLin, QualifiedLin, QCurve)` | `:636-776` | [`new_li_li_cu`](Geom2dGccCirc2d3TanIter::new_li_li_cu) |
//! | `(QualifiedLin, QCurve, QCurve)` | `:486-634` | [`new_li_cu_cu`](Geom2dGccCirc2d3TanIter::new_li_cu_cu) |
//! | `(QualifiedLin, QCurve, gp_Pnt2d)` | `:1037-1168` | [`new_li_cu_pt`](Geom2dGccCirc2d3TanIter::new_li_cu_pt) |
//! | `(QCurve, gp_Pnt2d, gp_Pnt2d)` | `:912-1034` | [`new_cu_pt_pt`](Geom2dGccCirc2d3TanIter::new_cu_pt_pt) |
//! | `(QCurve, QCurve, gp_Pnt2d)` | `:778-910` | [`new_cu_cu_pt`](Geom2dGccCirc2d3TanIter::new_cu_cu_pt) |
//! | `(QCurve, QCurve, QCurve)` | `:44-183` | [`new_cu_cu_cu`](Geom2dGccCirc2d3TanIter::new_cu_cu_cu) |
//! | `IsDone`..`IsTheSame3` | `:1455-1568` | this file |
//!
//! Translation notes:
//! - `Tol` is `|Tolerance|` for the analytic call and is then overwritten with
//!   `1.e-12` for the three `dot <= Tol` tests; the port keeps both uses,
//! - `Func.Value(Ufirst, Umin)` after `Root.Perform` only refills the vector
//!   OCCT hands in as output; the port's systems are stateless, so the call is
//!   dropped,
//! - `new_ci_li_cu` (`cxx:1228-1229`) calls `Func.Value` **before**
//!   `Root.Root`; harmless because `Func` is stateless,
//! - OCCT's `Geom2dGcc_FunctionTanCuCuCu::Values` returns `true` without
//!   writing the Jacobian (`cxx:275-278` of that file), so
//!   `math_FunctionSetRoot` runs on a matrix that is never refreshed; the port
//!   forwards the same empty `Values` and does not synthesize a Jacobian,
//! - `new_ci_ci_cu` (`cxx:455-459`) re-tests `Qualified1` against `R2` and never
//!   looks at `Qualified2`; `new_cu_cu_pt` (`cxx:887-889`) likewise re-tests
//!   `Qualified1` for the second operand. Both are OCCT bugs and are kept,
//! - `new_li_li_cu` (`cxx:707`) computes `pscal` from the **absolute** center
//!   coordinates and reuses it for `Qualified2`, unlike `new_li_cu_cu`
//!   (`cxx:556-560`) which measures from `L1.Location()`; kept,
//! - `new_li_cu_pt` (`cxx:1084`) and `new_ci_cu_pt` (`cxx:1381`) evaluate the
//!   curve operand at `Ufirst(2)` (the line / point-circle unknown) instead of
//!   `Ufirst(3)`; kept,
//! - `new_li_cu_pt` uses `tol(1) = 2.e-15` (`cxx:1084`) where the other circle
//!   operands use `2.e-15 * M_PI`,
//! - `new_ci_cu_pt` divides `Vec1.Dot(Tan1)` by `normevec1` alone and
//!   `Vec3.Dot(Tan3)` by `normevec3` alone (`cxx:1398-1418`), dropping the
//!   `normetan` factor; kept,
//! - `ThisSolution` does **not** check `WellDone`, although `hxx:147` says it
//!   raises, while `WhichQualifier` / `Tangency1` / `Tangency2` / `Tangency3`
//!   and the three `IsTheSame*` do,
//! - `IsTheSame1` returns `TheSame1`, `IsTheSame2` returns **`TheSame3`**
//!   (`cxx:1552`) and `IsTheSame3` returns `true` unconditionally
//!   (`cxx:1566`); all three are kept,
//! - `Tangency1` throws `StdFail_NotDone` both when `!WellDone` and when
//!   `TheSame1 != 0` (`cxx:1489-1498`),
//! - `GccEnt_BadQualifier` and `StdFail_NotDone` are panics, matching the rest
//!   of the `gcc` module.
//!
//! UNPORTED: none. The `.cxx` has no other member.

use std::f64::consts::PI;

use occt_core::gp::{GpAx22d, GpCirc2d, GpDir2d, GpPnt2d, GpVec2d};
use occt_core::precision::REAL_SMALL;
use occt_math::{MathFunctionSetRoot, MathVector};

use crate::gcc::function_tan_cu_cu_cu::Geom2dGccFunctionTanCuCuCu;
use crate::gcc::gcc_ent::GccEntPosition;

mod circ_curve_point;
mod circ_line_curve;
mod circ_two_circles;
mod circ_two_curves;
mod curve_two_points;
mod line_curve_point;
mod line_two_curves;
mod three_curves;
mod two_curves_point;
mod two_lines_curve;

/// `RealFirst()` (`Standard_Real.hxx:88-91`) = `-DBL_MAX`.
const REAL_FIRST: f64 = -f64::MAX;

/// `RealLast()` (`Standard_Real.hxx:94-97`) = `DBL_MAX`.
const REAL_LAST: f64 = f64::MAX;

/// `2.e-15 * M_PI`, the angular tolerance of every circle unknown except
/// `cxx:1085`.
const ANGULAR_TOL: f64 = 2.0e-15 * PI;

/// `1.e-15`, the linear tolerance of every line unknown.
const LINEAR_TOL: f64 = 1.0e-15;

/// `1.e-12`, the acceptance threshold assigned to `Tol` before the three dot
/// products are compared with it.
const ACCEPT_TOL: f64 = 1.0e-12;

/// `gp_Ax2d(P, gp_Dir2d(gp_Dir2d::D::X))`, the axis of every synthetic and
/// solution circle.
fn axis_at(center: GpPnt2d) -> GpAx22d {
    GpAx22d::from_xdir(center, GpDir2d::new(1.0, 0.0).expect("gp::DX2d()"))
}

/// `gp_Circ2d(gp_Ax2d(P, gp_Dir2d(gp_Dir2d::D::X)), 0.)`, the radius-0 circle
/// that turns a point into an operand (`cxx:809`, `:940-941`, `:1069`,
/// `:1352`).
fn point_circle(point: &GpPnt2d) -> GpCirc2d {
    GpCirc2d::new(axis_at(*point), 0.0)
}

/// `gp_Vec2d(P1, P2)` (`gp_Vec2d.cxx:83-90`).
fn vec_between(from: &GpPnt2d, to: &GpPnt2d) -> GpVec2d {
    GpVec2d::new(to.x() - from.x(), to.y() - from.y())
}

/// `Vec.Dot(Tan) / (normevec * normetan)` guarded by both
/// `>= gp::Resolution()` (`cxx:125-128` and its twins).
fn dot_ratio(vec: &GpVec2d, tan: &GpVec2d) -> f64 {
    let nv = vec.magnitude();
    let nt = tan.magnitude();
    if nv >= REAL_SMALL && nt >= REAL_SMALL {
        vec.dot(tan) / (nv * nt)
    } else {
        0.0
    }
}

/// `Vec.Angle(Tan)`; `gp_Vec2d::Angle` raises on a zero operand, which the
/// port's `angle` reports as `0.0`, the same convention as
/// `circ2d_2tan_on_iter::angle_or_zero`.
fn angle_or_zero(vec: &GpVec2d, tan: &GpVec2d) -> f64 {
    vec.angle(tan)
}

/// The `IsUnqualified() || (IsEnclosing() && Rsol >= R1 && dist <= Rsol) ||
/// (IsOutside() && dist >= Rsol) || (IsEnclosed() && Rsol <= R1 && dist <=
/// Rsol)` test of `cxx:452-454`, `:1284-1286`, `:1425-1427`.
fn accept_circ(q: GccEntPosition, r1: f64, r_sol: f64, dist: f64) -> bool {
    q.is_unqualified()
        || (q.is_enclosing() && r_sol >= r1 && dist <= r_sol)
        || (q.is_outside() && dist >= r_sol)
        || (q.is_enclosed() && r_sol <= r1 && dist <= r_sol)
}

/// The `IsUnqualified() || (IsOutside() && pscal <= 0.) || (IsEnclosed() &&
/// pscal >= 0.)` test of `cxx:600-602`, `:708-709`, `:749-750`, `:1144-1145`,
/// `:1289-1290`.
fn accept_lin(q: GccEntPosition, pscal: f64) -> bool {
    q.is_unqualified() || (q.is_outside() && pscal <= 0.0) || (q.is_enclosed() && pscal >= 0.0)
}

/// The `IsUnqualified() || (IsEnclosing() && angle <= 0.) || (IsOutside() &&
/// angle >= 0.) || (IsEnclosed() && angle <= 0.)` test of `cxx:152-167`,
/// `:461-463`, `:887-889`, `:1013-1015`, `:1147-1149`, `:1292-1294`,
/// `:1431-1433`.
fn accept_curve(q: GccEntPosition, angle: f64) -> bool {
    q.is_unqualified()
        || (q.is_enclosing() && angle <= 0.0)
        || (q.is_outside() && angle >= 0.0)
        || (q.is_enclosed() && angle <= 0.0)
}

/// `Geom2dGcc_Circ2d3TanIter` (`Geom2dGcc_Circ2d3TanIter.hxx:32-216`).
#[derive(Clone)]
pub struct Geom2dGccCirc2d3TanIter {
    /// `WellDone` (`hxx:196`).
    well_done: bool,
    /// `cirsol` (`hxx:197`).
    cirsol: GpCirc2d,
    /// `qualifier1` (`hxx:198`).
    qualifier1: GccEntPosition,
    /// `qualifier2` (`hxx:199`).
    qualifier2: GccEntPosition,
    /// `qualifier3` (`hxx:200`).
    qualifier3: GccEntPosition,
    /// `TheSame1` (`hxx:201`).
    the_same1: bool,
    /// `TheSame2` (`hxx:202`): written by every constructor and read nowhere,
    /// because OCCT's `IsTheSame2` returns `TheSame3` (`cxx:1552`).
    #[allow(dead_code)]
    the_same2: bool,
    /// `TheSame3` (`hxx:203`).
    the_same3: bool,
    /// `pnttg1sol` (`hxx:204`).
    pnttg1sol: GpPnt2d,
    /// `pnttg2sol` (`hxx:205`).
    pnttg2sol: GpPnt2d,
    /// `pnttg3sol` (`hxx:206`).
    pnttg3sol: GpPnt2d,
    /// `par1sol` (`hxx:207`).
    par1sol: f64,
    /// `par2sol` (`hxx:208`).
    par2sol: f64,
    /// `par3sol` (`hxx:209`).
    par3sol: f64,
    /// `pararg1` (`hxx:210`).
    pararg1: f64,
    /// `pararg2` (`hxx:211`).
    pararg2: f64,
    /// `pararg3` (`hxx:212`).
    pararg3: f64,
}

impl Geom2dGccCirc2d3TanIter {
    /// The member state every constructor starts from (`cxx:53-62` and twins).
    fn not_done() -> Self {
        Self {
            well_done: false,
            cirsol: GpCirc2d::default(),
            qualifier1: GccEntPosition::NoQualifier,
            qualifier2: GccEntPosition::NoQualifier,
            qualifier3: GccEntPosition::NoQualifier,
            the_same1: false,
            the_same2: false,
            the_same3: false,
            pnttg1sol: GpPnt2d::zero(),
            pnttg2sol: GpPnt2d::zero(),
            pnttg3sol: GpPnt2d::zero(),
            par1sol: 0.0,
            par2sol: 0.0,
            par3sol: 0.0,
            pararg1: 0.0,
            pararg2: 0.0,
            pararg3: 0.0,
        }
    }

    /// `math_FunctionSetRoot Root(Func, tol); Root.Perform(Func, Ufirst, Umin,
    /// Umax);` followed by `Root.Root(Ufirst)` (`cxx:98-104` and twins),
    /// answering `None` when `!Root.IsDone()`.
    fn solve(
        func: &mut Geom2dGccFunctionTanCuCuCu,
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

    /// The `qualifier1 = ...; ...; WellDone = true;` tails (`cxx:466-478` and
    /// twins).
    #[allow(clippy::too_many_arguments)]
    fn commit(
        &mut self,
        sol: GpCirc2d,
        quals: [GccEntPosition; 3],
        pnt1: GpPnt2d,
        pararg1: f64,
        par1sol: f64,
        pnt2: GpPnt2d,
        pararg2: f64,
        par2sol: f64,
        pnt3: GpPnt2d,
        pararg3: f64,
        par3sol: f64,
    ) {
        self.cirsol = sol;
        self.qualifier1 = quals[0];
        self.qualifier2 = quals[1];
        self.qualifier3 = quals[2];
        self.pnttg1sol = pnt1;
        self.pararg1 = pararg1;
        self.par1sol = par1sol;
        self.pnttg2sol = pnt2;
        self.pararg2 = pararg2;
        self.par2sol = par2sol;
        self.pnttg3sol = pnt3;
        self.pararg3 = pararg3;
        self.par3sol = par3sol;
        self.well_done = true;
    }

    /// `(QualifiedCirc, QualifiedCirc, QCurve)` (`cxx:335-484`).
    #[allow(clippy::too_many_arguments)]
    pub fn new_ci_ci_cu(
        qualified1: &crate::gcc::gcc_ent::GccEntQualifiedCirc,
        qualified2: &crate::gcc::gcc_ent::GccEntQualifiedCirc,
        qualified3: &crate::gcc::qualified_curve::Geom2dGccQCurve,
        param1: f64,
        param2: f64,
        param3: f64,
        tolerance: f64,
    ) -> Self {
        circ_two_circles::new(
            qualified1, qualified2, qualified3, param1, param2, param3, tolerance,
        )
    }

    /// `(QualifiedCirc, QCurve, QCurve)` (`cxx:187-332`).
    #[allow(clippy::too_many_arguments)]
    pub fn new_ci_cu_cu(
        qualified1: &crate::gcc::gcc_ent::GccEntQualifiedCirc,
        qualified2: &crate::gcc::qualified_curve::Geom2dGccQCurve,
        qualified3: &crate::gcc::qualified_curve::Geom2dGccQCurve,
        param1: f64,
        param2: f64,
        param3: f64,
        tolerance: f64,
    ) -> Self {
        circ_two_curves::new(
            qualified1, qualified2, qualified3, param1, param2, param3, tolerance,
        )
    }

    /// `(QualifiedCirc, QualifiedLin, QCurve)` (`cxx:1171-1315`).
    #[allow(clippy::too_many_arguments)]
    pub fn new_ci_li_cu(
        qualified1: &crate::gcc::gcc_ent::GccEntQualifiedCirc,
        qualified2: &crate::gcc::gcc_ent::GccEntQualifiedLin,
        qualified3: &crate::gcc::qualified_curve::Geom2dGccQCurve,
        param1: f64,
        param2: f64,
        param3: f64,
        tolerance: f64,
    ) -> Self {
        circ_line_curve::new(
            qualified1, qualified2, qualified3, param1, param2, param3, tolerance,
        )
    }

    /// `(QualifiedCirc, QCurve, gp_Pnt2d)` (`cxx:1318-1453`).
    #[allow(clippy::too_many_arguments)]
    pub fn new_ci_cu_pt(
        qualified1: &crate::gcc::gcc_ent::GccEntQualifiedCirc,
        qualified2: &crate::gcc::qualified_curve::Geom2dGccQCurve,
        point3: &GpPnt2d,
        param1: f64,
        param2: f64,
        tolerance: f64,
    ) -> Self {
        circ_curve_point::new(qualified1, qualified2, point3, param1, param2, tolerance)
    }

    /// `(QualifiedLin, QualifiedLin, QCurve)` (`cxx:636-776`).
    #[allow(clippy::too_many_arguments)]
    pub fn new_li_li_cu(
        qualified1: &crate::gcc::gcc_ent::GccEntQualifiedLin,
        qualified2: &crate::gcc::gcc_ent::GccEntQualifiedLin,
        qualified3: &crate::gcc::qualified_curve::Geom2dGccQCurve,
        param1: f64,
        param2: f64,
        param3: f64,
        tolerance: f64,
    ) -> Self {
        two_lines_curve::new(
            qualified1, qualified2, qualified3, param1, param2, param3, tolerance,
        )
    }

    /// `(QualifiedLin, QCurve, QCurve)` (`cxx:486-634`).
    #[allow(clippy::too_many_arguments)]
    pub fn new_li_cu_cu(
        qualified1: &crate::gcc::gcc_ent::GccEntQualifiedLin,
        qualified2: &crate::gcc::qualified_curve::Geom2dGccQCurve,
        qualified3: &crate::gcc::qualified_curve::Geom2dGccQCurve,
        param1: f64,
        param2: f64,
        param3: f64,
        tolerance: f64,
    ) -> Self {
        line_two_curves::new(
            qualified1, qualified2, qualified3, param1, param2, param3, tolerance,
        )
    }

    /// `(QualifiedLin, QCurve, gp_Pnt2d)` (`cxx:1037-1168`).
    #[allow(clippy::too_many_arguments)]
    pub fn new_li_cu_pt(
        qualified1: &crate::gcc::gcc_ent::GccEntQualifiedLin,
        qualified2: &crate::gcc::qualified_curve::Geom2dGccQCurve,
        point3: &GpPnt2d,
        param1: f64,
        param2: f64,
        tolerance: f64,
    ) -> Self {
        line_curve_point::new(qualified1, qualified2, point3, param1, param2, tolerance)
    }

    /// `(QCurve, gp_Pnt2d, gp_Pnt2d)` (`cxx:912-1034`).
    #[allow(clippy::too_many_arguments)]
    pub fn new_cu_pt_pt(
        qualified1: &crate::gcc::qualified_curve::Geom2dGccQCurve,
        point2: &GpPnt2d,
        point3: &GpPnt2d,
        param1: f64,
        tolerance: f64,
    ) -> Self {
        curve_two_points::new(qualified1, point2, point3, param1, tolerance)
    }

    /// `(QCurve, QCurve, gp_Pnt2d)` (`cxx:778-910`).
    #[allow(clippy::too_many_arguments)]
    pub fn new_cu_cu_pt(
        qualified1: &crate::gcc::qualified_curve::Geom2dGccQCurve,
        qualified2: &crate::gcc::qualified_curve::Geom2dGccQCurve,
        point3: &GpPnt2d,
        param1: f64,
        param2: f64,
        tolerance: f64,
    ) -> Self {
        two_curves_point::new(qualified1, qualified2, point3, param1, param2, tolerance)
    }

    /// `(QCurve, QCurve, QCurve)` (`cxx:44-183`).
    #[allow(clippy::too_many_arguments)]
    pub fn new_cu_cu_cu(
        qualified1: &crate::gcc::qualified_curve::Geom2dGccQCurve,
        qualified2: &crate::gcc::qualified_curve::Geom2dGccQCurve,
        qualified3: &crate::gcc::qualified_curve::Geom2dGccQCurve,
        param1: f64,
        param2: f64,
        param3: f64,
        tolerance: f64,
    ) -> Self {
        three_curves::new(
            qualified1, qualified2, qualified3, param1, param2, param3, tolerance,
        )
    }

    /// `IsDone()` (`cxx:1455-1458`).
    pub fn is_done(&self) -> bool {
        self.well_done
    }

    /// `ThisSolution()` (`cxx:1460-1463`): no `WellDone` check in the `.cxx`
    /// even though `hxx:147` documents `StdFail_NotDone`.
    pub fn this_solution(&self) -> GpCirc2d {
        self.cirsol
    }

    /// `WhichQualifier(Qualif1, Qualif2, Qualif3)` (`cxx:1465-1479`).
    pub fn which_qualifier(&self) -> (GccEntPosition, GccEntPosition, GccEntPosition) {
        if !self.well_done {
            panic!("StdFail_NotDone in Geom2dGcc_Circ2d3TanIter::WhichQualifier");
        }
        (self.qualifier1, self.qualifier2, self.qualifier3)
    }

    /// `Tangency1(ParSol, ParArg, PntSol)` (`cxx:1481-1500`).
    pub fn tangency1(&self) -> (f64, f64, GpPnt2d) {
        if !self.well_done || self.the_same1 {
            panic!("StdFail_NotDone in Geom2dGcc_Circ2d3TanIter::Tangency1");
        }
        (self.par1sol, self.pararg1, self.pnttg1sol)
    }

    /// `Tangency2(ParSol, ParArg, PntSol)` (`cxx:1502-1516`).
    pub fn tangency2(&self) -> (f64, f64, GpPnt2d) {
        if !self.well_done {
            panic!("StdFail_NotDone in Geom2dGcc_Circ2d3TanIter::Tangency2");
        }
        (self.par2sol, self.pararg2, self.pnttg2sol)
    }

    /// `Tangency3(ParSol, ParArg, PntSol)` (`cxx:1518-1528`).
    pub fn tangency3(&self) -> (f64, f64, GpPnt2d) {
        if !self.well_done {
            panic!("StdFail_NotDone in Geom2dGcc_Circ2d3TanIter::Tangency3");
        }
        (self.par3sol, self.pararg3, self.pnttg3sol)
    }

    /// `IsTheSame1()` (`cxx:1530-1542`).
    pub fn is_the_same1(&self) -> bool {
        if !self.well_done {
            panic!("StdFail_NotDone in Geom2dGcc_Circ2d3TanIter::IsTheSame1");
        }
        self.the_same1
    }

    /// `IsTheSame2()` (`cxx:1544-1557`): returns `TheSame3`, not `TheSame2`.
    pub fn is_the_same2(&self) -> bool {
        if !self.well_done {
            panic!("StdFail_NotDone in Geom2dGcc_Circ2d3TanIter::IsTheSame2");
        }
        self.the_same3
    }

    /// `IsTheSame3()` (`cxx:1559-1568`): `true` whenever `WellDone`.
    pub fn is_the_same3(&self) -> bool {
        if !self.well_done {
            panic!("StdFail_NotDone in Geom2dGcc_Circ2d3TanIter::IsTheSame3");
        }
        true
    }
}
