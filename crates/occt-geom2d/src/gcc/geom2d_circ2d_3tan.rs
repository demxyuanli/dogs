//! `Geom2dGcc_Circ2d3Tan` (`Geom2dGcc_Circ2d3Tan.{hxx:52-229,cxx:32-830}`):
//! circles tangent to three operands.
//!
//! The class is a front end: it dispatches on the three operand types between
//! the analytic [`GccAnaCirc2d3Tan`] (all three operands are lines or circles)
//! and the iterative [`Geom2dGccCirc2d3TanIter`] (at least one operand is a
//! general curve), then hands every solution back through the `cirsol` /
//! `qualifier1` / ... arrays of `hxx:210-229`.
//!
//! | OCCT | arguments | `cxx` |
//! |---|---|---|
//! | `(QualifiedCurve, QualifiedCurve, QualifiedCurve, Tol, Param1, Param2, Param3)` | `:32-264` |
//! | `(QualifiedCurve, QualifiedCurve, Geom2d_Point, Tol, Param1, Param2)` | `:266-415` |
//! | `(QualifiedCurve, Geom2d_Point, Geom2d_Point, Tol, Param1)` | `:417-517` |
//! | `(Geom2d_Point, Geom2d_Point, Geom2d_Point, Tol)` | `:519-556` |
//! | `Results` | `:558-650` |
//! | `IsDone`..`IsTheSame3` | `:652-830` |
//!
//! Translation notes:
//! - the iterative arms always build three `Geom2dGcc_QCurve` and call the
//!   `(QCurve, QCurve, QCurve)` / `(QCurve, QCurve, gp_Pnt2d)` /
//!   `(QCurve, gp_Pnt2d, gp_Pnt2d)` overloads (`cxx:226-229`, `:376-379`,
//!   `:480-483`); the mixed `QualifiedCirc` / `QualifiedLin` overloads of
//!   `Geom2dGcc_Circ2d3TanIter` are never reached from this front end and are
//!   kept on the solver for its own callers,
//! - the iterative arms set `NbrSol = 1` *before* `if (WellDone)` (`cxx:231`,
//!   `:380`, `:484`), so a failed solve still reports one solution while
//!   `IsDone()` is false; kept as written,
//! - `Results` (`cxx:558-650`) only writes a tangency triple when the selected
//!   `TheSame` flag is 0; the other arrays keep their initial value. OCCT
//!   leaves those `TColgp_Array1OfPnt2d` entries default-constructed, the port
//!   leaves the zeros of [`Geom2dGccCirc2d3Tan::empty`],
//! - `Results` asserts nothing but calls `Tangency1/2/3`, which throw
//!   `StdFail_NotDone` when `TheSame` is set; `rank_value` only calls them when
//!   the flag is 0, since [`GccAnaCirc2d3Tan::tangency1`] answers `None` in
//!   exactly that case,
//! - the `(Line, Circle)` arms build the solver with the operands swapped and
//!   map the results back: `WhichQualifier(i, qualifier2, qualifier1,
//!   qualifier3)` with `Results(Circ, 2, 1, 3)` (`cxx:184`, `:355`). One arm is
//!   inconsistent in OCCT itself: `(Line, Circle, Circle)` passes
//!   `WhichQualifier(i, qualifier3, qualifier1, qualifier2)` (`cxx:169`) while
//!   `Results` uses `(3, 1, 2)` (`cxx:171`); the two permutations disagree by a
//!   rotation and both are reproduced here,
//! - `Geom2dGcc_QualifiedCurve` only stores the curve in this port, so
//!   `Geom2dAdaptor_Curve C1 = Qualified1.Qualified()` (`cxx:56-58`) becomes
//!   [`Geom2dAdaptorCurve::from_curve`] on that curve, exactly like the
//!   `Geom2dAdaptor_Curve(const Handle(Geom2d_Curve)&)` constructor OCCT runs
//!   when the front end takes a `QualifiedCurve`,
//! - `StdFail_NotDone` / `Standard_OutOfRange` are panics, matching the rest of
//!   the `gcc` module.
//!
//! UNPORTED: none - all four constructors, `Results` and all ten accessors.

use occt_core::gp::{GpCirc2d, GpLin2d, GpPnt2d};
use occt_core::kernel::geomabs::CurveType;

use crate::adaptor_curve::{basis2d, type2d, Geom2dAdaptorCurve};
use crate::gcc::circ2d_3tan::{GccAnaCirc2d3Tan, MAX_SOL};
use crate::gcc::circ2d_3tan_iter::Geom2dGccCirc2d3TanIter;
use crate::gcc::gcc_ent::{GccEntPosition, GccEntQualifiedCirc, GccEntQualifiedLin};
use crate::gcc::qualified_curve::{Geom2dGccQCurve, Geom2dGccQualifiedCurve};

/// `Geom2dGcc_Circ2d3Tan` (`Geom2dGcc_Circ2d3Tan.hxx:52-229`).
#[derive(Clone)]
pub struct Geom2dGccCirc2d3Tan {
    /// `WellDone` (`hxx:213`).
    well_done: bool,
    /// `cirsol` (`hxx:211`).
    cirsol: [GpCirc2d; MAX_SOL],
    /// `NbrSol` (`hxx:212`).
    nbr_sol: usize,
    /// `qualifier1` (`hxx:214`).
    qualifier1: [GccEntPosition; MAX_SOL],
    /// `qualifier2` (`hxx:215`).
    qualifier2: [GccEntPosition; MAX_SOL],
    /// `qualifier3` (`hxx:216`).
    qualifier3: [GccEntPosition; MAX_SOL],
    /// `TheSame1` (`hxx:217`).
    the_same1: [i32; MAX_SOL],
    /// `TheSame2` (`hxx:218`).
    the_same2: [i32; MAX_SOL],
    /// `TheSame3` (`hxx:219`).
    the_same3: [i32; MAX_SOL],
    /// `pnttg1sol` (`hxx:220`).
    pnttg1sol: [GpPnt2d; MAX_SOL],
    /// `pnttg2sol` (`hxx:221`).
    pnttg2sol: [GpPnt2d; MAX_SOL],
    /// `pnttg3sol` (`hxx:222`).
    pnttg3sol: [GpPnt2d; MAX_SOL],
    /// `par1sol` (`hxx:223`).
    par1sol: [f64; MAX_SOL],
    /// `par2sol` (`hxx:224`).
    par2sol: [f64; MAX_SOL],
    /// `par3sol` (`hxx:225`).
    par3sol: [f64; MAX_SOL],
    /// `pararg1` (`hxx:226`).
    pararg1: [f64; MAX_SOL],
    /// `pararg2` (`hxx:227`).
    pararg2: [f64; MAX_SOL],
    /// `pararg3` (`hxx:228`).
    pararg3: [f64; MAX_SOL],
}

/// `Geom2d_Circle` / `Geom2d_Line` casts of the front end (`cxx:77`, `:83`,
/// `:88`, `:103`, `:118`, `:122`, `:136`, `:152`, `:158`, `:163`, `:175`,
/// `:193`, `:197`, `:307`, `:312`, `:324`, `:342`, `:348`, `:361`, `:452`,
/// `:466`).
fn as_circle(q: &Geom2dGccQualifiedCurve) -> GpCirc2d {
    basis2d(q.qualified().as_ref())
        .gp_circ2d()
        .expect("Geom2d_Circle down_cast")
}

fn as_line(q: &Geom2dGccQualifiedCurve) -> GpLin2d {
    basis2d(q.qualified().as_ref())
        .gp_lin2d()
        .expect("Geom2d_Line down_cast")
}

fn as_qualified_circ(q: &Geom2dGccQualifiedCurve) -> GccEntQualifiedCirc {
    GccEntQualifiedCirc::new(as_circle(q), q.qualifier())
}

fn as_qualified_lin(q: &Geom2dGccQualifiedCurve) -> GccEntQualifiedLin {
    GccEntQualifiedLin::new(as_line(q), q.qualifier())
}

/// `Geom2dGcc_QCurve Qc1(C1, Qualified1.Qualifier())` (`cxx:226-228`,
/// `:376-378`, `:480-481`).
fn as_qcurve(q: &Geom2dGccQualifiedCurve) -> Geom2dGccQCurve {
    Geom2dGccQCurve::new(
        Geom2dAdaptorCurve::from_curve(q.qualified().clone()),
        q.qualifier(),
    )
}

/// `C1.GetType()` (`cxx:64-66`, `:295-296`, `:440`).
fn curve_type(q: &Geom2dGccQualifiedCurve) -> CurveType {
    type2d(q.qualified().as_ref())
}

/// `TypeN == GeomAbs_Line || TypeN == GeomAbs_Circle` (`cxx:71-73`,
/// `:302-303`, `:447`).
fn line_or_circle(t: CurveType) -> bool {
    t == CurveType::Line || t == CurveType::Circle
}

/// `TheSameK(j) = 1` / `TheSameK(j) = 0` (`cxx:565-575`) on an `int` array.
fn flag(value: bool) -> i32 {
    if value {
        1
    } else {
        0
    }
}

/// One rank of `Results` (`cxx:577-647`): the `TheSameK` flag of the selected
/// solver tangency and that tangency (`None` where OCCT would throw).
fn rank_value(
    circ: &GccAnaCirc2d3Tan,
    index: usize,
    rank: usize,
) -> (i32, Option<(f64, f64, GpPnt2d)>) {
    match rank {
        1 => (flag(circ.is_the_same1(index)), circ.tangency1(index)),
        2 => (flag(circ.is_the_same2(index)), circ.tangency2(index)),
        _ => (flag(circ.is_the_same3(index)), circ.tangency3(index)),
    }
}

impl Geom2dGccCirc2d3Tan {
    /// The `: cirsol(1, 20), qualifier1(1, 20), ... :` initializer list of the
    /// widest overload (`cxx:274-288`); `WellDone` and `NbrSol` are not on it.
    fn empty() -> Self {
        Self {
            well_done: false,
            cirsol: [GpCirc2d::default(); MAX_SOL],
            nbr_sol: 0,
            qualifier1: [GccEntPosition::NoQualifier; MAX_SOL],
            qualifier2: [GccEntPosition::NoQualifier; MAX_SOL],
            qualifier3: [GccEntPosition::NoQualifier; MAX_SOL],
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

    /// The analytic arm shared by every overload (`cxx:90-98`, `:105-111`,
    /// `:125-131`, `:137-145`, `:164-171`, `:178-185`, `:197-205`, `:213-219`,
    /// `:314-323`, `:328-336`, `:349-356`, `:364-370`, `:455-461`, `:468-475`,
    /// `:546-553`): `WellDone`, `NbrSol`, the `WhichQualifier` loop and
    /// `Results`.
    ///
    /// `qual_order[k]` is the solver-argument index whose qualifier lands in
    /// facade slot `k + 1`; `rank1..3` are the `Results` ranks.
    fn copy_ana(
        &mut self,
        circ: &GccAnaCirc2d3Tan,
        qual_order: [usize; 3],
        rank1: usize,
        rank2: usize,
        rank3: usize,
    ) {
        self.well_done = circ.is_done();
        self.nbr_sol = circ.nb_solutions();
        for i in 1..=self.nbr_sol {
            let (s1, s2, s3) = circ.which_qualifier(i);
            let solver = [s1, s2, s3];
            for k in 0..3 {
                let q = solver[qual_order[k] - 1];
                match k {
                    0 => self.qualifier1[i - 1] = q,
                    1 => self.qualifier2[i - 1] = q,
                    _ => self.qualifier3[i - 1] = q,
                }
            }
        }
        self.results(circ, rank1, rank2, rank3);
    }

    /// `Results(Circ, Rank1, Rank2, Rank3)` (`cxx:558-650`).
    fn results(
        &mut self,
        circ: &GccAnaCirc2d3Tan,
        rank1: usize,
        rank2: usize,
        rank3: usize,
    ) {
        for j in 1..=self.nbr_sol {
            self.cirsol[j - 1] = circ.this_solution(j);
            let (same, tangency) = rank_value(circ, j, rank1);
            self.the_same1[j - 1] = same;
            if let Some((par_sol, par_arg, pnt)) = tangency {
                self.par1sol[j - 1] = par_sol;
                self.pararg1[j - 1] = par_arg;
                self.pnttg1sol[j - 1] = pnt;
            }
            let (same, tangency) = rank_value(circ, j, rank2);
            self.the_same2[j - 1] = same;
            if let Some((par_sol, par_arg, pnt)) = tangency {
                self.par2sol[j - 1] = par_sol;
                self.pararg2[j - 1] = par_arg;
                self.pnttg2sol[j - 1] = pnt;
            }
            let (same, tangency) = rank_value(circ, j, rank3);
            self.the_same3[j - 1] = same;
            if let Some((par_sol, par_arg, pnt)) = tangency {
                self.par3sol[j - 1] = par_sol;
                self.pararg3[j - 1] = par_arg;
                self.pnttg3sol[j - 1] = pnt;
            }
        }
    }

    /// The iterative arm shared by the three non-analytic overloads
    /// (`cxx:224-264`, `:374-414`, `:478-516`): one solution, the three
    /// `IsTheSame` flags, the three tangencies and `WhichQualifier`.
    fn copy_iter(&mut self, iter: &Geom2dGccCirc2d3TanIter) {
        self.well_done = iter.is_done();
        self.nbr_sol = 1;
        if self.well_done {
            self.cirsol[0] = iter.this_solution();
            self.the_same1[0] = flag(iter.is_the_same1());
            self.the_same2[0] = flag(iter.is_the_same2());
            self.the_same3[0] = flag(iter.is_the_same3());
            let (par_sol, par_arg, pnt) = iter.tangency1();
            self.par1sol[0] = par_sol;
            self.pararg1[0] = par_arg;
            self.pnttg1sol[0] = pnt;
            let (par_sol, par_arg, pnt) = iter.tangency2();
            self.par2sol[0] = par_sol;
            self.pararg2[0] = par_arg;
            self.pnttg2sol[0] = pnt;
            let (par_sol, par_arg, pnt) = iter.tangency3();
            self.par3sol[0] = par_sol;
            self.pararg3[0] = par_arg;
            self.pnttg3sol[0] = pnt;
            let (q1, q2, q3) = iter.which_qualifier();
            self.qualifier1[0] = q1;
            self.qualifier2[0] = q2;
            self.qualifier3[0] = q3;
        }
    }

    /// `Geom2dGcc_Circ2d3Tan(Qualified1, Qualified2, Qualified3, Tolerance,
    /// Param1, Param2, Param3)` (`cxx:32-264`).
    pub fn new_3curves(
        qualified1: &Geom2dGccQualifiedCurve,
        qualified2: &Geom2dGccQualifiedCurve,
        qualified3: &Geom2dGccQualifiedCurve,
        tolerance: f64,
        param1: f64,
        param2: f64,
        param3: f64,
    ) -> Self {
        let mut s = Self::empty();
        let type1 = curve_type(qualified1);
        let type2 = curve_type(qualified2);
        let type3 = curve_type(qualified3);
        s.nbr_sol = 0;
        if line_or_circle(type1) && line_or_circle(type2) && line_or_circle(type3) {
            // `cxx:75-223`: the analytic solver, in the eight type orders OCCT
            // writes them.
            match (type1, type2, type3) {
                (CurveType::Circle, CurveType::Circle, CurveType::Circle) => {
                    let solver = GccAnaCirc2d3Tan::new_circ_circ_circ(
                        &as_qualified_circ(qualified1),
                        &as_qualified_circ(qualified2),
                        &as_qualified_circ(qualified3),
                        tolerance,
                    );
                    s.copy_ana(&solver, [1, 2, 3], 1, 2, 3);
                }
                (CurveType::Circle, CurveType::Circle, CurveType::Line) => {
                    let solver = GccAnaCirc2d3Tan::new_circ_circ_lin(
                        &as_qualified_circ(qualified1),
                        &as_qualified_circ(qualified2),
                        &as_qualified_lin(qualified3),
                        tolerance,
                    );
                    s.copy_ana(&solver, [1, 2, 3], 1, 2, 3);
                }
                (CurveType::Circle, CurveType::Line, CurveType::Circle) => {
                    let solver = GccAnaCirc2d3Tan::new_circ_circ_lin(
                        &as_qualified_circ(qualified1),
                        &as_qualified_circ(qualified3),
                        &as_qualified_lin(qualified2),
                        tolerance,
                    );
                    s.copy_ana(&solver, [1, 3, 2], 1, 3, 2);
                }
                (CurveType::Circle, CurveType::Line, CurveType::Line) => {
                    let solver = GccAnaCirc2d3Tan::new_circ_lin_lin(
                        &as_qualified_circ(qualified1),
                        &as_qualified_lin(qualified2),
                        &as_qualified_lin(qualified3),
                        tolerance,
                    );
                    s.copy_ana(&solver, [1, 2, 3], 1, 2, 3);
                }
                (CurveType::Line, CurveType::Circle, CurveType::Circle) => {
                    let solver = GccAnaCirc2d3Tan::new_circ_circ_lin(
                        &as_qualified_circ(qualified2),
                        &as_qualified_circ(qualified3),
                        &as_qualified_lin(qualified1),
                        tolerance,
                    );
                    // `cxx:169` passes `(qualifier3, qualifier1, qualifier2)`
                    // while `cxx:171` passes the ranks `(3, 1, 2)`; both are
                    // kept, even though the permutations disagree.
                    s.copy_ana(&solver, [2, 3, 1], 3, 1, 2);
                }
                (CurveType::Line, CurveType::Circle, CurveType::Line) => {
                    let solver = GccAnaCirc2d3Tan::new_circ_lin_lin(
                        &as_qualified_circ(qualified2),
                        &as_qualified_lin(qualified1),
                        &as_qualified_lin(qualified3),
                        tolerance,
                    );
                    s.copy_ana(&solver, [2, 1, 3], 2, 1, 3);
                }
                (CurveType::Line, CurveType::Line, CurveType::Circle) => {
                    let solver = GccAnaCirc2d3Tan::new_circ_lin_lin(
                        &as_qualified_circ(qualified3),
                        &as_qualified_lin(qualified2),
                        &as_qualified_lin(qualified1),
                        tolerance,
                    );
                    s.copy_ana(&solver, [3, 2, 1], 3, 2, 1);
                }
                (CurveType::Line, CurveType::Line, CurveType::Line) => {
                    let solver = GccAnaCirc2d3Tan::new_lin_lin_lin(
                        &as_qualified_lin(qualified1),
                        &as_qualified_lin(qualified2),
                        &as_qualified_lin(qualified3),
                    );
                    s.copy_ana(&solver, [1, 2, 3], 1, 2, 3);
                }
                _ => unreachable!("the three types are line or circle"),
            }
        } else {
            // `cxx:224-264`: the geometric solver.
            let iter = Geom2dGccCirc2d3TanIter::new_cu_cu_cu(
                &as_qcurve(qualified1),
                &as_qcurve(qualified2),
                &as_qcurve(qualified3),
                param1,
                param2,
                param3,
                tolerance,
            );
            s.copy_iter(&iter);
        }
        s
    }

    /// `Geom2dGcc_Circ2d3Tan(Qualified1, Qualified2, Point, Tolerance, Param1,
    /// Param2)` (`cxx:266-415`).
    pub fn new_2curves_pnt(
        qualified1: &Geom2dGccQualifiedCurve,
        qualified2: &Geom2dGccQualifiedCurve,
        point: &GpPnt2d,
        tolerance: f64,
        param1: f64,
        param2: f64,
    ) -> Self {
        let mut s = Self::empty();
        let type1 = curve_type(qualified1);
        let type2 = curve_type(qualified2);
        s.nbr_sol = 0;
        if line_or_circle(type1) && line_or_circle(type2) {
            // `cxx:304-373`: the analytic solver.
            match (type1, type2) {
                (CurveType::Circle, CurveType::Circle) => {
                    let solver = GccAnaCirc2d3Tan::new_circ_circ_pnt(
                        &as_qualified_circ(qualified1),
                        &as_qualified_circ(qualified2),
                        point,
                        tolerance,
                    );
                    s.copy_ana(&solver, [1, 2, 3], 1, 2, 3);
                }
                (CurveType::Circle, CurveType::Line) => {
                    let solver = GccAnaCirc2d3Tan::new_circ_lin_pnt(
                        &as_qualified_circ(qualified1),
                        &as_qualified_lin(qualified2),
                        point,
                        tolerance,
                    );
                    s.copy_ana(&solver, [1, 2, 3], 1, 2, 3);
                }
                (CurveType::Line, CurveType::Circle) => {
                    let solver = GccAnaCirc2d3Tan::new_circ_lin_pnt(
                        &as_qualified_circ(qualified2),
                        &as_qualified_lin(qualified1),
                        point,
                        tolerance,
                    );
                    s.copy_ana(&solver, [2, 1, 3], 2, 1, 3);
                }
                (CurveType::Line, CurveType::Line) => {
                    let solver = GccAnaCirc2d3Tan::new_lin_lin_pnt(
                        &as_qualified_lin(qualified1),
                        &as_qualified_lin(qualified2),
                        point,
                        tolerance,
                    );
                    s.copy_ana(&solver, [1, 2, 3], 1, 2, 3);
                }
                _ => unreachable!("the two types are line or circle"),
            }
        } else {
            // `cxx:374-414`: the geometric solver.
            let iter = Geom2dGccCirc2d3TanIter::new_cu_cu_pt(
                &as_qcurve(qualified1),
                &as_qcurve(qualified2),
                point,
                param1,
                param2,
                tolerance,
            );
            s.copy_iter(&iter);
        }
        s
    }

    /// `Geom2dGcc_Circ2d3Tan(Qualified1, Point1, Point2, Tolerance, Param1)`
    /// (`cxx:417-517`).
    pub fn new_1curve_2pnts(
        qualified1: &Geom2dGccQualifiedCurve,
        point1: &GpPnt2d,
        point2: &GpPnt2d,
        tolerance: f64,
        param1: f64,
    ) -> Self {
        let mut s = Self::empty();
        let type1 = curve_type(qualified1);
        s.nbr_sol = 0;
        if line_or_circle(type1) {
            // `cxx:448-476`: the analytic solver.
            match type1 {
                CurveType::Circle => {
                    let solver = GccAnaCirc2d3Tan::new_circ_pnt_pnt(
                        &as_qualified_circ(qualified1),
                        point1,
                        point2,
                        tolerance,
                    );
                    s.copy_ana(&solver, [1, 2, 3], 1, 2, 3);
                }
                _ => {
                    let solver = GccAnaCirc2d3Tan::new_lin_pnt_pnt(
                        &as_qualified_lin(qualified1),
                        point1,
                        point2,
                        tolerance,
                    );
                    s.copy_ana(&solver, [1, 2, 3], 1, 2, 3);
                }
            }
        } else {
            // `cxx:478-516`: the geometric solver.
            let iter = Geom2dGccCirc2d3TanIter::new_cu_pt_pt(
                &as_qcurve(qualified1),
                point1,
                point2,
                param1,
                tolerance,
            );
            s.copy_iter(&iter);
        }
        s
    }

    /// `Geom2dGcc_Circ2d3Tan(Point1, Point2, Point3, Tolerance)`
    /// (`cxx:519-556`).
    pub fn new_3pnts(
        point1: &GpPnt2d,
        point2: &GpPnt2d,
        point3: &GpPnt2d,
        tolerance: f64,
    ) -> Self {
        let mut s = Self::empty();
        s.nbr_sol = 0;
        let solver = GccAnaCirc2d3Tan::new_pnt_pnt_pnt(point1, point2, point3, tolerance);
        s.copy_ana(&solver, [1, 2, 3], 1, 2, 3);
        s
    }

    /// `IsDone()` (`cxx:652-655`).
    pub fn is_done(&self) -> bool {
        self.well_done
    }

    /// `NbSolutions()` (`cxx:657-660`).
    pub fn nb_solutions(&self) -> usize {
        self.nbr_sol
    }

    /// `ThisSolution(Index)` (`cxx:662-674`).
    pub fn this_solution(&self, index: usize) -> GpCirc2d {
        self.check(index);
        self.cirsol[index - 1]
    }

    /// `WhichQualifier(Index, Qualif1, Qualif2, Qualif3)` (`cxx:676-694`).
    pub fn which_qualifier(
        &self,
        index: usize,
    ) -> (GccEntPosition, GccEntPosition, GccEntPosition) {
        self.check(index);
        (
            self.qualifier1[index - 1],
            self.qualifier2[index - 1],
            self.qualifier3[index - 1],
        )
    }

    /// `Tangency1(Index, ParSol, ParArg, PntSol)` (`cxx:696-723`).
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

    /// `Tangency2(Index, ParSol, ParArg, PntSol)` (`cxx:725-751`).
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

    /// `Tangency3(Index, ParSol, ParArg, PntSol)` (`cxx:753-779`).
    pub fn tangency3(&self, index: usize) -> Option<(f64, f64, GpPnt2d)> {
        self.check(index);
        if self.the_same3[index - 1] != 0 {
            return None;
        }
        Some((
            self.par3sol[index - 1],
            self.pararg3[index - 1],
            self.pnttg3sol[index - 1],
        ))
    }

    /// `IsTheSame1(Index)` (`cxx:781-796`).
    pub fn is_the_same1(&self, index: usize) -> bool {
        self.check(index);
        self.the_same1[index - 1] != 0
    }

    /// `IsTheSame2(Index)` (`cxx:798-813`).
    pub fn is_the_same2(&self, index: usize) -> bool {
        self.check(index);
        self.the_same2[index - 1] != 0
    }

    /// `IsTheSame3(Index)` (`cxx:815-830`).
    pub fn is_the_same3(&self, index: usize) -> bool {
        self.check(index);
        self.the_same3[index - 1] != 0
    }

    /// `if (!WellDone) throw StdFail_NotDone()` and
    /// `if (Index <= 0 || Index > NbrSol) throw Standard_OutOfRange()`, the
    /// guard of every accessor above (`cxx:664-672`, `:680-688`, ...).
    fn check(&self, index: usize) {
        assert!(
            self.well_done,
            "StdFail_NotDone in Geom2dGcc_Circ2d3Tan"
        );
        assert!(
            index >= 1 && index <= self.nbr_sol,
            "Standard_OutOfRange in Geom2dGcc_Circ2d3Tan"
        );
    }
}
