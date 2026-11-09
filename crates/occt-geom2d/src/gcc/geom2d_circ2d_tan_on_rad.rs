//! `Geom2dGcc_Circ2dTanOnRad`
//! (`Geom2dGcc_Circ2dTanOnRad.{hxx:44-190,cxx:33-356}`): circles of a given
//! radius tangent to one operand with the centre on a second one.
//!
//! The class is a front end: when both the operand and the carrier are lines or
//! circles it runs the analytic [`GccAnaCirc2dTanOnRad`], otherwise the
//! geometric [`Geom2dGccCirc2dTanOnRadGeo`], and it hands every solution back
//! through the `cirsol` / `qualifier1` / ... arrays of `hxx:176-190`.
//!
//! | OCCT | arguments | `cxx` |
//! |---|---|---|
//! | `(QualifiedCurve, Adaptor, Radius, Tol)` | `:33-153` |
//! | `(Geom2d_Point, Adaptor, Radius, Tol)` | `:155-217` |
//! | `Results(GccAna_Circ2dTanOnRad)` | `:219-236` |
//! | `Results(Geom2dGcc_Circ2dTanOnRadGeo)` | `:238-253` |
//! | `IsDone`..`IsTheSame1` | `:255-356` |
//!
//! Translation notes:
//! - `Radius < 0.` is tested before anything else and throws
//!   `Standard_NegativeValue` (`cxx:47-50`, `:168-172`), so the port panics
//!   before any type dispatch,
//! - the first constructor dispatches the geometric arm on `Type1` alone
//!   (`cxx:122-150`): a circle operand goes to
//!   `Circ2dTanOnRadGeo(QualifiedCirc, ...)`, a line operand to
//!   `(...QualifiedLin, ...)`, anything else to the `QCurve` overload,
//! - `WellDone` and `NbrSol` are assigned outside any `if (WellDone)` guard
//!   (`cxx:79-80`, `:86-87`, `:101-102`, `:108-109`, `:129-130`, `:138-139`,
//!   `:147-148`) and `Results` then indexes `1..NbrSol`, calling accessors that
//!   throw `StdFail_NotDone` when the solver failed; the port keeps that,
//! - `Results` calls `Tangency1` unconditionally (`cxx:232`, `:251`) and that
//!   accessor throws `StdFail_NotDone` when `TheSame1` is set (`cxx:311-320`),
//!   so such a solution aborts the constructor; the port panics accordingly,
//! - `Results(Geom2dGcc_Circ2dTanOnRadGeo)` takes the centre from the
//!   geometric solver's `CenterOn3`, which answers `pnttg1sol`
//!   (`Geom2dGcc_Circ2dTanOnRadGeo.cxx:890`), so `pntcen3` of the front end
//!   receives the tangency point; the port keeps that,
//! - `StdFail_NotDone` / `Standard_OutOfRange` are panics, matching the rest of
//!   the `gcc` module.
//!
//! UNPORTED: none - both constructors, both `Results` and all seven accessors.

use occt_core::gp::{GpCirc2d, GpLin2d, GpPnt2d};
use occt_core::kernel::geomabs::CurveType;

use crate::adaptor_curve::{basis2d, type2d, Geom2dAdaptorCurve};
use crate::gcc::circ2d_tan_on_rad::GccAnaCirc2dTanOnRad;
use crate::gcc::circ2d_tan_on_rad_geo::Geom2dGccCirc2dTanOnRadGeo;
use crate::gcc::gcc_ent::{GccEntPosition, GccEntQualifiedCirc, GccEntQualifiedLin};
use crate::gcc::qualified_curve::{Geom2dGccQCurve, Geom2dGccQualifiedCurve};

/// `cirsol(1, 8)` (`cxx:38`, `:158`): the front-end arrays hold eight slots
/// even though the analytic solver generates at most four.
pub const MAX_SOL: usize = 8;

/// `Geom2dGcc_Circ2dTanOnRad` (`Geom2dGcc_Circ2dTanOnRad.hxx:44-190`).
#[derive(Clone)]
pub struct Geom2dGccCirc2dTanOnRad {
    /// `WellDone` (`hxx:176`).
    well_done: bool,
    /// `cirsol` (`hxx:177`).
    cirsol: [GpCirc2d; MAX_SOL],
    /// `NbrSol` (`hxx:178`).
    nbr_sol: usize,
    /// `qualifier1` (`hxx:179`).
    qualifier1: [GccEntPosition; MAX_SOL],
    /// `TheSame1` (`hxx:180`).
    the_same1: [i32; MAX_SOL],
    /// `pnttg1sol` (`hxx:181`).
    pnttg1sol: [GpPnt2d; MAX_SOL],
    /// `par1sol` (`hxx:182`).
    par1sol: [f64; MAX_SOL],
    /// `pararg1` (`hxx:183`).
    pararg1: [f64; MAX_SOL],
    /// `pntcen3` (`hxx:184`).
    pntcen3: [GpPnt2d; MAX_SOL],
    /// `parcen3` (`hxx:185`).
    parcen3: [f64; MAX_SOL],
}

/// `Geom2d_Circle` cast of the operand (`cxx:69-70`, `:125-126`).
fn as_circle(q: &Geom2dGccQualifiedCurve) -> GpCirc2d {
    basis2d(q.qualified().as_ref())
        .gp_circ2d()
        .expect("Geom2d_Circle down_cast")
}

/// `Geom2d_Line` cast of the operand (`cxx:92-93`, `:135-136`).
fn as_line(q: &Geom2dGccQualifiedCurve) -> GpLin2d {
    basis2d(q.qualified().as_ref())
        .gp_lin2d()
        .expect("Geom2d_Line down_cast")
}

/// `GccEnt_QualifiedCirc(c1, Qualified1.Qualifier())` (`cxx:71`, `:127`).
fn as_qualified_circ(q: &Geom2dGccQualifiedCurve) -> GccEntQualifiedCirc {
    GccEntQualifiedCirc::new(as_circle(q), q.qualifier())
}

/// `GccEnt_QualifiedLin(l1, Qualified1.Qualifier())` (`cxx:95`, `:137`).
fn as_qualified_lin(q: &Geom2dGccQualifiedCurve) -> GccEntQualifiedLin {
    GccEntQualifiedLin::new(as_line(q), q.qualifier())
}

/// `Geom2dGcc_QCurve Qc1(C1, Qualified1.Qualifier())` (`cxx:146`).
fn as_qcurve(q: &Geom2dGccQualifiedCurve) -> Geom2dGccQCurve {
    Geom2dGccQCurve::new(
        Geom2dAdaptorCurve::from_curve(q.qualified().clone()),
        q.qualifier(),
    )
}

/// `C1.GetType()` (`cxx:54`).
fn curve_type(q: &Geom2dGccQualifiedCurve) -> CurveType {
    type2d(q.qualified().as_ref())
}

/// `Type == GeomAbs_Line || Type == GeomAbs_Circle`.
fn line_or_circle(t: CurveType) -> bool {
    t == CurveType::Line || t == CurveType::Circle
}

/// `TheSame1(j) = 1` / `TheSame1(j) = 0` (`cxx:226-233`, `:245-250`) on an
/// `int` array.
fn flag(value: bool) -> i32 {
    if value {
        1
    } else {
        0
    }
}

/// `CCon->Circ2d()` on the `OnCurve` adaptor (`cxx:75-76`, `:98-99`).
fn on_circle(on_curve: &Geom2dAdaptorCurve) -> GpCirc2d {
    on_curve
        .circle()
        .expect("Geom2d_Circle down_cast on OnCurve")
}

/// `LLon->Lin2d()` on the `OnCurve` adaptor (`cxx:83-84`, `:106-107`).
fn on_line(on_curve: &Geom2dAdaptorCurve) -> GpLin2d {
    on_curve
        .line()
        .expect("Geom2d_Line down_cast on OnCurve")
}

impl Geom2dGccCirc2dTanOnRad {
    /// The `: cirsol(1, 8), qualifier1(1, 8), ... :` initializer list both
    /// constructors share (`cxx:38-45`, `:158-166`).
    fn empty() -> Self {
        Self {
            well_done: false,
            cirsol: [GpCirc2d::default(); MAX_SOL],
            nbr_sol: 0,
            qualifier1: [GccEntPosition::NoQualifier; MAX_SOL],
            the_same1: [0; MAX_SOL],
            pnttg1sol: [GpPnt2d::zero(); MAX_SOL],
            par1sol: [0.0; MAX_SOL],
            pararg1: [0.0; MAX_SOL],
            pntcen3: [GpPnt2d::zero(); MAX_SOL],
            parcen3: [0.0; MAX_SOL],
        }
    }

    /// `Geom2dGcc_Circ2dTanOnRad(Qualified1, OnCurve, Radius, Tolerance)`
    /// (`cxx:33-153`).
    pub fn new(
        qualified1: &Geom2dGccQualifiedCurve,
        on_curve: &Geom2dAdaptorCurve,
        radius: f64,
        tolerance: f64,
    ) -> Self {
        if radius < 0.0 {
            panic!("Standard_NegativeValue in Geom2dGcc_Circ2dTanOnRad");
        }
        let mut s = Self::empty();
        let type1 = curve_type(qualified1);
        let type2 = on_curve.get_type().unwrap_or(CurveType::OtherCurve);
        if line_or_circle(type1) && line_or_circle(type2) {
            let solver = match (type1, type2) {
                (CurveType::Circle, CurveType::Circle) => {
                    GccAnaCirc2dTanOnRad::new_circ_on_circ(
                        &as_qualified_circ(qualified1),
                        &on_circle(on_curve),
                        radius,
                        tolerance,
                    )
                }
                (CurveType::Circle, CurveType::Line) => {
                    GccAnaCirc2dTanOnRad::new_circ_on_line(
                        &as_qualified_circ(qualified1),
                        &on_line(on_curve),
                        radius,
                        tolerance,
                    )
                }
                (CurveType::Line, CurveType::Circle) => {
                    GccAnaCirc2dTanOnRad::new_lin_on_circ(
                        &as_qualified_lin(qualified1),
                        &on_circle(on_curve),
                        radius,
                        tolerance,
                    )
                }
                (CurveType::Line, CurveType::Line) => {
                    GccAnaCirc2dTanOnRad::new_lin_on_line(
                        &as_qualified_lin(qualified1),
                        &on_line(on_curve),
                        radius,
                        tolerance,
                    )
                }
                _ => unreachable!("line_or_circle covers the four pairs"),
            };
            s.well_done = solver.is_done();
            s.nbr_sol = solver.nb_solutions();
            s.results_ana(&solver);
        } else {
            let geo = if type1 == CurveType::Circle {
                Geom2dGccCirc2dTanOnRadGeo::new_circ_on_curve(
                    &as_qualified_circ(qualified1),
                    on_curve,
                    radius,
                    tolerance,
                )
            } else if type1 == CurveType::Line {
                Geom2dGccCirc2dTanOnRadGeo::new_lin_on_curve(
                    &as_qualified_lin(qualified1),
                    on_curve,
                    radius,
                    tolerance,
                )
            } else {
                Geom2dGccCirc2dTanOnRadGeo::new_qcurve_on_curve(
                    &as_qcurve(qualified1),
                    on_curve,
                    radius,
                    tolerance,
                )
            };
            s.well_done = geo.is_done();
            s.nbr_sol = geo.nb_solutions();
            s.results_geo(&geo);
        }
        s
    }

    /// `Geom2dGcc_Circ2dTanOnRad(const Geom2d_Point&, OnCurve, Radius,
    /// Tolerance)` (`cxx:155-217`).
    pub fn new_with_point(
        point1: &GpPnt2d,
        on_curve: &Geom2dAdaptorCurve,
        radius: f64,
        tolerance: f64,
    ) -> Self {
        if radius < 0.0 {
            panic!("Standard_NegativeValue in Geom2dGcc_Circ2dTanOnRad");
        }
        let mut s = Self::empty();
        let type2 = on_curve.get_type().unwrap_or(CurveType::OtherCurve);
        if line_or_circle(type2) {
            let solver = if type2 == CurveType::Circle {
                GccAnaCirc2dTanOnRad::new_pnt_on_circ(
                    point1,
                    &on_circle(on_curve),
                    radius,
                    tolerance,
                )
            } else {
                GccAnaCirc2dTanOnRad::new_pnt_on_line(
                    point1,
                    &on_line(on_curve),
                    radius,
                    tolerance,
                )
            };
            s.well_done = solver.is_done();
            s.nbr_sol = solver.nb_solutions();
            s.results_ana(&solver);
        } else {
            let geo = Geom2dGccCirc2dTanOnRadGeo::new_pnt_on_curve(
                point1, on_curve, radius, tolerance,
            );
            s.well_done = geo.is_done();
            s.nbr_sol = geo.nb_solutions();
            s.results_geo(&geo);
        }
        s
    }

    /// `Results(const GccAna_Circ2dTanOnRad&)` (`cxx:219-236`).
    fn results_ana(&mut self, circ: &GccAnaCirc2dTanOnRad) {
        for j in 1..=self.nbr_sol {
            self.cirsol[j - 1] = circ.this_solution(j);
            self.the_same1[j - 1] = flag(circ.is_the_same1(j));
            let (par_sol, par_arg, pnt_tg) = circ.tangency1(j);
            self.par1sol[j - 1] = par_sol;
            self.pararg1[j - 1] = par_arg;
            self.pnttg1sol[j - 1] = pnt_tg;
            let (par_cen, pnt_cen) = circ.center_on3(j);
            self.parcen3[j - 1] = par_cen;
            self.pntcen3[j - 1] = pnt_cen;
            self.qualifier1[j - 1] = circ.which_qualifier(j);
        }
    }

    /// `Results(const Geom2dGcc_Circ2dTanOnRadGeo&)` (`cxx:238-253`).
    fn results_geo(&mut self, circ: &Geom2dGccCirc2dTanOnRadGeo) {
        for j in 1..=self.nbr_sol {
            self.cirsol[j - 1] = circ.this_solution(j);
            self.the_same1[j - 1] = flag(circ.is_the_same1(j));
            let (par_sol, par_arg, pnt_tg) = circ.tangency1(j);
            self.par1sol[j - 1] = par_sol;
            self.pararg1[j - 1] = par_arg;
            self.pnttg1sol[j - 1] = pnt_tg;
            let (par_cen, pnt_cen) = circ.center_on3(j);
            self.parcen3[j - 1] = par_cen;
            self.pntcen3[j - 1] = pnt_cen;
            self.qualifier1[j - 1] = circ.which_qualifier(j);
        }
    }

    /// `IsDone()` (`cxx:255-258`).
    pub fn is_done(&self) -> bool {
        self.well_done
    }

    /// `NbSolutions()` (`cxx:261-264`).
    pub fn nb_solutions(&self) -> usize {
        self.nbr_sol
    }

    /// `ThisSolution(Index)` (`cxx:267-277`).
    pub fn this_solution(&self, index: usize) -> GpCirc2d {
        self.check_done();
        self.check_index(index);
        self.cirsol[index - 1]
    }

    /// `WhichQualifier(Index, Qualif1)` (`cxx:280-295`).
    pub fn which_qualifier(&self, index: usize) -> GccEntPosition {
        self.check_done();
        self.check_index(index);
        self.qualifier1[index - 1]
    }

    /// `Tangency1(Index, ParSol, ParArg, PntSol)` (`cxx:298-322`): returns
    /// `None` where OCCT throws `StdFail_NotDone` because `TheSame1` is set.
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

    /// `CenterOn3(Index, ParArg, PntSol)` (`cxx:324-339`).
    pub fn center_on3(&self, index: usize) -> (f64, GpPnt2d) {
        self.check_done();
        self.check_index(index);
        (self.parcen3[index - 1], self.pntcen3[index - 1])
    }

    /// `IsTheSame1(Index)` (`cxx:341-356`).
    pub fn is_the_same1(&self, index: usize) -> bool {
        self.check_done();
        self.check_index(index);
        self.the_same1[index - 1] != 0
    }

    /// `if (!WellDone) throw StdFail_NotDone()`.
    fn check_done(&self) {
        assert!(
            self.well_done,
            "StdFail_NotDone in Geom2dGcc_Circ2dTanOnRad"
        );
    }

    /// `if (Index <= 0 || Index > NbrSol) throw Standard_OutOfRange()`.
    fn check_index(&self, index: usize) {
        assert!(
            index >= 1 && index <= self.nbr_sol,
            "Standard_OutOfRange in Geom2dGcc_Circ2dTanOnRad"
        );
    }
}
