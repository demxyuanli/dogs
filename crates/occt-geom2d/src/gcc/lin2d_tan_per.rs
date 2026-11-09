//! `GccAna_Lin2dTanPer` (`GccAna_Lin2dTanPer.cxx:40-366`): lines orthogonal
//! to a reference line (or circle), through a point or tangent to a circle.
//!
//! | OCCT | lines |
//! |---|---|
//! | `GccAna_Lin2dTanPer(gp_Pnt2d, gp_Lin2d)` | `cxx:40-79` |
//! | `GccAna_Lin2dTanPer(gp_Pnt2d, gp_Circ2d)` | `cxx:83-120` |
//! | `GccAna_Lin2dTanPer(GccEnt_QualifiedCirc, gp_Lin2d)` | `cxx:129-199` |
//! | `GccAna_Lin2dTanPer(GccEnt_QualifiedCirc, gp_Circ2d)` | `cxx:208-281` |
//! | `IsDone`..`Intersection2` | `cxx:283-366` |
//!
//! Translation notes:
//! - `IntAna2d_AnaIntersection(l1, l2)` / `(l, circ)` become
//!   [`IntAna2dAnaIntersection::new`] followed by `perform_lin_lin` /
//!   `perform_lin_circ`, which is what the OCCT constructors forward to,
//! - the "closest point" loops of `cxx:104-112` and `cxx:265-273` keep
//!   `maxdist = RealLast()` constant inside the loop in OCCT; the port does the
//!   same, so the surviving point is the last one below `f64::MAX`,
//! - `GccEnt_BadQualifier`, `StdFail_NotDone`, `Standard_RangeError` and
//!   `Standard_OutOfRange` have no Rust equivalent at the call sites of this
//!   package and panic with the OCCT exception name.
//!
//! UNPORTED: none - all four constructors and all six accessors.

use occt_core::elib::clib2d::{parameter_circ2d, parameter_lin2d};
use occt_core::gp::circ2d::GpCirc2d;
use occt_core::gp::dir2d::GpDir2d;
use occt_core::gp::lin2d::GpLin2d;
use occt_core::gp::pnt2d::GpPnt2d;
use occt_core::gp::xy::GpXY;
use occt_core::intana2d::ana_intersection::IntAna2dAnaIntersection;

use super::gcc_ent::{GccEntPosition, GccEntQualifiedCirc};

/// `NCollection_Array1<...>(1, 2)`: the largest solution count of the four
/// constructors.
pub const MAX_SOL: usize = 2;

/// `GccAna_Lin2dTanPer` (`GccAna_Lin2dTanPer.hxx`).
#[derive(Debug, Clone, Copy)]
pub struct GccAnaLin2dTanPer {
    pub(crate) linsol: [GpLin2d; MAX_SOL],
    pub(crate) qualifier1: [GccEntPosition; MAX_SOL],
    pub(crate) pnttg1sol: [GpPnt2d; MAX_SOL],
    pub(crate) pntint2sol: [GpPnt2d; MAX_SOL],
    pub(crate) par1sol: [f64; MAX_SOL],
    pub(crate) par2sol: [f64; MAX_SOL],
    pub(crate) pararg1: [f64; MAX_SOL],
    pub(crate) pararg2: [f64; MAX_SOL],
    pub(crate) well_done: bool,
    pub(crate) nbr_sol: usize,
}

impl GccAnaLin2dTanPer {
    /// The `: linsol(1, 1), qualifier1(1, 1), ...` initializer list every
    /// constructor shares (`cxx:41-48`, `cxx:84-91`, `cxx:131-138`,
    /// `cxx:210-217`).
    pub(crate) fn empty() -> Self {
        Self {
            linsol: [GpLin2d::default(); MAX_SOL],
            qualifier1: [GccEntPosition::Unqualified; MAX_SOL],
            pnttg1sol: [GpPnt2d::zero(); MAX_SOL],
            pntint2sol: [GpPnt2d::zero(); MAX_SOL],
            par1sol: [0.0; MAX_SOL],
            par2sol: [0.0; MAX_SOL],
            pararg1: [0.0; MAX_SOL],
            pararg2: [0.0; MAX_SOL],
            well_done: false,
            nbr_sol: 0,
        }
    }

    /// `GccAna_Lin2dTanPer(ThePnt, TheLin)` (`cxx:40-79`): the line through the
    /// point orthogonal to the line.
    pub fn from_point_and_line(the_pnt: &GpPnt2d, the_lin: &GpLin2d) -> Self {
        let mut s = Self::empty();
        let dir = GpDir2d::new(-the_lin.direction().y, the_lin.direction().x)
            .expect("gp_Dir2d of the normal of a line");
        s.linsol[0] = GpLin2d::from_pnt_dir(*the_pnt, dir);
        s.pnttg1sol[0] = *the_pnt;
        if let Some(p) = last_point(&intersect_lin_lin(&s.linsol[0], the_lin)) {
            s.pntint2sol[0] = p;
        }
        s.par1sol[0] = parameter_lin2d(&s.linsol[0], &s.pnttg1sol[0]);
        s.par2sol[0] = parameter_lin2d(&s.linsol[0], &s.pntint2sol[0]);
        s.pararg1[0] = 0.0;
        s.pararg2[0] = parameter_lin2d(the_lin, &s.pntint2sol[0]);
        s.nbr_sol = 1;
        s.well_done = true;
        s
    }

    /// `GccAna_Lin2dTanPer(ThePnt, TheCircle)` (`cxx:83-120`): the line through
    /// the point and the circle centre.
    pub fn from_point_and_circle(the_pnt: &GpPnt2d, the_circle: &GpCirc2d) -> Self {
        let mut s = Self::empty();
        let c = the_circle.location();
        let dir = GpDir2d::new(c.x() - the_pnt.x(), c.y() - the_pnt.y())
            .expect("gp_VectorWithNullMagnitude in GccAna_Lin2dTanPer");
        s.linsol[0] = GpLin2d::from_pnt_dir(*the_pnt, dir);
        s.pnttg1sol[0] = *the_pnt;
        if let Some(p) = closest_point(&intersect_lin_circ(&s.linsol[0], the_circle), the_pnt) {
            s.pntint2sol[0] = p;
        }
        s.par1sol[0] = parameter_lin2d(&s.linsol[0], &s.pnttg1sol[0]);
        s.par2sol[0] = parameter_lin2d(&s.linsol[0], &s.pntint2sol[0]);
        s.pararg1[0] = 0.0;
        s.pararg2[0] = parameter_circ2d(the_circle, &s.pntint2sol[0]);
        s.nbr_sol = 1;
        s.well_done = true;
        s
    }

    /// `GccAna_Lin2dTanPer(Qualified1, TheLin)` (`cxx:129-199`): the lines
    /// tangent to the circle and orthogonal to the line.
    ///
    /// `xy = R1 * Lin.Direction().XY()` (`cxx:169`) is *not* rotated here (the
    /// direction rotation moves into the solution line), unlike
    /// [`super::lin2d_tan_par`] where the radius vector is turned by 90
    /// degrees.
    pub fn tangent_to_circle_perpendicular_to_line(
        qualified1: &GccEntQualifiedCirc,
        the_lin: &GpLin2d,
    ) -> Self {
        let mut s = Self::empty();
        let nbsol: usize;
        let mut signe: i32;
        let c1 = qualified1.qualified();

        // `cxx:147-168`.
        if qualified1.is_enclosed() {
            panic!("GccEnt_BadQualifier in GccAna_Lin2dTanPer");
        } else if qualified1.is_enclosing() {
            nbsol = 1;
            signe = -1;
        } else if qualified1.is_outside() {
            nbsol = 1;
            signe = 1;
        } else {
            nbsol = 2;
            signe = -1;
        }
        let xy = GpXY::new(
            c1.radius() * the_lin.direction().x,
            c1.radius() * the_lin.direction().y,
        );
        for _ in 1..=nbsol {
            signe = -signe;
            s.nbr_sol += 1;
            let i = s.nbr_sol - 1;
            let pnt = GpPnt2d::from_xy(
                c1.location().xy().added(&xy.multiplied_scalar(signe as f64)),
            );
            let dir = GpDir2d::new(-the_lin.direction().y, the_lin.direction().x)
                .expect("gp_Dir2d of the normal of a line");
            s.linsol[i] = GpLin2d::from_pnt_dir(pnt, dir);
            s.pnttg1sol[i] = pnt;
            if let Some(p) = last_point(&intersect_lin_lin(&s.linsol[i], the_lin)) {
                s.pntint2sol[i] = p;
            }
            s.par1sol[i] = parameter_lin2d(&s.linsol[i], &s.pnttg1sol[i]);
            s.par2sol[i] = parameter_lin2d(&s.linsol[i], &s.pntint2sol[i]);
            s.pararg1[i] = parameter_circ2d(&c1, &s.pnttg1sol[i]);
            s.pararg2[i] = parameter_lin2d(the_lin, &s.pntint2sol[i]);
            s.well_done = true;
        }
        s
    }

    /// `GccAna_Lin2dTanPer(Qualified1, TheCircle)` (`cxx:208-281`): the lines
    /// tangent to the circle and orthogonal to the second circle's radius.
    ///
    /// The `for (j = 1; j <= 2; j++)` loop at `cxx:249` always runs twice, even
    /// when the qualifier only admits one position (`cxx:230-248`); the port
    /// keeps that, as OCCT does.
    pub fn tangent_to_circle_perpendicular_to_circle(
        qualified1: &GccEntQualifiedCirc,
        the_circle: &GpCirc2d,
    ) -> Self {
        let mut s = Self::empty();
        let mut signe: i32 = 0;
        let c1 = qualified1.qualified();

        // `cxx:225-248`.
        if qualified1.is_enclosed() {
            panic!("GccEnt_BadQualifier in GccAna_Lin2dTanPer");
        } else if qualified1.is_enclosing() {
            signe = -1;
            s.qualifier1[0] = GccEntPosition::Enclosing;
        } else if qualified1.is_outside() {
            signe = 1;
            s.qualifier1[0] = GccEntPosition::Outside;
        } else if qualified1.is_unqualified() {
            signe = -1;
            s.qualifier1[0] = GccEntPosition::Enclosing;
            s.qualifier1[1] = GccEntPosition::Outside;
        }
        for _ in 1..=2 {
            s.nbr_sol += 1;
            let i = s.nbr_sol - 1;
            signe = -signe;
            let c2 = the_circle.location();
            let d1 = GpDir2d::new(c2.x() - c1.location().x(), c2.y() - c1.location().y())
                .expect("gp_VectorWithNullMagnitude in GccAna_Lin2dTanPer");
            let pnt = GpPnt2d::from_xy(
                c1.location()
                    .xy()
                    .added(&GpXY::new(d1.x, d1.y).multiplied_scalar(signe as f64 * c1.radius())),
            );
            s.linsol[i] = GpLin2d::from_pnt_dir(pnt, GpDir2d::new(-d1.y, d1.x).expect("gp_Dir2d"));
            s.pnttg1sol[i] = pnt;
            if let Some(p) =
                closest_point(&intersect_lin_circ(&s.linsol[i], the_circle), &pnt)
            {
                s.pntint2sol[i] = p;
            }
            s.par1sol[i] = parameter_lin2d(&s.linsol[i], &s.pnttg1sol[i]);
            s.par2sol[i] = parameter_lin2d(&s.linsol[i], &s.pntint2sol[i]);
            s.pararg1[i] = parameter_circ2d(&c1, &s.pnttg1sol[i]);
            s.pararg2[i] = parameter_circ2d(the_circle, &s.pntint2sol[i]);
            s.well_done = true;
        }
        s
    }

    /// `IsDone()` (`cxx:283-286`).
    pub fn is_done(&self) -> bool {
        self.well_done
    }

    /// `NbSolutions()` (`cxx:288-295`).
    pub fn nb_solutions(&self) -> usize {
        self.check_done();
        self.nbr_sol
    }

    /// `ThisSolution(Index)` (`cxx:297-308`).
    pub fn this_solution(&self, index: usize) -> GpLin2d {
        self.check_done();
        // `cxx:303-306` raises `Standard_RangeError`, not `OutOfRange`.
        assert!(
            index >= 1 && index <= self.nbr_sol,
            "Standard_RangeError in GccAna_Lin2dTanPer"
        );
        self.linsol[index - 1]
    }

    /// `WhichQualifier(Index, Qualif1)` (`cxx:310-324`).
    pub fn which_qualifier(&self, index: usize) -> GccEntPosition {
        self.check_done();
        self.check_index(index);
        self.qualifier1[index - 1]
    }

    /// `Tangency1(Index, ParSol, ParArg, Pnt)` (`cxx:326-345`).
    pub fn tangency1(&self, index: usize) -> (f64, f64, GpPnt2d) {
        self.check_done();
        self.check_index(index);
        (self.par1sol[index - 1], self.pararg1[index - 1], self.pnttg1sol[index - 1])
    }

    /// `Intersection2(Index, ParSol, ParArg, PntSol)` (`cxx:347-366`).
    pub fn intersection2(&self, index: usize) -> (f64, f64, GpPnt2d) {
        self.check_done();
        self.check_index(index);
        (self.par2sol[index - 1], self.pararg2[index - 1], self.pntint2sol[index - 1])
    }

    /// `if (!WellDone) throw StdFail_NotDone()`.
    fn check_done(&self) {
        assert!(self.well_done, "StdFail_NotDone in GccAna_Lin2dTanPer");
    }

    /// `if (Index <= 0 || Index > NbrSol) throw Standard_OutOfRange()`.
    fn check_index(&self, index: usize) {
        assert!(
            index >= 1 && index <= self.nbr_sol,
            "Standard_OutOfRange in GccAna_Lin2dTanPer"
        );
    }
}

/// `IntAna2d_AnaIntersection(l1, l2)` (`IntAna2d_AnaIntersection_1.cxx`).
fn intersect_lin_lin(l1: &GpLin2d, l2: &GpLin2d) -> IntAna2dAnaIntersection {
    let mut intp = IntAna2dAnaIntersection::new();
    intp.perform_lin_lin(l1, l2);
    intp
}

/// `IntAna2d_AnaIntersection(l, circ)` (`IntAna2d_AnaIntersection_4.cxx`).
fn intersect_lin_circ(l: &GpLin2d, c: &GpCirc2d) -> IntAna2dAnaIntersection {
    let mut intp = IntAna2dAnaIntersection::new();
    intp.perform_lin_circ(l, c);
    intp
}

/// `if (Intp.IsDone()) { if (!Intp.IsEmpty()) { for (...) pntint = Point(i); } }`
/// (`cxx:58-67`, `cxx:183-192`): keeps the last intersection point.
fn last_point(intp: &IntAna2dAnaIntersection) -> Option<GpPnt2d> {
    if !intp.is_done() || intp.is_empty() {
        return None;
    }
    let mut p = None;
    for i in 1..=intp.nb_points() {
        p = Some(*intp.point(i).value());
    }
    p
}

/// Same loop with the `Point(i).Value().Distance(...) < maxdist` filter of
/// `cxx:104-112` / `cxx:265-273`; `maxdist` is `RealLast()` and never updated.
fn closest_point(intp: &IntAna2dAnaIntersection, target: &GpPnt2d) -> Option<GpPnt2d> {
    if !intp.is_done() || intp.is_empty() {
        return None;
    }
    let maxdist = f64::MAX;
    let mut p = None;
    for i in 1..=intp.nb_points() {
        let v = *intp.point(i).value();
        if v.distance(target) < maxdist {
            p = Some(v);
        }
    }
    p
}
