//! `GccAna_Lin2dTanObl` (`GccAna_Lin2dTanObl.cxx:42-311`): lines making a
//! given angle with a reference line, through a point or tangent to a circle.
//!
//! | OCCT | lines |
//! |---|---|
//! | `GccAna_Lin2dTanObl(gp_Pnt2d, gp_Lin2d, TheAngle)` | `cxx:42-88` |
//! | `GccAna_Lin2dTanObl(GccEnt_QualifiedCirc, gp_Lin2d, TheAngle)` | `cxx:99-225` |
//! | `IsDone`..`Intersection2` | `cxx:227-311` |
//!
//! Translation notes:
//! - the rotated direction is `xy = (cos(A) * Cosa - sin(A) * Sina,
//!   cos(A) * Sina + sin(A) * Cosa)` (`cxx:134-135`, `:157-158`, `:180-181`),
//!   evaluated once per branch as OCCT does,
//! - the `Unqualified` branch re-runs the intersection with `linsol(1)` (not
//!   `linsol(2)`) at `cxx:203`; the port keeps that asymmetry - `pntint2sol(2)`
//!   is the intersection of the *first* solution with `TheLine`,
//! - `GccEnt_BadQualifier`, `StdFail_NotDone` and `Standard_OutOfRange` panic
//!   with the OCCT exception name.
//!
//! UNPORTED: none - both constructors and all six accessors.

use occt_core::elib::clib2d::{parameter_circ2d, parameter_lin2d};
use occt_core::gp::dir2d::GpDir2d;
use occt_core::gp::lin2d::GpLin2d;
use occt_core::gp::pnt2d::GpPnt2d;
use occt_core::gp::xy::GpXY;
use occt_core::intana2d::ana_intersection::IntAna2dAnaIntersection;

use super::gcc_ent::{GccEntPosition, GccEntQualifiedCirc};

/// `NCollection_Array1<...>(1, 2)`: the largest solution count of the two
/// constructors.
pub const MAX_SOL: usize = 2;

/// `GccAna_Lin2dTanObl` (`GccAna_Lin2dTanObl.hxx`).
#[derive(Debug, Clone, Copy)]
pub struct GccAnaLin2dTanObl {
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

impl GccAnaLin2dTanObl {
    /// The `: linsol(1, 1), qualifier1(1, 1), ...` / `(1, 2)` initializer list
    /// both constructors share (`cxx:45-52`, `cxx:102-109`).
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

    /// `GccAna_Lin2dTanObl(ThePoint, TheLine, TheAngle)` (`cxx:42-88`).
    pub fn from_point(the_point: &GpPnt2d, the_line: &GpLin2d, the_angle: f64) -> Self {
        let mut s = Self::empty();
        let cosa = the_line.direction().x;
        let sina = the_line.direction().y;
        let dir = GpDir2d::new(
            cosa * the_angle.cos() - sina * the_angle.sin(),
            sina * the_angle.cos() + the_angle.sin() * cosa,
        )
        .expect("gp_Dir2d of the rotated direction");
        s.linsol[0] = GpLin2d::from_pnt_dir(*the_point, dir);
        s.qualifier1[0] = GccEntPosition::NoQualifier;
        s.pnttg1sol[0] = *the_point;
        let intp = intersect_lin_lin(&s.linsol[0], the_line);
        if intp.is_done() {
            if let Some(p) = last_point(&intp) {
                s.pntint2sol[0] = p;
            }
            s.par1sol[0] = parameter_lin2d(&s.linsol[0], &s.pnttg1sol[0]);
            s.par2sol[0] = parameter_lin2d(&s.linsol[0], &s.pntint2sol[0]);
            s.pararg1[0] = 0.0;
            s.pararg2[0] = parameter_lin2d(the_line, &s.pntint2sol[0]);
            s.nbr_sol = 1;
            s.well_done = true;
        }
        // `else { WellDone = false; NbrSol = 0; }` (`cxx:83-87`): `empty()`
        // already left both fields at their false/0 values.
        s
    }

    /// `GccAna_Lin2dTanObl(Qualified1, TheLine, TheAngle)` (`cxx:99-225`).
    pub fn tangent_to_circle(
        qualified1: &GccEntQualifiedCirc,
        the_line: &GpLin2d,
        the_angle: f64,
    ) -> Self {
        let mut s = Self::empty();
        // `cxx:114-119`.
        if !(qualified1.is_enclosed()
            || qualified1.is_enclosing()
            || qualified1.is_outside()
            || qualified1.is_unqualified())
        {
            panic!("GccEnt_BadQualifier in GccAna_Lin2dTanObl");
        }
        let cosa = the_line.direction().x;
        let sina = the_line.direction().y;

        if qualified1.is_enclosed() {
            panic!("GccEnt_BadQualifier in GccAna_Lin2dTanObl");
        }
        let c1 = qualified1.qualified();
        let r1 = c1.radius();
        // `gp_XY xy(cos(A) * Cosa - sin(A) * Sina, cos(A) * Sina + sin(A) * Cosa)`
        // (`cxx:134-135`).
        let xy = GpXY::new(
            the_angle.cos() * cosa - the_angle.sin() * sina,
            the_angle.cos() * sina + the_angle.sin() * cosa,
        );
        if qualified1.is_enclosing() {
            // `cxx:131-153`.
            let pnt = offset(&c1.location(), r1, GpXY::new(xy.y, -xy.x));
            s.pnttg1sol[0] = pnt;
            s.linsol[0] = GpLin2d::from_pnt_dir(pnt, GpDir2d::from_xy(&xy).expect("gp_Dir2d"));
            s.qualifier1[0] = qualified1.qualifier();
            let intp = intersect_lin_lin(&s.linsol[0], the_line);
            s.nbr_sol = 1;
            s.well_done = true;
            if intp.is_done() {
                if let Some(p) = last_point(&intp) {
                    s.pntint2sol[0] = p;
                }
            }
        } else if qualified1.is_outside() {
            // `cxx:154-176`.
            let pnt = offset(&c1.location(), r1, GpXY::new(-xy.y, xy.x));
            s.pnttg1sol[0] = pnt;
            s.linsol[0] = GpLin2d::from_pnt_dir(pnt, GpDir2d::from_xy(&xy).expect("gp_Dir2d"));
            s.qualifier1[0] = qualified1.qualifier();
            let intp = intersect_lin_lin(&s.linsol[0], the_line);
            s.well_done = true;
            s.nbr_sol = 1;
            if intp.is_done() {
                if let Some(p) = last_point(&intp) {
                    s.pntint2sol[0] = p;
                }
            }
        } else if qualified1.is_unqualified() {
            // `cxx:177-216`.
            let pnt1 = offset(&c1.location(), r1, GpXY::new(xy.y, -xy.x));
            s.pnttg1sol[0] = pnt1;
            s.linsol[0] = GpLin2d::from_pnt_dir(pnt1, GpDir2d::from_xy(&xy).expect("gp_Dir2d"));
            s.qualifier1[0] = GccEntPosition::Enclosing;
            let intp = intersect_lin_lin(&s.linsol[0], the_line);
            s.well_done = true;
            s.nbr_sol = 1;
            if intp.is_done() {
                if let Some(p) = last_point(&intp) {
                    s.pntint2sol[0] = p;
                }
            }
            let pnt2 = offset(&c1.location(), r1, GpXY::new(-xy.y, xy.x));
            s.pnttg1sol[1] = pnt2;
            s.linsol[1] = GpLin2d::from_pnt_dir(pnt2, GpDir2d::from_xy(&xy).expect("gp_Dir2d"));
            s.qualifier1[1] = GccEntPosition::Outside;
            // `Intp = IntAna2d_AnaIntersection(linsol(1), TheLine)` (`cxx:203`):
            // OCCT intersects the *first* solution again.
            let intp2 = intersect_lin_lin(&s.linsol[0], the_line);
            s.nbr_sol += 1;
            s.well_done = true;
            if intp2.is_done() {
                if let Some(p) = last_point(&intp2) {
                    s.pntint2sol[1] = p;
                }
            }
        }
        // `cxx:217-223`.
        for index in 0..s.nbr_sol {
            s.par1sol[index] = parameter_lin2d(&s.linsol[index], &s.pnttg1sol[index]);
            s.pararg1[index] = parameter_circ2d(&c1, &s.pnttg1sol[index]);
            s.par2sol[index] = parameter_lin2d(&s.linsol[index], &s.pntint2sol[index]);
            s.pararg2[index] = parameter_lin2d(the_line, &s.pntint2sol[index]);
        }
        s
    }

    /// `IsDone()` (`cxx:227-230`).
    pub fn is_done(&self) -> bool {
        self.well_done
    }

    /// `NbSolutions()` (`cxx:232-239`).
    pub fn nb_solutions(&self) -> usize {
        self.check_done();
        self.nbr_sol
    }

    /// `ThisSolution(Index)` (`cxx:241-253`).
    pub fn this_solution(&self, index: usize) -> GpLin2d {
        self.check_done();
        self.check_index(index);
        self.linsol[index - 1]
    }

    /// `WhichQualifier(Index, Qualif1)` (`cxx:255-269`).
    pub fn which_qualifier(&self, index: usize) -> GccEntPosition {
        self.check_done();
        self.check_index(index);
        self.qualifier1[index - 1]
    }

    /// `Tangency1(Index, ParSol, ParArg, PntSol)` (`cxx:271-290`).
    pub fn tangency1(&self, index: usize) -> (f64, f64, GpPnt2d) {
        self.check_done();
        self.check_index(index);
        (self.par1sol[index - 1], self.pararg1[index - 1], self.pnttg1sol[index - 1])
    }

    /// `Intersection2(Index, ParSol, ParArg, PntSol)` (`cxx:292-311`).
    pub fn intersection2(&self, index: usize) -> (f64, f64, GpPnt2d) {
        self.check_done();
        self.check_index(index);
        (self.par2sol[index - 1], self.pararg2[index - 1], self.pntint2sol[index - 1])
    }

    /// `if (!WellDone) throw StdFail_NotDone()`.
    fn check_done(&self) {
        assert!(self.well_done, "StdFail_NotDone in GccAna_Lin2dTanObl");
    }

    /// `if (Index <= 0 || Index > NbrSol) throw Standard_OutOfRange()`.
    fn check_index(&self, index: usize) {
        assert!(
            index >= 1 && index <= self.nbr_sol,
            "Standard_OutOfRange in GccAna_Lin2dTanObl"
        );
    }
}

/// `gp_Pnt2d(center.XY() + radius * dir)` (`cxx:136`, `:159`, `:182`, `:199`).
fn offset(center: &GpPnt2d, radius: f64, dir: GpXY) -> GpPnt2d {
    GpPnt2d::from_xy(center.xy().added(&dir.multiplied_scalar(radius)))
}

/// `IntAna2d_AnaIntersection(l1, l2)` (`IntAna2d_AnaIntersection_1.cxx`).
fn intersect_lin_lin(l1: &GpLin2d, l2: &GpLin2d) -> IntAna2dAnaIntersection {
    let mut intp = IntAna2dAnaIntersection::new();
    intp.perform_lin_lin(l1, l2);
    intp
}

/// `if (!Intp.IsEmpty()) { for (...) pntint2sol = Point(i); }` (`cxx:145-152`,
/// `:168-175`, `:189-198`, `:206-215`): keeps the last intersection point.
fn last_point(intp: &IntAna2dAnaIntersection) -> Option<GpPnt2d> {
    if intp.is_empty() {
        return None;
    }
    let mut p = None;
    for i in 1..=intp.nb_points() {
        p = Some(*intp.point(i).value());
    }
    p
}
