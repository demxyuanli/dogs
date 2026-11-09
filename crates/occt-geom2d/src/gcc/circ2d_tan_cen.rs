//! `GccAna_Circ2dTanCen` (`GccAna_Circ2dTanCen.cxx:42-356`): circles with a
//! given centre, tangent to a circle, to a line or through a point.
//!
//! | OCCT | lines |
//! |---|---|
//! | `GccAna_Circ2dTanCen(GccEnt_QualifiedCirc, gp_Pnt2d, Tolerance)` | `cxx:42-190` |
//! | `GccAna_Circ2dTanCen(gp_Lin2d, gp_Pnt2d)` | `cxx:198-238` |
//! | `GccAna_Circ2dTanCen(gp_Pnt2d, gp_Pnt2d)` | `cxx:246-272` |
//! | `IsDone`..`IsTheSame1` | `cxx:276-356` |
//!
//! Translation notes:
//! - the concentric case (`dist <= gp::Resolution()`, `cxx:118-121` and
//!   `cxx:180-188`) leaves `pnttg1sol`/`par1sol`/`pararg1` at their defaults and
//!   sets `TheSame1 = 1`; [`tangency1`](GccAnaCirc2dTanCen::tangency1) answers
//!   `None` there, where OCCT raises `StdFail_NotDone` (`cxx:334-335`),
//! - `gp::Resolution()` is `RealSmall()` (`DBL_MIN`), mirrored by
//!   [`occt_core::precision::REAL_SMALL`],
//! - `GccEnt_BadQualifier`, `StdFail_NotDone` and `Standard_OutOfRange` panic
//!   with the OCCT exception name.
//!
//! UNPORTED: none - all three constructors and all six accessors.

use occt_core::elib::clib2d::{parameter_circ2d, parameter_lin2d};
use occt_core::gp::circ2d::GpCirc2d;
use occt_core::gp::dir2d::GpDir2d;
use occt_core::gp::lin2d::GpLin2d;
use occt_core::gp::pnt2d::GpPnt2d;
use occt_core::precision::REAL_SMALL;

use super::circ2d_2tan_rad::{axis_at, offset_by};
use super::gcc_ent::{GccEntPosition, GccEntQualifiedCirc};

/// `NCollection_Array1<...>(1, 2)`: the largest solution count of the three
/// constructors.
pub const MAX_SOL: usize = 2;

/// `GccAna_Circ2dTanCen` (`GccAna_Circ2dTanCen.hxx`).
#[derive(Debug, Clone, Copy)]
pub struct GccAnaCirc2dTanCen {
    pub(crate) cirsol: [GpCirc2d; MAX_SOL],
    pub(crate) qualifier1: [GccEntPosition; MAX_SOL],
    pub(crate) the_same1: [i32; MAX_SOL],
    pub(crate) pnttg1sol: [GpPnt2d; MAX_SOL],
    pub(crate) par1sol: [f64; MAX_SOL],
    pub(crate) pararg1: [f64; MAX_SOL],
    pub(crate) well_done: bool,
    pub(crate) nbr_sol: usize,
}

impl GccAnaCirc2dTanCen {
    /// The `: cirsol(1, 2), qualifier1(1, 2), ...` initializer list every
    /// constructor shares (`cxx:51-56`, `cxx:205-210`, `cxx:253-258`).
    pub(crate) fn empty() -> Self {
        Self {
            cirsol: [GpCirc2d::default(); MAX_SOL],
            qualifier1: [GccEntPosition::Unqualified; MAX_SOL],
            the_same1: [0; MAX_SOL],
            pnttg1sol: [GpPnt2d::zero(); MAX_SOL],
            par1sol: [0.0; MAX_SOL],
            pararg1: [0.0; MAX_SOL],
            well_done: false,
            nbr_sol: 0,
        }
    }

    /// `GccAna_Circ2dTanCen(Qualified1, Pcenter, Tolerance)` (`cxx:42-190`).
    pub fn tangent_to_circle(
        qualified1: &GccEntQualifiedCirc,
        pcenter: &GpPnt2d,
        tolerance: f64,
    ) -> Self {
        let mut s = Self::empty();
        // `cxx:62-67`.
        if !(qualified1.is_enclosed()
            || qualified1.is_enclosing()
            || qualified1.is_outside()
            || qualified1.is_unqualified())
        {
            panic!("GccEnt_BadQualifier in GccAna_Circ2dTanCen");
        }
        let tol = tolerance.abs();
        let c1 = qualified1.qualified();
        let r1 = c1.radius();
        let center1 = c1.location();
        let mut radius = 0.0f64;
        let mut signe: i32 = 0;
        let dist = pcenter.distance(&center1);

        if !qualified1.is_unqualified() {
            // `cxx:79-131`.
            if qualified1.is_enclosed() {
                if dist - r1 <= tol {
                    radius = (r1 - dist).abs();
                    signe = 1;
                } else {
                    s.well_done = true;
                }
            } else if qualified1.is_enclosing() {
                radius = r1 + dist;
                signe = -1;
            } else if qualified1.is_outside() {
                if dist < r1 - tol {
                    s.well_done = true;
                } else {
                    radius = (r1 - dist).abs();
                    signe = -1;
                }
            }
            if signe != 0 {
                s.nbr_sol += 1;
                let i = s.nbr_sol - 1;
                s.cirsol[i] = GpCirc2d::new(axis_at(*pcenter), radius);
                s.qualifier1[i] = qualified1.qualifier();
                if dist <= REAL_SMALL {
                    s.the_same1[i] = 1;
                } else {
                    s.the_same1[i] = 0;
                    let d = GpDir2d::new(pcenter.x() - center1.x(), pcenter.y() - center1.y())
                        .expect("gp_VectorWithNullMagnitude in GccAna_Circ2dTanCen");
                    s.pnttg1sol[i] = offset_by(pcenter, signe as f64 * radius, &d);
                    s.par1sol[i] = parameter_circ2d(&s.cirsol[i], &s.pnttg1sol[i]);
                    s.pararg1[i] = parameter_circ2d(&c1, &s.pnttg1sol[i]);
                }
                s.well_done = true;
            }
        } else if dist >= REAL_SMALL {
            // `cxx:136-179`.
            signe = 1;
            for _ in 1..=2 {
                signe = -signe;
                let signe1: i32 = if r1 - dist <= 0.0 { -1 } else { -signe };
                radius = (r1 + signe as f64 * dist).abs();
                s.nbr_sol += 1;
                let i = s.nbr_sol - 1;
                s.cirsol[i] = GpCirc2d::new(axis_at(*pcenter), radius);
                let distcc1 = pcenter.distance(&center1);
                if !qualified1.is_unqualified() {
                    s.qualifier1[i] = qualified1.qualifier();
                } else if (distcc1 + radius - r1).abs() < tol {
                    s.qualifier1[i] = GccEntPosition::Enclosed;
                } else if (distcc1 - r1 - radius).abs() < tol {
                    s.qualifier1[i] = GccEntPosition::Outside;
                } else {
                    s.qualifier1[i] = GccEntPosition::Enclosing;
                }
                s.the_same1[i] = 0;
                s.well_done = true;
                let d = GpDir2d::new(pcenter.x() - center1.x(), pcenter.y() - center1.y())
                    .expect("gp_VectorWithNullMagnitude in GccAna_Circ2dTanCen");
                s.pnttg1sol[i] = offset_by(pcenter, signe1 as f64 * radius, &d);
                s.par1sol[i] = parameter_circ2d(&s.cirsol[i], &s.pnttg1sol[i]);
                s.pararg1[i] = parameter_circ2d(&c1, &s.pnttg1sol[i]);
            }
        } else {
            // `cxx:180-188`: the centre is the circle centre.
            s.nbr_sol += 1;
            s.cirsol[0] = c1;
            s.qualifier1[0] = qualified1.qualifier();
            s.the_same1[0] = 1;
            s.well_done = true;
        }
        s
    }

    /// `GccAna_Circ2dTanCen(Linetan, Pcenter)` (`cxx:198-238`).
    pub fn tangent_to_line(linetan: &GpLin2d, pcenter: &GpPnt2d) -> Self {
        let mut s = Self::empty();
        let rayon = linetan.distance(pcenter);
        s.cirsol[0] = GpCirc2d::new(axis_at(*pcenter), rayon);
        s.qualifier1[0] = GccEntPosition::NoQualifier;
        s.the_same1[0] = 0;
        let xloc = linetan.location().x();
        let yloc = linetan.location().y();
        let xdir = linetan.direction().x;
        let ydir = linetan.direction().y;
        let to_origin = GpDir2d::new(xloc - pcenter.x(), yloc - pcenter.y())
            .expect("gp_VectorWithNullMagnitude in GccAna_Circ2dTanCen");
        let normal = GpDir2d::new(-ydir, xdir).expect("gp_Dir2d of the line normal");
        if to_origin.dot(&normal) > 0.0 {
            s.pnttg1sol[0] = GpPnt2d::new(
                pcenter.x() + rayon * -ydir,
                pcenter.y() + rayon * xdir,
            );
        } else {
            s.pnttg1sol[0] = GpPnt2d::new(pcenter.x() + rayon * ydir, pcenter.y() + rayon * -xdir);
        }
        s.par1sol[0] = parameter_circ2d(&s.cirsol[0], &s.pnttg1sol[0]);
        s.pararg1[0] = parameter_lin2d(linetan, &s.pnttg1sol[0]);
        s.nbr_sol = 1;
        s.well_done = true;
        s
    }

    /// `GccAna_Circ2dTanCen(Point1, Pcenter)` (`cxx:246-272`).
    pub fn through_point(point1: &GpPnt2d, pcenter: &GpPnt2d) -> Self {
        let mut s = Self::empty();
        let rayon = point1.distance(pcenter);
        s.cirsol[0] = GpCirc2d::new(axis_at(*pcenter), rayon);
        s.qualifier1[0] = GccEntPosition::NoQualifier;
        s.the_same1[0] = 0;
        s.pnttg1sol[0] = *point1;
        s.par1sol[0] = parameter_circ2d(&s.cirsol[0], &s.pnttg1sol[0]);
        s.pararg1[0] = 0.0;
        s.nbr_sol = 1;
        s.well_done = true;
        s
    }

    /// `IsDone()` (`cxx:276-279`).
    pub fn is_done(&self) -> bool {
        self.well_done
    }

    /// `NbSolutions()` (`cxx:281-284`): no `WellDone` test in OCCT.
    pub fn nb_solutions(&self) -> usize {
        self.nbr_sol
    }

    /// `ThisSolution(Index)` (`cxx:286-293`): only the index is checked.
    pub fn this_solution(&self, index: usize) -> GpCirc2d {
        self.check_index(index);
        self.cirsol[index - 1]
    }

    /// `WhichQualifier(Index, Qualif1)` (`cxx:295-309`).
    pub fn which_qualifier(&self, index: usize) -> GccEntPosition {
        self.check_done();
        self.check_index(index);
        self.qualifier1[index - 1]
    }

    /// `Tangency1(Index, ParSol, ParArg, PntSol)` (`cxx:311-337`).
    ///
    /// Returns `None` where `TheSame1(Index) != 0` makes OCCT raise
    /// `StdFail_NotDone` (`cxx:332-335`).
    pub fn tangency1(&self, index: usize) -> Option<(f64, f64, GpPnt2d)> {
        self.check_done();
        self.check_index(index);
        if self.the_same1[index - 1] != 0 {
            return None;
        }
        Some((self.par1sol[index - 1], self.pararg1[index - 1], self.pnttg1sol[index - 1]))
    }

    /// `IsTheSame1(Index)` (`cxx:339-356`).
    pub fn is_the_same1(&self, index: usize) -> bool {
        self.check_done();
        self.check_index(index);
        self.the_same1[index - 1] != 0
    }

    /// `if (!WellDone) throw StdFail_NotDone()`.
    fn check_done(&self) {
        assert!(self.well_done, "StdFail_NotDone in GccAna_Circ2dTanCen");
    }

    /// `if (Index <= 0 || Index > NbrSol) throw Standard_OutOfRange()`.
    fn check_index(&self, index: usize) {
        assert!(
            index >= 1 && index <= self.nbr_sol,
            "Standard_OutOfRange in GccAna_Circ2dTanCen"
        );
    }
}
