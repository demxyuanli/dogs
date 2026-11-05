//! `GccAna_Circ2dBisec` (`GccAna_Circ2dBisec.cxx:36-414`): the bisector of two
//! circles, a line, an ellipse, a circle or a hyperbola branch depending on the
//! relative position of the two circles.
//!
//! The constructor only classifies the configuration (`intersection`,
//! `sameradius`, `NbrSol`) and `ThisSolution` builds the requested arc, exactly
//! as OCCT does - the constructor never allocates a solution.

use occt_core::gp::{GpAx2d, GpAx22d, GpCirc2d, GpDir2d, GpElips2d, GpHypr2d, GpLin2d, GpPnt2d};
use occt_core::precision::CONFUSION;

use super::{check_index, GccAnaBisecError};
use crate::gcc_int::GccIntBisec;

/// `class GccAna_Circ2dBisec`.
#[derive(Debug, Clone, Copy)]
pub struct GccAnaCirc2dBisec {
    circle1: GpCirc2d,
    circle2: GpCirc2d,
    /// Relative position of the two circles, `0..=4` (`cxx:65-129`).
    intersection: i32,
    sameradius: bool,
    nbr_sol: i32,
    well_done: bool,
}

impl GccAnaCirc2dBisec {
    /// `GccAna_Circ2dBisec::GccAna_Circ2dBisec` (`cxx:36-129`).
    pub fn new(circ1: &GpCirc2d, circ2: &GpCirc2d) -> Self {
        let tol = CONFUSION;
        let mut r1 = circ1.radius();
        let mut r2 = circ2.radius();
        let sameradius = (r1 - r2).abs() <= tol;
        let (circle1, circle2) = if r1 < r2 {
            r1 = circ2.radius();
            r2 = circ1.radius();
            (*circ2, *circ1)
        } else {
            (*circ1, *circ2)
        };
        let dist = circle2.location().distance(&circle1.location());
        let (intersection, nbr_sol) = if r1 - dist - r2 > tol {
            (0, 2)
        } else if (r1 - dist - r2).abs() <= tol {
            (1, if sameradius { 0 } else { 2 })
        } else if (dist + r2 - r1 > tol) && (r1 - dist + r2 > tol) {
            (2, if sameradius { 2 } else { 3 })
        } else if (r1 - dist + r2).abs() <= tol {
            (3, if sameradius { 2 } else { 3 })
        } else {
            (4, if sameradius { 3 } else { 4 })
        };
        Self {
            circle1,
            circle2,
            intersection,
            sameradius,
            nbr_sol,
            well_done: true,
        }
    }

    /// `GccAna_Circ2dBisec::IsDone` (`cxx:400-403`).
    pub fn is_done(&self) -> bool {
        self.well_done
    }

    /// `GccAna_Circ2dBisec::NbSolutions` (`cxx:405-414`).
    pub fn nb_solutions(&self) -> Result<i32, GccAnaBisecError> {
        if !self.well_done {
            return Err(GccAnaBisecError::NotDone);
        }
        Ok(self.nbr_sol)
    }

    /// `GccAna_Circ2dBisec::ThisSolution` (`cxx:139-397`).
    pub fn this_solution(&self, index: i32) -> Result<GccIntBisec, GccAnaBisecError> {
        check_index(self.well_done, self.nbr_sol, index)?;

        let tol = 1.0e-14;
        let xcencir1 = self.circle1.location().x();
        let ycencir1 = self.circle1.location().y();
        let xcencir2 = self.circle2.location().x();
        let ycencir2 = self.circle2.location().y();
        let dist = self.circle1.location().distance(&self.circle2.location());

        let pcen = GpPnt2d::new((xcencir1 + xcencir2) / 2.0, (ycencir1 + ycencir2) / 2.0);
        // `gp_Dir2d dircen, medcen;` (`cxx:157`) - both stay at the default
        // `(1, 0)` when the two centres coincide (`dist <= Tol`).
        let mut dircen = GpDir2d::default();
        let mut medcen = GpDir2d::default();
        if dist > tol {
            dircen = GpDir2d::new(xcencir2 - xcencir1, ycencir2 - ycencir1)
                .unwrap_or(GpDir2d::default());
            medcen = GpDir2d::new(ycencir2 - ycencir1, xcencir1 - xcencir2)
                .unwrap_or(GpDir2d::default());
        }
        let acenx = GpAx2d::new(pcen, GpDir2d::default());
        let acencen = GpAx2d::new(pcen, dircen);

        let r1 = self.circle1.radius();
        let r2 = self.circle2.radius();
        let same_center = (xcencir2 - xcencir1).abs() < tol && (ycencir2 - ycencir1).abs() < tol;

        let sol = if self.nbr_sol == 1 && self.intersection == 0 {
            // `cxx:170-183`: unreachable through the constructor, kept because
            // OCCT keeps it.
            let r = if index == 1 {
                (r1 + r2) / 2.0
            } else {
                (r1 - r2) / 2.0
            };
            GccIntBisec::Circle(GpCirc2d::new(GpAx22d::from_xdir(acenx.loc, acenx.vdir), r))
        } else if self.nbr_sol == 2 && self.intersection == 1 {
            if index == 1 {
                let e = GpElips2d::new(
                    GpAx22d::from_xdir(acencen.loc, acencen.vdir),
                    (r1 + r2) / 2.0,
                    ((r1 * r1 + r2 * r2 - dist * dist) / 4.0 + r1 * r2 / 2.0).sqrt(),
                );
                GccIntBisec::Ellipse(e)
            } else {
                GccIntBisec::Line(GpLin2d::from_pnt_dir(self.circle1.location(), dircen))
            }
        } else if self.nbr_sol == 2 && self.intersection == 0 {
            if index == 1 {
                if same_center {
                    GccIntBisec::Circle(GpCirc2d::new(
                        GpAx22d::from_xdir(acenx.loc, acenx.vdir),
                        (r1 + r2) / 2.0,
                    ))
                } else {
                    GccIntBisec::Ellipse(GpElips2d::new(
                        GpAx22d::from_xdir(acencen.loc, acencen.vdir),
                        (r1 + r2) / 2.0,
                        ((r1 * r1 + r2 * r2 - dist * dist) / 4.0 + r1 * r2 / 2.0).sqrt(),
                    ))
                }
            } else if same_center {
                GccIntBisec::Circle(GpCirc2d::new(
                    GpAx22d::from_xdir(acencen.loc, acencen.vdir),
                    (r1 - r2) / 2.0,
                ))
            } else {
                GccIntBisec::Ellipse(GpElips2d::new(
                    GpAx22d::from_xdir(acencen.loc, acencen.vdir),
                    (r1 - r2) / 2.0,
                    ((r1 * r1 + r2 * r2 - dist * dist) / 4.0 - r1 * r2 / 2.0).sqrt(),
                ))
            }
        } else if self.intersection == 2 {
            if self.sameradius {
                if index == 1 {
                    GccIntBisec::Line(GpLin2d::from_pnt_dir(pcen, medcen))
                } else {
                    GccIntBisec::Ellipse(GpElips2d::new(
                        GpAx22d::from_xdir(acencen.loc, acencen.vdir),
                        r1,
                        (r1 * r1 - dist * dist / 4.0).sqrt(),
                    ))
                }
            } else if index == 1 {
                let h = GpHypr2d::new(
                    GpAx22d::from_xdir(acencen.loc, acencen.vdir),
                    (r1 - r2) / 2.0,
                    (dist * dist - (r1 - r2) * (r1 - r2)).sqrt() / 2.0,
                );
                GccIntBisec::Hyperbola(h)
            } else if index == 2 {
                let h = GpHypr2d::new(
                    GpAx22d::from_xdir(acencen.loc, acencen.vdir),
                    (r1 - r2) / 2.0,
                    (dist * dist - (r1 - r2) * (r1 - r2)).sqrt() / 2.0,
                );
                GccIntBisec::Hyperbola(h.other_branch())
            } else {
                GccIntBisec::Ellipse(GpElips2d::new(
                    GpAx22d::from_xdir(acencen.loc, acencen.vdir),
                    (r1 + r2) / 2.0,
                    ((r1 * r1 + r2 * r2 - dist * dist) / 4.0 + r1 * r2 / 2.0).sqrt(),
                ))
            }
        } else if self.intersection == 3 {
            if self.sameradius {
                if index == 1 {
                    GccIntBisec::Line(GpLin2d::from_pnt_dir(pcen, dircen))
                } else if index == 2 {
                    GccIntBisec::Line(GpLin2d::from_pnt_dir(pcen, medcen))
                } else {
                    return Err(GccAnaBisecError::OutOfRange);
                }
            } else if index == 1 {
                GccIntBisec::Line(GpLin2d::from_pnt_dir(pcen, dircen))
            } else if index == 2 {
                let h = GpHypr2d::new(
                    GpAx22d::from_xdir(acencen.loc, acencen.vdir),
                    (r1 - r2) / 2.0,
                    (dist * dist - (r1 - r2) * (r1 - r2)).sqrt() / 2.0,
                );
                GccIntBisec::Hyperbola(h)
            } else {
                let h = GpHypr2d::new(
                    GpAx22d::from_xdir(acencen.loc, acencen.vdir),
                    (r1 - r2) / 2.0,
                    (dist * dist - (r1 - r2) * (r1 - r2)).sqrt() / 2.0,
                );
                GccIntBisec::Hyperbola(h.other_branch())
            }
        } else {
            // `intersection == 4` (`cxx:342-395`).
            if self.sameradius {
                if index == 1 {
                    GccIntBisec::Line(GpLin2d::from_pnt_dir(pcen, medcen))
                } else if index == 2 {
                    let h = GpHypr2d::new(
                        GpAx22d::from_xdir(acencen.loc, acencen.vdir),
                        r1,
                        (dist * dist - 4.0 * r1 * r1).sqrt() / 2.0,
                    );
                    GccIntBisec::Hyperbola(h)
                } else {
                    let h = GpHypr2d::new(
                        GpAx22d::from_xdir(acencen.loc, acencen.vdir),
                        r1,
                        (dist * dist - 4.0 * r1 * r1).sqrt() / 2.0,
                    );
                    GccIntBisec::Hyperbola(h.other_branch())
                }
            } else if index == 1 || index == 2 {
                let h = GpHypr2d::new(
                    GpAx22d::from_xdir(acencen.loc, acencen.vdir),
                    (r1 - r2) / 2.0,
                    (dist * dist - (r1 - r2) * (r1 - r2)).sqrt() / 2.0,
                );
                if index == 1 {
                    GccIntBisec::Hyperbola(h)
                } else {
                    GccIntBisec::Hyperbola(h.other_branch())
                }
            } else {
                let h = GpHypr2d::new(
                    GpAx22d::from_xdir(acencen.loc, acencen.vdir),
                    (r1 + r2) / 2.0,
                    (dist * dist - (r1 + r2) * (r1 + r2)).sqrt() / 2.0,
                );
                if index == 3 {
                    GccIntBisec::Hyperbola(h)
                } else {
                    GccIntBisec::Hyperbola(h.other_branch())
                }
            }
        };
        Ok(sol)
    }
}
