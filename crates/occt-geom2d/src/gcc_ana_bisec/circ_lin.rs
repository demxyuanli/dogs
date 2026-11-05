//! `GccAna_CircLin2dBisec` (`GccAna_CircLin2dBisec.cxx:35-160`): the bisector
//! of a circle and a straight line, a pair of parabolas (degenerating to a
//! straight line when the line is tangent to the circle).

use occt_core::gp::{GpAx2d, GpAx22d, GpCirc2d, GpDir2d, GpLin2d, GpParab2d, GpPnt2d};
use occt_core::precision::REAL_SMALL;

use super::{check_index, GccAnaBisecError};
use crate::gcc_int::GccIntBisec;

/// `class GccAna_CircLin2dBisec`.
#[derive(Debug, Clone, Copy)]
pub struct GccAnaCircLin2dBisec {
    circle: GpCirc2d,
    line: GpLin2d,
    nbr_sol: i32,
    well_done: bool,
}

impl GccAnaCircLin2dBisec {
    /// `GccAna_CircLin2dBisec::GccAna_CircLin2dBisec` (`cxx:35-56`): stores the
    /// arguments and fixes `NbrSol = 2`, `WellDone = true`.
    pub fn new(circle: &GpCirc2d, line: &GpLin2d) -> Self {
        Self {
            circle: *circle,
            line: *line,
            nbr_sol: 2,
            well_done: true,
        }
    }

    /// `GccAna_CircLin2dBisec::IsDone` (`cxx:134-137`).
    pub fn is_done(&self) -> bool {
        self.well_done
    }

    /// `GccAna_CircLin2dBisec::NbSolutions` (`cxx:141-144`).
    pub fn nb_solutions(&self) -> i32 {
        self.nbr_sol
    }

    /// `GccAna_CircLin2dBisec::ThisSolution` (`cxx:66-131`).
    ///
    /// The two parabolas are built lazily, exactly as OCCT does: the
    /// tangency branch (`cxx:81-86`) answers a straight line through the circle
    /// centre along the line direction rotated by +90 degrees.
    pub fn this_solution(&self, index: i32) -> Result<GccIntBisec, GccAnaBisecError> {
        check_index(self.well_done, self.nbr_sol, index)?;

        let xdir = self.line.direction().x();
        let ydir = self.line.direction().y();
        let xloc = self.line.location().x();
        let yloc = self.line.location().y();
        let xcencir = self.circle.location().x();
        let ycencir = self.circle.location().y();
        let r1 = self.circle.radius();
        let dist = self.line.distance(&self.circle.location());

        // Tangency test of `cxx:81`: `abs(dist - R1) <= gp::Resolution()`.
        if (dist - r1).abs() <= REAL_SMALL && index == 1 {
            let line = GpLin2d::from_pnt_dir(
                self.circle.location(),
                GpDir2d::new(-ydir, xdir).unwrap_or(GpDir2d::default()),
            );
            return Ok(GccIntBisec::Line(line));
        }

        let signe = if (-ydir * (xcencir - xloc) + xdir * (ycencir - yloc)) > 0.0 {
            1.0
        } else {
            -1.0
        };
        if dist != r1 {
            if index == 1 {
                let axe = GpAx2d::new(
                    GpPnt2d::new(
                        xcencir + signe * ydir * (dist + r1) / 2.0,
                        ycencir - signe * xdir * (dist + r1) / 2.0,
                    ),
                    GpDir2d::new(-signe * ydir, signe * xdir).unwrap_or(GpDir2d::default()),
                );
                Ok(GccIntBisec::Parabola(GpParab2d::new(
                    GpAx22d::from_xdir(axe.loc, axe.vdir),
                    (dist + r1) / 2.0,
                )))
            } else {
                let dir = if dist < r1 {
                    GpDir2d::new(signe * ydir, -signe * xdir)
                } else {
                    GpDir2d::new(-signe * ydir, signe * xdir)
                };
                let axe = GpAx2d::new(
                    GpPnt2d::new(
                        xcencir + signe * ydir * (dist - r1) / 2.0,
                        ycencir - signe * xdir * (dist - r1) / 2.0,
                    ),
                    dir.unwrap_or(GpDir2d::default()),
                );
                Ok(GccIntBisec::Parabola(GpParab2d::new(
                    GpAx22d::from_xdir(axe.loc, axe.vdir),
                    (dist - r1).abs() / 2.0,
                )))
            }
        } else {
            let axe = GpAx2d::new(
                GpPnt2d::new(
                    xcencir + signe * ydir * (dist + r1) / 2.0,
                    ycencir - signe * xdir * (dist + r1) / 2.0,
                ),
                GpDir2d::new(signe * (-ydir), signe * xdir).unwrap_or(GpDir2d::default()),
            );
            Ok(GccIntBisec::Parabola(GpParab2d::new(
                GpAx22d::from_xdir(axe.loc, axe.vdir),
                r1,
            )))
        }
    }
}
