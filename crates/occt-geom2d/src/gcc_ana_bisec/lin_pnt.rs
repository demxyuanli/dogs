//! `GccAna_LinPnt2dBisec` (`GccAna_LinPnt2dBisec.cxx:32-89`): the bisector of
//! a straight line and a point, a parabola (or a line when the point lies on
//! the line).

use occt_core::gp::{GpAx2d, GpAx22d, GpDir2d, GpLin2d, GpParab2d, GpPnt2d};

use super::GccAnaBisecError;
use crate::gcc_int::GccIntBisec;

/// `class GccAna_LinPnt2dBisec`.
#[derive(Debug, Clone, Copy)]
pub struct GccAnaLinPnt2dBisec {
    bissol: Option<GccIntBisec>,
    well_done: bool,
}

impl GccAnaLinPnt2dBisec {
    /// `GccAna_LinPnt2dBisec::GccAna_LinPnt2dBisec` (`cxx:36-72`).
    pub fn new(line1: &GpLin2d, point2: &GpPnt2d) -> Self {
        let xdir = line1.direction().x();
        let ydir = line1.direction().y();
        let xloc = line1.location().x();
        let yloc = line1.location().y();
        let dist = line1.distance(point2);
        let bissol;
        // OCCT: `if (dist > 1.0e-10)` (`cxx:44`).
        if dist > 1.0e-10 {
            let xpoint2 = point2.x();
            let ypoint2 = point2.y();
            if (-ydir * (xpoint2 - xloc) + xdir * (ypoint2 - yloc)) > 0.0 {
                let loc = GpPnt2d::new(
                    point2.x() - dist / 2.0 * -ydir,
                    point2.y() - dist / 2.0 * xdir,
                );
                let axe = GpAx2d::new(
                    loc,
                    GpDir2d::new(-ydir, xdir).unwrap_or(GpDir2d::default()),
                );
                bissol = GccIntBisec::Parabola(GpParab2d::new(
                    GpAx22d::from_xdir(axe.loc, axe.vdir),
                    dist / 2.0,
                ));
            } else {
                let loc = GpPnt2d::new(
                    point2.x() + dist / 2.0 * -ydir,
                    point2.y() + dist / 2.0 * xdir,
                );
                let axe = GpAx2d::new(
                    loc,
                    GpDir2d::new(ydir, -xdir).unwrap_or(GpDir2d::default()),
                );
                bissol = GccIntBisec::Parabola(GpParab2d::new(
                    GpAx22d::from_xdir(axe.loc, axe.vdir),
                    dist / 2.0,
                ));
            }
        } else {
            let line = GpLin2d::from_pnt_dir(
                *point2,
                GpDir2d::new(-ydir, xdir).unwrap_or(GpDir2d::default()),
            );
            bissol = GccIntBisec::Line(line);
        }
        Self {
            bissol: Some(bissol),
            well_done: true,
        }
    }

    /// `GccAna_LinPnt2dBisec::IsDone` (`cxx:76-79`).
    pub fn is_done(&self) -> bool {
        self.well_done
    }

    /// `GccAna_LinPnt2dBisec::ThisSolution` (`cxx:81-89`).
    pub fn this_solution(&self) -> Result<GccIntBisec, GccAnaBisecError> {
        if !self.well_done {
            return Err(GccAnaBisecError::NotDone);
        }
        self.bissol.ok_or(GccAnaBisecError::NotDone)
    }
}
