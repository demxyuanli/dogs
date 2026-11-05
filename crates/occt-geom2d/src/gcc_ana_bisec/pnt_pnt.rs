//! `GccAna_Pnt2dBisec` (`GccAna_Pnt2dBisec.cxx:28-77`): the bisector of two
//! points, a single straight line perpendicular to their joining segment at
//! its midpoint.

use occt_core::gp::{GpDir2d, GpLin2d, GpPnt2d};

use super::GccAnaBisecError;

/// `class GccAna_Pnt2dBisec`.
#[derive(Debug, Clone, Copy)]
pub struct GccAnaPnt2dBisec {
    linsol: GpLin2d,
    has_sol: bool,
    well_done: bool,
}

impl GccAnaPnt2dBisec {
    /// `GccAna_Pnt2dBisec::GccAna_Pnt2dBisec` (`GccAna_Pnt2dBisec.cxx:31-56`).
    pub fn new(point1: &GpPnt2d, point2: &GpPnt2d) -> Self {
        let mut linsol = GpLin2d::from_pnt_dir(GpPnt2d::zero(), GpDir2d::default());
        let mut has_sol = false;
        // OCCT tests `Point1.Distance(Point2) > 1.0e-10` (the `gp::Resolution()`
        // form is commented out at `GccAna_Pnt2dBisec.cxx:35-36`).
        if point1.distance(point2) > 1.0e-10 {
            let dir1 = GpDir2d::new(point2.x() - point1.x(), point2.y() - point1.y())
                .unwrap_or(GpDir2d::default());
            linsol = GpLin2d::from_pnt_dir(
                GpPnt2d::new(
                    (point2.x() + point1.x()) / 2.0,
                    (point2.y() + point1.y()) / 2.0,
                ),
                GpDir2d::new(-dir1.y(), dir1.x()).unwrap_or(GpDir2d::default()),
            );
            has_sol = true;
        }
        Self {
            linsol,
            has_sol,
            well_done: true,
        }
    }

    /// `GccAna_Pnt2dBisec::IsDone` (`cxx:60-63`).
    pub fn is_done(&self) -> bool {
        self.well_done
    }

    /// `GccAna_Pnt2dBisec::HasSolution` (`cxx:65-68`).
    pub fn has_solution(&self) -> bool {
        self.has_sol
    }

    /// `GccAna_Pnt2dBisec::ThisSolution` (`cxx:70-79`).
    pub fn this_solution(&self) -> Result<GpLin2d, GccAnaBisecError> {
        if !self.well_done {
            return Err(GccAnaBisecError::NotDone);
        }
        if !self.has_sol {
            return Err(GccAnaBisecError::NoSolution);
        }
        Ok(self.linsol)
    }
}
