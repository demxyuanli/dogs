//! `GccAna_Lin2dBisec` (`GccAna_Lin2dBisec.cxx:39-194`): the two bisecting
//! lines of two straight lines.
//!
//! | OCCT | lines |
//! |---|---|
//! | `GccAna_Lin2dBisec(gp_Lin2d, gp_Lin2d)` | `cxx:39-124` |
//! | `IsDone`..`Intersection2` | `cxx:126-194` |
//!
//! Translation notes:
//! - the parallel branch keeps OCCT's distance formula
//!   `Lin1.Distance(Lin2.Location()) / 2.0` with the `cross` sign flip
//!   (`cxx:68-73`), *not* `Lin1.Distance(Lin2)`,
//! - `linsol(NbrSol).Reverse()` (`cxx:102-105`) only flips the direction vector
//!   of the line, which is [`GpLin2d::pos`]` .vdir`,
//! - `StdFail_NotDone` / `Standard_OutOfRange` panic with the OCCT exception
//!   name, as the other `GccAna` ports of this crate do.

use occt_core::elib::clib2d::parameter_lin2d;
use occt_core::gp::{GpDir2d, GpLin2d, GpPnt2d};
use occt_core::intana2d::ana_intersection::IntAna2dAnaIntersection;

use super::GccAnaBisecError;

/// `NCollection_Array1<...>(1, 2)` of every field.
pub const MAX_SOL: usize = 2;

/// `class GccAna_Lin2dBisec`.
#[derive(Debug, Clone, Copy)]
pub struct GccAnaLin2dBisec {
    linsol: [GpLin2d; MAX_SOL],
    pntint1sol: [GpPnt2d; MAX_SOL],
    pntint2sol: [GpPnt2d; MAX_SOL],
    par1sol: [f64; MAX_SOL],
    par2sol: [f64; MAX_SOL],
    pararg1: [f64; MAX_SOL],
    pararg2: [f64; MAX_SOL],
    well_done: bool,
    nbr_sol: usize,
}

impl GccAnaLin2dBisec {
    /// `GccAna_Lin2dBisec::GccAna_Lin2dBisec` (`cxx:39-124`).
    pub fn new(lin1: &GpLin2d, lin2: &GpLin2d) -> Self {
        let mut me = Self {
            linsol: [GpLin2d::default(); MAX_SOL],
            pntint1sol: [GpPnt2d::zero(); MAX_SOL],
            pntint2sol: [GpPnt2d::zero(); MAX_SOL],
            par1sol: [0.0; MAX_SOL],
            par2sol: [0.0; MAX_SOL],
            pararg1: [0.0; MAX_SOL],
            pararg2: [0.0; MAX_SOL],
            well_done: false,
            nbr_sol: 0,
        };

        let mut intp = IntAna2dAnaIntersection::new();
        intp.perform_lin_lin(lin1, lin2);
        if intp.is_done() {
            if intp.parallel_elements() {
                if intp.identical_elements() {
                    me.nbr_sol = 1;
                    me.well_done = true;
                    me.linsol[0] = *lin1;
                } else {
                    // Attention: not `Lin1.Distance(Lin2)`, see `cxx:66-69`.
                    let mut dist = lin1.distance(&lin2.location()) / 2.0;
                    let cross = (-lin2.direction().y) * (lin1.location().x() - lin2.location().x())
                        + lin2.direction().x * (lin1.location().y() - lin2.location().y());
                    if cross < 0.0 {
                        dist = -dist;
                    }
                    me.nbr_sol += 1;
                    me.well_done = true;
                    me.linsol[me.nbr_sol - 1] = GpLin2d::from_pnt_dir(
                        GpPnt2d::new(
                            lin2.location().x() + (-lin2.direction().y * dist),
                            lin2.location().y() + (lin2.direction().x * dist),
                        ),
                        *lin2.direction(),
                    );
                }
            } else if !intp.is_empty() {
                for i in 1..=intp.nb_points() {
                    let p = *intp.point(i).value();
                    me.nbr_sol += 1;
                    me.linsol[me.nbr_sol - 1] = GpLin2d::from_pnt_dir(
                        p,
                        GpDir2d::new(
                            lin1.direction().x + lin2.direction().x,
                            lin1.direction().y + lin2.direction().y,
                        )
                        .expect("gp_Dir2d of two non-opposite directions"),
                    );
                    me.nbr_sol += 1;
                    me.linsol[me.nbr_sol - 1] = GpLin2d::from_pnt_dir(
                        p,
                        GpDir2d::new(
                            lin1.direction().x - lin2.direction().x,
                            lin1.direction().y - lin2.direction().y,
                        )
                        .expect("gp_Dir2d of two non-identical directions"),
                    );
                    if lin1.angle(lin2) >= 0.0 {
                        me.linsol[me.nbr_sol - 1].pos.vdir.reverse();
                    }
                    me.well_done = true;
                }
            }
        }

        for i in 1..=me.nbr_sol {
            me.pntint1sol[i - 1] = me.linsol[i - 1].location();
            me.pntint2sol[i - 1] = me.pntint1sol[i - 1];
            me.par1sol[i - 1] = parameter_lin2d(&me.linsol[i - 1], &me.pntint1sol[i - 1]);
            me.par2sol[i - 1] = parameter_lin2d(&me.linsol[i - 1], &me.pntint2sol[i - 1]);
            me.pararg1[i - 1] = parameter_lin2d(lin1, &me.pntint1sol[i - 1]);
            me.pararg2[i - 1] = parameter_lin2d(lin2, &me.pntint2sol[i - 1]);
        }
        me
    }

    /// `GccAna_Lin2dBisec::IsDone` (`cxx:126-129`).
    pub fn is_done(&self) -> bool {
        self.well_done
    }

    /// `GccAna_Lin2dBisec::NbSolutions` (`cxx:131-138`).
    pub fn nb_solutions(&self) -> Result<usize, GccAnaBisecError> {
        if !self.well_done {
            return Err(GccAnaBisecError::NotDone);
        }
        Ok(self.nbr_sol)
    }

    /// `GccAna_Lin2dBisec::ThisSolution` (`cxx:140-153`).
    pub fn this_solution(&self, index: usize) -> Result<GpLin2d, GccAnaBisecError> {
        check_index(self.well_done, self.nbr_sol, index)?;
        Ok(self.linsol[index - 1])
    }

    /// `GccAna_Lin2dBisec::Intersection1` (`cxx:155-174`).
    pub fn intersection1(
        &self,
        index: usize,
    ) -> Result<(f64, f64, GpPnt2d), GccAnaBisecError> {
        check_index(self.well_done, self.nbr_sol, index)?;
        Ok((
            self.par1sol[index - 1],
            self.pararg1[index - 1],
            self.pntint1sol[index - 1],
        ))
    }

    /// `GccAna_Lin2dBisec::Intersection2` (`cxx:176-194`).
    pub fn intersection2(
        &self,
        index: usize,
    ) -> Result<(f64, f64, GpPnt2d), GccAnaBisecError> {
        check_index(self.well_done, self.nbr_sol, index)?;
        Ok((
            self.par2sol[index - 1],
            self.pararg2[index - 1],
            self.pntint2sol[index - 1],
        ))
    }
}

/// `if (!WellDone) throw StdFail_NotDone(); if (Index <= 0 || Index > NbrSol)
/// throw Standard_OutOfRange();`.
fn check_index(
    done: bool,
    nbr_sol: usize,
    index: usize,
) -> Result<(), GccAnaBisecError> {
    if !done {
        return Err(GccAnaBisecError::NotDone);
    }
    if index < 1 || index > nbr_sol {
        return Err(GccAnaBisecError::OutOfRange);
    }
    Ok(())
}
