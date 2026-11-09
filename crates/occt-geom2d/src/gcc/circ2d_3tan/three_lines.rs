//! `GccAna_Circ2d3Tan(GccEnt_QualifiedLin, GccEnt_QualifiedLin,
//! GccEnt_QualifiedLin, Tolerance)` (`GccAna_Circ2d3Tan_3.cxx:38-547`): circles
//! tangent to three straight lines.
//!
//! The first pass (`cxx:117-207`) intersects `GccAna_Lin2dBisec(L1, L2)` with
//! `GccAna_Lin2dBisec(L1, L3)` and stores every center with its mean distance as
//! radius. The `ncote1..3` chain (`cxx:210-492`) then decides, per qualifier
//! combination, which of the eight in/exterior corners survive; the second pass
//! (`cxx:498-547`) appends the surviving corners.
//!
//! Translation notes:
//! - the tolerance argument is unnamed in OCCT (`const double)` at `cxx:40`), so
//!   this constructor ignores it; the port keeps the parameter unused,
//! - `cxx:323-325` writes `cote2` twice instead of `cote3` inside the
//!   all-unqualified branch, `cxx:352-353` writes `cote3(1)` twice, `cxx:408-409`
//!   writes `cote3(1) = -1` twice, and `cxx:427-428` sets `cote1(1) = 1` for an
//!   `IsOutside()` circle; all four slips are kept as written,
//! - the second pass reads and writes the same `cirsol` array in ascending index
//!   order (`cxx:519-521`), so a later copy can read a slot an earlier copy
//!   overwrote; the port performs the same in-place sequence,
//! - `NbrSol` can exceed the `(1, 4)` upper bound of `cxx:48-61` in that second
//!   pass; the port stops at [`super::MAX_SOL`] (memory-safety guard only).

use occt_core::elib::clib2d::{parameter_circ2d, parameter_lin2d};
use occt_core::gp::circ2d::GpCirc2d;
use occt_core::gp::dir2d::GpDir2d;
use occt_core::gp::pnt2d::GpPnt2d;
use occt_core::intana2d::IntAna2dAnaIntersection;

use crate::gcc_ana_bisec::GccAnaLin2dBisec;

use super::super::gcc_ent::{GccEntPosition, GccEntQualifiedLin};
use super::{axis_at, GccAnaCirc2d3Tan};

pub(crate) fn new(
    qualified1: &GccEntQualifiedLin,
    qualified2: &GccEntQualifiedLin,
    qualified3: &GccEntQualifiedLin,
) -> GccAnaCirc2d3Tan {
    let mut me = GccAnaCirc2d3Tan::empty();

    if !((qualified1.is_enclosed() || qualified1.is_outside() || qualified1.is_unqualified())
        && (qualified2.is_enclosed() || qualified2.is_outside() || qualified2.is_unqualified())
        && (qualified3.is_enclosed() || qualified3.is_outside() || qualified3.is_unqualified()))
    {
        panic!("GccEnt_BadQualifier in GccAna_Circ2d3Tan");
    }

    let l1 = qualified1.qualified();
    let l2 = qualified2.qualified();
    let l3 = qualified3.qualified();
    let origin1 = l1.location();
    let dir1 = *l1.direction();
    let norm_l1 = GpDir2d::new(-dir1.y, dir1.x).expect("gp_Dir2d(-dir.Y(), dir.X())");
    let origin2 = l2.location();
    let dir2 = *l2.direction();
    let norm_l2 = GpDir2d::new(-dir2.y, dir2.x).expect("gp_Dir2d(-dir.Y(), dir.X())");
    let origin3 = l3.location();
    let dir3 = *l3.direction();
    let norm_l3 = GpDir2d::new(-dir3.y, dir3.x).expect("gp_Dir2d(-dir.Y(), dir.X())");

    let bisec1 = GccAnaLin2dBisec::new(&l1, &l2);
    let bisec2 = GccAnaLin2dBisec::new(&l1, &l3);

    let mut ncote1 = 0i32;
    let mut ncote2 = 0i32;
    let mut ncote3 = 0i32;
    let mut cote1 = [0.0f64; 3];
    let mut cote2 = [0.0f64; 3];
    let mut cote3 = [0.0f64; 3];
    let mut nbsol = 0usize;

    if bisec1.is_done() && bisec2.is_done() {
        let nb1 = bisec1.nb_solutions().expect("GccAna_Lin2dBisec::NbSolutions");
        let nb2 = bisec2.nb_solutions().expect("GccAna_Lin2dBisec::NbSolutions");
        for i in 1..=nb1 {
            for j in 1..=nb2 {
                let mut intp = IntAna2dAnaIntersection::new();
                intp.perform_lin_lin(
                    &bisec1.this_solution(i).expect("GccAna_Lin2dBisec::ThisSolution"),
                    &bisec2.this_solution(j).expect("GccAna_Lin2dBisec::ThisSolution"),
                );
                if intp.is_done() && !intp.is_empty() {
                    for k in 1..=intp.nb_points() {
                        nbsol += 1;
                        let center = *intp.point(k).value();
                        let p = center;
                        let radius = (l1.distance(&p) + l2.distance(&p) + l3.distance(&p)) / 3.0;
                        let cx = center.x();
                        let cy = center.y();
                        me.cirsol[nbsol - 1] = GpCirc2d::new(axis_at(center), radius);

                        let dc1 = GpDir2d::new(origin1.x() - center.x(), origin1.y() - center.y())
                            .expect("gp_Dir2d(origin1.XY() - Center.XY())");
                        me.qualifier1[nbsol - 1] = if !qualified1.is_unqualified() {
                            qualified1.qualifier()
                        } else if dc1.dot(&norm_l1) > 0.0 {
                            GccEntPosition::Outside
                        } else {
                            GccEntPosition::Enclosed
                        };
                        let dc2 = GpDir2d::new(origin2.x() - center.x(), origin2.y() - center.y())
                            .expect("gp_Dir2d(origin2.XY() - Center.XY())");
                        me.qualifier2[nbsol - 1] = if !qualified2.is_unqualified() {
                            qualified2.qualifier()
                        } else if dc2.dot(&norm_l2) > 0.0 {
                            GccEntPosition::Outside
                        } else {
                            GccEntPosition::Enclosed
                        };
                        let dc3 = GpDir2d::new(origin3.x() - center.x(), origin3.y() - center.y())
                            .expect("gp_Dir2d(origin3.XY() - Center.XY())");
                        me.qualifier3[nbsol - 1] = if !qualified3.is_unqualified() {
                            qualified3.qualifier()
                        } else if dc3.dot(&norm_l3) > 0.0 {
                            GccEntPosition::Outside
                        } else {
                            GccEntPosition::Enclosed
                        };

                        let mut cross1 = (-dir1.y) * (origin1.x() - cx) + dir1.x * (origin1.y() - cy);
                        let mut cross2 = (-dir2.y) * (origin2.x() - cx) + dir2.x * (origin2.y() - cy);
                        let mut cross3 = (-dir3.y) * (origin3.x() - cx) + dir3.x * (origin3.y() - cy);
                        if cross1 != 0.0 {
                            cross1 /= cross1.abs();
                        }
                        me.pnttg1sol[nbsol - 1] = GpPnt2d::new(
                            cx + cross1 * radius * (-dir1.y),
                            cy + cross1 * radius * dir1.x,
                        );
                        if cross2 != 0.0 {
                            cross2 /= cross2.abs();
                        }
                        me.pnttg2sol[nbsol - 1] = GpPnt2d::new(
                            cx + cross2 * radius * (-dir2.y),
                            cy + cross2 * radius * dir2.x,
                        );
                        if cross3 != 0.0 {
                            cross3 /= cross3.abs();
                        }
                        me.pnttg3sol[nbsol - 1] = GpPnt2d::new(
                            cx + cross3 * radius * (-dir3.y),
                            cy + cross3 * radius * dir3.x,
                        );
                        me.par1sol[nbsol - 1] =
                            parameter_circ2d(&me.cirsol[nbsol - 1], &me.pnttg1sol[nbsol - 1]);
                        me.pararg1[nbsol - 1] = parameter_lin2d(&l1, &me.pnttg1sol[nbsol - 1]);
                        me.par2sol[nbsol - 1] =
                            parameter_circ2d(&me.cirsol[nbsol - 1], &me.pnttg2sol[nbsol - 1]);
                        me.pararg2[nbsol - 1] = parameter_lin2d(&l2, &me.pnttg2sol[nbsol - 1]);
                        me.par3sol[nbsol - 1] =
                            parameter_circ2d(&me.cirsol[nbsol - 1], &me.pnttg3sol[nbsol - 1]);
                        me.pararg3[nbsol - 1] = parameter_lin2d(&l3, &me.pnttg3sol[nbsol - 1]);
                    }
                    me.well_done = true;
                }
            }
        }
    }

    // `cxx:210-492`: the eight fixed corner combinations, then the unqualified
    // fallback tree.
    if qualified1.is_enclosed() && qualified2.is_enclosed() && qualified3.is_enclosed() {
        ncote1 = 1;
        ncote2 = 1;
        ncote3 = 1;
        cote1[0] = 1.0;
        cote2[0] = 1.0;
        cote3[0] = 1.0;
    } else if qualified1.is_enclosed() && qualified2.is_enclosed() && qualified3.is_outside() {
        ncote1 = 1;
        ncote2 = 1;
        ncote3 = 1;
        cote1[0] = 1.0;
        cote2[0] = 1.0;
        cote3[0] = -1.0;
    } else if qualified1.is_enclosed() && qualified2.is_outside() && qualified3.is_enclosed() {
        ncote1 = 1;
        ncote2 = 1;
        ncote3 = 1;
        cote1[0] = 1.0;
        cote2[0] = -1.0;
        cote3[0] = 1.0;
    } else if qualified1.is_enclosed() && qualified2.is_outside() && qualified3.is_outside() {
        ncote1 = 1;
        ncote2 = 1;
        ncote3 = 1;
        cote1[0] = 1.0;
        cote2[0] = -1.0;
        cote3[0] = -1.0;
    } else if qualified1.is_outside() && qualified2.is_enclosed() && qualified3.is_enclosed() {
        ncote1 = 1;
        ncote2 = 1;
        ncote3 = 1;
        cote1[0] = -1.0;
        cote2[0] = 1.0;
        cote3[0] = 1.0;
    } else if qualified1.is_outside() && qualified2.is_enclosed() && qualified3.is_outside() {
        ncote1 = 1;
        ncote2 = 1;
        ncote3 = 1;
        cote1[0] = -1.0;
        cote2[0] = 1.0;
        cote3[0] = -1.0;
    } else if qualified1.is_outside() && qualified2.is_outside() && qualified3.is_enclosed() {
        ncote1 = 1;
        ncote2 = 1;
        ncote3 = 1;
        cote1[0] = -1.0;
        cote2[0] = -1.0;
        cote3[0] = 1.0;
    } else if qualified1.is_outside() && qualified2.is_outside() && qualified3.is_outside() {
        ncote1 = 1;
        ncote2 = 1;
        ncote3 = 1;
        cote1[0] = -1.0;
        cote2[0] = -1.0;
        cote3[0] = -1.0;
    } else if qualified1.is_unqualified() {
        ncote1 = 2;
        cote1[0] = 1.0;
        cote1[1] = -1.0;
        if qualified2.is_unqualified() {
            ncote2 = 2;
            cote2[0] = 1.0;
            cote2[1] = -1.0;
            if qualified3.is_unqualified() {
                ncote3 = 2;
                // `cxx:323-325`: OCCT writes `cote2`, not `cote3`.
                cote2[0] = 1.0;
                cote2[1] = -1.0;
                me.nbr_sol = nbsol;
                me.well_done = true;
            } else if qualified3.is_enclosed() {
                ncote3 = 1;
                cote3[0] = 1.0;
            } else if qualified3.is_outside() {
                ncote3 = 1;
                cote3[0] = -1.0;
            }
        } else if qualified2.is_enclosed() {
            ncote2 = 1;
            cote2[0] = 1.0;
            if qualified3.is_unqualified() {
                ncote3 = 2;
                // `cxx:352-353`: both writes target `cote3(1)`.
                cote3[0] = 1.0;
                cote3[0] = -1.0;
            } else if qualified3.is_enclosed() {
                ncote3 = 1;
                cote3[0] = 1.0;
            } else if qualified3.is_outside() {
                ncote3 = 1;
                cote3[0] = -1.0;
            }
        } else if qualified2.is_outside() {
            ncote2 = 1;
            cote2[0] = -1.0;
            if qualified3.is_unqualified() {
                ncote3 = 2;
                cote3[0] = 1.0;
                cote3[1] = -1.0;
            } else if qualified3.is_enclosed() {
                ncote3 = 1;
                cote3[0] = 1.0;
            } else if qualified3.is_outside() {
                ncote3 = 1;
                cote3[0] = -1.0;
            }
        }
    } else if qualified2.is_unqualified() {
        ncote2 = 2;
        cote2[0] = 1.0;
        cote2[1] = -1.0;
        if qualified1.is_enclosed() {
            ncote1 = 1;
            cote1[0] = 1.0;
            if qualified3.is_unqualified() {
                ncote3 = 2;
                // `cxx:408-409`: both writes target `cote3(1)`.
                cote3[0] = -1.0;
                cote3[0] = -1.0;
            } else if qualified3.is_enclosed() {
                ncote3 = 1;
                cote3[0] = 1.0;
            } else if qualified3.is_outside() {
                ncote3 = 1;
                cote3[0] = -1.0;
            }
        } else if qualified1.is_outside() {
            ncote1 = 1;
            // `cxx:427-428`: `+1` despite `IsOutside()`.
            cote1[0] = 1.0;
            if qualified3.is_unqualified() {
                ncote3 = 2;
                cote3[0] = 1.0;
                cote3[1] = -1.0;
            } else if qualified3.is_enclosed() {
                ncote3 = 1;
                cote3[0] = 1.0;
            } else if qualified3.is_outside() {
                ncote3 = 1;
                cote3[0] = -1.0;
            }
        }
    } else if qualified3.is_unqualified() {
        ncote3 = 2;
        cote3[0] = 1.0;
        cote3[1] = -1.0;
        if qualified1.is_enclosed() {
            ncote1 = 1;
            cote1[0] = 1.0;
            if qualified2.is_enclosed() {
                ncote2 = 1;
                cote2[0] = 1.0;
            } else if qualified2.is_outside() {
                ncote2 = 1;
                cote2[0] = -1.0;
            }
        } else if qualified1.is_outside() {
            ncote1 = 1;
            cote1[0] = -1.0;
            if qualified2.is_enclosed() {
                ncote2 = 1;
                cote2[0] = 1.0;
            } else if qualified2.is_outside() {
                ncote2 = 1;
                cote2[0] = -1.0;
            }
        }
    }

    if me.nbr_sol > 0 {
        return me;
    }

    for i in 1..=nbsol {
        for j1 in 1..=ncote1 {
            for j2 in 1..=ncote2 {
                for j3 in 1..=ncote3 {
                    let loc = me.cirsol[i - 1].location();
                    if cote2[j2 as usize - 1]
                        * ((loc.x() - origin2.x()) * (-dir2.y)
                            + (loc.y() - origin2.y()) * (dir2.x))
                        > 0.0
                        && cote3[j3 as usize - 1]
                            * ((loc.x() - origin3.x()) * (-dir3.y)
                                + (loc.y() - origin3.y()) * (dir3.x))
                            > 0.0
                        && cote1[j1 as usize - 1]
                            * ((loc.x() - origin1.x()) * (-dir1.y)
                                + (loc.y() - origin1.y()) * (dir1.x))
                            > 0.0
                    {
                        if me.nbr_sol >= super::MAX_SOL {
                            break;
                        }
                        me.nbr_sol += 1;
                        let n = me.nbr_sol;
                        me.cirsol[n - 1] = me.cirsol[i - 1];
                        let radius = me.cirsol[n - 1].radius();
                        let center = me.cirsol[n - 1].location();

                        let dc = GpDir2d::new(origin1.x() - center.x(), origin1.y() - center.y())
                            .expect("gp_Dir2d(origin1.XY() - Center.XY())");
                        let sign = dc.dot(&norm_l1);
                        let dc = GpDir2d::new(sign * norm_l1.x, sign * norm_l1.y)
                            .expect("gp_Dir2d(sign * normL1.XY())");
                        me.pnttg1sol[n - 1] =
                            GpPnt2d::new(center.x() + radius * dc.x, center.y() + radius * dc.y);
                        me.par1sol[n - 1] =
                            parameter_circ2d(&me.cirsol[n - 1], &me.pnttg1sol[n - 1]);
                        me.pararg1[n - 1] = parameter_lin2d(&l1, &me.pnttg1sol[n - 1]);

                        let dc = GpDir2d::new(origin2.x() - center.x(), origin2.y() - center.y())
                            .expect("gp_Dir2d(origin2.XY() - Center.XY())");
                        let sign = dc.dot(&norm_l2);
                        let dc = GpDir2d::new(sign * norm_l2.x, sign * norm_l2.y)
                            .expect("gp_Dir2d(sign * normL2.XY())");
                        me.pnttg2sol[n - 1] =
                            GpPnt2d::new(center.x() + radius * dc.x, center.y() + radius * dc.y);
                        me.par2sol[n - 1] =
                            parameter_circ2d(&me.cirsol[n - 1], &me.pnttg2sol[n - 1]);
                        me.pararg2[n - 1] = parameter_lin2d(&l2, &me.pnttg2sol[n - 1]);

                        let dc = GpDir2d::new(origin3.x() - center.x(), origin3.y() - center.y())
                            .expect("gp_Dir2d(origin3.XY() - Center.XY())");
                        let sign = dc.dot(&norm_l3);
                        let dc = GpDir2d::new(sign * norm_l3.x, sign * norm_l3.y)
                            .expect("gp_Dir2d(sign * normL3.XY())");
                        me.pnttg3sol[n - 1] =
                            GpPnt2d::new(center.x() + radius * dc.x, center.y() + radius * dc.y);
                        me.par3sol[n - 1] =
                            parameter_circ2d(&me.cirsol[n - 1], &me.pnttg3sol[n - 1]);
                        me.pararg3[n - 1] = parameter_lin2d(&l3, &me.pnttg3sol[n - 1]);
                    }
                }
            }
        }
    }
    me
}
