//! `GccAna_Circ2d3Tan(GccEnt_QualifiedCirc, GccEnt_QualifiedCirc, gp_Pnt2d,
//! Tolerance)` (`GccAna_Circ2d3Tan_4.cxx:41-465`): circles tangent to two
//! circles and through a point.
//!
//! Centers come from intersecting `GccAna_Circ2dBisec(C1, C2)` with
//! `GccAna_CircPnt2dBisec(C1, Point3)` (`cxx:94-183`); the sixteen
//! `(type1, type2)` dispatch arms of `cxx:106-185` are transcribed one by one.
//!
//! Translation notes:
//! - this is the only constructor with `static int MaxSol = 20` (`cxx:35`): the
//!   `NbrSol == MaxSol` breaks of `cxx:397-418` are kept, and the port arrays are
//!   [`super::MAX_SOL`] = 20 wide, so no separate cap is needed,
//! - after the loops OCCT re-derives each radius from `Point3` (`cxx:424-433`)
//!   and removes centres closer than `Tol^2` (`cxx:436-463`, the "debug"
//!   passes); both are ported as written, including the fact that the removal
//!   pass leaves a stale `ll` after `NbrSol--`,
//! - the removal pass copies the fourteen arrays of `cxx:444-460` but not
//!   `TheSame1..3`; the port copies exactly the same set.

use occt_core::elib::clib2d::parameter_circ2d;
use occt_core::gp::circ2d::GpCirc2d;
use occt_core::gp::dir2d::GpDir2d;
use occt_core::gp::pnt2d::GpPnt2d;
use occt_core::intana2d::{IntAna2dAnaIntersection, IntAna2dConic};

use crate::gcc_ana_bisec::{GccAnaCirc2dBisec, GccAnaCircPnt2dBisec};
use crate::gcc_int::GccIntIType;

use super::super::gcc_ent::{GccEntPosition, GccEntQualifiedCirc};
use super::{axis_at, GccAnaCirc2d3Tan};

/// `static int MaxSol = 20` (`cxx:35`).
const MAX_SOL: usize = 20;

pub(crate) fn new(
    qualified1: &GccEntQualifiedCirc,
    qualified2: &GccEntQualifiedCirc,
    point3: &GpPnt2d,
    tolerance: f64,
) -> GccAnaCirc2d3Tan {
    let mut me = GccAnaCirc2d3Tan::empty();
    let tol = tolerance.abs();

    if !((qualified1.is_enclosed() || qualified1.is_enclosing() || qualified1.is_outside()
        || qualified1.is_unqualified())
        && (qualified2.is_enclosed() || qualified2.is_enclosing() || qualified2.is_outside()
            || qualified2.is_unqualified()))
    {
        panic!("GccEnt_BadQualifier in GccAna_Circ2d3Tan");
    }

    let c1 = qualified1.qualified();
    let c2 = qualified2.qualified();
    let r1 = c1.radius();
    let r2 = c2.radius();
    let center1 = c1.location();
    let center2 = c2.location();

    let mut radius = [0.0f64; 3];

    let bis1 = GccAnaCirc2dBisec::new(&c1, &c2);
    let bis2 = GccAnaCircPnt2dBisec::new(&c1, point3);
    if bis1.is_done() && bis2.is_done() {
        let nb1 = bis1.nb_solutions().expect("GccAna_Circ2dBisec::NbSolutions");
        let nb2 = bis2.nb_solutions();
        for i in 1..=nb1 {
            let sol1 = bis1.this_solution(i).expect("GccAna_Circ2dBisec::ThisSolution");
            let typ1 = sol1.arc_type();
            let mut intp = IntAna2dAnaIntersection::new();
            for k in 1..=nb2 {
                let sol2 = bis2.this_solution(k).expect("GccAna_CircPnt2dBisec::ThisSolution");
                let typ2 = sol2.arc_type();
                match typ1 {
                    GccIntIType::Cir => match typ2 {
                        GccIntIType::Cir => intp
                            .perform_circ_circ(&sol1.circle().unwrap(), &sol2.circle().unwrap()),
                        GccIntIType::Lin => intp
                            .perform_lin_circ(&sol2.line().unwrap(), &sol1.circle().unwrap()),
                        GccIntIType::Hpr => intp.perform_circ_conic(
                            &sol1.circle().unwrap(),
                            &IntAna2dConic::from_hypr2d(&sol2.hyperbola().unwrap()),
                        ),
                        GccIntIType::Ell => intp.perform_circ_conic(
                            &sol1.circle().unwrap(),
                            &IntAna2dConic::from_elips2d(&sol2.ellipse().unwrap()),
                        ),
                        _ => {}
                    },
                    GccIntIType::Ell => match typ2 {
                        GccIntIType::Cir => intp.perform_circ_conic(
                            &sol2.circle().unwrap(),
                            &IntAna2dConic::from_elips2d(&sol1.ellipse().unwrap()),
                        ),
                        GccIntIType::Lin => intp.perform_lin_conic(
                            &sol2.line().unwrap(),
                            &IntAna2dConic::from_elips2d(&sol1.ellipse().unwrap()),
                        ),
                        GccIntIType::Hpr => intp.perform_elips_conic(
                            &sol1.ellipse().unwrap(),
                            &IntAna2dConic::from_hypr2d(&sol2.hyperbola().unwrap()),
                        ),
                        GccIntIType::Ell => intp.perform_elips_conic(
                            &sol1.ellipse().unwrap(),
                            &IntAna2dConic::from_elips2d(&sol2.ellipse().unwrap()),
                        ),
                        _ => {}
                    },
                    GccIntIType::Lin => match typ2 {
                        GccIntIType::Cir => intp
                            .perform_lin_circ(&sol1.line().unwrap(), &sol2.circle().unwrap()),
                        GccIntIType::Lin => intp
                            .perform_lin_lin(&sol1.line().unwrap(), &sol2.line().unwrap()),
                        GccIntIType::Hpr => intp.perform_lin_conic(
                            &sol1.line().unwrap(),
                            &IntAna2dConic::from_hypr2d(&sol2.hyperbola().unwrap()),
                        ),
                        GccIntIType::Ell => intp.perform_lin_conic(
                            &sol1.line().unwrap(),
                            &IntAna2dConic::from_elips2d(&sol2.ellipse().unwrap()),
                        ),
                        _ => {}
                    },
                    GccIntIType::Hpr => match typ2 {
                        GccIntIType::Cir => intp.perform_circ_conic(
                            &sol2.circle().unwrap(),
                            &IntAna2dConic::from_hypr2d(&sol1.hyperbola().unwrap()),
                        ),
                        GccIntIType::Lin => intp.perform_lin_conic(
                            &sol2.line().unwrap(),
                            &IntAna2dConic::from_hypr2d(&sol1.hyperbola().unwrap()),
                        ),
                        GccIntIType::Hpr => intp.perform_hypr_conic(
                            &sol1.hyperbola().unwrap(),
                            &IntAna2dConic::from_hypr2d(&sol2.hyperbola().unwrap()),
                        ),
                        GccIntIType::Ell => intp.perform_elips_conic(
                            &sol2.ellipse().unwrap(),
                            &IntAna2dConic::from_hypr2d(&sol1.hyperbola().unwrap()),
                        ),
                        _ => {}
                    },
                    _ => {}
                }
                if intp.is_done() {
                    if !intp.is_empty() {
                        for j in 1..=intp.nb_points() {
                            let mut rradius = 0.0;
                            let center = *intp.point(j).value();
                            let dist1 = center.distance(&center1);
                            let dist2 = center.distance(&center2);
                            let dist3 = center.distance(point3);
                            let mut nbsol1 = 0i32;
                            let mut nbsol2 = 0i32;
                            let mut nbsol3 = 0i32;
                            let mut ok = false;
                            if qualified1.is_enclosed() {
                                if dist1 - r1 < tolerance {
                                    radius[0] = (r1 - dist1).abs();
                                    nbsol1 = 1;
                                    ok = true;
                                }
                            } else if qualified1.is_outside() {
                                if r1 - dist1 < tolerance {
                                    radius[0] = (r1 - dist1).abs();
                                    nbsol1 = 1;
                                    ok = true;
                                }
                            } else if qualified1.is_enclosing() {
                                ok = true;
                                nbsol1 = 1;
                                radius[0] = r1 + dist1;
                            } else if qualified1.is_unqualified() {
                                ok = true;
                                nbsol1 = 2;
                                radius[0] = (r1 - dist1).abs();
                                radius[1] = r1 + dist1;
                            }
                            if qualified2.is_enclosed() && ok {
                                if dist2 - r2 < tolerance {
                                    for ii in 1..=nbsol1 {
                                        if (radius[ii as usize - 1] - (r2 - dist2).abs()) < tol {
                                            radius[0] = (r2 - dist2).abs();
                                            ok = true;
                                            nbsol2 = 1;
                                        }
                                    }
                                }
                            } else if qualified2.is_outside() && ok {
                                if r2 - dist2 < tolerance {
                                    for ii in 1..=nbsol1 {
                                        if (radius[ii as usize - 1] - (r2 - dist2).abs()) < tol {
                                            radius[0] = (r2 - dist2).abs();
                                            ok = true;
                                            nbsol2 = 1;
                                        }
                                    }
                                }
                            } else if qualified2.is_enclosing() && ok {
                                for ii in 1..=nbsol1 {
                                    if (radius[ii as usize - 1] - r2 - dist2).abs() < tol {
                                        radius[0] = r2 + dist2;
                                        ok = true;
                                        nbsol2 = 1;
                                    }
                                }
                            } else if qualified2.is_unqualified() && ok {
                                for ii in 1..=nbsol1 {
                                    if (radius[ii as usize - 1] - (r2 - dist2).abs()) < tol {
                                        rradius = (r2 - dist2).abs();
                                        ok = true;
                                        nbsol2 += 1;
                                    } else if (radius[ii as usize - 1] - r2 - dist2).abs() < tol {
                                        rradius = r2 + dist2;
                                        ok = true;
                                        nbsol2 += 1;
                                    }
                                }
                                if nbsol2 == 1 {
                                    radius[0] = rradius;
                                } else if nbsol2 == 2 {
                                    radius[0] = (r2 - dist2).abs();
                                    radius[1] = r2 + dist2;
                                }
                            }
                            for ii in 1..=nbsol2 {
                                if (dist3 - radius[ii as usize - 1]).abs() <= tol {
                                    nbsol3 += 1;
                                    ok = true;
                                }
                            }
                            if ok {
                                for k1 in 1..=nbsol3 {
                                    me.nbr_sol += 1;
                                    let n = me.nbr_sol;
                                    let rad = radius[k1 as usize - 1];
                                    me.cirsol[n - 1] = GpCirc2d::new(axis_at(center), rad);

                                    let distcc1 = center.distance(&center1);
                                    me.qualifier1[n - 1] = if !qualified1.is_unqualified() {
                                        qualified1.qualifier()
                                    } else if (distcc1 + rad - r1).abs() < tol {
                                        GccEntPosition::Enclosed
                                    } else if (distcc1 - r1 - rad).abs() < tol {
                                        GccEntPosition::Outside
                                    } else {
                                        GccEntPosition::Enclosing
                                    };

                                    let distcc2 = center.distance(&center2);
                                    me.qualifier2[n - 1] = if !qualified2.is_unqualified() {
                                        qualified2.qualifier()
                                    } else if (distcc2 + rad - r2).abs() < tol {
                                        GccEntPosition::Enclosed
                                    } else if (distcc2 - r2 - rad).abs() < tol {
                                        GccEntPosition::Outside
                                    } else {
                                        GccEntPosition::Enclosing
                                    };
                                    me.qualifier3[n - 1] = GccEntPosition::NoQualifier;

                                    if center.distance(&center1) <= tolerance
                                        && (rad - r1).abs() <= tolerance
                                    {
                                        me.the_same1[n - 1] = 1;
                                    } else {
                                        me.the_same1[n - 1] = 0;
                                        let mut dc = GpDir2d::new(
                                            center1.x() - center.x(),
                                            center1.y() - center.y(),
                                        )
                                        .expect("gp_Dir2d(center1.XY() - Center.XY())");
                                        if me.qualifier1[n - 1] == GccEntPosition::Enclosed {
                                            dc.reverse();
                                        }
                                        me.pnttg1sol[n - 1] = GpPnt2d::new(
                                            center.x() + rad * dc.x,
                                            center.y() + rad * dc.y,
                                        );
                                        me.par1sol[n - 1] =
                                            parameter_circ2d(&me.cirsol[n - 1], &me.pnttg1sol[n - 1]);
                                        me.pararg1[n - 1] =
                                            parameter_circ2d(&c1, &me.pnttg1sol[n - 1]);
                                    }

                                    if center.distance(&center2) <= tolerance
                                        && (rad - r2).abs() <= tolerance
                                    {
                                        me.the_same2[n - 1] = 1;
                                    } else {
                                        me.the_same2[n - 1] = 0;
                                        let dc = GpDir2d::new(
                                            center2.x() - center.x(),
                                            center2.y() - center.y(),
                                        )
                                        .expect("gp_Dir2d(center2.XY() - Center.XY())");
                                        let alpha =
                                            if center1.distance(&center2) <= tolerance {
                                                -1.0
                                            } else {
                                                1.0
                                            };
                                        me.pnttg2sol[n - 1] = GpPnt2d::new(
                                            center.x() + alpha * rad * dc.x,
                                            center.y() + alpha * rad * dc.y,
                                        );
                                        me.par2sol[n - 1] =
                                            parameter_circ2d(&me.cirsol[n - 1], &me.pnttg2sol[n - 1]);
                                        me.pararg2[n - 1] =
                                            parameter_circ2d(&c2, &me.pnttg2sol[n - 1]);
                                    }

                                    me.the_same3[n - 1] = 0;
                                    me.pnttg3sol[n - 1] = *point3;
                                    me.par3sol[n - 1] = parameter_circ2d(
                                        &me.cirsol[n - 1],
                                        &me.pnttg3sol[n - 1],
                                    );
                                    me.pararg3[n - 1] = 0.0;
                                    me.well_done = true;
                                    if me.nbr_sol == MAX_SOL {
                                        break;
                                    }
                                }
                            }
                        }
                    }
                    me.well_done = true;
                    if me.nbr_sol == MAX_SOL {
                        break;
                    }
                }
                if me.nbr_sol == MAX_SOL {
                    break;
                }
            }
            if me.nbr_sol == MAX_SOL {
                break;
            }
        }
    }

    // `cxx:424-433`: re-derive the radius from `Point3`.
    for kk in 1..=me.nbr_sol {
        let radius = me.cirsol[kk - 1].location().distance(point3);
        if (radius - me.cirsol[kk - 1].radius()).abs() > tol {
            me.cirsol[kk - 1].set_radius(radius);
        }
    }

    // `cxx:436-463`: drop duplicated centers.
    let tol2 = tol * tol;
    let mut kk = 1usize;
    while kk < me.nbr_sol {
        let pk = me.cirsol[kk - 1].location();
        let mut ll = kk + 1;
        while ll <= me.nbr_sol {
            let pl = me.cirsol[ll - 1].location();
            if pk.square_distance(&pl) < tol2 {
                let mut mm = ll + 1;
                while mm <= me.nbr_sol {
                    me.cirsol[mm - 2] = me.cirsol[mm - 1];
                    me.pnttg1sol[mm - 2] = me.pnttg1sol[mm - 1];
                    me.pnttg2sol[mm - 2] = me.pnttg2sol[mm - 1];
                    me.pnttg3sol[mm - 2] = me.pnttg3sol[mm - 1];
                    me.par1sol[mm - 2] = me.par1sol[mm - 1];
                    me.par2sol[mm - 2] = me.par2sol[mm - 1];
                    me.par3sol[mm - 2] = me.par3sol[mm - 1];
                    me.pararg1[mm - 2] = me.pararg1[mm - 1];
                    me.pararg2[mm - 2] = me.pararg2[mm - 1];
                    me.pararg3[mm - 2] = me.pararg3[mm - 1];
                    me.qualifier1[mm - 2] = me.qualifier1[mm - 1];
                    me.qualifier2[mm - 2] = me.qualifier2[mm - 1];
                    me.qualifier3[mm - 2] = me.qualifier3[mm - 1];
                    mm += 1;
                }
                me.nbr_sol -= 1;
            }
            ll += 1;
        }
        kk += 1;
    }
    me
}
