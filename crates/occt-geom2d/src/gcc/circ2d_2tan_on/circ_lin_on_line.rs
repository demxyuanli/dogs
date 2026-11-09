//! `GccAna_Circ2d2TanOn(QualifiedCirc, QualifiedLin, gp_Lin2d, Tolerance)`
//! (`GccAna_Circ2d2TanOn_1.cxx:47-345`): circles tangent to a circle and a line
//! with the centre on a straight line.
//!
//! The boundary block (`cxx:92-200`) projects `center1` on `OnLine`, checks the
//! resulting `pinterm` against `L2` and, when the qualifiers agree, stores one
//! solution and returns. The general case (`cxx:205-345`) intersects each
//! `GccAna_CircLin2dBisec` branch of `(C1, L2)` with `OnLine`.
//!
//! Translation notes:
//! - `cxx:96-100` recomputes `pinterm` from the very same expression it was
//!   just built with, so the guarded branch has no effect; the port keeps the
//!   redundant assignment,
//! - the boundary block compares with `Tolerance` (`cxx:96`, `:104`, `:110`,
//!   `:129`, `:137`), the general case with `Tolerance` too (`cxx:237-249`,
//!   `:255`, `:277`); only `Tol = |Tolerance|` is used for the four-arm
//!   qualifier (`cxx:165-176`, `:294-306`),
//! - `cxx:269-277` sets `ok = true` again inside the `IsOutside` arm without
//!   clearing it first, as does `cxx:280-283` for `IsUnqualified`; the port
//!   keeps both no-ops,
//! - `dc2 = gp_Dir2d(sign * gp_XY(-dirL2.Y(), dirL2.X()))` (`cxx:190`, `:332`)
//!   raises `Standard_ConstructionError` in OCCT when `sign == 0`; the port's
//!   `expect` panics with the same message.

use occt_core::elib::clib2d::{parameter_circ2d, parameter_lin2d};
use occt_core::gp::circ2d::GpCirc2d;
use occt_core::gp::dir2d::GpDir2d;
use occt_core::gp::lin2d::GpLin2d;
use occt_core::intana2d::{IntAna2dAnaIntersection, IntAna2dConic};

use crate::gcc_ana_bisec::GccAnaCircLin2dBisec;
use crate::gcc_int::GccIntIType;

use super::super::gcc_ent::{GccEntQualifiedCirc, GccEntQualifiedLin};
use super::{
    axis_at, dir_between, left_normal, offset_xy, qualifier_of_circle, qualifier_of_line, side,
    GccAnaCirc2d2TanOn, MAX_SOL,
};

pub(crate) fn new(
    qualified1: &GccEntQualifiedCirc,
    qualified2: &GccEntQualifiedLin,
    on_line: &GpLin2d,
    tolerance: f64,
) -> GccAnaCirc2d2TanOn {
    let mut me = GccAnaCirc2d2TanOn::empty();

    if !(qualified1.is_enclosed()
        || qualified1.is_enclosing()
        || qualified1.is_outside()
        || qualified1.is_unqualified())
        || !(qualified2.is_enclosed()
            || qualified2.is_outside()
            || qualified2.is_unqualified())
    {
        panic!("GccEnt_BadQualifier in GccAna_Circ2d2TanOn");
    }

    let tol = tolerance.abs();
    let mut radius = 0.0f64;
    let mut ok = false;
    let c1 = qualified1.qualified();
    let l2 = qualified2.qualified();
    let r1 = c1.radius();
    let center1 = c1.location();
    let origin2 = l2.location();
    let dir_l2 = *l2.direction();
    let norm_l2 = left_normal(&dir_l2);

    // `cxx:92-200`: boundary case.
    let distcl = on_line.distance(&center1);
    let on_norm = left_normal(on_line.direction());
    let mut pinterm = offset_xy(&center1, distcl, &on_norm);
    if on_line.distance(&pinterm) > tolerance {
        pinterm = offset_xy(&center1, distcl, &on_norm);
    }
    let mut dist2 = l2.distance(&pinterm);
    if qualified1.is_enclosed() || qualified1.is_outside() {
        if (distcl - r1 - dist2).abs() <= tol {
            ok = true;
        }
    } else if qualified1.is_enclosing() {
        if (dist2 - distcl - r1).abs() <= tol {
            ok = true;
        }
    } else if qualified1.is_unqualified() {
        ok = true;
    } else {
        panic!("GccEnt_BadQualifier in GccAna_Circ2d2TanOn");
    }
    if ok {
        if qualified2.is_outside() {
            let pbid = offset_xy(&pinterm, dist2, &norm_l2);
            if l2.distance(&pbid) <= tol {
                me.well_done = true;
            }
        } else if qualified2.is_enclosed() {
            let pbid = offset_xy(&pinterm, -dist2, &norm_l2);
            if l2.distance(&pbid) <= tol {
                me.well_done = true;
            }
        } else if qualified2.is_unqualified() {
            me.well_done = false;
        } else {
            panic!("GccEnt_BadQualifier in GccAna_Circ2d2TanOn");
        }
    }
    if me.well_done {
        me.nbr_sol = 1;
        let n = me.nbr_sol;
        me.cirsol[n - 1] = GpCirc2d::new(axis_at(pinterm), dist2);

        let dc1 = dir_between(&pinterm, &center1);
        let dc2 = dir_between(&pinterm, &origin2);
        let distcc1 = pinterm.distance(&center1);
        me.qualifier1[n - 1] = qualifier_of_circle(qualified1, distcc1, dist2, r1, tol);
        me.qualifier2[n - 1] = qualifier_of_line(qualified2, dc2.dot(&norm_l2));

        let sign = dc2.dot(&norm_l2);
        let dc2 = signed_normal(sign, &norm_l2);
        me.pnttg1sol[n - 1] = offset_xy(&pinterm, dist2, &dc1);
        me.pnttg2sol[n - 1] = offset_xy(&pinterm, dist2, &dc2);
        me.par1sol[n - 1] = parameter_circ2d(&me.cirsol[n - 1], &me.pnttg1sol[n - 1]);
        me.pararg1[n - 1] = parameter_circ2d(&c1, &me.pnttg1sol[n - 1]);
        me.par2sol[n - 1] = parameter_circ2d(&me.cirsol[n - 1], &me.pnttg2sol[n - 1]);
        me.pararg2[n - 1] = parameter_lin2d(&l2, &me.pnttg2sol[n - 1]);
        me.pntcen[n - 1] = me.cirsol[n - 1].location();
        me.parcen3[n - 1] = parameter_lin2d(on_line, &me.pntcen[n - 1]);
        return me;
    }

    // `cxx:205-345`: general case.
    let bis = GccAnaCircLin2dBisec::new(&c1, &l2);
    if bis.is_done() {
        let nbsolution = bis.nb_solutions();
        for i in 1..=nbsolution {
            let sol = bis
                .this_solution(i)
                .expect("GccAna_CircLin2dBisec::ThisSolution");
            let typ = sol.arc_type();
            let mut intp = IntAna2dAnaIntersection::new();
            if typ == GccIntIType::Lin {
                intp.perform_lin_lin(on_line, &sol.line().unwrap());
            } else if typ == GccIntIType::Par {
                intp.perform_lin_conic(
                    on_line,
                    &IntAna2dConic::from_parab2d(&sol.parabola().unwrap()),
                );
            }
            if intp.is_done() {
                if !intp.is_empty() {
                    for j in 1..=intp.nb_points() {
                        let center = *intp.point(j).value();
                        let dist1 = center.distance(&center1);
                        dist2 = l2.distance(&center);
                        ok = false;
                        if qualified1.is_enclosed() {
                            if dist1 - r1 < tolerance {
                                if ((r1 - dist1).abs() - dist2).abs() < tolerance {
                                    ok = true;
                                }
                            }
                        } else if qualified1.is_outside() {
                            if r1 - dist1 < tolerance {
                                if ((r1 - dist1).abs() - dist2).abs() < tolerance {
                                    ok = true;
                                }
                            }
                        } else if qualified1.is_enclosing() || qualified1.is_unqualified() {
                            ok = true;
                        }
                        if qualified2.is_enclosed() && ok {
                            if side(&origin2, &center, &norm_l2) <= 0.0 {
                                ok = true;
                                radius = dist2;
                            }
                        } else if qualified2.is_outside() && ok {
                            if side(&origin2, &center, &norm_l2) >= 0.0 {
                                ok = true;
                                radius = dist2;
                            }
                        } else if qualified2.is_unqualified() && ok {
                            ok = true;
                            radius = dist2;
                        }
                        if ok {
                            if me.nbr_sol == MAX_SOL {
                                break;
                            }
                            me.nbr_sol += 1;
                            let n = me.nbr_sol;
                            me.cirsol[n - 1] = GpCirc2d::new(axis_at(center), radius);

                            let dc1 = dir_between(&center, &center1);
                            let dc2 = dir_between(&center, &origin2);
                            let distcc1 = center.distance(&center1);
                            me.qualifier1[n - 1] =
                                qualifier_of_circle(qualified1, distcc1, radius, r1, tol);
                            me.qualifier2[n - 1] =
                                qualifier_of_line(qualified2, dc2.dot(&norm_l2));

                            if center.distance(&center1) <= tolerance
                                && (radius - c1.radius()).abs() <= tolerance
                            {
                                me.the_same1[n - 1] = 1;
                            } else {
                                me.the_same1[n - 1] = 0;
                                me.pnttg1sol[n - 1] = offset_xy(&center, radius, &dc1);
                                me.par1sol[n - 1] =
                                    parameter_circ2d(&me.cirsol[n - 1], &me.pnttg1sol[n - 1]);
                                me.pararg1[n - 1] = parameter_circ2d(&c1, &me.pnttg1sol[n - 1]);
                            }
                            me.the_same2[n - 1] = 0;
                            let sign = dc2.dot(&norm_l2);
                            let dc2 = signed_normal(sign, &norm_l2);
                            me.pnttg2sol[n - 1] = offset_xy(&center, radius, &dc2);
                            me.par2sol[n - 1] =
                                parameter_circ2d(&me.cirsol[n - 1], &me.pnttg2sol[n - 1]);
                            me.pararg2[n - 1] = parameter_lin2d(&l2, &me.pnttg2sol[n - 1]);
                            me.pntcen[n - 1] = center;
                            me.parcen3[n - 1] = parameter_lin2d(on_line, &me.pntcen[n - 1]);
                        }
                    }
                }
                me.well_done = true;
            }
        }
    }
    me
}

/// `gp_Dir2d(sign * gp_XY(-dirL2.Y(), dirL2.X()))`.
fn signed_normal(sign: f64, normal: &GpDir2d) -> GpDir2d {
    GpDir2d::new(sign * normal.x, sign * normal.y).expect("gp_Dir2d from the signed normal")
}
