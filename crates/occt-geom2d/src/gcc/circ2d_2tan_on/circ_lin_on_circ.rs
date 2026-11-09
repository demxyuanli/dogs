//! `GccAna_Circ2d2TanOn(QualifiedCirc, QualifiedLin, gp_Circ2d, Tolerance)`
//! (`GccAna_Circ2d2TanOn_7.cxx:34-327`): circles tangent to a circle and a line
//! with the centre on a circle.
//!
//! The boundary block (`cxx:78-181`) checks the two centres
//! `center1 +/- (distcco - Ron) * dircc` against `L2` and the `C1` qualifier.
//! The general case (`cxx:186-326`) intersects each `GccAna_CircLin2dBisec`
//! branch of `(C1, L2)` with `OnCirc`.
//!
//! Translation notes:
//! - `cxx:92-96` recomputes `pinterm` with the very same expression before
//!   repeating the same test, so the block can only fail; the port keeps the
//!   redundant recomputation,
//! - `cxx:120` and `cxx:123` both test `dirbid.Dot(normal) < 0.` (the second was
//!   presumably meant to be `>`); the port keeps the copy-paste test,
//! - the boundary block never assigns `pntcen` (`cxx:173`) yet evaluates
//!   `ElCLib::Parameter(OnCirc, pntcen)` from it; `pntcen` stays at the
//!   default-constructed `(0, 0)`, which the port reproduces by leaving the
//!   field untouched,
//! - `cxx:312-313` dots and multiplies by `gp_Dir2d(normL2.XY())` when
//!   `normL2` already is that direction; the port keeps the round trip.

use occt_core::elib::clib2d::{parameter_circ2d, parameter_lin2d};
use occt_core::gp::circ2d::GpCirc2d;
use occt_core::gp::dir2d::GpDir2d;
use occt_core::intana2d::{IntAna2dAnaIntersection, IntAna2dConic};

use crate::gcc_ana_bisec::GccAnaCircLin2dBisec;
use crate::gcc_int::GccIntIType;

use super::super::gcc_ent::GccEntQualifiedCirc;
use super::super::gcc_ent::GccEntQualifiedLin;
use super::{
    axis_at, dir_between, left_normal, offset_xy, qualifier_of_circle, qualifier_of_line, side,
    GccAnaCirc2d2TanOn, MAX_SOL,
};

pub(crate) fn new(
    qualified1: &GccEntQualifiedCirc,
    qualified2: &GccEntQualifiedLin,
    on_circ: &GpCirc2d,
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

    let c1 = qualified1.qualified();
    let l2 = qualified2.qualified();
    let r1 = c1.radius();
    let center1 = c1.location();
    let origin2 = l2.location();
    let dir_l2 = *l2.direction();
    let norm_l2 = left_normal(&dir_l2);

    // `cxx:78-181`: boundary case.
    let tol = tolerance.abs();
    let mut rradius = [0.0f64; MAX_SOL];
    let mut nbsol1 = 1i32;
    let ron = on_circ.radius();
    let distcco = on_circ.location().distance(&center1);
    let dircc = dir_between(&center1, &on_circ.location());
    let mut pinterm = offset_xy(&center1, distcco - ron, &dircc);
    let distpl2 = l2.distance(&pinterm);
    let mut distcc1 = pinterm.distance(&center1);
    let mut d1 = (distpl2 - (distcc1 - r1).abs()).abs();
    let mut d2 = (distpl2 - (distcc1 + r1)).abs();
    if d1 > tol || d2 > tol {
        pinterm = offset_xy(&center1, distcco - ron, &dircc);
        d1 = (distpl2 - (distcc1 - r1).abs()).abs();
        d2 = (distpl2 - (distcc1 + r1)).abs();
        if d1 > tol || d2 > tol {
            nbsol1 = 0;
        }
    }
    if nbsol1 > 0 {
        if qualified1.is_enclosed() || qualified1.is_outside() {
            nbsol1 = 1;
            rradius[0] = (distcc1 - r1).abs();
        } else if qualified1.is_enclosing() {
            nbsol1 = 1;
            rradius[0] = r1 + distcc1;
        } else if qualified1.is_unqualified() {
            nbsol1 = 2;
            rradius[0] = (distcc1 - r1).abs();
            rradius[1] = r1 + distcc1;
        }
        let dirbid = dir_between(&pinterm, &origin2);
        let normal = left_normal(&dir_l2);
        if qualified1.is_enclosed() && dirbid.dot(&normal) < 0.0 {
            nbsol1 = 0;
        } else if qualified1.is_outside() && dirbid.dot(&normal) < 0.0 {
            nbsol1 = 0;
        }
        for i in 1..=nbsol1 as usize {
            if (rradius[i - 1] - distpl2).abs() <= tol {
                me.well_done = true;
                me.nbr_sol += 1;
                let n = me.nbr_sol;
                me.cirsol[n - 1] = GpCirc2d::new(axis_at(pinterm), rradius[i - 1]);

                let dc1 = dir_between(&pinterm, &center1);
                let dc2 = dir_between(&pinterm, &origin2);
                distcc1 = pinterm.distance(&center1);
                me.qualifier1[n - 1] =
                    qualifier_of_circle(qualified1, distcc1, rradius[i - 1], r1, tol);
                me.qualifier2[n - 1] = qualifier_of_line(qualified2, dc2.dot(&norm_l2));

                let sign = dc2.dot(&norm_l2);
                let dc2 = signed_normal(sign, &norm_l2);
                me.pnttg1sol[n - 1] = offset_xy(&pinterm, rradius[i - 1], &dc1);
                me.pnttg2sol[n - 1] = offset_xy(&pinterm, rradius[i - 1], &dc2);
                me.par1sol[n - 1] = parameter_circ2d(&me.cirsol[n - 1], &me.pnttg1sol[n - 1]);
                me.pararg1[n - 1] = parameter_circ2d(&c1, &me.pnttg1sol[n - 1]);
                me.par2sol[n - 1] = parameter_circ2d(&me.cirsol[n - 1], &me.pnttg2sol[n - 1]);
                me.pararg2[n - 1] = parameter_lin2d(&l2, &me.pnttg2sol[n - 1]);
                // `cxx:173`: `pntcen` is never assigned here.
                me.parcen3[n - 1] = parameter_circ2d(on_circ, &me.pntcen[n - 1]);
            }
        }
        if me.well_done {
            return me;
        }
    }

    // `cxx:186-326`: general case.
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
                intp.perform_lin_circ(&sol.line().unwrap(), on_circ);
            } else if typ == GccIntIType::Par {
                intp.perform_circ_conic(
                    on_circ,
                    &IntAna2dConic::from_parab2d(&sol.parabola().unwrap()),
                );
            }
            if intp.is_done() {
                if !intp.is_empty()
                    && !intp.parallel_elements()
                    && !intp.identical_elements()
                {
                    for j in 1..=intp.nb_points() {
                        let center = *intp.point(j).value();
                        let dist1 = center.distance(&center1);
                        let dist2 = l2.distance(&center);
                        let mut ok = false;
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
                        let mut radius = 0.0f64;
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

                            let dc2 = dir_between(&center, &origin2);
                            distcc1 = center.distance(&center1);
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
                                let dc1 = dir_between(&center, &center1);
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
                            me.parcen3[n - 1] = parameter_circ2d(on_circ, &me.pntcen[n - 1]);
                        }
                    }
                }
                me.well_done = true;
            }
        }
    }
    me
}

/// `gp_Dir2d(sign * gp_XY(normal.X(), normal.Y()))`.
fn signed_normal(sign: f64, normal: &GpDir2d) -> GpDir2d {
    GpDir2d::new(sign * normal.x, sign * normal.y).expect("gp_Dir2d from the signed normal")
}
