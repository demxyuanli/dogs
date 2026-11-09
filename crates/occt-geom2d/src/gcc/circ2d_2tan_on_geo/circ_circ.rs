//! `Geom2dGcc_Circ2d2TanOnGeo(QualifiedCirc, QualifiedCirc, OnCurv, Tolerance)`
//! (`Geom2dGcc_Circ2d2TanOnGeo.cxx:50-309`): the bisector of the two circles is
//! intersected with the carrier, and each intersection point yields the centres
//! of the candidate circles together with the radius pair
//! `(|dist - R|, dist + R)` selected by the qualifiers.
//!
//! Translation notes:
//! - `Rbid`, `RBid` and `Radius` are the function-level arrays of `cxx:74-76`;
//!   they are kept across the points of the `j` loop because OCCT declares them
//!   outside every loop,
//! - `cxx:199-200` writes `Rbid(1)` twice, so an unqualified first circle
//!   leaves `Rbid(2)` unwritten and the `isol` loop reads the zero-filled slot,
//! - `cxx:247-275` tests `Radius(i)` with the bisector index `i`, not the
//!   solution index `k`; kept as written,
//! - `double R1` is declared at `cxx:88` and reassigned inside the point loop
//!   (`cxx:172`), `double R2` likewise (`cxx:174`); both are the unchanged
//!   radii of the arguments.

use occt_core::elib::clib2d::parameter_circ2d;
use occt_core::gp::GpPnt2d;

use crate::adaptor_curve::Geom2dAdaptorCurve;
use crate::gcc::gcc_ent::{GccEntPosition, GccEntQualifiedCirc};
use crate::gcc_ana_bisec::GccAnaCirc2dBisec;
use crate::gcc_int::GccIntBisec;
use crate::geom2d_int::IntConicCurveGen;

use super::{
    carrier_domain, circle_at, circle_domain, dir_from, ellipse_domain, empty_domain,
    hyperbola_domain, par_on_circle, Geom2dGccCirc2d2TanOnGeo, MAX_SOL,
};

/// `Geom2dGcc_Circ2d2TanOnGeo::Geom2dGcc_Circ2d2TanOnGeo(QualifiedCirc,
/// QualifiedCirc, Geom2dAdaptor_Curve, double)` (`cxx:50-309`).
pub(super) fn new(
    qualified1: &GccEntQualifiedCirc,
    qualified2: &GccEntQualifiedCirc,
    on_curve: &Geom2dAdaptorCurve,
    tolerance: f64,
) -> Geom2dGccCirc2d2TanOnGeo {
    let mut s = Geom2dGccCirc2d2TanOnGeo::empty();
    let tol = tolerance.abs();
    // `NCollection_Array1<double> Rbid(1, 2); ... RBid(1, 2); ... Radius(1, 2);`
    // (`cxx:74-76`). `Radius` is written up to index 4.
    let mut rbid = [0.0_f64; MAX_SOL];
    let mut r_bid = [0.0_f64; MAX_SOL];
    let mut radius = [0.0_f64; MAX_SOL];
    if !(qualified1.is_enclosed()
        || qualified1.is_enclosing()
        || qualified1.is_outside()
        || qualified1.is_unqualified())
        || !(qualified2.is_enclosed()
            || qualified2.is_enclosing()
            || qualified2.is_outside()
            || qualified2.is_unqualified())
    {
        panic!("GccEnt_BadQualifier in Geom2dGcc_Circ2d2TanOnGeo");
    }
    let c1 = qualified1.qualified();
    let c2 = qualified2.qualified();
    let center1 = c1.location();
    let center2 = c2.location();
    let bis = GccAnaCirc2dBisec::new(&c1, &c2);
    if bis.is_done() {
        let mut inter = IntConicCurveGen::new();
        let nbsolution = bis
            .nb_solutions()
            .expect("GccAna_Circ2dBisec::NbSolutions");
        let (d2, curve) = carrier_domain(on_curve, tol);
        let tol1 = tol;
        let tol2 = tol1;
        for i in 1..=nbsolution {
            let sol = bis
                .this_solution(i)
                .expect("GccAna_Circ2dBisec::ThisSolution");
            match sol {
                GccIntBisec::Circle(circ) => {
                    let d1 = circle_domain(&circ, tol1, tol2);
                    inter.perform_circle(&circ, &d1, curve.as_ref(), &d2, tol1, tol2);
                }
                GccIntBisec::Ellipse(elips) => {
                    let d1 = ellipse_domain(&elips, tol1, tol2);
                    inter.perform_ellipse(&elips, &d1, curve.as_ref(), &d2, tol1, tol2);
                }
                GccIntBisec::Hyperbola(hypr) => {
                    let d1 = hyperbola_domain(&hypr, tol1, tol2);
                    inter.perform_hyperbola(&hypr, &d1, curve.as_ref(), &d2, tol1, tol2);
                }
                GccIntBisec::Line(line) => {
                    let d1 = empty_domain();
                    inter.perform_line(&line, &d1, curve.as_ref(), &d2, tol1, tol2);
                }
                GccIntBisec::Parabola(_) | GccIntBisec::Point(_) => {
                    panic!("Standard_ConstructionError in Geom2dGcc_Circ2d2TanOnGeo");
                }
            }
            if inter.is_done() && !inter.result().is_empty() {
                for j in 1..=inter.result().nb_points() {
                    let center = *inter.result().point(j).value();
                    let dist1 = center.distance(&center1);
                    let dist2 = center.distance(&center2);
                    let mut nbsol = 0_usize;
                    let mut nnsol = 0_usize;
                    let r1 = c1.radius();
                    let r2 = c2.radius();
                    if qualified1.is_enclosed() {
                        if dist1 - r1 < tol {
                            nbsol = 1;
                            rbid[0] = (r1 - dist1).abs();
                        }
                    } else if qualified1.is_outside() {
                        if r1 - dist1 < tol {
                            nbsol = 1;
                            rbid[0] = (dist1 - r1).abs();
                        }
                    } else if qualified1.is_enclosing() {
                        nbsol = 1;
                        rbid[0] = dist1 + r1;
                    } else if qualified1.is_unqualified() {
                        nbsol = 2;
                        rbid[0] = dist1 + r1;
                        rbid[0] = (dist1 - r1).abs();
                    }
                    if qualified2.is_enclosed() && nbsol != 0 {
                        if dist2 - r2 < tol {
                            r_bid[0] = (r2 - dist2).abs();
                        }
                    } else if qualified2.is_outside() && nbsol != 0 {
                        if r2 - dist2 < tol {
                            r_bid[0] = (r2 - dist2).abs();
                        }
                    } else if qualified2.is_enclosing() && nbsol != 0 {
                        r_bid[0] = dist2 + r2;
                    } else if qualified2.is_unqualified() && nbsol != 0 {
                        r_bid[0] = dist2 + r2;
                        r_bid[1] = (r2 - dist2).abs();
                    }
                    for isol in 1..=nbsol {
                        for jsol in 1..=nbsol {
                            if (rbid[isol - 1] - r_bid[jsol - 1]).abs() <= tol {
                                nnsol += 1;
                                radius[nnsol - 1] = (r_bid[jsol - 1] + rbid[isol - 1]) / 2.0;
                            }
                        }
                    }
                    if nnsol > 0 {
                        for k in 1..=nnsol {
                            if s.nbr_sol >= MAX_SOL {
                                break;
                            }
                            s.nbr_sol += 1;
                            let n = s.nbr_sol;
                            s.cirsol[n - 1] = circle_at(center, radius[k - 1]);
                            let distcc1 = center.distance(&center1);
                            let distcc2 = center.distance(&center2);
                            s.qualifier1[n - 1] = if !qualified1.is_unqualified() {
                                qualified1.qualifier()
                            } else if (distcc1 + radius[i as usize - 1] - r1).abs() < tol {
                                GccEntPosition::Enclosed
                            } else if (distcc1 - r1 - radius[i as usize - 1]).abs() < tol {
                                GccEntPosition::Outside
                            } else {
                                GccEntPosition::Enclosing
                            };
                            s.qualifier2[n - 1] = if !qualified2.is_unqualified() {
                                qualified2.qualifier()
                            } else if (distcc2 + radius[i as usize - 1] - r2).abs() < tol {
                                GccEntPosition::Enclosed
                            } else if (distcc2 - r2 - radius[i as usize - 1]).abs() < tol {
                                GccEntPosition::Outside
                            } else {
                                GccEntPosition::Enclosing
                            };
                            if dist1 <= tol && (radius[k - 1] - c1.radius()).abs() <= tol {
                                s.the_same1[n - 1] = 1;
                            } else {
                                s.the_same1[n - 1] = 0;
                                let dc1 = dir_from(&center1, &center);
                                let p = GpPnt2d::new(
                                    center.x() + radius[k - 1] * dc1.x(),
                                    center.y() + radius[k - 1] * dc1.y(),
                                );
                                let par1 = par_on_circle(&s.cirsol[n - 1], &p);
                                let arg1 = parameter_circ2d(&c1, &p);
                                s.pnttg1sol[n - 1] = p;
                                s.par1sol[n - 1] = par1;
                                s.pararg1[n - 1] = arg1;
                            }
                            if dist2 <= tol && (radius[k - 1] - c2.radius()).abs() <= tol {
                                s.the_same2[n - 1] = 1;
                            } else {
                                s.the_same2[n - 1] = 0;
                                let dc2 = dir_from(&center2, &center);
                                let p = GpPnt2d::new(
                                    center.x() + radius[k - 1] * dc2.x(),
                                    center.y() + radius[k - 1] * dc2.y(),
                                );
                                let par2 = par_on_circle(&s.cirsol[n - 1], &p);
                                let arg2 = parameter_circ2d(&c2, &p);
                                s.pnttg2sol[n - 1] = p;
                                s.par2sol[n - 1] = par2;
                                s.pararg2[n - 1] = arg2;
                            }
                            s.pntcen[n - 1] = center;
                            s.parcen3[n - 1] = inter.result().point(j).param_on_second();
                        }
                        s.well_done = true;
                    }
                }
            }
        }
    }
    s
}

