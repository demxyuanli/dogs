//! `Geom2dGcc_Circ2d2TanOnGeo(QualifiedLin, gp_Pnt2d, OnCurv, Tolerance)`
//! (`Geom2dGcc_Circ2d2TanOnGeo.cxx:912-1055`): the bisector of the line and the
//! point (a line or a parabola, `cxx:963-983`) is intersected with the carrier;
//! the radius is the distance from the centre to the line.
//!
//! Translation notes:
//! - `parabola_domain` keeps `Tol1` as the last tolerance (`cxx:972-978`),
//! - `pararg2` is zero and `pnttg2sol` is the argument point (`cxx:1046-1048`),
//! - `WellDone = true` (`cxx:1051`) is inside `if (Intp.IsDone())` and outside
//!   `if (!Intp.IsEmpty())`,
//! - the switch has no `GccInt_Cir` / `GccInt_Ell` / `GccInt_Hpr` arm, so those
//!   arc types reach `default:` and throw `Standard_ConstructionError`.

use occt_core::elib::clib2d::parameter_circ2d;
use occt_core::gp::{GpDir2d, GpPnt2d};

use crate::adaptor_curve::Geom2dAdaptorCurve;
use crate::gcc::gcc_ent::{GccEntPosition, GccEntQualifiedLin};
use crate::gcc_ana_bisec::GccAnaLinPnt2dBisec;
use crate::gcc_int::GccIntBisec;
use crate::geom2d_int::IntConicCurveGen;

use super::{
    carrier_domain, circle_at, dir_from, empty_domain, par_on_circle, par_on_line, parabola_domain,
    signed_normal, Geom2dGccCirc2d2TanOnGeo, MAX_SOL,
};

/// `Geom2dGcc_Circ2d2TanOnGeo::Geom2dGcc_Circ2d2TanOnGeo(QualifiedLin,
/// gp_Pnt2d, Geom2dAdaptor_Curve, double)` (`cxx:912-1055`).
pub(super) fn new(
    qualified1: &GccEntQualifiedLin,
    point2: &GpPnt2d,
    on_curve: &Geom2dAdaptorCurve,
    tolerance: f64,
) -> Geom2dGccCirc2d2TanOnGeo {
    let mut s = Geom2dGccCirc2d2TanOnGeo::empty();
    let tol = tolerance.abs();
    if !(qualified1.is_enclosed() || qualified1.is_outside() || qualified1.is_unqualified()) {
        panic!("GccEnt_BadQualifier in Geom2dGcc_Circ2d2TanOnGeo");
    }
    let l1 = qualified1.qualified();
    let origin1 = l1.location();
    let dir1 = *l1.direction();
    let normal = GpDir2d::new(-dir1.y(), dir1.x()).expect("gp_Dir2d of the line normal");
    let bis = GccAnaLinPnt2dBisec::new(&l1, point2);
    if bis.is_done() {
        let tol1 = tol;
        let tol2 = tol1;
        let mut inter = IntConicCurveGen::new();
        let (d2, curve) = carrier_domain(on_curve, tol);
        let sol = bis
            .this_solution()
            .expect("GccAna_LinPnt2dBisec::ThisSolution");
        match sol {
            GccIntBisec::Line(line) => {
                let d1 = empty_domain();
                inter.perform_line(&line, &d1, curve.as_ref(), &d2, tol1, tol2);
            }
            GccIntBisec::Parabola(parab) => {
                let d1 = parabola_domain(&parab, tol1);
                inter.perform_parabola(&parab, &d1, curve.as_ref(), &d2, tol1, tol2);
            }
            GccIntBisec::Circle(_)
            | GccIntBisec::Ellipse(_)
            | GccIntBisec::Hyperbola(_)
            | GccIntBisec::Point(_) => {
                panic!("Standard_ConstructionError in Geom2dGcc_Circ2d2TanOnGeo");
            }
        }
        if inter.is_done() {
            if !inter.result().is_empty() {
                for j in 1..=inter.result().nb_points() {
                    if s.nbr_sol >= MAX_SOL {
                        break;
                    }
                    let center = *inter.result().point(j).value();
                    let radius = l1.distance(&center);
                    let mut ok = false;
                    if qualified1.is_enclosed() {
                        if (origin1.x() - center.x()) * -dir1.y()
                            + (origin1.y() - center.y()) * dir1.x()
                            <= 0.0
                        {
                            ok = true;
                        }
                    } else if qualified1.is_outside() {
                        if (origin1.x() - center.x()) * -dir1.y()
                            + (origin1.y() - center.y()) * dir1.x()
                            >= 0.0
                        {
                            ok = true;
                        }
                    } else if qualified1.is_unqualified() {
                        ok = true;
                    }
                    if ok {
                        s.nbr_sol += 1;
                        let n = s.nbr_sol;
                        s.cirsol[n - 1] = circle_at(center, radius);
                        s.qualifier2[n - 1] = GccEntPosition::NoQualifier;
                        let dc2 = dir_from(&origin1, &center);
                        s.qualifier1[n - 1] = if !qualified1.is_unqualified() {
                            qualified1.qualifier()
                        } else if dc2.dot(&normal) > 0.0 {
                            GccEntPosition::Outside
                        } else {
                            GccEntPosition::Enclosed
                        };
                        s.the_same1[n - 1] = 0;
                        s.the_same2[n - 1] = 0;
                        let dc1 = dir_from(&origin1, &center);
                        let sign = dc1.dot(&normal);
                        let dc1 = signed_normal(sign, &dir1);
                        let p1 = GpPnt2d::new(
                            center.x() + radius * dc1.x(),
                            center.y() + radius * dc1.y(),
                        );
                        let par1 = par_on_circle(&s.cirsol[n - 1], &p1);
                        let arg1 = par_on_line(&l1, &p1);
                        s.pnttg1sol[n - 1] = p1;
                        s.par1sol[n - 1] = par1;
                        s.pararg1[n - 1] = arg1;
                        s.pnttg2sol[n - 1] = *point2;
                        s.par2sol[n - 1] = parameter_circ2d(&s.cirsol[n - 1], point2);
                        s.pararg2[n - 1] = 0.0;
                        s.pntcen[n - 1] = center;
                        s.parcen3[n - 1] = inter.result().point(j).param_on_second();
                    }
                }
            }
            s.well_done = true;
        }
    }
    s
}
