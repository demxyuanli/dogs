//! `Geom2dGcc_Circ2d2TanRadGeo(GccEnt_QualifiedLin, Geom2dGcc_QCurve, Radius,
//! Tolerance)` (`Geom2dGcc_Circ2d2TanRadGeo.cxx:51-249`): the line and the curve
//! are each displaced by their `cote` and intersected by
//! `Geom2dInt_TheIntConicCurveOfGInter`.
//!
//! Translation notes:
//! - `gp_Pnt2d Point(L1.Location().XY() + cote1(jcote1) * Dir.XY())`
//!   (`cxx:194`) moves the line along `Dir = (-y1dir, x1dir)`, the same
//!   direction as `normL1` (`cxx:105`),
//! - `IntRes2d_Domain D1;` (`cxx:195`) is the default (infinite) domain of the
//!   displaced line, while the curve side gets an explicit domain,
//! - the curve offset is `-cote2(jcote2)` (`cxx:200`), so the sign of the
//!   cascade is flipped once on each side,
//! - `par1sol` / `par2sol` are evaluated last (`cxx:240-241`), after both
//!   tangency points are stored.

use occt_core::elib::clib2d::{line_value_ax2d, parameter_circ2d};
use occt_core::gp::{GpLin2d, GpPnt2d};
use occt_core::intres2d::IntRes2dDomain;

use crate::gcc::curve_tool;
use crate::gcc::gcc_ent::{GccEntPosition, GccEntQualifiedLin};
use crate::gcc::qualified_curve::Geom2dGccQCurve;
use crate::gcc::tan_geo_common::{circle_at, curve_domain, offset_of};
use crate::geom2d_int::IntConicCurveGen;

use super::{
    check_qualifier, check_radius, cotes, Categories, Geom2dGccCirc2d2TanRadGeo, MAX_SOL,
};

const WHO: &str = "Geom2dGcc_Circ2d2TanRadGeo";

/// `Geom2dGcc_Circ2d2TanRadGeo(GccEnt_QualifiedLin, Geom2dGcc_QCurve, double,
/// double)` (`cxx:51-249`).
pub(super) fn new(
    qualified1: &GccEntQualifiedLin,
    qualified2: &Geom2dGccQCurve,
    radius: f64,
    tolerance: f64,
) -> Geom2dGccCirc2d2TanRadGeo {
    let mut s = Geom2dGccCirc2d2TanRadGeo::empty();
    let tol = tolerance.abs();
    let cat1 = Categories::of_lin(qualified1);
    let cat2 = Categories::of_qcurve(qualified2);
    check_qualifier(
        cat1.without_enclosing(),
        cat2.all_of(),
        WHO,
    );
    let l1 = qualified1.qualified();
    let x1dir = l1.direction().x();
    let y1dir = l1.direction().y();
    let lxloc = l1.location().x();
    let lyloc = l1.location().y();
    let origin1 = GpPnt2d::new(lxloc, lyloc);
    let norm_xy = occt_core::gp::xy::GpXY::new(-y1dir, x1dir);
    let cu2 = qualified2.qualified();
    check_radius(radius, WHO);
    let (nbrcote1, cote1, nbrcote2, cote2) = cotes(cat1, cat2, radius);
    let mut inter = IntConicCurveGen::new();
    for jcote1 in 1..=nbrcote1 {
        // `gp_Pnt2d Point(L1.Location().XY() + cote1(jcote1) * Dir.XY());`
        // `gp_Lin2d Line(Point, L1.Direction());` (`cxx:194-195`).
        let point = GpPnt2d::from_xy(
            l1.location()
                .xy()
                .added(&norm_xy.multiplied_scalar(cote1[jcote1 - 1])),
        );
        let line = GpLin2d::from_pnt_dir(point, *l1.direction());
        let d1 = IntRes2dDomain::new();
        for jcote2 in 1..=nbrcote2 {
            if s.nbr_sol >= MAX_SOL {
                break;
            }
            // `Adaptor2d_OffsetCurve C2(HCu2, -cote2(jcote2));` (`cxx:199-200`).
            let c2 = offset_of(cu2, -cote2[jcote2 - 1]);
            let d2 = curve_domain(c2.as_ref(), tol);
            inter.perform_line(&line, &d1, c2.as_ref(), &d2, tol, tol);
            if inter.is_done() {
                if !inter.result().is_empty() {
                    for i in 1..=inter.result().nb_points() {
                        if s.nbr_sol >= MAX_SOL {
                            break;
                        }
                        s.nbr_sol += 1;
                        let n = s.nbr_sol;
                        let pt = inter.result().point(i);
                        let center = *pt.value();
                        s.cirsol[n - 1] = circle_at(center, radius);
                        // `gp_Dir2d dc1(origin1.XY() - Center.XY());`
                        // (`cxx:221`).
                        let dc1 = origin1.xy().subtracted(center.xy());
                        s.qualifier2[n - 1] = qualified2.qualifier();
                        if !qualified1.is_unqualified() {
                            s.qualifier1[n - 1] = qualified1.qualifier();
                        } else if dc1.dot(&norm_xy) > 0.0 {
                            s.qualifier1[n - 1] = GccEntPosition::Outside;
                        } else {
                            s.qualifier1[n - 1] = GccEntPosition::Enclosed;
                        }
                        s.the_same1[n - 1] = 0;
                        s.the_same2[n - 1] = 0;
                        s.pararg1[n - 1] = pt.param_on_first();
                        s.pararg2[n - 1] = pt.param_on_second();
                        s.pnttg1sol[n - 1] = line_value_ax2d(s.pararg1[n - 1], l1.position());
                        s.pnttg2sol[n - 1] = curve_tool::value(cu2, s.pararg2[n - 1])
                            .expect("Geom2dAdaptor_Curve::Value on the curve operand");
                        s.par1sol[n - 1] = parameter_circ2d(&s.cirsol[n - 1], &s.pnttg1sol[n - 1]);
                        s.par2sol[n - 1] = parameter_circ2d(&s.cirsol[n - 1], &s.pnttg2sol[n - 1]);
                    }
                }
                s.well_done = true;
            }
        }
    }
    s
}
