//! `Geom2dGcc_Circ2d2TanRadGeo(GccEnt_QualifiedCirc, Geom2dGcc_QCurve, Radius,
//! Tolerance)` (`Geom2dGcc_Circ2d2TanRadGeo.cxx:265-474`): the operand circle is
//! enlarged by `cote1(j)` and intersected with the curve offset by
//! `-cote2(j)`.
//!
//! Translation notes:
//! - `gp_Circ2d Circ(C1.XAxis(), R1 + cote1(jcote1))` (`cxx:403`) reuses the
//!   operand X axis, so `pararg1` is a parameter on `C1`,
//! - `D1` is closed and equivalent on `[0, 2*M_PI]` (`cxx:404-410`), the same
//!   domain builder as the `Circ2dTanOnRadGeo` carriers,
//! - the outer sweep is bounded by `NbrSol < aNbSolMAX` but the inner one is not
//!   (`cxx:401`, `:411`); the point loop is bounded again (`cxx:429`),
//! - `par1sol` / `par2sol` are evaluated last (`cxx:466-467`).

use occt_core::elib::clib2d::{circle_value_ax22d, parameter_circ2d};
use occt_core::gp::{GpAx22d, GpCirc2d};
use occt_core::intres2d::IntRes2dDomain;

use crate::gcc::curve_tool;
use crate::gcc::gcc_ent::{GccEntPosition, GccEntQualifiedCirc};
use crate::gcc::qualified_curve::Geom2dGccQCurve;
use crate::gcc::tan_geo_common::{circle_at, circle_domain, curve_domain, offset_of};
use crate::geom2d_int::IntConicCurveGen;

use super::{
    check_qualifier, check_radius, cotes, Categories, Geom2dGccCirc2d2TanRadGeo, MAX_SOL,
};

const WHO: &str = "Geom2dGcc_Circ2d2TanRadGeo";

/// `Geom2dGcc_Circ2d2TanRadGeo(GccEnt_QualifiedCirc, Geom2dGcc_QCurve, double,
/// double)` (`cxx:265-474`).
pub(super) fn new(
    qualified1: &GccEntQualifiedCirc,
    qualified2: &Geom2dGccQCurve,
    radius: f64,
    tolerance: f64,
) -> Geom2dGccCirc2d2TanRadGeo {
    let mut s = Geom2dGccCirc2d2TanRadGeo::empty();
    let tol = tolerance.abs();
    let cat1 = Categories::of_circ(qualified1);
    let cat2 = Categories::of_qcurve(qualified2);
    check_qualifier(cat1.all_of(), cat2.all_of(), WHO);
    let c1 = qualified1.qualified();
    let center1 = c1.location();
    let cu2 = qualified2.qualified();
    check_radius(radius, WHO);
    let (nbrcote1, cote1, nbrcote2, cote2) = cotes(cat1, cat2, radius);
    let r1 = c1.radius();
    let mut inter = IntConicCurveGen::new();
    for jcote1 in 1..=nbrcote1 {
        if s.nbr_sol >= MAX_SOL {
            break;
        }
        // `gp_Circ2d Circ(C1.XAxis(), R1 + cote1(jcote1));` (`cxx:403`).
        let x_axis = c1.x_axis();
        let circ = GpCirc2d::new(
            GpAx22d::from_xdir(*x_axis.location(), *x_axis.direction()),
            r1 + cote1[jcote1 - 1],
        );
        let d1: IntRes2dDomain = circle_domain(&circ, tol);
        for jcote2 in 1..=nbrcote2 {
            // `Adaptor2d_OffsetCurve C2(HCu2, -cote2(jcote2));` (`cxx:417-418`).
            let c2 = offset_of(cu2, -cote2[jcote2 - 1]);
            let d2 = curve_domain(c2.as_ref(), tol);
            inter.perform_circle(&circ, &d1, c2.as_ref(), &d2, tol, tol);
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
                        let distcc1 = center.distance(&center1);
                        if !qualified1.is_unqualified() {
                            s.qualifier1[n - 1] = qualified1.qualifier();
                        } else if (distcc1 + radius - r1).abs() < tol {
                            s.qualifier1[n - 1] = GccEntPosition::Enclosed;
                        } else if (distcc1 - r1 - radius).abs() < tol {
                            s.qualifier1[n - 1] = GccEntPosition::Outside;
                        } else {
                            s.qualifier1[n - 1] = GccEntPosition::Enclosing;
                        }
                        s.qualifier2[n - 1] = qualified2.qualifier();
                        s.the_same1[n - 1] = 0;
                        s.the_same2[n - 1] = 0;
                        s.pararg1[n - 1] = pt.param_on_first();
                        s.pararg2[n - 1] = pt.param_on_second();
                        // `ElCLib::Value(pararg1, C1)` (`cxx:464`).
                        s.pnttg1sol[n - 1] = circle_value_ax22d(s.pararg1[n - 1], c1.position(), r1);
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
