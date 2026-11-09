//! `Geom2dGcc_Circ2d2TanRadGeo(Geom2dGcc_QCurve, const gp_Pnt2d&, Radius,
//! Tolerance)` (`Geom2dGcc_Circ2d2TanRadGeo.cxx:490-608`): the locus of the
//! centres at distance `Radius` from the point is the point's circle, which is
//! intersected with the curve offset by `-cote1(j)`.
//!
//! Translation notes:
//! - `gp_Circ2d Circ(gp_Ax2d(Point2, gp_Dir2d::D::X), Radius)` (`cxx:559`), and
//!   `D1` spans `[0, M_PI + M_PI]` with equivalent parameters (`cxx:560-564`),
//! - the intersection point is parameterised on the **curve** first
//!   (`ParamOnSecond`), the reversed side of the `Line`/`Circle` carriers, so
//!   `pararg1` takes `ParamOnSecond` and `pararg2` is fixed at `0.0`
//!   (`cxx:596-597`),
//! - `qualifier2` is `GccEnt_noqualifier` (`cxx:593`), the second operand being
//!   a point,
//! - `par1sol` / `par2sol` are evaluated last (`cxx:600-601`).

use occt_core::elib::clib2d::parameter_circ2d;
use occt_core::gp::GpPnt2d;

use crate::gcc::curve_tool;
use crate::gcc::gcc_ent::GccEntPosition;
use crate::gcc::qualified_curve::Geom2dGccQCurve;
use crate::gcc::tan_geo_common::{circle_at, circle_domain, curve_domain, offset_of};
use crate::geom2d_int::IntConicCurveGen;

use super::{
    check_qualifier, check_radius, Geom2dGccCirc2d2TanRadGeo, MAX_SOL,
};

const WHO: &str = "Geom2dGcc_Circ2d2TanRadGeo";

/// The `cote1` cascade of a single operand (`cxx:540-556`).
fn cotes(qualified1: &Geom2dGccQCurve, radius: f64) -> (usize, [f64; 2]) {
    let mut nbrcote1 = 0;
    let mut cote1 = [0.0f64; 2];
    if qualified1.is_enclosed() {
        nbrcote1 = 1;
        cote1[0] = radius;
    } else if qualified1.is_outside() {
        nbrcote1 = 1;
        cote1[0] = -radius;
    } else if qualified1.is_unqualified() {
        nbrcote1 = 2;
        cote1[0] = radius;
        cote1[1] = -radius;
    }
    (nbrcote1, cote1)
}

/// `Geom2dGcc_Circ2d2TanRadGeo(Geom2dGcc_QCurve, const gp_Pnt2d&, double,
/// double)` (`cxx:490-608`).
pub(super) fn new(
    qualified1: &Geom2dGccQCurve,
    point2: &GpPnt2d,
    radius: f64,
    tolerance: f64,
) -> Geom2dGccCirc2d2TanRadGeo {
    let mut s = Geom2dGccCirc2d2TanRadGeo::empty();
    let tol = tolerance.abs();
    let cat1 = super::Categories::of_qcurve(qualified1);
    check_qualifier(cat1.all_of(), true, WHO);
    let cu1 = qualified1.qualified();
    check_radius(radius, WHO);
    let (nbrcote1, cote1) = cotes(qualified1, radius);
    // `gp_Circ2d Circ(gp_Ax2d(Point2, gp_Dir2d(gp_Dir2d::D::X)), Radius);`
    // (`cxx:559`), and the closed equivalent domain of `cxx:560-564`.
    let circ = circle_at(*point2, radius);
    let d1 = circle_domain(&circ, tol);
    let mut inter = IntConicCurveGen::new();
    for jcote1 in 1..=nbrcote1 {
        if s.nbr_sol >= MAX_SOL {
            break;
        }
        // `Adaptor2d_OffsetCurve Cu2(HCu1, -cote1(jcote1));` (`cxx:571-572`).
        let cu2 = offset_of(cu1, -cote1[jcote1 - 1]);
        let d2 = curve_domain(cu2.as_ref(), tol);
        inter.perform_circle(&circ, &d1, cu2.as_ref(), &d2, tol, tol);
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
                    s.qualifier1[n - 1] = qualified1.qualifier();
                    s.qualifier2[n - 1] = GccEntPosition::NoQualifier;
                    s.the_same1[n - 1] = 0;
                    s.the_same2[n - 1] = 0;
                    s.pararg1[n - 1] = pt.param_on_second();
                    s.pararg2[n - 1] = 0.0;
                    s.pnttg1sol[n - 1] = curve_tool::value(cu1, s.pararg1[n - 1])
                        .expect("Geom2dAdaptor_Curve::Value on the curve operand");
                    s.pnttg2sol[n - 1] = *point2;
                    s.par1sol[n - 1] = parameter_circ2d(&s.cirsol[n - 1], &s.pnttg1sol[n - 1]);
                    s.par2sol[n - 1] = parameter_circ2d(&s.cirsol[n - 1], &s.pnttg2sol[n - 1]);
                }
            }
            s.well_done = true;
        }
    }
    s
}
