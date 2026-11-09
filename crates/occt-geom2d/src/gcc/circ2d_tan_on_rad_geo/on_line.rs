//! `Geom2dGcc_Circ2dTanOnRadGeo(Geom2dGcc_QCurve, gp_Lin2d, Radius, Tolerance)`
//! (`Geom2dGcc_Circ2dTanOnRadGeo.cxx:61-167`): the operand is offset by
//! `-Coef(j)` and that offset is intersected with the carrier line; every
//! intersection point is a centre.
//!
//! Translation notes:
//! - `par1sol` is computed from the previous `pnttg1sol` slot (`cxx:159-160`),
//! - `parcen3` is `ParamOnFirst` on the carrier line and `pararg1` is
//!   `ParamOnSecond` on the offset curve (`cxx:158-159`).

use occt_core::elib::clib2d::parameter_circ2d;
use occt_core::gp::GpLin2d;
use occt_core::intres2d::IntRes2dDomain;

use crate::gcc::curve_tool;
use crate::gcc::qualified_curve::Geom2dGccQCurve;
use crate::geom2d_int::IntConicCurveGen;

use super::{circle_at, cotes, curve_domain, offset_of, Geom2dGccCirc2dTanOnRadGeo, MAX_SOL};

/// `Geom2dGcc_Circ2dTanOnRadGeo(Geom2dGcc_QCurve, gp_Lin2d, double, double)`
/// (`cxx:61-167`).
pub(super) fn new(
    qualified1: &Geom2dGccQCurve,
    on_line: &GpLin2d,
    radius: f64,
    tolerance: f64,
) -> Geom2dGccCirc2dTanOnRadGeo {
    let mut s = Geom2dGccCirc2dTanOnRadGeo::empty();
    let tol = tolerance.abs();
    if !(qualified1.is_enclosed()
        || qualified1.is_enclosing()
        || qualified1.is_outside()
        || qualified1.is_unqualified())
    {
        panic!("GccEnt_BadQualifier in Geom2dGcc_Circ2dTanOnRadGeo");
    }
    let cu1 = qualified1.qualified();
    if radius < 0.0 {
        panic!("Standard_NegativeValue in Geom2dGcc_Circ2dTanOnRadGeo");
    }
    let (nbrcote1, coef) = cotes(
        qualified1.is_enclosed(),
        qualified1.is_outside(),
        qualified1.is_unqualified(),
        radius,
    );
    // `IntRes2d_Domain D1;` (`cxx:129`): the empty domain.
    let d1 = IntRes2dDomain::default();
    let mut inter = IntConicCurveGen::new();
    for jcote1 in 1..=nbrcote1 {
        // `Adaptor2d_OffsetCurve C2(HCu1, -Coef(jcote1));` (`cxx:136`).
        let c2 = offset_of(cu1, -coef[jcote1 - 1]);
        let d2 = curve_domain(c2.as_ref(), tol);
        inter.perform_line(on_line, &d1, c2.as_ref(), &d2, tol, tol);
        if inter.is_done() {
            if !inter.result().is_empty() {
                for i in 1..=inter.result().nb_points() {
                    if s.nbr_sol >= MAX_SOL {
                        break;
                    }
                    s.nbr_sol += 1;
                    let n = s.nbr_sol;
                    let point = inter.result().point(i);
                    let center = *point.value();
                    let par_on_second = point.param_on_second();
                    let par_on_first = point.param_on_first();
                    s.cirsol[n - 1] = circle_at(center, radius);
                    s.qualifier1[n - 1] = qualified1.qualifier();
                    s.the_same1[n - 1] = 0;
                    s.pararg1[n - 1] = par_on_second;
                    s.parcen3[n - 1] = par_on_first;
                    s.par1sol[n - 1] = parameter_circ2d(&s.cirsol[n - 1], &s.pnttg1sol[n - 1]);
                    s.pnttg1sol[n - 1] = curve_tool::value(cu1, s.pararg1[n - 1])
                        .expect("Geom2dAdaptor_Curve::Value on the operand");
                    s.pntcen3[n - 1] = center;
                }
            }
            s.well_done = true;
        }
    }
    s
}
