//! `Geom2dGcc_Circ2dTanOnRadGeo(Geom2dGcc_QCurve, Adaptor, Radius, Tolerance)`
//! (`Geom2dGcc_Circ2dTanOnRadGeo.cxx:604-717`): the operand offset by
//! `-cote1(j)` is intersected with the carrier by `Geom2dInt_GInter`.
//!
//! Translation notes:
//! - both sides are `Adaptor2d_OffsetCurve`, so both report `GeomAbs_OtherCurve`
//!   and `Geom2dInt_GInter` takes its generic `IntCurve_IntPolyPolyGen` branch;
//!   the port therefore wraps the carrier in a zero offset too (`cxx:683`) and
//!   does not hand it over as its native type,
//! - the carrier offset is built with the single-argument constructor, whose
//!   `myOffset` is `0.0` (`Adaptor2d_OffsetCurve.hxx:42-44`),
//! - `par1sol` is computed before `pnttg1sol` is written (`cxx:709-710`).

use occt_core::elib::clib2d::parameter_circ2d;

use crate::adaptor_curve::Geom2dAdaptorCurve;
use crate::gcc::curve_tool;
use crate::gcc::qualified_curve::Geom2dGccQCurve;
use crate::geom2d_int::Geom2dIntGInter;

use super::{circle_at, cotes, curve_domain, offset_of, Geom2dGccCirc2dTanOnRadGeo, MAX_SOL};

/// `Geom2dGcc_Circ2dTanOnRadGeo(Geom2dGcc_QCurve, Geom2dAdaptor_Curve, double,
/// double)` (`cxx:604-717`).
pub(super) fn new(
    qualified1: &Geom2dGccQCurve,
    on_curve: &Geom2dAdaptorCurve,
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
    let (nbrcote1, cote) = cotes(
        qualified1.is_enclosed(),
        qualified1.is_outside(),
        qualified1.is_unqualified(),
        radius,
    );
    // `Adaptor2d_OffsetCurve C2(HOnCurv);` (`cxx:683`): offset 0.
    let c2 = offset_of(on_curve, 0.0);
    let d2 = curve_domain(c2.as_ref(), tol);
    let mut inter = Geom2dIntGInter::new();
    for jcote1 in 1..=nbrcote1 {
        // `Adaptor2d_OffsetCurve C1(HCu1, -cote1(jcote1));` (`cxx:676-677`).
        let c1 = offset_of(cu1, -cote[jcote1 - 1]);
        let d1 = curve_domain(c1.as_ref(), tol);
        inter.perform(c1.as_ref(), &d1, c2.as_ref(), &d2, tol, tol);
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
                    let par_on_first = point.param_on_first();
                    let par_on_second = point.param_on_second();
                    s.cirsol[n - 1] = circle_at(center, radius);
                    s.qualifier1[n - 1] = qualified1.qualifier();
                    s.the_same1[n - 1] = 0;
                    s.pararg1[n - 1] = par_on_first;
                    s.parcen3[n - 1] = par_on_second;
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
