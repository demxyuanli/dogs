//! `Geom2dGcc_Circ2dTanOnRadGeo(GccEnt_QualifiedCirc, Adaptor, Radius, Tolerance)`
//! (`Geom2dGcc_Circ2dTanOnRadGeo.cxx:317-445`): the operand circle is
//! concentric-enlarged by `cote1(j)` and intersected with the carrier.
//!
//! Translation notes:
//! - `gp_Circ2d Circ(C1.XAxis(), R1 + cote1(jcote1))` (`cxx:391`) reuses the
//!   operand's X axis, so the parameters of `Circ` and `C1` agree,
//! - `pararg1` is `ParamOnFirst` on `Circ` and the tangency point is evaluated
//!   on `C1` with the same parameter (`cxx:436-437`),
//! - `par1sol` is computed before `pnttg1sol` is written (`cxx:436-437`).

use occt_core::elib::clib2d::parameter_circ2d;
use occt_core::gp::{GpAx22d, GpCirc2d};
use occt_core::intres2d::IntRes2dDomain;

use crate::adaptor_curve::Geom2dAdaptorCurve;
use crate::gcc::curve_tool;
use crate::gcc::gcc_ent::{GccEntPosition, GccEntQualifiedCirc};
use crate::geom2d_int::IntConicCurveGen;

use super::{
    circle_at, circle_domain, cotes, value_on_circ, Geom2dGccCirc2dTanOnRadGeo, MAX_SOL, THEFIRST,
    THELAST,
};

/// `Geom2dGcc_Circ2dTanOnRadGeo(GccEnt_QualifiedCirc, Geom2dAdaptor_Curve,
/// double, double)` (`cxx:317-445`).
pub(super) fn new(
    qualified1: &GccEntQualifiedCirc,
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
    let c1 = qualified1.qualified();
    let center1 = c1.location();
    let r1 = c1.radius();
    if radius < 0.0 {
        panic!("Standard_NegativeValue in Geom2dGcc_Circ2dTanOnRadGeo");
    }
    let (nbrcote1, cote) = cotes(
        qualified1.is_enclosed(),
        qualified1.is_outside(),
        qualified1.is_unqualified(),
        radius,
    );
    let carrier = on_curve
        .curve()
        .expect("Geom2dAdaptor_Curve::Curve on the carrier")
        .clone();
    let mut inter = IntConicCurveGen::new();
    for jcote1 in 1..=nbrcote1 {
        // `gp_Circ2d Circ(C1.XAxis(), R1 + cote1(jcote1));` (`cxx:391`).
        let x_axis = c1.x_axis();
        let circ = GpCirc2d::new(
            GpAx22d::from_xdir(*x_axis.location(), *x_axis.direction()),
            r1 + cote[jcote1 - 1],
        );
        let d1 = circle_domain(&circ, tol);
        let first = curve_tool::first_parameter(on_curve).max(THEFIRST);
        let last = curve_tool::last_parameter(on_curve).min(THELAST);
        let p_first = curve_tool::value(on_curve, first)
            .expect("Geom2dAdaptor_Curve::Value on the carrier");
        let p_last = curve_tool::value(on_curve, last)
            .expect("Geom2dAdaptor_Curve::Value on the carrier");
        let d2 = IntRes2dDomain::bounded(&p_first, first, tol, &p_last, last, tol);
        inter.perform_circle(&circ, &d1, carrier.as_ref(), &d2, tol, tol);
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
                    let distcc1 = center.distance(&center1);
                    s.qualifier1[n - 1] = if !qualified1.is_unqualified() {
                        qualified1.qualifier()
                    } else if (distcc1 + radius - r1).abs() < tol {
                        GccEntPosition::Enclosed
                    } else if (distcc1 - r1 - radius).abs() < tol {
                        GccEntPosition::Outside
                    } else {
                        GccEntPosition::Enclosing
                    };
                    s.the_same1[n - 1] = 0;
                    s.pararg1[n - 1] = par_on_first;
                    s.parcen3[n - 1] = par_on_second;
                    s.par1sol[n - 1] = parameter_circ2d(&s.cirsol[n - 1], &s.pnttg1sol[n - 1]);
                    s.pnttg1sol[n - 1] = value_on_circ(&c1, s.pararg1[n - 1]);
                    s.pntcen3[n - 1] = center;
                }
            }
            s.well_done = true;
        }
    }
    s
}
