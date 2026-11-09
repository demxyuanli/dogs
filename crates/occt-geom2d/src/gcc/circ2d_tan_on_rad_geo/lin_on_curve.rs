//! `Geom2dGcc_Circ2dTanOnRadGeo(GccEnt_QualifiedLin, Adaptor, Radius, Tolerance)`
//! (`Geom2dGcc_Circ2dTanOnRadGeo.cxx:466-583`): the operand line is displaced
//! by `cote1(j)` along its left normal and intersected with the carrier.
//!
//! Translation notes:
//! - OCCT builds the displaced line from `gp_Pnt2d(dir1.XY() + cote1(j) *
//!   norm1.XY())` (`cxx:537`), i.e. the **direction** is used where the
//!   operand's origin would be expected; the port copies that verbatim because
//!   `origin1` only feeds the qualifier test (`cxx:559-570`),
//! - the qualifier uses `gp_Dir2d dc1(origin1.XY() - Center.XY())`
//!   (`cxx:559`), which throws `Standard_ConstructionError` on a null vector
//!   and is ported as a panic,
//! - `par1sol` is computed before `pnttg1sol` is written (`cxx:575-576`).

use occt_core::elib::clib2d::parameter_circ2d;
use occt_core::gp::{GpAx2d, GpDir2d, GpLin2d, GpPnt2d};
use occt_core::intres2d::IntRes2dDomain;

use crate::adaptor_curve::Geom2dAdaptorCurve;
use crate::gcc::curve_tool;
use crate::gcc::gcc_ent::{GccEntPosition, GccEntQualifiedLin};
use crate::geom2d_int::IntConicCurveGen;

use super::{circle_at, cotes, value_on_line, Geom2dGccCirc2dTanOnRadGeo, MAX_SOL, THEFIRST, THELAST};

/// `Geom2dGcc_Circ2dTanOnRadGeo(GccEnt_QualifiedLin, Geom2dAdaptor_Curve,
/// double, double)` (`cxx:466-583`).
pub(super) fn new(
    qualified1: &GccEntQualifiedLin,
    on_curve: &Geom2dAdaptorCurve,
    radius: f64,
    tolerance: f64,
) -> Geom2dGccCirc2dTanOnRadGeo {
    let mut s = Geom2dGccCirc2dTanOnRadGeo::empty();
    let tol = tolerance.abs();
    // `IsEnclosing` is rejected here, unlike the other constructors
    // (`cxx:498-501`).
    if !(qualified1.is_enclosed() || qualified1.is_outside() || qualified1.is_unqualified()) {
        panic!("GccEnt_BadQualifier in Geom2dGcc_Circ2dTanOnRadGeo");
    }
    let l1 = qualified1.qualified();
    let origin1 = l1.location();
    let dir1 = *l1.direction();
    let norm1 = GpDir2d::new(-dir1.y(), dir1.x()).expect("normal of a unit direction");
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
        // `gp_Pnt2d Point(dir1.XY() + cote1(jcote1) * norm1.XY());` (`cxx:537`).
        let point = GpPnt2d::new(
            dir1.x() + cote[jcote1 - 1] * norm1.x(),
            dir1.y() + cote[jcote1 - 1] * norm1.y(),
        );
        // `gp_Lin2d Line(Point, dir1);` (`cxx:538`).
        let line = GpLin2d::new(GpAx2d::new(point, dir1));
        let d1 = IntRes2dDomain::default();
        let first = curve_tool::first_parameter(on_curve).max(THEFIRST);
        let last = curve_tool::last_parameter(on_curve).min(THELAST);
        let p_first = curve_tool::value(on_curve, first)
            .expect("Geom2dAdaptor_Curve::Value on the carrier");
        let p_last = curve_tool::value(on_curve, last)
            .expect("Geom2dAdaptor_Curve::Value on the carrier");
        let d2 = IntRes2dDomain::bounded(&p_first, first, tol, &p_last, last, tol);
        inter.perform_line(&line, &d1, carrier.as_ref(), &d2, tol, tol);
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
                    // `gp_Dir2d dc1(origin1.XY() - Center.XY());` (`cxx:559`).
                    let dc1 = GpDir2d::new(
                        origin1.x() - center.x(),
                        origin1.y() - center.y(),
                    )
                    .expect("Standard_ConstructionError: gp_Dir2d from the origin offset");
                    s.qualifier1[n - 1] = if !qualified1.is_unqualified() {
                        qualified1.qualifier()
                    } else if dc1.dot(&norm1) > 0.0 {
                        GccEntPosition::Outside
                    } else {
                        GccEntPosition::Enclosed
                    };
                    s.the_same1[n - 1] = 0;
                    s.pararg1[n - 1] = par_on_first;
                    s.parcen3[n - 1] = par_on_second;
                    s.par1sol[n - 1] = parameter_circ2d(&s.cirsol[n - 1], &s.pnttg1sol[n - 1]);
                    s.pnttg1sol[n - 1] = value_on_line(&l1, s.pararg1[n - 1]);
                    s.pntcen3[n - 1] = center;
                }
            }
            s.well_done = true;
        }
    }
    s
}
