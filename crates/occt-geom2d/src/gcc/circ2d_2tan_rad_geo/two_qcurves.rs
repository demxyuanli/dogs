//! `Geom2dGcc_Circ2d2TanRadGeo(Geom2dGcc_QCurve, Geom2dGcc_QCurve, Radius,
//! Tolerance)` (`Geom2dGcc_Circ2d2TanRadGeo.cxx:785-1006`) and the file-static
//! precision helper `PrecRoot` (`cxx:618-770`).
//!
//! Both operands are offset by `-cote(j)` and intersected by `Geom2dInt_GInter`,
//! which takes its generic `IntCurve_IntPolyPolyGen` branch. When the two
//! offsets meet in a tangent zone the 2D intersection point is refined by a
//! Newton iteration on the two "gap vector orthogonal to each tangent" equations.
//!
//! Translation notes:
//! - the outer sweep is unbounded and the inner one carries
//!   `NbrSol < aNbSolMAX` (`cxx:919`, `:935`); the point loop is bounded again
//!   (`cxx:953`),
//! - `PrecRoot(C1, C2, aU0, aV0, aU0, aV0)` (`cxx:980`) passes the coordinates
//!   both in and out, so the refinement overwrites `aU0` / `aV0` in place; the
//!   centre uses the refined value (`cxx:984`) while `pararg1` / `pnttg1sol`
//!   still use the unrefined intersection parameters (`cxx:991-994`),
//! - OCCT defect kept: `qualifier1(NbrSol) = Qualified1.Qualifier();` appears
//!   twice (`cxx:987-988`) and `qualifier2` is never written, so the second
//!   qualifier of every solution keeps its initial slot,
//! - `aTheta = aD1u * aD1v` (`cxx:721`) is the 2D cross product, `GpVec2d`'s
//!   `operator*`, not a dot product.

use occt_core::gp::GpVec2d;
use occt_core::precision::{epsilon, APPROXIMATION, PAPPROXIMATION, REAL_SMALL};

use crate::curve::Curve2d;
use crate::gcc::curve_tool;
use crate::gcc::qualified_curve::Geom2dGccQCurve;
use crate::gcc::tan_geo_common::{circle_at, offset_of};
use crate::geom2d_int::Geom2dIntGInter;

use super::{check_qualifier, check_radius, cotes, Categories, Geom2dGccCirc2d2TanRadGeo, MAX_SOL};

const WHO: &str = "Geom2dGcc_Circ2d2TanRadGeo";

/// `IsEqual(theValue1, theValue2)` (`Standard_Real.hxx:107-110`):
/// `Abs(a - b) < RealSmall()`.
fn is_equal(a: f64, b: f64) -> bool {
    (a - b).abs() < REAL_SMALL
}

/// `PrecRoot` (`cxx:618-770`): returns `theUfinal` / `theVfinal`.
///
/// `theC1` / `theC2` are the two `Adaptor2d_OffsetCurve`s. When the tangents at
/// `(theU0, theV0)` are not nearly parallel
/// (`|cross|^2 > 1e-6 |d1u|^2 |d1v|^2`, `cxx:678-683`) the initial point is
/// returned unchanged, because the generic intersector is already accurate
/// there.
fn prec_root(
    the_c1: &dyn Curve2d,
    the_c2: &dyn Curve2d,
    the_u0: f64,
    the_v0: f64,
) -> (f64, f64) {
    let mut the_u_final = the_u0;
    let mut the_v_final = the_v0;
    /// `aNbIterMax = 100` (`cxx:662`).
    const NB_ITER_MAX: i32 = 100;
    let mut a_u = the_u0;
    let mut a_v = the_v0;
    let (mut a_pu, mut a_d1u, mut a_d2u) = the_c1.d2(a_u);
    let (mut a_pv, mut a_d1v, mut a_d2v) = the_c2.d2(a_v);
    let a_cr_prod = a_d1u.crossed(&a_d1v).abs();
    if a_cr_prod * a_cr_prod > 1.0e-6 * a_d1u.square_magnitude() * a_d1v.square_magnitude() {
        return (the_u_final, the_v_final);
    }
    let mut a_step_u = 0.0;
    let mut a_step_v = 0.0;
    // `RealFirst()` (`Standard_Real.hxx:126-129`).
    let mut a_sq_dist_prev = f64::MIN_POSITIVE;
    let mut a_nb_iter = 0;
    loop {
        a_nb_iter += 1;
        // `gp_Vec2d aVuv(aPv, aPu);` (`cxx:690`).
        let a_vuv = GpVec2d::new(a_pu.x() - a_pv.x(), a_pu.y() - a_pv.y());
        let a_sq_dist = a_vuv.square_magnitude();
        if is_equal(a_sq_dist, 0.0) {
            break;
        }
        if a_nb_iter == 1 || a_sq_dist < a_sq_dist_prev {
            a_sq_dist_prev = a_sq_dist;
            the_u_final = a_u;
            the_v_final = a_v;
        }
        let a_g1 = a_d1u.magnitude();
        let a_g2 = a_d1v.magnitude();
        if is_equal(a_g1, 0.0) || is_equal(a_g2, 0.0) {
            break;
        }
        let a_f1 = a_vuv.dot(&a_d1u);
        let a_f2 = a_vuv.dot(&a_d1v);
        let a_fi_u = a_vuv.dot(&a_d2u);
        let a_fi_v = a_vuv.dot(&a_d2v);
        let a_psi_u = a_d1u.dot(&a_d2u);
        let a_psi_v = a_d1v.dot(&a_d2v);
        // `double aTheta = aD1u * aD1v;` (`cxx:721`): 2D cross product.
        let a_theta = a_d1u.crossed(&a_d1v);
        let a_s1 = a_f1 / a_g1;
        let a_s2 = a_f2 / a_g2;
        let a_ds1u = (a_g1 * a_g1 + a_fi_u) / a_g1 - (a_s1 * a_psi_u / (a_g1 * a_g1));
        let a_ds1v = -a_theta / a_g1;
        let a_ds2u = a_theta / a_g2;
        let a_ds2v = (a_fi_v - a_g2 * a_g2) / a_g2 - (a_s2 * a_psi_v / (a_g2 * a_g2));
        let a_det = a_ds1u * a_ds2v - a_ds1v * a_ds2u;
        if is_equal(a_det, 0.0) {
            if !is_equal(a_step_v, 0.0) && !is_equal(a_ds1u, 0.0) {
                a_v += a_step_v;
                a_u -= (a_ds1v * a_step_v - a_s1) / a_ds1u;
            } else if !is_equal(a_step_u, 0.0) && !is_equal(a_ds1v, 0.0) {
                a_u += a_step_u;
                a_v -= (a_ds1u * a_step_u - a_s1) / a_ds1v;
            } else {
                break;
            }
        } else {
            a_step_u = -(a_s1 * a_ds2v - a_s2 * a_ds1v) / a_det;
            a_step_v = -(a_s2 * a_ds1u - a_s1 * a_ds2u) / a_det;
            if a_step_u.abs() < epsilon(a_u.abs()) {
                if a_step_v.abs() < epsilon(a_v.abs()) {
                    break;
                }
            }
            a_u += a_step_u;
            a_v += a_step_v;
        }
        let (pu, d1u, d2u) = the_c1.d2(a_u);
        a_pu = pu;
        a_d1u = d1u;
        a_d2u = d2u;
        let (pv, d1v, d2v) = the_c2.d2(a_v);
        a_pv = pv;
        a_d1v = d1v;
        a_d2v = d2v;
        if a_nb_iter > NB_ITER_MAX {
            break;
        }
    }
    (the_u_final, the_v_final)
}

/// `Geom2dGcc_Circ2d2TanRadGeo(Geom2dGcc_QCurve, Geom2dGcc_QCurve, double,
/// double)` (`cxx:785-1006`).
pub(super) fn new(
    qualified1: &Geom2dGccQCurve,
    qualified2: &Geom2dGccQCurve,
    radius: f64,
    tolerance: f64,
) -> Geom2dGccCirc2d2TanRadGeo {
    let mut s = Geom2dGccCirc2d2TanRadGeo::empty();
    let tol = tolerance.abs();
    let cat1 = Categories::of_qcurve(qualified1);
    let cat2 = Categories::of_qcurve(qualified2);
    check_qualifier(cat1.all_of(), cat2.all_of(), WHO);
    let cu1 = qualified1.qualified();
    let cu2 = qualified2.qualified();
    check_radius(radius, WHO);
    let (nbrcote1, cote1, nbrcote2, cote2) = cotes(cat1, cat2, radius);
    let a_sq_approx_tol = APPROXIMATION * APPROXIMATION;
    let mut inter = Geom2dIntGInter::new();
    for jcote1 in 1..=nbrcote1 {
        // `Adaptor2d_OffsetCurve C1(HCu1, -cote1(jcote1));` (`cxx:921-923`).
        let c1 = offset_of(cu1, -cote1[jcote1 - 1]);
        for jcote2 in 1..=nbrcote2 {
            if s.nbr_sol >= MAX_SOL {
                break;
            }
            // `Adaptor2d_OffsetCurve C2(HCu2, -cote2(jcote2));` (`cxx:936-938`).
            let c2 = offset_of(cu2, -cote2[jcote2 - 1]);
            inter.perform_curves(c1.as_ref(), c2.as_ref(), tol, tol);
            if inter.is_done() {
                if !inter.result().is_empty() {
                    for i in 1..=inter.result().nb_points() {
                        if s.nbr_sol >= MAX_SOL {
                            break;
                        }
                        let pt = inter.result().point(i);
                        let mut a_u0 = pt.param_on_first();
                        let mut a_v0 = pt.param_on_second();
                        let a_u1 = a_u0 - PAPPROXIMATION;
                        let a_v1 = a_v0 - PAPPROXIMATION;
                        let a_u2 = a_u0 + PAPPROXIMATION;
                        let a_v2 = a_v0 + PAPPROXIMATION;
                        let p11 = c1.d0(a_u1);
                        let p12 = c2.d0(a_v1);
                        let p21 = c1.d0(a_u2);
                        let p22 = c2.d0(a_v2);
                        let a_dist1112 = p11.square_distance(&p12);
                        let a_dist1122 = p11.square_distance(&p22);
                        let a_dist1221 = p12.square_distance(&p21);
                        let a_dist2122 = p21.square_distance(&p22);
                        if a_dist1112.min(a_dist1122) <= a_sq_approx_tol
                            && a_dist1221.min(a_dist2122) <= a_sq_approx_tol
                        {
                            let (u, v) = prec_root(c1.as_ref(), c2.as_ref(), a_u0, a_v0);
                            a_u0 = u;
                            a_v0 = v;
                            // `aV0` is refined in place (`cxx:980`) but only
                            // `aU0` is read afterwards (`cxx:984`).
                            let _ = a_v0;
                        }
                        s.nbr_sol += 1;
                        let n = s.nbr_sol;
                        // `gp_Pnt2d Center(C1.Value(aU0));` (`cxx:984`).
                        let center = c1.d0(a_u0);
                        s.cirsol[n - 1] = circle_at(center, radius);
                        // `cxx:987-988` writes `qualifier1` twice; `qualifier2`
                        // is left untouched.
                        s.qualifier1[n - 1] = qualified1.qualifier();
                        s.qualifier1[n - 1] = qualified1.qualifier();
                        s.the_same1[n - 1] = 0;
                        s.the_same2[n - 1] = 0;
                        s.pararg1[n - 1] = pt.param_on_first();
                        s.pararg2[n - 1] = pt.param_on_second();
                        s.pnttg1sol[n - 1] = curve_tool::value(cu1, s.pararg1[n - 1])
                            .expect("Geom2dAdaptor_Curve::Value on the first operand");
                        s.pnttg2sol[n - 1] = curve_tool::value(cu2, s.pararg2[n - 1])
                            .expect("Geom2dAdaptor_Curve::Value on the second operand");
                        s.par1sol[n - 1] = occt_core::elib::clib2d::parameter_circ2d(
                            &s.cirsol[n - 1],
                            &s.pnttg1sol[n - 1],
                        );
                        s.par2sol[n - 1] = occt_core::elib::clib2d::parameter_circ2d(
                            &s.cirsol[n - 1],
                            &s.pnttg2sol[n - 1],
                        );
                    }
                }
                s.well_done = true;
            }
        }
    }
    s
}
