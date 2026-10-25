//! `ChFi3d_StripeEdgeInter` - intersect the p-curves of two stripes that share
//! a face, and fail when they meet.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_0.cxx:3346-3443`
//! (OCCT 8.0.0), declaration `ChFi3d_Builder_0.hxx`.
//!
//! `ChFi3d_Builder::Compute` calls it for every unordered pair of contours
//! (`ChFi3d_Builder.cxx:355-389`); the shared-corner-point test below is what
//! excludes the pairs that touch at a vertex. An intersection between the two
//! stripes' face interferences means the fillets overlap, and OCCT throws
//! `StdFail_NotDone("StripeEdgeInter : fillets have too big radiuses")`;
//! the port answers that as the `Err` string.
//!
//! Divergences that are not heuristics:
//! - The `TopOpeBRepDS_DataStructure&` parameter is commented out of the OCCT
//!   body (`cxx:3355`), so the port drops it: passing the builder's `myDS` as
//!   `&mut` would only borrow-conflict with the stripe list the caller walks.
//! - `Geom2dAdaptor_Curve(C, First, Last)` becomes
//!   [`Geom2dIntGInter::domain_over_range`], which builds the same bounded
//!   domain `ComputeDomain` would derive from the narrowed adaptor.

use occt_core::precision::{PCONFUSION, REAL_SMALL};
use occt_geom2d::curve::Curve2d;
use occt_geom2d::geom2d_int::Geom2dIntGInter;

use crate::chfi_ds::stripe::Stripe;
use crate::chfi_ds::surf_data::SurfData;

/// `IsEqual(A, B)` of `Standard_Real.hxx:148-151`.
fn is_equal_real(a: f64, b: f64) -> bool {
    (a - b).abs() < REAL_SMALL
}

/// `ChFi3d_StripeEdgeInter(S1, S2, DStr, Tol2d)` (`cxx:3346-3443`).
pub fn ch_fi3d_stripe_edge_inter(
    the_stripe1: &Stripe,
    the_stripe2: &Stripe,
    tol2d: f64,
) -> Result<(), String> {
    // `cxx:3358-3371`: stripes that share a corner point are not checked.
    for i_sur1 in 1..=2 {
        for i_sur2 in 1..=2 {
            if the_stripe1.index_point(false, i_sur1) == the_stripe2.index_point(false, i_sur2)
                || the_stripe1.index_point(false, i_sur1) == the_stripe2.index_point(true, i_sur2)
                || the_stripe1.index_point(true, i_sur1) == the_stripe2.index_point(false, i_sur2)
                || the_stripe1.index_point(true, i_sur1) == the_stripe2.index_point(true, i_sur2)
            {
                return Ok(());
            }
        }
    }

    let sur_dat1: &[std::sync::Arc<SurfData>] = the_stripe1.set_of_surf_data();
    let sur_dat2: &[std::sync::Arc<SurfData>] = the_stripe2.set_of_surf_data();

    let tol_domain = tol2d.max(PCONFUSION);
    let mut an_intersector = Geom2dIntGInter::new();

    // `cxx:3380-3441`: every part of the first stripe against every part of
    // the second one.
    for a_dat1 in sur_dat1.iter() {
        let ishape11 = a_dat1.index_of_s1();
        let ishape12 = a_dat1.index_of_s2();
        for a_dat2 in sur_dat2.iter() {
            let ishape21 = a_dat2.index_of_s1();
            let ishape22 = a_dat2.index_of_s2();

            // `cxx:3395-3417`: pick the two face interferences that lie on a
            // common face.
            let (afi1, afi2) = if ishape11 == ishape21 {
                (a_dat1.interference_on_s1(), a_dat2.interference_on_s1())
            } else if ishape11 == ishape22 {
                (a_dat1.interference_on_s1(), a_dat2.interference_on_s2())
            } else if ishape12 == ishape21 {
                (a_dat1.interference_on_s2(), a_dat2.interference_on_s1())
            } else if ishape12 == ishape22 {
                (a_dat1.interference_on_s2(), a_dat2.interference_on_s2())
            } else {
                // No common faces.
                continue;
            };

            // `cxx:3420-3426`: degenerates carry a single-parameter p-curve or
            // no p-curve at all; they are skipped.
            if is_equal_real(afi1.first_parameter(), afi1.last_parameter())
                || is_equal_real(afi2.first_parameter(), afi2.last_parameter())
            {
                continue;
            }
            let (pc2d1, pc2d2) = match (afi1.p_curve_on_face(), afi2.p_curve_on_face()) {
                (Some(p1), Some(p2)) => (p1, p2),
                _ => continue,
            };

            // `cxx:3428-3438`.
            let c1: &dyn Curve2d = pc2d1.as_ref();
            let c2: &dyn Curve2d = pc2d2.as_ref();
            let d1 = Geom2dIntGInter::domain_over_range(
                c1,
                afi1.first_parameter(),
                afi1.last_parameter(),
                tol_domain,
            );
            let d2 = Geom2dIntGInter::domain_over_range(
                c2,
                afi2.first_parameter(),
                afi2.last_parameter(),
                tol_domain,
            );
            an_intersector.perform(c1, &d1, c2, &d2, tol2d, PCONFUSION);
            if an_intersector.nb_segments() > 0 || an_intersector.nb_points() > 0 {
                return Err("StripeEdgeInter : fillets have too big radiuses".to_string());
            }
        }
    }
    Ok(())
}
