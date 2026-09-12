//! Bounding box of a 2D line.
//!
//! Source: `GeomBndLib_Line2d.hxx` (`Box(gp_Lin2d, U1, U2, Tol)`). Infinite
//! parameters open the box via [`crate::geom_bnd_lib_inf2d`]; finite
//! endpoints are added with `ElCLib::Value`. Analytic `BoxOptimal` equals
//! `Box`.

use occt_core::bnd::BndBox2d;
use occt_core::gp::GpLin2d;
use occt_core::precision::Precision;

use crate::geom_bnd_lib_elclib2d::line_value;
use crate::geom_bnd_lib_inf2d::{open_max, open_min, open_min_max};

/// `GeomBndLib_Line2d::Box(theLin, theU1, theU2, theTol)`.
pub fn box_lin(the_lin: &GpLin2d, the_u1: f64, the_u2: f64, the_tol: f64) -> Result<BndBox2d, String> {
    let mut a_box = BndBox2d::new();
    if Precision::is_negative_infinite(the_u1) {
        if Precision::is_negative_infinite(the_u2) {
            return Err("GeomBndLib_Line2d::Box - bad parameter".into());
        } else if Precision::is_positive_infinite(the_u2) {
            open_min_max(the_lin.direction(), &mut a_box);
            a_box.add_point(&line_value(0.0, the_lin));
        } else {
            open_min(the_lin.direction(), &mut a_box);
            a_box.add_point(&line_value(the_u2, the_lin));
        }
    } else if Precision::is_positive_infinite(the_u1) {
        if Precision::is_negative_infinite(the_u2) {
            open_min_max(the_lin.direction(), &mut a_box);
            a_box.add_point(&line_value(0.0, the_lin));
        } else if Precision::is_positive_infinite(the_u2) {
            return Err("GeomBndLib_Line2d::Box - bad parameter".into());
        } else {
            open_max(the_lin.direction(), &mut a_box);
            a_box.add_point(&line_value(the_u2, the_lin));
        }
    } else {
        a_box.add_point(&line_value(the_u1, the_lin));
        if Precision::is_negative_infinite(the_u2) {
            open_min(the_lin.direction(), &mut a_box);
        } else if Precision::is_positive_infinite(the_u2) {
            open_max(the_lin.direction(), &mut a_box);
        } else {
            a_box.add_point(&line_value(the_u2, the_lin));
        }
    }
    a_box.enlarge(the_tol);
    Ok(a_box)
}

/// `GeomBndLib_Line2d::BoxOptimal` — identical to [`box_lin`].
pub fn box_lin_optimal(
    the_lin: &GpLin2d,
    the_u1: f64,
    the_u2: f64,
    the_tol: f64,
) -> Result<BndBox2d, String> {
    box_lin(the_lin, the_u1, the_u2, the_tol)
}
