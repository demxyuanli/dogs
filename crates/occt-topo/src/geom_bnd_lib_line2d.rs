//! Bounding box of a 2D line.
//!
//! Source: `GeomBndLib_Line2d.hxx` (`Box(gp_Lin2d, U1, U2, Tol)`). Infinite
//! parameters open the box via [`crate::geom_bnd_lib_inf2d`]; finite
//! endpoints are added with `ElCLib::Value`. Analytic `BoxOptimal` equals
//! `Box`.
//! T-97: items below are faithful ports of the named OCCT source, but their
//! OCCT-side consumers are not all ported yet, so parts are not called from this
//! crate. The `dead_code` allowance is deliberate: **pending wiring**, not dead
//! code. Do not delete them to silence warnings (see
//! specs/_a3n00_gap_analysis.md §9.309/§9.310); wire the consumer instead.
#![allow(dead_code)]

use occt_core::bnd::BndBox2d;
use occt_core::gp::GpLin2d;
use occt_core::precision::Precision;

use crate::geom_bnd_lib_elclib2d::line_value_2d;
use crate::geom_bnd_lib_inf2d::{open_max_2d, open_min_2d, open_min_max_2d};

/// `GeomBndLib_Line2d::Box(theLin, theU1, theU2, theTol)`.
pub fn box_lin_2d(lin: &GpLin2d, u1: f64, u2: f64, tol: f64) -> Result<BndBox2d, String> {
    let mut bnd_box = BndBox2d::new();
    if Precision::is_negative_infinite(u1) {
        if Precision::is_negative_infinite(u2) {
            return Err("GeomBndLib_Line2d::Box - bad parameter".into());
        } else if Precision::is_positive_infinite(u2) {
            open_min_max_2d(lin.direction(), &mut bnd_box);
            bnd_box.add_point(&line_value_2d(0.0, lin));
        } else {
            open_min_2d(lin.direction(), &mut bnd_box);
            bnd_box.add_point(&line_value_2d(u2, lin));
        }
    } else if Precision::is_positive_infinite(u1) {
        if Precision::is_negative_infinite(u2) {
            open_min_max_2d(lin.direction(), &mut bnd_box);
            bnd_box.add_point(&line_value_2d(0.0, lin));
        } else if Precision::is_positive_infinite(u2) {
            return Err("GeomBndLib_Line2d::Box - bad parameter".into());
        } else {
            open_max_2d(lin.direction(), &mut bnd_box);
            bnd_box.add_point(&line_value_2d(u2, lin));
        }
    } else {
        bnd_box.add_point(&line_value_2d(u1, lin));
        if Precision::is_negative_infinite(u2) {
            open_min_2d(lin.direction(), &mut bnd_box);
        } else if Precision::is_positive_infinite(u2) {
            open_max_2d(lin.direction(), &mut bnd_box);
        } else {
            bnd_box.add_point(&line_value_2d(u2, lin));
        }
    }
    bnd_box.enlarge(tol);
    Ok(bnd_box)
}

/// `GeomBndLib_Line2d::BoxOptimal` — identical to [`box_lin`].
pub fn box_lin_optimal(
    lin: &GpLin2d,
    u1: f64,
    u2: f64,
    tol: f64,
) -> Result<BndBox2d, String> {
    box_lin_2d(lin, u1, u2, tol)
}
