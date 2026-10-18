//! Bounding box of a 2D parabola.
//!
//! Source: `GeomBndLib_Parabola2d.cxx`. Infinite parameters open X/Y min or
//! max together. Finite interval adds endpoints and the vertex at `t=0` when
//! the parameter range crosses zero.
//! T-97: items below are faithful ports of the named OCCT source, but their
//! OCCT-side consumers are not all ported yet, so parts are not called from this
//! crate. The `dead_code` allowance is deliberate: **pending wiring**, not dead
//! code. Do not delete them to silence warnings (see
//! specs/_a3n00_gap_analysis.md §9.309/§9.310); wire the consumer instead.
#![allow(dead_code)]

use occt_core::bnd::BndBox2d;
use occt_core::gp::GpParab2d;
use occt_core::precision::Precision;

use crate::geom_bnd_lib_elclib2d::parab_value;

/// `GeomBndLib_Parabola2d::Box(theParab, theU1, theU2, theTol)`.
pub fn box_parab(
    parab: &GpParab2d,
    u1: f64,
    u2: f64,
    tol: f64,
) -> Result<BndBox2d, String> {
    let mut bnd_box = BndBox2d::new();
    if Precision::is_negative_infinite(u1) {
        if Precision::is_negative_infinite(u2) {
            return Err("GeomBndLib_Parabola2d::Box - bad parameter".into());
        } else if Precision::is_positive_infinite(u2) {
            bnd_box.open_xmax();
            bnd_box.open_ymax();
        } else {
            bnd_box.add_point(&parab_value(u2, parab));
        }
        bnd_box.open_xmin();
        bnd_box.open_ymin();
    } else if Precision::is_positive_infinite(u1) {
        if Precision::is_negative_infinite(u2) {
            bnd_box.open_xmin();
            bnd_box.open_ymin();
        } else if Precision::is_positive_infinite(u2) {
            return Err("GeomBndLib_Parabola2d::Box - bad parameter".into());
        } else {
            bnd_box.add_point(&parab_value(u2, parab));
        }
        bnd_box.open_xmax();
        bnd_box.open_ymax();
    } else {
        bnd_box.add_point(&parab_value(u1, parab));
        if Precision::is_negative_infinite(u2) {
            bnd_box.open_xmin();
            bnd_box.open_ymin();
        } else if Precision::is_positive_infinite(u2) {
            bnd_box.open_xmax();
            bnd_box.open_ymax();
        } else {
            bnd_box.add_point(&parab_value(u2, parab));
            if u1 * u2 < 0.0 {
                bnd_box.add_point(&parab_value(0.0, parab));
            }
        }
    }
    bnd_box.enlarge(tol);
    Ok(bnd_box)
}
