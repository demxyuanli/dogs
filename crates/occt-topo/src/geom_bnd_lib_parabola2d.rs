//! Bounding box of a 2D parabola.
//!
//! Source: `GeomBndLib_Parabola2d.cxx`. Infinite parameters open X/Y min or
//! max together. Finite interval adds endpoints and the vertex at `t=0` when
//! the parameter range crosses zero.

use occt_core::bnd::BndBox2d;
use occt_core::gp::GpParab2d;
use occt_core::precision::Precision;

use crate::geom_bnd_lib_elclib2d::parab_value;

/// `GeomBndLib_Parabola2d::Box(theParab, theU1, theU2, theTol)`.
pub fn box_parab(
    the_parab: &GpParab2d,
    the_u1: f64,
    the_u2: f64,
    the_tol: f64,
) -> Result<BndBox2d, String> {
    let mut a_box = BndBox2d::new();
    if Precision::is_negative_infinite(the_u1) {
        if Precision::is_negative_infinite(the_u2) {
            return Err("GeomBndLib_Parabola2d::Box - bad parameter".into());
        } else if Precision::is_positive_infinite(the_u2) {
            a_box.open_xmax();
            a_box.open_ymax();
        } else {
            a_box.add_point(&parab_value(the_u2, the_parab));
        }
        a_box.open_xmin();
        a_box.open_ymin();
    } else if Precision::is_positive_infinite(the_u1) {
        if Precision::is_negative_infinite(the_u2) {
            a_box.open_xmin();
            a_box.open_ymin();
        } else if Precision::is_positive_infinite(the_u2) {
            return Err("GeomBndLib_Parabola2d::Box - bad parameter".into());
        } else {
            a_box.add_point(&parab_value(the_u2, the_parab));
        }
        a_box.open_xmax();
        a_box.open_ymax();
    } else {
        a_box.add_point(&parab_value(the_u1, the_parab));
        if Precision::is_negative_infinite(the_u2) {
            a_box.open_xmin();
            a_box.open_ymin();
        } else if Precision::is_positive_infinite(the_u2) {
            a_box.open_xmax();
            a_box.open_ymax();
        } else {
            a_box.add_point(&parab_value(the_u2, the_parab));
            if the_u1 * the_u2 < 0.0 {
                a_box.add_point(&parab_value(0.0, the_parab));
            }
        }
    }
    a_box.enlarge(the_tol);
    Ok(a_box)
}
