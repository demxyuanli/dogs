//! Bounding box of a 2D hyperbola.
//!
//! Source: `GeomBndLib_Hyperbola2d.cxx`. Infinite parameters open all four
//! box directions (OCCT opens Xmin/Ymin or Xmax/Ymax together). Finite arcs
//! add endpoints, the vertex at `t=0` when the interval crosses zero, and
//! per-coordinate `t = 0.5 * log(|B-A| / |A+B|)` extrema.

use occt_core::bnd::BndBox2d;
use occt_core::gp::GpHypr2d;
use occt_core::precision::Precision;

use crate::geom_bnd_lib_elclib2d::hypr_value;

/// Finite-arc extrema (`computeHyperbola2dBox`).
fn compute_hyperbola2d_box(the_hypr: &GpHypr2d, the_t1: f64, the_t2: f64, the_box: &mut BndBox2d) {
    the_box.add_point(&hypr_value(the_t1, the_hypr));
    the_box.add_point(&hypr_value(the_t2, the_hypr));
    if the_t1 * the_t2 < 0.0 {
        the_box.add_point(&hypr_value(0.0, the_hypr));
    }
    let a_xdir = the_hypr.pos.vxdir;
    let a_ydir = the_hypr.pos.vydir;
    let a_rmaj = the_hypr.major_radius;
    let a_rmin = the_hypr.minor_radius;
    let a_eps = f64::EPSILON;
    for i in 1..=2 {
        let a_a = a_rmin * if i == 1 { a_ydir.x() } else { a_ydir.y() };
        let a_b = a_rmaj * if i == 1 { a_xdir.x() } else { a_xdir.y() };
        let a_abp = (a_a + a_b).abs();
        let a_bam = (a_b - a_a).abs();
        if a_abp < a_eps || a_bam < a_eps {
            continue;
        }
        let a_cf = a_bam / a_abp;
        let a_t3 = 0.5 * a_cf.ln();
        if a_t3 < the_t1 || a_t3 > the_t2 {
            continue;
        }
        the_box.add_point(&hypr_value(a_t3, the_hypr));
    }
}

/// `GeomBndLib_Hyperbola2d::Box(theHypr, theU1, theU2, theTol)`.
pub fn box_hypr(the_hypr: &GpHypr2d, the_u1: f64, the_u2: f64, the_tol: f64) -> Result<BndBox2d, String> {
    let mut a_box = BndBox2d::new();
    if Precision::is_negative_infinite(the_u1) {
        if Precision::is_negative_infinite(the_u2) {
            return Err("GeomBndLib_Hyperbola2d::Box - bad parameter".into());
        } else if Precision::is_positive_infinite(the_u2) {
            a_box.open_xmax();
            a_box.open_ymax();
        } else {
            a_box.add_point(&hypr_value(the_u2, the_hypr));
        }
        a_box.open_xmin();
        a_box.open_ymin();
    } else if Precision::is_positive_infinite(the_u1) {
        if Precision::is_negative_infinite(the_u2) {
            a_box.open_xmin();
            a_box.open_ymin();
        } else if Precision::is_positive_infinite(the_u2) {
            return Err("GeomBndLib_Hyperbola2d::Box - bad parameter".into());
        } else {
            a_box.add_point(&hypr_value(the_u2, the_hypr));
        }
        a_box.open_xmax();
        a_box.open_ymax();
    } else {
        a_box.add_point(&hypr_value(the_u1, the_hypr));
        if Precision::is_negative_infinite(the_u2) {
            a_box.open_xmin();
            a_box.open_ymin();
        } else if Precision::is_positive_infinite(the_u2) {
            a_box.open_xmax();
            a_box.open_ymax();
        } else {
            compute_hyperbola2d_box(the_hypr, the_u1, the_u2, &mut a_box);
        }
    }
    a_box.enlarge(the_tol);
    Ok(a_box)
}
