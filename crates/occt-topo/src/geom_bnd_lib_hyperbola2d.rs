//! Bounding box of a 2D hyperbola.
//!
//! Source: `GeomBndLib_Hyperbola2d.cxx`. Infinite parameters open all four
//! box directions (OCCT opens Xmin/Ymin or Xmax/Ymax together). Finite arcs
//! add endpoints, the vertex at `t=0` when the interval crosses zero, and
//! per-coordinate `t = 0.5 * log(|B-A| / |A+B|)` extrema.
//! T-97: items below are faithful ports of the named OCCT source, but their
//! OCCT-side consumers are not all ported yet, so parts are not called from this
//! crate. The `dead_code` allowance is deliberate: **pending wiring**, not dead
//! code. Do not delete them to silence warnings (see
//! specs/_a3n00_gap_analysis.md §9.309/§9.310); wire the consumer instead.
#![allow(dead_code)]

use occt_core::bnd::BndBox2d;
use occt_core::gp::GpHypr2d;
use occt_core::precision::Precision;

use crate::geom_bnd_lib_elclib2d::hypr_value;

/// Finite-arc extrema (`computeHyperbola2dBox`).
fn compute_hyperbola2d_box(hypr: &GpHypr2d, t1: f64, t2: f64, bnd_box: &mut BndBox2d) {
    bnd_box.add_point(&hypr_value(t1, hypr));
    bnd_box.add_point(&hypr_value(t2, hypr));
    if t1 * t2 < 0.0 {
        bnd_box.add_point(&hypr_value(0.0, hypr));
    }
    let xdir = hypr.pos.vxdir;
    let ydir = hypr.pos.vydir;
    let rmaj = hypr.major_radius;
    let rmin = hypr.minor_radius;
    let eps = f64::EPSILON;
    for i in 1..=2 {
        let a = rmin * if i == 1 { ydir.x() } else { ydir.y() };
        let b = rmaj * if i == 1 { xdir.x() } else { xdir.y() };
        let abp = (a + b).abs();
        let bam = (b - a).abs();
        if abp < eps || bam < eps {
            continue;
        }
        let cf = bam / abp;
        let t3 = 0.5 * cf.ln();
        if t3 < t1 || t3 > t2 {
            continue;
        }
        bnd_box.add_point(&hypr_value(t3, hypr));
    }
}

/// `GeomBndLib_Hyperbola2d::Box(theHypr, theU1, theU2, theTol)`.
pub fn box_hypr(hypr: &GpHypr2d, u1: f64, u2: f64, tol: f64) -> Result<BndBox2d, String> {
    let mut a_box = BndBox2d::new();
    if Precision::is_negative_infinite(u1) {
        if Precision::is_negative_infinite(u2) {
            return Err("GeomBndLib_Hyperbola2d::Box - bad parameter".into());
        } else if Precision::is_positive_infinite(u2) {
            a_box.open_xmax();
            a_box.open_ymax();
        } else {
            a_box.add_point(&hypr_value(u2, hypr));
        }
        a_box.open_xmin();
        a_box.open_ymin();
    } else if Precision::is_positive_infinite(u1) {
        if Precision::is_negative_infinite(u2) {
            a_box.open_xmin();
            a_box.open_ymin();
        } else if Precision::is_positive_infinite(u2) {
            return Err("GeomBndLib_Hyperbola2d::Box - bad parameter".into());
        } else {
            a_box.add_point(&hypr_value(u2, hypr));
        }
        a_box.open_xmax();
        a_box.open_ymax();
    } else {
        a_box.add_point(&hypr_value(u1, hypr));
        if Precision::is_negative_infinite(u2) {
            a_box.open_xmin();
            a_box.open_ymin();
        } else if Precision::is_positive_infinite(u2) {
            a_box.open_xmax();
            a_box.open_ymax();
        } else {
            compute_hyperbola2d_box(hypr, u1, u2, &mut a_box);
        }
    }
    a_box.enlarge(tol);
    Ok(a_box)
}
