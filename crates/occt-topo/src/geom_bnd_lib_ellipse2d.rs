//! Bounding box of a 2D ellipse / elliptical arc.
//!
//! Source: `GeomBndLib_Ellipse2d.cxx`. Full ellipse: per-coordinate amplitude
//! `sqrt(Major^2 Xd.k^2 + Minor^2 Yd.k^2)`. Arc: endpoints plus
//! `atan((Minor Yk) / (Major Xk))` extrema wrapped with `ElCLib::InPeriod`.
//! T-97: items below are faithful ports of the named OCCT source, but their
//! OCCT-side consumers are not all ported yet, so parts are not called from this
//! crate. The `dead_code` allowance is deliberate: **pending wiring**, not dead
//! code. Do not delete them to silence warnings (see
//! specs/_a3n00_gap_analysis.md §9.309/§9.310); wire the consumer instead.
#![allow(dead_code)]

use occt_core::bnd::BndBox2d;
use occt_core::gp::GpElips2d;
use occt_core::precision::{PCONFUSION, RESOLUTION};

use crate::geom_bnd_lib_elclib2d::{elips_value, in_period};

const TWO_PI: f64 = 2.0 * std::f64::consts::PI;

/// `GeomBndLib_Ellipse2d::Box(theElips, theTol)` — full ellipse.
pub fn box_elips(the_elips: &GpElips2d, the_tol: f64) -> BndBox2d {
    let mut a_box = BndBox2d::new();
    let a_maj = the_elips.major_radius;
    let a_min_r = the_elips.minor_radius;
    let a_o = the_elips.pos.point;
    let a_xd = the_elips.pos.vxdir;
    let a_yd = the_elips.pos.vydir;
    let mut a_min = [0.0; 2];
    let mut a_max = [0.0; 2];
    for k in 1..=2 {
        let a_xk = if k == 1 { a_xd.x() } else { a_xd.y() };
        let a_yk = if k == 1 { a_yd.x() } else { a_yd.y() };
        let a_amp = (a_maj * a_maj * a_xk * a_xk + a_min_r * a_min_r * a_yk * a_yk).sqrt();
        let ok = if k == 1 { a_o.x() } else { a_o.y() };
        a_min[k - 1] = ok - a_amp;
        a_max[k - 1] = ok + a_amp;
    }
    a_box.update(a_min[0], a_min[1], a_max[0], a_max[1]);
    a_box.enlarge(the_tol);
    a_box
}

/// `GeomBndLib_Ellipse2d::Box(theElips, theU1, theU2, theTol)`.
pub fn box_elips_range(the_elips: &GpElips2d, the_u1: f64, the_u2: f64, the_tol: f64) -> BndBox2d {
    if the_u2 - the_u1 >= TWO_PI - PCONFUSION {
        return box_elips(the_elips, the_tol);
    }
    let mut a_box = BndBox2d::new();
    a_box.add_point(&elips_value(the_u1, the_elips));
    a_box.add_point(&elips_value(the_u2, the_elips));
    let a_maj = the_elips.major_radius;
    let a_min_r = the_elips.minor_radius;
    let a_xd = the_elips.pos.vxdir;
    let a_yd = the_elips.pos.vydir;
    for k in 1..=2 {
        let a_xk = if k == 1 { a_xd.x() } else { a_xd.y() };
        let a_yk = if k == 1 { a_yd.x() } else { a_yd.y() };
        let a_t_extr_min = if a_xk.abs() > RESOLUTION {
            in_period(((a_min_r * a_yk) / (a_maj * a_xk)).atan(), 0.0, TWO_PI)
        } else {
            std::f64::consts::FRAC_PI_2
        };
        let a_t_extr_max = if a_t_extr_min <= std::f64::consts::PI {
            a_t_extr_min + std::f64::consts::PI
        } else {
            a_t_extr_min - std::f64::consts::PI
        };
        let mut a_tk = in_period(a_t_extr_min, the_u1, the_u1 + TWO_PI);
        if a_tk >= the_u1 && a_tk <= the_u2 {
            a_box.add_point(&elips_value(a_t_extr_min, the_elips));
        }
        a_tk = in_period(a_t_extr_max, the_u1, the_u1 + TWO_PI);
        if a_tk >= the_u1 && a_tk <= the_u2 {
            a_box.add_point(&elips_value(a_t_extr_max, the_elips));
        }
    }
    a_box.enlarge(the_tol);
    a_box
}
