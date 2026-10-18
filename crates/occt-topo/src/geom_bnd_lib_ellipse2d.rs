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
pub fn box_elips(elips: &GpElips2d, tol: f64) -> BndBox2d {
    let mut bnd_box = BndBox2d::new();
    let maj = elips.major_radius;
    let min_r = elips.minor_radius;
    let o = elips.pos.point;
    let xd = elips.pos.vxdir;
    let yd = elips.pos.vydir;
    let mut min = [0.0; 2];
    let mut max = [0.0; 2];
    for k in 1..=2 {
        let xk = if k == 1 { xd.x() } else { xd.y() };
        let yk = if k == 1 { yd.x() } else { yd.y() };
        let amp = (maj * maj * xk * xk + min_r * min_r * yk * yk).sqrt();
        let ok = if k == 1 { o.x() } else { o.y() };
        min[k - 1] = ok - amp;
        max[k - 1] = ok + amp;
    }
    bnd_box.update(min[0], min[1], max[0], max[1]);
    bnd_box.enlarge(tol);
    bnd_box
}

/// `GeomBndLib_Ellipse2d::Box(theElips, theU1, theU2, theTol)`.
pub fn box_elips_range(elips: &GpElips2d, u1: f64, u2: f64, tol: f64) -> BndBox2d {
    if u2 - u1 >= TWO_PI - PCONFUSION {
        return box_elips(elips, tol);
    }
    let mut bnd_box = BndBox2d::new();
    bnd_box.add_point(&elips_value(u1, elips));
    bnd_box.add_point(&elips_value(u2, elips));
    let maj = elips.major_radius;
    let min_r = elips.minor_radius;
    let xd = elips.pos.vxdir;
    let yd = elips.pos.vydir;
    for k in 1..=2 {
        let xk = if k == 1 { xd.x() } else { xd.y() };
        let yk = if k == 1 { yd.x() } else { yd.y() };
        let t_extr_min = if xk.abs() > RESOLUTION {
            in_period(((min_r * yk) / (maj * xk)).atan(), 0.0, TWO_PI)
        } else {
            std::f64::consts::FRAC_PI_2
        };
        let t_extr_max = if t_extr_min <= std::f64::consts::PI {
            t_extr_min + std::f64::consts::PI
        } else {
            t_extr_min - std::f64::consts::PI
        };
        let mut tk = in_period(t_extr_min, u1, u1 + TWO_PI);
        if tk >= u1 && tk <= u2 {
            bnd_box.add_point(&elips_value(t_extr_min, elips));
        }
        tk = in_period(t_extr_max, u1, u1 + TWO_PI);
        if tk >= u1 && tk <= u2 {
            bnd_box.add_point(&elips_value(t_extr_max, elips));
        }
    }
    bnd_box.enlarge(tol);
    bnd_box
}
