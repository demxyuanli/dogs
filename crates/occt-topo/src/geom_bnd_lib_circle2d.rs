//! Bounding box of a 2D circle / circular arc.
//!
//! Source: `GeomBndLib_Circle2d.cxx`. Full circle: per-coordinate analytical
//! extrema `O.k ± sqrt(R^2 Xd.k^2 + R^2 Yd.k^2)`. Arc: endpoints plus
//! `atan(Yk/Xk)` extrema wrapped with `ElCLib::InPeriod`.

use occt_core::bnd::BndBox2d;
use occt_core::gp::GpCirc2d;
use occt_core::precision::{PCONFUSION, RESOLUTION};

use crate::geom_bnd_lib_elclib2d::{circ_value, in_period};

const TWO_PI: f64 = 2.0 * std::f64::consts::PI;

/// `GeomBndLib_Circle2d::Box(theCirc, theTol)` — full circle.
pub fn box_circ(circ: &GpCirc2d, tol: f64) -> BndBox2d {
    let mut bnd_box = BndBox2d::new();
    let r = circ.radius();
    let o = circ.location();
    let xd = *circ.position().x_direction();
    let yd = *circ.position().y_direction();
    let mut min = [0.0; 2];
    let mut max = [0.0; 2];
    for k in 1..=2 {
        let xk = if k == 1 { xd.x() } else { xd.y() };
        let yk = if k == 1 { yd.x() } else { yd.y() };
        let amp = (r * r * xk * xk + r * r * yk * yk).sqrt();
        let ok = if k == 1 { o.x() } else { o.y() };
        min[k - 1] = ok - amp;
        max[k - 1] = ok + amp;
    }
    bnd_box.update(min[0], min[1], max[0], max[1]);
    bnd_box.enlarge(tol);
    bnd_box
}

/// `GeomBndLib_Circle2d::Box(theCirc, theU1, theU2, theTol)` — arc or full.
pub fn box_circ_range_2d(circ: &GpCirc2d, u1: f64, u2: f64, tol: f64) -> BndBox2d {
    if u2 - u1 >= TWO_PI - PCONFUSION {
        return box_circ(circ, tol);
    }
    let mut bnd_box = BndBox2d::new();
    bnd_box.add_point(&circ_value(u1, circ));
    bnd_box.add_point(&circ_value(u2, circ));
    let xd = *circ.position().x_direction();
    let yd = *circ.position().y_direction();
    for k in 1..=2 {
        let xk = if k == 1 { xd.x() } else { xd.y() };
        let yk = if k == 1 { yd.x() } else { yd.y() };
        let t_extr_min = if xk.abs() > RESOLUTION {
            in_period((yk / xk).atan(), 0.0, TWO_PI)
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
            bnd_box.add_point(&circ_value(t_extr_min, circ));
        }
        tk = in_period(t_extr_max, u1, u1 + TWO_PI);
        if tk >= u1 && tk <= u2 {
            bnd_box.add_point(&circ_value(t_extr_max, circ));
        }
    }
    bnd_box.enlarge(tol);
    bnd_box
}
