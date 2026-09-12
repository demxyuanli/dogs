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
pub fn box_circ(the_circ: &GpCirc2d, the_tol: f64) -> BndBox2d {
    let mut a_box = BndBox2d::new();
    let a_r = the_circ.radius();
    let a_o = the_circ.location();
    let a_xd = *the_circ.position().x_direction();
    let a_yd = *the_circ.position().y_direction();
    let mut a_min = [0.0; 2];
    let mut a_max = [0.0; 2];
    for k in 1..=2 {
        let a_xk = if k == 1 { a_xd.x() } else { a_xd.y() };
        let a_yk = if k == 1 { a_yd.x() } else { a_yd.y() };
        let a_amp = (a_r * a_r * a_xk * a_xk + a_r * a_r * a_yk * a_yk).sqrt();
        let ok = if k == 1 { a_o.x() } else { a_o.y() };
        a_min[k - 1] = ok - a_amp;
        a_max[k - 1] = ok + a_amp;
    }
    a_box.update(a_min[0], a_min[1], a_max[0], a_max[1]);
    a_box.enlarge(the_tol);
    a_box
}

/// `GeomBndLib_Circle2d::Box(theCirc, theU1, theU2, theTol)` — arc or full.
pub fn box_circ_range(the_circ: &GpCirc2d, the_u1: f64, the_u2: f64, the_tol: f64) -> BndBox2d {
    if the_u2 - the_u1 >= TWO_PI - PCONFUSION {
        return box_circ(the_circ, the_tol);
    }
    let mut a_box = BndBox2d::new();
    a_box.add_point(&circ_value(the_u1, the_circ));
    a_box.add_point(&circ_value(the_u2, the_circ));
    let a_xd = *the_circ.position().x_direction();
    let a_yd = *the_circ.position().y_direction();
    for k in 1..=2 {
        let a_xk = if k == 1 { a_xd.x() } else { a_xd.y() };
        let a_yk = if k == 1 { a_yd.x() } else { a_yd.y() };
        let a_t_extr_min = if a_xk.abs() > RESOLUTION {
            in_period((a_yk / a_xk).atan(), 0.0, TWO_PI)
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
            a_box.add_point(&circ_value(a_t_extr_min, the_circ));
        }
        a_tk = in_period(a_t_extr_max, the_u1, the_u1 + TWO_PI);
        if a_tk >= the_u1 && a_tk <= the_u2 {
            a_box.add_point(&circ_value(a_t_extr_max, the_circ));
        }
    }
    a_box.enlarge(the_tol);
    a_box
}
