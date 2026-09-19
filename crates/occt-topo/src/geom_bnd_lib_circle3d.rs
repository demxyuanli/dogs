//! Port of `GeomBndLib_Circle.cxx` (`Box(gp_Circ, ...)`).
//!
//! Full circle: per-coordinate analytic extrema. Arc: endpoints plus any
//! extremal parameter that falls inside the arc.

use occt_core::bnd::BndBox;
use occt_core::elib::clib::circle_value;
use occt_core::gp::GpCirc;
use occt_core::precision::{epsilon, CONFUSION, PCONFUSION};

use crate::geom_bnd_lib_elclib2d::{adjust_periodic, in_period};

/// `GeomBndLib_Circle::Box(theCirc, theTol)` (`Circle.cxx:21-45`).
pub fn box_circ_full(the_circ: &GpCirc, the_tol: f64) -> BndBox {
    let mut a_box = BndBox::new();
    let a_r = the_circ.radius();
    let a_o = the_circ.location();
    let a_x_ax = the_circ.x_axis();
    let a_y_ax = the_circ.y_axis();
    let a_xd = a_x_ax.direction().xyz();
    let a_yd = a_y_ax.direction().xyz();
    let mut a_min = [0.0_f64; 3];
    let mut a_max = [0.0_f64; 3];
    for k in 0..3 {
        let a_xk = a_xd.coord(k);
        let a_yk = a_yd.coord(k);
        let a_amp = (a_r * a_r * a_xk * a_xk + a_r * a_r * a_yk * a_yk).sqrt();
        a_min[k] = a_o.xyz().coord(k) - a_amp;
        a_max[k] = a_o.xyz().coord(k) + a_amp;
    }
    a_box.update(a_min[0], a_min[1], a_min[2], a_max[0], a_max[1], a_max[2]);
    a_box.enlarge(the_tol);
    a_box
}

/// `GeomBndLib_Circle::Box(theCirc, theU1, theU2, theTol)` (`Circle.cxx:49-114`).
pub fn box_circ_range(the_circ: &GpCirc, the_u1: f64, the_u2: f64, the_tol: f64) -> BndBox {
    let a_period = 2.0 * std::f64::consts::PI - PCONFUSION;
    if the_u2 - the_u1 >= a_period {
        return box_circ_full(the_circ, the_tol);
    }

    let mut a_box = BndBox::new();
    let a_r = the_circ.radius();
    let a_o = the_circ.location();
    let a_x_ax = the_circ.x_axis();
    let a_y_ax = the_circ.y_axis();
    let a_xd = a_x_ax.direction().xyz();
    let a_yd = a_y_ax.direction().xyz();

    let mut a_u1 = the_u1;
    let mut a_u2 = the_u2;
    let a_tol = epsilon(1.0);
    adjust_periodic(0.0, 2.0 * std::f64::consts::PI, a_tol, &mut a_u1, &mut a_u2);

    a_box.add_point(&circle_value(the_circ, a_u1));
    a_box.add_point(&circle_value(the_circ, a_u2));
    for k in 0..3 {
        let a_xk = a_xd.coord(k);
        let a_yk = a_yd.coord(k);

        let mut a_t_extr_min;
        if a_xk.abs() > CONFUSION {
            a_t_extr_min = (a_yk / a_xk).atan();
            a_t_extr_min = in_period(a_t_extr_min, 0.0, 2.0 * std::f64::consts::PI);
        } else {
            a_t_extr_min = std::f64::consts::PI / 2.0;
        }
        let mut a_t_extr_max =
            if a_t_extr_min <= std::f64::consts::PI { a_t_extr_min + std::f64::consts::PI } else { a_t_extr_min - std::f64::consts::PI };

        let mut a_val_min = a_r * a_t_extr_min.cos() * a_xk
            + a_r * a_t_extr_min.sin() * a_yk
            + a_o.xyz().coord(k);
        let mut a_val_max = a_r * a_t_extr_max.cos() * a_xk
            + a_r * a_t_extr_max.sin() * a_yk
            + a_o.xyz().coord(k);
        if a_val_min > a_val_max {
            std::mem::swap(&mut a_val_min, &mut a_val_max);
            std::mem::swap(&mut a_t_extr_min, &mut a_t_extr_max);
        }

        let mut a_tk = in_period(a_t_extr_min, a_u1, a_u1 + 2.0 * std::f64::consts::PI);
        if a_tk >= a_u1 && a_tk <= a_u2 {
            a_box.add_point(&circle_value(the_circ, a_t_extr_min));
        }
        a_tk = in_period(a_t_extr_max, a_u1, a_u1 + 2.0 * std::f64::consts::PI);
        if a_tk >= a_u1 && a_tk <= a_u2 {
            a_box.add_point(&circle_value(the_circ, a_t_extr_max));
        }
    }

    a_box.enlarge(the_tol);
    a_box
}
