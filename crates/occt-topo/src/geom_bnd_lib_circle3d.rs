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
pub fn box_circ_full(circ: &GpCirc, tol: f64) -> BndBox {
    let mut bnd_box = BndBox::new();
    let r = circ.radius();
    let o = circ.location();
    let x_ax = circ.x_axis();
    let y_ax = circ.y_axis();
    let xd = x_ax.direction().xyz();
    let yd = y_ax.direction().xyz();
    let mut min = [0.0_f64; 3];
    let mut max = [0.0_f64; 3];
    for k in 0..3 {
        let xk = xd.coord(k);
        let yk = yd.coord(k);
        let amp = (r * r * xk * xk + r * r * yk * yk).sqrt();
        min[k] = o.xyz().coord(k) - amp;
        max[k] = o.xyz().coord(k) + amp;
    }
    bnd_box.update(min[0], min[1], min[2], max[0], max[1], max[2]);
    bnd_box.enlarge(tol);
    bnd_box
}

/// `GeomBndLib_Circle::Box(theCirc, theU1, theU2, theTol)` (`Circle.cxx:49-114`).
pub fn box_circ_range(circ: &GpCirc, u1: f64, u2: f64, tol: f64) -> BndBox {
    let period = 2.0 * std::f64::consts::PI - PCONFUSION;
    if u2 - u1 >= period {
        return box_circ_full(circ, tol);
    }

    let mut bnd_box = BndBox::new();
    let r = circ.radius();
    let o = circ.location();
    let x_ax = circ.x_axis();
    let y_ax = circ.y_axis();
    let xd = x_ax.direction().xyz();
    let yd = y_ax.direction().xyz();

    let mut a_u1 = u1;
    let mut a_u2 = u2;
    let a_tol = epsilon(1.0);
    adjust_periodic(0.0, 2.0 * std::f64::consts::PI, a_tol, &mut a_u1, &mut a_u2);

    bnd_box.add_point(&circle_value(circ, a_u1));
    bnd_box.add_point(&circle_value(circ, a_u2));
    for k in 0..3 {
        let xk = xd.coord(k);
        let yk = yd.coord(k);

        let mut t_extr_min;
        if xk.abs() > CONFUSION {
            t_extr_min = (yk / xk).atan();
            t_extr_min = in_period(t_extr_min, 0.0, 2.0 * std::f64::consts::PI);
        } else {
            t_extr_min = std::f64::consts::PI / 2.0;
        }
        let mut t_extr_max =
            if t_extr_min <= std::f64::consts::PI { t_extr_min + std::f64::consts::PI } else { t_extr_min - std::f64::consts::PI };

        let mut val_min = r * t_extr_min.cos() * xk
            + r * t_extr_min.sin() * yk
            + o.xyz().coord(k);
        let mut val_max = r * t_extr_max.cos() * xk
            + r * t_extr_max.sin() * yk
            + o.xyz().coord(k);
        if val_min > val_max {
            std::mem::swap(&mut val_min, &mut val_max);
            std::mem::swap(&mut t_extr_min, &mut t_extr_max);
        }

        let mut tk = in_period(t_extr_min, a_u1, a_u1 + 2.0 * std::f64::consts::PI);
        if tk >= a_u1 && tk <= a_u2 {
            bnd_box.add_point(&circle_value(circ, t_extr_min));
        }
        tk = in_period(t_extr_max, a_u1, a_u1 + 2.0 * std::f64::consts::PI);
        if tk >= a_u1 && tk <= a_u2 {
            bnd_box.add_point(&circle_value(circ, t_extr_max));
        }
    }

    bnd_box.enlarge(tol);
    bnd_box
}
