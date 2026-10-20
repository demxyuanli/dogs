//! `GetIntersection` (`gxx:1573-1783`).

use super::prelude::*;
use super::super::curve_tool;

/// `GetIntersection(...)` (`gxx:1573-1783`).
#[allow(clippy::too_many_arguments)]
pub(super) fn intersection(
    c1: &dyn Curve2d,
    t1f: f64,
    t1l: f64,
    c2: &dyn Curve2d,
    t2f: f64,
    t2l: f64,
    tol_conf: f64,
    max_count: i32,
    p_int: &mut IntRes2dIntersectionPoint,
    dist: &mut f64,
    count: &mut i32,
) {
    *count += 1;
    //
    let tol2 = tol_conf * tol_conf;
    let p_tol1 = (100.0 * epsilon(t1f.abs().max(t1l.abs()))).max(PCONFUSION);
    let p_tol2 = (100.0 * epsilon(t2f.abs().max(t2l.abs()))).max(PCONFUSION);
    //
    let p1f = curve_tool::geom2d_d0(c1, t1f);
    let p1l = curve_tool::geom2d_d0(c1, t1l);
    let mut b1 = BndBox2d::new();
    b1.add_point(&p1f);
    b1.add_point(&p1l);
    b1.enlarge(tol_conf);
    //
    let p2f = curve_tool::geom2d_d0(c2, t2f);
    let p2l = curve_tool::geom2d_d0(c2, t2l);
    let mut b2 = BndBox2d::new();
    b2.add_point(&p2f);
    b2.add_point(&p2l);
    b2.enlarge(tol_conf);
    //
    if b1.is_out_box(&b2) {
        *count -= 1;
        return;
    }
    //
    let is_small1 = (t1l - t1f) <= p_tol1
        || p1f.square_distance(&p1l) / 4.0 <= tol2;
    let is_small2 = (t2l - t2f) <= p_tol2
        || p2f.square_distance(&p2l) / 4.0 <= tol2;

    if (is_small1 && is_small2) || (*count > max_count) {
        // Seems to be intersection
        // Simple treatment of segment intersection
        let pnts1 = [
            *p1f.xy(),
            p1f.xy().added(&p1l.xy()).divided(2.0),
            *p1l.xy(),
        ];
        let pnts2 = [
            *p2f.xy(),
            p2f.xy().added(&p2l.xy()).divided(2.0),
            *p2l.xy(),
        ];
        let mut imin: i32 = -1;
        let mut jmin: i32 = -1;
        let mut dmin = f64::MAX;
        for i in 0..3 {
            for j in 0..3 {
                let d: f64 = pnts1[i].subtracted(&pnts2[j]).square_modulus();
                if d < dmin {
                    dmin = d;
                    imin = i as i32;
                    jmin = j as i32;
                }
            }
        }
        //
        dmin = dmin.sqrt();
        if *dist > dmin {
            *dist = dmin;
            //
            let t1 = match imin {
                0 => t1f,
                1 => (t1f + t1l) / 2.0,
                _ => t1l,
            };
            //
            let t2 = match jmin {
                0 => t2f,
                1 => (t2f + t2l) / 2.0,
                _ => t2l,
            };
            //
            let pint = GpPnt2d::from_xy(pnts1[imin as usize].added(&pnts2[jmin as usize]).divided(2.0));
            //
            let trans1 = IntRes2dTransition::new();
            let trans2 = IntRes2dTransition::new();
            p_int.set_values(&pint, t1, t2, &trans1, &trans2, false);
        }
        *count -= 1;
        return;
    }

    if is_small1 {
        let t2m = (t2l + t2f) / 2.0;
        intersection(
            c1, t1f, t1l, c2, t2f, t2m, tol_conf, max_count, p_int,
            dist, count,
        );
        intersection(
            c1, t1f, t1l, c2, t2m, t2l, tol_conf, max_count, p_int,
            dist, count,
        );
    } else if is_small2 {
        let t1m = (t1l + t1f) / 2.0;
        intersection(
            c1, t1f, t1m, c2, t2f, t2l, tol_conf, max_count, p_int,
            dist, count,
        );
        intersection(
            c1, t1m, t1l, c2, t2f, t2l, tol_conf, max_count, p_int,
            dist, count,
        );
    } else {
        let t1m = (t1l + t1f) / 2.0;
        let t2m = (t2l + t2f) / 2.0;
        intersection(
            c1, t1f, t1m, c2, t2f, t2m, tol_conf, max_count, p_int,
            dist, count,
        );
        intersection(
            c1, t1f, t1m, c2, t2m, t2l, tol_conf, max_count, p_int,
            dist, count,
        );
        intersection(
            c1, t1m, t1l, c2, t2f, t2m, tol_conf, max_count, p_int,
            dist, count,
        );
        intersection(
            c1, t1m, t1l, c2, t2m, t2l, tol_conf, max_count, p_int,
            dist, count,
        );
    }
}
