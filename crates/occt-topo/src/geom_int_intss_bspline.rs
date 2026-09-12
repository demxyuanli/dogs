//! `GeomInt_IntSS::MakeBSpline` / `MakeBSpline2d`. Source: `GeomInt_IntSS_1.cxx:1452`.

use std::sync::Arc;

use occt_geom::{Curve, GeomBSplineCurve};
use occt_geom2d::curve::Curve2d;
use occt_geom2d::Geom2dBSplineCurve;

use crate::int_tools_wline::WLine;

fn degree1_polyline_knots(n: usize) -> Vec<f64> {
    let mut k = Vec::with_capacity(n + 2);
    k.push(0.0);
    for i in 0..n {
        k.push(i as f64);
    }
    k.push((n - 1) as f64);
    k
}

/// `GeomInt_IntSS::MakeBSpline`.
pub fn make_bspline(wl: &WLine, ideb: i32, ifin: i32) -> Option<Arc<dyn Curve>> {
    if ifin < ideb {
        return None;
    }
    let nbpnt = (ifin - ideb + 1) as usize;
    if nbpnt < 2 {
        return None;
    }
    let mut poles = Vec::with_capacity(nbpnt);
    for i in 0..nbpnt {
        poles.push(wl.point(ideb + i as i32).value());
    }
    let knots = degree1_polyline_knots(nbpnt);
    GeomBSplineCurve::new(poles, knots, 1)
        .ok()
        .map(|c| Arc::new(c) as Arc<dyn Curve>)
}

/// `GeomInt_IntSS::MakeBSpline2d`.
pub fn make_bspline2d(wl: &WLine, ideb: i32, ifin: i32, on_first: bool) -> Option<Arc<dyn Curve2d>> {
    if ifin < ideb {
        return None;
    }
    let nbpnt = (ifin - ideb + 1) as usize;
    if nbpnt < 2 {
        return None;
    }
    let mut xs = Vec::with_capacity(nbpnt);
    let mut ys = Vec::with_capacity(nbpnt);
    for i in 0..nbpnt {
        let p = wl.point(ideb + i as i32);
        let (u, v) = if on_first {
            p.parameters_on_s1()
        } else {
            p.parameters_on_s2()
        };
        xs.push(u);
        ys.push(v);
    }
    let knots = degree1_polyline_knots(nbpnt);
    Geom2dBSplineCurve::new(xs, ys, knots, 1)
        .ok()
        .map(|c| Arc::new(c) as Arc<dyn Curve2d>)
}
