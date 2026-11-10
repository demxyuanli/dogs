//! The file-local helpers of `Bisector_BisecCC.cxx`: `Curvature`
//! (`:572-588`), `PointByInt` (`:1328-1481`), `ProjOnCurve` (`:1785-1833`),
//! `TestExtension` (`:1838-1876`) and `DiscretPar` (`:1935-1959`).

use occt_core::gp::{GpAx2d, GpDir2d, GpPnt2d, GpVec2d};
use occt_core::precision::{CONFUSION, INFINITE, PCONFUSION};

use crate::curve::Curve2d;
use crate::geom2d_api_project_point_on_curve::ProjectPointOnCurve2d;
use crate::geom2d_int::Geom2dIntGInter;
use crate::line::Geom2dLine;

use super::super::bisec_pc::BisectorBisecPC;
use super::super::is_convex;

/// `gp_Dir2d(X, Y)` with the OCCT throw replaced by a zero direction (see
/// `bisec_pc.rs`).
fn dir2d_or_zero(x: f64, y: f64) -> GpDir2d {
    GpDir2d::new(x, y).unwrap_or(GpDir2d { x: 0.0, y: 0.0 })
}

/// `Curvature(C, U, Tol)` (`cxx:572-588`).
pub(super) fn curvature(c: &dyn Curve2d, u: f64, tol: f64) -> f64 {
    let (_, d1, d2) = c.d2(u);
    let norm2 = d1.square_magnitude();
    if norm2 < tol {
        0.0
    } else {
        d1.crossed(&d2) / (norm2 * norm2.sqrt())
    }
}

/// `PointByInt(CA, CB, SignA, SignB, UOnA, UOnB, Dist)` (`cxx:1328-1481`).
pub(super) fn point_by_int(
    ca: &dyn Curve2d,
    cb: &dyn Curve2d,
    sign_a: f64,
    sign_b: f64,
    u_on_a: f64,
    u_on_b: &mut f64,
    dist: &mut f64,
) -> bool {
    let is_convex_a = is_convex(ca, sign_a);
    let is_convex_b = is_convex(cb, sign_b);

    let (p1, tan1) = ca.d1(u_on_a);
    let n1 = GpVec2d::new(tan1.y(), -tan1.x());

    if p1.distance(&cb.value(cb.first_parameter())) < CONFUSION {
        *u_on_b = cb.first_parameter();
        let (_, tan2) = cb.d1(*u_on_b);
        if is_convex_a && is_convex_b {
            *dist = 0.0;
            return true;
        }
        if !tan1.is_parallel(&tan2, occt_core::precision::ANGULAR) {
            *dist = 0.0;
            return false;
        }
    }
    if p1.distance(&cb.value(cb.last_parameter())) < CONFUSION {
        *u_on_b = cb.last_parameter();
        let (_, tan2) = cb.d1(*u_on_b);
        if is_convex_a && is_convex_b {
            *dist = 0.0;
            return true;
        }
        if !tan1.is_parallel(&tan2, occt_core::precision::ANGULAR) {
            *dist = 0.0;
            return false;
        }
    }

    let mut d_min = INFINITE;
    let mut upc;
    let mut ya_sol = false;
    let mut p_sol = GpPnt2d::zero();

    let bis_pc = BisectorBisecPC::with_dist_max(cb, &p1, sign_b, INFINITE);
    if bis_pc.is_empty() {
        *dist = INFINITE;
        return false;
    }

    let nor_li = Geom2dLine::new(GpAx2d::new(p1, dir2d_or_zero(n1.x(), n1.y())));

    let mut intersect = Geom2dIntGInter::new();
    intersect.perform_curves(&bis_pc, &nor_li, CONFUSION, CONFUSION);

    if intersect.is_done() && intersect.nb_points() > 0 {
        for i in 1..=intersect.nb_points() {
            let pt = intersect.point(i);
            if pt.param_on_second() * sign_a < PCONFUSION {
                let p = *pt.value();
                if p.square_distance(&p1) < d_min {
                    d_min = p.square_distance(&p1);
                    p_sol = p;
                    upc = pt.param_on_first();
                    *u_on_b = bis_pc.link_bis_curve(upc);
                    *dist = d_min;
                    ya_sol = true;
                }
            }
        }
    }
    if ya_sol {
        let p2 = cb.value(*u_on_b);
        if p1.square_distance(&p_sol) < 1.0e-32 {
            return false;
        }
        if p2.square_distance(&p_sol) < 1.0e-32 {
            return false;
        }

        let pp1_unit = dir2d_or_zero(p1.x() - p_sol.x(), p1.y() - p_sol.y());
        let pp2_unit = dir2d_or_zero(p2.x() - p_sol.x(), p2.y() - p_sol.y());

        if pp1_unit.dot(&pp2_unit) > 1.0 - occt_core::precision::ANGULAR {
            ya_sol = false;
        } else {
            *dist = dist.sqrt();
            if !is_convex_a {
                let k1 = curvature(ca, u_on_a, CONFUSION);
                if k1 != 0.0 && *dist > (1.0 / k1).abs() {
                    ya_sol = false;
                }
            }
            if ya_sol && !is_convex_b {
                let k2 = curvature(cb, *u_on_b, CONFUSION);
                if k2 != 0.0 && *dist > (1.0 / k2).abs() {
                    ya_sol = false;
                }
            }
        }
    }
    ya_sol
}

/// `ProjOnCurve(P, C, theParam)` (`cxx:1785-1833`).
pub(super) fn proj_on_curve(p: &GpPnt2d, c: &dyn Curve2d, the_param: &mut f64) -> bool {
    *the_param = 0.0;

    let (pf, mut tf) = c.d1(c.first_parameter());
    let (pl, mut tl) = c.d1(c.last_parameter());

    if p.is_equal(&pf, CONFUSION) {
        *the_param = c.first_parameter();
        return true;
    }
    if p.is_equal(&pl, CONFUSION) {
        *the_param = c.last_parameter();
        return true;
    }

    let ppf = GpVec2d::new(pf.x() - p.x(), pf.y() - p.y());
    let _ = tf.normalize();
    if ppf.dot(&tf).abs() < CONFUSION {
        *the_param = c.first_parameter();
        return true;
    }
    let ppl = GpVec2d::new(pl.x() - p.x(), pl.y() - p.y());
    let _ = tl.normalize();
    if ppl.dot(&tl).abs() < CONFUSION {
        *the_param = c.last_parameter();
        return true;
    }
    let proj = ProjectPointOnCurve2d::with_curve_range(p, &std::sync::Arc::from(c.clone_dyn()), c.first_parameter(), c.last_parameter());
    match proj.lower_distance_parameter() {
        Some(u) => {
            *the_param = u;
            true
        }
        None => false,
    }
}

/// `TestExtension(C1, C2, Start_End)` (`cxx:1838-1876`).
pub(super) fn test_extension(c1: &dyn Curve2d, c2: &dyn Curve2d, start_end: i32) -> bool {
    let (p1, mut t1) = if start_end == 1 {
        c1.d1(c1.first_parameter())
    } else {
        c1.d1(c1.last_parameter())
    };
    let (p2, mut t2) = c2.d1(c2.first_parameter());
    if p1.is_equal(&p2, CONFUSION) {
        let _ = t1.normalize();
        let _ = t2.normalize();
        return t1.dot(&t2) > 1.0 - CONFUSION;
    }
    let (p2, mut t2) = c2.d1(c2.last_parameter());
    if p1.is_equal(&p2, CONFUSION) {
        let _ = t2.normalize();
        return t1.dot(&t2) > 1.0 - CONFUSION;
    }
    false
}

/// `DiscretPar(DU, EpsMin, EpsMax, NbMin, NbMax, Eps, Nb)` (`cxx:1935-1959`).
pub(super) fn discret_par(
    du: f64,
    eps_min: f64,
    eps_max: f64,
    nb_min: i32,
    nb_max: i32,
) -> (f64, i32) {
    if du <= (nb_min as f64) * eps_min {
        return (du / ((nb_min + 1) as f64), nb_min);
    }

    let mut eps = eps_max.min(du / (nb_max as f64));
    let nb;
    if eps < eps_min {
        eps = eps_min;
        nb = (du / eps_min) as i32;
    } else {
        nb = nb_max;
    }
    (eps, nb)
}

/// `Epsilon(theValue)` (`Standard_Real.hxx:242-246`).
pub(super) fn epsilon_of(v: f64) -> f64 {
    if v == 0.0 {
        f64::from_bits(1)
    } else if v >= 0.0 {
        f64::from_bits(v.to_bits() + 1) - v
    } else {
        v - f64::from_bits(v.to_bits() - 1)
    }
}