//! `IntPatch_TheSearchInside` (`IntStart_SearchInside.gxx`).

use occt_core::gp::{GpDir2d, GpPnt, GpVec};
use occt_core::precision::CONFUSION;
use occt_geom::Surface;

use crate::fclass2d::FaceState;
use crate::geom_int::TopolTool;
use crate::int_tools_wline::{u_resolution, v_resolution};

use super::surf_func::{clamp_box, function_set_root, SurfFunction};

/// `IntSurf_InteriorPoint`.
#[derive(Clone, Copy)]
pub(crate) struct InteriorPoint {
    pub p: GpPnt,
    pub u: f64,
    pub v: f64,
    pub d3d: GpVec,
    pub d2d: GpDir2d,
}

/// `IntStart_SearchInside::Perform(Func, PS, T, Epsilon)`.
pub(crate) fn search_inside(
    func: &mut SurfFunction<'_>,
    surf: &dyn Surface,
    domain: &TopolTool,
    epsilon: f64,
) -> Vec<InteriorPoint> {
    let mut list = Vec::new();
    let (umin0, umax0) = surf.u_range();
    let (vmin0, vmax0) = surf.v_range();
    let (mut umin, mut umax) = finite_pair(umin0, umax0);
    let (mut vmin, mut vmax) = finite_pair(vmin0, vmax0);

    let nbs_u = domain.nb_samples_u(surf);
    let nbs_v = domain.nb_samples_v(surf);
    let nbs = domain.nb_samples(surf);
    if nbs_u <= 0 || nbs_v <= 0 || nbs <= 0 {
        return list;
    }

    let mut du = umax - umin;
    let mut dv = vmax - vmin;
    du /= nbs_u as f64 * 0.5;
    dv /= nbs_v as f64 * 0.5;

    let toler1 = u_resolution(surf, CONFUSION);
    let toler2 = v_resolution(surf, CONFUSION);
    let mut max_tol = toler1.max(toler2);
    max_tol *= 1000.0;
    if max_tol > du * 0.001 {
        max_tol = du * 0.001;
    }
    if max_tol > dv * 0.001 {
        max_tol = dv * 0.001;
    }

    umin += du * 0.01;
    vmin += dv * 0.01;
    umax -= du * 0.01;
    vmax -= dv * 0.01;

    let tol = func.tolerance();
    for i in 1..=nbs + 12 {
        let mut skip = false;
        let (uvap, binf, bsup) = if i <= nbs {
            let (s2d, _) = domain.sample_point(surf, i);
            let uv = (s2d.x(), s2d.y());
            let u1 = (uv.0 - du).max(umin);
            let v1 = (uv.1 - dv).max(vmin);
            let u2 = (uv.0 + du).min(umax);
            let v2 = (uv.1 + dv).min(vmax);
            let p1 = surf.d0(u1, v1);
            let p2 = surf.d0(u2, v2);
            let rvalf = func.value(uv.0, uv.1);
            let dist_pp = p1.square_distance(&p2);
            if rvalf * rvalf > 3.0 * dist_pp {
                skip = true;
            }
            (uv, (u1, v1), (u2, v2))
        } else {
            let s2d = extra_seed(i - nbs, umin, umax, vmin, vmax, du, dv);
            let uv = (s2d.0, s2d.1);
            let binf = ((uv.0 - du).max(umin), (uv.1 - dv).max(vmin));
            let bsup = ((uv.0 + du).min(umax), (uv.1 + dv).min(vmax));
            (uv, binf, bsup)
        };
        if skip {
            continue;
        }
        let (binf, bsup) = clamp_box(binf, bsup, umin, umax, vmin, vmax);
        let Some((u, v)) = function_set_root(func, uvap, binf, bsup, (toler1, toler2)) else {
            continue;
        };
        if func.root().abs() > tol || func.is_tangent() {
            continue;
        }
        let psol = func.point();
        let dup = list.iter().any(|ip: &InteriorPoint| {
            (ip.p.x() - psol.x()).abs() <= epsilon
                && (ip.p.y() - psol.y()).abs() <= epsilon
                && (ip.p.z() - psol.z()).abs() <= epsilon
                && (u - ip.u).abs() <= toler1
                && (v - ip.v).abs() <= toler2
        });
        if dup {
            continue;
        }
        if domain.classify(occt_core::gp::GpPnt2d::new(u, v), max_tol) != FaceState::In {
            continue;
        }
        list.push(InteriorPoint {
            p: psol,
            u,
            v,
            d3d: func.direction3d(),
            d2d: func.direction2d(),
        });
    }
    list
}

fn extra_seed(
    k: i32,
    umin: f64,
    umax: f64,
    vmin: f64,
    vmax: f64,
    du: f64,
    dv: f64,
) -> (f64, f64) {
    match k {
        1 | 5 => (umin + du * 0.02, vmin + dv * 0.02),
        2 | 6 => (umax - du * 0.02, vmin + dv * 0.02),
        3 | 7 => (umin + du * 0.02, vmax - dv * 0.02),
        4 | 8 => (umax - du * 0.02, vmax - dv * 0.02),
        9 => (umin + du * 0.005, vmin + dv * 0.005),
        10 => (umax - du * 0.005, vmin + dv * 0.005),
        11 => (umin + du * 0.005, vmax - dv * 0.005),
        _ => (umax - du * 0.005, vmax - dv * 0.005),
    }
}

fn finite_pair(a: f64, b: f64) -> (f64, f64) {
    let mut lo = a;
    let mut hi = b;
    if hi < lo {
        std::mem::swap(&mut lo, &mut hi);
    }
    if !lo.is_finite() && !hi.is_finite() {
        return (-1.0e5, 1.0e5);
    }
    if !lo.is_finite() {
        lo = hi - 2.0e5;
    }
    if !hi.is_finite() {
        hi = lo + 2.0e5;
    }
    (lo, hi)
}
