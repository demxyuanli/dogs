//! `IntPatch_WLineTool::ComputePurgedWLine`.

use occt_core::gp::{GpPnt, GpPnt2d, GpVec, GpVec2d};
use occt_core::precision::{CONFUSION, RESOLUTION};
use occt_geom::Surface;

use crate::brep_surface::{classify_surface, SurfaceKind};
use crate::fclass2d::FaceState;
use crate::geom_int::TopolTool;
use crate::int_tools_wline::{u_resolution, v_resolution, WLine};

/// `IntPatch_WLineTool::ComputePurgedWLine`.
pub fn compute_purged_wline(
    wl: &WLine,
    s1: &dyn Surface,
    s2: &dyn Surface,
    d1: &TopolTool,
    d2: &TopolTool,
) -> Option<WLine> {
    let nb = wl.nb_pnts();
    if nb == 2 {
        let p1 = wl.point(1).p;
        let p2 = wl.point(2).p;
        if p1.is_equal(&p2) {
            return None;
        }
    }
    let mut local = clone_wline(wl);
    delete_equal_points(&mut local);
    if local.nb_pnts() < 2 {
        return None;
    }
    local = delete_outer_points(&local, s1, s2, d1, d2);
    if local.nb_pnts() < 2 {
        return None;
    }
    Some(delete_by_tube(&local, s1, s2))
}

fn clone_wline(wl: &WLine) -> WLine {
    let mut out = WLine::new();
    out.set_creating_way(wl.creating_way());
    for i in 1..=wl.nb_pnts() {
        out.add(*wl.point(i));
    }
    out.vertices = wl.vertices.clone();
    out.has_first_point = wl.has_first_point;
    out.has_last_point = wl.has_last_point;
    out
}

fn delete_equal_points(wl: &mut WLine) {
    let mut i = 1i32;
    while i <= wl.nb_pnts() {
        let start = i + 1;
        let mut end = i + 5;
        let nb = wl.nb_pnts();
        if end > nb {
            end = nb;
        }
        if start > nb || end <= 1 {
            i += 1;
            continue;
        }
        let mut k = start;
        while k <= end {
            if i != k {
                let p1 = *wl.point(i);
                let p2 = *wl.point(k);
                let uv = [
                    p1.u1, p1.v1, p1.u2, p1.v2, p2.u1, p2.v1, p2.u2, p2.v2,
                ];
                let mut a_max = uv[0].abs();
                for x in &uv[1..] {
                    a_max = a_max.max(x.abs());
                }
                if p1.p.is_equal(&p2.p)
                    || (p1.u1 - p2.u1).abs() + (p1.v1 - p2.v1).abs() < 1.0e-16 * a_max
                    || (p1.u2 - p2.u2).abs() + (p1.v2 - p2.v2).abs() < 1.0e-16 * a_max
                {
                    wl.remove_point(k);
                    continue;
                }
            }
            k += 1;
        }
        i += 1;
    }
}

fn delete_outer_points(
    wl: &WLine,
    s1: &dyn Surface,
    s2: &dyn Surface,
    d1: &TopolTool,
    d2: &TopolTool,
) -> WLine {
    if s1.is_u_periodic() || s1.is_v_periodic() || s2.is_u_periodic() || s2.is_v_periodic() {
        return clone_wline(wl);
    }
    let n = wl.nb_pnts();
    let mut drop = vec![false; n as usize];
    let mut first = 1i32;
    for i in 1..=n {
        let (x1, y1, x2, y2) = wl.point(i).parameters();
        let st1 = d1.classify(GpPnt2d::new(x1, y1), CONFUSION);
        let st2 = d2.classify(GpPnt2d::new(x2, y2), CONFUSION);
        if st1 == FaceState::Out || st2 == FaceState::Out {
            drop[(i - 1) as usize] = true;
        } else {
            first = i;
            break;
        }
    }
    if drop.iter().all(|d| *d) {
        return clone_wline(wl);
    }
    let mut last = n;
    for i in (1..=n).rev() {
        let (x1, y1, x2, y2) = wl.point(i).parameters();
        let st1 = d1.classify(GpPnt2d::new(x1, y1), CONFUSION);
        let st2 = d2.classify(GpPnt2d::new(x2, y2), CONFUSION);
        if st1 == FaceState::Out || st2 == FaceState::Out {
            drop[(i - 1) as usize] = true;
        } else {
            last = i;
            break;
        }
    }
    let mut out = WLine::new();
    out.set_creating_way(wl.creating_way());
    for i in first.max(1)..=last.min(n) {
        if !drop[(i - 1) as usize] {
            out.add(*wl.point(i));
        }
    }
    if out.nb_pnts() < 2 {
        return clone_wline(wl);
    }
    out.ensure_end_vertices();
    out
}

fn delete_by_tube(wl: &WLine, s1: &dyn Surface, s2: &dyn Surface) -> WLine {
    if wl.nb_pnts() <= 2 {
        return clone_wline(wl);
    }
    let base_tol = 1.0e-3;
    let res1 = u_resolution(s1, base_tol).min(v_resolution(s1, base_tol));
    let res2 = u_resolution(s2, base_tol).min(v_resolution(s2, base_tol));
    let tol1 = res1 * res1;
    let tol2 = res2 * res2;
    let tol3d = base_tol * base_tol;
    let limit = 0.99 * 0.99;
    let plane_plane = is_plane_like(s1) && is_plane_like(s2);
    let n = wl.nb_pnts();
    let mut keep = vec![true; n as usize];
    let mut b1 = *wl.point(1);
    let mut b2 = *wl.point(2);
    for i in 3..=n {
        let cur = *wl.point(i);
        let v1 = GpVec2d::new(b2.u1 - b1.u1, b2.v1 - b1.v1);
        let v2 = GpVec2d::new(b2.u2 - b1.u2, b2.v2 - b1.v2);
        let v3 = GpVec::from_pnts(&b1.p, &b2.p);
        let inside = is_inside_2d(b1.u1, b1.v1, v1, cur.u1, cur.v1, tol1)
            && is_inside_2d(b1.u2, b1.v2, v2, cur.u2, cur.v2, tol2)
            && is_inside_3d(&b1.p, &v3, &cur.p, tol3d);
        if inside && keep[(i - 2) as usize] {
            let s1a = (b2.u1 - b1.u1).hypot(b2.v1 - b1.v1);
            let s1b = (cur.u1 - b2.u1).hypot(cur.v1 - b2.v1);
            let s2a = (b2.u2 - b1.u2).hypot(b2.v2 - b1.v2);
            let s2b = (cur.u2 - b2.u2).hypot(cur.v2 - b2.v2);
            let step1 = if s1b > RESOLUTION { (s1a * s1a) / (s1b * s1b) } else { 0.0 };
            let step2 = if s2b > RESOLUTION { (s2a * s2a) / (s2b * s2b) } else { 0.0 };
            if step1.min(step2) >= limit * step1.max(step2) && !plane_plane {
                keep[(i - 2) as usize] = false;
            }
        }
        if keep[(i - 2) as usize] {
            b1 = b2;
            b2 = cur;
        } else {
            b2 = cur;
        }
    }
    let mut out = WLine::new();
    out.set_creating_way(wl.creating_way());
    for i in 1..=n {
        if keep[(i - 1) as usize] {
            out.add(*wl.point(i));
        }
    }
    if out.nb_pnts() < 2 {
        return clone_wline(wl);
    }
    out.ensure_end_vertices();
    out
}

fn is_plane_like(s: &dyn Surface) -> bool {
    classify_surface(s) == SurfaceKind::Plane
}

fn is_inside_2d(px: f64, py: f64, v: GpVec2d, qx: f64, qy: f64, tol: f64) -> bool {
    let w = GpVec2d::new(qx - px, qy - py);
    let l2 = v.square_magnitude();
    if l2 <= f64::EPSILON {
        return w.square_magnitude() <= tol;
    }
    let t = w.dot(&v) / l2;
    if !(0.0..=1.0).contains(&t) {
        return false;
    }
    let hx = w.x() - t * v.x();
    let hy = w.y() - t * v.y();
    hx * hx + hy * hy <= tol
}

fn is_inside_3d(p: &GpPnt, v: &GpVec, q: &GpPnt, tol: f64) -> bool {
    let w = GpVec::from_pnts(p, q);
    let l2 = v.square_magnitude();
    if l2 <= f64::EPSILON {
        return w.square_magnitude() <= tol;
    }
    let t = w.dot(v) / l2;
    if !(0.0..=1.0).contains(&t) {
        return false;
    }
    let h = w.subtracted(&v.multiplied_scalar(t));
    h.square_magnitude() <= tol
}
