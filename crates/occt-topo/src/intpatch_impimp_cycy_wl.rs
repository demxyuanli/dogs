//! Domain inscription into a walking line, V-boundary 1D roots, densify.
//! Source: `AddPointIntoWL` / `AddBoundaryPoint` / `SeekAdditionalPoints`.

use occt_core::elib::slib;
use occt_core::gp::{GpCylinder, GpPnt};
use occt_core::precision::PCONFUSION;
use occt_math::brent;

use crate::int_tools_wline::{PntOn2S, WLine};

use super::{
    compute_params, inscribe_point, Coeffs, SurfDom,
};

fn mid_pnt(a: &GpPnt, b: &GpPnt) -> GpPnt {
    GpPnt::new(
        0.5 * (a.x() + b.x()),
        0.5 * (a.y() + b.y()),
        0.5 * (a.z() + b.z()),
    )
}

fn make_pnt(
    cyl1: &GpCylinder,
    cyl2: &GpCylinder,
    reversed: bool,
    u1: f64,
    v1: f64,
    u2: f64,
    v2: f64,
) -> PntOn2S {
    let p1 = slib::cylinder_value(cyl1, u1, v1);
    let p2 = slib::cylinder_value(cyl2, u2, v2);
    let p = mid_pnt(&p1, &p2);
    if reversed {
        PntOn2S {
            p,
            u1: u2,
            v1: v2,
            u2: u1,
            v2: v1,
        }
    } else {
        PntOn2S {
            p,
            u1,
            v1,
            u2,
            v2,
        }
    }
}

fn last_u1(line: &WLine, reversed: bool) -> Option<(f64, f64)> {
    let n = line.nb_pnts();
    if n < 1 {
        return None;
    }
    let p = line.point(n);
    Some(if reversed {
        p.parameters_on_s2()
    } else {
        p.parameters_on_s1()
    })
}

/// `AddPointIntoWL`.
pub(crate) fn add_point_into_wl(
    cyl1: &GpCylinder,
    cyl2: &GpCylinder,
    coeffs: &Coeffs,
    reversed: bool,
    precise: bool,
    u1: f64,
    v1: f64,
    u2: f64,
    v2: f64,
    dom: SurfDom,
    line: &mut WLine,
    wl: i32,
    fl_before: bool,
    only_check: bool,
) -> bool {
    let mid = 0.5 * (dom.u1f + dom.u1l);
    let mut u1par = u1;
    if !inscribe_point(dom.u1f, dom.u1l, &mut u1par, dom.tol2d, dom.period, u1 > mid) {
        return false;
    }
    if line.nb_pnts() > 0
        && (dom.u1l - dom.u1f) >= (dom.period - dom.tol2d)
        && (((u1par + dom.period - dom.u1l) <= dom.tol2d)
            || ((u1par - dom.period - dom.u1f) >= dom.tol2d))
    {
        if let Some((ul, _)) = last_u1(line, reversed) {
            let delta = ul - u1par;
            if 2.0 * delta.abs() > dom.period {
                u1par += delta.signum() * dom.period;
            }
        }
    }
    let mut u2par = u2;
    if !inscribe_point(dom.u2f, dom.u2l, &mut u2par, dom.tol2d, dom.period, false) {
        return false;
    }
    if (v1 - dom.v1l > dom.tol2d) || (dom.v1f - v1 > dom.tol2d) {
        return false;
    }
    if (v2 - dom.v2l > dom.tol2d) || (dom.v2f - v2 > dom.tol2d) {
        return false;
    }
    if let Some((ul, _)) = last_u1(line, reversed) {
        if !fl_before && u1par <= ul {
            u1par += dom.period;
            if (dom.u1f - u1par > dom.tol2d) || (u1par - dom.u1l > dom.tol2d) {
                return false;
            }
        }
        if only_check {
            return true;
        }
        let pnt = make_pnt(cyl1, cyl2, reversed, u1par, v1, u2par, v2);
        let d_tol = 1.0 - f64::EPSILON;
        let last = *line.point(line.nb_pnts());
        if pnt.is_same(&last, dom.tol3d * d_tol, dom.tol2d * d_tol) {
            line.remove_point(line.nb_pnts());
        }
        line.add(pnt);
    } else {
        if only_check {
            return true;
        }
        line.add(make_pnt(cyl1, cyl2, reversed, u1par, v1, u2par, v2));
    }
    if !precise {
        return true;
    }
    let n = line.nb_pnts();
    if n >= 3 {
        let (u_a, _) = if reversed {
            line.point(n - 2).parameters_on_s2()
        } else {
            line.point(n - 2).parameters_on_s1()
        };
        let (u_b, _) = if reversed {
            line.point(n - 1).parameters_on_s2()
        } else {
            line.point(n - 1).parameters_on_s1()
        };
        let (u_c, _) = if reversed {
            line.point(n).parameters_on_s2()
        } else {
            line.point(n).parameters_on_s1()
        };
        let step_prev = u_b - u_a;
        let step = u_c - u_b;
        if step.abs() > 0.0 {
            let delta_step = (step_prev / step) as i32;
            if (1 < delta_step) && (delta_step < 2000) {
                seek_additional_points(
                    cyl1,
                    cyl2,
                    line,
                    coeffs,
                    wl,
                    delta_step,
                    n - 2,
                    n - 1,
                    dom.tol2d,
                    dom.period,
                    reversed,
                );
            }
        }
    }
    true
}

#[derive(Clone, Copy)]
struct StPInfo {
    u1: f64,
    u2: f64,
    v1: f64,
    v2: f64,
    surf_id: i32,
}

fn find_v_bound(
    coeffs: &Coeffs,
    wl: i32,
    is_v1: bool,
    v_bound: f64,
    u_lo: f64,
    u_hi: f64,
) -> Option<f64> {
    if !(u_lo < u_hi) {
        return None;
    }
    let eval = |x: f64| -> Option<f64> {
        let (_u2, v1, v2) = compute_params(x, wl, coeffs)?;
        Some(if is_v1 { v1 } else { v2 } - v_bound)
    };
    let flo = eval(u_lo)?;
    let fhi = eval(u_hi)?;
    if flo == 0.0 {
        return Some(u_lo);
    }
    if fhi == 0.0 {
        return Some(u_hi);
    }
    if flo * fhi > 0.0 {
        return None;
    }
    brent(&|x| eval(x).unwrap_or(f64::MAX), u_lo, u_hi, PCONFUSION).ok()
}

/// `WorkWithBoundaries::AddBoundaryPoint`.
pub(crate) fn add_boundary_point(
    cyl1: &GpCylinder,
    cyl2: &GpCylinder,
    coeffs: &Coeffs,
    reversed: bool,
    line: &mut WLine,
    u1: f64,
    u1_prev: f64,
    u1_min: f64,
    u2: f64,
    v1: f64,
    v1_prev: f64,
    v2: f64,
    v2_prev: f64,
    wl: i32,
    fl_force: bool,
    dom: SurfDom,
) -> (bool, bool) {
    let vz = [dom.v1f, dom.v1l, dom.v2f, dom.v2l];
    let mut pts = [StPInfo {
        u1: f64::MAX,
        u2: 0.0,
        v1: 0.0,
        v2: 0.0,
        surf_id: 0,
    }; 4];
    for id_surf in (0..4).step_by(2) {
        let (vf, vl) = if id_surf == 0 {
            (v1, v1_prev)
        } else {
            (v2, v2_prev)
        };
        let is_v1 = id_surf == 0;
        for id_bound in 0..2 {
            let idx = id_surf + id_bound;
            pts[idx].surf_id = id_surf as i32;
            if (vf - vz[idx]).abs() > dom.tol2d && (vf - vz[idx]) * (vl - vz[idx]) > 0.0 {
                continue;
            }
            let u_lo = u1_prev.min(u1);
            let u_hi = u1_prev.max(u1);
            let Some(u_star) = find_v_bound(coeffs, wl, is_v1, vz[idx], u_lo, u_hi) else {
                pts[idx].u1 = f64::MAX;
                continue;
            };
            if u_star >= u1 || u_star < u1_min {
                pts[idx].u1 = f64::MAX;
                continue;
            }
            let Some((uu2, mut vv1, mut vv2)) = compute_params(u_star, wl, coeffs) else {
                pts[idx].u1 = f64::MAX;
                continue;
            };
            if is_v1 {
                vv1 = vz[idx];
            } else {
                vv2 = vz[idx];
            }
            pts[idx] = StPInfo {
                u1: u_star,
                u2: uu2,
                v1: vv1,
                v2: vv2,
                surf_id: id_surf as i32,
            };
            let _ = (u2, v1, v2);
        }
    }
    pts.sort_by(|a, b| a.u1.partial_cmp(&b.u1).unwrap_or(std::cmp::Ordering::Equal));
    let mut found1 = false;
    let mut found2 = false;
    for p in pts {
        if p.u1 == f64::MAX {
            break;
        }
        if !add_point_into_wl(
            cyl1, cyl2, coeffs, reversed, false, p.u1, p.v1, p.u2, p.v2, dom, line, wl, fl_force,
            false,
        ) {
            continue;
        }
        if p.surf_id == 0 {
            found1 = true;
        } else {
            found2 = true;
        }
    }
    (found1, found2)
}

/// `SeekAdditionalPoints`. Indices are 1-based inclusive.
pub(crate) fn seek_additional_points(
    cyl1: &GpCylinder,
    cyl2: &GpCylinder,
    line: &mut WLine,
    coeffs: &Coeffs,
    wl: i32,
    min_nb: i32,
    start: i32,
    end: i32,
    tol2d: f64,
    period: f64,
    reversed: bool,
) {
    if start < 1 || end < start {
        return;
    }
    let mut nb = end - start + 1;
    let (u1s, _) = if reversed {
        line.point(start).parameters_on_s2()
    } else {
        line.point(start).parameters_on_s1()
    };
    let (u1e, _) = if reversed {
        line.point(end).parameters_on_s2()
    } else {
        line.point(end).parameters_on_s1()
    };
    let min_du = ((u1e - u1s).abs() / min_nb as f64).max(tol2d);
    let mut last = end;
    loop {
        let prev = nb;
        let mut fp = start;
        while fp < last {
            let lp = fp + 1;
            let (u1f, _) = if reversed {
                line.point(fp).parameters_on_s2()
            } else {
                line.point(fp).parameters_on_s1()
            };
            let (u1l, _) = if reversed {
                line.point(lp).parameters_on_s2()
            } else {
                line.point(lp).parameters_on_s1()
            };
            let (mut u2f, _) = if reversed {
                line.point(fp).parameters_on_s1()
            } else {
                line.point(fp).parameters_on_s2()
            };
            let (mut u2l, _) = if reversed {
                line.point(lp).parameters_on_s1()
            } else {
                line.point(lp).parameters_on_s2()
            };
            if (u1l - u1f).abs() <= min_du {
                fp = lp;
                continue;
            }
            let u1p = 0.5 * (u1f + u1l);
            let Some((mut u2p, v1p, v2p)) = compute_params(u1p, wl, coeffs) else {
                fp = lp;
                continue;
            };
            if u2f > u2l {
                std::mem::swap(&mut u2f, &mut u2l);
            }
            if !inscribe_point(u2f, u2l, &mut u2p, tol2d, period, false) {
                fp = lp;
                continue;
            }
            let pnt = make_pnt(cyl1, cyl2, reversed, u1p, v1p, u2p, v2p);
            line.insert_before(lp, pnt);
            nb += 1;
            last += 1;
            fp = lp + 1;
        }
        if nb >= min_nb || nb == prev {
            return;
        }
    }
}
