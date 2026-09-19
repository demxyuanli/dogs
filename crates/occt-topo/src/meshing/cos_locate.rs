//! `Adaptor3d_CurveOnSurface` EvalFirstLastSurf / LocatePart stack.
//!
//! Source: `Adaptor3d_CurveOnSurface.cxx` Hunt/FindBounds/Locate1Coord/
//! Locate2Coord/LocatePart/LocatePart_RevExt/LocatePart_Offset/EvalFirstLastSurf
//! (`cxx:140-156`, `161-178`, `215-279`, `283-453`, `457-666`, `670-775`,
//! `779-869`, `1736-1924`).

use std::sync::Arc;

use occt_core::bspl::locate::{first_u_knot_index, last_u_knot_index};
use occt_core::gp::{GpPnt2d, GpVec2d};
use occt_core::precision::PCONFUSION;
use occt_geom::bspline_surface::GeomBSplineSurface;
use occt_geom::curve::Curve;
use occt_geom::rectangular_trimmed::GeomRectangularTrimmedSurface;
use occt_geom::Surface;
use occt_geom2d::curve::Curve2d;

const TOL: f64 = PCONFUSION / 10.0;

fn reverse_param_f(a: f64, b: f64) -> (f64, f64) {
    if a > b {
        (b, a)
    } else {
        (a, b)
    }
}

fn reverse_param_i(a: i32, b: i32) -> (i32, i32) {
    if a > b {
        (b, a)
    } else {
        (a, b)
    }
}

fn compare_bounds(p1: &mut GpPnt2d, p2: &mut GpPnt2d) {
    // `Adaptor3d_CurveOnSurface.cxx:140-156`.
    if p1.x() > p2.x() {
        let lx = p1.x();
        p1.set_x(p2.x());
        p2.set_x(lx);
    }
    if p1.y() > p2.y() {
        let ly = p1.y();
        p1.set_y(p2.y());
        p2.set_y(ly);
    }
}

/// Exact-knot Hunt used by LocatePart (`cxx:161-178`), not `BSplCLib::Hunt`.
fn hunt_exact(arr: &[f64], coord: f64) -> Option<i32> {
    let mut i = 1i32;
    let upper = arr.len() as i32;
    while i <= upper && (coord - arr[(i - 1) as usize]).abs() > TOL {
        i += 1;
    }
    if i <= upper && (coord - arr[(i - 1) as usize]).abs() < TOL {
        Some(i)
    } else {
        None
    }
}

fn knot1(arr: &[f64], i: i32) -> f64 {
    arr[(i - 1) as usize]
}

fn find_bounds(arr: &[f64], coord: f64, der: f64, bound1: &mut i32, bound2: &mut i32, der_null: &mut bool) {
    // `Adaptor3d_CurveOnSurface.cxx:215-279`.
    let n = match hunt_exact(arr, coord) {
        Some(n) => n,
        None => return,
    };
    *der_null = false;
    if n == *bound1 {
        *der_null = der.abs() <= TOL;
        *bound1 = n;
        *bound2 = n + 1;
        return;
    }
    if n == *bound2 {
        *der_null = der.abs() <= TOL;
        *bound1 = n - 1;
        *bound2 = n;
        return;
    }
    if der.abs() > TOL {
        if der > 0.0 {
            *bound1 = n;
            *bound2 = n + 1;
        } else {
            *bound1 = n - 1;
            *bound2 = n;
        }
        *der_null = false;
    } else {
        *der_null = true;
        *bound1 = n - 1;
        *bound2 = n + 1;
    }
}

struct BsplKnots {
    u_knots: Vec<f64>,
    v_knots: Vec<f64>,
    u_first: i32,
    u_last: i32,
    v_first: i32,
    v_last: i32,
}

fn bspl_knots(bs: &GeomBSplineSurface) -> Option<BsplKnots> {
    let (u_knots, u_mults) = GeomBSplineSurface::unique_knots_mults(&bs.knots_u);
    let (v_knots, v_mults) = GeomBSplineSurface::unique_knots_mults(&bs.knots_v);
    if u_knots.len() < 2 || v_knots.len() < 2 {
        return None;
    }
    let u_first = first_u_knot_index(bs.deg_u as i32, &u_mults);
    let u_last = last_u_knot_index(bs.deg_u as i32, &u_mults);
    let v_first = first_u_knot_index(bs.deg_v as i32, &v_mults);
    let v_last = last_u_knot_index(bs.deg_v as i32, &v_mults);
    if u_first < 1 || u_last > u_knots.len() as i32 || u_first >= u_last {
        return None;
    }
    if v_first < 1 || v_last > v_knots.len() as i32 || v_first >= v_last {
        return None;
    }
    Some(BsplKnots {
        u_knots,
        v_knots,
        u_first,
        u_last,
        v_first,
        v_last,
    })
}

fn locate1_coord_surf(
    index: i32,
    uv: GpPnt2d,
    duv: GpVec2d,
    k: &BsplKnots,
    d_is_null: &mut bool,
    left_bot: &mut GpPnt2d,
    right_top: &mut GpPnt2d,
) {
    // `Adaptor3d_CurveOnSurface.cxx:457-666`.
    *d_is_null = false;
    let (comp1, dcomp1, up, down, knots) = if index == 1 {
        (uv.x(), duv.x(), k.u_last, k.u_first, &k.u_knots)
    } else {
        (uv.y(), duv.y(), k.v_last, k.v_first, &k.v_knots)
    };

    let mut i = down;
    while (knot1(knots, i) - comp1).abs() > TOL && i != up {
        i += 1;
    }
    let cur = knot1(knots, i);

    if (comp1 - cur).abs() <= TOL {
        let mut bnd1 = down;
        let mut bnd2 = up;
        find_bounds(knots, cur, dcomp1, &mut bnd1, &mut bnd2, d_is_null);
        let (b1, b2) = reverse_param_i(bnd1, bnd2);
        bnd1 = b1;
        bnd2 = b2;
        if !*d_is_null {
            if index == 1 {
                left_bot.set_x(knot1(knots, bnd1));
                right_top.set_x(knot1(knots, bnd2));
            } else {
                left_bot.set_y(knot1(knots, bnd1));
                right_top.set_y(knot1(knots, bnd2));
            }
        }
        return;
    }

    if index == 1 && comp1 < knot1(knots, down) {
        left_bot.set_x(knot1(knots, down));
        right_top.set_x(knot1(knots, down + 1));
        return;
    }
    if index == 2 && comp1 < knot1(knots, down) {
        left_bot.set_y(knot1(knots, down));
        right_top.set_y(knot1(knots, down + 1));
        return;
    }
    if index == 1 && comp1 > knot1(knots, up) {
        right_top.set_x(knot1(knots, up - 1));
        left_bot.set_x(knot1(knots, up));
        return;
    }
    if index == 2 && comp1 > knot1(knots, up) {
        right_top.set_y(knot1(knots, up - 1));
        left_bot.set_y(knot1(knots, up));
        return;
    }

    let mut f = knot1(knots, down);
    let mut l = knot1(knots, up);
    i = down;
    if !(comp1 < f) && !(comp1 > l) {
        while i < up {
            f = knot1(knots, i);
            l = knot1(knots, i + 1);
            if !(f >= comp1 || l <= comp1) {
                break;
            }
            i += 1;
        }
    } else {
        let (a, b) = reverse_param_f(f, l);
        f = a;
        l = b;
    }

    if i != up {
        if dcomp1.abs() > TOL {
            if index == 1 {
                if dcomp1 > 0.0 {
                    left_bot.set_x(comp1);
                    right_top.set_x(l);
                } else {
                    left_bot.set_x(f);
                    right_top.set_x(comp1);
                }
            } else if dcomp1 > 0.0 {
                left_bot.set_y(comp1);
                right_top.set_y(l);
            } else {
                left_bot.set_y(f);
                right_top.set_y(comp1);
            }
        } else if index == 1 {
            left_bot.set_x(f);
            right_top.set_x(l);
        } else {
            left_bot.set_y(f);
            right_top.set_y(l);
        }
    } else if index == 1 {
        left_bot.set_x(comp1);
        right_top.set_x(knot1(knots, i));
    } else {
        left_bot.set_y(comp1);
        right_top.set_y(knot1(knots, i));
    }
}

/// Curve-knot `Locate1Coord` (`Adaptor3d_CurveOnSurface.cxx:283-453`).
fn locate1_coord_curve(
    index: i32,
    uv: GpPnt2d,
    duv: GpVec2d,
    knots: &[f64],
    lo: i32,
    up: i32,
    left_bot: &mut GpPnt2d,
    right_top: &mut GpPnt2d,
) {
    let (comp1, dcomp1) = if index == 1 {
        (uv.x(), duv.x())
    } else {
        (uv.y(), duv.y())
    };

    let mut i = lo;
    while (knot1(knots, i) - comp1).abs() > TOL && i != up {
        i += 1;
    }
    let cur = knot1(knots, i);

    if (comp1 - cur).abs() <= TOL {
        let mut bnd1 = lo;
        let mut bnd2 = up;
        let mut d_is_null = false;
        find_bounds(knots, cur, dcomp1, &mut bnd1, &mut bnd2, &mut d_is_null);
        let (b1, b2) = reverse_param_i(bnd1, bnd2);
        bnd1 = b1;
        bnd2 = b2;
        if !d_is_null {
            if index == 1 {
                left_bot.set_x(knot1(knots, bnd1));
                right_top.set_x(knot1(knots, bnd2));
            } else {
                left_bot.set_y(knot1(knots, bnd1));
                right_top.set_y(knot1(knots, bnd2));
            }
        } else if (comp1 - knot1(knots, lo)).abs() <= TOL {
            if index == 1 {
                left_bot.set_x(knot1(knots, lo));
                right_top.set_x(knot1(knots, lo + 1));
            } else {
                left_bot.set_y(knot1(knots, lo));
                right_top.set_y(knot1(knots, lo + 1));
            }
        } else if (comp1 - knot1(knots, up)).abs() <= TOL {
            if index == 1 {
                left_bot.set_x(knot1(knots, up - 1));
                right_top.set_x(knot1(knots, up));
            } else {
                left_bot.set_y(knot1(knots, up - 1));
                right_top.set_y(knot1(knots, up));
            }
        } else if index == 1 {
            left_bot.set_x(knot1(knots, bnd1));
            right_top.set_x(knot1(knots, bnd2));
        } else {
            left_bot.set_y(knot1(knots, bnd1));
            right_top.set_y(knot1(knots, bnd2));
        }
        return;
    }

    // Coord != Knot (`cxx:377-451`).
    let mut f = 0.0;
    let mut l = 0.0;
    i = lo;
    while i < up {
        f = knot1(knots, i);
        l = knot1(knots, i + 1);
        if f < comp1 && l > comp1 {
            break;
        }
        i += 1;
    }
    let (rf, rl) = reverse_param_f(f, l);
    f = rf;
    l = rl;

    if i != up {
        if dcomp1.abs() < TOL {
            if index == 1 {
                left_bot.set_x(f);
                right_top.set_x(l);
            } else {
                left_bot.set_y(f);
                right_top.set_y(l);
            }
        } else if dcomp1.abs() > TOL {
            if index == 1 {
                if dcomp1 > 0.0 {
                    left_bot.set_x(comp1);
                    right_top.set_x(l);
                } else {
                    left_bot.set_x(f);
                    right_top.set_x(comp1);
                }
            } else if dcomp1 > 0.0 {
                left_bot.set_y(comp1);
                right_top.set_y(l);
            } else {
                left_bot.set_y(f);
                right_top.set_y(comp1);
            }
        }
    } else if index == 1 {
        left_bot.set_x(comp1);
        right_top.set_x(knot1(knots, i));
    } else {
        left_bot.set_y(comp1);
        right_top.set_y(knot1(knots, i));
    }
}

/// Param-range `Locate2Coord` (`Adaptor3d_CurveOnSurface.cxx:670-775`).
fn locate2_coord_param(
    index: i32,
    uv: GpPnt2d,
    duv: GpVec2d,
    i1: f64,
    i2: f64,
    left_bot: &mut GpPnt2d,
    right_top: &mut GpPnt2d,
) {
    let (comp1, dcomp1) = if index == 1 {
        (uv.x(), duv.x())
    } else {
        (uv.y(), duv.y())
    };

    // Exact `!=` matches cxx:691 (not a tol test).
    if comp1 != i1 && comp1 != i2 {
        if dcomp1.abs() > TOL {
            if dcomp1 < 0.0 {
                if index == 1 {
                    left_bot.set_x(i1);
                    right_top.set_x(comp1);
                } else {
                    left_bot.set_y(i1);
                    right_top.set_y(comp1);
                }
            } else if dcomp1 > 0.0 {
                if index == 1 {
                    left_bot.set_x(comp1);
                    right_top.set_x(i2);
                } else {
                    left_bot.set_y(comp1);
                    right_top.set_y(i2);
                }
            } else if index == 1 {
                left_bot.set_x(i1);
                right_top.set_x(i2);
            } else {
                left_bot.set_y(i1);
                right_top.set_y(i2);
            }
        } else if index == 1 {
            left_bot.set_x(i1);
            right_top.set_x(i2);
        } else {
            left_bot.set_y(i1);
            right_top.set_y(i2);
        }
    } else if (comp1 - i1).abs() < TOL {
        if index == 1 {
            left_bot.set_x(i1);
            right_top.set_x(i2);
        } else {
            left_bot.set_y(i1);
            right_top.set_y(i2);
        }
    } else if (comp1 - i2).abs() < TOL {
        if index == 1 {
            left_bot.set_x(i1);
            right_top.set_x(i2);
        } else {
            left_bot.set_y(i1);
            right_top.set_y(i2);
        }
    }
}

fn curve_unique_knots(curve: &dyn Curve) -> Option<(Vec<f64>, i32, i32)> {
    let flat = curve.bspline_knots()?;
    let deg = curve.nurbs_degree()? as i32;
    let (knots, mults) = GeomBSplineSurface::unique_knots_mults(flat);
    if knots.len() < 2 {
        return None;
    }
    let lo = first_u_knot_index(deg, &mults);
    let up = last_u_knot_index(deg, &mults);
    if lo < 1 || up > knots.len() as i32 || lo >= up {
        return None;
    }
    Some((knots, lo, up))
}

/// `Adaptor3d_CurveOnSurface::LocatePart_RevExt` (`cxx:1833-1866`).
fn locate_part_rev_ext(
    uv: GpPnt2d,
    duv: GpVec2d,
    surf: &dyn Surface,
    left_bot: &mut GpPnt2d,
    right_top: &mut GpPnt2d,
) -> bool {
    let basis = if surf.is_surface_of_linear_extrusion() {
        surf.extrusion_basis_curve()
    } else if surf.is_surface_of_revolution() {
        surf.revolution_basis_curve()
    } else {
        None
    };
    let Some(basis) = basis else {
        return false;
    };
    let Some((knots, lo, up)) = curve_unique_knots(basis.as_ref()) else {
        return false;
    };
    let (u0, u1) = surf.u_range();
    let (v0, v1) = surf.v_range();
    if surf.is_surface_of_linear_extrusion() {
        locate1_coord_curve(1, uv, duv, &knots, lo, up, left_bot, right_top);
        locate2_coord_param(2, uv, duv, v0, v1, left_bot, right_top);
    } else {
        // SurfaceOfRevolution: V follows BSpline generatrix; U is the spin period.
        locate1_coord_curve(2, uv, duv, &knots, lo, up, left_bot, right_top);
        locate2_coord_param(1, uv, duv, u0, u1, left_bot, right_top);
    }
    let (u1b, u2b) = reverse_param_f(left_bot.x(), right_top.x());
    left_bot.set_x(u1b);
    right_top.set_x(u2b);
    let (v1b, v2b) = reverse_param_f(left_bot.y(), right_top.y());
    left_bot.set_y(v1b);
    right_top.set_y(v2b);
    true
}

fn locate2_coord_arr(
    index: i32,
    uv: GpPnt2d,
    duv: GpVec2d,
    k: &BsplKnots,
    arr: &[f64],
    left_bot: &mut GpPnt2d,
    right_top: &mut GpPnt2d,
) {
    // `Adaptor3d_CurveOnSurface.cxx:779-869`.
    let (comp, dcomp, n_up, n_lo) = if index == 1 {
        (uv.x(), duv.y(), k.u_last, k.u_first)
    } else {
        (uv.y(), duv.x(), k.v_last, k.v_first)
    };

    if dcomp > 0.0 && dcomp.abs() > TOL {
        let mut n = hunt_exact(arr, comp).unwrap_or(n_up);
        if n >= n_up {
            n = n_up - 1;
        }
        let (tmp1, tmp2) = reverse_param_f(knot1(arr, n), knot1(arr, n + 1));
        if index == 1 {
            left_bot.set_x(tmp1);
            right_top.set_x(tmp2);
        } else {
            left_bot.set_y(tmp1);
            right_top.set_y(tmp2);
        }
    } else if dcomp < 0.0 && dcomp.abs() > TOL {
        let mut n = hunt_exact(arr, comp).unwrap_or(n_lo);
        if n <= n_lo {
            n = n_lo + 1;
        }
        let (tmp1, tmp2) = reverse_param_f(knot1(arr, n - 1), knot1(arr, n));
        if index == 1 {
            left_bot.set_x(tmp1);
            right_top.set_x(tmp2);
        } else {
            left_bot.set_y(tmp1);
            right_top.set_y(tmp2);
        }
    }
}

fn locate_part(
    uv: GpPnt2d,
    duv: GpVec2d,
    k: &BsplKnots,
    left_bot: &mut GpPnt2d,
    right_top: &mut GpPnt2d,
) {
    // `Adaptor3d_CurveOnSurface.cxx:1901-1924`.
    let mut du_null = false;
    let mut dv_null = false;
    locate1_coord_surf(1, uv, duv, k, &mut du_null, left_bot, right_top);
    locate1_coord_surf(2, uv, duv, k, &mut dv_null, left_bot, right_top);
    if du_null && !dv_null {
        locate2_coord_arr(1, uv, duv, k, &k.u_knots, left_bot, right_top);
    } else if dv_null && !du_null {
        locate2_coord_arr(2, uv, duv, k, &k.v_knots, left_bot, right_top);
    }
}

fn locate_part_offset(
    uv: GpPnt2d,
    duv: GpVec2d,
    basis: &dyn Surface,
    left_bot: &mut GpPnt2d,
    right_top: &mut GpPnt2d,
) -> bool {
    // `Adaptor3d_CurveOnSurface.cxx:1871-1896`.
    if basis.is_surface_of_revolution() || basis.is_surface_of_linear_extrusion() {
        return locate_part_rev_ext(uv, duv, basis, left_bot, right_top);
    }
    let Some(bs) = basis.osculating_bspline() else {
        return false;
    };
    let Some(k) = bspl_knots(&bs) else {
        return false;
    };
    locate_part(uv, duv, &k, left_bot, right_top);
    true
}

fn needs_first_last(surf: &dyn Surface) -> bool {
    // `Adaptor3d_CurveOnSurface::Load` cxx:953-963 — Offset unwrap then type check.
    if let Some(b) = surf.offset_basis_surface() {
        return b.is_bspline_surface()
            || b.is_surface_of_linear_extrusion()
            || b.is_surface_of_revolution()
            || b.osculating_bspline().is_some();
    }
    surf.is_bspline_surface()
        || surf.is_surface_of_linear_extrusion()
        || surf.is_surface_of_revolution()
        || surf.osculating_bspline().is_some()
}

fn locate_for_surface(
    surf: &dyn Surface,
    uv: GpPnt2d,
    duv: GpVec2d,
    left_bot: &mut GpPnt2d,
    right_top: &mut GpPnt2d,
) -> bool {
    // `EvalFirstLastSurf` switch cxx:1755-1770 / 1799-1814.
    if surf.is_offset_surface() {
        let Some(basis) = surf.offset_basis_surface() else {
            return false;
        };
        return locate_part_offset(uv, duv, basis.as_ref(), left_bot, right_top);
    }
    if surf.is_surface_of_revolution() || surf.is_surface_of_linear_extrusion() {
        return locate_part_rev_ext(uv, duv, surf, left_bot, right_top);
    }
    let Some(bs) = surf.osculating_bspline() else {
        return false;
    };
    let Some(k) = bspl_knots(&bs) else {
        return false;
    };
    locate_part(uv, duv, &k, left_bot, right_top);
    true
}

fn trim_uv(surf: &Arc<dyn Surface>, left_bot: GpPnt2d, right_top: GpPnt2d) -> Arc<dyn Surface> {
    // `mySurface->UTrim` then `VTrim` cxx:1778-1779 / 1822-1823.
    Arc::new(GeomRectangularTrimmedSurface::uv(
        surf.clone(),
        left_bot.x(),
        right_top.x(),
        left_bot.y(),
        right_top.y(),
    ))
}

/// `Adaptor3d_CurveOnSurface::EvalFirstLastSurf` (`cxx:1736-1828`).
///
/// Returns `(myFirstSurf, myLastSurf)`. `None` means keep the full surface
/// (OCCT null / fallback to `mySurface`).
pub fn eval_first_last_surf(
    pcurve: &dyn Curve2d,
    surface: &Arc<dyn Surface>,
    first: f64,
    last: f64,
) -> (Option<Arc<dyn Surface>>, Option<Arc<dyn Surface>>) {
    if !needs_first_last(surface.as_ref()) {
        return (None, None);
    }

    let mut first_surf;
    let (uv, duv) = pcurve.d1(first);
    let mut ok = duv.magnitude() > TOL;
    if ok {
        let mut left_bot = GpPnt2d::new(0.0, 0.0);
        let mut right_top = GpPnt2d::new(0.0, 0.0);
        ok = locate_for_surface(surface.as_ref(), uv, duv, &mut left_bot, &mut right_top);
        if ok {
            compare_bounds(&mut left_bot, &mut right_top);
            first_surf = Some(trim_uv(surface, left_bot, right_top));
        } else {
            first_surf = Some(surface.clone());
        }
    } else {
        first_surf = Some(surface.clone());
    }

    let mut last_surf;
    let (uv, mut duv) = pcurve.d1(last);
    duv.reverse();
    let mut ok = duv.magnitude() > TOL;
    if ok {
        let mut left_bot = GpPnt2d::new(0.0, 0.0);
        let mut right_top = GpPnt2d::new(0.0, 0.0);
        ok = locate_for_surface(surface.as_ref(), uv, duv, &mut left_bot, &mut right_top);
        if ok {
            compare_bounds(&mut left_bot, &mut right_top);
            last_surf = Some(trim_uv(surface, left_bot, right_top));
        } else {
            last_surf = Some(surface.clone());
        }
    } else {
        last_surf = Some(surface.clone());
    }

    (first_surf, last_surf)
}
