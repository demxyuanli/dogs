//! `GeomAdaptor_Surface` LocalD1 / IfUVBound / Span for restricted UV domains.
//!
//! Source: `GeomAdaptor_Surface.cxx:1122-1205`, `2244-2350` and
//! `Geom_BSplineSurface::LocalD1/LocalD2` (`Geom_BSplineSurface_1.cxx:372-464`).

use occt_core::bspl::locate::{first_u_knot_index, last_u_knot_index};
use occt_core::gp::{GpPnt, GpVec};
use occt_core::precision::PCONFUSION;

use crate::bspline_surface::GeomBSplineSurface;

/// `PosTol` in `GeomAdaptor_Surface.cxx:72`.
const POS_TOL: f64 = PCONFUSION * 0.5;

/// `GeomAdaptor_Surface::Span` (`cxx:2283-2350`).
fn span(
    side: i32,
    ideb: i32,
    ifin: i32,
    fk_indx: i32,
    lk_indx: i32,
) -> (i32, i32) {
    if ideb != ifin {
        if ideb < fk_indx {
            (fk_indx, fk_indx + 1)
        } else if ifin > lk_indx {
            (lk_indx - 1, lk_indx)
        } else if ideb >= (lk_indx - 1) {
            (lk_indx - 1, lk_indx)
        } else if ifin <= fk_indx + 1 {
            (fk_indx, fk_indx + 1)
        } else if ideb > ifin {
            (ifin - 1, ifin)
        } else {
            (ideb, ifin)
        }
    } else if ideb <= fk_indx {
        (fk_indx, fk_indx + 1)
    } else if ifin >= lk_indx {
        (lk_indx - 1, lk_indx)
    } else if side == -1 {
        (ideb - 1, ifin)
    } else {
        (ideb, ifin + 1)
    }
}

/// `GeomAdaptor_Surface::IfUVBound` (`cxx:2244-2274`).
/// Returns `(Local, Ideb, Ifin, IVdeb, IVfin)`.
fn if_uv_bound(
    bs: &GeomBSplineSurface,
    u: f64,
    v: f64,
    u_side: i32,
    v_side: i32,
) -> Option<(bool, i32, i32, i32, i32)> {
    let (_, u_mults) = GeomBSplineSurface::unique_knots_mults(&bs.knots_u);
    let (_, v_mults) = GeomBSplineSurface::unique_knots_mults(&bs.knots_v);
    let u_fk = first_u_knot_index(bs.deg_u as i32, &u_mults);
    let u_lk = last_u_knot_index(bs.deg_u as i32, &u_mults);
    let v_fk = first_u_knot_index(bs.deg_v as i32, &v_mults);
    let v_lk = last_u_knot_index(bs.deg_v as i32, &v_mults);
    if u_fk >= u_lk || v_fk >= v_lk {
        return None;
    }
    let (mut ideb, mut ifin) = bs.locate_u(u, POS_TOL);
    let mut local = ideb == ifin;
    let (ou, ov) = span(u_side, ideb, ifin, u_fk, u_lk);
    ideb = ou;
    ifin = ov;
    let (mut ivdeb, mut ivfin) = bs.locate_v(v, POS_TOL);
    if ivdeb == ivfin {
        local = true;
    }
    let (vu, vv) = span(v_side, ivdeb, ivfin, v_fk, v_lk);
    ivdeb = vu;
    ivfin = vv;
    Some((local, ideb, ifin, ivdeb, ivfin))
}

/// Side flags from adaptor domain ends (`GeomAdaptor_Surface::EvalD1` cxx:1129-1148).
fn uv_sides(
    u: f64,
    v: f64,
    u_first: f64,
    u_last: f64,
    v_first: f64,
    v_last: f64,
    tol_u: f64,
    tol_v: f64,
) -> (i32, i32, f64, f64) {
    let mut u_side = 0i32;
    let mut v_side = 0i32;
    let mut uu = u;
    let mut vv = v;
    if (u - u_first).abs() <= tol_u {
        u_side = 1;
        uu = u_first;
    } else if (u - u_last).abs() <= tol_u {
        u_side = -1;
        uu = u_last;
    }
    if (v - v_first).abs() <= tol_v {
        v_side = 1;
        vv = v_first;
    } else if (v - v_last).abs() <= tol_v {
        v_side = -1;
        vv = v_last;
    }
    (u_side, v_side, uu, vv)
}

/// BSpline arm of `GeomAdaptor_Surface::EvalD1` LocalD1 (`cxx:1193-1195`).
pub fn try_local_d1(
    bs: &GeomBSplineSurface,
    u: f64,
    v: f64,
    u_first: f64,
    u_last: f64,
    v_first: f64,
    v_last: f64,
    tol_u: f64,
    tol_v: f64,
) -> Option<(GpPnt, GpVec, GpVec)> {
    let (u_side, v_side, uu, vv) =
        uv_sides(u, v, u_first, u_last, v_first, v_last, tol_u, tol_v);
    if u_side == 0 && v_side == 0 {
        return None;
    }
    let (local, ideb, ifin, ivdeb, ivfin) = if_uv_bound(bs, uu, vv, u_side, v_side)?;
    if !local {
        return None;
    }
    Some(bs.local_d1(uu, vv, ideb, ifin, ivdeb, ivfin))
}

/// BSpline arm of `GeomAdaptor_Surface::EvalD2` LocalD2 (`cxx:1273-1292` + LocalD2).
pub fn try_local_d2(
    bs: &GeomBSplineSurface,
    u: f64,
    v: f64,
    u_first: f64,
    u_last: f64,
    v_first: f64,
    v_last: f64,
    tol_u: f64,
    tol_v: f64,
) -> Option<(GpPnt, GpVec, GpVec, GpVec, GpVec, GpVec)> {
    let (u_side, v_side, uu, vv) =
        uv_sides(u, v, u_first, u_last, v_first, v_last, tol_u, tol_v);
    if u_side == 0 && v_side == 0 {
        return None;
    }
    let (local, ideb, ifin, ivdeb, ivfin) = if_uv_bound(bs, uu, vv, u_side, v_side)?;
    if !local {
        return None;
    }
    Some(bs.local_d2(uu, vv, ideb, ifin, ivdeb, ivfin))
}
