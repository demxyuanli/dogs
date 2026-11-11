//! `Geom_OffsetSurfaceUtils::EvaluateD3` (`pxx:1427-1617`) and
//! `EvaluateDN` (`pxx:1637-1794`).

use occt_core::cslib::{dn_normal, normal_d1_mag, normal_max_order, CSLibNormalStatus};
use occt_core::gp::{GpPnt, GpVec};

use crate::offset_surface_utils::{
    alloc_grid, compute_derivatives, grid_set, is_infinite_coord, shift_point, D1_TOL,
};
use crate::osculating_surface::OsculatingSurface;
use crate::surface::Surface;

struct OscPick {
    along_u: bool,
    along_v: bool,
    is_opposite: bool,
    osc: Option<crate::bspline_surface::GeomBSplineSurface>,
}

fn osc_pick(query: Option<&OsculatingSurface>, u: f64, v: f64) -> OscPick {
    let mut along_u = false;
    let mut along_v = false;
    let mut is_opposite = false;
    let mut osc = None;
    if let Some(q) = query {
        let (au, tu, lu) = q.u_osculating(u, v);
        along_u = au;
        if au {
            is_opposite = tu;
            osc = lu;
        }
        let (av, tv, lv) = q.v_osculating(u, v);
        along_v = av;
        if av {
            is_opposite = tv;
            osc = lv;
        }
    }
    OscPick {
        along_u,
        along_v,
        is_opposite,
        osc,
    }
}

fn osc_sign(pick: &OscPick) -> f64 {
    if (pick.along_u || pick.along_v) && pick.is_opposite {
        -1.0
    } else {
        1.0
    }
}

/// Result of `Geom_OffsetSurfaceUtils::EvaluateD2`
/// (`Geom_OffsetSurfaceUtils.pxx:1370-1404`).
pub struct OffsetSecond {
    pub point: GpPnt,
    pub d1u: GpVec,
    pub d1v: GpVec,
    pub d2u: GpVec,
    pub d2v: GpVec,
    pub d2uv: GpVec,
}

/// `EvaluateD2` with the same singular-point retry as the C++ wrapper.
pub fn evaluate_d2(
    u_in: f64,
    v_in: f64,
    basis: &dyn Surface,
    offset: f64,
    osc_query: Option<&OsculatingSurface>,
) -> Option<OffsetSecond> {
    let u_start = u_in;
    let v_start = v_in;
    let (u_min, u_max) = basis.u_range();
    let (v_min, v_max) = basis.v_range();
    let is_u_per = basis.is_u_periodic();
    let is_v_per = basis.is_v_periodic();
    let mut u = u_in;
    let mut v = v_in;

    loop {
        let jet = basis.d3(u, v)?;
        if is_infinite_coord(&jet.d1u) || is_infinite_coord(&jet.d1v) {
            return None;
        }

        let (n_status, _) = normal_d1_mag(&jet.d1u, &jet.d1v, D1_TOL);
        let max_order = if n_status == CSLibNormalStatus::Defined {
            0
        } else {
            3
        };
        let pick = if n_status != CSLibNormalStatus::Defined {
            osc_pick(osc_query, u, v)
        } else {
            OscPick {
                along_u: false,
                along_v: false,
                is_opposite: false,
                osc: None,
            }
        };

        let mut der_nuv = alloc_grid(max_order + 2, max_order + 2);
        let mut der_surf = alloc_grid(max_order + 3, max_order + 3);
        grid_set(&mut der_surf, 1, 0, jet.d1u);
        grid_set(&mut der_surf, 0, 1, jet.d1v);
        grid_set(&mut der_surf, 1, 1, jet.d2uv);
        grid_set(&mut der_surf, 2, 0, jet.d2u);
        grid_set(&mut der_surf, 0, 2, jet.d2v);
        grid_set(&mut der_surf, 3, 0, jet.d3u);
        grid_set(&mut der_surf, 2, 1, jet.d3_uuv);
        grid_set(&mut der_surf, 1, 2, jet.d3_uvv);
        grid_set(&mut der_surf, 0, 3, jet.d3v);

        let osc_ref = pick.osc.as_ref().map(|s| s as &dyn Surface);
        let has_osc = (pick.along_u || pick.along_v) && osc_ref.is_some();
        if !compute_derivatives(
            max_order,
            3,
            u,
            v,
            basis,
            2,
            2,
            has_osc && pick.along_u,
            has_osc && pick.along_v,
            osc_ref,
            &mut der_nuv,
            &mut der_surf,
        ) {
            return None;
        }

        let (status, _, order_u, order_v) = normal_max_order(
            max_order, &der_nuv, D1_TOL, u, v, u_min, u_max, v_min, v_max,
        );
        if status == CSLibNormalStatus::Defined {
            let sign = offset * osc_sign(&pick);
            let normal = dn_normal(0, 0, &der_nuv, order_u, order_v);
            let point = GpPnt::from_xyz(&jet.point.coord.added(&normal.xyz().multiplied(sign)));
            let add = |nu: i32, nv: i32| {
                basis
                    .eval_dn(u, v, nu, nv)
                    .added(&dn_normal(nu, nv, &der_nuv, order_u, order_v).multiplied_scalar(sign))
            };
            return Some(OffsetSecond {
                point,
                d1u: jet
                    .d1u
                    .added(&dn_normal(1, 0, &der_nuv, order_u, order_v).multiplied_scalar(sign)),
                d1v: jet
                    .d1v
                    .added(&dn_normal(0, 1, &der_nuv, order_u, order_v).multiplied_scalar(sign)),
                d2u: add(2, 0),
                d2v: add(0, 2),
                d2uv: add(1, 1),
            });
        }

        if !shift_point(
            u_start, v_start, &mut u, &mut v, u_min, u_max, v_min, v_max, is_u_per, is_v_per,
            &jet.d1u, &jet.d1v,
        ) {
            return None;
        }
    }
}

/// Third-order offset derivatives from `EvaluateD3`.
pub struct OffsetThird {
    pub d3u: GpVec,
    pub d3v: GpVec,
    pub d3uuv: GpVec,
    pub d3uvv: GpVec,
}

/// `EvaluateD3` (`pxx:1427-1617`). `None` is the `false` return.
pub fn evaluate_d3(
    u_in: f64,
    v_in: f64,
    basis: &dyn Surface,
    offset: f64,
    osc_query: Option<&OsculatingSurface>,
) -> Option<OffsetThird> {
    let u_start = u_in;
    let v_start = v_in;
    let (u_min, u_max) = basis.u_range();
    let (v_min, v_max) = basis.v_range();
    let is_u_per = basis.is_u_periodic();
    let is_v_per = basis.is_v_periodic();
    let mut u = u_in;
    let mut v = v_in;

    loop {
        let (_, d1u, d1v, d2u, d2v, d2uv) = basis.d2(u, v);
        let d3u = basis.eval_dn(u, v, 3, 0);
        let d3v = basis.eval_dn(u, v, 0, 3);
        let d3uuv = basis.eval_dn(u, v, 2, 1);
        let d3uvv = basis.eval_dn(u, v, 1, 2);
        if is_infinite_coord(&d1u) || is_infinite_coord(&d1v) {
            return None;
        }

        let (n_status, _) = normal_d1_mag(&d1u, &d1v, D1_TOL);
        let max_order = if n_status == CSLibNormalStatus::Defined {
            0
        } else {
            3
        };
        let pick = if n_status != CSLibNormalStatus::Defined {
            osc_pick(osc_query, u, v)
        } else {
            OscPick {
                along_u: false,
                along_v: false,
                is_opposite: false,
                osc: None,
            }
        };

        let mut der_nuv = alloc_grid(max_order + 3, max_order + 3);
        let mut der_surf = alloc_grid(max_order + 4, max_order + 4);
        grid_set(&mut der_surf, 1, 0, d1u);
        grid_set(&mut der_surf, 0, 1, d1v);
        grid_set(&mut der_surf, 1, 1, d2uv);
        grid_set(&mut der_surf, 2, 0, d2u);
        grid_set(&mut der_surf, 0, 2, d2v);
        grid_set(&mut der_surf, 3, 0, d3u);
        grid_set(&mut der_surf, 2, 1, d3uuv);
        grid_set(&mut der_surf, 1, 2, d3uvv);
        grid_set(&mut der_surf, 0, 3, d3v);

        let osc_ref = pick.osc.as_ref().map(|s| s as &dyn Surface);
        let has_osc = (pick.along_u || pick.along_v) && osc_ref.is_some();
        if !compute_derivatives(
            max_order,
            3,
            u,
            v,
            basis,
            3,
            3,
            has_osc && pick.along_u,
            has_osc && pick.along_v,
            osc_ref,
            &mut der_nuv,
            &mut der_surf,
        ) {
            return None;
        }

        let (status, _, order_u, order_v) = normal_max_order(
            max_order, &der_nuv, D1_TOL, u, v, u_min, u_max, v_min, v_max,
        );
        if status == CSLibNormalStatus::Defined {
            let sign = offset * osc_sign(&pick);
            let add = |nu: i32, nv: i32| {
                basis
                    .eval_dn(u, v, nu, nv)
                    .added(&dn_normal(nu, nv, &der_nuv, order_u, order_v).multiplied_scalar(sign))
            };
            return Some(OffsetThird {
                d3u: add(3, 0),
                d3v: add(0, 3),
                d3uuv: add(2, 1),
                d3uvv: add(1, 2),
            });
        }

        if !shift_point(
            u_start, v_start, &mut u, &mut v, u_min, u_max, v_min, v_max, is_u_per, is_v_per, &d1u,
            &d1v,
        ) {
            return None;
        }
    }
}

/// `EvaluateDN` (`pxx:1637-1794`, wrapper `pxx:1811-1826`). `None` is `false`.
pub fn evaluate_dn(
    u_in: f64,
    v_in: f64,
    nu: i32,
    nv: i32,
    basis: &dyn Surface,
    offset: f64,
    osc_query: Option<&OsculatingSurface>,
) -> Option<GpVec> {
    let u_start = u_in;
    let v_start = v_in;
    let (u_min, u_max) = basis.u_range();
    let (v_min, v_max) = basis.v_range();
    let is_u_per = basis.is_u_periodic();
    let is_v_per = basis.is_v_periodic();
    let mut u = u_in;
    let mut v = v_in;

    loop {
        let (_, d1u, d1v) = basis.d1(u, v);
        if is_infinite_coord(&d1u) || is_infinite_coord(&d1v) {
            return None;
        }
        let (n_status, _) = normal_d1_mag(&d1u, &d1v, D1_TOL);
        let max_order = if n_status == CSLibNormalStatus::Defined {
            0
        } else {
            3
        };
        let pick = if n_status != CSLibNormalStatus::Defined {
            osc_pick(osc_query, u, v)
        } else {
            OscPick {
                along_u: false,
                along_v: false,
                is_opposite: false,
                osc: None,
            }
        };

        let mut der_nuv = alloc_grid(max_order + nu, max_order + nv);
        let mut der_surf = alloc_grid(max_order + nu + 1, max_order + nv + 1);
        grid_set(&mut der_surf, 1, 0, d1u);
        grid_set(&mut der_surf, 0, 1, d1v);

        let osc_ref = pick.osc.as_ref().map(|s| s as &dyn Surface);
        let has_osc = (pick.along_u || pick.along_v) && osc_ref.is_some();
        if !compute_derivatives(
            max_order,
            1,
            u,
            v,
            basis,
            nu,
            nv,
            has_osc && pick.along_u,
            has_osc && pick.along_v,
            osc_ref,
            &mut der_nuv,
            &mut der_surf,
        ) {
            return None;
        }

        let (status, _, order_u, order_v) = normal_max_order(
            max_order, &der_nuv, D1_TOL, u, v, u_min, u_max, v_min, v_max,
        );
        if status == CSLibNormalStatus::Defined {
            let sign = offset * osc_sign(&pick);
            let base = basis.eval_dn(u, v, nu, nv);
            return Some(
                base.added(&dn_normal(nu, nv, &der_nuv, order_u, order_v).multiplied_scalar(sign)),
            );
        }

        if !shift_point(
            u_start, v_start, &mut u, &mut v, u_min, u_max, v_min, v_max, is_u_per, is_v_per, &d1u,
            &d1v,
        ) {
            return None;
        }
    }
}
