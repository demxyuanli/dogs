//! `Geom_OffsetSurfaceUtils`. Source: `Geom_OffsetSurfaceUtils.pxx`.
//! `EvaluateD1` queries `OsculatingSurface` (null query keeps the else arm).

use occt_core::cslib::{dn_normal, dnnuv, dnnuv2, normal_max_order, CSLibNormalStatus};
use occt_core::gp::{GpPnt, GpVec};
use occt_core::precision::{Precision, CONFUSION, PCONFUSION};

use crate::osculating_surface::OsculatingSurface;
use crate::surface::Surface;

/// `THE_D1_MAGNITUDE_TOL` (`pxx:45`).
const D1_TOL: f64 = 1.0e-9;

fn is_infinite_coord(v: &GpVec) -> bool {
    Precision::is_infinite(v.x()) || Precision::is_infinite(v.y()) || Precision::is_infinite(v.z())
}

fn grid_get(g: &[Vec<GpVec>], i: i32, j: i32) -> GpVec {
    if i < 0 || j < 0 {
        return GpVec::new(0.0, 0.0, 0.0);
    }
    g.get(i as usize)
        .and_then(|row| row.get(j as usize))
        .copied()
        .unwrap_or(GpVec::new(0.0, 0.0, 0.0))
}

fn grid_set(g: &mut [Vec<GpVec>], i: i32, j: i32, v: GpVec) {
    if i < 0 || j < 0 {
        return;
    }
    if let Some(slot) = g.get_mut(i as usize).and_then(|row| row.get_mut(j as usize)) {
        *slot = v;
    }
}

fn alloc_grid(u_upper: i32, v_upper: i32) -> Vec<Vec<GpVec>> {
    let rows = (u_upper + 1).max(0) as usize;
    let cols = (v_upper + 1).max(0) as usize;
    vec![vec![GpVec::new(0.0, 0.0, 0.0); cols]; rows]
}

fn upper_row(g: &[Vec<GpVec>]) -> i32 {
    g.len() as i32 - 1
}

fn upper_col(g: &[Vec<GpVec>]) -> i32 {
    g.first().map(|r| r.len() as i32 - 1).unwrap_or(-1)
}

/// `Geom_OffsetSurfaceUtils::ComputeDerivatives` (`pxx:267-414`).
pub fn compute_derivatives(
    mut max_order: i32,
    min_order: i32,
    u: f64,
    v: f64,
    basis: &dyn Surface,
    nu: i32,
    nv: i32,
    along_u: bool,
    along_v: bool,
    osc: Option<&dyn Surface>,
    der_nuv: &mut [Vec<GpVec>],
    der_surf: &mut [Vec<GpVec>],
) -> bool {
    if along_u || along_v {
        max_order = 0;
        let mut der_surf_l = alloc_grid(max_order + nu + 1, max_order + nv + 1);
        if let Some(osc) = osc {
            match min_order {
                1 => {
                    let (_, d1u, d1v) = osc.d1(u, v);
                    grid_set(&mut der_surf_l, 1, 0, d1u);
                    grid_set(&mut der_surf_l, 0, 1, d1v);
                }
                2 => {
                    let (_, d1u, d1v, d2u, d2v, d2uv) = osc.d2(u, v);
                    grid_set(&mut der_surf_l, 1, 0, d1u);
                    grid_set(&mut der_surf_l, 0, 1, d1v);
                    grid_set(&mut der_surf_l, 1, 1, d2uv);
                    grid_set(&mut der_surf_l, 2, 0, d2u);
                    grid_set(&mut der_surf_l, 0, 2, d2v);
                }
                3 => {
                    let (_, d1u, d1v, d2u, d2v, d2uv) = osc.d2(u, v);
                    grid_set(&mut der_surf_l, 1, 0, d1u);
                    grid_set(&mut der_surf_l, 0, 1, d1v);
                    grid_set(&mut der_surf_l, 1, 1, d2uv);
                    grid_set(&mut der_surf_l, 2, 0, d2u);
                    grid_set(&mut der_surf_l, 0, 2, d2v);
                    grid_set(&mut der_surf_l, 3, 0, osc.eval_dn(u, v, 3, 0));
                    grid_set(&mut der_surf_l, 2, 1, osc.eval_dn(u, v, 2, 1));
                    grid_set(&mut der_surf_l, 1, 2, osc.eval_dn(u, v, 1, 2));
                    grid_set(&mut der_surf_l, 0, 3, osc.eval_dn(u, v, 0, 3));
                }
                _ => {}
            }
        }
        if nu <= nv {
            let mut i = 0;
            while i <= max_order + 1 + nu {
                let mut j = i;
                while j <= max_order + nv + 1 {
                    if i + j > min_order {
                        if let Some(osc) = osc {
                            grid_set(&mut der_surf_l, i, j, osc.eval_dn(u, v, i, j));
                        }
                        grid_set(der_surf, i, j, basis.eval_dn(u, v, i, j));
                        if i != j && j <= nu + 1 {
                            grid_set(der_surf, j, i, basis.eval_dn(u, v, j, i));
                            if let Some(osc) = osc {
                                grid_set(&mut der_surf_l, j, i, osc.eval_dn(u, v, j, i));
                            }
                        }
                    }
                    j += 1;
                }
                i += 1;
            }
        } else {
            let mut j = 0;
            while j <= max_order + 1 + nv {
                let mut i = j;
                while i <= max_order + nu + 1 {
                    if i + j > min_order {
                        if let Some(osc) = osc {
                            grid_set(&mut der_surf_l, i, j, osc.eval_dn(u, v, i, j));
                        }
                        grid_set(der_surf, i, j, basis.eval_dn(u, v, i, j));
                        if i != j && i <= nv + 1 {
                            grid_set(der_surf, j, i, basis.eval_dn(u, v, j, i));
                            if let Some(osc) = osc {
                                grid_set(&mut der_surf_l, j, i, osc.eval_dn(u, v, j, i));
                            }
                        }
                    }
                    i += 1;
                }
                j += 1;
            }
        }
        let mut i = 0;
        while i <= max_order + nu {
            let mut j = 0;
            while j <= max_order + nv {
                if along_u {
                    grid_set(der_nuv, i, j, dnnuv2(i, j, &der_surf_l, der_surf));
                }
                if along_v {
                    grid_set(der_nuv, i, j, dnnuv2(i, j, der_surf, &der_surf_l));
                }
                j += 1;
            }
            i += 1;
        }
    } else {
        let mut i = 0;
        while i <= max_order + nu + 1 {
            let mut j = i;
            while j <= max_order + nv + 1 {
                if i + j > min_order {
                    grid_set(der_surf, i, j, basis.eval_dn(u, v, i, j));
                    if i != j && j <= upper_row(der_surf) && i <= upper_col(der_surf) {
                        grid_set(der_surf, j, i, basis.eval_dn(u, v, j, i));
                    }
                }
                j += 1;
            }
            i += 1;
        }
        let mut i = 0;
        while i <= max_order + nu {
            let mut j = 0;
            while j <= max_order + nv {
                grid_set(der_nuv, i, j, dnnuv(i, j, der_surf));
                j += 1;
            }
            i += 1;
        }
    }
    true
}

/// `Geom_OffsetSurfaceUtils::ReplaceDerivative` (`pxx:432-508`).
pub fn replace_derivative(
    u: f64,
    v: f64,
    u_min: f64,
    u_max: f64,
    v_min: f64,
    v_max: f64,
    du: &mut GpVec,
    dv: &mut GpVec,
    square_tol: f64,
    basis: &dyn Surface,
) -> bool {
    let is_replace_du = du.square_magnitude() < square_tol;
    let is_replace_dv = dv.square_magnitude() < square_tol;
    let mut is_replaced = false;
    if is_replace_du != is_replace_dv {
        let a_step = if is_replace_dv {
            let mut s = CONFUSION * du.magnitude();
            if s > u_max - u_min {
                s = (u_max - u_min) / 100.0;
            }
            s
        } else {
            let mut s = CONFUSION * dv.magnitude();
            if s > v_max - v_min {
                s = (v_max - v_min) / 100.0;
            }
            s
        };
        let mut step_sign = -1.0;
        while step_sign <= 1.0 && !is_replaced {
            let mut au = u;
            let mut av = v;
            if is_replace_dv {
                au = u + step_sign * a_step;
                if au < u_min || au > u_max {
                    step_sign += 2.0;
                    continue;
                }
            } else {
                av = v + step_sign * a_step;
                if av < v_min || av > v_max {
                    step_sign += 2.0;
                    continue;
                }
            }
            let (_, d1u, d1v) = basis.d1(au, av);
            if is_replace_du && d1u.square_magnitude() > square_tol {
                *du = d1u;
                is_replaced = true;
            }
            if is_replace_dv && d1v.square_magnitude() > square_tol {
                *dv = d1v;
                is_replaced = true;
            }
            step_sign += 2.0;
        }
    }
    is_replaced
}

/// `Geom_OffsetSurfaceUtils::ShiftPoint` (`pxx:527-568`).
pub fn shift_point(
    u_start: f64,
    v_start: f64,
    u: &mut f64,
    v: &mut f64,
    u_min: f64,
    u_max: f64,
    v_min: f64,
    v_max: f64,
    is_u_periodic: bool,
    is_v_periodic: bool,
    d1u: &GpVec,
    d1v: &GpVec,
) -> bool {
    let is_u_singular = d1u.square_magnitude() < D1_TOL * D1_TOL;
    let is_v_singular = d1v.square_magnitude() < D1_TOL * D1_TOL;
    let dir_u = if is_u_periodic || (is_u_singular && !is_v_singular) {
        0.0
    } else {
        0.5 * (u_min + u_max) - u_start
    };
    let dir_v = if is_v_periodic || (is_v_singular && !is_u_singular) {
        0.0
    } else {
        0.5 * (v_min + v_max) - v_start
    };
    let dist = (dir_u * dir_u + dir_v * dir_v).sqrt();
    let du = *u - u_start;
    let dv = *v - v_start;
    let mut step = (2.0 * (du * du + dv * dv).sqrt()).max(PCONFUSION);
    if step >= dist {
        return false;
    }
    step /= dist;
    *u += dir_u * step;
    *v += dir_v * step;
    true
}

fn fill_der_surf_d2(
    der_surf: &mut [Vec<GpVec>],
    d1u: GpVec,
    d1v: GpVec,
    d2u: GpVec,
    d2v: GpVec,
    d2uv: GpVec,
) {
    grid_set(der_surf, 1, 0, d1u);
    grid_set(der_surf, 0, 1, d1v);
    grid_set(der_surf, 1, 1, d2uv);
    grid_set(der_surf, 2, 0, d2u);
    grid_set(der_surf, 0, 2, d2v);
}

fn run_compute(
    max_order: i32,
    u: f64,
    v: f64,
    basis: &dyn Surface,
    along_u: bool,
    along_v: bool,
    osc: Option<&dyn Surface>,
    der_nuv: &mut [Vec<GpVec>],
    der_surf: &mut [Vec<GpVec>],
) -> bool {
    compute_derivatives(
        max_order, 2, u, v, basis, 1, 1, along_u, along_v, osc, der_nuv, der_surf,
    )
}

/// `EvaluateD1` with precomputed basis D2 (`pxx:804-1095`).
/// Returns `None` when cxx throws `Geom_UndefinedDerivative`.
pub fn evaluate_d1(
    u_in: f64,
    v_in: f64,
    basis: &dyn Surface,
    offset: f64,
    osc_query: Option<&OsculatingSurface>,
    value_in: GpPnt,
    d1u_in: GpVec,
    d1v_in: GpVec,
    d2u_in: GpVec,
    d2v_in: GpVec,
    d2uv_in: GpVec,
) -> Option<(GpPnt, GpVec, GpVec)> {
    let u_start = u_in;
    let v_start = v_in;
    let (u_min, u_max) = basis.u_range();
    let (v_min, v_max) = basis.v_range();
    let is_u_per = basis.is_u_periodic();
    let is_v_per = basis.is_v_periodic();

    let mut the_u = u_in;
    let mut the_v = v_in;
    let mut the_value = value_in;
    let mut the_d1u = d1u_in;
    let mut the_d1v = d1v_in;
    let mut a_d2u = d2u_in;
    let mut a_d2v = d2v_in;
    let mut a_d2uv = d2uv_in;
    let mut is_first = true;

    loop {
        if !is_first {
            let (p, d1u, d1v, d2u, d2v, d2uv) = basis.d2(the_u, the_v);
            the_value = p;
            the_d1u = d1u;
            the_d1v = d1v;
            a_d2u = d2u;
            a_d2v = d2v;
            a_d2uv = d2uv;
        }
        is_first = false;

        if is_infinite_coord(&the_d1u) || is_infinite_coord(&the_d1v) {
            return None;
        }

        let mut a_d1u = the_d1u;
        let mut a_d1v = the_d1v;
        let nu2 = a_d1u.square_magnitude();
        let nv2 = a_d1v.square_magnitude();
        if nu2 > 1.0 {
            a_d1u = a_d1u.divided(nu2.sqrt());
        }
        if nv2 > 1.0 {
            a_d1v = a_d1v.divided(nv2.sqrt());
        }

        let max_order = 3;
        let mut a_norm = a_d1u.crossed(&a_d1v);
        let is_singular = a_norm.square_magnitude() <= D1_TOL * D1_TOL;
        let mut along_u = false;
        let mut along_v = false;
        let mut is_opposite = false;
        let mut osc_surf: Option<crate::bspline_surface::GeomBSplineSurface> = None;
        if is_singular {
            if let Some(q) = osc_query {
                let (au, tu, lu) = q.u_osculating(the_u, the_v);
                along_u = au;
                if au {
                    is_opposite = tu;
                    osc_surf = lu;
                }
                let (av, tv, lv) = q.v_osculating(the_u, the_v);
                along_v = av;
                if av {
                    is_opposite = tv;
                    osc_surf = lv;
                }
            }
        }
        let a_sign = if (along_u || along_v) && is_opposite {
            -1.0
        } else {
            1.0
        };

        if !is_singular {
            a_norm.normalize();
            the_value = GpPnt::from_xyz(&the_value.coord.added(&a_norm.xyz().multiplied(offset * a_sign)));
            let a_n0 = a_norm;
            let a_scale = the_d1u.crossed(&the_d1v).dot(&a_n0);
            let mut a_n1u = GpVec::new(
                a_d2u.y() * the_d1v.z() + the_d1u.y() * a_d2uv.z()
                    - a_d2u.z() * the_d1v.y()
                    - the_d1u.z() * a_d2uv.y(),
                -(a_d2u.x() * the_d1v.z() + the_d1u.x() * a_d2uv.z()
                    - a_d2u.z() * the_d1v.x()
                    - the_d1u.z() * a_d2uv.x()),
                a_d2u.x() * the_d1v.y() + the_d1u.x() * a_d2uv.y()
                    - a_d2u.y() * the_d1v.x()
                    - the_d1u.y() * a_d2uv.x(),
            );
            let scale_u = a_n1u.dot(&a_n0);
            a_n1u = a_n1u
                .subtracted(&a_n0.multiplied_scalar(scale_u))
                .divided(a_scale);
            let mut a_n1v = GpVec::new(
                a_d2uv.y() * the_d1v.z() + a_d2v.z() * the_d1u.y()
                    - a_d2uv.z() * the_d1v.y()
                    - a_d2v.y() * the_d1u.z(),
                -(a_d2uv.x() * the_d1v.z() + a_d2v.z() * the_d1u.x()
                    - a_d2uv.z() * the_d1v.x()
                    - a_d2v.x() * the_d1u.z()),
                a_d2uv.x() * the_d1v.y() + a_d2v.y() * the_d1u.x()
                    - a_d2uv.y() * the_d1v.x()
                    - a_d2v.x() * the_d1u.y(),
            );
            let scale_v = a_n1v.dot(&a_n0);
            a_n1v = a_n1v
                .subtracted(&a_n0.multiplied_scalar(scale_v))
                .divided(a_scale);
            the_d1u = the_d1u.added(&a_n1u.multiplied_scalar(offset * a_sign));
            the_d1v = the_d1v.added(&a_n1v.multiplied_scalar(offset * a_sign));
            return Some((the_value, the_d1u, the_d1v));
        }

        let mut der_nuv = alloc_grid(max_order + 1, max_order + 1);
        let mut der_surf = alloc_grid(max_order + 2, max_order + 2);
        fill_der_surf_d2(&mut der_surf, the_d1u, the_d1v, a_d2u, a_d2v, a_d2uv);
        let osc_ref = osc_surf.as_ref().map(|s| s as &dyn Surface);
        let has_osc = (along_u || along_v) && osc_ref.is_some();
        if !run_compute(
            max_order,
            the_u,
            the_v,
            basis,
            has_osc && along_u,
            has_osc && along_v,
            osc_ref,
            &mut der_nuv,
            &mut der_surf,
        ) {
            return None;
        }

        let (mut n_status, mut normal, mut order_u, mut order_v) = normal_max_order(
            max_order, &der_nuv, D1_TOL, the_u, the_v, u_min, u_max, v_min, v_max,
        );

        if n_status == CSLibNormalStatus::InfinityOfSolutions {
            let mut new_du = the_d1u;
            let mut new_dv = the_d1v;
            if replace_derivative(
                the_u,
                the_v,
                u_min,
                u_max,
                v_min,
                v_max,
                &mut new_du,
                &mut new_dv,
                D1_TOL * D1_TOL,
                basis,
            ) {
                grid_set(&mut der_surf, 1, 0, new_du);
                grid_set(&mut der_surf, 0, 1, new_dv);
                let osc_ref = osc_surf.as_ref().map(|s| s as &dyn Surface);
                if !run_compute(
                    max_order,
                    the_u,
                    the_v,
                    basis,
                    has_osc && along_u,
                    has_osc && along_v,
                    osc_ref,
                    &mut der_nuv,
                    &mut der_surf,
                ) {
                    return None;
                }
                let r = normal_max_order(
                    max_order, &der_nuv, D1_TOL, the_u, the_v, u_min, u_max, v_min, v_max,
                );
                n_status = r.0;
                normal = r.1;
                order_u = r.2;
                order_v = r.3;
            }
        }

        if n_status == CSLibNormalStatus::Defined {
            if let Some(n) = normal {
                the_value =
                    GpPnt::from_xyz(&the_value.coord.added(&n.xyz().multiplied(offset * a_sign)));
                let mut out_d1u = dn_normal(1, 0, &der_nuv, order_u, order_v);
                let mut out_d1v = dn_normal(0, 1, &der_nuv, order_u, order_v);
                out_d1u = out_d1u.multiplied_scalar(offset * a_sign);
                out_d1u = out_d1u.added(&grid_get(&der_surf, 1, 0));
                out_d1v = out_d1v.multiplied_scalar(offset * a_sign);
                out_d1v = out_d1v.added(&grid_get(&der_surf, 0, 1));
                return Some((the_value, out_d1u, out_d1v));
            }
        }

        if !shift_point(
            u_start,
            v_start,
            &mut the_u,
            &mut the_v,
            u_min,
            u_max,
            v_min,
            v_max,
            is_u_per,
            is_v_per,
            &the_d1u,
            &the_d1v,
        ) {
            return None;
        }
    }
}
