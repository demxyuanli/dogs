//! Surface grid / box helpers for bean-face localization.
//!
//! Source: `IntTools_BeanFaceIntersector.cxx:2192-2543`
//! (`GetSurfaceBox`, `ComputeGridPoints`, `BuildBox`) plus
//! `BndLib_AddSurface::Add` for the non-BSpline path.

use occt_core::bnd::BndBox;
use occt_core::gp::{GpPnt, GpVec};
use occt_geom::Surface;

use crate::bean_face_sample::SurfaceRangeLocalizeData;
use crate::geom_bnd_lib_sample2d::{compute_nb_u_samples, compute_nb_v_samples, SampleSurfaceKind};

/// `BndLib_AddSurface::Add` for a general adaptor: sample the patch and
/// enlarge by `tol`. Sample counts follow `GeomBndLib` Other = 33, capped at 50.
pub(crate) fn add_surface_to_box(
    surf: &dyn Surface,
    first_u: f64,
    last_u: f64,
    first_v: f64,
    last_v: f64,
    tol: f64,
    box_: &mut BndBox,
) {
    let (su0, su1) = surf.u_range();
    let (sv0, sv1) = surf.v_range();
    let nu = compute_nb_u_samples(
        SampleSurfaceKind::Other,
        0,
        0,
        0,
        su0,
        su1,
        first_u,
        last_u,
    )
    .max(2) as usize;
    let nv = compute_nb_v_samples(
        SampleSurfaceKind::Other,
        0,
        0,
        0,
        sv0,
        sv1,
        first_v,
        last_v,
    )
    .max(2) as usize;
    let du = last_u - first_u;
    let dv = last_v - first_v;
    for i in 0..=nu {
        let u = first_u + du * i as f64 / nu as f64;
        for j in 0..=nv {
            let v = first_v + dv * j as f64 / nv as f64;
            box_.add_point(&surf.d0(u, v));
        }
    }
    box_.enlarge(tol);
}

/// `GetSurfaceBox`.
pub(crate) fn get_surface_box(
    surf: &dyn Surface,
    first_u: f64,
    last_u: f64,
    first_v: f64,
    last_v: f64,
    tolerance: f64,
    surface_data: &mut SurfaceRangeLocalizeData,
) -> BndBox {
    let mut total = BndBox::new();
    build_box(
        surf,
        first_u,
        last_u,
        first_v,
        last_v,
        surface_data,
        &mut total,
    );
    total.enlarge(tolerance);
    total
}

/// `ComputeGridPoints` — requires real BSpline knot vectors. Empty knots
/// return immediately (OCCT `iMax < iMin`). Do not invent knots.
#[allow(dead_code)]
pub(crate) fn compute_grid_points(
    surf: &dyn Surface,
    knots_u: &[f64],
    knots_v: &[f64],
    deg_u: i32,
    deg_v: i32,
    first_u: f64,
    last_u: f64,
    first_v: f64,
    last_v: f64,
    tolerance: f64,
    surface_data: &mut SurfaceRangeLocalizeData,
) {
    let a_nb_samples = [deg_u, deg_v];
    let a_nb_knots = [knots_u.len() as i32, knots_v.len() as i32];
    if a_nb_knots[0] < 1 || a_nb_knots[1] < 1 {
        return;
    }
    let a_f_par = [first_u, first_v];
    let a_l_par = [last_u, last_v];
    let a_fp_tol = [a_f_par[0] + tolerance, a_f_par[1] + tolerance];
    let a_fm_tol = [a_f_par[0] - tolerance, a_f_par[1] - tolerance];
    let a_lp_tol = [a_l_par[0] + tolerance, a_l_par[1] + tolerance];
    let a_lm_tol = [a_l_par[0] - tolerance, a_l_par[1] - tolerance];
    let mut i_min = [-1i32, -1];
    let mut i_max = [-1i32, -1];
    let mut a_nb_grid = [0i32, 0];

    for j in 0..2 {
        let a_knots: &[f64] = if j == 0 { knots_u } else { knots_v };
        let n = a_nb_knots[j];
        let knot = |i: i32| a_knots[(i - 1) as usize];
        let mut i = 1;
        while i <= n && (i_min[j] == -1 || i_max[j] == -1) {
            if i_min[j] == -1 && a_fp_tol[j] < knot(i) {
                i_min[j] = i - 1;
            }
            let i_lmi = n - i + 1;
            if i_max[j] == -1 && a_lm_tol[j] > knot(i_lmi) {
                i_max[j] = i_lmi + 1;
            }
            i += 1;
        }
        if i_min[j] == -1 {
            i_min[j] = 1;
        }
        if i_max[j] == -1 {
            i_max[j] = n;
        }
        if i_min[j] == 0 {
            i_min[j] = 1;
        }
        if i_max[j] > n {
            i_max[j] = n;
        }
        if i_max[j] < i_min[j] {
            return;
        }
        if i_max[j] == i_min[j] {
            i_max[j] += 1;
            i_min[j] -= 1;
            if i_min[j] == 0 {
                i_min[j] = 1;
            }
            if i_max[j] > n {
                i_max[j] = n;
            }
        }
        a_nb_grid[j] = (i_max[j] - i_min[j]) * a_nb_samples[j] + 1;
        if j == 0 {
            surface_data.set_range_u_grid(a_nb_grid[j]);
        } else {
            surface_data.set_range_v_grid(a_nb_grid[j]);
        }
        let mut i_abs = 1i32;
        let mut a_max_par = if j == 0 { last_u } else { last_v };
        for i in i_min[j]..i_max[j] {
            let a_min_par = if i == i_min[j] {
                if a_fm_tol[j] > knot(i_min[j]) {
                    a_f_par[j]
                } else {
                    knot(i_min[j])
                }
            } else {
                knot(i)
            };
            a_max_par = if i == i_max[j] - 1 {
                if a_lp_tol[j] < knot(i_max[j]) {
                    a_l_par[j]
                } else {
                    knot(i_max[j])
                }
            } else {
                knot(i + 1)
            };
            let a_delta = (a_max_par - a_min_par) / a_nb_samples[j] as f64;
            let mut a_par = a_min_par;
            for _k in 0..a_nb_samples[j] {
                if j == 0 {
                    surface_data.set_u_param(i_abs, a_par);
                } else {
                    surface_data.set_v_param(i_abs, a_par);
                }
                i_abs += 1;
                a_par += a_delta;
            }
        }
        if j == 0 {
            surface_data.set_u_param(i_abs, a_max_par);
        } else {
            surface_data.set_v_param(i_abs, a_max_par);
        }
    }

    let is_calc_defl = a_nb_grid[0] < 30 && a_nb_grid[1] < 30;
    let mut a_grid_box = BndBox::new();
    let mut an_ext_box = BndBox::new();
    for i in 1..=a_nb_grid[0] {
        let a_par_u = surface_data.u_param(i);
        let du = if is_calc_defl && i < a_nb_grid[0] {
            0.5 * (surface_data.u_param(i + 1) - a_par_u)
        } else {
            0.0
        };
        for j in 1..=a_nb_grid[1] {
            let a_par_v = surface_data.v_param(j);
            let (a_pnt, a_du, a_dv) = if is_calc_defl {
                surf.d1(a_par_u, a_par_v)
            } else {
                (surf.d0(a_par_u, a_par_v), GpVec::zero(), GpVec::zero())
            };
            surface_data.set_grid_point(i, j, a_pnt);
            if is_calc_defl {
                a_grid_box.add_point(&a_pnt);
                if i < a_nb_grid[0] && j < a_nb_grid[1] {
                    let dv = 0.5 * (surface_data.v_param(j + 1) - a_par_v);
                    let a_shift = a_du
                        .multiplied_scalar(du)
                        .added(&a_dv.multiplied_scalar(dv));
                    an_ext_box.add_point(&a_pnt.translated_vec(&a_shift));
                }
            }
        }
    }

    let mut a_def: f64 = 0.0;
    if is_calc_defl {
        if let (Some((xmin, ymin, zmin, xmax, ymax, zmax)), Some((xmin1, ymin1, zmin1, xmax1, ymax1, zmax1))) =
            (a_grid_box.get(), an_ext_box.get())
        {
            let mut ext_count = 0i32;
            if xmin1 < xmin {
                a_def = a_def.max(xmin - xmin1);
                ext_count += 1;
            }
            if ymin1 < ymin {
                a_def = a_def.max(ymin - ymin1);
                ext_count += 1;
            }
            if zmin1 < zmin {
                a_def = a_def.max(zmin - zmin1);
                ext_count += 1;
            }
            if xmax1 > xmax {
                a_def = a_def.max(xmax1 - xmax);
                ext_count += 1;
            }
            if ymax1 > ymax {
                a_def = a_def.max(ymax1 - ymax);
                ext_count += 1;
            }
            if zmax1 > zmax {
                a_def = a_def.max(zmax1 - zmax);
                ext_count += 1;
            }
            if ext_count < 3 {
                a_def /= 2.0;
            }
        }
    }
    if tolerance > a_def {
        a_def = 2.0 * tolerance;
    }
    surface_data.set_grid_deflection(a_def);
}

/// `BuildBox`.
pub(crate) fn build_box(
    surf: &dyn Surface,
    first_u: f64,
    last_u: f64,
    first_v: f64,
    last_v: f64,
    surface_data: &mut SurfaceRangeLocalizeData,
    the_box: &mut BndBox,
) {
    surface_data.set_frame(first_u, last_u, first_v, last_v);
    let nb_u = surface_data.nb_u_points_in_frame();
    let nb_v = surface_data.nb_v_points_in_frame();
    the_box.add_point(&surf.d0(first_u, first_v));
    the_box.add_point(&surf.d0(last_u, first_v));
    the_box.add_point(&surf.d0(first_u, last_v));
    the_box.add_point(&surf.d0(last_u, last_v));
    for i in 1..=nb_u {
        let a_param = surface_data.u_param_in_frame(i);
        the_box.add_point(&surf.d0(a_param, first_v));
        the_box.add_point(&surf.d0(a_param, last_v));
        for j in 1..=nb_v {
            the_box.add_point(&surface_data.point_in_frame(i, j));
        }
    }
    for j in 1..=nb_v {
        let a_param = surface_data.v_param_in_frame(j);
        the_box.add_point(&surf.d0(first_u, a_param));
        the_box.add_point(&surf.d0(last_u, a_param));
    }
    the_box.enlarge(surface_data.grid_deflection());
}
