//! `GeomConvert::SurfaceToBSplineSurface` elementary arms
//! (`GeomConvert_1.cxx:433-543`) through
//! `Convert_{Cylinder,Cone,Sphere,Torus}ToBSplineSurface`.
//!
//! A full-period seam that does not sit on the first period calls
//! `Geom_BSplineSurface::CheckAndSegment` (`Geom_BSplineSurface.cxx:880`).

use std::f64::consts::PI;

use occt_core::gp::{GpAx3, GpPnt, GpTrsf};
use occt_core::precision::{epsilon, PCONFUSION};
use occt_geom::bspline_surface::GeomBSplineSurface;
use occt_geom::Surface;

use super::gtransform::flat_knots;

/// Cylinder, cone, sphere, or torus trimmed to `[u1, u2] x [v1, v2]`.
pub(super) fn elementary_patch(
    surf: &dyn Surface,
    u1: f64,
    u2: f64,
    v1: f64,
    v2: f64,
) -> Result<Option<GeomBSplineSurface>, String> {
    let u_closed = (u2 - u1 - 2.0 * PI).abs() <= PCONFUSION;
    if let Some(cyl) = surf.gp_cylinder() {
        let built = if u_closed {
            cylinder_periodic(&cyl, v1, v2)?
        } else {
            cylinder_open(&cyl, u1, u2, v1, v2)?
        };
        return Ok(Some(shift_u(built, u_closed, u1, u2, v1, v2)?));
    }
    if let Some(cone) = surf.gp_cone() {
        let built = if u_closed {
            cone_periodic(&cone, v1, v2)?
        } else {
            cone_open(&cone, u1, u2, v1, v2)?
        };
        return Ok(Some(shift_u(built, u_closed, u1, u2, v1, v2)?));
    }
    if let Some(sph) = surf.gp_sphere() {
        let built = if u_closed {
            // `Convert_SphereToBSplineSurface(Sph, VFirst, VLast, false)`.
            sphere_param(&sph, v1, v2, false)?
        } else {
            sphere_open(&sph, u1, u2, v1, v2)?
        };
        return Ok(Some(shift_u(built, u_closed, u1, u2, v1, v2)?));
    }
    if let Some(tor) = surf.gp_torus() {
        let u_first_period = u1 >= 0.0 && u2 <= 2.0 * PI;
        let v_first_period = v1 >= 0.0 && v2 <= 2.0 * PI;
        let v_closed = (v2 - v1 - 2.0 * PI).abs() <= PCONFUSION;
        if u_closed && u_first_period {
            let built = torus_param(&tor, v1, v2, false)?;
            return Ok(Some(shift_u(built, true, u1, u2, v1, v2)?));
        }
        if v_closed && v_first_period {
            // `UTrim` defaults to true (`Convert_TorusToBSplineSurface.hxx:62`).
            let built = torus_param(&tor, u1, u2, true)?;
            return Ok(Some(shift_v(built, u1, u2, v1, v2)?));
        }
        return Ok(Some(torus_open(&tor, u1, u2, v1, v2)?));
    }
    Ok(None)
}

fn shift_u(
    surf: GeomBSplineSurface,
    closed: bool,
    u1: f64,
    u2: f64,
    v1: f64,
    v2: f64,
) -> Result<GeomBSplineSurface, String> {
    if !closed {
        return Ok(surf);
    }
    let (uk, _) = surf.distinct_knots_and_mults_u();
    let (vk, _) = surf.distinct_knots_and_mults_v();
    if knots_cover(&uk, u1, u2) && knots_cover(&vk, v1, v2) {
        return Ok(surf);
    }
    let mut surf = surf;
    surf.check_and_segment(u1, u2, v1, v2)?;
    Ok(surf)
}

fn shift_v(
    surf: GeomBSplineSurface,
    u1: f64,
    u2: f64,
    v1: f64,
    v2: f64,
) -> Result<GeomBSplineSurface, String> {
    let (uk, _) = surf.distinct_knots_and_mults_u();
    let (vk, _) = surf.distinct_knots_and_mults_v();
    if knots_cover(&uk, u1, u2) && knots_cover(&vk, v1, v2) {
        return Ok(surf);
    }
    let mut surf = surf;
    surf.check_and_segment(u1, u2, v1, v2)?;
    Ok(surf)
}

fn knots_cover(knots: &[f64], first: f64, last: f64) -> bool {
    let eps = 100.0 * epsilon(2.0 * PI);
    let Some(&k0) = knots.first() else {
        return false;
    };
    let Some(&k1) = knots.last() else {
        return false;
    };
    (k0 - first).abs() <= eps && (k1 - last).abs() <= eps
}

fn span_count(delta: f64) -> (usize, f64) {
    let spans = (1.2 * delta / PI).trunc() as usize + 1;
    let alfa = delta / (spans as f64 * 2.0);
    (spans, alfa)
}

fn collapsed(a: f64, b: f64) -> bool {
    (b - a).abs() <= epsilon(a).abs()
}

fn clamped_knots(start: f64, n_knots: usize, alfa: f64) -> (Vec<f64>, Vec<i32>) {
    let mut knots = Vec::with_capacity(n_knots);
    let mut mults = vec![2i32; n_knots];
    for i in 0..n_knots {
        knots.push(start + i as f64 * 2.0 * alfa);
    }
    if n_knots > 0 {
        mults[0] += 1;
        mults[n_knots - 1] += 1;
    }
    (knots, mults)
}

fn periodic_knots() -> (Vec<f64>, Vec<i32>) {
    let mut knots = Vec::with_capacity(4);
    for i in 0..4 {
        knots.push(i as f64 * 2.0 * PI / 3.0);
    }
    (knots, vec![2, 2, 2, 2])
}

fn place(poles: &mut [Vec<GpPnt>], pos: &GpAx3) {
    let mut tr = GpTrsf::identity();
    tr.set_transformation_from_to(pos, &GpAx3::standard());
    for row in poles {
        for pole in row {
            *pole = pole.transformed(&tr);
        }
    }
}

fn weights_uv(nu: usize, nv: usize, cos_u: f64, cos_v: f64) -> Vec<Vec<f64>> {
    let mut weights = vec![vec![1.0; nv]; nu];
    for i in 0..nu {
        let w1 = if (i + 1) % 2 == 0 { cos_u } else { 1.0 };
        for j in 0..nv {
            let w2 = if (j + 1) % 2 == 0 { cos_v } else { 1.0 };
            weights[i][j] = w1 * w2;
        }
    }
    weights
}

fn finish(
    mut poles: Vec<Vec<GpPnt>>,
    nu: usize,
    nv: usize,
    cos_u: f64,
    cos_v: f64,
    u_knots: Vec<f64>,
    u_mults: Vec<i32>,
    v_knots: Vec<f64>,
    v_mults: Vec<i32>,
    deg_u: usize,
    deg_v: usize,
    u_periodic: bool,
    v_periodic: bool,
    pos: &GpAx3,
) -> Result<GeomBSplineSurface, String> {
    poles.truncate(nu);
    for row in &mut poles {
        row.truncate(nv);
    }
    place(&mut poles, pos);
    let weights = weights_uv(poles.len(), poles.first().map_or(0, |r| r.len()), cos_u, cos_v);
    let period_u = u_knots.last().copied().unwrap_or(0.0) - u_knots.first().copied().unwrap_or(0.0);
    let period_v = v_knots.last().copied().unwrap_or(0.0) - v_knots.first().copied().unwrap_or(0.0);
    let flat_u = flat_knots(&u_knots, &u_mults, deg_u as i32, u_periodic, period_u);
    let flat_v = flat_knots(&v_knots, &v_mults, deg_v as i32, v_periodic, period_v);
    let mut built = GeomBSplineSurface::rational(poles, weights, flat_u, flat_v, deg_u, deg_v)
        .map_err(|e| format!("elementary BSpline: {e}"))?;
    built.set_u_periodic_flag(u_periodic);
    built.set_v_periodic_flag(v_periodic);
    Ok(built)
}

fn cylinder_local(r: f64, u1: f64, u2: f64, v1: f64, v2: f64) -> Vec<Vec<GpPnt>> {
    let (spans, alfa) = span_count(u2 - u1);
    let mut poles = Vec::new();
    let mut u = u1;
    poles.push(vec![
        GpPnt::new(r * u.cos(), r * u.sin(), v1),
        GpPnt::new(r * u.cos(), r * u.sin(), v2),
    ]);
    let c = alfa.cos();
    for _ in 0..spans {
        let mid = u + alfa;
        let end = u + 2.0 * alfa;
        poles.push(vec![
            GpPnt::new(r * mid.cos() / c, r * mid.sin() / c, v1),
            GpPnt::new(r * mid.cos() / c, r * mid.sin() / c, v2),
        ]);
        poles.push(vec![
            GpPnt::new(r * end.cos(), r * end.sin(), v1),
            GpPnt::new(r * end.cos(), r * end.sin(), v2),
        ]);
        u = end;
    }
    poles
}

fn cylinder_open(
    cyl: &occt_core::gp::GpCylinder,
    u1: f64,
    u2: f64,
    v1: f64,
    v2: f64,
) -> Result<GeomBSplineSurface, String> {
    let delta = u2 - u1;
    if collapsed(v1, v2) || delta > 2.0 * PI || delta < 0.0 {
        return Err("Convert_CylinderToBSplineSurface".into());
    }
    let (spans, alfa) = span_count(delta);
    let nu = 2 * spans + 1;
    let (u_knots, u_mults) = clamped_knots(u1, spans + 1, alfa);
    let poles = cylinder_local(cyl.radius(), u1, u2, v1, v2);
    let pos = cyl.position();
    finish(
        poles, nu, 2, alfa.cos(), 1.0, u_knots, u_mults, vec![v1, v2], vec![2, 2], 2, 1, false,
        false, &pos,
    )
}

fn cylinder_periodic(
    cyl: &occt_core::gp::GpCylinder,
    v1: f64,
    v2: f64,
) -> Result<GeomBSplineSurface, String> {
    if collapsed(v1, v2) {
        return Err("Convert_CylinderToBSplineSurface".into());
    }
    let poles = cylinder_local(cyl.radius(), 0.0, 2.0 * PI, v1, v2);
    let (u_knots, u_mults) = periodic_knots();
    let pos = cyl.position();
    finish(
        poles, 6, 2, 0.5, 1.0, u_knots, u_mults, vec![v1, v2], vec![2, 2], 2, 1, true, false,
        &pos,
    )
}

fn cone_local(r: f64, angle: f64, u1: f64, u2: f64, v1: f64, v2: f64) -> Vec<Vec<GpPnt>> {
    let (spans, alfa) = span_count(u2 - u1);
    let x0 = r + v1 * angle.sin();
    let z0 = v1 * angle.cos();
    let x1 = r + v2 * angle.sin();
    let z1 = v2 * angle.cos();
    let mut poles = Vec::new();
    let mut u = u1;
    poles.push(vec![
        GpPnt::new(x0 * u.cos(), x0 * u.sin(), z0),
        GpPnt::new(x1 * u.cos(), x1 * u.sin(), z1),
    ]);
    let c = alfa.cos();
    for _ in 0..spans {
        let mid = u + alfa;
        let end = u + 2.0 * alfa;
        poles.push(vec![
            GpPnt::new(x0 * mid.cos() / c, x0 * mid.sin() / c, z0),
            GpPnt::new(x1 * mid.cos() / c, x1 * mid.sin() / c, z1),
        ]);
        poles.push(vec![
            GpPnt::new(x0 * end.cos(), x0 * end.sin(), z0),
            GpPnt::new(x1 * end.cos(), x1 * end.sin(), z1),
        ]);
        u = end;
    }
    poles
}

fn cone_open(
    cone: &occt_core::gp::GpCone,
    u1: f64,
    u2: f64,
    v1: f64,
    v2: f64,
) -> Result<GeomBSplineSurface, String> {
    let delta = u2 - u1;
    if collapsed(v1, v2) || delta > 2.0 * PI || delta < 0.0 {
        return Err("Convert_ConeToBSplineSurface".into());
    }
    let (spans, alfa) = span_count(delta);
    let nu = 2 * spans + 1;
    let (u_knots, u_mults) = clamped_knots(u1, spans + 1, alfa);
    let poles = cone_local(cone.radius(), cone.semi_angle(), u1, u2, v1, v2);
    let pos = cone.position();
    finish(
        poles, nu, 2, alfa.cos(), 1.0, u_knots, u_mults, vec![v1, v2], vec![2, 2], 2, 1, false,
        false, &pos,
    )
}

fn cone_periodic(cone: &occt_core::gp::GpCone, v1: f64, v2: f64) -> Result<GeomBSplineSurface, String> {
    if collapsed(v1, v2) {
        return Err("Convert_ConeToBSplineSurface".into());
    }
    let poles = cone_local(cone.radius(), cone.semi_angle(), 0.0, 2.0 * PI, v1, v2);
    let (u_knots, u_mults) = periodic_knots();
    let pos = cone.position();
    finish(
        poles, 6, 2, 0.5, 1.0, u_knots, u_mults, vec![v1, v2], vec![2, 2], 2, 1, true, false,
        &pos,
    )
}

fn profile_poles(
    xs: &[f64],
    zs: &[f64],
    u1: f64,
    u2: f64,
) -> Vec<Vec<GpPnt>> {
    let (spans, alfa) = span_count(u2 - u1);
    let nv = xs.len();
    let mut poles = Vec::new();
    let mut u = u1;
    let mut row = Vec::with_capacity(nv);
    for j in 0..nv {
        row.push(GpPnt::new(xs[j] * u.cos(), xs[j] * u.sin(), zs[j]));
    }
    poles.push(row);
    let c = alfa.cos();
    for _ in 0..spans {
        let mid = u + alfa;
        let end = u + 2.0 * alfa;
        let mut mid_row = Vec::with_capacity(nv);
        let mut end_row = Vec::with_capacity(nv);
        for j in 0..nv {
            mid_row.push(GpPnt::new(xs[j] * mid.cos() / c, xs[j] * mid.sin() / c, zs[j]));
            end_row.push(GpPnt::new(xs[j] * end.cos(), xs[j] * end.sin(), zs[j]));
        }
        poles.push(mid_row);
        poles.push(end_row);
        u = end;
    }
    poles
}

fn sphere_profile(r: f64, v1: f64, v2: f64) -> (Vec<f64>, Vec<f64>, usize, f64) {
    let (spans, alfa) = span_count(v2 - v1);
    let mut x = Vec::new();
    let mut z = Vec::new();
    x.push(r * v1.cos());
    z.push(r * v1.sin());
    let mut v = v1;
    let c = alfa.cos();
    for _ in 0..spans {
        let mid = v + alfa;
        let end = v + 2.0 * alfa;
        x.push(r * mid.cos() / c);
        z.push(r * mid.sin() / c);
        x.push(r * end.cos());
        z.push(r * end.sin());
        v = end;
    }
    (x, z, spans, alfa)
}

fn sphere_open(
    sph: &occt_core::gp::GpSphere,
    u1: f64,
    u2: f64,
    v1: f64,
    v2: f64,
) -> Result<GeomBSplineSurface, String> {
    let delta_u = u2 - u1;
    let delta_v = v2 - v1;
    if delta_u > 2.0 * PI || delta_u < 0.0 || v1 < -PI / 2.0 || v2 > PI / 2.0 {
        return Err("Convert_SphereToBSplineSurface".into());
    }
    let (spans_u, alfa_u) = span_count(delta_u);
    let (x, z, spans_v, alfa_v) = sphere_profile(sph.radius(), v1, v2);
    let nu = 2 * spans_u + 1;
    let nv = 2 * spans_v + 1;
    let poles = profile_poles(&x, &z, u1, u2);
    let (u_knots, u_mults) = clamped_knots(u1, spans_u + 1, alfa_u);
    let (v_knots, v_mults) = clamped_knots(v1, spans_v + 1, alfa_v);
    let pos = *sph.position();
    finish(
        poles, nu, nv, alfa_u.cos(), alfa_v.cos(), u_knots, u_mults, v_knots, v_mults, 2, 2,
        false, false, &pos,
    )
}

fn sphere_param(
    sph: &occt_core::gp::GpSphere,
    param1: f64,
    param2: f64,
    u_trim: bool,
) -> Result<GeomBSplineSurface, String> {
    let delta = param2 - param1;
    if delta > 2.0 * PI || delta < 0.0 {
        return Err("Convert_SphereToBSplineSurface".into());
    }
    let r = sph.radius();
    let pos = *sph.position();
    if !u_trim {
        let (x, z, spans_v, alfa_v) = sphere_profile(r, param1, param2);
        let nv = 2 * spans_v + 1;
        let poles = profile_poles(&x, &z, 0.0, 2.0 * PI);
        let (u_knots, u_mults) = periodic_knots();
        let (v_knots, v_mults) = clamped_knots(param1, spans_v + 1, alfa_v);
        finish(
            poles, 6, nv, 0.5, alfa_v.cos(), u_knots, u_mults, v_knots, v_mults, 2, 2, true,
            false, &pos,
        )
    } else {
        let (x, z, _, _) = sphere_profile(r, -PI / 2.0, PI / 2.0);
        let (spans_u, alfa_u) = span_count(delta);
        let nu = 2 * spans_u + 1;
        let poles = profile_poles(&x, &z, param1, param2);
        let (u_knots, u_mults) = clamped_knots(param1, spans_u + 1, alfa_u);
        finish(
            poles,
            nu,
            5,
            alfa_u.cos(),
            0.5,
            u_knots,
            u_mults,
            vec![-PI / 2.0, 0.0, PI / 2.0],
            vec![3, 2, 3],
            2,
            2,
            false,
            false,
            &pos,
        )
    }
}

fn torus_profile(major: f64, minor: f64, v1: f64, v2: f64) -> (Vec<f64>, Vec<f64>, usize, f64) {
    let (spans, alfa) = span_count(v2 - v1);
    let mut x = Vec::new();
    let mut z = Vec::new();
    x.push(major + minor * v1.cos());
    z.push(minor * v1.sin());
    let mut v = v1;
    let c = alfa.cos();
    for _ in 0..spans {
        let mid = v + alfa;
        let end = v + 2.0 * alfa;
        x.push(major + minor * mid.cos() / c);
        z.push(minor * mid.sin() / c);
        x.push(major + minor * end.cos());
        z.push(minor * end.sin());
        v = end;
    }
    (x, z, spans, alfa)
}

fn torus_open(
    tor: &occt_core::gp::GpTorus,
    u1: f64,
    u2: f64,
    v1: f64,
    v2: f64,
) -> Result<GeomBSplineSurface, String> {
    let delta_u = u2 - u1;
    let delta_v = v2 - v1;
    if delta_u > 2.0 * PI || delta_u < 0.0 || delta_v > 2.0 * PI || delta_v < 0.0 {
        return Err("Convert_TorusToBSplineSurface".into());
    }
    let (spans_u, alfa_u) = span_count(delta_u);
    let (x, z, spans_v, alfa_v) = torus_profile(tor.major_radius(), tor.minor_radius(), v1, v2);
    let nu = 2 * spans_u + 1;
    let nv = 2 * spans_v + 1;
    let poles = profile_poles(&x, &z, u1, u2);
    let (u_knots, u_mults) = clamped_knots(u1, spans_u + 1, alfa_u);
    let (v_knots, v_mults) = clamped_knots(v1, spans_v + 1, alfa_v);
    let pos = *tor.position();
    finish(
        poles, nu, nv, alfa_u.cos(), alfa_v.cos(), u_knots, u_mults, v_knots, v_mults, 2, 2,
        false, false, &pos,
    )
}

fn torus_param(
    tor: &occt_core::gp::GpTorus,
    param1: f64,
    param2: f64,
    u_trim: bool,
) -> Result<GeomBSplineSurface, String> {
    let delta = param2 - param1;
    if delta > 2.0 * PI || delta < 0.0 {
        return Err("Convert_TorusToBSplineSurface".into());
    }
    let pos = *tor.position();
    let major = tor.major_radius();
    let minor = tor.minor_radius();
    if !u_trim {
        let (x, z, spans_v, alfa_v) = torus_profile(major, minor, param1, param2);
        let nv = 2 * spans_v + 1;
        let poles = profile_poles(&x, &z, 0.0, 2.0 * PI);
        let (u_knots, u_mults) = periodic_knots();
        let (v_knots, v_mults) = clamped_knots(param1, spans_v + 1, alfa_v);
        finish(
            poles, 6, nv, 0.5, alfa_v.cos(), u_knots, u_mults, v_knots, v_mults, 2, 2, true,
            false, &pos,
        )
    } else {
        let (x, z, _, _) = torus_profile(major, minor, 0.0, 2.0 * PI);
        let (spans_u, alfa_u) = span_count(delta);
        let nu = 2 * spans_u + 1;
        let poles = profile_poles(&x, &z, param1, param2);
        let (u_knots, u_mults) = clamped_knots(param1, spans_u + 1, alfa_u);
        let (v_knots, v_mults) = periodic_knots();
        finish(
            poles, nu, 6, alfa_u.cos(), 0.5, u_knots, u_mults, v_knots, v_mults, 2, 2, false,
            true, &pos,
        )
    }
}
