//! Entity 114. `IGESToBRep_BasicSurface::TransferSplineSurface`
//! (`IGESToBRep_BasicSurface.cxx:559-602`) and
//! `IGESConvGeom::SplineSurfaceFromIGES` (`IGESConvGeom.cxx:286-639`).
//!
//! `IncreaseSurfaceContinuity` (`IGESConvGeom.cxx:644-755`) runs at continuity 1.
//! `HasTransf` calls `ConvertLocation` with the default unit 1
//! (`IGESData_ToolLocation.hxx:148-151`), then `TransferBasicSurface` scales
//! the surface (`cxx:253`). `CompoundLocation` is the entity's own 124, the
//! same gap as `finish_plane`.

use occt_core::gp::{GpPnt, GpTrsf};
use occt_geom::bspline_surface::GeomBSplineSurface;
use occt_geom::surface::Surface;

use super::iges_data_iges_model::{EntityData, IgesModel};
use super::iges_data_param_reader::ParamReader;
use super::iges_to_brep_basic_curve::{convert_location, entity_location, eps_geom};
use super::iges_to_brep_curve_and_surface::Session;

/// `read.iges.bspline.continuity` default (`IGESData.cxx:66`).
const CONTINUITY: i32 = 1;

pub(in crate::iges_read) fn parse_spline_surface(
    pr: &mut ParamReader<'_>,
) -> Result<EntityData, String> {
    // `IGESGeom_ToolSplineSurface.cxx:61-227`.
    let boundary_type = pr.read_integer()?;
    let patch_type = pr.read_integer()?;
    let nb_u = pr.read_integer()?;
    let nb_v = pr.read_integer()?;
    if nb_u < 0 {
        return Err("XSTEP_141: IGESGeom_ToolSplineSurface.cxx:79".into());
    }
    if nb_v < 0 {
        return Err("XSTEP_142: IGESGeom_ToolSplineSurface.cxx:92".into());
    }
    let mut u_breaks = Vec::with_capacity((nb_u as usize) + 1);
    for _ in 0..=nb_u {
        u_breaks.push(pr.read_real()?);
    }
    let mut v_breaks = Vec::with_capacity((nb_v as usize) + 1);
    for _ in 0..=nb_v {
        v_breaks.push(pr.read_real()?);
    }
    let mut x_poly = Vec::with_capacity(nb_u as usize);
    let mut y_poly = Vec::with_capacity(nb_u as usize);
    let mut z_poly = Vec::with_capacity(nb_v as usize);
    for iu in 0..nb_u {
        let mut xr = Vec::with_capacity(nb_v as usize);
        let mut yr = Vec::with_capacity(nb_v as usize);
        let mut zr = Vec::with_capacity(nb_v as usize);
        for iv in 0..nb_v {
            xr.push(read16(pr)?);
            yr.push(read16(pr)?);
            let last = iu + 1 == nb_u && iv + 1 == nb_v;
            zr.push(read16_z(pr, last)?);
        }
        x_poly.push(xr);
        y_poly.push(yr);
        z_poly.push(zr);
        // Arbitrary values after each U row (`cxx:204-208`).
        pr.skip(48);
    }
    if nb_v >= 0 {
        pr.skip(48 * ((nb_v as usize) + 1));
    }
    Ok(EntityData::SplineSurface {
        boundary_type,
        patch_type,
        u_breaks,
        v_breaks,
        x_poly,
        y_poly,
        z_poly,
    })
}

pub(in crate::iges_read) fn transfer_spline_surface(
    session: &Session<'_>,
    index: usize,
) -> Result<GeomBSplineSurface, String> {
    let eps = eps_geom(session.model.resolution);
    let mut surf = spline_surface_from_iges(session, index, eps)?;
    increase_surface_continuity(&mut surf, eps, CONTINUITY)?;
    let unit = session.model.unit_factor;
    if unit != 1.0 {
        let mut scale = GpTrsf::identity();
        scale
            .set_scale(&GpPnt::new(0.0, 0.0, 0.0), unit)
            .map_err(|e| format!("IGES 114 scale: {e} (IGESToBRep_BasicSurface.cxx:253)"))?;
        surf.transform(&scale);
    }
    Ok(surf)
}

fn spline_surface_from_iges(
    session: &Session<'_>,
    index: usize,
    epsgeom: f64,
) -> Result<GeomBSplineSurface, String> {
    let (boundary, u_breaks, v_breaks, x_poly, y_poly, z_poly) =
        match &session.model.entity(index)?.data {
            EntityData::SplineSurface {
                boundary_type,
                u_breaks,
                v_breaks,
                x_poly,
                y_poly,
                z_poly,
                ..
            } => (
                *boundary_type,
                u_breaks.clone(),
                v_breaks.clone(),
                x_poly.clone(),
                y_poly.clone(),
                z_poly.clone(),
            ),
            _ => return Err("IGES 114 has no spline-surface parameters".into()),
        };
    let mut degree = boundary;
    if degree > 3 {
        degree = 3;
    }
    if degree < 0 {
        return Err("Standard_RangeError: IGESConvGeom.cxx:368".into());
    }
    let nb_u = u_breaks.len().saturating_sub(1);
    let nb_v = v_breaks.len().saturating_sub(1);
    if nb_u < 1 || nb_v < 1 || x_poly.len() < nb_u {
        return Err("IGES_1305: IGESToBRep_BasicSurface.cxx:577-580".into());
    }
    let degree_us = degree as usize;
    let nu = nb_u * degree_us + 1;
    let nv = nb_v * degree_us + 1;
    let mut poles = vec![vec![GpPnt::new(0.0, 0.0, 0.0); nv]; nu];
    let mut was_c0 = true;

    let bz = patch_poles(&x_poly, &y_poly, &z_poly, 1, 1, degree_us, &u_breaks, &v_breaks)?;
    copy_patch(&mut poles, 0, 0, &bz)?;

    for u_seg in 2..=nb_u {
        let bz = patch_poles(
            &x_poly, &y_poly, &z_poly, u_seg, 1, degree_us, &u_breaks, &v_breaks,
        )?;
        let mut u = (u_seg - 1) * degree_us;
        let mut v = 0usize;
        for j in 0..bz[0].len() {
            write_c0(&mut poles, u, v, bz[0][j], epsgeom, &mut was_c0)?;
            v += 1;
        }
        u += 1;
        for i in 1..bz.len() {
            v = 0;
            for j in 0..bz[i].len() {
                set_pole(&mut poles, u, v, bz[i][j])?;
                v += 1;
            }
            u += 1;
        }
    }

    for v_seg in 2..=nb_v {
        let bz = patch_poles(
            &x_poly, &y_poly, &z_poly, 1, v_seg, degree_us, &u_breaks, &v_breaks,
        )?;
        let mut u = 0usize;
        let v = (v_seg - 1) * degree_us;
        for i in 0..bz.len() {
            write_c0(&mut poles, u, v, bz[i][0], epsgeom, &mut was_c0)?;
            u += 1;
        }
        let mut v = v + 1;
        for j in 1..bz[0].len() {
            u = 0;
            for i in 0..bz.len() {
                set_pole(&mut poles, u, v, bz[i][j])?;
                u += 1;
            }
            v += 1;
        }
    }

    for v_seg in 2..=nb_v {
        for u_seg in 2..=nb_u {
            let bz = patch_poles(
                &x_poly, &y_poly, &z_poly, u_seg, v_seg, degree_us, &u_breaks, &v_breaks,
            )?;
            let u0 = (u_seg - 1) * degree_us;
            let v0 = (v_seg - 1) * degree_us;
            let mut u = u0;
            for i in 0..bz.len() {
                write_c0(&mut poles, u, v0, bz[i][0], epsgeom, &mut was_c0)?;
                u += 1;
            }
            u = u0;
            let mut v = v0;
            for j in 0..bz[0].len() {
                write_c0(&mut poles, u, v, bz[0][j], epsgeom, &mut was_c0)?;
                v += 1;
            }
            u = u0 + 1;
            for i in 1..bz.len() {
                v = v0 + 1;
                for j in 1..bz[i].len() {
                    set_pole(&mut poles, u, v, bz[i][j])?;
                    v += 1;
                }
                u += 1;
            }
        }
    }
    let _ = was_c0;

    if let Some(g) = entity_location(session.model, index)? {
        if let Some(tr) = convert_location(1.0e-4, &g, 1.0) {
            for row in &mut poles {
                for p in row.iter_mut() {
                    *p = p.transformed(&tr);
                }
            }
        }
    }

    let u_knots = u_breaks[..nb_u + 1].to_vec();
    let v_knots = v_breaks[..nb_v + 1].to_vec();
    GeomBSplineSurface::from_poles_knots_mults(
        poles,
        u_knots,
        v_knots,
        end_mults(nb_u + 1, degree),
        end_mults(nb_v + 1, degree),
        degree_us,
        degree_us,
    )
    .map_err(|e| format!("IGES 114 surface: {e} (IGESToBRep_BasicSurface.cxx:238-245)"))
}

fn patch_poles(
    x_poly: &[Vec<[f64; 16]>],
    y_poly: &[Vec<[f64; 16]>],
    z_poly: &[Vec<[f64; 16]>],
    u_seg: usize,
    v_seg: usize,
    degree: usize,
    u_breaks: &[f64],
    v_breaks: &[f64],
) -> Result<Vec<Vec<GpPnt>>, String> {
    let x = x_poly
        .get(u_seg - 1)
        .and_then(|row| row.get(v_seg - 1))
        .ok_or_else(|| "IGES 114 missing polynomial".to_string())?;
    let y = y_poly
        .get(u_seg - 1)
        .and_then(|row| row.get(v_seg - 1))
        .ok_or_else(|| "IGES 114 missing polynomial".to_string())?;
    let z = z_poly
        .get(u_seg - 1)
        .and_then(|row| row.get(v_seg - 1))
        .ok_or_else(|| "IGES 114 missing polynomial".to_string())?;
    let du = u_breaks[u_seg] - u_breaks[u_seg - 1];
    let dv = v_breaks[v_seg] - v_breaks[v_seg - 1];
    let n = degree + 1;
    let mut coef = vec![vec![GpPnt::new(0.0, 0.0, 0.0); n]; n];
    let mut param_u = 1.0;
    for i in 1..=n {
        let mut param_v = 1.0;
        for j in 1..=n {
            let k = i + 4 * (j - 1) - 1;
            let scale = param_u * param_v;
            coef[i - 1][j - 1] = GpPnt::new(x[k] * scale, y[k] * scale, z[k] * scale);
            param_v *= dv;
        }
        param_u *= du;
    }
    Ok(occt_core::bspl::plib::coefficients_poles_grid(&coef))
}

fn copy_patch(
    poles: &mut [Vec<GpPnt>],
    mut u: usize,
    v0: usize,
    bz: &[Vec<GpPnt>],
) -> Result<(), String> {
    for row in bz {
        let mut v = v0;
        for p in row {
            set_pole(poles, u, v, *p)?;
            v += 1;
        }
        u += 1;
    }
    Ok(())
}

fn write_c0(
    poles: &mut [Vec<GpPnt>],
    u: usize,
    v: usize,
    bz: GpPnt,
    eps: f64,
    was_c0: &mut bool,
) -> Result<(), String> {
    let prev = pole_at(poles, u, v)?;
    if prev.distance(&bz) > eps {
        *was_c0 = false;
        set_pole(
            poles,
            u,
            v,
            GpPnt::new(
                0.5 * (prev.x() + bz.x()),
                0.5 * (prev.y() + bz.y()),
                0.5 * (prev.z() + bz.z()),
            ),
        )
    } else {
        set_pole(poles, u, v, bz)
    }
}

fn set_pole(poles: &mut [Vec<GpPnt>], u: usize, v: usize, p: GpPnt) -> Result<(), String> {
    let row = poles.get_mut(u).ok_or_else(pole_err)?;
    let cell = row.get_mut(v).ok_or_else(pole_err)?;
    *cell = p;
    Ok(())
}

fn pole_at(poles: &[Vec<GpPnt>], u: usize, v: usize) -> Result<GpPnt, String> {
    poles
        .get(u)
        .and_then(|row| row.get(v))
        .copied()
        .ok_or_else(pole_err)
}

fn pole_err() -> String {
    "IGES 114 pole index: IGESToBRep_BasicSurface.cxx:238-245".into()
}

fn end_mults(n: usize, degree: i32) -> Vec<i32> {
    let mut mults = vec![degree; n];
    if let Some(first) = mults.first_mut() {
        *first = degree + 1;
    }
    if let Some(last) = mults.last_mut() {
        *last = degree + 1;
    }
    mults
}

/// `IGESConvGeom::IncreaseSurfaceContinuity` (`cxx:644-755`).
fn increase_surface_continuity(
    surf: &mut GeomBSplineSurface,
    eps: f64,
    continuity: i32,
) -> Result<(), String> {
    if continuity < 1 {
        return Ok(());
    }
    let degree_u = surf.deg_u as i32;
    loop {
        let (knots, mults) = surf.distinct_knots_and_mults_u();
        let last = occt_core::bspl::locate::last_u_knot_index(degree_u, &mults);
        let mut i = occt_core::bspl::locate::first_u_knot_index(degree_u, &mults) + 1;
        let mut modified = false;
        while i < last {
            let (_, mults) = surf.distinct_knots_and_mults_u();
            let mult = mults.get((i - 1) as usize).copied().unwrap_or(0);
            if degree_u - mult < continuity {
                let target = if continuity >= 2 {
                    degree_u - 2
                } else {
                    degree_u - 1
                };
                if remove_u_knot(surf, i, target, eps)? {
                    modified = true;
                } else if continuity >= 2 {
                    modified |= remove_u_knot(surf, i, degree_u - 1, eps)?;
                }
            }
            i += 1;
        }
        if !modified {
            break;
        }
    }
    let degree_v = surf.deg_v as i32;
    loop {
        let (_, mults) = surf.distinct_knots_and_mults_v();
        let last = occt_core::bspl::locate::last_u_knot_index(degree_v, &mults);
        let mut i = occt_core::bspl::locate::first_u_knot_index(degree_v, &mults) + 1;
        let mut modified = false;
        while i < last {
            let (_, mults) = surf.distinct_knots_and_mults_v();
            let mult = mults.get((i - 1) as usize).copied().unwrap_or(0);
            if degree_v - mult < continuity {
                let target = if continuity >= 2 {
                    degree_v - 2
                } else {
                    degree_v - 1
                };
                if remove_v_knot(surf, i, target, eps)? {
                    modified = true;
                } else if continuity >= 2 {
                    modified |= remove_v_knot(surf, i, degree_v - 1, eps)?;
                }
            }
            i += 1;
        }
        if !modified {
            break;
        }
    }
    Ok(())
}

fn remove_u_knot(
    surf: &mut GeomBSplineSurface,
    index: i32,
    m: i32,
    tolerance: f64,
) -> Result<bool, String> {
    if m < 0 {
        return Ok(true);
    }
    let (knots, mults) = surf.distinct_knots_and_mults_u();
    let degree = surf.deg_u as i32;
    let first = occt_core::bspl::locate::first_u_knot_index(degree, &mults);
    let last = occt_core::bspl::locate::last_u_knot_index(degree, &mults);
    if index <= first || index >= last {
        return Err("Geom_BSplineSurface::RemoveUKnot: invalid Index".into());
    }
    let step = mults[(index as usize) - 1] - m;
    if step <= 0 {
        return Ok(true);
    }
    let nu = surf.nb_poles_u();
    let nv = surf.nb_poles_v();
    let new_nu = nu as i32 - step;
    if new_nu < 1 {
        return Ok(false);
    }
    let new_nk = knots.len() as i32 - if m == 0 { 1 } else { 0 };
    let mut flat = Vec::with_capacity(3 * nu * nv);
    for row in &surf.poles {
        for p in row {
            flat.push(p.x());
            flat.push(p.y());
            flat.push(p.z());
        }
    }
    let dim = 3 * nv as i32;
    let mut new_flat = vec![0.0; 3 * new_nu as usize * nv];
    let mut new_knots = vec![0.0; new_nk as usize];
    let mut new_mults = vec![0i32; new_nk as usize];
    let ok = occt_core::bspl::remove_knot::remove_knot_dim(
        index,
        m,
        degree,
        false,
        dim,
        &flat,
        &knots,
        &mults,
        &mut new_flat,
        &mut new_knots,
        &mut new_mults,
        tolerance,
    );
    if !ok {
        return Ok(false);
    }
    let (v_knots, v_mults) = surf.distinct_knots_and_mults_v();
    let deg_v = surf.deg_v;
    let mut poles = Vec::with_capacity(new_nu as usize);
    let mut cursor = 0usize;
    for _ in 0..new_nu {
        let mut row = Vec::with_capacity(nv);
        for _ in 0..nv {
            row.push(GpPnt::new(
                new_flat[cursor],
                new_flat[cursor + 1],
                new_flat[cursor + 2],
            ));
            cursor += 3;
        }
        poles.push(row);
    }
    *surf = GeomBSplineSurface::from_poles_knots_mults(
        poles,
        new_knots,
        v_knots,
        new_mults,
        v_mults,
        degree as usize,
        deg_v,
    )
    .map_err(|e| format!("IGES 114 RemoveUKnot: {e}"))?;
    Ok(true)
}

fn remove_v_knot(
    surf: &mut GeomBSplineSurface,
    index: i32,
    m: i32,
    tolerance: f64,
) -> Result<bool, String> {
    if m < 0 {
        return Ok(true);
    }
    let (knots, mults) = surf.distinct_knots_and_mults_v();
    let degree = surf.deg_v as i32;
    let first = occt_core::bspl::locate::first_u_knot_index(degree, &mults);
    let last = occt_core::bspl::locate::last_u_knot_index(degree, &mults);
    if index <= first || index >= last {
        return Err("Geom_BSplineSurface::RemoveVKnot: invalid Index".into());
    }
    let step = mults[(index as usize) - 1] - m;
    if step <= 0 {
        return Ok(true);
    }
    let nu = surf.nb_poles_u();
    let nv = surf.nb_poles_v();
    let new_nv = nv as i32 - step;
    if new_nv < 1 {
        return Ok(false);
    }
    let new_nk = knots.len() as i32 - if m == 0 { 1 } else { 0 };
    // Column-major `(v, u)` (`BSplSLib::SetPoles`, `UDirection == false`).
    let mut flat = Vec::with_capacity(3 * nu * nv);
    for v in 0..nv {
        for u in 0..nu {
            let p = surf.poles[u][v];
            flat.push(p.x());
            flat.push(p.y());
            flat.push(p.z());
        }
    }
    let dim = 3 * nu as i32;
    let mut new_flat = vec![0.0; 3 * nu * new_nv as usize];
    let mut new_knots = vec![0.0; new_nk as usize];
    let mut new_mults = vec![0i32; new_nk as usize];
    let ok = occt_core::bspl::remove_knot::remove_knot_dim(
        index,
        m,
        degree,
        false,
        dim,
        &flat,
        &knots,
        &mults,
        &mut new_flat,
        &mut new_knots,
        &mut new_mults,
        tolerance,
    );
    if !ok {
        return Ok(false);
    }
    let (u_knots, u_mults) = surf.distinct_knots_and_mults_u();
    let deg_u = surf.deg_u;
    let mut poles = vec![Vec::with_capacity(new_nv as usize); nu];
    let mut cursor = 0usize;
    for _v in 0..new_nv {
        for u in 0..nu {
            poles[u].push(GpPnt::new(
                new_flat[cursor],
                new_flat[cursor + 1],
                new_flat[cursor + 2],
            ));
            cursor += 3;
        }
    }
    *surf = GeomBSplineSurface::from_poles_knots_mults(
        poles,
        u_knots,
        new_knots,
        u_mults,
        new_mults,
        deg_u,
        degree as usize,
    )
    .map_err(|e| format!("IGES 114 RemoveVKnot: {e}"))?;
    Ok(true)
}

fn read16(pr: &mut ParamReader<'_>) -> Result<[f64; 16], String> {
    let mut out = [0.0; 16];
    for slot in &mut out {
        *slot = pr.read_real()?;
    }
    Ok(out)
}

fn read16_z(pr: &mut ParamReader<'_>, last: bool) -> Result<[f64; 16], String> {
    if !last || pr.remaining() >= 16 {
        return read16(pr);
    }
    // Last patch may be short (`cxx:183-201`). Missing reals stay 0.
    let mut out = [0.0; 16];
    for slot in &mut out {
        match pr.read_real() {
            Ok(v) => *slot = v,
            Err(_) => break,
        }
    }
    Ok(out)
}
