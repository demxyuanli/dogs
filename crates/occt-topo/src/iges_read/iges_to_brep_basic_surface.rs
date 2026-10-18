//! `IGESToBRep_BasicSurface` for entities 128, 190, 192, 194, 196, and 198
//! (`IGESToBRep_BasicSurface.cxx:180-1076`).
//!
//! Entity 114 is `TransferSplineSurface` (`cxx:559-602`) in `spline_surface`.
//!
//! `GetModeTransfer` is `myModeIsTopo` (`IGESToBRep_CurveAndSurface.hxx:101-104`),
//! left true by `IGESToBRep_Actor.cxx:167`. The pole `gp_GTrsf` blocks
//! (`BasicSurface.cxx:656-678`, `:908-927`) are therefore not entered. The face
//! location is applied by `TransferTopoBasicSurface` (`cxx:312-329`).
//!
//! `read.iges.bspline.continuity` defaults to 1 (`IGESData.cxx:66`), so the
//! `RemoveUKnot` / `RemoveVKnot` loop (`cxx:1022-1056`, `icont >= 2`) is not
//! entered. `Geom_BSplineSurface::RemoveUKnot` is not ported.

use std::f64::consts::PI;
use std::sync::Arc;

use occt_core::bspl::banded_interp::knot_sequence;
use occt_core::gp::{GpAx1, GpAx3, GpCylinder, GpDir, GpPln, GpPnt, GpSphere, GpTorus, GpXyz};
use occt_core::precision::{CONFUSION, PCONFUSION};
use occt_core::toploc::TopLocLocation;
use occt_geom::bspline_curve::MAX_DEGREE;
use occt_geom::bspline_surface::GeomBSplineSurface;
use occt_geom::{GeomCone, GeomCylinder, GeomPlane, GeomSphere, GeomTorus, Surface};

use crate::brep_lib_make_face::make_face_from_surface;
use crate::builder::TopoBuilder;
use crate::shape::TopoShape;
use crate::shape_custom_surface::convert_to_periodic;

use super::iges_data_iges_model::EntityData;
use super::iges_to_brep_basic_curve::{convert_location, entity_location};
use super::iges_to_brep_curve_and_surface::Session;

pub(in crate::iges_read) fn transfer_topo_basic_surface(
    session: &Session<'_>,
    index: usize,
) -> Result<TopoShape, String> {
    // `IGESToBRep_TopoSurface::TransferTopoBasicSurface` (`cxx:260-331`).
    let type_number = session.model.entity(index)?.type_number;
    let surface = transfer_basic_surface(session, index)?;
    let builder = TopoBuilder::new();
    let mut face = if type_number == 190 {
        let pln = surface.gp_pln().ok_or_else(|| {
            "IGES 190 is not a plane: IGESToBRep_BasicSurface.cxx:305".to_string()
        })?;
        builder.make_face_plane(&pln)
    } else {
        make_face_from_surface(surface, CONFUSION)
    };
    if let Some(g) = entity_location(session.model, index)? {
        if let Some(tr) = convert_location(1.0e-4, &g, session.model.unit_factor) {
            face.0
                .move_location(&TopLocLocation::composed(&tr, TopLocLocation::identity()));
        }
    }
    Ok(face.0)
}

fn transfer_basic_surface(
    session: &Session<'_>,
    index: usize,
) -> Result<Arc<dyn Surface>, String> {
    // `TransferBasicSurface` (`cxx:180-256`). The unit scale is
    // `resurf->Scale(origin, GetUnitFactor())` (`cxx:253`).
    let type_number = session.model.entity(index)?.type_number;
    let surf: Arc<dyn Surface> = match type_number {
        128 => Arc::new(transfer_bspline_surface(session, index)?),
        114 => Arc::new(super::spline_surface::transfer_spline_surface(session, index)?),
        190 => Arc::new(GeomPlane::new(transfer_plane_surface(session, index)?)),
        192 => Arc::new(GeomCylinder::new(transfer_cylinder(session, index)?)),
        194 => Arc::new(GeomCone::new(transfer_cone(session, index)?)),
        196 => Arc::new(GeomSphere::new(transfer_sphere(session, index)?)),
        198 => Arc::new(GeomTorus::new(transfer_torus(session, index)?)),
        n => {
            return Err(format!(
                "UNPORTED basic surface {n}: IGESToBRep_BasicSurface.cxx:231-234"
            ));
        }
    };
    Ok(surf)
}

fn point_of(session: &Session<'_>, index: usize) -> Result<GpPnt, String> {
    match &session.model.entity(index)?.data {
        EntityData::Point { xyz } => Ok(GpPnt::new(xyz[0], xyz[1], xyz[2])),
        _ => Err(format!(
            "IGES entity {index} is not a point: IGESToBRep_BasicSurface.cxx:278-282"
        )),
    }
}

fn direction_of(session: &Session<'_>, index: usize) -> Result<GpDir, String> {
    match &session.model.entity(index)?.data {
        EntityData::Direction { xyz } => GpDir::from_xyz(&GpXyz::new(xyz[0], xyz[1], xyz[2]))
            .map_err(|_| format!("IGES entity {index} direction is null: IGESToBRep_BasicSurface.cxx:285")),
        _ => Err(format!(
            "IGES entity {index} is not a direction: IGESToBRep_BasicSurface.cxx:285-289"
        )),
    }
}

fn scale_pnt(p: GpPnt, unit: f64) -> GpPnt {
    GpPnt::new(p.x() * unit, p.y() * unit, p.z() * unit)
}

fn ax3(origin: GpPnt, normal: GpDir, reference: Option<GpDir>) -> Result<GpAx3, String> {
    let Some(vx) = reference else {
        return Ok(ax3_normal(origin, normal));
    };
    let cross = vx.xyz().crossed(normal.xyz());
    if cross.modulus() < CONFUSION {
        return Err(
            "IGES reference direction is parallel to the axis: IGESToBRep_BasicSurface.cxx:357-360"
                .into(),
        );
    }
    // `gp_Ax3(P, N, Vx)` (`gp_Ax3.hxx:84-90`): `vxdir.CrossCross(Vx, N)`.
    let x_xyz = cross_cross(normal.xyz(), vx.xyz(), normal.xyz());
    let xdir = GpDir::from_xyz(&x_xyz).map_err(|_| "IGES axis X direction is null".to_string())?;
    let ydir = normal
        .crossed(&xdir)
        .map_err(|_| "IGES axis Y direction is null".to_string())?;
    Ok(GpAx3 {
        axis: GpAx1::new(origin, normal),
        vxdir: xdir,
        vydir: ydir,
    })
}

fn ax3_normal(origin: GpPnt, normal: GpDir) -> GpAx3 {
    // `gp_Ax3(P, V)` (`gp_Ax3.cxx:29-79`).
    let a = normal.x();
    let b = normal.y();
    let c = normal.z();
    let aabs = a.abs();
    let babs = b.abs();
    let cabs = c.abs();
    let (dx, dy, dz) = if babs <= aabs && babs <= cabs {
        if aabs > cabs {
            (-c, 0.0, a)
        } else {
            (c, 0.0, -a)
        }
    } else if aabs <= babs && aabs <= cabs {
        if babs > cabs {
            (0.0, -c, b)
        } else {
            (0.0, c, -b)
        }
    } else if aabs > babs {
        (-b, a, 0.0)
    } else {
        (b, -a, 0.0)
    };
    let xdir = GpDir::new(dx, dy, dz).unwrap_or(normal);
    let ydir = normal.crossed(&xdir).unwrap_or(xdir);
    GpAx3 {
        axis: GpAx1::new(origin, normal),
        vxdir: xdir,
        vydir: ydir,
    }
}

fn cross_cross(this: &GpXyz, v1: &GpXyz, v2: &GpXyz) -> GpXyz {
    // `gp_XYZ::CrossCross` (`gp_XYZ.hxx:551-564`).
    let cx = v1.y() * v2.z() - v1.z() * v2.y();
    let cy = v1.z() * v2.x() - v1.x() * v2.z();
    let cz = v1.x() * v2.y() - v1.y() * v2.x();
    GpXyz::new(
        this.y() * cz - this.z() * cy,
        this.z() * cx - this.x() * cz,
        this.x() * cy - this.y() * cx,
    )
}

fn transfer_plane_surface(session: &Session<'_>, index: usize) -> Result<GpPln, String> {
    // `TransferPlaneSurface` (`cxx:261-305`).
    let (location, normal, reference) = match &session.model.entity(index)?.data {
        EntityData::PlaneSurface {
            location,
            normal,
            reference,
        } => (*location, *normal, *reference),
        _ => return Err("IGES 190 has no plane parameters".into()),
    };
    let origin = scale_pnt(point_of(session, location)?, session.model.unit_factor);
    let n = direction_of(session, normal)?;
    let reference = match reference {
        Some(r) => Some(direction_of(session, r)?),
        None => None,
    };
    Ok(GpPln::new(ax3(origin, n, reference)?))
}

fn transfer_cylinder(session: &Session<'_>, index: usize) -> Result<GpCylinder, String> {
    // `TransferRigthCylindricalSurface` (`cxx:310-365`).
    let (location, axis, radius, reference) = match &session.model.entity(index)?.data {
        EntityData::CylindricalSurface {
            location,
            axis,
            radius,
            reference,
        } => (*location, *axis, *radius, *reference),
        _ => return Err("IGES 192 has no cylinder parameters".into()),
    };
    if radius < CONFUSION {
        return Err("IGES 192 radius is below confusion: IGESToBRep_BasicSurface.cxx:340-342".into());
    }
    let origin = scale_pnt(point_of(session, location)?, session.model.unit_factor);
    let ax = direction_of(session, axis)?;
    let reference = match reference {
        Some(r) => Some(direction_of(session, r)?),
        None => None,
    };
    let pos = ax3(origin, ax, reference)?;
    GpCylinder::new(pos, radius * session.model.unit_factor)
        .map_err(|_| "IGES 192 cylinder: IGESToBRep_BasicSurface.cxx:364".to_string())
}

fn transfer_cone(session: &Session<'_>, index: usize) -> Result<occt_core::gp::GpCone, String> {
    // `TransferRigthConicalSurface` (`cxx:370-433`). Semi-angle is degrees.
    let (location, axis, radius, semi_angle, reference) = match &session.model.entity(index)?.data {
        EntityData::ConicalSurface {
            location,
            axis,
            radius,
            semi_angle,
            reference,
        } => (*location, *axis, *radius, *semi_angle, *reference),
        _ => return Err("IGES 194 has no cone parameters".into()),
    };
    let angle = semi_angle / 180.0 * PI;
    let mut radius = radius;
    if angle < CONFUSION || angle > PI / 2.0 {
        return Err("IGES 194 semi-angle is out of range: IGESToBRep_BasicSurface.cxx:401-403".into());
    }
    if radius < 0.0 {
        return Err("IGES 194 radius is negative: IGESToBRep_BasicSurface.cxx:405-407".into());
    }
    if radius < CONFUSION {
        radius = 0.0;
    }
    let origin = scale_pnt(point_of(session, location)?, session.model.unit_factor);
    let ax = direction_of(session, axis)?;
    let reference = match reference {
        Some(r) => Some(direction_of(session, r)?),
        None => None,
    };
    let pos = ax3(origin, ax, reference)?;
    occt_core::gp::GpCone::new(pos, radius * session.model.unit_factor, angle)
        .map_err(|_| "IGES 194 cone: IGESToBRep_BasicSurface.cxx:433".to_string())
}

fn transfer_sphere(session: &Session<'_>, index: usize) -> Result<GpSphere, String> {
    // `TransferSphericalSurface` (`cxx:438-493`). Form 0 leaves the axis null.
    let (center, radius, axis, reference) = match &session.model.entity(index)?.data {
        EntityData::SphericalSurface {
            center,
            radius,
            axis,
            reference,
        } => (*center, *radius, *axis, *reference),
        _ => return Err("IGES 196 has no sphere parameters".into()),
    };
    let Some(axis) = axis else {
        return Err("IGES 196 axis is null: IGESToBRep_BasicSurface.cxx:460-462".into());
    };
    if radius < CONFUSION {
        return Err("IGES 196 radius is below confusion: IGESToBRep_BasicSurface.cxx:468-471".into());
    }
    let origin = scale_pnt(point_of(session, center)?, session.model.unit_factor);
    let ax = direction_of(session, axis)?;
    let reference = match reference {
        Some(r) => Some(direction_of(session, r)?),
        None => None,
    };
    let pos = ax3(origin, ax, reference)?;
    GpSphere::new(pos, radius * session.model.unit_factor)
        .map_err(|_| "IGES 196 sphere: IGESToBRep_BasicSurface.cxx:493".to_string())
}

fn transfer_torus(session: &Session<'_>, index: usize) -> Result<GpTorus, String> {
    // `TransferToroidalSurface` (`cxx:498-554`).
    let (center, axis, major, minor, reference) = match &session.model.entity(index)?.data {
        EntityData::ToroidalSurface {
            center,
            axis,
            major,
            minor,
            reference,
        } => (*center, *axis, *major, *minor, *reference),
        _ => return Err("IGES 198 has no torus parameters".into()),
    };
    if major < CONFUSION || minor < CONFUSION {
        return Err(
            "IGES 198 radius is below confusion: IGESToBRep_BasicSurface.cxx:529-532".into(),
        );
    }
    let origin = scale_pnt(point_of(session, center)?, session.model.unit_factor);
    let ax = direction_of(session, axis)?;
    let reference = match reference {
        Some(r) => Some(direction_of(session, r)?),
        None => None,
    };
    let pos = ax3(origin, ax, reference)?;
    let unit = session.model.unit_factor;
    GpTorus::new(pos, major * unit, minor * unit)
        .map_err(|_| "IGES 198 torus: IGESToBRep_BasicSurface.cxx:554".to_string())
}

fn real_epsilon(x: f64) -> f64 {
    let a = x.abs();
    if a > 1.0 {
        f64::EPSILON * a
    } else {
        f64::EPSILON
    }
}

fn fix_knots(knots: &mut [f64]) {
    // `ShapeConstruct_Curve::FixKnots` (`ShapeConstruct_Curve.cxx:397-413`).
    if knots.is_empty() {
        return;
    }
    let mut knot_val = knots[0];
    for i in 1..knots.len() {
        let mut knot_next = knots[i];
        if knot_next - knot_val <= real_epsilon(knot_val) {
            knot_next = knot_val + 2.0 * real_epsilon(knot_val);
            knots[i] = knot_next;
        }
        knot_val = knot_next;
    }
}

struct PackedKnots {
    knots: Vec<f64>,
    mults: Vec<i32>,
    keep: Vec<usize>,
}

fn pack_knots(raw: &[f64], degree: usize, n_poles: usize) -> Result<PackedKnots, String> {
    // `TransferBSplineSurface` knot merge (`cxx:684-746` and the V copy).
    if raw.is_empty() {
        return Err("IGES 128 has no knots".into());
    }
    let mut temp_k = vec![raw[0]];
    let mut temp_m = vec![1i32];
    for i in 1..raw.len() {
        let k1 = raw[i];
        let k2 = raw[i - 1];
        if (k1 - k2).abs() <= real_epsilon(k2) {
            let last = temp_m.len() - 1;
            temp_m[last] += 1;
        } else {
            temp_k.push(k1);
            temp_m.push(1);
        }
    }
    let mut new_nb = n_poles;
    let mut seq: Vec<usize> = Vec::new();
    let mut old_sum = 0i32;
    let mut mults = Vec::with_capacity(temp_k.len());
    for i in 0..temp_k.len() {
        let mut mult = temp_m[i];
        let max_mult = if i == 0 || i + 1 == temp_k.len() {
            degree as i32 + 1
        } else {
            degree as i32
        };
        if mult > max_mult {
            let mut del = old_sum + 1;
            while mult > max_mult {
                new_nb = new_nb.saturating_sub(1);
                seq.push(del as usize);
                del += 1;
                mult -= 1;
            }
        }
        old_sum += temp_m[i];
        mults.push(mult);
    }
    let mut keep: Vec<usize> = (0..n_poles).collect();
    if new_nb < n_poles {
        let mut offset = 0usize;
        for del in &seq {
            let idx = del.saturating_sub(1 + offset);
            if idx < keep.len() {
                keep.remove(idx);
                offset += 1;
            }
        }
    }
    if keep.len() != new_nb {
        keep.clear();
    }
    for i in 1..temp_k.len() {
        if temp_k[i] < temp_k[i - 1] {
            return Err(
                "IGES 128 knots are not ascending: IGESToBRep_BasicSurface.cxx:132-136".into(),
            );
        }
    }
    fix_knots(&mut temp_k);
    Ok(PackedKnots {
        knots: temp_k,
        mults,
        keep,
    })
}

fn transfer_bspline_surface(
    session: &Session<'_>,
    index: usize,
) -> Result<GeomBSplineSurface, String> {
    // `TransferBSplineSurface` (`cxx:607-1076`).
    let data = match &session.model.entity(index)?.data {
        EntityData::BSplineSurface {
            degree_u,
            degree_v,
            closed_u,
            closed_v,
            polynomial,
            periodic_u,
            periodic_v,
            knots_u,
            knots_v,
            weights,
            poles,
            ..
        } => (
            *degree_u,
            *degree_v,
            *closed_u,
            *closed_v,
            *polynomial,
            *periodic_u,
            *periodic_v,
            knots_u.clone(),
            knots_v.clone(),
            weights.clone(),
            poles.clone(),
        ),
        _ => return Err("IGES 128 has no B-spline parameters".into()),
    };
    let (
        degree_u,
        degree_v,
        closed_u,
        closed_v,
        polynomial,
        periodic_u,
        periodic_v,
        knots_u,
        knots_v,
        weights,
        poles,
    ) = data;
    if degree_u <= 0
        || degree_v <= 0
        || degree_u as usize > MAX_DEGREE
        || degree_v as usize > MAX_DEGREE
    {
        return Err("IGES 128 degree is out of range: IGESToBRep_BasicSurface.cxx:620-626".into());
    }
    let nu = poles.len();
    let nv = poles.first().map(|r| r.len()).unwrap_or(0);
    if nu < 2 || nv < 2 {
        return Err("IGES 128 has fewer than 2 poles: IGESToBRep_BasicSurface.cxx:634-648".into());
    }
    let packed_u = pack_knots(&knots_u, degree_u as usize, nu)?;
    let packed_v = pack_knots(&knots_v, degree_v as usize, nv)?;
    let unit = session.model.unit_factor;
    let poles = if packed_u.keep.is_empty() || packed_v.keep.is_empty() {
        // Multiplicity correction did not line up (`cxx:885-889`). Poles stay at the origin.
        let nrow = packed_u.keep.len().max(1);
        let ncol = packed_v.keep.len().max(1);
        vec![vec![GpPnt::new(0.0, 0.0, 0.0); ncol]; nrow]
    } else {
        let mut grid = Vec::with_capacity(packed_u.keep.len());
        for &iu in &packed_u.keep {
            let mut row = Vec::with_capacity(packed_v.keep.len());
            for &iv in &packed_v.keep {
                let p = poles[iu][iv];
                row.push(scale_pnt(GpPnt::new(p[0], p[1], p[2]), unit));
            }
            grid.push(row);
        }
        grid
    };
    let flat_u = knot_sequence(&packed_u.knots, &packed_u.mults, degree_u);
    let flat_v = knot_sequence(&packed_v.knots, &packed_v.mults, degree_v);
    let mut surf = if polynomial {
        GeomBSplineSurface::new(poles, flat_u, flat_v, degree_u as usize, degree_v as usize)
    } else {
        let wgrid = weight_grid(&weights, &packed_u.keep, &packed_v.keep)?;
        GeomBSplineSurface::rational(
            poles,
            wgrid,
            flat_u,
            flat_v,
            degree_u as usize,
            degree_v as usize,
        )
    }
    .map_err(|e| format!("IGES 128 surface: {e} (IGESToBRep_BasicSurface.cxx:940-1013)"))?;
    // `cxx:1058-1074`. `ConvertToPeriodic(false)` uses preci -1 (`ShapeCustom_Surface.hxx:67`).
    let force_u = closed_u && (periodic_u || surf.is_u_closed());
    let force_v = closed_v && (periodic_v || surf.is_v_closed());
    if force_u || force_v {
        let wrapped: Arc<dyn Surface> = Arc::new(surf.clone());
        if let Some(periodic) = convert_to_periodic(&wrapped, -1.0) {
            if let Some(bs) = periodic.osculating_bspline() {
                surf = bs;
            }
        }
    }
    Ok(surf)
}

fn weight_grid(
    weights: &[Vec<f64>],
    keep_u: &[usize],
    keep_v: &[usize],
) -> Result<Vec<Vec<f64>>, String> {
    if keep_u.is_empty() || keep_v.is_empty() {
        return Ok(vec![vec![1.0]]);
    }
    let mut grid = Vec::with_capacity(keep_u.len());
    for &iu in keep_u {
        let mut row = Vec::with_capacity(keep_v.len());
        for &iv in keep_v {
            let w = weights
                .get(iu)
                .and_then(|r| r.get(iv))
                .copied()
                .unwrap_or(0.0);
            if w < PCONFUSION {
                return Err(
                    "IGES 128 weight is not positive: IGESToBRep_BasicSurface.cxx:959-963".into(),
                );
            }
            row.push(w);
        }
        grid.push(row);
    }
    Ok(grid)
}
