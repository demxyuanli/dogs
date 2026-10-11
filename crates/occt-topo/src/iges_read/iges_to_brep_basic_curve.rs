//! `IGESToBRep_BasicCurve` for entities 100, 110, 124, and 126.
//!
//! `TransferBasicCurve` (`IGESToBRep_BasicCurve.cxx:164-233`) scales the 3D
//! result by the unit factor. `Transfer2dBasicCurve` (`cxx:238-298`) does not.

use std::f64::consts::PI;
use std::sync::Arc;

use occt_core::elib::{clib, clib2d};
use occt_core::gp::{
    GpAx22d, GpAx2, GpAx3, GpCirc, GpCirc2d, GpDir, GpDir2d, GpGTrsf, GpPnt, GpPnt2d, GpTrsf,
    GpVec, GpXyz,
};
use occt_core::precision::{CONFUSION, PCONFUSION};
use occt_geom::curve::Curve;
use occt_geom::bspline_curve::MAX_DEGREE;
use occt_geom::{GeomBSplineCurve, GeomCircle, GeomTrimmedCurveBasis};
use occt_geom2d::curve::Curve2d;
use occt_geom2d::trimmed::Geom2dTrimmedCurve;
use occt_geom2d::{Geom2dBSplineCurve, Geom2dCircle, Geom2dLine};

use crate::builder::TopoBuilder;
use crate::shape::TopoShape;

use super::iges_data_iges_model::{de_pointer_to_entity, EntityData, IgesModel};

pub(in crate::iges_read) struct CurveSpan {
    pub curve: Arc<dyn Curve>,
    pub first: f64,
    pub last: f64,
}

pub(in crate::iges_read) struct Curve2dSpan {
    pub curve: Arc<dyn Curve2d>,
    pub first: f64,
    pub last: f64,
}

/// `IGESToBRep_Actor.cxx:171-190`: resolution above `1e-8` replaces the
/// `1e-4` default of `IGESToBRep_CurveAndSurface.cxx:58-62`.
pub(in crate::iges_read) fn eps_geom(resolution: f64) -> f64 {
    if resolution > 1.0e-8 {
        resolution
    } else {
        1.0e-4
    }
}

/// `IGESGeom_TransformationMatrix::Value` (`cxx:72-87`).
pub(in crate::iges_read) fn matrix_gtrsf(rows: &[[f64; 4]; 3]) -> GpGTrsf {
    let mut g = GpGTrsf::identity();
    for (i, row) in rows.iter().enumerate() {
        for (j, value) in row.iter().enumerate() {
            g.set_value(i + 1, j + 1, *value);
        }
    }
    g
}

pub(in crate::iges_read) fn entity_location(
    model: &IgesModel,
    index: usize,
) -> Result<Option<GpGTrsf>, String> {
    let ent = model.entity(index)?;
    if ent.transform == 0 {
        return Ok(None);
    }
    let ix = de_pointer_to_entity(ent.transform)?
        .ok_or_else(|| format!("IGES entity {index} transform pointer is null"))?;
    match &model.entity(ix)?.data {
        EntityData::TransformMatrix { rows } => Ok(Some(matrix_gtrsf(rows))),
        _ => Err(format!(
            "IGES entity {index} transform is not entity 124 (IGESGeom_TransformationMatrix.cxx:72)"
        )),
    }
}

pub(in crate::iges_read) fn apply_gtrsf(loc: &GpGTrsf, p: GpPnt) -> GpPnt {
    let mut xyz = *p.xyz();
    loc.transforms(&mut xyz);
    GpPnt::from_xyz(&xyz)
}

fn dir_of(xyz: GpXyz) -> Result<GpDir, String> {
    GpDir::from_xyz(&xyz).map_err(|_| "IGES direction is null".to_string())
}

/// `gp_Ax3` main direction (`gp_Ax3.hxx:89`): `N ^ (Vx ^ N)`, the part of `Vx`
/// perpendicular to `N`.
fn ax3_main_direction(normal: &GpDir, x_dir: &GpDir) -> Option<GpDir> {
    let n = normal.xyz();
    let proj = x_dir.xyz().subtracted(&n.multiplied(n.dot(x_dir.xyz())));
    GpDir::from_xyz(&proj).ok()
}

pub(in crate::iges_read) fn transformed_dir(loc: &GpGTrsf, xyz: GpXyz) -> Result<GpDir, String> {
    let mut rot = *loc;
    rot.set_translation_part(&GpXyz::zero());
    let mut v = xyz;
    rot.transforms(&mut v);
    dir_of(v)
}

/// `IGESData_ToolLocation::ConvertLocation` (`IGESData_ToolLocation.cxx:290-358`).
pub(in crate::iges_read) fn convert_location(
    prec: f64,
    loc: &GpGTrsf,
    unit: f64,
) -> Option<GpTrsf> {
    let mut result = GpTrsf::identity();
    let mut v1 = GpXyz::new(loc.value(1, 1), loc.value(1, 2), loc.value(1, 3));
    let mut v2 = GpXyz::new(loc.value(2, 1), loc.value(2, 2), loc.value(2, 3));
    let mut v3 = GpXyz::new(loc.value(3, 1), loc.value(3, 2), loc.value(3, 3));
    let m1 = v1.modulus();
    let m2 = v2.modulus();
    let m3 = v3.modulus();
    if m1 < prec || m2 < prec || m3 < prec {
        return None;
    }
    let mm = (m1 + m2 + m3) / 3.0;
    if (m1 - mm).abs() > prec * mm || (m2 - mm).abs() > prec * mm || (m3 - mm).abs() > prec * mm {
        return None;
    }
    v1 = v1.divide(m1);
    v2 = v2.divide(m2);
    v3 = v3.divide(m3);
    if v1.dot(&v2).abs() > prec || v2.dot(&v3).abs() > prec || v3.dot(&v1).abs() > prec {
        return None;
    }
    if (mm - 1.0).abs() > prec {
        result.set_scale(&GpPnt::new(0.0, 0.0, 0.0), mm).ok()?;
    }
    let mut tp = loc.translation_part();
    if unit != 1.0 {
        tp = tp.multiply(unit);
    }
    if tp.x() != 0.0 || tp.y() != 0.0 || tp.z() != 0.0 {
        result.set_translation_part(&GpVec::from_xyz(&tp));
    }
    if v1.x() != 1.0
        || v1.y() != 0.0
        || v1.z() != 0.0
        || v2.x() != 0.0
        || v2.y() != 1.0
        || v2.z() != 0.0
        || v3.x() != 0.0
        || v3.y() != 0.0
        || v3.z() != 1.0
    {
        let d1 = dir_of(v1).ok()?;
        let d2 = dir_of(v2).ok()?;
        let d3 = dir_of(v3).ok()?;
        // `gp_Ax3(P, N, Vx)` (`gp_Ax3.hxx:84-90`) sets X to `N ^ (Vx ^ N)`.
        // A 9-digit column can fail `GpAx3::new`'s `IsNormal` (1e-12) and still
        // pass the 1e-4 test above, which dropped the corner placement.
        let vx = ax3_main_direction(&d3, &d1)?;
        let mut axes = GpAx3::new(GpPnt::new(0.0, 0.0, 0.0), d3, &vx).ok()?;
        let cross = d3.crossed(&d1).ok()?;
        if cross.xyz().dot(d2.xyz()) < 0.0 {
            axes.y_reverse();
        }
        let mut transf = GpTrsf::identity();
        transf.set_transformation(&axes);
        result.multiply(&transf);
    }
    Some(result)
}

pub(in crate::iges_read) fn scale_curve(curve: Arc<dyn Curve>, unit: f64) -> Result<Arc<dyn Curve>, String> {
    // `TransferBasicCurve` (`cxx:231`): `res->Scale(origin, GetUnitFactor())`.
    if unit == 1.0 {
        return Ok(curve);
    }
    let mut t = GpTrsf::identity();
    t.set_scale(&GpPnt::new(0.0, 0.0, 0.0), unit)
        .map_err(|e| e.to_string())?;
    Ok(Arc::from(curve.transformed(&t)))
}

pub(in crate::iges_read) fn edge_from_span(span: &CurveSpan) -> TopoShape {
    let b = TopoBuilder::new();
    let p1 = span.curve.d0(span.first);
    let p2 = span.curve.d0(span.last);
    let mut edge = b.make_edge(span.curve.clone(), span.first, span.last);
    let v1 = b.make_vertex(p1, 0.0);
    let v2 = b.make_vertex(p2, 0.0);
    b.add_edge_vertices(&mut edge, &v1, &v2);
    edge.0
}

/// `IGESGeom_Line::TransformedStartPoint` / `TransformedEndPoint`
/// (`IGESGeom_Line.cxx:55-80`). `Location().Transforms` runs when the
/// directory transform pointer is non-zero. The points are then scaled from
/// the origin by the unit factor (`IGESToBRep_TopoSurface.cxx:760-761`).
pub(in crate::iges_read) fn line_transformed_ends(
    model: &IgesModel,
    index: usize,
) -> Result<(GpPnt, GpPnt), String> {
    let (start, end) = match &model.entity(index)?.data {
        EntityData::Line { start, end } => (*start, *end),
        _ => return Err("IGES 110 has no line parameters".into()),
    };
    let mut ps = GpPnt::new(start[0], start[1], start[2]);
    let mut pe = GpPnt::new(end[0], end[1], end[2]);
    if let Some(g) = entity_location(model, index)? {
        ps = apply_gtrsf(&g, ps);
        pe = apply_gtrsf(&g, pe);
    }
    let unit = model.unit_factor;
    ps = GpPnt::new(ps.x() * unit, ps.y() * unit, ps.z() * unit);
    pe = GpPnt::new(pe.x() * unit, pe.y() * unit, pe.z() * unit);
    Ok((ps, pe))
}

pub(in crate::iges_read) fn transfer_line_span(
    model: &IgesModel,
    index: usize,
) -> Result<CurveSpan, String> {
    // `TransferLine` (`IGESToBRep_BasicCurve.cxx:1372-1418`) then unit scale.
    let (ps, pe) = line_transformed_ends(model, index)?;
    if ps.distance(&pe) <= CONFUSION {
        return Err(
            "IGES 110 start and end coincide: IGESToBRep_BasicCurve.cxx:1413-1417".into(),
        );
    }
    let dir = GpDir::from_vec(&GpVec::from_pnts(&ps, &pe))
        .map_err(|_| "IGES 110 direction is null".to_string())?;
    let line = occt_geom::GeomLine::new(occt_core::gp::GpLin::from_pnt_dir(ps, dir));
    let t1 = clib::line_parameter(line.lin(), &ps);
    let t2 = clib::line_parameter(line.lin(), &pe);
    Ok(CurveSpan {
        curve: Arc::new(line),
        first: t1,
        last: t2,
    })
}

pub(in crate::iges_read) fn transfer_2d_line(
    model: &IgesModel,
    index: usize,
) -> Result<Curve2dSpan, String> {
    // `Transfer2dLine` (`cxx:1425-1472`). No unit scale.
    let (start, end) = match &model.entity(index)?.data {
        EntityData::Line { start, end } => (*start, *end),
        _ => return Err("IGES 110 has no line parameters".into()),
    };
    let loc = entity_location(model, index)?;
    let (beg, end_p) = if let Some(g) = loc.as_ref() {
        let ps = apply_gtrsf(g, GpPnt::new(start[0], start[1], start[2]));
        let pe = apply_gtrsf(g, GpPnt::new(end[0], end[1], end[2]));
        (GpPnt2d::new(ps.x(), ps.y()), GpPnt2d::new(pe.x(), pe.y()))
    } else {
        (
            GpPnt2d::new(start[0], start[1]),
            GpPnt2d::new(end[0], end[1]),
        )
    };
    if beg.distance(&end_p) <= PCONFUSION {
        return Err(
            "IGES 110 2d start and end coincide: IGESToBRep_BasicCurve.cxx:1467-1471".into(),
        );
    }
    let dir = GpDir2d::from_xy(&occt_core::gp::GpXY::new(end_p.x() - beg.x(), end_p.y() - beg.y()))
        .map_err(|_| "IGES 110 2d direction is null".to_string())?;
    let line = Geom2dLine::from_pnt_dir(beg, dir);
    let t1 = clib2d::line_parameter_ax2d(line.lin().position(), &beg);
    let t2 = clib2d::line_parameter_ax2d(line.lin().position(), &end_p);
    Ok(Curve2dSpan {
        curve: Arc::new(Geom2dTrimmedCurve::new(Arc::new(line), t1, t2)),
        first: t1,
        last: t2,
    })
}

fn arc_radius(center: [f64; 2], start: [f64; 2]) -> f64 {
    // `IGESGeom_CircularArc::Radius` (`cxx:100-109`).
    let dx = center[0] - start[0];
    let dy = center[1] - start[1];
    (dx * dx + dy * dy).sqrt()
}

fn arc_closed(start: [f64; 2], end: [f64; 2]) -> bool {
    // `IGESGeom_CircularArc::IsClosed` (`cxx:148-151`).
    (start[0] - end[0]).abs() < PCONFUSION && (start[1] - end[1]).abs() < PCONFUSION
}

fn trim_angles(closed: bool, t1: f64, mut t2: f64, eps: f64, chord: f64, radius: f64) -> Option<(f64, f64)> {
    // `TransferCircularArc` (`cxx:770-787`) and the 2d copy (`cxx:855-872`).
    if closed && t1 >= eps {
        t2 = t1 + 2.0 * PI;
    }
    if !closed && (t1 - t2).abs() <= PCONFUSION {
        t2 = t1 + chord / radius;
    }
    if !closed || t1 >= eps {
        if t2 < t1 {
            t2 += 2.0 * PI;
        }
        Some((t1, t2))
    } else {
        None
    }
}

pub(in crate::iges_read) fn transfer_circular_arc(
    model: &IgesModel,
    index: usize,
) -> Result<CurveSpan, String> {
    // `TransferCircularArc` (`cxx:717-789`) then unit scale (`cxx:231`).
    let (zt, center, start, end) = match &model.entity(index)?.data {
        EntityData::CircularArc {
            zt,
            center,
            start,
            end,
        } => (*zt, *center, *start, *end),
        _ => return Err("IGES 100 has no circular-arc parameters".into()),
    };
    let radius = arc_radius(center, start);
    if radius <= CONFUSION {
        return Err("IGES 100 radius is null: IGESToBRep_BasicCurve.cxx:761".into());
    }
    let loc = entity_location(model, index)?;
    let (frame_origin, norm, main, start_p, end_p) = if let Some(g) = loc.as_ref() {
        let norm = transformed_dir(g, GpXyz::new(0.0, 0.0, 1.0))?;
        let main = transformed_dir(g, GpXyz::new(1.0, 0.0, 0.0))?;
        let origin = apply_gtrsf(g, GpPnt::new(center[0], center[1], zt));
        let sp = apply_gtrsf(g, GpPnt::new(start[0], start[1], zt));
        let ep = apply_gtrsf(g, GpPnt::new(end[0], end[1], zt));
        (origin, norm, main, sp, ep)
    } else {
        (
            GpPnt::new(center[0], center[1], zt),
            GpDir::from_axis(occt_core::gp::dir::DirAxis::Z),
            GpDir::from_axis(occt_core::gp::dir::DirAxis::X),
            GpPnt::new(start[0], start[1], zt),
            GpPnt::new(end[0], end[1], zt),
        )
    };
    let frame = GpAx2::new(frame_origin, norm, main)
        .map_err(|_| "IGES 100 frame is singular: gp_Ax2".to_string())?;
    let circ = GpCirc::new(frame, radius);
    let t1 = clib::circle_parameter(&frame, &start_p);
    let t2 = clib::circle_parameter(&frame, &end_p);
    let eps = eps_geom(model.resolution);
    let closed = arc_closed(start, end);
    let chord = start_p.distance(&end_p);
    let (first, last) = match trim_angles(closed, t1, t2, eps, chord, radius) {
        Some(pair) => pair,
        None => (0.0, 2.0 * PI),
    };
    let curve = scale_curve(Arc::new(GeomCircle::new(circ)), model.unit_factor)?;
    Ok(CurveSpan { curve, first, last })
}

pub(in crate::iges_read) fn transfer_2d_circular_arc(
    model: &IgesModel,
    index: usize,
) -> Result<Curve2dSpan, String> {
    // `Transfer2dCircularArc` (`cxx:794-873`).
    let (zt, center, start, end) = match &model.entity(index)?.data {
        EntityData::CircularArc {
            zt,
            center,
            start,
            end,
        } => (*zt, *center, *start, *end),
        _ => return Err("IGES 100 has no circular-arc parameters".into()),
    };
    let radius = arc_radius(center, start);
    if radius <= PCONFUSION {
        return Err("IGES 100 2d radius is null: IGESToBRep_BasicCurve.cxx:829".into());
    }
    let axis_z = GpDir::from_axis(occt_core::gp::dir::DirAxis::Z);
    let loc = entity_location(model, index)?;
    let mut apply = false;
    let mut reversed = false;
    if let Some(g) = loc.as_ref() {
        let transformed = transformed_dir(g, GpXyz::new(0.0, 0.0, 1.0))?;
        // `SetEpsilon(1.E-03)` then `IsParallel` (`cxx:808-816`).
        if transformed.is_parallel_tol(&axis_z, 1.0e-3) {
            apply = true;
            if transformed.angle(&axis_z) > PI - 1.0e-3 {
                reversed = true;
            }
        }
    }
    let (cxy, main, sp, ep) = if apply {
        let g = loc.as_ref().expect("location");
        let c = apply_gtrsf(g, GpPnt::new(center[0], center[1], zt));
        let main = transformed_dir(g, GpXyz::new(1.0, 0.0, 0.0))?;
        let sp = apply_gtrsf(g, GpPnt::new(start[0], start[1], zt));
        let ep = apply_gtrsf(g, GpPnt::new(end[0], end[1], zt));
        (
            GpPnt2d::new(c.x(), c.y()),
            GpDir2d::new(main.x(), main.y()).map_err(|_| "IGES 100 2d axis is null".to_string())?,
            GpPnt2d::new(sp.x(), sp.y()),
            GpPnt2d::new(ep.x(), ep.y()),
        )
    } else {
        (
            GpPnt2d::new(center[0], center[1]),
            GpDir2d::new(1.0, 0.0).expect("x"),
            GpPnt2d::new(start[0], start[1]),
            GpPnt2d::new(end[0], end[1]),
        )
    };
    let ax = GpAx22d::from_xdir(cxy, main);
    let mut circle = Geom2dCircle::new(GpCirc2d::new(ax, radius));
    if reversed {
        circle.reverse();
    }
    let circ = *circle.circ();
    let t1 = clib2d::parameter_circ2d(&circ, &sp);
    let t2 = clib2d::parameter_circ2d(&circ, &ep);
    let eps = eps_geom(model.resolution);
    let closed = arc_closed(start, end);
    let chord = sp.distance(&ep);
    let basis: Arc<dyn Curve2d> = Arc::new(circle);
    let bounds = trim_angles(closed, t1, t2, eps, chord, radius);
    let (first, last) = bounds.unwrap_or((0.0, 2.0 * PI));
    let curve: Arc<dyn Curve2d> = if bounds.is_none() {
        basis
    } else {
        Arc::new(Geom2dTrimmedCurve::new(basis, first, last))
    };
    Ok(Curve2dSpan { curve, first, last })
}

fn real_epsilon(x: f64) -> f64 {
    let a = x.abs();
    if a > 1.0 {
        f64::EPSILON * a
    } else {
        f64::EPSILON
    }
}

struct PreparedSpline {
    poles: Vec<GpPnt>,
    weights: Option<Vec<f64>>,
    knots: Vec<f64>,
    mults: Vec<i32>,
    degree: usize,
    u_min: f64,
    u_max: f64,
}

fn prepare_bspline(model: &IgesModel, index: usize) -> Result<PreparedSpline, String> {
    // `TransferBSplineCurve` (`cxx:976-1203`).
    let data = match &model.entity(index)?.data {
        EntityData::BSplineCurve {
            degree,
            polynomial,
            knots,
            weights,
            poles,
            u_min,
            u_max,
            ..
        } => (
            *degree,
            *polynomial,
            knots.clone(),
            weights.clone(),
            poles.clone(),
            *u_min,
            *u_max,
        ),
        _ => return Err("IGES 126 has no B-spline parameters".into()),
    };
    let (degree_i, polynomial, knots, weights, poles, u_min, u_max) = data;
    if degree_i <= 0 || degree_i as usize > MAX_DEGREE {
        return Err(
            "IGES 126 degree is out of range: IGESToBRep_BasicCurve.cxx:992-997".into(),
        );
    }
    let degree = degree_i as usize;
    if poles.len() < 2 {
        return Err("IGES 126 has fewer than 2 poles: IGESToBRep_BasicCurve.cxx:1006-1011".into());
    }
    let loc = entity_location(model, index)?;
    let mut gp_poles: Vec<GpPnt> = Vec::with_capacity(poles.len());
    for p in &poles {
        let mut q = GpPnt::new(p[0], p[1], p[2]);
        if let Some(g) = loc.as_ref() {
            q = apply_gtrsf(g, q);
        }
        gp_poles.push(q);
    }
    if knots.is_empty() {
        return Err("IGES 126 has no knots".into());
    }
    let mut temp_k = vec![knots[0]];
    let mut temp_m = vec![1i32];
    for i in 1..knots.len() {
        let k1 = knots[i];
        if (k1 - knots[i - 1]).abs() <= real_epsilon(k1) {
            let last = temp_m.len() - 1;
            temp_m[last] += 1;
        } else {
            temp_k.push(k1);
            temp_m.push(1);
        }
    }
    let mut new_nb = gp_poles.len();
    let mut seq_index: Vec<usize> = Vec::new();
    let mut old_sum = 0i32;
    let mut mults = Vec::with_capacity(temp_k.len());
    for i in 0..temp_k.len() {
        let mut a_mult = temp_m[i];
        let max_mult = if i == 0 || i + 1 == temp_k.len() {
            degree as i32 + 1
        } else {
            degree as i32
        };
        if a_mult > max_mult {
            let mut del = old_sum + 1;
            while a_mult > max_mult {
                new_nb = new_nb.saturating_sub(1);
                seq_index.push(del as usize);
                del += 1;
                a_mult -= 1;
            }
        }
        old_sum += temp_m[i];
        mults.push(a_mult);
    }
    for i in 1..temp_k.len() {
        if temp_k[i] < temp_k[i - 1] {
            return Err(
                "IGES 126 knots are not ascending: IGESToBRep_BasicCurve.cxx:117-124".into(),
            );
        }
    }
    let mut keep: Vec<usize> = (0..gp_poles.len()).collect();
    if new_nb < gp_poles.len() {
        let mut offset = 0usize;
        for del in &seq_index {
            let idx = del.saturating_sub(1 + offset);
            if idx < keep.len() {
                keep.remove(idx);
                offset += 1;
            }
        }
    }
    if keep.len() < 2 {
        return Err("IGES 126 pole correction left fewer than 2 poles".into());
    }
    let poles_out: Vec<GpPnt> = keep.iter().map(|i| gp_poles[*i]).collect();
    let weights_out = if polynomial {
        None
    } else {
        let mut w = Vec::with_capacity(keep.len());
        for i in &keep {
            let weight = weights.get(*i).copied().unwrap_or(0.0);
            if weight < PCONFUSION {
                return Err(
                    "IGES 126 weight is not positive: IGESToBRep_BasicCurve.cxx:1169-1174".into(),
                );
            }
            w.push(weight);
        }
        Some(w)
    };
    Ok(PreparedSpline {
        poles: poles_out,
        weights: weights_out,
        knots: temp_k,
        mults,
        degree,
        u_min,
        u_max,
    })
}

fn segment_bounds(first: f64, last: f64, u_min: f64, u_max: f64) -> Option<(f64, f64)> {
    // `cxx:1251-1271`.
    if (u_min - first) > -PCONFUSION && (last - u_max) > -PCONFUSION && u_min <= u_max {
        Some((u_min, u_max))
    } else {
        None
    }
}

pub(in crate::iges_read) fn transfer_bspline(
    model: &IgesModel,
    index: usize,
) -> Result<CurveSpan, String> {
    let prep = prepare_bspline(model, index)?;
    let unit = model.unit_factor;
    let poles: Vec<GpPnt> = prep
        .poles
        .iter()
        .map(|p| GpPnt::new(p.x() * unit, p.y() * unit, p.z() * unit))
        .collect();
    let mut bs = match &prep.weights {
        Some(w) => {
            let flat = occt_core::bspl::banded_interp::knot_sequence(
                &prep.knots,
                &prep.mults,
                prep.degree as i32,
            );
            GeomBSplineCurve::rational(poles, w.clone(), flat, prep.degree)
        }
        None => GeomBSplineCurve::from_poles_knots_mults(
            poles,
            prep.knots.clone(),
            prep.mults.clone(),
            prep.degree,
        ),
    }
    .map_err(|e| format!("IGES 126 curve: {e} (IGESToBRep_BasicCurve.cxx:1215-1222)"))?;
    let first = bs.first_parameter();
    let last = bs.last_parameter();
    if let Some((u0, u1)) = segment_bounds(first, last, prep.u_min, prep.u_max) {
        if (u1 - u0).abs() > PCONFUSION && bs.segment(u0, u1, PCONFUSION).is_err() {
            let curve = Arc::new(bs);
            return Ok(CurveSpan {
                curve: Arc::new(GeomTrimmedCurveBasis::new(curve, u0, u1)),
                first: u0,
                last: u1,
            });
        }
        let f = bs.first_parameter();
        let l = bs.last_parameter();
        return Ok(CurveSpan {
            curve: Arc::new(bs),
            first: f,
            last: l,
        });
    }
    let f = bs.first_parameter();
    let l = bs.last_parameter();
    Ok(CurveSpan {
        curve: Arc::new(bs),
        first: f,
        last: l,
    })
}

pub(in crate::iges_read) fn transfer_2d_bspline(
    model: &IgesModel,
    index: usize,
) -> Result<Curve2dSpan, String> {
    // `Transfer2dBSplineCurve` (`cxx:1279-1367`): 3D transfer without the
    // outer unit scale, then XY poles.
    let prep = prepare_bspline(model, index)?;
    let xs: Vec<f64> = prep.poles.iter().map(|p| p.x()).collect();
    let ys: Vec<f64> = prep.poles.iter().map(|p| p.y()).collect();
    let flat = occt_core::bspl::banded_interp::knot_sequence(
        &prep.knots,
        &prep.mults,
        prep.degree as i32,
    );
    let mut bs = match &prep.weights {
        Some(w) => Geom2dBSplineCurve::rational(xs, ys, w.clone(), flat, prep.degree),
        None => Geom2dBSplineCurve::new(xs, ys, flat, prep.degree),
    }
    .map_err(|e| format!("IGES 126 2d curve: {e} (IGESToBRep_BasicCurve.cxx:1349-1356)"))?;
    let first = bs.first_parameter();
    let last = bs.last_parameter();
    if let Some((u0, u1)) = segment_bounds(first, last, prep.u_min, prep.u_max) {
        if (u1 - u0).abs() > PCONFUSION && bs.segment(u0, u1, PCONFUSION).is_err() {
            let curve: Arc<dyn Curve2d> =
                Arc::new(Geom2dTrimmedCurve::new(Arc::new(bs), u0, u1));
            return Ok(Curve2dSpan {
                curve,
                first: u0,
                last: u1,
            });
        }
        let f = bs.first_parameter();
        let l = bs.last_parameter();
        return Ok(Curve2dSpan {
            curve: Arc::new(bs),
            first: f,
            last: l,
        });
    }
    let f = bs.first_parameter();
    let l = bs.last_parameter();
    Ok(Curve2dSpan {
        curve: Arc::new(bs),
        first: f,
        last: l,
    })
}

pub(in crate::iges_read) fn basic_curve_3d(
    model: &IgesModel,
    index: usize,
) -> Result<CurveSpan, String> {
    match model.entity(index)?.type_number {
        110 => transfer_line_span(model, index),
        100 => transfer_circular_arc(model, index),
        104 => super::conic_arc::transfer_conic_arc(model, index),
        106 => super::copious_data::transfer_copious_data(model, index),
        112 => super::spline_curve::transfer_spline_curve(model, index),
        126 => transfer_bspline(model, index),
        n => Err(format!(
            "UNPORTED basic curve {n}: IGESToBRep_BasicCurve.cxx:164"
        )),
    }
}

pub(in crate::iges_read) fn basic_curve_2d(
    model: &IgesModel,
    index: usize,
) -> Result<Curve2dSpan, String> {
    match model.entity(index)?.type_number {
        110 => transfer_2d_line(model, index),
        100 => transfer_2d_circular_arc(model, index),
        104 => super::conic_arc::transfer_2d_conic_arc(model, index),
        106 => super::copious_data::transfer_2d_copious_data(model, index),
        112 => super::spline_curve::transfer_2d_spline_curve(model, index),
        126 => transfer_2d_bspline(model, index),
        n => Err(format!(
            "UNPORTED 2d basic curve {n}: IGESToBRep_BasicCurve.cxx:238"
        )),
    }
}
