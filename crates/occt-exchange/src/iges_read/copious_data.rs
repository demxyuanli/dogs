//! `IGESToBRep_BasicCurve::TransferCopiousData` (`cxx:1504-1736`) and
//! `IGESConvGeom::IncreaseCurveContinuity` (`IGESConvGeom.cxx:174-226`).
//!
//! `read.iges.bspline.continuity` defaults to 1 (`IGESData.cxx:66`). A form
//! other than 11, 12, or 63 warns `IGES_1240` and still transfers the points.
//! `TransferBasicCurve` scales the 3D curve afterwards (`cxx:231`).
//! `GetUVResolution` stays 0 when no surface is bound (`CurveAndSurface.cxx:678-686`).

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpPnt2d};
use occt_core::precision::{CONFUSION, REAL_SMALL};
use occt_geom::bspline_curve::GeomBSplineCurve;
use occt_geom::curve::Curve;
use occt_geom2d::bspline_curve::Geom2dBSplineCurve;
use occt_geom2d::curve::Curve2d;

use super::iges_data_iges_model::{EntityData, IgesModel};
use super::iges_to_brep_basic_curve::{
    apply_gtrsf, entity_location, eps_geom, scale_curve, Curve2dSpan, CurveSpan,
};

const CONTINUITY: i32 = 1;

pub(in crate::iges_read) fn transfer_copious_data(
    model: &IgesModel,
    index: usize,
) -> Result<CurveSpan, String> {
    let points = points_3d(model, index)?;
    let mut curve = polyline_3d(&points)?;
    let eps = (eps_geom(model.resolution) / 10.0).max(CONFUSION);
    increase_continuity_3d(&mut curve, eps, CONTINUITY)?;
    let curve = scale_curve(Arc::new(curve), model.unit_factor)?;
    let first = curve.first_parameter();
    let last = curve.last_parameter();
    Ok(CurveSpan { curve, first, last })
}

pub(in crate::iges_read) fn transfer_2d_copious_data(
    model: &IgesModel,
    index: usize,
) -> Result<Curve2dSpan, String> {
    let points = points_2d(model, index)?;
    let mut curve = polyline_2d(&points)?;
    increase_continuity_2d(&mut curve, CONFUSION, CONTINUITY)?;
    let first = curve.first_parameter();
    let last = curve.last_parameter();
    Ok(Curve2dSpan {
        curve: Arc::new(curve),
        first,
        last,
    })
}

fn points_3d(model: &IgesModel, index: usize) -> Result<Vec<GpPnt>, String> {
    let (data_type, z_plane, tuples) = load(model, index)?;
    let nb = nb_points(data_type, tuples.len())?;
    if nb < 2 {
        return Err("IGES_1195: IGESToBRep_BasicCurve.cxx:1526-1530".into());
    }
    let loc = entity_location(model, index)?;
    let mut poles = Vec::new();
    for i in 0..nb {
        let mut p = point_at(data_type, z_plane, &tuples, i)?;
        if let Some(g) = loc.as_ref() {
            p = apply_gtrsf(g, p);
        }
        if poles.last().is_none_or(|prev: &GpPnt| prev.distance(&p) > REAL_SMALL) {
            poles.push(p);
        }
    }
    if poles.len() == 1 {
        return Err("IGES_1235: IGESToBRep_BasicCurve.cxx:1574-1579".into());
    }
    Ok(poles)
}

fn points_2d(model: &IgesModel, index: usize) -> Result<Vec<GpPnt2d>, String> {
    let (data_type, z_plane, tuples) = load(model, index)?;
    let nb = nb_points(data_type, tuples.len())?;
    if nb < 2 {
        return Err("IGES_1195: IGESToBRep_BasicCurve.cxx:1644-1648".into());
    }
    let loc = entity_location(model, index)?;
    let mut poles = Vec::new();
    for i in 0..nb {
        let mut p = point_at(data_type, z_plane, &tuples, i)?;
        if let Some(g) = loc.as_ref() {
            p = apply_gtrsf(g, p);
        }
        let uv = GpPnt2d::new(p.x(), p.y());
        if poles.last().is_none_or(|prev: &GpPnt2d| prev.distance(&uv) > REAL_SMALL) {
            poles.push(uv);
        }
    }
    if poles.len() == 1 {
        return Err("IGES_1235: IGESToBRep_BasicCurve.cxx:1691-1696".into());
    }
    Ok(poles)
}

fn polyline_3d(poles: &[GpPnt]) -> Result<GeomBSplineCurve, String> {
    let (knots, mults) = chord_knots_3d(poles);
    GeomBSplineCurve::from_poles_knots_mults(poles.to_vec(), knots, mults, 1)
        .map_err(|e| format!("IGES 106 curve: {e} (IGESToBRep_BasicCurve.cxx:1611)"))
}

fn polyline_2d(poles: &[GpPnt2d]) -> Result<Geom2dBSplineCurve, String> {
    let (knots, mults) = chord_knots_2d(poles);
    let flat = occt_core::bspl::banded_interp::knot_sequence(&knots, &mults, 1);
    let xs: Vec<f64> = poles.iter().map(|p| p.x()).collect();
    let ys: Vec<f64> = poles.iter().map(|p| p.y()).collect();
    Geom2dBSplineCurve::new(xs, ys, flat, 1)
        .map_err(|e| format!("IGES 106 2d curve: {e} (IGESToBRep_BasicCurve.cxx:1728)"))
}

fn chord_knots_3d(poles: &[GpPnt]) -> (Vec<f64>, Vec<i32>) {
    let mut knots = Vec::with_capacity(poles.len());
    knots.push(0.0);
    for i in 1..poles.len() {
        let prev = knots[i - 1];
        knots.push(prev + poles[i].distance(&poles[i - 1]));
    }
    (knots, end_mults(poles.len()))
}

fn chord_knots_2d(poles: &[GpPnt2d]) -> (Vec<f64>, Vec<i32>) {
    let mut knots = Vec::with_capacity(poles.len());
    knots.push(0.0);
    for i in 1..poles.len() {
        let prev = knots[i - 1];
        knots.push(prev + poles[i].distance(&poles[i - 1]));
    }
    (knots, end_mults(poles.len()))
}

fn end_mults(n: usize) -> Vec<i32> {
    let mut mults = vec![1; n];
    if let Some(first) = mults.first_mut() {
        *first = 2;
    }
    if let Some(last) = mults.last_mut() {
        *last = 2;
    }
    mults
}

/// `IGESConvGeom::IncreaseCurveContinuity` (`cxx:174-226`).
pub(in crate::iges_read) fn increase_continuity_3d(
    curve: &mut GeomBSplineCurve,
    eps: f64,
    continuity: i32,
) -> Result<(), String> {
    if continuity < 1 {
        return Ok(());
    }
    loop {
        let degree = curve.degree() as i32;
        let (_, mults) = curve.distinct_knots_and_mults();
        let first = occt_core::bspl::locate::first_u_knot_index(degree, &mults).max(1);
        let mut i = first + 1;
        let mut modified = false;
        loop {
            let (_, mults) = curve.distinct_knots_and_mults();
            let last = occt_core::bspl::locate::last_u_knot_index(degree, &mults).max(1);
            if i >= last {
                break;
            }
            let mult = mults[(i - 1) as usize];
            if degree - mult < continuity {
                let target = if continuity >= 2 { degree - 2 } else { degree - 1 };
                let removed = curve.remove_knot(i, target, eps).map_err(|e| e.to_string())?;
                if removed {
                    modified = true;
                } else if continuity >= 2 {
                    let removed = curve
                        .remove_knot(i, degree - 1, eps)
                        .map_err(|e| e.to_string())?;
                    modified |= removed;
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

fn increase_continuity_2d(
    curve: &mut Geom2dBSplineCurve,
    eps: f64,
    continuity: i32,
) -> Result<(), String> {
    if continuity < 1 {
        return Ok(());
    }
    loop {
        let degree = curve.degree() as i32;
        let (_, mults) = curve.distinct_knots_and_mults();
        let first = occt_core::bspl::locate::first_u_knot_index(degree, &mults).max(1);
        let mut i = first + 1;
        let mut modified = false;
        loop {
            let (_, mults) = curve.distinct_knots_and_mults();
            let last = occt_core::bspl::locate::last_u_knot_index(degree, &mults).max(1);
            if i >= last {
                break;
            }
            let mult = mults[(i - 1) as usize];
            if degree - mult < continuity {
                let target = if continuity >= 2 { degree - 2 } else { degree - 1 };
                let removed = curve.remove_knot(i, target, eps).map_err(|e| e.to_string())?;
                if removed {
                    modified = true;
                } else if continuity >= 2 {
                    let removed = curve
                        .remove_knot(i, degree - 1, eps)
                        .map_err(|e| e.to_string())?;
                    modified |= removed;
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

fn load(model: &IgesModel, index: usize) -> Result<(i32, f64, Vec<f64>), String> {
    let ent = model.entity(index)?;
    match &ent.data {
        EntityData::Copious {
            data_type,
            z_plane,
            tuples,
        } => Ok((*data_type, *z_plane, tuples.clone())),
        _ => Err("IGES 106 has no copious-data parameters".into()),
    }
}

fn nb_points(data_type: i32, len: usize) -> Result<usize, String> {
    let stride = if data_type == 1 {
        2
    } else if data_type == 2 {
        3
    } else if data_type == 3 {
        6
    } else {
        return Err(format!(
            "IGES 106 data type {data_type}: IGESGeom_CopiousData.cxx:103-114"
        ));
    };
    if len % stride != 0 {
        return Err("IGES 106 tuple length: IGESGeom_CopiousData.cxx:118-133".into());
    }
    Ok(len / stride)
}

fn point_at(data_type: i32, z_plane: f64, tuples: &[f64], index: usize) -> Result<GpPnt, String> {
    if data_type == 1 {
        let i = 2 * index;
        Ok(GpPnt::new(tuples[i], tuples[i + 1], z_plane))
    } else if data_type == 2 || data_type == 3 {
        let stride = if data_type == 2 { 3 } else { 6 };
        let i = stride * index;
        Ok(GpPnt::new(tuples[i], tuples[i + 1], tuples[i + 2]))
    } else {
        Err(format!(
            "IGES 106 data type {data_type}: IGESGeom_CopiousData.cxx:141-171"
        ))
    }
}
