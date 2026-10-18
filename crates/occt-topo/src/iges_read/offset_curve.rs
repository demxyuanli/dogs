//! Entity 130. `IGESToBRep_TopoCurve::TransferOffsetCurve`
//! (`IGESToBRep_TopoCurve.cxx:583-751`) and `Transfer2dOffsetCurve`
//! (`cxx:755-855`). `IGESGeom_ToolOffsetCurve::ReadOwnParams`
//! (`IGESGeom_ToolOffsetCurve.cxx:47-209`).
//!
//! The offset distance is not multiplied by the unit factor: the basis curve
//! is already scaled by its own transfer. `TransformedNormalVector`
//! (`IGESGeom_OffsetCurve.cxx:107-118`) drops the translation, then the
//! direction is reversed (`cxx:615`). A directory 124 is the entity's own
//! matrix (`entity_location`), not `CompoundLocation` of the parents.

use std::sync::Arc;

use occt_core::gp::{GpDir, GpPnt, GpTrsf2d, GpXyz};
use occt_core::precision::{Precision, CONFUSION, PCONFUSION, REAL_SMALL};
use occt_core::toploc::TopLocLocation;
use occt_geom::curve::Curve;
use occt_geom::GeomOffsetCurve;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::offset::Geom2dOffsetCurve;

use crate::abs::ShapeType;
use crate::builder::TopoBuilder;
use crate::meshing::edge_discret::CurveOnSurface;
use crate::shape::{Edge, TopoShape};
use crate::tgeometry::GeometryRegistry;

use super::iges_boundary::connect_next_wire;
use super::iges_data_iges_model::EntityData;
use super::iges_data_param_reader::ParamReader;
use super::iges_to_brep_basic_curve::{
    convert_location, entity_location, eps_geom, transformed_dir,
};
use super::iges_to_brep_curve_and_surface::{is_basic_curve, is_topo_curve, Session};
use super::iges_to_brep_topo_curve::{
    boundary_max_tol, edges_from_2d_basic, edges_of, transfer_topo_curve,
};

struct OffsetRec {
    base: usize,
    offset_type: i32,
    distance1: f64,
    normal: [f64; 3],
    start_param: f64,
    end_param: f64,
}

pub(in crate::iges_read) fn parse_offset_curve(pr: &mut ParamReader<'_>) -> Result<EntityData, String> {
    // `IGESGeom_ToolOffsetCurve.cxx:66-209`. The second normal (`cxx:194`)
    // overwrites the first (`cxx:174`). A file that ends after the parameters
    // keeps the first normal.
    let base = pr
        .read_entity(false)?
        .ok_or_else(|| "IGES 130 null basis curve".to_string())?;
    let offset_type = pr.read_integer()?;
    let function = pr.read_entity(true)?;
    let function_coord = pr.read_integer()?;
    let tapered_type = pr.read_integer()?;
    let distance1 = pr.read_real()?;
    let arc_length1 = pr.read_real()?;
    let distance2 = pr.read_real()?;
    let arc_length2 = pr.read_real()?;
    let mut normal = pr.read_xyz()?;
    let start_param = pr.read_real()?;
    let end_param = pr.read_real()?;
    if pr.remaining() >= 3 {
        normal = pr.read_xyz()?;
    }
    Ok(EntityData::OffsetCurve {
        base,
        offset_type,
        function,
        function_coord,
        tapered_type,
        distance1,
        arc_length1,
        distance2,
        arc_length2,
        normal,
        start_param,
        end_param,
    })
}

fn offset_rec(session: &Session<'_>, index: usize) -> Result<OffsetRec, String> {
    match &session.model.entity(index)?.data {
        EntityData::OffsetCurve {
            base,
            offset_type,
            distance1,
            normal,
            start_param,
            end_param,
            ..
        } => Ok(OffsetRec {
            base: *base,
            offset_type: *offset_type,
            distance1: *distance1,
            normal: *normal,
            start_param: *start_param,
            end_param: *end_param,
        }),
        _ => Err("IGES 130 has no offset parameters".into()),
    }
}

/// `ElCLib::AdjustPeriodic` (`ElCLib.cxx:115`), used by `BRepLib_MakeEdge`
/// when the curve is periodic (`BRepLib_MakeEdge.cxx:628-632`).
fn adjust_periodic(u_first: f64, u_last: f64, preci: f64, u1: &mut f64, u2: &mut f64) {
    if Precision::is_infinite(u_first) || Precision::is_infinite(u_last) {
        *u1 = u_first;
        *u2 = u_last;
        return;
    }
    let period = u_last - u_first;
    if period < f64::EPSILON {
        *u1 = u_first;
        *u2 = u_last;
        return;
    }
    *u1 -= ((*u1 - u_first) / period).floor() * period;
    if u_last - *u1 < preci {
        *u1 -= period;
    }
    *u2 -= ((*u2 - *u1) / period).floor() * period;
    if *u2 - *u1 < preci {
        *u2 += period;
    }
}

/// `BRepLib_MakeEdge::Init` parameter check (`BRepLib_MakeEdge.cxx:618-665`).
/// `None` is `!IsDone()`, which the caller reports as `IGES_1005`.
fn ordered_range(curve: &dyn Curve, mut p1: f64, mut p2: f64) -> Option<(f64, f64)> {
    let cf = curve.first_parameter();
    let cl = curve.last_parameter();
    if curve.is_periodic() {
        adjust_periodic(cf, cl, PCONFUSION, &mut p1, &mut p2);
        return Some((p1, p2));
    }
    if p1 > p2 {
        std::mem::swap(&mut p1, &mut p2);
    }
    if (cf - p1) > PCONFUSION || (p2 - cl) > PCONFUSION {
        return None;
    }
    if (p2 - p1) <= REAL_SMALL {
        return None;
    }
    Some((p1, p2))
}

fn ordered_range_2d(curve: &dyn Curve2d, mut p1: f64, mut p2: f64) -> Option<(f64, f64)> {
    let cf = curve.first_parameter();
    let cl = curve.last_parameter();
    if curve.is_periodic() {
        adjust_periodic(cf, cl, PCONFUSION, &mut p1, &mut p2);
        return Some((p1, p2));
    }
    if p1 > p2 {
        std::mem::swap(&mut p1, &mut p2);
    }
    if (cf - p1) > PCONFUSION || (p2 - cl) > PCONFUSION {
        return None;
    }
    if (p2 - p1) <= REAL_SMALL {
        return None;
    }
    Some((p1, p2))
}

fn edge_on_curve(curve: Arc<dyn Curve>, first: f64, last: f64) -> Edge {
    // `BRepLib_MakeEdge` vertices use `BRepLib::Precision()` (`BRepLib.cxx:84`),
    // which is `Precision::Confusion()`. A closed curve shares one vertex
    // (`BRepLib_MakeEdge.cxx:691-696`).
    let b = TopoBuilder::new();
    let mut edge = b.make_edge(curve.clone(), first, last);
    let p1_inf = Precision::is_infinite(first);
    let p2_inf = Precision::is_infinite(last);
    if !p1_inf && !p2_inf {
        let p1 = curve.d0(first);
        let p2 = curve.d0(last);
        if p1.distance(&p2) <= CONFUSION {
            let v = b.make_vertex(p1, CONFUSION);
            b.add_edge_vertices(&mut edge, &v, &v);
        } else {
            let v1 = b.make_vertex(p1, CONFUSION);
            let v2 = b.make_vertex(p2, CONFUSION);
            b.add_edge_vertices(&mut edge, &v1, &v2);
        }
    }
    edge
}

fn attach_pcurve(edge: &Edge, face: &TopoShape, pc: Arc<dyn Curve2d>, first: f64, last: f64) {
    let reg = GeometryRegistry::global();
    let key = GeometryRegistry::shape_key(face);
    reg.set_edge_pcurve(&edge.0, key, pc);
    reg.set_pcurve_range(&edge.0, key, first, last);
    let (f3, l3) = reg.edge_parameters(&edge.0);
    if (f3 - first).abs() > PCONFUSION || (l3 - last).abs() > PCONFUSION {
        reg.set_same_parameter(&edge.0, false);
        reg.set_same_range(&edge.0, false);
    }
}

fn normal_dir(session: &Session<'_>, index: usize, normal: [f64; 3]) -> Result<GpDir, String> {
    // `cxx:606-615` then `TransformedNormalVector` (`IGESGeom_OffsetCurve.cxx:107-118`).
    let xyz = GpXyz::new(normal[0], normal[1], normal[2]);
    let mut dir = if let Some(g) = entity_location(session.model, index)? {
        transformed_dir(&g, xyz)?
    } else {
        GpDir::from_xyz(&xyz).map_err(|_| "IGES 130 normal is null".to_string())?
    };
    dir.reverse();
    Ok(dir)
}

fn move_own_transform(session: &Session<'_>, index: usize, shape: &mut TopoShape) -> Result<(), String> {
    // `cxx:730-748`. A failed conversion is warning `IGES_1035` and is not stored.
    if let Some(g) = entity_location(session.model, index)? {
        if let Some(tr) = convert_location(1.0e-4, &g, session.model.unit_factor) {
            shape.move_location(&TopLocLocation::composed(&tr, TopLocLocation::identity()));
        }
    }
    Ok(())
}

pub(in crate::iges_read) fn transfer_offset_curve(
    session: &mut Session<'_>,
    index: usize,
) -> Result<TopoShape, String> {
    // `TransferOffsetCurve` (`IGESToBRep_TopoCurve.cxx:583-751`).
    let start = offset_rec(session, index)?;
    if start.offset_type != 1 {
        return Err("IGES_1100: IGESToBRep_TopoCurve.cxx:596-601".into());
    }
    let nrm = normal_dir(session, index, start.normal)?;
    let mut offset = start.distance1;
    let mut base = start.base;
    while session.model.entity(base)?.type_number == 130 {
        let nested = offset_rec(session, base)?;
        if nested.offset_type != 1 {
            return Err("IGES_1100: IGESToBRep_TopoCurve.cxx:622-626".into());
        }
        offset += nested.distance1;
        base = nested.base;
    }
    if !is_topo_curve(session.model.entity(base)?.type_number) {
        return Err("XSTEP_110: IGESToBRep_TopoCurve.cxx:632-636".into());
    }
    let shape = transfer_topo_curve(session, base)?;
    let mut result = if shape.is_edge() {
        offset_edge(&shape, offset, nrm, start.start_param, start.end_param)?
    } else if shape.is_wire() {
        offset_wire(&shape, offset, nrm, start.start_param, start.end_param)?
    } else {
        return Err("IGES_1156: IGESToBRep_TopoCurve.cxx:645-653".into());
    };
    move_own_transform(session, index, &mut result)?;
    Ok(result)
}

fn offset_edge(
    shape: &TopoShape,
    offset: f64,
    nrm: GpDir,
    start_param: f64,
    end_param: f64,
) -> Result<TopoShape, String> {
    // `cxx:656-673`.
    let reg = GeometryRegistry::global();
    let curve = reg
        .edge_curve(shape)
        .ok_or_else(|| "IGES_1005: IGESToBRep_TopoCurve.cxx:663-668".to_string())?;
    let off: Arc<dyn Curve> = Arc::new(GeomOffsetCurve::new(curve, offset, nrm));
    let Some((first, last)) = ordered_range(off.as_ref(), start_param, end_param) else {
        return Err("IGES_1005: IGESToBRep_TopoCurve.cxx:663-668".into());
    };
    let mut edge = edge_on_curve(off, first, last);
    edge.0.move_location(shape.location());
    Ok(edge.0)
}

fn offset_wire(
    shape: &TopoShape,
    offset: f64,
    nrm: GpDir,
    sta_par: f64,
    end_par: f64,
) -> Result<TopoShape, String> {
    // `cxx:676-727`. Edges that end at or before `staPar` are skipped without
    // advancing `length`. `FixConnected()` uses the default precision, which
    // is `MaxTolerance()` (`ShapeFix_Wire.cxx:285`, `ShapeFix_Root.cxx:26`).
    let reg = GeometryRegistry::global();
    let mut built = Vec::new();
    let mut begin = true;
    let mut length = 0.0;
    for edge in edges_of(shape)? {
        let (first, last) = reg.edge_parameters(&edge.0);
        if length + last - first <= sta_par {
            continue;
        }
        if length >= end_par {
            if begin {
                return Err("IGES_1105: IGESToBRep_TopoCurve.cxx:696-701".into());
            }
            break;
        }
        let curve = reg.edge_curve(&edge.0).ok_or_else(|| {
            "IGES_1005: IGESToBRep_TopoCurve.cxx:709-714".to_string()
        })?;
        let off: Arc<dyn Curve> = Arc::new(GeomOffsetCurve::new(curve, offset, nrm));
        let Some((p1, p2)) = ordered_range(off.as_ref(), sta_par - length, end_par - length) else {
            return Err("IGES_1005: IGESToBRep_TopoCurve.cxx:709-714".into());
        };
        let mut another = edge_on_curve(off, p1, p2);
        another.0.move_location(edge.0.location());
        begin = false;
        length += last - first;
        built.push(another);
    }
    let mut wire = TopoBuilder::new().make_wire(&built);
    crate::shhealing::fix_connected_closed(&mut wire, CONFUSION, CONFUSION, true);
    Ok(wire.0)
}

pub(in crate::iges_read) fn transfer_2d_offset_curve(
    session: &mut Session<'_>,
    index: usize,
    face: &TopoShape,
    trans: &GpTrsf2d,
    u_fact: f64,
) -> Result<TopoShape, String> {
    // `Transfer2dOffsetCurve` (`cxx:755-855`). `OffsetType != 1` sends
    // `IGES_1100` and does not return (`cxx:771-775`). Nested offsets are not
    // walked. The distance is `FirstOffsetDistance * uFact` (`cxx:804`).
    let start = offset_rec(session, index)?;
    if !is_topo_curve(session.model.entity(start.base)?.type_number) {
        return Err("IGES 130 2d basis is not a topo curve: IGESToBRep_TopoCurve.cxx:783".into());
    }
    let Some(shape) = transfer_2d_shape(session, start.base, face, trans, u_fact)? else {
        return Err("IGES_1156: IGESToBRep_TopoCurve.cxx:788-795".into());
    };
    let distance = start.distance1 * u_fact;
    if shape.is_edge() {
        offset_edge_2d(&shape, face, distance, start.start_param, start.end_param)
    } else if shape.is_wire() {
        offset_wire_2d(&shape, face, distance, start.start_param, start.end_param)
    } else {
        Err("IGES_1156: IGESToBRep_TopoCurve.cxx:788-795".into())
    }
}

fn transfer_2d_shape(
    session: &mut Session<'_>,
    index: usize,
    face: &TopoShape,
    trans: &GpTrsf2d,
    u_fact: f64,
) -> Result<Option<TopoShape>, String> {
    // `Transfer2dTopoCurve` (`cxx:910-946`). A type with no arm returns null.
    // `SetModeTransfer(false)` (`cxx:786`) selects `TransformedValue` for a point.
    let ty = session.model.entity(index)?.type_number;
    if is_basic_curve(ty) {
        let Ok(mut edges) = edges_from_2d_basic(session, index, face, trans, u_fact) else {
            return Ok(None);
        };
        if edges.is_empty() {
            return Ok(None);
        }
        // `cxx:1294-1308`: a single edge is that edge, otherwise the wire.
        if edges.len() == 1 {
            return Ok(Some(edges.pop().expect("one edge").0));
        }
        let mut wire = TopoBuilder::new().make_wire(&edges);
        let _ = crate::shhealing::fix_connected_closed(&mut wire, CONFUSION, CONFUSION, true);
        return Ok(Some(wire.0));
    }
    match ty {
        130 => Ok(Some(transfer_2d_offset_curve(
            session, index, face, trans, u_fact,
        )?)),
        102 => transfer_2d_composite(session, index, face, trans, u_fact),
        116 => Ok(Some(transfer_2d_point(session, index)?)),
        _ => Ok(None),
    }
}

fn transfer_2d_point(session: &Session<'_>, index: usize) -> Result<TopoShape, String> {
    // `Transfer2dPoint` (`cxx:176-201`) with `GetModeTransfer` false.
    let xyz = match &session.model.entity(index)?.data {
        EntityData::Point { xyz } => *xyz,
        _ => return Err("IGES 116 has no point parameters".into()),
    };
    let point = if entity_location(session.model, index)?.is_some() {
        let g = entity_location(session.model, index)?.expect("location");
        let moved = super::iges_to_brep_basic_curve::apply_gtrsf(
            &g,
            GpPnt::new(xyz[0], xyz[1], xyz[2]),
        );
        GpPnt::new(moved.x(), moved.y(), 0.0)
    } else {
        GpPnt::new(xyz[0], xyz[1], 0.0)
    };
    Ok(TopoBuilder::new().make_vertex(point, CONFUSION).into())
}

fn transfer_2d_composite(
    session: &mut Session<'_>,
    index: usize,
    face: &TopoShape,
    trans: &GpTrsf2d,
    u_fact: f64,
) -> Result<Option<TopoShape>, String> {
    // `TransferCompositeCurveGeneral` with `is2d` (`cxx:210-347`).
    let curves = match &session.model.entity(index)?.data {
        EntityData::Composite { curves } => curves.clone(),
        _ => return Err("IGES 102 has no curve list".into()),
    };
    let precision = eps_geom(session.model.resolution) * session.model.unit_factor;
    let maxtol = boundary_max_tol(session);
    let mut wire_edges = Vec::new();
    let mut prev: Option<usize> = None;
    for (n, child) in curves.into_iter().enumerate() {
        if prev == Some(child) {
            prev = Some(child);
            continue;
        }
        prev = Some(child);
        let ty = session.model.entity(child)?.type_number;
        if ty == 141 || ty == 142 {
            return Err(format!(
                "IGES_1040: IGESToBRep_TopoCurve.cxx:260-268 component {}",
                n + 1
            ));
        }
        let Some(shape) = transfer_2d_shape(session, child, face, trans, u_fact)? else {
            continue;
        };
        if shape.shape_type() == ShapeType::Vertex {
            continue;
        }
        let mut next = edges_of(&shape)?;
        if !connect_next_wire(&mut wire_edges, &mut next, precision, maxtol) {
            return Err(format!(
                "IGES_1060: IGESToBRep_TopoCurve.cxx:329-333 component {}",
                n + 1
            ));
        }
    }
    let mut wire = TopoBuilder::new().make_wire(&wire_edges);
    // `ClosedWireMode` is false (`cxx:341`), then `FixConnected(maxtol)` and
    // `FixConnected(1, precision)` (`cxx:342-343`).
    let nb = wire_edges.len();
    if nb >= 2 {
        for i in (2..=nb).rev() {
            crate::shhealing::fix_connected_at(&mut wire, CONFUSION, maxtol, i);
        }
    }
    if nb >= 1 {
        crate::shhealing::fix_connected_at(&mut wire, CONFUSION, precision, 1);
    }
    Ok(Some(wire.0))
}

fn offset_edge_2d(
    edge: &TopoShape,
    face: &TopoShape,
    distance: f64,
    start_param: f64,
    end_param: f64,
) -> Result<TopoShape, String> {
    // `cxx:799-817` through `ShapeBuild_Edge::MakeEdge(pcurve, face, p1, p2)`
    // (`ShapeBuild_Edge.cxx:828-836`). The basis edge location is not applied.
    let built = make_offset_pcurve_edge(edge, face, distance, start_param, end_param)?;
    Ok(built.0)
}

fn offset_wire_2d(
    shape: &TopoShape,
    face: &TopoShape,
    distance: f64,
    start_param: f64,
    end_param: f64,
) -> Result<TopoShape, String> {
    // `cxx:819-850`. Every edge uses the offset entity's parameter window.
    let mut built = Vec::new();
    for edge in edges_of(shape)? {
        built.push(make_offset_pcurve_edge(
            &edge.0,
            face,
            distance,
            start_param,
            end_param,
        )?);
    }
    let mut wire = TopoBuilder::new().make_wire(&built);
    crate::shhealing::fix_connected_closed(&mut wire, CONFUSION, CONFUSION, true);
    Ok(wire.0)
}

fn make_offset_pcurve_edge(
    edge: &TopoShape,
    face: &TopoShape,
    distance: f64,
    start_param: f64,
    end_param: f64,
) -> Result<Edge, String> {
    let reg = GeometryRegistry::global();
    let key = GeometryRegistry::shape_key(face);
    let pc = reg.edge_pcurve(edge, key).ok_or_else(|| {
        "IGES_1005: IGESToBRep_TopoCurve.cxx:811-815".to_string()
    })?;
    let off: Arc<dyn Curve2d> = Arc::new(Geom2dOffsetCurve::new(pc, distance));
    let Some((first, last)) = ordered_range_2d(off.as_ref(), start_param, end_param) else {
        return Err("IGES_1005: IGESToBRep_TopoCurve.cxx:811-815".into());
    };
    let surf = reg
        .face_surface(face)
        .ok_or_else(|| "IGES face has no surface".to_string())?;
    let curve3d: Arc<dyn Curve> = Arc::new(CurveOnSurface::new(off.clone(), surf, first, last));
    let mut built = edge_on_curve(curve3d, first, last);
    attach_pcurve(&built, face, off, first, last);
    if !face.location().is_identity() {
        built.0.move_location(face.location());
    }
    Ok(built)
}
