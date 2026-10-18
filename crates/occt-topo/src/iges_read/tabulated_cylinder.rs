//! `IGESToBRep_TopoSurface::TransferTabulatedCylinder` (`cxx:884-1036`).
//!
//! One finite 3D edge builds `Geom_SurfaceOfLinearExtrusion` and
//! `BRepBuilderAPI_MakeFace` (`cxx:950-970`). Any other edge count takes
//! `BRepPrimAPI_MakePrism` (`cxx:989`) through `BRepSweep_Translation`.

use std::sync::Arc;

use occt_core::gp::{GpDir, GpPnt, GpVec};
use occt_core::precision::{CONFUSION, PCONFUSION};
use occt_core::toploc::TopLocLocation;
use occt_geom::curve::Curve;
use occt_geom::{GeomRectangularTrimmedSurface, GeomSurfaceOfLinearExtrusion, Surface};

use crate::abs::Orientation;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, TopoShape, Vertex, Wire};
use crate::topo_tools::{edge_vertices, edges_of_wire, is_same};

use super::iges_data_iges_model::EntityData;
use super::iges_to_brep_basic_curve::{convert_location, entity_location};
use super::iges_to_brep_curve_and_surface::{is_topo_curve, Session};
use super::iges_to_brep_topo_curve::transfer_topo_curve;
use super::iges_to_brep_topo_surface::extract_curve_3d;

pub(in crate::iges_read) fn transfer_tabulated_cylinder(
    session: &mut Session<'_>,
    index: usize,
) -> Result<TopoShape, String> {
    let (directrix, end) = match &session.model.entity(index)?.data {
        EntityData::Tabulated { directrix, end } => (*directrix, *end),
        _ => return Err("IGES 122 has no tabulated-cylinder parameters".into()),
    };
    if !is_topo_curve(session.model.entity(directrix)?.type_number) {
        return Err("XSTEP_153: IGESToBRep_TopoSurface.cxx:902-908".into());
    }
    let shape = transfer_topo_curve(session, directrix)?;
    let (first, _last) = find_bounds(&shape).ok_or_else(|| {
        "ShapeAnalysis::FindBounds: IGESToBRep_TopoSurface.cxx:929".to_string()
    })?;
    let pt1 = BRepTool::vertex_point_world(&first);
    let mut pt2 = GpPnt::new(end[0], end[1], end[2]);
    pt2.scale(&GpPnt::new(0.0, 0.0, 0.0), session.model.unit_factor);
    let extrusion = GpVec::from_pnts(&pt1, &pt2);
    session.u_length = pt1.distance(&pt2);
    if session.u_length < CONFUSION {
        return Err(
            "Tabulated cylinder with zero length: IGESToBRep_TopoSurface.cxx:935-939".into(),
        );
    }
    let Some(span) = extract_curve_3d(&shape) else {
        // `BRepPrimAPI_MakePrism(directrix, Vec(pt1, pt2))` (`cxx:989-995`).
        let prism = crate::brep_sweep::sweep_make_prism(&shape, &extrusion)?;
        return apply_location(session, index, prism);
    };
    // `BRep_Tool::Curve(E, First, Last)` (`BRep_Tool.cxx:205-218`) returns the
    // curve transformed by the edge location. `extract_curve_3d` keeps the
    // local curve for entity 120.
    let basis = curve_in_edge_frame(&shape, span.curve);
    let dir = GpDir::from_vec(&extrusion)
        .map_err(|_| "IGES 122 extrusion direction is null".to_string())?;
    let swept = Arc::new(GeomSurfaceOfLinearExtrusion::new(basis, dir));
    let surf: Arc<dyn Surface> = Arc::new(
        GeomRectangularTrimmedSurface::uv(swept, span.first, span.last, 0.0, extrusion.magnitude())
            .map_err(|e| format!("IGES 122: {e}"))?,
    );
    let mut face = TopoBuilder::new().make_face(surf, &[]);
    crate::tgeometry::GeometryRegistry::global().set_face_tolerance(&face.0, CONFUSION);
    let (_u0, _u1, v0, v1) = BRepTool::uv_bounds(&face);
    if v1 <= PCONFUSION && v0 < -PCONFUSION {
        session.u_length = -session.u_length;
        face.0.reverse();
    }
    apply_location(session, index, face.0)
}

/// `ShapeAnalysis::FindBounds` (`ShapeAnalysis.cxx:73-102`).
fn find_bounds(shape: &TopoShape) -> Option<(Vertex, Vertex)> {
    if shape.is_edge() {
        let edge = Edge::wrap(shape.clone())?;
        oriented_ends(&edge)
    } else if shape.is_wire() {
        let wire = Wire::wrap(shape.clone())?;
        wire_bounds(&wire)
    } else if shape.is_vertex() {
        let vertex = Vertex::wrap(shape.clone())?;
        Some((vertex.clone(), vertex))
    } else {
        None
    }
}

/// `ShapeAnalysis_Edge::FirstVertex` / `LastVertex` (`ShapeAnalysis_Edge.cxx:228-251`).
fn oriented_ends(edge: &Edge) -> Option<(Vertex, Vertex)> {
    let (stored_first, stored_last) = edge_vertices(edge);
    let (first, last) = if edge.0.orientation() == Orientation::Reversed {
        (stored_last, stored_first)
    } else {
        (stored_first, stored_last)
    };
    Some((first?, last?))
}

/// `TopExp::Vertices(Wire)` (`TopExp.cxx:255-318`).
fn wire_bounds(wire: &Wire) -> Option<(Vertex, Vertex)> {
    let mut open: Vec<(Vertex, Orientation)> = Vec::new();
    let mut last_end: Option<Vertex> = None;
    for edge in edges_of_wire(wire) {
        let (start, end) = oriented_ends(&edge)?;
        last_end = Some(end.clone());
        toggle(&mut open, start, Orientation::Forward);
        toggle(&mut open, end, Orientation::Reversed);
    }
    if open.is_empty() {
        let end = last_end?;
        return Some((end.clone(), end));
    }
    if open.len() != 2 {
        return None;
    }
    let first = open
        .iter()
        .find(|(_, ori)| *ori == Orientation::Forward)
        .map(|(v, _)| v.clone())?;
    let last = open
        .iter()
        .find(|(_, ori)| *ori == Orientation::Reversed)
        .map(|(v, _)| v.clone())?;
    Some((first, last))
}

/// `TopTools_Map` toggle: `IsSame` covers `TShape` and location
/// (`TopoDS_Shape.hxx:262-267`).
fn toggle(open: &mut Vec<(Vertex, Orientation)>, vertex: Vertex, orientation: Orientation) {
    if let Some(i) = open.iter().position(|(kept, _)| is_same(&kept.0, &vertex.0)) {
        open.remove(i);
    } else {
        open.push((vertex, orientation));
    }
}

fn curve_in_edge_frame(shape: &TopoShape, curve: Arc<dyn Curve>) -> Arc<dyn Curve> {
    let edge = if shape.is_edge() {
        Edge::wrap(shape.clone())
    } else if shape.is_wire() {
        Wire::wrap(shape.clone()).and_then(|wire| edges_of_wire(&wire).into_iter().next())
    } else {
        None
    };
    let Some(edge) = edge else {
        return curve;
    };
    if edge.0.location().is_identity() {
        return curve;
    }
    let mut moved = curve.clone_dyn();
    moved.transform(&edge.0.location().transformation());
    Arc::from(moved)
}

fn apply_location(
    session: &Session<'_>,
    index: usize,
    mut shape: TopoShape,
) -> Result<TopoShape, String> {
    if let Some(g) = entity_location(session.model, index)? {
        if let Some(tr) = convert_location(1.0e-4, &g, session.model.unit_factor) {
            if tr.form() != occt_core::gp::TrsfForm::Identity {
                shape.move_location(&TopLocLocation::composed(
                    &tr,
                    TopLocLocation::identity(),
                ));
            }
        }
    }
    Ok(shape)
}
