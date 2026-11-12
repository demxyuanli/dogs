//! Face, wire, and seam helpers used by `IntUnifyFaces`.
//!
//! Source: `ShapeUpgrade_UnifySameDomain.cxx:699-918`.

use std::collections::HashMap;

use occt_core::gp::GpPnt2d;
use occt_core::precision::CONFUSION;

use crate::abs::{Orientation, ShapeType};
use crate::boptools_2d::curve_on_surface_oriented;
use crate::brep_tool::BRepTool;
use crate::brep_top_adaptor_fclass2d::BrepTopAdaptorFClass2d;
use crate::fclass2d::FaceState;
use crate::shape::{Edge, Face, TopoShape, Wire};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{
    common_vertex, edge_vertices, edges_of_wire, explorer_shapes, is_same, shape_map_key,
    ShapeMapKey,
};

/// `InsertWiresIntoFaces` (`cxx:699-733`).
pub(super) fn insert_wires_into_faces(
    wires: &[Wire],
    faces: &mut [Face],
    reference_face: &Face,
) -> Result<(), String> {
    let builder = crate::builder::TopoBuilder::new();
    for wire in wires {
        let edge = edges_of_wire(wire)
            .into_iter()
            .next()
            .ok_or("InsertWiresIntoFaces: empty wire")?;
        let (pcurve, first, last) = curve_on_surface_oriented(&edge, reference_face, false)
            .ok_or("InsertWiresIntoFaces: edge has no pcurve")?;
        let point = pcurve.d0((first + last) * 0.5);
        let face = faces
            .iter_mut()
            .find(|face| {
                BrepTopAdaptorFClass2d::new(face, CONFUSION).perform(point, false) == FaceState::In
            })
            .ok_or("ShapeUpgrade_UnifySameDomain: wire remains unclassified")?;
        face.0.set_orientation(Orientation::Forward);
        builder.add_wire(face, wire);
    }
    Ok(())
}

/// `FindCommonFace` (`cxx:735-776`).
pub(super) fn find_common_face(
    first: &Edge,
    second: &Edge,
    vertex_faces: &HashMap<ShapeMapKey, Vec<Face>>,
) -> Option<(Face, Orientation, Orientation)> {
    let vertex = common_vertex(first, second)?;
    for candidate in vertex_faces.get(&shape_map_key(&vertex))? {
        let mut first_orientation = None;
        let mut second_orientation = None;
        for edge_shape in explorer_shapes(&candidate.0, ShapeType::Edge) {
            if is_same(&edge_shape, &first.0) {
                first_orientation = Some(edge_shape.orientation());
            }
            if is_same(&edge_shape, &second.0) {
                second_orientation = Some(edge_shape.orientation());
            }
            if let (Some(first_orientation), Some(second_orientation)) =
                (first_orientation, second_orientation)
            {
                let mut face = candidate.clone();
                face.0.set_orientation(Orientation::Forward);
                return Some((face, first_orientation, second_orientation));
            }
        }
    }
    None
}

pub(super) struct ClosestEdgeEnds {
    pub face: Face,
    pub squared_distance: f64,
    pub first_orientation: Orientation,
    pub second_orientation: Orientation,
    pub first_index: usize,
    pub second_index: usize,
    pub first_points: [GpPnt2d; 2],
    pub second_points: [GpPnt2d; 2],
}

/// `FindClosestPoints` (`cxx:778-824`).
pub(super) fn find_closest_points(
    first: &Edge,
    second: &Edge,
    vertex_faces: &HashMap<ShapeMapKey, Vec<Face>>,
) -> Option<ClosestEdgeEnds> {
    let (face, first_orientation, second_orientation) =
        find_common_face(first, second, vertex_faces)?;
    let (first_curve, first_parameter, last_parameter) =
        curve_on_surface_oriented(first, &face, false)?;
    let (second_curve, second_first, second_last) =
        curve_on_surface_oriented(second, &face, false)?;
    let first_points = [
        first_curve.d0(first_parameter),
        first_curve.d0(last_parameter),
    ];
    let second_points = [second_curve.d0(second_first), second_curve.d0(second_last)];
    let mut squared_distance = f64::MAX;
    let mut first_index = 0;
    let mut second_index = 0;
    for (i, first_point) in first_points.iter().enumerate() {
        for (j, second_point) in second_points.iter().enumerate() {
            let distance = first_point.square_distance(second_point);
            if distance < squared_distance {
                squared_distance = distance;
                first_index = i;
                second_index = j;
            }
        }
    }
    Some(ClosestEdgeEnds {
        face,
        squared_distance,
        first_orientation,
        second_orientation,
        first_index,
        second_index,
        first_points,
        second_points,
    })
}

fn oriented_vertices(edge: &Edge) -> (Option<TopoShape>, Option<TopoShape>) {
    let (first, last) = edge_vertices(edge);
    if edge.0.orientation() == Orientation::Reversed {
        (last.map(|vertex| vertex.0), first.map(|vertex| vertex.0))
    } else {
        (first.map(|vertex| vertex.0), last.map(|vertex| vertex.0))
    }
}

/// `ReconstructMissedSeam` (`cxx:829-918`).
pub(super) fn reconstruct_missed_seam(
    removed_edges: &[Edge],
    reference_face: &Face,
    current_edge: &Edge,
    current_vertex: &TopoShape,
    current_point: GpPnt2d,
    u_period: f64,
    v_period: f64,
) -> Option<(Edge, GpPnt2d)> {
    let registry = GeometryRegistry::global();
    let face_key = GeometryRegistry::shape_key(&reference_face.0);
    for source in removed_edges {
        if is_same(&source.0, &current_edge.0) {
            continue;
        }
        let mut edge = source.clone();
        let (curve, first, last) = curve_on_surface_oriented(&edge, reference_face, false)?;
        let (first_vertex, last_vertex) = oriented_vertices(&edge);
        if !first_vertex
            .as_ref()
            .is_some_and(|vertex| is_same(vertex, current_vertex))
            && !last_vertex
                .as_ref()
                .is_some_and(|vertex| is_same(vertex, current_vertex))
        {
            continue;
        }
        if !BRepTool::is_closed_edge_face(&edge, reference_face) {
            continue;
        }

        let parameter = if edge.0.orientation() == Orientation::Forward {
            first
        } else {
            last
        };
        let mut point = curve.d0(parameter);
        let mut u_difference = (point.x() - current_point.x()).abs();
        let mut v_difference = (point.y() - current_point.y()).abs();
        if (u_period != 0.0 && u_difference > u_period / 2.0)
            || (v_period != 0.0 && v_difference > v_period / 2.0)
        {
            if last_vertex
                .as_ref()
                .is_some_and(|vertex| is_same(vertex, current_vertex))
                || (u_period != 0.0 && v_period != 0.0)
            {
                edge.0.reverse();
            } else {
                let curves = registry.edge_pcurves(&edge.0, face_key);
                if curves.len() < 2 {
                    continue;
                }
                registry.set_edge_pcurves(
                    &edge.0,
                    face_key,
                    vec![curves[1].clone(), curves[0].clone()],
                );
            }
            let (curve, first, last) = curve_on_surface_oriented(&edge, reference_face, false)?;
            point = curve.d0(if edge.0.orientation() == Orientation::Forward {
                first
            } else {
                last
            });
            u_difference = (point.x() - current_point.x()).abs();
            v_difference = (point.y() - current_point.y()).abs();
        }
        if (u_period == 0.0 || u_difference < u_period / 2.0)
            && (v_period == 0.0 || v_difference < v_period / 2.0)
        {
            let (first_vertex, _) = oriented_vertices(&edge);
            if first_vertex
                .as_ref()
                .is_some_and(|vertex| is_same(vertex, current_vertex))
            {
                let (curve, first, last) = curve_on_surface_oriented(&edge, reference_face, false)?;
                let next = curve.d0(if edge.0.orientation() == Orientation::Forward {
                    last
                } else {
                    first
                });
                return Some((edge, next));
            }
        }
    }
    None
}
