//! Edge-chain construction and splitting.
//!
//! Source: `ShapeUpgrade_UnifySameDomain.cxx:2504-2877`.

use std::collections::{HashMap, HashSet};

use occt_core::gp::{GpDir, GpPnt, GpVec};
use occt_core::precision::{CONFUSION, SQUARE_CONFUSION};
use occt_geom::geom_abs_curve_type::{curve_kind, GeomAbsCurveType};
use occt_geom::Curve;

use crate::abs::Orientation;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::topo_tools::{edge_vertices, is_equal, is_same, shape_map_key, ShapeMapKey};

pub(super) type AncestorMap = HashMap<ShapeMapKey, Vec<TopoShape>>;

#[derive(Clone)]
pub(super) struct SubSequenceOfEdges {
    pub(super) edges: Vec<Edge>,
    pub(super) union_edge: Option<Edge>,
}

impl SubSequenceOfEdges {
    fn new(edge: Edge) -> Self {
        Self {
            edges: vec![edge],
            union_edge: None,
        }
    }
}

pub(super) fn oriented_vertices(edge: &Edge) -> (Option<Vertex>, Option<Vertex>) {
    let (first, last) = edge_vertices(edge);
    if edge.0.orientation() == Orientation::Reversed {
        (last, first)
    } else {
        (first, last)
    }
}

fn find_common_face(
    edge1: &Edge,
    edge2: &Edge,
    vertex_faces: &AncestorMap,
) -> Option<(Face, Orientation, Orientation)> {
    let common = crate::topo_tools::common_vertex(edge1, edge2)?;
    for face_shape in vertex_faces.get(&shape_map_key(&common))? {
        let mut face = Face::wrap(face_shape.clone())?;
        face.0.set_orientation(Orientation::Forward);
        let mut first_orientation = None;
        let mut second_orientation = None;
        for occurrence in crate::topo_tools::explorer_shapes(&face.0, crate::abs::ShapeType::Edge) {
            if is_same(&occurrence, &edge1.0) {
                first_orientation = Some(occurrence.orientation());
            }
            if is_same(&occurrence, &edge2.0) {
                second_orientation = Some(occurrence.orientation());
            }
            if let (Some(first), Some(second)) = (first_orientation, second_orientation) {
                return Some((face, first, second));
            }
        }
    }
    None
}

pub(super) struct ClosestPoints {
    pub(super) face: Face,
    pub(super) first_orientation: Orientation,
    pub(super) second_orientation: Orientation,
    pub(super) first_points: [occt_core::gp::GpPnt2d; 2],
    pub(super) second_points: [occt_core::gp::GpPnt2d; 2],
    pub(super) first_index: usize,
    pub(super) second_index: usize,
    pub(super) squared_distance: f64,
}

/// `FindClosestPoints` (`cxx:778-824`).
pub(super) fn find_closest_points(
    edge1: &Edge,
    edge2: &Edge,
    vertex_faces: &AncestorMap,
) -> Option<ClosestPoints> {
    let (face, first_orientation, second_orientation) =
        find_common_face(edge1, edge2, vertex_faces)?;
    let (first_curve, first, last) =
        crate::boptools_2d::curve_on_surface_oriented(edge1, &face, false)?;
    let (second_curve, second_first, second_last) =
        crate::boptools_2d::curve_on_surface_oriented(edge2, &face, false)?;
    let first_points = [first_curve.d0(first), first_curve.d0(last)];
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
    Some(ClosestPoints {
        face,
        first_orientation,
        second_orientation,
        first_points,
        second_points,
        first_index,
        second_index,
        squared_distance,
    })
}

fn line_edge_points(edge: &Edge) -> Option<(GpPnt, GpVec)> {
    let curve = BRepTool::edge_curve_adaptor_world(edge)?;
    if curve.gp_line().is_none() {
        return None;
    }
    let (mut first, mut last) = BRepTool::edge_parameters(edge);
    if edge.0.orientation() == Orientation::Reversed {
        std::mem::swap(&mut first, &mut last);
    }
    let first_point = curve.d0(first);
    let last_point = curve.d0(last);
    let mut direction = GpVec::from_pnts(&first_point, &last_point);
    if direction.normalize() <= 0.0 {
        return None;
    }
    Some((first_point, direction))
}

/// `IsLinear` (`cxx:121-140`).
pub(super) fn linear_direction(curve: &dyn Curve, first: f64, last: f64) -> Option<GpDir> {
    if let Some(line) = curve.gp_line() {
        return Some(line.direction());
    }
    let kind = curve_kind(curve);
    let has_two_poles = curve.bspline_poles().is_some_and(|poles| poles.len() == 2)
        || curve.bezier_poles().is_some_and(|poles| poles.len() == 2);
    if matches!(
        kind,
        GeomAbsCurveType::BezierCurve | GeomAbsCurveType::BSplineCurve
    ) && has_two_poles
    {
        return GpDir::from_vec(&GpVec::from_pnts(&curve.d0(first), &curve.d0(last))).ok();
    }
    None
}

/// `IsMergingPossible` (`cxx:2504-2638`).
fn is_merging_possible(
    edge1: &Edge,
    edge2: &Edge,
    angular_tolerance: f64,
    linear_tolerance: f64,
    avoid_vertices: &HashSet<ShapeMapKey>,
    line_reference: Option<(GpPnt, GpVec)>,
    vertex_faces: &AncestorMap,
) -> bool {
    let first_degenerated = BRepTool::is_degenerated(edge1);
    let second_degenerated = BRepTool::is_degenerated(edge2);
    if first_degenerated && second_degenerated {
        return find_closest_points(edge1, edge2, vertex_faces)
            .is_some_and(|points| points.squared_distance <= SQUARE_CONFUSION);
    }
    if first_degenerated || second_degenerated {
        return false;
    }

    let Some(common_vertex) = oriented_vertices(edge1).1 else {
        return false;
    };
    if avoid_vertices.contains(&shape_map_key(&common_vertex.0)) {
        return false;
    }

    let Some(first_curve) = BRepTool::edge_curve_adaptor_world(edge1) else {
        return false;
    };
    let Some(second_curve) = BRepTool::edge_curve_adaptor_world(edge2) else {
        return false;
    };
    let (first_start, first_end) = BRepTool::edge_parameters(edge1);
    let (second_start, second_end) = BRepTool::edge_parameters(edge2);
    let first_line = linear_direction(first_curve.as_ref(), first_start, first_end);
    let second_line = linear_direction(second_curve.as_ref(), second_start, second_end);
    let first_circle = first_curve.gp_circ();
    let second_circle = second_curve.gp_circ();
    if let (Some(first), Some(second)) = (&first_circle, &second_circle) {
        if first.location().distance(&second.location()) > CONFUSION {
            return false;
        }
    }
    let first_type = curve_kind(first_curve.as_ref());
    let second_type = curve_kind(second_curve.as_ref());
    let first_spline = matches!(
        first_type,
        GeomAbsCurveType::BezierCurve | GeomAbsCurveType::BSplineCurve
    );
    let second_spline = matches!(
        second_type,
        GeomAbsCurveType::BezierCurve | GeomAbsCurveType::BSplineCurve
    );
    if !(first_line.is_some() && second_line.is_some())
        && !(first_spline && second_spline)
        && first_type != second_type
    {
        return false;
    }

    let first_parameter = if edge1.0.orientation() == Orientation::Forward {
        first_end
    } else {
        first_start
    };
    let second_parameter = if edge2.0.orientation() == Orientation::Forward {
        second_start
    } else {
        second_end
    };
    let (_, mut first_derivative) = first_curve.d1(first_parameter);
    let (_, mut second_derivative) = second_curve.d1(second_parameter);
    if edge1.0.orientation() == Orientation::Reversed {
        first_derivative.reverse();
    }
    if edge2.0.orientation() == Orientation::Reversed {
        second_derivative.reverse();
    }
    if first_derivative.angle(&second_derivative) > angular_tolerance {
        return false;
    }

    if let (Some((first_point, direction)), Some(_)) = (line_reference, second_line) {
        let last = if edge2.0.orientation() == Orientation::Forward {
            second_end
        } else {
            second_start
        };
        let accumulated = GpVec::from_pnts(&first_point, &second_curve.d0(last));
        if direction.cross_square_magnitude(&accumulated) > linear_tolerance * linear_tolerance {
            return false;
        }
        if direction.angle(&accumulated) > angular_tolerance
            || second_derivative.angle(&accumulated) > angular_tolerance
        {
            return false;
        }
    }
    true
}

/// `generateSubSeq` (`cxx:2685-2751`).
pub(super) fn generate_subsequences(
    input: &[Edge],
    is_closed: bool,
    angular_tolerance: f64,
    linear_tolerance: f64,
    avoid_vertices: &HashSet<ShapeMapKey>,
    vertex_faces: &AncestorMap,
) -> Vec<SubSequenceOfEdges> {
    let Some(first) = input.first() else {
        return Vec::new();
    };
    let mut result = vec![SubSequenceOfEdges::new(first.clone())];
    let mut line_reference = line_edge_points(first);
    for pair in input.windows(2) {
        if is_merging_possible(
            &pair[0],
            &pair[1],
            angular_tolerance,
            linear_tolerance,
            avoid_vertices,
            line_reference,
            vertex_faces,
        ) {
            result.last_mut().unwrap().edges.push(pair[1].clone());
        } else {
            result.push(SubSequenceOfEdges::new(pair[1].clone()));
            line_reference = line_edge_points(&pair[1]);
        }
    }
    if is_closed
        && result.len() > 1
        && is_merging_possible(
            input.last().unwrap(),
            input.first().unwrap(),
            angular_tolerance,
            linear_tolerance,
            avoid_vertices,
            None,
            vertex_faces,
        )
    {
        let first_edges = result.remove(0).edges;
        result.last_mut().unwrap().edges.extend(first_edges);
    }
    result
}

fn vertex_edge_map(edges: &[Edge]) -> HashMap<ShapeMapKey, (TopoShape, Vec<Edge>)> {
    let mut result = HashMap::new();
    for edge in edges {
        let (first, last) = edge_vertices(edge);
        for vertex in [first, last].into_iter().flatten() {
            result
                .entry(shape_map_key(&vertex.0))
                .or_insert_with(|| (vertex.0.clone(), Vec::new()))
                .1
                .push(edge.clone());
        }
    }
    result
}

/// Chain-building half of `MergeEdges` (`cxx:2756-2862`).
pub(super) fn split_edge_chains(
    edges: &[Edge],
    vertex_faces: &AncestorMap,
    non_merging_vertices: &HashSet<ShapeMapKey>,
    angular_tolerance: f64,
    linear_tolerance: f64,
) -> Vec<SubSequenceOfEdges> {
    let vertex_edges = vertex_edge_map(edges);
    let mut used = HashSet::new();
    let mut result = Vec::new();
    for seed in edges {
        if !used.insert(shape_map_key(&seed.0)) {
            continue;
        }
        let mut chain = vec![seed.clone()];
        let (first, last) = edge_vertices(seed);
        let mut ends = [first, last];
        for side in 0..2 {
            loop {
                let Some(endpoint) = ends[side].clone() else {
                    break;
                };
                let Some((_, candidates)) = vertex_edges.get(&shape_map_key(&endpoint.0)) else {
                    break;
                };
                let mut selected = None;
                for candidate in candidates {
                    if used.contains(&shape_map_key(&candidate.0)) {
                        continue;
                    }
                    let candidate_ends = edge_vertices(candidate);
                    let opposite = if side == 0 {
                        candidate_ends.1.as_ref()
                    } else {
                        candidate_ends.0.as_ref()
                    };
                    let mut reversed_endpoint = endpoint.0.clone();
                    reversed_endpoint.reverse();
                    if opposite.is_some_and(|vertex| is_equal(&vertex.0, &reversed_endpoint)) {
                        selected = Some((candidate.clone(), candidate_ends));
                        break;
                    }
                }
                let Some((candidate, candidate_ends)) = selected else {
                    break;
                };
                used.insert(shape_map_key(&candidate.0));
                if side == 0 {
                    chain.insert(0, candidate);
                    ends[0] = candidate_ends.0;
                } else {
                    chain.push(candidate);
                    ends[1] = candidate_ends.1;
                }
            }
        }
        if chain.len() < 2 {
            continue;
        }
        let is_closed = matches!((&ends[0], &ends[1]), (Some(a), Some(b)) if is_same(&a.0, &b.0));
        result.extend(generate_subsequences(
            &chain,
            is_closed,
            angular_tolerance,
            linear_tolerance,
            non_merging_vertices,
            vertex_faces,
        ));
    }
    result
}

/// `CheckSharedVertices` (`cxx:2913-2948`).
pub(super) fn shared_vertices(
    edges: &[Edge],
    vertex_edges: &AncestorMap,
    keep_shapes: &HashSet<ShapeMapKey>,
) -> HashSet<ShapeMapKey> {
    let mut seen = HashSet::new();
    let mut repeated = Vec::new();
    for edge in edges {
        let (first, last) = oriented_vertices(edge);
        for vertex in [first, last].into_iter().flatten() {
            let key = shape_map_key(&vertex.0);
            if !seen.insert(key.clone()) {
                repeated.push(key);
            }
        }
    }
    repeated
        .into_iter()
        .filter(|key| {
            vertex_edges
                .get(key)
                .is_some_and(|ancestors| ancestors.len() > 2)
                || keep_shapes.contains(key)
        })
        .collect()
}
