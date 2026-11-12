//! UV-bound and pcurve relocation helpers used by face unification.
//!
//! Source: `ShapeUpgrade_UnifySameDomain.cxx:147-697`.

use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use occt_core::gp::{GpPnt2d, GpTrsf2d, GpVec2d};
use occt_core::precision::CONFUSION;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::trimmed::Geom2dTrimmedCurve;

use crate::abs::{Orientation, ShapeType};
use crate::boptools_2d::curve_on_surface_oriented;
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::{edge_vertices, explorer_shapes, is_same, shape_map_key, ShapeMapKey};

pub(super) type EdgeFaceMap = HashMap<ShapeMapKey, Vec<TopoShape>>;
pub(super) type VertexEdgeMap = HashMap<ShapeMapKey, (TopoShape, Vec<Edge>)>;

fn coordinate(point: &GpPnt2d, index: usize) -> f64 {
    if index == 1 {
        point.x()
    } else {
        point.y()
    }
}

/// `GetFaceFromSeq` (`cxx:158-236`).
fn face_from_sequence(
    faces: &[Face],
    reference_face: &Face,
    edge_faces: &EdgeFaceMap,
    used_faces: &mut Vec<Face>,
) -> Option<(Face, usize)> {
    if used_faces.is_empty() {
        let face = faces.first()?.clone();
        used_faces.push(face.clone());
        return Some((face, 0));
    }

    let mut target = None;
    for used in used_faces.iter() {
        let mut seam_target = None;
        for edge_shape in explorer_shapes(&used.0, ShapeType::Edge) {
            let adjacent = edge_faces.get(&shape_map_key(&edge_shape))?;
            if adjacent.len() < 2 {
                continue;
            }
            let candidate = if is_same(&adjacent[0], &used.0) {
                Face::wrap(adjacent[1].clone())?
            } else {
                Face::wrap(adjacent[0].clone())?
            };
            if used_faces.iter().any(|face| is_same(&face.0, &candidate.0)) {
                continue;
            }
            let edge = Edge::wrap(edge_shape)?;
            if crate::brep_tool::BRepTool::is_closed_edge_face(&edge, reference_face) {
                seam_target = Some(candidate);
                continue;
            }
            target = Some(candidate);
            break;
        }
        if target.is_none() {
            target = seam_target;
        }
        if target.is_some() {
            break;
        }
    }

    let target = target?;
    let index = faces.iter().position(|face| is_same(&face.0, &target.0))?;
    used_faces.push(target.clone());
    Some((target, index))
}

/// `UpdateBoundaries` (`cxx:241-266`).
pub(super) fn update_boundaries(
    pcurve: &dyn Curve2d,
    first: f64,
    last: f64,
    coordinate_index: usize,
    min_coordinate: &mut f64,
    max_coordinate: &mut f64,
) {
    let delta = (last - first) / 4.0;
    for sample in 0..=4 {
        let value = coordinate(&pcurve.d0(first + sample as f64 * delta), coordinate_index);
        *min_coordinate = min_coordinate.min(value);
        *max_coordinate = max_coordinate.max(value);
    }
}

/// `FindCoordBounds` (`cxx:384-507`).
pub(super) fn find_coord_bounds(
    faces: &[Face],
    reference_face: &Face,
    edge_faces: &EdgeFaceMap,
    boundary_edges: &HashSet<ShapeMapKey>,
    coordinate_index: usize,
    period: f64,
) -> Option<(f64, f64, usize, usize)> {
    let mut intervals: Vec<(f64, f64)> = Vec::new();
    let mut simple_max = f64::NEG_INFINITY;
    let mut max_face_index = 0;
    let mut used_faces = Vec::new();

    for _ in faces {
        let (face, face_index) =
            face_from_sequence(faces, reference_face, edge_faces, &mut used_faces)?;
        let mut local_min = f64::INFINITY;
        let mut local_max = f64::NEG_INFINITY;
        for edge_shape in explorer_shapes(&face.0, ShapeType::Edge) {
            if !boundary_edges.contains(&shape_map_key(&edge_shape)) {
                continue;
            }
            let edge = Edge::wrap(edge_shape)?;
            let (pcurve, first, last) = curve_on_surface_oriented(&edge, reference_face, false)?;
            update_boundaries(
                pcurve.as_ref(),
                first,
                last,
                coordinate_index,
                &mut local_min,
                &mut local_max,
            );
        }
        if !local_min.is_finite() || !local_max.is_finite() {
            continue;
        }
        if local_max > simple_max {
            simple_max = local_max;
            max_face_index = face_index;
        }

        let mut merged = false;
        for interval in &mut intervals {
            if local_min >= interval.0 && local_max <= interval.1 {
                merged = true;
                break;
            }
            if local_min <= interval.1 && local_max >= interval.0 {
                interval.0 = interval.0.min(local_min);
                interval.1 = interval.1.max(local_max);
                merged = true;
                break;
            }
        }
        if !merged {
            if intervals
                .first()
                .is_some_and(|interval| local_max < interval.0)
            {
                intervals.insert(0, (local_min, local_max));
            } else {
                intervals.push((local_min, local_max));
            }
        }
    }

    let count = intervals.len();
    let first = *intervals.first()?;
    let min_coordinate = if count == 2 {
        intervals[1].0 - period
    } else {
        first.0
    };
    Some((min_coordinate, first.1, count, max_face_index))
}

/// `getCurveParams` (`cxx:512-533`).
pub(super) fn curve_end_points(edge: &Edge, reference_face: &Face) -> Option<(GpPnt2d, GpPnt2d)> {
    let (pcurve, mut first, mut last) = curve_on_surface_oriented(edge, reference_face, false)?;
    if edge.0.orientation() != Orientation::Forward {
        std::mem::swap(&mut first, &mut last);
    }
    Some((pcurve.d0(first), pcurve.d0(last)))
}

fn remove_edge_from_vertex_map(edge: &Edge, vertex_edges: &mut VertexEdgeMap) -> bool {
    let mut removed = false;
    let (first, last) = edge_vertices(edge);
    for vertex in [first, last].into_iter().flatten() {
        if let Some((_, edges)) = vertex_edges.get_mut(&shape_map_key(&vertex.0)) {
            let old_len = edges.len();
            edges.retain(|candidate| !is_same(&candidate.0, &edge.0));
            removed |= old_len != edges.len();
        }
    }
    removed
}

fn translated_curve(
    curve: &Arc<dyn Curve2d>,
    coordinate_index: usize,
    offset: f64,
) -> Arc<dyn Curve2d> {
    let mut copy = curve.clone_dyn();
    let mut transformation = GpTrsf2d::identity();
    let vector = if coordinate_index == 1 {
        GpVec2d::new(offset, 0.0)
    } else {
        GpVec2d::new(0.0, offset)
    };
    transformation.set_translation_vec(&vector);
    copy.transform(&transformation);
    Arc::from(copy)
}

/// `RelocatePCurvesToNewUorigin` (`cxx:536-697`).
pub(super) fn relocate_pcurves_to_new_origin(
    edges: &[Edge],
    first_face: &Face,
    reference_face: &Face,
    coordinate_tolerance: f64,
    coordinate_index: usize,
    period: f64,
    vertex_edges: &mut VertexEdgeMap,
) -> Option<(HashMap<ShapeMapKey, Arc<dyn Curve2d>>, HashSet<ShapeMapKey>)> {
    let first_face_edges: HashSet<_> = explorer_shapes(&first_face.0, ShapeType::Edge)
        .iter()
        .map(shape_map_key)
        .collect();
    let mut relocated = HashMap::new();
    let mut used = HashSet::new();

    loop {
        let mut start: Option<(Edge, Orientation, f64)> = None;
        for edge in edges {
            let key = shape_map_key(&edge.0);
            if used.contains(&key) || !first_face_edges.contains(&key) {
                continue;
            }
            let (first, last) = curve_end_points(edge, reference_face)?;
            for (point, orientation) in
                [(first, Orientation::Forward), (last, Orientation::Reversed)]
            {
                let value = coordinate(&point, coordinate_index);
                if start.as_ref().is_none_or(|candidate| value < candidate.2) {
                    start = Some((edge.clone(), orientation, value));
                }
            }
        }
        let Some((mut current, walk_orientation, _)) = start else {
            break;
        };

        let (source, mut edge_start, mut edge_end) =
            curve_on_surface_oriented(&current, reference_face, false)?;
        let trimmed: Arc<dyn Curve2d> =
            Arc::new(Geom2dTrimmedCurve::new(source, edge_start, edge_end));
        relocated.insert(shape_map_key(&current.0), trimmed.clone());
        if current.0.orientation() == Orientation::Reversed {
            std::mem::swap(&mut edge_start, &mut edge_end);
        }
        let parameter = if walk_orientation == Orientation::Forward {
            edge_end
        } else {
            edge_start
        };
        let mut current_point = trimmed.d0(parameter);

        loop {
            if !remove_edge_from_vertex_map(&current, vertex_edges) {
                break;
            }
            used.insert(shape_map_key(&current.0));
            let (first_vertex, last_vertex) = edge_vertices(&current);
            let vertex = if walk_orientation == Orientation::Forward {
                if current.0.orientation() == Orientation::Reversed {
                    first_vertex?
                } else {
                    last_vertex?
                }
            } else if current.0.orientation() == Orientation::Reversed {
                last_vertex?
            } else {
                first_vertex?
            };
            let candidates = vertex_edges
                .get(&shape_map_key(&vertex.0))
                .map(|entry| entry.1.clone())
                .unwrap_or_default();
            let mut next = None;
            for edge in candidates {
                if is_same(&edge.0, &current.0) {
                    continue;
                }
                let (edge_first, edge_last) = edge_vertices(&edge);
                let expected = if walk_orientation == Orientation::Forward {
                    if edge.0.orientation() == Orientation::Reversed {
                        edge_last
                    } else {
                        edge_first
                    }
                } else if edge.0.orientation() == Orientation::Reversed {
                    edge_first
                } else {
                    edge_last
                };
                if !expected.is_some_and(|candidate| is_same(&candidate.0, &vertex.0)) {
                    continue;
                }
                let (curve, mut first, mut last) =
                    curve_on_surface_oriented(&edge, reference_face, false)?;
                let mut candidate: Arc<dyn Curve2d> =
                    Arc::new(Geom2dTrimmedCurve::new(curve, first, last));
                if edge.0.orientation() == Orientation::Reversed {
                    std::mem::swap(&mut first, &mut last);
                }
                let parameter = if walk_orientation == Orientation::Forward {
                    first
                } else {
                    last
                };
                let point = candidate.d0(parameter);
                let offset = coordinate(&current_point, coordinate_index)
                    - coordinate(&point, coordinate_index);
                if offset.abs() >= coordinate_tolerance
                    && (offset.abs() - period).abs() >= coordinate_tolerance
                {
                    continue;
                }
                if offset.abs() > period / 2.0 {
                    candidate = translated_curve(
                        &candidate,
                        coordinate_index,
                        if offset > 0.0 { period } else { -period },
                    );
                }
                relocated.insert(shape_map_key(&edge.0), candidate.clone());
                let composed_forward = walk_orientation == edge.0.orientation();
                current_point = candidate.d0(if composed_forward { last } else { first });
                next = Some(edge);
                break;
            }
            let Some(edge) = next else {
                break;
            };
            current = edge;
        }
    }
    Some((relocated, used))
}

/// `ComputeMinEdgeSize` (`cxx:344-380`).
pub(super) fn minimum_edge_size(
    edges: &[Edge],
    reference_face: &Face,
) -> Option<(f64, HashSet<ShapeMapKey>)> {
    let mut minimum = f64::MAX;
    let mut edge_map = HashSet::new();
    for edge in edges {
        edge_map.insert(shape_map_key(&edge.0));
        let (pcurve, first, last) = curve_on_surface_oriented(edge, reference_face, false)?;
        let first_point = pcurve.d0(first);
        let last_point = pcurve.d0(last);
        let (v1, v2) = edge_vertices(edge);
        let squared = if matches!((v1, v2), (Some(a), Some(b)) if is_same(&a.0, &b.0))
            && !crate::brep_tool::BRepTool::is_degenerated(edge)
        {
            first_point.square_distance(&pcurve.d0((first + last) * 0.5))
        } else {
            first_point.square_distance(&last_point)
        };
        minimum = minimum.min(squared);
    }
    Some((minimum.sqrt().max(CONFUSION), edge_map))
}
