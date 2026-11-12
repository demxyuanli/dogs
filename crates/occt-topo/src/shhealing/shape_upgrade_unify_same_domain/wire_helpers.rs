//! Wire traversal helpers for same-domain face unification.
//!
//! Source: `ShapeUpgrade_UnifySameDomain.cxx:93-119`, `:310-342`,
//! and `:4560-4688`.

use std::collections::{HashMap, HashSet};

use crate::abs::Orientation;
use crate::boptools_2d::curve_on_surface_oriented;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, TopoShape, Wire};
use crate::topo_tools::{edge_vertices, edges_of_wire, is_same, shape_map_key, ShapeMapKey};

use super::uv_bounds::VertexEdgeMap;

fn first_vertex(edge: &Edge) -> Option<TopoShape> {
    let (first, last) = edge_vertices(edge);
    if edge.0.orientation() == Orientation::Reversed {
        last.map(|vertex| vertex.0)
    } else {
        first.map(|vertex| vertex.0)
    }
}

fn last_vertex(edge: &Edge) -> Option<TopoShape> {
    let (first, last) = edge_vertices(edge);
    if edge.0.orientation() == Orientation::Reversed {
        first.map(|vertex| vertex.0)
    } else {
        last.map(|vertex| vertex.0)
    }
}

/// `IsOnSingularity` (`cxx:93-106`).
pub(super) fn is_on_singularity(edges: &[Edge]) -> bool {
    edges.iter().any(BRepTool::is_degenerated)
}

/// `IsUiso` (`cxx:109-119`).
pub(super) fn is_u_iso(edge: &Edge, face: &Face) -> Option<bool> {
    let (pcurve, first, _) = curve_on_surface_oriented(edge, face, false)?;
    let (_, derivative) = pcurve.d1(first);
    Some(derivative.y().abs() > derivative.x().abs())
}

/// `RemoveEdgeFromMap` (`cxx:305-342`).
pub(super) fn remove_edge_from_map(edge: &Edge, vertex_edges: &mut VertexEdgeMap) -> bool {
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

pub(super) fn vertex_edge_map(edges: &[Edge]) -> VertexEdgeMap {
    let mut map = HashMap::new();
    for edge in edges {
        let (first, last) = edge_vertices(edge);
        for vertex in [first, last].into_iter().flatten() {
            map.entry(shape_map_key(&vertex.0))
                .or_insert_with(|| (vertex.0.clone(), Vec::new()))
                .1
                .push(edge.clone());
        }
    }
    map
}

fn directed_angle(edge: &Edge, face: &Face, at_end: bool) -> Option<occt_core::gp::GpVec2d> {
    let (pcurve, first, last) = curve_on_surface_oriented(edge, face, false)?;
    let parameter = if at_end {
        if edge.0.orientation() == Orientation::Forward {
            last
        } else {
            first
        }
    } else if edge.0.orientation() == Orientation::Forward {
        first
    } else {
        last
    };
    let (_, mut direction) = pcurve.d1(parameter);
    if edge.0.orientation() == Orientation::Reversed {
        direction.reverse();
    }
    direction.normalize().ok()?;
    Some(direction)
}

/// `SplitWire` (`cxx:4560-4688`).
pub(super) fn split_wire(wire: &Wire, face: &Face, splitting_vertices: &[TopoShape]) -> Vec<Wire> {
    let mut outgoing: HashMap<ShapeMapKey, Vec<Edge>> = HashMap::new();
    let mut seen = HashSet::new();
    for edge in edges_of_wire(wire) {
        if !seen.insert(shape_map_key(&edge.0))
            || !matches!(
                edge.0.orientation(),
                Orientation::Forward | Orientation::Reversed
            )
        {
            continue;
        }
        let Some(vertex) = first_vertex(&edge) else {
            continue;
        };
        outgoing
            .entry(shape_map_key(&vertex))
            .or_default()
            .push(edge);
    }

    let builder = TopoBuilder::new();
    let mut result = Vec::new();
    for origin in splitting_vertices {
        loop {
            let key = shape_map_key(origin);
            let Some(current) = outgoing.get_mut(&key).and_then(|edges| {
                if edges.is_empty() {
                    None
                } else {
                    Some(edges.remove(0))
                }
            }) else {
                break;
            };
            let mut current = current;
            let mut chain = Vec::new();
            loop {
                chain.push(current.clone());
                let Some(vertex) = last_vertex(&current) else {
                    break;
                };
                if is_same(&vertex, origin) {
                    break;
                }
                let Some(candidates) = outgoing.get_mut(&shape_map_key(&vertex)) else {
                    break;
                };
                if candidates.is_empty() {
                    break;
                }
                let next_index = if candidates.len() == 1 {
                    0
                } else {
                    let Some(current_direction) = directed_angle(&current, face, true) else {
                        break;
                    };
                    let mut minimum = f64::MAX;
                    let mut selected = 0;
                    for (index, candidate) in candidates.iter().enumerate() {
                        let Some(direction) = directed_angle(candidate, face, false) else {
                            continue;
                        };
                        let angle = current_direction.angle(&direction);
                        if angle < minimum {
                            minimum = angle;
                            selected = index;
                        }
                    }
                    selected
                };
                current = candidates.remove(next_index);
            }
            result.push(builder.make_wire(&chain));
        }
    }
    result
}
