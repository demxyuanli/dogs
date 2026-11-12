//! Boundary and internal wire reconstruction for same-domain face unification.
//!
//! Source: `ShapeUpgrade_UnifySameDomain.cxx:3867-4228`.

use occt_core::precision::CONFUSION;

use crate::abs::Orientation;
use crate::boptools_2d::curve_on_surface_oriented;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, TopoShape, Wire};
use crate::topo_tools::{edge_vertices, is_same};

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

fn tangent(edge: &Edge, face: &Face, at_end: bool) -> Option<occt_core::gp::GpVec2d> {
    let (curve, first, last) = curve_on_surface_oriented(edge, face, false)?;
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
    let (_, mut direction) = curve.d1(parameter);
    direction.normalize().ok()?;
    if edge.0.orientation() == Orientation::Reversed {
        direction.reverse();
    }
    Some(direction)
}

pub(super) fn build_boundary_wires(
    edges: &[Edge],
    reference_face: &Face,
    removed_edges: &[Edge],
    u_period: f64,
    v_period: f64,
) -> Result<Vec<Wire>, String> {
    let builder = TopoBuilder::new();
    let mut remaining = edges.to_vec();
    let mut result = Vec::new();
    while !remaining.is_empty() {
        let start_index = remaining
            .iter()
            .position(|edge| !BRepTool::is_degenerated(edge))
            .unwrap_or(0);
        let start = remaining.remove(start_index);
        let start_vertex =
            first_vertex(&start).ok_or("IntUnifyFaces: boundary edge has no first vertex")?;
        let mut current = start.clone();
        let (start_curve, start_first, start_last) =
            curve_on_surface_oriented(&start, reference_face, false)
                .ok_or("IntUnifyFaces: boundary edge has no pcurve")?;
        let start_parameter = if start.0.orientation() == Orientation::Forward {
            start_first
        } else {
            start_last
        };
        let current_parameter = if start.0.orientation() == Orientation::Forward {
            start_last
        } else {
            start_first
        };
        let start_point = start_curve.d0(start_parameter);
        let mut current_point = start_curve.d0(current_parameter);
        let mut chain = vec![start];
        let mut splitting_vertices = Vec::new();
        loop {
            let current_vertex =
                last_vertex(&current).ok_or("IntUnifyFaces: boundary edge has no last vertex")?;
            if is_same(&current_vertex, &start_vertex) {
                let across_period = (u_period != 0.0
                    && (start_point.x() - current_point.x()).abs() > u_period / 2.0)
                    || (v_period != 0.0
                        && (start_point.y() - current_point.y()).abs() > v_period / 2.0);
                if !across_period {
                    break;
                }
                let Some((seam, next_point)) = super::face_helpers::reconstruct_missed_seam(
                    removed_edges,
                    reference_face,
                    &current,
                    &current_vertex,
                    current_point,
                    u_period,
                    v_period,
                ) else {
                    return Err("IntUnifyFaces: missed closing seam was not reconstructed".into());
                };
                current = seam;
                current_point = next_point;
                chain.push(current.clone());
                continue;
            }
            let candidates: Vec<usize> = remaining
                .iter()
                .enumerate()
                .filter_map(|(index, edge)| {
                    first_vertex(edge)
                        .filter(|vertex| is_same(vertex, &current_vertex))
                        .map(|_| index)
                })
                .collect();
            if candidates.len() > 1
                && !super::wire_helpers::is_on_singularity(
                    &candidates
                        .iter()
                        .map(|index| remaining[*index].clone())
                        .collect::<Vec<_>>(),
                )
            {
                splitting_vertices.push(current_vertex.clone());
            }
            let next_index = match candidates.as_slice() {
                [] => {
                    if u_period == 0.0 && v_period == 0.0 {
                        return Err(
                            "IntUnifyFaces: open non-periodic boundary; OCCT lines 4100-4108"
                                .into(),
                        );
                    }
                    let Some((seam, next_point)) = super::face_helpers::reconstruct_missed_seam(
                        removed_edges,
                        reference_face,
                        &current,
                        &current_vertex,
                        current_point,
                        u_period,
                        v_period,
                    ) else {
                        return Err("IntUnifyFaces: missed seam was not reconstructed".into());
                    };
                    current = seam;
                    current_point = next_point;
                    chain.push(current.clone());
                    continue;
                }
                [only] => *only,
                many if u_period != 0.0 || v_period != 0.0 => {
                    let mut selected = None;
                    for candidate in many {
                        let edge = &remaining[*candidate];
                        let (curve, first, last) =
                            curve_on_surface_oriented(edge, reference_face, false)
                                .ok_or("IntUnifyFaces: candidate edge has no pcurve")?;
                        let parameter = if edge.0.orientation() == Orientation::Forward {
                            first
                        } else {
                            last
                        };
                        let point = curve.d0(parameter);
                        let diff_u = (point.x() - current_point.x()).abs();
                        let diff_v = (point.y() - current_point.y()).abs();
                        let u_matches = u_period == 0.0
                            || diff_u <= CONFUSION
                            || (diff_u - u_period).abs() <= CONFUSION;
                        let v_matches = v_period == 0.0
                            || diff_v <= CONFUSION
                            || (diff_v - v_period).abs() <= CONFUSION;
                        if u_matches && v_matches {
                            selected = Some(*candidate);
                            break;
                        }
                    }
                    selected.ok_or("IntUnifyFaces: no periodic boundary candidate matches")?
                }
                many => {
                    let direction = tangent(&current, reference_face, true)
                        .ok_or("IntUnifyFaces: current edge has no pcurve tangent")?;
                    let mut maximum = f64::NEG_INFINITY;
                    let mut selected = many[0];
                    for candidate in many {
                        let next_direction = tangent(&remaining[*candidate], reference_face, false)
                            .ok_or("IntUnifyFaces: candidate edge has no pcurve tangent")?;
                        let angle = direction.angle(&next_direction);
                        if angle > maximum {
                            maximum = angle;
                            selected = *candidate;
                        }
                    }
                    selected
                }
            };
            let next = remaining[next_index].clone();
            let (next_curve, next_first, next_last) =
                curve_on_surface_oriented(&next, reference_face, false)
                    .ok_or("IntUnifyFaces: next edge has no pcurve")?;
            let first_parameter = if next.0.orientation() == Orientation::Forward {
                next_first
            } else {
                next_last
            };
            let first_point = next_curve.d0(first_parameter);
            let across_period = (u_period != 0.0
                && (first_point.x() - current_point.x()).abs() > u_period / 2.0)
                || (v_period != 0.0
                    && (first_point.y() - current_point.y()).abs() > v_period / 2.0);
            if across_period {
                let Some((seam, next_point)) = super::face_helpers::reconstruct_missed_seam(
                    removed_edges,
                    reference_face,
                    &current,
                    &current_vertex,
                    current_point,
                    u_period,
                    v_period,
                ) else {
                    return Err("IntUnifyFaces: periodic seam reconstruction failed".into());
                };
                current = seam;
                current_point = next_point;
                chain.push(current.clone());
                continue;
            }
            current = remaining.remove(next_index);
            let last_parameter = if current.0.orientation() == Orientation::Forward {
                next_last
            } else {
                next_first
            };
            current_point = next_curve.d0(last_parameter);
            chain.push(current.clone());
        }
        let wire = builder.make_wire(&chain);
        wire.0.set_closed(true);
        if splitting_vertices.is_empty() {
            result.push(wire);
        } else {
            result.extend(super::wire_helpers::split_wire(
                &wire,
                reference_face,
                &splitting_vertices,
            ));
        }
    }
    Ok(result)
}

pub(super) fn build_internal_wires(edges: &[Edge]) -> Vec<Wire> {
    let builder = TopoBuilder::new();
    let mut remaining = edges.to_vec();
    let mut result = Vec::new();
    while !remaining.is_empty() {
        let first = remaining.remove(0);
        let (mut left, mut right) = edge_vertices(&first);
        let mut component = vec![first];
        loop {
            let (Some(left_vertex), Some(right_vertex)) = (&left, &right) else {
                break;
            };
            if is_same(&left_vertex.0, &right_vertex.0) {
                break;
            }
            let mut found = false;
            for side in 0..2 {
                let endpoint = if side == 0 {
                    left.as_ref()
                } else {
                    right.as_ref()
                };
                let Some(endpoint) = endpoint else {
                    continue;
                };
                let Some(index) = remaining.iter().position(|edge| {
                    let (first, last) = edge_vertices(edge);
                    first
                        .as_ref()
                        .is_some_and(|vertex| is_same(&vertex.0, &endpoint.0))
                        || last
                            .as_ref()
                            .is_some_and(|vertex| is_same(&vertex.0, &endpoint.0))
                }) else {
                    continue;
                };
                let edge = remaining.remove(index);
                let (first_vertex, last_vertex) = edge_vertices(&edge);
                let other = match (first_vertex, last_vertex) {
                    (Some(first_vertex), Some(last_vertex))
                        if is_same(&first_vertex.0, &endpoint.0) =>
                    {
                        Some(last_vertex)
                    }
                    (Some(first_vertex), Some(_)) => Some(first_vertex),
                    _ => None,
                };
                if side == 0 {
                    left = other;
                } else {
                    right = other;
                }
                component.push(edge);
                found = true;
            }
            if !found {
                break;
            }
        }
        result.push(builder.make_wire(&component));
    }
    result
}
