//! `ShapeFix::FixVertexPosition` (`ShapeFix.cxx:310-787`).

use std::collections::HashMap;

use occt_core::gp::GpPnt;
use occt_core::precision::CONFUSION;

use super::prelude::*;
use super::{copy_replace_vertices_with, first_vertex, last_vertex};
use crate::brep_tool::BRepTool;
use crate::shape_fix_compose_shell::{ReShape, SharedReShape};
use crate::topo_tools::{is_same, shape_map_key, ShapeMapKey};

struct VertexEdges {
    vertex: Vertex,
    edges: Vec<Edge>,
}

fn midpoint(a: &GpPnt, b: &GpPnt) -> GpPnt {
    GpPnt::new(
        (a.x() + b.x()) * 0.5,
        (a.y() + b.y()) * 0.5,
        (a.z() + b.z()) * 0.5,
    )
}

fn average(a: &GpPnt, b: &GpPnt) -> GpPnt {
    midpoint(a, b)
}

fn get_near_point(seq1: &[GpPnt], seq2: &[GpPnt], center: &mut GpPnt) -> f64 {
    let mut min_dist = f64::MAX;
    let mut nearest = None;
    for p1 in seq1 {
        for p2 in seq2 {
            let distance = p1.distance(p2);
            if (distance - min_dist).abs() <= CONFUSION {
                continue;
            }
            if distance < min_dist {
                min_dist = distance;
                nearest = Some((*p1, *p2));
            }
        }
    }
    if let Some((p1, p2)) = nearest {
        *center = midpoint(&p1, &p2);
    }
    min_dist
}

fn endpoint_points(edge: &Edge, vertex: &Vertex) -> Option<(GpPnt, Option<GpPnt>)> {
    let first_vertex = first_vertex(edge)?;
    let last_vertex = last_vertex(edge)?;
    let curve = BRepTool::edge_curve(edge)?;
    let (first, last) = BRepTool::edge_parameters(edge);
    let is_first = is_same(&vertex.0, &first_vertex.0);
    let is_same_vertex = is_same(&first_vertex.0, &last_vertex.0);
    let point = if is_first {
        curve.d0(first)
    } else if !is_same_vertex {
        curve.d0(last)
    } else {
        return None;
    };
    let second = is_same_vertex.then(|| curve.d0(last));
    Some((point, second))
}

fn get_nearest_edges(
    edges: &mut Vec<Edge>,
    vertex: &Vertex,
    suitable: &mut Vec<Edge>,
    rejected: &mut Vec<Edge>,
    tolerance: f64,
    suitable_center: &mut GpPnt,
    rejected_center: &mut GpPnt,
) -> bool {
    if edges.is_empty() {
        return false;
    }
    let edge1 = edges[0].clone();
    let Some((p11, p12)) = endpoint_points(&edge1, vertex) else {
        return false;
    };
    let edge1_vertices = (first_vertex(&edge1), last_vertex(&edge1));
    let mut seen = vec![shape_map_key(&edge1.0)];
    let mut local_suitable = Vec::new();
    let mut local_rejected = Vec::new();
    let mut remaining = edges[1..].to_vec();
    let mut deferred_loops = 0usize;
    let index = 0usize;
    while !remaining.is_empty() {
        let edge = remaining[index].clone();
        let key = shape_map_key(&edge.0);
        if seen.contains(&key) {
            remaining.remove(index);
            continue;
        }
        let Some(v1) = first_vertex(&edge) else {
            remaining.remove(index);
            continue;
        };
        let Some(v2) = last_vertex(&edge) else {
            remaining.remove(index);
            continue;
        };
        let is_loop = match (&edge1_vertices.0, &edge1_vertices.1) {
            (Some(e11), Some(e12)) => {
                (is_same(&v1.0, &e11.0) && is_same(&v2.0, &e12.0))
                    || (is_same(&v1.0, &e12.0) && is_same(&v2.0, &e11.0))
            }
            _ => false,
        };
        if is_loop && remaining.len() > deferred_loops {
            let edge = remaining.remove(index);
            remaining.push(edge);
            deferred_loops += 1;
            continue;
        }
        seen.push(key);
        if let Some((p1, p2)) = endpoint_points(&edge, vertex) {
            let mut center = *suitable_center;
            let min_dist = if p12.is_none() && p2.is_none() {
                center = midpoint(&p1, &p11);
                p1.distance(&p11)
            } else {
                let seq1 = p12.map_or_else(|| vec![p11], |p| vec![p11, p]);
                let seq2 = p2.map_or_else(|| vec![p1], |p| vec![p1, p]);
                get_near_point(&seq1, &seq2, &mut center)
            };
            if min_dist > tolerance {
                if local_rejected.is_empty() {
                    *rejected_center = center;
                }
                local_rejected.push(edge);
            } else if local_suitable.is_empty() {
                *suitable_center = center;
                local_suitable.push(edge);
            } else if p12.is_none() || suitable_center.distance(&center) < tolerance {
                local_suitable.push(edge);
            } else {
                local_rejected.push(edge);
            }
        }
        remaining.remove(index);
    }

    let done = !local_suitable.is_empty() || !local_rejected.is_empty();
    if done {
        if local_suitable.is_empty() {
            rejected.push(edge1);
            edges.remove(0);
            let _ = get_nearest_edges(
                edges,
                vertex,
                suitable,
                rejected,
                tolerance,
                suitable_center,
                rejected_center,
            );
        } else {
            suitable.push(edge1);
            suitable.extend(local_suitable);
            rejected.extend(local_rejected);
        }
    } else {
        rejected.push(edge1);
    }
    done
}

fn replace_vertex(edge: &Edge, point: GpPnt, first: bool) -> Edge {
    let vertex = TopoBuilder::new().make_vertex(point, CONFUSION);
    if first {
        copy_replace_vertices_with(edge, Some(&vertex), None)
    } else {
        copy_replace_vertices_with(edge, None, Some(&vertex))
    }
}

pub fn fix_vertex_position(
    shape: &mut TopoShape,
    tolerance: f64,
    context: &mut SharedReShape,
) -> bool {
    let mut entries = Vec::<VertexEdges>::new();
    let mut index = HashMap::<ShapeMapKey, usize>::new();
    for edge in edges_of(shape) {
        let endpoints = [first_vertex(&edge), last_vertex(&edge)];
        for vertex in endpoints.into_iter().flatten() {
            let key = shape_map_key(&vertex.0);
            let slot = if let Some(slot) = index.get(&key).copied() {
                slot
            } else {
                let slot = entries.len();
                index.insert(key, slot);
                entries.push(VertexEdges {
                    vertex: vertex.clone(),
                    edges: Vec::new(),
                });
                slot
            };
            if entries[slot]
                .edges
                .last()
                .is_none_or(|previous| !is_same(&previous.0, &edge.0))
            {
                entries[slot].edges.push(edge.clone());
            }
        }
    }

    let mut done = false;
    for entry in entries {
        let vertex = entry.vertex;
        if BRepTool::vertex_tolerance(&vertex) <= tolerance {
            continue;
        }
        vertex.set_tolerance(tolerance);
        if entry.edges.len() == 1 {
            continue;
        }
        let old_point = BRepTool::vertex_point(&vertex);
        let mut suitable_center = old_point;
        let mut rejected_center = old_point;
        let mut suitable = Vec::new();
        let mut rejected = Vec::new();
        let mut work = entry.edges.clone();
        if !get_nearest_edges(
            &mut work,
            &vertex,
            &mut suitable,
            &mut rejected,
            tolerance,
            &mut suitable_center,
            &mut rejected_center,
        ) {
            continue;
        }

        let mut added = false;
        for (position, old_edge) in suitable.iter().enumerate() {
            let (Some(old_first), Some(old_last)) = (first_vertex(old_edge), last_vertex(old_edge))
            else {
                continue;
            };
            let is_first = is_same(&old_first.0, &vertex.0);
            let is_last = is_same(&old_last.0, &vertex.0);
            if !is_first && !is_last {
                continue;
            }
            let edge = Edge(context.apply(&old_edge.0));
            let (Some(new_first), Some(new_last), Some(curve)) = (
                first_vertex(&edge),
                last_vertex(&edge),
                BRepTool::edge_curve(&edge),
            ) else {
                continue;
            };
            let (first, last) = BRepTool::edge_parameters(&edge);
            let p1 = curve.d0(first);
            let p2 = curve.d0(last);
            let replace_same = is_same(&new_first.0, &new_last.0) && p1.distance(&p2) > tolerance;
            let endpoint = if is_first { p1 } else { p2 };
            if position + 1 > 2 {
                suitable_center = average(&suitable_center, &endpoint);
            }
            if replace_same {
                let replacement = if is_first {
                    if p1.distance(&suitable_center) < p2.distance(&suitable_center) {
                        replace_vertex(&edge, p2, false)
                    } else {
                        replace_vertex(&edge, p1, true)
                    }
                } else if p1.distance(&suitable_center) < p2.distance(&suitable_center) {
                    replace_vertex(&edge, p2, false)
                } else {
                    replace_vertex(&edge, p1, true)
                };
                context.replace(&edge.0, &replacement.0);
                done = true;
            }
            added = true;
        }

        if added && old_point.distance(&suitable_center) > tolerance {
            let mut new_vertex = TopoBuilder::new().make_vertex(suitable_center, CONFUSION);
            new_vertex.0.set_orientation(vertex.0.orientation());
            context.replace(&vertex.0, &new_vertex.0);
            done = true;
        }

        for old_edge in rejected {
            let (Some(old_first), Some(old_last)) =
                (first_vertex(&old_edge), last_vertex(&old_edge))
            else {
                continue;
            };
            let is_first = is_same(&old_first.0, &vertex.0);
            let is_last = is_same(&old_last.0, &vertex.0);
            if !is_first && !is_last {
                continue;
            }
            let same = is_same(&old_first.0, &old_last.0);
            let edge = Edge(context.apply(&old_edge.0));
            let Some(curve) = BRepTool::edge_curve(&edge) else {
                continue;
            };
            let (first, last) = BRepTool::edge_parameters(&edge);
            let p1 = curve.d0(first);
            let p2 = curve.d0(last);
            let mut replacement = if is_first {
                replace_vertex(&edge, p1, true)
            } else {
                replace_vertex(&edge, p2, false)
            };
            if same {
                replacement = if is_first {
                    replace_vertex(&replacement, p2, false)
                } else {
                    replace_vertex(&replacement, p1, true)
                };
            }
            context.replace(&edge.0, &replacement.0);
            done = true;
        }
    }
    if done {
        *shape = context.apply(shape);
    }
    done
}
