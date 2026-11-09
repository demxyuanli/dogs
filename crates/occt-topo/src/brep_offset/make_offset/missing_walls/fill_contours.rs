use super::prelude::*;
use super::statics::*;


/// `FillContours(aShape, Analyser, Contours, MapEF)` (`cxx:536-599`).
///
/// Answers the contours in the insertion order of the reference's
/// `NCollection_IndexedDataMap` (`cxx:597`) together with the edge -> face map
/// of `cxx:565`.
pub(crate) fn fill_contours(
    shape: &TopoShape,
    analyser: &Analyse,
    contours: &mut Vec<(TopoShape, Vec<TopoShape>)>,
    map_ef: &mut HashMap<ShapeMapKey, TopoShape>,
) {
    let mut edges: Vec<TopoShape> = Vec::new();

    // `cxx:548-570`: the free-boundary edges of every face wire, in wire order.
    for face_shape in explorer_shapes(shape, ShapeType::Face) {
        let Some(face) = Face::wrap(face_shape) else {
            continue;
        };
        // `TopoDS_Iterator itf(aFace)` (`cxx:551`): the face's own wires.
        for wire_shape in children_keep_orientation(&face.0) {
            if wire_shape.shape_type() != ShapeType::Wire {
                continue;
            }
            let Some(wire) = Wire::wrap(wire_shape) else {
                continue;
            };
            // `BRepTools_WireExplorer Wexp(wire, face)` (`cxx:555`).
            for edge in wire_edges_explorer(&wire, &face) {
                if BRepTool::is_degenerated(&edge) {
                    continue;
                }
                let intervals = analyser.type_of(&edge);
                let Some(first) = intervals.first() else {
                    continue;
                };
                if is_free_bound(first.kind()) {
                    map_ef.insert(shape_map_key(&edge.0), face.0.clone());
                    edges.push(edge.0);
                }
            }
        }
    }

    // `cxx:572-598`.
    while !edges.is_empty() {
        let start_edge = edges.remove(0);
        let Some(start) = Edge::wrap(start_edge) else {
            continue;
        };
        // `TopExp::Vertices(StartEdge, StartVertex, CurVertex, true)`.
        let (start_vertex, cur_vertex) = edge_vertices(&start);
        let (Some(start_vertex), Some(mut cur_vertex)) = (start_vertex, cur_vertex) else {
            continue;
        };
        let mut contour: Vec<TopoShape> = vec![start.0.clone()];
        // `cxx:581`: walk until the contour closes back on its start vertex.
        while !cur_vertex.0.same_tshape(&start_vertex.0) {
            let mut found: Option<usize> = None;
            let mut next_vertex: Option<Vertex> = None;
            for (i, e_shape) in edges.iter().enumerate() {
                let Some(candidate) = Edge::wrap(e_shape.clone()) else {
                    continue;
                };
                let (v1, v2) = edge_vertices(&candidate);
                let hit1 = v1
                    .as_ref()
                    .map(|v| v.0.same_tshape(&cur_vertex.0))
                    .unwrap_or(false);
                let hit2 = v2
                    .as_ref()
                    .map(|v| v.0.same_tshape(&cur_vertex.0))
                    .unwrap_or(false);
                if hit1 || hit2 {
                    found = Some(i);
                    next_vertex = if hit1 { v2 } else { v1 };
                    break;
                }
            }
            let Some(i) = found else {
                // The reference spins on this `while` when no edge touches the
                // current vertex (`cxx:581-596` has no exit); the port stops so
                // an open contour cannot hang the offset.
                break;
            };
            contour.push(edges.remove(i));
            match next_vertex {
                Some(v) => cur_vertex = v,
                None => break,
            }
        }
        contours.push((start_vertex.0.clone(), contour));
    }
}
