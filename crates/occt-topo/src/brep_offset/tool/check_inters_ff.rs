//! `CheckIntersFF` (`BRepOffset_Tool.cxx:1264-1343`): of all the section edges
//! stored in the pave filler's face/face interference curves, keep only those
//! of the connexity block nearest to `RefEdge`.
//!
//! Control flow, in the reference's order:
//! 1. every `BOPDS_PaveBlock::Edge()` of every `BOPDS_Curve` of every
//!    `BOPDS_InterfFF` is collected into one compound (`cxx:1275-1299`);
//! 2. an empty compound answers an empty `TrueEdges` (`cxx:1301-1304`);
//! 3. `BOPTools_AlgoTools::MakeConnexityBlocks(Edges, VERTEX, EDGE, ...)`
//!    splits the edges into blocks connected through shared vertices
//!    (`cxx:1307`);
//! 4. a single block is taken as is (`cxx:1310-1312`); otherwise the reference
//!    edge's midpoint becomes a vertex and each block's
//!    `BRepExtrema_DistShapeShape` distance selects the nearest one
//!    (`cxx:1314-1337`);
//! 5. `TopExp::MapShapes(NearestCompound, EDGE, TrueEdges)` fills the result
//!    (`cxx:1340`).
//!
//! Port guards (no OCCT counterpart): a pave-block edge index that does not
//! resolve to an edge shape is skipped instead of being cast; when no block
//! yields a distance the result stays empty.

use std::collections::HashMap;

use occt_core::precision::CONFUSION;
use occt_geom::extrema_surf::{ExtremaExtAlgo, ExtremaExtFlag};

use crate::abs::ShapeType;
use crate::bopds::BopdsDS;
use crate::brep_extrema_dist_shape_shape::DistShapeShape;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, TopoShape};
use crate::topo_tools::{explorer_shapes, is_same, shape_map_key, ShapeMapKey};

use super::BRepOffsetTool;

/// `BOPTools_AlgoTools::MakeConnexityBlocks(theS, VERTEX, EDGE, theLCB)`
/// (`BOPTools_AlgoTools.cxx:105-155`, compound variant at `:158-183`):
/// elements (`theElementType` = edges) connected when they share a connection
/// sub-shape (`theConnectionType` = vertices).
///
/// The block order follows the first element of each block, and an element is
/// appended to the growing block in the order of the connection map's lists,
/// exactly like the reference's `aLBlock` walk.
pub(crate) fn make_connexity_blocks(edges: Vec<Edge>) -> Vec<Vec<Edge>> {
    // `NCollection_Map<TopoDS_Shape>` fence (`cxx:114`): `IsSame` dedupe of the
    // elements added to the start compound.
    let mut unique: Vec<Edge> = Vec::new();
    for e in edges {
        if !unique.iter().any(|u| is_same(&u.0, &e.0)) {
            unique.push(e);
        }
    }

    // `TopExp::MapShapesAndAncestors(theS, VERTEX, EDGE, theConnectionMap)`
    // (`cxx:112`): per vertex (IsSame key) the edges that carry it, in element
    // traversal order.
    let vertex_lists: Vec<Vec<TopoShape>> = unique
        .iter()
        .map(|e| explorer_shapes(&e.0, ShapeType::Vertex))
        .collect();
    let mut connection: HashMap<ShapeMapKey, Vec<usize>> = HashMap::new();
    for (i, verts) in vertex_lists.iter().enumerate() {
        for v in verts {
            connection.entry(shape_map_key(v)).or_default().push(i);
        }
    }

    let mut fence = vec![false; unique.len()];
    let mut blocks: Vec<Vec<Edge>> = Vec::new();
    for start in 0..unique.len() {
        // `if (!aMFence.Add(aS)) continue;` (`cxx:121-124`).
        if fence[start] {
            continue;
        }
        fence[start] = true;
        let mut block = vec![start];
        // `for (aItB(aLBlock); aItB.More(); aItB.Next())` (`cxx:132-150`): the
        // list grows while it is scanned.
        let mut cursor = 0;
        while cursor < block.len() {
            for v in &vertex_lists[block[cursor]] {
                let Some(list) = connection.get(&shape_map_key(v)) else {
                    continue;
                };
                for &j in list {
                    if !fence[j] {
                        fence[j] = true;
                        block.push(j);
                    }
                }
            }
            cursor += 1;
        }
        blocks.push(block.into_iter().map(|i| unique[i].clone()).collect());
    }
    blocks
}

impl BRepOffsetTool {
    /// `CheckIntersFF(pDS, RefEdge, TrueEdges)` (`cxx:1264-1343`).
    ///
    /// `true_edges` is cleared first and receives the edges of the nearest
    /// connexity block (all section edges when there is a single block).
    pub fn check_inters_ff(pds: &BopdsDS, ref_edge: &Edge, true_edges: &mut Vec<Edge>) {
        true_edges.clear();

        // `cxx:1275-1299`: collect every pave-block edge into one compound.
        let mut all_edges: Vec<Edge> = Vec::new();
        for a_ff in pds.interf_ff() {
            for a_bc in a_ff.curves() {
                for a_pb in a_bc.pave_blocks() {
                    let Some(shape) = pds.shape(a_pb.edge()) else {
                        continue;
                    };
                    let Some(edge) = Edge::wrap(shape.clone()) else {
                        continue;
                    };
                    all_edges.push(edge);
                }
            }
        }
        if all_edges.is_empty() {
            return;
        }

        // `cxx:1307`.
        let blocks = make_connexity_blocks(all_edges);

        // `cxx:1310-1337`.
        let nearest: Vec<Edge> = if blocks.len() == 1 {
            blocks.into_iter().next().unwrap_or_default()
        } else {
            let Some(curve) = BRepTool::edge_curve_world(ref_edge) else {
                return;
            };
            let (first, last) = BRepTool::edge_parameters(ref_edge);
            // `BRepAdaptor_Curve::Value` at the abscissa midpoint of the
            // adaptor's own range (`cxx:1316-1318`).
            let p_ref = curve.d0(0.5 * (first + last));
            let v_ref = TopoBuilder::new().make_vertex(p_ref, CONFUSION);

            let mut min_dist = f64::MAX;
            let mut nearest: Option<Vec<Edge>> = None;
            for block in blocks {
                let shapes: Vec<TopoShape> = block.iter().map(|e| e.0.clone()).collect();
                let compound = TopoBuilder::new().make_compound_of(&shapes);
                // `BRepExtrema_DistShapeShape(Vref, aCompound)` (`cxx:1326`).
                let projector = DistShapeShape::new(
                    &v_ref.0,
                    &compound.0,
                    CONFUSION,
                    ExtremaExtFlag::MinMax,
                    ExtremaExtAlgo::Grad,
                );
                if !projector.is_done() || projector.nb_solution() == 0 {
                    continue;
                }
                let a_dist = projector.value();
                if a_dist < min_dist {
                    min_dist = a_dist;
                    nearest = Some(block);
                }
            }
            match nearest {
                Some(block) => block,
                None => return,
            }
        };

        // `cxx:1340`: `TopExp::MapShapes(NearestCompound, TopAbs_EDGE, ...)`,
        // an `IndexedMap` so the same edge is stored once.
        for edge in nearest {
            if !true_edges.iter().any(|t| is_same(&t.0, &edge.0)) {
                true_edges.push(edge);
            }
        }
    }
}
