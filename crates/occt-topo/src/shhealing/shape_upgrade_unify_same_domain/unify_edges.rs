//! `ShapeUpgrade_UnifySameDomain::UnifyEdges`.
//!
//! Source: `ShapeUpgrade_UnifySameDomain.cxx:4324-4450`.

use std::collections::{HashMap, HashSet};

use occt_core::precision::CONFUSION;

use crate::abs::ShapeType;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face};
use crate::topo_tools::{
    explorer_shapes, map_shapes_and_ancestors, map_shapes_and_unique_ancestors, shape_map_key,
    ShapeMapKey,
};

use super::edge_chain::AncestorMap;
use super::unify_faces::UnifySameDomain;

fn ordered_edges(shape: &crate::shape::TopoShape, edge_faces: &AncestorMap) -> Vec<Edge> {
    let mut seen = HashSet::new();
    explorer_shapes(shape, ShapeType::Edge)
        .into_iter()
        .filter_map(Edge::wrap)
        .filter(|edge| edge_faces.contains_key(&shape_map_key(&edge.0)))
        .filter(|edge| seen.insert(shape_map_key(&edge.0)))
        .collect()
}

fn update_keep_shapes(unifier: &mut UnifySameDomain) {
    let keep: Vec<_> = unifier.keep_shape_values.values().cloned().collect();
    for shape in keep {
        let applied = unifier.context.apply(&shape);
        let key = shape_map_key(&applied);
        unifier.keep_shapes.insert(key.clone());
        unifier.keep_shape_values.insert(key, applied);
    }
}

fn set_fix_wire_modes(fixer: &mut crate::shhealing::ShapeFixFace) {
    fixer.wire_modes.fix_self_intersection_mode = 0;
    fixer.wire_modes.fix_non_adjacent_intersecting_edges_mode = 0;
    fixer.wire_modes.fix_lacking_mode = 0;
    fixer.wire_modes.fix_notched_edges_mode = 0;
    fixer.wire_modes.topo_mode = false;
    fixer.wire_modes.remove_loop_mode = 0;
    fixer.wire_modes.fix_gaps_by_ranges = false;
    fixer.wire_modes.fix_small_mode = 0;
}

impl UnifySameDomain {
    /// `MergeEdges` and `MergeSeq` (`cxx:2756-2905`).
    fn merge_edge_sequence(
        &mut self,
        edges: &[Edge],
        vertex_faces: &AncestorMap,
        vertex_edges: &AncestorMap,
        edge_faces: &AncestorMap,
    ) -> Result<bool, String> {
        let shared = super::edge_chain::shared_vertices(edges, vertex_edges, &self.keep_shapes);
        let mut subsequences = super::edge_chain::split_edge_chains(
            edges,
            vertex_faces,
            &shared,
            self.angular_tolerance,
            self.linear_tolerance,
        );
        super::edge_merge::merge_sequences(
            &mut subsequences,
            vertex_faces,
            &mut self.context,
            self.safe_input_mode,
            self.concat_bsplines,
            self.angular_tolerance,
            edge_faces,
            &self.face_plane_map,
            &self.face_new_face,
        )
    }

    /// `UnifyEdges` (`cxx:4324-4450`).
    pub(crate) fn unify_edges(&mut self) -> Result<(), String> {
        let result = self.context.apply(&self.shape);
        let edge_faces = map_shapes_and_ancestors(&result, ShapeType::Edge, ShapeType::Face);
        let vertex_edges =
            map_shapes_and_unique_ancestors(&result, ShapeType::Vertex, ShapeType::Edge, false);
        let vertex_faces =
            map_shapes_and_unique_ancestors(&result, ShapeType::Vertex, ShapeType::Face, false);
        if self.safe_input_mode {
            update_keep_shapes(self);
        }

        let edges = ordered_edges(&result, &edge_faces);
        let initial_edge_faces = self.initial_edge_faces.clone();
        let merged =
            self.merge_edge_sequence(&edges, &vertex_faces, &vertex_edges, &initial_edge_faces)?;
        let mut changed_faces: HashMap<ShapeMapKey, crate::shape::TopoShape> = HashMap::new();
        if merged {
            for edge in &edges {
                if !self.context.is_recorded(&edge.0) {
                    continue;
                }
                if let Some(faces) = edge_faces.get(&shape_map_key(&edge.0)) {
                    for face in faces {
                        changed_faces
                            .entry(shape_map_key(face))
                            .or_insert_with(|| face.clone());
                    }
                }
            }
        }

        for source_face in changed_faces.values() {
            let applied = self.context.apply(source_face);
            let Some(face) = Face::wrap(applied) else {
                continue;
            };
            if !self.safe_input_mode
                && BRepTool::face_surface_world(&face)
                    .is_some_and(|surface| surface.gp_pln().is_some())
            {
                for edge_shape in explorer_shapes(&face.0, ShapeType::Edge) {
                    if let Some(edge) = Edge::wrap(edge_shape) {
                        crate::boptools_2d::build_pcurve_for_edge_on_face(&edge, &face)?;
                    }
                }
            }

            let mut fixer = crate::shhealing::ShapeFixFace::with_face(&face);
            fixer.precision = CONFUSION;
            fixer.max_tol = 1.0_f64.max(CONFUSION * 1000.0);
            fixer.modes.fix_orientation_mode = 0;
            fixer.modes.fix_add_natural_bound_mode = 0;
            fixer.modes.fix_intersecting_wires_mode = 0;
            fixer.modes.fix_loop_wires_mode = 0;
            fixer.modes.fix_split_face_mode = 0;
            fixer.modes.fix_periodic_degenerated = 0;
            set_fix_wire_modes(&mut fixer);
            fixer.perform();
            let fixed = fixer.result.unwrap_or_else(|| face.0.clone());
            self.context.replace(&face.0, &fixed)?;
        }

        if !changed_faces.is_empty() {
            let mut rebuilt = self.context.apply(&result);
            let mut shell_changed = false;
            for shell in explorer_shapes(&rebuilt, ShapeType::Shell) {
                let fixed = crate::shhealing::fix_shell_face_orientation(shell.clone());
                if !crate::topo_tools::is_same(&fixed, &shell) {
                    self.context.replace(&shell, &fixed)?;
                    shell_changed = true;
                }
            }
            if shell_changed {
                rebuilt = self.context.apply(&rebuilt);
            }
            self.context.replace(&self.shape, &rebuilt)?;
        }
        self.shape = self.context.apply(&self.shape);
        Ok(())
    }
}
