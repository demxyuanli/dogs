//! `ShapeUpgrade_UnifySameDomain::UnifyFaces` and `IntUnifyFaces`.
//!
//! Source: `ShapeUpgrade_UnifySameDomain.cxx:3047-4321`.

use std::collections::{HashMap, HashSet};

use occt_core::precision::{ANGULAR, CONFUSION};

use crate::abs::{Orientation, ShapeType};
use crate::bop_hist::BopHistory;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Compound, Edge, Face, TopoShape};
use crate::topo_tools::{
    children_keep_orientation, explorer_shapes, explorer_shapes_avoiding, is_same,
    map_shapes_and_ancestors, map_shapes_and_unique_ancestors_into, shape_map_key, ShapeMapKey,
};

use super::context::UnifyContext;
use super::same_domain::FacePlaneMap;
use super::uv_bounds::EdgeFaceMap;

type FaceShellMap = HashMap<ShapeMapKey, HashSet<ShapeMapKey>>;

/// Callable same-domain unification state.
pub(crate) struct UnifySameDomain {
    initial_shape: TopoShape,
    pub(super) shape: TopoShape,
    pub(super) context: UnifyContext,
    history: BopHistory,
    unify_faces_enabled: bool,
    unify_edges_enabled: bool,
    pub(super) linear_tolerance: f64,
    pub(super) angular_tolerance: f64,
    pub(super) concat_bsplines: bool,
    allow_internal: bool,
    pub(super) safe_input_mode: bool,
    pub(super) keep_shapes: HashSet<ShapeMapKey>,
    pub(super) keep_shape_values: HashMap<ShapeMapKey, TopoShape>,
    pub(super) face_plane_map: FacePlaneMap,
    pub(super) face_new_face: HashMap<ShapeMapKey, Face>,
    pub(super) initial_edge_faces: EdgeFaceMap,
}

impl UnifySameDomain {
    pub(crate) fn new(shape: TopoShape) -> Self {
        // `ShapeUpgrade_UnifySameDomain(shape)` defaults: unify edges+faces,
        // `ConcatBSplines = false` (`ShapeUpgrade_UnifySameDomain.cxx`).
        Self::with_flags(shape, true, true, false)
    }

    /// `ShapeUpgrade_UnifySameDomain(shape, UnifyEdges, UnifyFaces, ConcatBSplines)`
    /// (`ShapeUpgrade_UnifySameDomain.hxx`).
    pub(crate) fn with_flags(
        shape: TopoShape,
        unify_edges: bool,
        unify_faces: bool,
        concat_bsplines: bool,
    ) -> Self {
        Self {
            initial_shape: shape.clone(),
            shape,
            context: UnifyContext::new(),
            history: BopHistory::new(),
            unify_faces_enabled: unify_faces,
            unify_edges_enabled: unify_edges,
            linear_tolerance: CONFUSION,
            angular_tolerance: ANGULAR,
            concat_bsplines,
            allow_internal: false,
            safe_input_mode: true,
            keep_shapes: HashSet::new(),
            keep_shape_values: HashMap::new(),
            face_plane_map: FacePlaneMap::new(),
            face_new_face: HashMap::new(),
            initial_edge_faces: EdgeFaceMap::new(),
        }
    }

    pub(crate) fn shape(&self) -> &TopoShape {
        &self.shape
    }

    pub(crate) fn history(&self) -> &BopHistory {
        &self.history
    }

    /// `ShapeUpgrade_UnifySameDomain::Build` (`cxx:4454-4469`).
    pub(crate) fn build(&mut self) -> Result<(), String> {
        self.initial_edge_faces =
            map_shapes_and_ancestors(&self.initial_shape, ShapeType::Edge, ShapeType::Face);
        if self.unify_faces_enabled {
            self.unify_faces()?;
        }
        if self.unify_edges_enabled {
            self.unify_edges()?;
        }
        super::history::fill_history(
            &self.initial_shape,
            &self.shape,
            self.context.history(),
            &mut self.history,
        )
    }

    pub(crate) fn set_linear_tolerance(&mut self, tolerance: f64) {
        self.linear_tolerance = tolerance;
    }

    pub(crate) fn set_angular_tolerance(&mut self, tolerance: f64) {
        self.angular_tolerance = tolerance;
    }

    pub(crate) fn set_concat_bsplines(&mut self, concat: bool) {
        self.concat_bsplines = concat;
    }

    pub(crate) fn set_safe_input_mode(&mut self, safe: bool) {
        self.safe_input_mode = safe;
    }

    pub(crate) fn set_allow_internal(&mut self, allow: bool) {
        self.allow_internal = allow;
    }

    pub(crate) fn keep_shape(&mut self, shape: &TopoShape) {
        if !matches!(shape.shape_type(), ShapeType::Edge | ShapeType::Vertex) {
            return;
        }
        let key = shape_map_key(shape);
        self.keep_shapes.insert(key.clone());
        self.keep_shape_values.insert(key, shape.clone());
    }

    /// `UnifyFaces` (`cxx:3047-3129`).
    pub(crate) fn unify_faces(&mut self) -> Result<(), String> {
        let global_edge_faces =
            map_shapes_and_ancestors(&self.shape, ShapeType::Edge, ShapeType::Face);
        let mut face_shells: FaceShellMap = HashMap::new();
        for shell in explorer_shapes(&self.shape, ShapeType::Shell) {
            let shell_key = shape_map_key(&shell);
            for face in children_keep_orientation(&shell) {
                if face.shape_type() == ShapeType::Face {
                    face_shells
                        .entry(shape_map_key(&face))
                        .or_default()
                        .insert(shell_key.clone());
                }
            }
        }

        let mut free_boundaries = HashSet::new();
        for shell in explorer_shapes_avoiding(&self.shape, ShapeType::Shell, ShapeType::Solid) {
            let edge_faces = map_shapes_and_ancestors(&shell, ShapeType::Edge, ShapeType::Face);
            for (key, faces) in edge_faces {
                let edge = explorer_shapes(&shell, ShapeType::Edge)
                    .into_iter()
                    .find(|edge| shape_map_key(edge) == key)
                    .and_then(Edge::wrap);
                if edge.is_some_and(|edge| !BRepTool::is_degenerated(&edge) && faces.len() == 1) {
                    free_boundaries.insert(key);
                }
            }
        }

        let shells = explorer_shapes(&self.shape, ShapeType::Shell);
        for shell in shells {
            self.int_unify_faces(&shell, &global_edge_faces, &face_shells, &free_boundaries)?;
        }

        let builder = TopoBuilder::new();
        let mut loose = Compound(builder.make_shape(ShapeType::Compound));
        let loose_faces = explorer_shapes_avoiding(&self.shape, ShapeType::Face, ShapeType::Shell);
        for face in &loose_faces {
            builder.add_compound(&mut loose, face);
        }
        if !loose_faces.is_empty() {
            self.int_unify_faces(
                &loose.0,
                &global_edge_faces,
                &FaceShellMap::new(),
                &free_boundaries,
            )?;
        }
        self.shape = self.context.apply(&self.shape);
        Ok(())
    }

    /// `IntUnifyFaces` (`cxx:3184-4321`).
    fn int_unify_faces(
        &mut self,
        input: &TopoShape,
        global_edge_faces: &EdgeFaceMap,
        face_shells: &FaceShellMap,
        free_boundaries: &HashSet<ShapeMapKey>,
    ) -> Result<(), String> {
        let local_edge_faces = map_shapes_and_ancestors(input, ShapeType::Edge, ShapeType::Face);
        let mut processed = HashSet::new();
        for face_shape in explorer_shapes(input, ShapeType::Face) {
            let face = Face::wrap(face_shape).ok_or("IntUnifyFaces: invalid face")?;
            if processed.contains(&shape_map_key(&face.0)) {
                continue;
            }
            let Some(base_surface) = BRepTool::face_surface_world(&face) else {
                continue;
            };
            let mut base_surface = super::same_surface::clear_rectangular_trim(&base_surface);
            let mut boundary_shapes = Vec::new();
            let mut removed_shapes = Vec::new();
            super::boundaries::add_ordinary_edges(
                &mut boundary_shapes,
                &face.0,
                &mut removed_shapes,
            );
            let mut faces = vec![face.clone()];
            let reference_orientation = face.0.orientation();
            let mut index = 0;
            while index < boundary_shapes.len() {
                let edge = Edge::wrap(boundary_shapes[index].clone())
                    .ok_or("IntUnifyFaces: boundary is not an edge")?;
                if BRepTool::is_degenerated(&edge) {
                    index += 1;
                    continue;
                }
                let key = shape_map_key(&edge.0);
                let global_faces = global_edge_faces
                    .get(&key)
                    .map(Vec::as_slice)
                    .unwrap_or(&[]);
                if global_faces.len() != 2
                    || self.keep_shapes.contains(&key)
                    || free_boundaries.contains(&key)
                {
                    index += 1;
                    continue;
                }
                let local_faces = local_edge_faces.get(&key).map(Vec::as_slice).unwrap_or(&[]);
                if local_faces.len() < 2 {
                    index += 1;
                    continue;
                }

                let range = BRepTool::edge_parameters(&edge);
                let middle = (range.0 + range.1) * 0.5;
                let normal = super::same_domain::normal_to_surface_on_edge(&face, &edge, middle);
                let mut merged = false;
                for candidate_shape in local_faces {
                    if is_same(candidate_shape, &face.0)
                        || processed.contains(&shape_map_key(candidate_shape))
                    {
                        continue;
                    }
                    let candidate = Face::wrap(candidate_shape.clone())
                        .ok_or("IntUnifyFaces: ancestor is not a face")?;
                    if face_shells.get(&shape_map_key(&face.0))
                        != face_shells.get(&shape_map_key(&candidate.0))
                    {
                        continue;
                    }
                    if let (Some(first), Some(second)) = (
                        normal,
                        super::same_domain::normal_to_surface_on_edge(&candidate, &edge, middle),
                    ) {
                        if first.angle(&second) > self.angular_tolerance {
                            continue;
                        }
                    }
                    if !super::same_domain::is_same_domain(
                        &face,
                        &candidate,
                        self.linear_tolerance,
                        self.angular_tolerance,
                        &mut self.face_plane_map,
                    ) {
                        continue;
                    }
                    if let Some(changed_index) = super::boundaries::add_ordinary_edges(
                        &mut boundary_shapes,
                        &candidate.0,
                        &mut removed_shapes,
                    ) {
                        index = changed_index;
                    }
                    faces.push(candidate.clone());
                    processed.insert(shape_map_key(&candidate.0));
                    merged = true;
                    break;
                }
                if !merged {
                    index += 1;
                }
            }
            if faces.len() <= 1 {
                continue;
            }
            if let Some(common_plane) = self.face_plane_map.get(&shape_map_key(&faces[0].0)) {
                base_surface = common_plane.clone();
            }

            let mut selected_edge_faces = HashMap::new();
            let selected_compound = {
                let builder = TopoBuilder::new();
                let mut compound = Compound(builder.make_shape(ShapeType::Compound));
                for selected in &faces {
                    builder.add_compound(&mut compound, &selected.0);
                }
                compound
            };
            map_shapes_and_unique_ancestors_into(
                &selected_compound.0,
                ShapeType::Edge,
                ShapeType::Face,
                false,
                &mut selected_edge_faces,
            );
            let keep_edges: Vec<ShapeMapKey> = selected_edge_faces
                .iter()
                .filter_map(|(key, local)| {
                    (local.len() == 2
                        && (global_edge_faces
                            .get(key)
                            .is_some_and(|global| global.len() > 2)
                            || self.keep_shapes.contains(key)
                            || free_boundaries.contains(key)))
                    .then(|| key.clone())
                })
                .collect();
            if !keep_edges.is_empty() && !self.allow_internal {
                prune_faces_at_keep_edges(
                    &mut faces,
                    &mut boundary_shapes,
                    &mut removed_shapes,
                    &selected_edge_faces,
                    &keep_edges,
                );
                if faces.len() <= 1 {
                    continue;
                }
            }

            let builder = TopoBuilder::new();
            let mut reference_face =
                builder.make_face_with_tolerance(base_surface.clone(), &[], 0.0);
            reference_face.0.set_orientation(Orientation::Forward);
            let mut temporary_pcurves = HashSet::new();
            super::pcurves::add_pcurves(&faces, &reference_face, &mut temporary_pcurves)?;

            let mut occurrence = HashMap::new();
            for selected in &faces {
                for edge in explorer_shapes(&selected.0, ShapeType::Edge) {
                    occurrence.entry(shape_map_key(&edge)).or_insert(edge);
                }
            }
            for boundary in &mut boundary_shapes {
                let key = shape_map_key(boundary);
                if let Some(edge) = occurrence.get(&key) {
                    *boundary = edge.clone();
                }
            }
            if self.allow_internal {
                for key in keep_edges {
                    if let Some(edge) = occurrence.get(&key) {
                        let mut internal = edge.clone();
                        internal.set_orientation(Orientation::Internal);
                        boundary_shapes.push(internal);
                    }
                }
            }
            let mut internal_edges = Vec::new();
            boundary_shapes.retain(|boundary| {
                if boundary.orientation() == Orientation::Internal {
                    if let Some(edge) = Edge::wrap(boundary.clone()) {
                        internal_edges.push(edge);
                    }
                    false
                } else {
                    true
                }
            });
            if reference_orientation == Orientation::Reversed {
                for boundary in &mut boundary_shapes {
                    boundary.reverse();
                }
            }

            let edges: Vec<Edge> = boundary_shapes.into_iter().filter_map(Edge::wrap).collect();
            let boundary_edge_keys = edges.iter().map(|edge| shape_map_key(&edge.0)).collect();
            let periodic = super::periodic_faces::prepare_periodic_faces(
                &faces,
                &edges,
                base_surface,
                reference_face,
                &selected_edge_faces,
                &boundary_edge_keys,
                &temporary_pcurves,
                self.concat_bsplines,
                self.angular_tolerance,
            )?;
            base_surface = periodic.surface;
            reference_face = periodic.reference_face;
            let removed_edges: Vec<Edge> =
                removed_shapes.into_iter().filter_map(Edge::wrap).collect();
            let boundary_wires = super::boundary_rebuild::build_boundary_wires(
                &edges,
                &reference_face,
                &removed_edges,
                periodic.u_period,
                periodic.v_period,
            )?;
            let mut new_faces = Vec::new();
            let mut loose_wires = Vec::new();
            for wire in boundary_wires {
                let on_surface_bound = crate::topo_tools::edges_of_wire(&wire)
                    .iter()
                    .any(|edge| BRepTool::is_closed_edge_face(edge, &reference_face));
                if on_surface_bound {
                    let mut result =
                        builder.make_face_with_tolerance(base_surface.clone(), &[wire], 0.0);
                    result.0.set_orientation(reference_orientation);
                    new_faces.push(result);
                } else {
                    loose_wires.push(wire);
                }
            }
            let internal_wires = super::boundary_rebuild::build_internal_wires(&internal_edges);
            if new_faces.is_empty() {
                loose_wires.extend(internal_wires);
                let mut result = builder.make_face_with_tolerance(base_surface, &loose_wires, 0.0);
                result.0.set_orientation(reference_orientation);
                merge_face_result(&mut self.context, &mut self.face_new_face, &faces, &result)?;
            } else if new_faces.len() == 1 {
                for wire in loose_wires.iter().chain(internal_wires.iter()) {
                    builder.add_wire(&mut new_faces[0], wire);
                }
                merge_face_result(
                    &mut self.context,
                    &mut self.face_new_face,
                    &faces,
                    &new_faces[0],
                )?;
            } else {
                super::face_helpers::insert_wires_into_faces(
                    &loose_wires,
                    &mut new_faces,
                    &reference_face,
                )?;
                super::face_helpers::insert_wires_into_faces(
                    &internal_wires,
                    &mut new_faces,
                    &reference_face,
                )?;
                merge_split_face_results(
                    &mut self.context,
                    &mut self.face_new_face,
                    &faces,
                    &new_faces,
                    &reference_face,
                )?;
            }
        }
        Ok(())
    }
}

fn merge_face_result(
    context: &mut UnifyContext,
    face_new_face: &mut HashMap<ShapeMapKey, Face>,
    faces: &[Face],
    result: &Face,
) -> Result<(), String> {
    let parts: Vec<_> = faces.iter().map(|face| face.0.clone()).collect();
    context.merge(&parts, &result.0)?;
    for face in faces {
        face_new_face.insert(shape_map_key(&face.0), result.clone());
    }
    Ok(())
}

/// `cxx:4272-4317`: bind each seam-bounded result face to the old faces
/// contributing its ordinary, non-degenerate edges.
fn merge_split_face_results(
    context: &mut UnifyContext,
    face_new_face: &mut HashMap<ShapeMapKey, Face>,
    old_faces: &[Face],
    new_faces: &[Face],
    reference_face: &Face,
) -> Result<(), String> {
    let old_edge_sets: Vec<HashSet<ShapeMapKey>> = old_faces
        .iter()
        .map(|face| {
            explorer_shapes(&face.0, ShapeType::Edge)
                .iter()
                .map(shape_map_key)
                .collect()
        })
        .collect();
    for result in new_faces {
        let mut selected = Vec::new();
        let mut selected_keys = HashSet::new();
        for edge_shape in explorer_shapes(&result.0, ShapeType::Edge) {
            let Some(edge) = Edge::wrap(edge_shape) else {
                continue;
            };
            if BRepTool::is_degenerated(&edge)
                || BRepTool::is_closed_edge_face(&edge, reference_face)
            {
                continue;
            }
            let edge_key = shape_map_key(&edge.0);
            for (index, edge_set) in old_edge_sets.iter().enumerate() {
                if edge_set.contains(&edge_key)
                    && selected_keys.insert(shape_map_key(&old_faces[index].0))
                {
                    selected.push(old_faces[index].0.clone());
                }
            }
        }
        context.merge(&selected, &result.0)?;
        for old in old_faces {
            if selected_keys.contains(&shape_map_key(&old.0)) {
                face_new_face.insert(shape_map_key(&old.0), result.clone());
            }
        }
    }
    Ok(())
}

/// `cxx:3391-3468`: drop selected faces that cannot coexist with a keep or
/// multi-connected internal boundary, restoring their edges into the boundary.
fn prune_faces_at_keep_edges(
    faces: &mut Vec<Face>,
    boundaries: &mut Vec<TopoShape>,
    removed_edges: &mut Vec<TopoShape>,
    edge_faces: &EdgeFaceMap,
    keep_edges: &[ShapeMapKey],
) {
    let mut avoid_faces = HashSet::new();
    for edge in keep_edges {
        if let Some(adjacent) = edge_faces.get(edge) {
            if let (Some(first), Some(last)) = (adjacent.first(), adjacent.last()) {
                avoid_faces.insert(shape_map_key(first));
                avoid_faces.insert(shape_map_key(last));
            }
        }
    }

    let mut index = 0;
    while index < faces.len() {
        if !avoid_faces.contains(&shape_map_key(&faces[index].0)) {
            index += 1;
            continue;
        }
        let connected_elsewhere = explorer_shapes(&faces[index].0, ShapeType::Edge)
            .into_iter()
            .any(|edge| {
                edge_faces
                    .get(&shape_map_key(&edge))
                    .is_some_and(|adjacent| {
                        adjacent.len() > 1
                            && adjacent
                                .iter()
                                .any(|face| !avoid_faces.contains(&shape_map_key(face)))
                    })
            });
        if connected_elsewhere {
            index += 1;
            continue;
        }
        super::boundaries::add_ordinary_edges(boundaries, &faces[index].0, removed_edges);
        faces.remove(index);
    }

    let selected: HashSet<_> = faces.iter().map(|face| shape_map_key(&face.0)).collect();
    for edge in keep_edges {
        let Some(adjacent) = edge_faces.get(edge) else {
            continue;
        };
        if adjacent.len() < 2
            || !selected.contains(&shape_map_key(&adjacent[0]))
            || !selected.contains(&shape_map_key(&adjacent[1]))
        {
            continue;
        }
        let adjacent_keys = [shape_map_key(&adjacent[0]), shape_map_key(&adjacent[1])];
        let mut face_index = 0;
        while face_index < faces.len() {
            if adjacent_keys.contains(&shape_map_key(&faces[face_index].0)) {
                super::boundaries::add_ordinary_edges(
                    boundaries,
                    &faces[face_index].0,
                    removed_edges,
                );
                faces.remove(face_index);
            } else {
                face_index += 1;
            }
        }
    }
}
