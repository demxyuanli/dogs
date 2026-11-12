//! `DE_ShapeFixParameters` (`DE_ShapeFixParameters.hxx:20-85`).

use crate::abs::ShapeType;

use super::{FIX, FIX_OR_NOT, NOT_FIX};

#[derive(Debug, Clone, PartialEq)]
pub struct ShapeFixParameters {
    pub tolerance_3d: f64,
    pub max_tolerance_3d: f64,
    pub min_tolerance_3d: f64,
    pub detalization_level: ShapeType,
    pub non_manifold: bool,
    pub fix_free_shell_mode: i32,
    pub fix_free_face_mode: i32,
    pub fix_free_wire_mode: i32,
    pub fix_same_parameter_mode: i32,
    pub fix_solid_mode: i32,
    pub fix_shell_orientation_mode: i32,
    pub create_open_solid_mode: i32,
    pub fix_shell_mode: i32,
    pub fix_face_orientation_mode: i32,
    pub fix_face_mode: i32,
    pub fix_wire_mode: i32,
    pub fix_orientation_mode: i32,
    pub fix_add_natural_bound_mode: i32,
    pub fix_missing_seam_mode: i32,
    pub fix_small_area_wire_mode: i32,
    pub remove_small_area_face_mode: i32,
    pub fix_intersecting_wires_mode: i32,
    pub fix_loop_wires_mode: i32,
    pub fix_split_face_mode: i32,
    pub auto_correct_precision_mode: i32,
    pub modify_topology_mode: i32,
    pub modify_geometry_mode: i32,
    pub closed_wire_mode: i32,
    pub preference_pcurve_mode: i32,
    pub fix_reorder_mode: i32,
    pub fix_small_mode: i32,
    pub fix_connected_mode: i32,
    pub fix_edge_curves_mode: i32,
    pub fix_degenerated_mode: i32,
    pub fix_lacking_mode: i32,
    pub fix_self_intersection_mode: i32,
    pub remove_loop_mode: i32,
    pub fix_reversed_2d_mode: i32,
    pub fix_remove_pcurve_mode: i32,
    pub fix_remove_curve_3d_mode: i32,
    pub fix_add_pcurve_mode: i32,
    pub fix_add_curve_3d_mode: i32,
    pub fix_seam_mode: i32,
    pub fix_shifted_mode: i32,
    pub fix_edge_same_parameter_mode: i32,
    pub fix_notched_edges_mode: i32,
    pub fix_tail_mode: i32,
    pub max_tail_angle: i32,
    pub max_tail_width: i32,
    pub fix_self_intersecting_edge_mode: i32,
    pub fix_intersecting_edges_mode: i32,
    pub fix_non_adjacent_intersecting_edges_mode: i32,
    pub fix_vertex_position_mode: i32,
    pub fix_vertex_tolerance_mode: i32,
}

impl Default for ShapeFixParameters {
    fn default() -> Self {
        Self {
            tolerance_3d: 1.0e-6,
            max_tolerance_3d: 1.0,
            min_tolerance_3d: 1.0e-7,
            detalization_level: ShapeType::Vertex,
            non_manifold: false,
            fix_free_shell_mode: FIX_OR_NOT,
            fix_free_face_mode: FIX_OR_NOT,
            fix_free_wire_mode: FIX_OR_NOT,
            fix_same_parameter_mode: FIX_OR_NOT,
            fix_solid_mode: FIX_OR_NOT,
            fix_shell_orientation_mode: FIX_OR_NOT,
            create_open_solid_mode: NOT_FIX,
            fix_shell_mode: FIX_OR_NOT,
            fix_face_orientation_mode: FIX_OR_NOT,
            fix_face_mode: FIX_OR_NOT,
            fix_wire_mode: FIX_OR_NOT,
            fix_orientation_mode: FIX_OR_NOT,
            fix_add_natural_bound_mode: FIX_OR_NOT,
            fix_missing_seam_mode: FIX_OR_NOT,
            fix_small_area_wire_mode: FIX_OR_NOT,
            remove_small_area_face_mode: FIX_OR_NOT,
            fix_intersecting_wires_mode: FIX_OR_NOT,
            fix_loop_wires_mode: FIX_OR_NOT,
            fix_split_face_mode: FIX_OR_NOT,
            auto_correct_precision_mode: FIX,
            modify_topology_mode: NOT_FIX,
            modify_geometry_mode: FIX,
            closed_wire_mode: FIX,
            preference_pcurve_mode: FIX,
            fix_reorder_mode: FIX_OR_NOT,
            fix_small_mode: FIX_OR_NOT,
            fix_connected_mode: FIX_OR_NOT,
            fix_edge_curves_mode: FIX_OR_NOT,
            fix_degenerated_mode: FIX_OR_NOT,
            fix_lacking_mode: FIX_OR_NOT,
            fix_self_intersection_mode: FIX_OR_NOT,
            remove_loop_mode: FIX_OR_NOT,
            fix_reversed_2d_mode: FIX_OR_NOT,
            fix_remove_pcurve_mode: FIX_OR_NOT,
            fix_remove_curve_3d_mode: FIX_OR_NOT,
            fix_add_pcurve_mode: FIX_OR_NOT,
            fix_add_curve_3d_mode: FIX_OR_NOT,
            fix_seam_mode: FIX_OR_NOT,
            fix_shifted_mode: FIX_OR_NOT,
            fix_edge_same_parameter_mode: NOT_FIX,
            fix_notched_edges_mode: FIX_OR_NOT,
            fix_tail_mode: NOT_FIX,
            max_tail_angle: NOT_FIX,
            max_tail_width: FIX_OR_NOT,
            fix_self_intersecting_edge_mode: FIX_OR_NOT,
            fix_intersecting_edges_mode: FIX_OR_NOT,
            fix_non_adjacent_intersecting_edges_mode: FIX_OR_NOT,
            fix_vertex_position_mode: NOT_FIX,
            fix_vertex_tolerance_mode: FIX_OR_NOT,
        }
    }
}

impl ShapeFixParameters {
    pub fn step_reader() -> Self {
        Self::default()
    }

    pub fn iges_reader() -> Self {
        Self {
            detalization_level: ShapeType::Edge,
            create_open_solid_mode: FIX,
            fix_tail_mode: FIX_OR_NOT,
            max_tail_angle: FIX_OR_NOT,
            ..Self::default()
        }
    }
}
