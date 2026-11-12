//! `XSAlgo_ShapeProcessor` parameter-map bridge.

use std::collections::BTreeMap;

pub use occt_topo::shhealing::ShapeFixParameters;

pub type ShapeFixParameterMap = BTreeMap<String, String>;

fn set<T: ToString>(map: &mut ShapeFixParameterMap, key: &str, value: T, replace: bool) {
    let value = value.to_string();
    if let Some(existing) = map.get(key) {
        if !replace || existing == &value {
            return;
        }
    }
    map.insert(key.to_string(), value);
}

pub fn fill_parameter_map(
    parameters: &ShapeFixParameters,
    replace: bool,
    map: &mut ShapeFixParameterMap,
) {
    macro_rules! p {
        ($key:literal, $field:ident) => {
            set(map, concat!("FixShape.", $key), parameters.$field, replace)
        };
    }
    p!("Tolerance3d", tolerance_3d);
    p!("MaxTolerance3d", max_tolerance_3d);
    p!("MinTolerance3d", min_tolerance_3d);
    set(
        map,
        "DetalizationLevel",
        parameters.detalization_level as i32,
        replace,
    );
    set(
        map,
        "NonManifold",
        i32::from(parameters.non_manifold),
        replace,
    );
    p!("FixFreeShellMode", fix_free_shell_mode);
    p!("FixFreeFaceMode", fix_free_face_mode);
    p!("FixFreeWireMode", fix_free_wire_mode);
    p!("FixSameParameterMode", fix_same_parameter_mode);
    p!("FixSolidMode", fix_solid_mode);
    p!("FixShellOrientationMode", fix_shell_orientation_mode);
    p!("CreateOpenSolidMode", create_open_solid_mode);
    p!("FixShellMode", fix_shell_mode);
    p!("FixFaceOrientationMode", fix_face_orientation_mode);
    p!("FixFaceMode", fix_face_mode);
    p!("FixWireMode", fix_wire_mode);
    p!("FixOrientationMode", fix_orientation_mode);
    p!("FixAddNaturalBoundMode", fix_add_natural_bound_mode);
    p!("FixMissingSeamMode", fix_missing_seam_mode);
    p!("FixSmallAreaWireMode", fix_small_area_wire_mode);
    p!("RemoveSmallAreaFaceMode", remove_small_area_face_mode);
    p!("FixIntersectingWiresMode", fix_intersecting_wires_mode);
    p!("FixLoopWiresMode", fix_loop_wires_mode);
    p!("FixSplitFaceMode", fix_split_face_mode);
    p!("AutoCorrectPrecisionMode", auto_correct_precision_mode);
    p!("ModifyTopologyMode", modify_topology_mode);
    p!("ModifyGeometryMode", modify_geometry_mode);
    p!("ClosedWireMode", closed_wire_mode);
    p!("PreferencePCurveMode", preference_pcurve_mode);
    p!("FixReorderMode", fix_reorder_mode);
    p!("FixSmallMode", fix_small_mode);
    p!("FixConnectedMode", fix_connected_mode);
    p!("FixEdgeCurvesMode", fix_edge_curves_mode);
    p!("FixDegeneratedMode", fix_degenerated_mode);
    p!("FixLackingMode", fix_lacking_mode);
    p!("FixSelfIntersectionMode", fix_self_intersection_mode);
    p!("RemoveLoopMode", remove_loop_mode);
    p!("FixReversed2dMode", fix_reversed_2d_mode);
    p!("FixRemovePCurveMode", fix_remove_pcurve_mode);
    p!("FixRemoveCurve3dMode", fix_remove_curve_3d_mode);
    p!("FixAddPCurveMode", fix_add_pcurve_mode);
    p!("FixAddCurve3dMode", fix_add_curve_3d_mode);
    p!("FixSeamMode", fix_seam_mode);
    p!("FixShiftedMode", fix_shifted_mode);
    p!("FixEdgeSameParameterMode", fix_edge_same_parameter_mode);
    p!("FixNotchedEdgesMode", fix_notched_edges_mode);
    p!("FixTailMode", fix_tail_mode);
    p!("MaxTailAngle", max_tail_angle);
    p!("MaxTailWidth", max_tail_width);
    p!(
        "FixSelfIntersectingEdgeMode",
        fix_self_intersecting_edge_mode
    );
    p!("FixIntersectingEdgesMode", fix_intersecting_edges_mode);
    p!(
        "FixNonAdjacentIntersectingEdgesMode",
        fix_non_adjacent_intersecting_edges_mode
    );
    p!("FixVertexPositionMode", fix_vertex_position_mode);
    p!("FixVertexToleranceMode", fix_vertex_tolerance_mode);
}

pub fn set_shape_fix_parameters(
    parameters: &ShapeFixParameters,
    additional: &ShapeFixParameterMap,
    target: &mut ShapeFixParameterMap,
) {
    fill_parameter_map(parameters, true, target);
    for (key, value) in additional {
        set(target, key, value, false);
    }
}
