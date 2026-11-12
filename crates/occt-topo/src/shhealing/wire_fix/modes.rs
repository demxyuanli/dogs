//! `ShapeFix_Wire` mode fields and `ShapeFix_Root::NeedFix`.
//!
//! Source: `ShapeFix_Wire.hxx:471-499` (fields), `ShapeFix_Wire.cxx:167-201`
//! (`ClearModes` defaults), `ShapeFix_Root.lxx:101-104` (`NeedFix`).
//!
//! The reader-time passes used to carry their gates as ad-hoc booleans
//! (`reorder_ok && ...`, the `FIXEDGE_SAMEPARAM_ARM` const). This struct is the
//! same information in the shape OCCT keeps it, so the gates read
//! `need_fix(<mode>, <need>)` exactly as the `.cxx` does and the per-round
//! overrides (`ShapeFix_Face::Perform`, `ShapeFix_Shape::Perform`,
//! `ShapeProcess_OperLibrary::fixshape`) become presets instead of edits at
//! every call site.

/// `ShapeFix_Root::NeedFix` (`ShapeFix_Root.lxx:101-104`):
/// `Flag < 0 ? need : (Flag > 0)`.
pub fn need_fix(flag: i32, need: bool) -> bool {
    if flag < 0 {
        need
    } else {
        flag > 0
    }
}

/// `DE_ShapeFixParameters::FixMode` (`DE_ShapeFixParameters.hxx:24-29`).
pub const FIX_OR_NOT: i32 = -1;
pub const NOT_FIX: i32 = 0;
pub const FIX: i32 = 1;

/// The `ShapeFix_Wire` mode set (`ShapeFix_Wire.hxx:471-499`) plus the two
/// tail limits `ShapeFix_Wire::SetMaxTailAngle` / `SetMaxTailWidth` keep
/// (`ShapeFix_Wire.cxx:155-163`, `:4317-4324`).
#[derive(Debug, Clone, PartialEq)]
pub struct WireFixModes {
    pub geom_mode: bool,
    pub topo_mode: bool,
    pub closed_mode: bool,
    pub preference_2d: bool,
    pub fix_gaps_by_ranges: bool,

    pub remove_loop_mode: i32,

    pub fix_reversed2d_mode: i32,
    pub fix_remove_pcurve_mode: i32,
    pub fix_remove_curve3d_mode: i32,
    pub fix_add_pcurve_mode: i32,
    pub fix_add_curve3d_mode: i32,
    pub fix_seam_mode: i32,
    pub fix_shifted_mode: i32,
    pub fix_same_parameter_mode: i32,
    pub fix_vertex_tolerance_mode: i32,

    pub fix_notched_edges_mode: i32,
    pub fix_self_intersecting_edge_mode: i32,
    pub fix_intersecting_edges_mode: i32,
    pub fix_non_adjacent_intersecting_edges_mode: i32,
    pub fix_tail_mode: i32,

    pub fix_reorder_mode: i32,
    pub fix_small_mode: i32,
    pub fix_connected_mode: i32,
    pub fix_edge_curves_mode: i32,
    pub fix_degenerated_mode: i32,
    pub fix_self_intersection_mode: i32,
    pub fix_lacking_mode: i32,
    pub fix_gaps3d_mode: i32,
    pub fix_gaps2d_mode: i32,

    /// `myMaxTailAngleSine` (`ShapeFix_Wire.cxx:4324` stores `Sin(angle)`).
    pub max_tail_angle_sine: f64,
    /// `myMaxTailWidth`.
    pub max_tail_width: f64,
}

impl Default for WireFixModes {
    /// `ShapeFix_Wire::ClearModes` (`ShapeFix_Wire.cxx:167-201`).
    fn default() -> Self {
        WireFixModes {
            geom_mode: true,
            topo_mode: false,
            closed_mode: true,
            preference_2d: true,
            fix_gaps_by_ranges: false,
            remove_loop_mode: FIX_OR_NOT,
            fix_reversed2d_mode: FIX_OR_NOT,
            fix_remove_pcurve_mode: FIX_OR_NOT,
            fix_remove_curve3d_mode: FIX_OR_NOT,
            fix_add_pcurve_mode: FIX_OR_NOT,
            fix_add_curve3d_mode: FIX_OR_NOT,
            fix_seam_mode: FIX_OR_NOT,
            fix_shifted_mode: FIX_OR_NOT,
            fix_same_parameter_mode: FIX_OR_NOT,
            fix_vertex_tolerance_mode: FIX_OR_NOT,
            fix_notched_edges_mode: FIX_OR_NOT,
            fix_self_intersecting_edge_mode: FIX_OR_NOT,
            fix_intersecting_edges_mode: FIX_OR_NOT,
            fix_non_adjacent_intersecting_edges_mode: FIX_OR_NOT,
            fix_tail_mode: NOT_FIX,
            fix_reorder_mode: FIX_OR_NOT,
            fix_small_mode: FIX_OR_NOT,
            fix_connected_mode: FIX_OR_NOT,
            fix_edge_curves_mode: FIX_OR_NOT,
            fix_degenerated_mode: FIX_OR_NOT,
            fix_self_intersection_mode: FIX_OR_NOT,
            fix_lacking_mode: FIX_OR_NOT,
            fix_gaps3d_mode: FIX_OR_NOT,
            fix_gaps2d_mode: FIX_OR_NOT,
            max_tail_angle_sine: 0.0,
            max_tail_width: -1.0,
        }
    }
}

impl WireFixModes {
    pub fn from_parameters(parameters: &crate::shhealing::ShapeFixParameters) -> Self {
        let mut modes = Self::default();
        modes.geom_mode = parameters.modify_geometry_mode > 0;
        modes.topo_mode = parameters.modify_topology_mode > 0;
        modes.closed_mode = parameters.closed_wire_mode > 0;
        modes.preference_2d = parameters.preference_pcurve_mode > 0;
        modes.remove_loop_mode = parameters.remove_loop_mode;
        modes.fix_reversed2d_mode = parameters.fix_reversed_2d_mode;
        modes.fix_add_pcurve_mode = parameters.fix_add_pcurve_mode;
        modes.fix_add_curve3d_mode = parameters.fix_add_curve_3d_mode;
        modes.fix_shifted_mode = parameters.fix_shifted_mode;
        modes.fix_same_parameter_mode = parameters.fix_edge_same_parameter_mode;
        modes.fix_notched_edges_mode = parameters.fix_notched_edges_mode;
        modes.fix_self_intersecting_edge_mode = parameters.fix_self_intersecting_edge_mode;
        modes.fix_intersecting_edges_mode = parameters.fix_intersecting_edges_mode;
        modes.fix_non_adjacent_intersecting_edges_mode =
            parameters.fix_non_adjacent_intersecting_edges_mode;
        modes.fix_reorder_mode = parameters.fix_reorder_mode;
        modes.fix_small_mode = parameters.fix_small_mode;
        modes.fix_connected_mode = parameters.fix_connected_mode;
        modes.fix_edge_curves_mode = parameters.fix_edge_curves_mode;
        modes.fix_degenerated_mode = parameters.fix_degenerated_mode;
        modes.fix_self_intersection_mode = parameters.fix_self_intersection_mode;
        modes.fix_lacking_mode = parameters.fix_lacking_mode;
        // Not injected because the underlying fixers are absent:
        // FixRemovePCurve (ShapeFix_Wire.cxx:637), FixRemoveCurve3d
        // (cxx:866), FixSeam (cxx:1615-1670), and FixTails
        // (cxx:4316-4425). Their DE defaults keep all four paths inert.
        modes
    }

    /// `ClearModes` defaults (`ShapeFix_Wire.cxx:167-201`), i.e. a naked
    /// `ShapeFix_Wire` before any `ShapeProcess` operator touches it.
    pub fn naked() -> Self {
        Self::default()
    }

    /// `SetMaxTailAngle` (`ShapeFix_Wire.cxx:4317-4324`): the field keeps the
    /// sine of the degree value.
    pub fn set_max_tail_angle_deg(&mut self, deg: f64) {
        self.max_tail_angle_sine = deg.to_radians().sin();
    }

    /// The `ShapeFix_Wire` that `ShapeProcess_OperLibrary::fixshape`
    /// (`cxx:796-866`) hands to `ShapeFix_Wire::Perform` for a **face** on the
    /// `FromSTEP.FixShape` path.
    ///
    /// Values come from `STEPControl_Controller.cxx:229-262` (legacy map) and
    /// `DESTEP_Parameters::GetDefaultShapeFixParameters`
    /// (`DESTEP_Parameters.cxx:90-93`) -> `DE_ShapeFixParameters.hxx:31-85`
    /// (DE map); both agree on every field this struct carries. The two
    /// non-`ShapeProcess` overrides are
    /// `ShapeFix_Shape.cxx:199-200` (`ModifyTopologyMode() = true` before the
    /// face round) and `ShapeFix_Face.cxx:171-177` (the face tool forwards its
    /// precision / tolerances to the wire tool).
    pub fn step_face() -> Self {
        let mut m = Self::default();
        // `ShapeProcess_OperLibrary.cxx:849-853`: booleans.
        m.topo_mode = true; // `ShapeFix_Shape.cxx:200`
        m.geom_mode = true; // `ModifyGeometryMode` = `Fix`
        m.closed_mode = true; // `ClosedWireMode` = `Fix`
        m.preference_2d = true; // `PreferencePCurveMode` = `Fix`
                                // `ShapeProcess_OperLibrary.cxx:854-866`: integer modes.
                                // `DE_ShapeFixParameters.hxx:69-73` `FixMode::NotFix` for the three
                                // same-parameter / seam / tail keys below.
        m.fix_seam_mode = FIX_OR_NOT;
        m.fix_same_parameter_mode = NOT_FIX; // `FixEdgeSameParameterMode`, hxx:75
        m.fix_tail_mode = NOT_FIX; // hxx:77
                                   // `DE_ShapeFixParameters.hxx:29-33` / `STEPControl_Controller.cxx:262-263`.
        m.max_tail_angle_sine = 0.0_f64.to_radians().sin();
        m.max_tail_width = -1.0;
        m
    }

    /// The `ShapeFix_Wire` for a **face** on the `FromIGES.FixShape` path.
    ///
    /// `IGESData.cxx:136-166` (legacy map) and
    /// `DEIGES_Parameters::GetDefaultShapeFixParameters`
    /// (`DEIGES_Parameters.cxx:64-72`) agree on every wire key except
    /// `FixTailMode`, which the DE map sets to `FixOrNot` (`:69`) while the
    /// legacy default is `0`. `ShapeFix_Wire::FixTails` returns false either
    /// way because `MaxTailWidth` is `-1` (`ShapeFix_Wire.cxx:4317-4320`), so
    /// the distinction is not observable.
    pub fn iges_face() -> Self {
        let mut m = Self::default();
        m.topo_mode = true; // `ShapeFix_Shape.cxx:200`
        m.fix_same_parameter_mode = NOT_FIX; // `FixEdgeSameParameterMode`, IGESData.cxx:155
        m.fix_tail_mode = FIX_OR_NOT; // DEIGES_Parameters.cxx:69
        m.max_tail_angle_sine = 0.0_f64.to_radians().sin();
        m.max_tail_width = -1.0;
        m
    }

    /// `ShapeFix_Face::Perform`'s **first** wire round
    /// (`ShapeFix_Face.cxx:365-376`): `FixLacking` and `FixSelfIntersection`
    /// are saved and switched off, `FixNotchedEdges` is deliberately left
    /// alone (CR0024983, `cxx:370`). Both are restored at `cxx:457-458`
    /// before the second round.
    pub fn first_wire_round(&self) -> Self {
        let mut m = self.clone();
        m.fix_lacking_mode = NOT_FIX;
        m.fix_self_intersection_mode = NOT_FIX;
        m
    }

    /// `ShapeFix_Face::Perform`'s **second** wire round
    /// (`ShapeFix_Face.cxx:513-520`): `FixSmall`, `FixConnected`,
    /// `FixEdgeCurves` and `FixDegenerated` are switched off; the `Lacking`
    /// and `SelfIntersection` modes are the restored originals.
    pub fn second_wire_round(&self) -> Self {
        let mut m = self.clone();
        m.fix_small_mode = NOT_FIX;
        m.fix_connected_mode = NOT_FIX;
        m.fix_edge_curves_mode = NOT_FIX;
        m.fix_degenerated_mode = NOT_FIX;
        m
    }
}
