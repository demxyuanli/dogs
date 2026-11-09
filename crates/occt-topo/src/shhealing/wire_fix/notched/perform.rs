use super::super::*;
use super::statics::*;


/// First wire round of `ShapeFix_Face::Perform` (`ShapeFix_Face.cxx:365-480`).
///
/// `FixLackingMode` and `FixSelfIntersectionMode` are forced off for this
/// round (`cxx:369-374`) and restored only for the post-seam round.
///
/// `FixTails` (`ShapeFix_Wire.cxx:414`) makes no change on this path. The gate
/// is `myFixTailMode != 0`; the IGES DE default is
/// `DEIGES_Parameters.cxx:69` `FixTailMode = FixMode::FixOrNot = -1`, so the
/// call IS made, but `FixTails` opens with
/// `if (myMaxTailWidth < 0 || !IsReady()) return false;`
/// (`ShapeFix_Wire.cxx:4317-4320`) and `myMaxTailWidth` is
/// `ctx->RealVal("MaxTailWidth", -1)` (`ShapeProcess_OperLibrary.cxx:861`) from
/// `DE_ShapeFixParameters.hxx:79` `MaxTailWidth = FixMode::FixOrNot = -1`. The
/// STEP DE default is `FixTailMode = FixMode::NotFix = 0`
/// (`DE_ShapeFixParameters.hxx:77`, `DESTEP_Parameters.cxx:92-95`), so STEP
/// does not even call. Either way the pass is behaviourally absent, which is
/// why `FixTails` is intentionally not ported.
///
/// The singularity sub-branch (`ShapeFix_Wire.cxx:673-863`) runs here, after
/// `FixAddPCurve` reports `DONE2`. The STEP path `check_pcurves_and_shift`
/// does not call it.
pub fn trace_wire(tag: &str, wire: &Wire) {
    let edges = edges_of_wire(wire);
    eprintln!("TRACEWIRE {tag} n={}", edges.len());
    for (i, e) in edges.iter().enumerate() {
        let fmt = |v: &Option<Vertex>| match v {
            Some(v) => {
                let p = BRepTool::vertex_point(v);
                format!(
                    "({:.9},{:.9},{:.9})@{:#x}",
                    p.x(),
                    p.y(),
                    p.z(),
                    std::sync::Arc::as_ptr(&v.0.tshape) as usize
                )
            }
            None => "None".to_string(),
        };
        eprintln!(
            "TRACEWIRE   e{} first={} last={}",
            i + 1,
            fmt(&first_vertex(e)),
            fmt(&last_vertex(e))
        );
    }
}


pub fn perform_wire_first_round(
    wire: &mut Wire,
    face: &Face,
    preci: f64,
    ctx: &mut crate::shape_fix_compose_shell::SharedReShape,
    modes: &WireFixModes,
) {
    let mut reorder_ok = fix_reorder_wire(wire, face);
    let reorder_before_small = reorder_ok;
    if need_fix(modes.fix_small_mode, modes.topo_mode)
        && fix_small_all(
            wire,
            face,
            preci,
            SHAPE_FIX_MIN_TOLERANCE,
            !modes.topo_mode || !reorder_ok,
        )
    {
        if need_fix(modes.fix_reorder_mode, !reorder_before_small) {
            reorder_ok = fix_reorder_wire(wire, face);
        }
    }
    if need_fix(modes.fix_connected_mode, reorder_ok) {
        if crate::debug_flags::iges_trace_wire() {
            trace_wire("before-connected", wire);
        }
        // `ShapeFix_Wire::FixConnected(prec)` (`ShapeFix_Wire.cxx:557-596`)
        // runs with this tool's `Context()`, which the face round inherited
        // from `ShapeFix_Shape` (`ShapeFix_Face.cxx:379`).
        let changed = fix_connected_all_in_context(
            wire,
            Some(ctx),
            preci,
            SHAPE_FIX_MAX_TOLERANCE,
            modes.closed_mode,
        );
        if crate::debug_flags::iges_trace_wire() {
            eprintln!("TRACEWIRE connected-changed={changed}");
            trace_wire("after-connected", wire);
        }
    }
    let mut cache = crate::pcurve_full::ProjectorCache::new();
    let mut overdegen: Option<usize> = None;
    let mut i = 0usize;
    while i < edges_of_wire(wire).len() {
        let edges = edges_of_wire(wire);
        let edge = edges[i].clone();
        let is_seam = is_seam_use(&edges, i);
        cache.done4 = false;
        let _ = fix_add_pcurve(&edge, face, is_seam, preci, &mut cache);
        if !is_seam && cache.done4 {
            if super::super::singularity::split_edge_over_singularity(wire, face, &edge, i, ctx) {
                continue;
            }
            overdegen = Some(i);
        }
        i += 1;
    }
    if let Some(k) = overdegen {
        super::super::singularity::adjust_overdegen_pcurve(wire, face, k, preci, &mut cache);
    }
    // `ShapeFix_Wire::FixEdgeCurves` (`ShapeFix_Wire.cxx:881-918`):
    // `NeedFix(myFixAddCurve3dMode)` is true at the tool default `-1`
    // (`cxx:181`), so every edge without a 3D curve gets one built from its
    // pcurve before FixSeam / FixShifted / FixSameParameter run.
    if need_fix(modes.fix_edge_curves_mode, true)
        && need_fix(modes.fix_add_curve3d_mode, true)
    {
        fix_add_curve3d_or_drop(wire, face, preci);
    }
    // `ShapeFix_Wire::FixEdgeCurves` (`ShapeFix_Wire.cxx:920-933`): `FixSeam` is
    // NOT PORTED and is NOT REACHED. The gate is
    // `isReady && NeedFix(myFixSeamMode, false)` and `NeedFix(Flag, need)` is
    // `Flag < 0 ? need : Flag > 0` (`ShapeFix_Root.lxx:101-104`), so the
    // `false` argument is what a default `-1` mode resolves to - the arm is off
    // unless the mode is explicitly positive. `myFixSeamMode` is `-1`
    // (`ShapeFix_Wire.cxx:182`), the IGES operator leaves it at `-1`
    // (`FromIGES.FixShape.FixSeamMode`, `IGESData.cxx:150`) and so does the DE
    // map (`DE_ShapeFixParameters.hxx:73` `FixMode::FixOrNot = -1`,
    // `XSAlgo_ShapeProcessor.cxx:679`), and `ShapeProcess_OperLibrary.cxx:856`
    // passes that through. Both read paths therefore skip the arm, and the
    // `DONE6` / `FAIL6` accumulation with it. Not porting `FixSeam`
    // (`ShapeFix_Wire.cxx:1615-1670`) is the faithful result.
    //
    // The two arms in the same block that are off by the same rule are
    // `FixRemovePCurve` (`cxx:637`, `NeedFix(myFixRemovePCurveMode, false)`,
    // modes `-1` at `IGESData.cxx:146` / `DE_ShapeFixParameters.hxx:69`) and
    // `FixRemoveCurve3d` (`cxx:866`, `IGESData.cxx:147` / `hxx:70`), which are
    // likewise not ported and not reached.
    if need_fix(modes.fix_edge_curves_mode, true) {
        let shifted_mode = if modes.fix_shifted_mode == FIX_OR_NOT && !reorder_ok {
            NOT_FIX
        } else {
            modes.fix_shifted_mode
        };
        if need_fix(shifted_mode, true) {
            let _ = fix_shifted_wire(wire, face);
        }
    }
    // `ShapeFix_Wire::FixEdgeCurves` (`ShapeFix_Wire.cxx:953-1005`): the gate is
    // `isReady && NeedFix(myFixSameParameterMode)`, and this wire tool's mode is
    // the *separate* key `FixEdgeSameParameterMode`
    // (`ShapeProcess_OperLibrary.cxx:857`), whose IGES value is `0`
    // (`DE_ShapeFixParameters.hxx:75` `FixMode::NotFix`, `DEIGES_Parameters.cxx:64-72`,
    // `XSAlgo_ShapeProcessor.cxx:681`, and the legacy
    // `IGESData.cxx:152`), so `NeedFix(0)` is false and OCCT runs neither the
    // `CheckPCurveRange` re-projection below nor `FixSameParameter(Edge, Face())`.
    // The port still runs both, before `FixDegenerated` / `FixSelfIntersection`
    // / `FixLacking`. The whole-shape equivalent OCCT does run is
    // `ShapeFix_Shape.cxx:259` -> `ShapeFix::SameParameter(Result, false, 0.0)`
    // (`ShapeFix.cxx:70-258`), which the port currently covers by calling
    // `BRepLib::SameParameter` (`iges_process_shape.rs`, `same_parameter_shape`)
    // instead of that function. Bringing the two into line is tracked as the
    // next worker package; the oracle meshes are the check on it.
    let edges = edges_of_wire(wire);
    let reg = GeometryRegistry::global();
    let same_parameter_arm =
        need_fix(modes.fix_edge_curves_mode, true) && need_fix(modes.fix_same_parameter_mode, true);
    for (i, edge) in edges.iter().enumerate() {
        if same_parameter_arm {
            if let Some((c2d, fp2d, lp2d)) = curve_on_surface_oriented(edge, face, false) {
                let (first, last) = BRepTool::edge_parameters(edge);
                if first.is_finite() && last.is_finite() {
                    if (first - fp2d).abs() > PCONFUSION || (last - lp2d).abs() > PCONFUSION {
                        reg.set_same_range(&edge.0, false);
                    } else if !check_pcurve_range(first, last, c2d.as_ref()) {
                        let is_seam = is_seam_use(&edges, i);
                        reg.remove_pcurves_on_surface(&edge.0, &face.0);
                        let _ = fix_add_pcurve(edge, face, is_seam, preci, &mut cache);
                    }
                }
            }
            fix_same_parameter(edge, face);
        }
    }
    // `ShapeFix_Wire::FixEdgeCurves` (`ShapeFix_Wire.cxx:1009-1022`):
    // `NeedFix(myFixVertexToleranceMode)`.
    if need_fix(modes.fix_vertex_tolerance_mode, true) {
        for edge in edges_of_wire(wire) {
            fix_vertex_tolerance_edge_on_face(&edge, face);
        }
    }
    // `ShapeFix_Wire::Perform` (`cxx:386-392`): FixDegenerated.
    if need_fix(modes.fix_degenerated_mode, true) {
        fix_degenerated_all(wire, face, preci);
    }
    // `ShapeFix_Wire::Perform` (`cxx:400-406`): FixNotchedEdges, then FixShifted.
    if modes.fix_tail_mode <= 0
        && need_fix(modes.fix_notched_edges_mode, reorder_ok)
        && fix_notched_edges(wire, face, SHAPE_FIX_MIN_TOLERANCE, SHAPE_FIX_MAX_TOLERANCE)
    {
        let _ = fix_shifted_wire(wire, face);
    }
    // `ShapeFix_Wire::Perform` (`ShapeFix_Wire.cxx:462-470`): face overload,
    // not `CheckVertexTolerance(edge)` with `checkAll=true`.
    if need_fix(modes.fix_vertex_tolerance_mode, true) {
        for edge in edges_of_wire(wire) {
            fix_vertex_tolerance_edge_on_face(&edge, face);
        }
    }
    if crate::debug_flags::iges_trace_wire() {
        trace_wire("end-first", wire);
    }
}


/// Post-seam wire round of `ShapeFix_Face::Perform` (`ShapeFix_Face.cxx:508-627`).
///
/// `cxx:513-520` forces `FixSmall`, `FixConnected`, `FixEdgeCurves`, and
/// `FixDegenerated` off, then calls `ShapeFix_Wire::Perform`. With those
/// modes off, `Perform` (`ShapeFix_Wire.cxx:297-482`) still runs FixReorder,
/// FixNotchedEdges (`FixTailMode` is 0, `cxx:400`), FixSelfIntersection
/// (`myClosedMode` is true, `cxx:171` and `cxx:428`) followed by another
/// `FixReorder` (`cxx:441`), FixLacking (`cxx:448`), and
/// `FixVertexTolerance(edge, face)` (`cxx:462-470`).
/// `FixTails` stays off. The returned flag is `StatusRemovedSegment`
/// (`ShapeFix_Face.cxx:576-578`).
pub fn perform_wire_second_round(
    wire: &mut Wire,
    face: &Face,
    preci: f64,
    modes: &WireFixModes,
) -> bool {
    let reorder_ok = fix_reorder_wire(wire, face);
    if modes.fix_tail_mode <= 0
        && need_fix(modes.fix_notched_edges_mode, reorder_ok)
        && fix_notched_edges(wire, face, SHAPE_FIX_MIN_TOLERANCE, SHAPE_FIX_MAX_TOLERANCE)
    {
        let _ = fix_shifted_wire(wire, face);
    }
    // `cxx:428-443`: FixSelfIntersection under
    // `NeedFix(myFixSelfIntersectionMode, myClosedMode)`, with the
    // adjacent-edge arm forced off while `ReorderOK` is false.
    let mut inter = FixSelfIntersectionResult {
        done: false,
        removed_segment: false,
    };
    if need_fix(modes.fix_self_intersection_mode, modes.closed_mode) {
        let intersecting_edges = if modes.fix_intersecting_edges_mode == FIX_OR_NOT && !reorder_ok {
            NOT_FIX
        } else {
            modes.fix_intersecting_edges_mode
        };
        inter = fix_self_intersection(wire, face, preci, need_fix(intersecting_edges, true));
    }
    // `cxx:441`: `FixReorder()` runs again before `FixLacking`.
    let _ = fix_reorder_wire(wire, face);
    if need_fix(modes.fix_lacking_mode, reorder_ok) {
        let _ = fix_lacking_all(wire, face, false, preci, SHAPE_FIX_MAX_TOLERANCE);
    }
    if need_fix(modes.fix_vertex_tolerance_mode, true) {
        for edge in edges_of_wire(wire) {
            fix_vertex_tolerance_edge_on_face(&edge, face);
        }
    }
    if crate::debug_flags::iges_trace_wire() {
        trace_wire("end-second", wire);
    }
    inter.removed_segment
}


/// `ShapeFix_Wire::Perform` with `ClearModes` (`ShapeFix_Wire.cxx:167-201`,
/// `:297-482`). This is the naked constructor used by
/// `IGESControl_IGESBoundary::Transfer` (`cxx:484-488`), not the face round
/// that forces `FixSmall` on and `FixLacking` / `FixSelfIntersection` off.
///
/// `myTopoMode` is false, so `FixSmall` does not run. `myFixTailMode` is 0,
/// so `FixTails` does not run. `Context()` is null, so the singularity split
/// inside `FixEdgeCurves` (`cxx:676`) does not run. `precision` is both
/// `SetPrecision` and `SetMaxTolerance` (`cxx:137-139`).
pub fn perform_wire_constructor_defaults(wire: &mut Wire, face: &Face, precision: f64) {
    // `ShapeFix_Wire::ShapeFix_Wire()` -> `ClearModes()`
    // (`ShapeFix_Wire.cxx:167-201`): every mode at `-1` / `0`.
    let modes = WireFixModes::default();
    let mut fixed = false;
    let reorder_ok = fix_reorder_wire(wire, face);
    if need_fix(modes.fix_connected_mode, reorder_ok) {
        fixed |= fix_connected_all(wire, precision, precision);
    }
    let mut cache = crate::pcurve_full::ProjectorCache::new();
    let edges = edges_of_wire(wire);
    if need_fix(modes.fix_edge_curves_mode, true) && need_fix(modes.fix_add_pcurve_mode, true) {
        for (i, edge) in edges.iter().enumerate() {
            let is_seam = is_seam_use(&edges, i);
            let _ = fix_add_pcurve(edge, face, is_seam, precision, &mut cache);
        }
    }
    if reorder_ok {
        // `ShapeFix_Wire.cxx:939-951`: `FixEdgeCurves` runs `FixShifted()` and
        // folds `LastFixStatus(DONE)` into `myStatusEdgeCurves` as `DONE7`;
        // `FixEdgeCurves` then returns `StatusEdgeCurves(DONE)` (`cxx:1029`).
        // The `bool` returned by `FixShifted` itself is discarded there, so the
        // `fixed` flag follows the status bit, not the return value.
        fixed |= fix_shifted_outcome(wire, face).status.done();
    }
    let edges = edges_of_wire(wire);
    let reg = GeometryRegistry::global();
    for (i, edge) in edges.iter().enumerate() {
        if let Some((c2d, fp2d, lp2d)) = curve_on_surface_oriented(edge, face, false) {
            let (first, last) = BRepTool::edge_parameters(edge);
            if first.is_finite() && last.is_finite() {
                if (first - fp2d).abs() > PCONFUSION || (last - lp2d).abs() > PCONFUSION {
                    reg.set_same_range(&edge.0, false);
                } else if !check_pcurve_range(first, last, c2d.as_ref()) {
                    let is_seam = is_seam_use(&edges, i);
                    reg.remove_pcurves_on_surface(&edge.0, &face.0);
                    let _ = fix_add_pcurve(edge, face, is_seam, precision, &mut cache);
                }
            }
        }
        // `ShapeFix_Wire::FixEdgeCurves` (`ShapeFix_Wire.cxx:953-1005`). This
        // arm IS enabled here, unlike the two face-wire drivers whose preset
        // carries `fix_same_parameter_mode = 0`: this is the naked
        // `ShapeFix_Wire` of `IGESControl_IGESBoundary::Transfer`, built
        // through the `ClearModes` constructor (`ShapeFix_Wire.cxx:167-202`),
        // which leaves `myFixSameParameterMode = -1` (`cxx:185`), so
        // `NeedFix(-1)` is true (`ShapeFix_Root.lxx:100-104`). The face path
        // instead inherits the operator's `FixShape.FixEdgeSameParameterMode =
        // 0` (`ShapeProcess_OperLibrary.cxx:857`), see
        // `WireFixModes::step_face`.
        if need_fix(modes.fix_same_parameter_mode, true) {
            fix_same_parameter(edge, face);
        }
    }
    if need_fix(modes.fix_vertex_tolerance_mode, true) {
        for edge in edges_of_wire(wire) {
            fix_vertex_tolerance_edge_on_face(&edge, face);
        }
    }
    if need_fix(modes.fix_degenerated_mode, true) {
        fixed |= fix_degenerated_all(wire, face, precision);
    }
    if modes.fix_tail_mode <= 0
        && need_fix(modes.fix_notched_edges_mode, reorder_ok)
        && fix_notched_edges(wire, face, SHAPE_FIX_MIN_TOLERANCE, precision)
    {
        fixed = true;
    }
    if reorder_ok && fixed {
        let _ = fix_shifted_wire(wire, face);
    }
    if need_fix(modes.fix_self_intersection_mode, modes.closed_mode) {
        let intersecting_edges = if modes.fix_intersecting_edges_mode == FIX_OR_NOT && !reorder_ok {
            NOT_FIX
        } else {
            modes.fix_intersecting_edges_mode
        };
        let _ = fix_self_intersection(wire, face, precision, need_fix(intersecting_edges, true));
    }
    let _ = fix_reorder_wire(wire, face);
    if need_fix(modes.fix_lacking_mode, reorder_ok) {
        let _ = fix_lacking_all(wire, face, false, precision, precision);
    }
    if need_fix(modes.fix_vertex_tolerance_mode, true) {
        for edge in edges_of_wire(wire) {
            fix_vertex_tolerance_edge_on_face(&edge, face);
        }
    }
    if crate::debug_flags::iges_trace_wire() {
        trace_wire("end-ctor", wire);
    }
}
