use super::super::*;
use super::statics::*;


/// STEP post-pass after pcurve association.
/// Source: `TranslateEdgeLoop::CheckPCurves` first check, then
/// `ShapeFix_Wire::Perform` FixReorder (`cxx:317-325` / `cxx:487`) then
/// `FromSTEP.FixShape.FixShiftedMode` (default on) / `cxx:939` FixShifted,
/// then SameRange / CheckPCurveRange / FixAddPCurve / FixSameParameter
/// (`ShapeFix_Wire.cxx:938-995`), then `ShapeFix_Wire::FixDegenerated`
/// (`cxx:386-392`, `fix_degenerated_all` below).
///
/// `TranslateEdgeLoop.cxx:105-112`: for a plane face `CheckPCurves` calls
/// `RemovePCurves(aWire, aFace)` and returns, so the pcurve-range healing does
/// not run on a plane's edges and their tolerances keep the imported values
/// (linkrods FACE 35 edge 4 keeps 1.2432146551310009e-07). That removal is
/// ported in `check_pcurve_rep_range` below (`cxx:110-114`), together with the
/// non-planar `XSAlgo_ShapeProcessor::CheckPCurve` call at `cxx:175`. The
/// `ShapeFix_Wire` passes that follow are independent of the face type.
pub fn check_pcurves_and_shift(wire: &mut Wire, face: &Face, preci: f64, modes: &WireFixModes) {
    check_pcurve_rep_range(wire, face, preci);
    // `ShapeFix_Wire::Perform` (`cxx:317-325`): FixReorder before FixEdgeCurves.
    // `ReorderOK` starts as `myAnalyzer->CheckOrder(sawo, myClosedMode) == 0`
    // (`cxx:319`) and `FixReorder` only runs under
    // `NeedFix(myFixReorderMode, !ReorderOK)` (`cxx:320`); for an already
    // ordered wire the order status is `Same` and `FixReorder` applies nothing
    // (`ShapeFix_Wire.cxx:1360-1363`), so the combined helper is equivalent and
    // its return value is the post-`FixReorder` `ReorderOK` (`cxx:326`).
    let mut reorder_ok = fix_reorder_wire(wire, face);
    // `ShapeFix_Wire::Perform` (`cxx:333-347`): FixSmall under
    // `NeedFix(myFixSmallMode, myTopoMode)`. `FromSTEP.FixShape.FixSmallMode` is
    // "-1" (`STEPControl_Controller.cxx:233`) and `myTopoMode` is forced true by
    // `ShapeFix_Shape.cxx:200` for the face path, so the pass runs, with
    // `lockvtx = !myTopoMode || !ReorderOK` and `precsmall = MinTolerance()`.
    if need_fix(modes.fix_small_mode, modes.topo_mode)
        && fix_small_all(
            wire,
            face,
            preci,
            SHAPE_FIX_MIN_TOLERANCE,
            !modes.topo_mode || !reorder_ok,
        )
    {
        // `cxx:339-343`: after a removal the wire may reorder better; the retry
        // is gated by `NeedFix(myFixReorderMode = -1, !ReorderOK)` evaluated with
        // the ReorderOK from before FixSmall.
        if !reorder_ok {
            reorder_ok = fix_reorder_wire(wire, face);
        }
    }
    // `ShapeFix_Wire::Perform` (`cxx:353-357`): FixConnected under
    // `NeedFix(myFixConnectedMode, ReorderOK)`.
    if need_fix(modes.fix_connected_mode, reorder_ok) {
        let _ = fix_connected_all(wire, preci, SHAPE_FIX_MAX_TOLERANCE);
    }
    let edges = edges_of_wire(wire);
    // `ShapeFix_Wire::FixEdgeCurves` (`cxx:652-670`): FixAddPCurve runs for
    // every edge of the wire. `FromSTEP.exec.op = FixShape`
    // (`STEPControl_Controller.cxx:201`) routes every imported face through
    // `ShapeFix_Shape` -> `ShapeFix_Face` -> `ShapeFix_Wire::Perform`, so a
    // face whose STEP edges carry a 3D curve only (no `SURFACE_CURVE` /
    // `PCURVE`, as in ATU01038) is projected onto the surface here. With no
    // pcurve stored the face keeps the natural surface bounds
    // (`BRepTools::AddUVBounds`, `BRepTools.cxx:139-153`), which are infinite
    // for a cylinder/cone. The `cxx:673-863` singularity split is called from
    // the IGES first round (`perform_wire_first_round`), not from this STEP
    // pass.
    // `ShapeFix_Edge::myProjector` lives on the `ShapeFix_Edge` that
    // `ShapeFix_Wire` keeps for the whole wire (`ShapeFix_Wire.cxx:280-300`), so
    // one `myCache` covers every edge of this face's wire.
    let mut cache = crate::pcurve_full::ProjectorCache::new();
    // `ShapeFix_Wire::Perform` (`cxx:368-380`) -> `FixEdgeCurves`
    // (`ShapeFix_Wire.cxx:881-1029`) under `NeedFix(myFixEdgeCurvesMode)`;
    // inside it `FixShifted` is forced off while `ReorderOK` is false
    // (`cxx:369-380` `myFixShiftedMode == -1 && !ReorderOK -> 0`).
    let edge_curves = need_fix(modes.fix_edge_curves_mode, true);
    let shifted_mode = if modes.fix_shifted_mode == FIX_OR_NOT && !reorder_ok {
        NOT_FIX
    } else {
        modes.fix_shifted_mode
    };
    if edge_curves && need_fix(modes.fix_add_pcurve_mode, true) {
        for (i, e) in edges.iter().enumerate() {
            let is_seam = is_seam_use(&edges, i);
            let _ = fix_add_pcurve(e, face, is_seam, preci, &mut cache);
        }
    }
    // `ShapeFix_Wire.cxx:938-951`: `isReady && NeedFix(myFixShiftedMode)`.
    if edge_curves && need_fix(shifted_mode, true) {
        let _ = fix_shifted_wire(wire, face);
    }
    // `ShapeFix_Wire::FixEdgeCurves` (`ShapeFix_Wire.cxx:953-1005`) under
    // `isReady && NeedFix(myFixSameParameterMode)`; the field is set from the
    // *separate* key `FixEdgeSameParameterMode`
    // (`ShapeProcess_OperLibrary.cxx:857`), which is `0` for STEP
    // (`DE_ShapeFixParameters.hxx:75` `FixMode::NotFix`,
    // `DESTEP_Parameters.cxx:90-93`, `XSAlgo_ShapeProcessor.cxx:681`, legacy
    // `STEPControl_Controller.cxx:247`), so `NeedFix(0)` is false and OCCT runs
    // neither the `CheckPCurveRange` re-projection below nor
    // `FixSameParameter(Edge, Face())`. OCCT's shape-level equivalent
    // `ShapeFix_Shape.cxx:259` -> `ShapeFix::SameParameter(Result, false, 0.0)`
    // (`ShapeFix.cxx:70-258`) is what the port covers with
    // `same_parameter_shape` (`BRepLib::SameParameter`).
    let reg = GeometryRegistry::global();
    let same_parameter_arm = edge_curves && need_fix(modes.fix_same_parameter_mode, true);
    for (i, e) in edges.iter().enumerate() {
        if same_parameter_arm {
            if let Some((c2d, fp2d, lp2d)) = curve_on_surface_oriented(e, face, false) {
                let (first, last) = BRepTool::edge_parameters(e);
                if first.is_finite() && last.is_finite() {
                    let (_cfp, _clp) = (c2d.first_parameter(), c2d.last_parameter());
                    if (first - fp2d).abs() > PCONFUSION || (last - lp2d).abs() > PCONFUSION {
                        reg.set_same_range(&e.0, false);
                    } else if !check_pcurve_range(first, last, c2d.as_ref()) {
                        // `ShapeFix_Wire.cxx:970-985`: when the pcurve cannot
                        // cover the edge's 3D range OCCT does exactly two
                        // things — `ShapeBuild_Edge().RemovePCurve(E, S, L)`
                        // followed by `myFixEdge->FixAddPCurve(E, face, isSeam,
                        // myAnalyzer->Surface(), Precision())`, i.e. it
                        // *re-projects* the pcurve from the 3D curve. There is
                        // no range-remapping branch here, so neither is there
                        // one in the port: a stored pcurve whose 2D points fall
                        // outside the surface domain (e.g. the spherical pcurve
                        // of `data/Offset.step` `#1290`, whose line sits at
                        // v = -3*pi/4, outside v in [-pi/2, pi/2]) must be
                        // replaced by the projection, not re-parameterised.
                        let is_seam = is_seam_use(&edges, i);
                        reg.remove_pcurves_on_surface(&e.0, &face.0);
                        let _ = fix_add_pcurve(e, face, is_seam, preci, &mut cache);
                    }
                }
            }
            fix_same_parameter(e, face);
        }
    }
    // `ShapeFix_Wire::Perform` (`cxx:386-392`): FixDegenerated under
    // `NeedFix(myFixDegeneratedMode)`.
    if need_fix(modes.fix_degenerated_mode, true) {
        fix_degenerated_all(wire, face, preci);
    }
    // `ShapeFix_Wire::Perform` (`cxx:400`): FixNotchedEdges under
    // `myFixTailMode <= 0 && NeedFix(myFixNotchedEdgesMode, ReorderOK)`.
    // `FromSTEP.FixShape.FixTailMode` is "0" (`STEPControl_Controller.cxx:249`)
    // and `FixNotchedEdgesMode` is "-1" (`STEPControl_Controller.cxx:248`), so
    // the gate is ReorderOK. `cxx:402-406`: a fix is followed by FixShifted.
    if modes.fix_tail_mode <= 0
        && need_fix(modes.fix_notched_edges_mode, reorder_ok)
        && fix_notched_edges(wire, face, SHAPE_FIX_MIN_TOLERANCE, SHAPE_FIX_MAX_TOLERANCE)
    {
        let _ = fix_shifted_wire(wire, face);
    }
    // `ShapeFix_Wire::Perform` (`cxx:412`): FixTails under
    // `myFixTailMode != 0`. The STEP DE default is
    // `DE_ShapeFixParameters.hxx:77` `FixTailMode = FixMode::NotFix = 0`
    // (`DESTEP_Parameters.cxx:92-95`), so OCCT skips this pass in the STEP
    // read-in and porting it here would deviate. Intentionally absent. (The
    // IGES DE default is `-1`, so IGES does call `FixTails`, but its body
    // returns false at `ShapeFix_Wire.cxx:4317` because `MaxTailWidth` is `-1`;
    // see `perform_wire_first_round` below.)
    // `ShapeFix_Wire::Perform` (`cxx:428-443`): FixSelfIntersection under
    // `NeedFix(myFixSelfIntersectionMode, myClosedMode)`; `FixSelfIntersectionMode`
    // is "-1" (`STEPControl_Controller.cxx:238`), so with `myClosedMode = true`
    // (`ShapeFix_Wire.cxx:171`) the pass is requested, including its
    // non-adjacent tail (`ShapeFix_Wire.cxx:1200-1223`) with
    // `FromSTEP.FixShape.FixNonAdjacentIntersectingEdgesMode` = "-1"
    // (`STEPControl_Controller.cxx:254-257`) and the `maxtol = 1.0` default
    // (`ShapeFix_IntersectionTool.hxx:46`).
    //
    // `cxx:430-441`: `Perform` keeps the adjacent-edge pass (`FixIntersectingEdges`)
    // *off* while `ReorderOK` is false - `myFixIntersectingEdgesMode == -1 &&
    // !ReorderOK` is forced to `0`, so `NeedFix(0)` is false. Only the
    // single-edge pass (`FixSelfIntersectingEdge`) still runs.
    //
    // `ShapeFix_Face::Perform` (`cxx:369-374`) forces this mode off for its
    // first wire round, so on that path the whole block is skipped.
    if need_fix(modes.fix_self_intersection_mode, modes.closed_mode) {
        let intersecting_edges = if modes.fix_intersecting_edges_mode == FIX_OR_NOT && !reorder_ok {
            NOT_FIX
        } else {
            modes.fix_intersecting_edges_mode
        };
        let _ = fix_self_intersection(wire, face, preci, need_fix(intersecting_edges, true));
    }
    // `ShapeFix_Wire::Perform` (`cxx:448-456`): FixLacking under
    // `NeedFix(myFixLackingMode, ReorderOK)`; `NeedFix(flag, need)` is
    // `flag < 0 ? need : flag > 0` (`ShapeFix_Root.lxx:101-104`), and
    // `FromSTEP.FixShape.FixLackingMode` is "-1"
    // (`STEPControl_Controller.cxx:237`), so OCCT runs the pass exactly when
    // ReorderOK is true.
    // `ShapeFix_Face::Perform` (`ShapeFix_Face.cxx:369-374`) *disables*
    // `FixLacking` for its first wire round and restores it (`cxx:457`) for the
    // post-`FixMissingSeam` second round (`cxx:509+`), which is what
    // `WireFixModes::first_wire_round` / `second_wire_round` encode. On a
    // periodic face whose bounds are single closed seam edges the pass inserted
    // a spurious back-going edge, which made `FixMissingSeam`'s `check_wire`
    // sum to 0 and never merge the two wires (see
    // specs/_a3n00_gap_analysis.md 9.95).
    if need_fix(modes.fix_lacking_mode, reorder_ok) {
        let _ = fix_lacking_all(wire, face, false, preci, SHAPE_FIX_MAX_TOLERANCE);
    }
}
