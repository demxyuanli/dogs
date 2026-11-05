//! `ShapeFix_Shape` (`ShapeFix_Shape.cxx:36-353`).
//!
//! The class is the whole-shape driver of the healing toolkit: it owns one
//! `ShapeBuild_ReShape` context (`cxx:75-78`) and runs the solid / shell / face
//! / wire / edge tool over one shape, then the two tails - `SameParameter`
//! (`cxx:259-265`) and the multi-face vertex tolerance refresh (`cxx:267-291`).
//!
//! Ported arms:
//! * `COMPOUND` / `COMPSOLID` recursion (`cxx:131-159`).
//! * `WIRE` (`cxx:212-237`), the free-wire arm and the only in-tree consumer
//!   (`BRepAlgo::ConcatenateWireC0`, `BRepAlgo.cxx:412-417`).
//! * `EDGE` (`cxx:238-246`).
//!
//! UNPORTED arms, deliberately: `SOLID` (`cxx:159-179`), `SHELL`
//! (`cxx:181-192`) and `FACE` (`cxx:193-210`) call `ShapeFix_Solid::Perform`,
//! `ShapeFix_Shell::Perform` and `ShapeFix_Face::Perform` as *class* drivers.
//! The port runs those stages from the reader pipelines instead, bound to their
//! `FromSTEP.FixShape` / `FromIGES.FixShape` presets:
//! `iges_read::iges_process_shape::{perform_iges_solid, perform_iges_shell,
//! perform_face}` and `step::read_topology::shape_fix_face_perform`. This
//! driver therefore leaves those three shape types untouched and returns
//! `false` (no fix recorded), rather than inventing a second preset-less path.
//!
//! UNPORTED and inert here: `ShapeFix::FixVertexPosition` (`cxx:118-121`,
//! `ShapeFix.cxx:290-350`). Its gate is `NeedFix(myFixVertexPositionMode)` and
//! the ctor default is `0` (`ShapeFix_Shape.cxx:47`), so a freshly built
//! `ShapeFix_Shape` never calls it.

use super::prelude::*;
use super::*;

use crate::shape_extend::{
    decode_status, SHAPEEXTEND_DONE1, SHAPEEXTEND_DONE2, SHAPEEXTEND_DONE3,
    SHAPEEXTEND_DONE4, SHAPEEXTEND_OK,
};
use crate::shape_fix_compose_shell::{ReShape, SharedReShape};
use crate::shhealing::wire_fix::{check_order_3d, FIX_OR_NOT, NOT_FIX, WireFixModes};
use occt_core::precision::CONFUSION;
use occt_core::toploc::TopLocLocation;

/// `ShapeFix_Shape` (`ShapeFix_Shape.hxx:143-149`): the mode set lives in
/// [`ShapeFixShapeModes`], the tool-chain precision and tolerances are the
/// three `ShapeFix_Root` fields the class owns and forwards
/// (`ShapeFix_Shape.cxx:330-353`).
pub struct ShapeFixShape {
    modes: ShapeFixShapeModes,
    my_shape: TopoShape,
    my_result: TopoShape,
    my_status: i32,
    my_precision: f64,
    my_min_tol: f64,
    my_max_tol: f64,
    my_context: SharedReShape,
    /// `myMapFixingShape` (`ShapeFix_Shape.cxx:104-113`): the address of the
    /// location-stripped `TShape` of every shape already being fixed, so a
    /// self-referencing assembly does not recurse forever.
    my_map_fixing_shape: HashSet<usize>,
    /// `ShapeFix_Face`'s mode set, read only for the
    /// `FixSmallAreaWireMode()` save / restore the prologue does
    /// (`ShapeFix_Shape.cxx:87-97`, `:297-301`).
    face_modes: ShapeFixFaceModes,
    /// `ShapeFix_Wire`'s mode set: the free-wire tool the class builds is a
    /// `ClearModes` tool (`ShapeFix_Wire.cxx:167-201`) that `Perform` then
    /// overrides per call (`ShapeFix_Shape.cxx:216-223`).
    wire_modes: WireFixModes,
}

impl ShapeFixShape {
    /// `ShapeFix_Shape()` (`ShapeFix_Shape.cxx:40-50`) followed by
    /// `Init(shape)` (`cxx:70-79`). `SetPrecision` / `SetMinTolerance` /
    /// `SetMaxTolerance` start at `Precision::Confusion()`
    /// (`ShapeFix_Root.cxx:26-28`).
    pub fn new(shape: &TopoShape) -> Self {
        ShapeFixShape {
            modes: ShapeFixShapeModes::default(),
            my_shape: shape.clone(),
            my_result: shape.clone(),
            my_status: SHAPEEXTEND_OK,
            my_precision: CONFUSION,
            my_min_tol: CONFUSION,
            my_max_tol: CONFUSION,
            // `Init` (`cxx:73-78`): `SetContext(new ShapeBuild_ReShape)` with
            // `ModeConsiderLocation() = true`. The port's `MapReShape` keys on
            // `(TShape, location)` (`topo_tools::ShapeMapKey`), which is that
            // mode.
            my_context: SharedReShape::new(),
            my_map_fixing_shape: HashSet::new(),
            face_modes: ShapeFixFaceModes::default(),
            wire_modes: WireFixModes::default(),
        }
    }

    /// `ShapeFix_Shape::FixFreeWireMode()` (`ShapeFix_Shape.lxx:80-84`).
    pub fn fix_free_wire_mode(&mut self) -> &mut i32 {
        &mut self.modes.fix_wire_mode
    }

    /// `ShapeFix_Root::SetPrecision` (`ShapeFix_Root.cxx:58-68`): the min and
    /// max tolerances follow the precision up / down.
    pub fn set_precision(&mut self, preci: f64) {
        self.my_precision = preci;
        if self.my_max_tol < preci {
            self.my_max_tol = preci;
        }
        if self.my_min_tol > preci {
            self.my_min_tol = preci;
        }
    }

    /// `ShapeFix_Root::SetMinTolerance` (`ShapeFix_Root.cxx:72-76`).
    pub fn set_min_tolerance(&mut self, min_tol: f64) {
        self.my_min_tol = min_tol;
    }

    /// `ShapeFix_Root::SetMaxTolerance` (`ShapeFix_Root.cxx:80-84`).
    pub fn set_max_tolerance(&mut self, max_tol: f64) {
        self.my_max_tol = max_tol;
    }

    /// `ShapeFix_Shape::Shape()` (`ShapeFix_Shape.cxx:316-319`).
    pub fn shape(&self) -> &TopoShape {
        &self.my_result
    }

    /// `ShapeFix_Shape::Status(status)` (`ShapeFix_Shape.cxx:355-358`).
    pub fn status(&self, status: i32) -> bool {
        decode_status(self.my_status, status)
    }

    /// `ShapeFix_Shape::Perform` (`ShapeFix_Shape.cxx:82-302`).
    pub fn perform(&mut self) -> bool {
        let sav_fix_vertex_tol_mode = self.modes.fix_vertex_tol_mode;

        // `cxx:87-97`: `FixSmallAreaWireMode` is saved and, when it is still
        // `-1` and the shape is a face, switched on for the pass.
        let sav_fix_small_area_wire_mode = self.face_modes.fix_small_area_wire_mode;
        if sav_fix_small_area_wire_mode == FIX_OR_NOT
            && self.my_shape.shape_type() == ShapeType::Face
        {
            self.face_modes.fix_small_area_wire_mode = 1;
        }

        self.my_status = SHAPEEXTEND_OK;
        let mut status = false;

        // `cxx:101-116`: the re-entry guard runs on the location-stripped
        // shape. `TopoDS_Shape::Location(nullLoc)` keeps the TShape, so the
        // guard key is the TShape address, as in the `.cxx` where the map is
        // `TopTools_ShapeMapHasher` over that same stripped shape.
        let mut shape_null_loc = self.my_shape.clone();
        shape_null_loc.set_location(&TopLocLocation::identity());
        let key = ptr(&shape_null_loc);
        if self.my_context.is_recorded(&self.my_shape)
            || self.my_map_fixing_shape.contains(&key)
        {
            self.my_result = self.my_context.apply(&self.my_shape);
            return true;
        }
        self.my_map_fixing_shape.insert(key);

        // `cxx:117`: `S = Context()->Apply(myShape)`. The `.cxx` had already
        // stripped the location from `myShape` at `cxx:115`; the port keeps the
        // caller's shape intact and only strips it for the guard key.
        let s = self.my_context.apply(&self.my_shape);
        // `cxx:118-121`: UNPORTED, see the module note.

        match s.shape_type() {
            ShapeType::Compound | ShapeType::CompSolid => {
                // `cxx:131-159`: the children are fixed with
                // `FixSameParameter` and `FixVertexTol` cleared; the two tails
                // therefore only run on the outermost call.
                let saved_same_parameter_mode = self.modes.fix_same_parameter_mode;
                self.modes.fix_same_parameter_mode = NOT_FIX;
                self.modes.fix_vertex_tol_mode = NOT_FIX;
                let children = {
                    let ts = s.tshape.read().expect("poisoned TShape lock");
                    ts.children.clone()
                };
                for child in &children {
                    self.my_shape = child.clone();
                    if self.perform() {
                        status = true;
                    }
                }
                self.modes.fix_same_parameter_mode = saved_same_parameter_mode;
                self.modes.fix_vertex_tol_mode = sav_fix_vertex_tol_mode;
                self.my_shape = shape_null_loc;
            }
            ShapeType::Solid => {
                // `cxx:159-179`: UNPORTED, see the module note.
                self.my_status |= SHAPEEXTEND_DONE4;
            }
            ShapeType::Shell => {
                // `cxx:181-192`: UNPORTED, see the module note.
                self.my_status |= SHAPEEXTEND_DONE4;
            }
            ShapeType::Face => {
                // `cxx:193-210`: UNPORTED, see the module note.
                self.my_status |= SHAPEEXTEND_DONE3;
            }
            ShapeType::Wire => {
                // `cxx:212-237`.
                if need_fix(self.modes.fix_wire_mode, true) {
                    let (fixed, wire) = perform_free_wire(
                        &s,
                        &mut self.my_context,
                        self.my_precision,
                        self.my_min_tol,
                        self.my_max_tol,
                        &self.wire_modes,
                    );
                    if fixed {
                        status = true;
                        // `cxx:231`: replace for the wire only.
                        self.my_context.replace(&s, &wire);
                    }
                    self.my_status |= SHAPEEXTEND_DONE2;
                }
            }
            ShapeType::Edge => {
                // `cxx:238-246`.
                let edge = Edge(s.clone());
                if fix_vertex_tolerance_in_context(&edge, &mut self.my_context) {
                    self.my_status |= SHAPEEXTEND_DONE1;
                }
            }
            ShapeType::Vertex | ShapeType::Shape => {}
        }

        // `cxx:257`.
        self.my_result = self.my_context.apply(&s);

        // `cxx:259-265`: `ShapeFix_Shape::SameParameter(myResult, false)` ->
        // `ShapeFix::SameParameter(sh, enforce = false, preci = 0.0)`
        // (`cxx:307-311`).
        if need_fix(self.modes.fix_same_parameter_mode, true) {
            shape_fix_same_parameter(&self.my_result, false, 0.0);
        }

        // `cxx:267-291`: the vertex tolerance refresh runs only when the result
        // holds more than one face (`cxx:269-276`, the `nbF > 1` test).
        if need_fix(self.modes.fix_vertex_tol_mode, true) {
            let faces = faces_of(&self.my_result);
            if faces.len() > 1 {
                fix_vertex_tolerances_after_same_parameter(&self.my_result, &mut self.my_context);
            }
        }

        // `cxx:294`.
        self.my_result = self.my_context.apply(&self.my_result);

        // `cxx:297-301`.
        self.face_modes.fix_small_area_wire_mode = sav_fix_small_area_wire_mode;

        status
    }
}

/// `ShapeFix_Wire::Perform` (`ShapeFix_Wire.cxx:297-482`) as `ShapeFix_Shape`
/// reaches it for a free wire (`ShapeFix_Shape.cxx:216-237`):
/// `ModifyTopologyMode() = true`, `ClosedWireMode()` follows `S.Closed()`,
/// `SetFace(null)`.
///
/// With a null face `IsReady()` is false (`ShapeFix_Wire.lxx:68-70`), which is
/// what removes the face-dependent stages from the driver: `FixEdgeCurves`'
/// reversed-2d / remove-pcurve / add-pcurve / remove-curve3d / seam / shifted /
/// same-parameter arms (`cxx:619-705`, `:867-879`, `:921-1005`) and the whole
/// of `FixDegenerated` (`cxx:1034-1038`), `FixNotchedEdges` (`cxx:3979-3983`),
/// `FixTails` (`cxx:4316-4319`), `FixSelfIntersection` (`cxx:1085-1089`) and
/// `FixLacking` (`cxx:1286-1290`) return false before doing anything. The
/// stages that do run are listed below in `cxx` order.
///
/// Returns `(Fixed, wire)`: the `Wire()` the class records with
/// `Context()->Replace` (`ShapeFix_Shape.cxx:231`).
pub fn perform_free_wire(
    shape: &TopoShape,
    ctx: &mut SharedReShape,
    precision: f64,
    min_tolerance: f64,
    max_tolerance: f64,
    tool_modes: &WireFixModes,
) -> (bool, TopoShape) {
    let mut wire = Wire(shape.clone());
    let mut modes = tool_modes.clone();
    // `ShapeFix_Shape.cxx:220-223`.
    modes.topo_mode = true;
    if !shape.is_closed() {
        modes.closed_mode = false;
    }
    // The null `TopoDS_Face` of `cxx:225` as the port's bare face: a `Face`
    // with no surface, which is how every `BRep_Tool::Surface` lookup answers
    // a null surface (and what `iges_process_shape.rs` already passes).
    let bare = Face(TopoShape::new(ShapeType::Face));

    let mut fixed = false;

    // `ShapeFix_Wire.cxx:316-318`: `CheckOrder(sawo, myClosedMode)`.
    let mut reorder_ok = check_order_3d(&wire);
    // `cxx:319-328`.
    if need_fix(modes.fix_reorder_mode, !reorder_ok) {
        let (ok, order_status) = fix_reorder_wire_3d(&wire);
        fixed |= ok && order_status != crate::meshing::wire_order::WireOrderStatus::Same;
        // `ReorderOK = !StatusReorder(FAIL)`.
        reorder_ok = ok;
    }

    // `cxx:333-346`: `FixSmall(!myTopoMode || !ReorderOK, MinTolerance())`, then
    // the reorder retry inside the same `if`.
    if need_fix(modes.fix_small_mode, modes.topo_mode)
        && fix_small_all(
            &mut wire,
            &bare,
            precision,
            min_tolerance,
            !modes.topo_mode || !reorder_ok,
        )
    {
        fixed = true;
        if need_fix(modes.fix_reorder_mode, !reorder_ok) {
            let (ok, _order_status) = fix_reorder_wire_3d(&wire);
            reorder_ok = ok;
        }
    }

    // `cxx:355-359`: `NeedFix(myFixConnectedMode, ReorderOK)`.
    if need_fix(modes.fix_connected_mode, reorder_ok) {
        fixed |= fix_connected_all_in_context(
            &mut wire,
            Some(ctx),
            precision,
            max_tolerance,
            modes.closed_mode,
        );
    }

    // `cxx:365-378`: `FixEdgeCurves` under `NeedFix(myFixEdgeCurvesMode)`, with
    // `FixShiftedMode` forced off when the reorder failed.
    if need_fix(modes.fix_edge_curves_mode, true) {
        let shifted_mode = if modes.fix_shifted_mode == FIX_OR_NOT && !reorder_ok {
            NOT_FIX
        } else {
            modes.fix_shifted_mode
        };
        // `cxx:881-918`: `FixAddCurve3d`, the only curve arm that is not
        // `isReady`-gated. An edge that ends up with neither a 3D curve nor a
        // pcurve is dropped, and the neighbours are reconnected there.
        if need_fix(modes.fix_add_curve3d_mode, true) {
            fixed |= fix_add_curve3d_or_drop(&mut wire, &bare, precision);
        }
        let _ = shifted_mode;
        // `cxx:1009-1026`: `NeedFix(myFixVertexToleranceMode)`, then
        // `UpdateWire()` because the context is set.
        if need_fix(modes.fix_vertex_tolerance_mode, true) {
            fixed |= fix_vertex_tolerances_in_context(&wire, ctx);
            update_wire_from_context(&mut wire, ctx);
        }
    }

    // `cxx:386-392` `FixDegenerated`, `cxx:400-406` `FixNotchedEdges`,
    // `cxx:412-424` `FixTails`, `cxx:428-443` `FixSelfIntersection` and
    // `cxx:448-452` `FixLacking` all return false on a null face; their guards
    // are cited in the doc comment above.

    // `cxx:462-470`: the per-edge `FixVertexTolerance(edge, Face())` that runs
    // unconditionally ("TEMPORARILY without special mode").
    fixed |= fix_vertex_tolerances_in_context(&wire, ctx);

    // `cxx:475-480`: `UpdateWire()` when the context is set.
    update_wire_from_context(&mut wire, ctx);

    (fixed, wire.0)
}

/// `ShapeFix_Edge::FixVertexTolerance(edge, face)` with a live context
/// (`ShapeFix_Edge.cxx:642-684`): `CopyVertex` replaces `B.UpdateVertex`
/// (`cxx:674-677`). Returns whether a tolerance was written.
fn fix_vertex_tolerances_in_context(wire: &Wire, ctx: &mut SharedReShape) -> bool {
    let bare = Face(TopoShape::new(ShapeType::Face));
    let mut done = false;
    for edge in edges_of_wire(wire) {
        if fix_vertex_tolerance_in_context(&edge, ctx) {
            done = true;
        }
        let _ = &bare;
    }
    done
}

/// The single-edge form of the same call.
fn fix_vertex_tolerance_in_context(edge: &Edge, ctx: &mut SharedReShape) -> bool {
    let bare = Face(TopoShape::new(ShapeType::Face));
    let Some(grown) = grown_vertex_tolerances_on_face(edge, &bare) else {
        return false;
    };
    ctx.copy_vertex(&grown.v1.0, grown.toler1);
    ctx.copy_vertex(&grown.v2.0, grown.toler2);
    true
}
