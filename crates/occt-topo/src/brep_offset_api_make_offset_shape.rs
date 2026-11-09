//! `BRepOffsetAPI_MakeOffsetShape` - the public front end of the two 3D offset
//! algorithms.
//!
//! Source: `ModelingAlgorithms/TKOffset/BRepOffsetAPI/BRepOffsetAPI_MakeOffsetShape.cxx:22-146`
//! and `.hxx:21-149` (OCCT 8.0.0).
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | `OffsetAlgo_Type` | `hxx:137-141` | [`OffsetAlgoType`] |
//! | `BRepOffsetAPI_MakeOffsetShape()` | `cxx:22-25` | [`BRepOffsetApiMakeOffsetShape::new`] |
//! | `PerformByJoin(S, Offset, Tol, Mode, Inter, SelfInter, Join, RemoveIntEdges)` | `cxx:29-53` | [`BRepOffsetApiMakeOffsetShape::perform_by_join`] |
//! | `PerformBySimple(S, Offset)` | `cxx:57-72` | [`BRepOffsetApiMakeOffsetShape::perform_by_simple`] |
//! | `MakeOffset()` | `cxx:77-80` | [`BRepOffsetApiMakeOffsetShape::make_offset`] |
//! | `Build(Range)` | `cxx:84` | [`BRepOffsetApiMakeOffsetShape::build`] |
//! | `Generated(S)` | `cxx:88-105` | [`BRepOffsetApiMakeOffsetShape::generated`] |
//! | `Modified(S)` | `cxx:110-127` | [`BRepOffsetApiMakeOffsetShape::modified`] |
//! | `IsDeleted(S)` | `cxx:131-138` | [`BRepOffsetApiMakeOffsetShape::is_deleted`] |
//! | `GetJoinType()` | `cxx:142-146` | [`BRepOffsetApiMakeOffsetShape::get_join_type`] |
//!
//! # Translation notes
//!
//! * The class keeps OCCT's two member algorithms side by side (`hxx:146-147`)
//!   and `myLastUsedAlgo` (`hxx:144`) selects which one `Generated` /
//!   `Modified` / `IsDeleted` consult. Only one of them ever runs, because
//!   either `Perform*` sets the selector before it calls into its algorithm.
//! * `NotDone()` (`BRepBuilderAPI_Command.cxx`, `BRepBuilderAPI_Command.hxx:54-55`)
//!   only clears `myDone`; `myShape` keeps whatever the previous `Perform`
//!   left there. The port sets `my_done` and leaves `my_shape` untouched, so a
//!   failed re-run answers the old shape exactly like OCCT.
//! * `PerformByJoin` (`cxx:37-41`) passes `Thickening = false`; the thickening
//!   arm belongs to `BRepOffsetAPI_MakeThickSolid` (see
//!   [`crate::brep_offset_api_make_thick_solid`]).
//! * `NCollection_List` assignments of `Generated` / `Modified`
//!   (`cxx:93`, `cxx:116`) copy the algorithm's list; the port moves the `Vec`
//!   the algorithm returns, which already cleared `myGenerated` there.
//! * `TopoDS_Shape::IsNull` (`cxx:99`, `cxx:119`) is `Option::is_none` and
//!   `IsSame(S)` (`cxx:99`, `cxx:119`) is [`crate::topo_tools::is_same`].
//! * `Build` does nothing (`cxx:84`), so the port's is empty as well.

use occt_core::kernel::geomabs::JoinType;

use crate::brep_offset::{BRepOffsetMakeOffset, BRepOffsetMakeSimpleOffset, Mode};
use crate::shape::TopoShape;
use crate::topo_tools::is_same;

/// `OffsetAlgo_Type` (`BRepOffsetAPI_MakeOffsetShape.hxx:137-141`): which of
/// the two algorithms produced the current `myShape`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OffsetAlgoType {
    /// `OffsetAlgo_NONE`
    None,
    /// `OffsetAlgo_JOIN`: `BRepOffset_MakeOffset` ran.
    Join,
    /// `OffsetAlgo_SIMPLE`: `BRepOffset_MakeSimpleOffset` ran.
    Simple,
}

/// `class BRepOffsetAPI_MakeOffsetShape` (`hxx:39-148`).
pub struct BRepOffsetApiMakeOffsetShape {
    /// `myLastUsedAlgo` (`hxx:144`).
    pub(crate) my_last_used_algo: OffsetAlgoType,
    /// `myOffsetShape` (`hxx:146`).
    pub(crate) my_offset_shape: BRepOffsetMakeOffset,
    /// `mySimpleOffsetShape` (`hxx:147`).
    pub(crate) my_simple_offset_shape: BRepOffsetMakeSimpleOffset,
    /// `BRepBuilderAPI_MakeShape::myGenerated`, the list `Generated` /
    /// `Modified` answer.
    pub(crate) my_generated: Vec<TopoShape>,
    /// `BRepBuilderAPI_MakeShape::myShape`.
    pub(crate) my_shape: Option<TopoShape>,
    /// `BRepBuilderAPI_Command::myDone`.
    pub(crate) my_done: bool,
}

impl Default for BRepOffsetApiMakeOffsetShape {
    fn default() -> Self {
        Self::new()
    }
}

impl BRepOffsetApiMakeOffsetShape {
    /// `BRepOffsetAPI_MakeOffsetShape()` (`cxx:22-25`): only
    /// `myLastUsedAlgo = OffsetAlgo_NONE` is set, the two member algorithms are
    /// default-constructed.
    pub fn new() -> Self {
        Self {
            my_last_used_algo: OffsetAlgoType::None,
            my_offset_shape: BRepOffsetMakeOffset::new(),
            my_simple_offset_shape: BRepOffsetMakeSimpleOffset::new(),
            my_generated: Vec::new(),
            my_shape: None,
            my_done: false,
        }
    }

    /// `PerformByJoin(S, Offset, Tol, Mode, Intersection, SelfInter, Join,
    /// RemoveIntEdges)` (`cxx:29-53`).
    ///
    /// The defaults of `hxx:109-116` are the caller's business; the port takes
    /// every argument. `cxx:37-41` runs before the `IsDone` test, so a failed
    /// offset still records `OffsetAlgo_JOIN`.
    pub fn perform_by_join(
        &mut self,
        s: &TopoShape,
        offset: f64,
        tol: f64,
        mode: Mode,
        intersection: bool,
        self_inter: bool,
        join: JoinType,
        remove_int_edges: bool,
    ) {
        // `cxx:39`
        self.my_done = false;
        // `cxx:40`
        self.my_last_used_algo = OffsetAlgoType::Join;

        // `cxx:42-43`: `Initialize(..., Thickening = false, ...)`.
        self.my_offset_shape.initialize(
            s,
            offset,
            tol,
            mode,
            intersection,
            self_inter,
            join,
            false,
            remove_int_edges,
        );
        // `cxx:44`
        self.my_offset_shape.make_offset_shape();

        // `cxx:46-49`
        if !self.my_offset_shape.is_done() {
            return;
        }

        // `cxx:51-53`
        self.my_shape = self.my_offset_shape.shape().cloned();
        self.my_done = true;
    }

    /// `PerformBySimple(S, Offset)` (`cxx:57-72`).
    pub fn perform_by_simple(&mut self, s: &TopoShape, offset_value: f64) {
        // `cxx:59`
        self.my_done = false;
        // `cxx:60`
        self.my_last_used_algo = OffsetAlgoType::Simple;

        // `cxx:62-63`
        self.my_simple_offset_shape.initialize(s, offset_value);
        self.my_simple_offset_shape.perform();

        // `cxx:65-68`
        if !self.my_simple_offset_shape.is_done() {
            return;
        }

        // `cxx:70-72`
        self.my_shape = self.my_simple_offset_shape.result_shape().cloned();
        self.my_done = true;
    }

    /// `MakeOffset()` (`cxx:77-80`): the underlying
    /// `BRepOffset_MakeOffset`, even when the simple algorithm ran.
    pub fn make_offset(&self) -> &BRepOffsetMakeOffset {
        &self.my_offset_shape
    }

    /// `Build(Range)` (`cxx:84`): does nothing.
    pub fn build(&mut self) {}

    /// `IsDone()` (`BRepBuilderAPI_Command`).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// `Shape()` (`BRepBuilderAPI_MakeShape`).
    pub fn shape(&self) -> Option<&TopoShape> {
        self.my_shape.as_ref()
    }

    /// `Generated(S)` (`cxx:88-105`): `myGenerated` is cleared first, then the
    /// algorithm selected by `myLastUsedAlgo` fills it.
    pub fn generated(&mut self, s: &TopoShape) -> &[TopoShape] {
        // `cxx:90`
        self.my_generated.clear();
        match self.my_last_used_algo {
            // `cxx:91-94`
            OffsetAlgoType::Join => {
                self.my_generated = self.my_offset_shape.generated(s);
            }
            // `cxx:95-102`
            OffsetAlgoType::Simple => {
                // `cxx:97-98`
                if let Some(gen) = self.my_simple_offset_shape.generated(s) {
                    // `cxx:99-101`: `!IsNull() && !IsSame(S)`.
                    if !is_same(&gen, s) {
                        self.my_generated.push(gen);
                    }
                }
            }
            // `OffsetAlgo_NONE`: neither arm runs.
            OffsetAlgoType::None => {}
        }
        &self.my_generated
    }

    /// `Modified(S)` (`cxx:110-127`).
    pub fn modified(&mut self, s: &TopoShape) -> &[TopoShape] {
        // `cxx:112`
        self.my_generated.clear();
        match self.my_last_used_algo {
            // `cxx:113-116`
            OffsetAlgoType::Join => {
                self.my_generated = self.my_offset_shape.modified(s);
            }
            // `cxx:117-124`
            OffsetAlgoType::Simple => {
                // `cxx:119-120`
                if let Some(mod_shape) = self.my_simple_offset_shape.modified(s) {
                    // `cxx:121-123`: only the null test, there is no `IsSame`
                    // test on the `Modified` path.
                    self.my_generated.push(mod_shape);
                }
            }
            OffsetAlgoType::None => {}
        }
        &self.my_generated
    }

    /// `IsDeleted(S)` (`cxx:131-138`): the simple algorithm never deletes.
    pub fn is_deleted(&mut self, s: &TopoShape) -> bool {
        // `cxx:133-136`
        if self.my_last_used_algo == OffsetAlgoType::Join {
            return self.my_offset_shape.is_deleted(s);
        }
        // `cxx:137`
        false
    }

    /// `GetJoinType()` (`cxx:142-146`): always answered by `myOffsetShape`.
    pub fn get_join_type(&self) -> JoinType {
        self.my_offset_shape.get_join_type()
    }
}
