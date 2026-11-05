//! `BRepFill_AdvancedEvolved` (`BRepFill_AdvancedEvolved.cxx:1-1850`,
//! `BRepFill_AdvancedEvolved.hxx:1-126`).
//!
//! Volume-mode evolved sweep: sweep the profile along the spine with
//! [`PipeShell`](crate::brep_fill_pipe_shell::PipeShell), cut the pipe with its
//! own free edges and the two planar lids, and let
//! [`BopMakerVolume`](crate::bop_maker_volume::BopMakerVolume) build the solid
//! the (not necessarily closed) set of faces bounds.
//!
//! The class members are OCCT's (`hxx:114-126`); the performing half is split
//! over the sibling modules exactly as the reference's methods are:
//!
//! | module | reference |
//! | --- | --- |
//! | [`boolean`] | `PerformBoolean` `cxx:74-100`, `IsLid` `cxx:375-404` |
//! | [`spine_profile`] | `GetSpineAndProfile` `cxx:104-371` |
//! | [`sweep_lids`] | `Perform` `cxx:406-503`, `PerformSweep` `cxx:507-526`, `GetLids` `cxx:530-700` |
//! | [`build_solid`] | `BuildSolid` `cxx:704-831`, `UnifyShape` `cxx:835-848`, `ExtractOuterSolid` `cxx:852-906`, `RemoveExcessSolids` `cxx:910-1119` |
//! | [`singularity`] | `RebuildFaces` `cxx:1127-1155`, `MakeEdgeDegenerated` `cxx:1157-1214`, `InsertEDegenerated` `cxx:1216-1428`, `CheckSingularityAndAdd` `cxx:1430-1610` |
//! | [`reduce_tolerance`] | `FindInternals` `cxx:1630-1655`, `RemoveInternalWires` `cxx:1657-1679`, `ProcessVertex` `cxx:1682-1827`, `ReduceVertexTolerance` `cxx:1831-1850` |
//!
//! Translation notes:
//! - OCCT's anonymous error enum (`hxx:104-112`) is a `u32` status with the
//!   same values; `OK` is `UINT_MAX`, so [`AdvancedEvolved::is_done`] tests
//!   that one value and [`AdvancedEvolved::error_code`] reports the rest.
//! - A null `TopoDS_*` / `TopoDS_Compound` member becomes `Option<...>`, the
//!   crate's convention for a null `TShape` (`TopoShape::is_null` is always
//!   false here).
//! - `myDebugShapesPath` is only read by the `#ifdef
//!   BRepFill_AdvancedEvolved_DEBUG` blocks, which are not ported;
//!   [`AdvancedEvolved::set_temporary_directory`] keeps the member for API
//!   parity (`hxx:67`).
//! - `SetRunParallel` has no effect in this port (the algorithms run serially);
//!   the flag is stored so the reference's calls have a receiver.

use crate::shape::{TopoShape, Wire};

mod boolean;
mod build_solid;
mod maps;
mod reduce_tolerance;
mod singularity;
mod spine_profile;
mod sweep_lids;

/// `BRepFill_AdvancedEvolved_Empty` (`hxx:106`).
pub const STATUS_EMPTY: u32 = 0;
/// `BRepFill_AdvancedEvolved_NotPlanarSpine` (`hxx:107`).
pub const STATUS_NOT_PLANAR_SPINE: u32 = 1;
/// `BRepFill_AdvancedEvolved_SweepError` (`hxx:108`).
pub const STATUS_SWEEP_ERROR: u32 = 2;
/// `BRepFill_AdvancedEvolved_NoLids` (`hxx:109`).
pub const STATUS_NO_LIDS: u32 = 3;
/// `BRepFill_AdvancedEvolved_NotSolid` (`hxx:110`).
pub const STATUS_NOT_SOLID: u32 = 4;
/// `BRepFill_AdvancedEvolved_NotVolume` (`hxx:111`).
pub const STATUS_NOT_VOLUME: u32 = 5;
/// `BRepFill_AdvancedEvolved_OK` (`hxx:112`, `UINT_MAX`).
pub const STATUS_OK: u32 = u32::MAX;

/// `class BRepFill_AdvancedEvolved` (`hxx:36-126`).
#[derive(Clone)]
pub struct AdvancedEvolved {
    /// `myErrorStatus` (`hxx:113`).
    status: u32,
    /// `mySpine` (`hxx:115`).
    spine: Option<Wire>,
    /// `myProfile` (`hxx:116`).
    profile: Option<Wire>,
    /// `myPipeShell` (`hxx:117`).
    pipe_shell: Option<TopoShape>,
    /// `myTopBottom` (`hxx:118`): a compound, as the lids can be split on
    /// several faces.
    top_bottom: Option<TopoShape>,
    /// `myResult` (`hxx:120`).
    result: Option<TopoShape>,
    /// `myFuzzyValue` (`hxx:121`).
    fuzzy_value: f64,
    /// `myIsParallel` (`hxx:122`).
    is_parallel: bool,
    /// `myDebugShapesPath` (`hxx:123`); only the non-ported debug dumps read
    /// it.
    debug_shapes_path: String,
}

impl Default for AdvancedEvolved {
    fn default() -> Self {
        Self::new()
    }
}

impl AdvancedEvolved {
    /// `BRepFill_AdvancedEvolved()` (`hxx:41-47`).
    pub fn new() -> Self {
        Self {
            status: STATUS_EMPTY,
            spine: None,
            profile: None,
            pipe_shell: None,
            top_bottom: None,
            result: None,
            fuzzy_value: 0.0,
            is_parallel: true,
            debug_shapes_path: "C:\\Temp".to_string(),
        }
    }

    /// `Perform(theSpine, theProfile, theTolerance, theSolidReq = true)`
    /// (`cxx:406-503`).
    pub fn perform(
        &mut self,
        spine: &Wire,
        profile: &Wire,
        tolerance: f64,
        solid_req: bool,
    ) {
        // `cxx:412-424`.
        self.status = STATUS_EMPTY;
        if self.fuzzy_value < occt_core::precision::CONFUSION {
            self.fuzzy_value = tolerance;
        }

        self.get_spine_and_profile(spine, profile);

        // `cxx:444-446`: `myPipeShell.Nullify(); myTopBottom.Nullify();
        // myResult.Nullify()`.
        self.pipe_shell = None;
        self.top_bottom = None;
        self.result = None;

        self.perform_sweep();
        self.get_lids();

        // `cxx:464-467`.
        if self.status != STATUS_NOT_SOLID {
            return;
        }

        // `cxx:469-471`.
        self.result = self.pipe_shell.clone();
        self.build_solid();

        // `cxx:473-476`: only the `theSolidReq == false` arm continues.
        if self.status != STATUS_OK || solid_req {
            return;
        }

        // `cxx:478-501`: drop the lids from the result and answer the shell of
        // everything else. The shell is created inside the loop, before the
        // `IsLid` test, exactly as the reference does.
        let lids = match &self.top_bottom {
            Some(tb) => crate::topo_tools::map_shapes(tb, &[crate::abs::ShapeType::Face])
                .remove(&crate::abs::ShapeType::Face)
                .unwrap_or_default(),
            None => Vec::new(),
        };
        let builder = crate::builder::TopoBuilder::new();
        let mut shell: Option<crate::shape::Shell> = None;
        if let Some(result) = &self.result {
            for face in crate::topo_tools::explorer_shapes(result, crate::abs::ShapeType::Face) {
                if shell.is_none() {
                    shell = Some(builder.make_shell(&[]));
                }
                let face = crate::shape::Face(face);
                if Self::is_lid(&face, &lids) {
                    continue;
                }
                if let Some(s) = shell.as_mut() {
                    builder.add_face(s, &face);
                }
            }
        }
        if let Some(s) = shell {
            self.result = Some(s.0);
        }
    }

    /// `IsDone(unsigned int* theErrorCode = nullptr)` (`hxx:52-60`).
    pub fn is_done(&self) -> bool {
        self.status == STATUS_OK
    }

    /// The value `IsDone(theErrorCode)` writes out (`hxx:53-58`).
    pub fn error_code(&self) -> u32 {
        self.status
    }

    /// `Shape()` (`hxx:63`): the null `TopoDS_Shape` is `None`.
    pub fn shape(&self) -> Option<&TopoShape> {
        self.result.as_ref()
    }

    /// `SetTemporaryDirectory(thePath)` (`hxx:67`).
    pub fn set_temporary_directory(&mut self, path: &str) {
        self.debug_shapes_path = path.to_string();
    }

    /// `SetParallelMode(theVal)` (`hxx:69`).
    pub fn set_parallel_mode(&mut self, val: bool) {
        self.is_parallel = val;
    }

    /// `myFuzzyValue` as written by the last `Perform` (`cxx:412-417`).
    pub fn fuzzy_value(&self) -> f64 {
        self.fuzzy_value
    }
}
