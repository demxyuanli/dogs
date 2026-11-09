//! `BRepOffsetAPI_MakeThickSolid` - hollow a solid by offsetting it and
//! removing the closing faces.
//!
//! Source: `ModelingAlgorithms/TKOffset/BRepOffsetAPI/BRepOffsetAPI_MakeThickSolid.cxx:23-117`
//! and `.hxx:21-124` (OCCT 8.0.0).
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | `BRepOffsetAPI_MakeThickSolid()` | `cxx:23-27` | [`BRepOffsetApiMakeThickSolid::new`] |
//! | `MakeThickSolidByJoin(S, ClosingFaces, ...)` | `cxx:31-63` | [`BRepOffsetApiMakeThickSolid::make_thick_solid_by_join`] |
//! | `MakeThickSolidBySimple(S, Offset)` | `cxx:66-82` | [`BRepOffsetApiMakeThickSolid::make_thick_solid_by_simple`] |
//! | `Build(Range)` | `cxx:86` | [`BRepOffsetApiMakeThickSolid::build`] |
//! | `Modified(F)` | `cxx:90-117` | [`BRepOffsetApiMakeThickSolid::modified`] |
//!
//! # Translation notes
//!
//! * OCCT derives this class from `BRepOffsetAPI_MakeOffsetShape` (`hxx:47`)
//!   and only adds the two `MakeThickSolidBy*` drivers plus the `Modified`
//!   override; the port keeps the base object in [`base`](BRepOffsetApiMakeThickSolid::base)
//!   and forwards `IsDone` / `Shape` / `Generated` / `IsDeleted` /
//!   `MakeOffset` / `GetJoinType` to it.
//! * The constructor only flips the simple algorithm to the solid arm
//!   (`cxx:25-26`, `SetBuildSolidFlag(true)`), which makes
//!   [`BRepOffsetMakeSimpleOffset::perform`] run its `BuildMissingWalls`
//!   step.
//! * `TopoDS::Face(it.Value())` (`cxx:50`) throws `Standard_TypeMismatch` for a
//!   non-face closing shape. The exception is not caught here, so the port
//!   answers `Err` after the same state changes (`NotDone`,
//!   `myLastUsedAlgo = OffsetAlgo_JOIN`, `Initialize`, the faces already
//!   added).
//! * `MakeThickSolid` (`BRepOffset_MakeOffset.cxx:1083-1177`) is
//!   [`BRepOffsetMakeOffset::make_thick_solid`]; it reads
//!   `myThickening = false` and the closing faces registered by `AddFace`.
//! * `Modified` (`cxx:90-117`) clears `myGenerated` first. It then replaces
//!   the base `Modified` wholesale: the joining arm answers images only for a
//!   *closing* face, and reverses each of them (`cxx:100-107`).
//!   `OffsetFacesFromShapes()` is
//!   [`BRepOffsetMakeOffset::offset_faces_from_shapes`], `ClosingFaces()` is
//!   [`BRepOffsetMakeOffset::closing_faces`].
//! * `Build` does nothing (`cxx:86`).

use occt_core::kernel::geomabs::JoinType;

use crate::brep_offset::Mode;
use crate::brep_offset_api_make_offset_shape::{BRepOffsetApiMakeOffsetShape, OffsetAlgoType};
use crate::shape::{Face, TopoShape};

/// `class BRepOffsetAPI_MakeThickSolid` (`hxx:47-123`).
pub struct BRepOffsetApiMakeThickSolid {
    /// The `BRepOffsetAPI_MakeOffsetShape` base sub-object (`hxx:47`).
    pub base: BRepOffsetApiMakeOffsetShape,
}

impl Default for BRepOffsetApiMakeThickSolid {
    fn default() -> Self {
        Self::new()
    }
}

impl BRepOffsetApiMakeThickSolid {
    /// `BRepOffsetAPI_MakeThickSolid()` (`cxx:23-27`): the only action is
    /// `mySimpleOffsetShape.SetBuildSolidFlag(true)`.
    pub fn new() -> Self {
        let mut base = BRepOffsetApiMakeOffsetShape::new();
        // `cxx:25-26`
        base.my_simple_offset_shape.set_build_solid_flag(true);
        Self { base }
    }

    /// `MakeThickSolidByJoin(S, ClosingFaces, Offset, Tol, Mode, Intersection,
    /// SelfInter, Join, RemoveIntEdges)` (`cxx:31-63`).
    pub fn make_thick_solid_by_join(
        &mut self,
        s: &TopoShape,
        closing_faces: &[TopoShape],
        offset: f64,
        tol: f64,
        mode: Mode,
        intersection: bool,
        self_inter: bool,
        join: JoinType,
        remove_int_edges: bool,
    ) -> Result<(), String> {
        // `cxx:43`
        self.base.my_done = false;
        // `cxx:44`
        self.base.my_last_used_algo = OffsetAlgoType::Join;

        // `cxx:46-47`: `Initialize(..., Thickening = false, ...)`.
        self.base.my_offset_shape.initialize(
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
        // `cxx:48-51`
        for f in closing_faces {
            let Some(face) = Face::wrap(f.clone()) else {
                return Err("Standard_TypeMismatch: TopoDS::Face".into());
            };
            self.base.my_offset_shape.add_face(&face);
        }

        // `cxx:53`
        self.base.my_offset_shape.make_thick_solid();
        // `cxx:54-57`
        if !self.base.my_offset_shape.is_done() {
            return Ok(());
        }

        // `cxx:59-61`
        self.base.my_shape = self.base.my_offset_shape.shape().cloned();
        self.base.my_done = true;
        Ok(())
    }

    /// `MakeThickSolidBySimple(S, Offset)` (`cxx:66-82`).
    pub fn make_thick_solid_by_simple(&mut self, s: &TopoShape, offset_value: f64) {
        // `cxx:68`
        self.base.my_done = false;
        // `cxx:69`
        self.base.my_last_used_algo = OffsetAlgoType::Simple;

        // `cxx:71-72`
        self.base.my_simple_offset_shape.initialize(s, offset_value);
        self.base.my_simple_offset_shape.perform();

        // `cxx:74-77`
        if !self.base.my_simple_offset_shape.is_done() {
            return;
        }

        // `cxx:79-81`
        self.base.my_shape = self.base.my_simple_offset_shape.result_shape().cloned();
        self.base.my_done = true;
    }

    /// `Build(Range)` (`cxx:86`): does nothing (`cxx:86` overrides the also
    /// empty base `Build`).
    pub fn build(&mut self) {}

    /// `Modified(F)` (`cxx:90-117`): the closing-face images of the joining
    /// arm, reversed, or the simple arm's single modified shape.
    pub fn modified(&mut self, f: &TopoShape) -> &[TopoShape] {
        // `cxx:92`
        self.base.my_generated.clear();

        // `cxx:94-108`
        if self.base.my_last_used_algo == OffsetAlgoType::Join
            && self.base.my_offset_shape.offset_faces_from_shapes().has_image(f)
        {
            if self.base.my_offset_shape.closing_faces().contains(f) {
                // `cxx:98-99`
                let mut list = Vec::new();
                self.base
                    .my_offset_shape
                    .offset_faces_from_shapes()
                    .last_image(f, &mut list);
                // `cxx:101-107`: `it.ChangeValue().Reverse()`.
                for s in list.iter_mut() {
                    s.reverse();
                }
                self.base.my_generated = list;
            }
        } else if self.base.my_last_used_algo == OffsetAlgoType::Simple {
            // `cxx:109-114`
            if let Some(mod_shape) = self.base.my_simple_offset_shape.modified(f) {
                self.base.my_generated.push(mod_shape);
            }
        }

        // `cxx:116`
        &self.base.my_generated
    }

    /// `IsDone()` (`BRepBuilderAPI_Command`, inherited).
    pub fn is_done(&self) -> bool {
        self.base.is_done()
    }

    /// `Shape()` (`BRepBuilderAPI_MakeShape`, inherited).
    pub fn shape(&self) -> Option<&TopoShape> {
        self.base.shape()
    }

    /// `Generated(S)` (`BRepOffsetAPI_MakeOffsetShape::Generated`, inherited).
    pub fn generated(&mut self, s: &TopoShape) -> &[TopoShape] {
        self.base.generated(s)
    }

    /// `IsDeleted(S)` (`BRepOffsetAPI_MakeOffsetShape::IsDeleted`,
    /// inherited): the base's `IsDeleted` never consults the simple arm.
    pub fn is_deleted(&mut self, s: &TopoShape) -> bool {
        self.base.is_deleted(s)
    }

    /// `MakeOffset()` (`BRepOffsetAPI_MakeOffsetShape::MakeOffset`,
    /// inherited).
    pub fn make_offset(&self) -> &crate::brep_offset::BRepOffsetMakeOffset {
        self.base.make_offset()
    }

    /// `GetJoinType()` (`BRepOffsetAPI_MakeOffsetShape::GetJoinType`,
    /// inherited).
    pub fn get_join_type(&self) -> JoinType {
        self.base.get_join_type()
    }
}
