//! `BRepOffsetAPI_MakeDraft` - the draft surface along a wire, driven
//! through `BRepFill_Draft`.
//!
//! Source: `ModelingAlgorithms/TKOffset/BRepOffsetAPI/BRepOffsetAPI_MakeDraft.cxx:19-85`
//! and `.hxx:36-99` (OCCT 8.0.0).
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | `BRepOffsetAPI_MakeDraft(Shape, Dir, Angle)` | `cxx:19-24` | [`BRepOffsetApiMakeDraft::new`] |
//! | `SetOptions(Style, AngleMin, AngleMax)` | `cxx:26-36` | [`BRepOffsetApiMakeDraft::set_options`] |
//! | `SetDraft(IsInternal)` | `cxx:38-41` | [`BRepOffsetApiMakeDraft::set_draft`] |
//! | `Perform(LengthMax)` | `cxx:43-51` | [`BRepOffsetApiMakeDraft::perform`] |
//! | `Perform(Surface, KeepInsideSurface)` | `cxx:53-62` | [`BRepOffsetApiMakeDraft::perform_on_surface`] |
//! | `Perform(StopShape, KeepOutSide)` | `cxx:64-73` | [`BRepOffsetApiMakeDraft::perform_on_shape`] |
//! | `Shell()` | `cxx:75-78` | [`BRepOffsetApiMakeDraft::shell`] |
//! | `Generated(S)` | `cxx:80-85` | [`BRepOffsetApiMakeDraft::generated`] |
//!
//! # Translation notes
//!
//! * `myDraft(Shape, Dir, Angle)` (`cxx:23`) is [`BrepFillDraft::new`]; its
//!   `Standard_ConstructionError` arms (a shape that is neither a wire, a face
//!   nor a shell with free boundaries, `BRepFill_Draft.cxx:200-262`) become the
//!   `Err` of that constructor, and the reference's `NotDone()` (`cxx:24`)
//!   is `my_done = false`.
//! * `SetOptions` (`cxx:28-35`) does **not** forward `Style` verbatim: it
//!   starts from `BRepFill_Right` and switches to `BRepFill_Round` only for
//!   `BRepBuilderAPI_RoundCorner`. `Transformed` and `RightCorner` therefore
//!   both select `BRepFill_Right`, exactly as the reference's two-line body
//!   does.
//! * `Perform` (`cxx:43-51`, `:53-62`, `:64-73`) is `Done()` + `myShape =
//!   myDraft.Shape()` **only** when `myDraft.IsDone()`; a failed `BRepFill_Draft`
//!   leaves the API object not done and `myShape` untouched, and the port keeps
//!   that shape (an `Err` from the port's `BrepFillDraft::perform*` is the
//!   `Standard_ConstructionError` / `gp` exceptions of the reference and is
//!   returned before the `IsDone` test).
//! * `Shell()` returns a `TopoDS_Shell` by value and a null handle when the
//!   draft never ran (`BRepFill_Draft.cxx:963-968`); the port answers `None`.
//! * `Generated(S)` forwards to `BRepFill_Draft::Generated`, which raises
//!   `Standard_TypeMismatch` for a non-edge `S` (`BRepFill_Draft.cxx:928`) -
//!   the port's `Err`.

use std::sync::Arc;

use occt_core::gp::GpDir;
use occt_geom::Surface;

use crate::brep_builder_api_transition_mode::TransitionMode;
use crate::brep_fill_draft::BrepFillDraft;
use crate::brep_fill_sweep::TransitionStyle;
use crate::shape::{Shell, TopoShape};

/// `class BRepOffsetAPI_MakeDraft` (`hxx:36-99`).
pub struct BRepOffsetApiMakeDraft {
    /// `myDraft` (`hxx:97`).
    my_draft: BrepFillDraft,
    /// `BRepBuilderAPI_MakeShape::myDone`.
    my_done: bool,
    /// `BRepBuilderAPI_MakeShape::myShape`.
    my_shape: Option<TopoShape>,
}

impl BRepOffsetApiMakeDraft {
    /// `BRepOffsetAPI_MakeDraft(Shape, Dir, Angle)` (`cxx:19-24`).
    pub fn new(shape: &TopoShape, dir: &GpDir, angle: f64) -> Result<Self, String> {
        // `cxx:23`: `myDraft(Shape, Dir, Angle)`.
        let my_draft = BrepFillDraft::new(shape, dir, angle)?;
        Ok(Self {
            my_draft,
            // `cxx:24`: `NotDone()`.
            my_done: false,
            my_shape: None,
        })
    }

    /// `SetOptions(Style = BRepBuilderAPI_RightCorner, AngleMin = 0.01,
    /// AngleMax = 3.0)` (`cxx:26-36`).
    pub fn set_options(&mut self, style: TransitionMode, angle_min: f64, angle_max: f64) {
        // `cxx:28-32`: only `RoundCorner` maps to `BRepFill_Round`.
        let mut style_out = TransitionStyle::Right;
        if style == TransitionMode::RoundCorner {
            style_out = TransitionStyle::Round;
        }
        // `cxx:33`.
        self.my_draft.set_options(style_out, angle_min, angle_max);
    }

    /// `SetDraft(IsInternal = false)` (`cxx:38-41`).
    pub fn set_draft(&mut self, is_internal: bool) {
        self.my_draft.set_draft(is_internal);
    }

    /// `Perform(const double LengthMax)` (`cxx:43-51`).
    pub fn perform(&mut self, length_max: f64) -> Result<(), String> {
        // `cxx:45`.
        self.my_draft.perform(length_max)?;
        // `cxx:46-50`.
        if self.my_draft.is_done() {
            self.my_done = true;
            self.my_shape = self.my_draft.shape().cloned();
        }
        Ok(())
    }

    /// `Perform(Surface, KeepInsideSurface = true)` (`cxx:53-62`).
    pub fn perform_on_surface(
        &mut self,
        surface: &Arc<dyn Surface>,
        keep_inside_surface: bool,
    ) -> Result<(), String> {
        // `cxx:55`.
        self.my_draft
            .perform_on_surface(surface, keep_inside_surface)?;
        // `cxx:56-61`.
        if self.my_draft.is_done() {
            self.my_done = true;
            self.my_shape = self.my_draft.shape().cloned();
        }
        Ok(())
    }

    /// `Perform(StopShape, KeepOutSide = true)` (`cxx:64-73`).
    pub fn perform_on_shape(
        &mut self,
        stop_shape: &TopoShape,
        keep_out_side: bool,
    ) -> Result<(), String> {
        // `cxx:66`.
        self.my_draft.perform_on_shape(stop_shape, keep_out_side)?;
        // `cxx:67-72`.
        if self.my_draft.is_done() {
            self.my_done = true;
            self.my_shape = self.my_draft.shape().cloned();
        }
        Ok(())
    }

    /// `TopoDS_Shell Shell() const` (`cxx:75-78`).
    pub fn shell(&self) -> Option<&Shell> {
        self.my_draft.shell()
    }

    /// `const NCollection_List<TopoDS_Shape>& Generated(S)` (`cxx:80-85`).
    pub fn generated(&mut self, s: &TopoShape) -> Result<&[TopoShape], String> {
        self.my_draft.generated(s)
    }

    /// `IsDone()` (`BRepBuilderAPI_Command`, inherited).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// `Shape()` (`BRepBuilderAPI_MakeShape`, inherited).
    pub fn shape(&self) -> Option<&TopoShape> {
        self.my_shape.as_ref()
    }

    /// The `BRepFill_Draft` object the reference hides in `myDraft`
    /// (`hxx:97`).
    pub fn draft(&self) -> &BrepFillDraft {
        &self.my_draft
    }
}
