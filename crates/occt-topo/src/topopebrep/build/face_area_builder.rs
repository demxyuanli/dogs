//! `TopOpeBRepBuild_FaceAreaBuilder` — builds faces from wires and blocks of
//! edges.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_FaceAreaBuilder.hxx:38-72`,
//! `.cxx:22-36` (OCCT 8.0.0).
//!
//! `InitFaceAreaBuilder` (`cxx:32-36`) is `InitAreaBuilder` and nothing else, so
//! the class only carries the `Area2dBuilder` state.

use std::ops::{Deref, DerefMut};

use super::area2d_builder::Area2dBuilder;
use super::loop_classifier::LoopClassifier;
use super::loop_set::LoopSetSource;

/// `class TopOpeBRepBuild_FaceAreaBuilder` (`hxx:57`).
#[repr(transparent)]
pub struct FaceAreaBuilder {
    /// The `TopOpeBRepBuild_Area2dBuilder` state (`hxx:57`).
    pub base: Area2dBuilder,
}

impl FaceAreaBuilder {
    /// `TopOpeBRepBuild_FaceAreaBuilder()` (`cxx:24`).
    pub fn new() -> Self {
        Self {
            base: Area2dBuilder::new(),
        }
    }

    /// `TopOpeBRepBuild_FaceAreaBuilder(LS, LC, ForceClass)` (`cxx:28-31`).
    pub fn with_loop_set(
        ls: &mut dyn LoopSetSource,
        lc: &mut dyn LoopClassifier,
        force_class: bool,
    ) -> Self {
        let mut b = Self::new();
        b.init_face_area_builder(ls, lc, force_class);
        b
    }

    /// `InitFaceAreaBuilder(LS, LC, ForceClass)` (`cxx:32-36`).
    pub fn init_face_area_builder(
        &mut self,
        ls: &mut dyn LoopSetSource,
        lc: &mut dyn LoopClassifier,
        force_class: bool,
    ) {
        self.base.init_area_builder(ls, lc, force_class);
    }
}

impl Default for FaceAreaBuilder {
    fn default() -> Self {
        Self::new()
    }
}

impl Deref for FaceAreaBuilder {
    type Target = Area2dBuilder;

    fn deref(&self) -> &Area2dBuilder {
        &self.base
    }
}

impl DerefMut for FaceAreaBuilder {
    fn deref_mut(&mut self) -> &mut Area2dBuilder {
        &mut self.base
    }
}
