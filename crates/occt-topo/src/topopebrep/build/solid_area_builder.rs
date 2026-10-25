//! `TopOpeBRepBuild_SolidAreaBuilder` — builds solids from shells and blocks of
//! faces.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_SolidAreaBuilder.hxx:38-72`,
//! `.cxx:22-37` (OCCT 8.0.0).
//!
//! `InitSolidAreaBuilder` (`cxx:32-37`) is `InitAreaBuilder` and nothing else,
//! so the class only carries the `Area3dBuilder` state.

use std::ops::{Deref, DerefMut};

use super::area3d_builder::Area3dBuilder;
use super::loop_classifier::LoopClassifier;
use super::loop_set::LoopSetSource;

/// `class TopOpeBRepBuild_SolidAreaBuilder` (`hxx:57`).
#[repr(transparent)]
pub struct SolidAreaBuilder {
    /// The `TopOpeBRepBuild_Area3dBuilder` state (`hxx:57`).
    pub base: Area3dBuilder,
}

impl SolidAreaBuilder {
    /// `TopOpeBRepBuild_SolidAreaBuilder()` (`cxx:24`).
    pub fn new() -> Self {
        Self {
            base: Area3dBuilder::new(),
        }
    }

    /// `TopOpeBRepBuild_SolidAreaBuilder(LS, LC, ForceClass)` (`cxx:28-31`).
    pub fn with_loop_set(
        ls: &mut dyn LoopSetSource,
        lc: &mut dyn LoopClassifier,
        force_class: bool,
    ) -> Self {
        let mut b = Self::new();
        b.init_solid_area_builder(ls, lc, force_class);
        b
    }

    /// `InitSolidAreaBuilder(LS, LC, ForceClass)` (`cxx:32-37`).
    pub fn init_solid_area_builder(
        &mut self,
        ls: &mut dyn LoopSetSource,
        lc: &mut dyn LoopClassifier,
        force_class: bool,
    ) {
        self.base.init_area_builder(ls, lc, force_class);
    }
}

impl Default for SolidAreaBuilder {
    fn default() -> Self {
        Self::new()
    }
}

impl Deref for SolidAreaBuilder {
    type Target = Area3dBuilder;

    fn deref(&self) -> &Area3dBuilder {
        &self.base
    }
}

impl DerefMut for SolidAreaBuilder {
    fn deref_mut(&mut self) -> &mut Area3dBuilder {
        &mut self.base
    }
}
