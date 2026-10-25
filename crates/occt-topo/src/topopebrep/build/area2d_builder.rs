//! `TopOpeBRepBuild_Area2dBuilder` — builds the 2d areas (faces) of a
//! `LoopSet` of wires and blocks of edges.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_Area2dBuilder.hxx:38-74`,
//! `.cxx:24-262` (OCCT 8.0.0).
//!
//! `InitAreaBuilder` (`cxx:41-262`) repeats `AreaBuilder::InitAreaBuilder`
//! (`AreaBuilder.cxx:124-330`) verbatim apart from `OCCT_DEBUG` output, hence
//! the delegation.

use std::ops::{Deref, DerefMut};

use super::area_builder::AreaBuilder;
use super::loop_classifier::LoopClassifier;
use super::loop_set::LoopSetSource;

/// `class TopOpeBRepBuild_Area2dBuilder` (`hxx:57`).
#[repr(transparent)]
pub struct Area2dBuilder {
    /// The `TopOpeBRepBuild_AreaBuilder` state (`hxx:57`).
    pub base: AreaBuilder,
}

impl Area2dBuilder {
    /// `TopOpeBRepBuild_Area2dBuilder()` (`cxx:26`).
    pub fn new() -> Self {
        Self {
            base: AreaBuilder::new(),
        }
    }

    /// `TopOpeBRepBuild_Area2dBuilder(LS, LC, ForceClass)` (`cxx:30-35`).
    pub fn with_loop_set(
        ls: &mut dyn LoopSetSource,
        lc: &mut dyn LoopClassifier,
        force_class: bool,
    ) -> Self {
        let mut b = Self::new();
        b.init_area_builder(ls, lc, force_class);
        b
    }

    /// `InitAreaBuilder(LS, LC, ForceClass)` (`cxx:41-262`).
    pub fn init_area_builder(
        &mut self,
        ls: &mut dyn LoopSetSource,
        lc: &mut dyn LoopClassifier,
        force_class: bool,
    ) {
        self.base.init_area_builder(ls, lc, force_class);
    }
}

impl Default for Area2dBuilder {
    fn default() -> Self {
        Self::new()
    }
}

impl Deref for Area2dBuilder {
    type Target = AreaBuilder;

    fn deref(&self) -> &AreaBuilder {
        &self.base
    }
}

impl DerefMut for Area2dBuilder {
    fn deref_mut(&mut self) -> &mut AreaBuilder {
        &mut self.base
    }
}
