//! `TopOpeBRepBuild_Area1dBuilder` — builds the 1d areas (vertices/edges) of a
//! `PaveSet`.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_Area1dBuilder.hxx:38-75`,
//! `.cxx:52-413` (OCCT 8.0.0).
//!
//! `InitAreaBuilder` (`cxx:67-298`) repeats `AreaBuilder::InitAreaBuilder`
//! (`AreaBuilder.cxx:124-330`) verbatim apart from `OCCT_DEBUG` output, and
//! `ADD_Loop_TO_LISTOFLoop` / `REM_Loop_FROM_LISTOFLoop` /
//! `ADD_LISTOFLoop_TO_LISTOFLoop` (`cxx:302-413`) only `Append` / `Remove` like
//! the base ones (`AreaBuilder.cxx:401-430`); the class therefore holds the
//! base state and delegates.

use std::ops::{Deref, DerefMut};

use super::area_builder::AreaBuilder;
use super::loop_classifier::LoopClassifier;
use super::loop_set::LoopSetSource;

/// `class TopOpeBRepBuild_Area1dBuilder` (`hxx:58`).
#[repr(transparent)]
pub struct Area1dBuilder {
    /// The `TopOpeBRepBuild_AreaBuilder` state (`hxx:58`).
    pub base: AreaBuilder,
}

impl Area1dBuilder {
    /// `TopOpeBRepBuild_Area1dBuilder()` (`cxx:56`).
    pub fn new() -> Self {
        Self {
            base: AreaBuilder::new(),
        }
    }

    /// `TopOpeBRepBuild_Area1dBuilder(LS, LC, ForceClass)` (`cxx:62-65`).
    pub fn with_loop_set(
        ls: &mut dyn LoopSetSource,
        lc: &mut dyn LoopClassifier,
        force_class: bool,
    ) -> Self {
        let mut b = Self::new();
        b.init_area_builder(ls, lc, force_class);
        b
    }

    /// `InitAreaBuilder(LS, LC, ForceClass)` (`cxx:67-298`).
    pub fn init_area_builder(
        &mut self,
        ls: &mut dyn LoopSetSource,
        lc: &mut dyn LoopClassifier,
        force_class: bool,
    ) {
        self.base.init_area_builder(ls, lc, force_class);
    }
}

impl Default for Area1dBuilder {
    fn default() -> Self {
        Self::new()
    }
}

impl Deref for Area1dBuilder {
    type Target = AreaBuilder;

    fn deref(&self) -> &AreaBuilder {
        &self.base
    }
}

impl DerefMut for Area1dBuilder {
    fn deref_mut(&mut self) -> &mut AreaBuilder {
        &mut self.base
    }
}
