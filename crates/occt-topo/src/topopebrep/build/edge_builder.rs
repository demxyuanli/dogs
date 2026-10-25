//! `TopOpeBRepBuild_EdgeBuilder` — builds the edges out of a `PaveSet`.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_EdgeBuilder.hxx:31-83`,
//! `.cxx:22-82` (OCCT 8.0.0).
//!
//! `InitEdgeBuilder` (`cxx:34-38`) is `InitAreaBuilder` and nothing else; the
//! iteration methods are the `AreaBuilder` ones under the edge/vertex names
//! (`cxx:42-74`), and `Vertex` / `Parameter` down-cast the current loop to
//! `TopOpeBRepBuild_Pave`.

use std::ops::{Deref, DerefMut};

use crate::shape::TopoShape;

use super::area1d_builder::Area1dBuilder;
use super::loop_classifier::LoopClassifier;
use super::loop_set::LoopSetSource;

/// `class TopOpeBRepBuild_EdgeBuilder` (`hxx:39`).
#[repr(transparent)]
pub struct EdgeBuilder {
    /// The `TopOpeBRepBuild_Area1dBuilder` state (`hxx:39`).
    pub base: Area1dBuilder,
}

impl EdgeBuilder {
    /// `TopOpeBRepBuild_EdgeBuilder()` (`cxx:24`).
    pub fn new() -> Self {
        Self {
            base: Area1dBuilder::new(),
        }
    }

    /// `TopOpeBRepBuild_EdgeBuilder(LS, LC, ForceClass)` (`cxx:28-32`).
    pub fn with_loop_set(
        ls: &mut dyn LoopSetSource,
        lc: &mut dyn LoopClassifier,
        force_class: bool,
    ) -> Self {
        let mut b = Self::new();
        b.init_edge_builder(ls, lc, force_class);
        b
    }

    /// `InitEdgeBuilder(LS, LC, ForceClass)` (`cxx:34-38`).
    pub fn init_edge_builder(
        &mut self,
        ls: &mut dyn LoopSetSource,
        lc: &mut dyn LoopClassifier,
        force_class: bool,
    ) {
        self.base.init_area_builder(ls, lc, force_class);
    }

    /// `InitEdge()` (`cxx:42-45`).
    pub fn init_edge(&mut self) -> i32 {
        self.base.init_area()
    }

    /// `MoreEdge()` (`cxx:49-52`).
    pub fn more_edge(&self) -> bool {
        self.base.more_area()
    }

    /// `NextEdge()` (`cxx:56-59`).
    pub fn next_edge(&mut self) {
        self.base.next_area();
    }

    /// `InitVertex()` (`cxx:63-66`).
    pub fn init_vertex(&mut self) -> i32 {
        self.base.init_loop()
    }

    /// `MoreVertex()` (`cxx:70-73`).
    pub fn more_vertex(&self) -> bool {
        self.base.more_loop()
    }

    /// `NextVertex()` (`cxx:77-80`).
    pub fn next_vertex(&mut self) {
        self.base.next_loop();
    }

    /// `Vertex()` (`cxx:84-90`): `PV->Vertex()` of the current `Pave`.
    ///
    /// OCCT hands the shape back by reference; the port clones it out of the
    /// [`crate::topopebrep::build::LoopHandle`] lock.
    pub fn vertex(&self) -> Option<TopoShape> {
        let l = self.base.loop_()?;
        let guard = l.read().expect("poisoned Loop lock");
        guard.as_pave().map(|p| p.vertex().clone())
    }

    /// `Parameter()` (`cxx:94-100`): `PV->Parameter()` of the current `Pave`.
    pub fn parameter(&self) -> f64 {
        let Some(l) = self.base.loop_() else {
            return 0.0;
        };
        let guard = l.read().expect("poisoned Loop lock");
        guard.as_pave().map(|p| p.parameter()).unwrap_or(0.0)
    }
}

impl Default for EdgeBuilder {
    fn default() -> Self {
        Self::new()
    }
}

impl Deref for EdgeBuilder {
    type Target = Area1dBuilder;

    fn deref(&self) -> &Area1dBuilder {
        &self.base
    }
}

impl DerefMut for EdgeBuilder {
    fn deref_mut(&mut self) -> &mut Area1dBuilder {
        &mut self.base
    }
}
