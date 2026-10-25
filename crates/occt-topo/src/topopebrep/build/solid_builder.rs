//! `TopOpeBRepBuild_SolidBuilder` — rebuilds the solids of a `ShellFaceSet`.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_SolidBuilder.hxx:31-79`,
//! `.cxx:22-147` (OCCT 8.0.0).
//!
//! `InitSolidBuilder` (`cxx:34-42`) turns the `ShellFaceSet` into a `LoopSet` of
//! the old shells plus the blocks of faces (`MakeLoops`, `cxx:129-147`), then
//! runs the `SolidAreaBuilder` with the `ShellFaceClassifier`.

use std::sync::{Arc, RwLock};

use crate::shape::TopoShape;

use super::block_builder::BlockBuilder;
use super::block_iterator::BlockIterator;
use super::loop_set::{LoopSet, LoopSetSource};
use super::loop_shape::{Loop, LoopHandle};
use super::shape_set::ShapeSetSource;
use super::shell_face_classifier::ShellFaceClassifier;
use super::shell_face_set::ShellFaceSet;
use super::solid_area_builder::SolidAreaBuilder;

/// `class TopOpeBRepBuild_SolidBuilder` (`hxx:36`).
pub struct SolidBuilder {
    /// `myLoopSet` (`hxx:74`).
    my_loop_set: LoopSet,
    /// `myBlockIterator` (`hxx:75`).
    my_block_iterator: BlockIterator,
    /// `myBlockBuilder` (`hxx:76`).
    my_block_builder: BlockBuilder,
    /// `mySolidAreaBuilder` (`hxx:77`).
    my_solid_area_builder: SolidAreaBuilder,
}

impl SolidBuilder {
    /// `TopOpeBRepBuild_SolidBuilder()` (`cxx:26`).
    pub fn new() -> Self {
        Self {
            my_loop_set: LoopSet::new(),
            my_block_iterator: BlockIterator::new(),
            my_block_builder: BlockBuilder::new(),
            my_solid_area_builder: SolidAreaBuilder::new(),
        }
    }

    /// `TopOpeBRepBuild_SolidBuilder(SFS, ForceClass)` (`cxx:30-34`).
    pub fn with_shell_face_set(sfs: &mut ShellFaceSet, force_class: bool) -> Self {
        let mut b = Self::new();
        b.init_solid_builder(sfs, force_class);
        b
    }

    /// `InitSolidBuilder(SFS, ForceClass)` (`cxx:36-42`).
    pub fn init_solid_builder(&mut self, sfs: &mut ShellFaceSet, force_class: bool) {
        self.make_loops(sfs);
        let mut sfc = ShellFaceClassifier::new(&self.my_block_builder);
        self.my_solid_area_builder
            .init_solid_area_builder(&mut self.my_loop_set, &mut sfc, force_class);
    }

    /// `MakeLoops(SS)` (`cxx:129-147`).
    fn make_loops(&mut self, ss: &mut dyn ShapeSetSource) {
        // Build blocks on elements of SS
        self.my_block_builder.make_block(ss);

        let mut loops: Vec<LoopHandle> = Vec::new();

        // Add shapes of SS as shape loops
        ss.init_shapes();
        while ss.more_shapes() {
            let Some(s) = ss.shape().cloned() else {
                break;
            };
            loops.push(Arc::new(RwLock::new(Loop::new(&s))));
            ss.next_shape();
        }

        // Add blocks of BB as block loops
        self.my_block_builder.init_block();
        while self.my_block_builder.more_block() {
            let bi = self.my_block_builder.block_iterator();
            loops.push(Arc::new(RwLock::new(Loop::with_block_iterator(bi))));
            self.my_block_builder.next_block();
        }

        let ll = self.my_loop_set.change_list_of_loop();
        ll.clear();
        ll.extend(loops);
    }

    /// `InitSolid()` (`cxx:46-50`).
    pub fn init_solid(&mut self) -> i32 {
        self.my_solid_area_builder.init_area()
    }

    /// `MoreSolid()` (`cxx:54-58`).
    pub fn more_solid(&self) -> bool {
        self.my_solid_area_builder.more_area()
    }

    /// `NextSolid()` (`cxx:62-65`).
    pub fn next_solid(&mut self) {
        self.my_solid_area_builder.next_area();
    }

    /// `InitShell()` (`cxx:69-73`).
    pub fn init_shell(&mut self) -> i32 {
        self.my_solid_area_builder.init_loop()
    }

    /// `MoreShell()` (`cxx:77-81`).
    pub fn more_shell(&self) -> bool {
        self.my_solid_area_builder.more_loop()
    }

    /// `NextShell()` (`cxx:85-88`).
    pub fn next_shell(&mut self) {
        self.my_solid_area_builder.next_loop();
    }

    /// `IsOldShell()` (`cxx:92-96`).
    pub fn is_old_shell(&self) -> bool {
        match self.current_loop() {
            Some(l) => l.read().expect("poisoned Loop lock").is_shape(),
            None => false,
        }
    }

    /// `OldShell()` (`cxx:100-105`).
    ///
    /// OCCT returns a reference into the loop; the port clones it out of the
    /// [`LoopHandle`] lock.
    pub fn old_shell(&self) -> Option<TopoShape> {
        self.current_loop()
            .and_then(|l| l.read().expect("poisoned Loop lock").shape().cloned())
    }

    /// `InitFace()` (`cxx:109-121`).
    ///
    /// `Standard_DomainError("TopOpeBRepBuild_SolidBuilder:InitFace")` for a
    /// shell loop becomes the port's `panic!`.
    pub fn init_face(&mut self) -> i32 {
        let Some(l) = self.current_loop().cloned() else {
            panic!("TopOpeBRepBuild_SolidBuilder:InitFace");
        };
        if l.read().expect("poisoned Loop lock").is_shape() {
            panic!("TopOpeBRepBuild_SolidBuilder:InitFace");
        }
        let bi = *l.read().expect("poisoned Loop lock").block_iterator();
        self.my_block_iterator = bi;
        self.my_block_iterator.initialize();
        self.my_block_iterator.extent()
    }

    /// `MoreFace()` (`cxx:125-129`).
    pub fn more_face(&self) -> bool {
        self.my_block_iterator.more()
    }

    /// `NextFace()` (`cxx:133-136`).
    pub fn next_face(&mut self) {
        self.my_block_iterator.next();
    }

    /// `Face()` (`cxx:140-147`): `myBlockBuilder.Element(myBlockIterator)`.
    pub fn face(&self) -> Option<&TopoShape> {
        self.my_block_builder
            .element_of_iterator(&self.my_block_iterator)
    }

    /// The current loop of the area under `mySolidAreaBuilder`.
    fn current_loop(&self) -> Option<&LoopHandle> {
        self.my_solid_area_builder.loop_()
    }
}

impl Default for SolidBuilder {
    fn default() -> Self {
        Self::new()
    }
}
