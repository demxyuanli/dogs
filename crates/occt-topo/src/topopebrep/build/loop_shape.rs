//! `TopOpeBRepBuild_Loop` — a loop of the area builders.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_Loop.hxx:27-54`,
//! `.cxx:22-64`.
//!
//! OCCT spells the `TopOpeBRepBuild_Loop` / `TopOpeBRepBuild_Pave` hierarchy with
//! `Standard_Transient` RTTI and keeps `NCollection_List<handle<Loop>>`, i.e. one
//! homogeneous list of base handles. The Rust translation keeps one [`Loop`]
//! value with a [`LoopKind`] payload holding the [`PaveData`] subclass fields,
//! so a list stays homogeneous too: `occ::down_cast<TopOpeBRepBuild_Pave>` is
//! [`Loop::as_pave`] / [`Loop::as_pave_mut`], and `DynamicType()` is
//! [`Loop::dynamic_type`].

use std::sync::{Arc, RwLock};

use crate::shape::TopoShape;

use super::block_iterator::BlockIterator;
use super::pave::PaveData;

/// `occ::handle<TopOpeBRepBuild_Loop>`.
///
/// OCCT hands out the same object by reference (`Pave::ChangeVertex()` mutates
/// what every holder sees), so the Rust handle is a shared, lockable cell.
pub type LoopHandle = Arc<RwLock<Loop>>;

/// The subclass payload of a [`Loop`]. [`LoopKind::Base`] is
/// `TopOpeBRepBuild_Loop` itself; the subclasses of the package get one variant
/// each as they are ported.
#[derive(Debug, Clone, Default)]
pub enum LoopKind {
    /// `TopOpeBRepBuild_Loop`
    #[default]
    Base,
    /// `TopOpeBRepBuild_Pave`
    Pave(PaveData),
}

/// `class TopOpeBRepBuild_Loop` (`hxx:33`).
///
/// A Loop is an existing shape (Shell, Wire) or a set of shapes (Faces, Edges)
/// which are connex. A set of connex shapes is represented by a
/// [`BlockIterator`].
#[derive(Debug, Clone, Default)]
pub struct Loop {
    /// `myIsShape` (`hxx:49`).
    pub(super) my_is_shape: bool,
    /// `myShape` (`hxx:50`), left null by the `BlockIterator` constructor.
    pub(super) my_shape: Option<TopoShape>,
    /// `myBlockIterator` (`hxx:51`).
    pub(super) my_block_iterator: BlockIterator,
    /// The subclass payload.
    pub(super) kind: LoopKind,
}

impl Loop {
    /// `TopOpeBRepBuild_Loop(const TopoDS_Shape& S)` (`cxx:22-28`):
    /// `myIsShape = true`, `myBlockIterator(0, 0)`.
    pub fn new(s: &TopoShape) -> Self {
        Self {
            my_is_shape: true,
            my_shape: Some(s.clone()),
            my_block_iterator: BlockIterator::with_bounds(0, 0),
            kind: LoopKind::Base,
        }
    }

    /// `TopOpeBRepBuild_Loop(const BlockIterator& BI)` (`cxx:32-36`):
    /// `myIsShape = false`, `myShape` stays null.
    pub fn with_block_iterator(bi: BlockIterator) -> Self {
        Self {
            my_is_shape: false,
            my_shape: None,
            my_block_iterator: bi,
            kind: LoopKind::Base,
        }
    }

    /// `IsShape()` (`cxx:40-43`).
    ///
    /// `TopOpeBRepBuild_Pave` declares its own `myIsShape`
    /// (`Pave.hxx:74`) which hides the base field, and `Pave::IsShape()`
    /// (`Pave.cxx:76-79`) returns that one, so the answer is delegated to the
    /// payload when the loop is a `Pave`.
    pub fn is_shape(&self) -> bool {
        match &self.kind {
            LoopKind::Base => self.my_is_shape,
            LoopKind::Pave(p) => p.is_shape,
        }
    }

    /// `Shape()` (`cxx:47-50`). `Pave::Shape()` (`Pave.cxx:81-84`) returns its
    /// vertex instead of the base shape.
    pub fn shape(&self) -> Option<&TopoShape> {
        match &self.kind {
            LoopKind::Base => self.my_shape.as_ref(),
            LoopKind::Pave(p) => Some(&p.vertex),
        }
    }

    /// `BlockIterator()` (`cxx:54-57`).
    pub fn block_iterator(&self) -> &BlockIterator {
        &self.my_block_iterator
    }

    /// `Dump()` (`cxx:63-64`) — empty in `TopOpeBRepBuild_Loop`.
    pub fn dump(&self) {}

    /// `DynamicType()` for the RTTI queries of the package.
    pub fn dynamic_type(&self) -> &'static str {
        match &self.kind {
            LoopKind::Base => "TopOpeBRepBuild_Loop",
            LoopKind::Pave(_) => "TopOpeBRepBuild_Pave",
        }
    }

    /// `occ::down_cast<TopOpeBRepBuild_Pave>(L)`.
    pub fn as_pave(&self) -> Option<&PaveData> {
        match &self.kind {
            LoopKind::Pave(p) => Some(p),
            LoopKind::Base => None,
        }
    }

    /// `occ::down_cast<TopOpeBRepBuild_Pave>(L)`, mutable.
    pub fn as_pave_mut(&mut self) -> Option<&mut PaveData> {
        match &mut self.kind {
            LoopKind::Pave(p) => Some(p),
            LoopKind::Base => None,
        }
    }

    /// The subclass payload.
    pub fn kind(&self) -> &LoopKind {
        &self.kind
    }
}
