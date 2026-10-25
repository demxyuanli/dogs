//! `TopOpeBRepBuild_ShellFaceSet` — the face set of a solid: a bound is a
//! shell, a bound element is a face.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_ShellFaceSet.hxx:31-81`,
//! `.cxx:28-153` (OCCT 8.0.0).
//!
//! The `ShapeSet` stores the shells (shapes), the faces to start
//! reconstructions (start elements) and, for every edge, the faces incident
//! to it (`hxx:31-35`). `ShapeSet` is constructed on `TopAbs_EDGE`
//! (`cxx:29`, `:46`), so its sub-shapes are the edges shared by the faces.
//!
//! `AddShape` and `AddElement` forward to the base implementations
//! (`cxx:79-82`, `:93-96`); `AddStartElement` deliberately calls
//! `ProcessAddStartElement` directly (`cxx:88`), skipping the `CheckShape`
//! gate of `ShapeSet::AddStartElement`.

use crate::abs::ShapeType;
use crate::shape::{Solid, TopoShape};

use super::shape_set::{ShapeSet, ShapeSetSource};

/// `class TopOpeBRepBuild_ShellFaceSet` (`hxx:36`).
pub struct ShellFaceSet {
    /// The `TopOpeBRepBuild_ShapeSet` base (`TopAbs_EDGE`, `cxx:29`).
    base: ShapeSet,
    /// `mySolid` (`hxx:80`).
    my_solid: Option<Solid>,
}

impl ShellFaceSet {
    /// `TopOpeBRepBuild_ShellFaceSet()` (`cxx:28-34`).
    pub fn new() -> Self {
        Self {
            base: ShapeSet::new(ShapeType::Edge, false),
            my_solid: None,
        }
    }

    /// `TopOpeBRepBuild_ShellFaceSet(S, Addr)` (`cxx:38-75`).
    pub fn with_shape(s: &TopoShape) -> Self {
        Self {
            base: ShapeSet::new(ShapeType::Edge, false),
            my_solid: Solid::wrap(s.clone()),
        }
    }

    /// `Solid()` (`cxx:100-103`).
    pub fn solid(&self) -> Option<&Solid> {
        self.my_solid.as_ref()
    }

    /// `AddShape(S)` (`cxx:79-82`).
    pub fn add_shape(&mut self, s: &TopoShape) {
        ShapeSet::add_shape(&mut self.base, s);
    }

    /// `AddStartElement(S)` (`cxx:86-89`): `ProcessAddStartElement` directly.
    pub fn add_start_element(&mut self, s: &TopoShape) {
        self.base.process_add_start_element(s);
    }

    /// `AddElement(S)` (`cxx:93-96`).
    pub fn add_element(&mut self, s: &TopoShape) {
        ShapeSet::add_element(&mut self.base, s);
    }

    /// `DumpSS()` (`cxx:107-112`) — debug output only.
    pub fn dump_ss(&self) {}

    /// `SName(S, sb, sa)` (`cxx:116-122`).
    pub fn s_name(&self, _s: &TopoShape, sb: &str, _sa: &str) -> String {
        sb.to_string()
    }

    /// `SNameori(S, sb, sa)` (`cxx:126-132`).
    pub fn s_name_ori(&self, _s: &TopoShape, sb: &str, _sa: &str) -> String {
        sb.to_string()
    }

    /// `SName(L, sb, sa)` (`cxx:136-142`).
    pub fn s_name_list(&self, _l: &[TopoShape], _sb: &str, _sa: &str) -> String {
        String::new()
    }

    /// `SNameori(L, sb, sa)` (`cxx:146-153`).
    pub fn s_name_ori_list(&self, _l: &[TopoShape], _sb: &str, _sa: &str) -> String {
        String::new()
    }
}

impl Default for ShellFaceSet {
    fn default() -> Self {
        Self::new()
    }
}

impl ShapeSetSource for ShellFaceSet {
    fn base(&self) -> &ShapeSet {
        &self.base
    }

    fn base_mut(&mut self) -> &mut ShapeSet {
        &mut self.base
    }

    fn add_shape(&mut self, s: &TopoShape) {
        ShellFaceSet::add_shape(self, s);
    }

    fn add_start_element(&mut self, s: &TopoShape) {
        ShellFaceSet::add_start_element(self, s);
    }

    fn add_element(&mut self, s: &TopoShape) {
        ShellFaceSet::add_element(self, s);
    }

    fn init_neighbours(&mut self, s: &TopoShape) {
        self.base.init_neighbours(s);
    }

    fn find_neighbours(&mut self) {
        self.base.find_neighbours();
    }

    fn make_neighbours_list(&mut self, e: Option<&TopoShape>, v: &TopoShape) -> Vec<TopoShape> {
        ShapeSet::make_neighbours_list(&self.base, e, v)
    }
}
