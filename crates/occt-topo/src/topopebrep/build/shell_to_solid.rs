//! `TopOpeBRepBuild_ShellToSolid` — rebuild solids from a set of shells.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_ShellToSolid.cxx`
//! (OCCT 8.0.0, 61 lines) and `TopOpeBRepBuild_ShellToSolid.hxx`.
//!
//! `MakeSolids` (`cxx:43-61`) feeds the collected shells into a
//! `SolidBuilder` with `ForceClass = true` and hands the result to a *local*
//! `TopOpeBRepBuild_Builder::MakeSolids` — a fresh builder, not the caller's,
//! so nothing of the caller's split/merge bookkeeping is touched.

use crate::shape::TopoShape;
use crate::topopebrep::build::shell_face_set::ShellFaceSet;
use crate::topopebrep::build::solid_builder::SolidBuilder;
use crate::topopebrep::build::Builder;
use crate::topopebrep::ds::build_tool::BuildTool;

/// `class TopOpeBRepBuild_ShellToSolid` (`hxx`).
#[derive(Default)]
pub struct ShellToSolid {
    /// `myLSh` (`hxx`).
    my_l_sh: Vec<TopoShape>,
}

impl ShellToSolid {
    /// `TopOpeBRepBuild_ShellToSolid()` (`cxx:25`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `Init()` (`cxx:29-32`).
    pub fn init(&mut self) {
        self.my_l_sh.clear();
    }

    /// `AddShell(Sh)` (`cxx:36-39`).
    pub fn add_shell(&mut self, sh: &TopoShape) {
        self.my_l_sh.push(sh.clone());
    }

    /// `MakeSolids(So, LSo)` (`cxx:43-61`).
    pub fn make_solids(&self, so: &TopoShape, l_so: &mut Vec<TopoShape>) {
        l_so.clear();

        // `cxx:48-52`: the solid is the shape set, the shells are its shapes.
        let mut sfs = ShellFaceSet::with_shape(so);
        for sh in &self.my_l_sh {
            sfs.add_shape(sh);
        }

        // `cxx:54-56`.
        let force_class = true;
        let mut sb = SolidBuilder::with_shell_face_set(&mut sfs, force_class);

        // `cxx:58-60`: a local builder over a default build tool.
        let mut b = Builder::new(BuildTool::new());
        b.make_solids(&mut sb, l_so);
    }
}
