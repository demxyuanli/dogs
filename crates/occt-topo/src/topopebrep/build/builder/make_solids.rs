//! `TopOpeBRepBuild_makesolids.cxx` — solids out of a `ShellFaceSet`.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_makesolids.cxx`
//! (OCCT 8.0.0):
//!
//! | OCCT | lines |
//! |---|---|
//! | `GSFSMakeSolids(SOF, SFS, LOSO)` | `35-54` |
//! | `GSOBUMakeSolids(SOF, SOBU, LOSO)` | `58-172` |
//!
//! UNPORTED: the `OCCT_DEBUG` blocks (`cxx:39-46`, `:64-75`) — trace and
//! `debgsobu`/`GdumpSOBU` dumps only.

use crate::abs::{Orientation, ShapeType};
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, TopoShape};
use crate::topo_tools::{explorer_shapes, map_shapes_and_ancestors, shape_map_key};
use crate::topopebrep::build::shell_face_set::ShellFaceSet;
use crate::topopebrep::build::solid_builder::SolidBuilder;

use super::Builder;

impl Builder {
    /// `GSFSMakeSolids(SOF, SFS, LOSO)` (`cxx:35-54`): `ForceClass` is
    /// `true` (`cxx:49`).
    pub fn gsfs_make_solids(
        &mut self,
        sof: &TopoShape,
        sfs: &mut ShellFaceSet,
        loso: &mut Vec<TopoShape>,
    ) {
        let mut sobu = SolidBuilder::with_shell_face_set(sfs, true);
        self.gsobu_make_solids(sof, &mut sobu, loso);
    }

    /// `GSOBUMakeSolids(SOF, SOBU, LOSO)` (`cxx:58-172`).
    pub fn gsobu_make_solids(
        &mut self,
        sof: &TopoShape,
        sobu: &mut SolidBuilder,
        loso: &mut Vec<TopoShape>,
    ) {
        let mut new_solid = TopoShape::new(ShapeType::Solid);
        let mut new_shell = TopoShape::new(ShapeType::Shell);
        // `cxx:79-80`: `nfa`/`nsh` are re-assigned per solid/shell; the value
        // seen by `cxx:150` is the last `InitShell`/`InitFace` extent.
        let mut nfa: i32 = 0;
        let mut nsh: i32;
        sobu.init_solid();
        while sobu.more_solid() {
            self.my_build_tool.make_solid(&mut new_solid);
            nsh = sobu.init_shell();
            while sobu.more_shell() {
                let isold = sobu.is_old_shell();
                if isold {
                    new_shell = sobu
                        .old_shell()
                        .unwrap_or_else(|| TopoShape::new(ShapeType::Shell));
                } else {
                    self.my_build_tool.make_shell(&mut new_shell);
                    nfa = sobu.init_face();
                    while sobu.more_face() {
                        let f = sobu.face().cloned();
                        if let Some(f) = f {
                            self.my_build_tool.add_shell_face(&mut new_shell, &f);
                        }
                        sobu.next_face();
                    }
                }

                // `cxx:106-138`: closed character of the new shell.
                if !isold {
                    let mut closed = true;
                    let edgemap =
                        map_shapes_and_ancestors(&new_shell, ShapeType::Edge, ShapeType::Face);
                    for e in explorer_shapes(&new_shell, ShapeType::Edge) {
                        if matches!(
                            e.orientation(),
                            Orientation::Internal | Orientation::External
                        ) {
                            continue;
                        }
                        if BRepTool::is_degenerated(&Edge(e.clone())) {
                            continue;
                        }
                        let nbf = edgemap.get(&shape_map_key(&e)).map_or(0, Vec::len);
                        if nbf < 2 {
                            closed = false;
                            break;
                        }
                    }
                    self.my_build_tool.closed(&new_shell, closed);
                }

                self.my_build_tool.add_solid_shell(&mut new_solid, &new_shell);
                sobu.next_shell();
            }

            let has_vertex = !explorer_shapes(&new_solid, ShapeType::Vertex).is_empty();
            if has_vertex {
                sobu.next_solid();
                continue;
            }

            // `cxx:148-167`: a single shell made of a single face is kept only
            // when it carries a non-degenerated edge.
            let new_solid_ok;
            if nsh == 1 && nfa == 1 {
                let mut has_non_degenerated = false;
                for e in explorer_shapes(&new_solid, ShapeType::Edge) {
                    if !BRepTool::is_degenerated(&Edge(e.clone())) {
                        has_non_degenerated = true;
                        break;
                    }
                }
                new_solid_ok = has_non_degenerated;
                if !new_solid_ok {
                    sobu.next_solid();
                    continue;
                }
            }

            // `cxx:169-173`.
            let mut new_solid_los: Vec<TopoShape> = Vec::new();
            self.regularize_solid(sof, &new_solid, &mut new_solid_los);
            loso.extend(new_solid_los);
            sobu.next_solid();
        }
    }
}
