//! `BRepOffset_MakeOffset::MakeShells` (`BRepOffset_MakeOffset.cxx:3691-3792`)
//! and `SelectShells` (`:3896-3927`).
//!
//! `MakeShells` (`:3691-3792`) collects the last images of every
//! `myImageOffset` root (`:3703-3729`); for `myThickening` the root is reversed
//! first (`:3717-3720`) and the faces of `myShape` plus `myWalls` are appended
//! (`:3731-3747`). An empty list leaves `myOffsetShape` null (`:3749-3752`).
//! The Planar Complete arm calls `BuildShellsCompleteInter` (`:3754-3766`);
//! every other case glues with `BRepTools_Quilt` (`:3768-3777`) and then sets
//! the `Closed` flag of every shell whose `BRep_Tool::IsClosed` says so
//! (`:3780-3788`).
//!
//! `SelectShells` (`:3896-3927`) marks the edges with fewer than two ancestors
//! whose type is `ChFiDS_FreeBound`, returns untouched when there are free
//! borders but no caps (`:3918-3922`), and otherwise replaces `myOffsetShape`
//! with `BRepOffset_Tool::Deboucle3D(myOffsetShape, FreeEdges)` (`:3924`).
//!
//! The planar `Complete` shell builder `BuildShellsCompleteInter`
//! (`:5064-5227`) and its `UpdateHistory` (`:4649-4670`) live in
//! `shells_complete_inter.rs`.
//!
//! UNPORTED:
//!
//! * `Message_ProgressScope` / `!More()` -> `BRepOffset_UserBreak`
//!   (`:3710-3714`).
//! * the `OCCT_DEBUG` blocks (`:3693-3701`).
//!
use std::collections::HashSet;

use occt_core::kernel::geomabs::JoinType;

use crate::abs::ShapeType;
use crate::brep_offset::BRepOffsetTool;
use crate::brep_sweep::Quilt;
use crate::brep_tool::BRepTool;
use crate::chfi_ds::TypeOfConcavity;
use crate::shape::{Edge, TopoShape};
use crate::topo_tools::{explorer_shapes, shape_map_key};

use super::shape_prep::is_solid;
use super::shells_complete_inter::build_shells_complete_inter;
use super::BRepOffsetMakeOffset;

impl BRepOffsetMakeOffset {
    /// `MakeShells(theRange)` (`cxx:3691-3792`).
    pub fn make_shells(&mut self) {
        // `cxx:3703-3729`.
        let mut a_lsf: Vec<TopoShape> = Vec::new();
        let roots = self.my_image_offset.roots().to_vec();
        for mut a_f in roots {
            // `cxx:3717-3720`: offset faces change their orientation.
            if self.my_thickening {
                a_f.reverse();
            }
            let mut image: Vec<TopoShape> = Vec::new();
            self.my_image_offset.last_image(&a_f, &mut image);
            a_lsf.extend(image);
        }

        // `cxx:3731-3747`.
        if self.my_thickening {
            if let Some(shape) = self.my_shape.clone() {
                a_lsf.extend(explorer_shapes(&shape, ShapeType::Face));
            }
            a_lsf.extend(self.my_walls.iter().cloned());
        }

        // `cxx:3749-3752`.
        if a_lsf.is_empty() {
            return;
        }

        // `cxx:3753-3775`: the planar `Complete` arm, then the quilt fallback
        // whenever it answers `false`.
        let mut b_done = false;
        let planar_complete = self.my_join == JoinType::Intersection
            && self.my_inter
            && !self.my_thickening
            && self.my_faces.extent() == 0
            && self.my_shape.as_ref().is_some_and(is_solid)
            && self.my_is_planar;
        if planar_complete {
            if let Some(a_shells) =
                build_shells_complete_inter(&a_lsf, &mut self.my_image_offset)
            {
                self.my_offset_shape = Some(a_shells);
                b_done = true;
            }
        }

        // `cxx:3766-3775`: `BRepTools_Quilt`.
        if !b_done {
            let mut glue = Quilt::new();
            for a in &a_lsf {
                glue.add(a);
            }
            self.my_offset_shape = Some(glue.shells());
        }

        // `cxx:3780-3788`: set the correct value of the closed flag.
        let Some(shape) = self.my_offset_shape.clone() else {
            return;
        };
        for shell in explorer_shapes(&shape, ShapeType::Shell) {
            if !shell.is_closed() && BRepTool::is_closed(&shell) {
                shell.set_closed(true);
            }
        }
    }

    /// `SelectShells()` (`cxx:3896-3927`).
    pub fn select_shells(&mut self) {
        // `cxx:3898-3916`.
        let mut free_edges: HashSet<crate::topo_tools::ShapeMapKey> = HashSet::new();
        for e_shape in explorer_shapes(&self.my_face_comp.0, ShapeType::Edge) {
            let Some(e) = Edge::wrap(e_shape) else {
                continue;
            };
            let la = self.my_analyse.ancestors(&e.0);
            if la.len() < 2 {
                let first = self.my_analyse.type_of(&e);
                if first
                    .first()
                    .is_some_and(|i| i.kind() == TypeOfConcavity::FreeBound)
                {
                    free_edges.insert(shape_map_key(&e.0));
                }
            }
        }

        // `cxx:3918-3922`.
        if !free_edges.is_empty() && self.my_faces.extent() == 0 {
            return;
        }

        // `cxx:3924`.
        if let Some(shape) = self.my_offset_shape.clone() {
            self.my_offset_shape = BRepOffsetTool::deboucle_3d(&shape, &free_edges);
        }
    }
}
