//! `BRepOffset_MakeOffset::MakeThickSolid` (`BRepOffset_MakeOffset.cxx:1083-1177`):
//! the thickening driver on top of [`BRepOffsetMakeOffset::make_offset_shape`].
//!
//! Control flow, line by line:
//!
//! * `MakeOffsetShape` first (`:1088`); a failed run keeps its error and
//!   answers `!myDone` (`:1090-1094`).
//! * With closing faces (`:1100-1169`) the result solid is glued from the
//!   initial shape and the reversed offset faces: `BRepTools_Quilt::Add` over
//!   `myShape` faces (`:1110-1114`), then over `myOffsetShape` faces reversed
//!   (`:1116-1129`). No offset face at all is `BRepOffset_UnknownError`
//!   (`:1137-1142`). `Glue.Shells()` replaces `myOffsetShape` (`:1144`), its
//!   shells go into `Res` and `Res.Closed(true)` (`:1145-1150`).
//! * The face-count validity test (`:1152-1168`): fewer faces than the input
//!   solid is `BRepOffset_UnknownError`, an equal count zeroes `myOffset`.
//! * A positive offset reverses the final shape (`:1171-1174`), then
//!   `myDone = true` (`:1176`).
//!
//! Port notes:
//!
//! * `TopoDS_Shape::Reversed()` is `oriented(o.reversed())`, matching
//!   `TopAbs::Reverse(Or)` (`TopAbs.hxx:122-136`), not `oriented(Reversed)`.
//! * `TopoDS_Shape::Reverse()` (`:1173`) mutates in place -> [`TopoShape::reverse`].
//! * `B.MakeSolid(Res)` is `TopoBuilder::make_solid(&[])` (`builder.rs:399`).
//! * `Message_ProgressRange` is not threaded; the port has no progress
//!   machinery, exactly like the rest of the `BRepOffset` modules.

use crate::abs::ShapeType;
use crate::brep_sweep::Quilt;
use crate::builder::TopoBuilder;
use crate::shape::{Shell, TopoShape};
use crate::topo_tools::explorer_shapes;

use super::{BRepOffsetMakeOffset, Error};

/// `S.Reversed()` (`TopoDS_Shape::Reversed`, `TopoDS_Shape.lxx:58-63`).
fn reversed_copy(s: &TopoShape) -> TopoShape {
    s.oriented(s.orientation().reversed())
}

impl BRepOffsetMakeOffset {
    /// `MakeThickSolid(theRange)` (`cxx:1083-1177`).
    pub fn make_thick_solid(&mut self) {
        // `cxx:1088`: the shell parallel to the shell, without caps.
        self.make_offset_shape();

        // `cxx:1090-1094`: keep the saved return code and `myDone`.
        if !self.my_done {
            return;
        }

        // `cxx:1100-1169`: a solid from the initial shell plus the offset
        // shell limited by the caps.
        if self.my_faces.extent() != 0 {
            let builder = TopoBuilder::new();
            let mut nb_f = self.my_faces.extent();

            // `cxx:1107`.
            let mut res = builder.make_solid(&[]);

            // `cxx:1109-1114`.
            let mut glue = Quilt::new();
            if let Some(shape) = self.my_shape.clone() {
                for face in explorer_shapes(&shape, ShapeType::Face) {
                    nb_f += 1;
                    glue.add(&face);
                }
            }

            // `cxx:1115-1135`.
            let mut ya_result = false;
            if let Some(offset_shape) = self.my_offset_shape.clone() {
                for face in explorer_shapes(&offset_shape, ShapeType::Face) {
                    ya_result = true;
                    glue.add(&reversed_copy(&face));
                }
            }

            // `cxx:1137-1142`.
            if !ya_result {
                self.my_done = false;
                self.my_error = Error::UnknownError;
                return;
            }

            // `cxx:1144`: `Shells()` answers a compound of the glued shells.
            let glued = glue.shells();
            // `cxx:1145-1148`.
            for shell in explorer_shapes(&glued, ShapeType::Shell) {
                builder.add_shell(&mut res, &Shell(shell));
            }
            // `cxx:1149-1150`.
            res.0.set_closed(true);
            self.my_offset_shape = Some(res.0.clone());

            // `cxx:1152-1168`: the validity test of the thick solid.
            let mut nb_of = 0usize;
            if let Some(offset_shape) = self.my_offset_shape.clone() {
                for _ in explorer_shapes(&offset_shape, ShapeType::Face) {
                    nb_of += 1;
                }
            }
            if nb_of < nb_f {
                self.my_done = false;
                self.my_error = Error::UnknownError;
                return;
            }
            if nb_of == nb_f {
                self.my_offset = 0.0;
            }
        }

        // `cxx:1171-1174`.
        if self.my_offset > 0.0 {
            if let Some(offset_shape) = self.my_offset_shape.as_mut() {
                offset_shape.reverse();
            }
        }

        // `cxx:1176`.
        self.my_done = true;
    }
}
