//! `BRepOffset_MakeOffset::RemoveInternalEdges`
//! (`BRepOffset_MakeOffset.cxx:4512-4574`) and the file-local helper
//! `RemoveShapes` (`:4628-4642`).
//!
//! `RemoveInternalEdges` drops the edges the algorithm itself created twice:
//! an edge whose orientation is `TopAbs_INTERNAL` **and** that belongs to a
//! single face only (`:4547-4553`) is removed from its wire; a wire whose every
//! edge is like that is removed from its face (`:4524-4568`). The face passed
//! to `RemoveShapes` is a cast of `TopExp_Explorer::Current()`
//! (`:4525`, `*(TopoDS_Face*)&aExpF.Current()`), i.e. a copy sharing the
//! `TShape`, which is why the removal reaches the face stored inside
//! `myOffsetShape` - the port relies on the same sharing and passes a
//! `TopoShape` clone to `TopoBuilder::remove`.
//!
//! Port notes:
//!
//! * the file-local `UpdateHistory` (`:4649-4670`) and
//!   `BuildShellsCompleteInter` (`_1.cxx:5064-5229`) called from it
//!   (`_1.cxx:5084`) are `shells_complete_inter.rs`.
//!
//! * `aDMELF.FindFromKey(aE)` (`:4547`) throws `Standard_NoSuchObject` when the
//!   edge is absent; the map is built over the very shape being walked, so the
//!   port reads the entry as an empty list instead.
//! * `myOffsetShape` is null when `MakeOffsetShape` bailed out early
//!   (`:1003-1007`); OCCT explores the null shape and finds nothing, the port
//!   keeps `Option` and returns.

use crate::abs::{Orientation, ShapeType};
use crate::builder::TopoBuilder;
use crate::shape::TopoShape;
use crate::topo_tools::{explorer_shapes, map_shapes_and_ancestors, shape_map_key};

use super::BRepOffsetMakeOffset;

/// `RemoveShapes(theS, theLS)` (`cxx:4628-4642`): `theS` is made free, the
/// listed children are removed, then the previous free flag is restored.
fn remove_shapes(the_s: &mut TopoShape, the_ls: &[TopoShape]) {
    let builder = TopoBuilder::new();
    // `cxx:4631-4632`.
    let b_free = the_s.free();
    the_s.set_free(true);
    // `cxx:4634-4640`.
    for a_si in the_ls {
        builder.remove(the_s, a_si);
    }
    // `cxx:4642`.
    the_s.set_free(b_free);
}

impl BRepOffsetMakeOffset {
    /// `RemoveInternalEdges()` (`cxx:4512-4574`).
    pub fn remove_internal_edges(&mut self) {
        let Some(offset_shape) = self.my_offset_shape.clone() else {
            return;
        };

        // `cxx:4518-4519`: every edge of the result mapped to its faces.
        let ancestors = map_shapes_and_ancestors(&offset_shape, ShapeType::Edge, ShapeType::Face);

        // `cxx:4521-4572`.
        for a_f in explorer_shapes(&offset_shape, ShapeType::Face) {
            let mut a_liw: Vec<TopoShape> = Vec::new();

            // `cxx:4528-4567`: one wire at a time.
            for a_w in explorer_shapes(&a_f, ShapeType::Wire) {
                let mut b_remove_wire = true;
                let mut a_lie: Vec<TopoShape> = Vec::new();

                // `cxx:4537-4556`.
                for a_e in explorer_shapes(&a_w, ShapeType::Edge) {
                    if a_e.orientation() != Orientation::Internal {
                        b_remove_wire = false;
                        continue;
                    }
                    let a_lf = ancestors.get(&shape_map_key(&a_e));
                    let b_remove_edge = a_lf.map_or(0, |l| l.len()) == 1;
                    if b_remove_edge {
                        a_lie.push(a_e);
                    } else {
                        b_remove_wire = false;
                    }
                }

                // `cxx:4558-4565`.
                if b_remove_wire {
                    a_liw.push(a_w);
                } else if !a_lie.is_empty() {
                    remove_shapes(&mut a_w.clone(), &a_lie);
                }
            }

            // `cxx:4568-4571`.
            if !a_liw.is_empty() {
                remove_shapes(&mut a_f.clone(), &a_liw);
            }
        }
    }
}
