//! `Perform` (`TopOpeBRepTool_FuseEdges.cxx:384-468`).
//!
//! Every chain whose fused edge was built (`myMapEdg(iLst)`) is substituted
//! into `myShape` with `BRepTools_Substitution`:
//!
//! - the chain's first edge is replaced by the fused edge (`cxx:421-422`),
//! - every other edge of the chain is replaced by an **empty** list
//!   (`cxx:427-435`), i.e. it disappears from its ancestors,
//! - `Bsub.Build(myShape)` (`cxx:442`) rebuilds the parents,
//! - the faces that changed are recorded in `myMapFaces` **before** the shape
//!   itself is replaced (`cxx:446-457`), which is what
//!   `Builder::End`'s `UpdateSplitAndMerged` reads,
//! - `myShape` then becomes the rebuilt shape when it was copied
//!   (`cxx:458-460`).
//!
//! A chain can be missing from `myMapEdg` when `UpdatePCurve` refused it
//! (`cxx:373-377`); such a chain is skipped here (`cxx:411-414`) and its edges
//! stay as they are.

use crate::abs::ShapeType;
use crate::brep_tools_substitution::BRepToolsSubstitution;
use crate::topo_tools::{explorer_shapes, is_same};

use super::{FuseEdges, FuseEdgesError};

impl FuseEdges {
    /// `Perform()` (`cxx:384-468`).
    pub(super) fn perform(&mut self) -> Result<(), FuseEdgesError> {
        if !self.my_result_edges_done {
            self.build_list_result_edges()?;
        }

        let mut bsub = BRepToolsSubstitution::new();
        let empty: Vec<crate::shape::TopoShape> = Vec::new();

        // `cxx:405-437`
        let keys: Vec<i32> = self.my_map_lst_edg.keys().copied().collect();
        for i_lst in keys {
            // `cxx:411-414`
            let Some(new_edge) = self.my_map_edg.get(&i_lst).cloned() else {
                continue;
            };
            let Some(lmap_edg) = self.my_map_lst_edg.get(&i_lst).cloned() else {
                continue;
            };
            let Some(old_edge) = lmap_edg.first().cloned() else {
                continue;
            };

            // `cxx:421-422`
            let edge_to_subs = vec![new_edge];
            bsub.substitute(&old_edge, &edge_to_subs);

            // `cxx:427-435`: the other edges of the list are removed.
            for e in &lmap_edg {
                if !is_same(&old_edge, e) {
                    bsub.substitute(e, &empty);
                }
            }
        }

        // `cxx:442`
        bsub.build(&self.my_shape);

        // `cxx:446-457`: the modified faces, read off the **old** shape.
        self.my_map_faces.clear();
        for face in explorer_shapes(&self.my_shape, ShapeType::Face) {
            if bsub.is_copied(&face) {
                if let Some(first) = bsub.copy(&face).and_then(|l| l.first()).cloned() {
                    self.my_map_faces.push((face, first));
                }
            }
        }

        // `cxx:458-460`
        if bsub.is_copied(&self.my_shape) {
            if let Some(first) = bsub.copy(&self.my_shape).and_then(|l| l.first()).cloned() {
                self.my_shape = first;
            }
        }

        self.my_shape_done = true;
        Ok(())
    }
}
