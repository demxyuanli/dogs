//! `ChFi3d_Builder::ExtentAnalyse` — the pre-analysis that evaluates the
//! extensions at every vertex carrying stripes.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder.cxx:144-174`
//! (OCCT 8.0.0).
//!
//! Three branches, by the number of stripes loaded on the vertex
//! (`myVDataMap(iv).Extent()`):
//!
//! | `nbs` | condition | call |
//! |---|---|---|
//! | 1 | — | `ExtentOneCorner(Vtx, myVDataMap.FindFromIndex(iv).First())` |
//! | 2 | `nbedges <= 3` | `ExtentTwoCorner(Vtx, myVDataMap.FindFromIndex(iv))` |
//! | 3 | `nbedges <= 3` | `ExtentThreeCorner(Vtx, myVDataMap.FindFromIndex(iv))` |
//!
//! `nbedges` is `ChFi3d_NumberOfSharpEdges(Vtx, myVEMap, myEFMap)`; the
//! `ChFi3d_NumberOfEdges` call above it is commented out in OCCT (`:151`) and
//! is not translated.

use super::{ChFi3dBuilder, ChFi3dCornerOps};
use crate::chfi3d::builder_helpers::number_of_sharp_edges;
use crate::shape::Vertex;

impl ChFi3dBuilder {
    /// `ExtentAnalyse()` (`ChFi3d_Builder.cxx:144-174`).
    ///
    /// `iv` is OCCT's 1-based index over `myVDataMap`; the port walks the same
    /// order because [`crate::chfi_ds::StripeMap`] keeps insertion order.
    pub fn extent_analyse(&mut self, ops: &mut dyn ChFi3dCornerOps) -> Result<(), String> {
        let extent = self.my_v_data_map.extent();
        for iv in 1..=extent {
            let nbs = self
                .my_v_data_map
                .find_from_index(iv)
                .map(|l| l.len())
                .ok_or_else(|| "ExtentAnalyse: myVDataMap index out of range".to_string())?;
            let vtx = self
                .my_v_data_map
                .find_key(iv)
                .cloned()
                .ok_or_else(|| "ExtentAnalyse: myVDataMap index out of range".to_string())?;
            let vertex = Vertex::wrap(vtx.clone())
                .ok_or_else(|| "ExtentAnalyse: myVDataMap key is not a vertex".to_string())?;
            let nbedges = number_of_sharp_edges(&vertex, &self.my_ve_map, &self.my_ef_map)?;
            match nbs {
                1 => {
                    // `myVDataMap.FindFromIndex(iv).First()`.
                    let first = self
                        .my_v_data_map
                        .find_from_index(iv)
                        .and_then(|l| l.first())
                        .cloned()
                        .ok_or_else(|| {
                            "ExtentAnalyse: empty stripe list for ExtentOneCorner".to_string()
                        })?;
                    ops.extent_one_corner(&vtx, first)?;
                }
                2 => {
                    if nbedges <= 3 {
                        let list = self
                            .my_v_data_map
                            .find_from_index(iv)
                            .map(|l| l.to_vec())
                            .unwrap_or_default();
                        ops.extent_two_corner(&vtx, &list)?;
                    }
                }
                3 => {
                    if nbedges <= 3 {
                        let list = self
                            .my_v_data_map
                            .find_from_index(iv)
                            .map(|l| l.to_vec())
                            .unwrap_or_default();
                        ops.extent_three_corner(&vtx, &list)?;
                    }
                }
                // `default: break;` (`:170-171`).
                _ => {}
            }
        }
        Ok(())
    }
}
