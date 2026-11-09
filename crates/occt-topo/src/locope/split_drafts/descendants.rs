//! Descendant fix-ups and assembly (`LocOpe_SplitDrafts.cxx:1360-1435`).
//!
//! `myMap` is rewritten so that every entry points at the shapes the draft left
//! behind: edge entries through `MapW`, the split face `F` through the drafted
//! faces, everything else through `BRepTools_Substitution`. The distinct faces
//! of the whole map are then closed into a shape by `LocOpe_BuildShape`.

use crate::abs::ShapeType;
use crate::locope::build_shape::BuildShape;
use crate::locope::split_drafts::state::PerformState;
use crate::locope::split_drafts::SplitDrafts;
use crate::shape::{Face, TopoShape};
use crate::topo_tools::is_same;

impl SplitDrafts {
    /// `cxx:1360-1435`. `Err` is the reference's `return` on a substitution that
    /// does not have exactly one descendant.
    pub(crate) fn build_descendants(
        &mut self,
        f: &Face,
        st: &mut PerformState,
    ) -> Result<(), String> {
        // `cxx:1364-1418`: the reference assigns into `myMap` while iterating
        // it; the port snapshots the entries first.
        let pairs: Vec<(TopoShape, Vec<TopoShape>)> = self
            .my_map
            .iter()
            .map(|(k, v)| (k.clone(), v.clone()))
            .collect();

        for (key, value) in pairs {
            if key.shape_type() == ShapeType::Edge {
                // `cxx:1367-1379`.
                let mut thedesc: Vec<TopoShape> = Vec::new();
                st.the_map.clear();
                for v in &value {
                    // `MapW(v)` raises in the reference when `v` was never
                    // bound; the port skips it.
                    let Some(mapped) = st.map_w.find(v).cloned().flatten() else {
                        continue;
                    };
                    if st.the_map.add(&mapped) {
                        thedesc.push(mapped);
                    }
                }
                self.my_map.bind(&key, thedesc);
            } else if is_same(&key, &f.0) {
                // `cxx:1380-1388`.
                let faces = st.the_faces.clone();
                self.my_map.bind(&key, faces);
            } else {
                // `cxx:1389-1417`.
                let mut thedesc: Vec<TopoShape> = Vec::new();
                st.the_map.clear();
                for v in &value {
                    if st.subs.is_copied(v) {
                        let copied: Vec<TopoShape> = st.subs.copy(v).map(|s| s.to_vec()).unwrap_or_default();
                        if copied.len() != 1 {
                            return Err("Invalid number of descendant".to_string());
                        }
                        let c = copied[0].clone();
                        if st.the_map.add(&c) {
                            thedesc.push(c);
                        }
                    } else if st.the_map.add(v) {
                        thedesc.push(v.clone());
                    }
                }
                self.my_map.bind(&key, thedesc);
            }
        }

        // `cxx:1420-1431`.
        st.the_map.clear();
        st.the_faces.clear();
        for (_, value) in self.my_map.iter() {
            for v in value {
                if v.shape_type() == ShapeType::Face && st.the_map.add(v) {
                    st.the_faces.push(v.clone());
                }
            }
        }

        // `cxx:1432-1434`.
        let mut bs = BuildShape::new();
        bs.perform(&st.the_faces)?;
        self.my_result = bs.shape().cloned();
        Ok(())
    }
}
