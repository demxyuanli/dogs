//! The face-split stage of `LocOpe_SplitDrafts::Perform`
//! (`LocOpe_SplitDrafts.cxx:353-403`).
//!
//! The wire produced by the previous stage is handed to `LocOpe_Spliter`, which
//! cuts the face `F` (and, through the ancestor tracking, the whole solid) along
//! it. The stage records the split result, which of the pieces ended up on the
//! left of the wire, the descendant faces of every face of `myShape`, the set of
//! edges/vertices that still sit on the original wire (`MapW`), and the
//! edge/face ancestor index of the result (`theMapEF`).

use crate::abs::ShapeType;
use crate::locope::split_drafts::perform::{CONSTRUCTION_ERROR, NULL_OBJECT};
use crate::locope::split_drafts::state::PerformState;
use crate::locope::split_drafts::SplitDrafts;
use crate::locope::{Spliter, WiresOnShape};
use crate::shape::{Face, TopoShape};
use crate::topo_tools::{self, explorer_shapes, AncestorMap};

impl SplitDrafts {
    /// `cxx:353-371`. `Ok(false)` is `if (!Spls.IsDone()) return;`, which leaves
    /// `myResult` null.
    pub(crate) fn split_face(&mut self, f: &Face, st: &mut PerformState) -> Result<bool, String> {
        // `cxx:355-357`.
        let my_shape = self
            .my_shape
            .clone()
            .ok_or_else(|| NULL_OBJECT.to_string())?;
        let mut wons = WiresOnShape::new(my_shape.clone());
        let mut spls = Spliter::with_shape(&my_shape);
        wons.bind_wire(&st.the_w, f)?;

        // `cxx:367-370`.
        if spls.perform(&mut wons).is_err() || !spls.is_done() {
            return Ok(false);
        }
        let res: TopoShape = spls
            .resulting_shape()
            .cloned()
            .ok_or_else(|| CONSTRUCTION_ERROR.to_string())?;
        st.the_left = spls.direct_left().map(|l| l.to_vec()).unwrap_or_default();

        // `cxx:375-383`: descendants of every face of `myShape`.
        for fac in explorer_shapes(&my_shape, ShapeType::Face) {
            let desc = spls.descendant_shapes(&fac);
            if let Some(slot) = self.my_map.change(&fac) {
                slot.extend(desc);
            }
        }

        // `cxx:385-399`: `MapW` marks the edges of `theW` and their vertices; the
        // null value is the reference's `TopoDS_Shape()`.
        for edg in topo_tools::explorer_shapes(&st.the_w.0, ShapeType::Edge) {
            if st.map_w.is_bound(&edg) {
                continue;
            }
            st.map_w.bind(&edg, None);
            for v in topo_tools::vertices_of(&edg) {
                if !st.map_w.is_bound(&v.0) {
                    st.map_w.bind(&v.0, None);
                }
            }
        }

        // `cxx:400-403`.
        st.the_map_ef =
            AncestorMap::map_shapes_and_ancestors(&res, ShapeType::Edge, ShapeType::Face);
        st.res = Some(res);
        Ok(true)
    }
}
