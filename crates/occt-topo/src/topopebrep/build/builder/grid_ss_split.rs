//! `TopOpeBRepBuild_GridSS.cxx` — splitting a face of the solid/solid grid.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_GridSS.cxx`
//! (OCCT 8.0.0):
//!
//! | OCCT | lines |
//! |---|---|
//! | `GSplitFace(FOR, GG1, LSclass)` | `1348-1453` |
//!
//! `FillOnPatches(anEdgesON, aBaseFace, avoidMap)` (`cxx:1784-2036`, with
//! `AreFacesCoincideInArea` `:1595-1782` and `FindFacesTouchingEdge`
//! `:2037-2085`) lives in `build/builder/on_patches.rs` and is called at
//! `cxx:1441`. `SplitEvisoONperiodicF` and its helpers live in
//! `build/builder/split_eviso_on_periodic.rs`.
//!
//! UNPORTED: the `OCCT_DEBUG` blocks (`cxx:1366-1375`, `cxx:1382-1386`,
//! `:1400-1405`, `:1411-1415`) — trace and dump hooks only.

use crate::abs::{Orientation, ShapeType};
use crate::shape::TopoShape;
use crate::topo_tools::explorer_shapes;
use crate::topopebrep::ds::data_structure::IndexedShapeSet;
use crate::topopebrep::build::g_topo::GTopo;
use crate::topopebrep::build::wire_edge_set::WireEdgeSet;

use super::Builder;

impl Builder {
    /// `GSplitFace(FOR, GG1, LSclass)` (`cxx:1348-1453`).
    pub fn g_split_face(&mut self, for_: &TopoShape, gg1: &GTopo, ls_class: &[TopoShape]) {
        // `cxx:1351-1353`: the grid is forced to `Reverse = false`.
        let mut g1 = gg1.clone();
        g1.set_reverse(false);

        let (tb1, _tb2) = g1.states_on();
        let _ = g1.type_();

        // `cxx:1361-1362`: work on a `FORWARD` face.
        let mut ff = for_.clone();
        ff.set_orientation(Orientation::Forward);

        // `cxx:1379-1380`: a WireEdgeSet on `FF`.
        let mut wes = WireEdgeSet::new(&ff);

        // `cxx:1382`: add the ON parts (edges ON solid).
        self.g_fill_on_parts_wes(for_, &g1, ls_class, &mut wes);

        // `cxx:1388-1399`: save these edges.
        let mut an_edges_on: Vec<TopoShape> = Vec::new();
        if self.my_process_on {
            let to_rev_ori = self.opefus();
            for e in wes.start_elements() {
                if to_rev_ori {
                    let mut e = e.clone();
                    e.reverse();
                    an_edges_on.push(e);
                } else {
                    an_edges_on.push(e.clone());
                }
            }
            self.my_on_elem_map.clear();
        }

        // `cxx:1402`: split the edges of `FF`, adding the split edges to `WES`.
        self.g_fill_face_wes(&ff, ls_class, &g1, &mut wes);
        let n1 = wes.start_elements().len() as i32;

        // `cxx:1410`: add the edges built on the curves supported by `FF`.
        self.g_fill_curve_topology_wes(&ff, &g1, &mut wes);
        let n2 = wes.start_elements().len() as i32;

        // `cxx:1416-1419`: `myEdgeAvoid` = the start elements added by
        // `GFillCurveTopologyWES`.
        self.my_edge_avoid.clear();
        let mut avoid: Vec<TopoShape> = Vec::new();
        Self::copy_list_range(wes.start_elements(), n1 + 1, n2, &mut avoid);
        self.my_edge_avoid.extend(avoid);

        // `cxx:1421`: mark `FF` as split on `TB1`.
        self.mark_split(&ff, tb1, true);

        // `cxx:1424-1425`: build the new faces of `FF` from `WES`.
        let mut lof: Vec<TopoShape> = Vec::new();
        self.gwes_make_faces(&ff, &mut wes, &mut lof);

        if self.my_process_on && (!an_edges_on.is_empty() || self.my_on_elem_map.extent() != 0) {
            // `cxx:1431-1439`: the map of edges already used by `LOF`.
            let mut a_map_oe = IndexedShapeSet::new();
            for f in &lof {
                for e in explorer_shapes(f, ShapeType::Edge) {
                    a_map_oe.add(&e);
                }
            }

            // `cxx:1441`
            self.fill_on_patches(&an_edges_on, for_, &a_map_oe);

            self.my_on_elem_map.clear();
        }

        // `cxx:1447-1449`: `LOFS` = the split faces of state `TB1` of `FF`.
        let empty: Vec<TopoShape> = Vec::new();
        let mut lofs: Vec<TopoShape> = Vec::new();
        self.g_keep_shapes(&ff, &empty, tb1, &lof, &mut lofs);
        let lo = self.change_split(&ff, tb1);
        lo.clear();
        lo.extend(lofs);
    }
}
