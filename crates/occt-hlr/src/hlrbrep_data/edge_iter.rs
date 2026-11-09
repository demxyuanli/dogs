//! Edge iteration of `HLRBRep_Data` - `InitEdge` (`HLRBRep_Data.cxx:1022-1075`),
//! `MoreEdge` (`:1077-1112`), `NextEdge` (`:1114-1207`) and `Edge`
//! (`:1209-1221`).
//!
//! `InitEdge(FI, MST)` starts the pass of one hiding face over the edges it may
//! hide. Two modes, selected by `iFaceTest = !FaceData.Simple()`:
//!
//! * `iFaceTest` - the face is not simple, so **every** edge of the face is a
//!   hiding-boundary candidate: `myFaceItr2` walks the face wire blocks
//!   (`HLRBRep_FaceIterator`).
//! * otherwise - a simple face only tests the edges `InitBoundSort` kept, in
//!   that order: `myEdgeIndices(myCurSortEd)`, `myCurSortEd = 1..myNbrSortEd`.
//!
//! [`HlrBRepData::next_edge`] is the rejection cascade (`cxx:1163-1201`):
//! vertical edges, edges already hidden at this level (`HideCount`), edges the
//! status says are `AllHidden`, edges whose projected box is separated from the
//! face box in one of the sixteen directions, and edges the face
//! `IsAbove`. Each rejection re-enters the loop with `skip = true`, which is
//! the OCCT recursion `NextEdge()` (default `skip = true`).
//!
//! Deviations from the OCCT members:
//!
//! * `myLEData` / `myLEGeom` / `myLEMinMax` / `iFaceData` / `iFaceGeom` are raw
//!   pointers into the arrays. The port keeps the OCCT 1-based **index**
//!   (`my_le`, `i_face`) and caches the values the algorithms read into the
//!   same-named scalar fields (`my_le_min_max`, `my_le_tol`, `my_le_type`,
//!   `i_face_min_max`, `i_face_type`), so the borrow checker never needs two
//!   mutable views of `my_e_data`.
//! * `myClassifier` is a `handle<BRepTopAdaptor_TopolTool>` obtained from the
//!   `MST` map, which OCCT binds on demand; the port stores `Arc<TopolTool>`
//!   (see [`occt_topo::brep_top_adaptor_tool`]).
//! * `mySLProps.SetSurface` / `myIntersector.Load` receive the face geometry by
//!   reference. OCCT stores the pointer and dereferences it later; the port
//!   clones the surface inside those two calls.

use occt_topo::brep_top_adaptor_tool::BrepTopAdaptorTool;
use crate::hlr_algo::MinMaxIndices;
use crate::hlrtopobrep::out_liner::ShapeToTopAdaptorToolMap;
use occt_topo::topo_tools::shape_map_key;

use super::HlrBRepData;
use super::MINMAX_SIGN_MASK;

/// The sixteen `((A.Max[i] - B.Min[i]) & 0x80008000) != 0` tests of
/// `cxx:1181-1198` (edge against face), `cxx:1245-1266` (wire against edge) and
/// `cxx:1302-1316` (face edge against edge): `true` as soon as one direction
/// separates the two boxes.
pub(crate) fn separated(a: &MinMaxIndices, b: &MinMaxIndices) -> bool {
    for i in 0..8 {
        if ((a.max[i] - b.min[i]) & MINMAX_SIGN_MASK) != 0 {
            return true;
        }
        if ((b.max[i] - a.min[i]) & MINMAX_SIGN_MASK) != 0 {
            return true;
        }
    }
    false
}

impl HlrBRepData {
    /// `InitEdge(FI, MST)` (`cxx:1022-1075`).
    ///
    /// Raises `myHideCount` by two (one level for the hiding face, one for its
    /// own edges), loads the face geometry into `mySLProps` and
    /// `myIntersector`, binds the `MST` classifier entry for the face, then
    /// either initialises the wire iterator (`iFaceTest`) or pre-hides every
    /// edge of a simple face and starts the sorted-edge cursor (`cxx:1065-1071`).
    pub fn init_edge(&mut self, fi: i32, mst: &mut ShapeToTopAdaptorToolMap) {
        self.my_hide_count += 1;
        self.my_hide_count += 1;

        self.i_face = fi;
        let idx = fi as usize;
        self.i_face_geom = Some(fi);
        self.i_face_back = self.my_f_data[idx].back();
        self.i_face_simp = self.my_f_data[idx].simple();
        self.i_face_min_max = self.my_f_data[idx].wires().map(|w| *w.min_max());
        self.i_face_type = self.my_f_data[idx].geometry().get_type();
        self.i_face_test = !self.i_face_simp;
        self.my_sl_props.set_surface(self.my_f_data[idx].geometry());
        self.my_intersector.load(self.my_f_data[idx].geometry());

        // `p1->Surface()` then `BRepAdaptor_Surface::Face()`; the classifier
        // map is keyed by the face (`cxx:1045-1053`).
        let face = self.my_f_data[idx]
            .geometry()
            .surface_adaptor()
            .face()
            .cloned();
        self.my_classifier = None;
        if let Some(face) = face {
            let key = shape_map_key(face.shape());
            if !mst.contains_key(&key) {
                // `BRepTopAdaptor_Tool BRT(topodsface, Precision::PConfusion())`
                // (`cxx:1051`), i.e. the restricted adaptor constructor.
                mst.insert(key.clone(), BrepTopAdaptorTool::of_face(&face));
            }
            self.my_classifier = mst.get(&key).map(|brt| brt.topol_tool_shared());
        }

        if self.i_face_test {
            self.i_face_smpl = !self.my_f_data[idx].cut();
            self.my_face_itr2.init_edge(&self.my_f_data[idx]);
        } else {
            self.my_face_itr1.init_edge(&self.my_f_data[idx]);
            while self.my_face_itr1.more_edge() {
                self.my_fe = self.my_face_itr1.edge();
                // `myEData(myFE).HideCount(myHideCount - 1)`: the face must be
                // jumped by its own edges (`cxx:1067-1070`).
                self.my_e_data[self.my_fe as usize].set_hide_count(self.my_hide_count - 1);
                self.my_face_itr1.next_edge();
            }
            self.my_cur_sort_ed = 1;
        }
        self.next_edge(false);
    }

    /// `MoreEdge()` (`cxx:1077-1112`).
    ///
    /// In test mode it loads the next wire edge (and pre-hides it unless it is
    /// a double edge); at the end of the wires it **becomes** a simple face
    /// (`iFaceTest = false`, `iFaceSimp = iFaceSmpl`, `FaceData.Simple`), starts
    /// the sorted cursor and falls through to the sort test.
    pub fn more_edge(&mut self) -> bool {
        if self.i_face_test {
            if self.my_face_itr2.more_edge() {
                self.load_le_from_face_iterator();
                if !self.my_le_double {
                    self.my_e_data[self.my_le as usize].set_hide_count(self.my_hide_count - 1);
                }
                return true;
            } else {
                // `cxx:1101-1107`: at the end of the test the face is known to
                // be simple iff it was not cut.
                self.i_face_test = false;
                self.i_face_simp = self.i_face_smpl;
                let idx = self.i_face as usize;
                self.my_f_data[idx].set_simple(self.i_face_simp);
                self.my_cur_sort_ed = 1;
                self.next_edge(false);
            }
        }
        self.my_cur_sort_ed <= self.my_nbr_sort_ed
    }

    /// `NextEdge(skip)` (`cxx:1114-1207`).
    ///
    /// `skip` advances the current iterator first (the wire one in test mode,
    /// the sort cursor otherwise). The four OCCT self-recursions
    /// (`NextEdge()` = `NextEdge(true)`) are the `continue`s of the loop.
    pub fn next_edge(&mut self, skip: bool) {
        let mut skip = skip;
        loop {
            if skip {
                if self.i_face_test {
                    self.my_face_itr2.next_edge();
                } else {
                    self.my_cur_sort_ed += 1;
                }
            }
            skip = true;

            if !self.more_edge() {
                return;
            }

            if self.i_face_test {
                // `cxx:1133-1150`: the test-mode edge needs no other filter
                // than "vertical, or already hidden at this level". The
                // rejection re-enters `NextEdge()`, which re-assigns
                // `myLEData` to the next edge, so OCCT's trailing
                // `myLEData->HideCount(myHideCount - 1)` only ever lands on the
                // accepted edge - the rejected one is left untouched.
                self.load_le_from_face_iterator();
                let i = self.my_le as usize;
                if self.my_e_data[i].vertical()
                    || (self.my_le_double && self.my_e_data[i].hide_count() == self.my_hide_count - 1)
                {
                    continue;
                }
                self.my_e_data[i].set_hide_count(self.my_hide_count - 1);
                return;
            }

            // `cxx:1152-1161`: the sorted edge is not a face boundary, so the
            // four flags are cleared.
            self.my_le = self.edge();
            self.my_le_out_line = false;
            self.my_le_internal = false;
            self.my_le_double = false;
            self.my_le_iso_line = false;
            self.load_le_data();

            let i = self.my_le as usize;
            if self.my_e_data[i].vertical() {
                continue;
            }
            if self.my_e_data[i].hide_count() > self.my_hide_count - 2 {
                continue;
            }
            if self.my_e_data[i].status().all_hidden() {
                continue;
            }
            let Some(face_min_max) = self.i_face_min_max else {
                // `iFaceMinMax` is a pointer set by `InitEdge`; OCCT would
                // dereference a null pointer here.
                return;
            };
            let Some(le_min_max) = self.my_le_min_max else {
                return;
            };
            if separated(&face_min_max, &le_min_max) {
                continue;
            }
            let fi = self.i_face as usize;
            if self.my_f_data[fi].geometry().is_above(
                self.i_face_back,
                self.my_e_data[i].geometry(),
                self.my_le_tol as f64,
            ) {
                continue;
            }
            // edge is OK
            return;
        }
    }

    /// `Edge()` (`cxx:1209-1221`): the wire-iterator edge in test mode, the
    /// sorted edge index otherwise.
    pub fn edge(&self) -> i32 {
        if self.i_face_test {
            self.my_face_itr2.edge()
        } else {
            self.my_edge_indices[self.my_cur_sort_ed as usize]
        }
    }

    /// `myLE = myFaceItr2.Edge()` and the seven fields copied with it
    /// (`cxx:1083-1093`, `:1133-1143`).
    fn load_le_from_face_iterator(&mut self) {
        self.my_le = self.my_face_itr2.edge();
        self.my_le_out_line = self.my_face_itr2.out_line();
        self.my_le_internal = self.my_face_itr2.internal();
        self.my_le_double = self.my_face_itr2.double();
        self.my_le_iso_line = self.my_face_itr2.iso_line();
        self.load_le_data();
    }

    /// `myLEData = &myEData(myLE)` and the three derived fields
    /// (`cxx:1089-1093`).
    fn load_le_data(&mut self) {
        let i = self.my_le as usize;
        self.my_le_min_max = Some(*self.my_e_data[i].min_max());
        self.my_le_tol = self.my_e_data[i].tolerance();
        self.my_le_type = self.my_e_data[i].geometry().get_type();
    }

    /// `myFEData = &myEData(myFE)` and the two derived fields
    /// (`cxx:1272-1277`).
    pub(crate) fn load_fe_data(&mut self) {
        let i = self.my_fe as usize;
        self.my_fe_tol = self.my_e_data[i].tolerance();
        self.my_fe_type = self.my_e_data[i].geometry().get_type();
    }
}
