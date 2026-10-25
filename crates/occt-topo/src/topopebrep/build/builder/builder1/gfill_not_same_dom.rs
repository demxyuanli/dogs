//! `TopOpeBRepBuild_Builder1` — the faces, wires and edges *without* a
//! same-domain twin.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_Builder1.cxx`
//! (OCCT 8.0.0):
//!
//! | OCCT | lines |
//! |---|---|
//! | `GFillSolidSFS(SO1, LSO2, G1, SFS)` (virtual, `hxx:69-73`) | `186-241` |
//! | `GFillShellSFS(SH, LSO2, G1, SFS)` (virtual, `hxx:76-79`) | `243-301` |
//! | `GFillFaceNotSameDomSFS(FOR, LSO2, Gin, SFS)` | `303-378` |
//! | `GFillFaceNotSameDomWES(FOR1, LFclass, G1, WES)` | `380-430` |
//! | `GFillWireNotSameDomWES(W, LSclass, G1, WES)` | `431-477` |
//! | `GFillEdgeNotSameDomWES(EOR, LSclass, G1, WES)` | `478-532` |
//!
//! The two virtual members are reached through
//! `Builder::g_fill_solid_sfs` / `Builder::g_fill_shell_sfs`, which branch on
//! [`Builder::is_builder1`] exactly where OCCT's vtable does; the bodies here
//! fall back to the base members when `GLOBAL_USE_NEW_BUILDER` is false
//! (`cxx:192-196`, `:250-254`), as OCCT's overrides do.
//!
//! The non-virtual members are called only by the two above, so they need no
//! dispatch: OCCT declares them in `Builder1` alone (`hxx:85-107`).
//!
//! `GFillSplitsPVS` (`hxx:83-85`) has no definition anywhere in OCCT 8.0.0 —
//! declared only, nothing to translate.

use crate::abs::{ShapeType, State};
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, TopoShape};
use crate::topo_tools::explorer_shapes;
use crate::topopebrep::build::g_topo::GTopo;
use crate::topopebrep::build::shell_face_set::ShellFaceSet;
use crate::topopebrep::build::wire_edge_set::WireEdgeSet;
use crate::topopebrep::ds::data_structure::{IndexedShapeData, IndexedShapeSet};

use super::global_use_new_builder;
use super::super::Builder;

/// `TopExp::MapShapesAndAncestors(S, TopAbs_EDGE, TopAbs_FACE, M)`
/// (`TopExp.cxx:80-120`) into an `NCollection_IndexedDataMap`: every edge of
/// `shape` mapped to the faces it sits under, appended in the ancestor walk
/// order.
fn add_edge_face_ancestors(shape: &TopoShape, map: &mut IndexedShapeData<Vec<TopoShape>>) {
    for f in explorer_shapes(shape, ShapeType::Face) {
        for e in explorer_shapes(&f, ShapeType::Edge) {
            if map.find_from_key(&e).is_none() {
                map.add(&e, vec![f.clone()]);
            } else if let Some(faces) = map.change_from_key(&e) {
                faces.push(f.clone());
            }
        }
    }
}

impl Builder {
    /// `DS().GetShapeWithState(S).State()` (`ShapeWithStateMap.cxx:1327-1341`):
    /// the dummy map OCCT returns for an absent shape answers `UNKNOWN`
    /// (`ShapeWithState.cxx:22-26`).
    fn state_with_state(&self, s: &TopoShape) -> State {
        self.data_structure()
            .and_then(|hds| hds.ds().get_shape_with_state(s))
            .map(|sws| sws.state())
            .unwrap_or(State::Unknown)
    }

    /// The split parts of a shape in one state, cloned so the data structure
    /// is not borrowed while the caller mutates the builder.
    pub(crate) fn parts_in_state(&self, s: &TopoShape, state: State) -> Vec<TopoShape> {
        self.data_structure()
            .and_then(|hds| hds.ds().get_shape_with_state(s))
            .map(|sws| sws.part(state).to_vec())
            .unwrap_or_default()
    }

    /// `Builder1::GFillSolidSFS(SO1, LSO2, G1, SFS)` (`cxx:186-241`).
    ///
    /// `myMapOfEdgeFaces` is filled with the edges of *both* operands
    /// (`cxx:199-203`); `PerformPieceOn2D` and `GFillEdgeSameDomWES` read it.
    pub(crate) fn b1_g_fill_solid_sfs(
        &mut self,
        so1: &TopoShape,
        lso2: &[TopoShape],
        g1: &GTopo,
        sfs: &mut ShellFaceSet,
    ) {
        // `cxx:192-196`
        if !global_use_new_builder() {
            self.base_g_fill_solid_sfs(so1, lso2, g1, sfs);
            return;
        }

        // `cxx:199-203`
        self.b1.my_map_of_edge_faces.clear();
        let shape1 = self.my_shape1.clone();
        let shape2 = self.my_shape2.clone();
        if let Some(s) = shape1.as_ref() {
            add_edge_face_ancestors(s, &mut self.b1.my_map_of_edge_faces);
        }
        if let Some(s) = shape2.as_ref() {
            add_edge_face_ancestors(s, &mut self.b1.my_map_of_edge_faces);
        }

        let (tb1, _tb2) = g1.states_on();
        let rev_ori1 = g1.is_to_reverse1();

        // `cxx:211-213`: `SOF = SO1; mySolidToFill = TopoDS::Solid(SOF)`.
        let sof = so1.clone();
        self.my_solid_to_fill = Some(sof.clone());

        // `cxx:215-238`
        for sh in explorer_shapes(&sof, ShapeType::Shell) {
            let hasshape = self
                .data_structure()
                .map(|hds| hds.has_shape(&sh, true))
                .unwrap_or(false);

            if !hasshape {
                // `cxx:221-233`: the shell is not in the DS, its state is.
                let keep = self.state_with_state(&sh) == tb1;
                if keep {
                    let mut sh = sh;
                    let newori = Builder::orient(sh.orientation(), rev_ori1);
                    sh.set_orientation(newori);
                    sfs.add_shape(&sh);
                }
            } else {
                // `cxx:237`: the shell has faces with geometry.
                self.g_fill_shell_sfs(&sh, lso2, g1, sfs);
            }
        }
    }

    /// `Builder1::GFillShellSFS(SH, LSO2, G1, SFS)` (`cxx:243-301`).
    pub(crate) fn b1_g_fill_shell_sfs(
        &mut self,
        sh: &TopoShape,
        lso2: &[TopoShape],
        g1: &GTopo,
        sfs: &mut ShellFaceSet,
    ) {
        // `cxx:250-254`
        if !global_use_new_builder() {
            self.base_g_fill_shell_sfs(sh, lso2, g1, sfs);
            return;
        }

        let (tb1, _tb2) = g1.states_on();
        let rev_ori1 = g1.is_to_reverse1();

        // 1) same-domain faces and faces absent from the DS (`cxx:265-287`)
        for for_ in explorer_shapes(sh, ShapeType::Face) {
            let hasshape = self
                .data_structure()
                .map(|hds| hds.has_shape(&for_, true))
                .unwrap_or(false);

            if !hasshape {
                // `cxx:269-280`
                let keep = self.state_with_state(&for_) == tb1;
                if keep {
                    let mut for_ = for_;
                    let newori = Builder::orient(for_.orientation(), rev_ori1);
                    for_.set_orientation(newori);
                    sfs.add_element(&for_);
                }
                continue;
            }

            // `cxx:281-285`
            let hsd = self
                .data_structure()
                .map(|hds| hds.has_same_domain(&for_, true))
                .unwrap_or(false);
            if hsd && !self.b1.my_same_dom_map.contains(&for_) {
                self.b1_g_fill_face_same_dom_sfs(&for_, lso2, g1, sfs);
            }
        }

        // 2) all other faces (`cxx:289-297`)
        for for_ in explorer_shapes(sh, ShapeType::Face) {
            let in_ds = self
                .data_structure()
                .map(|hds| hds.has_shape(&for_, true))
                .unwrap_or(false);
            let hsd = self
                .data_structure()
                .map(|hds| hds.has_same_domain(&for_, true))
                .unwrap_or(false);
            if !in_ds || hsd {
                continue;
            }
            self.b1_g_fill_face_not_same_dom_sfs(&for_, lso2, g1, sfs);
        }
    }

    /// `GFillFaceNotSameDomSFS(FOR, LSO2, Gin, SFS)` (`cxx:303-378`).
    pub(crate) fn b1_g_fill_face_not_same_dom_sfs(
        &mut self,
        for_: &TopoShape,
        lso2: &[TopoShape],
        gin: &GTopo,
        sfs: &mut ShellFaceSet,
    ) {
        // `cxx:308-310`: the local grid always has `RevOri = false`.
        let mut g1 = gin.clone();
        g1.set_reverse(false);

        let (tb1, _tb2) = g1.states_on();

        // `cxx:315-317`: work on a FORWARD face.
        let mut ff = for_.clone();
        ff.set_orientation(crate::abs::Orientation::Forward);

        // `cxx:320`
        let mut wes = WireEdgeSet::new(&ff);

        // `cxx:322-323`: the ON parts of the solid.
        self.g_fill_on_parts_wes(for_, &g1, lso2, &mut wes);

        // `cxx:326-336`: save these edges, ON matter reversed for a fuse.
        let mut an_edges_on: Vec<TopoShape> = Vec::new();
        if self.my_process_on {
            let to_rev_ori = self.opefus();
            for e in wes.start_elements() {
                let mut e = e.clone();
                if to_rev_ori {
                    e.reverse();
                }
                an_edges_on.push(e);
            }
            self.my_on_elem_map.clear();
        }

        // `cxx:339` then `:342`
        self.b1_g_fill_face_not_same_dom_wes(&ff, lso2, &g1, &mut wes);
        self.g_fill_curve_topology_wes(&ff, &g1, &mut wes);

        // `cxx:344-347`
        self.my_edge_avoid.clear();
        self.mark_split(&ff, tb1, true);

        // `cxx:350-351`
        let mut lof: Vec<TopoShape> = Vec::new();
        self.gwes_make_faces(&ff, &mut wes, &mut lof);

        // `cxx:353-369`: patches made only of ON matter, avoiding the edges
        // already used by `LOF`.
        if self.my_process_on && (!an_edges_on.is_empty() || self.my_on_elem_map.extent() != 0) {
            let mut a_map_oe = IndexedShapeSet::new();
            for f in lof.iter() {
                for e in explorer_shapes(f, ShapeType::Edge) {
                    a_map_oe.add(&e);
                }
            }
            self.fill_on_patches(&an_edges_on, for_, &a_map_oe);
            self.my_on_elem_map.clear();
        }

        // `cxx:372-376`
        let mut lofs: Vec<TopoShape> = Vec::new();
        self.g_keep_shapes(&ff, &[], tb1, &lof, &mut lofs);
        let lo = self.change_split(&ff, tb1);
        lo.clear();
        lo.extend(lofs);

        // `cxx:377`
        self.g_split_face_sfs(for_, lso2, gin, sfs);
    }

    /// `GFillFaceNotSameDomWES(FOR1, LFclass, G1, WES)` (`cxx:380-430`).
    pub(crate) fn b1_g_fill_face_not_same_dom_wes(
        &mut self,
        for1: &TopoShape,
        lf_class: &[TopoShape],
        g1: &GTopo,
        wes: &mut WireEdgeSet,
    ) {
        let (tb1, _tb2) = g1.states_on();
        let rev_ori1 = g1.is_to_reverse1();

        // `cxx:389`
        self.b1.my_source_shapes.clear();

        // `cxx:392-393`
        let mut ff = for1.clone();
        ff.set_orientation(crate::abs::Orientation::Forward);

        // `cxx:395-427`
        for w in explorer_shapes(&ff, ShapeType::Wire) {
            let hasshape = self
                .data_structure()
                .map(|hds| hds.has_shape(&w, true))
                .unwrap_or(false);

            if !hasshape {
                // `cxx:402-422`
                let sh_st = self.state_with_state(&w);
                let keep = sh_st == tb1;
                if keep || (self.my_process_on && sh_st == State::On) {
                    let mut w = w;
                    let newori = Builder::orient(w.orientation(), rev_ori1);
                    w.set_orientation(newori);
                    if keep {
                        wes.add_shape(&w);
                    } else {
                        self.my_on_elem_map.add(&w);
                    }
                    self.b1.my_source_shapes.add(&w);
                }
            } else {
                // `cxx:423-426`
                self.b1_g_fill_wire_not_same_dom_wes(&w, lf_class, g1, wes);
            }
        }
    }

    /// `GFillWireNotSameDomWES(W, LSclass, G1, WES)` (`cxx:431-477`).
    pub(crate) fn b1_g_fill_wire_not_same_dom_wes(
        &mut self,
        w: &TopoShape,
        ls_class: &[TopoShape],
        g1: &GTopo,
        wes: &mut WireEdgeSet,
    ) {
        let (tb1, _tb2) = g1.states_on();
        let rev_ori1 = g1.is_to_reverse1();

        // `cxx:440`: the wire keeps its own orientation (`WW = W`).
        for eor in explorer_shapes(w, ShapeType::Edge) {
            let hasshape = self
                .data_structure()
                .map(|hds| hds.has_shape(&eor, true))
                .unwrap_or(false);

            if !hasshape {
                // `cxx:449-467`
                let sh_st = self.state_with_state(&eor);
                let keep = sh_st == tb1;
                if keep || (self.my_process_on && sh_st == State::On) {
                    let mut eor = eor;
                    let newori = Builder::orient(eor.orientation(), rev_ori1);
                    eor.set_orientation(newori);
                    if keep {
                        wes.add_element(&eor);
                    } else {
                        self.my_on_elem_map.add(&eor);
                    }
                    self.b1.my_source_shapes.add(&eor);
                }
            } else {
                // `cxx:468-471`
                self.b1_g_fill_edge_not_same_dom_wes(&eor, ls_class, g1, wes);
            }
        }
    }

    /// `GFillEdgeNotSameDomWES(EOR, LSclass, G1, WES)` (`cxx:478-532`).
    pub(crate) fn b1_g_fill_edge_not_same_dom_wes(
        &mut self,
        eor: &TopoShape,
        _ls_class: &[TopoShape],
        g1: &GTopo,
        wes: &mut WireEdgeSet,
    ) {
        let (tb1, _tb2) = g1.states_on();
        let rev_ori1 = g1.is_to_reverse1();

        let ori_e = eor.orientation();
        let newori_e = Builder::orient(ori_e, rev_ori1);

        // 1) the split parts of the edge in state `TB1` (`cxx:493-502`).
        for new_e in self.parts_in_state(eor, tb1) {
            let mut new_e = new_e;
            new_e.set_orientation(newori_e);
            wes.add_start_element(&new_e);
            self.b1.my_source_shapes.add(&new_e);
        }

        // 2) the ON parts, kept or not (`cxx:504-527`).
        for new_e in self.parts_in_state(eor, State::On) {
            let mut new_e = new_e;
            new_e.set_orientation(newori_e);

            if self.b1.my_splits_on_to_keep.contains(&new_e) {
                wes.add_start_element(&new_e);
                continue;
            }

            // Degenerated edges are kept here: `FillONPartsWES` cannot
            // process them (`cxx:517-521`).
            let degenerated = Edge::wrap(new_e.clone())
                .map(|e| BRepTool::is_degenerated(&e))
                .unwrap_or(false);
            if degenerated {
                wes.add_start_element(&new_e);
                self.b1.my_source_shapes.add(&new_e);
            }
            if self.my_process_on {
                self.my_on_elem_map.add(&new_e);
                self.b1.my_source_shapes.add(&new_e);
            }
        }
    }
}
