//! The `TopOpeBRepBuild_Builder1` state block.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_Builder1.hxx:186-203`
//! (OCCT 8.0.0).
//!
//! OCCT keeps these maps in the subclass. The three members OCCT declares
//! `virtual` in `TopOpeBRepBuild_Builder.hxx` and overrides in
//! `TopOpeBRepBuild_Builder1` (`GFillSolidSFS` `Builder1.hxx:69-73`,
//! `GFillShellSFS` `:76-79`, `GWESMakeFaces` `:80-82`) are reached from the
//! base class bodies (`TopOpeBRepBuild_GridSS.cxx:433-494`,
//! `:596-645`, `TopOpeBRepBuild_makefaces.cxx:85-246`). Rust has no virtual
//! dispatch here because `Builder1` holds the base as a sub-object, so the
//! subclass state lives in [`super::Builder::b1`] and the base member checks
//! [`super::Builder::is_builder1`] before running the override — the same
//! branch OCCT resolves in the vtable.
//!
//! `myFSplits` / `myESplits` are declared a second time in `Builder1.hxx:186-189`,
//! shadowing the base members; they are kept here as well so a shape split
//! through `Builder1` and through `Builder` can never alias. Nothing in the
//! ported subclass code reads them (`Builder1.cxx` does not either), hence
//! the `dead_code` allowance on those two.

use std::cell::RefCell;
use std::collections::HashMap;
use std::rc::Rc;

use crate::shape::TopoShape;
use crate::topo_tools::ShapeMapKey;
use crate::topopebrep::ds::data_structure::{IndexedShapeData, IndexedShapeSet};

/// `class TopOpeBRepBuild_Builder1` private block (`hxx:186-203`).
#[derive(Clone, Default)]
pub(crate) struct Builder1State {
    /// `myFSplits` (`hxx:186-187`).
    #[allow(dead_code)]
    pub(crate) my_f_splits: HashMap<ShapeMapKey, Vec<TopoShape>>,
    /// `myESplits` (`hxx:188-189`).
    #[allow(dead_code)]
    pub(crate) my_e_splits: HashMap<ShapeMapKey, Vec<TopoShape>>,
    /// `mySameDomMap` (`hxx:192`): the same-domain faces already filled.
    pub(crate) my_same_dom_map: IndexedShapeSet,
    /// `mySDFaceToFill` (`hxx:193`).
    pub(crate) my_sd_face_to_fill: Option<TopoShape>,
    /// `myBaseFaceToFill` (`hxx:194`).
    pub(crate) my_base_face_to_fill: Option<TopoShape>,
    /// `myMapOfEdgeFaces` (`hxx:195-196`): every face of `myShape1` /
    /// `myShape2` under each edge, filled by `GFillSolidSFS` (`cxx:201-202`).
    pub(crate) my_map_of_edge_faces: IndexedShapeData<Vec<TopoShape>>,
    /// `myMapOfEdgeWithFaceState` (`hxx:197`).
    pub(crate) my_map_of_edge_with_face_state: HashMap<ShapeMapKey, bool>,
    /// `myProcessedPartsOut2d` (`hxx:198`).
    pub(crate) my_processed_parts_out_2d: IndexedShapeSet,
    /// `myProcessedPartsON2d` (`hxx:199`).
    pub(crate) my_processed_parts_on_2d: IndexedShapeSet,
    /// `mySplitsONtoKeep` (`hxx:200`).
    pub(crate) my_splits_on_to_keep: IndexedShapeSet,
    /// `mySourceShapes` (`hxx:201`).
    pub(crate) my_source_shapes: IndexedShapeSet,
    /// `myMapOfCorrect2dEdges` (`hxx:202-203`): the corrected edge for each
    /// edge `CorrectFace2d` moved.
    ///
    /// `TopOpeBRepBuild_CorrectFace2d` receives this map by reference
    /// (`CorrectFace2d.cxx:56-68`) and writes the copies into it while it
    /// walks the face, so the port shares it through `Rc<RefCell<_>>` and the
    /// caller's copy sees the moves.
    pub(crate) my_map_of_correct_2d_edges:
        Rc<RefCell<IndexedShapeData<TopoShape>>>,
}

impl Builder1State {
    /// The maps start empty (`Builder1.cxx:68-73`).
    pub(crate) fn new() -> Self {
        Self::default()
    }
}
