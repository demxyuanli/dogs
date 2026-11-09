//! `HLRBRep_ShapeBounds` (`ModelingAlgorithms/TKHLR/HLRBRep/HLRBRep_ShapeBounds.hxx:31-91`,
//! `.cxx:23-98`, `.lxx:19-68`): one shape of the exact HLR data structure
//! together with the index ranges its vertices, edges and faces occupy inside
//! `HLRBRep_Data`.
//!
//! `myShapeData` is an `occ::handle<Standard_Transient>` in OCCT. The slot is
//! written by `HLRBRep_InternalAlgo::ShapeData(I, SData)`
//! (`HLRBRep_InternalAlgo.cxx:236-242`) and only ever forwarded by this class
//! (`cxx:23-44`, `cxx:46-65`, `HLRBRep_InternalAlgo.cxx:115`); the object stored
//! in practice is an `HLRBRep_Data`, which is not ported yet. The field is
//! therefore a type-erased `Arc<dyn Any + Send + Sync>` - the
//! `Standard_Transient` contract without inventing a concrete type.
//!
//! The `int&` out-parameters of `Sizes` / `Bounds` (`cxx:79-97`) become a
//! returned tuple; the OCCT defaults (`lxx:19-30`) are preserved.

use std::any::Any;
use std::sync::Arc;

use crate::hlr_algo::MinMaxIndices;
use crate::hlrtopobrep::HlrTopoBRepOutLiner;

/// `HLRBRep_ShapeBounds`.
#[derive(Clone)]
pub struct HlrBRepShapeBounds {
    /// `occ::handle<HLRTopoBRep_OutLiner> myShape` (`hxx:83`).
    my_shape: Option<Arc<HlrTopoBRepOutLiner>>,
    /// `occ::handle<Standard_Transient> myShapeData` (`hxx:84`).
    my_shape_data: Option<Arc<dyn Any + Send + Sync>>,
    /// `int myNbIso` (`hxx:85`).
    my_nb_iso: i32,
    /// `int myVertStart` (`hxx:86`).
    my_vert_start: i32,
    /// `int myVertEnd` (`hxx:87`).
    my_vert_end: i32,
    /// `int myEdgeStart` (`hxx:88`).
    my_edge_start: i32,
    /// `int myEdgeEnd` (`hxx:89`).
    my_edge_end: i32,
    /// `int myFaceStart` (`hxx:90`).
    my_face_start: i32,
    /// `int myFaceEnd` (`hxx:91`).
    my_face_end: i32,
    /// `HLRAlgo_EdgesBlock::MinMaxIndices myMinMax` (`hxx:92`).
    my_min_max: MinMaxIndices,
}

impl Default for HlrBRepShapeBounds {
    /// `HLRBRep_ShapeBounds()` (`lxx:19-30`): a null shape / null data handle
    /// and every bound at zero. `myMinMax` is not initialised by OCCT.
    fn default() -> Self {
        Self::new()
    }
}

impl HlrBRepShapeBounds {
    /// `HLRBRep_ShapeBounds()` (`lxx:19-30`).
    pub fn new() -> Self {
        Self {
            my_shape: None,
            my_shape_data: None,
            my_nb_iso: 0,
            my_vert_start: 0,
            my_vert_end: 0,
            my_edge_start: 0,
            my_edge_end: 0,
            my_face_start: 0,
            my_face_end: 0,
            my_min_max: MinMaxIndices::default(),
        }
    }

    /// `HLRBRep_ShapeBounds(S, SData, nbIso, V1, V2, E1, E2, F1, F2)`
    /// (`cxx:23-44`).
    #[allow(clippy::too_many_arguments)]
    pub fn with_shape_data(
        s: Arc<HlrTopoBRepOutLiner>,
        s_data: Option<Arc<dyn Any + Send + Sync>>,
        nb_iso: i32,
        v1: i32,
        v2: i32,
        e1: i32,
        e2: i32,
        f1: i32,
        f2: i32,
    ) -> Self {
        Self {
            my_shape: Some(s),
            my_shape_data: s_data,
            my_nb_iso: nb_iso,
            my_vert_start: v1,
            my_vert_end: v2,
            my_edge_start: e1,
            my_edge_end: e2,
            my_face_start: f1,
            my_face_end: f2,
            my_min_max: MinMaxIndices::default(),
        }
    }

    /// `HLRBRep_ShapeBounds(S, nbIso, V1, V2, E1, E2, F1, F2)` (`cxx:46-65`):
    /// the same construction leaving `myShapeData` null.
    #[allow(clippy::too_many_arguments)]
    pub fn new_bounded(
        s: Arc<HlrTopoBRepOutLiner>,
        nb_iso: i32,
        v1: i32,
        v2: i32,
        e1: i32,
        e2: i32,
        f1: i32,
        f2: i32,
    ) -> Self {
        Self::with_shape_data(s, None, nb_iso, v1, v2, e1, e2, f1, f2)
    }

    /// `Translate(NV, NE, NF)` (`cxx:67-77`): shifts all six bounds.
    pub fn translate(&mut self, nv: i32, ne: i32, nf: i32) {
        self.my_vert_start += nv;
        self.my_vert_end += nv;
        self.my_edge_start += ne;
        self.my_edge_end += ne;
        self.my_face_start += nf;
        self.my_face_end += nf;
    }

    /// `Sizes(NV, NE, NF)` (`cxx:79-86`): `end + 1 - start` per entity kind.
    pub fn sizes(&self) -> (i32, i32, i32) {
        (
            self.my_vert_end + 1 - self.my_vert_start,
            self.my_edge_end + 1 - self.my_edge_start,
            self.my_face_end + 1 - self.my_face_start,
        )
    }

    /// `Bounds(V1, V2, E1, E2, F1, F2)` (`cxx:88-97`).
    pub fn bounds(&self) -> (i32, i32, i32, i32, i32, i32) {
        (
            self.my_vert_start,
            self.my_vert_end,
            self.my_edge_start,
            self.my_edge_end,
            self.my_face_start,
            self.my_face_end,
        )
    }

    /// `Shape(S)` (`lxx:32-35`).
    pub fn set_shape(&mut self, s: Arc<HlrTopoBRepOutLiner>) {
        self.my_shape = Some(s);
    }

    /// `Shape()` (`lxx:37-42`).
    pub fn shape(&self) -> Option<&Arc<HlrTopoBRepOutLiner>> {
        self.my_shape.as_ref()
    }

    /// The non-const `Shape()` of `lxx:37-42`, used by
    /// `HLRBRep_InternalAlgo::Update` (`cxx:95`) to fill the outliner in place
    /// through its handle.
    pub fn take_shape(&mut self) -> Option<Arc<HlrTopoBRepOutLiner>> {
        self.my_shape.take()
    }

    /// `ShapeData(SD)` (`lxx:44-49`).
    pub fn set_shape_data(&mut self, sd: Option<Arc<dyn Any + Send + Sync>>) {
        self.my_shape_data = sd;
    }

    /// `ShapeData()` (`lxx:51-56`).
    pub fn shape_data(&self) -> Option<&Arc<dyn Any + Send + Sync>> {
        self.my_shape_data.as_ref()
    }

    /// `NbOfIso(nbIso)` (`lxx:58-63`).
    pub fn set_nb_of_iso(&mut self, nb_iso: i32) {
        self.my_nb_iso = nb_iso;
    }

    /// `NbOfIso()` (`lxx:65-70`).
    pub fn nb_of_iso(&self) -> i32 {
        self.my_nb_iso
    }

    /// `UpdateMinMax(theTotMinMax)` (`hxx:73-76`).
    pub fn update_min_max(&mut self, tot_min_max: &MinMaxIndices) {
        self.my_min_max = *tot_min_max;
    }

    /// `MinMax()` (`hxx:78`).
    pub fn min_max(&self) -> &MinMaxIndices {
        &self.my_min_max
    }

    /// The non-const `MinMax()` of `hxx:78`.
    pub fn min_max_mut(&mut self) -> &mut MinMaxIndices {
        &mut self.my_min_max
    }
}
