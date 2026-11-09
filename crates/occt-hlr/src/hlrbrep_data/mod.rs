//! `HLRBRep_Data` (`ModelingAlgorithms/TKHLR/HLRBRep/HLRBRep_Data.hxx:56-289`,
//! `.cxx:522-2682`, `.lxx:19-120`): the exact-HLR data structure - one
//! `HLRBRep_EdgeData` / `HLRBRep_FaceData` per topological edge / face plus the
//! projection bookkeeping the hiding pass reads.
//!
//! The class is split by OCCT responsibility:
//!
//! * this file - storage, the two constructors, the accessors of `.lxx`;
//! * [`init_bound_sort`] - `InitBoundSort` (`cxx:984-1019`);
//! * [`update`] - `Update(P)` (`cxx:612-980`);
//! * [`write`] - `Write(DS, dv, de, df)` (`cxx:548-608`);
//! * [`reject`] - the `TableauRejection` cache (`cxx:78-491`);
//! * [`edge_iter`] - `InitEdge` / `MoreEdge` / `NextEdge` / `Edge`
//!   (`cxx:1022-1221`);
//! * [`interference`] - `InitInterference` / `MoreInterference` /
//!   `NextInterference` / `SameVertex` / `LocalLEGeometry2D` /
//!   `LocalFEGeometry2D` (`cxx:1223-1578`, `:2597-2653`);
//! * [`rejected_point`] - `RejectedInterference` / `AboveInterference` /
//!   `RejectedPoint` (`cxx:1486-1530`, `:2330-2595`);
//! * [`edge_state`] - `EdgeState` / `HidingStartLevel` / `Compare`
//!   (`cxx:1580-1744`);
//! * [`classify`] - `Classify` / `SimplClassify` / `IsBadFace`
//!   (`cxx:1993-2328`, `:2655-2682`).
//!
//! # What this batch ports
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | `HLRBRep_Data(NV, NE, NF)` | `cxx:522-537` | [`HlrBRepData::new`] |
//! | `Destroy()` | `cxx:539-544` | [`HlrBRepData::destroy`] |
//! | `Write(DS, dv, de, df)` | `cxx:548-608` | [`HlrBRepData::write`] |
//! | `EDataArray` / `FDataArray` | `lxx:19-28` | [`HlrBRepData::e_data`] / [`HlrBRepData::f_data`] (+ `_mut`) |
//! | `Tolerance(float)` | `lxx:33-37` | [`HlrBRepData::set_tolerance`] |
//! | `Tolerance()` | `lxx:40-44` | [`HlrBRepData::tolerance`] |
//! | `Projector()` | `lxx:47-51` | [`HlrBRepData::projector`] (+ `_mut`) |
//! | `NbVertices` / `NbEdges` / `NbFaces` | `lxx:54-72` | the same names |
//! | `EdgeMap` / `FaceMap` | `lxx:75-86` | [`HlrBRepData::edge_map`] / [`HlrBRepData::face_map`] |
//! | `HidingTheFace` / `SimpleHidingFace` | `lxx:89-98` | the same names |
//! | `EdgeOfTheHidingFace` | `lxx:115-120` | [`HlrBRepData::edge_of_the_hiding_face`] |
//! | `InitBoundSort` | `cxx:984-1019` | [`HlrBRepData::init_bound_sort`] |
//! | `Update(P)` | `cxx:612-980` | [`HlrBRepData::update`] |
//! | `TableauRejection` | `cxx:78-491` | [`reject::TableauRejection`] |
//! | `InitEdge` / `MoreEdge` / `NextEdge` / `Edge` | `cxx:1022-1221` | [`HlrBRepData::init_edge`] / `more_edge` / `next_edge` / `edge` |
//! | `InitInterference` / `MoreInterference` / `NextInterference` | `cxx:1223-1484`, `lxx:103-106` | the same names |
//! | `SameVertex` | `cxx:2597-2653` | `same_vertex` (private) |
//! | `LocalLEGeometry2D` / `LocalFEGeometry2D` | `cxx:1532-1578` | the same names |
//! | `RejectedInterference` / `AboveInterference` / `RejectedPoint` | `cxx:1486-1530`, `:2330-2595` | the same names |
//! | `EdgeState` / `HidingStartLevel` / `Compare` | `cxx:1580-1744` | `edge_state` / the same names |
//! | `Classify` / `SimplClassify` / `IsBadFace` | `cxx:1993-2328`, `:2655-2682` | the same names |
//!
//! # UNPORTED (each needs an algorithm that has no port yet)
//!
//! * nothing: the batch after T-419 filled `RejectedInterference`,
//!   `RejectedPoint`, `AboveInterference`, `EdgeState`, `HidingStartLevel`,
//!   `Compare`, `Classify`, `SimplClassify` and `IsBadFace`. `HLRBRep_Hider`
//!   (`HLRBRep_Hider.cxx`) is the remaining consumer of all of them.
//! * `OrientOutLine` (`cxx:1746-1879`) and `OrientOthEdge` (`cxx:1881-1991`)
//!   are ported in [`orient`]; `Update` reaches them only through
//!   `FD.WithOutL() && !FD.Side()`.
//!
//! # Translation notes
//!
//! * `NCollection_Array1<T>(0, NE)` is a `Vec` of `NE + 1` slots with the
//!   OCCT 1-based indices kept: slot 0 is the spare `Write` uses through its
//!   `ChangeValue(0)` + `ed++` idiom (`cxx:557-563`), slots `1..=NE` are the
//!   entities. [`HlrBRepData::e_data`] takes the OCCT index.
//! * `myEMap` / `myFMap` are `NCollection_IndexedMap` (`FindKey(i)` is a
//!   1-based lookup, `Add` appends). The port keeps a `Vec<TopoShape>` in
//!   insertion order; the key set is not queried anywhere in this batch.
//! * `HLRBRep_CLProps` (`cxx:530-531`, `HLRBRep_CLProps.hxx:24-31`) is the
//!   curve-side `GeomLProp_CLProps` specialisation; the only use inside
//!   `Update` is the surface-side holder, and `HLRBRep_SLProps`
//!   (`HLRBRep_SLProps.hxx:29-33`) is ported in [`crate::hlrbrep_sl_props`].
//!   `mySLProps` therefore lives on the struct (`hxx:238`) instead of being
//!   recomputed inline.
//! * `OCC_CATCH_SIGNALS` / `OCCT_DEBUG` printing has no port counterpart.

use std::sync::Arc;

use occt_topo::abs::Orientation;
use occt_topo::brep_top_adaptor_topol_tool::TopolTool;
use occt_topo::geom_abs_curve_type::GeomAbsCurveType;
use occt_topo::geom_bnd_lib_surface3d::GeomAbsSurfaceType;
use crate::hlr_algo::{HlrAlgoInterference, HlrAlgoProjector, MinMaxIndices};
use crate::hlrbrep_edge_data::HlrBRepEdgeData;
use crate::hlrbrep_face_data::HlrBRepFaceData;
use crate::hlrbrep_face_iterator::HlrBRepFaceIterator;
use crate::hlrbrep_intersector::HlrBRepIntersector;
use crate::hlrbrep_sl_props::HlrBRepSlProps;
use occt_topo::shape::TopoShape;

pub(crate) mod classify;
pub(crate) mod edge_iter;
pub(crate) mod edge_state;
pub(crate) mod interference;
mod init_bound_sort;
mod orient;
pub(crate) mod reject;
pub(crate) mod rejected_point;
mod update;
mod write;

use reject::TableauRejection;

/// `MASK` used by OCCT's 16-bit packed min/max sign test (`cxx:996`, the
/// `& 0x80008000` pairs of `InitBoundSort` and `HLRBRep_InternalAlgo::Hide`).
pub(crate) const MINMAX_SIGN_MASK: i32 = 0x8000_8000u32 as i32;

/// `HLRBRep_Data`.
#[derive(Clone)]
pub struct HlrBRepData {
    /// `int myNbVertices` (`hxx:225`).
    my_nb_vertices: i32,
    /// `int myNbEdges` (`hxx:226`).
    my_nb_edges: i32,
    /// `int myNbFaces` (`hxx:227`).
    my_nb_faces: i32,
    /// `NCollection_IndexedMap<TopoDS_Shape, ...> myEMap` (`hxx:228`).
    my_e_map: Vec<TopoShape>,
    /// `NCollection_IndexedMap<TopoDS_Shape, ...> myFMap` (`hxx:229`).
    my_f_map: Vec<TopoShape>,
    /// `NCollection_Array1<HLRBRep_EdgeData> myEData` (`hxx:230`), index 0 spare.
    my_e_data: Vec<HlrBRepEdgeData>,
    /// `NCollection_Array1<HLRBRep_FaceData> myFData` (`hxx:231`), index 0 spare.
    my_f_data: Vec<HlrBRepFaceData>,
    /// `NCollection_Array1<int> myEdgeIndices` (`hxx:232`), index 0 unused.
    my_edge_indices: Vec<i32>,
    /// `float myToler` (`hxx:233`).
    my_toler: f32,
    /// `HLRAlgo_Projector myProj` (`hxx:234`).
    my_proj: Arc<HlrAlgoProjector>,
    /// `double myBigSize` (`hxx:239`).
    my_big_size: f64,
    /// `double myDeca[16]` (`hxx:253`).
    my_deca: [f64; 16],
    /// `double mySurD[16]` (`hxx:254`).
    my_sur_d: [f64; 16],
    /// `int myCurSortEd` (`hxx:255`).
    my_cur_sort_ed: i32,
    /// `int myNbrSortEd` (`hxx:256`).
    my_nbr_sort_ed: i32,
    /// `int myHideCount` (`hxx:252`).
    my_hide_count: i32,
    /// `bool iFaceTest` (`hxx:250`): "the current face hides the edge", set by
    /// `InitEdge` (`cxx:1035`).
    i_face_test: bool,
    /// `bool iFaceSimp` (`hxx:248`): "the current hiding face is simple", set
    /// by `InitEdge` (`cxx:1036`) and cleared by `RejectedPoint` (`cxx:2402`).
    i_face_simp: bool,
    /// `HLRBRep_Surface* iFaceGeom` (`hxx:245`): the face `Update` is working
    /// on, as an index into `my_f_data` (OCCT keeps the pointer).
    i_face_geom: Option<i32>,
    /// `HLRBRep_SLProps mySLProps` (`hxx:238`), constructed as
    /// `mySLProps(2, Epsilon(1.))` (`cxx:532`).
    my_sl_props: HlrBRepSlProps,
    /// `HLRBRep_FaceIterator myFaceItr1` (`hxx:240`): the wire iterator of the
    /// `InitInterference` / `NextInterference` pass.
    my_face_itr1: HlrBRepFaceIterator,
    /// `HLRBRep_FaceIterator myFaceItr2` (`hxx:241`): the wire iterator of the
    /// `InitEdge` / `MoreEdge` / `NextEdge` pass over a non-simple face.
    my_face_itr2: HlrBRepFaceIterator,
    /// `int iFace` (`hxx:242`): the hiding face, 1-based into `my_f_data`.
    i_face: i32,
    /// `HLRAlgo_EdgesBlock::MinMaxIndices* iFaceMinMax` (`hxx:246`): a copy of
    /// `FaceData.Wires().MinMax()` (`cxx:1034`).
    i_face_min_max: Option<MinMaxIndices>,
    /// `GeomAbs_SurfaceType iFaceType` (`hxx:247`): OCCT writes it in
    /// `InitEdge` (`cxx:1036`) and never reads it.
    #[allow(dead_code)]
    i_face_type: GeomAbsSurfaceType,
    /// `bool iFaceBack` (`hxx:248`), the `IsAbove` orientation flag.
    i_face_back: bool,
    /// `bool iFaceSmpl` (`hxx:249`): `!FaceData.Cut()`, decided at the start of
    /// a non-simple face and stored back when its wires are exhausted.
    i_face_smpl: bool,
    /// `int myLE` (`hxx:257`): the current hiding edge, 1-based.
    my_le: i32,
    /// `bool myLEOutLine` (`hxx:258`).
    my_le_out_line: bool,
    /// `bool myLEInternal` (`hxx:259`), read by `SameVertex`.
    my_le_internal: bool,
    /// `bool myLEDouble` (`hxx:260`).
    my_le_double: bool,
    /// `bool myLEIsoLine` (`hxx:261`), written by `InitEdge`; only the OCCT
    /// `OCCT_DEBUG` output and `HLRBRep_Hider` read it.
    #[allow(dead_code)]
    my_le_iso_line: bool,
    /// `HLRAlgo_EdgesBlock::MinMaxIndices* myLEMinMax` (`hxx:264`).
    my_le_min_max: Option<MinMaxIndices>,
    /// `GeomAbs_CurveType myLEType` (`hxx:265`).
    my_le_type: GeomAbsCurveType,
    /// `float myLETol` (`hxx:266`).
    my_le_tol: f32,
    /// `int myFE` (`hxx:267`): the current edge of the hiding face, 1-based.
    my_fe: i32,
    /// `TopAbs_Orientation myFEOri` (`hxx:268`).
    my_fe_ori: Orientation,
    /// `bool myFEOutLine` (`hxx:269`).
    my_fe_out_line: bool,
    /// `bool myFEInternal` (`hxx:270`): a face edge marked internal, which
    /// forces `Level = 2` in `RejectedPoint` (`cxx:2569-2572`).
    my_fe_internal: bool,
    /// `bool myFEDouble` (`hxx:271`).
    my_fe_double: bool,
    /// `GeomAbs_CurveType myFEType` (`hxx:274`).
    my_fe_type: GeomAbsCurveType,
    /// `float myFETol` (`hxx:275`): OCCT declares it next to `myLETol` but no
    /// algorithm reads it (`RejectedPoint` publishes `myLETol` only).
    #[allow(dead_code)]
    my_fe_tol: f32,
    /// `HLRBRep_Intersector myIntersector` (`hxx:276`).
    my_intersector: HlrBRepIntersector,
    /// `occ::handle<BRepTopAdaptor_TopolTool> myClassifier` (`hxx:277`), bound
    /// by `InitEdge` from the `MST` map and read by `Classify`
    /// (`cxx:2254`) through [`HlrBRepData::classify`].
    my_classifier: Option<Arc<TopolTool>>,
    /// `bool mySameVertex` (`hxx:278`).
    my_same_vertex: bool,
    /// `bool myIntersected` (`hxx:279`).
    my_intersected: bool,
    /// `int myNbPoints` (`hxx:280`): the points of the current pair.
    my_nb_points: i32,
    /// `int myNbSegments` (`hxx:281`): the segments of the current pair.
    my_nb_segments: i32,
    /// `int iInterf` (`hxx:282`): 1-based interference cursor.
    i_interf: i32,
    /// `TableauRejection* myReject` (`hxx:285`), built as
    /// `new TableauRejection` + `SetDim(myNbEdges)` (`cxx:535-536`).
    my_reject: TableauRejection,
    /// `HLRAlgo_Interference myIntf` (`hxx:282`): the interference
    /// `RejectedPoint` fills and the hider reads back through
    /// [`HlrBRepData::interference`].
    my_intf: HlrAlgoInterference,
    /// `bool myAboveIntf` (`hxx:283`): the last `RejectedPoint` found the point
    /// above the face (`cxx:2363-2368`).
    my_above_intf: bool,
}

impl HlrBRepData {
    /// `HLRBRep_Data(NV, NE, NF)` (`cxx:522-537`): the three array sizes, the
    /// `1e-5` rejection tolerance, the rejection table sized to `NE`
    /// (`myReject = new TableauRejection; myReject->SetDim(myNbEdges)`,
    /// `cxx:535-536`) and the two property holders.
    pub fn new(nv: i32, ne: i32, nf: i32) -> Self {
        let ne_len = if ne < 0 { 0 } else { ne as usize + 1 };
        let nf_len = if nf < 0 { 0 } else { nf as usize + 1 };
        let mut my_reject = TableauRejection::new();
        my_reject.set_dim(if ne < 0 { 0 } else { ne as usize });
        Self {
            my_nb_vertices: nv,
            my_nb_edges: ne,
            my_nb_faces: nf,
            my_e_map: Vec::new(),
            my_f_map: Vec::new(),
            my_e_data: vec![HlrBRepEdgeData::new(); ne_len],
            my_f_data: vec![HlrBRepFaceData::new(); nf_len],
            my_edge_indices: vec![0; ne_len],
            // `myToler((float)1e-5)` (`cxx:529`).
            my_toler: 1e-5f32,
            my_proj: Arc::new(HlrAlgoProjector::new()),
            my_big_size: 0.0,
            my_deca: [0.0; 16],
            my_sur_d: [0.0; 16],
            my_cur_sort_ed: 0,
            my_nbr_sort_ed: 0,
            my_hide_count: 0,
            i_face_test: false,
            i_face_simp: false,
            i_face_geom: None,
            my_sl_props: HlrBRepSlProps::for_hlr_data(),
            my_face_itr1: HlrBRepFaceIterator::new(),
            my_face_itr2: HlrBRepFaceIterator::new(),
            i_face: 0,
            i_face_min_max: None,
            i_face_type: GeomAbsSurfaceType::Plane,
            i_face_back: false,
            i_face_smpl: false,
            my_le: 0,
            my_le_out_line: false,
            my_le_internal: false,
            my_le_double: false,
            my_le_iso_line: false,
            my_le_min_max: None,
            my_le_type: GeomAbsCurveType::Line,
            my_le_tol: 0.0,
            my_fe: 0,
            my_fe_ori: Orientation::Forward,
            my_fe_out_line: false,
            my_fe_internal: false,
            my_fe_double: false,
            my_fe_type: GeomAbsCurveType::Line,
            my_fe_tol: 0.0,
            my_intersector: HlrBRepIntersector::new(),
            my_classifier: None,
            my_same_vertex: false,
            my_intersected: false,
            my_nb_points: 0,
            my_nb_segments: 0,
            i_interf: 0,
            my_reject,
            my_intf: HlrAlgoInterference::new(),
            my_above_intf: false,
        }
    }

    /// `Destroy()` (`cxx:539-544`): `myReject->Destroy()` and the deletion of
    /// the table; the owned arrays are dropped with the struct.
    pub fn destroy(&mut self) {
        self.my_cur_sort_ed = 0;
        self.my_nbr_sort_ed = 0;
        self.my_reject.destroy();
    }

    /// `myBigSize` (`hxx:239`): the largest of the sixteen global extents,
    /// written by [`HlrBRepData::update`].
    pub fn big_size(&self) -> f64 {
        self.my_big_size
    }

    /// `myDeca[16]` (`hxx:253`), index 0..15.
    pub fn deca(&self) -> &[f64; 16] {
        &self.my_deca
    }

    /// `mySurD[16]` (`hxx:254`), index 0..15.
    pub fn sur_d(&self) -> &[f64; 16] {
        &self.my_sur_d
    }

    /// `myHideCount` (`hxx:252`), the level counter `InitEdge` raises by two per
    /// hiding face (`cxx:1039-1040`).
    pub fn hide_count(&self) -> i32 {
        self.my_hide_count
    }

    /// The mutable `myHideCount`.
    pub fn set_hide_count(&mut self, c: i32) {
        self.my_hide_count = c;
    }

    /// `myNbrSortEd` (`hxx:256`): the number of edges `InitBoundSort` kept.
    pub fn nb_sort_edges(&self) -> i32 {
        self.my_nbr_sort_ed
    }

    /// `myEdgeIndices(i)` (`cxx:1015` / `:1150`), OCCT 1-based.
    pub fn sort_edge_index(&self, i: i32) -> i32 {
        self.my_edge_indices[i as usize]
    }

    /// `myCurSortEd` (`hxx:255`), the cursor of `NextEdge` (`cxx:1150-1160`).
    pub fn cur_sort_edge(&self) -> i32 {
        self.my_cur_sort_ed
    }

    /// The mutable `myCurSortEd`.
    pub fn set_cur_sort_edge(&mut self, i: i32) {
        self.my_cur_sort_ed = i;
    }

    /// `EDataArray()` (`lxx:19-22`) element access with the OCCT 1-based index.
    pub fn e_data(&self, i: i32) -> &HlrBRepEdgeData {
        &self.my_e_data[i as usize]
    }

    /// The mutable `EDataArray()` (`cxx:550`, `:589`).
    pub fn e_data_mut(&mut self, i: i32) -> &mut HlrBRepEdgeData {
        &mut self.my_e_data[i as usize]
    }

    /// `FDataArray()` (`lxx:24-28`) element access with the OCCT 1-based index.
    pub fn f_data(&self, i: i32) -> &HlrBRepFaceData {
        &self.my_f_data[i as usize]
    }

    /// The mutable `FDataArray()` (`cxx:546`, `:588`).
    pub fn f_data_mut(&mut self, i: i32) -> &mut HlrBRepFaceData {
        &mut self.my_f_data[i as usize]
    }

    /// `Tolerance(tol)` (`lxx:33-37`).
    pub fn set_tolerance(&mut self, tol: f32) {
        self.my_toler = tol;
    }

    /// `Tolerance()` (`lxx:40-44`).
    pub fn tolerance(&self) -> f32 {
        self.my_toler
    }

    /// `Projector()` (`lxx:47-51`).
    pub fn projector(&self) -> &Arc<HlrAlgoProjector> {
        &self.my_proj
    }

    /// The non-const `Projector()` of `lxx:47-51`.
    pub fn projector_mut(&mut self) -> &mut Arc<HlrAlgoProjector> {
        &mut self.my_proj
    }

    /// `mySLProps` (`hxx:238`).
    pub(crate) fn my_sl_props(&self) -> &HlrBRepSlProps {
        &self.my_sl_props
    }

    /// The mutable `mySLProps`.
    pub(crate) fn my_sl_props_mut(&mut self) -> &mut HlrBRepSlProps {
        &mut self.my_sl_props
    }

    /// `NbVertices()` (`lxx:54-57`).
    pub fn nb_vertices(&self) -> i32 {
        self.my_nb_vertices
    }

    /// `NbEdges()` (`lxx:60-63`).
    pub fn nb_edges(&self) -> i32 {
        self.my_nb_edges
    }

    /// `NbFaces()` (`lxx:66-70`).
    pub fn nb_faces(&self) -> i32 {
        self.my_nb_faces
    }

    /// `EdgeMap()` (`lxx:75-78`), insertion order.
    pub fn edge_map(&self) -> &[TopoShape] {
        &self.my_e_map
    }

    /// `EdgeMap().Add(E)` (`HLRBRep_ShapeToHLR.cxx:168`).
    pub fn edge_map_push(&mut self, e: &TopoShape) {
        self.my_e_map.push(e.clone());
    }

    /// `FaceMap()` (`lxx:81-86`), insertion order.
    pub fn face_map(&self) -> &[TopoShape] {
        &self.my_f_map
    }

    /// `FaceMap().Add(F)` (`HLRBRep_ShapeToHLR.cxx:239`).
    pub fn face_map_push(&mut self, f: &TopoShape) {
        self.my_f_map.push(f.clone());
    }

    /// `HidingTheFace()` (`lxx:95-98`): `iFaceTest`, set by `InitEdge`
    /// (`cxx:1035`).
    pub fn hiding_the_face(&self) -> bool {
        self.i_face_test
    }

    /// `SimpleHidingFace()` (`lxx:89-93`): `iFaceSimp`, set by `InitEdge`
    /// (`cxx:1036`).
    pub fn simple_hiding_face(&self) -> bool {
        self.i_face_simp
    }

    /// `EdgeOfTheHidingFace(E, ED)` (`lxx:115-120`): `ED.HideCount() ==
    /// myHideCount - 1`.
    pub fn edge_of_the_hiding_face(&self, _e: i32, ed: &HlrBRepEdgeData) -> bool {
        ed.hide_count() == self.my_hide_count - 1
    }
}

impl Default for HlrBRepData {
    /// `HLRBRep_Data(0, 0, 0)`: the empty structure `HLRBRep_InternalAlgo`
    /// builds when loading a shape fails (`HLRBRep_InternalAlgo.cxx:109`).
    fn default() -> Self {
        Self::new(0, 0, 0)
    }
}
