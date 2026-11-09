//! `BRepOffset_MakeOffset` (`BRepOffset_MakeOffset.hxx` 288 lines /
//! `.cxx` 5,660 lines / `_1.cxx` 9,49x lines): the driver that turns an initial
//! shape into its offset shape.
//!
//! Ported here so far (control flow taken from
//! `src/ModelingAlgorithms/TKOffset/BRepOffset/BRepOffset_MakeOffset.cxx`):
//!
//! | OCCT | lines | Rust |
//! |---|---|---|
//! | `BRepOffset_MakeOffset()` | `:606-609` | [`BRepOffsetMakeOffset::new`] |
//! | `Initialize(...)` | `:643-669` | [`BRepOffsetMakeOffset::initialize`] |
//! | `Clear()` | `:672-688` | [`BRepOffsetMakeOffset::clear`] |
//! | `AllowLinearization(theIsAllowed)` | `:691-694` | [`BRepOffsetMakeOffset::allow_linearization`] |
//! | `AddFace(F)` | `:698-703` | [`BRepOffsetMakeOffset::add_face`] |
//! | `SetOffsetOnFace(F, Off)` | `:706-708` | [`BRepOffsetMakeOffset::set_offset_on_face`] |
//! | `SetFaces()` | `:795-815` | `face_comp.rs` |
//! | `SetFacesWithOffset()` | `:818-833` | `face_comp.rs` |
//! | `BuildFaceComp()` | `:1873-1890` | `face_comp.rs` |
//! | `ReplaceRoots()` | `:1825-1867` | `face_comp.rs` |
//! | `IsPlanar()` | `:5447-5552` | `is_planar.rs` |
//! | file-local `RemoveSeamAndDegeneratedEdges(theFace, theOldFace)` | `:5556-5633` | `is_planar.rs` |
//! | `MakeOffsetFaces(theMapSF, theRange)` | `:1202-1275` | `make_offset_faces.rs` |
//! | `UpdateFaceOffset()` | `:2342-2394` | `make_offset_faces.rs` |
//! | `BuildOffsetByArc(theRange)` | `:1893-2152` | `build_offset_by_arc.rs` |
//! | `ToContext(MapSF)` | `:2171-2340` | `build_offset_by_arc.rs` |
//! | `BuildOffsetByInter(theRange)` | `:1277-1823` | `build_offset_by_inter.rs` |
//! | file-local `TrimEdges(...)` | `:4752-4918` | `intersect_trim.rs` |
//! | file-local `TrimEdge(NE, ...)` | `:4920-5033` | `intersect_trim.rs` |
//! | `GetEnlargedFaces(...)` | `:5036-5056` | `intersect_trim.rs` |
//! | file-local `RemoveCorks(S, Faces)` | `:712-738` | `shape_prep.rs` |
//! | file-local `IsConnectedShell(S)` | `:742-750` | `shape_prep.rs` |
//! | file-local `MakeList(...)` | `:754-775` | `shape_prep.rs` |
//! | file-local `EvalMax(S, Tol)` | `:779-791` | `shape_prep.rs` |
//! | file-local `IsSolid(S)` / `AppendToList(...)` | `:5637-5659` | `shape_prep.rs` |
//! | `CheckInputData(theRange)` | `:4384-4498` | `input_checks.rs` |
//! | file-local `checkSinglePoint(...)` | `:4583-4622` | `input_checks.rs` |
//! | `IsDone()` / `Shape()` / `Error()` | `:1181-1197` | accessors |
//! | `OffsetFacesFromShapes()` / `OffsetEdgesFromShapes()` / `GetJoinType()` / `ClosingFaces()` | `:3929-3956` | accessors |
//! | `GetBadShape()` | `:4505-4507` | [`BRepOffsetMakeOffset::get_bad_shape`] |
//! | `InitShape()` | `hxx:98` | [`BRepOffsetMakeOffset::init_shape`] |
//! | `Intersection3D(Inter, ...)` | `:2855-2939` | `intersection.rs` |
//! | `Intersection2D(Modif, NewEdges, ...)` | `:2941-2989` | `intersection.rs` |
//! | `MakeLoops(Modif, ...)` | `:2991-3057` | `loop_rebuild.rs` |
//! | `MakeFaces(Modif, ...)` | `:3059-3146` | `loop_rebuild.rs` |
//! | `MakeShells(theRange)` | `:3691-3792` | `shells.rs` |
//! | `SelectShells()` | `:3896-3927` | `shells.rs` |
//! | `MakeSolid(theRange)` | `:3794-3894` | `solids.rs` |
//! | file-local `UpdateInitOffset(...)` | `:3110-3147` | `solids.rs` |
//! | file-local `CorrectSolid(...)` | `:4304-4382` | `solids.rs` |
//! | `EncodeRegularity()` | `:3962-4194` | `regularity.rs` |
//! | `RemoveInternalEdges()` | `:4512-4574` | `internal_edges.rs` |
//! | file-local `RemoveShapes(S, LS)` | `:4628-4642` | `internal_edges.rs` |
//! | file-local `ComputeMaxDist(Pln, C, f, l)` | `:4165-4192` | `tolerance_update.rs` |
//! | `UpdateTolerance(S, Faces, InitShape)` | `:4196-4302` | `tolerance_update.rs` |
//! | `MakeOffsetShape(theRange)` | `:837-1080` | `make_offset_shape.rs` |
//! | `CorrectConicalFaces()` | `:2396-2853` | `conical_faces.rs` (collection `:2415-2466`, rebuild `:2467-2760`, re-assembly `:2762-2851`) |
//! | `GetEdgePoints(E, F, f, m, l)` | `:515-530` | `conical_faces.rs` |
//! | `MakeMissingWalls(theRange)` | `:3148-3687` | `missing_walls.rs` (contours `:3151-3160`, per-edge wall `:3167-3685`, arc corner `:3619-3669`) |
//! | `MakeThickSolid(theRange)` | `:1083-1177` | `thick_solid.rs` |
//! | `Generated(theS)` | `:5231-5330` | `history.rs` (vertex `:5237-5266`, edge `:5269-5274`, face `:5276-5299`, solid `:5301-5306`, `myResMap` filter `:5312-5327`) |
//! | `Modified(theShape)` | `:5334-5365` | `history.rs` |
//! | `IsDeleted(theS)` | `:5369-5382` | `history.rs` |
//! | `SelfInter(Modif)` | `:2154-2167` | `self_inter.rs` (throws `Standard_NotImplemented`) |
//! | file-local `FillContours(...)` | `:536-599` | `missing_walls.rs` |
//! | file-local `OrientationOfEdgeInFace(E, F)` | `:326-343` | `missing_walls.rs` |
//! | `MakeShells(theRange)` | `:3691-3790` | `shells.rs` |
//! | `SelectShells()` | `:3896-3927` | `shells.rs` |
//! | file-local `BuildShellsCompleteInter(...)` | `:5064-5227` | `shells_complete_inter.rs` |
//! | file-local `UpdateHistory(...)` | `:4649-4670` | `shells_complete_inter.rs` |
//! | file-local `GetSubShapes(...)` | `:5386-5406` | `shape_prep.rs` |
//! | `gce_MakeCone(P1, P2, R1, R2)` (`gce_MakeCone.cxx:234-283`) | - | `missing_walls.rs` |
//!
//! UNPORTED (with the OCCT lines they live in):
//!
//! * the 10-argument constructor (`:613-640`) calls `MakeOffsetShape` itself;
//!   the port keeps `new` + [`BRepOffsetMakeOffset::initialize`] +
//!   [`BRepOffsetMakeOffset::make_offset_shape`].
//! * `analyzeProgress` (`:5410-5441`) and `Message_ProgressScope` /
//!   `PIOperation_*`: the port has no progress machinery, so no `!More()` arm
//!   can fire.
//! * the `OCCT_DEBUG` blocks of the ported methods.
//!
//! Every other method of the class - including the two file-local shell
//! builders `BuildShellsCompleteInter` (`:5064-5227`) and `UpdateHistory`
//! (`:4649-4670`, in `shells_complete_inter.rs`), the file-local
//! `GetSubShapes` (`:5386-5406`, in `shape_prep.rs`), and the two
//! `BRepOffset_MakeOffset` wrappers around the `BRepOffset_BuildOffsetFaces`
//! tool (`_1.cxx:9497-9506`, `:9514-9533`, in `build_offset_faces/mod.rs`) -
//! is ported, together with the whole tool (`_1.cxx:529-9496`, see
//! `build_offset_faces/`).
//!
//! Field notes, all "OCCT leaves it undefined, the port initializes it":
//!
//! * `myInitialShape` / `myShape` / `myOffsetShape` / `myBadShape` are
//!   `TopoDS_Shape` handles that OCCT leaves null via `Nullify()`; the port
//!   carries them as `Option<TopoShape>` ([`TopoShape::is_null`] is always
//!   false here by construction).
//! * `myFaceComp` (`TopoDS_Compound`): `BuildFaceComp` does
//!   `aBB.MakeCompound(myFaceComp)`, which replaces the TShape wholesale; the
//!   port keeps an empty `TopoDS_Compound` before that call, which both
//!   `TopExp_Explorer` sites (`:1829`, `:1877`) enumerate identically to a
//!   null compound.
//! * `myError` is not in any OCCT initializer list (both constructors leave it
//!   unset); the port starts at [`Error::NoError`].
//! * `myFacePlanfaceMap`, `myEdgeIntEdges`, `myResMap` are
//!   `NCollection_DataMap` / `NCollection_Map` keyed by `TopTools_ShapeMapHasher`
//!   (`IsSame`), carried as `HashMap<ShapeMapKey, _>` / `HashSet<ShapeMapKey>`.
//!   `SetFacesWithOffset` (`:822-833`) walks the map and rebinds independent
//!   keys, so the walk order does not change the outcome.
//! * `myFaceOffset` cannot be a bare `HashMap<ShapeMapKey, f64>`: `:2356`
//!   casts `it.Key()` back to a face, so [`FaceOffsetMap`] keeps the keys and
//!   the insertion order like the `NCollection_DataMap` it replaces.

#![allow(dead_code)]

use std::collections::{HashMap, HashSet};

use occt_core::kernel::geomabs::JoinType;

use crate::brep_algo::AsDes;
use crate::brep_algo_image::BRepAlgoImage;
use crate::shape::{Compound, Face, TopoShape};
use crate::topo_tools::{shape_map_key, ShapeMapKey};

use super::inter3d::ShapeIndexedMap;
use super::{Analyse, BRepOffsetMakeLoops, Error, Mode};

mod build_offset_by_arc;
mod build_offset_by_inter;
mod build_offset_faces;
mod conical_faces;
mod input_checks;
mod missing_walls;
mod internal_edges;
mod intersection;
mod intersect_trim;
mod is_planar;
mod loop_rebuild;
mod make_offset_shape;
mod regularity;
mod tolerance_update;
mod face_comp;
mod history;
mod self_inter;
mod make_offset_faces;
mod shape_prep;
mod shells;
mod shells_complete_inter;
mod solids;
mod thick_solid;

pub use build_offset_by_arc::{build_offset_by_arc, to_context};
pub use build_offset_by_inter::build_offset_by_inter;
pub use build_offset_faces::*;
pub use input_checks::check_single_point;
pub use shape_prep::{
    append_to_list, eval_max, is_connected_shell, is_solid, make_list, remove_corks,
};
pub use intersect_trim::{
    get_enlarged_faces, intersect_edges, trim_edge, trim_edges, EdgeIntEdges, EdgesOrigins,
};
pub use face_comp::{build_face_comp, replace_roots, set_faces, set_faces_with_offset};
pub use is_planar::is_planar;
pub use make_offset_faces::{make_offset_faces, update_face_offset};

/// `NCollection_DataMap<TopoDS_Shape, double, TopTools_ShapeMapHasher>`, the
/// type of `myFaceOffset` (`hxx:262`).
///
/// It has to keep the **keys themselves**, not only their
/// [`ShapeMapKey`]: `UpdateFaceOffset` (`BRepOffset_MakeOffset.cxx:2356`)
/// iterates the map and casts `it.Key()` back to a face, which a bare
/// key -> `f64` map cannot answer. `Bind` replaces the item of an existing key
/// and keeps its position (`NCollection_DataMap.hxx:463-472`), `UnBind` removes
/// it, and `Assign` is a plain clone.
#[derive(Debug, Default, Clone)]
pub struct FaceOffsetMap {
    order: Vec<TopoShape>,
    index: HashMap<ShapeMapKey, usize>,
    values: Vec<f64>,
}

impl FaceOffsetMap {
    pub fn new() -> Self {
        Self::default()
    }

    /// `Bind(K, I)`: insert, or overwrite the item in place.
    pub fn bind(&mut self, key: &TopoShape, value: f64) {
        let k = crate::topo_tools::shape_map_key(key);
        if let Some(&i) = self.index.get(&k) {
            self.values[i] = value;
            return;
        }
        self.order.push(key.clone());
        self.values.push(value);
        self.index.insert(k, self.order.len() - 1);
    }

    /// `IsBound(K)`.
    pub fn is_bound(&self, key: &TopoShape) -> bool {
        self.index
            .contains_key(&crate::topo_tools::shape_map_key(key))
    }

    /// `operator()(K)`; `None` where OCCT raises `Standard_NoSuchObject`.
    pub fn get(&self, key: &TopoShape) -> Option<f64> {
        self.index
            .get(&crate::topo_tools::shape_map_key(key))
            .map(|&i| self.values[i])
    }

    /// `UnBind(K)`.
    pub fn un_bind(&mut self, key: &TopoShape) {
        let k = crate::topo_tools::shape_map_key(key);
        let Some(i) = self.index.remove(&k) else {
            return;
        };
        self.order.remove(i);
        self.values.remove(i);
        for idx in self.index.values_mut() {
            if *idx > i {
                *idx -= 1;
            }
        }
    }

    /// `Clear()`.
    pub fn clear(&mut self) {
        self.order.clear();
        self.index.clear();
        self.values.clear();
    }

    /// `NCollection_DataMap::Iterator` / `Extent()`: the key/item pairs in
    /// insertion order.
    pub fn iter(&self) -> impl Iterator<Item = (&TopoShape, f64)> {
        self.order.iter().zip(self.values.iter().copied())
    }

    pub fn extent(&self) -> usize {
        self.order.len()
    }

    /// The plain key -> item view `BRepOffset_Analyse::SetFaceOffsetMap`
    /// (`BRepOffset_Analyse.hxx:140-146`) asks for.
    pub fn to_hash_map(&self) -> HashMap<ShapeMapKey, f64> {
        self.index
            .iter()
            .map(|(k, &i)| (k.clone(), self.values[i]))
            .collect()
    }
}

/// `class BRepOffset_MakeOffset` (`hxx:47-284`).
pub struct BRepOffsetMakeOffset {
    /// `myOffset` (`hxx:249`).
    my_offset: f64,
    /// `myTol` (`hxx:250`).
    my_tol: f64,
    /// `myInitialShape` (`hxx:251`).
    my_initial_shape: Option<TopoShape>,
    /// `myShape` (`hxx:252`): the working copy, `RemoveCorks` shrinks it.
    my_shape: Option<TopoShape>,
    /// `myFaceComp` (`hxx:253`).
    my_face_comp: Compound,
    /// `myMode` (`hxx:255`).
    my_mode: Mode,
    /// `myIsLinearizationAllowed` (`hxx:256`).
    my_is_linearization_allowed: bool,
    /// `myInter` (`hxx:257`).
    my_inter: bool,
    /// `mySelfInter` (`hxx:258`).
    my_self_inter: bool,
    /// `myJoin` (`hxx:259`).
    my_join: JoinType,
    /// `myThickening` (`hxx:260`).
    my_thickening: bool,
    /// `myRemoveIntEdges` (`hxx:261`).
    my_remove_int_edges: bool,
    /// `myFaceOffset` (`hxx:262`): per-face offset override, `SetOffsetOnFace`.
    my_face_offset: FaceOffsetMap,
    /// `myFaces` (`hxx:263`): the faces to offset, caps excluded.
    my_faces: ShapeIndexedMap,
    /// `myOriginalFaces` (`hxx:264`): the closing faces added by `AddFace`.
    my_original_faces: ShapeIndexedMap,
    /// `myAnalyse` (`hxx:265`).
    my_analyse: Analyse,
    /// `myOffsetShape` (`hxx:266`).
    my_offset_shape: Option<TopoShape>,
    /// `myInitOffsetFace` (`hxx:267`).
    my_init_offset_face: BRepAlgoImage,
    /// `myInitOffsetEdge` (`hxx:268`).
    my_init_offset_edge: BRepAlgoImage,
    /// `myImageOffset` (`hxx:269`).
    my_image_offset: BRepAlgoImage,
    /// `myImageVV` (`hxx:270`).
    my_image_vv: BRepAlgoImage,
    /// `myWalls` (`hxx:271`).
    my_walls: Vec<TopoShape>,
    /// `myAsDes` (`hxx:272`).
    my_as_des: AsDes,
    /// `myEdgeIntEdges` (`hxx:273-274`).
    my_edge_int_edges: HashMap<ShapeMapKey, Vec<TopoShape>>,
    /// `myDone` (`hxx:275`).
    my_done: bool,
    /// `myError` (`hxx:276`).
    my_error: Error,
    /// `myMakeLoops` (`hxx:277`).
    my_make_loops: BRepOffsetMakeLoops,
    /// `myIsPerformSewing` (`hxx:278`).
    my_is_perform_sewing: bool,
    /// `myIsPlanar` (`hxx:279`).
    my_is_planar: bool,
    /// `myBadShape` (`hxx:280`).
    my_bad_shape: Option<TopoShape>,
    /// `myFacePlanfaceMap` (`hxx:281`): face -> the plane that replaces it when
    /// `IsPlanar` linearizes it.
    my_face_planface_map: HashMap<ShapeMapKey, TopoShape>,
    /// `myGenerated` (`hxx:282`).
    my_generated: Vec<TopoShape>,
    /// `myResMap` (`hxx:283`).
    my_res_map: HashSet<ShapeMapKey>,
}

impl Default for BRepOffsetMakeOffset {
    fn default() -> Self {
        Self::new()
    }
}

impl BRepOffsetMakeOffset {
    /// `BRepOffset_MakeOffset()` (`cxx:606-609`): only `myAsDes` is built.
    ///
    /// The other counters start at their neutral value (see the module notes);
    /// in particular `myIsLinearizationAllowed` is `false` here whereas the
    /// 10-argument constructor sets it to `true` (`cxx:636`), so a caller that
    /// only runs `new` + `is_planar` must call
    /// [`BRepOffsetMakeOffset::allow_linearization`] first, exactly like OCCT.
    pub fn new() -> Self {
        Self {
            my_offset: 0.0,
            my_tol: 0.0,
            my_initial_shape: None,
            my_shape: None,
            my_face_comp: Compound::new(),
            my_mode: Mode::Skin,
            my_is_linearization_allowed: false,
            my_inter: false,
            my_self_inter: false,
            my_join: JoinType::Arc,
            my_thickening: false,
            my_remove_int_edges: false,
            my_face_offset: FaceOffsetMap::new(),
            my_faces: ShapeIndexedMap::new(),
            my_original_faces: ShapeIndexedMap::new(),
            my_analyse: Analyse::new(),
            my_offset_shape: None,
            my_init_offset_face: BRepAlgoImage::new(),
            my_init_offset_edge: BRepAlgoImage::new(),
            my_image_offset: BRepAlgoImage::new(),
            my_image_vv: BRepAlgoImage::new(),
            my_walls: Vec::new(),
            my_as_des: AsDes::new(),
            my_edge_int_edges: HashMap::new(),
            my_done: false,
            my_error: Error::NoError,
            my_make_loops: BRepOffsetMakeLoops::new(),
            my_is_perform_sewing: false,
            my_is_planar: false,
            my_bad_shape: None,
            my_face_planface_map: HashMap::new(),
            my_generated: Vec::new(),
            my_res_map: HashSet::new(),
        }
    }

    /// `Initialize(S, Offset, Tol, Mode, Inter, SelfInter, Join, Thickening,
    /// RemoveIntEdges)` (`cxx:643-669`).
    ///
    /// The tail order matters: `myDone`, `myIsPerformSewing` and `myIsPlanar`
    /// are reset first and `Clear()` runs last, so the caches it wipes never
    /// carry a stale `myDone`.
    pub fn initialize(
        &mut self,
        s: &TopoShape,
        offset: f64,
        tol: f64,
        mode: Mode,
        inter: bool,
        self_inter: bool,
        join: JoinType,
        thickening: bool,
        remove_int_edges: bool,
    ) {
        self.my_offset = offset;
        self.my_initial_shape = Some(s.clone());
        self.my_shape = Some(s.clone());
        self.my_tol = tol;
        self.my_mode = mode;
        self.my_inter = inter;
        self.my_self_inter = self_inter;
        self.my_join = join;
        self.my_thickening = thickening;
        self.my_remove_int_edges = remove_int_edges;
        self.my_is_linearization_allowed = true;
        self.my_done = false;
        self.my_is_perform_sewing = false;
        self.my_is_planar = false;
        self.clear();
    }

    /// `Clear()` (`cxx:672-688`).
    pub fn clear(&mut self) {
        self.my_offset_shape = None;
        self.my_init_offset_face.clear();
        self.my_init_offset_edge.clear();
        self.my_image_offset.clear();
        self.my_image_vv.clear();
        self.my_faces = ShapeIndexedMap::new();
        self.my_original_faces = ShapeIndexedMap::new();
        self.my_face_offset.clear();
        self.my_edge_int_edges.clear();
        self.my_as_des.clear();
        self.my_done = false;
        self.my_generated.clear();
        self.my_res_map.clear();
    }

    /// `AllowLinearization(theIsAllowed)` (`cxx:691-694`).
    pub fn allow_linearization(&mut self, is_allowed: bool) {
        self.my_is_linearization_allowed = is_allowed;
    }

    /// `AddFace(F)` (`cxx:698-703`): register a closing face.
    pub fn add_face(&mut self, f: &crate::shape::Face) {
        self.my_original_faces.add(&f.0);
    }

    /// `SetOffsetOnFace(F, Off)` (`cxx:706-708`).
    pub fn set_offset_on_face(&mut self, f: &crate::shape::Face, off: f64) {
        self.my_face_offset.bind(&f.0, off);
    }

    /// `IsDone()` (`cxx:1181-1183`).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// `Error()` (`cxx:1186-1191`).
    pub fn error(&self) -> Error {
        self.my_error
    }

    /// `Shape()` (`cxx:1195-1197`).
    pub fn shape(&self) -> Option<&TopoShape> {
        self.my_offset_shape.as_ref()
    }

    /// `InitShape()` (`hxx:98`).
    pub fn init_shape(&self) -> Option<&TopoShape> {
        self.my_initial_shape.as_ref()
    }

    /// `GetJoinType()` (`cxx:3936-3941`).
    pub fn get_join_type(&self) -> JoinType {
        self.my_join
    }

    /// `OffsetFacesFromShapes()` (`cxx:3929-3932`).
    pub fn offset_faces_from_shapes(&self) -> &BRepAlgoImage {
        &self.my_init_offset_face
    }

    /// `OffsetEdgesFromShapes()` (`cxx:3945-3948`).
    pub fn offset_edges_from_shapes(&self) -> &BRepAlgoImage {
        &self.my_init_offset_edge
    }

    /// `ClosingFaces()` (`cxx:3954-3957`).
    pub fn closing_faces(&self) -> &ShapeIndexedMap {
        &self.my_original_faces
    }

    /// `GetBadShape()` (`cxx:4505-4507`).
    pub fn get_bad_shape(&self) -> Option<&TopoShape> {
        self.my_bad_shape.as_ref()
    }

    /// Read-only view of `myAsDes`.
    pub fn as_des(&self) -> &AsDes {
        &self.my_as_des
    }

    /// Mutable view of the owned `myAsDes` handle.
    pub fn as_des_mut(&mut self) -> &mut AsDes {
        &mut self.my_as_des
    }

    /// Read-only view of `myAnalyse`.
    pub fn analyse(&self) -> &Analyse {
        &self.my_analyse
    }

    /// `myIsPlanar` / `myDone` / `myError` accessors used by the rest of the
    /// class once `MakeOffsetShape` lands.
    pub fn is_planar_flag(&self) -> bool {
        self.my_is_planar
    }

    /// `myMode`.
    pub fn mode(&self) -> Mode {
        self.my_mode
    }

    /// `myOffset`.
    pub fn offset(&self) -> f64 {
        self.my_offset
    }

    /// `myTol`.
    pub fn tol(&self) -> f64 {
        self.my_tol
    }

    /// `myFaces` (`hxx:263`): insertion-ordered, `IsSame`-unique.
    pub fn faces(&self) -> &ShapeIndexedMap {
        &self.my_faces
    }

    /// `myFaceOffset` (`hxx:262`).
    pub fn face_offsets(&self) -> &FaceOffsetMap {
        &self.my_face_offset
    }

    /// `TopoDS::Face(it.Key())` on a `myFaceOffset` key.
    ///
    /// The map stores the keys themselves, so this is just the key walk that
    /// `UpdateFaceOffset` (`cxx:2356`) and `SetFacesWithOffset` (`cxx:825`) do;
    /// it stays here so both call sites share the one lookup.
    fn face_of_key(&self, key: &ShapeMapKey) -> Option<Face> {
        self.my_face_offset
            .iter()
            .find(|(s, _)| shape_map_key(s) == *key)
            .and_then(|(s, _)| Face::wrap(s.clone()))
    }

    /// `myFaceComp` (`hxx:253`).
    pub fn face_comp(&self) -> &Compound {
        &self.my_face_comp
    }
}
