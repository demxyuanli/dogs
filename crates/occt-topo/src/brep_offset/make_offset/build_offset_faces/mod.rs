//! `BRepOffset_BuildOffsetFaces` - the file-local class of
//! `BRepOffset_MakeOffset_1.cxx` (`cxx:529-1032`).
//!
//! `BRepOffset_MakeOffset` instantiates it in mode `Complete` with join type
//! `Intersection` to build the splits of the offset faces and, from them, the
//! solids (`BRepOffset_MakeOffset.cxx:1508`, `:3020`, `:3089`). The whole file
//! is that one class plus its anonymous-namespace helpers, which is why the
//! class has no header of its own.
//!
//! Ported so far (T-461 to T-475):
//!
//! | OCCT | Rust |
//! |------|------|
//! | `AddToContainer` x5 (`cxx:78-114`) | [`AddToContainer`], [`add_to_container_of_map`] |
//! | `TakeModified` x2 (`cxx:124-164`) | [`take_modified`] |
//! | `hasData` (`cxx:170-175`) | [`has_data`] |
//! | `AppendToList` (`cxx:180-191`) | [`append_to_list`] |
//! | `ProcessMicroEdge` (`cxx:197-216`) | [`process_micro_edge`] |
//! | `UpdateOrigins` (`cxx:220-261`) | [`update_origins`] |
//! | `UpdateImages` (`cxx:265-305`) | [`update_images`] |
//! | `FindCommonParts` (`cxx:312-346`) | [`find_common_parts`] |
//! | `NbPoints` (`cxx:351-362`) | [`nb_points`] |
//! | `FindShape` (`cxx:367-405`) | [`find_shape`] |
//! | `GetAverageTangent` (`cxx:491-521`) | [`get_average_tangent`] |
//! | `BuildSplitsOfTrimmedFace` (`cxx:411-433`) | `splits::build_splits_of_trimmed_face` |
//! | `BuildSplitsOfFace` (`cxx:438-486`) | `splits::build_splits_of_face` |
//! | `BuildSplitsOfTrimmedFaces` (`cxx:1036-1086`) | `BRepOffsetBuildOffsetFaces::build_splits_of_trimmed_faces` |
//! | `IntersectTrimmedEdges` (`cxx:1207-1292`) | `BRepOffsetBuildOffsetFaces::intersect_trimmed_edges` |
//! | `UpdateIntersectedEdges` (`cxx:1172-1203`) | `BRepOffsetBuildOffsetFaces::update_intersected_edges` |
//! | `GetEdges` (`cxx:1839-1953`) | `BRepOffsetBuildOffsetFaces::get_edges` |
//! | `FillHistory` (`cxx:9413-9496`) | `BRepOffsetBuildOffsetFaces::fill_history` |
//! | `checkConnectionsOfFace` (`cxx:1301-1341`) | [`check_connections_of_face`] |
//! | `CheckIfArtificial` (`cxx:1960-2056`) | `BRepOffsetBuildOffsetFaces::check_if_artificial` |
//! | `CheckInverted` (`cxx:3316-3516`) | `BRepOffsetBuildOffsetFaces::check_inverted` |
//! | `FindInvalidEdges` (`cxx:2057-2528`, per face) | `BRepOffsetBuildOffsetFaces::find_invalid_edges` |
//! | `addAsNeutral` (`cxx:2534-2561`) | [`add_as_neutral`] |
//! | `FindInvalidEdges` (`cxx:2569-2774`, additional pass) | `BRepOffsetBuildOffsetFaces::find_unclassified_invalid_edges` |
//! | `MakeInvertedEdgesInvalid` (`cxx:2778-2853`) | `BRepOffsetBuildOffsetFaces::make_inverted_edges_invalid` |
//! | `FindInvalidFaces` (`cxx:2856-3085`) | `BRepOffsetBuildOffsetFaces::find_invalid_faces` |
//! | `FindFacesInsideHoleWires` (`cxx:3089-3310`) | `BRepOffsetBuildOffsetFaces::find_faces_inside_hole_wires` |
//! | `FilterEdgesImages` (`cxx:5320-5363`) | `BRepOffsetBuildOffsetFaces::filter_edges_images` |
//! | `FilterInvalidFaces` (`cxx:5369-5549`) | `BRepOffsetBuildOffsetFaces::filter_invalid_faces` |
//! | `CheckEdgesCreatedByVertex` (`cxx:5555-5601`) | `BRepOffsetBuildOffsetFaces::check_edges_created_by_vertex` |
//! | `FilterInvalidEdges` (`cxx:5606-5750`) | `BRepOffsetBuildOffsetFaces::filter_invalid_edges` |
//! | `GetVerticesOnEdges` (`cxx:3522-3542`) | `inverted_blocks::get_vertices_on_edges` |
//! | `CheckInvertedBlock` (`cxx:3545-3681`) | `BRepOffsetBuildOffsetFaces::check_inverted_block` |
//! | `RemoveInvalidSplitsByInvertedEdges` (`cxx:3685-3886`) | `BRepOffsetBuildOffsetFaces::remove_invalid_splits_by_inverted_edges` |
//! | `RemoveInvalidSplitsFromValid` (`cxx:3893-4070`) | `BRepOffsetBuildOffsetFaces::remove_invalid_splits_from_valid` |
//! | `buildPairs` (`cxx:4076-4106`) | `inside_faces::build_pairs` |
//! | `buildIntersectionPairs` (`cxx:4110-4287`) | `inside_faces::build_intersection_pairs` |
//! | `RemoveInsideFaces` (`cxx:4295-4636`) | `BRepOffsetBuildOffsetFaces::remove_inside_faces` |
//! | `ShapesConnections` (`cxx:4643-4964`) | `BRepOffsetBuildOffsetFaces::shapes_connections` |
//! | `RemoveHangingParts` (`cxx:4968-5164`) | `BRepOffsetBuildOffsetFaces::remove_hanging_parts` |
//! | `RemoveValidSplits` (`cxx:5168-5223`) | `BRepOffsetBuildOffsetFaces::remove_valid_splits` |
//! | `RemoveInvalidSplits` (`cxx:5228-5315`) | `BRepOffsetBuildOffsetFaces::remove_invalid_splits` |
//! | `FindFacesToRebuild` (`cxx:5758-5897`) | `BRepOffsetBuildOffsetFaces::find_faces_to_rebuild` |
//! | `mapShapes` (`cxx:5904-5917`) | `rebuild_prep::map_shapes` |
//! | `PrepareFacesForIntersection` (`cxx:6571-6674`) | `BRepOffsetBuildOffsetFaces::prepare_faces_for_intersection` |
//! | `FindVerticesToAvoid` (`cxx:6680-6782`) | `BRepOffsetBuildOffsetFaces::find_vertices_to_avoid` |
//! | `FindFacesForIntersection` (`cxx:6788-6969`) | `BRepOffsetBuildOffsetFaces::find_faces_for_intersection` |
//! | `ProcessCommonEdges` (`cxx:6975-7106`) | `BRepOffsetBuildOffsetFaces::process_common_edges` |
//! | `FindOrigins` (`cxx:7113-7152`) | `common_edges::find_origins` |
//! | `UpdateIntersectedFaces` (`cxx:7159-7227`) | `BRepOffsetBuildOffsetFaces::update_intersected_faces` |
//! | `IntersectFaces` driver (`cxx:5924-6565`) | `BRepOffsetBuildOffsetFaces::intersect_faces_driver` |
//! | `IntersectFaces` pair (`cxx:7233-7311`) | `BRepOffsetBuildOffsetFaces::intersect_faces_of_pair` |
//! | `IntersectAndTrimEdges` (`cxx:7319-7604`) | `BRepOffsetBuildOffsetFaces::intersect_and_trim_edges` |
//! | `GetInvalidEdges` (`cxx:7607-7676`) | `BRepOffsetBuildOffsetFaces::get_invalid_edges` |
//! | `UpdateValidEdges` (`cxx:7679-8236`) | `BRepOffsetBuildOffsetFaces::update_valid_edges` |
//! | `TrimNewIntersectionEdges` (`cxx:8239-8401`) | `BRepOffsetBuildOffsetFaces::trim_new_intersection_edges` |
//! | `IntersectEdges` (`cxx:8408-8545`) | `BRepOffsetBuildOffsetFaces::intersect_edges` |
//! | `GetBounds` (`cxx:8551-8587`) | `BRepOffsetBuildOffsetFaces::get_bounds` |
//! | `GetBoundsToUpdate` (`cxx:8591-8667`) | `BRepOffsetBuildOffsetFaces::get_bounds_to_update` |
//! | `GetInvalidEdgesByBounds` (`cxx:8672-8990`) | `BRepOffsetBuildOffsetFaces::get_invalid_edges_by_bounds` |
//! | `FilterSplits` (`cxx:8994-9042`) | `BRepOffsetBuildOffsetFaces::filter_splits` |
//! | `UpdateNewIntersectionEdges` (`cxx:9048-9292`) | `BRepOffsetBuildOffsetFaces::update_new_intersection_edges` |
//! | `FillGaps` (`cxx:9297-9407`) | `BRepOffsetBuildOffsetFaces::fill_gaps` |
//! | `BuildSplitsOfFaces` (`cxx:1342-1833`) | `BRepOffsetBuildOffsetFaces::build_splits_of_faces` |
//! | `BuildSplitsOfExtendedFaces` (`cxx:1088-1169`) | `BRepOffsetBuildOffsetFaces::build_splits_of_extended_faces` |
//! | the `BuildSplitsOfTrimmedFaces` wrapper (`cxx:9497-9506`) | [`build_splits_of_trimmed_faces`] |
//! | the `BuildSplitsOfExtendedFaces` wrapper (`cxx:9514-9533`) | [`build_splits_of_extended_faces`] |
//! | `UpdateImages` over `myOEImages` (`cxx:265-305`) | [`update_images_dmvv`] |
//! | `NCollection_DataMap::UnBind` (`cxx:8289`, `:8371`, `:9038`) | [`Dmvv::unbind`], `HashMap::remove` |
//! | `IntTools_Context::ComputeVF` UV outputs (`IntTools_Context.cxx:545-592`) | `IntToolsContext::compute_vf_uv` |
//! | `IntTools_Context::IsPointInOnFace` (`IntTools_Context.cxx:590-608`) | `IntToolsContext::is_point_in_on_face_2d` |
//! | `BOPAlgo_BOP` kept for `IsDeleted` (`cxx:8210-8228`) | `bop_builder2::builder_bop_algo` |
//! | `MakeConnexityBlocks(S, VERTEX, EDGE, LCB)` (`BOPTools_AlgoTools.cxx:158-183`) | `AlgoTools::make_connexity_blocks_in` |
//! | constructor + setters (`cxx:534-576`) | [`BRepOffsetBuildOffsetFaces`] |
//! | the field list (`cxx:984-1032`) | [`BRepOffsetBuildOffsetFaces`] fields |
//!
//! UNPORTED (with the OCCT lines they live in):
//!
//! * `analyzeProgress` (`BRepOffset_MakeOffset.cxx:5410-5441`): it only fills
//!   the per-operation share table `aSteps` that slices the
//!   `Message_ProgressRange` of `MakeOffset` (`:1286-1326`). The port has no
//!   progress machinery, so no caller and no consumer exists; every other
//!   method of the class is ported.
//!
//! Container notes (the reference's pointers and OCCT collections):
//!
//! * `myFaces`, `myEdgesOrigins`, `myFacesOrigins` and `myETrimEInf` are
//!   pointers to containers owned by `BRepOffset_MakeOffset`. The port owns
//!   copies (a `Vec` / three `HashMap`s) and offers read-back accessors, so a
//!   caller that must observe the updates reads them back through
//!   [`BRepOffsetBuildOffsetFaces::faces_origins`] and
//!   [`BRepOffsetBuildOffsetFaces::edges_origins`]. `hasData` therefore cannot
//!   see a null pointer (see [`has_data`]).
//! * `myImage` (`cxx:1031`) is *not* stored: `BRepAlgoImage` is a value type
//!   here, so a stored copy would swallow the history writes OCCT makes
//!   through the pointer. The methods that write the history take
//!   `&mut BRepAlgoImage` explicitly instead.
//! * `NCollection_IndexedDataMap<Shape, List>` fields (`myOFImages`,
//!   `myInvalidFaces`, `myFacesToRebuild`) are carried by [`Dmvv`], which keeps
//!   the key order OCCT's `IndexedDataMap` iterates. `NCollection_Map` fields
//!   are carried by [`ShapeIndexedMap`] (IsSame-unique; the extra insertion
//!   order does not change membership).
//! * `myArtInvalidFaces` (`DataMap<Shape, IndexedMap<Shape>>`) is carried by
//!   `HashMap<ShapeMapKey, ShapeIndexedMap>`; `myIntersectionPairs`
//!   (`DataMap<Shape, DataMap<Shape, Map>>`) nests one more level of the same
//!   two carriers.

#![allow(dead_code)]

use std::collections::HashMap;

use crate::brep_algo::AsDes;
use crate::brep_algo_image::BRepAlgoImage;
use crate::brep_offset::inter3d::ShapeIndexedMap;
use crate::brep_offset::{Analyse, Dmvv};
use crate::int_tools_full::IntToolsContext;
use crate::shape::TopoShape;
use crate::topo_tools::{shape_map_key, ShapeMapKey};

mod edges;
mod hole_wires;
mod invalid_edges;
mod invalid_edges_by_bounds;
mod invalid_edges_filter;
mod invalid_faces;
mod inverted_blocks;
mod splits;
mod statics;
mod validity;

mod bounds;
mod build_splits;
mod common_edges;
mod connections;
mod fill_gaps;
mod hanging_splits;
mod inside_faces;
mod intersect_and_trim;
mod intersect_faces;
mod rebuild_prep;
mod rebuild_driver;
mod trim_edges;
mod update_valid_edges;
pub use invalid_edges::add_as_neutral;
pub use statics::*;

/// `class BRepOffset_BuildOffsetFaces` (`cxx:529-1032`).
pub struct BRepOffsetBuildOffsetFaces {
    // Input data (`cxx:984-992`).
    /// `myFaces` (`cxx:985`): the faces to split.
    my_faces: Vec<TopoShape>,
    /// `myAsDes` (`cxx:986`): ascendants/descendants of the edges and faces.
    my_as_des: AsDes,
    /// `myAnalyzer` (`cxx:987`): analysis of the input shape.
    ///
    /// The reference keeps a pointer initialised to null (`cxx:535`) and only
    /// `SetAnalysis` (`cxx:552`) binds it, so the two entry points differ:
    /// `BuildSplitsOfExtendedFaces` sets it (`cxx:9527`) while
    /// `BuildSplitsOfTrimmedFaces` does not (`cxx:9501-9505`). `None` is that
    /// null pointer; `FindInvalidEdges` tests it at `cxx:2150`, `:2271`,
    /// `:2361`, `:2591`.
    my_analyzer: Option<Analyse>,
    /// `myEdgesOrigins` (`cxx:990`): offset edge -> original edges.
    my_edges_origins: HashMap<ShapeMapKey, Vec<TopoShape>>,
    /// `myFacesOrigins` (`cxx:991`): offset face -> original face.
    my_faces_origins: HashMap<ShapeMapKey, TopoShape>,
    /// `myETrimEInf` (`cxx:992`): trimmed edge -> infinite edge.
    my_e_trim_e_inf: HashMap<ShapeMapKey, TopoShape>,

    // Intermediate data (`cxx:995-1027`).
    /// `myOEImages` (`cxx:995`): offset edge -> its splits.
    ///
    /// The reference type is `NCollection_DataMap<Shape, List>` (unordered).
    /// It is carried by [`Dmvv`] because `FillHistory` (`cxx:9464-9490`)
    /// iterates the *keys* and hands each one to the history; a
    /// `HashMap<ShapeMapKey, _>` would keep only the identity key and lose the
    /// shape. The extra insertion order `Dmvv` keeps is inert here - the
    /// reference's comment at `cxx:9462` states the order is irrelevant.
    my_oe_images: Dmvv,
    /// `myOEOrigins` (`cxx:996`): split of an offset edge -> its origins.
    my_oe_origins: HashMap<ShapeMapKey, Vec<TopoShape>>,
    /// `myOFImages` (`cxx:997`): offset edge -> its images, key order kept.
    my_of_images: Dmvv,
    /// `myInvalidEdges` (`cxx:999`): edges whose orientation changed.
    my_invalid_edges: ShapeIndexedMap,
    /// `myValidEdges` (`cxx:1000`): edges whose orientation was kept.
    my_valid_edges: ShapeIndexedMap,
    /// `myInvertedEdges` (`cxx:1001`): edges with swapped vertices.
    my_inverted_edges: ShapeIndexedMap,
    /// `myEdgesToAvoid` (`cxx:1002`): splits avoided while splitting faces.
    my_edges_to_avoid: ShapeIndexedMap,
    /// `myLastInvEdges` (`cxx:1003`).
    my_last_inv_edges: ShapeIndexedMap,
    /// `myModifiedEdges` (`cxx:1004`): edges used to build the splits.
    my_modified_edges: ShapeIndexedMap,
    /// `myInsideEdges` (`cxx:1005`): edges inside the solids.
    my_inside_edges: ShapeIndexedMap,
    /// `myInvalidFaces` (`cxx:1007`): splits made only of invalid edges.
    my_invalid_faces: Dmvv,
    /// `myArtInvalidFaces` (`cxx:1008`): valid faces marked invalid.
    my_art_invalid_faces: HashMap<ShapeMapKey, ShapeIndexedMap>,
    /// The faces `myArtInvalidFaces` is keyed on, in bind order.
    ///
    /// `NCollection_DataMap<Shape, IndexedMap>` iterates its keys
    /// (`cxx:6002-6011`), which the identity key alone cannot reproduce, so the
    /// port keeps the shapes alongside. `BuildSplitsOfFaces` binds through
    /// [`BRepOffsetBuildOffsetFaces::bind_art_invalid_face`] (`cxx:1516`) and
    /// `FilterInvalidFaces` unbinds through
    /// [`BRepOffsetBuildOffsetFaces::unbind_art_invalid_face`] (`cxx:5397`).
    my_art_invalid_face_order: Vec<TopoShape>,
    /// `myAlreadyInvFaces` (`cxx:1011`): rebuild count per face.
    my_already_inv_faces: HashMap<ShapeMapKey, i32>,
    /// `myFNewHoles` (`cxx:1012`): hole images of the original face.
    my_f_new_holes: HashMap<ShapeMapKey, Vec<TopoShape>>,
    /// `mySSInterfs` (`cxx:1014`): intersection pairs collected on rebuild.
    my_ss_interfs: HashMap<ShapeMapKey, Vec<TopoShape>>,
    /// `mySSInterfsArt` (`cxx:1015`).
    my_ss_interfs_art: HashMap<ShapeMapKey, Vec<TopoShape>>,
    /// `myIntersectionPairs` (`cxx:1016-1019`).
    my_intersection_pairs: HashMap<ShapeMapKey, HashMap<ShapeMapKey, ShapeIndexedMap>>,
    /// `myFacesToRebuild` (`cxx:1022`): invalid and close-to-invalid faces.
    my_faces_to_rebuild: Dmvv,
    /// `myFSelfRebAvoid` (`cxx:1023`): faces avoided while rebuilding.
    my_f_self_reb_avoid: ShapeIndexedMap,
    /// `mySolids` (`cxx:1025`): solids built from the face splits.
    my_solids: Option<TopoShape>,

    // Auxiliary tools (`cxx:1028`).
    /// `myContext` (`cxx:1028`).
    my_context: IntToolsContext,
}

impl Default for BRepOffsetBuildOffsetFaces {
    fn default() -> Self {
        Self::new()
    }
}

impl BRepOffsetBuildOffsetFaces {
    /// `BRepOffset_BuildOffsetFaces(BRepAlgo_Image&)` (`cxx:534-541`).
    ///
    /// `myImage` is not stored (see the module note); `myContext` is built the
    /// way `cxx:541` does (`new IntTools_Context()`), and every other field
    /// starts empty, which is the reference's initializer list plus the
    /// default-constructed members it leaves out.
    pub fn new() -> Self {
        Self {
            my_faces: Vec::new(),
            my_as_des: AsDes::new(),
            my_analyzer: None,
            my_edges_origins: HashMap::new(),
            my_faces_origins: HashMap::new(),
            my_e_trim_e_inf: HashMap::new(),
            my_oe_images: Dmvv::new(),
            my_oe_origins: HashMap::new(),
            my_of_images: Dmvv::new(),
            my_invalid_edges: ShapeIndexedMap::new(),
            my_valid_edges: ShapeIndexedMap::new(),
            my_inverted_edges: ShapeIndexedMap::new(),
            my_edges_to_avoid: ShapeIndexedMap::new(),
            my_last_inv_edges: ShapeIndexedMap::new(),
            my_modified_edges: ShapeIndexedMap::new(),
            my_inside_edges: ShapeIndexedMap::new(),
            my_invalid_faces: Dmvv::new(),
            my_art_invalid_faces: HashMap::new(),
            my_art_invalid_face_order: Vec::new(),
            my_already_inv_faces: HashMap::new(),
            my_f_new_holes: HashMap::new(),
            my_ss_interfs: HashMap::new(),
            my_ss_interfs_art: HashMap::new(),
            my_intersection_pairs: HashMap::new(),
            my_faces_to_rebuild: Dmvv::new(),
            my_f_self_reb_avoid: ShapeIndexedMap::new(),
            my_solids: None,
            my_context: IntToolsContext::new(),
        }
    }

    /// `SetFaces` (`cxx:546`).
    pub fn set_faces(&mut self, the_faces: Vec<TopoShape>) {
        self.my_faces = the_faces;
    }

    /// `SetAsDesInfo` (`cxx:549`).
    pub fn set_as_des_info(&mut self, the_as_des: AsDes) {
        self.my_as_des = the_as_des;
    }

    /// `SetAnalysis` (`cxx:552`): binds `myAnalyzer` (see the field note).
    pub fn set_analysis(&mut self, the_analyse: Analyse) {
        self.my_analyzer = Some(the_analyse);
    }

    /// `myAnalyzer` (`cxx:987`), `None` when the reference's pointer is null.
    pub fn analysis(&self) -> Option<&Analyse> {
        self.my_analyzer.as_ref()
    }

    /// `SetEdgesOrigins` (`cxx:555-560`).
    pub fn set_edges_origins(&mut self, the_edges_origins: HashMap<ShapeMapKey, Vec<TopoShape>>) {
        self.my_edges_origins = the_edges_origins;
    }

    /// `SetFacesOrigins` (`cxx:562-566`).
    pub fn set_faces_origins(&mut self, the_faces_origins: HashMap<ShapeMapKey, TopoShape>) {
        self.my_faces_origins = the_faces_origins;
    }

    /// `SetInfEdges` (`cxx:570-575`).
    pub fn set_inf_edges(&mut self, the_e_trim_e_inf: HashMap<ShapeMapKey, TopoShape>) {
        self.my_e_trim_e_inf = the_e_trim_e_inf;
    }

    /// `myFaces` (`cxx:985`), read back for `hasData` at `cxx:1039`.
    pub fn faces(&self) -> &[TopoShape] {
        &self.my_faces
    }

    /// `myEdgesOrigins` (`cxx:990`), read back by the caller that must observe
    /// the updates the reference makes through the pointer.
    pub fn edges_origins(&self) -> &HashMap<ShapeMapKey, Vec<TopoShape>> {
        &self.my_edges_origins
    }

    /// `myFacesOrigins` (`cxx:991`), read back like
    /// [`BRepOffsetBuildOffsetFaces::edges_origins`].
    pub fn faces_origins(&self) -> &HashMap<ShapeMapKey, TopoShape> {
        &self.my_faces_origins
    }

    /// `myETrimEInf` (`cxx:992`).
    pub fn inf_edges(&self) -> &HashMap<ShapeMapKey, TopoShape> {
        &self.my_e_trim_e_inf
    }

    /// `mySolids` (`cxx:1025`).
    pub fn solids(&self) -> Option<&TopoShape> {
        self.my_solids.as_ref()
    }

    /// `myArtInvalidFaces.Bind(theF, theMIE)` (`cxx:1516`).
    pub fn bind_art_invalid_face(&mut self, the_f: &TopoShape, the_mie: ShapeIndexedMap) {
        let a_key = shape_map_key(the_f);
        if !self.my_art_invalid_faces.contains_key(&a_key) {
            self.my_art_invalid_face_order.push(the_f.clone());
        }
        self.my_art_invalid_faces.insert(a_key, the_mie);
    }

    /// `myArtInvalidFaces.UnBind(theF)` (`cxx:5397`).
    pub fn unbind_art_invalid_face(&mut self, the_f: &TopoShape) {
        let a_key = shape_map_key(the_f);
        if self.my_art_invalid_faces.remove(&a_key).is_some() {
            self.my_art_invalid_face_order
                .retain(|a| shape_map_key(a) != a_key);
        }
    }

    /// `myArtInvalidFaces.Seek(theF)` (`cxx:5674`, `:6242`).
    pub fn art_invalid_face(&self, the_f: &TopoShape) -> Option<&ShapeIndexedMap> {
        self.my_art_invalid_faces.get(&shape_map_key(the_f))
    }
}

/// `void BRepOffset_MakeOffset::BuildSplitsOfTrimmedFaces(const List& theLF,
/// const handle<BRepAlgo_AsDes>& theAsDes, BRepAlgo_Image& theImage, const
/// Message_ProgressRange&)` (`_1.cxx:9497-9506`).
///
/// The reference's wrapper builds the tool, hands it the faces and the `AsDes`
/// and runs the trimmed-face splits; `myImage` is the caller's `theImage` here
/// (see the module note), so the tool receives it explicitly.
pub fn build_splits_of_trimmed_faces(
    the_lf: &[TopoShape],
    the_as_des: &AsDes,
    the_image: &mut BRepAlgoImage,
) {
    let mut a_bf_tool = BRepOffsetBuildOffsetFaces::new();
    a_bf_tool.set_faces(the_lf.to_vec());
    a_bf_tool.set_as_des_info(the_as_des.clone());
    a_bf_tool.build_splits_of_trimmed_faces(the_image);
}

/// `void BRepOffset_MakeOffset::BuildSplitsOfExtendedFaces(const List& theLF,
/// const BRepOffset_Analyse& theAnalyse, const handle<BRepAlgo_AsDes>&
/// theAsDes, DataMap<Shape, List>& theEdgesOrigins, DataMap<Shape, Shape>&
/// theFacesOrigins, DataMap<Shape, Shape>& theETrimEInf, BRepAlgo_Image&
/// theImage, const Message_ProgressRange&)` (`_1.cxx:9514-9533`).
///
/// The reference binds `myEdgesOrigins` / `myFacesOrigins` / `myETrimEInf` as
/// *pointers*, so the tool writes straight into the caller's maps. The port
/// owns copies, so the two maps the tool updates are written back after the
/// run - `myETrimEInf` is only read and needs no write-back.
#[allow(clippy::too_many_arguments)]
pub fn build_splits_of_extended_faces(
    the_lf: &[TopoShape],
    the_analyse: &Analyse,
    the_as_des: &AsDes,
    the_edges_origins: &mut HashMap<ShapeMapKey, Vec<TopoShape>>,
    the_faces_origins: &mut HashMap<ShapeMapKey, TopoShape>,
    the_e_trim_e_inf: &HashMap<ShapeMapKey, TopoShape>,
    the_image: &mut BRepAlgoImage,
) {
    let mut a_bf_tool = BRepOffsetBuildOffsetFaces::new();
    a_bf_tool.set_faces(the_lf.to_vec());
    a_bf_tool.set_as_des_info(the_as_des.clone());
    a_bf_tool.set_analysis(the_analyse.clone());
    a_bf_tool.set_edges_origins(the_edges_origins.clone());
    a_bf_tool.set_faces_origins(the_faces_origins.clone());
    a_bf_tool.set_inf_edges(the_e_trim_e_inf.clone());

    a_bf_tool.build_splits_of_extended_faces(the_image);

    // `myEdgesOrigins` / `myFacesOrigins` are the caller's maps in the
    // reference; write the tool's updates back.
    *the_edges_origins = a_bf_tool.edges_origins().clone();
    *the_faces_origins = a_bf_tool.faces_origins().clone();
}
