//! `TopOpeBRepBuild_Builder` — state, split/merged maps and the bookkeeping
//! helpers.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_Builder.hxx:78-869`,
//! `.cxx:77-924` (OCCT 8.0.0), `TopOpeBRepBuild_Merge.cxx:174-758`,
//! `TopOpeBRepBuild_Section.cxx:248-254`.
//!
//! Ported in this module:
//! - the state and every accessor of the class (`hxx:801-869`),
//! - `Clear` (`.cxx:182-246`), `InitSection` (`Section.cxx:248-254`),
//! - `NewVertex` / `NewEdges` / `NewFaces` and their `Change*` forms
//!   (`.cxx:249-306`),
//! - the split / merged map protocol (`.cxx:308-494`, `Merge.cxx:626-731`),
//! - `ShapePosition` / `KeepShape` (`.cxx:499-556`),
//! - `TopType` / `Reverse` / `Orient` / `Contains` (`.cxx:560-645`, `:836-849`),
//! - `FindSameDomain` / `FindSameDomainSameOrientation` (`.cxx:650-762`),
//! - `MapShapes` / `ClearMaps` / `FindSameRank` / `ShapeRank` / `IsShapeOf`
//!   (`.cxx:766-830`),
//! - `Opec12` / `Opec21` / `Opecom` / `Opefus` (`.cxx:853-882`),
//! - `Classify` / `ChangeClassify` (`Merge.cxx:341-354`),
//! - `MergeEdges` / `MergeFaces` (`Merge.cxx:733-758`), both empty in OCCT,
//! - the split drivers and the fill family: see `splits.rs` and `fill.rs`,
//! - `MergeShapes` / `MergeSolids` / `MergeSolid` (`Merge.cxx:174-336`,
//!   `:356-361`, `:368-371`) and `IsKPart` (`KPart.cxx:169-172`): see
//!   `merge.rs`,
//! - `Perform` both overloads (`Builder.cxx:116-145`) and `BuildVertices`
//!   (`BuildVertices.cxx:25-34`): see `perform.rs`,
//! - `Section` / `SectionCurves` / `SectionEdges` /
//!   `FillSecEdgeAncestorMap` (`Section.cxx:1234-1392`): see `section.rs`,
//! - `SplitSectionEdge` (`Section.cxx:928-1230`) and `SplitSectionEdges`
//!   (`Section.cxx:304-915`): see `section_edge_split.rs` and `section_edge.rs`,
//! - `BuildFaces` / `BuildFaces(HDS)` (`BuildFaces.cxx:40-107`): see
//!   `build_faces.rs`,
//! - `BuildEdges` / `BuildEdges(HDS)` (`BuildEdges.cxx:38-218`): see
//!   `build_edges.rs`.
//!
//! UNPORTED, with their OCCT location:
//! - Every `KPart*` / `KP*` member (`TopOpeBRepBuild_KPart.cxx:1-2670`,
//!   including `FindIsKPart` `:111-166` and `MergeKPart` `:176-220`).
//! - Every `G*` member (`TopOpeBRepBuild_GTool.cxx` and the `G*` files),
//!   `End` (`End.cxx:45-125`), `Section*` (`Section.cxx:1-247`),
//!   `SplitSectionEdge` (`Section.cxx:928-1248`).

use std::cell::RefCell;
use std::collections::{BTreeMap, HashMap};
use std::rc::Rc;

mod build_edges;
mod build_faces;
mod builder1;
mod curve_topology_wes;
mod end;
mod end_fuse_edges;
mod fill;
mod fregu;
mod grid;
mod grid_ee;
mod gwes_make_faces;
mod k_part;
mod make_edges;
mod make_faces;
mod make_shapes;
mod make_solids;
mod merge;
mod on;
mod on_patches;
mod perform;
mod section;
mod section_edge;
mod section_edge_split;
mod splits;
mod sregu;
mod ffwesk;
mod grid_ff;
mod grid_ff_edges;
mod keep;
mod ffsfs;
mod grid_ss;
mod grid_ss_faces;
mod grid_ss_globals;
mod grid_ss_on;
mod grid_ss_split;
mod split_eviso_on_periodic;
mod merge_states;

pub(crate) use merge_states::map_shapes_into;

pub use builder1::{global_use_new_builder, set_global_use_new_builder, Builder1};
use builder1::Builder1State;

use crate::abs::{Orientation, ShapeType, State};
use crate::shape::TopoShape;
use crate::topo_tools::{explorer_shapes, is_same, shape_map_key, ShapeMapKey};
use crate::topopebrep::ds::build_tool::BuildTool;
use crate::topopebrep::ds::data_structure::{IndexedShapeData, IndexedShapeSet};
use crate::topopebrep::ds::h_data_structure::HDataStructure;
use crate::topopebrep::ds::tool::DShapeOn1StateMap;
use crate::topopebrep::tool::shape_classifier::ShapeClassifier;

/// `class TopOpeBRepBuild_Builder` (`hxx:78`).
pub struct Builder {
    /// `myState1` (`hxx:801`).
    my_state1: State,
    /// `myState2` (`hxx:802`).
    my_state2: State,
    /// `myShape1` (`hxx:803`).
    my_shape1: Option<TopoShape>,
    /// `myShape2` (`hxx:804`).
    my_shape2: Option<TopoShape>,
    /// `myDataStructure` (`hxx:805`).
    my_data_structure: Option<HDataStructure>,
    /// `myBuildTool` (`hxx:806`).
    my_build_tool: BuildTool,
    /// `myNewVertices` (`hxx:807`), a 1-based `HArray1<TopoDS_Shape>`.
    my_new_vertices: Vec<Option<TopoShape>>,
    /// `myNewEdges` (`hxx:808`).
    my_new_edges: BTreeMap<i32, Vec<TopoShape>>,
    /// `myNewFaces` (`hxx:809`), a 1-based `HArray1<List<Shape>>`.
    my_new_faces: Vec<Vec<TopoShape>>,
    /// `mySplitIN` (`hxx:811-812`).
    my_split_in: DShapeOn1StateMap,
    /// `mySplitON` (`hxx:813-814`).
    my_split_on: DShapeOn1StateMap,
    /// `mySplitOUT` (`hxx:815-816`).
    my_split_out: DShapeOn1StateMap,
    /// `myMergedIN` (`hxx:817-818`).
    my_merged_in: DShapeOn1StateMap,
    /// `myMergedON` (`hxx:819-820`).
    my_merged_on: DShapeOn1StateMap,
    /// `myMergedOUT` (`hxx:821-822`).
    my_merged_out: DShapeOn1StateMap,
    /// `myEmptyShapeList` (`hxx:823`).
    my_empty_shape_list: Vec<TopoShape>,
    /// `myListOfSolid` (`hxx:824`).
    my_list_of_solid: Vec<TopoShape>,
    /// `myListOfFace` (`hxx:825`).
    my_list_of_face: Vec<TopoShape>,
    /// `myListOfEdge` (`hxx:826`).
    my_list_of_edge: Vec<TopoShape>,
    /// `myFSplits` (`hxx:827-828`).
    my_f_splits: HashMap<ShapeMapKey, Vec<TopoShape>>,
    /// `myESplits` (`hxx:829-830`).
    my_e_splits: HashMap<ShapeMapKey, Vec<TopoShape>>,
    /// `mySectionDone` (`hxx:831`).
    my_section_done: bool,
    /// `mySplitSectionEdgesDone` (`hxx:832`).
    my_split_section_edges_done: bool,
    /// `mySection` (`hxx:833`).
    my_section: Vec<TopoShape>,
    /// `mySolidReference` (`hxx:834`).
    my_solid_reference: Option<TopoShape>,
    /// `mySolidToFill` (`hxx:835`).
    my_solid_to_fill: Option<TopoShape>,
    /// `myFaceAvoid` (`hxx:836`).
    my_face_avoid: Vec<TopoShape>,
    /// `myFaceReference` (`hxx:837`).
    my_face_reference: Option<TopoShape>,
    /// `myFaceToFill` (`hxx:838`).
    my_face_to_fill: Option<TopoShape>,
    /// `myEdgeAvoid` (`hxx:839`).
    my_edge_avoid: Vec<TopoShape>,
    /// `myEdgeReference` (`hxx:840`).
    my_edge_reference: Option<TopoShape>,
    /// `myEdgeToFill` (`hxx:841`).
    my_edge_to_fill: Option<TopoShape>,
    /// `myVertexAvoid` (`hxx:842`).
    my_vertex_avoid: Vec<TopoShape>,
    /// `myMAP1` (`hxx:843`).
    my_map1: IndexedShapeSet,
    /// `myMAP2` (`hxx:844`).
    my_map2: IndexedShapeSet,
    /// `myIsKPart` (`hxx:845`).
    my_is_k_part: i32,
    /// `myKPMAPf1f2` (`hxx:846-847`).
    ///
    /// `NCollection_DataMap<TopoDS_Shape, List<TopoDS_Shape>>`; the port keeps
    /// a `Vec` of `(key, values)` so the iteration order stays the insertion
    /// order `MergeKPartiskole` relies on when picking `pfOU`
    /// (`kpkole.cxx:135-152`).
    my_kp_map_f1f2: Vec<(TopoShape, Vec<TopoShape>)>,
    /// `mySHASETindex` (`hxx:848`).
    my_sha_set_index: i32,
    /// `myClassifyDef` (`hxx:849`).
    my_classify_def: bool,
    /// `myClassifyVal` (`hxx:850`).
    my_classify_val: bool,
    /// `myShapeClassifier` (`hxx:851`).
    ///
    /// OCCT hands `&myShapeClassifier` to `TopOpeBRepDS_Filter`
    /// (`Builder.cxx:131`), which keeps the raw pointer past the call, so the
    /// classification state (`TopOpeBRepTool_ShapeClassifier`'s cached solid /
    /// shell) is shared rather than copied; the port shares one cell for the
    /// same reason.
    my_shape_classifier: Rc<RefCell<ShapeClassifier>>,
    /// `myMemoSplit` (`hxx:852`).
    my_memo_split: IndexedShapeSet,
    /// `myEmptyAS` (`hxx:853`).
    my_empty_as: String,
    /// `myProcessON` (`hxx:854`).
    my_process_on: bool,
    /// `myONFacesMap` (`hxx:855`).
    my_on_faces_map: IndexedShapeData<TopoShape>,
    /// `myONElemMap` (`hxx:856`).
    my_on_elem_map: IndexedShapeSet,
    /// The `TopOpeBRepBuild_Builder1` private block (`Builder1.hxx:186-203`).
    ///
    /// OCCT stores these in the subclass. The three base members it overrides
    /// (`GFillSolidSFS`, `GFillShellSFS`, `GWESMakeFaces`) are reached from
    /// base bodies that only have `this`, so the port keeps the block in the
    /// base object and lets those members branch on [`Builder::is_builder1`]
    /// exactly where the vtable would. Empty for a plain `Builder`.
    pub(crate) b1: Builder1State,
    /// The dynamic type is `TopOpeBRepBuild_Builder1` (`Builder1.hxx:45`):
    /// set by `Builder1::new` only.
    pub(crate) is_builder1: bool,
}

impl Builder {
    /// `TopOpeBRepBuild_Builder(BT)` (`.cxx:77-85`).
    pub fn new(build_tool: BuildTool) -> Self {
        let mut b = Self {
            my_state1: State::Unknown,
            my_state2: State::Unknown,
            my_shape1: None,
            my_shape2: None,
            my_data_structure: None,
            my_build_tool: build_tool,
            my_new_vertices: Vec::new(),
            my_new_edges: BTreeMap::new(),
            my_new_faces: Vec::new(),
            my_split_in: DShapeOn1StateMap::new(),
            my_split_on: DShapeOn1StateMap::new(),
            my_split_out: DShapeOn1StateMap::new(),
            my_merged_in: DShapeOn1StateMap::new(),
            my_merged_on: DShapeOn1StateMap::new(),
            my_merged_out: DShapeOn1StateMap::new(),
            my_empty_shape_list: Vec::new(),
            my_list_of_solid: Vec::new(),
            my_list_of_face: Vec::new(),
            my_list_of_edge: Vec::new(),
            my_f_splits: HashMap::new(),
            my_e_splits: HashMap::new(),
            my_section_done: false,
            my_split_section_edges_done: false,
            my_section: Vec::new(),
            my_solid_reference: None,
            my_solid_to_fill: None,
            my_face_avoid: Vec::new(),
            my_face_reference: None,
            my_face_to_fill: None,
            my_edge_avoid: Vec::new(),
            my_edge_reference: None,
            my_edge_to_fill: None,
            my_vertex_avoid: Vec::new(),
            my_map1: IndexedShapeSet::new(),
            my_map2: IndexedShapeSet::new(),
            my_is_k_part: 0,
            my_kp_map_f1f2: Vec::new(),
            my_sha_set_index: 0,
            my_classify_def: false,
            my_classify_val: true,
            my_shape_classifier: Rc::new(RefCell::new(ShapeClassifier::new())),
            my_memo_split: IndexedShapeSet::new(),
            my_empty_as: String::new(),
            my_process_on: false,
            my_on_faces_map: IndexedShapeData::new(),
            my_on_elem_map: IndexedShapeSet::new(),
            b1: Builder1State::new(),
            is_builder1: false,
        };
        b.init_section();
        b
    }

    /// `ChangeBuildTool()` (`.cxx:96-99`).
    pub fn change_build_tool(&mut self) -> &mut BuildTool {
        &mut self.my_build_tool
    }

    /// `BuildTool()` (`.cxx:102-105`).
    pub fn build_tool(&self) -> &BuildTool {
        &self.my_build_tool
    }

    /// `DataStructure()` (`.cxx:109-112`).
    pub fn data_structure(&self) -> Option<&HDataStructure> {
        self.my_data_structure.as_ref()
    }

    /// `DataStructure()`, mutable view.
    pub fn change_data_structure(&mut self) -> Option<&mut HDataStructure> {
        self.my_data_structure.as_mut()
    }

    /// Store the DS (the first statement of `Perform`, `.cxx:122`).
    pub fn set_data_structure(&mut self, hds: HDataStructure) {
        self.my_data_structure = Some(hds);
    }

    /// `Clear()` (`.cxx:182-246`): section edges keep their split list, every
    /// other edge is cleared, and the three merged maps are emptied.
    pub fn clear(&mut self) {
        let Some(hds) = self.my_data_structure.as_ref() else {
            self.my_merged_out.clear();
            self.my_merged_in.clear();
            self.my_merged_on.clear();
            return;
        };

        let bds = hds.ds();
        for map in [
            &mut self.my_split_out,
            &mut self.my_split_in,
            &mut self.my_split_on,
        ] {
            for (key, losos) in map.iter_mut() {
                // The map key is the `TopTools_ShapeMapHasher` value, so the
                // key shape is resolved through the DS to run the
                // `ShapeType() == TopAbs_EDGE` and `IsSectionEdge` tests of
                // `cxx:200-205` / `:216-221` / `:231-236`. A key the DS does
                // not hold is not a section edge either, so it is cleared.
                let mut is_edge = false;
                let mut isse = false;
                let nb = bds.nb_shapes();
                for i in 1..=nb {
                    let Some(s) = bds.shape(i, false) else {
                        continue;
                    };
                    if shape_map_key(&s) == *key {
                        is_edge = s.shape_type() == ShapeType::Edge;
                        isse = is_edge && bds.is_section_edge(&s, false);
                        break;
                    }
                }
                if is_edge && !isse {
                    losos.clear();
                }
            }
        }

        self.my_merged_out.clear();
        self.my_merged_in.clear();
        self.my_merged_on.clear();
    }

    /// The shape a `mySplit*` key stands for.
    ///
    /// `TopOpeBRepDS_ListOfShapeOn1State`'s map is keyed by
    /// `TopTools_ShapeMapHasher`, which does not keep the shape itself, while
    /// `Clear` (`cxx:207-213`) and `TopOpeBRepBuild_HBuilder::MakeEdgeAncestorMap`
    /// (`HBuilder.cxx:319-327`) need it back. Both resolve it through the DS;
    /// `clear` does the same lookup inline because it iterates the map it is
    /// mutating. `None` means the DS does not hold that key, which no lookup
    /// should have produced.
    pub fn shape_of_key(&self, key: &ShapeMapKey) -> Option<TopoShape> {
        let hds = self.my_data_structure.as_ref()?;
        let bds = hds.ds();
        let nb = bds.nb_shapes();
        for i in 1..=nb {
            let Some(s) = bds.shape(i, false) else {
                continue;
            };
            if shape_map_key(&s) == *key {
                return Some(s);
            }
        }
        None
    }

    /// `InitSection()` (`Section.cxx:248-254`).
    pub fn init_section(&mut self) {
        self.my_section_done = false;
        self.my_section.clear();
        self.my_split_section_edges_done = false;
        self.my_split_on.clear();
    }

    /// `NewFaces(I)` (`.cxx:249-253`).
    pub fn new_faces(&self, i: i32) -> Option<&[TopoShape]> {
        self.my_new_faces.get((i - 1) as usize).map(|v| v.as_slice())
    }

    /// `ChangeNewFaces(I)` (`.cxx:257-262`).
    pub fn change_new_faces(&mut self, i: i32) -> Option<&mut Vec<TopoShape>> {
        self.my_new_faces.get_mut((i - 1) as usize)
    }

    /// `NewEdges(I)` (`.cxx:265-275`): the empty list when `I` is unbound.
    pub fn new_edges(&self, i: i32) -> &[TopoShape] {
        match self.my_new_edges.get(&i) {
            Some(l) => l.as_slice(),
            None => &self.my_empty_shape_list,
        }
    }

    /// `ChangeNewEdges(I)` (`.cxx:279-287`): binds an empty list first.
    pub fn change_new_edges(&mut self, i: i32) -> &mut Vec<TopoShape> {
        self.my_new_edges.entry(i).or_default()
    }

    /// `NewVertex(I)` (`.cxx:292-296`).
    pub fn new_vertex(&self, i: i32) -> Option<&TopoShape> {
        self.my_new_vertices.get((i - 1) as usize).and_then(|v| v.as_ref())
    }

    /// `ChangeNewVertex(I)` (`.cxx:300-304`).
    pub fn change_new_vertex(&mut self, i: i32) -> Option<&mut Option<TopoShape>> {
        self.my_new_vertices.get_mut((i - 1) as usize)
    }

    /// Size of the `myNewVertices` array (`HArray1::Length`).
    pub fn set_new_vertices_len(&mut self, n: usize) {
        self.my_new_vertices = vec![None; n];
    }

    /// Resize and clear `myNewFaces` (`new HArray1<List<Shape>>(1, N)`).
    pub fn set_new_faces_len(&mut self, n: usize) {
        self.my_new_faces = vec![Vec::new(); n];
    }

    /// `ToSplit(S, ToBuild)` (`.cxx:308-326`).
    pub fn to_split(&self, s: &TopoShape, to_build: State) -> bool {
        let issplit = self.is_split(s, to_build);
        let Some(hds) = self.my_data_structure.as_ref() else {
            return false;
        };
        let hasgeom = hds.has_geometry(s);
        let hassame = hds.has_same_domain(s, true);
        (!issplit) && (hasgeom || hassame)
    }

    /// `MarkSplit(S, ToBuild, Bval)` (`.cxx:330-374`).
    pub fn mark_split(&mut self, s: &TopoShape, to_build: State, bval: bool) {
        let map = match to_build {
            State::Out => &mut self.my_split_out,
            State::In => &mut self.my_split_in,
            State::On => &mut self.my_split_on,
            _ => return,
        };
        map.entry(shape_map_key(s)).or_default().split(bval);
    }

    /// `IsSplit(S, ToBuild)` (`.cxx:378-410`).
    pub fn is_split(&self, s: &TopoShape, to_build: State) -> bool {
        let map = match to_build {
            State::Out => &self.my_split_out,
            State::In => &self.my_split_in,
            State::On => &self.my_split_on,
            _ => return false,
        };
        match map.get(&shape_map_key(s)) {
            Some(losos) => losos.is_split(),
            None => false,
        }
    }

    /// `Splits(S, ToBuild)` (`.cxx:414-446`).
    pub fn splits(&self, s: &TopoShape, to_build: State) -> &[TopoShape] {
        let map = match to_build {
            State::Out => &self.my_split_out,
            State::In => &self.my_split_in,
            State::On => &self.my_split_on,
            _ => return &self.my_empty_shape_list,
        };
        match map.get(&shape_map_key(s)) {
            Some(losos) => losos.list_on_state(),
            None => &self.my_empty_shape_list,
        }
    }

    /// `ChangeSplit(S, ToBuild)` (`.cxx:449-491`).
    ///
    /// `TopAbs_UNKNOWN` answers `myEmptyShapeList` exactly as OCCT does, i.e. a
    /// member list shared by every such call.
    pub fn change_split(&mut self, s: &TopoShape, to_build: State) -> &mut Vec<TopoShape> {
        let key = shape_map_key(s);
        let map = match to_build {
            State::Out => &mut self.my_split_out,
            State::In => &mut self.my_split_in,
            State::On => &mut self.my_split_on,
            _ => return &mut self.my_empty_shape_list,
        };
        map.entry(key).or_default().change_list_on_state()
    }

    /// `IsMerged(S, ToBuild)` (`Merge.cxx:626-657`).
    pub fn is_merged(&self, s: &TopoShape, to_build: State) -> bool {
        let map = match to_build {
            State::Out => &self.my_merged_out,
            State::In => &self.my_merged_in,
            State::On => &self.my_merged_on,
            _ => return false,
        };
        match map.get(&shape_map_key(s)) {
            Some(losos) => !losos.list_on_state().is_empty(),
            None => false,
        }
    }

    /// `Merged(S, ToBuild)` (`Merge.cxx:663-693`).
    pub fn merged(&self, s: &TopoShape, to_build: State) -> &[TopoShape] {
        let map = match to_build {
            State::Out => &self.my_merged_out,
            State::In => &self.my_merged_in,
            State::On => &self.my_merged_on,
            _ => return &self.my_empty_shape_list,
        };
        match map.get(&shape_map_key(s)) {
            Some(losos) => losos.list_on_state(),
            None => &self.my_empty_shape_list,
        }
    }

    /// `ChangeMerged(S, ToBuild)` (`Merge.cxx:700-731`).
    pub fn change_merged(&mut self, s: &TopoShape, to_build: State) -> &mut Vec<TopoShape> {
        let key = shape_map_key(s);
        let map = match to_build {
            State::Out => &mut self.my_merged_out,
            State::In => &mut self.my_merged_in,
            State::On => &mut self.my_merged_on,
            _ => return &mut self.my_empty_shape_list,
        };
        map.entry(key).or_default().change_list_on_state()
    }

    /// `MSplit(s)` (`.cxx:885-901`): the default for an unknown state is
    /// `mySplitIN`.
    pub fn m_split(&self, s: State) -> &DShapeOn1StateMap {
        match s {
            State::Out => &self.my_split_out,
            State::On => &self.my_split_on,
            _ => &self.my_split_in,
        }
    }

    /// `ChangeMSplit(s)` (`.cxx:905-921`): the default for an unknown state is
    /// `mySplitIN`.
    pub fn change_m_split(&mut self, s: State) -> &mut DShapeOn1StateMap {
        match s {
            State::Out => &mut self.my_split_out,
            State::On => &mut self.my_split_on,
            _ => &mut self.my_split_in,
        }
    }

    /// `ShapePosition(S, LS)` (`.cxx:499-531`).
    pub fn shape_position(&mut self, s: &TopoShape, ls: &[TopoShape]) -> State {
        let mut state = State::Unknown;

        // `myEdgeAvoid` is the avoid list for a face classification.
        let p_los: &[TopoShape] = if s.shape_type() == ShapeType::Face {
            &self.my_edge_avoid
        } else {
            &self.my_empty_shape_list
        };

        for sls in ls {
            let st = self
                .my_shape_classifier
                .borrow_mut()
                .state_shape_shape_av_ls(s, p_los, sls);
            match st {
                Ok(st) => state = st,
                // `StateShapeShape` throws `Standard_ConstructionError` on a
                // null reference; the port keeps the previous state.
                Err(_) => state = State::Unknown,
            }
            if state != State::Out && state != State::Unknown {
                return state;
            }
        }
        if state == State::Unknown {
            return state;
        }
        State::Out
    }

    /// `KeepShape(S1, LS2, ToBuild1)` (`.cxx:540-554`).
    pub fn keep_shape(&mut self, s1: &TopoShape, ls2: &[TopoShape], to_build1: State) -> bool {
        if !ls2.is_empty() {
            let pos2 = self.shape_position(s1, ls2);
            if pos2 != to_build1 {
                return false;
            }
        }
        true
    }

    /// `TopType(S)` (`.cxx:560-615`): the first non-empty type, `COMPOUND`
    /// first and `SHAPE` last.
    pub fn top_type(s: &TopoShape) -> ShapeType {
        for t in [
            ShapeType::Compound,
            ShapeType::CompSolid,
            ShapeType::Solid,
            ShapeType::Shell,
            ShapeType::Face,
            ShapeType::Wire,
            ShapeType::Edge,
            ShapeType::Vertex,
        ] {
            if !explorer_shapes(s, t).is_empty() {
                return t;
            }
        }
        ShapeType::Shape
    }

    /// `Reverse(ToBuild1, ToBuild2)` (`.cxx:621-633`): `IN`/`IN` does not
    /// reverse, any other pair reverses when `ToBuild1` is `IN`.
    pub fn reverse(to_build1: State, to_build2: State) -> bool {
        if to_build1 == State::In && to_build2 == State::In {
            false
        } else {
            to_build1 == State::In
        }
    }

    /// `Orient(Ori, Reverse)` (`.cxx:637-640`).
    pub fn orient(o: Orientation, reverse: bool) -> Orientation {
        if !reverse {
            o
        } else {
            o.complement()
        }
    }

    /// `Contains(S, L)` (`.cxx:836-849`): a `TopoDS_Shape::IsSame` scan.
    pub fn contains(s: &TopoShape, l: &[TopoShape]) -> bool {
        l.iter().any(|sl| is_same(sl, s))
    }

    /// `FindSameDomain(L1, L2)` (`.cxx:650-708`).
    ///
    /// Both lists grow until no new same-domain shape appears; new shapes are
    /// prepended, as OCCT does.
    pub fn find_same_domain(&self, l1: &mut Vec<TopoShape>, l2: &mut Vec<TopoShape>) {
        let Some(hds) = self.my_data_structure.as_ref() else {
            return;
        };
        let mut nl1 = l1.len();
        let mut nl2 = l2.len();

        while nl1 > 0 || nl2 > 0 {
            for i in 0..nl1 {
                let s1 = l1[i].clone();
                for s2 in hds.same_domain(&s1) {
                    if !Self::contains(s2, l2) {
                        l2.insert(0, s2.clone());
                        nl2 += 1;
                    }
                }
            }
            nl1 = 0;

            for i in 0..nl2 {
                let s2 = l2[i].clone();
                for s1 in hds.same_domain(&s2) {
                    if !Self::contains(s1, l1) {
                        l1.insert(0, s1.clone());
                        nl1 += 1;
                    }
                }
            }
            nl2 = 0;
        }
    }

    /// `FindSameDomainSameOrientation(L1, L2)` (`.cxx:713-762`).
    ///
    /// The reference orientation is that of the first `L1` element, taken
    /// after `FindSameDomain` prepended to it.
    pub fn find_same_domain_same_orientation(
        &self,
        l1: &mut Vec<TopoShape>,
        l2: &mut Vec<TopoShape>,
    ) {
        self.find_same_domain(l1, l2);
        let Some(sref) = l1.first().cloned() else {
            return;
        };
        let Some(hds) = self.my_data_structure.as_ref() else {
            return;
        };
        let oref = hds.same_domain_orientation(&sref);

        let mut ll1: Vec<TopoShape> = Vec::new();
        let mut ll2: Vec<TopoShape> = Vec::new();

        for s in l1.iter().chain(l2.iter()) {
            let o = hds.same_domain_orientation(s);
            if o == oref && !Self::contains(s, &ll1) {
                ll1.push(s.clone());
            } else if o != oref && !Self::contains(s, &ll2) {
                ll2.push(s.clone());
            }
        }

        *l1 = ll1;
        *l2 = ll2;
    }

    /// `MapShapes(S1, S2)` (`.cxx:766-779`).
    pub fn map_shapes(&mut self, s1: Option<&TopoShape>, s2: Option<&TopoShape>) {
        self.clear_maps();
        if let Some(s1) = s1 {
            map_shapes_into(s1, &mut self.my_map1);
        }
        if let Some(s2) = s2 {
            map_shapes_into(s2, &mut self.my_map2);
        }
    }

    /// `ClearMaps()` (`.cxx:783-787`).
    pub fn clear_maps(&mut self) {
        self.my_map1.clear();
        self.my_map2.clear();
    }

    /// `FindSameRank(L1, rank, L2)` (`.cxx:791-804`).
    pub fn find_same_rank(&self, l1: &[TopoShape], rank: i32, l2: &mut Vec<TopoShape>) {
        for s in l1 {
            let r = self.shape_rank(s);
            if r == rank && !Self::contains(s, l2) {
                l2.push(s.clone());
            }
        }
    }

    /// `ShapeRank(s)` (`.cxx:808-814`).
    pub fn shape_rank(&self, s: &TopoShape) -> i32 {
        let isof1 = self.is_shape_of(s, 1);
        let isof2 = self.is_shape_of(s, 2);
        if isof1 || isof2 {
            if isof1 {
                1
            } else {
                2
            }
        } else {
            0
        }
    }

    /// `IsShapeOf(s, i)` (`.cxx:818-830`).
    pub fn is_shape_of(&self, s: &TopoShape, i: i32) -> bool {
        if i == 1 {
            self.my_map1.contains(s)
        } else if i == 2 {
            self.my_map2.contains(s)
        } else {
            false
        }
    }
}
