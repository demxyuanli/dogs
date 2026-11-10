//! `TopOpeBRepBuild_Builder1` — the SOLID/SOLID rebuild that keeps the state
//! of every sub-shape.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_Builder1.hxx:43-203`,
//! `TopOpeBRepBuild_Builder1.cxx:63-182` (OCCT 8.0.0).
//!
//! `Builder1` is the subclass `TopOpeBRepBuild_HBuilder` owns
//! (`HBuilder.hxx:170`); the port keeps it as a wrapper around [`Builder`]
//! because OCCT's inheritance only ever means "call the base member with the
//! same name, or extend it". `my_base` is that base sub-object and
//! [`Builder1::base`] / [`Builder1::base_mut`] expose it.
//!
//! Ported here:
//! - the constructor (`cxx:65-73`), `Clear` (`cxx:86-97`), `Perform` both
//!   overloads (`cxx:101-104`, `cxx:108-141`), `MergeKPart` both overloads
//!   (`cxx:145-148`, `cxx:152-182`),
//! - [`shape_states`]: `PerformShapeWithStates` both overloads and
//!   `PerformFacesWithStates` (`Builder1_1.cxx:183-552`),
//! - [`status_edges`]: `StatusEdgesToSplit`, `SplitEdge` and the
//!   `ClassifyEdgeToSolidByOnePoint` helper (`Builder1_1.cxx:556-1144`),
//! - [`correct_result_2d`]: `CorrectResult2d` (`Builder1_2.cxx:33-224`),
//! - [`on_parts`]: `PerformONParts` (`cxx:1151-1372`),
//! - [`fill_2d`], [`pieces_on`], [`is_same_2d`]: `PerformPieceIn2D`
//!   (`cxx:1415-2226`), `PerformPieceOn2D`, `TwoPiecesON`, `IsSame2d` and
//!   `OrientateEdgeOnFace` (`cxx:2227-2290`),
//! - [`gfill_not_same_dom`]: `GFillSolidSFS` (`cxx:186-241`), `GFillShellSFS`
//!   (`cxx:243-301`), `GFillFaceNotSameDomSFS` (`:303-378`),
//!   `GFillFaceNotSameDomWES` (`:380-430`), `GFillWireNotSameDomWES`
//!   (`:431-477`), `GFillEdgeNotSameDomWES` (`:478-532`),
//! - [`gfill_same_dom`]: `GFillFaceSameDomSFS` (`:533-613`),
//!   `GFillFaceSameDomWES` (`:614-746`), `GFillWireSameDomWES` (`:747-817`),
//!   `GFillEdgeSameDomWES` (`:818-1150`),
//! - [`gwes_make_faces`]: the `GWESMakeFaces` override (`:1373-1412`).
//!
//! Virtual dispatch: three members of `TopOpeBRepBuild_Builder` are virtual and
//! overridden here — `GFillSolidSFS` (`hxx:69-73`), `GFillShellSFS`
//! (`hxx:76-79`) and `GWESMakeFaces` (`hxx:80-82`). OCCT reaches the override
//! through the vtable; the port reaches it through
//! `Builder::g_fill_solid_sfs`, `Builder::g_fill_shell_sfs` and
//! `Builder::gwes_make_faces`, which branch on [`Builder::is_builder1`] and call
//! the `b1_*` member or the base one, exactly where the vtable would land.
//! `MergeShapes` (`Merge.cxx:203-208`) calls `MergeKPart()` virtually as well,
//! so every boolean built through `Builder1` lands in the override below.
//!
//! UNPORTED, with their OCCT location:
//! - `GFillSplitsPVS` (`hxx:83-85`): declared in `Builder1.hxx` with no
//!   definition anywhere in OCCT 8.0.0.
//!
//! `myFSplits` / `myESplits`: `Builder1.hxx:186-188` declares them a second
//! time, shadowing the base members of the same name. The port keeps two
//! separate maps so a shape split through `Builder1` and through `Builder` can
//! never alias.

use std::cell::{Cell, RefCell};

use crate::abs::State;
use crate::shape::TopoShape;
use crate::topo_tools::ShapeMapKey;
use crate::topopebrep::ds::build_tool::BuildTool;
use crate::topopebrep::ds::data_structure::IndexedShapeSet;
use crate::topopebrep::ds::h_data_structure::HDataStructure;
use crate::topopebrep::ds::tool::DShapeOn1StateMap;

use super::k_part::{KP_ISAFA, KP_ISDISJ, KP_ISKOLE, KP_ISSOSO, KP_ISKOLETGE};
use super::Builder;

mod correct_result_2d;
mod fill_2d;
mod gfill_not_same_dom;
mod gfill_same_dom;
mod gwes_make_faces;
mod is_same_2d;
mod on_parts;
mod pieces_on;
mod shape_states;
mod state;
mod status_edges;

pub(crate) use state::Builder1State;

thread_local! {
    /// `bool GLOBAL_USE_NEW_BUILDER` (`Builder1.cxx:61`), false at start-up:
    /// `BRepAlgoAPI_BooleanOperation` sets it, LocOpe / features leave it
    /// alone and get the classic `Builder` behaviour.
    static GLOBAL_USE_NEW_BUILDER: Cell<bool> = const { Cell::new(false) };

    /// `static IndexedMap processedEdges` (`Builder1_1.cxx:55`).
    static PROCESSED_EDGES: RefCell<IndexedShapeSet> =
        RefCell::new(IndexedShapeSet::new());
    /// `static IndexedMap theUsedVertexMap` (`Builder1_1.cxx:56`).
    static THE_USED_VERTEX_MAP: RefCell<IndexedShapeSet> =
        RefCell::new(IndexedShapeSet::new());
    /// `static IndexedMap theUnkStateVer` (`Builder1_1.cxx:57`), a
    /// `TopTools_MapOfShape` in OCCT; a set here as well.
    static THE_UNK_STATE_VER: RefCell<IndexedShapeSet> =
        RefCell::new(IndexedShapeSet::new());

    /// `static IndexedMap mySDEdgeMap` (`Builder1.cxx:51`): every edge of the
    /// same-domain faces currently being filled. Filled by
    /// `GFillFaceSameDomSFS` (`cxx:541`, `:638`) and read by `PerformONParts`
    /// (`cxx:1199`).
    static MY_SD_EDGE_MAP: RefCell<IndexedShapeSet> =
        RefCell::new(IndexedShapeSet::new());
}

/// `mySDEdgeMap.Clear()` (`Builder1.cxx:541`).
fn sd_edge_map_clear() {
    MY_SD_EDGE_MAP.with(|m| m.borrow_mut().clear());
}

/// `TopExp::MapShapes(SDF, TopAbs_EDGE, mySDEdgeMap)` (`Builder1.cxx:638`).
fn sd_edge_map_add_shapes(sdf: &TopoShape) {
    MY_SD_EDGE_MAP.with(|m| {
        let mut m = m.borrow_mut();
        for e in crate::topo_tools::explorer_shapes(sdf, crate::abs::ShapeType::Edge) {
            m.add(&e);
        }
    });
}

/// `mySDEdgeMap.Contains(EG)` (`Builder1.cxx:1199`).
fn sd_edge_map_contains(e: &TopoShape) -> bool {
    MY_SD_EDGE_MAP.with(|m| m.borrow().contains(e))
}

/// Current `GLOBAL_USE_NEW_BUILDER` (`Builder1.cxx:61`).
pub fn global_use_new_builder() -> bool {
    GLOBAL_USE_NEW_BUILDER.with(|c| c.get())
}

/// Assign `GLOBAL_USE_NEW_BUILDER` (`Builder1.cxx:61`).
pub fn set_global_use_new_builder(v: bool) {
    GLOBAL_USE_NEW_BUILDER.with(|c| c.set(v));
}

/// `processedEdges.Contains(S)` (`Builder1_1.cxx:587`).
fn processed_edges_contains(s: &TopoShape) -> bool {
    PROCESSED_EDGES.with(|m| m.borrow().contains(s))
}

/// `processedEdges.Add(S)` (`Builder1_1.cxx:620`).
fn processed_edges_add(s: &TopoShape) {
    PROCESSED_EDGES.with(|m| m.borrow_mut().add(s));
}

/// `processedEdges.Clear()` (`Builder1_1.cxx:337`, `:339`).
fn processed_edges_clear() {
    PROCESSED_EDGES.with(|m| m.borrow_mut().clear());
}

/// `theUsedVertexMap.Add(V)` (`Builder1_1.cxx:202` and the other sites).
fn used_vertex_add(v: &TopoShape) {
    THE_USED_VERTEX_MAP.with(|m| m.borrow_mut().add(v));
}

/// `theUsedVertexMap.Contains(V)` (`Builder1_1.cxx:932`).
fn used_vertex_contains(v: &TopoShape) -> bool {
    THE_USED_VERTEX_MAP.with(|m| m.borrow().contains(v))
}

/// `theUsedVertexMap.Clear()` (`Builder1_1.cxx:185`).
fn used_vertex_clear() {
    THE_USED_VERTEX_MAP.with(|m| m.borrow_mut().clear());
}

/// `theUnkStateVer.Add(V)` (`Builder1_1.cxx:749`, `:756`, `:797`, `:804`,
/// `:991`, `:993`).
fn unk_state_ver_add(v: &TopoShape) {
    THE_UNK_STATE_VER.with(|m| m.borrow_mut().add(v));
}

/// `theUnkStateVer.Clear()` (`Builder1_1.cxx:186`).
fn unk_state_ver_clear() {
    THE_UNK_STATE_VER.with(|m| m.borrow_mut().clear());
}

/// `const double PAR_T = 0.43213918` (`Builder1_1.cxx:45`): the parameter
/// division of `ClassifyEdgeToSolidByOnePoint`.
const PAR_T: f64 = 0.43213918;

/// `class TopOpeBRepBuild_Builder1 : public TopOpeBRepBuild_Builder`
/// (`hxx:43`).
///
/// The subclass state lives in [`Builder::b1`]; this wrapper only adds the
/// `my_base` sub-object OCCT gets from inheritance.
pub struct Builder1 {
    /// The `TopOpeBRepBuild_Builder` sub-object.
    my_base: Builder,
}

impl Builder1 {
    /// `TopOpeBRepBuild_Builder1(BT)` (`cxx:65-73`).
    pub fn new(bt: BuildTool) -> Self {
        let mut b = Self {
            my_base: Builder::new(bt),
        };
        b.my_base.is_builder1 = true;
        // `cxx:68-73`
        b.my_base.b1 = Builder1State::new();
        b.my_base.b1.my_same_dom_map.clear();
        b.my_base.b1.my_map_of_edge_faces.clear();
        b.my_base.b1.my_splits_on_to_keep.clear();
        b.my_base.b1.my_processed_parts_out_2d.clear();
        b.my_base.b1.my_processed_parts_on_2d.clear();
        b
    }

    /// The `TopOpeBRepBuild_Builder` sub-object.
    pub fn base(&self) -> &Builder {
        &self.my_base
    }

    /// The `TopOpeBRepBuild_Builder` sub-object, mutable.
    pub fn base_mut(&mut self) -> &mut Builder {
        &mut self.my_base
    }

    /// `Clear()` (`cxx:86-97`): only the base members are cleared; every
    /// `Builder1` map is commented out in OCCT, so they survive.
    pub fn clear(&mut self) {
        self.my_base.clear();
    }

    /// `Perform(HDS)` (`cxx:101-104`).
    pub fn perform(&mut self, hds: HDataStructure) {
        self.my_base.perform(hds);
    }

    /// `Perform(HDS, S1, S2)` (`cxx:108-141`).
    pub fn perform_with_shapes(
        &mut self,
        hds: HDataStructure,
        s1: Option<&TopoShape>,
        s2: Option<&TopoShape>,
    ) {
        // `cxx:114-118`: LocOpe and the mechanical features keep the classic
        // builder.
        if !global_use_new_builder() {
            self.my_base.perform_with_shapes(hds, s1, s2);
            return;
        }

        // `cxx:120-124`
        self.my_base.b1.my_same_dom_map.clear();
        self.my_base.b1.my_map_of_edge_faces.clear();
        self.my_base.b1.my_splits_on_to_keep.clear();
        self.my_base.b1.my_processed_parts_out_2d.clear();
        self.my_base.b1.my_processed_parts_on_2d.clear();

        // `cxx:126-128`
        self.my_base.my_shape1 = s1.cloned();
        self.my_base.my_shape2 = s2.cloned();
        self.perform(hds);

        // `cxx:130-134`
        let mut kp = self.my_base.find_is_k_part();
        if kp == KP_ISKOLE || kp == KP_ISKOLETGE {
            kp = KP_ISSOSO;
        }
        self.my_base.my_is_k_part = kp;

        // `cxx:136-140`: only SOLID/SOLID gets the per-shape states.
        if kp == KP_ISSOSO {
            self.my_base.perform_shape_with_states();
        }
    }

    /// `MergeKPart(TB1, TB2)` (`cxx:145-148`).
    pub fn merge_k_part_states(&mut self, tb1: State, tb2: State) {
        self.my_base.merge_k_part_states(tb1, tb2);
    }

    /// `MergeKPart()` (`cxx:152-182`): the five k-part cases, with the
    /// `issoso` case additionally running `CorrectResult2d` on every merged
    /// solid, then `End()`.
    ///
    /// `Merged` answers a `const` list in OCCT and the `Iterator` copies it,
    /// so `its.ChangeValue()` edits a temporary; the port's `to_vec()` keeps
    /// that same "the list stored in the builder is untouched" behaviour.
    pub fn merge_k_part(&mut self) {
        let kp = self.my_base.my_is_k_part;
        if kp == KP_ISKOLE {
            self.my_base.merge_k_part_iskole();
        } else if kp == KP_ISKOLETGE {
            self.my_base.merge_k_part_iskoletge();
        } else if kp == KP_ISDISJ {
            self.my_base.merge_k_part_isdisj();
        } else if kp == KP_ISAFA {
            self.my_base.merge_k_part_isfafa();
        } else if kp == KP_ISSOSO {
            self.my_base.merge_k_part_issoso();

            // `cxx:174-178`
            let shape1 = self.my_base.my_shape1.clone();
            let state1 = self.my_base.my_state1;
            if let Some(shape1) = shape1 {
                let merged = self.my_base.merged(&shape1, state1).to_vec();
                for mut r in merged {
                    self.my_base.correct_result_2d(&mut r);
                }
            }
        }
        // `cxx:181`
        self.my_base.end();
    }

    // -----------------------------------------------------------------
    // The inherited members `TopOpeBRepBuild_HBuilder` calls. Each one is a
    // `TopOpeBRepBuild_Builder` member with the same name and signature.
    // -----------------------------------------------------------------

    /// `DataStructure()` (base).
    pub fn data_structure(&self) -> Option<&HDataStructure> {
        self.my_base.data_structure()
    }

    /// `ChangeBuildTool()` (base).
    pub fn change_build_tool(&mut self) -> &mut BuildTool {
        self.my_base.change_build_tool()
    }

    /// `BuildTool()` (base).
    pub fn build_tool(&self) -> &BuildTool {
        self.my_base.build_tool()
    }

    /// `MergeShapes(...)` (base).
    pub fn merge_shapes(
        &mut self,
        s1: Option<&TopoShape>,
        to_build1: State,
        s2: Option<&TopoShape>,
        to_build2: State,
    ) {
        self.my_base.merge_shapes(s1, to_build1, s2, to_build2);
    }

    /// `MergeSolids(...)` (base).
    pub fn merge_solids(
        &mut self,
        s1: &TopoShape,
        to_build1: State,
        s2: &TopoShape,
        to_build2: State,
    ) {
        self.my_base.merge_solids(s1, to_build1, s2, to_build2);
    }

    /// `MergeSolid(...)` (base).
    pub fn merge_solid(&mut self, s: &TopoShape, to_build: State) {
        self.my_base.merge_solid(s, to_build);
    }

    /// `IsSplit(...)` (base).
    pub fn is_split(&self, s: &TopoShape, to_build: State) -> bool {
        self.my_base.is_split(s, to_build)
    }

    /// `Splits(...)` (base).
    pub fn splits(&self, s: &TopoShape, to_build: State) -> &[TopoShape] {
        self.my_base.splits(s, to_build)
    }

    /// `IsMerged(...)` (base).
    pub fn is_merged(&self, s: &TopoShape, to_build: State) -> bool {
        self.my_base.is_merged(s, to_build)
    }

    /// `Merged(...)` (base).
    pub fn merged(&self, s: &TopoShape, to_build: State) -> &[TopoShape] {
        self.my_base.merged(s, to_build)
    }

    /// `NewVertex(I)` (base).
    pub fn new_vertex(&self, i: i32) -> Option<&TopoShape> {
        self.my_base.new_vertex(i)
    }

    /// `NewEdges(I)` (base).
    pub fn new_edges(&self, i: i32) -> &[TopoShape] {
        self.my_base.new_edges(i)
    }

    /// `ChangeNewEdges(I)` (base).
    pub fn change_new_edges(&mut self, i: i32) -> &mut Vec<TopoShape> {
        self.my_base.change_new_edges(i)
    }

    /// `NewFaces(I)` (base).
    pub fn new_faces(&self, i: i32) -> Option<&[TopoShape]> {
        self.my_base.new_faces(i)
    }

    /// `Section()` (base).
    pub fn section(&mut self) -> &[TopoShape] {
        self.my_base.section()
    }

    /// `SectionCurves(lse)` (base).
    pub fn section_curves_into(&mut self, lse: &mut Vec<TopoShape>) {
        self.my_base.section_curves_into(lse);
    }

    /// `SectionEdges(lse)` (base).
    pub fn section_edges_into(&mut self, lse: &mut Vec<TopoShape>) {
        self.my_base.section_edges_into(lse);
    }

    /// `Section(l)` (base).
    pub fn section_list(&mut self, l: &mut Vec<TopoShape>) {
        self.my_base.section_list(l);
    }

    /// `MSplit(s)` (base).
    pub fn m_split(&self, s: State) -> &DShapeOn1StateMap {
        self.my_base.m_split(s)
    }

    /// The shape a `mySplit*` key stands for (base helper).
    pub fn shape_of_key(&self, key: &ShapeMapKey) -> Option<TopoShape> {
        self.my_base.shape_of_key(key)
    }

    /// `IsKPart()` (base, `KPart.cxx:169-172`).
    pub fn is_k_part(&self) -> i32 {
        self.my_base.is_k_part()
    }

    /// `End()` (base, `End.cxx:45-...`).
    pub fn end(&mut self) {
        self.my_base.end();
    }
}
