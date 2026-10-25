//! `TopOpeBRepBuild_FaceBuilder` — rebuilds the faces of a `WireEdgeSet`.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_FaceBuilder.hxx:31-140`,
//! `.cxx:58-661` (OCCT 8.0.0).
//!
//! `InitFaceBuilder` (`cxx:72-82`) turns the `WireEdgeSet` into a `LoopSet` of
//! the old wires plus the blocks of edges (`MakeLoops`, `cxx:634-661`), then
//! runs the `FaceAreaBuilder` with the `WireEdgeClassifier`.
//!
//! The three repair passes (`DetectUnclosedWire` `cxx:205-346`,
//! `CorrectGclosedWire` `cxx:348-384`, `DetectPseudoInternalEdge` `cxx:388-448`)
//! invalidate edges of the `BlockBuilder`; only the non-`RM_HANGING` arm is
//! translated, exactly as OCCT builds it (`cxx:44`, `:47-49`).

use std::collections::HashMap;
use std::sync::{Arc, RwLock};

use occt_core::precision::CONFUSION;

use crate::abs::{Orientation, ShapeType};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, TopoShape, Vertex};
use crate::topo_tools::{
    explorer_shapes, explorer_shapes_avoiding, is_same, shape_map_key, ShapeMapKey,
};
use crate::topopebrep::ds::build_tool::BuildTool;
use crate::topopebrep::ds::data_structure::{IndexedShapeData, IndexedShapeSet};

use super::block_builder::BlockBuilder;
use super::block_iterator::BlockIterator;
use super::face_area_builder::FaceAreaBuilder;
use super::loop_set::{LoopSet, LoopSetSource};
use super::loop_shape::{Loop, LoopHandle};
use super::shape_set::ShapeSetSource;
use super::wire_edge_classifier::WireEdgeClassifier;
use super::wire_edge_set::WireEdgeSet;

/// `#define ISUNKNOWN` (`cxx:118`).
const ISUNKNOWN: i32 = -1;
/// `#define ISVERTEX` (`cxx:119`).
const ISVERTEX: i32 = 0;
/// `#define GCLOSEDW` (`cxx:120`).
const GCLOSEDW: i32 = 1;
/// `#define UNCLOSEDW` (`cxx:121`).
const UNCLOSEDW: i32 = 2;
/// `#define CLOSEDW` (`cxx:122`).
const CLOSEDW: i32 = 10;

/// `class TopOpeBRepBuild_FaceBuilder` (`hxx:37`).
pub struct FaceBuilder {
    /// `myFace` (`hxx:131`).
    my_face: Option<TopoShape>,
    /// `myLoopSet` (`hxx:132`).
    my_loop_set: LoopSet,
    /// `myBlockIterator` (`hxx:133`).
    my_block_iterator: BlockIterator,
    /// `myBlockBuilder` (`hxx:134`).
    my_block_builder: BlockBuilder,
    /// `myFaceAreaBuilder` (`hxx:135`).
    my_face_area_builder: FaceAreaBuilder,
    /// `myMOSI` (`hxx:136`).
    my_mosi: HashMap<ShapeMapKey, i32>,
}

impl FaceBuilder {
    /// `TopOpeBRepBuild_FaceBuilder()` (`cxx:58`).
    pub fn new() -> Self {
        Self {
            my_face: None,
            my_loop_set: LoopSet::new(),
            my_block_iterator: BlockIterator::new(),
            my_block_builder: BlockBuilder::new(),
            my_face_area_builder: FaceAreaBuilder::new(),
            my_mosi: HashMap::new(),
        }
    }

    /// `TopOpeBRepBuild_FaceBuilder(WES, F, ForceClass)` (`cxx:63-68`).
    pub fn with_wire_edge_set(
        wes: &mut WireEdgeSet,
        f: &TopoShape,
        force_class: bool,
    ) -> Self {
        let mut b = Self::new();
        b.init_face_builder(wes, f, force_class);
        b
    }

    /// `InitFaceBuilder(WES, F, ForceClass)` (`cxx:72-82`).
    pub fn init_face_builder(
        &mut self,
        wes: &mut WireEdgeSet,
        f: &TopoShape,
        force_class: bool,
    ) {
        // `cxx:74`: `myFace = TopoDS::Face(F)`
        assert!(f.shape_type() == ShapeType::Face, "TopOpeBRepBuild_FaceBuilder: not a face");
        self.my_face = Some(f.clone());
        self.make_loops(wes);
        let mut wec = WireEdgeClassifier::new(f, &self.my_block_builder);
        self.my_face_area_builder
            .init_face_area_builder(&mut self.my_loop_set, &mut wec, force_class);
    }

    /// `MakeLoops(SS)` (`cxx:634-661`).
    fn make_loops(&mut self, ss: &mut dyn ShapeSetSource) {
        // Build blocks on elements of SS
        self.my_block_builder.make_block(ss);

        let mut loops: Vec<LoopHandle> = Vec::new();

        // Add shapes of SS as shape loops
        ss.init_shapes();
        while ss.more_shapes() {
            let Some(s) = ss.shape().cloned() else {
                break;
            };
            loops.push(Arc::new(RwLock::new(Loop::new(&s))));
            ss.next_shape();
        }

        // Add blocks of BB as block loops
        self.my_block_builder.init_block();
        while self.my_block_builder.more_block() {
            let bi = self.my_block_builder.block_iterator();
            loops.push(Arc::new(RwLock::new(Loop::with_block_iterator(bi))));
            self.my_block_builder.next_block();
        }

        let ll = self.my_loop_set.change_list_of_loop();
        ll.clear();
        ll.extend(loops);
    }

    /// `DetectUnclosedWire(mapVVsameG, mapVon1Edge)` (`cxx:205-346`).
    pub fn detect_unclosed_wire(
        &mut self,
        map_v_v_same_g: &mut IndexedShapeData<TopoShape>,
        map_v_on1_edge: &mut IndexedShapeData<TopoShape>,
    ) {
        map_v_v_same_g.clear();
        map_v_on1_edge.clear();

        self.init_face();
        while self.more_face() {
            self.init_wire();
            while self.more_wire() {
                if self.is_old_wire() {
                    self.next_wire();
                    continue;
                }

                let builder = TopoBuilder::new();
                let mut cmp = builder.make_shape(ShapeType::Compound);
                self.init_edge();
                while self.more_edge() {
                    let e = self.edge().clone();
                    self.add_edge_wire(&e, &mut cmp);
                    self.next_edge();
                }
                let w = cmp;

                // <mapVon1E> binds vertices of connexity 1 attached to one
                // non-closed, non-degenerated edge.
                let mut map_von1e: IndexedShapeData<TopoShape> = IndexedShapeData::new();
                fun_detect_vertices_on1_edge(&w, &mut map_von1e);

                let mut map_vv: IndexedShapeData<TopoShape> = IndexedShapeData::new();
                let res = fun_analyze_map_von1e(&map_von1e, &mut map_vv);

                if res == ISVERTEX {
                    self.next_wire();
                    continue;
                } else if res == CLOSEDW {
                    self.next_wire();
                    continue;
                } else if res == GCLOSEDW {
                    for i in 1..=map_vv.extent() {
                        let key = map_vv.find_key(i).cloned();
                        let val = map_vv.find_from_index(i).cloned();
                        if let (Some(key), Some(val)) = (key, val) {
                            map_v_v_same_g.add(&key, val);
                        }
                    }
                    for i in 1..=map_von1e.extent() {
                        let key = map_von1e.find_key(i).cloned();
                        let val = map_von1e.find_from_index(i).cloned();
                        if let (Some(key), Some(val)) = (key, val) {
                            map_v_on1_edge.add(&key, val);
                        }
                    }
                } else if res == UNCLOSEDW {
                    // `cxx:333-341` (the `#else` arm of `RM_HANGING`).
                    for e in explorer_shapes(&w, ShapeType::Edge) {
                        let i = self.my_block_builder.element_index(&e);
                        self.my_block_builder.set_valid(i, false);
                    }
                }
                self.next_wire();
            }
            self.next_face();
        }
    }

    /// `CorrectGclosedWire(mapVVref, mapVon1Edge)` (`cxx:348-384`).
    ///
    /// The `mapVVsameG` / `mapVon1Edge` pair filled by
    /// [`Self::detect_unclosed_wire`] is the input.
    pub fn correct_gclosed_wire(
        &self,
        map_v_v_ref: &IndexedShapeData<TopoShape>,
        map_v_on1_edge: &IndexedShapeData<TopoShape>,
    ) {
        let builder = TopoBuilder::new();
        let n_vv = map_v_v_ref.extent();
        for i in 1..=n_vv {
            let (Some(v), Some(v_ref)) = (
                map_v_v_ref.find_key(i).cloned(),
                map_v_v_ref.find_from_index(i).cloned(),
            ) else {
                continue;
            };

            if is_same(&v, &v_ref) {
                continue;
            }

            let Some(e) = map_v_on1_edge.find_from_key(&v).cloned() else {
                continue;
            };
            let paron_e = BRepTool::parameter_of_vertex(&e, &v);

            // `cxx:373-378`: `E.Free(true); BB.Remove(E, V);`
            e.set_free(true);
            let mut e_rem = e.clone();
            builder.remove(&mut e_rem, &v);

            // `cxx:377-378`: `aLocalShape = Vref.Oriented(V.Orientation())`
            let a_local_shape = v_ref.oriented(v.orientation());
            let new_v_ref = a_local_shape;
            // `cxx:381`: `BB.Add(E, newVref)`; `cxx:382-383`:
            // `TopOpeBRepDS_BuildTool BT; BT.Parameter(E, newVref, paronE);`
            builder.add(&mut e_rem, &new_v_ref);
            BuildTool::new().parameter(&e_rem, &new_v_ref, paron_e);
        }
    }

    /// `DetectPseudoInternalEdge(MapE)` (`cxx:388-448`).
    pub fn detect_pseudo_internal_edge(&mut self, map_e: &mut IndexedShapeSet) {
        let builder = TopoBuilder::new();
        let mut cmp = builder.make_shape(ShapeType::Compound);
        self.init_face();
        while self.more_face() {
            self.init_wire();
            while self.more_wire() {
                if self.is_old_wire() {
                    self.next_wire();
                    continue;
                }
                self.init_edge();
                while self.more_edge() {
                    let e = self.edge().clone();
                    self.add_edge_wire(&e, &mut cmp);
                    self.next_edge();
                }
                self.next_wire();
            }
            self.next_face();
        }

        let map_voe = map_shapes_and_ancestors_ordered(&cmp, ShapeType::Vertex, ShapeType::Edge);
        let nv = map_voe.extent();

        map_e.clear();
        for i in 1..=nv {
            let Some(le) = map_voe.find_from_index(i) else {
                continue;
            };
            if le.len() == 2 {
                let e1 = &le[0];
                let e2 = &le[1];
                let same = is_same(e1, e2);
                let o1 = e1.orientation();
                let o2 = e2.orientation();
                let o1co2 = o1 == o2.complement();

                if same && o1co2 {
                    map_e.add(e1);

                    let ie1 = self.my_block_builder.element_index(e1);
                    self.my_block_builder.set_valid(ie1, false);

                    let ie2 = self.my_block_builder.element_index(e2);
                    self.my_block_builder.set_valid(ie2, false);
                }
            }
        }
    }

    /// `Face()` (`cxx:452-455`).
    pub fn face(&self) -> Option<&TopoShape> {
        self.my_face.as_ref()
    }

    /// `InitFace()` (`cxx:459-463`).
    pub fn init_face(&mut self) -> i32 {
        self.my_face_area_builder.init_area()
    }

    /// `MoreFace()` (`cxx:467-471`).
    pub fn more_face(&self) -> bool {
        self.my_face_area_builder.more_area()
    }

    /// `NextFace()` (`cxx:475-478`).
    pub fn next_face(&mut self) {
        self.my_face_area_builder.next_area();
    }

    /// `InitWire()` (`cxx:482-486`).
    pub fn init_wire(&mut self) -> i32 {
        self.my_face_area_builder.init_loop()
    }

    /// `MoreWire()` (`cxx:490-494`).
    pub fn more_wire(&self) -> bool {
        self.my_face_area_builder.more_loop()
    }

    /// `NextWire()` (`cxx:498-501`).
    pub fn next_wire(&mut self) {
        self.my_face_area_builder.next_loop();
    }

    /// `IsOldWire()` (`cxx:505-509`).
    pub fn is_old_wire(&self) -> bool {
        match self.current_loop() {
            Some(l) => l.read().expect("poisoned Loop lock").is_shape(),
            None => false,
        }
    }

    /// `OldWire()` (`cxx:513-518`).
    ///
    /// OCCT returns a reference into the loop; the port hands out a clone
    /// because the loop lives behind the [`LoopHandle`] lock.
    pub fn old_wire(&self) -> Option<TopoShape> {
        self.current_loop()
            .and_then(|l| l.read().expect("poisoned Loop lock").shape().cloned())
    }

    /// `FindNextValidElement()` (`cxx:524-542`).
    pub fn find_next_valid_element(&mut self) {
        // prerequisites : myBlockIterator.Initialize
        let _ = self.current_loop();
        while self.my_block_iterator.more() {
            let i = self.my_block_iterator.value();
            if self.my_block_builder.element_is_valid(i) {
                break;
            } else {
                self.my_block_iterator.next();
            }
        }
    }

    /// `InitEdge()` (`cxx:546-560`).
    ///
    /// `Standard_DomainError("TopOpeBRepBuild_FaceBuilder:InitEdge")` for a wire
    /// loop becomes the port's `panic!`.
    pub fn init_edge(&mut self) -> i32 {
        let Some(l) = self.current_loop().cloned() else {
            panic!("TopOpeBRepBuild_FaceBuilder:InitEdge");
        };
        if l.read().expect("poisoned Loop lock").is_shape() {
            panic!("TopOpeBRepBuild_FaceBuilder:InitEdge");
        }
        let bi = *l.read().expect("poisoned Loop lock").block_iterator();
        self.my_block_iterator = bi;
        self.my_block_iterator.initialize();
        self.find_next_valid_element();
        self.my_block_iterator.extent()
    }

    /// `MoreEdge()` (`cxx:564-568`).
    pub fn more_edge(&self) -> bool {
        self.my_block_iterator.more()
    }

    /// `NextEdge()` (`cxx:572-576`).
    pub fn next_edge(&mut self) {
        self.my_block_iterator.next();
        self.find_next_valid_element();
    }

    /// `Edge()` (`cxx:580-596`).
    pub fn edge(&self) -> &TopoShape {
        if !self.my_block_iterator.more() {
            panic!("OutOfRange");
        }

        let i = self.my_block_iterator.value();
        let isvalid = self.my_block_builder.element_is_valid(i);
        if !isvalid {
            panic!("Edge not Valid");
        }

        self.my_block_builder
            .element_at(i)
            .expect("TopOpeBRepBuild_FaceBuilder::Edge")
    }

    /// `EdgeConnexity(E)` (`cxx:600-610`).
    ///
    /// The `OCCT_DEBUG` build throws
    /// `Standard_ProgramError("FaceBuilder::EdgeConnexity management
    /// disactivated")`; the release body answers 0, which is what the port
    /// keeps (`myMOSI` is never filled, `cxx:609-610`).
    pub fn edge_connexity(&self, e: &TopoShape) -> i32 {
        *self.my_mosi.get(&shape_map_key(e)).unwrap_or(&0)
    }

    /// `AddEdgeWire(E, W)` (`cxx:614-629`).
    pub fn add_edge_wire(&self, e: &TopoShape, w: &mut TopoShape) -> i32 {
        let builder = TopoBuilder::new();
        builder.add(w, e);
        1
    }

    /// The current loop of the area under `myFaceAreaBuilder`.
    fn current_loop(&self) -> Option<&LoopHandle> {
        self.my_face_area_builder.loop_()
    }
}

impl Default for FaceBuilder {
    fn default() -> Self {
        Self::new()
    }
}

/// `FUN_DetectVerticesOn1Edge(W, mapVon1E)` (`cxx:84-116`).
///
/// Fills `mapVon1E` with the vertices of `W` of 3d connexity 1.
fn fun_detect_vertices_on1_edge(w: &TopoShape, map_von1e: &mut IndexedShapeData<TopoShape>) {
    let map_vedges = map_shapes_and_ancestors_ordered(w, ShapeType::Vertex, ShapeType::Edge);
    let n_v = map_vedges.extent();

    for i in 1..=n_v {
        let Some(v) = map_vedges.find_key(i).cloned() else {
            continue;
        };
        if v.orientation() == Orientation::Internal {
            continue;
        }

        let Some(lo_e) = map_vedges.find_from_index(i) else {
            continue;
        };
        if lo_e.len() < 2 {
            // Keeping INTERNAL or EXTERNAL edges
            let Some(e) = lo_e.first().cloned() else {
                continue;
            };
            let ori_e = e.orientation();
            if ori_e == Orientation::Internal || ori_e == Orientation::External {
                continue;
            }
            map_von1e.add(&v, e);
        }
    }
}

/// `FUN_AnalyzemapVon1E(mapVon1E, mapVV)` (`cxx:124-200`).
// The `ISUNKNOWN` initializer (`cxx:127`) is kept although every branch of the
// `if`/`else if`/`else` chain reassigns `res`.
#[allow(unused_assignments)]
fn fun_analyze_map_von1e(
    map_von1e: &IndexedShapeData<TopoShape>,
    map_vv: &mut IndexedShapeData<TopoShape>,
) -> i32 {
    let mut res = ISUNKNOWN;

    let n_v = map_von1e.extent();
    if n_v == 0 {
        res = CLOSEDW;
    } else if n_v == 1 {
        let e = map_von1e
            .find_from_index(1)
            .cloned()
            .expect("FUN_AnalyzemapVon1E: one edge");
        let e_closed = BRepTool::is_closed(&e);
        let dg_e = BRepTool::is_degenerated(
            &Edge::wrap(e.clone()).expect("FUN_AnalyzemapVon1E: edge"),
        );
        if dg_e {
            res = ISVERTEX;
        } else if e_closed {
            res = CLOSEDW;
        } else {
            res = UNCLOSEDW;
        }
    } else {
        // Finding among all vertices, couple of vertices falling on same
        // geometry. Filling up map <mapVV>, with (vi,vj), vi and vj are on
        // same point.
        let tol = CONFUSION;
        for i in 1..=n_v {
            let Some(vi) = map_von1e.find_key(i).cloned() else {
                continue;
            };
            let vi_p = Vertex::wrap(vi.clone())
                .map(|v| BRepTool::vertex_point(&v))
                .expect("FUN_AnalyzemapVon1E: vertex");
            for j in (i + 1)..=n_v {
                let Some(vj) = map_von1e.find_key(j).cloned() else {
                    continue;
                };
                let vj_p = Vertex::wrap(vj.clone())
                    .map(|v| BRepTool::vertex_point(&v))
                    .expect("FUN_AnalyzemapVon1E: vertex");
                let same = vi_p.is_equal_with_tolerance(&vj_p, tol);
                if same {
                    map_vv.add(&vi, vj.clone());
                    map_vv.add(&vj, vi);
                    break;
                }            }
        }
        let n_vv = map_vv.extent();
        if n_vv == n_v {
            res = GCLOSEDW;
        } else {
            res = UNCLOSEDW;
        }
    }

    res
}

/// `TopExp::MapShapesAndAncestors(S, TS, TA, M)` (`TopExp.cxx:80-120`) with the
/// index order OCCT's `NCollection_IndexedDataMap` gives the keys (1..Extent in
/// insertion order). [`crate::topo_tools::map_shapes_and_ancestors`] is the
/// same pass over a `HashMap`, which loses that order; the
/// `FindKey(i)` walks of `FUN_DetectVerticesOn1Edge` (`cxx:132-174`) and
/// `DetectPseudoInternalEdge` (`cxx:420-446`) iterate by index, so the ordered
/// form is what they need.
pub(crate) fn map_shapes_and_ancestors_ordered(
    shape: &TopoShape,
    ts: ShapeType,
    ta: ShapeType,
) -> IndexedShapeData<Vec<TopoShape>> {
    let mut map: IndexedShapeData<Vec<TopoShape>> = IndexedShapeData::new();
    for anc in explorer_shapes(shape, ta) {
        for sub in explorer_shapes(&anc, ts) {
            let index = map.add(&sub, Vec::new());
            if let Some(list) = map.change_from_index(index) {
                list.push(anc.clone());
            }
        }
    }
    for sub in explorer_shapes_avoiding(shape, ts, ta) {
        map.add(&sub, Vec::new());
    }
    map
}
