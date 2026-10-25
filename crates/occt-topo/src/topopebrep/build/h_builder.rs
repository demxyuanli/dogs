//! `TopOpeBRepBuild_HBuilder` — the public face of the classic builder.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_HBuilder.hxx`
//! and `.cxx:39-738` (OCCT 8.0.0).
//!
//! `HBuilder` owns a [`Builder1`] — `TopOpeBRepBuild_Builder1`
//! (`HBuilder.hxx:170`: `TopOpeBRepBuild_Builder1 myBuilder`), the subclass
//! that adds the per-shape states of a SOLID/SOLID operation — and exposes it
//! through the API the `TopOpeBRep` package calls: the `Perform` / `Merge*`
//! drivers, the `IsSplit` / `Splits` / `IsMerged` / `Merged` queries, a
//! section cursor, and the four ancestor maps that answer "which DS edge /
//! curve / point produced this result sub-shape".
//!
//! UNPORTED (same file): nothing. `MergeKPart(TB1, TB2)` (`cxx:725-731`)
//! reaches `Builder1::MergeKPart(TB1, TB2)` (`Builder1.cxx:144-147`), which
//! sets the two states and runs the five-case `MergeKPart()` merge.
//!
//! The section cursor is kept as members rather than as the `thread_local`
//! `PLE` / `PITLE` statics of `cxx:197-198`: OCCT effectively keeps one cursor
//! per thread, the port one per `HBuilder`, which is the same thing for the
//! single-builder-per-thread use the statics already force.

use std::collections::HashMap;

use crate::abs::State;
use crate::shape::TopoShape;
use crate::topo_tools::{shape_map_key, ShapeMapKey};
use crate::topopebrep::ds::curve_explorer::CurveExplorer;
use crate::topopebrep::ds::export::fds_data;
use crate::topopebrep::ds::h_data_structure::HDataStructure;
use crate::topopebrep::ds::kind::Kind;

use super::builder::{Builder, Builder1};
use crate::topopebrep::ds::build_tool::BuildTool;

/// `class TopOpeBRepBuild_HBuilder` (`hxx:30-...`).
pub struct HBuilder {
    my_builder: Builder1,
    my_make_edge_ancestor_is_done: bool,
    my_make_curve_ancestor_is_done: bool,
    my_make_point_ancestor_is_done: bool,
    /// `mySectEdgeDSEdges1` / `:2` (`cxx:325`): section edge -> DS edge index,
    /// per ancestor rank.
    my_sect_edge_ds_edges1: HashMap<ShapeMapKey, i32>,
    my_sect_edge_ds_edges2: HashMap<ShapeMapKey, i32>,
    /// `myDSEdgesDSFaces1` / `:2` (`cxx:399`): DS edge index -> DS face
    /// indexes, per ancestor rank.
    my_ds_edges_ds_faces1: HashMap<i32, Vec<i32>>,
    my_ds_edges_ds_faces2: HashMap<i32, Vec<i32>>,
    /// `mySectEdgeDSCurve` (`cxx:503`): new edge -> DS curve index.
    my_sect_edge_ds_curve: HashMap<ShapeMapKey, i32>,
    /// `myNewVertexDSPoint` (`cxx:566`): new vertex -> DS point index.
    my_new_vertex_ds_point: HashMap<ShapeMapKey, i32>,
    /// `myEmptyIntegerList` (`cxx:476`): the list
    /// `GetDSFaceFromDSEdge` returns when nothing is bound.
    my_empty_integer_list: Vec<i32>,
    /// `PLE` / `PITLE` (`cxx:197-198`): the list `InitSection` filled and the
    /// read position in it.
    my_section_list: Vec<TopoShape>,
    my_section_index: usize,
    /// The `PLE == nullptr` state of `cxx:197-198`: `MoreSection` /
    /// `NextSection` are no-ops and `CurrentSection` raises until
    /// `InitSection` has run.
    my_section_init: bool,
}

impl HBuilder {
    /// `TopOpeBRepBuild_HBuilder(BT)` (`cxx:41-47`).
    pub fn new(bt: BuildTool) -> Self {
        Self {
            my_builder: Builder1::new(bt),
            my_make_edge_ancestor_is_done: false,
            my_make_curve_ancestor_is_done: false,
            my_make_point_ancestor_is_done: false,
            my_sect_edge_ds_edges1: HashMap::new(),
            my_sect_edge_ds_edges2: HashMap::new(),
            my_ds_edges_ds_faces1: HashMap::new(),
            my_ds_edges_ds_faces2: HashMap::new(),
            my_sect_edge_ds_curve: HashMap::new(),
            my_new_vertex_ds_point: HashMap::new(),
            my_empty_integer_list: Vec::new(),
            my_section_list: Vec::new(),
            my_section_index: 0,
            my_section_init: false,
        }
    }

    /// `Perform(HDS)` (`cxx:50-53`).
    pub fn perform(&mut self, hds: HDataStructure) {
        self.my_builder.perform(hds);
    }

    /// `Perform(HDS, S1, S2)` (`cxx:57-62`).
    pub fn perform_with_shapes(
        &mut self,
        hds: HDataStructure,
        s1: Option<&TopoShape>,
        s2: Option<&TopoShape>,
    ) {
        self.my_builder.perform_with_shapes(hds, s1, s2);
    }

    /// `Clear()` (`cxx:67-71`): `myBuilder.Clear()` then
    /// `InitExtendedSectionDS()` with the default `k = 3` (`hxx`).
    pub fn clear(&mut self) {
        self.my_builder.clear();
        self.init_extended_section_ds(3);
    }

    /// `DataStructure()` (`cxx:74-77`).
    pub fn data_structure(&self) -> Option<&HDataStructure> {
        self.my_builder.data_structure()
    }

    /// `ChangeBuildTool()` (`cxx:82-86`).
    pub fn change_build_tool(&mut self) -> &mut BuildTool {
        self.my_builder.change_build_tool()
    }

    /// `BuildTool()` (`cxx:89-92`).
    pub fn build_tool(&self) -> &BuildTool {
        self.my_builder.build_tool()
    }

    /// `MergeShapes(S1, ToBuild1, S2, ToBuild2)` (`cxx:97-102`).
    pub fn merge_shapes(
        &mut self,
        s1: Option<&TopoShape>,
        to_build1: State,
        s2: Option<&TopoShape>,
        to_build2: State,
    ) {
        self.my_builder
            .merge_shapes(s1, to_build1, s2, to_build2);
    }

    /// `MergeSolids(S1, ToBuild1, S2, ToBuild2)` (`cxx:107-113`).
    pub fn merge_solids(&mut self, s1: &TopoShape, to_build1: State, s2: &TopoShape, to_build2: State) {
        self.my_builder.merge_solids(s1, to_build1, s2, to_build2);
    }

    /// `MergeSolid(S, ToBuild)` (`cxx:116-119`).
    pub fn merge_solid(&mut self, s: &TopoShape, to_build: State) {
        self.my_builder.merge_solid(s, to_build);
    }

    /// `IsSplit(S, ToBuild)` (`cxx:124-128`).
    pub fn is_split(&self, s: &TopoShape, to_build: State) -> bool {
        self.my_builder.is_split(s, to_build)
    }

    /// `Splits(S, ToBuild)` (`cxx:132-138`).
    pub fn splits(&self, s: &TopoShape, to_build: State) -> &[TopoShape] {
        self.my_builder.splits(s, to_build)
    }

    /// `IsMerged(S, ToBuild)` (`cxx:142-146`).
    pub fn is_merged(&self, s: &TopoShape, to_build: State) -> bool {
        self.my_builder.is_merged(s, to_build)
    }

    /// `Merged(S, ToBuild)` (`cxx:149-155`).
    pub fn merged(&self, s: &TopoShape, to_build: State) -> &[TopoShape] {
        self.my_builder.merged(s, to_build)
    }

    /// `NewVertex(I)` (`cxx:159-163`).
    pub fn new_vertex(&self, i: i32) -> Option<&TopoShape> {
        self.my_builder.new_vertex(i)
    }

    /// `NewEdges(I)` (`cxx:167-172`).
    pub fn new_edges(&self, i: i32) -> &[TopoShape] {
        self.my_builder.new_edges(i)
    }

    /// `ChangeNewEdges(I)` (`cxx:176-180`).
    pub fn change_new_edges(&mut self, i: i32) -> &mut Vec<TopoShape> {
        self.my_builder.change_new_edges(i)
    }

    /// `NewFaces(I)` (`cxx:184-188`).
    pub fn new_faces(&self, i: i32) -> Option<&[TopoShape]> {
        self.my_builder.new_faces(i)
    }

    /// `Section()` (`cxx:191-195`).
    pub fn section(&mut self) -> &[TopoShape] {
        self.my_builder.section()
    }

    /// `InitExtendedSectionDS(k)` (`cxx:203-222`).
    pub fn init_extended_section_ds(&mut self, k: i32) {
        match k {
            1 => self.my_make_curve_ancestor_is_done = false,
            2 => self.my_make_edge_ancestor_is_done = false,
            3 => {
                self.my_make_edge_ancestor_is_done = false;
                self.my_make_curve_ancestor_is_done = false;
            }
            _ => {}
        }
    }

    /// `InitSection(k)` (`cxx:226-255`).
    pub fn init_section(&mut self, k: i32) {
        self.my_section_list.clear();
        self.my_section_index = 0;
        self.my_section_init = true;
        self.init_extended_section_ds(k);
        match k {
            1 => self.my_builder.section_curves_into(&mut self.my_section_list),
            2 => self.my_builder.section_edges_into(&mut self.my_section_list),
            3 => {
                let l = std::mem::take(&mut self.my_section_list);
                let mut l = l;
                self.my_builder.section_list(&mut l);
                self.my_section_list = l;
            }
            _ => return,
        }
    }

    /// `MoreSection()` (`cxx:259-266`).
    pub fn more_section(&self) -> bool {
        self.my_section_init && self.my_section_index < self.my_section_list.len()
    }

    /// `NextSection()` (`cxx:272-283`).
    pub fn next_section(&mut self) {
        if !self.my_section_init {
            return;
        }
        if self.my_section_index < self.my_section_list.len() {
            self.my_section_index += 1;
        }
    }

    /// `CurrentSection()` (`cxx:287-296`): panics on the two
    /// `Standard_ProgramError("no more CurrentSection")` raises.
    pub fn current_section(&self) -> &TopoShape {
        if !self.my_section_init {
            panic!("no more CurrentSection");
        }
        if self.my_section_index >= self.my_section_list.len() {
            panic!("no more CurrentSection");
        }
        &self.my_section_list[self.my_section_index]
    }

    /// `MakeEdgeAncestorMap()` (`cxx:301-416`).
    pub fn make_edge_ancestor_map(&mut self) {
        if self.my_make_edge_ancestor_is_done {
            return;
        }
        self.my_sect_edge_ds_edges1.clear();
        self.my_sect_edge_ds_edges2.clear();
        self.my_ds_edges_ds_faces1.clear();
        self.my_ds_edges_ds_faces2.clear();
        self.my_make_edge_ancestor_is_done = true;

        let mut sect1: Vec<(ShapeMapKey, i32)> = Vec::new();
        let mut sect2: Vec<(ShapeMapKey, i32)> = Vec::new();
        let mut faces1: Vec<(i32, i32)> = Vec::new();
        let mut faces2: Vec<(i32, i32)> = Vec::new();

        if let Some(hds) = self.my_builder.data_structure() {
            let ds = hds.ds();

            // `cxx:319-357`: the split(ON) parts of the rank-1 / rank-2
            // section edges name the DS edge that produced them.
            for (key, losos) in self.my_builder.m_split(State::On).iter() {
                let Some(sha_spl) = self.my_builder.shape_of_key(key) else {
                    continue;
                };
                let ei = ds.shape_index(&sha_spl, true);
                let re = ds.ancestor_rank_of(&sha_spl);
                if re == 0 {
                    continue;
                }
                // `if (!mySectEdgeDSEdges1.IsBound(SecEdg)) Bind(...)`: first
                // binding wins.
                for sec_edg in losos.list_on_state() {
                    let k = shape_map_key(sec_edg);
                    if re == 1 {
                        if !sect1.iter().any(|(ek, _)| *ek == k) {
                            sect1.push((k, ei));
                        }
                    } else if re == 2 && !sect2.iter().any(|(ek, _)| *ek == k) {
                        sect2.push((k, ei));
                    }
                }
            }

            // `cxx:361-414`
            let ns = ds.nb_shapes();
            for fi in 1..=ns {
                let Some(fds) = ds.shape(fi, false) else {
                    continue;
                };
                if fds.shape_type() != crate::abs::ShapeType::Face {
                    continue;
                }
                let rf = ds.ancestor_rank_of(&fds);
                if rf == 0 {
                    continue;
                }
                for ssi in ds.shape_interferences(&fds, true) {
                    let i = ssi.read().unwrap();
                    // `down_cast<TopOpeBRepDS_ShapeShapeInterference>`:
                    // `cxx:377-382` skips a null cast.
                    if !i.is_shape_shape() {
                        continue;
                    }
                    let (gk, gi, _st, _si) = fds_data(&i);
                    if gk != Kind::Edge {
                        continue;
                    }
                    if rf == 1 {
                        faces1.push((gi, fi));
                    } else if rf == 2 {
                        faces2.push((gi, fi));
                    }
                }
            }
        }

        for (k, ei) in sect1 {
            self.my_sect_edge_ds_edges1.entry(k).or_insert(ei);
        }
        for (k, ei) in sect2 {
            self.my_sect_edge_ds_edges2.entry(k).or_insert(ei);
        }
        for (gi, fi) in faces1 {
            self.my_ds_edges_ds_faces1.entry(gi).or_default().push(fi);
        }
        for (gi, fi) in faces2 {
            self.my_ds_edges_ds_faces2.entry(gi).or_default().push(fi);
        }
    }

    /// `GetDSEdgeFromSectEdge(E, rank)` (`cxx:420-446`).
    pub fn get_ds_edge_from_sect_edge(&mut self, e: &TopoShape, rank: i32) -> i32 {
        if !self.my_make_edge_ancestor_is_done {
            self.make_edge_ancestor_map();
        }
        let key = shape_map_key(e);
        let map = if rank == 1 {
            &self.my_sect_edge_ds_edges1
        } else if rank == 2 {
            &self.my_sect_edge_ds_edges2
        } else {
            return 0;
        };
        map.get(&key).copied().unwrap_or(0)
    }

    /// `GetDSFaceFromDSEdge(indexEdg, rank)` (`cxx:449-477`).
    pub fn get_ds_face_from_ds_edge(&mut self, index_edg: i32, rank: i32) -> &[i32] {
        if !self.my_make_edge_ancestor_is_done {
            self.make_edge_ancestor_map();
        }
        let map = if rank == 1 {
            &self.my_ds_edges_ds_faces1
        } else if rank == 2 {
            &self.my_ds_edges_ds_faces2
        } else {
            return &self.my_empty_integer_list;
        };
        match map.get(&index_edg) {
            Some(loi) => loi,
            None => &self.my_empty_integer_list,
        }
    }

    /// `MakeCurveAncestorMap()` (`cxx:480-504`).
    pub fn make_curve_ancestor_map(&mut self) {
        if self.my_make_curve_ancestor_is_done {
            return;
        }
        self.my_sect_edge_ds_curve.clear();
        self.my_make_curve_ancestor_is_done = true;

        let mut binds: Vec<(ShapeMapKey, i32)> = Vec::new();
        if let Some(hds) = self.my_builder.data_structure() {
            // `TopOpeBRepDS_CurveExplorer cex(DS, true)` (`cxx:490`).
            let mut cex = CurveExplorer::with_ds(hds.ds(), true);
            while cex.more() {
                let ic = cex.index();
                for e in self.my_builder.new_edges(ic) {
                    binds.push((shape_map_key(e), ic));
                }
                cex.next();
            }
        }
        for (k, ic) in binds {
            self.my_sect_edge_ds_curve.insert(k, ic);
        }
    }

    /// `GetDSCurveFromSectEdge(SectEdge)` (`cxx:508-522`).
    pub fn get_ds_curve_from_sect_edge(&mut self, sect_edge: &TopoShape) -> i32 {
        if !self.my_make_curve_ancestor_is_done {
            self.make_curve_ancestor_map();
        }
        self.my_sect_edge_ds_curve
            .get(&shape_map_key(sect_edge))
            .copied()
            .unwrap_or(0)
    }

    /// `GetDSFaceFromDSCurve(indexCur, rank)` (`cxx:526-549`).
    pub fn get_ds_face_from_ds_curve(&mut self, index_cur: i32, rank: i32) -> i32 {
        if !self.my_make_curve_ancestor_is_done {
            self.make_curve_ancestor_map();
        }
        let Some(hds) = self.my_builder.data_structure() else {
            return 0;
        };
        let ds = hds.ds();
        let dsc = ds.curve(index_cur);
        let shape = if rank == 1 {
            dsc.shape1()
        } else if rank == 2 {
            dsc.shape2()
        } else {
            None
        };
        match shape {
            Some(f) => ds.shape_index(f, true),
            None => 0,
        }
    }

    /// `GetDSPointFromNewVertex(NewVert)` (`cxx:553-577`).
    pub fn get_ds_point_from_new_vertex(&mut self, new_vert: &TopoShape) -> i32 {
        if !self.my_make_point_ancestor_is_done {
            self.my_make_point_ancestor_is_done = true;
            let nb_point = match self.my_builder.data_structure() {
                Some(hds) => hds.ds().nb_points(),
                None => 0,
            };
            let mut binds: Vec<(ShapeMapKey, i32)> = Vec::new();
            for i in 1..=nb_point {
                if let Some(vertex) = self.my_builder.new_vertex(i) {
                    if !vertex.is_null() {
                        binds.push((shape_map_key(vertex), i));
                    }
                }
            }
            for (k, i) in binds {
                self.my_new_vertex_ds_point.insert(k, i);
            }
        }
        self.my_new_vertex_ds_point
            .get(&shape_map_key(new_vert))
            .copied()
            .unwrap_or(0)
    }

    /// `EdgeCurveAncestors(E, F1, F2, IC)` (`cxx:580-607`).
    pub fn edge_curve_ancestors(
        &mut self,
        e: &TopoShape,
        f1: &mut Option<TopoShape>,
        f2: &mut Option<TopoShape>,
        ic: &mut i32,
    ) -> bool {
        *f1 = None;
        *f2 = None;
        *ic = 0;

        *ic = self.get_ds_curve_from_sect_edge(e);
        if *ic == 0 {
            return false;
        }

        let if1 = self.get_ds_face_from_ds_curve(*ic, 1);
        let if2 = self.get_ds_face_from_ds_curve(*ic, 2);

        let Some(hds) = self.my_builder.data_structure() else {
            return true;
        };
        let ds = hds.ds();
        *f1 = ds.shape(if1, true);
        *f2 = ds.shape(if2, true);
        true
    }

    /// `EdgeSectionAncestors(E, LF1, LF2, LE1, LE2)` (`cxx:612-714`).
    pub fn edge_section_ancestors(
        &mut self,
        e: &TopoShape,
        lf1: &mut Vec<TopoShape>,
        lf2: &mut Vec<TopoShape>,
        le1: &mut Vec<TopoShape>,
        le2: &mut Vec<TopoShape>,
    ) -> bool {
        if e.shape_type() != crate::abs::ShapeType::Edge {
            return false;
        }

        lf1.clear();
        lf2.clear();
        le1.clear();
        le2.clear();

        let ie1 = self.get_ds_edge_from_sect_edge(e, 1);
        let ie2 = self.get_ds_edge_from_sect_edge(e, 2);

        // `cxx:636-687`
        let (mut f1, mut f2): (Vec<i32>, Vec<i32>) = (Vec::new(), Vec::new());
        if ie1 != 0 && ie2 != 0 {
            let mut moi: Vec<i32> = Vec::new();
            f1 = self.get_ds_face_from_ds_edge(ie1, 1).to_vec();
            moi.clear();
            for v in &f1 {
                moi.push(*v);
            }
            for v in self.get_ds_face_from_ds_edge(ie2, 1).to_vec() {
                if !moi.contains(&v) {
                    moi.push(v);
                    f1.push(v);
                }
            }
            f2 = self.get_ds_face_from_ds_edge(ie1, 2).to_vec();
            moi.clear();
            for v in &f2 {
                moi.push(*v);
            }
            for v in self.get_ds_face_from_ds_edge(ie2, 2).to_vec() {
                if !moi.contains(&v) {
                    moi.push(v);
                    f2.push(v);
                }
            }
        } else if ie1 != 0 {
            f1 = self.get_ds_face_from_ds_edge(ie1, 1).to_vec();
            f2 = self.get_ds_face_from_ds_edge(ie1, 2).to_vec();
        } else if ie2 != 0 {
            f1 = self.get_ds_face_from_ds_edge(ie2, 1).to_vec();
            f2 = self.get_ds_face_from_ds_edge(ie2, 2).to_vec();
        }

        let Some(hds) = self.my_builder.data_structure() else {
            return false;
        };
        let ds = hds.ds();

        // `cxx:689-696`
        if ie1 != 0 {
            if let Some(s) = ds.shape(ie1, true) {
                le1.push(s);
            }
        }
        if ie2 != 0 {
            if let Some(s) = ds.shape(ie2, true) {
                le2.push(s);
            }
        }

        // `cxx:698-707`
        for curr in f1 {
            if let Some(s) = ds.shape(curr, true) {
                lf1.push(s);
            }
        }
        for curr in f2 {
            if let Some(s) = ds.shape(curr, true) {
                lf2.push(s);
            }
        }

        // `cxx:709-711`
        let mut r = !lf1.is_empty() && !lf2.is_empty();
        r = r && (!le1.is_empty() || !le2.is_empty());
        r
    }

    /// `IsKPart()` (`cxx:718-721`).
    pub fn is_k_part(&self) -> i32 {
        self.my_builder.is_k_part()
    }

    /// `MergeKPart(TB1, TB2)` (`cxx:725-731`).
    pub fn merge_k_part(&mut self, tb1: State, tb2: State) {
        self.my_builder.merge_k_part_states(tb1, tb2);
    }

    /// `ChangeBuilder()` (`cxx:735-738`).
    pub fn change_builder(&mut self) -> &mut Builder {
        self.my_builder.base_mut()
    }

    /// The `TopOpeBRepBuild_Builder1` this builder owns.
    pub fn builder1(&self) -> &Builder1 {
        &self.my_builder
    }
}
