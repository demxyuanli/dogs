//! `TopOpeBRepTool_FuseEdges` — fuse the chains of edges that meet at a
//! useless vertex.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepTool/TopOpeBRepTool_FuseEdges.cxx`
//! (OCCT 8.0.0, 1,075 lines) and `.hxx` (121 lines). Split by OCCT function:
//!
//! | OCCT | lines | file |
//! |---|---|---|
//! | constructor, `AvoidEdges`, `Edges`, `ResultEdges`, `Faces`, `NbVertices`, `Shape` | `cxx:66-206` | here |
//! | `BuildListEdges` | `cxx:214-267` | here |
//! | `BuildListConnexEdge` | `cxx:480-537` | [`connex`] |
//! | `NextConnexEdge` | `cxx:549-618` | [`connex`] |
//! | `SameSupport` | `cxx:624-880` | [`same_support`] |
//! | `BuildListResultEdges` | `cxx:272-380` | [`result_edges`] |
//! | `UpdatePCurve` | `cxx:885-1038` | [`update_pcurve`] |
//! | `Perform` | `cxx:384-468` | [`perform`] |
//!
//! A "useless" vertex is one with exactly two connex edges that (a) share
//! exactly the same faces (`cxx:576-600`), (b) lie on the same geometric
//! support ([`same_support::same_support`]), and (c) are `FORWARD` or
//! `REVERSED` (`cxx:248-249`), so the chain can be collapsed into a single edge
//! built on the support between the two outer vertices.
//!
//! The three maps that drive it are `TopExp::MapShapesAndUniqueAncestors`
//! (`VERTEX -> EDGE`) and `TopExp::MapShapesAndAncestors` (`EDGE -> FACE`);
//! both are `NCollection_IndexedDataMap`, so **insertion order is part of the
//! behaviour**: `BuildListEdges` walks `myMapEdgLstFac.FindKey(1..Extent)`
//! (`cxx:242-262`) and that order decides the keys of `myMapLstEdg`, which
//! `Builder::End`'s `UpdateSplitAndMerged` then reads. [`IndexedAncestors`]
//! keeps that order.
//!
//! UNPORTED: `GeomLib::ExtendCurveToPoint` (`GeomLib.cxx:1268-1411`), the
//! fallback of `BuildListResultEdges` (`cxx:334-363`) for a `MakeEdge` that
//! failed because a new vertex lies outside a bounded support curve. OCCT
//! extends the curve by a C1 Hermite patch and retries; here the fusion is
//! reported as [`FuseEdgesError::FusionFailed`], which is the
//! `Standard_ConstructionError` OCCT throws at `cxx:358` / `cxx:362`.

pub mod connex;
pub mod perform;
pub mod result_edges;
pub mod same_support;
pub mod update_pcurve;

use std::collections::BTreeMap;

use crate::abs::{Orientation, ShapeType};
use crate::shape::{Edge, TopoShape};
use crate::topo_tools::{
    explorer_shapes, explorer_shapes_avoiding, shape_map_key, ShapeMapKey,
};

/// `FuseEdges : Fusion failed` (`cxx:358`, `cxx:362`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FuseEdgesError {
    /// `Standard_ConstructionError("FuseEdges : Fusion failed")` — the new
    /// edge could not be built on the support curve (`cxx:352-363`), including
    /// the UNPORTED `ExtendCurveToPoint` retry.
    FusionFailed,
}

/// `NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>,
/// TopTools_ShapeMapHasher>` (`hxx:112-117`): keys in insertion order, each
/// with an appended list.
#[derive(Debug, Default, Clone)]
pub struct IndexedAncestors {
    keys: Vec<TopoShape>,
    index: std::collections::HashMap<ShapeMapKey, usize>,
    lists: Vec<Vec<TopoShape>>,
}

impl IndexedAncestors {
    pub fn new() -> Self {
        Self::default()
    }

    /// `M.FindIndex(S)` followed by `M.Add(S, empty)` when absent
    /// (`TopExp.cxx:98-103`).
    fn index_of(&mut self, s: &TopoShape) -> usize {
        let key = shape_map_key(s);
        if let Some(&i) = self.index.get(&key) {
            return i;
        }
        self.keys.push(s.clone());
        self.lists.push(Vec::new());
        let i = self.keys.len() - 1;
        self.index.insert(key, i);
        i
    }

    /// `M.FindFromKey(S)` — the list bound to `S`, `None` when absent (OCCT
    /// raises `Standard_NoSuchObject` there).
    pub fn find_from_key(&self, s: &TopoShape) -> Option<&[TopoShape]> {
        let i = *self.index.get(&shape_map_key(s))?;
        Some(&self.lists[i])
    }

    /// `M.FindKey(i)` (1-based, `cxx:245`).
    pub fn find_key(&self, i: usize) -> Option<&TopoShape> {
        self.keys.get(i.checked_sub(1)?)
    }

    /// `M.Extent()`.
    pub fn extent(&self) -> usize {
        self.keys.len()
    }

    /// `TopExp::MapShapesAndAncestors` (`TopExp.cxx:80-120`).
    pub fn map_shapes_and_ancestors(shape: &TopoShape, ts: ShapeType, ta: ShapeType) -> Self {
        let mut m = Self::new();
        for anc in explorer_shapes(shape, ta) {
            for sub in explorer_shapes(&anc, ts) {
                let i = m.index_of(&sub);
                m.lists[i].push(anc.clone());
            }
        }
        // `TopExp.cxx:110-119`: sub-shapes under no ancestor.
        for sub in explorer_shapes_avoiding(shape, ts, ta) {
            m.index_of(&sub);
        }
        m
    }

    /// `TopExp::MapShapesAndUniqueAncestors` (`TopExp.cxx:124-177`) with
    /// `useOrientation == false`: `IsSame` decides (`cxx:156`).
    pub fn map_shapes_and_unique_ancestors(shape: &TopoShape, ts: ShapeType, ta: ShapeType) -> Self {
        let mut m = Self::new();
        for anc in explorer_shapes(shape, ta) {
            for sub in explorer_shapes(&anc, ts) {
                let i = m.index_of(&sub);
                let already = m.lists[i].iter().any(|x| crate::topo_tools::is_same(x, &anc));
                if !already {
                    m.lists[i].push(anc.clone());
                }
            }
        }
        for sub in explorer_shapes_avoiding(shape, ts, ta) {
            m.index_of(&sub);
        }
        m
    }
}

/// `class TopOpeBRepTool_FuseEdges` (`hxx:43`).
pub struct FuseEdges {
    my_shape: TopoShape,
    my_shape_done: bool,
    my_edges_done: bool,
    my_result_edges_done: bool,
    /// `myMapVerLstEdg` (`hxx:112`): `VERTEX -> EDGE`, unique ancestors.
    my_map_ver_lst_edg: IndexedAncestors,
    /// `myMapEdgLstFac` (`hxx:115`): `EDGE -> FACE`, non-unique ancestors.
    my_map_edg_lst_fac: IndexedAncestors,
    /// `myMapLstEdg` (`hxx:117`).
    my_map_lst_edg: BTreeMap<i32, Vec<TopoShape>>,
    /// `myMapEdg` (`hxx:118`).
    my_map_edg: BTreeMap<i32, TopoShape>,
    /// `myMapFaces` (`hxx:119`).
    my_map_faces: Vec<(TopoShape, TopoShape)>,
    /// `myNbConnexEdge` (`hxx:120`).
    my_nb_connex_edge: i32,
    /// `myAvoidEdg` (`hxx:121`).
    my_avoid_edg: Vec<ShapeMapKey>,
}

impl FuseEdges {
    /// Constructor (`cxx:66-81`). `Standard_NullObject` on a null shape is
    /// unreachable in the port (no null `TopoShape`).
    pub fn new(shape: &TopoShape) -> Self {
        Self {
            my_shape: shape.clone(),
            my_shape_done: false,
            my_edges_done: false,
            my_result_edges_done: false,
            my_map_ver_lst_edg: IndexedAncestors::new(),
            my_map_edg_lst_fac: IndexedAncestors::new(),
            my_map_lst_edg: BTreeMap::new(),
            my_map_edg: BTreeMap::new(),
            my_map_faces: Vec::new(),
            my_nb_connex_edge: 0,
            my_avoid_edg: Vec::new(),
        }
    }

    /// `AvoidEdges(theMapEdg)` (`cxx:88-91`): an `IndexedMap`, read only for
    /// membership.
    pub fn avoid_edges(&mut self, edges: impl IntoIterator<Item = TopoShape>) {
        self.my_avoid_edg = edges.into_iter().map(|e| shape_map_key(&e)).collect();
    }

    fn is_avoided(&self, e: &TopoShape) -> bool {
        let key = shape_map_key(e);
        self.my_avoid_edg.contains(&key)
    }

    /// `Edges(theMapLstEdg)` (`cxx:102-111`).
    pub fn edges(&mut self) -> Vec<(i32, Vec<TopoShape>)> {
        if !self.my_edges_done {
            self.build_list_edges();
        }
        self.my_map_lst_edg.iter().map(|(k, v)| (*k, v.clone())).collect()
    }

    /// `ResultEdges(theMapEdg)` (`cxx:115-128`).
    pub fn result_edges(&mut self) -> Result<BTreeMap<i32, TopoShape>, FuseEdgesError> {
        if !self.my_edges_done {
            self.build_list_edges();
        }
        if !self.my_result_edges_done {
            self.build_list_result_edges()?;
        }
        Ok(self.my_map_edg.clone())
    }

    /// `Faces(theMapFac)` (`cxx:136-157`): the modified faces, original ->
    /// rebuilt. Runs `Perform` when it has not run yet.
    pub fn faces(&mut self) -> Result<Vec<(TopoShape, TopoShape)>, FuseEdgesError> {
        if !self.my_edges_done {
            self.build_list_edges();
        }
        if !self.my_result_edges_done {
            self.build_list_result_edges()?;
        }
        if !self.my_shape_done {
            self.perform()?;
        }
        Ok(self.my_map_faces.clone())
    }

    /// `NbVertices()` (`cxx:158-184`): candidate vertices, one less than the
    /// size of every chain.
    pub fn nb_vertices(&mut self) -> i32 {
        if !self.my_edges_done {
            self.build_list_edges();
        }
        self.my_map_lst_edg.values().map(|l| l.len() as i32 - 1).sum()
    }

    /// `Shape()` (`cxx:187-206`).
    pub fn shape(&mut self) -> Result<TopoShape, FuseEdgesError> {
        if !self.my_edges_done {
            self.build_list_edges();
        }
        if !self.my_result_edges_done {
            self.build_list_result_edges()?;
        }
        if !self.my_shape_done {
            self.perform()?;
        }
        Ok(self.my_shape.clone())
    }

    /// `BuildListEdges()` (`cxx:214-267`).
    fn build_list_edges(&mut self) {
        self.my_map_lst_edg.clear();
        self.my_map_ver_lst_edg =
            IndexedAncestors::map_shapes_and_unique_ancestors(&self.my_shape, ShapeType::Vertex, ShapeType::Edge);
        self.my_map_edg_lst_fac =
            IndexedAncestors::map_shapes_and_ancestors(&self.my_shape, ShapeType::Edge, ShapeType::Face);

        // `cxx:239`: the `mapUniqEdg` of already-treated edges.
        let mut map_uniq_edg: Vec<ShapeMapKey> = Vec::new();

        // `cxx:242-262`: `for (iEdg = 1; iEdg <= myMapEdgLstFac.Extent(); iEdg++)`.
        for i_edg in 1..=self.my_map_edg_lst_fac.extent() {
            let Some(edgecur) = self.my_map_edg_lst_fac.find_key(i_edg).cloned() else {
                continue;
            };
            // `cxx:248-249`: only a FORWARD / REVERSED occurrence takes part.
            if map_uniq_edg.contains(&shape_map_key(&edgecur))
                || !(edgecur.orientation() == Orientation::Forward
                    || edgecur.orientation() == Orientation::Reversed)
            {
                continue;
            }
            if self.is_avoided(&edgecur) {
                continue;
            }
            let mut lst_edg: Vec<TopoShape> = Vec::new();
            self.build_list_connex_edge(&edgecur, &mut map_uniq_edg, &mut lst_edg);
            if lst_edg.len() > 1 {
                self.my_nb_connex_edge += 1;
                self.my_map_lst_edg.insert(self.my_nb_connex_edge, lst_edg);
            }
        }

        self.my_edges_done = true;
        self.my_result_edges_done = false;
    }


}

/// `TopExp::FirstVertex(E, CumOri = true)` (`TopExp.cxx:182-194`) /
/// `LastVertex(E, true)` (`cxx:198-209`).
///
/// `TopExp.cxx` iterates `TopoDS_Iterator(E, true)`, so the edge's own
/// orientation is composed into the children before the `FORWARD` (first) or
/// `REVERSED` (last) vertex is picked; the port picks among the stored
/// children and swaps for a `REVERSED` edge, which selects the same occurrence.
/// The returned vertex keeps the stored orientation, not the composed one.
pub(super) fn cumulated_vertex(e: &TopoShape, last: bool) -> Option<TopoShape> {
    let (v1, v2) = crate::topo_tools::edge_vertices(&Edge(e.clone()));
    let (v1, v2) = (v1.map(|v| v.0), v2.map(|v| v.0));
    if e.orientation() == Orientation::Reversed {
        if last { v1 } else { v2 }
    } else if last {
        v2
    } else {
        v1
    }
}
