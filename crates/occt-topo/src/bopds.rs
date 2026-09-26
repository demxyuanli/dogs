//! BOPDS — the data structure hub of the Boolean Component (Phase 15, wave 1).
//!
//! Port of the OCCT `BOPDS_*` classes:
//!
//! | OCCT class                  | Rust port                 | Purpose                              |
//! |-----------------------------|---------------------------|--------------------------------------|
//! | `BOPDS_ShapeInfo`           | [`BopdsShapeInfo`]        | shape + type + sub-shape indices     |
//! | `BOPDS_Pave`                | [`BopdsPave`]             | a vertex on an edge (index+param)    |
//! | `BOPDS_PaveBlock`           | [`BopdsPaveBlock`]        | an interval on an edge + splitting   |
//! | `BOPDS_CommonBlock`         | [`BopdsCommonBlock`]      | edges sharing a geometric interval   |
//! | `BOPDS_FaceInfo`            | [`BopdsFaceInfo`]         | state of a face (IN/ON paves)        |
//! | `BOPDS_IndexRange`          | [`BopdsIndexRange`]       | shape parameter range (rank)         |
//! | `BOPDS_Interf`              | [`BopdsInterf`]           | an interference between two shapes   |
//! | `BOPDS_DS`                  | [`BopdsDS`]               | the DS: shape index, pave-block pool |
//! | `BOPDS_Iterator`            | [`BopdsIterator`]         | inter-argument interference pairs    |
//! | `BOPDS_SubIterator`         | [`BopdsSubIterator`]      | interference pairs of two sub-sets   |
//! | `BOPDS_IteratorSI`          | [`BopdsIteratorSI`]       | self-intersection candidates         |
//! | `BOPDS_Tools`               | [`bopds_tools`]           | type helpers (static functions)      |
//!
//! The DS stores every participating shape (arguments + all their sub-shapes)
//! once, keyed by `TShape` identity, and answers the queries the boolean
//! builder needs: shape→index round-trips, per-shape boundary sub-shape
//! indices, per-edge pave-block pools, and the interference-candidate
//! iterators.
//!
//! This module is **self-contained**: it depends only on already-ported
//! modules (`crate::topo_tools_full`, `crate::brep_tool`,
//! `crate::bbox_from_geometry`, `occt_core::kernel::containers`) and not on
//! the other Phase-15 wave-1 modules.

use std::collections::{HashMap, HashSet};


use occt_core::bnd::BndBox;
use occt_core::kernel::containers::IndexedMap;

use crate::abs::ShapeType;
use crate::bbox_from_geometry::shape_bbox;
use crate::bopds_ff::BopdsInterfFf;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, TopoShape};
use crate::topo_tools_full::shapes_of;
use crate::bopds_pave::ShapeKey;
use crate::bopds_tree::{prepared_sub_shapes, structural_children};

pub use crate::bopds_cb::{BopdsCommonBlock, BopdsFaceInfo, BopdsInterf};
pub use crate::bopds_iter::{intersect_ext_pairs, BopdsIterator, BopdsIteratorSI, BopdsSubIterator};
pub use crate::bopds_pave::{BopdsIndexRange, BopdsPave, BopdsPaveBlock, BopdsShapeInfo};
pub use crate::bopds_tools;

// ---------------------------------------------------------------------------
// BopdsDS
// ---------------------------------------------------------------------------

/// The data structure of the Boolean Component.
/// Source: `BOPDS_DS.hxx/.cxx`.
///
/// Contents:
/// 1. the arguments of an operation;
/// 2. one [`BopdsShapeInfo`] per argument and sub-shape, indexed by `TShape`
///    identity through an [`IndexedMap`];
/// 3. per-argument index ranges (ranks);
/// 4. the pave-block pool for source edges;
/// 5. the face-info pool;
/// 6. same-domain shape links and interference tracking.
#[derive(Debug, Clone)]
pub struct BopdsDS {
    /// Shape → index (identity keyed by TShape pointer).
    pub(crate) map: IndexedMap<ShapeKey>,
    /// Parallel array: per-shape information, index-aligned with `map`.
    pub(crate) shape_infos: Vec<BopdsShapeInfo>,
    /// Parallel array: cached bounding boxes, index-aligned with `map`.
    pub(crate) boxes: Vec<BndBox>,
    /// The arguments of the operation.
    pub(crate) arguments: Vec<TopoShape>,
    /// Per-argument index ranges.
    pub(crate) ranges: Vec<BopdsIndexRange>,
    /// Number of source shapes (arguments + sub-shapes present after `init`).
    pub(crate) nb_source_shapes: usize,
    /// Pool of pave-block lists, referenced by `BopdsShapeInfo::pb_reference`.
    pub(crate) pave_blocks_pool: Vec<Vec<BopdsPaveBlock>>,
    /// Pool of face states.
    pub(crate) face_info_pool: Vec<BopdsFaceInfo>,
    /// Common blocks collected by `update_common_block`.
    pub(crate) common_blocks: Vec<BopdsCommonBlock>,
    /// Same-domain shape map `index → same-domain index`.
    pub(crate) shapes_sd: HashMap<usize, usize>,
    /// Indices of the vertices whose tolerance was increased during the
    /// intersection (`BOPAlgo_PaveFiller::myIncreasedSS`). The repeat
    /// intersection stage re-runs V/V, V/E and V/F for these vertices.
    pub(crate) increased_ss: HashSet<usize>,
    /// Flat set of interfering pairs `(min, max)`. Mirror of OCCT's
    /// `myInterfTB` — membership queries only (`has_interf_pair`); the full
    /// per-type records (with the new-vertex index) live in `interf_vv..`.
    pub(crate) interferences: HashSet<(usize, usize)>,
    /// Set of interfered shape indices.
    pub(crate) interfered: HashSet<usize>,
    /// Typed interference records V/V. Source: `BOPDS_DS::InterfVV`.
    pub(crate) interf_vv: Vec<BopdsInterf>,
    /// Typed interference records V/E.
    pub(crate) interf_ve: Vec<BopdsInterf>,
    /// Typed interference records V/F.
    pub(crate) interf_vf: Vec<BopdsInterf>,
    /// Typed interference records E/E.
    pub(crate) interf_ee: Vec<BopdsInterf>,
    /// Typed interference records E/F.
    pub(crate) interf_ef: Vec<BopdsInterf>,
    /// Typed F/F interference records with section curves/points.
    /// Source: `BOPDS_DS::InterfFF`.
    pub(crate) interf_ff: Vec<BopdsInterfFf>,
}

impl BopdsDS {
    /// Empty constructor.
    pub fn new() -> Self {
        Self {
            map: IndexedMap::new(),
            shape_infos: Vec::new(),
            boxes: Vec::new(),
            arguments: Vec::new(),
            ranges: Vec::new(),
            nb_source_shapes: 0,
            pave_blocks_pool: Vec::new(),
            face_info_pool: Vec::new(),
            common_blocks: Vec::new(),
            shapes_sd: HashMap::new(),
            increased_ss: HashSet::new(),
            interferences: HashSet::new(),
            interfered: HashSet::new(),
            interf_vv: Vec::new(),
            interf_ve: Vec::new(),
            interf_vf: Vec::new(),
            interf_ee: Vec::new(),
            interf_ef: Vec::new(),
            interf_ff: Vec::new(),
        }
    }

    /// Clears the contents.
    pub fn clear(&mut self) {
        self.map = IndexedMap::new();
        self.shape_infos.clear();
        self.boxes.clear();
        self.arguments.clear();
        self.ranges.clear();
        self.nb_source_shapes = 0;
        self.pave_blocks_pool.clear();
        self.face_info_pool.clear();
        self.common_blocks.clear();
        self.shapes_sd.clear();
        self.increased_ss.clear();
        self.interferences.clear();
        self.interfered.clear();
        self.interf_vv.clear();
        self.interf_ve.clear();
        self.interf_vf.clear();
        self.interf_ee.clear();
        self.interf_ef.clear();
        self.interf_ff.clear();
    }

    /// Sets the arguments of the operation (they are appended by [`BopdsDS::init`]).
    pub fn set_arguments(&mut self, arguments: Vec<TopoShape>) {
        self.arguments = arguments;
    }

    /// Returns the arguments of the operation.
    pub fn arguments(&self) -> &[TopoShape] {
        &self.arguments
    }

    /// Initializes the data structure for the arguments: appends each argument
    /// and its sub-shapes, and builds the per-argument index ranges.
    pub fn init(&mut self, arguments: &[TopoShape]) {
        self.clear();
        self.arguments = arguments.to_vec();
        self.ranges.clear();
        let mut start = 0usize;
        for a in arguments {
            if self.map.contains(&ShapeKey(a.clone())) {
                continue;
            }
            self.append(a.clone()).expect("BopdsDS::init: append failed");
            let end = self.nb_shapes() - 1;
            self.ranges.push(BopdsIndexRange::new(start, end));
            start = end + 1;
        }
        self.nb_source_shapes = self.nb_shapes();
        self.mark_degenerated_edge_flags();
    }

    /// Returns the total number of shapes stored.
    pub fn nb_shapes(&self) -> usize {
        self.shape_infos.len()
    }

    /// Returns the number of source shapes stored (arguments + sub-shapes).
    pub fn nb_source_shapes(&self) -> usize {
        self.nb_source_shapes
    }

    /// Returns the number of index ranges.
    pub fn nb_ranges(&self) -> usize {
        self.ranges.len()
    }

    /// Returns the index range `i`.
    pub fn range(&self, i: usize) -> Option<BopdsIndexRange> {
        self.ranges.get(i).copied()
    }

    /// Returns the rank (argument index) of the shape of index `i`.
    ///
    /// Returns `0` when the shape is not covered by any range (a new shape).
    pub fn rank(&self, i: usize) -> usize {
        for (r, range) in self.ranges.iter().enumerate() {
            if range.contains(i) {
                return r;
            }
        }
        0
    }

    /// True if the shape of index `i` is not a source shape/sub-shape.
    pub fn is_new_shape(&self, i: usize) -> bool {
        i >= self.nb_source_shapes
    }

    /// Appends the shape and its whole sub-shape subtree to the data
    /// structure. Returns the index of the shape (existing if already present).
    pub fn append(&mut self, shape: TopoShape) -> Result<usize, String> {
        if let Some(i) = self.index(&shape) {
            return Ok(i);
        }
        let kind = shape.shape_type();
        let index = self.map.add(ShapeKey(shape.clone()));
        let bbox = shape_bbox(&shape);
        self.shape_infos
            .push(BopdsShapeInfo {
                shape: shape.clone(),
                kind,
                sub_indices: Vec::new(),
                pb_reference: -1,
                flag: None,
            });
        self.boxes.push(bbox);

        // Walk the structural children so the whole subtree is present.
        let children = structural_children(&shape);
        for c in children {
            self.append(c)?;
        }

        // Store the prepared boundary sub-shape indices.
        let prepared = prepared_sub_shapes(&shape);
        let sub: Vec<usize> = prepared
            .iter()
            .map(|s| self.index(s).expect("BopdsDS::append: sub-shape must be indexed"))
            .collect();
        self.shape_infos[index].sub_indices = sub;
        Ok(index)
    }

    /// `BOPDS_DS::prepareFaces`: degenerated edges store the parent-face index
    /// in `ShapeInfo::flag` (`SetFlag(aFaceIndex)`).
    fn mark_degenerated_edge_flags(&mut self) {
        let n = self.nb_shapes();
        for n_f in 0..n {
            if self.shape_infos.get(n_f).map(|s| s.kind) != Some(ShapeType::Face) {
                continue;
            }
            let subs = self.shape_infos[n_f].sub_indices.clone();
            for n_e in subs {
                if self.shape_infos.get(n_e).map(|s| s.kind) != Some(ShapeType::Edge) {
                    continue;
                }
                let degenerated = self.shape_infos[n_e]
                    .shape
                    .clone();
                if BRepTool::is_degenerated(&Edge(degenerated)) {
                    self.shape_infos[n_e].set_flag(n_f);
                }
            }
        }
    }

    /// Appends a pre-built shape info. Returns the index of the shape
    /// (existing if already present).
    pub fn append_info(&mut self, info: BopdsShapeInfo) -> usize {
        if let Some(i) = self.map.find_index(&ShapeKey(info.shape.clone())) {
            return i;
        }
        let index = self.map.add(ShapeKey(info.shape.clone()));
        let bbox = shape_bbox(&info.shape);
        self.shape_infos.push(info);
        self.boxes.push(bbox);
        index
    }

    /// Returns the information about the shape with index `i`.
    pub fn shape_info(&self, i: usize) -> Option<&BopdsShapeInfo> {
        self.shape_infos.get(i)
    }

    /// Mutable access to the information about the shape with index `i`.
    pub fn change_shape_info(&mut self, i: usize) -> Option<&mut BopdsShapeInfo> {
        self.shape_infos.get_mut(i)
    }

    /// Returns the shape with index `i`.
    pub fn shape(&self, i: usize) -> Option<&TopoShape> {
        self.shape_infos.get(i).map(|s| &s.shape)
    }

    /// Returns the index of the shape (by TShape identity).
    pub fn index(&self, shape: &TopoShape) -> Option<usize> {
        self.map.find_index(&ShapeKey(shape.clone()))
    }

    /// Returns the cached bounding box of the shape with index `i`.
    pub fn box_of(&self, i: usize) -> Option<&BndBox> {
        self.boxes.get(i)
    }

    /// Rebuilds the bounding box of the vertex `index` so it covers the point
    /// and the tolerance sphere.
    ///
    /// Port of the `BRepBndLib::Add` + `SetGap` part of
    /// `BOPAlgo_PaveFiller::UpdateVertex`: a vertex whose tolerance grows must
    /// have its DS box enlarged, or the repeat intersection (which re-selects
    /// pairs by box overlap) would miss its newly-touching neighbours.
    pub fn refresh_vertex_box(&mut self, index: usize, tol: f64) {
        let Some(shape) = self.shape(index) else { return };
        let mut b = shape_bbox(shape);
        b.enlarge(tol);
        if let Some(slot) = self.boxes.get_mut(index) {
            *slot = b;
        }
    }

    /// Rebuild the DS bounding box of any shape (`BRepBndLib::Add` + `SetGap`).
    pub fn refresh_shape_box(&mut self, index: usize, extra_gap: f64) {
        let Some(shape) = self.shape(index).cloned() else {
            return;
        };
        let mut b = shape_bbox(&shape);
        if extra_gap > 0.0 {
            b.enlarge(extra_gap);
        }
        if let Some(slot) = self.boxes.get_mut(index) {
            *slot = b;
        }
    }

    /// True if the shape with index `i` has pave-block information.
    pub fn has_pave_blocks(&self, i: usize) -> bool {
        self.shape_infos.get(i).map_or(false, |s| s.pb_reference >= 0)
    }

    /// Returns the pave blocks of the shape with index `i`.
    pub fn pave_blocks(&self, i: usize) -> &[BopdsPaveBlock] {
        match self.shape_infos.get(i).and_then(|s| {
            if s.pb_reference >= 0 {
                Some(s.pb_reference as usize)
            } else {
                None
            }
        }) {
            Some(r) => &self.pave_blocks_pool[r],
            None => &[],
        }
    }

    /// Returns the mutable list of pave blocks of the shape with index `i`,
    /// lazily allocating a pool slot on first access.
    pub fn change_pave_blocks_mut(&mut self, i: usize) -> &mut Vec<BopdsPaveBlock> {
        let r = self.ensure_pave_blocks(i);
        &mut self.pave_blocks_pool[r]
    }

    /// Returns the full pave-block pool.
    pub fn pave_blocks_pool(&self) -> &[Vec<BopdsPaveBlock>] {
        &self.pave_blocks_pool
    }

    /// Mutable pave-block pool (`BOPDS_DS::ChangePaveBlocksPool`).
    pub fn change_pave_blocks_pool(&mut self) -> &mut Vec<Vec<BopdsPaveBlock>> {
        &mut self.pave_blocks_pool
    }

    /// Initializes the default pave block of the edge with index `edge_index`
    /// from its boundary vertices (when the edge has finite parameters).
    pub fn init_pave_blocks_for_edge(&mut self, edge_index: usize) {
        if self.has_pave_blocks(edge_index) {
            return;
        }
        if self.shape_infos.get(edge_index).map(|s| s.kind) != Some(ShapeType::Edge) {
            return;
        }
        let shape = self.shape_infos[edge_index].shape.clone();
        let (first, last) = BRepTool::edge_parameters(&Edge(shape));
        if !first.is_finite() || !last.is_finite() {
            return;
        }
        let verts = self.shape_infos[edge_index].sub_indices.clone();
        if verts.len() < 2 {
            return;
        }
        let mut pb = BopdsPaveBlock::new();
        pb.edge_index = edge_index;
        pb.original_edge = edge_index;
        pb.index1 = verts[0];
        pb.index2 = verts[1];
        pb.first = first;
        pb.last = last;
        let r = self.ensure_pave_blocks(edge_index);
        self.pave_blocks_pool[r].push(pb);
    }

    fn ensure_pave_blocks(&mut self, i: usize) -> usize {
        if self.shape_infos[i].pb_reference >= 0 {
            return self.shape_infos[i].pb_reference as usize;
        }
        let r = self.pave_blocks_pool.len();
        self.pave_blocks_pool.push(Vec::new());
        self.shape_infos[i].pb_reference = r as i64;
        r
    }

    /// Updates the pave blocks for all shapes in the data structure: every
    /// block carrying extra paves is split into elementary blocks.
    pub fn update_pave_blocks(&mut self) {
        for list in &mut self.pave_blocks_pool {
            let mut new_list = Vec::new();
            for mut pb in list.drain(..) {
                if !pb.is_to_update() {
                    new_list.push(pb);
                } else {
                    pb.update(&mut new_list, true);
                }
            }
            *list = new_list;
        }
    }

    /// Replaces the bound vertex indices of every pave block in the pool by
    /// their same-domain (SD) representatives.
    ///
    /// Source: `BOPDS_DS::UpdatePaveBlocksWithSDVertices` +
    /// `BOPDS_DS::UpdatePaveBlockWithSDVertices`. The counterpart of
    /// [`BopdsDS::update_pave_blocks`] (which *splits* blocks carrying extra
    /// paves); this step only *redirects* the two bound paves of each block to
    /// the SD vertices created by the V/V, E/E, V/F and E/F stages.
    pub fn update_pave_blocks_with_sd_vertices(&mut self) {
        let sd = &self.shapes_sd;
        for list in &mut self.pave_blocks_pool {
            for pb in list.iter_mut() {
                let (n1, n2) = pb.indices();
                let mut c1 = n1;
                while let Some(&nx) = sd.get(&c1) {
                    c1 = nx;
                }
                let mut c2 = n2;
                while let Some(&nx) = sd.get(&c2) {
                    c2 = nx;
                }
                pb.set_indices(c1, c2);
            }
        }
    }

    /// Updates the pave blocks of the edges of the common block and records
    /// the block. Source: `BOPDS_DS::UpdateCommonBlock` (simplified).
    pub fn update_common_block(&mut self, cb: &BopdsCommonBlock) {
        if cb.indices.is_empty() {
            return;
        }
        for &e in &cb.indices {
            let r = self.ensure_pave_blocks(e);
            let list = &mut self.pave_blocks_pool[r];
            let mut new_list = Vec::new();
            for mut pb in list.drain(..) {
                if !pb.is_to_update() {
                    new_list.push(pb);
                } else {
                    pb.update(&mut new_list, true);
                }
            }
            *list = new_list;
        }
        self.common_blocks.push(cb.clone());
    }

    /// Returns the common blocks recorded by [`BopdsDS::update_common_block`].
    pub fn common_blocks(&self) -> &[BopdsCommonBlock] {
        &self.common_blocks
    }

    /// Records that pave blocks of edge `n_e` on `[t1, t2]` lie on face `n_f`.
    ///
    /// Port of `BOPAlgo_Tools::FillMap(aPB, nF)` + `PerformCommonBlocks(aMPBLI)`
    /// (`BOPAlgo_Tools.cxx`): the common block of that pave range gains the
    /// face, so a later `UpdateFaceInfoIn` can recover IN blocks via
    /// `CommonBlock::Contains(face)`.
    pub fn add_face_to_common_block(&mut self, n_e: usize, t1: f64, t2: f64, n_f: usize) {
        if let Some(cb) = self.common_blocks.iter_mut().find(|cb| {
            cb.contains_index(n_e) && cb.contains_range(t1, t2, 1e-7)
        }) {
            cb.add_face(n_f);
            return;
        }
        let mut cb = BopdsCommonBlock::new();
        cb.add_index(n_e);
        cb.add_range(t1, t2);
        cb.add_face(n_f);
        self.common_blocks.push(cb);
    }
}

impl Default for BopdsDS {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
#[path = "bopds_tests.rs"]
mod tests;
