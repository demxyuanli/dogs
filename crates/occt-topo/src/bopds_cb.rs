//! BOPDS common-block, face-info and interference types.
use crate::bopds_pave::BopdsPaveBlock;

/// Collection of pave blocks that have geometrical coincidence (within a
/// tolerance) with each other and/or with faces. Source: `BOPDS_CommonBlock`.
///
/// In this port the block is stored as the shared parameter ranges and the
/// indices of the edges whose pave blocks lie on those ranges.
#[derive(Debug, Clone)]
pub struct BopdsCommonBlock {
    /// Shared parameter ranges.
    pub ranges: Vec<(f64, f64)>,
    /// Indices of the edges sharing the ranges.
    pub indices: Vec<usize>,
    /// Faces this block lies on (`BOPDS_CommonBlock::myFaces`).
    pub faces: Vec<usize>,
    /// Tolerance of the common block.
    pub tolerance: f64,
    /// Pave blocks of this common block (`BOPDS_CommonBlock::myPaveBlocks`).
    pub pave_blocks: Vec<BopdsPaveBlock>,
}

impl BopdsCommonBlock {
    /// Empty constructor.
    pub fn new() -> Self {
        Self {
            ranges: Vec::new(),
            indices: Vec::new(),
            faces: Vec::new(),
            tolerance: 0.0,
            pave_blocks: Vec::new(),
        }
    }

    /// Adds a shared range.
    pub fn add_range(&mut self, first: f64, last: f64) {
        self.ranges.push((first, last));
    }

    /// Adds the index of an edge to the block.
    pub fn add_index(&mut self, index: usize) {
        if !self.indices.contains(&index) {
            self.indices.push(index);
        }
    }

    /// Adds the index of a face to the block (`BOPDS_CommonBlock::AddFace`).
    pub fn add_face(&mut self, face: usize) {
        if !self.contains_face(face) {
            self.faces.push(face);
        }
    }

    /// Returns the shared ranges.
    pub fn ranges(&self) -> &[(f64, f64)] {
        &self.ranges
    }

    /// Returns the indices of the edges of the block.
    pub fn indices(&self) -> &[usize] {
        &self.indices
    }

    /// Faces this block lies on.
    pub fn faces(&self) -> &[usize] {
        &self.faces
    }

    /// True if the block contains the edge index `i`.
    pub fn contains_index(&self, i: usize) -> bool {
        self.indices.contains(&i)
    }

    /// True if the block contains a range equal to `(first, last)` within `tol`.
    pub fn contains_range(&self, first: f64, last: f64, tol: f64) -> bool {
        self.ranges
            .iter()
            .any(|&(a, b)| (a - first).abs() <= tol && (b - last).abs() <= tol)
    }

    /// True if the block contains the face index `f`.
    pub fn contains_face(&self, f: usize) -> bool {
        self.faces.contains(&f)
    }

    /// Sets the tolerance of the block.
    pub fn set_tolerance(&mut self, tol: f64) {
        self.tolerance = tol;
    }

    /// Returns the tolerance of the block.
    pub fn tolerance(&self) -> f64 {
        self.tolerance
    }

    /// `BOPDS_CommonBlock::AddPaveBlock`.
    pub fn add_pave_block(&mut self, pb: BopdsPaveBlock) {
        let e = if pb.has_edge() {
            pb.edge()
        } else {
            pb.original_edge()
        };
        let (t1, t2) = pb.range();
        self.add_index(e);
        self.add_range(t1, t2);
        self.pave_blocks.push(pb);
    }

    /// `BOPDS_CommonBlock::SetPaveBlocks`.
    pub fn set_pave_blocks(&mut self, pbs: Vec<BopdsPaveBlock>) {
        self.pave_blocks.clear();
        self.indices.clear();
        self.ranges.clear();
        for pb in pbs {
            self.add_pave_block(pb);
        }
    }

    /// `BOPDS_CommonBlock::PaveBlocks`.
    pub fn pave_blocks(&self) -> &[BopdsPaveBlock] {
        &self.pave_blocks
    }

    /// `BOPDS_CommonBlock::PaveBlock1`.
    pub fn pave_block1(&self) -> Option<&BopdsPaveBlock> {
        self.pave_blocks.first()
    }

    /// `BOPDS_CommonBlock::SetFaces`.
    pub fn set_faces(&mut self, faces: Vec<usize>) {
        self.faces = faces;
    }

    /// `BOPDS_CommonBlock::AppendFaces`.
    pub fn append_faces(&mut self, faces: &[usize]) {
        for &n_f in faces {
            self.add_face(n_f);
        }
    }

    /// `BOPDS_CommonBlock::SetEdge` — assign the same edge index to every
    /// stored pave block.
    pub fn set_edge(&mut self, n_e: usize) {
        for pb in &mut self.pave_blocks {
            pb.set_edge(n_e);
        }
    }

    /// `BOPDS_CommonBlock::SetRealPaveBlock` — move `pb` to the front of the
    /// stored pave-block list so `PaveBlock1` returns it.
    pub fn set_real_pave_block(&mut self, pb: &BopdsPaveBlock) {
        if let Some(i) = self.pave_blocks.iter().position(|p| pb_same(p, pb)) {
            if i != 0 {
                self.pave_blocks.swap(0, i);
            }
        } else {
            self.pave_blocks.insert(0, pb.clone());
        }
    }
}

impl Default for BopdsCommonBlock {
    fn default() -> Self {
        Self::new()
    }
}

pub(crate) fn pb_same(a: &BopdsPaveBlock, b: &BopdsPaveBlock) -> bool {
    a.edge() == b.edge()
        && a.original_edge() == b.original_edge()
        && (a.first - b.first).abs() <= 1e-7
        && (a.last - b.last).abs() <= 1e-7
}

pub(crate) fn cb_contains_pb(cb: &BopdsCommonBlock, pb: &BopdsPaveBlock) -> bool {
    if cb.pave_blocks.iter().any(|p| pb_same(p, pb)) {
        return true;
    }
    // Face-only common blocks store indices+ranges without a pave-block list.
    // A CB that already lists pave blocks is one geometric interval; do not
    // match a different range of the same original via the index fallback
    // (`BOPDS_DS::CommonBlock` is keyed by the pave-block handle).
    if !cb.pave_blocks.is_empty() {
        return false;
    }
    let on_edge = cb.contains_index(pb.edge()) || cb.contains_index(pb.original_edge());
    on_edge && cb.contains_range(pb.first, pb.last, 1e-7)
}

// ---------------------------------------------------------------------------
// BopdsFaceInfo
// ---------------------------------------------------------------------------

/// Handy information about the state of a face. Source: `BOPDS_FaceInfo.hxx`.
///
/// The pave blocks on the face are split by kind, mirroring
/// `BOPDS_FaceInfo::PaveBlocksSc` / `PaveBlocksIn` / `PaveBlocksOn`:
/// - `paves` (`Sc`): section edges from face/face intersection;
/// - `paves_in` (`In`): coincident edges of the other operand lying on the
///   face (e.g. a cylinder rim on a box top). `BuildSplitFaces` appends them
///   as 1.2 In edges (`BOPAlgo_Builder_2.cxx`);
/// - `paves_on` (`On`): the face's own boundary-edge pave blocks, filled by
///   [`BopdsDS::update_face_info_on`]. `BOPDS_DS::RefineFaceInfoIn` drops the
///   IN blocks that are also ON.
/// - `paves_on` (`On`): pave blocks of the face's own boundary edges, built by
///   [`BopdsDS::update_face_info_on`]. `BOPDS_DS::RefineFaceInfoIn` drops the
///   IN blocks that are also ON (boundary) blocks.
///
/// Each pave-block entry is `(edge index, first parameter, last parameter)`.
/// Vertices lying on the face are kept separately in `verts` — the OCCT
/// `VerticesSc` / `VerticesOn` / `VerticesIn` maps. Mixing them into `paves`
/// would make the section-edge list un-usable as an edge list (a vertex index
/// is not an edge).
#[derive(Debug, Clone, Default)]
pub struct BopdsFaceInfo {
    /// Index of the face in the DS.
    pub face_index: usize,
    /// Section pave blocks lying on the face (F/F intersection edges).
    pub paves: Vec<(usize, f64, f64)>,
    /// IN pave blocks lying on the face (coincident edges of the other operand).
    pub paves_in: Vec<(usize, f64, f64)>,
    /// ON pave blocks lying on the face (the face's own boundary edges).
    pub paves_on: Vec<(usize, f64, f64)>,
    /// Vertices lying on the face (V/F and E/F hits), `(vertex index, u, v)`.
    pub verts: Vec<(usize, f64, f64)>,
    /// IN vertex indices (`BOPDS_FaceInfo::VerticesIn`): internal vertices,
    /// V/F hits, and E/F new vertices. Separate from [`verts`], which also
    /// holds UV data for On/Sc consumers.
    pub verts_in: Vec<usize>,
    /// ON vertex indices (`BOPDS_FaceInfo::VerticesOn`).
    pub verts_on: Vec<usize>,
    /// Section vertex indices (`BOPDS_FaceInfo::VerticesSc`).
    pub verts_sc: Vec<usize>,
}

impl BopdsFaceInfo {
    /// Constructor with the face index.
    pub fn new(face_index: usize) -> Self {
        Self {
            face_index,
            paves: Vec::new(),
            paves_in: Vec::new(),
            paves_on: Vec::new(),
            verts: Vec::new(),
            verts_in: Vec::new(),
            verts_on: Vec::new(),
            verts_sc: Vec::new(),
        }
    }

    /// Set the index of the face.
    pub fn set_index(&mut self, index: usize) {
        self.face_index = index;
    }

    /// Returns the index of the face.
    pub fn index(&self) -> usize {
        self.face_index
    }

    /// Adds a section pave block to the face.
    pub fn add_pave(&mut self, edge: usize, first: f64, last: f64) {
        self.paves.push((edge, first, last));
    }

    /// Adds an IN pave block (a coincident edge of the other operand) to the face.
    pub fn add_pave_in(&mut self, edge: usize, first: f64, last: f64) {
        self.paves_in.push((edge, first, last));
    }

    /// Returns the section pave blocks of the face.
    pub fn paves(&self) -> &[(usize, f64, f64)] {
        &self.paves
    }

    /// Returns the IN pave blocks of the face.
    pub fn paves_in(&self) -> &[(usize, f64, f64)] {
        &self.paves_in
    }

    /// Adds an IN pave block from a `BOPDS_PaveBlock` handle
    /// (`ChangePaveBlocksIn().Add(aPB)`).
    pub fn add_pave_block_in(&mut self, pb: &BopdsPaveBlock) {
        let e = if pb.has_edge() {
            pb.edge()
        } else {
            pb.original_edge()
        };
        let (first, last) = pb.range();
        if !self.paves_in.iter().any(|&(x, f, l)| {
            x == e && (f - first).abs() <= 1e-7 && (l - last).abs() <= 1e-7
        }) {
            self.add_pave_in(e, first, last);
        }
    }

    /// Adds an ON pave block (a pave block of the face's boundary edge).
    pub fn add_pave_on(&mut self, edge: usize, first: f64, last: f64) {
        self.paves_on.push((edge, first, last));
    }

    /// Returns the ON pave blocks of the face (its boundary-edge pave blocks).
    pub fn paves_on(&self) -> &[(usize, f64, f64)] {
        &self.paves_on
    }

    /// Adds a vertex lying on the face (`BOPDS_FaceInfo::VerticesSc`/`On`/`In`).
    pub fn add_vert(&mut self, vertex: usize, u: f64, v: f64) {
        self.verts.push((vertex, u, v));
    }

    /// Returns the vertices lying on the face.
    pub fn verts(&self) -> &[(usize, f64, f64)] {
        &self.verts
    }

    /// Adds an IN vertex index (`BOPDS_FaceInfo::ChangeVerticesIn`).
    pub fn add_vert_in(&mut self, vertex: usize) {
        if !self.verts_in.contains(&vertex) {
            self.verts_in.push(vertex);
        }
    }

    /// `BOPDS_FaceInfo::ChangeVerticesOn().Add`.
    pub fn add_vert_on(&mut self, vertex: usize) {
        if !self.verts_on.contains(&vertex) {
            self.verts_on.push(vertex);
        }
    }

    /// `BOPDS_FaceInfo::ChangeVerticesSc().Add`.
    pub fn add_vert_sc(&mut self, vertex: usize) {
        if !self.verts_sc.contains(&vertex) {
            self.verts_sc.push(vertex);
        }
    }

    /// Returns the IN vertex indices of the face.
    pub fn verts_in(&self) -> &[usize] {
        &self.verts_in
    }

    /// Returns the ON vertex indices of the face.
    pub fn verts_on(&self) -> &[usize] {
        &self.verts_on
    }

    /// Returns the section vertex indices of the face.
    pub fn verts_sc(&self) -> &[usize] {
        &self.verts_sc
    }
}

// ---------------------------------------------------------------------------
// BopdsInterf
// ---------------------------------------------------------------------------

/// Information about an interference between two shapes.
/// Source: `BOPDS_Interf.hxx` (root class).
#[derive(Debug, Clone, Copy)]
pub struct BopdsInterf {
    /// Index of the first interfered shape.
    pub index1: usize,
    /// Index of the second interfered shape.
    pub index2: usize,
    /// Index of a new shape created by the interference, or `-1`.
    pub index_new: i64,
    /// Common-part range on the first shape (`IntTools_CommonPrt::Range1`).
    pub common_first: f64,
    /// Common-part last parameter on the first shape.
    pub common_last: f64,
}

impl BopdsInterf {
    /// Constructor with the indices of the interfered shapes.
    pub fn new(index1: usize, index2: usize) -> Self {
        Self {
            index1,
            index2,
            index_new: -1,
            common_first: 0.0,
            common_last: 0.0,
        }
    }

    /// Set both indices.
    pub fn set_indices(&mut self, index1: usize, index2: usize) {
        self.index1 = index1;
        self.index2 = index2;
    }

    /// Returns both indices.
    pub fn indices(&self) -> (usize, usize) {
        (self.index1, self.index2)
    }

    /// Set the first index.
    pub fn set_index1(&mut self, index: usize) {
        self.index1 = index;
    }

    /// Set the second index.
    pub fn set_index2(&mut self, index: usize) {
        self.index2 = index;
    }

    /// Returns the first index.
    pub fn index1(&self) -> usize {
        self.index1
    }

    /// Returns the second index.
    pub fn index2(&self) -> usize {
        self.index2
    }

    /// Returns the index opposite to `i` (`None` when `i` is not a member).
    pub fn opposite_index(&self, i: usize) -> Option<usize> {
        if i == self.index1 {
            Some(self.index2)
        } else if i == self.index2 {
            Some(self.index1)
        } else {
            None
        }
    }

    /// True if the interference contains the index `i`.
    pub fn contains(&self, i: usize) -> bool {
        self.index1 == i || self.index2 == i
    }

    /// Set the index of a new shape.
    pub fn set_index_new(&mut self, index: usize) {
        self.index_new = index as i64;
    }

    /// Returns the index of a new shape (`-1` when unset).
    pub fn index_new(&self) -> i64 {
        self.index_new
    }

    /// True if a new-shape index is set.
    pub fn has_index_new(&self) -> bool {
        self.index_new >= 0
    }

    /// Returns the new-shape index, if set.
    pub fn get_index_new(&self) -> Option<usize> {
        if self.index_new >= 0 {
            Some(self.index_new as usize)
        } else {
            None
        }
    }

    /// Store `IntTools_CommonPrt::Range1` on this interference.
    pub fn set_common_range(&mut self, first: f64, last: f64) {
        self.common_first = first;
        self.common_last = last;
    }

    /// `IntTools_CommonPrt::Range1` stored with this interference.
    pub fn common_range(&self) -> (f64, f64) {
        (self.common_first, self.common_last)
    }
}
