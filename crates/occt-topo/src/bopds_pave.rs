//! BOPDS pave, pave-block, index-range and shape-info types.
use std::cmp::Ordering;
use std::collections::HashSet;
use std::hash::{Hash, Hasher};
use std::sync::Arc;

use crate::abs::ShapeType;
use crate::shape::TopoShape;

// Shape-key wrapper
// ---------------------------------------------------------------------------

/// A `TopoShape` usable as a hash/equality key.
///
/// `TopoShape` does not implement `Hash`/`Eq`, so this newtype hashes and
/// compares by the shared `TShape` pointer — the same identity OCCT's
/// `TopTools_ShapeMapHasher` uses (topologically equal, not geometrically).
#[derive(Clone)]
pub(crate) struct ShapeKey(pub(crate) TopoShape);

impl PartialEq for ShapeKey {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.0.tshape, &other.0.tshape)
    }
}

impl Eq for ShapeKey {}

impl Hash for ShapeKey {
    fn hash<H: Hasher>(&self, state: &mut H) {
        (Arc::as_ptr(&self.0.tshape) as usize).hash(state);
    }
}

impl std::fmt::Debug for ShapeKey {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "ShapeKey({:p})", self.0.tshape)
    }
}

// ---------------------------------------------------------------------------
// BopdsPave
// ---------------------------------------------------------------------------

/// Information about a vertex on an edge. Source: `BOPDS_Pave.hxx`.
#[derive(Debug, Clone, Copy)]
pub struct BopdsPave {
    /// Index of the vertex in the DS.
    pub index: usize,
    /// Parameter of the vertex on the edge's curve.
    pub param: f64,
}

impl BopdsPave {
    /// Constructor with index and parameter.
    pub fn new(index: usize, param: f64) -> Self {
        Self { index, param }
    }

    /// Set the index of the vertex.
    pub fn set_index(&mut self, index: usize) {
        self.index = index;
    }

    /// Returns the index of the vertex.
    pub fn index(&self) -> usize {
        self.index
    }

    /// Set the parameter of the vertex.
    pub fn set_parameter(&mut self, param: f64) {
        self.param = param;
    }

    /// Returns the parameter of the vertex.
    pub fn parameter(&self) -> f64 {
        self.param
    }

    /// Returns the index and parameter.
    pub fn contents(&self) -> (usize, f64) {
        (self.index, self.param)
    }

    /// True if `self.param < other.param` (ordering by parameter only).
    pub fn is_less(&self, other: &Self) -> bool {
        self.param < other.param
    }

    /// True if both index and parameter are equal.
    pub fn is_equal(&self, other: &Self) -> bool {
        self.index == other.index && self.param == other.param
    }
}

impl Default for BopdsPave {
    fn default() -> Self {
        Self::new(0, 99.0)
    }
}

impl PartialEq for BopdsPave {
    fn eq(&self, other: &Self) -> bool {
        self.is_equal(other)
    }
}

impl Eq for BopdsPave {}

impl PartialOrd for BopdsPave {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.param.partial_cmp(&other.param)
    }
}

impl Ord for BopdsPave {
    fn cmp(&self, other: &Self) -> Ordering {
        self.param.partial_cmp(&other.param).unwrap_or(Ordering::Equal)
    }
}

impl Hash for BopdsPave {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.index.hash(state);
        self.param.to_bits().hash(state);
    }
}

// ---------------------------------------------------------------------------
// BopdsPaveBlock
// ---------------------------------------------------------------------------

/// Information about an interval (pave block) on an edge.
/// Source: `BOPDS_PaveBlock.hxx/.cxx`.
///
/// Two adjacent paves on an edge make up a pave block. A block can carry
/// extra paves ([`BopdsPaveBlock::append_ext_pave`]); [`BopdsPaveBlock::update`]
/// then splits it into the elementary blocks delimited by those paves.
#[derive(Debug, Clone)]
pub struct BopdsPaveBlock {
    /// Index of the edge the block currently lies on (0 if unset).
    pub edge_index: usize,
    /// Index of the original edge the block was created for.
    pub original_edge: usize,
    /// First parameter of the block.
    pub first: f64,
    /// Last parameter of the block.
    pub last: f64,
    /// DS index of the vertex at the first bound.
    pub index1: usize,
    /// DS index of the vertex at the last bound.
    pub index2: usize,
    /// Extra paves (potential split points) on the block.
    pub(crate) ext_paves: Vec<BopdsPave>,
    /// Fence of vertex indices already present in `ext_paves`.
    pub(crate) fence: HashSet<usize>,
    /// Shrunk range data.
    pub(crate) ts1: f64,
    pub(crate) ts2: f64,
    pub(crate) has_shrunk: bool,
    pub(crate) is_splittable: bool,
}

impl BopdsPaveBlock {
    /// Empty constructor.
    pub fn new() -> Self {
        Self {
            edge_index: 0,
            original_edge: 0,
            first: 0.0,
            last: 0.0,
            index1: 0,
            index2: 0,
            ext_paves: Vec::new(),
            fence: HashSet::new(),
            ts1: -99.0,
            ts2: -99.0,
            has_shrunk: false,
            is_splittable: false,
        }
    }

    /// Set the index of the edge of the pave block.
    pub fn set_edge(&mut self, edge: usize) {
        self.edge_index = edge;
    }

    /// Returns the index of the edge of the pave block.
    pub fn edge(&self) -> usize {
        self.edge_index
    }

    /// True if the pave block has an edge assigned.
    pub fn has_edge(&self) -> bool {
        self.edge_index != usize::MAX
    }

    /// Set the index of the original edge.
    pub fn set_original_edge(&mut self, edge: usize) {
        self.original_edge = edge;
    }

    /// Returns the index of the original edge.
    pub fn original_edge(&self) -> usize {
        self.original_edge
    }

    /// True if the block lies on a split edge (edge differs from original).
    pub fn is_split_edge(&self) -> bool {
        self.edge_index != self.original_edge
    }

    /// Set the first bound as a pave.
    pub fn set_pave1(&mut self, pave: BopdsPave) {
        self.index1 = pave.index;
        self.first = pave.param;
    }

    /// Set the last bound as a pave.
    pub fn set_pave2(&mut self, pave: BopdsPave) {
        self.index2 = pave.index;
        self.last = pave.param;
    }

    /// Returns the first bound as a pave.
    pub fn pave1(&self) -> BopdsPave {
        BopdsPave::new(self.index1, self.first)
    }

    /// Returns the last bound as a pave.
    pub fn pave2(&self) -> BopdsPave {
        BopdsPave::new(self.index2, self.last)
    }

    /// Returns the parametric range `(t1, t2)` of the pave block.
    pub fn range(&self) -> (f64, f64) {
        (self.first, self.last)
    }

    /// Set the parametric range.
    pub fn set_range(&mut self, first: f64, last: f64) {
        self.first = first;
        self.last = last;
    }

    /// Returns the pave indices `(index1, index2)` of the block.
    pub fn indices(&self) -> (usize, usize) {
        (self.index1, self.index2)
    }

    /// Set the pave indices.
    pub fn set_indices(&mut self, index1: usize, index2: usize) {
        self.index1 = index1;
        self.index2 = index2;
    }

    /// True if this block has the same bound indices as `other`
    /// (in either orientation).
    pub fn has_same_bounds(&self, other: &Self) -> bool {
        let (n11, n12) = self.indices();
        let (n21, n22) = other.indices();
        (n11 == n21 && n12 == n22) || (n11 == n22 && n12 == n21)
    }

    /// True if the block contains extra paves and must be updated (split).
    pub fn is_to_update(&self) -> bool {
        !self.ext_paves.is_empty()
    }

    /// Append an extra pave, skipping duplicates by vertex index.
    pub fn append_ext_pave(&mut self, pave: BopdsPave) {
        if self.fence.insert(pave.index) {
            self.ext_paves.push(pave);
        }
    }

    /// Append an extra pave unconditionally.
    pub fn append_ext_pave1(&mut self, pave: BopdsPave) {
        self.ext_paves.push(pave);
    }

    /// Remove all extra paves referencing the given vertex index.
    pub fn remove_ext_pave(&mut self, vert_num: usize) {
        if self.fence.remove(&vert_num) {
            self.ext_paves.retain(|p| p.index != vert_num);
        }
    }

    /// Returns the extra paves.
    pub fn ext_paves(&self) -> &[BopdsPave] {
        &self.ext_paves
    }

    /// Mutable access to the extra paves.
    pub fn change_ext_paves(&mut self) -> &mut Vec<BopdsPave> {
        &mut self.ext_paves
    }

    /// True if an extra pave has a parameter within `tol` of `t`.
    /// Returns the index of the found pave.
    pub fn contains_parameter(&self, t: f64, tol: f64) -> Option<usize> {
        self.ext_paves
            .iter()
            .find(|p| (p.param - t).abs() < tol)
            .map(|p| p.index)
    }

    /// Update the pave block. The extra paves (plus, when `flag`, the two
    /// bound paves) are used to create the new elementary pave blocks `out`.
    /// Source: `BOPDS_PaveBlock::Update`.
    pub fn update(&mut self, out: &mut Vec<BopdsPaveBlock>, flag: bool) {
        let mut paves: Vec<BopdsPave> = Vec::new();
        if flag {
            paves.push(self.pave1());
            paves.push(self.pave2());
        }
        paves.append(&mut self.ext_paves);
        self.fence.clear();

        if paves.len() <= 1 {
            return;
        }

        paves.sort();
        let original = self.original_edge;
        for w in paves.windows(2) {
            let p1 = w[0];
            let p2 = w[1];
            let mut pb = BopdsPaveBlock::new();
            pb.original_edge = original;
            pb.edge_index = original;
            pb.index1 = p1.index;
            pb.index2 = p2.index;
            pb.first = p1.param;
            pb.last = p2.param;
            out.push(pb);
        }
    }

    /// Set the shrunk range data.
    pub fn set_shrunk_data(&mut self, ts1: f64, ts2: f64, is_splittable: bool) {
        self.ts1 = ts1;
        self.ts2 = ts2;
        self.has_shrunk = true;
        self.is_splittable = is_splittable;
    }

    /// Drop shrunk data so `HasShrunkData` is false.
    ///
    /// OCCT `HasShrunkData` is `!myShrunkBox.IsVoid()`. `AnalyzeShrunkData`
    /// stores an empty box when `!IsDone()`, which makes `HasShrunkData`
    /// return false even though `SetShrunkData` was called.
    pub fn clear_shrunk_data(&mut self) {
        self.ts1 = 0.0;
        self.ts2 = 0.0;
        self.has_shrunk = false;
        self.is_splittable = false;
    }

    /// Returns the shrunk range data `(ts1, ts2, is_splittable)`.
    pub fn shrunk_data(&self) -> (f64, f64, bool) {
        (self.ts1, self.ts2, self.is_splittable)
    }

    /// True if the block has shrunk data.
    pub fn has_shrunk_data(&self) -> bool {
        self.has_shrunk
    }

    /// True if the shrunk range is long enough for the edge to be split.
    pub fn is_splittable(&self) -> bool {
        self.is_splittable
    }
}

impl Default for BopdsPaveBlock {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// BopdsIndexRange
// ---------------------------------------------------------------------------

/// A range of two indices `[first, last]`. Source: `BOPDS_IndexRange.hxx`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BopdsIndexRange {
    /// First index of the range.
    pub first: usize,
    /// Last index of the range.
    pub last: usize,
}

impl BopdsIndexRange {
    /// Constructor with initial indices.
    pub fn new(first: usize, last: usize) -> Self {
        Self { first, last }
    }

    /// Set the first index.
    pub fn set_first(&mut self, first: usize) {
        self.first = first;
    }

    /// Set the last index.
    pub fn set_last(&mut self, last: usize) {
        self.last = last;
    }

    /// Returns the first index.
    pub fn first(&self) -> usize {
        self.first
    }

    /// Returns the last index.
    pub fn last(&self) -> usize {
        self.last
    }

    /// Set both indices.
    pub fn set_indices(&mut self, first: usize, last: usize) {
        self.first = first;
        self.last = last;
    }

    /// Returns both indices.
    pub fn indices(&self) -> (usize, usize) {
        (self.first, self.last)
    }

    /// True if `index` lies inside the range.
    pub fn contains(&self, index: usize) -> bool {
        index >= self.first && index <= self.last
    }
}

impl Default for BopdsIndexRange {
    fn default() -> Self {
        Self::new(0, 0)
    }
}

// ---------------------------------------------------------------------------
// BopdsShapeInfo
// ---------------------------------------------------------------------------

/// Handy information about a shape stored in the DS.
/// Source: `BOPDS_ShapeInfo.hxx/.lxx`.
#[derive(Debug, Clone)]
pub struct BopdsShapeInfo {
    /// The shape itself.
    pub shape: TopoShape,
    /// The type of the shape.
    pub kind: ShapeType,
    /// Indices of the boundary sub-shapes (prepared form: a solid lists its
    /// faces and edges, a face its edges and vertices, an edge its vertices).
    pub sub_indices: Vec<usize>,
    /// Reference into the pave-block pool; `-1` when none. (In OCCT this is
    /// the single `Reference` field shared with the face-info pool; here it is
    /// used for pave blocks, and the face-info pool is indexed separately.)
    pub(crate) pb_reference: i64,
    /// Optional flag (`BOPDS_ShapeInfo::myFlag`). Degenerated edges store the
    /// parent-face index; infinite-curve vertices store `1`.
    pub flag: Option<usize>,
}

impl BopdsShapeInfo {
    /// Constructor computing the type from the shape. Sub-shape indices are
    /// empty until filled by the DS.
    pub fn new(shape: TopoShape) -> Self {
        let kind = shape.shape_type();
        Self {
            shape,
            kind,
            sub_indices: Vec::new(),
            pb_reference: -1,
            flag: None,
        }
    }

    /// Returns the shape.
    pub fn shape(&self) -> &TopoShape {
        &self.shape
    }

    /// Returns the type of the shape.
    pub fn shape_type(&self) -> ShapeType {
        self.kind
    }

    /// Returns the indices of the boundary sub-shapes.
    pub fn sub_shapes(&self) -> &[usize] {
        &self.sub_indices
    }

    /// Mutable access to the sub-shape indices.
    pub fn change_sub_shapes(&mut self) -> &mut Vec<usize> {
        &mut self.sub_indices
    }

    /// True if the shape has a sub-shape with index `i`.
    pub fn has_subshape(&self, i: usize) -> bool {
        self.sub_indices.contains(&i)
    }

    /// True if the shape has a boundary representation (vertex/edge/face).
    pub fn has_brep(&self) -> bool {
        crate::bopds_tools::has_brep(self.kind)
    }

    /// True if the shape can participate in an interference
    /// (boundary representation, or a solid).
    pub fn is_interfering(&self) -> bool {
        crate::bopds_tools::is_interfering(self.kind)
    }

    /// `BOPDS_ShapeInfo::HasFlag`.
    pub fn has_flag(&self) -> bool {
        self.flag.is_some()
    }

    /// `BOPDS_ShapeInfo::HasFlag(int&)`.
    pub fn flag(&self) -> Option<usize> {
        self.flag
    }

    /// `BOPDS_ShapeInfo::SetFlag`.
    pub fn set_flag(&mut self, flag: usize) {
        self.flag = Some(flag);
    }
}
