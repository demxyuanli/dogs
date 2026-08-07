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

use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};
use std::hash::{Hash, Hasher};
use std::sync::Arc;

use occt_core::bnd::BndBox;
use occt_core::kernel::containers::IndexedMap;

use crate::abs::ShapeType;
use crate::bbox_from_geometry::shape_bbox;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, Shell, Solid, TopoShape, Vertex, Wire};
use crate::topo_tools_full::{edges_of, edges_of_wire, faces_of, shapes_of, vertices_of, wires_of_face};

// ---------------------------------------------------------------------------
// Shape-key wrapper
// ---------------------------------------------------------------------------

/// A `TopoShape` usable as a hash/equality key.
///
/// `TopoShape` does not implement `Hash`/`Eq`, so this newtype hashes and
/// compares by the shared `TShape` pointer — the same identity OCCT's
/// `TopTools_ShapeMapHasher` uses (topologically equal, not geometrically).
#[derive(Clone)]
struct ShapeKey(TopoShape);

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
}

impl BopdsShapeInfo {
    /// Constructor computing the type from the shape. Sub-shape indices are
    /// empty until filled by the DS.
    pub fn new(shape: TopoShape) -> Self {
        let kind = shape.shape_type();
        Self { shape, kind, sub_indices: Vec::new(), pb_reference: -1 }
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
        bopds_tools::has_brep(self.kind)
    }

    /// True if the shape can participate in an interference
    /// (boundary representation, or a solid).
    pub fn is_interfering(&self) -> bool {
        bopds_tools::is_interfering(self.kind)
    }
}

// ---------------------------------------------------------------------------
// BopdsCommonBlock
// ---------------------------------------------------------------------------

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
    /// Tolerance of the common block.
    pub tolerance: f64,
}

impl BopdsCommonBlock {
    /// Empty constructor.
    pub fn new() -> Self {
        Self { ranges: Vec::new(), indices: Vec::new(), tolerance: 0.0 }
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

    /// Adds the index of a face to the block.
    pub fn add_face(&mut self, face: usize) {
        if !self.contains_face(face) {
            self.indices.push(face);
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
        self.indices.contains(&f)
    }

    /// Sets the tolerance of the block.
    pub fn set_tolerance(&mut self, tol: f64) {
        self.tolerance = tol;
    }

    /// Returns the tolerance of the block.
    pub fn tolerance(&self) -> f64 {
        self.tolerance
    }
}

impl Default for BopdsCommonBlock {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// BopdsFaceInfo
// ---------------------------------------------------------------------------

/// Handy information about the state of a face. Source: `BOPDS_FaceInfo.hxx`.
///
/// The pave blocks on the face are split by kind, mirroring
/// `BOPDS_FaceInfo::PaveBlocksSc` / `PaveBlocksIn` / `PaveBlocksOn`:
/// - `paves` (`Sc`): section edges — the edges produced by the face/face
///   intersection. `BOPAlgo_Builder::BuildSplitFaces` uses only these to cut
///   the face into split pieces (plus the IN edges as internal shapes);
/// - `paves_in` (`In`): edges lying on the face (an edge of the other operand
///   coincident with the face, e.g. the triangulated base of a cylinder
///   standing on a planar face). They do not split the face.
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
    /// Vertices lying on the face (V/F and E/F hits), `(vertex index, u, v)`.
    pub verts: Vec<(usize, f64, f64)>,
}

impl BopdsFaceInfo {
    /// Constructor with the face index.
    pub fn new(face_index: usize) -> Self {
        Self { face_index, paves: Vec::new(), paves_in: Vec::new(), verts: Vec::new() }
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

    /// Adds a vertex lying on the face (`BOPDS_FaceInfo::VerticesSc`/`On`/`In`).
    pub fn add_vert(&mut self, vertex: usize, u: f64, v: f64) {
        self.verts.push((vertex, u, v));
    }

    /// Returns the vertices lying on the face.
    pub fn verts(&self) -> &[(usize, f64, f64)] {
        &self.verts
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
}

impl BopdsInterf {
    /// Constructor with the indices of the interfered shapes.
    pub fn new(index1: usize, index2: usize) -> Self {
        Self { index1, index2, index_new: -1 }
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
}

// ---------------------------------------------------------------------------
// bopds_tools
// ---------------------------------------------------------------------------

/// Type helper functions of the BOPDS package. Source: `BOPDS_Tools.hxx`.
pub mod bopds_tools {
    use occt_core::gp::GpPnt;

    use crate::abs::ShapeType;
    use crate::brep_tool::BRepTool;
    use crate::shape::{TopoShape, Vertex};

    /// Returns the type of the shape.
    pub fn shape_type(s: &TopoShape) -> ShapeType {
        s.shape_type()
    }

    /// True if the shape is a vertex.
    pub fn is_vertex(s: &TopoShape) -> bool {
        s.shape_type() == ShapeType::Vertex
    }

    /// True if the shape is an edge.
    pub fn is_edge(s: &TopoShape) -> bool {
        s.shape_type() == ShapeType::Edge
    }

    /// True if the shape is a face.
    pub fn is_face(s: &TopoShape) -> bool {
        s.shape_type() == ShapeType::Face
    }

    /// True if the shape is a wire.
    pub fn is_wire(s: &TopoShape) -> bool {
        s.shape_type() == ShapeType::Wire
    }

    /// True if the shape is a shell.
    pub fn is_shell(s: &TopoShape) -> bool {
        s.shape_type() == ShapeType::Shell
    }

    /// True if the shape is a solid.
    pub fn is_solid(s: &TopoShape) -> bool {
        s.shape_type() == ShapeType::Solid
    }

    /// The 3D point of a vertex (`BRep_Tool::Pnt`); `None` for non-vertices.
    pub fn vertex_point(s: &TopoShape) -> Option<GpPnt> {
        if is_vertex(s) {
            Some(BRepTool::vertex_point(&Vertex(s.clone())))
        } else {
            None
        }
    }

    /// True if the type corresponds to a shape having a boundary
    /// representation (vertex / edge / face).
    pub fn has_brep(t: ShapeType) -> bool {
        matches!(t, ShapeType::Vertex | ShapeType::Edge | ShapeType::Face)
    }

    /// True if the type can be a participant of an interference.
    pub fn is_interfering(t: ShapeType) -> bool {
        has_brep(t) || t == ShapeType::Solid
    }

    /// Converts a shape type to an integer (the OCCT `TopAbs_ShapeEnum` value).
    pub fn type_to_integer(t: ShapeType) -> i32 {
        t as i32
    }

    /// Converts the combination of two shape types to the index of the
    /// corresponding interference type (VV=0, VE=1, EE=2, VF=3, EF=4, FF=5,
    /// VZ=6, EZ=7, FZ=8, ZZ=9), or `-1` for an incompatible combination.
    pub fn type_to_integer2(t1: ShapeType, t2: ShapeType) -> i32 {
        let x = type_to_integer(t2) * 10 + type_to_integer(t1);
        match x {
            77 => 0,             // VV
            76 | 67 => 1,        // VE
            66 => 2,             // EE
            74 | 47 => 3,        // VF
            64 | 46 => 4,        // EF
            44 => 5,             // FF
            72 | 27 => 6,        // VZ
            62 | 26 => 7,        // EZ
            42 | 24 => 8,        // FZ
            22 => 9,             // ZZ
            _ => -1,
        }
    }
}

// ---------------------------------------------------------------------------
// Sub-shape extraction helpers
// ---------------------------------------------------------------------------

/// Direct vertex children of a shape (used for edges).
fn direct_vertex_children(shape: &TopoShape) -> Vec<TopoShape> {
    shape
        .tshape
        .read()
        .unwrap()
        .children
        .iter()
        .filter(|h| h.read().unwrap().shape_type() == ShapeType::Vertex)
        .map(|h| TopoShape::from_handle(h.clone()))
        .collect()
}

/// All direct children of a shape.
fn direct_children(shape: &TopoShape) -> Vec<TopoShape> {
    shape
        .tshape
        .read()
        .unwrap()
        .children
        .iter()
        .map(|h| TopoShape::from_handle(h.clone()))
        .collect()
}

/// The two boundary vertices of an edge.
fn edge_vertex_shapes(edge: &Edge) -> Vec<TopoShape> {
    direct_vertex_children(&edge.0)
}

/// The edges of a wire.
fn wire_edge_shapes(wire: &Wire) -> Vec<TopoShape> {
    edges_of_wire(wire).into_iter().map(|e| e.0).collect()
}

/// The boundary sub-shapes of a face: the distinct edges of its wires plus
/// its distinct vertices (wires themselves are not boundary sub-shapes).
fn face_sub_shapes(face: &Face) -> Vec<TopoShape> {
    let mut seen = HashSet::new();
    let mut out = Vec::new();
    for w in wires_of_face(face) {
        for e in edges_of_wire(&w) {
            if seen.insert(Arc::as_ptr(&e.0.tshape) as usize) {
                out.push(e.0);
            }
        }
    }
    for v in vertices_of(&face.0) {
        if seen.insert(Arc::as_ptr(&v.0.tshape) as usize) {
            out.push(v.0);
        }
    }
    out
}

/// The faces of a shell.
fn shell_face_shapes(shell: &Shell) -> Vec<TopoShape> {
    faces_of(&shell.0).into_iter().map(|f| f.0).collect()
}

/// The shells of a solid (structural children).
fn solid_shell_shapes(solid: &Solid) -> Vec<TopoShape> {
    shapes_of(&solid.0, ShapeType::Shell)
}

/// The prepared boundary sub-shapes of a solid: its faces and edges.
fn solid_face_edge_shapes(solid: &Solid) -> Vec<TopoShape> {
    let mut seen = HashSet::new();
    let mut out = Vec::new();
    for f in faces_of(&solid.0) {
        if seen.insert(Arc::as_ptr(&f.0.tshape) as usize) {
            out.push(f.0);
        }
    }
    for e in edges_of(&solid.0) {
        if seen.insert(Arc::as_ptr(&e.0.tshape) as usize) {
            out.push(e.0);
        }
    }
    out
}

/// The structural children used to walk the whole subtree of a shape while
/// appending it (a solid walks its shells, a face its edges and vertices).
fn structural_children(shape: &TopoShape) -> Vec<TopoShape> {
    match shape.shape_type() {
        ShapeType::Vertex => Vec::new(),
        ShapeType::Edge => edge_vertex_shapes(&Edge(shape.clone())),
        ShapeType::Wire => wire_edge_shapes(&Wire(shape.clone())),
        ShapeType::Face => face_sub_shapes(&Face(shape.clone())),
        ShapeType::Shell => shell_face_shapes(&Shell(shape.clone())),
        ShapeType::Solid => solid_shell_shapes(&Solid(shape.clone())),
        _ => direct_children(shape),
    }
}

/// The prepared boundary sub-shapes of a shape, in the form the DS stores
/// (matching OCCT after `BOPDS_DS::Init` + `prepareFaces`/`prepareSolids`):
/// a solid lists faces and edges, a shell its faces, a face its edges and
/// vertices, an edge its vertices.
fn prepared_sub_shapes(shape: &TopoShape) -> Vec<TopoShape> {
    match shape.shape_type() {
        ShapeType::Vertex => Vec::new(),
        ShapeType::Edge => edge_vertex_shapes(&Edge(shape.clone())),
        ShapeType::Wire => wire_edge_shapes(&Wire(shape.clone())),
        ShapeType::Face => face_sub_shapes(&Face(shape.clone())),
        ShapeType::Shell => shell_face_shapes(&Shell(shape.clone())),
        ShapeType::Solid => solid_face_edge_shapes(&Solid(shape.clone())),
        _ => direct_children(shape),
    }
}

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
    map: IndexedMap<ShapeKey>,
    /// Parallel array: per-shape information, index-aligned with `map`.
    shape_infos: Vec<BopdsShapeInfo>,
    /// Parallel array: cached bounding boxes, index-aligned with `map`.
    boxes: Vec<BndBox>,
    /// The arguments of the operation.
    arguments: Vec<TopoShape>,
    /// Per-argument index ranges.
    ranges: Vec<BopdsIndexRange>,
    /// Number of source shapes (arguments + sub-shapes present after `init`).
    nb_source_shapes: usize,
    /// Pool of pave-block lists, referenced by `BopdsShapeInfo::pb_reference`.
    pave_blocks_pool: Vec<Vec<BopdsPaveBlock>>,
    /// Pool of face states.
    face_info_pool: Vec<BopdsFaceInfo>,
    /// Common blocks collected by `update_common_block`.
    common_blocks: Vec<BopdsCommonBlock>,
    /// Same-domain shape map `index → same-domain index`.
    shapes_sd: HashMap<usize, usize>,
    /// Indices of the vertices whose tolerance was increased during the
    /// intersection (`BOPAlgo_PaveFiller::myIncreasedSS`). The repeat
    /// intersection stage re-runs V/V, V/E and V/F for these vertices.
    increased_ss: HashSet<usize>,
    /// Flat set of interfering pairs `(min, max)`. Mirror of OCCT's
    /// `myInterfTB` — membership queries only (`has_interf_pair`); the full
    /// per-type records (with the new-vertex index) live in `interf_vv..`.
    interferences: HashSet<(usize, usize)>,
    /// Set of interfered shape indices.
    interfered: HashSet<usize>,
    /// Typed interference records V/V. Source: `BOPDS_DS::InterfVV`.
    interf_vv: Vec<BopdsInterf>,
    /// Typed interference records V/E.
    interf_ve: Vec<BopdsInterf>,
    /// Typed interference records V/F.
    interf_vf: Vec<BopdsInterf>,
    /// Typed interference records E/E.
    interf_ee: Vec<BopdsInterf>,
    /// Typed interference records E/F.
    interf_ef: Vec<BopdsInterf>,
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
            .push(BopdsShapeInfo { shape: shape.clone(), kind, sub_indices: Vec::new(), pb_reference: -1 });
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

    /// Returns the face-info pool.
    pub fn face_info_pool(&self) -> &[BopdsFaceInfo] {
        &self.face_info_pool
    }

    /// Mutable access to the face-info pool.
    pub fn change_face_info_pool(&mut self) -> &mut Vec<BopdsFaceInfo> {
        &mut self.face_info_pool
    }

    /// Adds same-domain shape information (`index` and `index_sd`).
    pub fn add_shape_sd(&mut self, index: usize, index_sd: usize) {
        if index != index_sd {
            self.shapes_sd.insert(index, index_sd);
        }
    }

    /// Returns the same-domain index of `index` following the chain, or the
    /// index itself when no same-domain shape exists.
    pub fn get_same_domain_index(&self, index: usize) -> usize {
        let mut cur = index;
        while let Some(&next) = self.shapes_sd.get(&cur) {
            cur = next;
        }
        cur
    }

    /// Returns the final same-domain index of `index`, if any.
    pub fn has_shape_sd(&self, index: usize) -> Option<usize> {
        let result = self.get_same_domain_index(index);
        if result == index {
            None
        } else {
            Some(result)
        }
    }

    /// Returns the raw same-domain shape map (`index → same-domain index`).
    pub fn shapes_sd(&self) -> &HashMap<usize, usize> {
        &self.shapes_sd
    }

    /// Mutable access to the raw same-domain shape map.
    pub fn shapes_sd_mut(&mut self) -> &mut HashMap<usize, usize> {
        &mut self.shapes_sd
    }

    /// Returns the indices of the vertices whose tolerance was increased.
    pub fn increased_ss(&self) -> &HashSet<usize> {
        &self.increased_ss
    }

    /// Mutable access to the set of vertices with increased tolerance.
    pub fn increased_ss_mut(&mut self) -> &mut HashSet<usize> {
        &mut self.increased_ss
    }

    /// Adds an interference between shapes with indices `i1` and `i2`.
    /// Returns true when the pair was not already present.
    pub fn add_interf(&mut self, i1: usize, i2: usize) -> bool {
        let (a, b) = (i1.min(i2), i1.max(i2));
        if self.interferences.insert((a, b)) {
            self.interfered.insert(i1);
            self.interfered.insert(i2);
            true
        } else {
            false
        }
    }

    /// True if the shape with index `i` is interfered with any other shape.
    pub fn has_interf(&self, i: usize) -> bool {
        self.interfered.contains(&i)
    }

    /// True if the shapes with indices `i1` and `i2` are interfered.
    pub fn has_interf_pair(&self, i1: usize, i2: usize) -> bool {
        let (a, b) = (i1.min(i2), i1.max(i2));
        self.interferences.contains(&(a, b))
    }

    /// Returns the set of interfering pairs.
    pub fn interferences(&self) -> &HashSet<(usize, usize)> {
        &self.interferences
    }

    /// Returns the typed V/V interference records.
    pub fn interf_vv(&self) -> &[BopdsInterf] {
        &self.interf_vv
    }

    /// Returns the typed V/E interference records.
    pub fn interf_ve(&self) -> &[BopdsInterf] {
        &self.interf_ve
    }

    /// Returns the typed V/F interference records.
    pub fn interf_vf(&self) -> &[BopdsInterf] {
        &self.interf_vf
    }

    /// Returns the typed E/E interference records.
    pub fn interf_ee(&self) -> &[BopdsInterf] {
        &self.interf_ee
    }

    /// Returns the typed E/F interference records.
    pub fn interf_ef(&self) -> &[BopdsInterf] {
        &self.interf_ef
    }

    /// Appends a typed V/V interference record, registering the flat pair as
    /// well. The typed record is appended only when the pair was not interfered
    /// before (matching OCCT's `if (AddInterf(n1, n2)) { InterfVV().Appended() }`
    /// guard in `MakeSDVertices`). Returns true when a record was appended.
    pub fn add_interf_vv(&mut self, i1: usize, i2: usize, index_new: Option<usize>) -> bool {
        if !self.add_interf(i1, i2) {
            return false;
        }
        Self::append_interf(&mut self.interf_vv, i1, i2, index_new);
        true
    }

    /// Appends a typed V/E interference record (see [`BopdsDS::add_interf_vv`]).
    pub fn add_interf_ve(&mut self, i1: usize, i2: usize, index_new: Option<usize>) -> bool {
        if !self.add_interf(i1, i2) {
            return false;
        }
        Self::append_interf(&mut self.interf_ve, i1, i2, index_new);
        true
    }

    /// Appends a typed V/F interference record (see [`BopdsDS::add_interf_vv`]).
    pub fn add_interf_vf(&mut self, i1: usize, i2: usize, index_new: Option<usize>) -> bool {
        if !self.add_interf(i1, i2) {
            return false;
        }
        Self::append_interf(&mut self.interf_vf, i1, i2, index_new);
        true
    }

    /// Appends a typed E/E interference record (see [`BopdsDS::add_interf_vv`]).
    pub fn add_interf_ee(&mut self, i1: usize, i2: usize, index_new: Option<usize>) -> bool {
        if !self.add_interf(i1, i2) {
            return false;
        }
        Self::append_interf(&mut self.interf_ee, i1, i2, index_new);
        true
    }

    /// Appends a typed E/F interference record (see [`BopdsDS::add_interf_vv`]).
    pub fn add_interf_ef(&mut self, i1: usize, i2: usize, index_new: Option<usize>) -> bool {
        if !self.add_interf(i1, i2) {
            return false;
        }
        Self::append_interf(&mut self.interf_ef, i1, i2, index_new);
        true
    }

    fn append_interf(arr: &mut Vec<BopdsInterf>, i1: usize, i2: usize, index_new: Option<usize>) {
        let mut it = BopdsInterf::new(i1, i2);
        if let Some(n) = index_new {
            it.set_index_new(n);
        }
        arr.push(it);
    }

    /// Redirects the new-vertex index of every typed interference to its
    /// same-domain (SD) representative.
    ///
    /// Source: `BOPAlgo_PaveFiller::UpdateInterfsWithSDVertices` +
    /// `UpdateIntfsWithSDVertices` (`BOPAlgo_PaveFiller_10.cxx`): each typed
    /// interference whose `index_new` names a vertex that has since been merged
    /// into an SD cluster is re-pointed at the cluster representative. Called
    /// in `PerformInternal` after the E/F stage and again after `MakeBlocks`.
    pub fn update_interfs_with_sd_vertices(&mut self) {
        Self::update_intfs_with_sd_vertices(&self.shapes_sd, &mut self.interf_vv);
        Self::update_intfs_with_sd_vertices(&self.shapes_sd, &mut self.interf_ve);
        Self::update_intfs_with_sd_vertices(&self.shapes_sd, &mut self.interf_vf);
        Self::update_intfs_with_sd_vertices(&self.shapes_sd, &mut self.interf_ee);
        Self::update_intfs_with_sd_vertices(&self.shapes_sd, &mut self.interf_ef);
    }

    fn update_intfs_with_sd_vertices(sd: &HashMap<usize, usize>, interfs: &mut [BopdsInterf]) {
        for it in interfs.iter_mut() {
            let Some(n) = it.get_index_new() else { continue };
            let mut cur = n;
            while let Some(&next) = sd.get(&cur) {
                cur = next;
            }
            if cur != n {
                it.set_index_new(cur);
            }
        }
    }
}

impl Default for BopdsDS {
    fn default() -> Self {
        Self::new()
    }
}

/// Extended interference pairs for the vertices of `the_indices` whose
/// tolerance was increased.
///
/// Port of `BOPDS_Iterator::IntersectExt` (`BOPDS_Iterator.cxx`): every vertex
/// of `the_indices` (using its same-domain representative's box) is tested
/// against *every* interfering shape of the DS (using the representative box
/// for map members, the own box otherwise). The resulting pairs are returned
/// bucketed by interference type (VV..ZZ) — `BOPDS_Iterator::Initialize` then
/// serves these buckets *instead of* the regular pair lists, so the repeated
/// V/V, V/E and V/F stages process only the pairs involving the increased
/// vertices.
///
/// `the_indices` holds *source* vertex indices (as built by
/// `BOPAlgo_PaveFiller::RepeatIntersection`). The box-overlap test is done
/// brute-force (mirroring [`BopdsIterator::prepare`]) rather than with the
/// OCCT BVH tree — the semantics are identical.
pub fn intersect_ext_pairs(ds: &BopdsDS, the_indices: &HashSet<usize>) -> Vec<Vec<(usize, usize)>> {
    let mut buckets: Vec<Vec<(usize, usize)>> = vec![Vec::new(); NB_INTERF_TYPES];
    if the_indices.is_empty() {
        return buckets;
    }
    let n = ds.nb_source_shapes();
    // Resolved box of a shape: the same-domain representative's box for the
    // map members (OCCT `ShapeInfo(nVSD).Box()`), the own box otherwise
    // (OCCT `aSI.Box()` for the tree).
    let box_of = |k: usize, in_map: bool| -> BndBox {
        if in_map {
            ds.boxes[ds.get_same_domain_index(k)].clone()
        } else {
            ds.boxes[k].clone()
        }
    };
    let mut fence: HashSet<(usize, usize)> = HashSet::new();
    for i in 0..n {
        let si = &ds.shape_infos[i];
        // `IntersectExt` skips shapes without a boundary representation and
        // solids; `has_brep` is exactly the effective set.
        if !si.is_interfering() || si.kind == ShapeType::Solid {
            continue;
        }
        if !the_indices.contains(&i) {
            // Only map members are queries (every shape still lives in the tree).
            continue;
        }
        let box_i = box_of(i, true);
        for j in 0..n {
            if i == j {
                continue;
            }
            let sj = &ds.shape_infos[j];
            if !sj.is_interfering() || sj.kind == ShapeType::Solid {
                continue;
            }
            if ds.rank(i) == ds.rank(j) {
                continue;
            }
            let box_j = box_of(j, the_indices.contains(&j));
            if box_i.is_out_box(&box_j) {
                continue;
            }
            // Avoid interfering a shape with its own sub-shapes.
            let t1 = bopds_tools::type_to_integer(si.kind);
            let t2 = bopds_tools::type_to_integer(sj.kind);
            if (t1 < t2 && si.has_subshape(j)) || (t1 > t2 && sj.has_subshape(i)) {
                continue;
            }
            let ix = bopds_tools::type_to_integer2(si.kind, sj.kind);
            if (0..NB_INTERF_TYPES as i32).contains(&ix) {
                let (a, b) = (i.min(j), i.max(j));
                if fence.insert((a, b)) {
                    buckets[ix as usize].push((i, j));
                }
            }
        }
    }
    buckets
}

// ---------------------------------------------------------------------------
// BopdsIterator
// ---------------------------------------------------------------------------

/// Iterator over pairs of sub-shapes of the arguments whose bounding boxes
/// interfere. Source: `BOPDS_Iterator.hxx/.cxx`.
///
/// The regular iterator only reports *inter-argument* pairs (sub-shapes of
/// different arguments). Use [`BopdsIteratorSI`] for self-intersections.
#[derive(Debug)]
pub struct BopdsIterator<'a> {
    ds: Option<&'a BopdsDS>,
    /// Pair lists indexed by interference type (VV..ZZ).
    lists: Vec<Vec<(usize, usize)>>,
    /// The active (initialized) pair list.
    current: Vec<(usize, usize)>,
    /// Position in `current`.
    pos: usize,
    /// When true, self-intersection mode (no inter-argument skip, solids
    /// participate) — used by [`BopdsIteratorSI`].
    self_intersection: bool,
}

const NB_INTERF_TYPES: usize = 10;

impl<'a> BopdsIterator<'a> {
    /// Empty constructor (regular, inter-argument mode).
    pub fn new() -> Self {
        Self {
            ds: None,
            lists: vec![Vec::new(); NB_INTERF_TYPES],
            current: Vec::new(),
            pos: 0,
            self_intersection: false,
        }
    }

    /// Empty constructor in self-intersection mode.
    pub fn new_self_intersection() -> Self {
        Self {
            ds: None,
            lists: vec![Vec::new(); NB_INTERF_TYPES],
            current: Vec::new(),
            pos: 0,
            self_intersection: true,
        }
    }

    /// Sets the data structure to process.
    pub fn set_ds(&mut self, ds: &'a BopdsDS) {
        self.ds = Some(ds);
    }

    /// Returns the data structure.
    pub fn ds(&self) -> &'a BopdsDS {
        self.ds.expect("BopdsIterator: DS not set")
    }

    /// Initializes the iterator for the given pair of shape types.
    pub fn initialize(&mut self, t1: ShapeType, t2: ShapeType) {
        self.pos = 0;
        let ix = bopds_tools::type_to_integer2(t1, t2);
        if ix >= 0 {
            let mut list = self.lists[ix as usize].clone();
            list.sort_unstable();
            self.current = list;
        } else {
            self.current.clear();
        }
    }

    /// True if there are still pairs to iterate.
    pub fn more(&self) -> bool {
        self.pos < self.current.len()
    }

    /// Moves the iteration ahead.
    pub fn next(&mut self) {
        self.pos += 1;
    }

    /// Returns the indices (DS) of the interfering shapes, ordered so the
    /// shape of the higher type comes first.
    pub fn value(&self) -> (usize, usize) {
        let (i, j) = self.current[self.pos];
        let ds = self.ds();
        let t1 = bopds_tools::type_to_integer(ds.shape_infos[i].kind);
        let t2 = bopds_tools::type_to_integer(ds.shape_infos[j].kind);
        if t1 < t2 {
            (j, i)
        } else {
            (i, j)
        }
    }

    /// Performs the interference (bounding-box) computation and prepares the
    /// pair lists to be used.
    pub fn prepare(&mut self) {
        for list in &mut self.lists {
            list.clear();
        }
        self.current.clear();
        self.pos = 0;
        let Some(ds) = self.ds else { return };
        let n = ds.nb_source_shapes();
        for i in 0..n {
            let si = &ds.shape_infos[i];
            let ok = if self.self_intersection {
                bopds_tools::is_interfering(si.kind)
            } else {
                bopds_tools::has_brep(si.kind)
            };
            if !ok {
                continue;
            }
            for j in (i + 1)..n {
                let sj = &ds.shape_infos[j];
                let okj = if self.self_intersection {
                    bopds_tools::is_interfering(sj.kind)
                } else {
                    bopds_tools::has_brep(sj.kind)
                };
                if !okj {
                    continue;
                }
                // Inter-argument mode: skip pairs within the same argument.
                if !self.self_intersection && ds.rank(i) == ds.rank(j) {
                    continue;
                }
                // Bounding-box rejection.
                if ds.boxes[i].is_out_box(&ds.boxes[j]) {
                    continue;
                }
                // Avoid interfering a shape with its own sub-shapes.
                let t1 = bopds_tools::type_to_integer(si.kind);
                let t2 = bopds_tools::type_to_integer(sj.kind);
                if (t1 < t2 && si.has_subshape(j)) || (t1 > t2 && sj.has_subshape(i)) {
                    continue;
                }
                let ix = bopds_tools::type_to_integer2(si.kind, sj.kind);
                if (0..NB_INTERF_TYPES as i32).contains(&ix) {
                    self.lists[ix as usize].push((i, j));
                }
            }
        }
    }

    /// Returns the number of pairs in the currently initialized list.
    pub fn expected_length(&self) -> usize {
        self.current.len()
    }

    /// Returns the block length (half the expected length, at least 1).
    pub fn block_length(&self) -> usize {
        let n = self.expected_length();
        if n <= 1 {
            1
        } else {
            (0.5 * n as f64) as usize
        }
    }

    /// Clears the pair lists for interference types above `level`.
    /// Level semantics: 0 — only V/V; 1 — +V/E; 2 — +E/E; 3 — +V/F;
    /// 4 — +E/F; anything else — all interferences.
    pub fn update_by_level_of_check(&mut self, level: i32) {
        for k in (level + 1).max(0) as usize..self.lists.len() {
            self.lists[k].clear();
        }
    }
}

impl Default for BopdsIterator<'_> {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// BopdsSubIterator
// ---------------------------------------------------------------------------

/// Iterator over pairs of interfering sub-shapes of two given sub-sets of DS
/// indices. Source: `BOPDS_SubIterator.hxx/.cxx`.
#[derive(Debug)]
pub struct BopdsSubIterator<'a> {
    ds: Option<&'a BopdsDS>,
    /// First set of indices.
    subset1: Vec<usize>,
    /// Second set of indices.
    subset2: Vec<usize>,
    /// The interfering pairs found by [`BopdsSubIterator::prepare`].
    pairs: Vec<(usize, usize)>,
    /// Position in `pairs`.
    pos: usize,
}

impl<'a> BopdsSubIterator<'a> {
    /// Empty constructor.
    pub fn new() -> Self {
        Self { ds: None, subset1: Vec::new(), subset2: Vec::new(), pairs: Vec::new(), pos: 0 }
    }

    /// Sets the data structure to process.
    pub fn set_ds(&mut self, ds: &'a BopdsDS) {
        self.ds = Some(ds);
    }

    /// Returns the data structure.
    pub fn ds(&self) -> &'a BopdsDS {
        self.ds.expect("BopdsSubIterator: DS not set")
    }

    /// Sets the first set of indices to process.
    pub fn set_subset1(&mut self, subset: Vec<usize>) {
        self.subset1 = subset;
    }

    /// Returns the first set of indices.
    pub fn subset1(&self) -> &[usize] {
        &self.subset1
    }

    /// Sets the second set of indices to process.
    pub fn set_subset2(&mut self, subset: Vec<usize>) {
        self.subset2 = subset;
    }

    /// Returns the second set of indices.
    pub fn subset2(&self) -> &[usize] {
        &self.subset2
    }

    /// Initializes the iterator (sorts the pairs).
    pub fn initialize(&mut self) {
        self.pairs.sort_unstable();
        self.pos = 0;
    }

    /// True if there are more pairs of interfering shapes.
    pub fn more(&self) -> bool {
        self.pos < self.pairs.len()
    }

    /// Moves the iteration ahead.
    pub fn next(&mut self) {
        self.pos += 1;
    }

    /// Returns the indices (DS) of the interfering shapes.
    pub fn value(&self) -> (usize, usize) {
        let (i, j) = self.pairs[self.pos];
        let ds = self.ds();
        let t1 = bopds_tools::type_to_integer(ds.shape_infos[i].kind);
        let t2 = bopds_tools::type_to_integer(ds.shape_infos[j].kind);
        if t1 < t2 {
            (j, i)
        } else {
            (i, j)
        }
    }

    /// Performs the intersection of bounding boxes of the two sub-sets and
    /// prepares the results.
    pub fn prepare(&mut self) {
        self.pairs.clear();
        let Some(ds) = self.ds else { return };
        if self.subset1.is_empty() || self.subset2.is_empty() {
            return;
        }
        let mut fence = HashSet::new();
        for &i in &self.subset1 {
            if i >= ds.shape_infos.len() {
                continue;
            }
            let si = &ds.shape_infos[i];
            for &j in &self.subset2 {
                if j >= ds.shape_infos.len() || i == j {
                    continue;
                }
                let (a, b) = (i.min(j), i.max(j));
                if !fence.insert((a, b)) {
                    continue;
                }
                if ds.boxes[i].is_out_box(&ds.boxes[j]) {
                    continue;
                }
                let sj = &ds.shape_infos[j];
                let t1 = bopds_tools::type_to_integer(si.kind);
                let t2 = bopds_tools::type_to_integer(sj.kind);
                if (t1 < t2 && si.has_subshape(j)) || (t1 > t2 && sj.has_subshape(i)) {
                    continue;
                }
                self.pairs.push((a, b));
            }
        }
        self.pairs.sort_unstable();
        self.pos = 0;
    }

    /// Returns the number of interfering pairs.
    pub fn expected_length(&self) -> usize {
        self.pairs.len()
    }
}

impl Default for BopdsSubIterator<'_> {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// BopdsIteratorSI
// ---------------------------------------------------------------------------

/// Iterator computing *self-intersections* between the sub-shapes of each
/// argument. Source: `BOPDS_IteratorSI.hxx/.cxx`.
///
/// Unlike [`BopdsIterator`], it considers pairs inside the same argument and
/// lets solids participate, so it can detect an argument that intersects
/// itself (e.g. two overlapping boxes).
#[derive(Debug)]
pub struct BopdsIteratorSI<'a> {
    inner: BopdsIterator<'a>,
}

impl<'a> BopdsIteratorSI<'a> {
    /// Empty constructor (self-intersection mode).
    pub fn new() -> Self {
        Self { inner: BopdsIterator::new_self_intersection() }
    }

    /// Sets the data structure to process.
    pub fn set_ds(&mut self, ds: &'a BopdsDS) {
        self.inner.set_ds(ds);
    }

    /// Clears the pair lists for interference types above `level`
    /// (see [`BopdsIterator::update_by_level_of_check`]).
    pub fn update_by_level_of_check(&mut self, level: i32) {
        self.inner.update_by_level_of_check(level);
    }
}

impl<'a> std::ops::Deref for BopdsIteratorSI<'a> {
    type Target = BopdsIterator<'a>;
    fn deref(&self) -> &Self::Target {
        &self.inner
    }
}

impl<'a> std::ops::DerefMut for BopdsIteratorSI<'a> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.inner
    }
}

impl Default for BopdsIteratorSI<'_> {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::brep_extrema::test_box::unit_box;

    fn box_counts(ds: &BopdsDS) -> (usize, usize, usize, usize, usize) {
        let (mut v, mut e, mut f, mut sh, mut so) = (0, 0, 0, 0, 0);
        for si in &ds.shape_infos {
            match si.kind {
                ShapeType::Vertex => v += 1,
                ShapeType::Edge => e += 1,
                ShapeType::Face => f += 1,
                ShapeType::Shell => sh += 1,
                ShapeType::Solid => so += 1,
                _ => {}
            }
        }
        (v, e, f, sh, so)
    }

    #[test]
    fn append_box_counts_all_subshapes() {
        let mut ds = BopdsDS::new();
        let b = unit_box();
        let idx = ds.append(b.solid.0.clone()).unwrap();
        assert_eq!(idx, 0);
        // 8 vertices + 12 edges + 6 faces + 1 shell + 1 solid = 28.
        assert_eq!(ds.nb_shapes(), 28);
        assert_eq!(box_counts(&ds), (8, 12, 6, 1, 1));
    }

    #[test]
    fn append_is_idempotent_and_index_roundtrips() {
        let mut ds = BopdsDS::new();
        let b = unit_box();
        ds.append(b.solid.0.clone()).unwrap();
        let again = ds.append(b.solid.0.clone()).unwrap();
        assert_eq!(again, 0);
        assert_eq!(ds.nb_shapes(), 28);

        // Every distinct sub-shape round-trips through index()/shape().
        for v in &b.vertices {
            let i = ds.index(&v.0).expect("vertex indexed");
            assert!(ds.shape(i).unwrap().same_tshape(&v.0));
        }
        for e in &b.edges {
            let i = ds.index(&e.0).expect("edge indexed");
            assert!(ds.shape(i).unwrap().same_tshape(&e.0));
        }
        for f in &b.faces {
            let i = ds.index(&f.0).expect("face indexed");
            assert!(ds.shape(i).unwrap().same_tshape(&f.0));
        }
        assert_eq!(ds.index(&b.solid.0), Some(0));
        assert_eq!(ds.index(&b.solid.0), ds.index(&b.solid.0));
    }

    #[test]
    fn sub_shape_indices_are_prepared() {
        let mut ds = BopdsDS::new();
        let b = unit_box();
        ds.append(b.solid.0.clone()).unwrap();

        // A face's sub-shapes are its 4 edges + 4 vertices (no wires).
        let face_idx = ds.index(&b.faces[0].0).unwrap();
        let face_info = ds.shape_info(face_idx).unwrap();
        let n_edges = face_info
            .sub_indices
            .iter()
            .filter(|&&i| ds.shape_info(i).unwrap().kind == ShapeType::Edge)
            .count();
        let n_verts = face_info
            .sub_indices
            .iter()
            .filter(|&&i| ds.shape_info(i).unwrap().kind == ShapeType::Vertex)
            .count();
        assert_eq!(n_edges, 4);
        assert_eq!(n_verts, 4);
        // An edge's sub-shapes are its two vertices.
        let edge_idx = ds.index(&b.edges[0].0).unwrap();
        let edge_info = ds.shape_info(edge_idx).unwrap();
        assert_eq!(edge_info.sub_indices.len(), 2);
        assert!(edge_info.sub_indices.iter().all(|&i| ds.shape_info(i).unwrap().kind == ShapeType::Vertex));
    }

    #[test]
    fn ranks_follow_arguments() {
        let mut ds = BopdsDS::new();
        let a = unit_box();
        let c = unit_box();
        ds.init(&[a.solid.0.clone(), c.solid.0.clone()]);
        assert_eq!(ds.nb_shapes(), 56);
        assert_eq!(ds.nb_ranges(), 2);
        assert_eq!(ds.nb_source_shapes(), 56);
        // Shapes of the first argument have rank 0, second rank 1.
        let e0 = ds.index(&a.edges[0].0).unwrap();
        let e1 = ds.index(&c.edges[0].0).unwrap();
        assert_eq!(ds.rank(e0), 0);
        assert_eq!(ds.rank(e1), 1);
        // New shapes (beyond source shapes) are not source shapes.
        assert!(ds.is_new_shape(ds.nb_source_shapes()));
        assert!(!ds.is_new_shape(e0));
    }

    #[test]
    fn pave_blocks_split_on_update() {
        let mut ds = BopdsDS::new();
        let b = unit_box();
        ds.append(b.solid.0.clone()).unwrap();
        let e_idx = ds.index(&b.edges[0].0).unwrap();
        let v0 = ds.index(&b.vertices[0].0).unwrap();
        let v1 = ds.index(&b.vertices[1].0).unwrap();
        let v4 = ds.index(&b.vertices[4].0).unwrap();
        let v5 = ds.index(&b.vertices[5].0).unwrap();

        assert!(!ds.has_pave_blocks(e_idx));
        assert!(ds.pave_blocks(e_idx).is_empty());

        // Add two pave blocks; the second carries an extra pave at t = 0.5.
        {
            let pbs = ds.change_pave_blocks_mut(e_idx);
            let mut pb1 = BopdsPaveBlock::new();
            pb1.edge_index = e_idx;
            pb1.original_edge = e_idx;
            pb1.index1 = v0;
            pb1.index2 = v1;
            pb1.first = 0.0;
            pb1.last = 1.0;

            let mut pb2 = BopdsPaveBlock::new();
            pb2.edge_index = e_idx;
            pb2.original_edge = e_idx;
            pb2.index1 = v4;
            pb2.index2 = v5;
            pb2.first = 0.0;
            pb2.last = 1.0;
            pb2.append_ext_pave(BopdsPave::new(v0, 0.5));

            pbs.push(pb1);
            pbs.push(pb2);
        }

        assert!(ds.has_pave_blocks(e_idx));
        assert_eq!(ds.pave_blocks(e_idx).len(), 2);
        assert!(ds.pave_blocks(e_idx)[1].is_to_update());

        ds.update_pave_blocks();

        // pb1 untouched (1 block) + pb2 split into [0, 0.5] and [0.5, 1] = 3.
        let blocks = ds.pave_blocks(e_idx);
        assert_eq!(blocks.len(), 3);
        assert_eq!(blocks[1].range(), (0.0, 0.5));
        assert_eq!(blocks[2].range(), (0.5, 1.0));
        assert!(!ds.pave_blocks(e_idx)[0].is_to_update());
    }

    #[test]
    fn common_block_updates_and_pool_grows() {
        let mut ds = BopdsDS::new();
        let b = unit_box();
        ds.append(b.solid.0.clone()).unwrap();
        let e_idx = ds.index(&b.edges[0].0).unwrap();

        let mut cb = BopdsCommonBlock::new();
        cb.add_range(0.0, 0.5);
        cb.add_index(e_idx);
        assert!(cb.contains_index(e_idx));
        assert!(cb.contains_range(0.0, 0.5, 1e-9));
        assert!(!cb.contains_range(0.2, 0.8, 1e-9));

        ds.update_common_block(&cb);
        assert_eq!(ds.common_blocks().len(), 1);
        // update_common_block lazily allocated a pave-block slot for the edge.
        assert!(ds.has_pave_blocks(e_idx));
    }

    #[test]
    fn tools_report_shape_types() {
        let b = unit_box();
        assert_eq!(bopds_tools::shape_type(&b.solid.0), ShapeType::Solid);
        assert!(bopds_tools::is_solid(&b.solid.0));
        assert!(bopds_tools::is_vertex(&b.vertices[0].0));
        assert!(bopds_tools::is_edge(&b.edges[0].0));
        assert!(bopds_tools::is_face(&b.faces[0].0));

        let shell = shapes_of(&b.solid.0, ShapeType::Shell);
        assert_eq!(shell.len(), 1);
        assert!(bopds_tools::is_shell(&shell[0]));

        let p = bopds_tools::vertex_point(&b.vertices[0].0);
        assert!(p.is_some());
        assert!((p.unwrap().x() - 0.0).abs() < 1e-12);
        assert!(bopds_tools::vertex_point(&b.solid.0).is_none());

        // Interference-type encoding.
        assert_eq!(bopds_tools::type_to_integer(ShapeType::Vertex), 7);
        assert_eq!(bopds_tools::type_to_integer2(ShapeType::Vertex, ShapeType::Vertex), 0);
        assert_eq!(bopds_tools::type_to_integer2(ShapeType::Edge, ShapeType::Face), 4);
        assert_eq!(bopds_tools::type_to_integer2(ShapeType::Face, ShapeType::Face), 5);
        assert_eq!(bopds_tools::type_to_integer2(ShapeType::Solid, ShapeType::Solid), 9);
    }

    #[test]
    fn pave_and_range_value_semantics() {
        let a = BopdsPave::new(1, 0.5);
        let b = BopdsPave::new(2, 1.0);
        assert!(a.is_less(&b));
        assert!(a < b);
        assert_eq!(a, BopdsPave::new(1, 0.5));
        assert_ne!(a, BopdsPave::new(1, 0.6));
        assert_ne!(a, BopdsPave::new(3, 0.5));

        let r = BopdsIndexRange::new(2, 7);
        assert!(r.contains(5));
        assert!(!r.contains(8));
        assert_eq!(r.indices(), (2, 7));

        let info = BopdsShapeInfo::new(crate::shape::TopoShape::new(ShapeType::Vertex));
        assert_eq!(info.shape_type(), ShapeType::Vertex);
        assert!(info.has_brep());
        // A vertex can participate in V/V, V/E and V/F interferences.
        assert!(info.is_interfering());

        let mut fi = BopdsFaceInfo::new(3);
        fi.add_pave(1, 0.0, 0.5);
        assert_eq!(fi.index(), 3);
        assert_eq!(fi.paves(), &[(1, 0.0, 0.5)]);

        let mut interf = BopdsInterf::new(4, 9);
        assert!(interf.contains(4));
        assert_eq!(interf.opposite_index(4), Some(9));
        assert_eq!(interf.opposite_index(8), None);
        interf.set_index_new(12);
        assert_eq!(interf.get_index_new(), Some(12));
    }

    #[test]
    fn append_info_registers_shape() {
        let mut ds = BopdsDS::new();
        let b = unit_box();
        ds.append(b.solid.0.clone()).unwrap();
        let before = ds.nb_shapes();
        let info = BopdsShapeInfo::new(b.vertices[0].0.clone());
        let i = ds.append_info(info);
        // Already present → returns the existing index.
        assert_eq!(i, ds.index(&b.vertices[0].0).unwrap());
        assert_eq!(ds.nb_shapes(), before);
    }

    #[test]
    fn iterator_finds_inter_argument_candidates() {
        let mut ds = BopdsDS::new();
        let a = unit_box();
        let c = unit_box();
        ds.init(&[a.solid.0.clone(), c.solid.0.clone()]);

        let mut it = BopdsIterator::new();
        it.set_ds(&ds);
        it.prepare();
        it.initialize(ShapeType::Face, ShapeType::Face);
        assert!(it.more(), "overlapping faces of two identical boxes must interfere");
        let (i, j) = it.value();
        // Inter-argument mode never pairs shapes of the same rank.
        assert_ne!(ds.rank(i), ds.rank(j));
        // Each of the 6 faces overlaps the 5 non-parallel faces of the other
        // coincident box → 6 × 5 = 30 face/face candidate pairs.
        assert_eq!(it.expected_length(), 30);
    }

    #[test]
    fn si_iterator_reports_overlapping_cubes() {
        let mut ds = BopdsDS::new();
        let a = unit_box();
        let c = unit_box();
        ds.init(&[a.solid.0.clone(), c.solid.0.clone()]);

        let mut it = BopdsIteratorSI::new();
        it.set_ds(&ds);
        it.prepare();
        it.initialize(ShapeType::Vertex, ShapeType::Vertex);
        assert!(it.more(), "self-intersection VV candidates must be non-empty");
        let (i, j) = it.value();
        assert_ne!(i, j);

        // Face/face candidates are non-empty as well.
        it.initialize(ShapeType::Face, ShapeType::Face);
        assert!(it.more());
        // And solid/solid.
        it.initialize(ShapeType::Solid, ShapeType::Solid);
        assert!(it.more());
        let (s1, s2) = it.value();
        assert_eq!(ds.shape_info(s1).unwrap().kind, ShapeType::Solid);
        assert_eq!(ds.shape_info(s2).unwrap().kind, ShapeType::Solid);
    }

    #[test]
    fn sub_iterator_limited_to_given_subsets() {
        let mut ds = BopdsDS::new();
        let a = unit_box();
        let c = unit_box();
        ds.init(&[a.solid.0.clone(), c.solid.0.clone()]);

        let vertices: Vec<usize> = ds
            .shape_infos
            .iter()
            .enumerate()
            .filter(|(_, si)| si.kind == ShapeType::Vertex)
            .map(|(i, _)| i)
            .collect();
        assert_eq!(vertices.len(), 16);

        let mut it = BopdsSubIterator::new();
        it.set_ds(&ds);
        it.set_subset1(vertices.clone());
        it.set_subset2(vertices.clone());
        it.prepare();
        assert!(it.more());
        let (i, j) = it.value();
        assert!(ds.shape_info(i).unwrap().kind == ShapeType::Vertex);
        assert!(ds.shape_info(j).unwrap().kind == ShapeType::Vertex);

        // Restricting to the vertices of a single box yields no pairs
        // (distinct corners do not overlap).
        let first_box_vertices: Vec<usize> = vertices[..8].to_vec();
        let mut it2 = BopdsSubIterator::new();
        it2.set_ds(&ds);
        it2.set_subset1(first_box_vertices.clone());
        it2.set_subset2(first_box_vertices);
        it2.prepare();
        assert!(!it2.more());
    }

    #[test]
    fn interference_tracking() {
        let mut ds = BopdsDS::new();
        assert!(ds.add_interf(3, 7));
        assert!(!ds.add_interf(7, 3)); // unordered duplicate
        assert!(ds.has_interf(3));
        assert!(ds.has_interf(7));
        assert!(ds.has_interf_pair(7, 3));
        assert!(!ds.has_interf(5));
        assert!(!ds.has_interf_pair(3, 5));

        ds.add_shape_sd(1, 4);
        ds.add_shape_sd(4, 9);
        assert_eq!(ds.get_same_domain_index(1), 9);
        assert_eq!(ds.has_shape_sd(1), Some(9));
        assert_eq!(ds.get_same_domain_index(2), 2);
    }

    #[test]
    fn typed_interferences_and_sd_redirection() {
        let mut ds = BopdsDS::new();
        // A V/V record with a new-vertex index.
        assert!(ds.add_interf_vv(3, 7, Some(5)));
        // Unordered duplicate: no second flat pair, no second typed record.
        assert!(!ds.add_interf_vv(7, 3, Some(5)));
        assert_eq!(ds.interf_vv().len(), 1);
        assert_eq!(ds.interf_vv()[0].get_index_new(), Some(5));
        // The flat table still tracks membership.
        assert!(ds.has_interf_pair(3, 7));
        assert!(ds.has_interf(7));

        // A V/E record whose new vertex is later merged into an SD cluster.
        assert!(ds.add_interf_ve(1, 8, Some(6)));
        ds.add_shape_sd(6, 12);
        ds.add_shape_sd(12, 20);
        ds.update_interfs_with_sd_vertices();
        assert_eq!(
            ds.interf_ve()[0].get_index_new(),
            Some(20),
            "redirected to the final SD representative"
        );
        // The V/V new vertex has no SD partner — left untouched.
        assert_eq!(ds.interf_vv()[0].get_index_new(), Some(5));

        // An interference with no new-vertex index stays untouched.
        assert!(ds.add_interf_ee(2, 9, None));
        ds.update_interfs_with_sd_vertices();
        assert_eq!(ds.interf_ee()[0].get_index_new(), None);
    }
}
