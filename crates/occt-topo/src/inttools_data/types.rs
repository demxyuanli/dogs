use super::prelude::*;


/// Tolerance used when merging adjacent intervals (`MarkedRangeSet`,
/// `LocalizeData`). Two intervals are considered adjacent / touching when their
/// boundaries are within this distance.

pub(super) const MERGE_EPS: f64 = 1e-9;

// ---------------------------------------------------------------------------
// IntRange / IntToolsRange
// ---------------------------------------------------------------------------

/// A closed 1D range `[first, last]`.
///
/// Mirrors `IntTools_Range` (a pair of `double` bounds). Ranges are inclusive
/// on both ends: `contains` accepts both boundary parameters. No ordering of
/// `first`/`last` is implied beyond `first <= last`, which [`IntRange::new`]
/// enforces.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IntRange {
    pub first: f64,
    pub last: f64,
}

/// OCCT-named alias for [`IntRange`].
pub type IntToolsRange = IntRange;

impl IntRange {
    /// Validated constructor. `Err` when either bound is non-finite or
    /// `first > last`.
    pub fn new(first: f64, last: f64) -> Result<Self, String> {
        if !first.is_finite() || !last.is_finite() {
            return Err(format!(
                "IntRange: non-finite bound ({first}, {last})"
            ));
        }
        if first > last {
            return Err(format!(
                "IntRange: first {first} is greater than last {last}"
            ));
        }
        Ok(Self { first, last })
    }

    /// Unchecked constructor for internally-known-valid ranges. Prefer [`new`].
    ///
    /// [`new`]: Self::new
    pub const fn new_unchecked(first: f64, last: f64) -> Self {
        Self { first, last }
    }

    /// Degenerate `[0, 0]` range.
    pub const fn zero() -> Self {
        Self { first: 0.0, last: 0.0 }
    }

    pub fn first(&self) -> f64 {
        self.first
    }

    pub fn last(&self) -> f64 {
        self.last
    }

    /// `last - first`.
    pub fn length(&self) -> f64 {
        self.last - self.first
    }

    /// True when both bounds are finite and `first <= last`.
    pub fn is_valid(&self) -> bool {
        self.first.is_finite() && self.last.is_finite() && self.first <= self.last
    }

    /// Whether `t` lies inside the closed interval (both ends inclusive).
    pub fn contains(&self, t: f64) -> bool {
        self.first <= t && t <= self.last
    }

    /// Whether two closed ranges share at least one point. Touching ranges
    /// (`[0,1]` and `[1,2]`) overlap at the shared boundary.
    pub fn overlaps(&self, other: &Self) -> bool {
        self.first <= other.last && other.first <= self.last
    }

    /// Whether the ranges are disjoint (no shared point, boundaries excluded).
    pub fn disjoint(&self, other: &Self) -> bool {
        !self.overlaps(other)
    }

    /// The smallest range spanning both inputs (`min first`, `max last`).
    /// Meaningful for overlapping/adjacent ranges; for disjoint ranges it
    /// returns the bounding hull (the two gaps are not covered).
    pub fn merge(&self, other: &Self) -> Self {
        Self {
            first: self.first.min(other.first),
            last: self.last.max(other.last),
        }
    }
}

impl Default for IntRange {
    fn default() -> Self {
        Self::zero()
    }
}

// ---------------------------------------------------------------------------
// CommonPartType / CommonPrt
// ---------------------------------------------------------------------------

/// The kind of a common part between two shapes.
///
/// OCCT stores this as a `TopAbs_ShapeEnum` (`VERTEX`/`EDGE`/`FACE`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum CommonPartType {
    /// Uninitialized / not yet classified.
    #[default]
    Unknown,
    /// `TopAbs_VERTEX`: a point common part (touch or piercing).
    Vertex,
    /// `TopAbs_EDGE`: a coincident sub-range of an edge on a face (or on
    /// another edge).
    Edge,
    /// The common part is a face region.
    Face,
}

/// A common part between two edges (or an edge and a face), described by a
/// range and the shapes it touches.
///
/// Mirrors `IntTools_CommonPrt`. The C++ class stores two edges, a primary
/// range, a sequence of secondary ranges, vertex parameters and bounding
/// points; this port keeps the payload the Phase 16 boolean layer consumes:
/// the part type, the primary range, the face the part lies on (when any), and
/// the vertex shapes at the part ends (when any).
#[derive(Debug, Clone)]
pub struct CommonPrt {
    pub part_type: CommonPartType,
    pub range: IntRange,
    pub face: Option<TopoShape>,
    pub vertices: Vec<TopoShape>,
    /// `IntTools_CommonPrt::VertexParameter1`.
    pub vertex_parameter1: Option<f64>,
    /// `IntTools_CommonPrt::AllNullFlag`.
    pub all_null_flag: bool,
    /// `IntTools_CommonPrt::SetBoundingPoints` / `BoundingPoints`.
    pub bounding_p1: GpPnt,
    pub bounding_p2: GpPnt,
}

impl CommonPrt {
    /// Empty common part of [`CommonPartType::Unknown`].
    pub fn new() -> Self {
        Self {
            part_type: CommonPartType::Unknown,
            range: IntRange::zero(),
            face: None,
            vertices: Vec::new(),
            vertex_parameter1: None,
            all_null_flag: false,
            bounding_p1: GpPnt::zero(),
            bounding_p2: GpPnt::zero(),
        }
    }

    /// Fully-specified constructor.
    pub fn with(
        part_type: CommonPartType,
        range: IntRange,
        face: Option<TopoShape>,
        vertices: Vec<TopoShape>,
    ) -> Self {
        Self {
            part_type,
            range,
            face,
            vertices,
            vertex_parameter1: None,
            all_null_flag: false,
            bounding_p1: GpPnt::zero(),
            bounding_p2: GpPnt::zero(),
        }
    }

    pub fn part_type(&self) -> CommonPartType {
        self.part_type
    }

    pub fn set_part_type(&mut self, t: CommonPartType) {
        self.part_type = t;
    }

    pub fn range(&self) -> IntRange {
        self.range
    }

    pub fn set_range(&mut self, r: IntRange) {
        self.range = r;
    }

    /// `IntTools_CommonPrt::SetBoundingPoints`.
    pub fn set_bounding_points(&mut self, p1: GpPnt, p2: GpPnt) {
        self.bounding_p1 = p1;
        self.bounding_p2 = p2;
    }

    /// `IntTools_CommonPrt::BoundingPoints`.
    pub fn bounding_points(&self) -> (GpPnt, GpPnt) {
        (self.bounding_p1, self.bounding_p2)
    }

    pub fn face(&self) -> Option<&TopoShape> {
        self.face.as_ref()
    }

    pub fn set_face(&mut self, f: Option<TopoShape>) {
        self.face = f;
    }

    pub fn vertices(&self) -> &[TopoShape] {
        &self.vertices
    }

    /// Appends a vertex shape.
    pub fn add_vertex(&mut self, v: TopoShape) {
        self.vertices.push(v);
    }

    /// `IntTools_CommonPrt::SetVertexParameter1`.
    pub fn set_vertex_parameter1(&mut self, t: f64) {
        self.vertex_parameter1 = Some(t);
    }

    /// `IntTools_CommonPrt::VertexParameter1`.
    pub fn vertex_parameter1(&self) -> Option<f64> {
        self.vertex_parameter1
    }

    /// `IntTools_CommonPrt::SetAllNullFlag`.
    pub fn set_all_null_flag(&mut self, f: bool) {
        self.all_null_flag = f;
    }

    /// `IntTools_CommonPrt::AllNullFlag`.
    pub fn all_null_flag(&self) -> bool {
        self.all_null_flag
    }

    /// True when the common part is classified as a point (`TopAbs_VERTEX`).
    pub fn is_vertex_part(&self) -> bool {
        self.part_type == CommonPartType::Vertex
    }

    /// True when the common part is classified as an edge sub-range.
    pub fn is_edge_part(&self) -> bool {
        self.part_type == CommonPartType::Edge
    }

    /// True when the common part is classified as a face region.
    pub fn is_face_part(&self) -> bool {
        self.part_type == CommonPartType::Face
    }

    /// True when the part has no face and no vertices yet (only a type/range).
    pub fn is_empty(&self) -> bool {
        self.face.is_none() && self.vertices.is_empty()
    }
}

impl Default for CommonPrt {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// RootType / IntRoot
// ---------------------------------------------------------------------------

/// The nature of a root of a 1D intersection function.
///
/// OCCT's `IntTools_Root` stores an integer `Type()` (0 = simple/bisection,
/// 1 = pure zero interval, 2 = smart/Fibonacci). This port classifies the root
/// geometrically instead.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum RootType {
    /// The function has a transversal root here.
    IsRoot,
    /// The function touches zero without crossing (tangency).
    IsTangent,
    /// The classification could not be determined.
    #[default]
    IsUnknown,
}

/// A root of the edge/edge or edge/surface function.
///
/// Mirrors `IntTools_Root`. Instead of a scalar root value, this port carries
/// the interval `[first, last]` the root was found in (the bisection window),
/// a diagnostic `root_index`, the [`RootType`] and a conflict flag. Keeping the
/// interval rather than the collapsed value lets the caller split edges exactly
/// at the root's uncertainty window.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IntRoot {
    pub root_index: i32,
    pub root_type: RootType,
    pub range: IntRange,
    pub is_conflict: bool,
}

impl IntRoot {
    pub fn new(root_index: i32, root_type: RootType, range: IntRange) -> Self {
        Self { root_index, root_type, range, is_conflict: false }
    }

    pub fn with_conflict(
        root_index: i32,
        root_type: RootType,
        range: IntRange,
        is_conflict: bool,
    ) -> Self {
        Self { root_index, root_type, range, is_conflict }
    }

    pub fn root_index(&self) -> i32 {
        self.root_index
    }

    pub fn root_type(&self) -> RootType {
        self.root_type
    }

    pub fn range(&self) -> IntRange {
        self.range
    }

    pub fn is_conflict(&self) -> bool {
        self.is_conflict
    }

    pub fn set_conflict(&mut self, v: bool) {
        self.is_conflict = v;
    }

    /// The mid-point of the root interval — the best single-parameter estimate
    /// of the root.
    pub fn root_value(&self) -> f64 {
        0.5 * (self.range.first + self.range.last)
    }

    pub fn is_root(&self) -> bool {
        self.root_type == RootType::IsRoot
    }

    pub fn is_tangent(&self) -> bool {
        self.root_type == RootType::IsTangent
    }
}

impl Default for IntRoot {
    fn default() -> Self {
        Self::new(0, RootType::IsUnknown, IntRange::zero())
    }
}

// ---------------------------------------------------------------------------
// PntOnFace / PntOn2Faces
// ---------------------------------------------------------------------------

/// A 3D point on a face together with its UV parameters.
///
/// Mirrors `IntTools_PntOnFace`, which stores a `TopoDS_Face` handle; this port
/// stores a `face_index` into the caller's face array instead, avoiding shape
/// handles in hot geometric data. `uv` is the `(U, V)` surface parameter pair.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PntOnFace {
    pub face_index: usize,
    pub uv: (f64, f64),
    pub pnt: GpPnt,
}

impl PntOnFace {
    pub fn new(face_index: usize, uv: (f64, f64), pnt: GpPnt) -> Self {
        Self { face_index, uv, pnt }
    }

    pub fn face_index(&self) -> usize {
        self.face_index
    }

    pub fn u(&self) -> f64 {
        self.uv.0
    }

    pub fn v(&self) -> f64 {
        self.uv.1
    }

    pub fn uv(&self) -> (f64, f64) {
        self.uv
    }

    pub fn pnt(&self) -> &GpPnt {
        &self.pnt
    }
}

impl Default for PntOnFace {
    fn default() -> Self {
        Self::new(0, (0.0, 0.0), GpPnt::zero())
    }
}

/// A pair of points on two faces, used to represent a touching/intersection
/// location sampled on both surfaces.
///
/// Mirrors `IntTools_PntOn2Faces` (two `IntTools_PntOnFace` values). Each side
/// keeps its face index, 3D point and UV parameters.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PntOn2Faces {
    pub face1_index: usize,
    pub face2_index: usize,
    pub pnt1: GpPnt,
    pub pnt2: GpPnt,
    pub uv1: (f64, f64),
    pub uv2: (f64, f64),
}

impl PntOn2Faces {
    pub fn new(
        face1_index: usize,
        face2_index: usize,
        pnt1: GpPnt,
        pnt2: GpPnt,
        uv1: (f64, f64),
        uv2: (f64, f64),
    ) -> Self {
        Self { face1_index, face2_index, pnt1, pnt2, uv1, uv2 }
    }

    pub fn face1_index(&self) -> usize {
        self.face1_index
    }

    pub fn face2_index(&self) -> usize {
        self.face2_index
    }

    pub fn pnt1(&self) -> &GpPnt {
        &self.pnt1
    }

    pub fn pnt2(&self) -> &GpPnt {
        &self.pnt2
    }

    pub fn uv1(&self) -> (f64, f64) {
        self.uv1
    }

    pub fn uv2(&self) -> (f64, f64) {
        self.uv2
    }
}

impl Default for PntOn2Faces {
    fn default() -> Self {
        Self::new(0, 0, GpPnt::zero(), GpPnt::zero(), (0.0, 0.0), (0.0, 0.0))
    }
}

// ---------------------------------------------------------------------------
// CurveKind / IntCurve
// ---------------------------------------------------------------------------

/// The analytic family of an intersection curve.
///
/// Mirrors OCCT's `GeomAbs_CurveType` (as returned by `IntTools_Curve::Type`)
/// with the `BezierCurve` case folded into `BSpline`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum CurveKind {
    Line,
    Circle,
    Ellipse,
    Parabola,
    Hyperbola,
    BSpline,
    #[default]
    Other,
}

/// A face/face intersection curve container.
///
/// Mirrors `IntTools_Curve`. The C++ class holds a 3D curve handle plus two 2D
/// pcurve handles and tolerance values; this port keeps the classification the
/// boolean layer needs: the curve family, the two faces it came from and the
/// parameter range over which it is valid.
#[derive(Debug, Clone)]
pub struct IntCurve {
    pub kind: CurveKind,
    pub face1: Option<TopoShape>,
    pub face2: Option<TopoShape>,
    pub range: IntRange,
}

impl IntCurve {
    pub fn new(kind: CurveKind, face1: Option<TopoShape>, face2: Option<TopoShape>, range: IntRange) -> Self {
        Self { kind, face1, face2, range }
    }

    pub fn kind(&self) -> CurveKind {
        self.kind
    }

    pub fn face1(&self) -> Option<&TopoShape> {
        self.face1.as_ref()
    }

    pub fn face2(&self) -> Option<&TopoShape> {
        self.face2.as_ref()
    }

    pub fn range(&self) -> IntRange {
        self.range
    }

    pub fn set_range(&mut self, r: IntRange) {
        self.range = r;
    }
}

impl Default for IntCurve {
    fn default() -> Self {
        Self::new(CurveKind::Other, None, None, IntRange::zero())
    }
}

// ---------------------------------------------------------------------------
// MarkedRangeSet
// ---------------------------------------------------------------------------

/// A sorted, non-overlapping set of 1D ranges, each carrying a mark flag.
///
/// Mirrors `IntTools_MarkedRangeSet`, which keeps a sorted boundary list plus a
/// parallel flag sequence. In OCCT the set is used to track which parameter
/// sub-ranges of an edge are already accounted for by intersection results
/// (`InsertRange` with a flag). This port exposes the same idea with a simpler
/// surface:
///
/// - [`insert`](Self::insert) adds an *unmarked* range, replacing (and
///   unmarking) any existing coverage it overlaps;
/// - [`mark`](Self::mark) sets the flag of the range containing a parameter;
/// - [`unite`](Self::unite) / [`intersect`](Self::intersect) /
///   [`subtract`](Self::subtract) combine two sets, propagating marks (a
///   location is marked in the result when the corresponding point is marked in
///   the source(s), per set operation);
/// - [`sorted_ranges`](Self::sorted_ranges) returns the (sorted, adjacent-
///   merged) intervals.
///
/// The internal invariant is always "sorted, non-overlapping, adjacent intervals
/// merged when their flags agree", restored after every mutation.
#[derive(Debug, Clone, PartialEq)]
pub struct MarkedRangeSet {
    pub(super) ranges: Vec<(IntRange, bool)>,
}

/// A `(lo, hi, flag)` triple used by the internal interval algebra.
pub(super) type Flagged = (f64, f64, bool);

/// Sort `v`, drop degenerate (≤ `MERGE_EPS`) intervals and merge adjacent
/// intervals that carry the same flag.
pub(super) fn merge_flagged(mut v: Vec<Flagged>) -> Vec<Flagged> {
    v.retain(|(lo, hi, _)| hi - lo > MERGE_EPS);
    if v.is_empty() {
        return v;
    }
    v.sort_by(|a, b| {
        a.0.partial_cmp(&b.0)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
    });
    let mut out: Vec<Flagged> = vec![v[0]];
    for (lo, hi, f) in v.into_iter().skip(1) {
        let last = out.last_mut().unwrap();
        if lo <= last.1 + MERGE_EPS && f == last.2 {
            last.1 = last.1.max(hi);
        } else {
            out.push((lo, hi, f));
        }
    }
    out
}

#[derive(Clone, Copy)]
pub(super) enum CombineOp {
    Unite,
    Intersect,
    Subtract,
}

/// Boolean combination of two flagged interval sets via breakpoint sweep.
///
/// Every candidate sub-interval `[x, y]` (bounded by breakpoints from either
/// set) is classified by sampling its midpoint in both inputs, then the result
/// flag is derived from the operation:
///
/// - `Unite`: covered by either; marked when marked in either.
/// - `Intersect`: covered by both; marked when marked in both.
/// - `Subtract`: covered by `a` and not `b`; keeps `a`'s mark.
pub(super) fn combine(a: &[Flagged], b: &[Flagged], op: CombineOp) -> Vec<Flagged> {
    let mut pts: Vec<f64> = Vec::with_capacity(2 * (a.len() + b.len()));
    for (lo, hi, _) in a.iter().chain(b.iter()) {
        pts.push(*lo);
        pts.push(*hi);
    }
    pts.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
    pts.dedup_by(|x, y| (*x - *y).abs() <= MERGE_EPS);
    let mut out = Vec::new();
    for w in pts.windows(2) {
        let (x, y) = (w[0], w[1]);
        if y - x <= MERGE_EPS {
            continue;
        }
        let mid = 0.5 * (x + y);
        let fa = a.iter().find(|(lo, hi, _)| *lo <= mid && mid <= *hi);
        let fb = b.iter().find(|(lo, hi, _)| *lo <= mid && mid <= *hi);
        let (ca, fa_f) = fa.map_or((false, false), |(_, _, f)| (true, *f));
        let (cb, fb_f) = fb.map_or((false, false), |(_, _, f)| (true, *f));
        let (cover, flag) = match op {
            CombineOp::Unite => (ca || cb, fa_f || fb_f),
            CombineOp::Intersect => (ca && cb, fa_f && fb_f),
            CombineOp::Subtract => (ca && !cb, fa_f),
        };
        if cover {
            out.push((x, y, flag));
        }
    }
    merge_flagged(out)
}
