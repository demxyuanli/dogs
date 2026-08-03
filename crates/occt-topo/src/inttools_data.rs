//! IntTools data classes — Phase 16 Wave A1.
//!
//! Ports the `IntTools_*` data containers from TKBO
//! (`ModelingAlgorithms/TKBO/IntTools`) into self-contained Rust value types:
//!
//! - [`IntRange`] / [`IntToolsRange`] — a 1D range `[first, last]`
//!   (`IntTools_Range`).
//! - [`CommonPrt`] / [`CommonPartType`] — a common part between shapes
//!   (`IntTools_CommonPrt`).
//! - [`IntRoot`] / [`RootType`] — a root of a 1D function, kept as the
//!   interval the root was found in (`IntTools_Root`).
//! - [`PntOnFace`] / [`PntOn2Faces`] — a 3D point plus UV parameters on one /
//!   two faces (`IntTools_PntOnFace`, `IntTools_PntOn2Faces`).
//! - [`IntCurve`] / [`CurveKind`] — an intersection-curve container
//!   (`IntTools_Curve`).
//! - [`MarkedRangeSet`] — a flagged interval set with union / intersection /
//!   difference (`IntTools_MarkedRangeSet`).
//! - [`LocalizeData`] / [`LocalizeData2`] — local edge/face intersection
//!   bookkeeping (face range, edge ranges, UV point lists).
//!
//! These are pure data classes: they carry geometry-adjacent numbers and shape
//! references but perform no surface/curve evaluation themselves. The module is
//! self-contained — it depends only on `GpPnt` (coordinates) and `TopoShape`
//! (shape handles), never on sibling Phase 16 modules.

use crate::shape::TopoShape;
use occt_core::gp::GpPnt;

/// Tolerance used when merging adjacent intervals (`MarkedRangeSet`,
/// `LocalizeData`). Two intervals are considered adjacent / touching when their
/// boundaries are within this distance.
const MERGE_EPS: f64 = 1e-9;

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
/// OCCT stores this as a `TopAbs_ShapeEnum` (`VERTEX`/`EDGE`/`FACE`); this port
/// collapses it to the three cases the exact boolean layer needs.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum CommonPartType {
    /// Uninitialized / not yet classified.
    #[default]
    Unknown,
    /// The common part is a sub-range of an edge (a vertex or coincident arc).
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
}

impl CommonPrt {
    /// Empty common part of [`CommonPartType::Unknown`].
    pub fn new() -> Self {
        Self {
            part_type: CommonPartType::Unknown,
            range: IntRange::zero(),
            face: None,
            vertices: Vec::new(),
        }
    }

    /// Fully-specified constructor.
    pub fn with(
        part_type: CommonPartType,
        range: IntRange,
        face: Option<TopoShape>,
        vertices: Vec<TopoShape>,
    ) -> Self {
        Self { part_type, range, face, vertices }
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
    ranges: Vec<(IntRange, bool)>,
}

/// A `(lo, hi, flag)` triple used by the internal interval algebra.
type Flagged = (f64, f64, bool);

/// Sort `v`, drop degenerate (≤ `MERGE_EPS`) intervals and merge adjacent
/// intervals that carry the same flag.
fn merge_flagged(mut v: Vec<Flagged>) -> Vec<Flagged> {
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
enum CombineOp {
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
fn combine(a: &[Flagged], b: &[Flagged], op: CombineOp) -> Vec<Flagged> {
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

impl MarkedRangeSet {
    pub fn new() -> Self {
        Self { ranges: Vec::new() }
    }

    pub fn len(&self) -> usize {
        self.ranges.len()
    }

    pub fn is_empty(&self) -> bool {
        self.ranges.is_empty()
    }

    /// Flatten to `(lo, hi, flag)` triples.
    fn flagged(&self) -> Vec<Flagged> {
        self.ranges.iter().map(|(r, f)| (r.first, r.last, *f)).collect()
    }

    fn from_flagged(v: Vec<Flagged>) -> Self {
        Self {
            ranges: v.into_iter().map(|(lo, hi, f)| (IntRange::new_unchecked(lo, hi), f)).collect(),
        }
    }

    /// Inserts `range` as an *unmarked* region.
    ///
    /// Any existing coverage overlapping `range` is clipped away (its parts
    /// outside `range` keep their marks); the region `range` itself is added
    /// unmarked. This matches `IntTools_MarkedRangeSet::InsertRange` when the
    /// inserted flag is "not used". Adjacent same-flag intervals are merged.
    pub fn insert(&mut self, range: IntRange) {
        let (r0, r1) = (range.first, range.last);
        let mut out: Vec<Flagged> = Vec::new();
        for (lo, hi, f) in self.flagged() {
            if hi <= r0 || lo >= r1 {
                out.push((lo, hi, f));
                continue;
            }
            if lo < r0 {
                out.push((lo, r0, f));
            }
            if hi > r1 {
                out.push((r1, hi, f));
            }
        }
        out.push((r0, r1, false));
        self.ranges = Self::from_flagged(merge_flagged(out)).ranges;
    }

    /// Marks the entire range that contains `param` (both boundaries
    /// inclusive). `Err` when `param` is not covered by any range.
    pub fn mark(&mut self, param: f64) -> Result<(), String> {
        match self.ranges.iter().position(|(r, _)| r.contains(param)) {
            Some(i) => {
                self.ranges[i].1 = true;
                Ok(())
            }
            None => Err(format!("MarkedRangeSet::mark: parameter {param} not covered")),
        }
    }

    /// Whether the range containing `param` is marked. Returns `false` when
    /// `param` is not covered by any range.
    pub fn is_marked(&self, param: f64) -> bool {
        self.ranges
            .iter()
            .find(|(r, _)| r.contains(param))
            .map_or(false, |(_, f)| *f)
    }

    /// The ranges (marked and unmarked), sorted and adjacent-merged.
    pub fn sorted_ranges(&self) -> Vec<IntRange> {
        self.ranges.iter().map(|(r, _)| *r).collect()
    }

    /// Only the marked ranges, sorted.
    pub fn marked_ranges(&self) -> Vec<IntRange> {
        self.ranges.iter().filter(|(_, f)| *f).map(|(r, _)| *r).collect()
    }

    /// Union of coverage. A point of the result is marked when it is marked in
    /// either operand.
    pub fn unite(&self, other: &MarkedRangeSet) -> MarkedRangeSet {
        Self::from_flagged(combine(&self.flagged(), &other.flagged(), CombineOp::Unite))
    }

    /// Intersection of coverage. A point of the result is marked when it is
    /// marked in both operands.
    pub fn intersect(&self, other: &MarkedRangeSet) -> MarkedRangeSet {
        Self::from_flagged(combine(&self.flagged(), &other.flagged(), CombineOp::Intersect))
    }

    /// Difference of coverage (`self` minus `other`). Marks are carried from
    /// `self` where the coverage survives.
    pub fn subtract(&self, other: &MarkedRangeSet) -> MarkedRangeSet {
        Self::from_flagged(combine(&self.flagged(), &other.flagged(), CombineOp::Subtract))
    }
}

impl Default for MarkedRangeSet {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// LocalizeData / LocalizeData2
// ---------------------------------------------------------------------------

/// Localization bookkeeping for one edge/face intersection query.
///
/// OCCT's older `IntTools_LocalizeData` recorded, for an edge–face pass, the
/// sub-range of the face under consideration, the candidate sub-ranges of the
/// edge, and the UV points computed on the face. This port keeps exactly that
/// payload. `edge_ranges` are expected to be sorted and non-overlapping; use
/// [`sort_edge_ranges`](Self::sort_edge_ranges) to normalize.
#[derive(Debug, Clone, PartialEq)]
pub struct LocalizeData {
    pub face_range: IntRange,
    pub edge_ranges: Vec<IntRange>,
    pub uv_points: Vec<(f64, f64)>,
}

impl LocalizeData {
    pub fn new(face_range: IntRange) -> Self {
        Self { face_range, edge_ranges: Vec::new(), uv_points: Vec::new() }
    }

    pub fn face_range(&self) -> IntRange {
        self.face_range
    }

    pub fn edge_ranges(&self) -> &[IntRange] {
        &self.edge_ranges
    }

    pub fn uv_points(&self) -> &[(f64, f64)] {
        &self.uv_points
    }

    pub fn add_edge_range(&mut self, r: IntRange) {
        self.edge_ranges.push(r);
    }

    pub fn add_uv_point(&mut self, u: f64, v: f64) {
        self.uv_points.push((u, v));
    }

    /// Sorts `edge_ranges` and merges adjacent/overlapping intervals, so the
    /// list is normalized for consumption by downstream splitting code.
    pub fn sort_edge_ranges(&mut self) {
        self.edge_ranges = merge_plain(std::mem::take(&mut self.edge_ranges));
    }

    /// True when the edge sub-ranges are already sorted and pairwise disjoint
    /// (with touching intervals merged).
    pub fn is_sorted(&self) -> bool {
        self.edge_ranges.windows(2).all(|w| {
            w[0].first <= w[1].first && (w[1].first - w[0].last) > MERGE_EPS
        })
    }
}

impl Default for LocalizeData {
    fn default() -> Self {
        Self::new(IntRange::zero())
    }
}

/// Same bookkeeping as [`LocalizeData`], but for a two-face pass: one UV list
/// per face, plus the per-face ranges and the shared edge sub-ranges.
#[derive(Debug, Clone, PartialEq)]
pub struct LocalizeData2 {
    pub face1_range: IntRange,
    pub face2_range: IntRange,
    pub edge_ranges: Vec<IntRange>,
    pub uv1_points: Vec<(f64, f64)>,
    pub uv2_points: Vec<(f64, f64)>,
}

impl LocalizeData2 {
    pub fn new(face1_range: IntRange, face2_range: IntRange) -> Self {
        Self { face1_range, face2_range, edge_ranges: Vec::new(), uv1_points: Vec::new(), uv2_points: Vec::new() }
    }

    pub fn face1_range(&self) -> IntRange {
        self.face1_range
    }

    pub fn face2_range(&self) -> IntRange {
        self.face2_range
    }

    pub fn edge_ranges(&self) -> &[IntRange] {
        &self.edge_ranges
    }

    pub fn uv1_points(&self) -> &[(f64, f64)] {
        &self.uv1_points
    }

    pub fn uv2_points(&self) -> &[(f64, f64)] {
        &self.uv2_points
    }

    pub fn add_edge_range(&mut self, r: IntRange) {
        self.edge_ranges.push(r);
    }

    pub fn add_uv1(&mut self, u: f64, v: f64) {
        self.uv1_points.push((u, v));
    }

    pub fn add_uv2(&mut self, u: f64, v: f64) {
        self.uv2_points.push((u, v));
    }

    /// Sorts and merges `edge_ranges` (see [`LocalizeData::sort_edge_ranges`]).
    pub fn sort_edge_ranges(&mut self) {
        self.edge_ranges = merge_plain(std::mem::take(&mut self.edge_ranges));
    }
}

impl Default for LocalizeData2 {
    fn default() -> Self {
        Self::new(IntRange::zero(), IntRange::zero())
    }
}

/// Sort a list of ranges and merge adjacent/overlapping entries.
fn merge_plain(mut v: Vec<IntRange>) -> Vec<IntRange> {
    if v.is_empty() {
        return v;
    }
    v.sort_by(|a, b| {
        a.first
            .partial_cmp(&b.first)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(a.last.partial_cmp(&b.last).unwrap_or(std::cmp::Ordering::Equal))
    });
    let mut out: Vec<IntRange> = vec![v[0]];
    for r in v.into_iter().skip(1) {
        let last = out.last_mut().unwrap();
        if r.first <= last.last + MERGE_EPS {
            last.last = last.last.max(r.last);
        } else {
            out.push(r);
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::abs::ShapeType;

    /// Test helper: a validated range.
    fn r(a: f64, b: f64) -> IntRange {
        IntRange::new(a, b).unwrap()
    }

    fn shape(t: ShapeType) -> TopoShape {
        TopoShape::new(t)
    }

    // ---- IntRange ----

    #[test]
    fn intrange_contains_inclusive_boundaries() {
        let range = r(0.0, 1.0);
        assert!(range.contains(0.0), "lower boundary");
        assert!(range.contains(1.0), "upper boundary");
        assert!(range.contains(0.5));
        assert!(!range.contains(-0.1));
        assert!(!range.contains(1.1));
        assert_eq!(range.length(), 1.0);
    }

    #[test]
    fn intrange_overlaps_touching_and_disjoint() {
        assert!(r(0.0, 1.0).overlaps(&r(1.0, 2.0)), "touching ranges overlap");
        assert!(r(0.0, 2.0).overlaps(&r(1.0, 3.0)));
        assert!(r(2.0, 3.0).overlaps(&r(0.0, 2.5)));
        assert!(!r(0.0, 1.0).overlaps(&r(1.5, 2.0)), "gap between");
        assert!(r(0.0, 1.0).disjoint(&r(1.5, 2.0)));
        assert!(!r(0.0, 1.0).disjoint(&r(1.0, 2.0)));
    }

    #[test]
    fn intrange_merge_spans() {
        assert_eq!(r(1.0, 3.0).merge(&r(0.0, 2.0)), r(0.0, 3.0));
        // Disjoint ranges merge to the bounding hull.
        assert_eq!(r(0.0, 1.0).merge(&r(2.0, 3.0)), r(0.0, 3.0));
    }

    #[test]
    fn intrange_new_rejects_inverted_and_nan() {
        assert!(IntRange::new(2.0, 1.0).is_err(), "inverted bounds rejected");
        assert!(IntRange::new(f64::NAN, 1.0).is_err());
        assert!(IntRange::new(0.0, f64::INFINITY).is_err());
        assert!(IntRange::new(0.0, 0.0).is_ok(), "degenerate point range allowed");
        assert!(IntRange::new(-1.0, 1.0).is_ok());
        assert!(r(0.0, 1.0).is_valid());
    }

    #[test]
    fn inttools_range_alias() {
        let a: IntToolsRange = r(0.0, 2.0);
        let b: IntRange = a;
        assert_eq!(b.first(), 0.0);
        assert_eq!(b.last(), 2.0);
    }

    // ---- CommonPrt ----

    #[test]
    fn commonprt_default_is_unknown_empty() {
        let cp = CommonPrt::new();
        assert_eq!(cp.part_type(), CommonPartType::Unknown);
        assert!(cp.is_empty());
        assert!(!cp.is_edge_part());
        assert!(!cp.is_face_part());
    }

    #[test]
    fn commonprt_construct_and_accessors() {
        let face = shape(ShapeType::Face);
        let v1 = shape(ShapeType::Vertex);
        let v2 = shape(ShapeType::Vertex);
        let mut cp = CommonPrt::with(CommonPartType::Edge, r(0.5, 2.5), Some(face), vec![v1, v2]);
        assert!(cp.is_edge_part());
        assert_eq!(cp.range(), r(0.5, 2.5));
        assert!(cp.face().is_some());
        assert_eq!(cp.vertices().len(), 2);

        cp.set_part_type(CommonPartType::Face);
        assert!(cp.is_face_part());
        cp.set_range(r(1.0, 1.0));
        assert_eq!(cp.range(), r(1.0, 1.0));
        cp.set_face(None);
        assert!(cp.face().is_none());
        cp.add_vertex(shape(ShapeType::Vertex));
        assert_eq!(cp.vertices().len(), 3);
    }

    // ---- IntRoot ----

    #[test]
    fn introot_construct_and_accessors() {
        let root = IntRoot::new(3, RootType::IsRoot, r(0.9, 1.1));
        assert_eq!(root.root_index(), 3);
        assert_eq!(root.root_type(), RootType::IsRoot);
        assert!(root.is_root());
        assert!(!root.is_tangent());
        assert!(!root.is_conflict());
        assert!((root.root_value() - 1.0).abs() < 1e-12);

        let tan = IntRoot::with_conflict(1, RootType::IsTangent, r(0.0, 1.0), true);
        assert!(tan.is_tangent());
        assert!(tan.is_conflict());
        assert_eq!(tan.root_index(), 1);
    }

    #[test]
    fn introot_default_unknown() {
        let root = IntRoot::default();
        assert_eq!(root.root_type(), RootType::IsUnknown);
        assert_eq!(root.root_index(), 0);
        assert_eq!(root.range(), r(0.0, 0.0));
        root_ok(root);
    }

    fn root_ok(root: IntRoot) {
        assert!(!root.is_root());
    }

    // ---- PntOnFace / PntOn2Faces ----

    #[test]
    fn pntonface_construct() {
        let pnt = GpPnt::new(1.0, 2.0, 3.0);
        let p = PntOnFace::new(7, (0.25, 0.75), pnt);
        assert_eq!(p.face_index(), 7);
        assert_eq!(p.uv(), (0.25, 0.75));
        assert_eq!(p.u(), 0.25);
        assert_eq!(p.v(), 0.75);
        assert_eq!(*p.pnt(), pnt);
        assert_eq!(PntOnFace::default().face_index(), 0);
    }

    #[test]
    fn pnton2faces_construct() {
        let p1 = GpPnt::new(1.0, 0.0, 0.0);
        let p2 = GpPnt::new(0.0, 1.0, 0.0);
        let hit = PntOn2Faces::new(0, 1, p1, p2, (0.0, 1.0), (1.0, 0.0));
        assert_eq!(hit.face1_index(), 0);
        assert_eq!(hit.face2_index(), 1);
        assert_eq!(*hit.pnt1(), p1);
        assert_eq!(*hit.pnt2(), p2);
        assert_eq!(hit.uv1(), (0.0, 1.0));
        assert_eq!(hit.uv2(), (1.0, 0.0));
    }

    // ---- IntCurve ----

    #[test]
    fn intcurve_construct_and_kinds() {
        let f1 = shape(ShapeType::Face);
        let f2 = shape(ShapeType::Face);
        let curve = IntCurve::new(CurveKind::Circle, Some(f1), Some(f2), r(-1.0, 1.0));
        assert_eq!(curve.kind(), CurveKind::Circle);
        assert!(curve.face1().is_some());
        assert!(curve.face2().is_some());
        assert_eq!(curve.range(), r(-1.0, 1.0));

        let line = IntCurve::new(CurveKind::Line, None, None, r(0.0, 2.0));
        assert_eq!(line.kind(), CurveKind::Line);
        assert!(line.face1().is_none());

        let mut b = IntCurve::default();
        assert_eq!(b.kind(), CurveKind::Other);
        b.set_range(r(0.0, 5.0));
        assert_eq!(b.range(), r(0.0, 5.0));
        // All curve kinds are representable.
        for k in [CurveKind::Ellipse, CurveKind::Parabola, CurveKind::Hyperbola, CurveKind::BSpline] {
            assert_eq!(IntCurve::new(k, None, None, r(0.0, 1.0)).kind(), k);
        }
    }

    // ---- MarkedRangeSet ----

    #[test]
    fn marked_insert_disjoint_sorted() {
        let mut s = MarkedRangeSet::new();
        s.insert(r(4.0, 6.0));
        s.insert(r(0.0, 2.0));
        assert_eq!(s.sorted_ranges(), vec![r(0.0, 2.0), r(4.0, 6.0)]);
        assert_eq!(s.len(), 2);
        assert!(s.is_marked(1.0) == false);
    }

    #[test]
    fn marked_insert_overlap_merges() {
        let mut s = MarkedRangeSet::new();
        s.insert(r(0.0, 5.0));
        s.insert(r(2.0, 3.0)); // interior, same (unmarked) flag -> merged
        assert_eq!(s.sorted_ranges(), vec![r(0.0, 5.0)]);
        assert_eq!(s.len(), 1);
    }

    #[test]
    fn marked_mark_whole_interval() {
        let mut s = MarkedRangeSet::new();
        s.insert(r(0.0, 1.0));
        s.insert(r(3.0, 4.0));
        s.mark(0.5).unwrap();
        assert!(s.is_marked(0.5));
        assert!(!s.is_marked(3.5), "other interval stays unmarked");
        assert_eq!(s.marked_ranges(), vec![r(0.0, 1.0)]);
    }

    #[test]
    fn marked_is_marked_uncovered_false() {
        let mut s = MarkedRangeSet::new();
        s.insert(r(0.0, 1.0));
        assert!(!s.is_marked(10.0));
        assert!(!s.is_marked(-1.0));
    }

    #[test]
    fn marked_mark_outside_errors() {
        let mut s = MarkedRangeSet::new();
        s.insert(r(0.0, 1.0));
        assert!(s.mark(2.0).is_err());
        assert!(s.mark(f64::NAN).is_err());
    }

    #[test]
    fn marked_insert_inside_marked_unmarks_overlap() {
        let mut s = MarkedRangeSet::new();
        s.insert(r(0.0, 5.0));
        s.mark(2.0).unwrap();
        assert!(s.is_marked(4.0));
        // Inserting an unmarked interior range splits and unmarks the overlap.
        s.insert(r(1.0, 3.0));
        assert!(s.is_marked(0.5), "left fragment keeps mark");
        assert!(!s.is_marked(2.0), "overlap becomes unmarked");
        assert!(s.is_marked(4.0), "right fragment keeps mark");
        assert_eq!(s.sorted_ranges(), vec![r(0.0, 1.0), r(1.0, 3.0), r(3.0, 5.0)]);
        assert_eq!(s.marked_ranges(), vec![r(0.0, 1.0), r(3.0, 5.0)]);
    }

    #[test]
    fn marked_unite() {
        let mut a = MarkedRangeSet::new();
        a.insert(r(0.0, 2.0));
        a.insert(r(4.0, 6.0));
        a.mark(1.0).unwrap();
        let mut b = MarkedRangeSet::new();
        b.insert(r(1.0, 5.0));
        b.mark(3.0).unwrap();

        let u = a.unite(&b);
        assert_eq!(u.sorted_ranges(), vec![r(0.0, 5.0), r(5.0, 6.0)]);
        assert!(u.is_marked(0.5), "covered only by marked A");
        assert!(u.is_marked(2.5), "covered only by marked B");
        assert!(!u.is_marked(5.5), "covered by unmarked A only");
    }

    #[test]
    fn marked_intersect() {
        let mut a = MarkedRangeSet::new();
        a.insert(r(0.0, 2.0));
        a.insert(r(4.0, 6.0));
        a.mark(1.0).unwrap();
        let mut b = MarkedRangeSet::new();
        b.insert(r(1.0, 5.0));
        b.mark(3.0).unwrap();

        let i = a.intersect(&b);
        assert_eq!(i.sorted_ranges(), vec![r(1.0, 2.0), r(4.0, 5.0)]);
        assert!(i.is_marked(1.5), "marked in both");
        assert!(!i.is_marked(4.5), "marked only in B -> not marked in intersection");
    }

    #[test]
    fn marked_subtract() {
        let mut a = MarkedRangeSet::new();
        a.insert(r(0.0, 2.0));
        a.insert(r(4.0, 6.0));
        a.mark(1.0).unwrap();
        let mut b = MarkedRangeSet::new();
        b.insert(r(1.0, 5.0));
        b.mark(3.0).unwrap();

        let d = a.subtract(&b);
        assert_eq!(d.sorted_ranges(), vec![r(0.0, 1.0), r(5.0, 6.0)]);
        assert!(d.is_marked(0.5), "surviving marked part of A");
        assert!(!d.is_marked(5.5), "surviving unmarked part of A");
    }

    #[test]
    fn marked_unite_with_empty() {
        let mut a = MarkedRangeSet::new();
        a.insert(r(0.0, 2.0));
        a.mark(1.0).unwrap();
        let empty = MarkedRangeSet::new();
        let u = a.unite(&empty);
        assert_eq!(u.sorted_ranges(), a.sorted_ranges());
        assert!(u.is_marked(1.0));
        assert!(empty.is_empty());
        assert!(a.subtract(&a).is_empty(), "self-subtraction clears everything");
    }

    // ---- LocalizeData ----

    #[test]
    fn localizedata_construct_and_uv() {
        let mut ld = LocalizeData::new(r(0.0, 3.0));
        assert_eq!(ld.face_range(), r(0.0, 3.0));
        ld.add_edge_range(r(1.0, 2.0));
        ld.add_uv_point(0.5, 1.5);
        assert_eq!(ld.edge_ranges(), &[r(1.0, 2.0)]);
        assert_eq!(ld.uv_points(), &[(0.5, 1.5)]);
        assert_eq!(LocalizeData::default().face_range(), r(0.0, 0.0));
    }

    #[test]
    fn localizedata_sort_edge_ranges_merges() {
        let mut ld = LocalizeData::new(r(0.0, 10.0));
        ld.add_edge_range(r(6.0, 8.0));
        ld.add_edge_range(r(0.0, 2.0));
        ld.add_edge_range(r(2.0, 4.0)); // adjacent to [0,2] -> merged
        ld.add_edge_range(r(9.0, 10.0));
        assert!(!ld.is_sorted());
        ld.sort_edge_ranges();
        assert_eq!(ld.edge_ranges(), &[r(0.0, 4.0), r(6.0, 8.0), r(9.0, 10.0)]);
        assert!(ld.is_sorted());
    }

    #[test]
    fn localizedata2_construct() {
        let mut ld2 = LocalizeData2::new(r(0.0, 1.0), r(2.0, 3.0));
        ld2.add_edge_range(r(0.5, 1.5));
        ld2.add_uv1(0.1, 0.2);
        ld2.add_uv2(0.3, 0.4);
        assert_eq!(ld2.face1_range(), r(0.0, 1.0));
        assert_eq!(ld2.face2_range(), r(2.0, 3.0));
        assert_eq!(ld2.edge_ranges(), &[r(0.5, 1.5)]);
        assert_eq!(ld2.uv1_points(), &[(0.1, 0.2)]);
        assert_eq!(ld2.uv2_points(), &[(0.3, 0.4)]);
        ld2.sort_edge_ranges();
        assert_eq!(ld2.edge_ranges(), &[r(0.5, 1.5)]);
    }

    #[test]
    fn merge_plain_dedup_and_touch() {
        let v = vec![r(2.0, 3.0), r(0.0, 1.0), r(1.0, 2.0), r(2.0, 3.0)];
        assert_eq!(merge_plain(v), vec![r(0.0, 3.0)]);
    }
}
