
use super::*;

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
    pub(super) fn flagged(&self) -> Vec<Flagged> {
        self.ranges.iter().map(|(r, f)| (r.first, r.last, *f)).collect()
    }

    pub(super) fn from_flagged(v: Vec<Flagged>) -> Self {
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
pub(super) fn merge_plain(mut v: Vec<IntRange>) -> Vec<IntRange> {
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
