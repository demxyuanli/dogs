use super::prelude::*;


/// Deflection used to derive the angular sampling step of analytic surfaces.
/// Source: `IntTools_TopolTool::ComputeSamplePoints` (1.e-02).

pub(super) const SAMPLE_DEFLECTION: f64 = 1e-2;
/// Fallback span for unbounded parametric directions. Source: the 1.e5 sentinel
/// in `IntTools_TopolTool::ComputeSamplePoints`.
pub(super) const BIG_RANGE: f64 = 1e5;
/// Cap on the number of samples per direction. Source:
/// `aMaxNbSample = 50` in `IntTools_TopolTool::ComputeSamplePoints`.
pub(super) const MAX_NB_SAMPLE: usize = 50;

// ---------------------------------------------------------------------------
// BaseRangeSample
// ---------------------------------------------------------------------------

/// Base class for range-index management.
///
/// A range sample is identified by a non-negative `depth`: depth 0 denotes the
/// whole parameter domain; a deeper sample at depth `d` addresses one of
/// `nb_sample^d` equal sub-intervals. Mirrors `IntTools_BaseRangeSample`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct BaseRangeSample {
    pub(super) depth: usize,
}

impl BaseRangeSample {
    /// Empty sample at depth 0 (the whole domain).
    pub fn new() -> Self {
        Self { depth: 0 }
    }

    /// Sample at the given depth.
    pub fn with_depth(depth: usize) -> Self {
        Self { depth }
    }

    /// The subdivision depth of this sample.
    pub fn get_depth(&self) -> usize {
        self.depth
    }

    /// Sets the subdivision depth.
    pub fn set_depth(&mut self, depth: usize) {
        self.depth = depth;
    }
}

// ---------------------------------------------------------------------------
// CurveRangeSample
// ---------------------------------------------------------------------------

/// A curve parameter sub-range addressed by `(depth, index)`.
///
/// At a given [`depth`](Self::get_depth), `index` selects one of the
/// `nb_sample^depth` equal sub-intervals into which the curve parameter domain
/// `[first, last]` is split by [`get_range`](Self::get_range). Depth 0 means
/// the sample covers the whole domain.
///
/// Mirrors `IntTools_CurveRangeSample`. OCCT does **not** store a materialised
/// range — `GetRange(theFirst, theLast, theNbSample)` computes it from the
/// depth and index, so this port keeps the same on-demand `get_range`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CurveRangeSample {
    pub(super) depth: usize,
    pub(super) index: usize,
}

impl CurveRangeSample {
    /// Sample covering the whole domain: depth 0, index 0.
    pub fn new() -> Self {
        Self { depth: 0, index: 0 }
    }

    /// Sample at depth 0 with the given index.
    pub fn with_index(index: usize) -> Self {
        Self { depth: 0, index }
    }

    /// Sample with the given depth and index.
    pub fn with_index_depth(index: usize, depth: usize) -> Self {
        Self { depth, index }
    }

    /// The subdivision depth of this sample.
    pub fn get_depth(&self) -> usize {
        self.depth
    }

    /// Sets the subdivision depth.
    pub fn set_depth(&mut self, depth: usize) {
        self.depth = depth;
    }

    /// The sub-interval index within the current depth.
    pub fn get_index(&self) -> usize {
        self.index
    }

    /// Sets the sub-interval index.
    pub fn set_index(&mut self, index: usize) {
        self.index = index;
    }

    /// Equal when both depth and index match. Source: `IsEqual`.
    pub fn is_equal(&self, other: &Self) -> bool {
        self.index == other.index && self.depth == other.depth
    }

    /// The materialised `[first, last]` sub-range of the curve.
    ///
    /// Depth 0 returns the whole `[first, last]`; otherwise the domain is
    /// split into `nb_sample^depth` equal intervals and this sample's `index`
    /// selects one. Source: `IntTools_CurveRangeSample::GetRange`.
    pub fn get_range(&self, first: f64, last: f64, nb_sample: usize) -> IntRange {
        if self.depth == 0 {
            return IntRange::new_unchecked(first, last);
        }
        let tmp = (nb_sample as f64).powi(self.depth as i32);
        let local_diff = (last - first) / tmp;
        let sample_first = first + self.index as f64 * local_diff;
        IntRange::new_unchecked(sample_first, sample_first + local_diff)
    }

    /// The index of the deeper sub-interval `[index * nb_sample,
    /// (index + 1) * nb_sample)` at depth + 1. Source:
    /// `GetRangeIndexDeeper`.
    pub fn get_range_index_deeper(&self, nb_sample: usize) -> usize {
        self.index * nb_sample
    }
}

impl Default for CurveRangeSample {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// SurfaceRangeSample
// ---------------------------------------------------------------------------

/// A surface parameter cell addressed by `(depth_u, index_u, depth_v, index_v)`.
///
/// The `U` and `V` directions are each a [`CurveRangeSample`]; the cell bounds
/// are materialised on demand by [`get_range_u`](Self::get_range_u) /
/// [`get_range_v`](Self::get_range_v). Mirrors `IntTools_SurfaceRangeSample`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SurfaceRangeSample {
    pub(super) range_u: CurveRangeSample,
    pub(super) range_v: CurveRangeSample,
}

impl SurfaceRangeSample {
    /// Empty cell: depth 0 / index 0 in both directions.
    pub fn new() -> Self {
        Self { range_u: CurveRangeSample::new(), range_v: CurveRangeSample::new() }
    }

    /// Cell from explicit U/V indexes and depths.
    pub fn with_indexes_depths(
        index_u: usize,
        depth_u: usize,
        index_v: usize,
        depth_v: usize,
    ) -> Self {
        Self {
            range_u: CurveRangeSample::with_index_depth(index_u, depth_u),
            range_v: CurveRangeSample::with_index_depth(index_v, depth_v),
        }
    }

    /// Cell from two [`CurveRangeSample`]s (one per direction). Source:
    /// the `(theRangeU, theRangeV)` constructor.
    pub fn from_ranges(range_u: CurveRangeSample, range_v: CurveRangeSample) -> Self {
        Self { range_u, range_v }
    }

    /// Sets both direction samples at once. Source: `SetRanges`.
    pub fn set_ranges(&mut self, range_u: CurveRangeSample, range_v: CurveRangeSample) {
        self.range_u = range_u;
        self.range_v = range_v;
    }

    /// Returns the two direction samples. Source: `GetRanges`.
    pub fn get_ranges(&self) -> (CurveRangeSample, CurveRangeSample) {
        (self.range_u, self.range_v)
    }

    /// Sets the U/V indexes (depths unchanged). Source: `SetIndexes`.
    pub fn set_indexes(&mut self, index_u: usize, index_v: usize) {
        self.range_u.set_index(index_u);
        self.range_v.set_index(index_v);
    }

    /// Returns the U/V indexes. Source: `GetIndexes`.
    pub fn get_indexes(&self) -> (usize, usize) {
        (self.range_u.get_index(), self.range_v.get_index())
    }

    /// Returns the U/V depths. Source: `GetDepths`.
    pub fn get_depths(&self) -> (usize, usize) {
        (self.range_u.get_depth(), self.range_v.get_depth())
    }

    /// Sets the U-direction sample. Source: `SetSampleRangeU`.
    pub fn set_sample_range_u(&mut self, range_u: CurveRangeSample) {
        self.range_u = range_u;
    }

    /// Returns the U-direction sample. Source: `GetSampleRangeU`.
    pub fn get_sample_range_u(&self) -> CurveRangeSample {
        self.range_u
    }

    /// Sets the V-direction sample. Source: `SetSampleRangeV`.
    pub fn set_sample_range_v(&mut self, range_v: CurveRangeSample) {
        self.range_v = range_v;
    }

    /// Returns the V-direction sample. Source: `GetSampleRangeV`.
    pub fn get_sample_range_v(&self) -> CurveRangeSample {
        self.range_v
    }

    /// Sets the U index. Source: `SetIndexU`.
    pub fn set_index_u(&mut self, index_u: usize) {
        self.range_u.set_index(index_u);
    }

    /// Returns the U index. Source: `GetIndexU`.
    pub fn get_index_u(&self) -> usize {
        self.range_u.get_index()
    }

    /// Sets the V index. Source: `SetIndexV`.
    pub fn set_index_v(&mut self, index_v: usize) {
        self.range_v.set_index(index_v);
    }

    /// Returns the V index. Source: `GetIndexV`.
    pub fn get_index_v(&self) -> usize {
        self.range_v.get_index()
    }

    /// Sets the U depth. Source: `SetDepthU`.
    pub fn set_depth_u(&mut self, depth_u: usize) {
        self.range_u.set_depth(depth_u);
    }

    /// Returns the U depth. Source: `GetDepthU`.
    pub fn get_depth_u(&self) -> usize {
        self.range_u.get_depth()
    }

    /// Sets the V depth. Source: `SetDepthV`.
    pub fn set_depth_v(&mut self, depth_v: usize) {
        self.range_v.set_depth(depth_v);
    }

    /// Returns the V depth. Source: `GetDepthV`.
    pub fn get_depth_v(&self) -> usize {
        self.range_v.get_depth()
    }

    /// The materialised `[first_u, last_u]` sub-range. Source: `GetRangeU`.
    pub fn get_range_u(&self, first_u: f64, last_u: f64, nb_sample_u: usize) -> IntRange {
        self.range_u.get_range(first_u, last_u, nb_sample_u)
    }

    /// The materialised `[first_v, last_v]` sub-range. Source: `GetRangeV`.
    pub fn get_range_v(&self, first_v: f64, last_v: f64, nb_sample_v: usize) -> IntRange {
        self.range_v.get_range(first_v, last_v, nb_sample_v)
    }

    /// The `(U, V)` cell bounds materialised from both domains at once.
    /// Equivalent to calling [`get_range_u`](Self::get_range_u) and
    /// [`get_range_v`](Self::get_range_v).
    pub fn get_range(
        &self,
        first_u: f64,
        last_u: f64,
        nb_sample_u: usize,
        first_v: f64,
        last_v: f64,
        nb_sample_v: usize,
    ) -> (IntRange, IntRange) {
        (
            self.range_u.get_range(first_u, last_u, nb_sample_u),
            self.range_v.get_range(first_v, last_v, nb_sample_v),
        )
    }

    /// Equal when both direction samples are equal. Source: `IsEqual`.
    pub fn is_equal(&self, other: &Self) -> bool {
        self.range_u.is_equal(&other.range_u) && self.range_v.is_equal(&other.range_v)
    }

    /// The deeper U index at depth + 1. Source: `GetRangeIndexUDeeper`.
    pub fn get_range_index_u_deeper(&self, nb_sample_u: usize) -> usize {
        self.range_u.get_range_index_deeper(nb_sample_u)
    }

    /// The deeper V index at depth + 1. Source: `GetRangeIndexVDeeper`.
    pub fn get_range_index_v_deeper(&self, nb_sample_v: usize) -> usize {
        self.range_v.get_range_index_deeper(nb_sample_v)
    }
}

impl Default for SurfaceRangeSample {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// CurveRangeLocalizeData
// ---------------------------------------------------------------------------

/// Localization bookkeeping for one curve parameter domain.
///
/// Splits a curve parameter domain `[first, last]` into `nb_samples` equal
/// sub-ranges and maps each sub-range to a curve index — the range→curve
/// mapping the intersection localizer walks. Cells that are already known to
/// be outside the intersection region are tracked as "out" indices.
///
/// Mirrors `IntTools_CurveRangeLocalizeData`. OCCT stores `(CurveRangeSample →
/// Bnd_Box)` and out-range sets keyed by sample; this port keeps the same idea
/// with integer cell indices over a flattened [`ranges`](Self::ranges) list.
#[derive(Debug, Clone, PartialEq)]
pub struct CurveRangeLocalizeData {
    /// Seed curve index used when building the range mapping.
    pub root_index: usize,
    /// Number of sub-intervals the domain is split into.
    pub(super) nb_samples: usize,
    /// Minimal acceptable range length (informational; the caller checks it).
    pub(super) min_range: f64,
    /// `(sub-range, curve index)` cells covering the domain.
    pub(super) ranges: Vec<(IntRange, usize)>,
    /// Indices of cells already known to be outside the intersection region.
    pub(super) out_indices: Vec<usize>,
}

impl CurveRangeLocalizeData {
    /// Empty localizer. Source: the `(theNbSample, theMinRange)` constructor.
    pub fn new(nb_samples: usize, min_range: f64) -> Self {
        Self {
            root_index: 0,
            nb_samples,
            min_range,
            ranges: Vec::new(),
            out_indices: Vec::new(),
        }
    }

    /// The configured number of sub-intervals.
    pub fn get_nb_sample(&self) -> usize {
        self.nb_samples
    }

    /// The configured minimal range length.
    pub fn get_min_range(&self) -> f64 {
        self.min_range
    }

    /// The seed curve index. Source: OCCT's `myRootIndex`-style seed field.
    pub fn get_root_index(&self) -> usize {
        self.root_index
    }

    /// Sets the seed curve index.
    pub fn set_root_index(&mut self, root_index: usize) {
        self.root_index = root_index;
    }

    /// Splits `domain` into `nb_samples` equal sub-ranges, mapping every cell
    /// to `curve_index`. Replaces any previous mapping; returns the number of
    /// cells created. The curve parameter domain is covered by adjacent,
    /// non-overlapping intervals, so a parameter can be located with
    /// [`find_index`](Self::find_index).
    pub fn build(&mut self, domain: IntRange, curve_index: usize) -> usize {
        let n = self.nb_samples.max(1);
        self.ranges.clear();
        self.out_indices.clear();
        if !domain.is_valid() {
            return 0;
        }
        let step = domain.length() / n as f64;
        for i in 0..n {
            let first = domain.first + i as f64 * step;
            let last = if i + 1 == n { domain.last } else { first + step };
            self.ranges.push((IntRange::new_unchecked(first, last), curve_index));
        }
        self.ranges.len()
    }

    /// All `(sub-range, curve index)` cells.
    pub fn ranges(&self) -> &[(IntRange, usize)] {
        &self.ranges
    }

    /// The sub-range at cell `index`.
    pub fn range(&self, index: usize) -> Option<IntRange> {
        self.ranges.get(index).map(|(r, _)| *r)
    }

    /// The curve index mapped at cell `index`.
    pub fn curve_index(&self, index: usize) -> Option<usize> {
        self.ranges.get(index).map(|(_, c)| *c)
    }

    /// Re-maps cell `index` to `curve_index`. `Err` when `index` is out of
    /// range.
    pub fn set_curve_index(&mut self, index: usize, curve_index: usize) -> Result<(), String> {
        let cell = self
            .ranges
            .get_mut(index)
            .ok_or_else(|| format!("CurveRangeLocalizeData::set_curve_index: index {index} out of range"))?;
        cell.1 = curve_index;
        Ok(())
    }

    /// Locates the cell containing parameter `t` (both bounds inclusive),
    /// returning its cell index.
    pub fn find_index(&self, t: f64) -> Option<usize> {
        self.ranges
            .iter()
            .position(|(r, _)| r.contains(t))
    }

    /// Marks cell `index` as out. `Err` when `index` is out of range.
    pub fn add_out_range(&mut self, index: usize) -> Result<(), String> {
        if index >= self.ranges.len() {
            return Err(format!(
                "CurveRangeLocalizeData::add_out_range: index {index} out of range"
            ));
        }
        if !self.out_indices.contains(&index) {
            self.out_indices.push(index);
        }
        Ok(())
    }

    /// Whether cell `index` is known out.
    pub fn is_range_out(&self, index: usize) -> bool {
        self.out_indices.contains(&index)
    }

    /// The indices of all out cells.
    pub fn list_range_out(&self) -> &[usize] {
        &self.out_indices
    }

    /// Clears the out-cell set.
    pub fn remove_range_out_all(&mut self) {
        self.out_indices.clear();
    }
}

impl Default for CurveRangeLocalizeData {
    fn default() -> Self {
        Self::new(1, 0.0)
    }
}

// ---------------------------------------------------------------------------
// SurfaceRangeLocalizeData
// ---------------------------------------------------------------------------

/// Localization bookkeeping for a `(U, V)` surface parameter grid.
///
/// Splits a `U × V` domain into `nb_samples_u × nb_samples_v` cells and maps
/// each cell to a surface index. Cells known to be outside the intersection
/// region are tracked as out. The [`get_depth_range`](Self::get_depth_range)
/// helper materialises the `(U, V)` bounds of one cell.
///
/// Mirrors `IntTools_SurfaceRangeLocalizeData`. OCCT also keeps an optimised
/// grid-point/frame layer (U/V parameter arrays, frame sub-range); this port
/// keeps the essential cell-mapping layer the localizer walks.
#[derive(Debug, Clone, PartialEq)]
pub struct SurfaceRangeLocalizeData {
    /// Seed surface index used when building the grid mapping.
    pub root_index: usize,
    /// Number of U sub-intervals.
    pub(super) nb_samples_u: usize,
    /// Number of V sub-intervals.
    pub(super) nb_samples_v: usize,
    /// Minimal acceptable U range length.
    pub(super) min_range_u: f64,
    /// Minimal acceptable V range length.
    pub(super) min_range_v: f64,
    /// `(U sub-range, V sub-range, surface index)` cells, row-major
    /// (`iv * nb_samples_u + iu`).
    pub(super) ranges: Vec<(IntRange, IntRange, usize)>,
    /// Flattened cell indices known to be out.
    pub(super) out_indices: Vec<usize>,
}

impl SurfaceRangeLocalizeData {
    /// Empty localizer. Source: the `(theNbSampleU, theNbSampleV,
    /// theMinRangeU, theMinRangeV)` constructor.
    pub fn new(nb_samples_u: usize, nb_samples_v: usize, min_range_u: f64, min_range_v: f64) -> Self {
        Self {
            root_index: 0,
            nb_samples_u,
            nb_samples_v,
            min_range_u,
            min_range_v,
            ranges: Vec::new(),
            out_indices: Vec::new(),
        }
    }

    /// Number of U sub-intervals.
    pub fn get_nb_samples_u(&self) -> usize {
        self.nb_samples_u
    }

    /// Number of V sub-intervals.
    pub fn get_nb_samples_v(&self) -> usize {
        self.nb_samples_v
    }

    /// Minimal acceptable U range length.
    pub fn get_min_range_u(&self) -> f64 {
        self.min_range_u
    }

    /// Minimal acceptable V range length.
    pub fn get_min_range_v(&self) -> f64 {
        self.min_range_v
    }

    /// The seed surface index.
    pub fn get_root_index(&self) -> usize {
        self.root_index
    }

    /// Sets the seed surface index.
    pub fn set_root_index(&mut self, root_index: usize) {
        self.root_index = root_index;
    }

    /// Splits `u_domain × v_domain` into the configured grid, mapping every
    /// cell to `surface_index`. Replaces any previous mapping; returns the
    /// number of cells (`nb_samples_u * nb_samples_v`).
    pub fn build(&mut self, u_domain: IntRange, v_domain: IntRange, surface_index: usize) -> usize {
        let nu = self.nb_samples_u.max(1);
        let nv = self.nb_samples_v.max(1);
        self.ranges.clear();
        self.out_indices.clear();
        if !u_domain.is_valid() || !v_domain.is_valid() {
            return 0;
        }
        let du = u_domain.length() / nu as f64;
        let dv = v_domain.length() / nv as f64;
        for j in 0..nv {
            let v_first = v_domain.first + j as f64 * dv;
            let v_last = if j + 1 == nv { v_domain.last } else { v_first + dv };
            for i in 0..nu {
                let u_first = u_domain.first + i as f64 * du;
                let u_last = if i + 1 == nu { u_domain.last } else { u_first + du };
                self.ranges.push((
                    IntRange::new_unchecked(u_first, u_last),
                    IntRange::new_unchecked(v_first, v_last),
                    surface_index,
                ));
            }
        }
        self.ranges.len()
    }

    /// All `(U sub-range, V sub-range, surface index)` cells, row-major.
    pub fn ranges(&self) -> &[(IntRange, IntRange, usize)] {
        &self.ranges
    }

    /// The `(U, V)` bounds of the cell at grid position `(iu, iv)`.
    pub fn get_depth_range(&self, iu: usize, iv: usize) -> Option<(IntRange, IntRange)> {
        if iu >= self.nb_samples_u.max(1) || iv >= self.nb_samples_v.max(1) {
            return None;
        }
        self.ranges
            .get(iv * self.nb_samples_u.max(1) + iu)
            .map(|(u, v, _)| (*u, *v))
    }

    /// The surface index mapped at grid position `(iu, iv)`.
    pub fn surface_index(&self, iu: usize, iv: usize) -> Option<usize> {
        if iu >= self.nb_samples_u.max(1) || iv >= self.nb_samples_v.max(1) {
            return None;
        }
        self.ranges.get(iv * self.nb_samples_u.max(1) + iu).map(|(_, _, s)| *s)
    }

    /// Locates the cell containing parameter `(u, v)`, returning its grid
    /// position `(iu, iv)`.
    pub fn find_cell(&self, u: f64, v: f64) -> Option<(usize, usize)> {
        for iv in 0..self.nb_samples_v.max(1) {
            for iu in 0..self.nb_samples_u.max(1) {
                if let Some((ur, vr, _)) = self.ranges.get(iv * self.nb_samples_u.max(1) + iu) {
                    if ur.contains(u) && vr.contains(v) {
                        return Some((iu, iv));
                    }
                }
            }
        }
        None
    }

    /// Marks the cell at `(iu, iv)` as out.
    pub fn add_out_range(&mut self, iu: usize, iv: usize) -> Result<(), String> {
        if iu >= self.nb_samples_u.max(1) || iv >= self.nb_samples_v.max(1) {
            return Err(format!(
                "SurfaceRangeLocalizeData::add_out_range: cell ({iu}, {iv}) out of grid"
            ));
        }
        let key = iv * self.nb_samples_u.max(1) + iu;
        if !self.out_indices.contains(&key) {
            self.out_indices.push(key);
        }
        Ok(())
    }

    /// Whether the cell at `(iu, iv)` is known out.
    pub fn is_range_out(&self, iu: usize, iv: usize) -> bool {
        let key = iv * self.nb_samples_u.max(1) + iu;
        self.out_indices.contains(&key)
    }

    /// Clears the out-cell set.
    pub fn remove_range_out_all(&mut self) {
        self.out_indices.clear();
    }
}

impl Default for SurfaceRangeLocalizeData {
    fn default() -> Self {
        Self::new(1, 1, 0.0, 0.0)
    }
}

// ---------------------------------------------------------------------------
// TopolTool
// ---------------------------------------------------------------------------

/// Sample-point generator for intersection algorithms.
///
/// Given a [`Surface`], computes a uniform `U × V` grid of sample parameters
/// (`compute_sample_points`, ported from `IntTools_TopolTool`), evaluates
/// individual samples (`sample_point`, returning the UV pair plus the 3D
/// surface point), and builds deflection-adaptive grids for BSpline surfaces
/// (`sample_pnts`, ported from `Adaptor3d_TopolTool::SamplePnts` semantics).
///
/// Unbounded parametric directions are clamped to `[-1e5, 1e5]` (or
/// `[last - 2e5, last]` / `[first, first + 2e5]` when only one end is
/// unbounded), exactly as OCCT's `ComputeSamplePoints` does.
///
/// `ponytail:` `dyn Surface` implements neither `Debug` nor `Clone`, so the
/// struct is not `Clone`/`Debug`; surfaces are re-attached with
/// [`initialize`](Self::initialize).
pub struct TopolTool {
    pub(super) surface: Option<Box<dyn Surface>>,
    /// Computed sample count along U.
    pub(super) nb_samples_u: usize,
    /// Computed sample count along V.
    pub(super) nb_samples_v: usize,
    /// Lower U bound of the (clamped) parameter domain.
    pub(super) u0: f64,
    /// Lower V bound of the (clamped) parameter domain.
    pub(super) v0: f64,
    /// U step between interior samples: `(u_last - u0) / (nb_samples_u + 1)`.
    pub(super) du: f64,
    /// V step between interior samples: `(v_last - v0) / (nb_samples_v + 1)`.
    pub(super) dv: f64,
    /// Clamped `[u_first, u_last]` parameter domain.
    pub(super) u_range: (f64, f64),
    /// Clamped `[v_first, v_last]` parameter domain.
    pub(super) v_range: (f64, f64),
    /// U sample parameters produced by the last [`sample_pnts`](Self::sample_pnts)
    /// call (adaptive for BSpline, uniform otherwise). `None` before that.
    pub(super) u_pars: Option<Vec<f64>>,
    /// V sample parameters produced by the last [`sample_pnts`](Self::sample_pnts)
    /// call.
    pub(super) v_pars: Option<Vec<f64>>,
}
