//! IntTools interval sampling and topology tools — Phase 16b.
//!
//! Ports the TKBO `IntTools` sample/range-localization classes into
//! self-contained Rust value types:
//!
//! - [`BaseRangeSample`] — base class for range-index management
//!   (`IntTools_BaseRangeSample`).
//! - [`CurveRangeSample`] — a curve parameter sub-range addressed by
//!   `(depth, index)`; `get_range` materialises the `[first, last]` bounds
//!   (`IntTools_CurveRangeSample`).
//! - [`SurfaceRangeSample`] — the 2D analogue: a `(U, V)` cell addressed by
//!   `(depth_u, index_u, depth_v, index_v)` (`IntTools_SurfaceRangeSample`).
//! - [`CurveRangeLocalizeData`] / [`SurfaceRangeLocalizeData`] — split a
//!   curve/UV domain by sample points and map each cell to a curve/surface
//!   index, tracking which cells are already known "out"
//!   (`IntTools_CurveRangeLocalizeData`, `IntTools_SurfaceRangeLocalizeData`).
//! - [`TopolTool`] — the sample-point generator for intersection algorithms:
//!   computes a uniform `U × V` grid (`ComputeSamplePoints`), answers
//!   `sample_point` (UV + evaluated 3D point), and builds an
//!   deflection-adaptive grid for BSpline surfaces (`SamplePnts`)
//!   (`IntTools_TopolTool`).
//!
//! The range samples reuse [`IntRange`](crate::inttools_data::IntRange) from
//! the Phase 16a data port. The module depends only on `occt_geom::Surface`
//! (via `dyn Surface`) and the surface-type classifier reused from the BRepMesh
//! range-splitter port — never on sibling Phase 16 modules.
//!
//! `ponytail:` `IntTools_TopolTool` in OCCT reads the concrete `Geom_*` surface
//! type, radii and pole/knot counts through `Adaptor3d_Surface` accessors. The
//! Rust port only sees `dyn Surface`, so the analytic type is recovered with
//! [`classify_surface`], radii are measured by sampling diametrically opposite
//! iso-parameter points, and BSpline/Bezier pole counts (unreachable through the
//! trait object) fall back to a fixed 10×10 base grid that `sample_pnts` then
//! refines adaptively.

use occt_core::gp::{GpPnt, GpPnt2d};
use occt_core::precision::ANGULAR;
use occt_geom::Surface;

use crate::inttools_data::IntRange;
use crate::meshing::range_splitter::{classify_surface, SurfaceType};

/// Deflection used to derive the angular sampling step of analytic surfaces.
/// Source: `IntTools_TopolTool::ComputeSamplePoints` (1.e-02).
const SAMPLE_DEFLECTION: f64 = 1e-2;
/// Fallback span for unbounded parametric directions. Source: the 1.e5 sentinel
/// in `IntTools_TopolTool::ComputeSamplePoints`.
const BIG_RANGE: f64 = 1e5;
/// Cap on the number of samples per direction. Source:
/// `aMaxNbSample = 50` in `IntTools_TopolTool::ComputeSamplePoints`.
const MAX_NB_SAMPLE: usize = 50;

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
    depth: usize,
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
    depth: usize,
    index: usize,
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
    range_u: CurveRangeSample,
    range_v: CurveRangeSample,
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
    nb_samples: usize,
    /// Minimal acceptable range length (informational; the caller checks it).
    min_range: f64,
    /// `(sub-range, curve index)` cells covering the domain.
    ranges: Vec<(IntRange, usize)>,
    /// Indices of cells already known to be outside the intersection region.
    out_indices: Vec<usize>,
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
    nb_samples_u: usize,
    /// Number of V sub-intervals.
    nb_samples_v: usize,
    /// Minimal acceptable U range length.
    min_range_u: f64,
    /// Minimal acceptable V range length.
    min_range_v: f64,
    /// `(U sub-range, V sub-range, surface index)` cells, row-major
    /// (`iv * nb_samples_u + iu`).
    ranges: Vec<(IntRange, IntRange, usize)>,
    /// Flattened cell indices known to be out.
    out_indices: Vec<usize>,
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
    surface: Option<Box<dyn Surface>>,
    /// Computed sample count along U.
    nb_samples_u: usize,
    /// Computed sample count along V.
    nb_samples_v: usize,
    /// Lower U bound of the (clamped) parameter domain.
    u0: f64,
    /// Lower V bound of the (clamped) parameter domain.
    v0: f64,
    /// U step between interior samples: `(u_last - u0) / (nb_samples_u + 1)`.
    du: f64,
    /// V step between interior samples: `(v_last - v0) / (nb_samples_v + 1)`.
    dv: f64,
    /// Clamped `[u_first, u_last]` parameter domain.
    u_range: (f64, f64),
    /// Clamped `[v_first, v_last]` parameter domain.
    v_range: (f64, f64),
    /// U sample parameters produced by the last [`sample_pnts`](Self::sample_pnts)
    /// call (adaptive for BSpline, uniform otherwise). `None` before that.
    u_pars: Option<Vec<f64>>,
    /// V sample parameters produced by the last [`sample_pnts`](Self::sample_pnts)
    /// call.
    v_pars: Option<Vec<f64>>,
}

impl TopolTool {
    /// Empty tool — no surface. Call [`initialize`](Self::initialize) before
    /// sampling.
    pub fn new() -> Self {
        Self {
            surface: None,
            nb_samples_u: 0,
            nb_samples_v: 0,
            u0: 0.0,
            v0: 0.0,
            du: 1.0,
            dv: 1.0,
            u_range: (0.0, 1.0),
            v_range: (0.0, 1.0),
            u_pars: None,
            v_pars: None,
        }
    }

    /// Whether a surface has been attached.
    pub fn is_initialized(&self) -> bool {
        self.surface.is_some()
    }

    /// Binds `surface` (cloned) and computes the sample grid. Any adaptive
    /// parameters from a previous [`sample_pnts`](Self::sample_pnts) call are
    /// cleared. Source: `IntTools_TopolTool::Initialize`.
    pub fn initialize(&mut self, surface: &dyn Surface) {
        self.surface = Some(surface.clone_dyn());
        self.nb_samples_u = 0;
        self.nb_samples_v = 0;
        self.u0 = 0.0;
        self.v0 = 0.0;
        self.du = 1.0;
        self.dv = 1.0;
        self.u_pars = None;
        self.v_pars = None;
        self.compute_sample_points();
    }

    /// Computes `nb_samples_u`/`nb_samples_v` (and the `du`/`dv` steps) from
    /// the surface type and its parameter domain. Source:
    /// `IntTools_TopolTool::ComputeSamplePoints`.
    ///
    /// `ponytail:` analytic type is recovered via [`classify_surface`]; radii
    /// for the angular-step count are measured by sampling the surface, and
    /// BSpline/Bezier pole/knot counts (not reachable through `dyn Surface`)
    /// fall back to a 10×10 base grid.
    pub fn compute_sample_points(&mut self) {
        let Some(surface) = self.surface.as_ref() else { return };
        let s: &dyn Surface = surface.as_ref();
        let (mut uinf, mut usup) = s.u_range();
        let (mut vinf, mut vsup) = s.v_range();

        if usup < uinf {
            std::mem::swap(&mut uinf, &mut usup);
        }
        if vsup < vinf {
            std::mem::swap(&mut vinf, &mut vsup);
        }

        // Clamp unbounded directions to a big-but-finite span (OCCT sentinel).
        let is_big_uinf = !uinf.is_finite() && uinf < 0.0;
        let is_big_usup = !usup.is_finite() && usup > 0.0;
        let is_big_vinf = !vinf.is_finite() && vinf < 0.0;
        let is_big_vsup = !vsup.is_finite() && vsup > 0.0;
        if is_big_uinf && is_big_usup {
            uinf = -BIG_RANGE;
            usup = BIG_RANGE;
        } else if is_big_uinf {
            uinf = usup - 2.0 * BIG_RANGE;
        } else if is_big_usup {
            usup = uinf + 2.0 * BIG_RANGE;
        }
        if is_big_vinf && is_big_vsup {
            vinf = -BIG_RANGE;
            vsup = BIG_RANGE;
        } else if is_big_vinf {
            vinf = vsup - 2.0 * BIG_RANGE;
        } else if is_big_vsup {
            vsup = vinf + 2.0 * BIG_RANGE;
        }

        self.u0 = uinf;
        self.v0 = vinf;
        self.u_range = (uinf, usup);
        self.v_range = (vinf, vsup);
        self.u_pars = None;
        self.v_pars = None;

        let typ = classify_surface(s);
        let (mut nbsu, mut nbsv): (usize, usize) = match typ {
            SurfaceType::Plane => (10, 10),
            SurfaceType::Cylinder => {
                let radius = cylinder_radius(s);
                let max_angle = max_angle_for_radius(radius);
                let nbsu = if max_angle > ANGULAR {
                    ((usup - uinf) / max_angle) as usize
                } else {
                    0
                };
                let nbsv = ((vsup - vinf) / 10.0) as usize;
                (nbsu.max(2).min(MAX_NB_SAMPLE), nbsv.max(2).min(MAX_NB_SAMPLE))
            }
            SurfaceType::Cone => {
                let radius = cone_radius_at(s, vinf).max(cone_radius_at(s, vsup));
                let max_angle = max_angle_for_radius(radius);
                let nbsu = if max_angle > ANGULAR {
                    ((usup - uinf) / max_angle) as usize
                } else {
                    0
                };
                let nbsv = ((vsup - vinf) / 10.0) as usize;
                (nbsu.max(10).min(MAX_NB_SAMPLE), nbsv.max(10).min(MAX_NB_SAMPLE))
            }
            SurfaceType::Sphere => {
                let radius = sphere_radius(s);
                let max_angle = max_angle_for_radius(radius);
                let nbsu = if max_angle > ANGULAR {
                    ((usup - uinf) / max_angle) as usize
                } else {
                    0
                };
                let nbsv = if max_angle > ANGULAR {
                    ((vsup - vinf) / max_angle) as usize
                } else {
                    0
                };
                (nbsu.max(10).min(MAX_NB_SAMPLE), nbsv.max(10).min(MAX_NB_SAMPLE))
            }
            SurfaceType::Torus => {
                let (_major, minor) = torus_radii(s);
                let max_angle = max_angle_for_radius(minor);
                let nbsu = if max_angle > ANGULAR {
                    ((usup - uinf) / max_angle) as usize
                } else {
                    0
                };
                let nbsv = if max_angle > ANGULAR {
                    ((vsup - vinf) / max_angle) as usize
                } else {
                    0
                };
                (nbsu.max(10).min(MAX_NB_SAMPLE), nbsv.max(10).min(MAX_NB_SAMPLE))
            }
            SurfaceType::BezierSurface | SurfaceType::BSplineSurface => {
                // Pole/knot counts are unreachable through `dyn Surface`; the
                // base grid is refined adaptively by `sample_pnts`.
                (10, 10)
            }
            SurfaceType::SurfaceOfExtrusion => {
                let nbsv = ((vsup - vinf) / 10.0) as usize;
                (15, nbsv.max(15).min(MAX_NB_SAMPLE))
            }
            SurfaceType::SurfaceOfRevolution => (15, 15),
            SurfaceType::OffsetSurface | SurfaceType::OtherSurface => (10, 10),
        };
        if nbsu == 0 {
            nbsu = 10;
        }
        if nbsv == 0 {
            nbsv = 10;
        }

        self.nb_samples_u = nbsu;
        self.nb_samples_v = nbsv;
        self.du = (usup - uinf) / (nbsu + 1) as f64;
        self.dv = (vsup - vinf) / (nbsv + 1) as f64;
    }

    /// Number of sample points along U. Returns 0 when no surface is attached.
    pub fn nb_samples_u(&self) -> usize {
        self.nb_samples_u
    }

    /// Number of sample points along V. Returns 0 when no surface is attached.
    pub fn nb_samples_v(&self) -> usize {
        self.nb_samples_v
    }

    /// Total number of sample points: `nb_samples_u * nb_samples_v`.
    pub fn nb_samples(&self) -> usize {
        self.nb_samples_u * self.nb_samples_v
    }

    /// The clamped `[u_first, u_last]` parameter domain.
    pub fn u_range(&self) -> (f64, f64) {
        self.u_range
    }

    /// The clamped `[v_first, v_last]` parameter domain.
    pub fn v_range(&self) -> (f64, f64) {
        self.v_range
    }

    /// The U step between interior samples.
    pub fn u_step(&self) -> f64 {
        self.du
    }

    /// The V step between interior samples.
    pub fn v_step(&self) -> f64 {
        self.dv
    }

    /// Returns the `index`-th sample: the `(u, v)` parameter pair plus the 3D
    /// surface point `surface.d0(u, v)`.
    ///
    /// `index` is 1-based, from `1` to [`nb_samples`](Self::nb_samples), in
    /// row-major order (`U` fastest). When [`sample_pnts`](Self::sample_pnts)
    /// has been called, the adaptive parameter grid is used instead of the
    /// uniform one. Source: `IntTools_TopolTool::SamplePoint`.
    pub fn sample_point(&self, index: usize) -> Result<(GpPnt2d, GpPnt), String> {
        let s = self
            .surface
            .as_ref()
            .ok_or("TopolTool::sample_point: no surface initialized")?;
        let s: &dyn Surface = s.as_ref();
        if self.nb_samples_u == 0 || self.nb_samples_v == 0 {
            return Err("TopolTool::sample_point: sample grid not computed".into());
        }
        let (nu, nv) = match (&self.u_pars, &self.v_pars) {
            (Some(u), Some(v)) => (u.len(), v.len()),
            _ => (self.nb_samples_u, self.nb_samples_v),
        };
        if index == 0 || index > nu * nv {
            return Err(format!(
                "TopolTool::sample_point: index {index} out of range [1, {}]",
                nu * nv
            ));
        }
        let (iu, iv) = ((index - 1) % nu, (index - 1) / nu);
        let (u, v) = match (&self.u_pars, &self.v_pars) {
            (Some(u), Some(v)) => (u[iu], v[iv]),
            _ => (
                self.u0 + (iu + 1) as f64 * self.du,
                self.v0 + (iv + 1) as f64 * self.dv,
            ),
        };
        let p3d = s.d0(u, v);
        Ok((GpPnt2d::new(u, v), p3d))
    }

    /// Builds a full grid of sample points.
    ///
    /// For BSpline/Bezier surfaces the `U` and `V` parameter sequences are
    /// refined adaptively until the chord deviation of each interval stays
    /// under `deflection`; every other surface uses a uniform grid of
    /// `max(nb_samples, nu_min) × max(nb_samples, nv_min)` interior samples.
    /// The `u_pars`/`v_pars` are stored so subsequent
    /// [`sample_point`](Self::sample_point) calls walk the refined grid.
    ///
    /// Returns `(u, v, surface.d0(u, v))` for every parameter pair, in
    /// row-major order. Source: `IntTools_TopolTool::SamplePnts` +
    /// `Adaptor3d_TopolTool::SamplePnts`.
    pub fn sample_pnts(
        &mut self,
        deflection: f64,
        nu_min: usize,
        nv_min: usize,
    ) -> Result<Vec<(GpPnt2d, GpPnt)>, String> {
        if self.surface.is_none() {
            return Err("TopolTool::sample_pnts: no surface initialized".into());
        }
        // Recompute the analytic base grid (idempotent; also clears any stale
        // adaptive parameters from a previous call).
        self.compute_sample_points();
        let s = self.surface.as_ref().expect("checked above");
        let s: &dyn Surface = s.as_ref();
        let (u0, u1) = self.u_range;
        let (v0, v1) = self.v_range;
        let nbsu = self.nb_samples_u.max(nu_min).max(1);
        let nbsv = self.nb_samples_v.max(nv_min).max(1);

        let uniform_u = || -> Vec<f64> {
            let du = (u1 - u0) / (nbsu + 1) as f64;
            (1..=nbsu).map(|i| u0 + i as f64 * du).collect()
        };
        let uniform_v = || -> Vec<f64> {
            let dv = (v1 - v0) / (nbsv + 1) as f64;
            (1..=nbsv).map(|i| v0 + i as f64 * dv).collect()
        };

        let typ = classify_surface(s);
        let (u_pars, v_pars) = if matches!(
            typ,
            SurfaceType::BSplineSurface | SurfaceType::BezierSurface
        ) {
            let u_base = uniform_u();
            let v_base = uniform_v();
            let v_fixed = 0.5 * (v0 + v1);
            let u_fixed = 0.5 * (u0 + u1);
            (
                refine_params(s, deflection, u_base, v_fixed, true, (u0, u1)),
                refine_params(s, deflection, v_base, u_fixed, false, (v0, v1)),
            )
        } else {
            (uniform_u(), uniform_v())
        };

        self.u_pars = Some(u_pars.clone());
        self.v_pars = Some(v_pars.clone());
        self.nb_samples_u = u_pars.len();
        self.nb_samples_v = v_pars.len();

        let mut out = Vec::with_capacity(u_pars.len() * v_pars.len());
        for &u in &u_pars {
            for &v in &v_pars {
                out.push((GpPnt2d::new(u, v), s.d0(u, v)));
            }
        }
        Ok(out)
    }
}

impl Default for TopolTool {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Surface-classification helpers (measured through `dyn Surface`)
// ---------------------------------------------------------------------------

/// Angular step for a circle of `radius` under the 1e-2 sampling deflection.
/// Source: `IntTools_TopolTool::ComputeSamplePoints`
/// (`acos(1 - deflection / radius) * 2`, floored by `π/2`).
fn max_angle_for_radius(radius: f64) -> f64 {
    let mut max_angle = std::f64::consts::PI * 0.5;
    if radius > SAMPLE_DEFLECTION {
        max_angle = (1.0 - SAMPLE_DEFLECTION / radius).acos() * 2.0;
    }
    max_angle
}

/// Radius of a cylinder measured from the surface: half the distance between
/// the diametrically opposite points `(0, v)` and `(π, v)`.
fn cylinder_radius(s: &dyn Surface) -> f64 {
    s.d0(0.0, 0.0).distance(&s.d0(std::f64::consts::PI, 0.0)) * 0.5
}

/// Radius of a sphere measured from the surface.
fn sphere_radius(s: &dyn Surface) -> f64 {
    s.d0(0.0, 0.0).distance(&s.d0(std::f64::consts::PI, 0.0)) * 0.5
}

/// Radius of the U-circle of a cone at parameter `v`, measured from the
/// surface.
fn cone_radius_at(s: &dyn Surface, v: f64) -> f64 {
    s.d0(0.0, v).distance(&s.d0(std::f64::consts::PI, v)) * 0.5
}

/// `(major, minor)` radii of a torus measured from the surface.
fn torus_radii(s: &dyn Surface) -> (f64, f64) {
    let minor = s.d0(0.0, 0.0).distance(&s.d0(0.0, std::f64::consts::PI)) * 0.5;
    let major_plus_minor = s.d0(0.0, 0.0).distance(&s.d0(std::f64::consts::PI, 0.0)) * 0.5;
    (major_plus_minor - minor, minor)
}

/// Chord deviation of the surface point at parameter `m` against the linear
/// interpolation of the points at `a` and `b`, along the `is_u` iso-line
/// (the orthogonal parameter held at `fixed`).
fn slice_deviation(s: &dyn Surface, a: f64, b: f64, m: f64, fixed: f64, is_u: bool) -> f64 {
    let pa = if is_u { s.d0(a, fixed) } else { s.d0(fixed, a) };
    let pb = if is_u { s.d0(b, fixed) } else { s.d0(fixed, b) };
    let pm = if is_u { s.d0(m, fixed) } else { s.d0(fixed, m) };
    let t = (m - a) / (b - a);
    let interp = GpPnt::new(
        pa.x() + t * (pb.x() - pa.x()),
        pa.y() + t * (pb.y() - pa.y()),
        pa.z() + t * (pb.z() - pa.z()),
    );
    pm.distance(&interp)
}

/// Deflection-adaptive 1D refinement of a parameter sequence.
///
/// Starts from the interior `base` samples plus the two domain endpoints;
/// every interval whose midpoint deviates from its chord by more than
/// `deflection` is split at the midpoint. Repeats until no interval needs
/// splitting (bounded by a safety cap). The returned sequence drops the domain
/// endpoints, keeping the interior-only semantics of the uniform grid (and of
/// OCCT's `SamplePoint` `iu`/`iv` in `1..=nb`).
fn refine_params(
    s: &dyn Surface,
    deflection: f64,
    base: Vec<f64>,
    fixed: f64,
    is_u: bool,
    domain: (f64, f64),
) -> Vec<f64> {
    const MAX_PARAMS: usize = 4096;
    let deflection = deflection.max(1e-9);
    let mut params: Vec<f64> = base;
    params.push(domain.0);
    params.push(domain.1);
    params.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    params.dedup_by(|a, b| (*a - *b).abs() < 1e-12);

    loop {
        let mut changed = false;
        let mut next: Vec<f64> = Vec::with_capacity(params.len().min(MAX_PARAMS) + 16);
        for w in params.windows(2) {
            let (a, b) = (w[0], w[1]);
            if b - a <= 1e-12 {
                next.push(a);
                continue;
            }
            let m = 0.5 * (a + b);
            let dev = slice_deviation(s, a, b, m, fixed, is_u);
            if dev > deflection && next.len() < MAX_PARAMS {
                next.push(a);
                next.push(m);
                changed = true;
            } else {
                next.push(a);
            }
        }
        if let Some(&last) = params.last() {
            next.push(last);
        }
        next.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        next.dedup_by(|a, b| (*a - *b).abs() < 1e-12);
        params = next;
        if !changed || params.len() >= MAX_PARAMS {
            break;
        }
    }

    params
        .into_iter()
        .filter(|&p| p > domain.0 + 1e-9 && p < domain.1 - 1e-9)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;
    use std::sync::Arc;

    use occt_core::gp::{GpAx3, GpCylinder, GpPln};
    use occt_geom::bspline_surface::{bspline_surface_uniform_knots, GeomBSplineSurface};
    use occt_geom::{GeomCylinder, GeomPlane};

    use crate::brep_extrema::test_box::unit_box;
    use crate::tgeometry::GeometryRegistry;

    /// Test helper: an unchecked range.
    fn r(a: f64, b: f64) -> IntRange {
        IntRange::new_unchecked(a, b)
    }

    fn plane_surface() -> Arc<dyn Surface> {
        Arc::new(GeomPlane::new(GpPln::new(GpAx3::standard())))
    }

    fn cylinder_surface(radius: f64) -> Arc<dyn Surface> {
        Arc::new(GeomCylinder::new(GpCylinder::new(GpAx3::standard(), radius).unwrap()))
    }

    /// A degree-3 BSpline patch over `[0,3]²` that is curved along U (a cubic
    /// in `u`) and constant along V — a ruled "extrusion" of a cubic curve.
    fn curved_bspline() -> Arc<dyn Surface> {
        let (ku, kv) = bspline_surface_uniform_knots(4, 4, 3, 3);
        let poles: Vec<Vec<GpPnt>> = (0..4)
            .map(|i| {
                let u = i as f64;
                (0..4).map(|j| GpPnt::new(u, j as f64, u * u)).collect()
            })
            .collect();
        Arc::new(GeomBSplineSurface::new(poles, ku, kv, 3, 3).unwrap())
    }

    // ---- BaseRangeSample ----

    #[test]
    fn base_range_sample_depth_access() {
        let mut s = BaseRangeSample::new();
        assert_eq!(s.get_depth(), 0);
        s.set_depth(3);
        assert_eq!(s.get_depth(), 3);
        let s2 = BaseRangeSample::with_depth(5);
        assert_eq!(s2.get_depth(), 5);
        assert_eq!(BaseRangeSample::default().get_depth(), 0);
    }

    // ---- CurveRangeSample ----

    #[test]
    fn curve_range_sample_construct_depth_index_equal() {
        let a = CurveRangeSample::with_index_depth(7, 2);
        assert_eq!(a.get_depth(), 2);
        assert_eq!(a.get_index(), 7);
        let b = CurveRangeSample::with_index_depth(7, 2);
        assert!(a.is_equal(&b));
        assert_eq!(a, b);
        let c = CurveRangeSample::with_index_depth(7, 3);
        assert!(!a.is_equal(&c), "different depth");
        let d = CurveRangeSample::with_index_depth(8, 2);
        assert!(!a.is_equal(&d), "different index");
        let e = CurveRangeSample::with_index(4);
        assert_eq!(e.get_depth(), 0);
        assert_eq!(e.get_index(), 4);
    }

    #[test]
    fn curve_range_sample_get_range_depth_zero_covers_domain() {
        let s = CurveRangeSample::new();
        assert_eq!(s.get_range(-2.0, 5.0, 10), r(-2.0, 5.0));
    }

    #[test]
    fn curve_range_sample_get_range_subdivides() {
        // Depth 1 with 4 samples splits [0, 1] into 4 quarters.
        let s = CurveRangeSample::with_index_depth(2, 1);
        assert_eq!(s.get_range(0.0, 1.0, 4), r(0.5, 0.75));
        // Index 0 -> first quarter.
        let s0 = CurveRangeSample::with_index_depth(0, 1);
        assert_eq!(s0.get_range(0.0, 1.0, 4), r(0.0, 0.25));
        // Last index -> last quarter (upper bound inclusive).
        let s3 = CurveRangeSample::with_index_depth(3, 1);
        assert_eq!(s3.get_range(0.0, 1.0, 4), r(0.75, 1.0));
        // Depth 2 squares the sample count: 4^2 = 16 intervals of length 1/16.
        let sd = CurveRangeSample::with_index_depth(5, 2);
        assert_eq!(sd.get_range(0.0, 1.0, 4), r(5.0 / 16.0, 6.0 / 16.0));
    }

    #[test]
    fn curve_range_sample_index_deeper() {
        let s = CurveRangeSample::with_index_depth(3, 1);
        assert_eq!(s.get_range_index_deeper(4), 12);
        assert_eq!(CurveRangeSample::with_index(0).get_range_index_deeper(10), 0);
    }

    #[test]
    fn curve_range_sample_setters() {
        let mut s = CurveRangeSample::new();
        s.set_index(2);
        s.set_depth(1);
        assert_eq!(s, CurveRangeSample::with_index_depth(2, 1));
    }

    // ---- SurfaceRangeSample ----

    #[test]
    fn surface_range_sample_construct_and_accessors() {
        let s = SurfaceRangeSample::with_indexes_depths(2, 1, 3, 2);
        assert_eq!(s.get_index_u(), 2);
        assert_eq!(s.get_depth_u(), 1);
        assert_eq!(s.get_index_v(), 3);
        assert_eq!(s.get_depth_v(), 2);
        assert_eq!(s.get_indexes(), (2, 3));
        assert_eq!(s.get_depths(), (1, 2));

        let mut m = SurfaceRangeSample::new();
        m.set_index_u(5);
        m.set_index_v(6);
        m.set_depth_u(1);
        m.set_depth_v(1);
        assert_eq!(m.get_indexes(), (5, 6));
        assert_eq!(m.get_depths(), (1, 1));

        let (ru, rv) = m.get_ranges();
        assert_eq!(ru, CurveRangeSample::with_index_depth(5, 1));
        assert_eq!(rv, CurveRangeSample::with_index_depth(6, 1));
    }

    #[test]
    fn surface_range_sample_equality() {
        let a = SurfaceRangeSample::with_indexes_depths(1, 1, 2, 1);
        let b = SurfaceRangeSample::with_indexes_depths(1, 1, 2, 1);
        let c = SurfaceRangeSample::with_indexes_depths(1, 1, 2, 2);
        assert!(a.is_equal(&b));
        assert_eq!(a, b);
        assert!(!a.is_equal(&c), "V depth differs");
    }

    #[test]
    fn surface_range_sample_get_ranges_and_deeper() {
        let s = SurfaceRangeSample::with_indexes_depths(1, 1, 2, 1);
        assert_eq!(s.get_range_u(0.0, 1.0, 4), r(0.25, 0.5));
        assert_eq!(s.get_range_v(0.0, 2.0, 4), r(1.0, 1.5));
        let (ur, vr) = s.get_range(0.0, 1.0, 4, 0.0, 2.0, 4);
        assert_eq!(ur, r(0.25, 0.5));
        assert_eq!(vr, r(1.0, 1.5));
        assert_eq!(s.get_range_index_u_deeper(4), 4);
        assert_eq!(s.get_range_index_v_deeper(4), 8);
    }

    #[test]
    fn surface_range_sample_from_ranges() {
        let ru = CurveRangeSample::with_index_depth(1, 1);
        let rv = CurveRangeSample::with_index_depth(2, 1);
        let s = SurfaceRangeSample::from_ranges(ru, rv);
        assert_eq!(s.get_sample_range_u(), ru);
        assert_eq!(s.get_sample_range_v(), rv);
        s.get_ranges();
    }

    // ---- CurveRangeLocalizeData ----

    #[test]
    fn curve_localize_data_build_maps_intervals() {
        let mut ld = CurveRangeLocalizeData::new(4, 1e-6);
        ld.set_root_index(2);
        let n = ld.build(r(0.0, 1.0), 7);
        assert_eq!(n, 4);
        assert_eq!(ld.get_nb_sample(), 4);
        assert_eq!(ld.get_root_index(), 2);
        assert_eq!(ld.ranges().len(), 4);
        assert_eq!(ld.range(0), Some(r(0.0, 0.25)));
        assert_eq!(ld.range(1), Some(r(0.25, 0.5)));
        assert_eq!(ld.range(2), Some(r(0.5, 0.75)));
        assert_eq!(ld.range(3), Some(r(0.75, 1.0)));
        for i in 0..4 {
            assert_eq!(ld.curve_index(i), Some(7), "all cells mapped to curve 7");
        }
        // Every parameter is located in exactly one cell.
        assert_eq!(ld.find_index(0.0), Some(0));
        assert_eq!(ld.find_index(0.3), Some(1));
        assert_eq!(ld.find_index(0.99), Some(3));
        assert_eq!(ld.find_index(1.0), Some(3), "upper bound inclusive");
        assert_eq!(ld.find_index(-0.1), None);
    }

    #[test]
    fn curve_localize_data_reassign_and_out() {
        let mut ld = CurveRangeLocalizeData::new(3, 0.0);
        ld.build(r(0.0, 3.0), 0);
        ld.set_curve_index(1, 9).unwrap();
        assert_eq!(ld.curve_index(1), Some(9));
        assert_eq!(ld.curve_index(0), Some(0));
        assert!(ld.set_curve_index(9, 1).is_err(), "out of range rejected");

        assert!(!ld.is_range_out(0));
        ld.add_out_range(0).unwrap();
        assert!(ld.is_range_out(0));
        assert!(!ld.is_range_out(1));
        ld.add_out_range(0).unwrap(); // idempotent
        assert_eq!(ld.list_range_out(), &[0]);
        ld.remove_range_out_all();
        assert!(!ld.is_range_out(0));
        assert!(ld.add_out_range(99).is_err());
    }

    // ---- SurfaceRangeLocalizeData ----

    #[test]
    fn surface_localize_data_build_grid() {
        let mut ld = SurfaceRangeLocalizeData::new(2, 3, 1e-6, 1e-6);
        ld.set_root_index(5);
        let n = ld.build(r(0.0, 1.0), r(0.0, 3.0), 11);
        assert_eq!(n, 6, "2 × 3 cells");
        assert_eq!(ld.get_nb_samples_u(), 2);
        assert_eq!(ld.get_nb_samples_v(), 3);
        assert_eq!(ld.get_root_index(), 5);
        assert_eq!(ld.ranges().len(), 6);

        // Cell (0, 0): U in [0, 0.5], V in [0, 1].
        let (ur, vr) = ld.get_depth_range(0, 0).unwrap();
        assert_eq!(ur, r(0.0, 0.5));
        assert_eq!(vr, r(0.0, 1.0));
        // Cell (1, 1): U in [0.5, 1], V in [1, 2].
        let (ur, vr) = ld.get_depth_range(1, 1).unwrap();
        assert_eq!(ur, r(0.5, 1.0));
        assert_eq!(vr, r(1.0, 2.0));
        // Last cell (1, 2): V upper bound inclusive.
        let (ur, vr) = ld.get_depth_range(1, 2).unwrap();
        assert_eq!(ur, r(0.5, 1.0));
        assert_eq!(vr, r(2.0, 3.0));

        for iu in 0..2 {
            for iv in 0..3 {
                assert_eq!(ld.surface_index(iu, iv), Some(11));
            }
        }
        assert_eq!(ld.get_depth_range(2, 0), None, "out of grid");

        // Parametric lookup.
        assert_eq!(ld.find_cell(0.25, 0.5), Some((0, 0)));
        assert_eq!(ld.find_cell(0.75, 2.5), Some((1, 2)));
        assert_eq!(ld.find_cell(-1.0, 0.0), None);
    }

    #[test]
    fn surface_localize_data_out_tracking() {
        let mut ld = SurfaceRangeLocalizeData::new(2, 2, 0.0, 0.0);
        ld.build(r(0.0, 1.0), r(0.0, 1.0), 0);
        assert!(!ld.is_range_out(1, 1));
        ld.add_out_range(1, 1).unwrap();
        assert!(ld.is_range_out(1, 1));
        assert!(!ld.is_range_out(0, 1));
        assert!(ld.add_out_range(5, 5).is_err());
        ld.remove_range_out_all();
        assert!(!ld.is_range_out(1, 1));
    }

    // ---- TopolTool ----

    #[test]
    fn topol_tool_plane_nb_samples_and_domain() {
        let mut tool = TopolTool::new();
        assert!(!tool.is_initialized());
        tool.initialize(plane_surface().as_ref());
        assert!(tool.is_initialized());
        // A plane gets a 10 × 10 grid, domain clamped to [-1e5, 1e5]².
        assert_eq!(tool.nb_samples_u(), 10);
        assert_eq!(tool.nb_samples_v(), 10);
        assert_eq!(tool.nb_samples(), 100);
        let (u0, u1) = tool.u_range();
        let (v0, v1) = tool.v_range();
        assert!((u0 + 1e5).abs() < 1e-9 && (u1 - 1e5).abs() < 1e-9);
        assert!((v0 + 1e5).abs() < 1e-9 && (v1 - 1e5).abs() < 1e-9);
        // The step spans the clamped domain.
        assert!((tool.u_step() - (2e5 / 11.0)).abs() < 1e-9);
    }

    #[test]
    fn topol_tool_plane_sample_point_matches_d0() {
        let surf = plane_surface();
        let mut tool = TopolTool::new();
        tool.initialize(surf.as_ref());
        // The standard plane maps (u, v) -> (u, v, 0); every sample's 3D point
        // must agree with surface.d0(u, v) and lie in the plane.
        for i in 1..=tool.nb_samples() {
            let (p2d, p3d) = tool.sample_point(i).unwrap();
            let expected = surf.d0(p2d.x(), p2d.y());
            assert!(p3d.distance(&expected) < 1e-12, "sample {i}: 3D != d0(UV)");
            assert!(p3d.z().abs() < 1e-9, "plane sample {i} off the plane");
            assert!(p2d.x() >= -1e5 && p2d.x() <= 1e5);
            assert!(p2d.y() >= -1e5 && p2d.y() <= 1e5);
        }
    }

    #[test]
    fn topol_tool_plane_sample_point_row_major() {
        let mut tool = TopolTool::new();
        tool.initialize(plane_surface().as_ref());
        let (p1, _) = tool.sample_point(1).unwrap();
        let (p2, _) = tool.sample_point(2).unwrap();
        // Row-major, U fastest: consecutive samples differ by one U step.
        assert!((p2.x() - p1.x()).abs() - tool.u_step() < 1e-9);
        assert!((p2.y() - p1.y()).abs() < 1e-9);
        // Sample 11 is the start of the second V row.
        let (p11, _) = tool.sample_point(11).unwrap();
        assert!((p11.y() - p1.y()).abs() - tool.v_step() < 1e-9);
        assert!((p11.x() - p1.x()).abs() < 1e-9);
    }

    #[test]
    fn topol_tool_cylinder_nb_samples_match_domain() {
        let surf = cylinder_surface(1.0);
        let mut tool = TopolTool::new();
        tool.initialize(surf.as_ref());
        // Radius 1 -> max_angle = 2*acos(1 - 0.01); U count from the 2π span,
        // V count clamped to the 50 cap (2e5/10 = 2e4).
        let max_angle: f64 = 2.0 * (1.0 - 0.01f64).acos();
        let expected_u = (2.0 * PI / max_angle) as usize;
        assert_eq!(tool.nb_samples_u(), expected_u.max(2));
        assert_eq!(tool.nb_samples_v(), 50);
        assert_eq!(tool.nb_samples(), tool.nb_samples_u() * tool.nb_samples_v());
    }

    #[test]
    fn topol_tool_cylinder_sample_point_on_surface() {
        let surf = cylinder_surface(1.0);
        let mut tool = TopolTool::new();
        tool.initialize(surf.as_ref());
        for i in 1..=tool.nb_samples() {
            let (p2d, p3d) = tool.sample_point(i).unwrap();
            let expected = surf.d0(p2d.x(), p2d.y());
            assert!(p3d.distance(&expected) < 1e-12, "sample {i}: 3D != d0(UV)");
            // Cylinder about the Z axis: distance to the axis == radius.
            let radial = (p3d.x() * p3d.x() + p3d.y() * p3d.y()).sqrt();
            assert!((radial - 1.0).abs() < 1e-9, "sample {i} radial {radial} != 1");
        }
    }

    #[test]
    fn topol_tool_sample_point_out_of_range_errors() {
        let mut tool = TopolTool::new();
        assert!(tool.sample_point(1).is_err(), "no surface");
        tool.initialize(plane_surface().as_ref());
        assert!(tool.sample_point(0).is_err());
        assert!(tool.sample_point(101).is_err(), "100 samples, 101 is out");
        assert!(tool.sample_point(100).is_ok());
    }

    #[test]
    fn topol_tool_sample_pnts_plane_uniform_grid_count() {
        let mut tool = TopolTool::new();
        tool.initialize(plane_surface().as_ref());
        let pts = tool.sample_pnts(0.001, 4, 4).unwrap();
        // Plane is not a BSpline: uniform grid of max(10, 4)² = 100 points.
        assert_eq!(pts.len(), 100);
        for (p2d, p3d) in &pts {
            let expected = plane_surface().d0(p2d.x(), p2d.y());
            assert!(p3d.distance(&expected) < 1e-12);
            assert!(p3d.z().abs() < 1e-9);
        }
    }

    #[test]
    fn topol_tool_sample_pnts_plane_floor_applies() {
        let mut tool = TopolTool::new();
        tool.initialize(plane_surface().as_ref());
        let pts = tool.sample_pnts(0.001, 20, 20).unwrap();
        // nu_min/nv_min are a floor: max(10, 20)² = 400.
        assert_eq!(pts.len(), 400);
    }

    #[test]
    fn topol_tool_sample_pnts_cylinder_grid_on_surface() {
        let surf = cylinder_surface(1.0);
        let mut tool = TopolTool::new();
        tool.initialize(surf.as_ref());
        let pts = tool.sample_pnts(0.001, 4, 4).unwrap();
        let expected_len = tool.nb_samples_u() * tool.nb_samples_v();
        assert_eq!(pts.len(), expected_len, "cylinder stays uniform, count from base grid");
        for (p2d, p3d) in &pts {
            let expected = surf.d0(p2d.x(), p2d.y());
            assert!(p3d.distance(&expected) < 1e-12);
            let radial = (p3d.x() * p3d.x() + p3d.y() * p3d.y()).sqrt();
            assert!((radial - 1.0).abs() < 1e-9, "point off cylinder");
        }
    }

    #[test]
    fn topol_tool_sample_pnts_bspline_adaptive_refines() {
        let surf = curved_bspline();
        let mut coarse = TopolTool::new();
        coarse.initialize(surf.as_ref());
        let coarse_pts = coarse.sample_pnts(10.0, 4, 4).unwrap();

        let mut fine = TopolTool::new();
        fine.initialize(surf.as_ref());
        let fine_pts = fine.sample_pnts(1e-4, 4, 4).unwrap();

        // Small deflection must refine the curved (U) direction, producing a
        // strictly larger sample set than the coarse grid.
        assert!(fine_pts.len() > coarse_pts.len());
        assert!(coarse_pts.len() >= 16, "at least the 4×4 floor");
        for (p2d, p3d) in &fine_pts {
            let expected = surf.d0(p2d.x(), p2d.y());
            assert!(p3d.distance(&expected) < 1e-12, "fine sample not on surface");
        }
    }

    #[test]
    fn topol_tool_unit_box_face_sample_points() {
        let box_ = unit_box();
        let face = &box_.faces[0];
        let surf = GeometryRegistry::global()
            .face_surface(&face.0)
            .expect("unit-box face has a registered surface");
        let mut tool = TopolTool::new();
        tool.initialize(surf.as_ref());
        assert_eq!(tool.nb_samples_u(), 10);
        assert_eq!(tool.nb_samples_v(), 10);
        for i in 1..=tool.nb_samples() {
            let (p2d, p3d) = tool.sample_point(i).unwrap();
            let expected = surf.d0(p2d.x(), p2d.y());
            assert!(p3d.distance(&expected) < 1e-12, "box face sample {i}");
        }
        let pts = tool.sample_pnts(0.001, 4, 4).unwrap();
        assert_eq!(pts.len(), 100, "box faces are planes -> uniform 10×10");
    }
}
