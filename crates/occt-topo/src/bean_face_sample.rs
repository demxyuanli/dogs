//! Range-sample and localize-data types for bean/face localization.
//!
//! Ports `IntTools_BaseRangeSample`, `IntTools_CurveRangeSample`,
//! `IntTools_SurfaceRangeSample`, `IntTools_CurveRangeLocalizeData`, and
//! `IntTools_SurfaceRangeLocalizeData` used by
//! `IntTools_BeanFaceIntersector::LocalizeSolutions`.

use std::collections::{HashMap, HashSet};

use occt_core::bnd::BndBox;
use occt_core::gp::GpPnt;

use crate::inttools_data::IntRange;

/// `IntTools_BaseRangeSample` — subdivision depth of a parameter range.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub(crate) struct BaseRangeSample {
    pub depth: i32,
}

impl BaseRangeSample {
    pub(crate) fn new() -> Self {
        Self { depth: 0 }
    }

    pub(crate) fn with_depth(depth: i32) -> Self {
        Self { depth }
    }
}

/// `IntTools_CurveRangeSample` — (depth, index) of a curve parameter cell.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub(crate) struct CurveRangeSample {
    pub base: BaseRangeSample,
    pub index: i32,
}

impl CurveRangeSample {
    pub(crate) fn new() -> Self {
        Self {
            base: BaseRangeSample::new(),
            index: 0,
        }
    }

    pub(crate) fn with_index(index: i32) -> Self {
        Self {
            base: BaseRangeSample::new(),
            index,
        }
    }

    pub(crate) fn set_depth(&mut self, depth: i32) {
        self.base.depth = depth;
    }

    pub(crate) fn depth(&self) -> i32 {
        self.base.depth
    }

    pub(crate) fn set_range_index(&mut self, index: i32) {
        self.index = index;
    }

    pub(crate) fn range_index(&self) -> i32 {
        self.index
    }

    /// `GetRangeIndexDeeper`.
    pub(crate) fn range_index_deeper(&self, nb_sample: i32) -> i32 {
        self.index * nb_sample
    }

    /// `IntTools_CurveRangeSample::GetRange`.
    pub(crate) fn get_range(&self, first: f64, last: f64, nb_sample: i32) -> IntRange {
        if self.depth() <= 0 {
            return IntRange::new_unchecked(first, last);
        }
        let n = nb_sample.max(1) as f64;
        let tmp = n.powi(self.depth());
        let local = (last - first) / tmp;
        let a_first = first + self.index as f64 * local;
        IntRange::new_unchecked(a_first, a_first + local)
    }
}

/// `IntTools_SurfaceRangeSample` — independent U/V curve-range samples.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub(crate) struct SurfaceRangeSample {
    pub range_u: CurveRangeSample,
    pub range_v: CurveRangeSample,
}

impl SurfaceRangeSample {
    pub(crate) fn new() -> Self {
        Self {
            range_u: CurveRangeSample::new(),
            range_v: CurveRangeSample::new(),
        }
    }

    pub(crate) fn with_indices(index_u: i32, depth_u: i32, index_v: i32, depth_v: i32) -> Self {
        let mut u = CurveRangeSample::with_index(index_u);
        u.set_depth(depth_u);
        let mut v = CurveRangeSample::with_index(index_v);
        v.set_depth(depth_v);
        Self {
            range_u: u,
            range_v: v,
        }
    }

    pub(crate) fn set_depth_u(&mut self, d: i32) {
        self.range_u.set_depth(d);
    }

    pub(crate) fn set_depth_v(&mut self, d: i32) {
        self.range_v.set_depth(d);
    }

    pub(crate) fn depth_u(&self) -> i32 {
        self.range_u.depth()
    }

    pub(crate) fn depth_v(&self) -> i32 {
        self.range_v.depth()
    }

    pub(crate) fn set_index_u(&mut self, i: i32) {
        self.range_u.set_range_index(i);
    }

    pub(crate) fn set_index_v(&mut self, i: i32) {
        self.range_v.set_range_index(i);
    }

    pub(crate) fn index_u(&self) -> i32 {
        self.range_u.range_index()
    }

    pub(crate) fn index_v(&self) -> i32 {
        self.range_v.range_index()
    }

    pub(crate) fn range_index_u_deeper(&self, nb_sample_u: i32) -> i32 {
        self.range_u.range_index_deeper(nb_sample_u)
    }

    pub(crate) fn range_index_v_deeper(&self, nb_sample_v: i32) -> i32 {
        self.range_v.range_index_deeper(nb_sample_v)
    }

    pub(crate) fn get_range_u(&self, first: f64, last: f64, nb_sample: i32) -> IntRange {
        self.range_u.get_range(first, last, nb_sample)
    }

    pub(crate) fn get_range_v(&self, first: f64, last: f64, nb_sample: i32) -> IntRange {
        self.range_v.get_range(first, last, nb_sample)
    }

    pub(crate) fn is_equal(&self, other: &Self) -> bool {
        self == other
    }
}

/// `IntTools_CurveRangeLocalizeData`.
#[derive(Debug, Clone)]
pub(crate) struct CurveRangeLocalizeData {
    nb_sample: i32,
    min_range: f64,
    out_ranges: HashSet<CurveRangeSample>,
    boxes: HashMap<CurveRangeSample, BndBox>,
}

impl CurveRangeLocalizeData {
    pub(crate) fn new(nb_sample: i32, min_range: f64) -> Self {
        Self {
            nb_sample,
            min_range,
            out_ranges: HashSet::new(),
            boxes: HashMap::new(),
        }
    }

    pub(crate) fn nb_sample(&self) -> i32 {
        self.nb_sample
    }

    pub(crate) fn min_range(&self) -> f64 {
        self.min_range
    }

    pub(crate) fn add_out_range(&mut self, range: CurveRangeSample) {
        self.out_ranges.insert(range);
    }

    pub(crate) fn add_box(&mut self, range: CurveRangeSample, box_: BndBox) {
        self.boxes.insert(range, box_);
    }

    pub(crate) fn find_box(&self, range: &CurveRangeSample) -> Option<BndBox> {
        self.boxes.get(range).copied()
    }

    pub(crate) fn is_range_out(&self, range: &CurveRangeSample) -> bool {
        self.out_ranges.contains(range)
    }

    pub(crate) fn list_range_out(&self) -> Vec<CurveRangeSample> {
        self.out_ranges.iter().copied().collect()
    }
}

/// `IntTools_SurfaceRangeLocalizeData` including the optional BSpline grid.
#[derive(Debug, Clone)]
pub(crate) struct SurfaceRangeLocalizeData {
    nb_sample_u: i32,
    nb_sample_v: i32,
    min_range_u: f64,
    min_range_v: f64,
    out_ranges: HashSet<SurfaceRangeSample>,
    boxes: HashMap<SurfaceRangeSample, BndBox>,
    u_params: Vec<f64>,
    v_params: Vec<f64>,
    grid_points: Vec<Vec<GpPnt>>,
    u_ind_min: i32,
    u_ind_max: i32,
    v_ind_min: i32,
    v_ind_max: i32,
    deflection: f64,
}

impl SurfaceRangeLocalizeData {
    pub(crate) fn new() -> Self {
        Self {
            nb_sample_u: 1,
            nb_sample_v: 1,
            min_range_u: 0.0,
            min_range_v: 0.0,
            out_ranges: HashSet::new(),
            boxes: HashMap::new(),
            u_params: Vec::new(),
            v_params: Vec::new(),
            grid_points: Vec::new(),
            u_ind_min: 0,
            u_ind_max: 0,
            v_ind_min: 0,
            v_ind_max: 0,
            deflection: 0.0,
        }
    }

    pub(crate) fn with_samples(
        nb_sample_u: i32,
        nb_sample_v: i32,
        min_range_u: f64,
        min_range_v: f64,
    ) -> Self {
        let mut s = Self::new();
        s.nb_sample_u = nb_sample_u;
        s.nb_sample_v = nb_sample_v;
        s.min_range_u = min_range_u;
        s.min_range_v = min_range_v;
        s
    }

    pub(crate) fn nb_sample_u(&self) -> i32 {
        self.nb_sample_u
    }

    pub(crate) fn nb_sample_v(&self) -> i32 {
        self.nb_sample_v
    }

    pub(crate) fn min_range_u(&self) -> f64 {
        self.min_range_u
    }

    pub(crate) fn min_range_v(&self) -> f64 {
        self.min_range_v
    }

    pub(crate) fn add_out_range(&mut self, range: SurfaceRangeSample) {
        self.out_ranges.insert(range);
    }

    pub(crate) fn add_box(&mut self, range: SurfaceRangeSample, box_: BndBox) {
        self.boxes.insert(range, box_);
    }

    pub(crate) fn find_box(&self, range: &SurfaceRangeSample) -> Option<BndBox> {
        self.boxes.get(range).copied()
    }

    pub(crate) fn is_range_out(&self, range: &SurfaceRangeSample) -> bool {
        self.out_ranges.contains(range)
    }

    pub(crate) fn remove_range_out_all(&mut self) {
        self.out_ranges.clear();
    }

    pub(crate) fn set_grid_deflection(&mut self, d: f64) {
        self.deflection = d;
    }

    pub(crate) fn grid_deflection(&self) -> f64 {
        self.deflection
    }

    /// `SetRangeUGrid` — 1-based OCCT arrays are stored 0-based here.
    pub(crate) fn set_range_u_grid(&mut self, nb: i32) {
        self.u_ind_min = 0;
        self.u_ind_max = 0;
        self.v_ind_min = 0;
        self.v_ind_max = 0;
        let n = nb.max(0) as usize;
        if self.u_params.len() != n {
            self.u_params = vec![0.0; n];
            if !self.v_params.is_empty() {
                self.grid_points = vec![vec![GpPnt::zero(); self.v_params.len()]; n];
            }
        }
    }

    pub(crate) fn set_range_v_grid(&mut self, nb: i32) {
        self.u_ind_min = 0;
        self.u_ind_max = 0;
        self.v_ind_min = 0;
        self.v_ind_max = 0;
        let n = nb.max(0) as usize;
        if self.v_params.len() != n {
            self.v_params = vec![0.0; n];
            if !self.u_params.is_empty() {
                self.grid_points = vec![vec![GpPnt::zero(); n]; self.u_params.len()];
            }
        }
    }

    pub(crate) fn range_u_grid(&self) -> i32 {
        self.u_params.len() as i32
    }

    pub(crate) fn range_v_grid(&self) -> i32 {
        self.v_params.len() as i32
    }

    /// OCCT 1-based `SetUParam`.
    pub(crate) fn set_u_param(&mut self, index: i32, u: f64) {
        let i = (index - 1) as usize;
        if i < self.u_params.len() {
            self.u_params[i] = u;
        }
    }

    pub(crate) fn u_param(&self, index: i32) -> f64 {
        self.u_params
            .get((index - 1) as usize)
            .copied()
            .unwrap_or(0.0)
    }

    pub(crate) fn set_v_param(&mut self, index: i32, v: f64) {
        let i = (index - 1) as usize;
        if i < self.v_params.len() {
            self.v_params[i] = v;
        }
    }

    pub(crate) fn v_param(&self, index: i32) -> f64 {
        self.v_params
            .get((index - 1) as usize)
            .copied()
            .unwrap_or(0.0)
    }

    pub(crate) fn set_grid_point(&mut self, u_index: i32, v_index: i32, p: GpPnt) {
        let i = (u_index - 1) as usize;
        let j = (v_index - 1) as usize;
        if i < self.grid_points.len() && j < self.grid_points[i].len() {
            self.grid_points[i][j] = p;
        }
    }

    pub(crate) fn grid_point(&self, u_index: i32, v_index: i32) -> GpPnt {
        let i = (u_index - 1) as usize;
        let j = (v_index - 1) as usize;
        self.grid_points
            .get(i)
            .and_then(|row| row.get(j))
            .copied()
            .unwrap_or_else(GpPnt::zero)
    }

    /// `SetFrame`.
    pub(crate) fn set_frame(&mut self, u_min: f64, u_max: f64, v_min: f64, v_max: f64) {
        self.u_ind_min = 0;
        self.u_ind_max = 0;
        self.v_ind_min = 0;
        self.v_ind_max = 0;
        if self.u_params.is_empty() || self.v_params.is_empty() {
            return;
        }
        let a_len = self.u_params.len() as i32;
        for i in 1..=a_len {
            if self.u_ind_min == 0 && u_min < self.u_param(i) {
                self.u_ind_min = i;
            }
            let a_lmi = a_len - i + 1;
            if self.u_ind_max == 0 && u_max > self.u_param(a_lmi) {
                self.u_ind_max = a_lmi;
            }
        }
        if self.u_ind_min == 0 {
            self.u_ind_min = a_len + 1;
        }
        let a_len = self.v_params.len() as i32;
        for i in 1..=a_len {
            if self.v_ind_min == 0 && v_min < self.v_param(i) {
                self.v_ind_min = i;
            }
            let a_lmi = a_len - i + 1;
            if self.v_ind_max == 0 && v_max > self.v_param(a_lmi) {
                self.v_ind_max = a_lmi;
            }
        }
        if self.v_ind_min == 0 {
            self.v_ind_min = a_len + 1;
        }
    }

    pub(crate) fn nb_u_points_in_frame(&self) -> i32 {
        self.u_ind_max - self.u_ind_min + 1
    }

    pub(crate) fn nb_v_points_in_frame(&self) -> i32 {
        self.v_ind_max - self.v_ind_min + 1
    }

    pub(crate) fn point_in_frame(&self, u_index: i32, v_index: i32) -> GpPnt {
        let a_frm_u = u_index + self.u_ind_min - 1;
        let a_frm_v = v_index + self.v_ind_min - 1;
        if self.grid_points.is_empty() || a_frm_u > self.u_ind_max || a_frm_v > self.v_ind_max {
            return GpPnt::zero();
        }
        self.grid_point(a_frm_u, a_frm_v)
    }

    pub(crate) fn u_param_in_frame(&self, index: i32) -> f64 {
        let a_frm = index + self.u_ind_min - 1;
        if self.u_params.is_empty() || a_frm > self.u_ind_max {
            return f64::INFINITY;
        }
        self.u_param(a_frm)
    }

    pub(crate) fn v_param_in_frame(&self, index: i32) -> f64 {
        let a_frm = index + self.v_ind_min - 1;
        if self.v_params.is_empty() || a_frm > self.v_ind_max {
            return f64::INFINITY;
        }
        self.v_param(a_frm)
    }

    pub(crate) fn clear_grid(&mut self) {
        self.deflection = 0.0;
        self.u_ind_min = 0;
        self.u_ind_max = 0;
        self.v_ind_min = 0;
        self.v_ind_max = 0;
        self.u_params.clear();
        self.v_params.clear();
        self.grid_points.clear();
    }

    pub(crate) fn has_grid(&self) -> bool {
        !self.u_params.is_empty() && !self.v_params.is_empty()
    }
}

/// `CheckSampling` (`IntTools_BeanFaceIntersector.cxx:2596`).
pub(crate) fn check_sampling(
    curve_range: &CurveRangeSample,
    surface_range: &SurfaceRangeSample,
    curve_data: &CurveRangeLocalizeData,
    surface_data: &SurfaceRangeLocalizeData,
    diff_c: f64,
    diff_u: f64,
    diff_v: f64,
) -> (bool, bool, bool) {
    const D_LIMIT: f64 = 1000.0;
    let mut allow_c = true;
    let mut allow_u = true;
    let mut allow_v = true;

    let mut samples_nb = if curve_range.depth() == 0 {
        1
    } else {
        curve_data.nb_sample()
    };
    if (curve_data.nb_sample() as f64).powi(curve_range.depth() + 1) > D_LIMIT
        || (diff_c / samples_nb as f64) < curve_data.min_range()
    {
        allow_c = false;
    }

    samples_nb = if surface_range.depth_u() == 0 {
        1
    } else {
        surface_data.nb_sample_u()
    };
    if (surface_data.nb_sample_u() as f64).powi(surface_range.depth_u() + 1) > D_LIMIT
        || (diff_u / samples_nb as f64) < surface_data.min_range_u()
    {
        allow_u = false;
    }

    samples_nb = if surface_range.depth_v() == 0 {
        1
    } else {
        surface_data.nb_sample_v()
    };
    if (surface_data.nb_sample_v() as f64).powi(surface_range.depth_v() + 1) > D_LIMIT
        || (diff_v / samples_nb as f64) < surface_data.min_range_v()
    {
        allow_v = false;
    }

    (allow_c, allow_u, allow_v)
}
