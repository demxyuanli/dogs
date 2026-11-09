//! `HLRBRep_EdgeData` (`ModelingAlgorithms/TKHLR/HLRBRep/HLRBRep_EdgeData.hxx:37-172`,
//! `HLRBRep_EdgeData.cxx:22-54`, `HLRBRep_EdgeData.lxx:19-304`): one edge of the
//! exact HLR data structure - the 2D projected curve plus the visibility state
//! that the hidden-line removal mutates.
//!
//! Deviations:
//! * `myVSta`, `myVEnd` and `myTolerance` are read by OCCT only after `Set`
//!   wrote them; the C++ default constructor leaves them indeterminate and the
//!   port zero-initialises them instead of exposing uninitialised memory.

use crate::hlr_algo::{HlrAlgoEdgeStatus, MinMaxIndices};
use crate::hlrbrep_curve::HlrBRepCurve;
use occt_topo::shape::Edge;
use occt_topo::topopebrep::tool::shape_tool;

/// `HLRBRep_EdgeData::EMaskFlags` (`hxx:139-154`).
const EMASK_SELECTED: i32 = 1;
const EMASK_USED: i32 = 2;
const EMASK_RG1_LINE: i32 = 4;
const EMASK_VERTICAL: i32 = 8;
const EMASK_SIMPLE: i32 = 16;
const EMASK_OUT_LV_STA: i32 = 32;
const EMASK_OUT_LV_END: i32 = 64;
const EMASK_INT_DONE: i32 = 128;
const EMASK_CUT_AT_STA: i32 = 256;
const EMASK_CUT_AT_END: i32 = 512;
const EMASK_VER_AT_STA: i32 = 1024;
const EMASK_VER_AT_END: i32 = 2048;
const EMASK_RGN_LINE: i32 = 4096;

/// `HLRBRep_EdgeData`.
#[derive(Clone)]
pub struct HlrBRepEdgeData {
    /// `int myFlags` (`hxx:158`).
    my_flags: i32,
    /// `int myHideCount` (`hxx:159`).
    my_hide_count: i32,
    /// `int myVSta` (`hxx:160`).
    my_v_sta: i32,
    /// `int myVEnd` (`hxx:161`).
    my_v_end: i32,
    /// `HLRAlgo_EdgesBlock::MinMaxIndices myMinMax` (`hxx:162`).
    my_min_max: MinMaxIndices,
    /// `HLRAlgo_EdgeStatus myStatus` (`hxx:163`).
    my_status: HlrAlgoEdgeStatus,
    /// `HLRBRep_Curve myGeometry` (`hxx:164`).
    my_geometry: HlrBRepCurve,
    /// `float myTolerance` (`hxx:165`).
    my_tolerance: f32,
}

impl Default for HlrBRepEdgeData {
    /// `HLRBRep_EdgeData()` (`hxx:45-51`): `myFlags = 0`, `myHideCount = 0`,
    /// then `Selected(true)`.
    fn default() -> Self {
        let mut self_ = Self {
            my_flags: 0,
            my_hide_count: 0,
            my_v_sta: 0,
            my_v_end: 0,
            my_min_max: MinMaxIndices::default(),
            my_status: HlrAlgoEdgeStatus::new(),
            my_geometry: HlrBRepCurve::new(),
            my_tolerance: 0.0,
        };
        self_.set_selected(true);
        self_
    }
}

impl HlrBRepEdgeData {
    /// `HLRBRep_EdgeData()` (`hxx:45-51`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `Set(Reg1, RegN, EG, V1, V2, Out1, Out2, Cut1, Cut2, Start, TolStart,
    /// End, TolEnd)` (`cxx:25-54`).
    #[allow(clippy::too_many_arguments)]
    pub fn set(
        &mut self,
        rg1_line: bool,
        rgn_line: bool,
        edge: &Edge,
        v1: i32,
        v2: i32,
        out1: bool,
        out2: bool,
        cut1: bool,
        cut2: bool,
        start: f64,
        tol_start: f32,
        end: f64,
        tol_end: f32,
    ) {
        self.set_rg1_line(rg1_line);
        self.set_rgn_line(rgn_line);
        self.set_used(false);
        self.my_geometry.set_curve(edge);
        self.my_tolerance = shape_tool::tolerance_of(edge.shape())
            .map(|t| t as f32)
            .unwrap_or(0.0);
        self.set_v_sta(v1);
        self.set_v_end(v2);
        self.set_out_lv_sta(out1);
        self.set_out_lv_end(out2);
        self.set_cut_at_sta(cut1);
        self.set_cut_at_end(cut2);
        let res_start = self.my_geometry.resolution(tol_start as f64) as f32;
        let res_end = self.my_geometry.resolution(tol_end as f64) as f32;
        self.my_status.initialize(start, res_start, end, res_end);
    }

    fn mask(&self, m: i32) -> bool {
        (self.my_flags & m) != 0
    }

    fn set_mask(&mut self, m: i32, b: bool) {
        if b {
            self.my_flags |= m;
        } else {
            self.my_flags &= !m;
        }
    }

    /// `Selected()` / `Selected(B)` (`lxx:19-31`).
    pub fn selected(&self) -> bool {
        self.mask(EMASK_SELECTED)
    }

    pub fn set_selected(&mut self, b: bool) {
        self.set_mask(EMASK_SELECTED, b);
    }

    /// `Rg1Line()` / `Rg1Line(B)` (`lxx:35-47`).
    pub fn rg1_line(&self) -> bool {
        self.mask(EMASK_RG1_LINE)
    }

    pub fn set_rg1_line(&mut self, b: bool) {
        self.set_mask(EMASK_RG1_LINE, b);
    }

    /// `RgNLine()` / `RgNLine(B)` (`lxx:51-63`).
    pub fn rgn_line(&self) -> bool {
        self.mask(EMASK_RGN_LINE)
    }

    pub fn set_rgn_line(&mut self, b: bool) {
        self.set_mask(EMASK_RGN_LINE, b);
    }

    /// `Vertical()` / `Vertical(B)` (`lxx:67-79`).
    pub fn vertical(&self) -> bool {
        self.mask(EMASK_VERTICAL)
    }

    pub fn set_vertical(&mut self, b: bool) {
        self.set_mask(EMASK_VERTICAL, b);
    }

    /// `Simple()` / `Simple(B)` (`lxx:83-95`).
    pub fn simple(&self) -> bool {
        self.mask(EMASK_SIMPLE)
    }

    pub fn set_simple(&mut self, b: bool) {
        self.set_mask(EMASK_SIMPLE, b);
    }

    /// `OutLVSta()` / `OutLVSta(B)` (`lxx:99-111`).
    pub fn out_lv_sta(&self) -> bool {
        self.mask(EMASK_OUT_LV_STA)
    }

    pub fn set_out_lv_sta(&mut self, b: bool) {
        self.set_mask(EMASK_OUT_LV_STA, b);
    }

    /// `OutLVEnd()` / `OutLVEnd(B)` (`lxx:115-127`).
    pub fn out_lv_end(&self) -> bool {
        self.mask(EMASK_OUT_LV_END)
    }

    pub fn set_out_lv_end(&mut self, b: bool) {
        self.set_mask(EMASK_OUT_LV_END, b);
    }

    /// `CutAtSta()` / `CutAtSta(B)` (`lxx:131-143`).
    pub fn cut_at_sta(&self) -> bool {
        self.mask(EMASK_CUT_AT_STA)
    }

    pub fn set_cut_at_sta(&mut self, b: bool) {
        self.set_mask(EMASK_CUT_AT_STA, b);
    }

    /// `CutAtEnd()` / `CutAtEnd(B)` (`lxx:147-159`).
    pub fn cut_at_end(&self) -> bool {
        self.mask(EMASK_CUT_AT_END)
    }

    pub fn set_cut_at_end(&mut self, b: bool) {
        self.set_mask(EMASK_CUT_AT_END, b);
    }

    /// `VerAtSta()` / `VerAtSta(B)` (`lxx:163-175`).
    pub fn ver_at_sta(&self) -> bool {
        self.mask(EMASK_VER_AT_STA)
    }

    pub fn set_ver_at_sta(&mut self, b: bool) {
        self.set_mask(EMASK_VER_AT_STA, b);
    }

    /// `VerAtEnd()` / `VerAtEnd(B)` (`lxx:179-191`).
    pub fn ver_at_end(&self) -> bool {
        self.mask(EMASK_VER_AT_END)
    }

    pub fn set_ver_at_end(&mut self, b: bool) {
        self.set_mask(EMASK_VER_AT_END, b);
    }

    /// `AutoIntersectionDone()` / `AutoIntersectionDone(B)` (`lxx:195-207`).
    pub fn auto_intersection_done(&self) -> bool {
        self.mask(EMASK_INT_DONE)
    }

    pub fn set_auto_intersection_done(&mut self, b: bool) {
        self.set_mask(EMASK_INT_DONE, b);
    }

    /// `Used()` / `Used(B)` (`lxx:211-223`).
    pub fn used(&self) -> bool {
        self.mask(EMASK_USED)
    }

    pub fn set_used(&mut self, b: bool) {
        self.set_mask(EMASK_USED, b);
    }

    /// `HideCount()` / `HideCount(I)` (`lxx:227-237`).
    pub fn hide_count(&self) -> i32 {
        self.my_hide_count
    }

    pub fn set_hide_count(&mut self, i: i32) {
        self.my_hide_count = i;
    }

    /// `VSta()` / `VSta(I)` (`lxx:241-251`).
    pub fn v_sta(&self) -> i32 {
        self.my_v_sta
    }

    pub fn set_v_sta(&mut self, i: i32) {
        self.my_v_sta = i;
    }

    /// `VEnd()` / `VEnd(I)` (`lxx:255-265`).
    pub fn v_end(&self) -> i32 {
        self.my_v_end
    }

    pub fn set_v_end(&mut self, i: i32) {
        self.my_v_end = i;
    }

    /// `UpdateMinMax(theTotMinMax)` (`hxx:129-132`).
    pub fn update_min_max(&mut self, tot_min_max: &MinMaxIndices) {
        self.my_min_max = *tot_min_max;
    }

    /// `MinMax()` (`hxx:134`).
    pub fn min_max(&self) -> &MinMaxIndices {
        &self.my_min_max
    }

    /// `MinMax()` mutable overload.
    pub fn min_max_mut(&mut self) -> &mut MinMaxIndices {
        &mut self.my_min_max
    }

    /// `Status()` (`lxx:269-273`).
    pub fn status(&self) -> &HlrAlgoEdgeStatus {
        &self.my_status
    }

    /// `Status()` mutable overload.
    pub fn status_mut(&mut self) -> &mut HlrAlgoEdgeStatus {
        &mut self.my_status
    }

    /// `ChangeGeometry()` (`lxx:277-281`).
    pub fn change_geometry(&mut self) -> &mut HlrBRepCurve {
        &mut self.my_geometry
    }

    /// `Geometry()` (`lxx:285-289`).
    pub fn geometry(&self) -> &HlrBRepCurve {
        &self.my_geometry
    }

    /// `Curve()` (`hxx:136-137`).
    pub fn curve(&mut self) -> &mut HlrBRepCurve {
        &mut self.my_geometry
    }

    /// `Tolerance()` (`lxx:293-297`).
    pub fn tolerance(&self) -> f32 {
        self.my_tolerance
    }
}
