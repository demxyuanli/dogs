//! `Bisector_BisecCC` (`Bisector_BisecCC.cxx:60-1960`, `.hxx:49-241`): the
//! algorithmic bisector between two curves.
//!
//! The bisector is the parametric curve
//!
//! ```text
//!                1   ||P2(v)P1(u)||**2
//! F(u,v) = P1(u) - --- * ---------------- * N(u)
//!                2    (N(u).P2(v)P1(u))
//! ```
//!
//! where `v` is the solution of `H(u,v) = 0`
//!
//! ```text
//! H(u,v) = (Tu.P1P2(v))**2 * ||Tv||**2 - (Tv.P1P2(v))**2 * ||Tu||**2
//! ```
//!
//! (`cxx:596-608`). The class keeps a polygon of sampled points
//! (`Bisector_PolyBis`) to seed both the parameterisation and the iterative
//! evaluation.
//!
//! Split into
//! * `mod.rs` - the object, its accessors and the curve trait impls,
//! * `perform.rs` - `Perform` and the definition-range bookkeeping,
//! * `evaluate.rs` - `ValueAndDist` / `ValueByInt` / `Values` / `Extension`,
//! * `support.rs` - the file-local static helpers.

use std::sync::Arc;

use occt_core::gp::{GpPnt2d, GpTrsf2d, GpVec2d};

use crate::curve::Curve2d;

use super::{curve_is_closed, null_curve, BisectorBisecAnaError, BisectorCurve, BisectorPolyBis};

mod change_guide;
mod evaluate;
mod perform;
mod support;

/// `class Bisector_BisecCC` (`Bisector_BisecCC.hxx:49-241`).
pub struct BisectorBisecCC {
    /// `curve1` (`.hxx:224`).
    curve1: Arc<dyn Curve2d>,
    /// `curve2` (`.hxx:225`).
    curve2: Arc<dyn Curve2d>,
    /// `sign1` (`.hxx:226`).
    sign1: f64,
    /// `sign2` (`.hxx:227`).
    sign2: f64,
    /// `myPolygon` (`.hxx:228`).
    my_polygon: BisectorPolyBis,
    /// `startIntervals` (`.hxx:218`), 1-based like OCCT.
    start_intervals: Vec<f64>,
    /// `endIntervals` (`.hxx:219`).
    end_intervals: Vec<f64>,
    /// `currentInterval` (`.hxx:220`).
    current_interval: i32,
    /// `shiftParameter` (`.hxx:221`).
    shift_parameter: f64,
    /// `distMax` (`.hxx:222`).
    dist_max: f64,
    /// `isEmpty` (`.hxx:223`).
    is_empty: bool,
    /// `isConvex1` (`.hxx:229`).
    is_convex1: bool,
    /// `isConvex2` (`.hxx:230`).
    is_convex2: bool,
    /// `extensionStart` (`.hxx:231`).
    extension_start: bool,
    /// `extensionEnd` (`.hxx:232`).
    extension_end: bool,
    /// `pointStart` (`.hxx:233`).
    point_start: GpPnt2d,
    /// `pointEnd` (`.hxx:234`).
    point_end: GpPnt2d,
}

impl Default for BisectorBisecCC {
    fn default() -> Self {
        Self::new()
    }
}

impl BisectorBisecCC {
    /// `Bisector_BisecCC()` (`cxx:64-77`).
    pub fn new() -> Self {
        Self {
            curve1: null_curve(),
            curve2: null_curve(),
            sign1: 0.0,
            sign2: 0.0,
            my_polygon: BisectorPolyBis::new(),
            start_intervals: Vec::new(),
            end_intervals: Vec::new(),
            current_interval: 0,
            shift_parameter: 0.0,
            dist_max: 0.0,
            is_empty: true,
            is_convex1: false,
            is_convex2: false,
            extension_start: false,
            extension_end: false,
            point_start: GpPnt2d::zero(),
            point_end: GpPnt2d::zero(),
        }
    }

    /// `Bisector_BisecCC(Cu1, Cu2, Side1, Side2, Origin, DistMax)` (`cxx:81-88`).
    pub fn with_curves(
        cu1: &dyn Curve2d,
        cu2: &dyn Curve2d,
        side1: f64,
        side2: f64,
        origin: &GpPnt2d,
        dist_max: f64,
    ) -> Self {
        let mut b = Self::new();
        b.perform(cu1, cu2, side1, side2, origin, dist_max);
        b
    }

    // --- 1-based interval sequences (NCollection_Sequence<double>) -----------

    fn start_first(&self) -> f64 {
        *self.start_intervals.first().unwrap_or(&0.0)
    }

    fn start_value(&self, k: usize) -> f64 {
        self.start_intervals[k - 1]
    }

    fn end_last(&self) -> f64 {
        *self.end_intervals.last().unwrap_or(&0.0)
    }

    fn end_value(&self, k: usize) -> f64 {
        self.end_intervals[k - 1]
    }

    /// `IsExtendAtStart()` (`cxx:332-335`).
    pub fn is_extend_at_start(&self) -> bool {
        self.extension_start
    }

    /// `IsExtendAtEnd()` (`cxx:339-342`).
    pub fn is_extend_at_end(&self) -> bool {
        self.extension_end
    }

    /// `IsEmpty()` (`cxx:346-349`).
    pub fn is_empty(&self) -> bool {
        self.is_empty
    }

    /// `ReversedParameter(U)` (`cxx:358-361`).
    pub fn reversed_parameter(&self, u: f64) -> f64 {
        self.last_parameter() + self.first_parameter() - u
    }

    /// `Copy()` (`cxx:366-390`).
    pub fn copy_bisector(&self) -> Self {
        let mut c = Self::new();
        c.curve1 = Arc::from(self.curve1.clone_dyn());
        c.curve2 = Arc::from(self.curve2.clone_dyn());
        c.sign1 = self.sign1;
        c.sign2 = self.sign2;
        c.is_convex1 = self.is_convex1;
        c.is_convex2 = self.is_convex2;
        c.my_polygon = self.my_polygon.clone();
        c.is_empty = self.is_empty;
        c.dist_max = self.dist_max;
        c.start_intervals = self.start_intervals.clone();
        c.end_intervals = self.end_intervals.clone();
        c.current_interval = self.current_interval;
        c.shift_parameter = self.shift_parameter;
        c.extension_start = self.extension_start;
        c.extension_end = self.extension_end;
        c.point_start = self.point_start;
        c.point_end = self.point_end;
        c
    }

    /// `Transform(T)` (`cxx:453-460`).
    pub fn transform_curve(&mut self, t: &GpTrsf2d) {
        let mut c1 = self.curve1.clone_dyn();
        c1.transform(t);
        self.curve1 = Arc::from(c1);
        let mut c2 = self.curve2.clone_dyn();
        c2.transform(t);
        self.curve2 = Arc::from(c2);
        self.my_polygon.transform(t);
        self.point_start.transform(t);
        self.point_end.transform(t);
    }

    /// `IsCN(N)` (`cxx:464-467`).
    pub fn is_cn(&self, n: i32) -> bool {
        self.curve1.continuity() as i32 >= n + 1 && self.curve2.continuity() as i32 >= n + 1
    }

    /// `FirstParameter()` (`cxx:472-475`); `0` when the sequence is empty.
    pub fn first_parameter(&self) -> f64 {
        self.start_first()
    }

    /// `LastParameter()` (`cxx:479-482`); `0` when the sequence is empty.
    pub fn last_parameter(&self) -> f64 {
        self.end_last()
    }

    /// `Continuity()` (`cxx:486-502`).
    pub fn continuity(&self) -> u8 {
        match self.curve1.continuity() {
            2 => 0,
            4 => 2,
            5 => 4,
            6 => 6,
            _ => 0,
        }
    }

    /// `NbIntervals()` (`cxx:506-509`).
    pub fn nb_intervals(&self) -> i32 {
        self.start_intervals.len() as i32
    }

    /// `IntervalFirst(Index)` (`cxx:513-516`).
    pub fn interval_first(&self, index: i32) -> Result<f64, BisectorBisecAnaError> {
        if index < 1 || (index as usize) > self.start_intervals.len() {
            return Err(BisectorBisecAnaError::OutOfRange);
        }
        Ok(self.start_value(index as usize))
    }

    /// `IntervalLast(Index)` (`cxx:520-523`).
    pub fn interval_last(&self, index: i32) -> Result<f64, BisectorBisecAnaError> {
        if index < 1 || (index as usize) > self.end_intervals.len() {
            return Err(BisectorBisecAnaError::OutOfRange);
        }
        Ok(self.end_value(index as usize))
    }

    /// `IntervalContinuity()` (`cxx:527-543`).
    pub fn interval_continuity(&self) -> u8 {
        self.continuity()
    }

    /// `IsClosed()` (`cxx:547-559`).
    pub fn is_closed(&self) -> bool {
        if curve_is_closed(&*self.curve1) {
            if self.start_first() == self.curve1.first_parameter()
                && self.end_last() == self.curve1.last_parameter()
            {
                return true;
            }
        }
        false
    }

    /// `Curve(I)` (`cxx:1514-1528`).
    pub fn curve(&self, i: i32) -> Result<&Arc<dyn Curve2d>, BisectorBisecAnaError> {
        match i {
            1 => Ok(&self.curve1),
            2 => Ok(&self.curve2),
            _ => Err(BisectorBisecAnaError::OutOfRange),
        }
    }

    /// `Curve(I, C)` (`cxx:1611-1623`).
    pub fn set_curve(&mut self, i: i32, c: &dyn Curve2d) -> Result<(), BisectorBisecAnaError> {
        match i {
            1 => self.curve1 = Arc::from(c.clone_dyn()),
            2 => self.curve2 = Arc::from(c.clone_dyn()),
            _ => return Err(BisectorBisecAnaError::OutOfRange),
        }
        Ok(())
    }

    /// `Sign(I)` (`Bisector_BisecCC.hxx:74`).
    pub fn sign(&self, i: i32) -> Result<f64, BisectorBisecAnaError> {
        match i {
            1 => Ok(self.sign1),
            2 => Ok(self.sign2),
            _ => Err(BisectorBisecAnaError::OutOfRange),
        }
    }

    /// `Sign(I, S)` (`cxx:1627-1643`).
    pub fn set_sign(&mut self, i: i32, s: f64) -> Result<(), BisectorBisecAnaError> {
        match i {
            1 => self.sign1 = s,
            2 => self.sign2 = s,
            _ => return Err(BisectorBisecAnaError::OutOfRange),
        }
        Ok(())
    }

    /// `IsConvex(I)` (`Bisector_BisecCC.hxx:70`).
    pub fn is_convex_of(&self, i: i32) -> Result<bool, BisectorBisecAnaError> {
        match i {
            1 => Ok(self.is_convex1),
            2 => Ok(self.is_convex2),
            _ => Err(BisectorBisecAnaError::OutOfRange),
        }
    }

    /// `IsConvex(I, IsConvex)` (`cxx:1662-1674`).
    pub fn set_is_convex(&mut self, i: i32, value: bool) -> Result<(), BisectorBisecAnaError> {
        match i {
            1 => self.is_convex1 = value,
            2 => self.is_convex2 = value,
            _ => return Err(BisectorBisecAnaError::OutOfRange),
        }
        Ok(())
    }

    /// `Polygon()` (`cxx:1559-1562`).
    pub fn polygon(&self) -> &BisectorPolyBis {
        &self.my_polygon
    }

    /// `Polygon(P)` (`cxx:1647-1650`).
    pub fn set_polygon(&mut self, p: &BisectorPolyBis) {
        self.my_polygon = p.clone();
    }

    /// `DistMax()` (`Bisector_BisecCC.hxx:90`).
    pub fn dist_max(&self) -> f64 {
        self.dist_max
    }

    /// `DistMax(D)` (`cxx:1653-1656`).
    pub fn set_dist_max(&mut self, d: f64) {
        self.dist_max = d;
    }

    /// `IsEmpty(IsEmpty)` (`cxx:1678-1681`).
    pub fn set_is_empty(&mut self, value: bool) {
        self.is_empty = value;
    }

    /// `ExtensionStart(ExtensionStart)` (`cxx:1685-1688`).
    pub fn set_extension_start(&mut self, value: bool) {
        self.extension_start = value;
    }

    /// `ExtensionEnd(ExtensionEnd)` (`cxx:1692-1695`).
    pub fn set_extension_end(&mut self, value: bool) {
        self.extension_end = value;
    }

    /// `PointStart(Point)` (`cxx:1699-1702`).
    pub fn set_point_start(&mut self, p: &GpPnt2d) {
        self.point_start = *p;
    }

    /// `PointEnd(Point)` (`cxx:1706-1709`).
    pub fn set_point_end(&mut self, p: &GpPnt2d) {
        self.point_end = *p;
    }

    /// `StartIntervals(StartIntervals)` (`cxx:1713-1716`).
    pub fn set_start_intervals(&mut self, values: &[f64]) {
        self.start_intervals = values.to_vec();
    }

    /// `EndIntervals(EndIntervals)` (`cxx:1720-1723`).
    pub fn set_end_intervals(&mut self, values: &[f64]) {
        self.end_intervals = values.to_vec();
    }

    /// `FirstParameter(U)` (`cxx:1727-1730`): appends to `startIntervals`.
    pub fn append_first_parameter(&mut self, u: f64) {
        self.start_intervals.push(u);
    }

    /// `LastParameter(U)` (`cxx:1734-1737`): appends to `endIntervals`.
    pub fn append_last_parameter(&mut self, u: f64) {
        self.end_intervals.push(u);
    }

    /// `LinkBisCurve(U)` (`cxx:1532-1535`).
    pub fn link_bis_curve(&self, u: f64) -> f64 {
        u - self.shift_parameter
    }

    /// `LinkCurveBis(U)` (`cxx:1539-1542`).
    pub fn link_curve_bis(&self, u: f64) -> f64 {
        u + self.shift_parameter
    }

    /// `Parameter(P)` (`cxx:1565-1583`).
    pub fn parameter(&self, p: &GpPnt2d) -> f64 {
        if p.is_equal(&self.eval_d0(self.first_parameter()), occt_core::precision::CONFUSION) {
            self.first_parameter()
        } else if p
            .is_equal(&self.eval_d0(self.last_parameter()), occt_core::precision::CONFUSION)
        {
            self.last_parameter()
        } else {
            let mut u_on_curve = 0.0;
            support::proj_on_curve(p, &*self.curve1, &mut u_on_curve);
            u_on_curve
        }
    }
}

impl Curve2d for BisectorBisecCC {
    fn d0(&self, u: f64) -> GpPnt2d {
        self.eval_d0(u)
    }

    fn d1(&self, u: f64) -> (GpPnt2d, GpVec2d) {
        self.eval_d1(u)
    }

    fn d2(&self, u: f64) -> (GpPnt2d, GpVec2d, GpVec2d) {
        self.eval_d2(u)
    }

    fn d3(&self, u: f64) -> (GpPnt2d, GpVec2d, GpVec2d, GpVec2d) {
        self.eval_d3(u)
    }

    fn eval_dn(&self, u: f64, n: i32) -> GpVec2d {
        BisectorBisecCC::eval_dn(self, u, n)
    }

    fn first_parameter(&self) -> f64 {
        BisectorBisecCC::first_parameter(self)
    }

    fn last_parameter(&self) -> f64 {
        BisectorBisecCC::last_parameter(self)
    }

    /// `Bisector_BisecCC::IsPeriodic` (`cxx:563-566`): always false.
    fn is_periodic(&self) -> bool {
        false
    }

    fn continuity(&self) -> u8 {
        BisectorBisecCC::continuity(self)
    }

    fn transform(&mut self, t: &GpTrsf2d) {
        self.transform_curve(t);
    }

    /// `Bisector_BisecCC::Reverse` (`cxx:353-356`) raises
    /// `Standard_NotImplemented`; the port panics for the same reason.
    fn reverse(&mut self) {
        panic!("Standard_NotImplemented: Bisector_BisecCC::Reverse");
    }

    fn reversed_parameter(&self, u: f64) -> f64 {
        BisectorBisecCC::reversed_parameter(self, u)
    }

    fn clone_dyn(&self) -> Box<dyn Curve2d> {
        Box::new(self.copy_bisector())
    }

    /// `STANDARD_TYPE(Bisector_BisecCC)` downcast (`Bisector_Inter.cxx:290`).
    fn bisector_bisec_cc(&self) -> Option<&BisectorBisecCC> {
        Some(self)
    }
}

impl BisectorCurve for BisectorBisecCC {
    fn is_extend_at_start(&self) -> bool {
        BisectorBisecCC::is_extend_at_start(self)
    }

    fn is_extend_at_end(&self) -> bool {
        BisectorBisecCC::is_extend_at_end(self)
    }

    fn nb_intervals(&self) -> i32 {
        BisectorBisecCC::nb_intervals(self)
    }

    fn interval_first(&self, index: i32) -> Result<f64, BisectorBisecAnaError> {
        BisectorBisecCC::interval_first(self, index)
    }

    fn interval_last(&self, index: i32) -> Result<f64, BisectorBisecAnaError> {
        BisectorBisecCC::interval_last(self, index)
    }
}
