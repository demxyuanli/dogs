//! Interval with start / end tolerances. Source: `Intrv_Interval`
//! (`ModelingAlgorithms/TKHLR/Intrv/Intrv_Interval.hxx`,
//! `Intrv_Interval.lxx`, `Intrv_Interval.cxx:38-151`) and `Intrv_Position`
//! (`Intrv_Position.hxx`).

use occt_core::precision::epsilon;

/// `RealFirst()` (`Standard_Real.hxx:222-225`).
const REAL_FIRST: f64 = -f64::MAX;
/// `RealLast()` (`Standard_Real.hxx:234-237`).
const REAL_LAST: f64 = f64::MAX;

/// Relative position of one interval against another. Source:
/// `Intrv_Position.hxx`; the 13 cases are the ones classified by
/// `Intrv_Interval::Position` (`cxx:78-151`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum IntrvPosition {
    Before,
    JustBefore,
    OverlappingAtStart,
    JustEnclosingAtEnd,
    Enclosing,
    JustOverlappingAtStart,
    Similar,
    JustEnclosingAtStart,
    Inside,
    JustOverlappingAtEnd,
    OverlappingAtEnd,
    JustAfter,
    After,
}

/// `AreFused` (`Intrv_Interval.lxx:170-173`).
fn are_fused(c1: f64, t1: f32, c2: f64, t2: f32) -> bool {
    t1 as f64 + t2 as f64 >= (c1 - c2).abs()
}

/// Interval of the real line with a tolerance on each bound. Source:
/// `Intrv_Interval`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IntrvInterval {
    start: f64,
    end: f64,
    tol_start: f32,
    tol_end: f32,
}

impl Default for IntrvInterval {
    /// `Intrv_Interval::Intrv_Interval()` (`cxx:36-42`): the whole real line,
    /// tolerances `(float)Epsilon(RealFirst())` / `(float)Epsilon(RealLast())`,
    /// both `0.0f` because `RealFirst()` / `RealLast()` are already the extreme
    /// doubles.
    fn default() -> Self {
        Self {
            start: REAL_FIRST,
            end: REAL_LAST,
            tol_start: 0.0,
            tol_end: 0.0,
        }
    }
}

impl IntrvInterval {
    /// `Intrv_Interval()` (`cxx:36-42`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `Intrv_Interval(Start, End)` (`cxx:44-52`): tolerances are
    /// `Epsilon(Start)` / `Epsilon(End)`.
    pub fn from_bounds(start: f64, end: f64) -> Self {
        Self {
            start,
            end,
            tol_start: epsilon(start) as f32,
            tol_end: epsilon(end) as f32,
        }
    }

    /// `Intrv_Interval(Start, TolStart, End, TolEnd)` (`cxx:54-76`).
    pub fn new_with_tol(start: f64, tol_start: f32, end: f64, tol_end: f32) -> Self {
        let eps_start = epsilon(start) as f32;
        let eps_end = epsilon(end) as f32;
        Self {
            start,
            end,
            tol_start: if tol_start < eps_start { eps_start } else { tol_start },
            tol_end: if tol_end < eps_end { eps_end } else { tol_end },
        }
    }

    pub fn start(&self) -> f64 {
        self.start
    }

    pub fn end(&self) -> f64 {
        self.end
    }

    pub fn tol_start(&self) -> f32 {
        self.tol_start
    }

    pub fn tol_end(&self) -> f32 {
        self.tol_end
    }

    /// `Bounds` (`lxx:63-68`): `(Start, TolStart, End, TolEnd)`.
    pub fn bounds(&self) -> (f64, f32, f64, f32) {
        (self.start, self.tol_start, self.end, self.tol_end)
    }

    /// `SetStart` (`lxx:73-77`).
    pub fn set_start(&mut self, start: f64, tol_start: f32) {
        self.start = start;
        self.tol_start = tol_start;
    }

    /// `FuseAtStart` (`lxx:88-97`). No-op on an interval that starts at
    /// `RealFirst()`.
    pub fn fuse_at_start(&mut self, start: f64, tol_start: f32) {
        if self.start != REAL_FIRST {
            let a = (self.start - self.tol_start as f64).min(start - tol_start as f64);
            let b = (self.start + self.tol_start as f64).min(start + tol_start as f64);
            self.start = (a + b) / 2.0;
            self.tol_start = ((b - a) / 2.0) as f32;
        }
    }

    /// `CutAtStart` (`lxx:108-118`). No-op on an interval that starts at
    /// `RealFirst()`.
    pub fn cut_at_start(&mut self, start: f64, tol_start: f32) {
        if self.start != REAL_FIRST {
            let a = (self.start - self.tol_start as f64).max(start - tol_start as f64);
            let b = (self.start + self.tol_start as f64).max(start + tol_start as f64);
            self.start = (a + b) / 2.0;
            self.tol_start = ((b - a) / 2.0) as f32;
        }
    }

    /// `SetEnd` (`lxx:122-126`).
    pub fn set_end(&mut self, end: f64, tol_end: f32) {
        self.end = end;
        self.tol_end = tol_end;
    }

    /// `FuseAtEnd` (`lxx:138-147`). No-op on an interval that ends at
    /// `RealLast()`.
    pub fn fuse_at_end(&mut self, end: f64, tol_end: f32) {
        if self.end != REAL_LAST {
            let a = (self.end - self.tol_end as f64).max(end - tol_end as f64);
            let b = (self.end + self.tol_end as f64).max(end + tol_end as f64);
            self.end = (a + b) / 2.0;
            self.tol_end = ((b - a) / 2.0) as f32;
        }
    }

    /// `CutAtEnd` (`lxx:158-167`). No-op on an interval that ends at
    /// `RealLast()`.
    pub fn cut_at_end(&mut self, end: f64, tol_end: f32) {
        if self.end != REAL_LAST {
            let a = (self.end - self.tol_end as f64).min(end - tol_end as f64);
            let b = (self.end + self.tol_end as f64).min(end + tol_end as f64);
            self.end = (a + b) / 2.0;
            self.tol_end = ((b - a) / 2.0) as f32;
        }
    }

    /// `IsProbablyEmpty` (`lxx:179-182`).
    pub fn is_probably_empty(&self) -> bool {
        are_fused(self.start, self.tol_start, self.end, self.tol_end)
    }

    /// `Position` (`cxx:78-151`): the 13-case classification of `self`
    /// against `other`.
    pub fn position(&self, other: &Self) -> IntrvPosition {
        let my_s_min = self.start - self.tol_start as f64;
        let my_s_max = self.start + self.tol_start as f64;
        let my_e_min = self.end - self.tol_end as f64;
        let my_e_max = self.end + self.tol_end as f64;
        let ot_s_min = other.start - other.tol_start as f64;
        let ot_s_max = other.start + other.tol_start as f64;
        let ot_e_min = other.end - other.tol_end as f64;
        let ot_e_max = other.end + other.tol_end as f64;

        if my_s_max < ot_s_min {
            if my_e_max < ot_s_min {
                IntrvPosition::Before
            } else if ot_s_max >= my_e_min {
                IntrvPosition::JustBefore
            } else if my_e_max < ot_e_min {
                IntrvPosition::OverlappingAtStart
            } else if ot_e_max >= my_e_min {
                IntrvPosition::JustEnclosingAtEnd
            } else {
                IntrvPosition::Enclosing
            }
        } else if ot_s_max >= my_s_min {
            if my_e_max < ot_e_min {
                IntrvPosition::JustOverlappingAtStart
            } else if ot_e_max >= my_e_min {
                IntrvPosition::Similar
            } else {
                IntrvPosition::JustEnclosingAtStart
            }
        } else if my_s_max < ot_e_min {
            if my_e_max < ot_e_min {
                IntrvPosition::Inside
            } else if ot_e_max >= my_e_min {
                IntrvPosition::JustOverlappingAtEnd
            } else {
                IntrvPosition::OverlappingAtEnd
            }
        } else if ot_e_max >= my_s_min {
            IntrvPosition::JustAfter
        } else {
            IntrvPosition::After
        }
    }

    /// `IsBefore` (`lxx:189-192`).
    pub fn is_before(&self, other: &Self) -> bool {
        (self.tol_end as f64 + other.tol_start as f64) < other.start - self.end
    }

    /// `IsAfter` (`lxx:199-202`).
    pub fn is_after(&self, other: &Self) -> bool {
        (self.tol_start as f64 + other.tol_end as f64) < self.start - other.end
    }

    /// `IsInside` (`lxx:209-213`).
    pub fn is_inside(&self, other: &Self) -> bool {
        (self.tol_start as f64 + other.tol_start as f64) < self.start - other.start
            && (self.tol_end as f64 + other.tol_end as f64) < other.end - self.end
    }

    /// `IsEnclosing` (`lxx:220-224`).
    pub fn is_enclosing(&self, other: &Self) -> bool {
        (self.tol_start as f64 + other.tol_start as f64) < other.start - self.start
            && (self.tol_end as f64 + other.tol_end as f64) < self.end - other.end
    }

    /// `IsJustEnclosingAtStart` (`lxx:232-236`).
    pub fn is_just_enclosing_at_start(&self, other: &Self) -> bool {
        are_fused(self.start, self.tol_start, other.start, other.tol_start)
            && (self.tol_end as f64 + other.tol_end as f64) < self.end - other.end
    }

    /// `IsJustEnclosingAtEnd` (`lxx:243-247`).
    pub fn is_just_enclosing_at_end(&self, other: &Self) -> bool {
        (self.tol_start as f64 + other.tol_start as f64) < other.start - self.start
            && are_fused(other.end, other.tol_end, self.end, self.tol_end)
    }

    /// `IsJustBefore` (`lxx:254-257`).
    pub fn is_just_before(&self, other: &Self) -> bool {
        are_fused(self.end, self.tol_end, other.start, other.tol_start)
    }

    /// `IsJustAfter` (`lxx:264-267`).
    pub fn is_just_after(&self, other: &Self) -> bool {
        are_fused(other.end, other.tol_end, self.start, self.tol_start)
    }

    /// `IsOverlappingAtStart` (`lxx:274-278`).
    pub fn is_overlapping_at_start(&self, other: &Self) -> bool {
        (self.tol_start as f64 + other.tol_start as f64) < other.start - self.start
            && (self.tol_end as f64 + other.tol_start as f64) < self.end - other.start
            && (self.tol_end as f64 + other.tol_end as f64) < other.end - self.end
    }

    /// `IsOverlappingAtEnd` (`lxx:285-289`).
    pub fn is_overlapping_at_end(&self, other: &Self) -> bool {
        (self.tol_start as f64 + other.tol_start as f64) < self.start - other.start
            && (self.tol_start as f64 + other.tol_end as f64) < other.end - self.start
            && (self.tol_end as f64 + other.tol_end as f64) < self.end - other.end
    }

    /// `IsJustOverlappingAtStart` (`lxx:297-301`).
    pub fn is_just_overlapping_at_start(&self, other: &Self) -> bool {
        are_fused(self.start, self.tol_start, other.start, other.tol_start)
            && (self.tol_end as f64 + other.tol_end as f64) < other.end - self.end
    }

    /// `IsJustOverlappingAtEnd` (`lxx:309-313`).
    pub fn is_just_overlapping_at_end(&self, other: &Self) -> bool {
        (self.tol_start as f64 + other.tol_start as f64) < self.start - other.start
            && are_fused(other.end, other.tol_end, self.end, self.tol_end)
    }

    /// `IsSimilar` (`lxx:319-325`).
    pub fn is_similar(&self, other: &Self) -> bool {
        let b1 = are_fused(self.start, self.tol_start, other.start, other.tol_start);
        let b2 = are_fused(self.end, self.tol_end, other.end, other.tol_end);
        b1 && b2
    }
}
