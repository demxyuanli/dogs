//! `IntTools_MarkedRangeSet` — flagged parameter intervals for bean/face work.
//!
//! Ranges are stored as consecutive `[boundaries[i], boundaries[i+1]]`
//! intervals; inserting a range splits existing intervals and marks the overlap.

use crate::inttools_data::IntRange;

/// A sorted set of parameter ranges, each carrying an integer flag.
#[derive(Debug, Clone)]
pub(crate) struct MarkedRangeSet {
    boundaries: Vec<f64>,
    flags: Vec<i32>,
}

impl MarkedRangeSet {
    pub(crate) fn new() -> Self {
        Self {
            boundaries: Vec::new(),
            flags: Vec::new(),
        }
    }

    /// `[first, last]` covered by one range carrying `init_flag`.
    pub(crate) fn set_boundaries(&mut self, first: f64, last: f64, init_flag: i32) {
        self.boundaries = vec![first, last];
        self.flags = vec![init_flag];
    }

    pub(crate) fn len(&self) -> usize {
        self.flags.len()
    }

    pub(crate) fn is_empty(&self) -> bool {
        self.flags.is_empty()
    }

    pub(crate) fn range(&self, i: usize) -> IntRange {
        IntRange::new_unchecked(self.boundaries[i], self.boundaries[i + 1])
    }

    pub(crate) fn flag(&self, i: usize) -> i32 {
        self.flags[i]
    }

    pub(crate) fn set_flag(&mut self, i: usize, flag: i32) {
        self.flags[i] = flag;
    }

    /// Index (0-based) of the range containing `value`, using the OCCT
    /// `GetIndex(value, UseLower)` semantics; `-1` when `value` is outside the
    /// set (or exactly at the last boundary with `UseLower`).
    pub(crate) fn get_index(&self, value: f64, use_lower: bool) -> isize {
        if self.boundaries.is_empty() {
            return -1;
        }
        if (use_lower && value < self.boundaries[0])
            || (!use_lower && value <= self.boundaries[0])
        {
            return -1;
        }
        for i in 1..self.boundaries.len() {
            if (use_lower && value < self.boundaries[i])
                || (!use_lower && value <= self.boundaries[i])
            {
                return (i - 1) as isize;
            }
        }
        -1
    }

    /// Indices (0-based) of every range containing `value` (a boundary value
    /// belongs to both adjacent ranges). Port of `GetIndices`.
    pub(crate) fn get_indices(&self, value: f64) -> Vec<usize> {
        let mut out = Vec::new();
        if self.boundaries.is_empty() || value < self.boundaries[0] {
            return out;
        }
        let mut found = false;
        for i in 1..self.boundaries.len() {
            if found {
                if value >= self.boundaries[i - 1] {
                    out.push(i - 1);
                } else {
                    break;
                }
            } else if value <= self.boundaries[i] {
                out.push(i - 1);
                found = true;
            }
        }
        out
    }

    /// Insert `[first, last]` with `flag`, splitting covered ranges. Returns
    /// `false` when the boundaries do not resolve (OCCT `InsertRange`).
    pub(crate) fn insert_range(&mut self, first: f64, last: f64, flag: i32) -> bool {
        if self.boundaries.is_empty() {
            return false;
        }
        let mut idx1 = self.get_index(first, true);
        if idx1 < 0 {
            return false;
        }
        let mut idx2 = self.get_index(last, false);
        if idx2 < 0 {
            return false;
        }
        if idx2 < idx1 {
            std::mem::swap(&mut idx1, &mut idx2);
            if last < first {
                return false;
            }
        }
        let idx1 = idx1 as usize;
        let mut idx2 = idx2 as usize;
        let are_equal = idx1 == idx2;
        let prev_flag = self.flags[idx1];

        self.boundaries.insert(idx1 + 1, first);
        self.flags.insert(idx1 + 1, flag);
        idx2 += 1;
        self.boundaries.insert(idx2 + 1, last);
        if are_equal {
            self.flags.insert(idx2 + 1, prev_flag);
        } else {
            self.flags.insert(idx2, flag);
        }
        if !are_equal {
            for i in (idx1 + 1)..idx2 {
                self.flags[i] = flag;
            }
        }
        true
    }
}
