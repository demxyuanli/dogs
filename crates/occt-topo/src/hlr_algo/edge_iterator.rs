//! Iterator over the visible or hidden parts of an edge. Source:
//! `HLRAlgo_EdgeIterator`
//! (`ModelingAlgorithms/TKHLR/HLRAlgo/HLRAlgo_EdgeIterator.hxx:26-66`,
//! `HLRAlgo_EdgeIterator.cxx:26-95`, `HLRAlgo_EdgeIterator.lxx:22-71`).

use super::edge_status::HlrAlgoEdgeStatus;

/// Iterator on the visible or hidden parts of an edge. Source:
/// `HLRAlgo_EdgeIterator`.
///
/// OCCT keeps a raw `HLRAlgo_EdgeStatus*` (`hxx:57-58`); the port keeps a
/// borrowed reference, so the status must outlive the iterator.
///
/// The hidden / visible intervals are 1-based in OCCT; `i_hid` / `i_vis` keep
/// that numbering.
#[derive(Debug, Clone)]
pub struct HlrAlgoEdgeIterator<'a> {
    nb_vis: usize,
    nb_hid: usize,
    e_vis: Option<&'a HlrAlgoEdgeStatus>,
    e_hid: Option<&'a HlrAlgoEdgeStatus>,
    i_vis: usize,
    i_hid: usize,
    hid_start: f64,
    hid_end: f64,
    hid_tol_start: f32,
    hid_tol_end: f32,
}

impl Default for HlrAlgoEdgeIterator<'_> {
    /// `HLRAlgo_EdgeIterator()` (`cxx:26-38`).
    fn default() -> Self {
        Self {
            nb_vis: 0,
            nb_hid: 0,
            e_vis: None,
            e_hid: None,
            i_vis: 0,
            i_hid: 0,
            hid_start: 0.0,
            hid_end: 0.0,
            hid_tol_start: 0.0,
            hid_tol_end: 0.0,
        }
    }
}

/// OCCT would raise `Standard_OutOfRange` here; the stored zero is only
/// reachable for the `VisiblePart(1)` call of `InitHidden` on a status whose
/// `NbVisiblePart()` is `0` while `AllHidden()` is false (`cxx:52-57`).
const NO_PART: (f64, f32, f64, f32) = (0.0, 0.0, 0.0, 0.0);

impl<'a> HlrAlgoEdgeIterator<'a> {
    /// `HLRAlgo_EdgeIterator()` (`cxx:26-38`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `InitHidden(status)` (`cxx:42-63`).
    pub fn init_hidden(&mut self, status: &'a HlrAlgoEdgeStatus) {
        self.e_hid = Some(status);
        self.i_hid = 1;
        if status.all_hidden() {
            let (s, ts, e, te) = status.bounds();
            self.hid_start = s;
            self.hid_tol_start = ts;
            self.hid_end = e;
            self.hid_tol_end = te;
            self.nb_hid = 0;
        } else {
            self.nb_hid = status.nb_visible_part();
            let (s, ts, _, _) = status.bounds();
            self.hid_start = s;
            self.hid_tol_start = ts;
            let (e, te, _, _) = self.visible_part_or_none(1);
            self.hid_end = e;
            self.hid_tol_end = te;
        }
        if (self.hid_start + self.hid_tol_start as f64) >= (self.hid_end - self.hid_tol_end as f64)
            && (self.hid_end + self.hid_tol_end as f64) >= (self.hid_start - self.hid_tol_start as f64)
        {
            self.next_hidden();
        }
    }

    /// `MoreHidden()` (`lxx:22-25`).
    pub fn more_hidden(&self) -> bool {
        self.i_hid <= self.nb_hid + 1
    }

    /// `NextHidden()` (`cxx:67-94`).
    pub fn next_hidden(&mut self) {
        if self.i_hid >= self.nb_hid + 1 {
            self.i_hid += 1;
        } else {
            let (s, ts, _, _) = self.visible_part_or_none(self.i_hid);
            self.hid_start = s;
            self.hid_tol_start = ts;
            self.i_hid += 1;
            if self.i_hid == self.nb_hid + 1 {
                let status = self.e_hid.expect("next_hidden before init_hidden");
                let (_, _, e, te) = status.bounds();
                self.hid_end = e;
                self.hid_tol_end = te;
                if (self.hid_start + self.hid_tol_start as f64)
                    >= (self.hid_end - self.hid_tol_end as f64)
                    && (self.hid_end + self.hid_tol_end as f64)
                        >= (self.hid_start - self.hid_tol_start as f64)
                {
                    self.i_hid += 1;
                }
            } else {
                let (e, te, _, _) = self.visible_part_or_none(self.i_hid);
                self.hid_end = e;
                self.hid_tol_end = te;
            }
        }
    }

    /// `Hidden(Start, TolStart, End, TolEnd)` (`lxx:27-36`):
    /// `(Start, TolStart, End, TolEnd)`.
    pub fn hidden(&self) -> (f64, f32, f64, f32) {
        (self.hid_start, self.hid_tol_start, self.hid_end, self.hid_tol_end)
    }

    /// `InitVisible(status)` (`lxx:41-47`).
    pub fn init_visible(&mut self, status: &'a HlrAlgoEdgeStatus) {
        self.e_vis = Some(status);
        self.i_vis = 1;
        self.nb_vis = status.nb_visible_part();
    }

    /// `MoreVisible()` (`lxx:51-54`).
    pub fn more_visible(&self) -> bool {
        self.i_vis <= self.nb_vis
    }

    /// `NextVisible()` (`lxx:57-60`).
    pub fn next_visible(&mut self) {
        self.i_vis += 1;
    }

    /// `Visible(Start, TolStart, End, TolEnd)` (`lxx:63-71`):
    /// `(Start, TolStart, End, TolEnd)` of the current visible interval.
    pub fn visible(&self) -> (f64, f32, f64, f32) {
        match self.e_vis {
            Some(status) => status.visible_part(self.i_vis).unwrap_or(NO_PART),
            None => NO_PART,
        }
    }

    /// `VisiblePart(i, ...)` of the hidden status; the OCCT out-of-range
    /// exception becomes [`NO_PART`].
    fn visible_part_or_none(&self, index: usize) -> (f64, f32, f64, f32) {
        match self.e_hid {
            Some(status) => status.visible_part(index).unwrap_or(NO_PART),
            None => NO_PART,
        }
    }
}
