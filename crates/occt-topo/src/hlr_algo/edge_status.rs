//! Hidden Line status of an edge. Source: `HLRAlgo_EdgeStatus`
//! (`ModelingAlgorithms/TKHLR/HLRAlgo/HLRAlgo_EdgeStatus.hxx`,
//! `HLRAlgo_EdgeStatus.cxx:22-123`).

use super::intrv_interval::IntrvInterval;
use super::intrv_intervals::IntrvIntervals;

/// Describes the Hidden Line status of an edge: the bounds of the edge and
/// their tolerances, the `AllVisible` / `AllHidden` flags, and the sequence of
/// visible intervals. Source: `HLRAlgo_EdgeStatus`.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct HlrAlgoEdgeStatus {
    start: f64,
    end: f64,
    tol_start: f32,
    tol_end: f32,
    all_hidden: bool,
    all_visible: bool,
    visibles: IntrvIntervals,
}

impl HlrAlgoEdgeStatus {
    /// `HLRAlgo_EdgeStatus()` (`cxx:24-34`): both flags false.
    pub fn new() -> Self {
        Self::default()
    }

    /// `HLRAlgo_EdgeStatus(Start, TolStart, End, TolEnd)` (`cxx:38-49`):
    /// default visible.
    pub fn from_bounds(start: f64, tol_start: f32, end: f64, tol_end: f32) -> Self {
        let mut s = Self {
            start,
            end,
            tol_start,
            tol_end,
            all_hidden: false,
            all_visible: false,
            visibles: IntrvIntervals::new(),
        };
        s.show_all();
        s
    }

    /// `Initialize` (`cxx:53-64`): default visible.
    pub fn initialize(&mut self, start: f64, tol_start: f32, end: f64, tol_end: f32) {
        self.start = start;
        self.tol_start = tol_start;
        self.end = end;
        self.tol_end = tol_end;
        self.show_all();
    }

    /// `Bounds` (`hxx:55-62`): `(Start, TolStart, End, TolEnd)`.
    pub fn bounds(&self) -> (f64, f32, f64, f32) {
        (self.start, self.tol_start, self.end, self.tol_end)
    }

    /// `NbVisiblePart` (`cxx:66-81`).
    pub fn nb_visible_part(&self) -> usize {
        if self.all_hidden {
            0
        } else if self.all_visible {
            1
        } else {
            self.visibles.nb_intervals()
        }
    }

    /// `VisiblePart` (`cxx:85-98`): `(Start, TolStart, End, TolEnd)` of the
    /// `index`-th visible part. `index` keeps the OCCT 1-based numbering of
    /// `HLRAlgo_EdgeStatus::VisiblePart(Index)`.
    pub fn visible_part(&self, index: usize) -> Option<(f64, f32, f64, f32)> {
        if self.all_visible {
            Some(self.bounds())
        } else if index == 0 {
            None
        } else {
            self.visibles.value(index - 1).map(|i| i.bounds())
        }
    }

    /// `Hide` (`cxx:102-122`): subtracts the interval from the visible parts.
    /// `on_boundary` is unused by OCCT too (`cxx:108`).
    pub fn hide(
        &mut self,
        start: f64,
        tol_start: f32,
        end: f64,
        tol_end: f32,
        on_face: bool,
        _on_boundary: bool,
    ) {
        if !on_face {
            if self.all_visible {
                self.visibles = IntrvIntervals::from_interval(&IntrvInterval::new_with_tol(
                    self.start, self.tol_start, self.end, self.tol_end,
                ));
                self.all_visible = false;
            }
            self.visibles
                .subtract_interval(&IntrvInterval::new_with_tol(start, tol_start, end, tol_end));
            if !self.all_hidden {
                self.all_hidden = self.visibles.nb_intervals() == 0;
            }
        }
    }

    /// `HideAll` (`hxx:88-92`).
    pub fn hide_all(&mut self) {
        self.all_visible = false;
        self.all_hidden = true;
    }

    /// `ShowAll` (`hxx:94-99`).
    pub fn show_all(&mut self) {
        self.all_visible = true;
        self.all_hidden = false;
    }

    /// `AllHidden()` (`hxx:103`).
    pub fn all_hidden(&self) -> bool {
        self.all_hidden
    }

    /// `AllHidden(bool)` (`hxx:105`).
    pub fn set_all_hidden(&mut self, b: bool) {
        self.all_hidden = b;
    }

    /// `AllVisible()` (`hxx:107`).
    pub fn all_visible(&self) -> bool {
        self.all_visible
    }

    /// `AllVisible(bool)` (`hxx:109`).
    pub fn set_all_visible(&mut self, b: bool) {
        self.all_visible = b;
    }
}
