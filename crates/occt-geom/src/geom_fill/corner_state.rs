//! `GeomFill_CornerState` — the per-corner bookkeeping of
//! `GeomFill_ConstrainedFilling`.
//!
//! Source: `GeomFill_CornerState.hxx:23-65`, `GeomFill_CornerState.cxx:22-107`
//! (OCCT 8.0.0).
//!
//! Translation notes:
//! - `Gap()` starts at `RealLast()` (`cxx:23`), i.e. `DBL_MAX`; the port has no
//!   `RealLast` constant in `occt_core::precision`, so `f64::MAX` is used and
//!   the source is noted here.
//! - `IsToKill(Scal)` writes the scale through its out parameter; the port
//!   returns `(killed, scale)`.

/// `GeomFill_CornerState` (`GeomFill_CornerState.hxx:32-65`).
#[derive(Debug, Clone, Copy)]
pub struct CornerState {
    /// `gap` (`hxx:59`).
    gap: f64,
    /// `tgtang` (`hxx:60`).
    tgtang: f64,
    /// `isconstrained` (`hxx:61`).
    isconstrained: bool,
    /// `norang` (`hxx:62`).
    norang: f64,
    /// `scal` (`hxx:63`).
    scal: f64,
    /// `coonscnd` (`hxx:64`).
    coonscnd: bool,
}

impl Default for CornerState {
    /// `GeomFill_CornerState::GeomFill_CornerState` (`cxx:22-29`).
    fn default() -> Self {
        Self::new()
    }
}

impl CornerState {
    /// `GeomFill_CornerState::GeomFill_CornerState` (`cxx:22-29`).
    pub fn new() -> Self {
        Self {
            // `cxx:23`: `gap = RealLast()`.
            gap: f64::MAX,
            tgtang: 0.0,
            isconstrained: false,
            norang: 0.0,
            scal: 1.0,
            coonscnd: true,
        }
    }

    /// `GeomFill_CornerState::Gap` (`cxx:33-36`).
    pub fn gap(&self) -> f64 {
        self.gap
    }

    /// `GeomFill_CornerState::Gap(G)` (`cxx:40-43`).
    pub fn set_gap(&mut self, g: f64) {
        self.gap = g;
    }

    /// `GeomFill_CornerState::TgtAng` (`cxx:47-50`).
    pub fn tgt_ang(&self) -> f64 {
        self.tgtang
    }

    /// `GeomFill_CornerState::TgtAng(Ang)` (`cxx:54-57`).
    pub fn set_tgt_ang(&mut self, ang: f64) {
        self.tgtang = ang;
    }

    /// `GeomFill_CornerState::HasConstraint` (`cxx:63-66`).
    pub fn has_constraint(&self) -> bool {
        self.isconstrained
    }

    /// `GeomFill_CornerState::Constraint` (`cxx:70-73`).
    pub fn constraint(&mut self) {
        self.isconstrained = true;
    }

    /// `GeomFill_CornerState::NorAng` (`cxx:77-80`).
    pub fn nor_ang(&self) -> f64 {
        self.norang
    }

    /// `GeomFill_CornerState::NorAng(Ang)` (`cxx:84-87`).
    pub fn set_nor_ang(&mut self, ang: f64) {
        self.norang = ang;
    }

    /// `GeomFill_CornerState::IsToKill(Scal)` (`cxx:91-97`): always writes the
    /// corner's scale, but only reports a kill for a constrained corner whose
    /// Coons condition already failed.
    pub fn is_to_kill(&self) -> (bool, f64) {
        // `cxx:92`.
        let scal = self.scal;
        if !self.isconstrained {
            return (false, scal);
        }
        (!self.coonscnd, scal)
    }

    /// `GeomFill_CornerState::DoKill(Scal)` (`cxx:101-105`).
    pub fn do_kill(&mut self, scal: f64) {
        self.scal = scal;
        self.coonscnd = false;
    }
}
