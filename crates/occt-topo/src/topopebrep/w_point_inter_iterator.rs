//! `TopOpeBRep_WPointInterIterator` — walks the `WPoint`s of a `LineInter`.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRep/TopOpeBRep_WPointInterIterator.hxx:27-52`,
//! `.cxx:23-82`.
//!
//! Translation notes:
//! - `TopOpeBRep_PLineInter myLineInter` (`hxx:49`) becomes
//!   `Option<&'a mut LineInter>`: `CurrentWP` (`cxx:69-77`) calls the non-const
//!   `LineInter::WPoint(IW)` (`LineInter.cxx:462-476`), which writes
//!   `myCurrentWP`.
//! - `CurrentWP` raises `Standard_ProgramError` past the end (`cxx:71-74`); the
//!   port returns `Err` with the same message.

use crate::topopebrep::interf::WPointInter;
use crate::topopebrep::line_inter::LineInter;

/// `class TopOpeBRep_WPointInterIterator` (`hxx:27-52`).
pub struct WPointInterIterator<'a> {
    /// `myLineInter` (`hxx:49`).
    line: Option<&'a mut LineInter>,
    /// `myWPointIndex` (`hxx:50`).
    index: i32,
    /// `myWPointNb` (`hxx:51`).
    nb: i32,
}

impl<'a> Default for WPointInterIterator<'a> {
    /// `TopOpeBRep_WPointInterIterator()` (`cxx:23-28`).
    fn default() -> Self {
        Self::new()
    }
}

impl<'a> WPointInterIterator<'a> {
    /// `TopOpeBRep_WPointInterIterator()` (`cxx:23-28`).
    pub fn new() -> Self {
        Self {
            line: None,
            index: 0,
            nb: 0,
        }
    }

    /// `TopOpeBRep_WPointInterIterator(LI)` (`cxx:32-35`).
    pub fn of_line(line: &'a mut LineInter) -> Self {
        let mut it = Self::new();
        it.init_line(line);
        it
    }

    /// `Init(LI)` (`cxx:39-43`).
    pub fn init_line(&mut self, line: &'a mut LineInter) {
        self.line = Some(line);
        self.init();
    }

    /// `Init()` (`cxx:47-51`).
    pub fn init(&mut self) {
        self.index = 1;
        self.nb = self.line.as_deref().map(|l| l.nb_w_point()).unwrap_or(0);
    }

    /// `More()` (`cxx:55-58`).
    pub fn more(&self) -> bool {
        self.index <= self.nb
    }

    /// `Next()` (`cxx:62-65`).
    pub fn next(&mut self) {
        self.index += 1;
    }

    /// `CurrentWP()` (`cxx:69-77`).
    pub fn current_wp(&mut self) -> Result<&WPointInter, String> {
        if !self.more() {
            return Err("TopOpeBRep_WPointInterIterator::Current".to_string());
        }
        let i = self.index;
        match self.line.as_deref_mut() {
            Some(l) => Ok(l.w_point(i)),
            None => Err("TopOpeBRep_WPointInterIterator::Current".to_string()),
        }
    }

    /// `PLineInterDummy()` (`cxx:79-82`).
    pub fn p_line_inter_dummy(&self) -> Option<&LineInter> {
        self.line.as_deref()
    }
}
