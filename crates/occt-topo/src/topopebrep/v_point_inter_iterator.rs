//! `TopOpeBRep_VPointInterIterator` — walks the `VPoint`s of a `LineInter`.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRep/TopOpeBRep_VPointInterIterator.hxx:27-57`,
//! `.cxx:23-137`.
//!
//! Translation notes:
//! - `TopOpeBRep_PLineInter myLineInter` is a `TopOpeBRep_LineInter*`
//!   (`TopOpeBRep_PLineInter.hxx`); the iterator reads and writes it, so the
//!   port stores `Option<&'a mut LineInter>`. `ChangeCurrentVP`
//!   (`cxx:113-121`) is the write path; it is the reason the borrow is mutable.
//! - The index is 1-based over `myVPointIndex <= myVPointNb` (`cxx:53`, `:72-76`),
//!   exactly as `LineInter::VPoint(I)` expects.
//! - `CurrentVP` / `ChangeCurrentVP` / `CurrentVPIndex` raise
//!   `Standard_ProgramError` past the end (`cxx:101-105`, `:113-117`, `:125-129`);
//!   the port returns `Err` with the same message.

use crate::topopebrep::line_inter::LineInter;
use crate::topopebrep::interf::VPointInter;

/// `class TopOpeBRep_VPointInterIterator` (`hxx:27-57`).
pub struct VPointInterIterator<'a> {
    /// `myLineInter` (`hxx:53`).
    line: Option<&'a mut LineInter>,
    /// `myVPointIndex` (`hxx:54`).
    index: i32,
    /// `myVPointNb` (`hxx:55`).
    nb: i32,
    /// `mycheckkeep` (`hxx:56`): skip the `VPoint`s whose `Keep()` is false.
    check_keep: bool,
}

impl<'a> Default for VPointInterIterator<'a> {
    /// `TopOpeBRep_VPointInterIterator()` (`cxx:23-26`).
    fn default() -> Self {
        Self::new()
    }
}

impl<'a> VPointInterIterator<'a> {
    /// `TopOpeBRep_VPointInterIterator()` (`cxx:23-26`).
    pub fn new() -> Self {
        Self {
            line: None,
            index: 0,
            nb: 0,
            check_keep: false,
        }
    }

    /// `TopOpeBRep_VPointInterIterator(LI)` (`cxx:33-36`): `Init(LI)` with
    /// `checkkeep = false`.
    pub fn of_line(line: &'a mut LineInter) -> Self {
        let mut it = Self::new();
        it.init_line(line, false);
        it
    }

    /// `Init(LI, checkkeep = false)` (`cxx:40-45`).
    pub fn init_line(&mut self, line: &'a mut LineInter, check_keep: bool) {
        self.line = Some(line);
        self.check_keep = check_keep;
        self.init();
    }

    /// `Init()` (`cxx:49-70`): the `mycheckkeep` arm walks forward until the
    /// first kept `VPoint`.
    pub fn init(&mut self) {
        self.index = 1;
        self.nb = self.line.as_deref().map(|l| l.nb_v_point()).unwrap_or(0);
        if self.check_keep {
            while self.more() {
                let keep = self
                    .current_vp()
                    .map(|vp| vp.keep())
                    .unwrap_or(false);
                if keep {
                    break;
                } else {
                    self.index += 1;
                }
            }
        }
    }

    /// `More()` (`cxx:72-75`).
    pub fn more(&self) -> bool {
        self.index <= self.nb
    }

    /// `Next()` (`cxx:79-94`).
    pub fn next(&mut self) {
        self.index += 1;
        if self.check_keep {
            while self.more() {
                let keep = self
                    .current_vp()
                    .map(|vp| vp.keep())
                    .unwrap_or(false);
                if keep {
                    break;
                } else {
                    self.index += 1;
                }
            }
        }
    }

    /// `CurrentVP()` (`cxx:101-110`).
    pub fn current_vp(&self) -> Result<&VPointInter, String> {
        if !self.more() {
            return Err("TopOpeBRep_VPointInterIterator::CurrentVP".to_string());
        }
        self.line
            .as_deref()
            .and_then(|l| l.v_point(self.index))
            .ok_or_else(|| "TopOpeBRep_VPointInterIterator::CurrentVP".to_string())
    }

    /// `ChangeCurrentVP()` (`cxx:113-122`).
    pub fn change_current_vp(&mut self) -> Result<&mut VPointInter, String> {
        if !self.more() {
            return Err("TopOpeBRep_VPointInterIterator::ChangeCurrentVP".to_string());
        }
        let i = self.index;
        self.line
            .as_deref_mut()
            .and_then(|l| l.change_v_point(i))
            .ok_or_else(|| "TopOpeBRep_VPointInterIterator::ChangeCurrentVP".to_string())
    }

    /// `CurrentVPIndex()` (`cxx:125-132`).
    pub fn current_vp_index(&self) -> Result<i32, String> {
        if !self.more() {
            return Err("TopOpeBRep_VPointInterIterator::CurrentVPIndex".to_string());
        }
        Ok(self.index)
    }

    /// `PLineInterDummy()` (`cxx:134-137`).
    pub fn p_line_inter_dummy(&self) -> Option<&LineInter> {
        self.line.as_deref()
    }
}
