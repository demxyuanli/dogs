//! Sorted sequence of non-overlapping intervals with the boolean operations of
//! `Intrv_Intervals` (`ModelingAlgorithms/TKHLR/Intrv/Intrv_Intervals.hxx`,
//! `Intrv_Intervals.lxx`, `Intrv_Intervals.cxx:38-248`).

use super::intrv_interval::{IntrvInterval, IntrvPosition};

/// Sorted sequence of non-overlapping intervals. Source: `Intrv_Intervals`.
///
/// `Value` keeps the OCCT 1-based order in `value`: `value(i)` is
/// `Intrv_Intervals::Value(i + 1)`.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct IntrvIntervals {
    inter: Vec<IntrvInterval>,
}

impl IntrvIntervals {
    /// `Intrv_Intervals()` (`cxx:38`): a void sequence.
    pub fn new() -> Self {
        Self { inter: Vec::new() }
    }

    /// `Intrv_Intervals(const Intrv_Interval&)` (`cxx:42-45`).
    pub fn from_interval(int: &IntrvInterval) -> Self {
        Self { inter: vec![*int] }
    }

    /// `NbIntervals` (`lxx`).
    pub fn nb_intervals(&self) -> usize {
        self.inter.len()
    }

    /// `Value(Index)` with a 0-based `index`.
    pub fn value(&self, index: usize) -> Option<&IntrvInterval> {
        self.inter.get(index)
    }

    /// `Intersect(const Intrv_Interval&)` (`cxx:50-54`).
    pub fn intersect_interval(&mut self, tool: &IntrvInterval) {
        let inter = IntrvIntervals::from_interval(tool);
        self.intersect(&inter);
    }

    /// `Intersect(const Intrv_Intervals&)` (`cxx:58-63`):
    /// `XUni = this ^ Tool`, `this = this | Tool`, `this = this - XUni`.
    pub fn intersect(&mut self, tool: &IntrvIntervals) {
        let mut x_uni = self.clone();
        x_uni.x_unite(tool);
        self.unite(tool);
        self.subtract(&x_uni);
    }

    /// `Subtract(const Intrv_Interval&)` (`cxx:68-121`).
    pub fn subtract_interval(&mut self, tool: &IntrvInterval) {
        let mut index = 0usize;
        while index < self.inter.len() {
            match tool.position(&self.inter[index]) {
                // sortir
                IntrvPosition::Before => index = self.inter.len(),
                // modifier le debut, then sortir
                IntrvPosition::JustBefore => {
                    self.inter[index].cut_at_start(tool.end(), tool.tol_end());
                    index = self.inter.len();
                }
                // garder la fin, then sortir
                IntrvPosition::OverlappingAtStart | IntrvPosition::JustOverlappingAtStart => {
                    self.inter[index].set_start(tool.end(), tool.tol_end());
                    index = self.inter.len();
                }
                // detruire et continuer
                IntrvPosition::JustEnclosingAtEnd
                | IntrvPosition::Enclosing
                | IntrvPosition::Similar
                | IntrvPosition::JustEnclosingAtStart => {
                    self.inter.remove(index);
                    continue;
                }
                // garder le debut et la fin, then sortir
                IntrvPosition::Inside => {
                    let dup = self.inter[index];
                    self.inter.insert(index + 1, dup);
                    self.inter[index].set_end(tool.start(), tool.tol_start());
                    self.inter[index + 1].set_start(tool.end(), tool.tol_end());
                    index = self.inter.len();
                }
                // garder le debut, then continuer
                IntrvPosition::JustOverlappingAtEnd | IntrvPosition::OverlappingAtEnd => {
                    self.inter[index].set_end(tool.start(), tool.tol_start());
                }
                // modifier la fin, then continuer
                IntrvPosition::JustAfter => {
                    self.inter[index].cut_at_end(tool.start(), tool.tol_start());
                }
                // continuer
                IntrvPosition::After => {}
            }
            index += 1;
        }
    }

    /// `Subtract(const Intrv_Intervals&)` (`cxx:125-133`).
    pub fn subtract(&mut self, tool: &IntrvIntervals) {
        for i in 0..tool.inter.len() {
            let t = tool.inter[i];
            self.subtract_interval(&t);
        }
    }

    /// `Unite(const Intrv_Interval&)` (`cxx:137-218`).
    pub fn unite_interval(&mut self, tool: &IntrvInterval) {
        let mut inserted = false;
        let mut tins = *tool;
        let mut index = 0usize;
        while index < self.inter.len() {
            match tins.position(&self.inter[index]) {
                // inserer avant et sortir
                IntrvPosition::Before => {
                    inserted = true;
                    self.inter.insert(index, tins);
                    index = self.inter.len();
                }
                // changer le debut, then sortir
                IntrvPosition::JustBefore | IntrvPosition::OverlappingAtStart => {
                    inserted = true;
                    self.inter[index].set_start(tins.start(), tins.tol_start());
                    index = self.inter.len();
                }
                // modifier le debut et la fin, then detruire et continuer
                IntrvPosition::Similar => {
                    tins.fuse_at_start(self.inter[index].start(), self.inter[index].tol_start());
                    tins.fuse_at_end(self.inter[index].end(), self.inter[index].tol_end());
                    self.inter.remove(index);
                    continue;
                }
                // modifier la fin, then detruire et continuer
                IntrvPosition::JustEnclosingAtEnd => {
                    tins.fuse_at_end(self.inter[index].end(), self.inter[index].tol_end());
                    self.inter.remove(index);
                    continue;
                }
                // detruire et continuer
                IntrvPosition::Enclosing => {
                    self.inter.remove(index);
                    continue;
                }
                // changer le debut, modifier la fin, then detruire et continuer
                IntrvPosition::JustOverlappingAtEnd => {
                    tins.set_start(self.inter[index].start(), self.inter[index].tol_start());
                    tins.fuse_at_end(self.inter[index].end(), self.inter[index].tol_end());
                    self.inter.remove(index);
                    continue;
                }
                // modifier le debut, then sortir
                IntrvPosition::JustOverlappingAtStart => {
                    inserted = true;
                    self.inter[index].fuse_at_start(tins.start(), tins.tol_start());
                    index = self.inter.len();
                }
                // modifier le debut, then detruire et continuer
                IntrvPosition::JustEnclosingAtStart => {
                    tins.fuse_at_start(self.inter[index].start(), self.inter[index].tol_start());
                    self.inter.remove(index);
                    continue;
                }
                // sortir
                IntrvPosition::Inside => {
                    inserted = true;
                    index = self.inter.len();
                }
                // changer le debut, then detruire et continuer
                IntrvPosition::OverlappingAtEnd | IntrvPosition::JustAfter => {
                    tins.set_start(self.inter[index].start(), self.inter[index].tol_start());
                    self.inter.remove(index);
                    continue;
                }
                // continuer
                IntrvPosition::After => {}
            }
            index += 1;
        }
        if !inserted {
            self.inter.push(tins);
        }
    }

    /// `Unite(const Intrv_Intervals&)` (`cxx:222-230`).
    pub fn unite(&mut self, tool: &IntrvIntervals) {
        for i in 0..tool.inter.len() {
            let t = tool.inter[i];
            self.unite_interval(&t);
        }
    }

    /// `XUnite(const Intrv_Interval&)` (`cxx:234-238`).
    pub fn x_unite_interval(&mut self, tool: &IntrvInterval) {
        let inter = IntrvIntervals::from_interval(tool);
        self.x_unite(&inter);
    }

    /// `XUnite(const Intrv_Intervals&)` (`cxx:242-248`): symmetric difference.
    pub fn x_unite(&mut self, tool: &IntrvIntervals) {
        let mut sub2 = tool.clone();
        sub2.subtract(self);
        self.subtract(tool);
        self.unite(&sub2);
    }
}
