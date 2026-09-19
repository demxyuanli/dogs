//! Port of `Intf_SectionLine`
//! (`src/ModelingAlgorithms/TKGeomAlgo/Intf/Intf_SectionLine.hxx/.cxx/.lxx`).

use super::section_point::IntfSectionPoint;

/// `Intf_SectionLine` (`Intf_SectionLine.hxx:31-98`).
#[derive(Clone, Debug)]
pub struct IntfSectionLine {
    /// `myPoints` (`hxx:95`)
    my_points: Vec<IntfSectionPoint>,
    /// `closed` (`hxx:96`)
    closed: bool,
}

impl Default for IntfSectionLine {
    /// `Intf_SectionLine()` (`Intf_SectionLine.cxx:22-25`).
    fn default() -> Self {
        Self::new()
    }
}

impl IntfSectionLine {
    /// `Intf_SectionLine()` (`Intf_SectionLine.cxx:22-25`).
    pub fn new() -> Self {
        Self {
            my_points: Vec::new(),
            closed: false,
        }
    }

    /// `Intf_SectionLine(const Intf_SectionLine&)` (`Intf_SectionLine.cxx:29-33`).
    pub fn from_other(other: &Self) -> Self {
        Self {
            my_points: other.my_points.clone(),
            closed: false,
        }
    }

    /// `NumberOfPoints()` (`SectionLine.lxx:22-25`).
    pub fn number_of_points(&self) -> usize {
        self.my_points.len()
    }

    /// `GetPoint(Index)` (`Intf_SectionLine.cxx:79-82`). OCCT uses 1-based
    /// indexing.
    pub fn get_point(&self, index: usize) -> &IntfSectionPoint {
        &self.my_points[index - 1]
    }

    /// `IsClosed()` (`Intf_SectionLine.cxx:85-95`).
    pub fn is_closed(&self) -> bool {
        match (self.my_points.first(), self.my_points.last()) {
            (Some(f), Some(l)) => f == l,
            _ => false,
        }
    }

    /// `Contains(ThePI)` (`Intf_SectionLine.cxx:98-109`).
    pub fn contains(&self, the_pi: &IntfSectionPoint) -> bool {
        self.my_points.iter().any(|p| the_pi.is_equal(p))
    }

    /// `IsEnd(ThePI)` (`Intf_SectionLine.cxx:113-124`).
    pub fn is_end(&self, the_pi: &IntfSectionPoint) -> usize {
        if self.my_points.first().map(|f| f.is_equal(the_pi)) == Some(true) {
            return 1;
        }
        if self.my_points.last().map(|l| l.is_equal(the_pi)) == Some(true) {
            return self.my_points.len();
        }
        0
    }

    /// `IsEqual(Other)` (`Intf_SectionLine.cxx:128-142`).
    pub fn is_equal(&self, other: &Self) -> bool {
        if self.my_points.len() != other.my_points.len() {
            return false;
        }
        for i in 0..self.my_points.len() {
            if !self.my_points[i].is_equal(&other.my_points[i]) {
                return false;
            }
        }
        true
    }

    /// `Append(Pi)` (`Intf_SectionLine.cxx:37-40`).
    pub fn append_point(&mut self, pi: &IntfSectionPoint) {
        self.my_points.push(*pi);
    }

    /// `Append(LS)` (`Intf_SectionLine.cxx:44-47`).
    pub fn append_line(&mut self, ls: &Self) {
        self.my_points.extend_from_slice(&ls.my_points);
    }

    /// `Prepend(Pi)` (`Intf_SectionLine.cxx:52-55`).
    pub fn prepend_point(&mut self, pi: &IntfSectionPoint) {
        self.my_points.insert(0, *pi);
    }

    /// `Prepend(LS)` (`Intf_SectionLine.cxx:59-62`).
    pub fn prepend_line(&mut self, ls: &Self) {
        let mut new_points = ls.my_points.clone();
        new_points.extend_from_slice(&self.my_points);
        self.my_points = new_points;
    }

    /// `Reverse()` (`Intf_SectionLine.cxx:66-69`).
    pub fn reverse(&mut self) {
        self.my_points.reverse();
    }

    /// `Close()` (`Intf_SectionLine.cxx:73-76`).
    pub fn close(&mut self) {
        self.closed = true;
    }

    /// Access to the stored points; the port uses this where OCCT exposes the
    /// protected `myPoints` sequence indirectly.
    pub fn points(&self) -> &[IntfSectionPoint] {
        &self.my_points
    }
}

impl PartialEq for IntfSectionLine {
    fn eq(&self, other: &Self) -> bool {
        self.is_equal(other)
    }
}
