//! `IntPolyh_SectionLine` (`IntPolyh_SectionLine.hxx:26-70`, `.cxx:26-139`).
//!
//! A sequence of start points forming one intersection polyline. The `.cxx`
//! keeps a leading dummy element (`Init` appends one, `NbStartPoints` is
//! `Length() - 1`), which `Prepend` relies on.

use super::start_point::IntPolyhStartPoint;

/// `IntPolyh_SectionLine` (`hxx:26`).
#[derive(Debug, Clone, Default)]
pub struct IntPolyhSectionLine {
    points: Vec<IntPolyhStartPoint>,
}

impl IntPolyhSectionLine {
    /// `Init(N)` (`cxx:38-46`): append the leading dummy only when empty.
    pub fn init(&mut self, _n: i32) {
        if self.points.is_empty() {
            self.increment_nb_start_points();
        }
    }

    /// `GetN()` (`cxx:51-55`), the raw sequence length.
    pub fn get_n(&self) -> i32 {
        self.points.len() as i32
    }

    /// `NbStartPoints()` (`cxx:58-62`).
    pub fn nb_start_points(&self) -> i32 {
        self.points.len() as i32 - 1
    }

    /// `IncrementNbStartPoints()` (`cxx:65-71`).
    pub fn increment_nb_start_points(&mut self) {
        self.points.push(IntPolyhStartPoint::default());
    }

    /// `Value(Index)` (`cxx:76-79`): the sequence is 1-based in the `.cxx`, the
    /// dummy occupies slot 0.
    pub fn value(&self, index: i32) -> IntPolyhStartPoint {
        self.points[index as usize]
    }

    /// `ChangeValue(Index)` (`cxx:84-87`).
    pub fn change_value(&mut self, index: i32) -> &mut IntPolyhStartPoint {
        &mut self.points[index as usize]
    }

    /// `Copy(Other)` (`cxx:105-114`).
    pub fn copy_from(&mut self, other: &Self) {
        self.points = other.points.clone();
    }

    /// `Prepend(SP)` (`cxx:133-136`).
    pub fn prepend(&mut self, sp: &IntPolyhStartPoint) {
        self.points.insert(0, *sp);
    }
}
