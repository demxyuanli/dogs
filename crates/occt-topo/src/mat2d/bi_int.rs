//! `MAT2d_BiInt` (`MAT2d_BiInt.cxx:19-42`, `.hxx:...`): the pair of indices that
//! names a bisector by its two sides in the `MAT2d` graph.

/// `class MAT2d_BiInt` (`MAT2d_BiInt.hxx`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Mat2dBiInt {
    /// `i1` (`.hxx`).
    i1: i32,
    /// `i2` (`.hxx`).
    i2: i32,
}

impl Mat2dBiInt {
    /// `MAT2d_BiInt(I1, I2)` (`cxx:19-23`).
    pub fn new(i1: i32, i2: i32) -> Self {
        Self { i1, i2 }
    }

    /// `FirstIndex()` (`cxx:25-28`).
    pub fn first_index(&self) -> i32 {
        self.i1
    }

    /// `SecondIndex()` (`cxx:30-33`).
    pub fn second_index(&self) -> i32 {
        self.i2
    }

    /// `FirstIndex(I1)` (`cxx:35-38`).
    pub fn set_first_index(&mut self, i1: i32) {
        self.i1 = i1;
    }

    /// `SecondIndex(I2)` (`cxx:40-43`).
    pub fn set_second_index(&mut self, i2: i32) {
        self.i2 = i2;
    }

    /// `IsEqual(B)` (`cxx:45-49`).
    pub fn is_equal(&self, b: &Self) -> bool {
        self.i1 == b.first_index() && self.i2 == b.second_index()
    }
}
