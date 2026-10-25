//! `TopOpeBRepBuild_BlockIterator` — iterator on the elements of a block.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_BlockIterator.hxx:25-50`,
//! `.cxx:21-36`, `.lxx:22-60`.

/// `class TopOpeBRepBuild_BlockIterator` (`hxx:28`).
///
/// The default constructor leaves a cursor that is already past its end:
/// `myLower = myUpper = 0` with `myValue = 1` (`cxx:21-26`), so `More()` is
/// false until `Initialize()` resets `myValue` to `myLower` (`lxx:24-27`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BlockIterator {
    /// `myLower` (`hxx:47`).
    my_lower: i32,
    /// `myUpper` (`hxx:48`).
    my_upper: i32,
    /// `myValue` (`hxx:49`).
    my_value: i32,
}

impl BlockIterator {
    /// `TopOpeBRepBuild_BlockIterator()` (`cxx:21-26`).
    pub fn new() -> Self {
        Self {
            my_lower: 0,
            my_upper: 0,
            my_value: 1,
        }
    }

    /// `TopOpeBRepBuild_BlockIterator(Lower, Upper)` (`cxx:30-35`).
    pub fn with_bounds(lower: i32, upper: i32) -> Self {
        Self {
            my_lower: lower,
            my_upper: upper,
            my_value: lower,
        }
    }

    /// `Initialize()` (`lxx:24-27`).
    pub fn initialize(&mut self) {
        self.my_value = self.my_lower;
    }

    /// `More()` (`lxx:31-35`).
    pub fn more(&self) -> bool {
        self.my_value <= self.my_upper
    }

    /// `Next()` (`lxx:39-42`).
    pub fn next(&mut self) {
        self.my_value += 1;
    }

    /// `Value()` (`lxx:46-49`).
    pub fn value(&self) -> i32 {
        self.my_value
    }

    /// `Extent()` (`lxx:53-60`): 0 when `myLower` is 0, the block is then
    /// considered empty whatever the bounds are.
    pub fn extent(&self) -> i32 {
        if self.my_lower != 0 {
            self.my_upper - self.my_lower + 1
        } else {
            0
        }
    }
}

impl Default for BlockIterator {
    fn default() -> Self {
        Self::new()
    }
}
