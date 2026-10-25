//! `TopOpeBRep_Bipoint` — a pair of integer indices, used as the
//! `(EdgesIntersector edge, point)` key of `TopOpeBRep_EdgesFiller`.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRep/TopOpeBRep_Bipoint.hxx:26-43`,
//! `.cxx:20-43`.

/// `class TopOpeBRep_Bipoint` (`hxx:26-43`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Bipoint {
    /// `myI1`
    i1: i32,
    /// `myI2`
    i2: i32,
}

impl Bipoint {
    /// `TopOpeBRep_Bipoint()` (`cxx:20-25`): both indices `0`.
    pub fn new() -> Self {
        Self { i1: 0, i2: 0 }
    }

    /// `TopOpeBRep_Bipoint(I1, I2)` (`cxx:27-32`).
    pub fn with_indices(i1: i32, i2: i32) -> Self {
        Self { i1, i2 }
    }

    /// `I1()` (`cxx:34-41`): raises `Standard_ProgramError` when unset.
    pub fn i1(&self) -> Result<i32, String> {
        if self.i1 <= 0 {
            return Err("TopOpeBRep_Bipoint I1=0".into());
        }
        Ok(self.i1)
    }

    /// `I2()` (`cxx:43-50`): raises `Standard_ProgramError` when unset.
    pub fn i2(&self) -> Result<i32, String> {
        if self.i2 <= 0 {
            return Err("TopOpeBRep_Bipoint I2=0".into());
        }
        Ok(self.i2)
    }
}
