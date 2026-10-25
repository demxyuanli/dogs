//! `ChFiDS_Regul` — the storage of a curve with its two faces or surfaces of
//! support.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFiDS/ChFiDS_Regul.hxx:27-54`,
//! `ChFiDS_Regul.cxx:20-96` (OCCT 8.0.0).
//!
//! The sign of `is1` / `is2` carries the "face or surface" bit: `SetS1` stores
//! `+|IS1|` for a face and `-|IS1|` for a surface (`cxx:36-48`), `IsSurface1`
//! reads `is1 < 0` (`cxx:64-67`) and `S1` returns `|is1|` (`cxx:86-89`).

/// Source: `ChFiDS_Regul.hxx:28-54`.
#[derive(Debug, Clone, Copy, Default)]
pub struct Regul {
    /// `icurv`
    icurv: i32,
    /// `is1` — signed index; negative means "surface".
    is1: i32,
    /// `is2` — signed index; negative means "surface".
    is2: i32,
}

impl Regul {
    /// `ChFiDS_Regul()` (`cxx:20-24`): all three indices start at `0`.
    pub fn new() -> Self {
        Self::default()
    }

    /// `SetCurve(IC)` (`cxx:30-33`): `std::abs`.
    pub fn set_curve(&mut self, ic: i32) {
        self.icurv = ic.abs();
    }

    /// `SetS1(IS1, IsFace = true)` (`cxx:35-49`).
    pub fn set_s1(&mut self, is1: i32, is_face: bool) {
        self.is1 = if is_face { is1.abs() } else { -is1.abs() };
    }

    /// `SetS2(IS2, IsFace = true)` (`cxx:51-64`).
    pub fn set_s2(&mut self, is2: i32, is_face: bool) {
        self.is2 = if is_face { is2.abs() } else { -is2.abs() };
    }

    /// `IsSurface1()` (`cxx:66-69`).
    pub fn is_surface1(&self) -> bool {
        self.is1 < 0
    }

    /// `IsSurface2()` (`cxx:71-74`).
    pub fn is_surface2(&self) -> bool {
        self.is2 < 0
    }

    /// `Curve()` (`cxx:78-81`).
    pub fn curve(&self) -> i32 {
        self.icurv
    }

    /// `S1()` (`cxx:84-87`): `std::abs`, so the sign never leaks to callers.
    pub fn s1(&self) -> i32 {
        self.is1.abs()
    }

    /// `S2()` (`cxx:91-94`).
    pub fn s2(&self) -> i32 {
        self.is2.abs()
    }
}
