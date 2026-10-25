//! `GeomFill_Line` — the row index list `AppBlend_AppSurf` walks.
//!
//! Source: `TKGeomAlgo/GeomFill/GeomFill_Line.hxx`, `GeomFill_Line.cxx:22-35`
//! and `GeomFill_Line.lxx:20-33` (OCCT 8.0.0). `Point(Index)` is the identity,
//! so `AppBlend_AppSurf` walks the sections by index.

/// `GeomFill_Line` (`GeomFill_Line.hxx:30-50`).
#[derive(Debug, Clone, Copy, Default)]
pub struct Line {
    /// `myNbPoints` (`hxx:46`).
    nb_points: i32,
}

impl Line {
    /// `GeomFill_Line::GeomFill_Line()` (`cxx:24-27`): `myNbPoints = 0`.
    pub fn new() -> Self {
        Self { nb_points: 0 }
    }

    /// `GeomFill_Line::GeomFill_Line(const int)` (`cxx:31-34`).
    pub fn with_nb_points(nb_points: i32) -> Self {
        Self { nb_points }
    }

    /// `GeomFill_Line::NbPoints()` (`lxx:20-23`).
    pub fn nb_points(&self) -> i32 {
        self.nb_points
    }

    /// `GeomFill_Line::Point(Index)` (`lxx:27-31`): the identity, so a section
    /// index maps to itself.
    pub fn point(&self, index: i32) -> i32 {
        index
    }
}
