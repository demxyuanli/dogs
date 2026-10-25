//! `IntPatch_LineConstructor` — splits the raw intersection lines between two
//! surfaces into segments that stay inside the two faces' domains.
//!
//! Source: `ModelingAlgorithms/TKGeomAlgo/IntPatch/IntPatch_LineConstructor.hxx:33-52`,
//! `.cxx:1-2085`.
//!
//! Translation notes:
//! - `slin` (`hxx:50`) is a `Vec<GeomIntLine>`.
//! - `Perform` (`cxx:1433-2085`) is ported in `perform.rs`; the walking-line
//!   to restriction-line conversion it drives (`TestWLineIsARLine`,
//!   `TestIfWLineIsRestriction`, `ProjectOnArc`, `TestWLineToRLine`,
//!   `cxx:949-1432`) lives in `wline_to_rline.rs`.
//! - `Line(index)` (`cxx:882-886`) is 1-based in OCCT; the port answers
//!   `None` outside `1..=nb_lines`.
//! - The domains are `&mut` because `Classify` / `IsThePointOn` /
//!   `Init` / `Next` mutate OCCT's handle targets (the lazy `myFClass2d` and the
//!   arc iterator).

pub mod add_line;
pub mod helpers;
pub mod perform;
pub mod wline_to_rline;

use crate::geom_int::GeomIntLine;

/// `class IntPatch_LineConstructor` (`hxx:33-52`).
pub struct LineConstructor {
    /// `slin` (`hxx:50`).
    slin: Vec<GeomIntLine>,
}

impl LineConstructor {
    /// `IntPatch_LineConstructor(const int mode)` (`cxx:412`): `mode` is
    /// ignored by OCCT.
    pub fn new(_mode: i32) -> Self {
        Self { slin: Vec::new() }
    }

    /// `Line(index)` (`cxx:882-886`), 1-based.
    pub fn line(&self, index: i32) -> Option<&GeomIntLine> {
        if index < 1 {
            return None;
        }
        self.slin.get((index - 1) as usize)
    }

    /// `NbLines()` (`cxx:889-892`).
    pub fn nb_lines(&self) -> i32 {
        self.slin.len() as i32
    }

    /// The constructed lines, in `slin` order.
    pub fn lines(&self) -> &[GeomIntLine] {
        &self.slin
    }

    /// Mutable access to `slin` for [`perform`].
    pub(crate) fn slin_mut(&mut self) -> &mut Vec<GeomIntLine> {
        &mut self.slin
    }

    /// `Perform(SL, L, S1, D1, S2, D2, Tol)` (`cxx:1433-2085`).
    #[allow(clippy::too_many_arguments)]
    pub fn perform(
        &mut self,
        sl: &[GeomIntLine],
        l: &GeomIntLine,
        s1: &dyn occt_geom::Surface,
        d1: &mut crate::brep_top_adaptor_topol_tool::TopolTool,
        s2: &dyn occt_geom::Surface,
        d2: &mut crate::brep_top_adaptor_topol_tool::TopolTool,
        tol: f64,
    ) {
        perform::perform(self, sl, l, s1, d1, s2, d2, tol);
    }
}
