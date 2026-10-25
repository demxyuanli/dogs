//! `SearchIndex` — the first point of a line whose parameter is not below a
//! given value.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_6.cxx:94-103`
//! (OCCT 8.0.0, a file-static helper of `ChFi3d_Builder`).
//!
//! ```text
//! int NbPnt = Lin->NbPoints(), Ind;
//! for (Ind = 1; (Ind < NbPnt) && (Lin->Point(Ind).Parameter() < Value);) Ind++;
//! return Ind;
//! ```
//!
//! The loop stops at `NbPnt - 1`, so the last point is returned when no
//! earlier one qualifies, and `Ind = 1` when `Lin->NbPoints() <= 1`. The port
//! answers the 1-based index of OCCT.
//!
//! Divergences that are not heuristics:
//! - `Lin->Point(Ind)` is `Line::point(Ind)`, which answers `None` past the end;
//!   the loop bound already keeps `Ind < NbPnt`, so the port reads the point it
//!   guarded on.

use crate::brep_blend::line::Line;

/// `SearchIndex(Value, Lin)` (`cxx:94-103`): the 1-based index of the first
/// point whose parameter is `>= value`, capped at `nb_points - 1`.
pub fn search_index(value: f64, lin: &Line) -> usize {
    let nb_pnt = lin.nb_points();
    let mut ind = 1usize;
    while ind < nb_pnt && lin.point(ind).is_some_and(|p| p.parameter() < value) {
        ind += 1;
    }
    ind
}
