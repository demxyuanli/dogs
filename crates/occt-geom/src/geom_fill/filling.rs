//! `GeomFill_Filling` — the abstract pole grid shared by the three filling
//! styles.
//!
//! Source: `GeomFill_Filling.hxx`/`GeomFill_Filling.cxx:22-61`.
//!
//! Translation notes:
//! - OCCT stores the grid in a `NCollection_HArray2<gp_Pnt>` indexed
//!   `(Row = u, Col = v)`; the port stores `Vec<Vec<GpPnt>>` with
//!   `poles[i][j]` (0-based `u = i`, `v = j`), which is the grid
//!   `GeomBSplineSurface::from_poles_knots_mults` and `GeomBezierSurface`
//!   take.
//! - `NbUPoles()` / `NbVPoles()` read `ColLength()` / `RowLength()`, whose
//!   names are swapped in `NCollection_Array2` (`RowLength()` is the number of
//!   columns, `ColLength()` the number of rows,
//!   `NCollection_Array2.hxx:198-201`); the port exposes the u / v counts
//!   directly.

use occt_core::gp::GpPnt;

/// `GeomFill_Filling` (`GeomFill_Filling.hxx:26-60`).
pub trait Filling {
    /// `GeomFill_Filling::NbUPoles`.
    fn nb_u_poles(&self) -> usize;
    /// `GeomFill_Filling::NbVPoles`.
    fn nb_v_poles(&self) -> usize;
    /// `GeomFill_Filling::isRational`.
    fn is_rational(&self) -> bool;
    /// `GeomFill_Filling::Poles`: the `[u][v]` pole grid.
    fn poles(&self) -> &[Vec<GpPnt>];
    /// `GeomFill_Filling::Weights`: `None` where OCCT dereferences a null
    /// weight handle (a polynomial filling).
    fn weights(&self) -> Option<&[Vec<f64>]>;
}
