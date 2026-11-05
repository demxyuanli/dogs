//! `Plate_LinearXYZConstraint` — one or several constraints expressed as a
//! linear combination of `PinpointConstraint`s, keeping X, Y and Z uncoupled.
//!
//! Source: `Plate_LinearXYZConstraint.hxx:24-64`,
//! `Plate_LinearXYZConstraint.cxx:19-77`, `Plate_LinearXYZConstraint.lxx`.
//!
//! Storage note: OCCT keeps `myPPC` as `HArray1` (1-based) and `myCoef` as
//! `HArray2(1, ColLen, 1, RowLen)`. The port keeps a flat row-major `Vec<f64>`
//! with the same `(row, col)` 1-based accessors (`row` in `1..=row_len`,
//! `col` in `1..=col_len`), plus `ppc` in 1-based order.

use crate::plate::PinpointConstraint;

/// `Plate_LinearXYZConstraint` (`hxx:29-64`).
#[derive(Debug, Clone, Default)]
pub struct LinearXyzConstraint {
    ppc: Vec<PinpointConstraint>,
    /// Row-major, `row_len * col_len` entries.
    coeff: Vec<f64>,
    row_len: usize,
    col_len: usize,
}

impl LinearXyzConstraint {
    /// `Plate_LinearXYZConstraint()` (`cxx:19`): both handles null, so both
    /// lengths read as zero.
    pub fn new() -> Self {
        Self::default()
    }

    /// `Plate_LinearXYZConstraint(thePPC, theCoeff)` with
    /// `Array1<double> coeff` (`cxx:21-35`): one row, `theCoeff.Length()`
    /// columns. `Standard_DimensionMismatch` becomes `Err`.
    pub fn from_coeff_row(
        the_ppc: Vec<PinpointConstraint>,
        the_coeff: &[f64],
    ) -> Result<Self, String> {
        if the_coeff.len() != the_ppc.len() {
            return Err("Plate_LinearXYZConstraint: dimension mismatch".to_string());
        }
        Ok(Self {
            ppc: the_ppc,
            coeff: the_coeff.to_vec(),
            row_len: 1,
            col_len: the_coeff.len(),
        })
    }

    /// `Plate_LinearXYZConstraint(thePPC, Array2<double> theCoeff)`
    /// (`cxx:37-49`): `RowLength` must equal `thePPC.Length()`.
    ///
    /// `the_coeff` is given row-major with `row_len = theCoeff.RowLength()` and
    /// `col_len = theCoeff.ColLength()`; `the_coeff[(row - 1) * row_len + col - 1]`.
    pub fn from_coeff_matrix(
        the_ppc: Vec<PinpointConstraint>,
        the_coeff: &[f64],
        row_len: usize,
        col_len: usize,
    ) -> Result<Self, String> {
        if row_len != the_ppc.len() {
            return Err("Plate_LinearXYZConstraint: dimension mismatch".to_string());
        }
        if the_coeff.len() != row_len * col_len {
            return Err("Plate_LinearXYZConstraint: coefficient array size".to_string());
        }
        Ok(Self {
            ppc: the_ppc,
            coeff: the_coeff.to_vec(),
            row_len,
            col_len,
        })
    }

    /// `Plate_LinearXYZConstraint(ColLen, RowLen)` (`cxx:51-56`): `RowLen`
    /// point constraints and a zero-filled `ColLen` x `RowLen` coefficient
    /// array.
    pub fn with_size(col_len: usize, row_len: usize) -> Self {
        Self {
            ppc: vec![PinpointConstraint::default(); row_len],
            coeff: vec![0.0; row_len * col_len],
            row_len,
            col_len,
        }
    }

    /// `GetPPC()` (`hxx:44`).
    pub fn ppc(&self) -> &[PinpointConstraint] {
        &self.ppc
    }

    /// `RowLength()` of `Coeff()`.
    pub fn row_length(&self) -> usize {
        self.row_len
    }

    /// `ColLength()` of `Coeff()`.
    pub fn col_length(&self) -> usize {
        self.col_len
    }

    /// `Coeff()(row, col)` with 1-based `row` / `col`.
    pub fn coeff(&self, row: usize, col: usize) -> f64 {
        self.coeff[(row - 1) * self.row_len + (col - 1)]
    }

    /// `SetPPC(Index, Value)` (`cxx:59-62`), 1-based.
    pub fn set_ppc(&mut self, index: usize, value: PinpointConstraint) {
        self.ppc[index - 1] = value;
    }

    /// `SetCoeff(Row, Col, Value)` (`cxx:64-67`), 1-based.
    pub fn set_coeff(&mut self, row: usize, col: usize, value: f64) {
        self.coeff[(row - 1) * self.row_len + (col - 1)] = value;
    }
}
