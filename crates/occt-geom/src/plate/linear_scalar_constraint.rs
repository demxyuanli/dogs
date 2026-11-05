//! `Plate_LinearScalarConstraint` — one or several constraints expressed as a
//! linear combination of the X, Y and Z components of a set of
//! `PinpointConstraint`s.
//!
//! Source: `Plate_LinearScalarConstraint.hxx:24-67`,
//! `Plate_LinearScalarConstraint.cxx:19-80`, `Plate_LinearScalarConstraint.lxx`.
//!
//! Storage note: OCCT keeps `myCoef` as `HArray2(1, ColLen, 1, RowLen)` of
//! `gp_XYZ`. The port keeps a flat row-major `Vec<GpXyz>` with the same
//! 1-based `(row, col)` accessors.

use occt_core::gp::GpXyz;

use crate::plate::PinpointConstraint;

/// `Plate_LinearScalarConstraint` (`hxx:29-67`).
#[derive(Debug, Clone, Default)]
pub struct LinearScalarConstraint {
    ppc: Vec<PinpointConstraint>,
    /// Row-major, `row_len * col_len` entries.
    coeff: Vec<GpXyz>,
    row_len: usize,
    col_len: usize,
}

impl LinearScalarConstraint {
    /// `Plate_LinearScalarConstraint()` (`cxx:19`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `Plate_LinearScalarConstraint(PPC1, coeff)` (`cxx:21-30`): a 1x1
    /// coefficient array.
    pub fn from_ppc_coeff(ppc1: PinpointConstraint, coeff: GpXyz) -> Self {
        Self { ppc: vec![ppc1], coeff: vec![coeff], row_len: 1, col_len: 1 }
    }

    /// `Plate_LinearScalarConstraint(thePPC, Array1<gp_XYZ> theCoeff)`
    /// (`cxx:32-47`): one row, `theCoeff.Length()` columns.
    /// `Standard_DimensionMismatch` becomes `Err`.
    pub fn from_coeff_row(
        the_ppc: Vec<PinpointConstraint>,
        the_coeff: &[GpXyz],
    ) -> Result<Self, String> {
        if the_coeff.len() != the_ppc.len() {
            return Err("Plate_LinearScalarConstraint: dimension mismatch".to_string());
        }
        Ok(Self {
            ppc: the_ppc,
            coeff: the_coeff.to_vec(),
            row_len: 1,
            col_len: the_coeff.len(),
        })
    }

    /// `Plate_LinearScalarConstraint(thePPC, Array2<gp_XYZ> theCoeff)`
    /// (`cxx:49-63`): `RowLength` must equal `thePPC.Length()`.
    pub fn from_coeff_matrix(
        the_ppc: Vec<PinpointConstraint>,
        the_coeff: &[GpXyz],
        row_len: usize,
        col_len: usize,
    ) -> Result<Self, String> {
        if row_len != the_ppc.len() {
            return Err("Plate_LinearScalarConstraint: dimension mismatch".to_string());
        }
        if the_coeff.len() != row_len * col_len {
            return Err("Plate_LinearScalarConstraint: coefficient array size".to_string());
        }
        Ok(Self { ppc: the_ppc, coeff: the_coeff.to_vec(), row_len, col_len })
    }

    /// `Plate_LinearScalarConstraint(ColLen, RowLen)` (`cxx:65-70`):
    /// `RowLen` point constraints and a zero-filled `ColLen` x `RowLen`
    /// coefficient array.
    pub fn with_size(col_len: usize, row_len: usize) -> Self {
        Self {
            ppc: vec![PinpointConstraint::default(); row_len],
            coeff: vec![GpXyz::zero(); row_len * col_len],
            row_len,
            col_len,
        }
    }

    /// `GetPPC()` (`hxx:45`).
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
    pub fn coeff(&self, row: usize, col: usize) -> GpXyz {
        self.coeff[(row - 1) * self.row_len + (col - 1)]
    }

    /// `SetPPC(Index, Value)` (`cxx:72-75`), 1-based.
    pub fn set_ppc(&mut self, index: usize, value: PinpointConstraint) {
        self.ppc[index - 1] = value;
    }

    /// `SetCoeff(Row, Col, Value)` (`cxx:77-80`), 1-based.
    pub fn set_coeff(&mut self, row: usize, col: usize, value: GpXyz) {
        self.coeff[(row - 1) * self.row_len + (col - 1)] = value;
    }
}
