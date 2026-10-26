//! `NCollection_Array2<T>` with OCCT lower/upper bounds.
//!
//! Source: `FoundationClasses/TKernel/NCollection/NCollection_Array2.hxx`.
//!
//! `FEmTool` stores element coefficients with an OCCT row range that does not
//! start at 1 (`AppDef_LinearCriteria.cxx:98-101` builds `(0, 0, 1, NbDim)`,
//! `AppDef_Variational.cxx:330` rebuilds `(0, deg, 1, NbDim)`), and the
//! reference-matrix criteria read `myCoeff->LowerRow()` / `LowerCol()`. The port
//! keeps the bounds instead of shifting the indices.
//!
//! The type is generic because the same container also carries the `int`
//! dependence / assembly tables (`AppDef_SmoothCriterion::AssemblyTable`,
//! `FEmTool_ElementsOfRefMatrix`) and the per-constraint global index lists
//! (`Vec<i32>`).

/// `NCollection_Array2<T>` (`NCollection_Array2.hxx`).
#[derive(Clone, Debug)]
pub struct Array2<T> {
    lower_row: usize,
    upper_row: usize,
    lower_col: usize,
    upper_col: usize,
    data: Vec<T>,
}

impl<T: Clone> Array2<T> {
    /// `NCollection_Array2(LowerRow, UpperRow, LowerCol, UpperCol)` filled with
    /// `init`.
    pub fn new(lower_row: i32, upper_row: i32, lower_col: i32, upper_col: i32, init: T) -> Self {
        let lr = lower_row.max(0) as usize;
        let lc = lower_col.max(0) as usize;
        let ur = (upper_row.max(0) as usize).max(lr);
        let uc = (upper_col.max(0) as usize).max(lc);
        let nr = ur - lr + 1;
        let nc = uc - lc + 1;
        Self {
            lower_row: lr,
            upper_row: ur,
            lower_col: lc,
            upper_col: uc,
            data: vec![init; nr * nc],
        }
    }

    /// `Init(Value)`.
    pub fn init(&mut self, value: T) {
        for v in self.data.iter_mut() {
            *v = value.clone();
        }
    }

    /// `LowerRow()`.
    pub fn lower_row(&self) -> usize {
        self.lower_row
    }

    /// `UpperRow()`.
    pub fn upper_row(&self) -> usize {
        self.upper_row
    }

    /// `LowerCol()`.
    pub fn lower_col(&self) -> usize {
        self.lower_col
    }

    /// `UpperCol()`.
    pub fn upper_col(&self) -> usize {
        self.upper_col
    }

    /// `NbRows()`.
    pub fn nb_rows(&self) -> usize {
        self.upper_row - self.lower_row + 1
    }

    /// `NbColumns()`.
    pub fn nb_columns(&self) -> usize {
        self.upper_col - self.lower_col + 1
    }

    /// `RowLength()` — OCCT returns the **number of columns**
    /// (`NCollection_Array2.hxx:198-199`). The FEmTool criteria rely on this:
    /// `NbDim = myCoeff->RowLength()` where the coefficient rows of `myCoeff`
    /// are the polynomial degree and its columns are the dimensions.
    pub fn row_length(&self) -> usize {
        self.nb_columns()
    }

    /// `ColLength()` — OCCT returns the **number of rows**
    /// (`NCollection_Array2.hxx:201`), i.e. the number of coefficients per
    /// dimension for the `FEmTool` criterion arrays.
    pub fn col_length(&self) -> usize {
        self.nb_rows()
    }

    /// `Value(Row, Col)` / `ChangeValue(Row, Col)` address the element at the
    /// given (row, col) offsets from the lower bounds.
    fn offset(&self, row: usize, col: usize) -> usize {
        (col - self.lower_col) * self.nb_rows() + (row - self.lower_row)
    }

    /// `Value(Row, Col)`.
    pub fn value(&self, row: usize, col: usize) -> &T {
        &self.data[self.offset(row, col)]
    }

    /// `SetValue(Row, Col, Value)`.
    pub fn set_value(&mut self, row: usize, col: usize, value: T) {
        let i = self.offset(row, col);
        self.data[i] = value;
    }

    /// `ChangeValue(Row, Col)`.
    pub fn change_value(&mut self, row: usize, col: usize) -> &mut T {
        let i = self.offset(row, col);
        &mut self.data[i]
    }
}

impl Array2<f64> {
    /// `Value(Row, Col)` for the coefficient arrays.
    pub fn at(&self, row: usize, col: usize) -> f64 {
        *self.value(row, col)
    }

    /// `SetValue(Row, Col, Value)` for the coefficient arrays.
    pub fn set(&mut self, row: usize, col: usize, value: f64) {
        self.set_value(row, col, value);
    }

    /// `ChangeValue(Row, Col) += Value`.
    pub fn add(&mut self, row: usize, col: usize, value: f64) {
        *self.change_value(row, col) += value;
    }

    /// `ChangeValue(Row, Col) *= Value`.
    pub fn scale(&mut self, row: usize, col: usize, factor: f64) {
        *self.change_value(row, col) *= factor;
    }
}
