//! `NCollection_Array1<double>` with OCCT bounds.
//!
//! Source: `FoundationClasses/TKernel/NCollection/NCollection_Array1.hxx`.
//!
//! `FEmTool_Assembly` (constraint coefficient segments, `cxx:532-547`) and
//! `AppDef_Variational` (`myTabPoints`, `myTypConstraints`, `myTabConstraints`,
//! `myTtheta`, `myTfthet`) address 1-based arrays whose lower bound is not
//! always 1, so the port keeps the bound instead of shifting the indices.

/// `NCollection_Array1<double>` (`NCollection_Array1.hxx`).
#[derive(Clone, Debug)]
pub struct Array1d {
    lower: i32,
    data: Vec<f64>,
}

impl Array1d {
    /// `NCollection_Array1(Lower, Upper)` — OCCT leaves the values undefined;
    /// the port zero-fills them.
    pub fn zeros(lower: i32, upper: i32) -> Self {
        Self {
            lower,
            data: vec![0.0; (upper - lower + 1).max(0) as usize],
        }
    }

    /// Array holding `values` starting at `lower`.
    pub fn from_slice(lower: i32, values: &[f64]) -> Self {
        Self {
            lower,
            data: values.to_vec(),
        }
    }

    /// `Lower()`.
    #[inline]
    pub fn lower(&self) -> i32 {
        self.lower
    }

    /// `Upper()`.
    #[inline]
    pub fn upper(&self) -> i32 {
        self.lower + self.data.len() as i32 - 1
    }

    /// `Length()`.
    #[inline]
    pub fn len(&self) -> usize {
        self.data.len()
    }

    /// `IsEmpty()`.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// `Value(i)`.
    #[inline]
    pub fn value(&self, i: i32) -> f64 {
        self.data[(i - self.lower) as usize]
    }

    /// `ChangeValue(i) = v` / `SetValue(i, v)`.
    #[inline]
    pub fn set_value(&mut self, i: i32, v: f64) {
        self.data[(i - self.lower) as usize] = v;
    }

    /// `ChangeValue(i) += v`.
    #[inline]
    pub fn add_at(&mut self, i: i32, v: f64) {
        self.data[(i - self.lower) as usize] += v;
    }

    /// `ChangeValue(i) *= f`.
    #[inline]
    pub fn scale_at(&mut self, i: i32, f: f64) {
        self.data[(i - self.lower) as usize] *= f;
    }

    /// `Init(v)`.
    pub fn init(&mut self, v: f64) {
        self.data.fill(v);
    }

    /// Mutable view of the whole array (OCCT `ChangeArray1()`).
    #[inline]
    pub fn as_mut_slice(&mut self) -> &mut [f64] {
        &mut self.data
    }

    /// Read-only view of the whole array (OCCT `Array1()`).
    #[inline]
    pub fn as_slice(&self) -> &[f64] {
        &self.data
    }
}
