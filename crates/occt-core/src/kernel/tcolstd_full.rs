//! Collection helpers mirroring Open CASCADE's `TColStd` package.
//!
//! Provides small `Array1`-style wrappers over `Vec` with a configurable
//! lower bound (OCCT arrays are indexable from an arbitrary integer, most
//! commonly 1), a handle-style wrapper, and a handful of numeric utilities.

/// A 1-D array of doubles with a configurable lower index bound.
pub struct Array1OfReal {
    /// Backing storage, indexed `[0..len)` corresponding to `[lower..upper]`.
    pub data: Vec<f64>,
    /// Index of the first element.
    pub lower: i32,
}

impl Array1OfReal {
    /// Creates an array spanning `[lower, upper]`, zero-filled.
    ///
    /// If `upper < lower` the array is empty (length 0).
    pub fn new(lower: i32, upper: i32) -> Self {
        let len = if upper < lower {
            0
        } else {
            (upper - lower + 1) as usize
        };
        Self {
            data: vec![0.0; len],
            lower,
        }
    }

    /// Creates a 1-based array from a slice.
    pub fn from_slice(values: &[f64]) -> Self {
        Self {
            data: values.to_vec(),
            lower: 1,
        }
    }

    /// Returns the element at index `i`. Panics if out of range.
    pub fn value(&self, i: i32) -> f64 {
        self.data[(i - self.lower) as usize]
    }

    /// Writes `v` at index `i`. Panics if out of range.
    pub fn set_value(&mut self, i: i32, v: f64) {
        self.data[(i - self.lower) as usize] = v;
    }

    /// Lower index bound.
    pub fn lower(&self) -> i32 {
        self.lower
    }

    /// Upper index bound (`lower - 1` when empty).
    pub fn upper(&self) -> i32 {
        if self.data.is_empty() {
            self.lower - 1
        } else {
            self.lower + self.data.len() as i32 - 1
        }
    }

    /// Number of elements.
    pub fn len(&self) -> i32 {
        self.data.len() as i32
    }

    /// Returns `true` if the array holds no elements.
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// Fills every element with `v`.
    pub fn init(&mut self, v: f64) {
        for x in &mut self.data {
            *x = v;
        }
    }

    /// Borrows the backing storage as a plain slice.
    pub fn as_slice(&self) -> &[f64] {
        &self.data
    }
}

/// A 1-D array of 32-bit integers with a configurable lower index bound.
pub struct Array1OfInteger {
    /// Backing storage, indexed `[0..len)` corresponding to `[lower..upper]`.
    pub data: Vec<i32>,
    /// Index of the first element.
    pub lower: i32,
}

impl Array1OfInteger {
    /// Creates an array spanning `[lower, upper]`, zero-filled.
    pub fn new(lower: i32, upper: i32) -> Self {
        let len = if upper < lower {
            0
        } else {
            (upper - lower + 1) as usize
        };
        Self {
            data: vec![0; len],
            lower,
        }
    }

    /// Creates a 1-based array from a slice.
    pub fn from_slice(values: &[i32]) -> Self {
        Self {
            data: values.to_vec(),
            lower: 1,
        }
    }

    /// Returns the element at index `i`. Panics if out of range.
    pub fn value(&self, i: i32) -> i32 {
        self.data[(i - self.lower) as usize]
    }

    /// Writes `v` at index `i`. Panics if out of range.
    pub fn set_value(&mut self, i: i32, v: i32) {
        self.data[(i - self.lower) as usize] = v;
    }

    /// Lower index bound.
    pub fn lower(&self) -> i32 {
        self.lower
    }

    /// Upper index bound (`lower - 1` when empty).
    pub fn upper(&self) -> i32 {
        if self.data.is_empty() {
            self.lower - 1
        } else {
            self.lower + self.data.len() as i32 - 1
        }
    }

    /// Number of elements.
    pub fn len(&self) -> i32 {
        self.data.len() as i32
    }

    /// Returns `true` if the array holds no elements.
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// Fills every element with `v`.
    pub fn init(&mut self, v: i32) {
        for x in &mut self.data {
            *x = v;
        }
    }

    /// Borrows the backing storage as a plain slice.
    pub fn as_slice(&self) -> &[i32] {
        &self.data
    }
}

/// Handle-style wrapper around [`Array1OfReal`] (analogous to `Handle<TColStd_Array1OfReal>`).
pub struct HArray1OfReal {
    /// The wrapped array.
    pub array: Array1OfReal,
}

impl HArray1OfReal {
    /// Creates a handle wrapping an array spanning `[lower, upper]`.
    pub fn new(lower: i32, upper: i32) -> Self {
        Self {
            array: Array1OfReal::new(lower, upper),
        }
    }

    /// Creates a handle wrapping a 1-based array from a slice.
    pub fn from_slice(values: &[f64]) -> Self {
        Self {
            array: Array1OfReal::from_slice(values),
        }
    }

    /// Returns the element at index `i`. Panics if out of range.
    pub fn value(&self, i: i32) -> f64 {
        self.array.value(i)
    }

    /// Writes `v` at index `i`. Panics if out of range.
    pub fn change_value(&mut self, i: i32, v: f64) {
        self.array.set_value(i, v);
    }

    /// Number of elements.
    pub fn len(&self) -> i32 {
        self.array.len()
    }

    /// Returns `true` if the array holds no elements.
    pub fn is_empty(&self) -> bool {
        self.array.is_empty()
    }
}

/// Smallest element in `a`. Returns `INFINITY` for an empty slice.
pub fn min_of_array(a: &[f64]) -> f64 {
    a.iter().copied().fold(f64::INFINITY, f64::min)
}

/// Largest element in `a`. Returns `-INFINITY` for an empty slice.
pub fn max_of_array(a: &[f64]) -> f64 {
    a.iter().copied().fold(f64::NEG_INFINITY, f64::max)
}

/// Sum of all elements in `a`.
pub fn sum_array(a: &[f64]) -> f64 {
    a.iter().sum()
}

/// Arithmetic mean of `a`; `0.0` for an empty slice.
pub fn mean_array(a: &[f64]) -> f64 {
    if a.is_empty() {
        0.0
    } else {
        sum_array(a) / a.len() as f64
    }
}

/// Sorts `a` in ascending order (`NaN` sorted to the end).
pub fn sort_f64(a: &mut [f64]) {
    a.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Greater));
}

/// Sorts `a` ascending and removes duplicate values.
pub fn unique_f64(a: &mut Vec<f64>) {
    sort_f64(a);
    a.dedup();
}

/// Linear interpolation/extrapolation through `(x0, y0)` and `(x1, y1)` at `x`.
///
/// If `x0 == x1` the segment is degenerate and `y0` is returned.
pub fn interpolate_linear(x0: f64, y0: f64, x1: f64, y1: f64, x: f64) -> f64 {
    if x1 == x0 {
        return y0;
    }
    y0 + (y1 - y0) * (x - x0) / (x1 - x0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn array1_of_real_access_and_bounds() {
        let mut a = Array1OfReal::new(1, 3);
        a.set_value(1, 1.0);
        a.set_value(2, 2.0);
        a.set_value(3, 3.0);
        assert_eq!(a.value(1), 1.0);
        assert_eq!(a.value(3), 3.0);
        assert_eq!(a.lower(), 1);
        assert_eq!(a.upper(), 3);
        assert_eq!(a.len(), 3);
        assert!(!a.is_empty());

        a.init(7.0);
        assert_eq!(a.as_slice(), &[7.0, 7.0, 7.0]);
    }

    #[test]
    #[should_panic]
    fn array1_of_real_out_of_bounds_panics() {
        let a = Array1OfReal::new(1, 3);
        let _ = a.value(4);
    }

    #[test]
    fn array1_offsets_and_empty() {
        let b = Array1OfReal::new(0, 2);
        assert_eq!(b.len(), 3);
        assert_eq!(b.lower(), 0);
        assert_eq!(b.upper(), 2);

        let e = Array1OfReal::new(5, 2);
        assert_eq!(e.len(), 0);
        assert!(e.is_empty());
        assert_eq!(e.upper(), 4);
    }

    #[test]
    fn array1_of_integer_from_slice() {
        let mut a = Array1OfInteger::from_slice(&[10, 20, 30]);
        assert_eq!(a.value(2), 20);
        a.set_value(3, 99);
        assert_eq!(a.value(3), 99);
        assert_eq!(a.lower(), 1);
        assert_eq!(a.upper(), 3);
    }

    #[test]
    fn harray1_of_real_wraps() {
        let mut h = HArray1OfReal::new(1, 2);
        h.change_value(1, 0.5);
        h.change_value(2, 1.5);
        assert_eq!(h.value(1), 0.5);
        assert_eq!(h.value(2), 1.5);
        assert_eq!(h.len(), 2);
    }

    #[test]
    fn min_max_sum_mean() {
        let a = [1.0, 2.0, 3.0, 4.0];
        assert_eq!(min_of_array(&a), 1.0);
        assert_eq!(max_of_array(&a), 4.0);
        assert_eq!(sum_array(&a), 10.0);
        assert_eq!(mean_array(&a), 2.5);
        assert!(mean_array(&[]).is_finite());
    }

    #[test]
    fn sort_and_unique() {
        let mut v = vec![3.0, 1.0, 2.0];
        sort_f64(&mut v);
        assert_eq!(v, vec![1.0, 2.0, 3.0]);

        let mut u = vec![3.0, 1.0, 2.0, 1.0, 3.0];
        unique_f64(&mut u);
        assert_eq!(u, vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn interpolate_linear_works() {
        assert_eq!(interpolate_linear(0.0, 0.0, 10.0, 100.0, 5.0), 50.0);
        assert_eq!(interpolate_linear(0.0, 0.0, 10.0, 100.0, 20.0), 200.0);
        // Degenerate segment returns y0.
        assert_eq!(interpolate_linear(2.0, 5.0, 2.0, 9.0, 2.0), 5.0);
    }
}
