//! `math_FunctionSample`.
//! Source: `math_FunctionSample.hxx:29-52`, `math_FunctionSample.cxx:25-46`.

/// `math_FunctionSample` (`math_FunctionSample.hxx:29-52`): a default sample
/// with a constant parameter step between `A` and `B` over `N` points.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FunctionSample {
    a: f64,
    b: f64,
    n: i32,
}

impl FunctionSample {
    /// `math_FunctionSample(A, B, N)` (`math_FunctionSample.cxx:25-31`).
    pub fn new(a: f64, b: f64, n: i32) -> Self {
        Self { a, b, n }
    }

    /// `Bounds(A, B)` (`math_FunctionSample.cxx:33-37`).
    pub fn bounds(&self) -> (f64, f64) {
        (self.a, self.b)
    }

    /// `NbPoints()` (`math_FunctionSample.cxx:39-42`).
    pub fn nb_points(&self) -> i32 {
        self.n
    }

    /// `GetParameter(Index)` (`math_FunctionSample.cxx:44-47`), 1-based.
    /// Raises `Standard_OutOfRange` when `Index <= 0 || Index > N`.
    pub fn get_parameter(&self, index: i32) -> f64 {
        assert!(
            index > 0 && index <= self.n,
            "Standard_OutOfRange in math_FunctionSample::GetParameter"
        );
        ((self.n - index) as f64 * self.a + (index - 1) as f64 * self.b) / (self.n - 1) as f64
    }
}
