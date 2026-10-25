//! `GeomLib_LogSample` — a `math_FunctionSample` whose inner parameters grow
//! geometrically.
//! Source: `ModelingData/TKGeomBase/GeomLib/GeomLib_LogSample.cxx:19-48`
//! (OCCT 8.0.0), declaration `GeomLib_LogSample.hxx:26-44`.
//!
//! `GetParameter(1)` and `GetParameter(NbPoints)` are the declared bounds;
//! every other index `i` answers `myF + exp(myexp * i)` with
//! `myF = A - 1` and `myexp = log(B - A) / N` (`cxx:21-27`, `cxx:40-42`).

use occt_core::math_function_sample::SampleSource;

/// `GeomLib_LogSample` (`hxx:26-44`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LogSample {
    /// `myF` (`hxx:40`): `A - 1`.
    f: f64,
    /// `myexp` (`hxx:42`): `log(B - A) / N`.
    exp: f64,
    a: f64,
    b: f64,
    n: i32,
}

impl LogSample {
    /// `GeomLib_LogSample(A, B, N)` (`cxx:21-27`); the base
    /// `math_FunctionSample(A, B, N)` only stores the bounds.
    pub fn new(a: f64, b: f64, n: i32) -> Self {
        Self {
            f: a - 1.0,
            exp: (b - a).ln() / n as f64,
            a,
            b,
            n,
        }
    }
}

impl SampleSource for LogSample {
    /// `NbPoints()` (`math_FunctionSample.cxx:39-42`).
    fn nb_points(&self) -> i32 {
        self.n
    }

    /// `GetParameter(Index)` (`cxx:30-46`): the `Standard_OutOfRange` of
    /// `cxx:44` becomes an assertion, as in
    /// [`FunctionSample::parameter`](occt_core::math_function_sample::FunctionSample::parameter).
    fn parameter(&self, index: i32) -> f64 {
        if index < self.n && index > 1 {
            return self.f + (self.exp * index as f64).exp();
        }
        if index == 1 {
            return self.a;
        }
        if index == self.n {
            return self.b;
        }
        panic!("Standard_OutOfRange in GeomLib_LogSample::GetParameter");
    }

    /// `Bounds(A, B)` (`math_FunctionSample.cxx:33-37`).
    fn bounds(&self) -> (f64, f64) {
        (self.a, self.b)
    }
}
