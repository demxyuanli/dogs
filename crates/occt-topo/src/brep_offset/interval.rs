//! `BRepOffset_Interval` (`BRepOffset_Interval.hxx`, `.cxx`, `.lxx`): one
//! parameter interval of an edge together with the connectivity type its two
//! faces have along that interval.

use crate::chfi_ds::TypeOfConcavity;

/// `BRepOffset_Interval` (`hxx:27-54`).
///
/// OCCT's default constructor is `= default` (`cxx:22`), which leaves `f`, `l`
/// and `type` uninitialised; every call site writes all three before the value
/// is read (`BRepOffset_Analyse.cxx:89-91`, `:143`, `:344-348`). The port
/// starts them at zero / `ChFiDS_Other` instead of leaving them undefined.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Interval {
    /// `f`
    f: f64,
    /// `l`
    l: f64,
    /// `type`
    kind: TypeOfConcavity,
}

impl Default for Interval {
    fn default() -> Self {
        Self {
            f: 0.0,
            l: 0.0,
            kind: TypeOfConcavity::Other,
        }
    }
}

impl Interval {
    /// `BRepOffset_Interval(U1, U2, Type)` (`cxx:26-33`).
    pub fn new(u1: f64, u2: f64, kind: TypeOfConcavity) -> Self {
        Self {
            f: u1,
            l: u2,
            kind,
        }
    }

    /// `First(U)` (`lxx:20-23`).
    pub fn set_first(&mut self, u: f64) {
        self.f = u;
    }

    /// `Last(U)` (`lxx:28-31`).
    pub fn set_last(&mut self, u: f64) {
        self.l = u;
    }

    /// `Type(T)` (`lxx:36-39`).
    pub fn set_type(&mut self, t: TypeOfConcavity) {
        self.kind = t;
    }

    /// `First()` (`lxx:44-47`).
    pub fn first(&self) -> f64 {
        self.f
    }

    /// `Last()` (`lxx:52-55`).
    pub fn last(&self) -> f64 {
        self.l
    }

    /// `Type()` (`lxx:60-63`).
    pub fn kind(&self) -> TypeOfConcavity {
        self.kind
    }
}
