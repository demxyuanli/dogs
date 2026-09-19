//! Intersection point between two 2d elements.
//! Source: `IntAna2d_IntPoint.hxx` + `IntAna2d_IntPoint.cxx` + `IntAna2d_IntPoint.lxx`.
use crate::gp::GpPnt2d;

/// `RealLast()` (`Standard_Real.hxx:179-182`).
const REAL_LAST: f64 = f64::MAX;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IntAna2dIntPoint {
    myu1: f64,
    myu2: f64,
    myp: GpPnt2d,
    myimplicit: bool,
}

impl IntAna2dIntPoint {
    /// Ctor `(X, Y, U1, U2)` (`cxx:19-28`).
    pub fn new4(x: f64, y: f64, u1: f64, u2: f64) -> Self {
        Self {
            myu1: u1,
            myu2: u2,
            myp: GpPnt2d::new(x, y),
            myimplicit: false,
        }
    }

    /// Ctor `(X, Y, U1)` — point on an implicit curve (`cxx:30-36`).
    pub fn new3(x: f64, y: f64, u1: f64) -> Self {
        Self {
            myu1: u1,
            myu2: REAL_LAST,
            myp: GpPnt2d::new(x, y),
            myimplicit: true,
        }
    }

    /// Empty ctor (`cxx:38-45`).
    pub fn new() -> Self {
        Self {
            myu1: REAL_LAST,
            myu2: REAL_LAST,
            myp: GpPnt2d::new(REAL_LAST, REAL_LAST),
            myimplicit: false,
        }
    }

    /// `SetValue(X, Y, U1, U2)` (`cxx:47-55`).
    pub fn set_value4(&mut self, x: f64, y: f64, u1: f64, u2: f64) {
        self.myimplicit = false;
        self.myp.set_coord(x, y);
        self.myu1 = u1;
        self.myu2 = u2;
    }

    /// `SetValue(X, Y, U1)` — implicit second curve (`cxx:57-64`).
    pub fn set_value3(&mut self, x: f64, y: f64, u1: f64) {
        self.myimplicit = true;
        self.myp.set_coord(x, y);
        self.myu1 = u1;
        self.myu2 = REAL_LAST;
    }

    /// `Value` (`lxx:18-21`).
    pub fn value(&self) -> &GpPnt2d {
        &self.myp
    }

    /// `SecondIsImplicit` (`lxx:38-41`).
    pub fn second_is_implicit(&self) -> bool {
        self.myimplicit
    }

    /// `ParamOnFirst` (`lxx:23-26`).
    pub fn param_on_first(&self) -> f64 {
        self.myu1
    }

    /// `ParamOnSecond` (`lxx:28-36`); raises `Standard_DomainError` when the
    /// second curve is implicit.
    pub fn param_on_second(&self) -> f64 {
        assert!(
            !self.myimplicit,
            "Standard_DomainError in IntAna2d_IntPoint::ParamOnSecond"
        );
        self.myu2
    }
}

impl Default for IntAna2dIntPoint {
    fn default() -> Self {
        Self::new()
    }
}
