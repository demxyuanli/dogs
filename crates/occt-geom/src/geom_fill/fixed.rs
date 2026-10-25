//! `GeomFill_Fixed` — a constant trihedron.
//!
//! Source: `TKGeomAlgo/GeomFill/GeomFill_Fixed.hxx` and
//! `GeomFill_Fixed.cxx:27-118` (OCCT 8.0.0).
//!
//! Translation notes: `Standard_ConstructionError` from the constructor
//! (`cxx:29-32`) becomes `Err`; `Intervals` fills a one-element array with
//! `-Precision::Infinite()` / `Precision::Infinite()` (`cxx:102-106`), which in
//! the port is the two-element vector `[-Infinite, Infinite]`.

use std::sync::Arc;

use occt_core::gp::GpVec;
use occt_core::precision::INFINITE;

use crate::curve::Curve;
use crate::geom_fill::trihedron_law::{TrihedronD1, TrihedronD2, TrihedronLaw, TrihedronLawCurve};

/// `GeomFill_Fixed` (`GeomFill_Fixed.hxx:33-72`).
pub struct Fixed {
    /// `T` (`hxx:68`).
    t: GpVec,
    /// `N` (`hxx:69`).
    n: GpVec,
    /// `B` (`hxx:70`).
    b: GpVec,
    /// `myCurve` / `myTrimmed` (`GeomFill_TrihedronLaw.hxx:111-112`).
    data: TrihedronLawCurve,
}

impl Fixed {
    /// `GeomFill_Fixed::GeomFill_Fixed` (`cxx:27-38`).
    pub fn new(tangent: &GpVec, normal: &GpVec) -> Result<Self, String> {
        if tangent.is_parallel_ang(normal, 0.01) {
            return Err("Standard_ConstructionError: GeomFill_Fixed : Two parallel vectors !".into());
        }
        let mut t = *tangent;
        t.normalize();
        let mut n = *normal;
        n.normalize();
        let mut b = t.crossed(&n);
        b.normalize();
        Ok(Self {
            t,
            n,
            b,
            data: TrihedronLawCurve::new(),
        })
    }

    /// `T` (`hxx:68`).
    pub fn tangent(&self) -> &GpVec {
        &self.t
    }

    /// `N` (`hxx:69`).
    pub fn normal(&self) -> &GpVec {
        &self.n
    }

    /// `B` (`hxx:70`).
    pub fn bi_normal(&self) -> &GpVec {
        &self.b
    }
}

impl TrihedronLaw for Fixed {
    /// `GeomFill_Fixed` does not override `SetCurve`; the base stores the
    /// handles and returns true (`GeomFill_TrihedronLaw.cxx:25-30`).
    fn set_curve(&self, curve: Arc<dyn Curve>) -> bool {
        self.data.set_curve(curve)
    }

    /// `GeomFill_Fixed::Copy` (`cxx:40-46`): rebuild from `(T, N)` and reattach
    /// the spine.
    fn copy(&self) -> Arc<dyn TrihedronLaw> {
        let copy = Fixed::new(&self.t, &self.n)
            .expect("GeomFill_Fixed::Copy: the stored frame is not parallel");
        if let Some(curve) = self.data.curve() {
            copy.set_curve(curve);
        }
        Arc::new(copy)
    }

    /// `GeomFill_Fixed::D0` (`cxx:48-54`).
    fn d0(&self, _param: f64) -> Result<(GpVec, GpVec, GpVec), String> {
        Ok((self.t, self.n, self.b))
    }

    /// `GeomFill_Fixed::D1` (`cxx:56-72`).
    fn d1(&self, _param: f64) -> Result<TrihedronD1, String> {
        let zero = GpVec::zero();
        Ok(TrihedronD1 {
            tangent: self.t,
            dtangent: zero,
            normal: self.n,
            dnormal: zero,
            binormal: self.b,
            dbinormal: zero,
        })
    }

    /// `GeomFill_Fixed::D2` (`cxx:74-96`).
    fn d2(&self, _param: f64) -> Result<TrihedronD2, String> {
        let zero = GpVec::zero();
        Ok(TrihedronD2 {
            tangent: self.t,
            dtangent: zero,
            d2tangent: zero,
            normal: self.n,
            dnormal: zero,
            d2normal: zero,
            binormal: self.b,
            dbinormal: zero,
            d2binormal: zero,
        })
    }

    /// `GeomFill_Fixed::NbIntervals` (`cxx:98-101`).
    fn nb_intervals(&self, _s: u8) -> i32 {
        1
    }

    /// `GeomFill_Fixed::Intervals` (`cxx:103-107`).
    fn intervals(&self, _s: u8) -> Vec<f64> {
        vec![-INFINITE, INFINITE]
    }

    /// `GeomFill_TrihedronLaw::SetInterval` (`cxx:57-60`).
    fn set_interval(&self, first: f64, last: f64) -> Result<(), String> {
        self.data.set_interval(first, last)
    }

    /// `GeomFill_TrihedronLaw::GetInterval` (`cxx:62-67`).
    fn get_interval(&self) -> (f64, f64) {
        self.data.interval()
    }

    /// `GeomFill_Fixed::GetAverageLaw` (`cxx:109-114`).
    fn get_average_law(&self) -> (GpVec, GpVec, GpVec) {
        (self.t, self.n, self.b)
    }

    /// `GeomFill_Fixed::IsConstant` (`cxx:116-119`).
    fn is_constant(&self) -> bool {
        true
    }
}
