//! `GeomFill_LocFunction` — the `GeomFill_LocationLaw` view `GeomFill_Sweep`
//! and `Approx_SweepApproximation` evaluate as a single vector function: the
//! frame origin followed by the three frame columns.
//!
//! Source: `TKGeomAlgo/GeomFill/GeomFill_LocFunction.{hxx,cxx:20-150}`
//! (OCCT 8.0.0).
//!
//! Translation notes:
//! - `NCollection_Array1<gp_Vec> V(1,4)` becomes `[GpVec; 4]`; `V(1)` is the
//!   law translation, `V(2..4)` the matrix columns (`cxx:44-46`).
//! - The `First` / `Last` arguments are unused by every `D0`/`D1`/`D2` body in
//!   OCCT (`cxx:33-34`, `cxx:49-50`, `cxx:71-72`) but stay in the signatures
//!   because `AdvApprox_EvaluatorFunction` passes them.
//! - `myLaw->D0/D1/D2` return a `bool` in OCCT; the port's [`LocationLaw`]
//!   returns `Result`, and `Err` is exactly that `false`.
//! - `DN` (`cxx:107-150`) copies the twelve doubles of `V` / `DV` / `D2V`
//!   into the caller's array. OCCT leaves `LocalResult` null for `Order > 2`
//!   and then dereferences it (undefined behaviour); the port returns the
//!   `Ier` code and an all-zero array instead.
//! - `M` / `DM` / `D2M` (`hxx:63-65`) are never read in OCCT; they are kept
//!   and filled here so the struct layout matches the class.

use std::sync::Arc;

use occt_core::gp::{GpMat, GpVec};

use crate::geom_fill::location_law::LocationLaw;

/// `GeomFill_LocFunction` (`GeomFill_LocFunction.hxx:33-67`).
pub struct LocFunction {
    /// `myLaw` (`hxx:61`).
    law: Arc<dyn LocationLaw>,
    /// `V` (`hxx:62`): origin, then the three matrix columns.
    v: [GpVec; 4],
    /// `DV` (`hxx:63`).
    dv: [GpVec; 4],
    /// `D2V` (`hxx:64`).
    d2v: [GpVec; 4],
    /// `M` (`hxx:65`).
    m: GpMat,
    /// `DM` (`hxx:66`).
    dm: GpMat,
    /// `D2M` (`hxx:67`).
    d2m: GpMat,
}

impl LocFunction {
    /// `GeomFill_LocFunction(Law)` (`cxx:20-29`).
    pub fn new(law: Arc<dyn LocationLaw>) -> Self {
        Self {
            law,
            v: [GpVec::zero(); 4],
            dv: [GpVec::zero(); 4],
            d2v: [GpVec::zero(); 4],
            m: GpMat::identity(),
            dm: GpMat::identity(),
            d2m: GpMat::identity(),
        }
    }

    /// `V` (`hxx:62`): `[origin, col1, col2, col3]`.
    pub fn v(&self) -> &[GpVec; 4] {
        &self.v
    }

    /// `DV` (`hxx:63`).
    pub fn dv(&self) -> &[GpVec; 4] {
        &self.dv
    }

    /// `D2V` (`hxx:64`).
    pub fn d2v(&self) -> &[GpVec; 4] {
        &self.d2v
    }

    fn store_columns(m: &GpMat, out: &mut [GpVec; 4]) {
        out[1] = GpVec::from_xyz(&m.column(1));
        out[2] = GpVec::from_xyz(&m.column(2));
        out[3] = GpVec::from_xyz(&m.column(3));
    }

    /// `D0(Param, First, Last)` (`cxx:31-47`).
    pub fn d0(&mut self, param: f64, _first: f64, _last: f64) -> bool {
        let (m, position) = match self.law.d0(param) {
            Ok(r) => r,
            Err(_) => return false,
        };
        self.v[0] = position;
        Self::store_columns(&m, &mut self.v);
        self.m = m;
        true
    }

    /// `D1(Param, First, Last)` (`cxx:49-69`).
    pub fn d1(&mut self, param: f64, _first: f64, _last: f64) -> bool {
        let r = match self.law.d1(param) {
            Ok(r) => r,
            Err(_) => return false,
        };
        self.v[0] = r.v;
        self.dv[0] = r.dv;
        Self::store_columns(&r.m, &mut self.v);
        Self::store_columns(&r.dm, &mut self.dv);
        self.m = r.m;
        self.dm = r.dm;
        true
    }

    /// `D2(Param, First, Last)` (`cxx:71-105`).
    pub fn d2(&mut self, param: f64, _first: f64, _last: f64) -> bool {
        let r = match self.law.d2(param) {
            Ok(r) => r,
            Err(_) => return false,
        };
        self.v[0] = r.v;
        self.dv[0] = r.dv;
        self.d2v[0] = r.d2v;
        Self::store_columns(&r.m, &mut self.v);
        Self::store_columns(&r.dm, &mut self.dv);
        Self::store_columns(&r.d2m, &mut self.d2v);
        self.m = r.m;
        self.dm = r.dm;
        self.d2m = r.d2m;
        true
    }

    /// `DN(Param, First, Last, Order, Result, Ier)` (`cxx:107-150`): copies
    /// the twelve doubles of the requested derivative level.
    pub fn dn(
        &mut self,
        param: f64,
        first: f64,
        last: f64,
        order: i32,
    ) -> ([f64; 12], i32) {
        let mut ier = 0;
        let mut out = [0.0; 12];
        let ok = match order {
            0 => self.d0(param, first, last),
            1 => self.d1(param, first, last),
            2 => self.d2(param, first, last),
            _ => false,
        };
        if !ok {
            ier = order + 1;
        }
        let src: Option<&[GpVec; 4]> = match order {
            0 => Some(&self.v),
            1 => Some(&self.dv),
            2 => Some(&self.d2v),
            _ => None,
        };
        if let Some(src) = src {
            for (k, vec) in src.iter().enumerate() {
                out[3 * k] = vec.coord.x;
                out[3 * k + 1] = vec.coord.y;
                out[3 * k + 2] = vec.coord.z;
            }
        }
        (out, ier)
    }
}
