//! `GeomFill_Boundary` and its unattached subclasses — the contour pieces a
//! `GeomFill_ConstrainedFilling` fits to.
//!
//! Source: `GeomFill_Boundary.hxx/.cxx:22-90`
//! (`GeomFill_Boundary::HasNormals` `cxx:31-34`, `Norm` `cxx:36-46`,
//! `D1Norm` `cxx:48-53`, `Points` `cxx:55-60`), `GeomFill_SimpleBound.cxx:24-100`
//! and `GeomFill_DegeneratedBound.cxx:24-84` (OCCT 8.0.0).
//!
//! Translation notes:
//! - `GeomFill_Boundary` is an abstract `Standard_Transient` root with virtual
//!   `Value`/`D1`/`HasNormals`/`Norm`/`D1Norm`/`Reparametrize`/`Bounds`/
//!   `IsDegenerated`; the port keeps it as the [`Boundary`] trait. The virtual
//!   defaults (`HasNormals` false, `Norm`/`D1Norm` raising
//!   `Standard_Failure`, `Points` from `Bounds`/`Value`) live in the trait.
//! - OCCT shares one object behind several handles and calls the non-const
//!   `Reparametrize` through any of them (`GeomFill_ConstrainedFilling` keeps a
//!   `bound[]` array and a `GeomFill_CoonsAlgPatch` over the same objects). The
//!   port therefore stores the mutated state (`myPar`, the degenerate bounds)
//!   behind a `Mutex` and takes `&self` in `reparametrize`, which reproduces
//!   the handle semantics.
//! - `myT3d`/`myTang` are plain `f64` with `&mut` setters: OCCT's
//!   `Tol3d(Tol)`/`Tolang(Tol)` are non-const, and nothing in the filling
//!   algorithm calls them.
//! - `occ::handle<Adaptor3d_Curve>` is `Arc<dyn Curve>`; `Law_Function` is
//!   [`LawBSpFunc`], null until `Reparametrize` runs.

use std::sync::{Arc, Mutex};

use occt_core::gp::{GpPnt, GpVec};

use crate::curve::Curve;
use crate::law::bspline::LawBSpFunc;
use crate::law::function::LawFunction;
use crate::law::services::law_reparametrize;

/// `GeomFill_Boundary` (`GeomFill_Boundary.hxx:37-81`).
pub trait Boundary: Send + Sync {
    /// `GeomFill_Boundary::Value` (`hxx:40`).
    fn value(&self, u: f64) -> GpPnt;

    /// `GeomFill_Boundary::D1` (`hxx:42`).
    fn d1(&self, u: f64) -> (GpPnt, GpVec);

    /// `GeomFill_Boundary::HasNormals` (`cxx:31-34`): false on the root class.
    fn has_normals(&self) -> bool {
        false
    }

    /// `GeomFill_Boundary::Norm` (`cxx:36-46`): `Standard_Failure` on the root.
    fn norm(&self, _u: f64) -> Result<GpVec, String> {
        Err("Standard_Failure: GeomFill_Boundary::Norm : Undefined normals".into())
    }

    /// `GeomFill_Boundary::D1Norm` (`cxx:48-53`): `Standard_Failure` on the root.
    fn d1_norm(&self, _u: f64) -> Result<(GpVec, GpVec), String> {
        Err("Standard_Failure: GeomFill_Boundary::Norm : Undefined normals".into())
    }

    /// `GeomFill_Boundary::Reparametrize` (`hxx:50-57`).
    #[allow(clippy::too_many_arguments)]
    fn reparametrize(
        &self,
        first: f64,
        last: f64,
        has_df: bool,
        has_dl: bool,
        df: f64,
        dl: f64,
        rev: bool,
    ) -> Result<(), String>;

    /// `GeomFill_Boundary::Bounds` (`hxx:61`).
    fn bounds(&self) -> (f64, f64);

    /// `GeomFill_Boundary::IsDegenerated` (`hxx:63`).
    fn is_degenerated(&self) -> bool;

    /// `GeomFill_Boundary::Tol3d` (`cxx:61-65`).
    fn tol3d(&self) -> f64;

    /// `GeomFill_Boundary::Tol3d(Tol)` (`cxx:67-71`).
    fn set_tol3d(&mut self, tol: f64);

    /// `GeomFill_Boundary::Tolang` (`cxx:73-77`).
    fn tolang(&self) -> f64;

    /// `GeomFill_Boundary::Tolang(Tol)` (`cxx:79-83`).
    fn set_tolang(&mut self, tol: f64);

    /// `GeomFill_Boundary::Points` (`cxx:55-60`): `Value` at both bounds.
    fn points(&self) -> (GpPnt, GpPnt) {
        let (first, last) = self.bounds();
        (self.value(first), self.value(last))
    }
}

/// `GeomFill_SimpleBound` (`GeomFill_SimpleBound.hxx:30-78`).
///
/// A boundary given by an unattached 3D curve, optionally reparametrized by a
/// law.
pub struct SimpleBound {
    /// `myC3d` (`hxx:73`).
    curve: Arc<dyn Curve>,
    /// `myPar` (`hxx:74`), set by `Reparametrize`.
    par: Mutex<Option<LawBSpFunc>>,
    /// `GeomFill_Boundary::myT3d`.
    tol3d: f64,
    /// `GeomFill_Boundary::myTang`.
    tolang: f64,
}

impl SimpleBound {
    /// `GeomFill_SimpleBound::GeomFill_SimpleBound` (`cxx:30-36`).
    pub fn new(curve: Arc<dyn Curve>, tol3d: f64, tolang: f64) -> Self {
        Self {
            curve,
            par: Mutex::new(None),
            tol3d,
            tolang,
        }
    }

    /// The boundary curve (`hxx:66`, `myC3d`).
    pub fn curve(&self) -> &Arc<dyn Curve> {
        &self.curve
    }

    /// A snapshot of the reparametrization law (`myPar`), if it has run.
    pub fn law(&self) -> Option<LawBSpFunc> {
        self.par.lock().unwrap_or_else(|e| e.into_inner()).clone()
    }

    /// The law's value at `u`, or `u` itself when `myPar` is null (`cxx:41-44`).
    fn mapped_value(&self, u: f64) -> f64 {
        let par = self.par.lock().unwrap_or_else(|e| e.into_inner());
        match par.as_ref() {
            Some(law) => law.value(u),
            None => u,
        }
    }

    /// `D1`'s `(x, dx)` pair (`cxx:58-62`).
    fn mapped_d1(&self, u: f64) -> (f64, f64) {
        let par = self.par.lock().unwrap_or_else(|e| e.into_inner());
        match par.as_ref() {
            Some(law) => law.d1(u),
            None => (u, 1.0),
        }
    }
}

impl Boundary for SimpleBound {
    /// `GeomFill_SimpleBound::Value` (`cxx:40-48`).
    fn value(&self, u: f64) -> GpPnt {
        self.curve.d0(self.mapped_value(u))
    }

    /// `GeomFill_SimpleBound::D1` (`cxx:52-62`).
    fn d1(&self, u: f64) -> (GpPnt, GpVec) {
        let (x, dx) = self.mapped_d1(u);
        let (p, v) = self.curve.d1(x);
        (p, v.multiplied_scalar(dx))
    }

    /// `GeomFill_SimpleBound::Reparametrize` (`cxx:66-77`).
    fn reparametrize(
        &self,
        first: f64,
        last: f64,
        has_df: bool,
        has_dl: bool,
        df: f64,
        dl: f64,
        rev: bool,
    ) -> Result<(), String> {
        let curve = law_reparametrize(
            self.curve.as_ref(),
            first,
            last,
            has_df,
            has_dl,
            df,
            dl,
            rev,
            30,
        )?;
        let mut par = LawBSpFunc::new();
        par.set_curve(curve);
        *self.par.lock().unwrap_or_else(|e| e.into_inner()) = Some(par);
        Ok(())
    }

    /// `GeomFill_SimpleBound::Bounds` (`cxx:81-93`).
    fn bounds(&self) -> (f64, f64) {
        let par = self.par.lock().unwrap_or_else(|e| e.into_inner());
        match par.as_ref() {
            Some(law) => LawBSpFunc::bounds(law),
            None => (self.curve.first_parameter(), self.curve.last_parameter()),
        }
    }

    /// `GeomFill_SimpleBound::IsDegenerated` (`cxx:97-100`).
    fn is_degenerated(&self) -> bool {
        false
    }

    fn tol3d(&self) -> f64 {
        self.tol3d
    }

    fn set_tol3d(&mut self, tol: f64) {
        self.tol3d = tol;
    }

    fn tolang(&self) -> f64 {
        self.tolang
    }

    fn set_tolang(&mut self, tol: f64) {
        self.tolang = tol;
    }
}

/// `GeomFill_DegeneratedBound` (`GeomFill_DegeneratedBound.hxx:30-78`).
///
/// A boundary collapsed to a single point.
pub struct DegeneratedBound {
    /// `myPoint` (`hxx:73`).
    point: GpPnt,
    /// `myFirst`/`myLast` (`hxx:74-75`), rewritten by `Reparametrize`.
    range: Mutex<(f64, f64)>,
    tol3d: f64,
    tolang: f64,
}

impl DegeneratedBound {
    /// `GeomFill_DegeneratedBound::GeomFill_DegeneratedBound` (`cxx:26-37`).
    pub fn new(point: GpPnt, first: f64, last: f64, tol3d: f64, tolang: f64) -> Self {
        Self {
            point,
            range: Mutex::new((first, last)),
            tol3d,
            tolang,
        }
    }

    /// The collapsed point (`hxx:70`).
    pub fn point(&self) -> &GpPnt {
        &self.point
    }
}

impl Boundary for DegeneratedBound {
    /// `GeomFill_DegeneratedBound::Value` (`cxx:41-44`).
    fn value(&self, _u: f64) -> GpPnt {
        self.point
    }

    /// `GeomFill_DegeneratedBound::D1` (`cxx:48-53`).
    fn d1(&self, _u: f64) -> (GpPnt, GpVec) {
        (self.point, GpVec::zero())
    }

    /// `GeomFill_DegeneratedBound::Reparametrize` (`cxx:57-67`): only the
    /// bounds move.
    fn reparametrize(
        &self,
        first: f64,
        last: f64,
        _has_df: bool,
        _has_dl: bool,
        _df: f64,
        _dl: f64,
        _rev: bool,
    ) -> Result<(), String> {
        *self.range.lock().unwrap_or_else(|e| e.into_inner()) = (first, last);
        Ok(())
    }

    /// `GeomFill_DegeneratedBound::Bounds` (`cxx:72-77`).
    fn bounds(&self) -> (f64, f64) {
        *self.range.lock().unwrap_or_else(|e| e.into_inner())
    }

    /// `GeomFill_DegeneratedBound::IsDegenerated` (`cxx:80-83`).
    fn is_degenerated(&self) -> bool {
        true
    }

    fn tol3d(&self) -> f64 {
        self.tol3d
    }

    fn set_tol3d(&mut self, tol: f64) {
        self.tol3d = tol;
    }

    fn tolang(&self) -> f64 {
        self.tolang
    }

    fn set_tolang(&mut self, tol: f64) {
        self.tolang = tol;
    }
}

/// Shared reparametrization step of [`SimpleBound`] and
/// [`crate::geom_fill::bound_with_surf::BoundWithSurf`] (`Law::Reparametrize`
/// with 30 samples, `Law_BSpFunc::SetCurve`).
pub(crate) fn make_par_law(
    curve: &dyn Curve,
    first: f64,
    last: f64,
    has_df: bool,
    has_dl: bool,
    df: f64,
    dl: f64,
    rev: bool,
) -> Result<LawBSpFunc, String> {
    let bs = law_reparametrize(curve, first, last, has_df, has_dl, df, dl, rev, 30)?;
    let mut par = LawBSpFunc::new();
    par.set_curve(bs);
    Ok(par)
}
