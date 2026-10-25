//! `GeomFill_BoundWithSurf` — a boundary attached to an existing surface.
//!
//! Source: `GeomFill_BoundWithSurf.cxx:33-193` (OCCT 8.0.0):
//! constructor `cxx:36-40`, `Value` `cxx:44-54`, `D1` `cxx:58-66`,
//! `HasNormals` `cxx:70-73`, `Norm` `cxx:77-100`, `D1Norm` `cxx:104-158`,
//! `Reparametrize` `cxx:162-175`, `Bounds` `cxx:179-187`,
//! `IsDegenerated` `cxx:191-194`.
//!
//! Translation notes:
//! - OCCT stores one `Adaptor3d_CurveOnSurface myConS`, which is the bundle
//!   `(3D curve, pcurve, surface)`; the port takes the three pieces apart:
//!   `curve` (what `myConS.Value`/`D1` evaluate), `pcurve` + `surface` (what
//!   `myConS.GetCurve()`/`GetSurface()` expose to `Norm`/`D1Norm`). Callers
//!   pass the same pcurve/surface pair the adaptor was built from, e.g. the
//!   port's `CurveOnSurface` in `occt-topo::meshing`.
//! - `myPar` is a `Law_BSpFunc` handle, null until `Reparametrize`; it lives
//!   behind a `Mutex` so `Reparametrize` can keep taking `&self` like OCCT's
//!   non-const method on a shared handle.
//! - `Norm`/`D1Norm` reproduce the OCCT algebra term for term, including the
//!   `deno < 1e-16` fallback `(Norm(U + 1e-12) - N) * 1e-12`.

use std::sync::{Arc, Mutex};

use occt_core::gp::{GpPnt, GpVec};
use occt_geom2d::curve::Curve2d;

use crate::curve::Curve;
use crate::geom_fill::boundary::{make_par_law, Boundary};
use crate::law::bspline::LawBSpFunc;
use crate::law::function::LawFunction;
use crate::surface::Surface;

/// `GeomFill_BoundWithSurf` (`GeomFill_BoundWithSurf.hxx:37-103`).
pub struct BoundWithSurf {
    /// `myConS` (`hxx:101`) — its 3D evaluation.
    curve: Arc<dyn Curve>,
    /// `myConS.GetCurve()` (`hxx:101`).
    pcurve: Arc<dyn Curve2d>,
    /// `myConS.GetSurface()` (`hxx:101`).
    surface: Arc<dyn Surface>,
    /// `myPar` (`hxx:102`).
    par: Mutex<Option<LawBSpFunc>>,
    tol3d: f64,
    tolang: f64,
}

impl BoundWithSurf {
    /// `GeomFill_BoundWithSurf(CurveOnSurf, Tol3d, Tolang)` (`cxx:36-40`).
    pub fn new(
        curve: Arc<dyn Curve>,
        pcurve: Arc<dyn Curve2d>,
        surface: Arc<dyn Surface>,
        tol3d: f64,
        tolang: f64,
    ) -> Self {
        Self {
            curve,
            pcurve,
            surface,
            par: Mutex::new(None),
            tol3d,
            tolang,
        }
    }

    /// A snapshot of the reparametrization law (`myPar`), if it has run.
    pub fn law(&self) -> Option<LawBSpFunc> {
        self.par.lock().unwrap_or_else(|e| e.into_inner()).clone()
    }

    /// The attached surface (`myConS.GetSurface()`).
    pub fn surface(&self) -> &Arc<dyn Surface> {
        &self.surface
    }

    /// The pcurve (`myConS.GetCurve()`).
    pub fn pcurve(&self) -> &Arc<dyn Curve2d> {
        &self.pcurve
    }

    /// The 3D evaluation curve (`myConS` as an `Adaptor3d_Curve`).
    pub fn curve(&self) -> &Arc<dyn Curve> {
        &self.curve
    }
}

impl Boundary for BoundWithSurf {
    /// `GeomFill_BoundWithSurf::Value` (`cxx:44-54`).
    fn value(&self, u: f64) -> GpPnt {
        let par = self.par.lock().unwrap_or_else(|e| e.into_inner());
        match par.as_ref() {
            Some(law) => self.curve.d0(law.value(u)),
            None => self.curve.d0(u),
        }
    }

    /// `GeomFill_BoundWithSurf::D1` (`cxx:58-66`).
    fn d1(&self, u: f64) -> (GpPnt, GpVec) {
        let par = self.par.lock().unwrap_or_else(|e| e.into_inner());
        let (x, dx) = match par.as_ref() {
            Some(law) => law.d1(u),
            None => (u, 1.0),
        };
        drop(par);
        let (p, v) = self.curve.d1(x);
        (p, v.multiplied_scalar(dx))
    }

    /// `GeomFill_BoundWithSurf::HasNormals` (`cxx:70-73`): always true.
    fn has_normals(&self) -> bool {
        true
    }

    /// `GeomFill_BoundWithSurf::Norm` (`cxx:77-100`).
    fn norm(&self, u: f64) -> Result<GpVec, String> {
        // `cxx:88-91`.
        let w = match self.par.lock().unwrap_or_else(|e| e.into_inner()).as_ref() {
            Some(law) => law.value(u),
            None => u,
        };
        let uv = self.pcurve.d0(w);
        let (_, su, sv) = self.surface.d1(uv.x(), uv.y());
        let mut n = su.crossed(&sv);
        n.normalize();
        Ok(n)
    }

    /// `GeomFill_BoundWithSurf::D1Norm` (`cxx:104-158`).
    fn d1_norm(&self, u: f64) -> Result<(GpVec, GpVec), String> {
        // `cxx:114-118`.
        let (w, dw) = match self.par.lock().unwrap_or_else(|e| e.into_inner()).as_ref() {
            Some(law) => law.d1(u),
            None => (u, 1.0),
        };
        // `cxx:120-124`.
        let (p2d, v2d) = self.pcurve.d1(w);
        let scaled = v2d.multiplied_scalar(dw);
        let (dx, dy) = (scaled.x(), scaled.y());
        let (x, y) = (p2d.x(), p2d.y());
        let (_, su, sv, suu, svv, suv) = self.surface.d2(x, y);
        // `cxx:128-132`.
        let mut n = su.crossed(&sv);
        n.normalize();
        let nsuu = n.dot(&suu);
        let nsuv = n.dot(&suv);
        let nsvv = n.dot(&svv);
        let susu = su.dot(&su);
        let susv = su.dot(&sv);
        let svsv = sv.dot(&sv);
        let deno = (susu * svsv) - (susv * susv);
        if deno.abs() < 1.0e-16 {
            // `cxx:135-139`: the approximate fallback.
            let temp = self.norm(u + 1.0e-12)?;
            let dn = n
                .multiplied_scalar(-1.0)
                .added(&temp)
                .multiplied_scalar(1.0e-12);
            Ok((n, dn))
        } else {
            // `cxx:142-156`.
            let a = (-nsuu * svsv + nsuv * susv) / deno;
            let b = (nsuu * susv - nsuv * susu) / deno;
            let c = (-nsuv * svsv + nsvv * susv) / deno;
            let d = (nsuv * susv - nsvv * susu) / deno;
            let temp1 = su.multiplied_scalar(a).added(&sv.multiplied_scalar(b));
            let temp2 = su.multiplied_scalar(c).added(&sv.multiplied_scalar(d));
            let dn = temp1
                .multiplied_scalar(dx)
                .added(&temp2.multiplied_scalar(dy));
            Ok((n, dn))
        }
    }

    /// `GeomFill_BoundWithSurf::Reparametrize` (`cxx:162-175`).
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
        let par = make_par_law(
            self.curve.as_ref(),
            first,
            last,
            has_df,
            has_dl,
            df,
            dl,
            rev,
        )?;
        *self.par.lock().unwrap_or_else(|e| e.into_inner()) = Some(par);
        Ok(())
    }

    /// `GeomFill_BoundWithSurf::Bounds` (`cxx:179-187`).
    fn bounds(&self) -> (f64, f64) {
        match self.par.lock().unwrap_or_else(|e| e.into_inner()).as_ref() {
            Some(law) => LawBSpFunc::bounds(law),
            None => (self.curve.first_parameter(), self.curve.last_parameter()),
        }
    }

    /// `GeomFill_BoundWithSurf::IsDegenerated` (`cxx:191-194`).
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
