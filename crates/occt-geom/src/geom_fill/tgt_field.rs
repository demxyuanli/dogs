//! `GeomFill_TgtField` and `GeomFill_TgtOnCoons` — the algorithmic tangency
//! fields of a Coons patch.
//!
//! Source: `GeomFill_TgtField.hxx:27-53`, `GeomFill_TgtField.cxx:24-32`
//! (default `IsScalable` false, `Scale` empty) and
//! `GeomFill_TgtOnCoons.cxx:26-185` (OCCT 8.0.0):
//! `Value` `cxx:34-74`, `D1(W)` `cxx:78-132`, `D1(W, T, DT)` `cxx:136-185`.
//!
//! Translation notes:
//! - `GeomFill_TgtField` becomes a trait; `IsScalable`/`Scale` keep their
//!   virtual defaults (`false` / no-op).
//! - `Scale` takes `handle<Law_BSpline>` in OCCT (`hxx:37`), here
//!   [`LawBSpline`] by reference; the method is `&mut self` because a scalable
//!   field would rewrite its own law.
//! - `myK->Bound(I)->Norm(W)` can raise on a boundary without normals; the port
//!   returns `Result` from the trait so the failure surfaces as OCCT's
//!   `Standard_Failure`.

use std::sync::Arc;

use occt_core::gp::GpVec;

use crate::geom_fill::coons_alg_patch::CoonsAlgPatch;
use crate::law::bspline::LawBSpline;

/// `GeomFill_TgtField` (`GeomFill_TgtField.hxx:30-51`).
pub trait TgtField: Send + Sync {
    /// `GeomFill_TgtField::IsScalable` (`cxx:28-31`): false on the root class.
    fn is_scalable(&self) -> bool {
        false
    }

    /// `GeomFill_TgtField::Scale` (`cxx:35-36`): empty on the root class.
    ///
    /// OCCT mutates the field through its handle; the port keeps `&self` and
    /// leaves the mutation to implementors (none of them is scalable yet).
    fn scale(&self, _func: &LawBSpline) {}

    /// `GeomFill_TgtField::Value(W)` (`hxx:40`).
    fn value(&self, w: f64) -> Result<GpVec, String>;

    /// `GeomFill_TgtField::D1(W)` (`hxx:44`).
    fn d1(&self, w: f64) -> Result<GpVec, String>;

    /// `GeomFill_TgtField::D1(W, V, DV)` (`hxx:48`).
    fn d1_with(&self, w: f64) -> Result<(GpVec, GpVec), String>;
}

/// `GeomFill_TgtOnCoons` (`GeomFill_TgtOnCoons.hxx:30-54`).
#[derive(Clone)]
pub struct TgtOnCoons {
    /// `myK` (`hxx:52`).
    patch: Arc<CoonsAlgPatch>,
    /// `ibound` (`hxx:53`).
    ibound: usize,
}

impl TgtOnCoons {
    /// `GeomFill_TgtOnCoons(K, I)` (`cxx:27-31`).
    pub fn new(patch: Arc<CoonsAlgPatch>, ibound: usize) -> Self {
        Self { patch, ibound }
    }

    /// The `(U, V)` pair `Value`/`D1` evaluate the patch at, after the
    /// `cxx:39-55` bound switch. `w` replaces the free parameter in the second
    /// switch (`cxx:57-68`).
    fn parameters(&self, w: f64) -> (f64, f64) {
        let mut u = 0.0;
        let mut v = 0.0;
        match self.ibound {
            // case 0: `myK->Bound(1)->Bounds(V, bid)`.
            0 => v = self.patch.bound(1).bounds().0,
            // case 1: `myK->Bound(0)->Bounds(bid, U)`.
            1 => u = self.patch.bound(0).bounds().1,
            // case 2: `myK->Bound(1)->Bounds(bid, V)`.
            2 => v = self.patch.bound(1).bounds().1,
            // case 3: `myK->Bound(0)->Bounds(U, bid)`.
            3 => u = self.patch.bound(0).bounds().0,
            _ => {}
        }
        // `cxx:57-68`.
        match self.ibound {
            0 | 2 => u = w,
            _ => v = w,
        }
        (u, v)
    }
}

impl TgtField for TgtOnCoons {
    /// `GeomFill_TgtOnCoons::Value` (`cxx:34-74`): the patch derivative
    /// projected onto the plane normal to the boundary's normal.
    fn value(&self, w: f64) -> Result<GpVec, String> {
        let (u, v) = self.parameters(w);
        let tgk = match self.ibound {
            0 | 2 => self.patch.d1_v(u, v),
            _ => self.patch.d1_u(u, v),
        };
        let mut n = self.patch.bound(self.ibound).norm(w)?;
        let scal = tgk.dot(&n);
        n = n.multiplied_scalar(-scal);
        Ok(tgk.added(&n))
    }

    /// `GeomFill_TgtOnCoons::D1(W)` (`cxx:78-132`).
    fn d1(&self, w: f64) -> Result<GpVec, String> {
        let (u, v) = self.parameters(w);
        let tgsc = match self.ibound {
            0 | 2 => self.patch.d1_v(u, v),
            _ => self.patch.d1_u(u, v),
        };
        let dtgsc = self.patch.duv(u, v);
        let (n, dn) = self.patch.bound(self.ibound).d1_norm(w)?;
        // `cxx:117-130`.
        let scal = tgsc.dot(&n);
        let scaln = n.multiplied_scalar(-scal);
        let mut tgsc = tgsc.added(&scaln);
        let scaldn = dn.multiplied_scalar(-scal);
        let scal2 = -dtgsc.dot(&n) - tgsc.dot(&dn);
        let mut temp = n.multiplied_scalar(scal2);
        temp = temp.added(&scaldn);
        tgsc = dtgsc.added(&temp);
        Ok(tgsc)
    }

    /// `GeomFill_TgtOnCoons::D1(W, T, DT)` (`cxx:136-185`).
    fn d1_with(&self, w: f64) -> Result<(GpVec, GpVec), String> {
        let (u, v) = self.parameters(w);
        let tgsc = match self.ibound {
            0 | 2 => self.patch.d1_v(u, v),
            _ => self.patch.d1_u(u, v),
        };
        let dtgsc = self.patch.duv(u, v);
        let (n, dn) = self.patch.bound(self.ibound).d1_norm(w)?;
        // `cxx:174-184`.
        let scal = tgsc.dot(&n);
        let scaln = n.multiplied_scalar(-scal);
        let t = tgsc.added(&scaln);
        let scaldn = dn.multiplied_scalar(-scal);
        let scal2 = -dtgsc.dot(&n) - tgsc.dot(&dn);
        let mut temp = n.multiplied_scalar(scal2);
        temp = temp.added(&scaldn);
        let dt = dtgsc.added(&temp);
        Ok((t, dt))
    }
}
