//! `GeomFill_CoonsAlgPatch` — the algorithmic Coons patch over four
//! [`Boundary`] handles.
//!
//! Source: `GeomFill_CoonsAlgPatch.cxx:28-301` (OCCT 8.0.0):
//! constructor `cxx:28-65`, `SetFunc` `cxx:69-73`, `Value` `cxx:87-126`,
//! `D1U` `cxx:130-176`, `D1V` `cxx:180-226`, `DUV` `cxx:230-278`.
//!
//! Translation notes:
//! - `occ::handle<GeomFill_Boundary> bound[4]` is `[Arc<dyn Boundary>; 4]`;
//!   `occ::handle<Law_Function> a[2]` is `[Arc<dyn LawFunction>; 2]`.
//! - The default blending functions are `Law_Linear`: `a[0]` goes 1 → 0 over
//!   `B2`'s bounds, `a[1]` goes 0 → 1 over `B1`'s (`cxx:40-48`).
//! - The argument each boundary is evaluated at is reproduced exactly as in
//!   OCCT, including `D1U` evaluating `bound[0]`/`bound[2]` at `U` while
//!   `Value` evaluates every boundary at `V`, and `D1V`/`DUV` following suit.
//!   These calls are the OCCT control flow; they are not "fixed" here.

use std::sync::{Arc, Mutex};

use occt_core::gp::{GpPnt, GpVec};

use crate::geom_fill::boundary::Boundary;
use crate::law::function::{LawFunction, LawLinear};

/// `GeomFill_CoonsAlgPatch` (`GeomFill_CoonsAlgPatch.hxx:35-82`).
pub struct CoonsAlgPatch {
    /// `bound[4]` (`hxx:79`).
    bound: [Arc<dyn Boundary>; 4],
    /// `c[4]` (`hxx:80`): the averaged corners.
    c: [GpPnt; 4],
    /// `a[2]` (`hxx:81`).
    ///
    /// `SetFunc` is called through a shared handle after the patch has already
    /// been handed to the `GeomFill_TgtOnCoons` fields
    /// (`GeomFill_ConstrainedFilling.cxx:478-479`, `:570-571`), so the port
    /// stores the pair behind a `Mutex` and keeps `&self`.
    a: Mutex<[Arc<dyn LawFunction>; 2]>,
}

impl Clone for CoonsAlgPatch {
    fn clone(&self) -> Self {
        Self {
            bound: self.bound.clone(),
            c: self.c,
            a: Mutex::new(self.laws()),
        }
    }
}

impl CoonsAlgPatch {
    /// `GeomFill_CoonsAlgPatch(B1, B2, B3, B4)` (`cxx:28-65`).
    ///
    /// OCCT performs no bound check: `B1`/`B3` and `B2`/`B4` must have the same
    /// range and be oriented consistently.
    pub fn new(
        b1: Arc<dyn Boundary>,
        b2: Arc<dyn Boundary>,
        b3: Arc<dyn Boundary>,
        b4: Arc<dyn Boundary>,
    ) -> Self {
        let bound = [b1, b2, b3, b4];
        // `cxx:40-48`.
        let (deb1, fin1) = bound[1].bounds();
        let mut law0 = LawLinear::new();
        law0.set(deb1, 1.0, fin1, 0.0);
        let a0: Arc<dyn LawFunction> = Arc::new(law0);

        let (deb0, fin0) = bound[0].bounds();
        let mut law1 = LawLinear::new();
        law1.set(deb0, 0.0, fin0, 1.0);
        let a1: Arc<dyn LawFunction> = Arc::new(law1);

        // `cxx:50-64`: the four averaged corners.
        let c0 = midpoint(&bound[3].value(deb1), &bound[0].value(deb0));
        let c1 = midpoint(&bound[0].value(fin0), &bound[1].value(deb1));
        let c2 = midpoint(&bound[1].value(fin1), &bound[2].value(fin0));
        let c3 = midpoint(&bound[2].value(deb0), &bound[3].value(fin1));
        Self {
            bound,
            c: [c0, c1, c2, c3],
            a: Mutex::new([a0, a1]),
        }
    }

    /// A snapshot of `a[2]`, so no lock is held while the laws evaluate.
    fn laws(&self) -> [Arc<dyn LawFunction>; 2] {
        let a = self.a.lock().unwrap_or_else(|e| e.into_inner());
        [a[0].clone(), a[1].clone()]
    }

    /// `GeomFill_CoonsAlgPatch::SetFunc(f1, f2)` (`cxx:69-73`).
    pub fn set_func(&self, f1: Arc<dyn LawFunction>, f2: Arc<dyn LawFunction>) {
        *self.a.lock().unwrap_or_else(|e| e.into_inner()) = [f1, f2];
    }

    /// `GeomFill_CoonsAlgPatch::Func(f1, f2)` (`cxx:75-81`).
    pub fn funcs(&self) -> (Arc<dyn LawFunction>, Arc<dyn LawFunction>) {
        let a = self.laws();
        (a[0].clone(), a[1].clone())
    }

    /// `GeomFill_CoonsAlgPatch::Func(I)` (`cxx:296-301`).
    pub fn func(&self, i: usize) -> Arc<dyn LawFunction> {
        self.laws()[i].clone()
    }

    /// `GeomFill_CoonsAlgPatch::Value(U, V)` (`cxx:87-126`).
    pub fn value(&self, u: f64, v: f64) -> GpPnt {
        let a = self.laws();
        let a0 = a[0].value(v);
        let a1 = a[1].value(u);
        let a2 = 1.0 - a0;
        let a3 = 1.0 - a1;

        // `cxx:93-104`.
        let mut cor = self.bound[0].value(v).xyz().multiplied(a0);
        cor = cor.added(&self.bound[1].value(v).xyz().multiplied(a1));
        cor = cor.added(&self.bound[2].value(v).xyz().multiplied(a2));
        cor = cor.added(&self.bound[3].value(v).xyz().multiplied(a3));
        // `cxx:109-123`.
        cor = cor.added(&self.c[0].xyz().multiplied(-a0 * a3));
        cor = cor.added(&self.c[1].xyz().multiplied(-a0 * a1));
        cor = cor.added(&self.c[2].xyz().multiplied(-a1 * a2));
        cor = cor.added(&self.c[3].xyz().multiplied(-a2 * a3));
        GpPnt::from_xyz(&cor)
    }

    /// `GeomFill_CoonsAlgPatch::D1U(U, V)` (`cxx:130-176`).
    pub fn d1_u(&self, u: f64, v: f64) -> GpVec {
        let a = self.laws();
        let a0 = a[0].value(v);
        let (_, a1) = a[1].d1(u);
        let a2 = 1.0 - a0;
        let a3 = -a1;

        // `cxx:141-158`.
        let mut cor = self.bound[0].d1(u).1.xyz().multiplied(a0);
        cor = cor.added(&self.bound[1].value(v).xyz().multiplied(a1));
        cor = cor.added(&self.bound[2].d1(u).1.xyz().multiplied(a2));
        cor = cor.added(&self.bound[3].value(v).xyz().multiplied(a3));
        // `cxx:159-173`.
        cor = cor.added(&self.c[0].xyz().multiplied(-a0 * a3));
        cor = cor.added(&self.c[1].xyz().multiplied(-a0 * a1));
        cor = cor.added(&self.c[2].xyz().multiplied(-a1 * a2));
        cor = cor.added(&self.c[3].xyz().multiplied(-a2 * a3));
        GpVec::from_xyz(&cor)
    }

    /// `GeomFill_CoonsAlgPatch::D1V(U, V)` (`cxx:180-226`).
    pub fn d1_v(&self, u: f64, v: f64) -> GpVec {
        let a = self.laws();
        let (_, a0) = a[0].d1(v);
        let a1 = a[1].value(u);
        let a2 = -a0;
        let a3 = 1.0 - a1;

        // `cxx:191-207`: the two boundary *values* use `U`, matching OCCT.
        let mut cor = self.bound[0].value(u).xyz().multiplied(a0);
        cor = cor.added(&self.bound[1].d1(v).1.xyz().multiplied(a1));
        cor = cor.added(&self.bound[2].value(u).xyz().multiplied(a2));
        cor = cor.added(&self.bound[3].d1(v).1.xyz().multiplied(a3));
        // `cxx:208-223`.
        cor = cor.added(&self.c[0].xyz().multiplied(-a0 * a3));
        cor = cor.added(&self.c[1].xyz().multiplied(-a0 * a1));
        cor = cor.added(&self.c[2].xyz().multiplied(-a1 * a2));
        cor = cor.added(&self.c[3].xyz().multiplied(-a2 * a3));
        GpVec::from_xyz(&cor)
    }

    /// `GeomFill_CoonsAlgPatch::DUV(U, V)` (`cxx:230-278`).
    pub fn duv(&self, u: f64, v: f64) -> GpVec {
        let a = self.laws();
        let (_, a0) = a[0].d1(v);
        let (_, a1) = a[1].d1(u);
        let a2 = -a0;
        let a3 = -a1;

        // `cxx:243-258`.
        let mut cor = self.bound[0].d1(u).1.xyz().multiplied(a0);
        cor = cor.added(&self.bound[1].d1(v).1.xyz().multiplied(a1));
        cor = cor.added(&self.bound[2].d1(u).1.xyz().multiplied(a2));
        cor = cor.added(&self.bound[3].d1(v).1.xyz().multiplied(a3));
        // `cxx:259-274`.
        cor = cor.added(&self.c[0].xyz().multiplied(-a0 * a3));
        cor = cor.added(&self.c[1].xyz().multiplied(-a0 * a1));
        cor = cor.added(&self.c[2].xyz().multiplied(-a1 * a2));
        cor = cor.added(&self.c[3].xyz().multiplied(-a2 * a3));
        GpVec::from_xyz(&cor)
    }

    /// `GeomFill_CoonsAlgPatch::Corner(I)` (`cxx:290-294`).
    pub fn corner(&self, i: usize) -> &GpPnt {
        &self.c[i]
    }

    /// `GeomFill_CoonsAlgPatch::Bound(I)` (`cxx:284-288`).
    pub fn bound(&self, i: usize) -> &Arc<dyn Boundary> {
        &self.bound[i]
    }
}

/// `gp_XYZ::Added` + `Multiply(0.5)` as OCCT's corner average.
fn midpoint(a: &GpPnt, b: &GpPnt) -> GpPnt {
    let xyz = a.xyz().added(b.xyz()).multiplied(0.5);
    GpPnt::from_xyz(&xyz)
}
