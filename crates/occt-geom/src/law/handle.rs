//! A shared, mutable `Law_Composite` seen through [`LawFunction`].
//!
//! Source: `TKGeomAlgo/Law/Law_Composite.hxx`/`.cxx` (OCCT 8.0.0) — more
//! precisely, the way OCCT hands a law to its consumers.
//!
//! In OCCT every law is an `occ::handle<Law_Function>`: `ChFiDS_FilSpine` owns
//! the composite of a guideline (`ChFiDS_FilSpine.cxx:858-890`) and
//! `ChFi3d_FilBuilder::SimulSurf` passes that very handle into
//! `BRepBlend_EvolRad` (`ChFi3d_FilBuilder.cxx:666-667`). Both objects then see
//! every later `ChangeElementaryLaw` — the handle is shared, not copied.
//!
//! The port's [`LawFunction`] is a trait object, and its concrete owners store a
//! `LawComposite` behind `Arc<RwLock<_>>` (see
//! [`crate::law::LawComposite`]). This adapter is the handle: it borrows
//! nothing, forwards every `LawFunction` method to the current composite, and
//! therefore gives the consumer the same live view OCCT's handle gives.
//!
//! `Trim` answers a *new* law in OCCT (`Law_Composite::Trim` builds a trimmed
//! composite), so the adapter trims the composite and wraps the result — the
//! returned law is no longer tied to the owner, exactly like OCCT.

use std::sync::{Arc, RwLock};

use occt_core::kernel::geomabs::Shape;

use super::composite::LawComposite;
use super::function::LawFunction;

/// `occ::handle<Law_Composite>` — a shared, interior-mutable law handle.
pub struct LawHandle {
    inner: Arc<RwLock<LawComposite>>,
}

impl LawHandle {
    /// Wrap an already shared composite (the one `ChFiDS_FilSpine::Law`
    /// answers).
    pub fn new(inner: Arc<RwLock<LawComposite>>) -> Self {
        Self { inner }
    }

    /// The underlying handle, for the callers that need to change the law
    /// through it (`Law_Composite::ChangeElementaryLaw`).
    pub fn composite(&self) -> Arc<RwLock<LawComposite>> {
        Arc::clone(&self.inner)
    }
}

impl LawFunction for LawHandle {
    /// `Law_Function::Continuity` (`Law_Composite.cxx:56-112`): the weakest
    /// continuity of the components.
    fn continuity(&self) -> Shape {
        let guard = self.inner.read().unwrap_or_else(|e| e.into_inner());
        guard.continuity()
    }

    /// `Law_Function::NbIntervals(S)` (`Law_Composite.cxx:114-156`).
    fn nb_intervals(&self, s: Shape) -> Result<i32, String> {
        let guard = self.inner.read().unwrap_or_else(|e| e.into_inner());
        guard.nb_intervals(s)
    }

    /// `Law_Function::Intervals(T, S)` (`Law_Composite.cxx:158-204`).
    fn intervals(&self, s: Shape) -> Result<Vec<f64>, String> {
        let guard = self.inner.read().unwrap_or_else(|e| e.into_inner());
        guard.intervals(s)
    }

    /// `Law_Function::Value(X)` (`Law_Composite.cxx:206-212`).
    fn value(&self, x: f64) -> f64 {
        let guard = self.inner.read().unwrap_or_else(|e| e.into_inner());
        guard.value(x)
    }

    /// `Law_Function::D1(X, F, D)` (`Law_Composite.cxx:228-240`).
    fn d1(&self, x: f64) -> (f64, f64) {
        let guard = self.inner.read().unwrap_or_else(|e| e.into_inner());
        guard.d1(x)
    }

    /// `Law_Function::D2(X, F, D, D2)` (`Law_Composite.cxx:242-254`).
    fn d2(&self, x: f64) -> (f64, f64, f64) {
        let guard = self.inner.read().unwrap_or_else(|e| e.into_inner());
        guard.d2(x)
    }

    /// `Law_Function::Trim(PFirst, PLast, Tol)` (`Law_Composite.cxx:256-296`):
    /// a new law, no longer sharing the owner's composite.
    fn trim(&self, p_first: f64, p_last: f64, tol: f64) -> Result<Arc<dyn LawFunction>, String> {
        let guard = self.inner.read().unwrap_or_else(|e| e.into_inner());
        guard.trim(p_first, p_last, tol)
    }

    /// `Law_Function::Bounds(PFirst, PLast)` (`Law_Composite.cxx:298-306`).
    fn bounds(&self) -> (f64, f64) {
        let guard = self.inner.read().unwrap_or_else(|e| e.into_inner());
        guard.bounds()
    }
}
