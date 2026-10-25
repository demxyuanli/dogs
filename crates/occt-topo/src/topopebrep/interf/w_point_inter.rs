//! `TopOpeBRep_WPointInter` — a non-owning view of an `IntSurf_PntOn2S`.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRep/TopOpeBRep_WPointInter.hxx:26-53`,
//! `.cxx:20-82`.
//!
//! OCCT stores `myPP2S` as a raw `IntSurf_PntOn2S*` handed out by
//! `Set(const IntSurf_PntOn2S&)`; the point is a value type that the owner never
//! mutates while the view is alive, so the port stores the value (`Copy`). The
//! accessors keep OCCT's `double&` outputs as tuples.

use occt_core::gp::{GpPnt, GpPnt2d};

use crate::int_tools_wline::PntOn2S;

/// `class TopOpeBRep_WPointInter` (`hxx:26-53`).
#[derive(Debug, Clone, Copy, Default)]
pub struct WPointInter {
    /// `myPP2S` — `None` is the OCCT null pointer.
    pp2s: Option<PntOn2S>,
}

impl WPointInter {
    /// `TopOpeBRep_WPointInter()` (`cxx:22`).
    pub fn new() -> Self {
        Self { pp2s: None }
    }

    /// `Set(const IntSurf_PntOn2S& P)` (`cxx:26-30`).
    pub fn set(&mut self, p: &PntOn2S) {
        self.pp2s = Some(*p);
    }

    /// `ParametersOnS1(U1, V1)` (`cxx:34-37`).
    pub fn parameters_on_s1(&self) -> (f64, f64) {
        match self.pp2s {
            Some(p) => p.parameters_on_s1(),
            // OCCT dereferences a null `myPP2S`; the empty state is not
            // reachable in OCCT either (every reader runs after `Set`).
            None => (0.0, 0.0),
        }
    }

    /// `ParametersOnS2(U2, V2)` (`cxx:41-44`).
    pub fn parameters_on_s2(&self) -> (f64, f64) {
        match self.pp2s {
            Some(p) => p.parameters_on_s2(),
            None => (0.0, 0.0),
        }
    }

    /// `Parameters(U1, V1, U2, V2)` (`cxx:48-51`).
    pub fn parameters(&self) -> (f64, f64, f64, f64) {
        match self.pp2s {
            Some(p) => p.parameters(),
            None => (0.0, 0.0, 0.0, 0.0),
        }
    }

    /// `ValueOnS1()` (`cxx:55-60`): the UV on surface 1 as a 2d point.
    pub fn value_on_s1(&self) -> GpPnt2d {
        let (u, v) = self.parameters_on_s1();
        GpPnt2d::new(u, v)
    }

    /// `ValueOnS2()` (`cxx:64-69`): the UV on surface 2 as a 2d point.
    pub fn value_on_s2(&self) -> GpPnt2d {
        let (u, v) = self.parameters_on_s2();
        GpPnt2d::new(u, v)
    }

    /// `Value()` (`cxx:73-76`).
    pub fn value(&self) -> GpPnt {
        match self.pp2s {
            Some(p) => p.value(),
            None => GpPnt::zero(),
        }
    }

    /// `PPntOn2SDummy()` (`cxx:78-81`).
    pub fn ppnt_on2s(&self) -> Option<&PntOn2S> {
        self.pp2s.as_ref()
    }
}
