//! `ApproxInt_SvSurfaces` — the root class of the surface-pair evaluators.
//!
//! Source: `ModelingAlgorithms/TKGeomAlgo/ApproxInt/ApproxInt_SvSurfaces.hxx:22-101`,
//! `.cxx:5-19` (OCCT 8.0.0).
//!
//! The class carries the `myUseSolver` flag (`hxx:96-99`) and five pure virtual
//! evaluators. OCCT's `myUseSolver` documentation (`hxx:70-78`) is the contract
//! the implementations honour: with `true` the four parameters are a first
//! approximation refined by the intersection solver; with `false` they are
//! taken as exact.
//!
//! `myUseSolver` is a plain field with inline `SetUseSolver` / `GetUseSolver`
//! (`hxx:94-99`). The callers reach the evaluator through an untyped pointer
//! (`ApproxInt_MultiLine.gxx:233`, `:467`, `:632`), so the port stores it as
//! `Rc<RefCell<dyn SvSurfaces>>` and the two accessors stay on the trait.

use occt_core::gp::{GpPnt, GpVec, GpVec2d};

use crate::int_tools_wline::PntOn2S;

/// `class ApproxInt_SvSurfaces` (`hxx:38`).
pub trait SvSurfaces {
    /// `Compute(u1, v1, u2, v2, Pt, Tg, Tguv1, Tguv2)` (`hxx:52-59`): the point
    /// and the 3D / 2D tangents of the intersection line, refining the four
    /// parameters in place. `false` when the point is tangent or the solver
    /// fails.
    fn compute(
        &mut self,
        u1: &mut f64,
        v1: &mut f64,
        u2: &mut f64,
        v2: &mut f64,
        pt: &mut GpPnt,
        tg: &mut GpVec,
        tguv1: &mut GpVec2d,
        tguv2: &mut GpVec2d,
    ) -> bool;

    /// `Pnt(u1, v1, u2, v2, P)` (`hxx:61-66`).
    fn pnt(&mut self, u1: f64, v1: f64, u2: f64, v2: f64, p: &mut GpPnt);

    /// `SeekPoint(u1, v1, u2, v2, Point)` (`hxx:69-74`).
    fn seek_point(&mut self, u1: f64, v1: f64, u2: f64, v2: f64, point: &mut PntOn2S) -> bool;

    /// `Tangency(u1, v1, u2, v2, Tg)` (`hxx:76-81`).
    fn tangency(&mut self, u1: f64, v1: f64, u2: f64, v2: f64, tg: &mut GpVec) -> bool;

    /// `TangencyOnSurf1(u1, v1, u2, v2, Tg)` (`hxx:83-88`).
    fn tangency_on_surf1(&mut self, u1: f64, v1: f64, u2: f64, v2: f64, tg: &mut GpVec2d) -> bool;

    /// `TangencyOnSurf2(u1, v1, u2, v2, Tg)` (`hxx:90-95`).
    fn tangency_on_surf2(&mut self, u1: f64, v1: f64, u2: f64, v2: f64, tg: &mut GpVec2d) -> bool;

    /// `SetUseSolver(theUseSol)` (`hxx:96`).
    fn set_use_solver(&mut self, use_solver: bool);

    /// `GetUseSolver()` (`hxx:98`).
    fn get_use_solver(&self) -> bool;
}

/// The `myUseSolver` field of `ApproxInt_SvSurfaces` (`hxx:100`), whose
/// initial value is `false` (`hxx:50-52`).
///
/// OCCT declares it in the base class, so every implementation shares this one
/// field shape; implementations embed [`SolverFlag`] and forward the two
/// accessors to it.
#[derive(Debug, Clone, Copy, Default)]
pub struct SolverFlag {
    use_solver: bool,
}

impl SolverFlag {
    /// `SetUseSolver(theUseSol)` (`hxx:96`).
    pub fn set_use_solver(&mut self, use_solver: bool) {
        self.use_solver = use_solver;
    }

    /// `GetUseSolver()` (`hxx:98`).
    pub fn get_use_solver(&self) -> bool {
        self.use_solver
    }
}
