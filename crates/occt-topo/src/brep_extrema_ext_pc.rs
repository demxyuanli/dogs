//! `BRepExtrema_ExtPC` — the extrema between a vertex and an edge.
//!
//! Source: `ModelingAlgorithms/TKTopAlgo/BRepExtrema/BRepExtrema_ExtPC.cxx:24-53`
//! and `.hxx:24-69` (OCCT 8.0.0).
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | `(V, E)` | `cxx:27-30` | [`ExtPC::new`] |
//! | `Initialize(E)` | `cxx:34-45` | [`ExtPC::initialize`] |
//! | `Perform(V)` | `cxx:49-55` | [`ExtPC::perform`] |
//! | `IsDone` / `NbExt` / `IsMin` / `SquareDistance` / `Parameter` / `Point` | `hxx:39-55` | the same names |
//!
//! # Translation notes
//!
//! * `Initialize` stores the edge's curve with the edge's own parameter range
//!   and `Tol = max(Resolution(min(Tol(E), Confusion)), PConfusion)`
//!   (`cxx:39-43`); `BRepAdaptor_Curve::Resolution` is
//!   [`occt_geom::extrema_pc::resolution`] over the same edge curve.
//! * `Perform` is a no-op when the edge was not geometric (`cxx:51`); the port
//!   keeps the empty solution list.
//! * `Extrema_ExtPC` itself is already ported: [`extrema_ext_pc_range_done`]
//!   returns the solutions together with `IsDone()` (see its doc comment for
//!   why a non-empty list can coexist with `IsDone() == false`).

use occt_core::precision::Precision;
use occt_geom::curve::Curve;
use occt_geom::extrema_pc::{extrema_ext_pc_range_done, ExtPcSolution};

use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Vertex};

/// `BRepExtrema_ExtPC` (`BRepExtrema_ExtPC.hxx:24-69`).
pub struct ExtPC {
    /// `myHC` (`hxx:68`): null while the edge is not geometric.
    my_hc: Option<std::sync::Arc<dyn Curve>>,
    /// The `(U1, U2)` range `myExtPC` was initialised with (`cxx:42-43`).
    my_range: (f64, f64),
    /// `myExtPC`'s solutions (`cxx:50-54`).
    my_sols: Vec<ExtPcSolution>,
    /// `myExtPC.IsDone()`.
    my_done: bool,
}

impl ExtPC {
    /// `BRepExtrema_ExtPC(V, E)` (`cxx:27-30`): `Initialize(E)` then
    /// `Perform(V)`.
    pub fn new(v: &Vertex, e: &Edge) -> Self {
        let mut this = Self {
            my_hc: None,
            my_range: (0.0, 0.0),
            my_sols: Vec::new(),
            my_done: false,
        };
        this.initialize(e);
        this.perform(v);
        this
    }

    /// `Initialize(E)` (`cxx:34-45`).
    pub fn initialize(&mut self, e: &Edge) {
        // `cxx:37-40`: `if (!BRep_Tool::IsGeometric(E)) return;`. The port has
        // no `IsGeometric`; the equivalent test is "the edge has a 3D curve",
        // which is what every later `BRep_Tool::Curve` call would need.
        let curve = match BRepTool::edge_curve(e) {
            Some(c) => c,
            None => return,
        };
        // `cxx:41-42`.
        // `Tol = min(Tol(E), Confusion)` raised to `Resolution(Tol)` and
        // clamped at `PConfusion`. `Extrema_ExtPC::Initialize` stores it as
        // `myTol` (`Extrema_ExtPC.cxx:60-71`); the ported entry point
        // `extrema_ext_pc_range_done` takes the window only, so the value is
        // not consumed (see the module deviations).
        let tol = BRepTool::edge_tolerance(e).min(Precision::CONFUSION);
        let _tol = curve.resolution(tol)
            .max(Precision::PCONFUSION);
        // `cxx:43`.
        let range = BRepTool::edge_parameters(e);
        self.my_hc = Some(curve);
        self.my_range = range;
        // `myExtPC.Initialize(*myHC, U1, U2, Tol)` (`cxx:44`): the range is
        // replayed by `Perform` below.
    }

    /// `Perform(V)` (`cxx:49-55`): solve `Extrema_ExtPC` for the vertex point.
    pub fn perform(&mut self, v: &Vertex) {
        let curve = match &self.my_hc {
            Some(c) => c.clone(),
            None => return,
        };
        // `cxx:52-53`.
        let p = BRepTool::vertex_point(v);
        let (sols, done) = extrema_ext_pc_range_done(curve.as_ref(), &p, self.my_range.0, self.my_range.1);
        self.my_sols = sols;
        self.my_done = done;
    }

    /// `IsDone()` (`hxx:39`).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// `NbExt()` (`hxx:42`).
    pub fn nb_ext(&self) -> usize {
        self.my_sols.len()
    }

    /// `IsMin(N)` (`hxx:45`).
    pub fn is_min(&self, n: usize) -> bool {
        self.my_sols
            .get(n.wrapping_sub(1))
            .map(|s| s.is_min)
            .unwrap_or(false)
    }

    /// `SquareDistance(N)` (`hxx:48`).
    pub fn square_distance(&self, n: usize) -> f64 {
        self.my_sols
            .get(n.wrapping_sub(1))
            .map(|s| s.sq_dist)
            .unwrap_or(0.0)
    }

    /// `Parameter(N)` (`hxx:51`).
    pub fn parameter(&self, n: usize) -> f64 {
        self.my_sols
            .get(n.wrapping_sub(1))
            .map(|s| s.u)
            .unwrap_or(0.0)
    }

    /// `Point(N)` (`hxx:54`).
    pub fn point(&self, n: usize) -> occt_core::gp::GpPnt {
        self.my_sols
            .get(n.wrapping_sub(1))
            .map(|s| s.point)
            .unwrap_or_else(occt_core::gp::GpPnt::zero)
    }
}
