//! `BRepExtrema_ExtCC` — the extrema between two edges.
//!
//! Source: `ModelingAlgorithms/TKTopAlgo/BRepExtrema/BRepExtrema_ExtCC.cxx:26-119`
//! and `.hxx:24-75` (OCCT 8.0.0).
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | `(E1, E2)` | `cxx:26-29` | [`ExtCC::new`] |
//! | `Initialize(E2)` | `cxx:33-47` | [`ExtCC::initialize`] |
//! | `Perform(E1)` | `cxx:52-69` | [`ExtCC::perform`] |
//! | `IsDone` / `IsParallel` / `NbExt` / `SquareDistance` | `hxx:43-60` | the same names |
//! | `ParameterOnE1/2` / `PointOnE1/2` | `cxx:73-105` | the same names |
//!
//! # Translation notes
//!
//! * `myExtCC` here is OCCT's `Extrema_ExtCC`, which dispatches `Extrema_ECC`
//!   (analytic `Extrema_ExtElC` arms + `math_FunctionSetRoot` from defined
//!   starts) and exposes `SetCurve` / `SetTolerance` / `Perform`. The port's
//!   equivalent engine is [`occt_geom::extrema_cc::GGenExtCC`]
//!   (`Extrema_GenExtCC`), which carries the same `IsDone` / `IsParallel` /
//!   `NbExt` / `SquareDistance` / `Points` surface. See `extrema_cc`'s module
//!   doc: the general-path **seed set** there is an UNPORTED stand-in for
//!   `math_GlobOptMin`, so this arm inherits that deviation.
//! * `BRepAdaptor_Curve::Resolution(Tol)` (`cxx:43` / `cxx:60`) is
//!   [`occt_geom::extrema_pc::resolution`] over the edge's own curve.
//! * `IsGeometric` (`cxx:36` / `cxx:54`) is the "edge has a 3D curve" test, as
//!   in [`super::brep_extrema_ext_pc`].
//! * **UNPORTED**: `TrimmedSquareDistances` (`cxx:108-119`). The ported
//!   `GGenExtCC` does not keep `Extrema_ECC`'s `myTrDist` / `myP1` / `myP2`
//!   members, so there is no equal branch to read them from. No ported caller
//!   of `BRepExtrema_ExtCC` uses them
//!   (`BRepExtrema_DistanceSS.cxx:967-990` reads only the `Points` /
//!   `SquareDistance` pair).

use std::sync::Arc;

use occt_core::gp::GpPnt;
use occt_core::precision::Precision;
use occt_geom::curve::Curve;
use occt_geom::extrema_cc::GGenExtCC;

use crate::brep_tool::BRepTool;
use crate::shape::Edge;

/// `BRepExtrema_ExtCC` (`BRepExtrema_ExtCC.hxx:24-75`).
pub struct ExtCC {
    /// `myHC` (`hxx:73`): E2's curve, `None` for the `IsGeometric` guard.
    my_hc_e2: Option<Arc<dyn Curve>>,
    /// E2's `(V1, V2)` range (`cxx:45-46`).
    my_range_e2: (f64, f64),
    /// E2's `Tol` (`cxx:43-44`).
    my_tol_e2: f64,
    /// E1's `(U1, U2)` range (`cxx:63-64`).
    my_range_e1: (f64, f64),
    /// E1's `Tol` (`cxx:60-61`).
    my_tol_e1: f64,
    /// `myExtCC`'s `mySqDist` (`cxx:68`), kept as a sequence as in
    /// `Extrema_ExtCC`.
    my_sq_dist: Vec<f64>,
    /// `myExtCC::Points(N)` as `(u1, p1, u2, p2)`.
    my_points: Vec<(f64, GpPnt, f64, GpPnt)>,
    /// `myExtCC.IsDone()`.
    my_done: bool,
    /// `myExtCC.IsParallel()`.
    my_parallel: bool,
}

impl ExtCC {
    /// `BRepExtrema_ExtCC(E1, E2)` (`cxx:26-29`).
    pub fn new(e1: &Edge, e2: &Edge) -> Self {
        let mut this = Self {
            my_hc_e2: None,
            my_range_e2: (0.0, 0.0),
            my_tol_e2: 0.0,
            my_range_e1: (0.0, 0.0),
            my_tol_e1: 0.0,
            my_sq_dist: Vec::new(),
            my_points: Vec::new(),
            my_done: false,
            my_parallel: false,
        };
        this.initialize(e2);
        this.perform(e1);
        this
    }

    /// `Initialize(E2)` (`cxx:33-47`).
    pub fn initialize(&mut self, e2: &Edge) {
        // `cxx:36-39`.
        let curve = match BRepTool::edge_curve(e2) {
            Some(c) => c,
            None => return,
        };
        // `cxx:42-44`.
        let tol = BRepTool::edge_tolerance(e2).min(Precision::CONFUSION);
        let tol = curve.resolution(tol).max(Precision::PCONFUSION);
        // `cxx:45-47`.
        self.my_range_e2 = BRepTool::edge_parameters(e2);
        self.my_hc_e2 = Some(curve);
        self.my_tol_e2 = tol;
    }

    /// `Perform(E1)` (`cxx:52-69`).
    pub fn perform(&mut self, e1: &Edge) {
        // `cxx:54-57`.
        let curve1 = match BRepTool::edge_curve(e1) {
            Some(c) => c,
            None => return,
        };
        // `cxx:59-61`.
        let tol = BRepTool::edge_tolerance(e1).min(Precision::CONFUSION);
        let tol = curve1.resolution(tol).max(Precision::PCONFUSION);
        // `cxx:62-65`.
        self.my_range_e1 = BRepTool::edge_parameters(e1);
        self.my_tol_e1 = tol;

        // `cxx:66-68`: `SetCurve(2, *myHC, V1, V2)`, `SetTolerance(2, Tol)`,
        // `SetCurve(1, *HC, U1, U2)`, `SetTolerance(1, Tol)`, `Perform()`.
        let curve2 = match &self.my_hc_e2 {
            Some(c) => c.clone(),
            None => return,
        };
        let mut ext = GGenExtCC::new(
            curve1.as_ref(),
            curve2.as_ref(),
            self.my_range_e1,
            self.my_range_e2,
        );
        ext.set_tolerance(self.my_tol_e1.max(self.my_tol_e2));
        ext.perform().ok();

        self.my_done = ext.is_done();
        self.my_parallel = ext.is_parallel();
        self.my_sq_dist.clear();
        self.my_points.clear();
        for i in 1..=ext.nb_ext() {
            self.my_sq_dist.push(ext.square_distance(i));
            self.my_points.push(ext.points(i));
        }
    }

    /// `IsDone()` (`hxx:43`).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// `IsParallel()` (`hxx:46`).
    pub fn is_parallel(&self) -> bool {
        self.my_parallel
    }

    /// `NbExt()` (`hxx:49`).
    pub fn nb_ext(&self) -> usize {
        self.my_sq_dist.len()
    }

    /// `SquareDistance(N)` (`hxx:52`).
    pub fn square_distance(&self, n: usize) -> f64 {
        self.my_sq_dist.get(n.wrapping_sub(1)).copied().unwrap_or(0.0)
    }

    /// `ParameterOnE1(N)` (`cxx:73-77`).
    pub fn parameter_on_e1(&self, n: usize) -> f64 {
        self.my_points
            .get(n.wrapping_sub(1))
            .map(|p| p.0)
            .unwrap_or(0.0)
    }

    /// `PointOnE1(N)` (`cxx:82-86`).
    pub fn point_on_e1(&self, n: usize) -> GpPnt {
        self.my_points
            .get(n.wrapping_sub(1))
            .map(|p| p.1)
            .unwrap_or_else(GpPnt::zero)
    }

    /// `ParameterOnE2(N)` (`cxx:91-95`).
    pub fn parameter_on_e2(&self, n: usize) -> f64 {
        self.my_points
            .get(n.wrapping_sub(1))
            .map(|p| p.2)
            .unwrap_or(0.0)
    }

    /// `PointOnE2(N)` (`cxx:100-104`).
    pub fn point_on_e2(&self, n: usize) -> GpPnt {
        self.my_points
            .get(n.wrapping_sub(1))
            .map(|p| p.3)
            .unwrap_or_else(GpPnt::zero)
    }
}
