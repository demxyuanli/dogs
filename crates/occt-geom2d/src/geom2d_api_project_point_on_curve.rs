//! `Geom2dAPI_ProjectPointOnCurve` — all extremal projections of a 2D point
//! onto a 2D curve.
//!
//! Source: `ModelingAlgorithms/TKGeomAlgo/Geom2dAPI/Geom2dAPI_ProjectPointOnCurve.cxx:22-186`
//! and `.hxx`. The body is a thin wrapper: `Init` loads a `Geom2dAdaptor_Curve`
//! over `[Umin, Usup]` (`cxx:59-62`), runs `Extrema_ExtPC2d(P, myC)`
//! (`cxx:64-68`), sets `myIsDone = IsDone() && NbExt() > 0` (`cxx:70`) and
//! records the index of the smallest `SquareDistance` (`cxx:73-88`). Every
//! accessor is one of those two sequences.
//!
//! [`crate::extrema2d::ExtremaExtPC2d`] is the ported `Extrema_ExtPC2d`, so the
//! range-capable constructor [`ProjectPointOnCurve2d::new_range`] mirrors
//! `Init(P, C, Umin, Usup)` exactly.
//!
//! Divergences that are not heuristics:
//! - OCCT's `myIndex = -1` sentinel becomes `Option<usize>` (`myIndex` is only
//!   read when `myIsDone`, where it is >= 1),
//! - `Point` / `Parameter` / `Distance` raise `Standard_OutOfRange` on an index
//!   outside `1..=NbPoints()`; the port's 1-based accessors return `None`.

use std::sync::Arc;

use occt_core::gp::GpPnt2d;

use crate::curve::Curve2d;
use crate::extrema2d::ExtremaExtPC2d;

/// One `myExtPC` entry: `Point(i).Parameter()`, `Point(i).Value()` and
/// `SquareDistance(i)` (`Geom2dAPI_ProjectPointOnCurve.cxx:64-68`).
#[derive(Debug, Clone, Copy)]
struct Solution {
    parameter: f64,
    point: GpPnt2d,
    square_distance: f64,
}

/// `Geom2dAPI_ProjectPointOnCurve` (`Geom2dAPI_ProjectPointOnCurve.hxx`).
#[derive(Clone, Default)]
pub struct ProjectPointOnCurve2d {
    /// `myC` — the adaptor over `[Umin, Usup]`. OCCT keeps the loaded adaptor so
    /// that a later `Init` can reload it; the port only needs the curve handle
    /// for its parameter range, kept for the record.
    curve: Option<Arc<dyn Curve2d>>,
    /// `myExtPC` — snapshotted, since the extrema object borrows the curve.
    solutions: Vec<Solution>,
    /// `myIsDone` (`cxx:70`).
    is_done: bool,
    /// `myIndex`, 1-based (`cxx:75-87`).
    index: Option<usize>,
}

impl ProjectPointOnCurve2d {
    /// `Geom2dAPI_ProjectPointOnCurve()` (`cxx:26-30`): `myIndex = -1`,
    /// `myIsDone = false`, empty `myExtPC`.
    pub fn new() -> Self {
        Self::default()
    }

    /// `Geom2dAPI_ProjectPointOnCurve(P, C)` (`cxx:34-38`).
    pub fn with_curve(p: &GpPnt2d, curve: &Arc<dyn Curve2d>) -> Self {
        let mut proj = Self::new();
        proj.init(p, curve);
        proj
    }

    /// `Geom2dAPI_ProjectPointOnCurve(P, C, Umin, Usup)` (`cxx:42-47`).
    pub fn with_curve_range(
        p: &GpPnt2d,
        curve: &Arc<dyn Curve2d>,
        u_min: f64,
        u_max: f64,
    ) -> Self {
        let mut proj = Self::new();
        proj.init_range(p, curve, u_min, u_max);
        proj
    }

    /// `Init(P, C)` (`cxx:51-54`): `Init(P, C, C->FirstParameter(),
    /// C->LastParameter())`.
    pub fn init(&mut self, p: &GpPnt2d, curve: &Arc<dyn Curve2d>) {
        let (first, last) = (curve.first_parameter(), curve.last_parameter());
        self.init_range(p, curve, first, last);
    }

    /// `Init(P, C, Umin, Usup)` (`cxx:57-88`).
    pub fn init_range(&mut self, p: &GpPnt2d, curve: &Arc<dyn Curve2d>, u_min: f64, u_max: f64) {
        let ext = ExtremaExtPC2d::new_range(curve.as_ref(), p, u_min, u_max);
        // `myIsDone = myExtPC.IsDone() && (myExtPC.NbExt() > 0);` (`cxx:70`).
        self.is_done = ext.is_done() && ext.nb_ext() > 0;
        self.solutions = if self.is_done {
            ext.solutions()
                .iter()
                .map(|s| Solution {
                    parameter: s.u,
                    point: s.point,
                    square_distance: s.sq_dist,
                })
                .collect()
        } else {
            Vec::new()
        };
        // `cxx:73-88`: the smallest `SquareDistance`, first index winning ties.
        self.index = if self.is_done {
            let mut best = 1usize;
            for i in 2..=self.solutions.len() {
                if self.solutions[i - 1].square_distance
                    < self.solutions[best - 1].square_distance
                {
                    best = i;
                }
            }
            Some(best)
        } else {
            None
        };
        self.curve = Some(Arc::clone(curve));
    }

    /// `IsDone()` — OCCT exposes it through `operator int()` / `NbPoints()`;
    /// kept so callers can distinguish "not run" from "no extremum".
    pub fn is_done(&self) -> bool {
        self.is_done
    }

    /// `NbPoints()` (`cxx:92-102`): `NbExt()` when done, else `0`.
    pub fn nb_points(&self) -> usize {
        if self.is_done {
            self.solutions.len()
        } else {
            0
        }
    }

    /// `Point(Index)` (`cxx:106-111`), `index` 1-based.
    pub fn point(&self, index: usize) -> Option<GpPnt2d> {
        self.solution(index).map(|s| s.point)
    }

    /// `Parameter(Index)` (`cxx:115-121`), `index` 1-based.
    pub fn parameter(&self, index: usize) -> Option<f64> {
        self.solution(index).map(|s| s.parameter)
    }

    /// `Distance(Index)` (`cxx:135-141`): `sqrt(SquareDistance(Index))`.
    pub fn distance(&self, index: usize) -> Option<f64> {
        self.solution(index).map(|s| s.square_distance.sqrt())
    }

    /// `NearestPoint()` (`cxx:145-149`): raises `StdFail_NotDone` when
    /// `!myIsDone`; the port answers `None` in that case.
    pub fn nearest_point(&self) -> Option<GpPnt2d> {
        self.solution(self.index?).map(|s| s.point)
    }

    /// `LowerDistanceParameter()` (`cxx:165-170`).
    pub fn lower_distance_parameter(&self) -> Option<f64> {
        self.solution(self.index?).map(|s| s.parameter)
    }

    /// `LowerDistance()` (`cxx:174-179`).
    pub fn lower_distance(&self) -> Option<f64> {
        self.solution(self.index?)
            .map(|s| s.square_distance.sqrt())
    }

    /// The curve the adaptor was loaded with (OCCT's `myC`).
    pub fn curve(&self) -> Option<&Arc<dyn Curve2d>> {
        self.curve.as_ref()
    }

    fn solution(&self, index: usize) -> Option<&Solution> {
        if index < 1 || index > self.nb_points() {
            return None;
        }
        self.solutions.get(index - 1)
    }
}
