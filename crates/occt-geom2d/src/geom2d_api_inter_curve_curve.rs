//! `Geom2dAPI_InterCurveCurve` — the point/segment intersection of two 2D
//! curves.
//!
//! Source: `ModelingAlgorithms/TKGeomAlgo/Geom2dAPI/Geom2dAPI_InterCurveCurve.cxx:30-210`.
//! `Init(C1, C2, Tol)` (`cxx:56-68`) copies both curves, builds a
//! `Geom2dAdaptor_Curve` over each and runs `Geom2dInt_GInter(AC1, AC2, Tol,
//! Tol)`; `NbPoints` / `Point` / `NbSegments` (`cxx:84-118`) are its
//! `NbPoints` / `Point(i).Value()` / `NbSegments`.
//!
//! [`crate::geom2d_int::Geom2dIntGInter::perform_curves`] is the ported
//! `Geom2dInt_GInter` overload with the same signature, so the body is a
//! straight transcription.
//!
//! UNPORTED: `Segment(theIndex, theCurve1, theCurve2)` (`cxx:147-210`), the only
//! consumer of the `IntRes2d_IntersectionSegment` payload. Its sole OCCT caller
//! (`ChFi2d_FilletAlgo.cxx:210`) has the call commented out with
//! `//not implemented (bug in OCC)`, so no ported code path reaches it.
//!
//! Divergences that are not heuristics:
//! - `Standard_NullObject` on a null curve handle is not reachable in the port
//!   (`Arc<dyn Curve2d>` is never null),
//! - `Point(Index)` raises `Standard_OutOfRange`; the port's 1-based accessor
//!   returns `None`.

use std::sync::Arc;

use occt_core::gp::GpPnt2d;
use occt_core::intres2d::IntRes2dIntersectionPoint;

use crate::curve::Curve2d;
use crate::geom2d_int::Geom2dIntGInter;

/// `Geom2dAPI_InterCurveCurve` (`Geom2dAPI_InterCurveCurve.hxx`).
#[derive(Clone, Default)]
pub struct InterCurveCurve2d {
    /// `myCurve1` (`cxx:62`): `C1->Copy()`.
    curve1: Option<Arc<dyn Curve2d>>,
    /// `myCurve2` (`cxx:63`). `None` is the single-curve `Init(C1, Tol)`.
    curve2: Option<Arc<dyn Curve2d>>,
    /// `myIntersector`.
    intersector: Geom2dIntGInter,
    /// `myIsDone` (`cxx:67`).
    is_done: bool,
}

impl InterCurveCurve2d {
    /// `Geom2dAPI_InterCurveCurve()` (`cxx:32-35`): `myIsDone = false`.
    pub fn new() -> Self {
        Self::default()
    }

    /// `Geom2dAPI_InterCurveCurve(C1, C2, Tol)` (`cxx:37-43`).
    pub fn with_curves(c1: &Arc<dyn Curve2d>, c2: &Arc<dyn Curve2d>, tol: f64) -> Self {
        let mut inter = Self::new();
        inter.init(c1, c2, tol);
        inter
    }

    /// `Init(C1, C2, Tol)` (`cxx:56-68`).
    ///
    /// `Geom2dAdaptor_Curve AC1(C1)` in the port is the curve itself: a
    /// `Curve2d` carries its own `FirstParameter()` / `LastParameter()`, which is
    /// exactly the range the adaptor reports.
    pub fn init(&mut self, c1: &Arc<dyn Curve2d>, c2: &Arc<dyn Curve2d>, tol: f64) {
        self.curve1 = Some(Arc::clone(c1));
        self.curve2 = Some(Arc::clone(c2));
        // `myIntersector = Geom2dInt_GInter(AC1, AC2, Tol, Tol)` (`cxx:66`):
        // `TolConf = Tol` and `Tol = Tol`.
        self.intersector
            .perform_curves(c1.as_ref(), c2.as_ref(), tol, tol);
        self.is_done = self.intersector.is_done();
    }

    /// `Init(C1, Tol)` (`cxx:72-80`), the self-intersection arm:
    /// `Geom2dInt_GInter(AC1, Tol, Tol)` is
    /// `perform_curve(AC1, ComputeDomain(AC1, max(Tol, Tol)), Tol, Tol)`
    /// (`gxx:91-117`). No ported caller needs it yet; it is here so the class is
    /// complete.
    pub fn init_self(&mut self, c1: &Arc<dyn Curve2d>, tol: f64) {
        self.curve1 = Some(Arc::clone(c1));
        self.curve2 = None;
        let domain = self.intersector.compute_domain(c1.as_ref(), tol);
        self.intersector
            .perform_curve(c1.as_ref(), &domain, tol, tol);
        self.is_done = self.intersector.is_done();
    }

    /// `NbPoints()` (`cxx:84-96`).
    pub fn nb_points(&self) -> usize {
        if self.is_done {
            self.intersector.nb_points()
        } else {
            0
        }
    }

    /// `Point(Index)` (`cxx:100-105`), `index` 1-based.
    pub fn point(&self, index: usize) -> Option<GpPnt2d> {
        if index < 1 || index > self.nb_points() {
            return None;
        }
        Some(*self.intersector.point(index).value())
    }

    /// `NbSegments()` (`cxx:109-118`).
    pub fn nb_segments(&self) -> usize {
        if self.is_done {
            self.intersector.nb_segments()
        } else {
            0
        }
    }

    /// The 1-based `IntRes2d_IntersectionPoint` of `Point(Index)`.
    pub fn intersection_point(&self, index: usize) -> Option<IntRes2dIntersectionPoint> {
        if index < 1 || index > self.nb_points() {
            return None;
        }
        Some(self.intersector.point(index))
    }

    /// `myIsDone`.
    pub fn is_done(&self) -> bool {
        self.is_done
    }

    /// `myCurve1` / `myCurve2` (`cxx:85-86`), the copies `Segment` trims.
    pub fn curves(&self) -> (&Option<Arc<dyn Curve2d>>, &Option<Arc<dyn Curve2d>>) {
        (&self.curve1, &self.curve2)
    }
}
