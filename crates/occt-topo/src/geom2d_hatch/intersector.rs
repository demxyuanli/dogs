//! `Geom2dHatch_Intersector`: the curve/line intersector the hatcher trims
//! hatchings with, plus the local-geometry helper its transition code needs.
//!
//! Sources: `ModelingAlgorithms/TKGeomAlgo/Geom2dHatch/Geom2dHatch_Intersector.hxx:27-85`,
//! `.lxx:18-78`, `.cxx:26-104` (OCCT 8.0.0).
//!
//! `Geom2dHatch_Intersector` derives from `Geom2dInt_GInter`, so the port keeps
//! a [`Geom2dIntGInter`] member and re-exports the `IntRes2d_Intersection`
//! accessors the hatcher reads from it.

use occt_core::gp::{GpDir2d, GpLin2d};
use occt_core::intres2d::{
    IntRes2dDomain, IntRes2dIntersectionPoint, IntRes2dIntersectionSegment,
};
use occt_core::precision::{PCONFUSION, PINTERSECTION};
use occt_geom2d::cl_props2d::GeomLPropCLProps2d;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::geom2d_int::Geom2dIntGInter;
use occt_geom2d::line::Geom2dLine;

use crate::geom_bnd_lib_elclib2d::line_value_2d;

/// `RealLast()` (`Standard_Real.hxx:179-182`).
const REAL_LAST: f64 = f64::MAX;

/// `Geom2dHatch_Intersector` (`hxx:27-85`).
#[derive(Clone, Debug)]
pub struct Intersector {
    /// The inherited `Geom2dInt_GInter` part.
    base: Geom2dIntGInter,
    /// `myConfusionTolerance`
    my_confusion_tolerance: f64,
    /// `myTangencyTolerance`
    my_tangency_tolerance: f64,
}

impl Default for Intersector {
    /// `Geom2dHatch_Intersector()` (`cxx:26-31`): both tolerances are null.
    fn default() -> Self {
        Self::new(0.0, 0.0)
    }
}

impl Intersector {
    /// `Geom2dHatch_Intersector(Confusion, Tangency)` (`.lxx:19-23`).
    pub fn new(confusion: f64, tangency: f64) -> Self {
        Self {
            base: Geom2dIntGInter::new(),
            my_confusion_tolerance: confusion,
            my_tangency_tolerance: tangency,
        }
    }

    /// `ConfusionTolerance()` (`.lxx:31-34`).
    pub fn confusion_tolerance(&self) -> f64 {
        self.my_confusion_tolerance
    }

    /// `SetConfusionTolerance(Confusion)` (`.lxx:42-45`).
    pub fn set_confusion_tolerance(&mut self, confusion: f64) {
        self.my_confusion_tolerance = confusion;
    }

    /// `TangencyTolerance()` (`.lxx:52-55`).
    pub fn tangency_tolerance(&self) -> f64 {
        self.my_tangency_tolerance
    }

    /// `SetTangencyTolerance(Tangency)` (`.lxx:62-65`).
    pub fn set_tangency_tolerance(&mut self, tangency: f64) {
        self.my_tangency_tolerance = tangency;
    }

    /// `Intersect(C1, C2)` (`.lxx:71-78`):
    /// `Geom2dInt_GInter::Perform(C1, C2, myConfusionTolerance,
    /// myTangencyTolerance)`.
    pub fn intersect(&mut self, c1: &dyn Curve2d, c2: &dyn Curve2d) {
        self.base.perform_curves(
            c1,
            c2,
            self.my_confusion_tolerance,
            self.my_tangency_tolerance,
        );
    }

    /// `Perform(L, P, Tol, E)` (`cxx:35-70`): intersects the `[0, P]` part of
    /// `L` with the curve `E`. A null `P` domain is semi-infinite on `[0, +inf[`
    /// only when `P == RealLast()`; otherwise the domain ends at
    /// `ElCLib::Value(P, L)`. The edge domain always runs over the whole
    /// parameter range with `Precision::PIntersection()`, and the inter-sector
    /// is run with `Precision::PConfusion()` / `PIntersection()`.
    pub fn perform(&mut self, l: &GpLin2d, p: f64, tol: f64, c: &dyn Curve2d) {
        let dl = if p != REAL_LAST {
            IntRes2dDomain::bounded(&l.location(), 0., tol, &line_value_2d(p, l), p, tol)
        } else {
            IntRes2dDomain::semi_infinite(&l.location(), 0., tol, true)
        };

        let de = IntRes2dDomain::bounded(
            &c.d0(c.first_parameter()),
            c.first_parameter(),
            PINTERSECTION,
            &c.d0(c.last_parameter()),
            c.last_parameter(),
            PINTERSECTION,
        );

        let gl = Geom2dLine::from_pnt_dir(l.location(), *l.direction());
        // `cxx:60-69`: OCCT builds a temporary `Geom2dInt_GInter` and calls
        // `SetValues(Inter)` on itself; running the same `Perform` in place on
        // the inherited base is that same sequence.
        self.base
            .perform(&gl, &dl, c, &de, PCONFUSION, PINTERSECTION);
    }

    /// `LocalGeometry(E, U, T, N, C)` (`cxx:74-104`): tangent, normal and
    /// curvature of `E` at `U` from a `GeomLProp_CLProps2d` with order 2 and
    /// `Precision::PConfusion()`. When the tangent is undefined it falls back to
    /// `(1, 0)` and a null curvature; the normal is the analytic one only when
    /// the curvature is in `]PConfusion(), RealLast()[`, else the tangent
    /// rotated by -90 degrees.
    pub fn local_geometry(
        &self,
        e: &dyn Curve2d,
        u: f64,
        tang: &mut GpDir2d,
        norm: &mut GpDir2d,
        c: &mut f64,
    ) {
        let mut prop = GeomLPropCLProps2d::new(e, u, 2, PCONFUSION);

        *c = 0.;
        if prop.is_tangent_defined() {
            *tang = prop.tangent();
            *c = prop.curvature();
        } else {
            *tang = GpDir2d::new(1., 0.).expect("Geom2dHatch_Intersector::LocalGeometry: X dir");
        }

        if *c > PCONFUSION && *c < REAL_LAST {
            *norm = prop.normal();
        } else {
            *norm = GpDir2d::new(tang.y(), -tang.x())
                .expect("Geom2dHatch_Intersector::LocalGeometry: null normal");
        }
    }

    /// `Geom2dInt_GInter::IsDone()` (inherited), read by the hatcher.
    pub fn is_done(&self) -> bool {
        self.base.is_done()
    }

    /// `Geom2dInt_GInter::IsEmpty()` (inherited).
    pub fn is_empty(&self) -> bool {
        self.base.result().is_empty()
    }

    /// `Geom2dInt_GInter::NbPoints()` (inherited).
    pub fn nb_points(&self) -> usize {
        self.base.nb_points()
    }

    /// `Geom2dInt_GInter::Point(Index)` (inherited), 1-based.
    pub fn point(&self, index: usize) -> IntRes2dIntersectionPoint {
        self.base.point(index - 1)
    }

    /// `Geom2dInt_GInter::NbSegments()` (inherited).
    pub fn nb_segments(&self) -> usize {
        self.base.nb_segments()
    }

    /// `Geom2dInt_GInter::Segment(Index)` (inherited), 1-based.
    pub fn segment(&self, index: usize) -> IntRes2dIntersectionSegment {
        self.base.segment(index - 1)
    }
}
