//! Port of `IntCurve_IntPolyPolyGen` instantiated as
//! `Geom2dInt_TheIntPCurvePCurveOfGInter`
//! (`IntCurve_IntPolyPolyGen.gxx:1-1797`,
//! `Geom2dInt_TheIntPCurvePCurveOfGInter.hxx:35-102`,
//! `Geom2dInt_TheIntPCurvePCurveOfGInter_0.cxx:27-...`).
//!
//! The general polygon-approximation curve/curve intersector: it discretizes
//! both curves into `Geom2dIntPolygon2d`, runs `Intf_InterferencePolygon2d` on
//! them, refines every interference point with
//! `Geom2dIntExactIntersectionPoint` and turns tangency zones into
//! intersection segments (recursively refining the zone's sub-domain while the
//! polygons are still coarser than `TolConf`).
//!
//! `Perform(C1, D1, C2, D2, TolConf, Tol)` (`gxx:94-292`),
//! `Perform(C1, D1, TolConf, Tol)` (`gxx:296-387`),
//! `Perform(..., NbIter, DeltaU, DeltaV)` (`gxx:390-870` and `:1040-1146`),
//! `findIntersect` (`gxx:1152-1568`), `GetIntersection` (`gxx:1573-1783`),
//! `GetMinNbSamples` / `SetMinNbSamples` (`gxx:1787-1797`).

mod prelude {
    pub(crate) use occt_core::bnd::box2d::BndBox2d;
    pub(crate) use occt_core::gp::GpPnt2d;
    pub(crate) use occt_core::intf::{IntfInterferencePolygon2d, IntfPolygon2d};
    pub(crate) use occt_core::intimpargen::gen::{determine_transition_simple, determine_transition_touch};
    pub(crate) use occt_core::intres2d::{
        IntRes2dDomain, IntRes2dIntersection, IntRes2dIntersectionPoint,
        IntRes2dIntersectionSegment, IntRes2dPosition, IntRes2dTransition,
    };
    pub(crate) use occt_core::precision::{epsilon, PCONFUSION};
    pub(crate) use crate::curve::Curve2d;
    pub(crate) use super::super::exact_intersection_point::Geom2dIntExactIntersectionPoint;
    pub(crate) use super::super::polygon2d::Geom2dIntPolygon2d;
}

use prelude::*;


/// `NBITER_MAX_POLYGON` (`gxx:50`).
const NBITER_MAX_POLYGON: i32 = 10;
/// `TOL_CONF_MINI` (`gxx:51`).
const TOL_CONF_MINI: f64 = 0.0000000001;
/// `TOL_MINI` (`gxx:52`).
const TOL_MINI: f64 = 0.0000000001;

/// `IntCurve_IntPolyPolyGen` (`Geom2dInt_TheIntPCurvePCurveOfGInter.hxx:35-102`),
/// carrying the `IntRes2d_Intersection` result.
#[derive(Clone, Debug)]
pub struct Geom2dIntIntPolyPolyGen {
    /// The `IntRes2d_Intersection` base class.
    pub(crate) base: IntRes2dIntersection,
    /// `DomainOnCurve1` (`hxx:96`).
    domain_on_curve1: IntRes2dDomain,
    /// `DomainOnCurve2` (`hxx:97`).
    domain_on_curve2: IntRes2dDomain,
    /// `myMinPntNb` (`hxx:100`).
    min_pnt_nb: usize,
}

impl Default for Geom2dIntIntPolyPolyGen {
    fn default() -> Self {
        Self::new()
    }
}
impl Geom2dIntIntPolyPolyGen {
    /// `IntCurve_IntPolyPolyGen()` (`gxx:86-90`).
    pub fn new() -> Self {
        // Minimum number of samples (`gxx:87`).
        let min_pnt_nb_val = 20;
        Self {
            base: IntRes2dIntersection::new(),
            domain_on_curve1: IntRes2dDomain::new(),
            domain_on_curve2: IntRes2dDomain::new(),
            min_pnt_nb: min_pnt_nb_val,
        }
    }

    /// The `IntRes2d_Intersection` base (`Geom2dInt_TheIntPCurvePCurveOfGInter`
    /// derives from it).
    pub fn result(&self) -> &IntRes2dIntersection {
        &self.base
    }

    /// `IntRes2d_Intersection::ResetFields` (inherited, `gxx:100`).
    fn reset_fields(&mut self) {
        self.base.reset_fields();
    }

    /// `GetMinNbSamples()` (`gxx:1787-1790`).
    pub fn min_nb_samples(&self) -> usize {
        self.min_pnt_nb
    }

    /// `SetMinNbSamples(theMinNbSamples)` (`gxx:1794-1797`).
    pub fn set_min_nb_samples(&mut self, min_nb_samples: usize) {
        self.min_pnt_nb = min_nb_samples;
    }
}

mod find_intersect;
mod geometry;
mod get_intersection;
mod perform;
mod perform_poly;
