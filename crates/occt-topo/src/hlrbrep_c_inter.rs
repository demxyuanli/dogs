//! `HLRBRep_CInter` (`ModelingAlgorithms/TKHLR/HLRBRep/HLRBRep_CInter.hxx:43-179`).
//!
//! `HLRBRep_CInter` is `IntCurve_IntCurveCurveGen` instantiated on
//! `HLRBRep_CurveTool` + `Geom2dInt_GInter` (`hxx:232-233` defines
//! `IntCurve_IntCurveCurveGen` as `HLRBRep_CInter`), and it adds nothing to the
//! generic class: the whole header is the six `Perform` overloads, the six
//! constructor twins, `ComputeDomain` and `Set/GetMinNbSamples`
//! (`int_curve_gen` / `IntCurve_IntCurveCurveGen.lxx:91-151`,
//! `gxx:1023-1031`).
//!
//! The port's `Geom2dIntGInter`
//! (`crates/occt-geom2d/src/geom2d_int/ginter.rs`) is that same generic class
//! instantiated on the `Geom2d` curve tool, already working against
//! `&dyn Curve2d`. Two differences separate the two instantiations and both are
//! handled outside this file:
//!
//! * `TheCurveTool::GetType` - `Geom2dIntGInter` asks the `Curve2d` for
//!   `is_line` / `gp_circ2d` / ..., and `impl Curve2d for HlrBRepCurve`
//!   (`crate::hlrbrep_curve_2d`) answers from `HLRBRep_Curve::GetType()`, which
//!   is `HLRBRep_CurveTool::GetType` (`HLRBRep_CurveTool.lxx:173-176`);
//! * `TheCurveTool::NbSamples` - `Geom2dIntGInter` routes the polygon sample
//!   counts through `Curve2d::nb_samples_2d` / `nb_samples_2d_full`, which the
//!   `HlrBRepCurve` impl overrides with
//!   `HLRBRep_CurveTool::NbSamples` (`HLRBRep_CurveTool.cxx:22-83`).
//!
//! `SetMinNbSamples` is the one piece the generic class keeps in its own
//! `intcurvcurv` member (`gxx:1023-1026`); the port keeps it on
//! `Geom2dIntGInter` because that class builds the poly-poly intersector per
//! `Perform`.

use occt_core::intres2d::{
    IntRes2dDomain, IntRes2dIntersectionPoint, IntRes2dIntersectionSegment,
};
use occt_geom2d::geom2d_int::Geom2dIntGInter;

use crate::hlrbrep_curve::HlrBRepCurve;

/// `HLRBRep_CInter` (`HLRBRep_CInter.hxx:43-179`).
#[derive(Clone)]
pub struct HlrBRepCInter {
    /// `IntCurve_IntPolyPolyGen intcurvcurv`'s owning instance: the generic
    /// base is `IntRes2d_Intersection` (`hxx:45`), which the port keeps inside
    /// `Geom2dIntGInter`.
    base: Geom2dIntGInter,
}

impl Default for HlrBRepCInter {
    fn default() -> Self {
        Self::new()
    }
}

impl HlrBRepCInter {
    /// `HLRBRep_CInter()` (`IntCurve_IntCurveCurveGen.lxx:40-53`): the empty
    /// constructor; `HLRBRep_Intersector` then calls `SetMinNbSamples(4)`
    /// (`HLRBRep_Intersector.cxx:86-87`).
    pub fn new() -> Self {
        Self {
            base: Geom2dIntGInter::new(),
        }
    }

    /// `SetMinNbSamples(theMinNbSamples)` (`gxx:1023-1026`).
    pub fn set_min_nb_samples(&mut self, the_min_nb_samples: usize) {
        self.base.set_min_nb_samples(the_min_nb_samples);
    }

    /// `GetMinNbSamples()` (`gxx:1030-1033`).
    pub fn min_nb_samples(&self) -> usize {
        self.base.min_nb_samples()
    }

    /// `ComputeDomain(C1, TolDomain)` (`gxx:120-177`).
    pub fn compute_domain(&self, c: &HlrBRepCurve, tol_domain: f64) -> IntRes2dDomain {
        self.base.compute_domain(c, tol_domain)
    }

    /// `Perform(C1, C2, TolConf, Tol)` (`lxx:110-121`).
    pub fn perform_curves(
        &mut self,
        c1: &HlrBRepCurve,
        c2: &HlrBRepCurve,
        tol_conf: f64,
        tol: f64,
    ) {
        self.base.perform_curves(c1, c2, tol_conf, tol);
    }

    /// `Perform(C1, D1, C2, D2, TolConf, Tol)` (`gxx:182-225`).
    pub fn perform(
        &mut self,
        c1: &HlrBRepCurve,
        d1: &IntRes2dDomain,
        c2: &HlrBRepCurve,
        d2: &IntRes2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        self.base.perform(c1, d1, c2, d2, tol_conf, tol);
    }

    /// `Perform(C1, D1, C2, TolConf, Tol)` (`lxx:153-166`).
    pub fn perform_with_d1(
        &mut self,
        c1: &HlrBRepCurve,
        d1: &IntRes2dDomain,
        c2: &HlrBRepCurve,
        tol_conf: f64,
        tol: f64,
    ) {
        let d2 = self.base.compute_domain(c2, tol_conf.max(tol));
        self.base.perform(c1, d1, c2, &d2, tol_conf, tol);
    }

    /// `Perform(C1, C2, D2, TolConf, Tol)` (`lxx:139-151`).
    pub fn perform_with_d2(
        &mut self,
        c1: &HlrBRepCurve,
        c2: &HlrBRepCurve,
        d2: &IntRes2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        self.base.perform_with_d2(c1, c2, d2, tol_conf, tol);
    }

    /// `Perform(C1, D1, TolConf, Tol)` (`gxx:44-88`).
    pub fn perform_curve(
        &mut self,
        c1: &HlrBRepCurve,
        d1: &IntRes2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        self.base.perform_curve(c1, d1, tol_conf, tol);
    }

    /// `Perform(C1, TolConf, Tol)` (`gxx:91-117`).
    pub fn perform_auto(&mut self, c1: &HlrBRepCurve, tol_conf: f64, tol: f64) {
        let d1 = self.base.compute_domain(c1, tol_conf.max(tol));
        self.base.perform_curve(c1, &d1, tol_conf, tol);
    }

    /// `IsDone()` (`IntRes2d_Intersection::IsDone`, inherited).
    pub fn is_done(&self) -> bool {
        self.base.is_done()
    }

    /// `NbPoints()` (`IntRes2d_Intersection::NbPoints`, inherited).
    pub fn nb_points(&self) -> usize {
        self.base.nb_points()
    }

    /// `Point(N)` (`IntRes2d_Intersection::Point`, inherited). OCCT answers a
    /// `const&`; the port's `Geom2dIntGInter` answers a copy.
    pub fn point(&self, n: usize) -> IntRes2dIntersectionPoint {
        self.base.point(n)
    }

    /// `NbSegments()` (`IntRes2d_Intersection::NbSegments`, inherited).
    pub fn nb_segments(&self) -> usize {
        self.base.nb_segments()
    }

    /// `Segment(N)` (`IntRes2d_Intersection::Segment`, inherited).
    pub fn segment(&self, n: usize) -> IntRes2dIntersectionSegment {
        self.base.segment(n)
    }
}
