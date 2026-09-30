//! Port of `Geom2dInt_GInter` (`Geom2dInt_GInter_0.cxx:29-89`) =
//! `IntCurve_IntCurveCurveGen` instantiated on `Geom2dInt_Geom2dCurveTool`
//! (`IntCurve_IntCurveCurveGen.gxx` / `.lxx`).
//!
//! Only the `typ1 == Line` arms of `InternalPerform` are ported, because
//! `ShapeFix_ComposeShell::SplitByLine` always passes the cutting line as the
//! first curve. The other arms are marked UNPORTED with their OCCT lines.

use std::f64::consts::PI;

use occt_core::intres2d::{
    IntRes2dDomain, IntRes2dIntersection, IntRes2dIntersectionPoint,
    IntRes2dIntersectionSegment,
};
use occt_core::precision::INFINITE;

use crate::curve::Curve2d;
use super::curve_tool::{self, GeomAbsCurveType};
use super::int_conic_conic::IntCurveIntConicConic;
use super::int_conic_curve::IntConicCurveGen;

/// `Geom2dInt_GInter` (`Geom2dInt_GInter.hxx:43-...`), holding the
/// `IntRes2d_Intersection` result.
#[derive(Clone, Debug, Default)]
pub struct Geom2dIntGInter {
    base: IntRes2dIntersection,
}

impl Geom2dIntGInter {
    pub fn new() -> Self {
        Self { base: IntRes2dIntersection::new() }
    }

    /// `IntRes2d_Intersection::SetReversedParameters` (inherited).
    pub fn set_reversed_parameters(&mut self, flag: bool) {
        self.base.set_reversed_parameters(flag);
    }

    pub fn reversed_parameters(&self) -> bool {
        self.base.reversed_parameters()
    }

    pub fn is_done(&self) -> bool {
        self.base.is_done()
    }

    pub fn result(&self) -> &IntRes2dIntersection {
        &self.base
    }

    pub fn nb_points(&self) -> usize {
        self.base.nb_points()
    }

    pub fn point(&self, n: usize) -> IntRes2dIntersectionPoint {
        *self.base.point(n)
    }

    pub fn nb_segments(&self) -> usize {
        self.base.nb_segments()
    }

    pub fn segment(&self, n: usize) -> IntRes2dIntersectionSegment {
        *self.base.segment(n)
    }

    /// `Perform(C1, D1, C2, D2, TolConf, Tol)` (`gxx:182-225` ->
    /// `InternalPerform`).
    pub fn perform(
        &mut self,
        c1: &dyn Curve2d,
        d1: &IntRes2dDomain,
        c2: &dyn Curve2d,
        d2: &IntRes2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        self.internal_perform(c1, d1, c2, d2, tol_conf, tol);
    }

    /// `Perform(C1, C2, D2, TolConf, Tol)` (`lxx:139-151`), the overload
    /// `SplitByLine` uses: `D1 = ComputeDomain(C1, max(TolConf, Tol))`.
    pub fn perform_with_d2(
        &mut self,
        c1: &dyn Curve2d,
        c2: &dyn Curve2d,
        d2: &IntRes2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        let tol_domain = tol.max(tol_conf);
        let d1 = self.compute_domain(c1, tol_domain);
        self.perform(c1, &d1, c2, d2, tol_conf, tol);
    }

    /// `Perform(C1, C2, TolConf, Tol)` (`lxx:110-121`).
    pub fn perform_curves(
        &mut self,
        c1: &dyn Curve2d,
        c2: &dyn Curve2d,
        tol_conf: f64,
        tol: f64,
    ) {
        let tol_domain = tol.max(tol_conf);
        let d1 = self.compute_domain(c1, tol_domain);
        let d2 = self.compute_domain(c2, tol_domain);
        self.perform(c1, &d1, c2, &d2, tol_conf, tol);
    }

    /// `ComputeDomain(C1, TolDomain)` (`gxx:120-177`).
    pub fn compute_domain(&self, c: &dyn Curve2d, tol_domain: f64) -> IntRes2dDomain {
        match curve_tool::get_type(c) {
            GeomAbsCurveType::Ellipse | GeomAbsCurveType::Circle => {
                let first = curve_tool::first_parameter(c);
                let last = curve_tool::last_parameter(c);
                let p1 = curve_tool::value(c, first);
                let p2 = curve_tool::value(c, last);
                let mut d = IntRes2dDomain::bounded(&p1, first, tol_domain, &p2, last, tol_domain);
                d.set_equivalent_parameters(first, first + PI + PI);
                d
            }
            _ => {
                let param_inf = curve_tool::first_parameter(c);
                let param_sup = curve_tool::last_parameter(c);
                let mut d = IntRes2dDomain::new();
                if param_inf > -INFINITE {
                    if param_sup < INFINITE {
                        d.set_bounded(
                            &curve_tool::value(c, param_inf),
                            param_inf,
                            tol_domain,
                            &curve_tool::value(c, param_sup),
                            param_sup,
                            tol_domain,
                        );
                    } else {
                        d.set_semi_infinite(&curve_tool::value(c, param_inf), param_inf, tol_domain, true);
                    }
                } else if param_sup < INFINITE {
                    d.set_semi_infinite(&curve_tool::value(c, param_sup), param_sup, tol_domain, false);
                }
                d
            }
        }
    }

    /// `IntCurve_IntCurveCurveGen::InternalPerform` (`gxx:235-...`).
    fn internal_perform(
        &mut self,
        c1: &dyn Curve2d,
        d1: &IntRes2dDomain,
        c2: &dyn Curve2d,
        d2: &IntRes2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        let typ1 = curve_tool::get_type(c1);
        let typ2 = curve_tool::get_type(c2);
        if typ1 == GeomAbsCurveType::Line {
            match typ2 {
                GeomAbsCurveType::Line => {
                    // `gxx:249-264`: Line/Line -> IntConicConic.
                    let l1 = curve_tool::line(c1).expect("line");
                    let l2 = curve_tool::line(c2).expect("line");
                    let mut icc = IntCurveIntConicConic::new();
                    icc.set_reversed_parameters(false);
                    icc.perform_line_line(&l1, d1, &l2, d2, tol_conf, tol);
                    self.base.set_values(icc.result());
                }
                // T-100: OCCT routes Line/Circle, Line/Ellipse and the
                // parabolic/hyperbolic pairs to the dedicated `IntConicConic`
                // overloads (`IntCurve_IntConicConic_1.cxx:2236`, `:2861`, ...),
                // of which this port only has Line/Line. Those overloads are a
                // **specialisation**, not new capability: the generic arm below
                // builds `IntCurveIConicTool::from_lin2d` and calls
                // `IntConicCurveGen::perform` (`gxx:245-779`), whose
                // `MyImpParTool` takes *any* `Curve2d`, and
                // `IntCurveIConicTool` already represents every conic
                // (`from_circ2d` / `from_elips2d` / `from_parab2d` /
                // `from_hypr2d`). So a Line/Conic pair is handled correctly here;
                // the dedicated overloads remain unported as an optimisation
                // only. See specs/_a3n00_gap_analysis.md §9.305.
                _ => {
                    // `gxx:327-339`: Line/default -> IntConicCurve.
                    let l1 = curve_tool::line(c1).expect("line");
                    let mut cc = IntConicCurveGen::new();
                    cc.set_reversed_parameters(false);
                    cc.perform_line(&l1, d1, c2, d2, tol_conf, tol);
                    self.base.set_values(cc.result());
                }
            }
        } else {
            // UNPORTED (IntCurve_IntCurveCurveGen.gxx:339-...): the typ1 != Line
            // arms (conic/curve and curve/curve). **Unreachable from this
            // codebase**: the only caller is
            // `shape_fix_compose_shell/split_by_line.rs:203-204`, whose first
            // curve is always `Geom2dLine` (`:42`), so `typ1` is always `Line`.
            // Left unimplemented rather than guessed; see §9.305.
        }
    }
}
