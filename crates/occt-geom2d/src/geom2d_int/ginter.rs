//! Port of `Geom2dInt_GInter` (`Geom2dInt_GInter_0.cxx:29-89`) =
//! `IntCurve_IntCurveCurveGen` instantiated on `Geom2dInt_Geom2dCurveTool`
//! (`IntCurve_IntCurveCurveGen.gxx` / `.lxx`).
//!
//! `typ1 == Line` row of `InternalPerform` (`gxx:249-339`) is fully ported
//! (Line/Line, Line/Circle and Line/Ellipse through the dedicated
//! `IntCurve_IntConicConic` overloads of `IntCurve_IntConicConic_1.cxx`, the
//! other conic pairs through the `IntCurve_IntConicConic.cxx` overloads, and
//! everything else through `IntConicCurveGen`), and so are the `Circle`
//! (`gxx:342-433`), `Ellipse` (`gxx:438-529`), `Parabola` (`gxx:532-624`) and
//! `Hyperbola` (`gxx:627-723`) rows, plus the conic-first arms of the `default:`
//! row (`gxx:724-...`).
//!
//! The two `default: typ2` arms of the `default: typ1` row (`gxx:799-812`)
//! call `intcurvcurv.Perform(C1, D1, C2, D2, TolConf, Tol)` - the
//! `IntCurve_IntPolyPolyGen` polygon-approximation intersector
//! (`Geom2dInt_TheIntPCurvePCurveOfGInter` + `Intf_InterferencePolygon2d`),
//! ported as `Geom2dIntIntPolyPolyGen` and wired into `perform_curve` below.
//! See specs/_a3n00_gap_analysis.md 9.305 / 9.628 / 9.629.

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
use super::int_poly_poly_gen::Geom2dIntIntPolyPolyGen;

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

    /// `Perform(C, D, TolConf, Tol)` (`gxx:91-117`), the single-curve overload
    /// `ShapeAnalysis_Wire::CheckSelfIntersectingEdge`
    /// (`ShapeAnalysis_Wire.cxx:1300-1301`) uses.
    ///
    /// `gxx:99-106`: when the adaptor type is `Ellipse` / `Circle` /
    /// `Parabola` / `Hyperbola` / `Line` OCCT **resets the fields and marks the
    /// result done without intersecting anything** - a conic is never tested
    /// against itself. The `default:` arm goes to
    /// `intcurvcurv.Perform(C, D, TolConf, Tol)` (`gxx:108-113`), the
    /// polygon-approximation intersector `IntCurve_IntPolyPolyGen`
    /// (`IntCurve_IntPolyPolyGen.gxx`) instantiated as
    /// `Geom2dInt_TheIntPCurvePCurveOfGInter` (`Geom2dInt_TheIntPCurvePCurveOfGInter_0.cxx`
    /// + `Intf_InterferencePolygon2d.cxx`), ported as `Geom2dIntIntPolyPolyGen`.
    pub fn perform_curve(
        &mut self,
        c: &dyn Curve2d,
        d: &IntRes2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        match curve_tool::get_type(c) {
            GeomAbsCurveType::Ellipse
            | GeomAbsCurveType::Circle
            | GeomAbsCurveType::Parabola
            | GeomAbsCurveType::Hyperbola
            | GeomAbsCurveType::Line => {
                // `gxx:97-106`.
                self.base.reset_fields();
                self.base.done = true;
            }
            _ => {
                // `gxx:108-113`.
                let mut intcurvcurv = Geom2dIntIntPolyPolyGen::new();
                intcurvcurv.base.set_reversed_parameters(false);
                intcurvcurv.perform_self(c, d, tol_conf, tol);
                self.base.set_values(intcurvcurv.result());
                self.base.done = true;
            }
        }
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

        use GeomAbsCurveType::*;

        match typ1 {
            // `gxx:249-339`: the `typ1 == Line` row.
            Line => match typ2 {
                Line => {
                    // `gxx:249-264`.
                    let l1 = curve_tool::line(c1).expect("line");
                    let l2 = curve_tool::line(c2).expect("line");
                    self.run_conic_conic(false, |icc| {
                        icc.perform_line_line(&l1, d1, &l2, d2, tol_conf, tol)
                    });
                }
                Circle => {
                    // `gxx:266-277`.
                    let l1 = curve_tool::line(c1).expect("line");
                    let c2c = curve_tool::circle(c2).expect("circle");
                    self.run_conic_conic(false, |icc| {
                        icc.perform_line_circle(&l1, d1, &c2c, d2, tol_conf, tol)
                    });
                }
                Ellipse => {
                    // `gxx:279-293`.
                    let l1 = curve_tool::line(c1).expect("line");
                    let e2 = curve_tool::ellipse(c2).expect("ellipse");
                    self.run_conic_conic(false, |icc| {
                        icc.perform_line_ellipse(&l1, d1, &e2, d2, tol_conf, tol)
                    });
                }
                Parabola => {
                    // `gxx:295-309`: `Perform(Lin2d, ..., Parab2d, ...)`
                    // (`IntCurve_IntConicConic.cxx:109-225`).
                    let l1 = curve_tool::line(c1).expect("line");
                    let p2 = curve_tool::parabola(c2).expect("parabola");
                    self.run_conic_conic(false, |icc| {
                        icc.perform_line_parabola(&l1, d1, &p2, d2, tol_conf, tol)
                    });
                }
                Hyperbola => {
                    // `gxx:311-325`: `Perform(Lin2d, ..., Hypr2d, ...)`
                    // (`IntCurve_IntConicConic.cxx:229-333`).
                    let l1 = curve_tool::line(c1).expect("line");
                    let h2 = curve_tool::hyperbola(c2).expect("hyperbola");
                    self.run_conic_conic(false, |icc| {
                        icc.perform_line_hyperbola(&l1, d1, &h2, d2, tol_conf, tol)
                    });
                }
                _ => {
                    // `gxx:327-339`: Line/default -> IntConicCurve.
                    let l1 = curve_tool::line(c1).expect("line");
                    self.run_conic_curve(false, |cc| cc.perform_line(&l1, d1, c2, d2, tol_conf, tol));
                }
            },

            // `gxx:342-433`: the `typ1 == Circle` row.
            Circle => {
                let c1c = curve_tool::circle(c1).expect("circle");
                match typ2 {
                    Line => {
                        // `gxx:346-358`: reversed into
                        // `Perform(Line(C2), D2, Circ(C1), D1, ...)`.
                        let l2 = curve_tool::line(c2).expect("line");
                        self.run_conic_conic(true, |icc| {
                            icc.perform_line_circle(&l2, d2, &c1c, d1, tol_conf, tol)
                        });
                    }
                    Circle => {
                        // `gxx:360-374`.
                        let c2c = curve_tool::circle(c2).expect("circle");
                        self.run_conic_conic(false, |icc| {
                            icc.perform_circle_circle(&c1c, d1, &c2c, d2, tol_conf, tol)
                        });
                    }
                    Ellipse => {
                        // `gxx:377-388`: `Perform(Circ2d, ..., Elips2d, ...)`
                        // (`IntCurve_IntConicConic.cxx:438-483`).
                        let e2 = curve_tool::ellipse(c2).expect("ellipse");
                        self.run_conic_conic(false, |icc| {
                            icc.perform_circle_ellipse(&c1c, d1, &e2, d2, tol_conf, tol)
                        });
                    }
                    Parabola => {
                        // `gxx:390-404`: `Perform(Circ2d, ..., Parab2d, ...)`
                        // (`IntCurve_IntConicConic.cxx:337-435`).
                        let p2 = curve_tool::parabola(c2).expect("parabola");
                        self.run_conic_conic(false, |icc| {
                            icc.perform_circle_parabola(&c1c, d1, &p2, d2, tol_conf, tol)
                        });
                    }
                    Hyperbola => {
                        // `gxx:406-419`: `Perform(Circ2d, ..., Hypr2d, ...)`
                        // (`IntCurve_IntConicConic.cxx:486-581`).
                        let h2 = curve_tool::hyperbola(c2).expect("hyperbola");
                        self.run_conic_conic(false, |icc| {
                            icc.perform_circle_hyperbola(&c1c, d1, &h2, d2, tol_conf, tol)
                        });
                    }
                    _ => {
                        // `gxx:421-433`: Circle/default -> IntConicCurve.
                        self.run_conic_curve(false, |cc| {
                            cc.perform_circle(&c1c, d1, c2, d2, tol_conf, tol)
                        });
                    }
                }
            }

            // `gxx:438-529`: the `typ1 == Ellipse` row.
            Ellipse => {
                let e1 = curve_tool::ellipse(c1).expect("ellipse");
                match typ2 {
                    Line => {
                        // `gxx:441-455`: reversed into
                        // `Perform(Line(C2), D2, Elips(C1), D1, ...)`.
                        let l2 = curve_tool::line(c2).expect("line");
                        self.run_conic_conic(true, |icc| {
                            icc.perform_line_ellipse(&l2, d2, &e1, d1, tol_conf, tol)
                        });
                    }
                    Circle => {
                        // `gxx:458-471`: reversed into
                        // `Perform(Circ(C2), D2, Elips(C1), D1, ...)`.
                        let c2c = curve_tool::circle(c2).expect("circle");
                        self.run_conic_conic(true, |icc| {
                            icc.perform_circle_ellipse(&c2c, d2, &e1, d1, tol_conf, tol)
                        });
                    }
                    Ellipse => {
                        // `gxx:473-486`: `Perform(Elips2d, ..., Elips2d, ...)`
                        // (`IntCurve_IntConicConic.cxx:913-957`).
                        let e2 = curve_tool::ellipse(c2).expect("ellipse");
                        self.run_conic_conic(false, |icc| {
                            icc.perform_ellipse_ellipse(&e1, d1, &e2, d2, tol_conf, tol)
                        });
                    }
                    Parabola => {
                        // `gxx:488-502`: `Perform(Elips2d, ..., Parab2d, ...)`
                        // (`IntCurve_IntConicConic.cxx:690-807`).
                        let p2 = curve_tool::parabola(c2).expect("parabola");
                        self.run_conic_conic(false, |icc| {
                            icc.perform_ellipse_parabola(&e1, d1, &p2, d2, tol_conf, tol)
                        });
                    }
                    Hyperbola => {
                        // `gxx:504-518`: `Perform(Elips2d, ..., Hypr2d, ...)`
                        // (`IntCurve_IntConicConic.cxx:959-1062`).
                        let h2 = curve_tool::hyperbola(c2).expect("hyperbola");
                        self.run_conic_conic(false, |icc| {
                            icc.perform_ellipse_hyperbola(&e1, d1, &h2, d2, tol_conf, tol)
                        });
                    }
                    _ => {
                        // `gxx:517-529`: Ellipse/default -> IntConicCurve.
                        self.run_conic_curve(false, |cc| {
                            cc.perform_ellipse(&e1, d1, c2, d2, tol_conf, tol)
                        });
                    }
                }
            }

            // `gxx:532-624`: the `typ1 == Parabola` row.
            Parabola => {
                let p1 = curve_tool::parabola(c1).expect("parabola");
                match typ2 {
                    Line => {
                        // `gxx:537-549`: reversed into
                        // `Perform(Line(C2), D2, Parab(C1), D1, ...)`.
                        let l2 = curve_tool::line(c2).expect("line");
                        self.run_conic_conic(true, |icc| {
                            icc.perform_line_parabola(&l2, d2, &p1, d1, tol_conf, tol)
                        });
                    }
                    Circle => {
                        // `gxx:552-564`: reversed into
                        // `Perform(Circ(C2), D2, Parab(C1), D1, ...)`.
                        let c2c = curve_tool::circle(c2).expect("circle");
                        self.run_conic_conic(true, |icc| {
                            icc.perform_circle_parabola(&c2c, d2, &p1, d1, tol_conf, tol)
                        });
                    }
                    Ellipse => {
                        // `gxx:567-579`: reversed into
                        // `Perform(Elips(C2), D2, Parab(C1), D1, ...)`.
                        let e2 = curve_tool::ellipse(c2).expect("ellipse");
                        self.run_conic_conic(true, |icc| {
                            icc.perform_ellipse_parabola(&e2, d2, &p1, d1, tol_conf, tol)
                        });
                    }
                    Parabola => {
                        // `gxx:582-594`: `Perform(Parab2d, ..., Parab2d, ...)`
                        // (`IntCurve_IntConicConic.cxx:584-688`).
                        let p2 = curve_tool::parabola(c2).expect("parabola");
                        self.run_conic_conic(false, |icc| {
                            icc.perform_parabola_parabola(&p1, d1, &p2, d2, tol_conf, tol)
                        });
                    }
                    Hyperbola => {
                        // `gxx:597-609`: `Perform(Parab2d, ..., Hypr2d, ...)`
                        // (`IntCurve_IntConicConic.cxx:809-911`).
                        let h2 = curve_tool::hyperbola(c2).expect("hyperbola");
                        self.run_conic_conic(false, |icc| {
                            icc.perform_parabola_hyperbola(&p1, d1, &h2, d2, tol_conf, tol)
                        });
                    }
                    _ => {
                        // `gxx:612-624`: Parabola/default -> IntConicCurve.
                        self.run_conic_curve(false, |cc| {
                            cc.perform_parabola(&p1, d1, c2, d2, tol_conf, tol)
                        });
                    }
                }
            }

            // `gxx:627-723`: the `typ1 == Hyperbola` row.
            Hyperbola => {
                let h1 = curve_tool::hyperbola(c1).expect("hyperbola");
                match typ2 {
                    Line => {
                        // `gxx:631-643`: reversed into
                        // `Perform(Line(C2), D2, Hypr(C1), D1, ...)`.
                        let l2 = curve_tool::line(c2).expect("line");
                        self.run_conic_conic(true, |icc| {
                            icc.perform_line_hyperbola(&l2, d2, &h1, d1, tol_conf, tol)
                        });
                    }
                    Circle => {
                        // `gxx:646-658`: reversed into
                        // `Perform(Circ(C2), D2, Hypr(C1), D1, ...)`.
                        let c2c = curve_tool::circle(c2).expect("circle");
                        self.run_conic_conic(true, |icc| {
                            icc.perform_circle_hyperbola(&c2c, d2, &h1, d1, tol_conf, tol)
                        });
                    }
                    Ellipse => {
                        // `gxx:661-673`: reversed into
                        // `Perform(Elips(C2), D2, Hypr(C1), D1, ...)`.
                        let e2 = curve_tool::ellipse(c2).expect("ellipse");
                        self.run_conic_conic(true, |icc| {
                            icc.perform_ellipse_hyperbola(&e2, d2, &h1, d1, tol_conf, tol)
                        });
                    }
                    Parabola => {
                        // `gxx:676-688`: reversed into
                        // `Perform(Parab(C2), D2, Hypr(C1), D1, ...)`.
                        let p2 = curve_tool::parabola(c2).expect("parabola");
                        self.run_conic_conic(true, |icc| {
                            icc.perform_parabola_hyperbola(&p2, d2, &h1, d1, tol_conf, tol)
                        });
                    }
                    Hyperbola => {
                        // `gxx:691-707`: `Perform(Hypr2d, ..., Hypr2d, ...)`
                        // (`IntCurve_IntConicConic.cxx:1064-1166`).
                        let h2 = curve_tool::hyperbola(c2).expect("hyperbola");
                        self.run_conic_conic(false, |icc| {
                            icc.perform_hyperbola_hyperbola(&h1, d1, &h2, d2, tol_conf, tol)
                        });
                    }
                    _ => {
                        // `gxx:710-722`: Hyperbola/default -> IntConicCurve.
                        self.run_conic_curve(false, |cc| {
                            cc.perform_hyperbola(&h1, d1, c2, d2, tol_conf, tol)
                        });
                    }
                }
            }

            // `gxx:726-813`: the `default:` row (`typ1` is a Bezier, BSpline,
            // offset or other curve).
            _ => match typ2 {
                Line => {
                    // `gxx:729-741`.
                    let l2 = curve_tool::line(c2).expect("line");
                    self.run_conic_curve(true, |cc| {
                        cc.perform_line(&l2, d2, c1, d1, tol_conf, tol)
                    });
                }
                Circle => {
                    // `gxx:744-756`.
                    let c2c = curve_tool::circle(c2).expect("circle");
                    self.run_conic_curve(true, |cc| {
                        cc.perform_circle(&c2c, d2, c1, d1, tol_conf, tol)
                    });
                }
                Ellipse => {
                    // `gxx:759-771`.
                    let e2 = curve_tool::ellipse(c2).expect("ellipse");
                    self.run_conic_curve(true, |cc| {
                        cc.perform_ellipse(&e2, d2, c1, d1, tol_conf, tol)
                    });
                }
                Parabola => {
                    // `gxx:774-786`.
                    let p2 = curve_tool::parabola(c2).expect("parabola");
                    self.run_conic_curve(true, |cc| {
                        cc.perform_parabola(&p2, d2, c1, d1, tol_conf, tol)
                    });
                }
                Hyperbola => {
                    // `gxx:788-798`.
                    let h2 = curve_tool::hyperbola(c2).expect("hyperbola");
                    self.run_conic_curve(true, |cc| {
                        cc.perform_hyperbola(&h2, d2, c1, d1, tol_conf, tol)
                    });
                }
                _ => {
                    // `gxx:799-812`: `intcurvcurv.SetReversedParameters(false);
                    // intcurvcurv.Perform(C1, D1, C2, D2, TolConf, Tol)`, then
                    // `SetValues` (this instantiation has `Composite == false`),
                    // `done = true`. This is the `IntCurve_IntPolyPolyGen`
                    // polygon-approximation intersector ported as
                    // `Geom2dIntIntPolyPolyGen`.
                    let mut intcurvcurv = Geom2dIntIntPolyPolyGen::new();
                    intcurvcurv.base.set_reversed_parameters(false);
                    intcurvcurv.perform(c1, d1, c2, d2, tol_conf, tol);
                    self.base.set_values(intcurvcurv.result());
                    self.base.done = true;
                }
            },
        }
    }

    /// One `intconiconi` arm of `IntCurve_IntCurveCurveGen::InternalPerform`
    /// (`gxx:253-260` and its twins): build the dedicated
    /// `IntCurve_IntConicConic`, run the arm's overload on it, and copy the
    /// result into the `IntRes2d_Intersection` base (`Composite` is always false
    /// for this instantiation, so `SetValues` is the only tail).
    fn run_conic_conic<F>(&mut self, reversed: bool, f: F)
    where
        F: FnOnce(&mut IntCurveIntConicConic),
    {
        let mut icc = IntCurveIntConicConic::new();
        icc.set_reversed_parameters(reversed);
        f(&mut icc);
        self.base.set_values(icc.result());
    }

    /// One `intconicurv` arm of `IntCurve_IntCurveCurveGen::InternalPerform`
    /// (`gxx:328-336` and its twins).
    fn run_conic_curve<F>(&mut self, reversed: bool, f: F)
    where
        F: FnOnce(&mut IntConicCurveGen),
    {
        let mut cc = IntConicCurveGen::new();
        cc.set_reversed_parameters(reversed);
        f(&mut cc);
        self.base.set_values(cc.result());
    }
}
