//! The curve/point, point/curve and point/point overloads of
//! `Bisector_BisecAna::Perform` (`Bisector_BisecAna.cxx:966-1305`).

use std::sync::Arc;

use occt_core::elib::clib2d;
use occt_core::gp::{GpAx2d, GpCirc2d, GpLin2d, GpPnt2d, GpVec2d};
use occt_core::precision::INFINITE;

use crate::circle::Geom2dCircle;
use crate::curve::Curve2d;
use crate::ellipse::Geom2dEllipse;
use crate::gcc_ana_bisec::{GccAnaCircPnt2dBisec, GccAnaLinPnt2dBisec, GccAnaPnt2dBisec};
use crate::gcc_int::GccIntBisec;
use crate::hyperbola::Geom2dHyperbola;
use crate::line::Geom2dLine;
use crate::parabola::Geom2dParabola;
use crate::trimmed::Geom2dTrimmedCurve;

use super::degenerate::degenerate;
use super::BisectorBisecAna;

impl BisectorBisecAna {
    /// `Bisector_BisecAna::Perform(const occ::handle<Geom2d_Curve>&,
    /// const occ::handle<Geom2d_Point>&, ...)` (`cxx:977-1224`).
    ///
    /// Curve first, point second. The `switch (cas) default` arm only prints
    /// "Not yet implemented" (`cxx:1218-1221`) and leaves the bisector unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn perform_curve_point(
        &mut self,
        first_curve: &dyn Curve2d,
        second_point: &GpPnt2d,
        apoint: &GpPnt2d,
        first_vector: &GpVec2d,
        second_vector: &GpVec2d,
        adirection: f64,
        tolerance: f64,
        oncurve: bool,
    ) {
        let mut firstparameter = 0.0;
        // `gp_Vec2d VecRef(0., 0.);` (`cxx:989`): the reference vector is never
        // set in this overload.
        let vec_ref = GpVec2d::zero();

        let curve: &dyn Curve2d = first_curve.trimmed_basis().unwrap_or(first_curve);

        let mut circle = GpCirc2d::default();
        let mut line = GpLin2d::new(GpAx2d::new(GpPnt2d::zero(), Default::default()));
        let cas: u8;
        if curve.is_geom2d_circle() {
            cas = 1;
            circle = curve.gp_circ2d().unwrap_or_default();
        } else if curve.is_geom2d_line() {
            cas = 2;
            line = curve.gp_lin2d().unwrap_or(line);
        } else {
            cas = 0;
        }

        match cas {
            1 => {
                // `cxx:1035-1163`: bisector point - circle.
                let bisector =
                    GccAnaCircPnt2dBisec::with_tolerance(&circle, second_point, tolerance);
                let mut distancemini = INFINITE;
                let mut the_sol: Option<GccIntBisec> = None;
                let mut thesense = false;
                if bisector.is_done() {
                    let nbsolution = bisector.nb_solutions();
                    for i in 1..=nbsolution {
                        let mut solution = match bisector.this_solution(i) {
                            Ok(s) => s,
                            Err(_) => continue,
                        };
                        degenerate(&mut solution, tolerance);
                        let r = Self::distance(
                            apoint,
                            &solution,
                            first_vector,
                            second_vector,
                            &vec_ref,
                            adirection,
                            false,
                            false,
                        );
                        if r.distance <= distancemini {
                            the_sol = Some(solution);
                            firstparameter = r.parameter;
                            thesense = r.sense;
                            distancemini = r.distance;
                        }
                    }
                    if let Some(the_sol) = the_sol {
                        match the_sol {
                            GccIntBisec::Line(l) => {
                                let circlecenter = circle.location();
                                let mut secondparameter =
                                    clib2d::parameter_lin2d(&l, &circlecenter);
                                let bisectorcurve: Arc<dyn Curve2d> = Arc::from(
                                    Geom2dLine::new(GpAx2d::new(l.location(), *l.direction()))
                                        .clone_dyn(),
                                );
                                if !thesense {
                                    if secondparameter > firstparameter {
                                        secondparameter = -INFINITE;
                                    } else {
                                        secondparameter -= 1.0e-8;
                                    }
                                } else if secondparameter < firstparameter {
                                    secondparameter = INFINITE;
                                } else {
                                    secondparameter += 1.0e-8;
                                }
                                self.thebisector = Some(Geom2dTrimmedCurve::new(
                                    bisectorcurve,
                                    firstparameter,
                                    secondparameter,
                                ));
                            }
                            GccIntBisec::Circle(c) => {
                                let bisectorcurve: Arc<dyn Curve2d> =
                                    Arc::from(Geom2dCircle::new(c).clone_dyn());
                                self.thebisector = Some(if !thesense {
                                    Geom2dTrimmedCurve::new_sense(
                                        bisectorcurve,
                                        firstparameter - 2.0 * std::f64::consts::PI,
                                        firstparameter,
                                        thesense,
                                        true,
                                    )
                                } else {
                                    Geom2dTrimmedCurve::new_sense(
                                        bisectorcurve,
                                        firstparameter,
                                        firstparameter + 2.0 * std::f64::consts::PI,
                                        thesense,
                                        true,
                                    )
                                });
                            }
                            GccIntBisec::Hyperbola(h) => {
                                let bisectorcurve: Arc<dyn Curve2d> =
                                    Arc::from(Geom2dHyperbola::new(h).clone_dyn());
                                self.thebisector = Some(if !thesense {
                                    Geom2dTrimmedCurve::new(
                                        bisectorcurve,
                                        firstparameter,
                                        -INFINITE,
                                    )
                                } else {
                                    Geom2dTrimmedCurve::new(
                                        bisectorcurve,
                                        firstparameter,
                                        INFINITE,
                                    )
                                });
                            }
                            GccIntBisec::Ellipse(e) => {
                                let bisectorcurve: Arc<dyn Curve2d> =
                                    Arc::from(Geom2dEllipse::new(e).clone_dyn());
                                self.thebisector = Some(if !thesense {
                                    Geom2dTrimmedCurve::new_sense(
                                        bisectorcurve,
                                        firstparameter - 2.0 * std::f64::consts::PI,
                                        firstparameter,
                                        thesense,
                                        true,
                                    )
                                } else {
                                    Geom2dTrimmedCurve::new_sense(
                                        bisectorcurve,
                                        firstparameter,
                                        firstparameter + 2.0 * std::f64::consts::PI,
                                        thesense,
                                        true,
                                    )
                                });
                            }
                            GccIntBisec::Parabola(_) | GccIntBisec::Point(_) => {}
                        }
                    }
                }
            }
            2 => {
                // `cxx:1166-1216`: bisector point - straight line.
                let bisector = GccAnaLinPnt2dBisec::new(&line, second_point);
                let mut solution = match bisector.this_solution() {
                    Ok(s) => s,
                    Err(_) => return,
                };
                degenerate(&mut solution, tolerance);
                // `if (type == GccInt_Lin) ... else if (type == GccInt_Par) ...`
                // - those are the only two arcs `GccAnaLinPnt2dBisec` produces.
                let bisectorcurve: Arc<dyn Curve2d> = match &solution {
                    GccIntBisec::Line(l) => Arc::from(
                        Geom2dLine::new(GpAx2d::new(l.location(), *l.direction())).clone_dyn(),
                    ),
                    GccIntBisec::Parabola(pb) => Arc::from(Geom2dParabola::new(*pb).clone_dyn()),
                    _ => return,
                };
                let mut thesense = false;
                let r = Self::distance(
                    apoint,
                    &solution,
                    first_vector,
                    second_vector,
                    &vec_ref,
                    adirection,
                    thesense,
                    false,
                );
                if r.status || !oncurve {
                    firstparameter = r.parameter;
                    thesense = r.sense;
                }
                self.thebisector = Some(if !thesense {
                    Geom2dTrimmedCurve::new(bisectorcurve, firstparameter, -INFINITE)
                } else {
                    Geom2dTrimmedCurve::new(bisectorcurve, firstparameter, INFINITE)
                });
            }
            _ => {}
        }
    }

    /// `Bisector_BisecAna::Perform(const occ::handle<Geom2d_Point>&,
    /// const occ::handle<Geom2d_Curve>&, ...)` (`cxx:1237-1259`): the point/curve
    /// overload simply calls the curve/point one with the two sides swapped and
    /// `-adirection`, dropping `tolerance` (`0.` is passed at `cxx:1256`).
    #[allow(clippy::too_many_arguments)]
    pub fn perform_point_curve(
        &mut self,
        first_point: &GpPnt2d,
        second_curve: &dyn Curve2d,
        apoint: &GpPnt2d,
        first_vector: &GpVec2d,
        second_vector: &GpVec2d,
        adirection: f64,
        _tolerance: f64,
        oncurve: bool,
    ) {
        let adirectionreverse = -adirection;
        self.perform_curve_point(
            second_curve,
            first_point,
            apoint,
            second_vector,
            first_vector,
            adirectionreverse,
            0.0,
            oncurve,
        );
    }

    /// `Bisector_BisecAna::Perform(const occ::handle<Geom2d_Point>&,
    /// const occ::handle<Geom2d_Point>&, ...)` (`cxx:1273-1305`).
    #[allow(clippy::too_many_arguments)]
    pub fn perform_point_point(
        &mut self,
        first_point: &GpPnt2d,
        second_point: &GpPnt2d,
        apoint: &GpPnt2d,
        first_vector: &GpVec2d,
        second_vector: &GpVec2d,
        adirection: f64,
        _tolerance: f64,
        oncurve: bool,
    ) {
        let vec_ref = GpVec2d::zero();
        let bisector = GccAnaPnt2dBisec::new(first_point, second_point);
        let line = match bisector.this_solution() {
            Ok(l) => l,
            Err(_) => return,
        };
        let solution = GccIntBisec::Line(line);
        let r = Self::distance(
            apoint,
            &solution,
            first_vector,
            second_vector,
            &vec_ref,
            adirection,
            false,
            false,
        );
        if r.status || !oncurve {
            let bisectorcurve: Arc<dyn Curve2d> =
                Arc::from(Geom2dLine::new(GpAx2d::new(line.location(), *line.direction())).clone_dyn());
            self.thebisector = Some(if !r.sense {
                Geom2dTrimmedCurve::new(bisectorcurve, r.parameter, -INFINITE)
            } else {
                Geom2dTrimmedCurve::new(bisectorcurve, r.parameter, INFINITE)
            });
        }
    }
}
