//! `Bisector_BisecAna::Perform(Cu1, Cu2, ...)` (`Bisector_BisecAna.cxx:232-963`):
//! the bisector between two curves, each of which must reduce to a circle or a
//! straight line.

use std::sync::Arc;

use occt_core::elib::clib2d;
use occt_core::gp::{GpAx2d, GpCirc2d, GpDir2d, GpLin2d, GpPnt2d, GpVec2d};
use occt_core::kernel::geomabs::JoinType;
use occt_core::precision::{CONFUSION, INFINITE};

use crate::circle::Geom2dCircle;
use crate::curve::Curve2d;
use crate::ellipse::Geom2dEllipse;
use crate::gcc_ana_bisec::{GccAnaCirc2dBisec, GccAnaCircLin2dBisec};
use crate::gcc_int::GccIntBisec;
use crate::hyperbola::Geom2dHyperbola;
use crate::line::Geom2dLine;
use crate::parabola::Geom2dParabola;
use crate::trimmed::Geom2dTrimmedCurve;

use super::degenerate::degenerate;
use super::{BisectorBisecAna, BisectorBisecAnaError};

impl BisectorBisecAna {
    /// `Bisector_BisecAna::Perform(const occ::handle<Geom2d_Curve>&,
    /// const occ::handle<Geom2d_Curve>&, ...)` (`cxx:232-963`).
    ///
    /// The nature of both arguments is read from the **basis** curve, so a
    /// `Geom2d_TrimmedCurve` is unwrapped first (`cxx:272-288`). Any other
    /// combination prints "Not yet implemented" and leaves `cas == 0`, which the
    /// `switch` answers with `StdFail_NotDone` (`cxx:957-961`).
    #[allow(clippy::too_many_arguments)]
    pub fn perform_curve_curve(
        &mut self,
        first_curve: &dyn Curve2d,
        second_curve: &dyn Curve2d,
        apoint: &GpPnt2d,
        first_vector: &GpVec2d,
        second_vector: &GpVec2d,
        adirection: f64,
        ajointype: JoinType,
        tolerance: f64,
        oncurve: bool,
    ) -> Result<(), BisectorBisecAnaError> {
        thread_local! {
            static DEPTH: std::cell::Cell<u32> = const { std::cell::Cell::new(0) };
        }
        struct Guard;
        impl Drop for Guard {
            fn drop(&mut self) {
                DEPTH.with(|d| d.set(d.get() - 1));
            }
        }
        DEPTH.with(|d| d.set(d.get() + 1));
        let _guard = Guard;
        if DEPTH.with(|d| d.get()) > 40 {
            panic!("perform_curve_curve recursion depth");
        }
        // jgv: for OCC26296 (`cxx:258-267`).
        let mut line_bis_vec = GpVec2d::zero();
        let (pnt1, mut tan1) = first_curve.d1(first_curve.last_parameter());
        let (pnt2, tan2) = second_curve.d1(second_curve.first_parameter());
        if !oncurve {
            line_bis_vec = GpVec2d::new(pnt2.x() - pnt1.x(), pnt2.y() - pnt1.y());
            line_bis_vec.rotate(std::f64::consts::FRAC_PI_2);
        }
        tan1.reverse();

        // `cxx:272-288`: a trimmed argument contributes its basis curve.
        let curve_f: &dyn Curve2d = first_curve.trimmed_basis().unwrap_or(first_curve);
        let curve_e: &dyn Curve2d = second_curve.trimmed_basis().unwrap_or(second_curve);

        // `cxx:296-350`: nature of the arguments.
        let mut circle1 = GpCirc2d::default();
        let mut circle2 = GpCirc2d::default();
        let mut line1 = GpLin2d::new(GpAx2d::new(GpPnt2d::zero(), GpDir2d::default()));
        let mut line2 = line1;
        let cas: u8;
        if curve_f.is_geom2d_circle() {
            if curve_e.is_geom2d_circle() {
                cas = 1;
                circle1 = curve_f.gp_circ2d().unwrap_or_default();
                circle2 = curve_e.gp_circ2d().unwrap_or_default();
            } else if curve_e.is_geom2d_line() {
                cas = 2;
                circle1 = curve_f.gp_circ2d().unwrap_or_default();
                line2 = curve_e.gp_lin2d().unwrap_or(line2);
            } else {
                cas = 0;
            }
        } else if curve_f.is_geom2d_line() {
            if curve_e.is_geom2d_circle() {
                cas = 2;
                circle1 = curve_e.gp_circ2d().unwrap_or_default();
                line2 = curve_f.gp_lin2d().unwrap_or(line2);
            } else if curve_e.is_geom2d_line() {
                cas = 3;
                line1 = curve_f.gp_lin2d().unwrap_or(line1);
                line2 = curve_e.gp_lin2d().unwrap_or(line2);
            } else {
                cas = 0;
            }
        } else {
            cas = 0;
        }

        match cas {
            1 => self.perform_circle_circle(                first_curve,
                second_curve,
                apoint,
                first_vector,
                second_vector,
                adirection,
                tolerance,
                oncurve,
                &tan1,
                &tan2,
                &line_bis_vec,
                circle1,
                circle2,
            )?,
            2 => {
                self.perform_circle_line(
                    apoint,
                    first_vector,
                    second_vector,
                    adirection,
                    ajointype,
                    tolerance,
                    oncurve,
                    &tan1,
                    &tan2,
                    &line_bis_vec,
                    circle1,
                    line2,
                )?;
            }
            3 => {
                self.perform_line_line(
                    apoint,
                    first_vector,
                    second_vector,
                    adirection,
                    oncurve,
                    &tan2,
                    &tan1,
                    &line_bis_vec,
                    line1,
                    line2,
                );
            }
            _ => return Err(BisectorBisecAnaError::NotDone),
        }
        Ok(())
    }

    /// `case 1` of `Perform` (`cxx:355-678`): circle - circle.
    #[allow(clippy::too_many_arguments)]
    fn perform_circle_circle(
        &mut self,
        first_curve: &dyn Curve2d,
        second_curve: &dyn Curve2d,
        apoint: &GpPnt2d,
        first_vector: &GpVec2d,
        second_vector: &GpVec2d,
        adirection: f64,
        tolerance: f64,
        oncurve: bool,
        tan1: &GpVec2d,
        tan2: &GpVec2d,
        line_bis_vec: &GpVec2d,
        mut circle1: GpCirc2d,
        mut circle2: GpCirc2d,
    ) -> Result<(), BisectorBisecAnaError> {
        let pre_conf = CONFUSION;
        let mut radius1 = circle1.radius();
        let mut radius2 = circle2.radius();
        let mut thesense = false;
        // `cxx:362-433`: the two circles are mixed (same centre, same radius).
        if circle1
            .location()
            .is_equal(&circle2.location(), pre_conf)
            && (radius1 - radius2).abs() <= pre_conf
        {
            let p1 = first_curve.value(first_curve.last_parameter());
            let p2 = second_curve.value(second_curve.first_parameter());
            let p_mil = GpPnt2d::new((p1.x() + p2.x()) / 2.0, (p1.y() + p2.y()) / 2.0);
            let line = if !circle1.location().is_equal(&p_mil, pre_conf) {
                GpLin2d::from_pnt_dir(
                    p_mil,
                    GpDir2d::new(
                        circle1.location().x() - p_mil.x(),
                        circle1.location().y() - p_mil.y(),
                    )
                    .unwrap_or_default(),
                )
            } else if radius1 >= pre_conf {
                GpLin2d::from_pnt_dir(
                    circle1.location(),
                    GpDir2d::new(
                        p1.y() - circle1.location().y(),
                        circle1.location().x() - p1.x(),
                    )
                    .unwrap_or_default(),
                )
            } else {
                GpLin2d::from_pnt_dir(circle1.location(), GpDir2d::default())
            };
            let solution = GccIntBisec::Line(line);
            let res = if oncurve {
                Self::distance(
                    apoint,
                    &solution,
                    tan2,
                    tan1,
                    line_bis_vec,
                    adirection,
                    false,
                    false,
                )
            } else {
                Self::distance(
                    apoint,
                    &solution,
                    first_vector,
                    second_vector,
                    line_bis_vec,
                    adirection,
                    false,
                    false,
                )
            };
            let parameter = res.parameter;
            let sense = res.sense;
            let bisectorcurve: Arc<dyn Curve2d> = Arc::from(
                Geom2dLine::new(GpAx2d::new(line.location(), *line.direction())).clone_dyn(),
            );
            self.thebisector = Some(if !sense {
                Geom2dTrimmedCurve::new(bisectorcurve, parameter, -INFINITE)
            } else {
                let mut parameter2 = clib2d::parameter_lin2d(&line, &circle1.location());
                parameter2 += 1.0e-8;
                Geom2dTrimmedCurve::new(bisectorcurve, parameter, parameter2)
            });
            return Ok(());
        }

        if radius1 < radius2 {
            std::mem::swap(&mut circle1, &mut circle2);
            std::mem::swap(&mut radius1, &mut radius2);
        }

        // `cxx:449-477`: small reframing of circles. When the circles are
        // OnCurve and almost tangent they become tangent.
        let entre_axe = circle1.location().distance(&circle2.location());
        let mut d1 = 0.5 * (radius1 - entre_axe - radius2);
        let mut circles_tangent = false;
        if oncurve && d1.abs() < pre_conf && tan1.is_parallel(tan2, 1.0e-8) {
            circle1.set_radius(radius1 - d1);
            circle2.set_radius(radius2 + d1);
            circles_tangent = true;
        } else {
            d1 = 0.5 * (radius1 - entre_axe + radius2);
            if oncurve && d1.abs() < pre_conf && tan1.is_parallel(tan2, 1.0e-8) {
                circle1.set_radius(radius1 - d1);
                circle2.set_radius(radius2 - d1);
                circles_tangent = true;
            }
        }

        let bisector = GccAnaCirc2dBisec::new(&circle1, &circle2);
        let mut distancemini = INFINITE;
        let mut the_sol: Option<GccIntBisec> = None;
        let mut firstparameter = 0.0;

        if bisector.is_done() {
            let nbsolution = bisector.nb_solutions()?;
            for i in 1..=nbsolution {
                let mut solution = bisector.this_solution(i)?;
                degenerate(&mut solution, tolerance);
                // `sense` is `false` for the decisive call below (`cxx:489-502`);
                // the first pass only contributes its status.
                let _sense = false;
                let ok = if oncurve {
                    Self::distance(
                        apoint,
                        &solution,
                        tan1,
                        tan2,
                        line_bis_vec,
                        adirection,
                        true,
                        false,
                    )
                    .status
                } else {
                    true
                };

                if ok {
                    // OCCT sets `sense = false` before the decisive call
                    // (`cxx:502`); the value is overwritten by `Distance` unless
                    // both `asense` branches are skipped.
                    let sense = false;
                    let r = if oncurve {
                        Self::distance(
                            apoint,
                            &solution,
                            tan2,
                            tan1,
                            line_bis_vec,
                            adirection,
                            sense,
                            false,
                        )
                    } else {
                        Self::distance(
                            apoint,
                            &solution,
                            first_vector,
                            second_vector,
                            line_bis_vec,
                            adirection,
                            sense,
                            false,
                        )
                    };
                    if r.distance <= distancemini {
                        the_sol = Some(solution);
                        firstparameter = r.parameter;
                        thesense = r.sense;
                        distancemini = r.distance;
                    }
                }
            }
            if let Some(the_sol) = the_sol {
                match the_sol {
                    GccIntBisec::Line(gpline) => {
                        let mut bisectorcurve: Arc<dyn Curve2d> = Arc::from(
                            Geom2dLine::new(GpAx2d::new(gpline.location(), *gpline.direction()))
                                .clone_dyn(),
                        );
                        let mut secondparameter = INFINITE;
                        if !thesense {
                            secondparameter = -INFINITE;
                        }
                        if oncurve && circles_tangent {
                            // `cxx:565-621`: the bisector is a segment between
                            // the two circle centres of the same radius when the
                            // circles are tangent.
                            if distancemini > tolerance {
                                let a_ploc = gpline.location();
                                let a_newdir = GpDir2d::new(
                                    apoint.x() - a_ploc.x(),
                                    apoint.y() - a_ploc.y(),
                                )
                                .unwrap_or_default();
                                let a_newlin = GpLin2d::from_pnt_dir(a_ploc, a_newdir);
                                let a_cc2 = circle2.location();
                                let a_tolconf = 1.0e-3;
                                if a_newlin.distance(&a_cc2) <= a_tolconf {
                                    // OCCT also stores `aNewDMin` into
                                    // `distancemini`, which is never read again.
                                    firstparameter =
                                        clib2d::parameter_lin2d(&a_newlin, apoint);
                                    bisectorcurve = Arc::from(
                                        Geom2dLine::new(GpAx2d::new(
                                            a_newlin.location(),
                                            *a_newlin.direction(),
                                        ))
                                        .clone_dyn(),
                                    );
                                }
                            }
                            if tan1.dot(tan2) < 0.0 {
                                // Flat and not turned back.
                                let par1 =
                                    clib2d::parameter_lin2d(&gpline, &circle1.location());
                                let par2 =
                                    clib2d::parameter_lin2d(&gpline, &circle2.location());
                                let min_par = par1.min(par2);
                                let max_par = par1.max(par2);
                                if !thesense {
                                    if max_par < firstparameter {
                                        secondparameter = max_par - 1.0e-8;
                                    } else if min_par < firstparameter {
                                        secondparameter = min_par - 1.0e-8;
                                    }
                                } else if min_par > firstparameter {
                                    secondparameter = min_par + 1.0e-8;
                                } else if max_par > firstparameter {
                                    secondparameter = max_par + 1.0e-8;
                                }
                            }
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
                            Geom2dTrimmedCurve::new(bisectorcurve, firstparameter, -INFINITE)
                        } else {
                            Geom2dTrimmedCurve::new(bisectorcurve, firstparameter, INFINITE)
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
        Ok(())
    }

    /// `case 2` of `Perform` (`cxx:684-834`): circle - straight line.
    #[allow(clippy::too_many_arguments)]
    fn perform_circle_line(
        &mut self,
        apoint: &GpPnt2d,
        first_vector: &GpVec2d,
        second_vector: &GpVec2d,
        adirection: f64,
        ajointype: JoinType,
        tolerance: f64,
        oncurve: bool,
        tan1: &GpVec2d,
        tan2: &GpVec2d,
        line_bis_vec: &GpVec2d,
        mut circle1: GpCirc2d,
        line2: GpLin2d,
    ) -> Result<(), BisectorBisecAnaError> {
        let pre_conf = CONFUSION;
        if oncurve {
            let radius1 = circle1.radius();
            let d1 = line2.distance(&circle1.location()) - radius1;
            if d1.abs() < pre_conf && tan1.is_parallel(tan2, 1.0e-8) {
                circle1.set_radius(radius1 + d1);
            }
        }

        let bisector = GccAnaCircLin2dBisec::new(&circle1, &line2);
        let mut distancemini = INFINITE;
        let mut the_sol: Option<GccIntBisec> = None;
        let mut firstparameter = 0.0;
        let mut thesense = false;

        if bisector.is_done() {
            let nbsolution = bisector.nb_solutions();
            for i in 1..=nbsolution {
                let mut solution = bisector.this_solution(i)?;
                let degen = degenerate(&mut solution, tolerance);
                let kind = match &solution {
                    GccIntBisec::Line(_) => "line".to_string(),
                    GccIntBisec::Circle(_) => "circle".to_string(),
                    GccIntBisec::Hyperbola(_) => "hyperbola".to_string(),
                    GccIntBisec::Ellipse(_) => "ellipse".to_string(),
                    GccIntBisec::Parabola(p) => format!("parabola(focal={:?})", p.focal),
                    GccIntBisec::Point(_) => "point".to_string(),
                };
                eprintln!("PCL i={i} tol={tolerance} degen={degen} kind={kind}");
                // `cxx:711-718`: the first pass only contributes `ok`; the second
                // one starts from `sense = false`.
                let ok = Self::distance(
                    apoint,
                    &solution,
                    tan1,
                    tan2,
                    line_bis_vec,
                    adirection,
                    true,
                    false,
                )
                .status;
                if ok || !oncurve {
                    let sense = false;
                    let r = if oncurve {
                        Self::distance(
                            apoint,
                            &solution,
                            tan2,
                            tan1,
                            line_bis_vec,
                            adirection,
                            sense,
                            false,
                        )
                    } else {
                        Self::distance(
                            apoint,
                            &solution,
                            first_vector,
                            second_vector,
                            line_bis_vec,
                            adirection,
                            sense,
                            false,
                        )
                    };
                    if r.distance <= distancemini {
                        the_sol = Some(solution);
                        firstparameter = r.parameter;
                        thesense = r.sense;
                        distancemini = r.distance + 1.0e-8;
                    }
                }
            }
            if let Some(the_sol) = the_sol {
                match the_sol {
                    GccIntBisec::Line(gpline) => {
                        let circlecenter = circle1.location();
                        let mut secondparameter = clib2d::parameter_lin2d(&gpline, &circlecenter);
                        let bisectorcurve: Arc<dyn Curve2d> = Arc::from(
                            Geom2dLine::new(GpAx2d::new(gpline.location(), *gpline.direction()))
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
                    GccIntBisec::Parabola(pb) => {
                        let bisectorcurve: Arc<dyn Curve2d> =
                            Arc::from(Geom2dParabola::new(pb).clone_dyn());
                        let apex = clib2d::parabola_value_ax22d(0.0, &pb.pos, pb.focal);
                        let firstpnt =
                            clib2d::parabola_value_ax22d(firstparameter, &pb.pos, pb.focal);
                        let chordlen = apex.distance(&firstpnt);
                        let tol_par = 1.0e-5;
                        let mut secondparameter = INFINITE;
                        if !thesense {
                            if ajointype == JoinType::Intersection
                                && tol_par < firstparameter
                                && chordlen >= circle1.radius()
                            {
                                // The first parameter is too far from the apex.
                                secondparameter = 0.0;
                            }
                            self.thebisector = Some(Geom2dTrimmedCurve::new(
                                bisectorcurve,
                                firstparameter,
                                -secondparameter,
                            ));
                        } else {
                            if ajointype == JoinType::Intersection
                                && firstparameter < -tol_par
                                && chordlen >= circle1.radius()
                            {
                                secondparameter = 0.0;
                            }
                            self.thebisector = Some(Geom2dTrimmedCurve::new(
                                bisectorcurve,
                                firstparameter,
                                secondparameter,
                            ));
                        }
                    }
                    _ => {}
                }
            }
        }
        Ok(())
    }

    /// `case 3` of `Perform` (`cxx:838-957`): straight line - straight line.
    #[allow(clippy::too_many_arguments)]
    fn perform_line_line(
        &mut self,
        apoint: &GpPnt2d,
        first_vector: &GpVec2d,
        second_vector: &GpVec2d,
        adirection: f64,
        oncurve: bool,
        tan2: &GpVec2d,
        tan1: &GpVec2d,
        line_bis_vec: &GpVec2d,
        line1: GpLin2d,
        line2: GpLin2d,
    ) {
        let pre_conf = CONFUSION;
        let direc1 = *line1.direction();
        let direc2 = *line2.direction();
        // `distancemini` / `distanceptsol` of OCCT (`cxx:843-844`) are local
        // bookkeeping that no branch reads back.
        let mut thesense = false;

        if direc1.is_parallel(&direc2, 1.0e-8) {
            // Same criterion as `MAT2d_Circuit::InitOpen` (`cxx:845-850`).
            let line = if line1.distance(&line2.location()) / 2.0 <= pre_conf {
                GpLin2d::from_pnt_dir(
                    *apoint,
                    GpDir2d::new(-line1.direction().y(), line1.direction().x())
                        .unwrap_or_default(),
                )
            } else {
                GpLin2d::from_pnt_dir(*apoint, direc2)
            };
            let solution = GccIntBisec::Line(line);
            let sense = false;
            let r = if oncurve {
                Self::distance(
                    apoint,
                    &solution,
                    tan2,
                    tan1,
                    line_bis_vec,
                    adirection,
                    sense,
                    false,
                )
            } else {
                Self::distance(
                    apoint,
                    &solution,
                    first_vector,
                    second_vector,
                    line_bis_vec,
                    adirection,
                    sense,
                    false,
                )
            };
            let firstparameter = r.parameter;
            let bisectorcurve: Arc<dyn Curve2d> = Arc::from(
                Geom2dLine::new(GpAx2d::new(line.location(), *line.direction())).clone_dyn(),
            );
            self.thebisector = Some(if !r.sense {
                Geom2dTrimmedCurve::new(bisectorcurve, firstparameter, -INFINITE)
            } else {
                Geom2dTrimmedCurve::new(bisectorcurve, firstparameter, INFINITE)
            });
        } else {
            let l = GpLin2d::from_pnt_dir(
                *apoint,
                GpDir2d::new(direc2.x() - direc1.x(), direc2.y() - direc1.y())
                    .unwrap_or_default(),
            );
            let solution = GccIntBisec::Line(l);
            let sense = false;
            let r = if oncurve {
                Self::distance(
                    apoint,
                    &solution,
                    tan2,
                    tan1,
                    line_bis_vec,
                    adirection,
                    sense,
                    false,
                )
            } else {
                // `cxx:926-935`: the two-line variant, which skips the sector
                // block of `Distance`.
                Self::distance(
                    apoint,
                    &solution,
                    first_vector,
                    second_vector,
                    line_bis_vec,
                    adirection,
                    sense,
                    true,
                )
            };
            if r.status || !oncurve {
                thesense = r.sense;
            }
            let bisectorcurve: Arc<dyn Curve2d> = Arc::from(
                Geom2dLine::new(GpAx2d::new(l.location(), *l.direction())).clone_dyn(),
            );
            self.thebisector = Some(if !thesense {
                Geom2dTrimmedCurve::new(bisectorcurve, 0.0, -INFINITE)
            } else {
                Geom2dTrimmedCurve::new(bisectorcurve, 0.0, INFINITE)
            });
        }
    }
}
