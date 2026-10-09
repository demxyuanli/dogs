//! `Extrema_ExtCC2d` -- all extremal distances between two 2D curves.
//! Source: `Extrema_ExtCC2d.cxx` (684 lines), `Extrema_ExtCC2d.hxx:36-143`.
//!
//! `Perform` (`Extrema_ExtCC2d.cxx:101-453`) dispatches on
//! `Extrema_Curve2dTool::GetType` of both curves. The elementary pairs build an
//! `Extrema_ExtElC2d` (`:125-430`); every other pair goes through
//! `Extrema_ECC2d` (=`Extrema_GGenExtCC`, `:167-178`, `:229-240`, `:295-306`,
//! `:360-371`, `:416-427`, `:435-451`). The trimming `[U1,U2]x[V1,V2]` is
//! applied in `Results` (`:609-659`), *not* inside the solver: the solver is
//! constructed from the two adaptors, whose range is the curves' own.
//!
//! The elementary arms build the faithful `Extrema_ExtElC2d`
//! (`super::ext_el_c2d`, source `Extrema_ExtElC2d.cxx`, 542 lines) and go
//! through the `Results(const Extrema_ExtElC2d&, ...)` overload (`:527-605`);
//! `elementary_ext_el_c2d` reproduces the `(inverse, Period1, Period2)` triple
//! of every such arm. The `Tol`/`AngTol` arguments those arms pass are unused by
//! the `Extrema_ExtElC2d` constructors they call (`Extrema_ExtElC2d.cxx:46`,
//! `:105` take an unnamed `double`), so `TolC1`/`TolC2` are not threaded
//! through.
//!
//! This is inert for `curve_ops::curve2d_intersections`: a pair of elementary
//! curves is solved by the analytic `IntAna2d_AnaIntersection` route before this
//! function is reached. It is exercised by
//! `int_conic_conic_line_ellipse::line_ellipse_geometric_intersection`'s
//! `D < 0` tangent branch (`IntCurve_IntConicConic_1.cxx:2702-2727`), which
//! builds an `Extrema_ExtElC2d(Line, Ellipse)` directly.

use super::prelude::*;

use super::analytic_solvers::PI;
use super::curve2d_tool::{
    adaptor_basis2d, adaptor_is_periodic2d, adaptor_period2d, adaptor_type2d,
};
use super::curve_curve::in_period2d;
use super::ext_el_c2d::ExtremaExtElC2d;
use super::general_extrema::GGenExtCC;

use occt_core::kernel::geomabs::CurveType;

const TWO_PI: f64 = 2.0 * PI;

/// The result of `Extrema_ExtCC2d` for one curve pair: the extrema
/// (`mypoints`/`mySqDist`) plus the four trimmed endpoint distances.
#[derive(Debug, Clone)]
pub struct ExtremaExtCC2d {
    done: bool,
    is_par: bool,
    /// `myIsFindSingleSolution` (`Extrema_ExtCC2d.hxx:121`, default `false` at
    /// `Extrema_ExtCC2d.cxx:35`, `:58`, `:78`); every solver arm installs it with
    /// `SetSingleSolutionFlag(GetSingleSolutionFlag())` (`:169`, `:206`, `:215`,
    /// `:224`, `:231`, `:271`, `:280`, `:290`, `:297`, `:337`, `:346`, `:355`,
    /// `:362`, `:418`, `:437`).
    single_solution: bool,
    /// One entry per extremum: `(u on C1, P1, u on C2, P2)`, i.e. the
    /// `mypoints` pairs of `Extrema_ExtCC2d.cxx:500-501`.
    points: Vec<(f64, GpPnt2d, f64, GpPnt2d)>,
    /// `mySqDist`; entry `n-1` is `SquareDistance(n)`.
    sq_dist: Vec<f64>,
    dist11: f64,
    dist12: f64,
    dist21: f64,
    dist22: f64,
}

impl ExtremaExtCC2d {
    /// `Extrema_ExtCC2d(C1, C2, TolC1, TolC2)`
    /// (`Extrema_ExtCC2d.cxx:54-66`): the whole range of both curves.
    pub fn new(c1: &dyn Curve2d, c2: &dyn Curve2d) -> Self {
        Self::new_range(
            c1,
            c2,
            c1.first_parameter(),
            c1.last_parameter(),
            c2.first_parameter(),
            c2.last_parameter(),
        )
    }

    /// `Extrema_ExtCC2d(C1, C2, U1, U2, V1, V2, TolC1, TolC2)`
    /// (`Extrema_ExtCC2d.cxx:70-82` + `Initialize` `:86-97`).
    ///
    /// `TolC1`/`TolC2` are only consumed by the `Extrema_ExtElC2d` arms
    /// (`:107`, `:136-139`, `:384-399`), whose constructors take an unnamed
    /// `double` and ignore it (`Extrema_ExtElC2d.cxx:46`, `:105`), so they are
    /// not parameters here.
    pub fn new_range(
        c1: &dyn Curve2d,
        c2: &dyn Curve2d,
        u1: f64,
        u2: f64,
        v1: f64,
        v2: f64,
    ) -> Self {
        let mut out = Self {
            done: false,
            is_par: false,
            single_solution: false,
            points: Vec::new(),
            sq_dist: Vec::new(),
            dist11: 0.0,
            dist12: 0.0,
            dist21: 0.0,
            dist22: 0.0,
        };
        out.perform(c1, c2, u1, u2, v1, v2);
        out
    }

    /// `IsDone()` (`Extrema_ExtCC2d.cxx:457-460`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `IsParallel()` (`:663-670`); OCCT throws when not done, the port
    /// reports `false`.
    pub fn is_parallel(&self) -> bool {
        self.done && self.is_par
    }

    /// `NbExt()` (`:479-486`).
    pub fn nb_ext(&self) -> usize {
        if !self.done {
            return 0;
        }
        self.sq_dist.len()
    }

    /// `SquareDistance(N)` (`:464-475`), `N` 1-based.
    pub fn square_distance(&self, n: usize) -> f64 {
        self.sq_dist[n - 1]
    }

    /// `Points(N, P1, P2)` (`:490-502`): `(u1, P1, u2, P2)`.
    pub fn points(&self, n: usize) -> (f64, GpPnt2d, f64, GpPnt2d) {
        self.points[n - 1]
    }

    /// `TrimmedSquareDistances` (`:506-523`), for the two ends of each curve.
    pub fn trimmed_square_distances(&self) -> (f64, f64, f64, f64) {
        (self.dist11, self.dist12, self.dist21, self.dist22)
    }

    /// `SetSingleSolutionFlag` (`:674-677`).
    pub fn set_single_solution_flag(&mut self, flag: bool) {
        self.single_solution = flag;
    }

    /// `GetSingleSolutionFlag` (`:681-684`).
    pub fn single_solution_flag(&self) -> bool {
        self.single_solution
    }

    /// `Extrema_ExtCC2d::Perform` (`Extrema_ExtCC2d.cxx:101-453`).
    fn perform(&mut self, c1: &dyn Curve2d, c2: &dyn Curve2d, u1: f64, u2: f64, v1: f64, v2: f64) {
        self.points.clear();
        self.sq_dist.clear();
        let type1 = adaptor_type2d(c1);
        let type2 = adaptor_type2d(c2);

        let p1f = c1.d0(u1); // `:117`
        let p1l = c1.d0(u2); // `:118`
        let p2f = c2.d0(v1); // `:119`
        let p2l = c2.d0(v2); // `:120`

        // `:125-430`: an elementary pair builds an `Extrema_ExtElC2d` and goes
        // through the `Results(const Extrema_ExtElC2d&, ...)` overload
        // (`:527-605`); every other pair -- including the arms OCCT leaves
        // commented out -- runs `Extrema_ECC2d` below.
        if let Some((ext, inverse, period1, period2)) =
            elementary_ext_el_c2d(type1, type2, c1, c2)
        {
            self.results_el_c2d(&ext, inverse, u1, u2, v1, v2, period1, period2);
            if self.done {
                // `:600-603`
                self.dist11 = p1f.square_distance(&p2f);
                self.dist12 = p1f.square_distance(&p2l);
                self.dist21 = p1l.square_distance(&p2f);
                self.dist22 = p1l.square_distance(&p2l);
            }
            return;
        }

        // `Results(*aParamSolver, U11, U12, U21, U22, Period1, Period2)`
        // (`Extrema_ExtCC2d.cxx:609-659`); the periods are the branch values of
        // the OCCT switch below.
        let (period1, period2) = match type1 {
            CurveType::Circle => match type2 {
                CurveType::Line => (TWO_PI, 0.0),                     // `:134-140`
                CurveType::Circle => (TWO_PI, TWO_PI),                // `:142-146`
                CurveType::Ellipse => (TWO_PI, TWO_PI),               // `:148-152`
                CurveType::Parabola => (TWO_PI, 0.0),                 // `:154-158`
                CurveType::Hyperbola => (TWO_PI, 0.0),                // `:160-165`
                _ => (TWO_PI, periodic_period2d(c2)),                 // `:167-178`
            },
            CurveType::Ellipse => match type2 {
                CurveType::Line => (TWO_PI, 0.0),                     // `:190-194`
                CurveType::Circle => (TWO_PI, TWO_PI),                // `:197-202`
                CurveType::Ellipse => (TWO_PI, TWO_PI),               // `:204-208`
                CurveType::Parabola => (TWO_PI, 0.0),                 // `:211-218`
                CurveType::Hyperbola => (TWO_PI, 0.0),                // `:220-227`
                _ => (TWO_PI, periodic_period2d(c2)),                 // `:229-240`
            },
            CurveType::Parabola => match type2 {
                CurveType::Line => (0.0, 0.0),                        // `:252-256`
                CurveType::Circle => (0.0, TWO_PI),                   // `:259-263`
                CurveType::Ellipse => (0.0, TWO_PI),                  // `:266-274`
                CurveType::Parabola => (0.0, 0.0),                    // `:276-283`
                CurveType::Hyperbola => (0.0, 0.0),                   // `:285-293`
                _ => (0.0, periodic_period2d(c2)),                    // `:295-306`
            },
            CurveType::Hyperbola => match type2 {
                CurveType::Line => (0.0, 0.0),                        // `:318-322`
                CurveType::Circle => (0.0, TWO_PI),                   // `:325-329`
                CurveType::Ellipse => (0.0, TWO_PI),                  // `:332-339`
                CurveType::Parabola => (0.0, 0.0),                    // `:342-348`
                CurveType::Hyperbola => (0.0, 0.0),                   // `:351-357`
                _ => (0.0, periodic_period2d(c2)),                    // `:360-369`
            },
            CurveType::Line => match type2 {
                CurveType::Line => (0.0, 0.0),                        // `:383-388`
                CurveType::Circle => (0.0, TWO_PI),                   // `:390-394`
                CurveType::Ellipse => (0.0, TWO_PI),                  // `:397-401`
                CurveType::Parabola => (0.0, 0.0),                    // `:403-407`
                CurveType::Hyperbola => (0.0, 0.0),                   // `:409-414`
                _ => (0.0, periodic_period2d(c2)),                    // `:416-427`
            },
            // `default:` -- Bezier / BSpline / Offset / Other (`:432-451`).
            _ => (periodic_period2d(c1), periodic_period2d(c2)),        // `:439-448`
        };

        // `aParamSolver = Extrema_ECC2d(C1, *myC)` -- the two adaptors, so the
        // search box is each curve's own range, not `[U1,U2]x[V1,V2]`.
        let mut alg = GGenExtCC::new(
            c1,
            c2,
            (c1.first_parameter(), c2.first_parameter()),
            (c1.last_parameter(), c2.last_parameter()),
        );
        alg.set_single_solution_flag(self.single_solution_flag());
        if alg.perform().is_err() {
            self.done = false;
            return;
        }

        // `Results(const Extrema_ECC2d& AlgExt, ...)` (`:609-659`).
        self.done = alg.is_done();
        if !self.done {
            return;
        }
        self.is_par = alg.is_parallel();
        for i in 1..=alg.nb_ext() {
            // `:629-639`
            let (u, pt1, v, pt2) = alg.points(i);
            let mut uu = u;
            if period1 != 0.0 {
                uu = in_period2d(uu, u1, period1); // `ElCLib::InPeriod(U, U11, U11 + Period1)`
            }
            let mut uu2 = v;
            if period2 != 0.0 {
                uu2 = in_period2d(uu2, v1, period2); // `ElCLib::InPeriod(U, U21, U21 + Period2)`
            }
            if uu >= u1 - PCONFUSION
                && uu <= u2 + PCONFUSION
                && uu2 >= v1 - PCONFUSION
                && uu2 <= v2 + PCONFUSION
            {
                // `:641-651`
                self.sq_dist.push(alg.square_distance(i));
                self.points.push((uu, pt1, uu2, pt2));
            }
        }

        // `:654-657`
        self.dist11 = p1f.square_distance(&p2f);
        self.dist12 = p1f.square_distance(&p2l);
        self.dist21 = p1l.square_distance(&p2f);
        self.dist22 = p1l.square_distance(&p2l);
    }

    /// `Results(const Extrema_ExtElC2d& AlgExt, Ut11, Ut12, Ut21, Ut22, Period1,
    /// Period2)` (`Extrema_ExtCC2d.cxx:527-605`).
    ///
    /// `inverse` mirrors the OCCT flag of the type switch (`:135`, `:191`,
    /// `:198`, `:253`, `:260`, `:319`, `:326`): the `Extrema_ExtElC2d` was then
    /// built with the two operands swapped, so `P1` belongs to the second curve
    /// and `P2` to the first.
    #[allow(clippy::too_many_arguments)]
    fn results_el_c2d(
        &mut self,
        ext: &ExtremaExtElC2d,
        inverse: bool,
        ut11: f64,
        ut12: f64,
        ut21: f64,
        ut22: f64,
        period1: f64,
        period2: f64,
    ) {
        self.done = ext.is_done(); // `:539`
        self.is_par = ext.is_parallel(); // `:540`
        if !self.done || self.is_par {
            // `:541-544`
            return;
        }
        for i in 1..=ext.nb_ext() {
            let (p1, p2) = ext.points(i); // `:549`
            let (u, u2) = if !inverse {
                // `:551-563`
                let mut u = p1.parameter();
                if period1 != 0.0 {
                    u = in_period2d(u, ut11, period1);
                }
                let mut u2 = p2.parameter();
                if period2 != 0.0 {
                    u2 = in_period2d(u2, ut21, period2);
                }
                (u, u2)
            } else {
                // `:564-576`
                let mut u2 = p1.parameter();
                if period2 != 0.0 {
                    u2 = in_period2d(u2, ut21, period2);
                }
                let mut u = p2.parameter();
                if period1 != 0.0 {
                    u = in_period2d(u, ut11, period1);
                }
                (u, u2)
            };
            // `:577-578`
            if u >= ut11 - PCONFUSION
                && u <= ut12 + PCONFUSION
                && u2 >= ut21 - PCONFUSION
                && u2 <= ut22 + PCONFUSION
            {
                // `:579-594`: `mypoints` always stores the first curve's
                // parameter/point first.
                self.sq_dist.push(ext.square_distance(i));
                if !inverse {
                    self.points.push((u, p1.value(), u2, p2.value()));
                } else {
                    self.points.push((u, p2.value(), u2, p1.value()));
                }
            }
        }
    }
}

/// The `Extrema_ExtElC2d` arms of `Extrema_ExtCC2d::Perform`'s type switch
/// (`Extrema_ExtCC2d.cxx:125-430`), as `(solver, inverse, Period1, Period2)`.
/// `None` selects the general `Extrema_ECC2d` engine, which is also what the
/// arms OCCT leaves commented out (`:210-227`, `:266-293`, `:332-357`) do.
fn elementary_ext_el_c2d(
    type1: CurveType,
    type2: CurveType,
    c1: &dyn Curve2d,
    c2: &dyn Curve2d,
) -> Option<(ExtremaExtElC2d, bool, f64, f64)> {
    let b1 = adaptor_basis2d(c1);
    let b2 = adaptor_basis2d(c2);
    match type1 {
        // The first curve is a circle (`:130-180`).
        CurveType::Circle => {
            let c = b1.gp_circ2d()?;
            match type2 {
                // `:134-140`, `inverse = true`.
                CurveType::Line => {
                    let l = b2.gp_lin2d()?;
                    Some((ExtremaExtElC2d::new_line_circle(&l, &c), true, TWO_PI, 0.0))
                }
                // `:142-146`.
                CurveType::Circle => {
                    let o = b2.gp_circ2d()?;
                    Some((ExtremaExtElC2d::new_circle_circle(&c, &o), false, TWO_PI, TWO_PI))
                }
                // `:148-152`.
                CurveType::Ellipse => {
                    let e = b2.gp_elips2d()?;
                    Some((ExtremaExtElC2d::new_circle_ellipse(&c, &e), false, TWO_PI, TWO_PI))
                }
                // `:154-158`.
                CurveType::Parabola => {
                    let p = b2.gp_parab2d()?;
                    Some((ExtremaExtElC2d::new_circle_parabola(&c, &p), false, TWO_PI, 0.0))
                }
                // `:160-165`.
                CurveType::Hyperbola => {
                    let h = b2.gp_hypr2d()?;
                    Some((ExtremaExtElC2d::new_circle_hyperbola(&c, &h), false, TWO_PI, 0.0))
                }
                _ => None, // `:167-178`
            }
        }
        // The first curve is an ellipse (`:186-241`).
        CurveType::Ellipse => {
            let e = b1.gp_elips2d()?;
            match type2 {
                // `:190-194`, `inverse = true`.
                CurveType::Line => {
                    let l = b2.gp_lin2d()?;
                    Some((ExtremaExtElC2d::new_line_ellipse(&l, &e), true, TWO_PI, 0.0))
                }
                // `:197-202`, `inverse = true`.
                CurveType::Circle => {
                    let c = b2.gp_circ2d()?;
                    Some((ExtremaExtElC2d::new_circle_ellipse(&c, &e), true, TWO_PI, TWO_PI))
                }
                // `:204-208` and `:211-240` use `Extrema_ECC2d`.
                _ => None,
            }
        }
        // The first curve is a parabola (`:247-308`).
        CurveType::Parabola => {
            let p = b1.gp_parab2d()?;
            match type2 {
                // `:252-256`, `inverse = true`.
                CurveType::Line => {
                    let l = b2.gp_lin2d()?;
                    Some((ExtremaExtElC2d::new_line_parabola(&l, &p), true, 0.0, 0.0))
                }
                // `:259-263`, `inverse = true`.
                CurveType::Circle => {
                    let c = b2.gp_circ2d()?;
                    Some((ExtremaExtElC2d::new_circle_parabola(&c, &p), true, 0.0, TWO_PI))
                }
                _ => None, // `:266-306`
            }
        }
        // The first curve is a hyperbola (`:314-372`).
        CurveType::Hyperbola => {
            let h = b1.gp_hypr2d()?;
            match type2 {
                // `:318-322`, `inverse = true`.
                CurveType::Line => {
                    let l = b2.gp_lin2d()?;
                    Some((ExtremaExtElC2d::new_line_hyperbola(&l, &h), true, 0.0, 0.0))
                }
                // `:325-329`, `inverse = true`.
                CurveType::Circle => {
                    let c = b2.gp_circ2d()?;
                    Some((ExtremaExtElC2d::new_circle_hyperbola(&c, &h), true, 0.0, TWO_PI))
                }
                _ => None, // `:332-369`
            }
        }
        // The first curve is a line (`:379-428`).
        CurveType::Line => {
            let l = b1.gp_lin2d()?;
            match type2 {
                // `:383-388`.
                CurveType::Line => {
                    let o = b2.gp_lin2d()?;
                    Some((ExtremaExtElC2d::new_line_line(&l, &o), false, 0.0, 0.0))
                }
                // `:390-394`.
                CurveType::Circle => {
                    let c = b2.gp_circ2d()?;
                    Some((ExtremaExtElC2d::new_line_circle(&l, &c), false, 0.0, TWO_PI))
                }
                // `:397-401`.
                CurveType::Ellipse => {
                    let e = b2.gp_elips2d()?;
                    Some((ExtremaExtElC2d::new_line_ellipse(&l, &e), false, 0.0, TWO_PI))
                }
                // `:403-407`.
                CurveType::Parabola => {
                    let p = b2.gp_parab2d()?;
                    Some((ExtremaExtElC2d::new_line_parabola(&l, &p), false, 0.0, 0.0))
                }
                // `:409-414`.
                CurveType::Hyperbola => {
                    let h = b2.gp_hypr2d()?;
                    Some((ExtremaExtElC2d::new_line_hyperbola(&l, &h), false, 0.0, 0.0))
                }
                _ => None, // `:416-427`
            }
        }
        // `default:` -- Bezier / BSpline / Offset / Other (`:432-451`).
        _ => None,
    }
}

/// `Extrema_Curve2dTool::Period` guarded by `IsPeriodic`
/// (`Extrema_ExtCC2d.cxx:172-175`, `:233-237`, `...`): the adaptor period, or 0
/// for a non-periodic curve.
fn periodic_period2d(c: &dyn Curve2d) -> f64 {
    if adaptor_is_periodic2d(c) {
        adaptor_period2d(c)
    } else {
        0.0
    }
}