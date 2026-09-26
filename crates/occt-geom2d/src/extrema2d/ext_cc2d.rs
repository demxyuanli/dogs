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
//! **UNPORTED**: `Extrema_ExtElC2d` (`Extrema_ExtElC2d.cxx`, 542 lines) has no
//! port, so the elementary branches below call the general `Extrema_ECC2d`
//! engine with the `Period1`/`Period2` values OCCT's branch would have passed to
//! `Results`. This is inert for `curve_ops::curve2d_intersections`: a pair of
//! elementary curves is solved by the analytic `IntAna2d_AnaIntersection`
//! route before this function is reached.

use super::prelude::*;

use super::analytic_solvers::PI;
use super::curve2d_tool::{adaptor_is_periodic2d, adaptor_period2d, adaptor_type2d};
use super::curve_curve::in_period2d;
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
    /// (`:107`, `:136-139`, `:384-399`), which are UNPORTED (see the module
    /// header), so they are not parameters here.
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