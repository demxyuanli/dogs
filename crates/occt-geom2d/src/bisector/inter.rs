//! `Bisector_Inter` (`Bisector_Inter.cxx:45-463`, `.hxx:32-90`): the
//! intersection between two bisectors.
//!
//! `Perform` cuts both bisectors at their continuity intervals (clipped by the
//! caller's `IntRes2d_Domain`s), turns the extension parts — `IsExtendAtStart`
//! on the first interval, `IsExtendAtEnd` on the last (`cxx:136-149`) — into
//! line segments (`ConstructSegment`, `cxx:61-72`), and intersects the parts
//! pairwise.
//!
//! `SinglePerform` picks the intersector from the dynamic types
//! (`cxx:234-322`):
//! * a `Bisector_BisecAna` is replaced by its `Geom2dCurve()` basis and the
//!   pair goes through `Geom2dInt_GInter`, after `TestBound` for a line /
//!   non-line pair;
//! * a `Bisector_BisecPC` pair goes straight to `Geom2dInt_GInter`;
//! * two `Bisector_BisecCC` with `ComunElement` go to `NeighbourPerform`, which
//!   re-parameterises the second one along the shared guide (`ChangeGuide`) and
//!   solves `Bisector_FunctionInter` with `math_BissecNewton`;
//! * the remaining (segment) cases run `TestBound` then `Geom2dInt_GInter`.
//!
//! The results land in the inherited `IntRes2d_Intersection` through
//! `Append(Intersect, ...)` / `Append(Point)`, which is why `IsDone` only turns
//! true once a pair produced something (`IntRes2d_Intersection.cxx:188-384`).

use occt_core::elib::clib2d;
use occt_core::gp::{GpAx2d, GpDir2d, GpLin2d, GpPnt2d};
use occt_core::intres2d::{
    IntRes2dDomain, IntRes2dIntersection, IntRes2dIntersectionPoint, IntRes2dTransition,
};
use occt_core::precision::{INFINITE, PCONFUSION};
use occt_math::bissec_newton::MathBissecNewton;

use crate::curve::Curve2d;
use crate::geom2d_int::Geom2dIntGInter;
use crate::line::Geom2dLine;

use std::sync::Arc;

use super::function_inter::BisectorFunctionInter;
use super::{BisectorBisec, BisectorBisecCC};

/// `class Bisector_Inter : public IntRes2d_Intersection`
/// (`Bisector_Inter.hxx:32-90`).
#[derive(Clone, Debug, Default)]
pub struct BisectorInter {
    /// The `IntRes2d_Intersection` base.
    base: IntRes2dIntersection,
}

/// `occ::handle<Geom2d_Curve> SBis1[IB1]` (`cxx:86-88`): either the bisector
/// itself or the `Geom2d_Line` built by `ConstructSegment`.
#[derive(Clone)]
enum PartCurve<'a> {
    /// `SBis1[IB1] = Bis1` (`cxx:146`).
    Bisector(&'a dyn Curve2d),
    /// `SBis1[IB1] = ConstructSegment(...)` (`cxx:141`).
    Segment(Geom2dLine),
}

impl PartCurve<'_> {
    fn as_curve(&self) -> &dyn Curve2d {
        match self {
            Self::Bisector(c) => *c,
            Self::Segment(l) => l,
        }
    }
}

/// `Bisector_Curve::IsExtendAtStart` / `NbIntervals` / ... for a curve OCCT has
/// already `down_cast` to `Bisector_Curve` (`cxx:83-84`).
fn bisector_nb_intervals(c: &dyn Curve2d) -> i32 {
    if let Some(b) = c.bisector_bisec_ana() {
        return b.nb_intervals();
    }
    if let Some(b) = c.bisector_bisec_cc() {
        return b.nb_intervals();
    }
    if let Some(b) = c.bisector_bisec_pc() {
        return b.nb_intervals();
    }
    unreachable!("down_cast<Bisector_Curve> failed")
}

fn bisector_interval_first(c: &dyn Curve2d, index: i32) -> f64 {
    let r = if let Some(b) = c.bisector_bisec_ana() {
        b.interval_first(index)
    } else if let Some(b) = c.bisector_bisec_cc() {
        b.interval_first(index)
    } else if let Some(b) = c.bisector_bisec_pc() {
        b.interval_first(index)
    } else {
        unreachable!("down_cast<Bisector_Curve> failed")
    };
    r.expect("Bisector_Curve::IntervalFirst: Standard_OutOfRange")
}

fn bisector_interval_last(c: &dyn Curve2d, index: i32) -> f64 {
    let r = if let Some(b) = c.bisector_bisec_ana() {
        b.interval_last(index)
    } else if let Some(b) = c.bisector_bisec_cc() {
        b.interval_last(index)
    } else if let Some(b) = c.bisector_bisec_pc() {
        b.interval_last(index)
    } else {
        unreachable!("down_cast<Bisector_Curve> failed")
    };
    r.expect("Bisector_Curve::IntervalLast: Standard_OutOfRange")
}

fn bisector_is_extend_at_start(c: &dyn Curve2d) -> bool {
    if let Some(b) = c.bisector_bisec_ana() {
        return b.is_extend_at_start();
    }
    if let Some(b) = c.bisector_bisec_cc() {
        return b.is_extend_at_start();
    }
    if let Some(b) = c.bisector_bisec_pc() {
        return b.is_extend_at_start();
    }
    unreachable!("down_cast<Bisector_Curve> failed")
}

fn bisector_is_extend_at_end(c: &dyn Curve2d) -> bool {
    if let Some(b) = c.bisector_bisec_ana() {
        return b.is_extend_at_end();
    }
    if let Some(b) = c.bisector_bisec_cc() {
        return b.is_extend_at_end();
    }
    if let Some(b) = c.bisector_bisec_pc() {
        return b.is_extend_at_end();
    }
    unreachable!("down_cast<Bisector_Curve> failed")
}

/// `ConstructSegment(PMin, PMax, UMin, UMax)` (`cxx:61-72`): the line through
/// `PMin` with direction `PMax - PMin`, shifted back so that its parameter
/// `UMin` evaluates to `PMin`.
fn construct_segment(pmin: &GpPnt2d, pmax: &GpPnt2d, umin: f64) -> Geom2dLine {
    let dir = GpDir2d::new(pmax.x() - pmin.x(), pmax.y() - pmin.y()).unwrap_or_default();
    let loc = GpPnt2d::new(
        pmin.x() - umin * dir.x(),
        pmin.y() - umin * dir.y(),
    );
    Geom2dLine::new(GpAx2d::new(loc, dir))
}

impl BisectorInter {
    /// `Bisector_Inter()` (`cxx:44`): the defaulted constructor.
    pub fn new() -> Self {
        Self { base: IntRes2dIntersection::new() }
    }

    /// `Bisector_Inter(C1, D1, C2, D2, TolConf, Tol, ComunElement)`
    /// (`cxx:48-58`).
    #[allow(clippy::too_many_arguments)]
    pub fn with_bisectors(
        c1: &BisectorBisec,
        d1: &IntRes2dDomain,
        c2: &BisectorBisec,
        d2: &IntRes2dDomain,
        tol_conf: f64,
        tol: f64,
        comun_element: bool,
    ) -> Self {
        let mut inter = Self::new();
        inter.perform(c1, d1, c2, d2, tol_conf, tol, comun_element);
        inter
    }

    /// The inherited `IntRes2d_Intersection` result.
    pub fn result(&self) -> &IntRes2dIntersection {
        &self.base
    }

    /// `IntRes2d_Intersection::IsDone`.
    pub fn is_done(&self) -> bool {
        self.base.is_done()
    }

    /// `IntRes2d_Intersection::IsEmpty`; panics where OCCT raises
    /// `StdFail_NotDone`.
    pub fn is_empty(&self) -> bool {
        self.base.is_empty()
    }

    /// `IntRes2d_Intersection::NbPoints`.
    pub fn nb_points(&self) -> usize {
        self.base.nb_points()
    }

    /// `IntRes2d_Intersection::Point`.
    pub fn point(&self, n: usize) -> IntRes2dIntersectionPoint {
        *self.base.point(n)
    }

    /// `IntRes2d_Intersection::NbSegments`.
    pub fn nb_segments(&self) -> usize {
        self.base.nb_segments()
    }

    /// `IntRes2d_Intersection::Segment`.
    pub fn segment(&self, n: usize) -> occt_core::intres2d::IntRes2dIntersectionSegment {
        *self.base.segment(n)
    }

    /// `Perform(C1, D1, C2, D2, TolConf, Tol, ComunElement)` (`cxx:76-222`).
    #[allow(clippy::too_many_arguments)]
    pub fn perform(
        &mut self,
        c1: &BisectorBisec,
        d1: &IntRes2dDomain,
        c2: &BisectorBisec,
        d2: &IntRes2dDomain,
        tol_conf: f64,
        tol: f64,
        comun_element: bool,
    ) {
        // `cxx:83-84`: `down_cast<Bisector_Curve>(C1.Value()->BasisCurve())`.
        let bis1: &Arc<dyn Curve2d> = c1
            .value()
            .expect("Bisector_Inter: null thebisector (C1)")
            .basis_curve();
        let bis2: &Arc<dyn Curve2d> = c2
            .value()
            .expect("Bisector_Inter: null thebisector (C2)")
            .basis_curve();
        let bis1: &dyn Curve2d = &**bis1;
        let bis2: &dyn Curve2d = &**bis2;
        if !is_bisector_curve(bis1) || !is_bisector_curve(bis2) {
            panic!("Bisector_Inter: down_cast<Bisector_Curve> failed");
        }

        // `cxx:94-114`: the bounds of D1 (RealFirst / RealLast when unbounded).
        let (min_domain1, max_domain1) = domain_bounds(d1);

        // `cxx:119-151`: cut the first bisector by its continuity intervals.
        let nb_intervals1 = bisector_nb_intervals(bis1);
        // `SD1` / `SBis1` are 1-based in OCCT (`cxx:86-88`).
        let mut parts1: Vec<Option<(PartCurve<'_>, IntRes2dDomain)>> =
            vec![None; nb_intervals1 as usize + 1];
        let mut nb1 = 0usize;
        for ib1 in 1..=nb_intervals1 {
            let mut umin = bisector_interval_first(bis1, ib1);
            let mut umax = bisector_interval_last(bis1, ib1);
            if umax > min_domain1 && umin < max_domain1 {
                umin = umin.max(min_domain1);
                umax = umax.min(max_domain1);
                let pmin = bis1.d0(umin);
                let pmax = bis1.d0(umax);
                // `cxx:133`: `SD1[IB1].SetValues(PMin, UMin, D1.FirstTolerance(),
                // PMax, UMax, D1.LastTolerance())`.
                let sd = IntRes2dDomain::bounded(
                    &pmin,
                    umin,
                    d1.first_tolerance(),
                    &pmax,
                    umax,
                    d1.last_tolerance(),
                );
                let seg = if (ib1 == 1 && bisector_is_extend_at_start(bis1))
                    || (ib1 == nb_intervals1 && bisector_is_extend_at_end(bis1))
                {
                    // `cxx:136-142`.
                    PartCurve::Segment(construct_segment(&pmin, &pmax, umin))
                } else {
                    // `cxx:144-147`.
                    PartCurve::Bisector(bis1)
                };
                parts1[ib1 as usize] = Some((seg, sd));
                nb1 += 1;
            }
        }

        // `cxx:157-177`: the bounds of D2.
        let (min_domain2, max_domain2) = domain_bounds(d2);

        // `cxx:179-206`: cut the second bisector by its continuity intervals.
        let nb_intervals2 = bisector_nb_intervals(bis2);
        let mut parts2: Vec<Option<(PartCurve<'_>, IntRes2dDomain)>> =
            vec![None; nb_intervals2 as usize + 1];
        let mut nb2 = 0usize;
        for ib2 in 1..=nb_intervals2 {
            let mut umin = bisector_interval_first(bis2, ib2);
            let mut umax = bisector_interval_last(bis2, ib2);
            if umax > min_domain2 && umin < max_domain2 {
                umin = umin.max(min_domain2);
                umax = umax.min(max_domain2);
                let pmin = bis2.d0(umin);
                let pmax = bis2.d0(umax);
                let sd = IntRes2dDomain::bounded(
                    &pmin,
                    umin,
                    d2.first_tolerance(),
                    &pmax,
                    umax,
                    d2.last_tolerance(),
                );
                // `cxx:192-193` reads `Bis1->NbIntervals()` in the
                // `IsExtendAtEnd` test; the OCCT source uses `Bis1`, not
                // `Bis2`, and the port keeps that exact predicate.
                let seg = if (ib2 == 1 && bisector_is_extend_at_start(bis2))
                    || (ib2 == nb_intervals1 && bisector_is_extend_at_end(bis2))
                {
                    PartCurve::Segment(construct_segment(&pmin, &pmax, umin))
                } else {
                    PartCurve::Bisector(bis2)
                };
                parts2[ib2 as usize] = Some((seg, sd));
                nb2 += 1;
            }
        }

        // `cxx:209-217`: OCCT walks `IB1 in 1..=NB1` and `IB2 in 1..=NB2`, i.e.
        // it indexes the slots by the *count* of accepted intervals, not by the
        // interval indices. A rejected leading interval therefore leaves that
        // slot unset (`SBis1[IB1]` is the null handle, `SD1[IB1]` the default
        // infinite domain) and OCCT dereferences a null handle there. The port
        // skips such an unset slot.
        for ib1 in 1..=nb1 {
            let Some((part1, sd1)) = parts1[ib1].as_ref() else {
                continue;
            };
            for ib2 in 1..=nb2 {
                let Some((part2, sd2)) = parts2[ib2].as_ref() else {
                    continue;
                };
                self.single_perform(
                    part1.as_curve(),
                    sd1,
                    part2.as_curve(),
                    sd2,
                    tol_conf,
                    tol,
                    comun_element,
                );
            }
        }
    }

    /// `SinglePerform(C1, D1, C2, D2, TolConf, Tol, ComunElement)`
    /// (`cxx:226-324`).
    #[allow(clippy::too_many_arguments)]
    fn single_perform(
        &mut self,
        cbis1: &dyn Curve2d,
        d1: &IntRes2dDomain,
        cbis2: &dyn Curve2d,
        d2: &IntRes2dDomain,
        tol_conf: f64,
        tol: f64,
        comun_element: bool,
    ) {
        // `cxx:238-239`.
        let type1_is_ana = cbis1.bisector_bisec_ana().is_some();
        let type2_is_ana = cbis2.bisector_bisec_ana().is_some();

        if type1_is_ana || type2_is_ana {
            // `cxx:240-258`: replace an analytic bisector by its basis curve.
            let c2bis1: &dyn Curve2d = if type1_is_ana {
                match cbis1.bisector_bisec_ana().unwrap().geom2d_curve() {
                    Some(c) => &**c,
                    None => panic!("Bisector_BisecAna: null thebisector"),
                }
            } else {
                cbis1
            };
            let c2bis2: &dyn Curve2d = if type2_is_ana {
                match cbis2.bisector_bisec_ana().unwrap().geom2d_curve() {
                    Some(c) => &**c,
                    None => panic!("Bisector_BisecAna: null thebisector"),
                }
            } else {
                cbis2
            };
            let type1_is_line = c2bis1.is_geom2d_line();
            let type2_is_line = c2bis2.is_geom2d_line();
            // `cxx:260-268`.
            if type1_is_line && !type2_is_line {
                let l1 = c2bis1.gp_lin2d().expect("Geom2d_Line");
                self.test_bound(&l1, d1, c2bis2, d2, tol_conf, false);
            } else if type2_is_line && !type1_is_line {
                let l1 = c2bis2.gp_lin2d().expect("Geom2d_Line");
                self.test_bound(&l1, d2, c2bis1, d1, tol_conf, true);
            }
            // `cxx:269-277`.
            let mut intersect = Geom2dIntGInter::new();
            intersect.perform(c2bis1, d1, c2bis2, d2, tol_conf, tol);
            self.append_intersection(&intersect, d1, d2);
        } else if cbis1.bisector_bisec_pc().is_some() || cbis2.bisector_bisec_pc().is_some() {
            // `cxx:279-289`.
            let mut intersect = Geom2dIntGInter::new();
            intersect.perform(cbis1, d1, cbis2, d2, tol_conf, tol);
            self.append_intersection(&intersect, d1, d2);
        } else if comun_element
            && cbis1.bisector_bisec_cc().is_some()
            && cbis2.bisector_bisec_cc().is_some()
        {
            // `cxx:290-298`.
            let bis1 = cbis1.bisector_bisec_cc().unwrap();
            let bis2 = cbis2.bisector_bisec_cc().unwrap();
            self.neighbour_perform(bis1, d1, bis2, d2, tol);
        } else {
            // `cxx:299-322`: one of the two is a segment; test the extremities
            // of the other one against it before the generic intersection.
            let type1_is_line = cbis1.is_geom2d_line();
            let type2_is_line = cbis2.is_geom2d_line();
            if type1_is_line && !type2_is_line {
                let l1 = cbis1.gp_lin2d().expect("Geom2d_Line");
                self.test_bound(&l1, d1, cbis2, d2, tol_conf, false);
            } else if type2_is_line && !type1_is_line {
                let l1 = cbis2.gp_lin2d().expect("Geom2d_Line");
                self.test_bound(&l1, d2, cbis1, d1, tol_conf, true);
            }
            let mut intersect = Geom2dIntGInter::new();
            intersect.perform(cbis1, d1, cbis2, d2, tol_conf, tol);
            self.append_intersection(&intersect, d1, d2);
        }
    }

    /// `Append(Intersect, D1.FirstParameter(), D1.LastParameter(),
    /// D2.FirstParameter(), D2.LastParameter())` (`cxx:273-276`).
    fn append_intersection(
        &mut self,
        intersect: &Geom2dIntGInter,
        d1: &IntRes2dDomain,
        d2: &IntRes2dDomain,
    ) {
        self.base.append_intersection(
            intersect.result(),
            d1.first_parameter(),
            d1.last_parameter(),
            d2.first_parameter(),
            d2.last_parameter(),
        );
    }

    /// `NeighbourPerform(Bis1, D1, Bis2, D2, Tol)` (`cxx:338-391`).
    fn neighbour_perform(
        &mut self,
        bis1: &BisectorBisecCC,
        d1: &IntRes2dDomain,
        bis2: &BisectorBisecCC,
        d2: &IntRes2dDomain,
        tol: f64,
    ) {
        const EPS: f64 = PCONFUSION;

        // `cxx:350-351`: change the guide line of Bis2, then the guide is
        // Bis2's own second curve.
        let bis_temp = bis2.change_guide();
        let guide = bis2
            .curve(2)
            .expect("Bisector_BisecCC::Curve: Standard_OutOfRange");

        // `cxx:355-358`: the returned points are not used in OCCT (kept for the
        // debugger); only `UMax` / `UMin` matter.
        let mut u1 = 0.0;
        let mut u_max = 0.0;
        let mut dist = 0.0;
        let _p2s = bis2.value_and_dist(d2.first_parameter(), &mut u1, &mut u_max, &mut dist);
        let mut u_min = 0.0;
        let _p2e = bis2.value_and_dist(d2.last_parameter(), &mut u1, &mut u_min, &mut dist);

        // `cxx:360-361`: the domain of intersection on the guideline.
        let umin = d1.first_parameter().max(u_min);
        let umax = d1.last_parameter().min(u_max);

        // `cxx:363`.
        self.base.done = true;

        // `cxx:365-368`.
        if umin - EPS > umax + EPS {
            return;
        }

        // `cxx:371`: solution of F = 0.
        let mut fint = BisectorFunctionInter::with_curves(&**guide, bis1, &bis_temp);
        let mut a_solution = MathBissecNewton::new(tol);
        a_solution.perform(&mut fint, umin, umax, 20);
        if !a_solution.is_done() {
            // `cxx:380-383`.
            return;
        }
        let u_sol = a_solution.root();

        // `cxx:386-389`.
        let mut u1 = 0.0;
        let mut u2 = 0.0;
        let mut dist = 0.0;
        let p_sol = bis_temp.value_and_dist(u_sol, &mut u1, &mut u2, &mut dist);

        let trans1 = IntRes2dTransition::new();
        let trans2 = IntRes2dTransition::new();
        let point_inter_sol =
            IntRes2dIntersectionPoint::with_transitions(&p_sol, u_sol, u2, &trans1, &trans2, false);
        // `cxx:390`: `Append(PointInterSol)` -> `IntRes2d_Intersection::Insert`.
        self.base.insert(&point_inter_sol);
    }

    /// `TestBound(Bis1, D1, Bis2, D2, TolConf, Reverse)` (`cxx:398-462`): test
    /// whether the extremities of `Bis2` lie on the line `Bis1`.
    #[allow(clippy::too_many_arguments)]
    fn test_bound(
        &mut self,
        l1: &GpLin2d,
        d1: &IntRes2dDomain,
        bis2: &dyn Curve2d,
        d2: &IntRes2dDomain,
        tol_conf: f64,
        reverse: bool,
    ) {
        let trans1 = IntRes2dTransition::new();
        let trans2 = IntRes2dTransition::new();
        let mut point_inter_sol = IntRes2dIntersectionPoint::new();

        let mut pf = bis2.d0(d2.first_parameter());
        let mut pl = bis2.d0(d2.last_parameter());
        // `cxx:411-418`: `Tol = TolConf` since OCC616, plus the `BisecAlgo`
        // flag that decides whether the point is projected onto the line.
        let tol = tol_conf;
        let bisec_algo = bis2.bisector_bisec_cc().is_some();

        // `cxx:425-443`.
        if l1.distance(&pf) < tol {
            let u1 = clib2d::parameter_lin2d(l1, &pf);
            if d1.first_parameter() - d1.first_tolerance() < u1
                && d1.last_parameter() + d1.last_tolerance() > u1
            {
                if bisec_algo {
                    pf = clib2d::line_value_ax2d(u1, l1.position());
                }
                point_inter_sol.set_values(
                    &pf,
                    u1,
                    d2.first_parameter(),
                    &trans1,
                    &trans2,
                    reverse,
                );
                self.base.insert(&point_inter_sol);
            }
        }

        // `cxx:446-461`.
        if l1.distance(&pl) < tol {
            let u1 = clib2d::parameter_lin2d(l1, &pl);
            if d1.first_parameter() - d1.first_tolerance() < u1
                && d1.last_parameter() + d1.last_tolerance() > u1
            {
                if bisec_algo {
                    pl = clib2d::line_value_ax2d(u1, l1.position());
                }
                point_inter_sol.set_values(
                    &pl,
                    u1,
                    d2.last_parameter(),
                    &trans1,
                    &trans2,
                    reverse,
                );
                self.base.insert(&point_inter_sol);
            }
        }
    }
}

/// `occ::down_cast<Bisector_Curve>(C).IsNull()` (`cxx:83-84`).
fn is_bisector_curve(c: &dyn Curve2d) -> bool {
    c.bisector_bisec_ana().is_some()
        || c.bisector_bisec_cc().is_some()
        || c.bisector_bisec_pc().is_some()
}

/// `D1.HasFirstPoint() ? D1.FirstParameter() : RealFirst()` and its `RealLast`
/// twin (`cxx:100-114`).
fn domain_bounds(d: &IntRes2dDomain) -> (f64, f64) {
    let min = if d.has_first_point() {
        d.first_parameter()
    } else {
        -INFINITE
    };
    let max = if d.has_last_point() {
        d.last_parameter()
    } else {
        INFINITE
    };
    (min, max)
}
