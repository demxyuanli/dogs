//! `Geom2dConvert_ApproxArcsSegments` — approximation of a 2D adaptor curve by
//! arcs of circles and line segments.
//!
//! Source: `ModelingData/TKGeomBase/Geom2dConvert/Geom2dConvert_ApproxArcsSegments.cxx:37-919`.
//!
//! The OCCT class stores `const Adaptor2d_Curve2d& myCurve`; the port stores the
//! same adaptor as `&dyn Curve2d` (the crate's `Adaptor2d_Curve2d` equivalent,
//! see `gcpnts_quasi_uniform_deflection.rs`). The constructor and the
//! single-piece builders live here; the free-form chain (`makeFreeform` and the
//! inflection-point helpers) is in [`make_freeform`].
//!
//! Fidelity notes:
//! - `aLastPnt` (`cxx:403`, `:410`, `:421`) is dead in OCCT: it is assigned in
//!   both sampling loops and never read. It is not carried here.
//! - A null `occ::handle<Geom2d_Curve>` becomes `None`: OCCT appends nulls to
//!   `mySeqCurves` on one path (`cxx:563`, when the line is not built but the
//!   distance is below `Precision::Confusion()`), and `GetResult()` hands them
//!   back, so the result vector keeps `Option` slots at the same positions.

use std::sync::Arc;

use occt_core::elib::clib2d::{circle_value_ax22d, parameter_circ2d, parameter_lin2d};
use occt_core::gp::{GpAx2d, GpAx22d, GpCirc2d, GpDir2d, GpPnt2d, GpVec2d};
use occt_core::precision::{ANGULAR, CONFUSION, PCONFUSION, REAL_SMALL};

use crate::circle::Geom2dCircle;
use crate::curve::Curve2d;
use crate::gce2d::make_segment;
use crate::gce_make_circ2d::gce_make_circ2d_3p;
use crate::geom2d_convert_ppoint::PPoint;
use crate::line::Geom2dLine;
use crate::trimmed::Geom2dTrimmedCurve;

mod make_freeform;

/// `MAXPOINTS` (`cxx:34`), shared with the free-form chain.
pub(super) const MAX_POINTS: usize = 100;
/// `MyCurvatureTolerance` (`cxx:35`), shared with the free-form chain.
pub(super) const MY_CURVATURE_TOLERANCE: f64 = 0.0001;

/// `Geom2dConvert_ApproxArcsSegments::Status`
/// (`Geom2dConvert_ApproxArcsSegments.hxx`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Status {
    NotDone,
    Ok,
    Error,
}

/// `GC_MakeArcOfCircle2d(P1, P2, P3).Value()`
/// (`GC_MakeArcOfCircle2d.cxx:31-47`): the circle through the three points,
/// trimmed from `ElCLib::Parameter(C, P1)` to `ElCLib::Parameter(C, P3)` with
/// `Sense = true`.
fn gc_make_arc_of_circle2d(
    p1: &GpPnt2d,
    p2: &GpPnt2d,
    p3: &GpPnt2d,
) -> Option<Geom2dTrimmedCurve> {
    let circ = gce_make_circ2d_3p(p1, p2, p3).ok()?;
    let alpha1 = parameter_circ2d(&circ, p1);
    let alpha2 = parameter_circ2d(&circ, p3);
    let basis: Arc<dyn Curve2d> = Arc::new(Geom2dCircle::new(circ));
    Some(Geom2dTrimmedCurve::new_sense(
        basis, alpha1, alpha2, true, true,
    ))
}

/// `Geom2dConvert_ApproxArcsSegments` (`cxx:55-117` and members).
pub struct ApproxArcsSegments<'a> {
    pub(super) curve: &'a dyn Curve2d,
    pub(super) tolerance: f64,
    pub(super) angle_tolerance: f64,
    pub(super) seq_params: Vec<PPoint>,
    pub(super) seq_curves: Vec<Option<Arc<dyn Curve2d>>>,
    pub(super) ext: [PPoint; 2],
    pub(super) status: Status,
}

impl<'a> ApproxArcsSegments<'a> {
    /// `Geom2dConvert_ApproxArcsSegments(theCurve, theTolerance, theAngleTol)`
    /// (`cxx:55-115`).
    pub fn new(curve: &'a dyn Curve2d, tolerance: f64, angle_tol: f64) -> Self {
        let mut me = Self {
            curve,
            tolerance,
            angle_tolerance: angle_tol,
            seq_params: Vec::new(),
            seq_curves: Vec::new(),
            ext: [
                PPoint::on_curve(curve.first_parameter(), curve),
                PPoint::on_curve(curve.last_parameter(), curve),
            ],
            status: Status::NotDone,
        };

        if curve.is_line() {
            // `cxx:71-82`: a single line segment.
            let dist = me.ext[0].dist(&me.ext[1]);
            if dist > CONFUSION {
                if let Ok(dir) = GpDir2d::new(
                    me.ext[1].point().x() - me.ext[0].point().x(),
                    me.ext[1].point().y() - me.ext[0].point().y(),
                ) {
                    let line = Geom2dLine::new(GpAx2d::new(me.ext[0].point(), dir));
                    let basis: Arc<dyn Curve2d> = Arc::new(line);
                    let seg = Geom2dTrimmedCurve::new(basis, 0.0, dist);
                    me.seq_curves.push(Some(Arc::new(seg)));
                    me.status = Status::Ok;
                }
            }
        } else if curve.gp_circ2d().is_some() {
            // `cxx:83-97`: a couple of arcs of equal size.
            let mid_param = 0.5 * (me.ext[0].parameter() + me.ext[1].parameter());
            let mid = PPoint::on_curve(mid_param, curve);
            let first = me.ext[0];
            let last = me.ext[1];
            if let Some(arc) = me.make_circle(&first, &mid) {
                me.seq_curves.push(Some(arc));
                if let Some(arc2) = me.make_circle(&mid, &last) {
                    me.seq_curves.push(Some(arc2));
                }
            }
        } else {
            me.make_freeform();
        }

        // `cxx:102-116`.
        if me.status == Status::NotDone {
            me.status = if me.seq_curves.is_empty() {
                Status::Error
            } else {
                Status::Ok
            };
        }
        me
    }

    /// `Status()` / `IsDone()` (`cxx:807-820`).
    pub fn status(&self) -> Status {
        self.status
    }

    /// `GetResult()` (`cxx:787-790`): the approximated curves, in order of
    /// appearance along the adaptor curve.
    pub fn result(&self) -> &[Option<Arc<dyn Curve2d>>] {
        &self.seq_curves
    }

    /// `makeCircle` (`cxx:119-137`).
    fn make_circle(&self, first: &PPoint, last: &PPoint) -> Option<Arc<dyn Curve2d>> {
        let para_m = (first.parameter() + last.parameter()) * 0.5;
        let point_m = self.curve.d0(para_m);
        let arc = gc_make_arc_of_circle2d(&first.point(), &point_m, &last.point())?;
        Some(Arc::new(arc))
    }

    /// `makeArc` (`cxx:140-217`). `the_curve` is only assigned on success,
    /// which the `true` return marks.
    pub(super) fn make_arc(
        &self,
        param1: &PPoint,
        param2: &mut PPoint,
        is_first: bool,
        the_curve: &mut Option<Geom2dTrimmedCurve>,
    ) -> bool {
        let a_p1 = param1.point().coord;
        let a_p2 = param2.point().coord;
        let a_vec = if is_first {
            param1.d1().xy()
        } else {
            param1.d1().xy().reversed()
        };

        // `cxx:150-164`: the sense, from the sign of `aVec ^ aDelta`.
        let a_delta = a_p2.subtracted(&a_p1);
        let mut a_sense = a_vec.crossed(&a_delta);
        if a_sense > ANGULAR {
            a_sense = 1.0;
        } else if a_sense < -ANGULAR {
            a_sense = -1.0;
        } else {
            return false;
        }

        // `cxx:167-178`: the centre, from the two tangent directions.
        let a_middle = a_p1.added(&a_p2).multiplied_scalar(0.5);
        let prod_p1v = a_p1.dot(&a_vec);
        let prod_de_m = a_delta.dot(&a_middle);
        let vprod_vd = a_vec.crossed(&a_delta);
        let a_resolution = REAL_SMALL;
        if vprod_vd < -a_resolution || vprod_vd > a_resolution {
            let center = GpPnt2d::new(
                (prod_p1v * a_delta.y() - prod_de_m * a_vec.y()) / vprod_vd,
                (prod_de_m * a_vec.x() - prod_p1v * a_delta.x()) / vprod_vd,
            );
            let a_rad =
                (center.distance(&param1.point()) + center.distance(&param2.point())) * 0.5;
            // `gp_Ax22d(aCenter, gp::DX2d(), gp::DY2d())`.
            let ax22d = GpAx22d {
                point: center,
                ..GpAx22d::standard()
            };
            let a_cir = GpCirc2d::new(ax22d, a_rad);
            let mut circ = Geom2dCircle::new(a_cir);

            // `cxx:184-195`: the arc parameters, swapped when the arc runs
            // backwards.
            let (mut alpha1, mut alpha2) = if is_first {
                (
                    parameter_circ2d(&a_cir, &param1.point()),
                    parameter_circ2d(&a_cir, &param2.point()),
                )
            } else {
                a_sense = -a_sense;
                (
                    parameter_circ2d(&a_cir, &param2.point()),
                    parameter_circ2d(&a_cir, &param1.point()),
                )
            };

            if (alpha1 - alpha2).abs() < 1e-100 {
                return false;
            }

            // `cxx:203-208`.
            if a_sense < 0.0 {
                alpha1 = circ.reversed_parameter(alpha1);
                alpha2 = circ.reversed_parameter(alpha2);
                circ.reverse();
            }
            let basis: Arc<dyn Curve2d> = Arc::new(circ);
            *the_curve = Some(Geom2dTrimmedCurve::new(basis, alpha1, alpha2));

            // `cxx:211-213`: the tangent at the second point.
            let a_rad_v = param2.point().xy().subtracted(&center.xy());
            param2.set_d1(GpVec2d::new(-a_rad_v.y() * a_sense, a_rad_v.x() * a_sense));
            return true;
        }
        false
    }

    /// `makeLine` (`cxx:222-286`).
    pub(super) fn make_line(
        &self,
        first: &mut PPoint,
        last: &mut PPoint,
        is_check: bool,
    ) -> Option<Arc<dyn Curve2d>> {
        let slope = last.point().xy().subtracted(&first.point().xy());
        if slope.square_modulus() < REAL_SMALL {
            return None;
        }
        let dir_line = GpDir2d::new(slope.x(), slope.y()).ok()?;

        if is_check {
            if first.d1().square_magnitude() < REAL_SMALL
                || last.d1().square_magnitude() < REAL_SMALL
            {
                return None;
            }

            // `cxx:246-256`: G1 is only required away from the curve ends.
            let mut abs_angle = [0.0f64, 0.0];
            if *first != self.ext[0] {
                let d = GpDir2d::from_xy(&first.d1().xy()).ok()?;
                abs_angle[0] = dir_line.angle(&d).abs();
            }
            if *last != self.ext[1] {
                let d = GpDir2d::from_xy(&last.d1().xy()).ok()?;
                abs_angle[1] = dir_line.angle(&d).abs();
            }

            let cont_tolerance = self.angle_tolerance.max(0.01);
            if abs_angle[0] > cont_tolerance || abs_angle[1] > cont_tolerance {
                return None;
            }
        }

        // `cxx:270-284`: `GC_MakeSegment2d`.
        let seg = make_segment(&first.point(), &last.point());
        if self.check_curve(&seg, first.parameter(), last.parameter()) {
            let new_d1 = last.point().xy().subtracted(&first.point().xy());
            first.set_d1(GpVec2d::from_xy(new_d1));
            last.set_d1(GpVec2d::from_xy(new_d1));
            Some(Arc::new(seg))
        } else {
            None
        }
    }

    /// `checkCurve` (`cxx:744-822`): the maximum deflection between the adaptor
    /// curve and the line/circle, sampled at `mySeqParams` when available and
    /// at `MAXPOINTS` uniform parameters otherwise.
    pub(super) fn check_curve(
        &self,
        a_curve: &dyn Curve2d,
        first_param: f64,
        last_param: f64,
    ) -> bool {
        let is_uniform_done = !self.seq_params.is_empty();
        let nb_pnts = if is_uniform_done {
            self.seq_params.len() as i32
        } else {
            MAX_POINTS as i32
        };
        let param_step = (last_param - first_param) / MAX_POINTS as f64;

        // `cxx:757-783`: the trimmed curve is replaced by its basis, then
        // down-cast to `Geom2d_Line` / `Geom2d_Circle`.
        let curve1: &dyn Curve2d = a_curve.trimmed_basis().unwrap_or(a_curve);
        let (is_line, a_lin2d, a_circ2d) = if curve1.is_geom2d_line() {
            match curve1.gp_lin2d() {
                Some(l) => (true, l, GpCirc2d::default()),
                None => return false,
            }
        } else if curve1.is_geom2d_circle() {
            match curve1.gp_circ2d() {
                Some(c) => (false, Default::default(), c),
                None => return false,
            }
        } else {
            return false;
        };

        let mut is_less = true;
        let mut i = 1i32;
        while i <= nb_pnts && is_less {
            let a_param = if is_uniform_done {
                self.seq_params[(i - 1) as usize].parameter()
            } else {
                first_param + i as f64 * param_step
            };
            i += 1;
            if a_param < (first_param - PCONFUSION) || a_param > (last_param + PCONFUSION) {
                continue;
            }

            let point_adaptor = self.curve.d0(a_param);
            let point_proj = if is_line {
                let t = parameter_lin2d(&a_lin2d, &point_adaptor);
                let loc = a_lin2d.location();
                GpPnt2d::new(
                    loc.x() + t * a_lin2d.direction().x,
                    loc.y() + t * a_lin2d.direction().y,
                )
            } else {
                let t = parameter_circ2d(&a_circ2d, &point_adaptor);
                circle_value_ax22d(t, &a_circ2d.pos, a_circ2d.radius)
            };

            is_less = point_proj.distance(&point_adaptor) < self.tolerance + CONFUSION;
        }
        is_less
    }
}
