//! `Bisector_Bisec` (`Bisector_Bisec.cxx:49-681`, `.hxx:40-100`): the dispatcher
//! that builds the trimmed bisector between two curves, a curve and a point, or
//! two points.
//!
//! The evaluator is picked from the dynamic type of the arguments, exactly as in
//! OCCT : analytic pairs (`Geom2d_Circle` / `Geom2d_Line`) go to
//! `Bisector_BisecAna`, the curve/curve remainder goes to `Bisector_BisecCC`, and
//! the curve/point remainder to `Bisector_BisecPC`. The same `Perform` also
//! covers the degenerate half-straight cases that OCCT special-cases in place.

use std::sync::Arc;

use occt_core::gp::{GpAx2d, GpDir2d, GpPnt2d, GpVec2d};
use occt_core::kernel::geomabs::JoinType;
use occt_core::precision::{ANGULAR, CONFUSION, INFINITE, PCONFUSION, REAL_EPSILON, RESOLUTION};

use crate::curve::Curve2d;
use crate::line::Geom2dLine;
use crate::trimmed::Geom2dTrimmedCurve;

use super::{is_convex, BisectorBisecAna, BisectorBisecAnaError, BisectorBisecCC, BisectorBisecPC};

/// Error of the `Bisector_Bisec` dispatcher.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BisectorBisecError {
    /// The analytic/algorithmic construction raised `StdFail_NotDone`.
    NotDone,
    /// `thebisector` is the null handle (`Value` / `ChangeValue`).
    NullBisector,
}

impl std::fmt::Display for BisectorBisecError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = match self {
            Self::NotDone => "StdFail_NotDone",
            Self::NullBisector => "null thebisector",
        };
        write!(f, "{s}")
    }
}

impl std::error::Error for BisectorBisecError {}

impl From<BisectorBisecAnaError> for BisectorBisecError {
    fn from(_: BisectorBisecAnaError) -> Self {
        Self::NotDone
    }
}

/// `class Bisector_Bisec` (`Bisector_Bisec.hxx:40-100`).
#[derive(Clone, Default)]
pub struct BisectorBisec {
    /// `occ::handle<Geom2d_TrimmedCurve> thebisector` (`.hxx`); the handle is
    /// null after the default constructor.
    thebisector: Option<Geom2dTrimmedCurve>,
}

impl BisectorBisec {
    /// `Bisector_Bisec()` (`cxx:49`): the defaulted constructor, which leaves
    /// `thebisector` null.
    pub fn new() -> Self {
        Self { thebisector: None }
    }

    /// `Value()` (`cxx:567-570`).
    pub fn value(&self) -> Option<&Geom2dTrimmedCurve> {
        self.thebisector.as_ref()
    }

    /// `ChangeValue()` (`cxx:574-577`).
    pub fn change_value(&mut self) -> Result<&mut Geom2dTrimmedCurve, BisectorBisecError> {
        self.thebisector
            .as_mut()
            .ok_or(BisectorBisecError::NullBisector)
    }

    /// Replaces `thebisector`; used by the `Perform` overloads of this file and
    /// by the `BRepFill` callers that store a bisector by value.
    pub fn set_value(&mut self, bis: Geom2dTrimmedCurve) {
        self.thebisector = Some(bis);
    }

    /// `Perform(afirstcurve, asecondcurve, ...)` (`cxx:63-246`).
    #[allow(clippy::too_many_arguments)]
    pub fn perform_curve_curve(
        &mut self,
        afirstcurve: &dyn Curve2d,
        asecondcurve: &dyn Curve2d,
        apoint: &GpPnt2d,
        afirstvector: &GpVec2d,
        asecondvector: &GpVec2d,
        adirection: f64,
        ajointype: JoinType,
        tolerance: f64,
        oncurve: bool,
    ) -> Result<(), BisectorBisecError> {
        let mut first_circle_or_line = is_circle_or_line(afirstcurve);
        let mut second_circle_or_line = is_circle_or_line(asecondcurve);

        // `cxx:88-132`: a degree-1 two-pole B-spline whose poles coincide is
        // replaced by the segment between them, which makes it analytic.
        let mut first: Box<dyn Curve2d> = afirstcurve.clone_dyn();
        let mut second: Box<dyn Curve2d> = asecondcurve.clone_dyn();
        if let Some(segment) = collapsed_segment(afirstcurve) {
            first = segment;
            first_circle_or_line = true;
        }
        if let Some(segment) = collapsed_segment(asecondcurve) {
            second = segment;
            second_circle_or_line = true;
        }

        let mut u_first;
        let mut u_last;
        let bis: Arc<dyn Curve2d>;

        if first_circle_or_line && second_circle_or_line {
            // `cxx:137-153`: analytic bisector.
            let mut bis_ana = BisectorBisecAna::new();
            bis_ana.perform_curve_curve(
                &*first,
                &*second,
                apoint,
                afirstvector,
                asecondvector,
                adirection,
                ajointype,
                tolerance,
                oncurve,
            )?;
            u_first = bis_ana.parameter_of_start_point();
            u_last = bis_ana.parameter_of_end_point();
            bis = Arc::from(bis_ana.clone_dyn());
        } else {
            let fd = GpDir2d::from_vec2d(afirstvector).unwrap_or_default();
            let sd = GpDir2d::from_vec2d(asecondvector).unwrap_or_default();
            let is_line = oncurve && fd.dot(&sd) < (2.0 * ANGULAR).sqrt() - 1.0;
            if is_line {
                // `cxx:169-181`: half-straight.
                let bis_ana = half_straight(apoint, afirstvector, adirection, INFINITE);
                u_first = bis_ana.parameter_of_start_point();
                u_last = bis_ana.parameter_of_end_point();
                bis = Arc::from(bis_ana.clone_dyn());
            } else {
                // `cxx:187-240`: algorithmic bisector.
                let mut bis_cc = BisectorBisecCC::new();
                bis_cc.perform(&*second, &*first, adirection, adirection, apoint, INFINITE);
                if bis_cc.is_empty() {
                    // `cxx:194-232`: the algorithmic bisector is empty, build a
                    // false bisector through the origin.
                    let bis_ana = false_bisector(
                        apoint,
                        afirstvector,
                        asecondvector,
                        adirection,
                        INFINITE,
                    );
                    u_first = bis_ana.parameter_of_start_point();
                    u_last = bis_ana.parameter_of_end_point();
                    bis = Arc::from(bis_ana.clone_dyn());
                } else {
                    u_first = bis_cc.first_parameter();
                    u_last = bis_cc.last_parameter();
                    let raw: Arc<dyn Curve2d> = Arc::from(bis_cc.clone_dyn());
                    bis = replace_by_line_if_is_to_small(&raw, &mut u_first, &mut u_last);
                }
            }
        }

        u_first = u_first.max(bis.first_parameter());
        u_last = u_last.min(bis.last_parameter());
        self.thebisector = Some(Geom2dTrimmedCurve::new(bis, u_first, u_last));
        Ok(())
    }

    /// `Perform(afirstcurve, asecondpoint, ...)` (`cxx:260-392`): bisector
    /// between a curve and a point.
    #[allow(clippy::too_many_arguments)]
    pub fn perform_curve_point(
        &mut self,
        afirstcurve: &dyn Curve2d,
        asecondpoint: &GpPnt2d,
        apoint: &GpPnt2d,
        afirstvector: &GpVec2d,
        asecondvector: &GpVec2d,
        adirection: f64,
        tolerance: f64,
        oncurve: bool,
    ) -> Result<(), BisectorBisecError> {
        let mut u_first;
        let mut u_last;
        let bis: Arc<dyn Curve2d>;

        if is_circle_or_line(afirstcurve) {
            // `cxx:280-296`: analytic bisector.
            let mut bis_ana = BisectorBisecAna::new();
            bis_ana.perform_curve_point(
                afirstcurve,
                asecondpoint,
                apoint,
                afirstvector,
                asecondvector,
                adirection,
                tolerance,
                oncurve,
            );
            u_first = bis_ana.parameter_of_start_point();
            u_last = bis_ana.parameter_of_end_point();
            bis = Arc::from(bis_ana.clone_dyn());
        } else {
            let mut is_line = false;
            let mut rc = INFINITE;
            if oncurve {
                if is_convex(afirstcurve, adirection)
                    || is_max_rc(afirstcurve, afirstcurve.last_parameter(), &mut rc)
                {
                    is_line = true;
                }
            }
            if is_line {
                // `cxx:312-325`: half-right.
                let bis_ana = half_straight(apoint, afirstvector, adirection, rc);
                u_first = bis_ana.parameter_of_start_point();
                u_last = bis_ana.parameter_of_end_point();
                bis = Arc::from(bis_ana.clone_dyn());
            } else {
                // `cxx:330-380`: algorithmic bisector on the reversed curve.
                let reversed = afirstcurve.reversed();
                let mut bis_pc = BisectorBisecPC::new();
                bis_pc.perform(&*reversed, asecondpoint, -adirection, INFINITE);
                if bis_pc.is_empty() {
                    // `cxx:334-354`: false bisector.
                    let bis_ana = false_bisector(
                        apoint,
                        afirstvector,
                        asecondvector,
                        adirection,
                        rc,
                    );
                    u_first = bis_ana.parameter_of_start_point();
                    u_last = bis_ana.parameter_of_end_point();
                    bis = Arc::from(bis_ana.clone_dyn());
                } else {
                    u_first = bis_pc.parameter(apoint);
                    u_last = bis_pc.last_parameter();
                    let raw = if u_first >= u_last {
                        // `cxx:360-375`: the curve parameter is beyond the last
                        // one, extrapolate by a line along the bisector.
                        let v = vec_between(
                            &bis_pc.value(bis_pc.first_parameter()),
                            &bis_pc.value(u_last),
                        );
                        let line = line_through(apoint, &v, rc);
                        u_first = line.parameter_of_start_point();
                        u_last = line.parameter_of_end_point();
                        Arc::from(line.clone_dyn())
                    } else {
                        Arc::from(bis_pc.clone_dyn())
                    };
                    bis = raw;
                }
            }
        }

        if u_first < bis.first_parameter() {
            u_first = bis.first_parameter();
        }
        if u_last > bis.last_parameter() {
            u_last = bis.last_parameter();
        }
        self.thebisector = Some(Geom2dTrimmedCurve::new(bis, u_first, u_last));
        Ok(())
    }

    /// `Perform(afirstpoint, asecondcurve, ...)` (`cxx:405-528`): bisector
    /// between a point and a curve.
    #[allow(clippy::too_many_arguments)]
    pub fn perform_point_curve(
        &mut self,
        afirstpoint: &GpPnt2d,
        asecondcurve: &dyn Curve2d,
        apoint: &GpPnt2d,
        afirstvector: &GpVec2d,
        asecondvector: &GpVec2d,
        adirection: f64,
        tolerance: f64,
        oncurve: bool,
    ) -> Result<(), BisectorBisecError> {
        let mut u_first;
        let mut u_last;
        let bis: Arc<dyn Curve2d>;

        if is_circle_or_line(asecondcurve) {
            // `cxx:426-444`: analytic bisector.
            let mut bis_ana = BisectorBisecAna::new();
            bis_ana.perform_point_curve(
                afirstpoint,
                asecondcurve,
                apoint,
                afirstvector,
                asecondvector,
                adirection,
                tolerance,
                oncurve,
            );
            u_first = bis_ana.parameter_of_start_point();
            u_last = bis_ana.parameter_of_end_point();
            bis = Arc::from(bis_ana.clone_dyn());
        } else {
            let mut is_line = false;
            let mut rc = INFINITE;
            if oncurve {
                if is_convex(asecondcurve, adirection)
                    || is_max_rc(asecondcurve, asecondcurve.first_parameter(), &mut rc)
                {
                    is_line = true;
                }
            }
            if is_line {
                // `cxx:459-471`: half-straight.
                let bis_ana = half_straight(apoint, afirstvector, adirection, rc);
                u_first = bis_ana.parameter_of_start_point();
                u_last = bis_ana.parameter_of_end_point();
                bis = Arc::from(bis_ana.clone_dyn());
            } else {
                // `cxx:477-522`: algorithmic bisector.
                let mut bis_pc = BisectorBisecPC::new();
                bis_pc.perform(asecondcurve, afirstpoint, adirection, INFINITE);
                if bis_pc.is_empty() {
                    // `cxx:480-498`: false bisector.
                    let bis_ana = false_bisector(
                        apoint,
                        afirstvector,
                        asecondvector,
                        adirection,
                        rc,
                    );
                    u_first = bis_ana.parameter_of_start_point();
                    u_last = bis_ana.parameter_of_end_point();
                    bis = Arc::from(bis_ana.clone_dyn());
                } else {
                    u_first = bis_pc.parameter(apoint);
                    u_last = bis_pc.last_parameter();
                    if u_first >= u_last {
                        // `cxx:507-517`: extrapolate by line.
                        let v = vec_between(
                            &bis_pc.value(bis_pc.first_parameter()),
                            &bis_pc.value(u_last),
                        );
                        let line = line_through(apoint, &v, rc);
                        u_first = line.parameter_of_start_point();
                        u_last = line.parameter_of_end_point();
                        bis = Arc::from(line.clone_dyn());
                    } else {
                        bis = Arc::from(bis_pc.clone_dyn());
                    }
                }
            }
        }

        if u_first < bis.first_parameter() {
            u_first = bis.first_parameter();
        }
        if u_last > bis.last_parameter() {
            u_last = bis.last_parameter();
        }
        self.thebisector = Some(Geom2dTrimmedCurve::new(bis, u_first, u_last));
        Ok(())
    }

    /// `Perform(afirstpoint, asecondpoint, ...)` (`cxx:543-564`): the bisector
    /// between two points is always analytic.
    #[allow(clippy::too_many_arguments)]
    pub fn perform_point_point(
        &mut self,
        afirstpoint: &GpPnt2d,
        asecondpoint: &GpPnt2d,
        apoint: &GpPnt2d,
        afirstvector: &GpVec2d,
        asecondvector: &GpVec2d,
        adirection: f64,
        tolerance: f64,
        oncurve: bool,
    ) -> Result<(), BisectorBisecError> {
        let mut bis = BisectorBisecAna::new();
        bis.perform_point_point(
            afirstpoint,
            asecondpoint,
            apoint,
            afirstvector,
            asecondvector,
            adirection,
            tolerance,
            oncurve,
        );
        let u_first = bis.parameter_of_start_point();
        let u_last = bis.parameter_of_end_point();
        let curve: Arc<dyn Curve2d> = Arc::from(bis.clone_dyn());
        self.thebisector = Some(Geom2dTrimmedCurve::new(curve, u_first, u_last));
        Ok(())
    }
}

/// `Type == Geom2d_Circle || Type == Geom2d_Line` (`cxx:78-86`, `:112-123`): the
/// dynamic type is read on the basis of a `Geom2d_TrimmedCurve`.
fn is_circle_or_line(c: &dyn Curve2d) -> bool {
    let t = c.trimmed_basis().unwrap_or(c);
    t.is_geom2d_circle() || t.is_geom2d_line()
}

/// `cxx:90-132`: a degree-1 two-pole B-spline whose poles are closer than
/// `1.e-4` is replaced by the segment between those poles.
fn collapsed_segment(c: &dyn Curve2d) -> Option<Box<dyn Curve2d>> {
    let basis = c.trimmed_basis().unwrap_or(c);
    if basis.bspline_degree() != Some(1) {
        return None;
    }
    let poles = basis.poles2d()?;
    if poles.len() != 2 || poles[0].distance(&poles[1]) >= 1.0e-4 {
        return None;
    }
    Some(make_segment2d(&poles[0], &poles[1]))
}

/// `GC_MakeSegment2d(P1, P2).Value()` (`GCE2d_MakeSegment.cxx`): the line through
/// `P1` towards `P2`, trimmed from `0` to `P1.Distance(P2)`.
fn make_segment2d(p1: &GpPnt2d, p2: &GpPnt2d) -> Box<dyn Curve2d> {
    Box::new(crate::gce2d::make_segment(p1, p2))
}

/// `gp_Vec2d(P1, P2)` (`gp_Vec2d.hxx`): the translation `P2 - P1`.
fn vec_between(p1: &GpPnt2d, p2: &GpPnt2d) -> GpVec2d {
    GpVec2d::new(p2.x() - p1.x(), p2.y() - p1.y())
}

/// `new Geom2d_Line(apoint, N)` + `Bisector_BisecAna::Init` over `[0, u_last]`
/// (`cxx:172-180`, `:313-322`, `:458-469`, `:366-373`, `:508-515`).
fn half_straight(
    apoint: &GpPnt2d,
    afirstvector: &GpVec2d,
    adirection: f64,
    u_last: f64,
) -> BisectorBisecAna {
    let n = GpDir2d::new(
        -adirection * afirstvector.y(),
        adirection * afirstvector.x(),
    )
    .unwrap_or_default();
    let line = Geom2dLine::new(GpAx2d::new(*apoint, n));
    let basis: Arc<dyn Curve2d> = Arc::from(line.clone_dyn());
    let mut bis = BisectorBisecAna::new();
    bis.init(Geom2dTrimmedCurve::new(basis, 0.0, u_last));
    bis
}

/// `cxx:196-232` / `:336-352` / `:482-497`: the false bisector built when the
/// algorithmic construction comes back empty. Its direction is `-(dir1+dir2)`
/// unless that cancels, in which case it is the normal of the first vector.
fn false_bisector(
    apoint: &GpPnt2d,
    afirstvector: &GpVec2d,
    asecondvector: &GpVec2d,
    adirection: f64,
    u_last: f64,
) -> BisectorBisecAna {
    let dir1 = GpDir2d::from_vec2d(afirstvector).unwrap_or_default();
    let dir2 = GpDir2d::from_vec2d(asecondvector).unwrap_or_default();
    let mut nx = -dir1.x() - dir2.x();
    let mut ny = -dir1.y() - dir2.y();
    if nx.abs() <= RESOLUTION && ny.abs() <= RESOLUTION {
        nx = -afirstvector.y();
        ny = afirstvector.x();
    }
    let n = GpDir2d::new(adirection * nx, adirection * ny).unwrap_or_default();
    let line = Geom2dLine::new(GpAx2d::new(*apoint, n));
    let basis: Arc<dyn Curve2d> = Arc::from(line.clone_dyn());
    let mut bis = BisectorBisecAna::new();
    bis.init(Geom2dTrimmedCurve::new(basis, 0.0, u_last));
    bis
}

/// `new Geom2d_Line(apoint, gp_Dir2d(V))` + `Bisector_BisecAna::Init` over
/// `[0, RC]` (`cxx:360-372`, `:505-517`).
fn line_through(apoint: &GpPnt2d, v: &GpVec2d, u_last: f64) -> BisectorBisecAna {
    let n = GpDir2d::from_vec2d(v).unwrap_or_default();
    let line = Geom2dLine::new(GpAx2d::new(*apoint, n));
    let basis: Arc<dyn Curve2d> = Arc::from(line.clone_dyn());
    let mut bis = BisectorBisecAna::new();
    bis.init(Geom2dTrimmedCurve::new(basis, 0.0, u_last));
    bis
}

/// `ReplaceByLineIfIsToSmall(Bis, UFirst, ULast)` (`cxx:586-612`): an
/// algorithmic bisector whose size is negligible is replaced by a half-straight.
fn replace_by_line_if_is_to_small(
    bis: &Arc<dyn Curve2d>,
    u_first: &mut f64,
    u_last: &mut f64,
) -> Arc<dyn Curve2d> {
    if (*u_last - *u_first).abs() > 2.0 * PCONFUSION * 10.0 {
        return bis.clone();
    }
    let pf = bis.value(*u_first);
    let pl = bis.value(*u_last);
    if pf.distance(&pl) > CONFUSION * 10.0 {
        return bis.clone();
    }
    let t1 = bis.eval_dn(*u_first, 1);
    let replaced = line_through(&pf, &t1, INFINITE);
    *u_first = replaced.parameter_of_start_point();
    *u_last = replaced.parameter_of_end_point();
    Arc::from(replaced.clone_dyn())
}

/// Curvature at one end, the `KF` / `KL` computation of `IsMaxRC`
/// (`cxx:625-646`); `Norm2 < gp::Resolution()` answers `0`, and
/// `gp::Resolution()` is `RealEpsilon()`.
fn end_curvature(c: &dyn Curve2d, u: f64) -> f64 {
    let (_, d1, d2) = c.d2(u);
    let norm2 = d1.square_magnitude();
    if norm2 < REAL_EPSILON {
        0.0
    } else {
        d1.crossed(&d2).abs() / (norm2 * norm2.sqrt())
    }
}

/// `IsMaxRC(C, U, R)` (`cxx:615-681`): true when the curvature at `U` is the
/// maximum over the two ends of `C`, in which case `R` is the corresponding
/// radius of curvature.
fn is_max_rc(c: &dyn Curve2d, u: f64, r: &mut f64) -> bool {
    let us = c.first_parameter();
    let ul = c.last_parameter();
    let kf = end_curvature(c, us);
    let kl = end_curvature(c, ul);
    if u == ul {
        if kl < kf {
            *r = if kl == 0.0 { INFINITE } else { 1.0 / kl };
            return true;
        }
    } else if kf < kl {
        *r = if kf == 0.0 { INFINITE } else { 1.0 / kf };
        return true;
    }
    false
}
