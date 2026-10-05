//! `GeomLib` — the arms needed by `BRepLib::BuildCurve3d`.
//!
//! Ported:
//! * `GeomLib::To3d(const gp_Ax2&, const Handle(Geom2d_Curve)&)`
//!   (`GeomLib.cxx:559-679`);
//! * `GeomLib::isIsoLine` (`GeomLib.cxx:2991-3077`);
//! * `GeomLib::buildC3dOnIsoLine` (`GeomLib.cxx:3079-3221`);
//! * `GeomLib::BuildCurve3d` (`GeomLib.cxx:1051-1165`).
//!
//! `BuildCurve3d`'s OCCT parameter is an `Adaptor3d_CurveOnSurface&`. The port
//! passes the same object as `&dyn Curve` (the port's `CurveOnSurface`
//! implements `Curve`) together with its `GetSurface()` / `GetCurve()`
//! handles, which is what the three arms of the OCCT body actually read.

use std::sync::Arc;

use occt_core::elib::clib;
use occt_core::gp::{GpAx2, GpDir2d, GpPnt2d, GpVec2d};
use occt_core::kernel::geomabs::Shape;
use occt_core::precision::{Precision, ANGULAR, CONFUSION, PCONFUSION};
use occt_geom2d::curve::Curve2d;

use crate::adv_approx::{ApproxAFunction3d, PrefAndRec};
use crate::bezier_curve::GeomBezierCurve;
use crate::bspline_curve::GeomBSplineCurve;
use crate::circle::GeomCircle;
use crate::curve::Curve;
use crate::ellipse::GeomEllipse;
use crate::hyperbola::GeomHyperbola;
use crate::line::GeomLine;
use crate::offset::GeomOffsetCurve;
use crate::parabola::GeomParabola;
use crate::surface::Surface;
use crate::trimmed::GeomTrimmedCurveBasis;

/// `gp::DX2d()`, the reference direction `ElCLib`/`GeomLib` compare against.
fn dx2d() -> GpDir2d {
    GpDir2d::new(1.0, 0.0).expect("DX2d")
}

/// `gp::DY2d()`.
fn dy2d() -> GpDir2d {
    GpDir2d::new(0.0, 1.0).expect("DY2d")
}

/// `GeomLib::To3d(Position, Curve2d)` (`GeomLib.cxx:559-679`).
///
/// The dispatch is on `Curve2d->DynamicType()`; the port's `Curve2d` hooks
/// (`trimmed_basis` / `offset_basis` / `is_bezier2d` / `is_bspline2d` /
/// `gp_lin2d` / `gp_circ2d` / `gp_elips2d` / `gp_parab2d` / `gp_hypr2d`) are
/// each `Some` for exactly one `Geom2d_*` class.
///
/// `None` stands for the `Standard_NotImplemented` thrown at `cxx:674-676`
/// for an unrecognised curve type.
///
/// UNPORTED arms, both unreachable through the port's 2D types:
/// the rational `Geom_BezierCurve(Poles, Weights)` (`cxx:592-596`) because
/// `crate::bezier_curve::GeomBezierCurve` carries no weights, and the periodic
/// `Geom_BSplineCurve(Poles, Knots, Mults, Degree, IsPeriodic)` (`cxx:619-629`)
/// because `Geom_BSplineCurve::from_poles_knots_mults` builds the
/// non-periodic flat sequence only. Both report `None` rather than a curve
/// that is not the OCCT one.
pub fn to_3d(position: &GpAx2, curve2d: &dyn Curve2d) -> Option<Arc<dyn Curve>> {
    // `Geom2d_TrimmedCurve` (`cxx:565-573`).
    if let Some(basis) = curve2d.trimmed_basis() {
        let u1 = curve2d.first_parameter();
        let u2 = curve2d.last_parameter();
        let cc = to_3d(position, basis)?;
        // `new Geom_TrimmedCurve(CC, U1, U2)`: both parameters live on the
        // basis curve, which is exactly `GeomTrimmedCurveBasis`.
        return Some(Arc::new(GeomTrimmedCurveBasis::new(cc, u1, u2)));
    }

    // `Geom2d_OffsetCurve` (`cxx:575-581`).
    if let Some(basis) = curve2d.offset_basis() {
        let offset = curve2d.offset_value()?;
        let cc = to_3d(position, basis)?;
        return Some(Arc::new(GeomOffsetCurve::new(
            cc,
            offset,
            position.direction(),
        )));
    }

    // `Geom2d_BezierCurve` (`cxx:582-602`).
    if curve2d.is_bezier2d() {
        let poles2d = curve2d.poles2d()?;
        let poles3d = poles2d.iter().map(|p| clib::to_3d_pnt(position, p)).collect();
        if curve2d.is_rational() {
            // `new Geom_BezierCurve(Poles3d, CBSpl2d->WeightsArray())`.
            // UNPORTED: no rational 3D Bezier in the port.
            return None;
        }
        // `new Geom_BezierCurve(Poles3d)` (`cxx:598`).
        return GeomBezierCurve::new(poles3d)
            .ok()
            .map(|c| Arc::new(c) as Arc<dyn Curve>);
    }

    // `Geom2d_BSplineCurve` (`cxx:603-632`).
    if curve2d.is_bspline2d() {
        let poles2d = curve2d.poles2d()?;
        let poles3d = poles2d.iter().map(|p| clib::to_3d_pnt(position, p)).collect();
        let degree = curve2d.bspline_degree()?;
        let (knots, mults) = curve2d.bspline_distinct_knots_mults()?;
        // `IsPeriodic = CBSpl2d->IsPeriodic()` (`cxx:607`).
        // UNPORTED: `Geom_BSplineCurve(Poles, Knots, Mults, Degree, IsPeriodic)`
        // for `IsPeriodic == true` (`cxx:619-629`).
        if curve2d.is_periodic() {
            return None;
        }
        if let Some(weights) = curve2d.bspline_weights2d() {
            return GeomBSplineCurve::from_poles_knots_mults(poles3d, knots, mults, degree)
                .ok()
                .and_then(|c| {
                    GeomBSplineCurve::rational(
                        c.poles.clone(),
                        weights.to_vec(),
                        c.knots.clone(),
                        degree,
                    )
                    .ok()
                })
                .map(|c| Arc::new(c) as Arc<dyn Curve>);
        }
        return GeomBSplineCurve::from_poles_knots_mults(poles3d, knots, mults, degree)
            .ok()
            .map(|c| Arc::new(c) as Arc<dyn Curve>);
    }

    // `Geom2d_Line` (`cxx:633-640`).
    if let Some(l) = curve2d.gp_lin2d() {
        return Some(Arc::new(GeomLine::new(clib::to_3d_lin(position, &l))));
    }

    // `Geom2d_Circle` (`cxx:641-648`).
    if let Some(c) = curve2d.gp_circ2d() {
        return clib::to_3d_circ(position, &c)
            .ok()
            .map(|c3| Arc::new(GeomCircle::new(c3)) as Arc<dyn Curve>);
    }

    // `Geom2d_Ellipse` (`cxx:649-656`).
    if let Some(e) = curve2d.gp_elips2d() {
        return clib::to_3d_elips(position, &e)
            .ok()
            .map(|e3| Arc::new(GeomEllipse::new(e3)) as Arc<dyn Curve>);
    }

    // `Geom2d_Parabola` (`cxx:657-664`).
    if let Some(p) = curve2d.gp_parab2d() {
        return clib::to_3d_parab(position, &p)
            .ok()
            .map(|p3| Arc::new(GeomParabola::new(p3)) as Arc<dyn Curve>);
    }

    // `Geom2d_Hyperbola` (`cxx:665-672`).
    if let Some(h) = curve2d.gp_hypr2d() {
        return clib::to_3d_hypr(position, &h)
            .ok()
            .map(|h3| Arc::new(GeomHyperbola::new(h3)) as Arc<dyn Curve>);
    }

    // `throw Standard_NotImplemented()` (`cxx:674-676`).
    None
}

/// Result of `GeomLib::isIsoLine`'s three out-parameters.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IsoLine {
    /// `theIsU`: `true` when the iso is `U = theParam` (a vertical line).
    pub is_u: bool,
    /// `theParam`: the constant parameter.
    pub param: f64,
    /// `theIsForward`: `theDir . gp::DX2d()/DY2d() > 0`.
    pub is_forward: bool,
}

/// `GeomLib::isIsoLine(theC2D, theIsU, theParam, theIsForward)`
/// (`GeomLib.cxx:2991-3077`).
///
/// `Geom2dAdaptor_Curve::load` / `Load` unwraps a `Geom2d_TrimmedCurve` onto
/// its basis (`Geom2dAdaptor_Curve.cxx:285-288`), so the type test runs on the
/// basis curve; that unwrapping is done here as well.
pub fn is_iso_line(c2d: &dyn Curve2d) -> Option<IsoLine> {
    // `load` recursion onto the basis (`cxx:285-288`).
    let mut curve = c2d;
    while let Some(basis) = curve.trimmed_basis() {
        curve = basis;
    }

    // These variables are used to check line state (vertical or horizontal).
    let loc2d: GpPnt2d;
    let dir2d: Option<GpDir2d>;

    if let Some(l) = curve.gp_lin2d() {
        // `aType == GeomAbs_Line` (`cxx:3003-3009`).
        loc2d = l.location();
        dir2d = Some(*l.direction());
    } else if curve.is_bspline2d() {
        // `aType == GeomAbs_BSplineCurve` (`cxx:3010-3029`).
        let poles = curve.poles2d()?;
        if curve.bspline_degree() != Some(1) || poles.len() != 2 {
            return None; // Not a line or uneven parameterization.
        }
        loc2d = poles[0];
        // Vector should be non-degenerated.
        let v = GpVec2d::new(poles[1].x() - poles[0].x(), poles[1].y() - poles[0].y());
        if v.square_magnitude() < CONFUSION {
            return None; // Degenerated spline.
        }
        dir2d = GpDir2d::from_vec2d(&v).ok();
    } else if curve.is_bezier2d() {
        // `aType == GeomAbs_BezierCurve` (`cxx:3030-3049`).
        let poles = curve.poles2d()?;
        if curve.bezier_nb_poles() != Some(2) {
            return None;
        }
        loc2d = poles[0];
        let v = GpVec2d::new(poles[1].x() - poles[0].x(), poles[1].y() - poles[0].y());
        if v.square_magnitude() < CONFUSION {
            return None;
        }
        dir2d = GpDir2d::from_vec2d(&v).ok();
    } else {
        // `!isAppropriateType` (`cxx:3050-3053`).
        return None;
    }

    let dir2d = dir2d?;

    // Check line to be vertical or horizontal (`cxx:3055-3075`).
    if dir2d.is_parallel(&dx2d(), ANGULAR) {
        // Horizontal line. V = const.
        Some(IsoLine {
            is_u: false,
            param: loc2d.y(),
            is_forward: dir2d.dot(&dx2d()) > 0.0,
        })
    } else if dir2d.is_parallel(&dy2d(), ANGULAR) {
        // Vertical line. U = const.
        Some(IsoLine {
            is_u: true,
            param: loc2d.x(),
            is_forward: dir2d.dot(&dy2d()) > 0.0,
        })
    } else {
        None
    }
}

/// `GeomLib::buildC3dOnIsoLine` (`GeomLib.cxx:3079-3221`).
///
/// `theFirst` / `theLast` are the `Adaptor3d_CurveOnSurface` range and
/// `theParam` the constant iso parameter from [`is_iso_line`]. `None` is the
/// OCCT `return occ::handle<Geom_Curve>()` (sphere, out-of-bounds iso,
/// degenerate span, or the deviation check at `cxx:3216-3219`).
pub fn build_c3d_on_iso_line(
    c2d: &dyn Curve2d,
    surface: &Arc<dyn Surface>,
    the_first: f64,
    the_last: f64,
    tolerance: f64,
    is_u: bool,
    param: f64,
    is_forward: bool,
) -> Option<Arc<dyn Curve>> {
    // `theSurf->GetType() == GeomAbs_Sphere` (`cxx:3095-3098`).
    if surface.gp_sphere().is_some() {
        return None;
    }

    // Extract isoline (`cxx:3100-3105`).
    let mut a_surf: Arc<dyn Surface> = surface.clone();

    // `aF2d = theC2D->Value(FirstParameter())`, `aL2d = ... Value(LastParameter())`.
    let f2d = c2d.d0(the_first);
    let l2d = c2d.d0(the_last);

    let (u1, u2) = a_surf.u_range();
    let (v1, v2) = a_surf.v_range();
    let mut is_to_trim = true;
    let a_iso: Arc<dyn Curve>;
    // `aV1Param` / `aV2Param` (`theIsU`) or `aU1Param` / `aU2Param` — the trim
    // range used by `new Geom_TrimmedCurve` below.
    let trim_first: f64;
    let trim_last: f64;

    if is_u {
        // `cxx:3111-3145`.
        let mut a_v1 = f2d.y().min(l2d.y());
        let mut a_v2 = f2d.y().max(l2d.y());
        if a_v2 < v1 - tolerance || a_v1 > v2 + tolerance {
            return None;
        } else if Precision::is_infinite(v1) || Precision::is_infinite(v2) {
            if (a_v2 - a_v1).abs() < PCONFUSION {
                return None;
            }
            a_surf = Arc::new(
                crate::rectangular_trimmed::GeomRectangularTrimmedSurface::uv(
                    a_surf, u1, u2, a_v1, a_v2,
                ),
            );
            is_to_trim = false;
        } else {
            a_v1 = a_v1.max(v1);
            a_v2 = a_v2.min(v2);
            if (a_v2 - a_v1).abs() < PCONFUSION {
                return None;
            }
        }
        trim_first = a_v1;
        trim_last = a_v2;
        a_iso = a_surf.u_iso_curve(param)?;
    } else {
        // `cxx:3146-3173`.
        let mut a_u1 = f2d.x().min(l2d.x());
        let mut a_u2 = f2d.x().max(l2d.x());
        if a_u2 < u1 - tolerance || a_u1 > u2 + tolerance {
            return None;
        } else if Precision::is_infinite(u1) || Precision::is_infinite(u2) {
            if (a_u2 - a_u1).abs() < PCONFUSION {
                return None;
            }
            a_surf = Arc::new(
                crate::rectangular_trimmed::GeomRectangularTrimmedSurface::uv(
                    a_surf, a_u1, a_u2, v1, v2,
                ),
            );
            is_to_trim = false;
        } else {
            a_u1 = a_u1.max(u1);
            a_u2 = a_u2.min(u2);
            if (a_u2 - a_u1).abs() < PCONFUSION {
                return None;
            }
        }
        trim_first = a_u1;
        trim_last = a_u2;
        a_iso = a_surf.v_iso_curve(param)?;
    }

    // `if (isToTrim) aC3d = new Geom_TrimmedCurve(aC3d, aV1Param, aV2Param)`
    // (`cxx:3137-3142` / `cxx:3169-3172`): the trim keeps the basis
    // parameters, which is `GeomTrimmedCurveBasis`.
    let a_c3d: Arc<dyn Curve> = if is_to_trim {
        Arc::new(GeomTrimmedCurveBasis::new(a_iso, trim_first, trim_last))
    } else {
        a_iso
    };

    // Convert arbitrary curve type to the b-spline (`cxx:3175-3178`).
    let mut curve3d =
        crate::convert_bspl::curve_to_bspline_curve(
            a_c3d.as_ref(),
            occt_core::convert::ParameterisationType::QuasiAngular,
        )
        .ok()?;
    if !is_forward {
        // `aCurve3d->Reverse()` (`cxx:3179-3182`).
        curve3d.reverse();
    }

    // Rebuild parameterization for the 3d curve to have the same
    // parameterization with a two-dimensional curve (`cxx:3184-3188`).
    let (mut knots, _mults) = curve3d.distinct_knots_and_mults();
    knots_reparametrize(c2d.first_parameter(), c2d.last_parameter(), &mut knots);
    curve3d.set_knots(&knots).ok()?;

    // Evaluate error (`cxx:3190-3213`).
    let mut error3d = 0.0f64;
    let par_f = the_first;
    let par_l = the_last;
    const NB_PNT: i32 = 23;
    for idx in 0..=NB_PNT {
        let par = par_f + ((par_l - par_f) * idx as f64) / NB_PNT as f64;
        let pnt2d = c2d.d0(par);
        let pnt_c3d = curve3d.d0(par);
        let pnt_c2d = surface.d0(pnt2d.x(), pnt2d.y());
        let sq_deviation = pnt_c3d.square_distance(&pnt_c2d);
        error3d = error3d.max(sq_deviation);
    }
    let error3d = error3d.sqrt();

    // Target tolerance is not obtained (`cxx:3216-3219`).
    if error3d > tolerance {
        return None;
    }

    Some(Arc::new(curve3d))
}

/// `BSplCLib::Reparametrize(U1, U2, Knots)`; a thin re-export so the arm above
/// reads like the OCCT one.
pub fn knots_reparametrize(u1: f64, u2: f64, knots: &mut [f64]) {
    occt_core::bspl::knots::reparametrize(u1, u2, knots);
}

/// Result of [`build_curve3d`]: the OCCT out-parameters
/// `NewCurvePtr`, `MaxDeviation` and `AverageDeviation`.
pub struct BuildCurve3d {
    /// `NewCurvePtr`; `None` where OCCT leaves the handle null (and, for the
    /// plane arm, where `GeomLib::To3d` throws `Standard_NotImplemented`).
    pub curve: Option<Arc<dyn Curve>>,
    /// `MaxDeviation` — non-zero only on the approximation arm.
    pub max_deviation: f64,
    /// `AverageDeviation` — non-zero only on the approximation arm.
    pub average_deviation: f64,
}

/// `GeomLib::BuildCurve3d(Tolerance, Curve, FirstParameter, LastParameter,
/// NewCurvePtr, MaxDeviation, AverageDeviation, Continuity, MaxDegree,
/// MaxSegment)` (`GeomLib.cxx:1051-1165`).
///
/// `curve` is the `Adaptor3d_CurveOnSurface` (`GetSurface()` / `GetCurve()` are
/// passed separately because the `Curve` trait does not expose them). The plane
/// arm (`cxx:1077-1093`) and the iso-line arm (`cxx:1096-1114`) come from
/// [`to_3d`] and [`build_c3d_on_iso_line`]; the general arm (`cxx:1116-1164`)
/// runs `AdvApprox_ApproxAFunction` over the composed `D0`/`D1`/`D2`.
///
/// UNPORTED (general arm): `GeomLib_CurveOnSurfaceEvaluator` re-trims the
/// `CurveOnSurface` to each interval handed to `Evaluate`
/// (`GeomLib.cxx:1009-1016`). The port evaluates `curve` directly, which gives
/// the same point/tangent/curvature for interior parameters; only the
/// `EvalFirstLastSurf` end patch of a sub-interval differs.
///
/// UNPORTED (cut points): OCCT reads them from
/// `Adaptor3d_CurveOnSurface::NbIntervals` (`Adaptor3d_CurveOnSurface.cxx:1045-1114`),
/// which merges the 2D curve's spans with the surface U/V discontinuity
/// parameters crossed by the curve. The port's `CurveOnSurface` keeps the
/// `Curve` trait default (one span `[first, last]`); `AdvApprox` then cuts
/// purely adaptively.
#[allow(clippy::too_many_arguments)]
pub fn build_curve3d(
    tolerance: f64,
    curve: &dyn Curve,
    surface: &Arc<dyn Surface>,
    pcurve: &Arc<dyn Curve2d>,
    first: f64,
    last: f64,
    continuity: Shape,
    max_degree: i32,
    max_segment: i32,
) -> BuildCurve3d {
    let mut max_deviation = 0.0e0;
    let mut average_deviation = 0.0e0;

    // `RT = down_cast<Geom_RectangularTrimmedSurface>(geom_surface.Surface())`;
    // `P = RT.IsNull() ? Geom_Plane : RT->BasisSurface()` (`cxx:1077-1086`).
    let plane = match surface.rectangular_trimmed_basis() {
        Some(basis) => basis,
        None => surface.clone(),
    };

    if let Some(pln) = plane.gp_pln() {
        // Compute the 3d curve (`cxx:1088-1092`). `GeomLib::To3d` throws for an
        // unknown 2D type; the handle stays null here and `BRepLib::BuildCurve3d`
        // turns that into `false` (`BRepLib.cxx:361-364`).
        let axes = pln.position().ax2();
        return BuildCurve3d {
            curve: to_3d(&axes, pcurve.as_ref()),
            max_deviation,
            average_deviation,
        };
    }

    // `TrimmedC2D = geom_adaptor_curve_ptr->Trim(First, Last, PConfusion())`;
    // `Geom2dAdaptor_Curve::Trim` only rewrites the kept range
    // (`Geom2dAdaptor_Curve.cxx:577-584`), so the type/geometry queries below
    // run on the same 2D curve (`cxx:1096-1114`).
    if let Some(iso) = is_iso_line(pcurve.as_ref()) {
        let c3d = build_c3d_on_iso_line(
            pcurve.as_ref(),
            surface,
            first,
            last,
            tolerance,
            iso.is_u,
            iso.param,
            iso.is_forward,
        );
        if c3d.is_some() {
            return BuildCurve3d {
                curve: c3d,
                max_deviation,
                average_deviation,
            };
        }
    }

    //
    // Entree
    //
    // Search for discontinuities (`cxx:1123-1133`).
    let cut_pnts_c2 = curve.parameter_intervals(4);
    let cut_pnts_c3 = curve.parameter_intervals(5);
    let cut_tool = match PrefAndRec::new(cut_pnts_c2, cut_pnts_c3, 5.0) {
        Ok(cut_tool) => cut_tool,
        Err(()) => {
            return BuildCurve3d {
                curve: None,
                max_deviation,
                average_deviation,
            }
        }
    };

    // `GeomLib_CurveOnSurfaceEvaluator::Evaluate` (`cxx:1001-1047`): order 0 is
    // the point, 1 the first derivative, 2 the second; any other order returns
    // error code 3 with a zero result.
    let eval = |u: f64, order: i32, out: &mut [f64]| -> i32 {
        match order {
            0 => {
                let p = curve.d0(u);
                out[0] = p.x();
                out[1] = p.y();
                out[2] = p.z();
            }
            1 => {
                let (_, v) = curve.d1(u);
                out[0] = v.x();
                out[1] = v.y();
                out[2] = v.z();
            }
            2 => {
                let (_, _, v) = curve.d2(u);
                out[0] = v.x();
                out[1] = v.y();
                out[2] = v.z();
            }
            _ => {
                out[0] = 0.0;
                out[1] = 0.0;
                out[2] = 0.0;
                return 3;
            }
        }
        0
    };

    // `AdvApprox_ApproxAFunction` (`cxx:1139-1152`): Num1DSS = Num2DSS = 0,
    // Num3DSS = 1, Tolerance3D = {Tolerance}.
    match ApproxAFunction3d::approx(
        first,
        last,
        continuity,
        max_degree,
        max_segment,
        tolerance,
        &cut_tool,
        &eval,
    ) {
        Ok(approx) => {
            if approx.has_result {
                // `GeomLib_MakeCurvefromApprox(anApproximator).Curve(1)`
                // (`cxx:1156-1162`).
                let curve3d = GeomBSplineCurve::from_poles_knots_mults(
                    approx.poles.clone(),
                    approx.knots.clone(),
                    approx.mults.clone(),
                    approx.degree.max(0) as usize,
                )
                .ok();
                max_deviation = approx.max_error;
                average_deviation = approx.average_error;
                BuildCurve3d {
                    curve: curve3d.map(|c| Arc::new(c) as Arc<dyn Curve>),
                    max_deviation,
                    average_deviation,
                }
            } else {
                BuildCurve3d {
                    curve: None,
                    max_deviation,
                    average_deviation,
                }
            }
        }
        Err(()) => BuildCurve3d {
            curve: None,
            max_deviation,
            average_deviation,
        },
    }
}
