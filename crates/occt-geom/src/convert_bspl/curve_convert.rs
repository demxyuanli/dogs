// ---------------------------------------------------------------------------
// CurveToBSplineCurve.
// ---------------------------------------------------------------------------

use occt_core::bspl::banded_interp::knot_sequence;
use occt_core::bspl::knots;
use occt_core::convert::{
    circle_to_bspline_curve, circle_to_bspline_curve_range, ellipse_to_bspline_curve,
    ellipse_to_bspline_curve_range, hyperbola_to_bspline_curve, parabola_to_bspline_curve,
    ConicToBSplineCurve, ConvertError, ParameterisationType,
};
use occt_core::gp::{
    GpAx2, GpAx22d, GpAx3, GpCirc2d, GpElips2d, GpHypr2d, GpParab2d, GpPnt, GpTrsf,
};
use occt_core::kernel::geomabs::Shape;

use crate::bezier_curve::GeomBezierCurve;
use crate::bspline_curve::GeomBSplineCurve;
use crate::convert_approx_curve::GeomConvertApproxCurve;
use crate::curve::Curve;

use super::comp_curve::CompCurveToBSplineCurve;

/// `Geom_BSplineCurve`'s file-static `Rational(Weights)`
/// (`Geom_BSplineCurve.cxx:97-108`): the curve is rational as soon as two
/// consecutive weights differ by more than `gp::Resolution()` (`REAL_SMALL`).
pub(super) fn weights_are_rational(weights: &[f64]) -> bool {
    for w in weights.windows(2) {
        if (w[0] - w[1]).abs() > occt_core::precision::REAL_SMALL {
            return true;
        }
    }
    false
}

/// `BSplineCurveBuilder` (`GeomConvert.cxx:60-83`): lift the converter's 2D
/// poles to `z = 0`, build the B-spline from its knot/multiplicity tables, then
/// move it from the conic's own frame into the world by
/// `gp_Trsf::SetTransformation(TheConic->Position(), gp::XOY())`.
///
/// OCCT builds the **rational** constructor and lets `CheckRational` decide
/// (`Geom_BSplineCurve.cxx:205-223`); [`weights_are_rational`] reproduces that
/// decision so `IsRational()` matches (a parabola's unit weights end up
/// non-rational, a circular arc's do not).
fn bspline_curve_builder(
    conic_position: &GpAx2,
    convert: &ConicToBSplineCurve,
) -> Result<GeomBSplineCurve, ConvertError> {
    // `CheckCurveData` (`Geom_BSplineCurve.cxx:91-94`):
    // `Poles.Length() == BSplCLib::NbPoles(Degree, Periodic, Mults)`.
    if knots::nb_poles(
        convert.degree() as i32,
        convert.is_periodic(),
        convert.multiplicities(),
    ) as usize
        != convert.poles().len()
    {
        return Err(ConvertError::ConstructionError);
    }
    let poles: Vec<GpPnt> = convert
        .poles()
        .iter()
        .map(|p| GpPnt::new(p.x(), p.y(), 0.0))
        .collect();
    // `Geom_BSplineCurve`'s rational/non-rational constructors call
    // `updateKnots()`, i.e. `BSplCLib::KnotSequence(…, Periodic)`; for a
    // periodic result that is the sequence extended by one period each side.
    let flat = if convert.is_periodic() {
        knots::knot_sequence_periodic(
            convert.knots(),
            convert.multiplicities(),
            convert.degree() as i32,
        )
    } else {
        knot_sequence(
            convert.knots(),
            convert.multiplicities(),
            convert.degree() as i32,
        )
    };
    let weights = convert.weights().to_vec();
    let weights = if weights_are_rational(&weights) { Some(weights) } else { None };
    let mut curve = GeomBSplineCurve::from_parts(
        poles,
        weights,
        flat,
        convert.degree(),
        convert.is_periodic(),
    );
    let mut trsf = GpTrsf::identity();
    trsf.set_transformation_from_to(&conic_position.to_ax3(), &GpAx3::standard());
    Curve::transform(&mut curve, &trsf);
    Ok(curve)
}

/// Dynamic type test of OCCT's `IsKind(STANDARD_TYPE(Geom_BSplineCurve))`
/// down_cast. The port's `Curve` trait is not `Any`, so the B-spline query is
/// the discriminator; a `Geom_TrimmedCurve` delegates that query to its basis
/// (`trimmed.rs:127`) and is excluded here exactly as the OCCT down_cast
/// excludes it.
///
/// Shared with `ShapeUpgrade::C0BSplineToSequenceOfC1BSplineCurve`, which the
/// IGES read path reaches at `IGESToBRep_TopoCurve.cxx:1016`.
pub fn is_bspline(c: &dyn Curve) -> bool {
    !c.is_geom_trimmed() && c.bspline_poles().is_some()
}

/// `Geom_BSplineCurve::Copy()` (`Geom_BSplineCurve.cxx:112-115`): a verbatim
/// copy of the `Geom_BSplineCurve` fields, recovered from the `Curve`
/// queries (`bspline_poles`/`bspline_weights`/`bspline_knots`/
/// `nurbs_degree`/`is_periodic`). Equal data, not a re-approximation.
///
/// Callers must have taken the `is_bspline` test first; the result carries no
/// guard of its own because OCCT's counterpart is an unchecked down_cast.
pub fn bspline_copy(c: &dyn Curve) -> Result<GeomBSplineCurve, ConvertError> {
    let poles = c.bspline_poles().ok_or(ConvertError::DomainError)?;
    let knots = c.bspline_knots().ok_or(ConvertError::DomainError)?;
    let degree = c.nurbs_degree().ok_or(ConvertError::DomainError)?;
    Ok(GeomBSplineCurve::from_parts(
        poles.to_vec(),
        c.bspline_weights().map(|w| w.to_vec()),
        knots.to_vec(),
        degree,
        c.is_periodic(),
    ))
}

/// `GeomConvert::CurveToBSplineCurve(C, Parameterisation)` (`GeomConvert.cxx:163-430`);
/// `Parameterisation` defaults to `Convert_TgtThetaOver2`
/// (`GeomConvert.hxx:270-272`).
///
/// The two `Geom_OffsetCurve` arms (`:340-354`, `:436-450`) approximate the
/// (possibly trimmed) curve through `GeomConvert_ApproxCurve(C, 1e-4, C2, 16,
/// 14)` and raise `Standard_ConstructionError` when there is no result.
///
/// The other arms are ported: the trimmed line/conic arms (including the
/// `RationalC1` + `U2 - U1 >= 6` split-and-stitch, `:224-242`, `:262-280`),
/// the trimmed Bezier arm (`:300-321`, its rational sub-branch is noted
/// in place) and the trimmed/non-trimmed B-spline copy + `Segment` arms
/// (`:322-339`, `:431-434`).
///
/// The non-trimmed `Geom_Circle`/`Geom_Ellipse` arms (`:363-408`) end with
/// `TheCurve->SetPeriodic()` (`:378`, `:383`, `:406`); the periodic
/// `Geom_BSplineCurve` representation they need is ported
/// (`GeomBSplineCurve::set_periodic`, `Geom_BSplineCurve.cxx:777-815`) —
/// board card **R2-21**.
///
/// An unrecognised curve type throws `Standard_DomainError("No such curve")`
/// (`:355-358`, `:451-454`) — returned as [`ConvertError::DomainError`].
pub fn curve_to_bspline_curve(
    c: &dyn Curve,
    parameterisation: ParameterisationType,
) -> Result<GeomBSplineCurve, ConvertError> {
    if c.is_geom_trimmed() {
        // `Curv = Ctrim->BasisCurve()`, `U1 = FirstParameter()`,
        // `U2 = LastParameter()` (`GeomConvert.cxx:171-175`).
        let (basis, _b_first, _b_last) = c.untrimmed_basis().ok_or(ConvertError::DomainError)?;
        let mut u1 = c.first_parameter();
        let mut u2 = c.last_parameter();
        // `cxx:177-189`: for a non-periodic basis the range is clamped, so that
        // `BS->Segment` cannot raise.
        if !basis.is_periodic() {
            if u1 < basis.first_parameter() {
                u1 = basis.first_parameter();
            }
            if u2 > basis.last_parameter() {
                u2 = basis.last_parameter();
            }
        }

        if basis.gp_line().is_some() {
            // `cxx:191-206`: `Poles = {StartPoint, EndPoint}` (the trim's own
            // values), `Knots = {FirstParameter, LastParameter}`, both with
            // multiplicity 2, `Degree = 1`.
            let poles = vec![c.d0(c.first_parameter()), c.d0(c.last_parameter())];
            let knots = vec![c.first_parameter(), c.last_parameter()];
            let flat = knot_sequence(&knots, &[2, 2], 1);
            return GeomBSplineCurve::new(poles, flat, 1)
                .map_err(|_| ConvertError::ConstructionError);
        }

        if let Some(circ) = basis.gp_circ() {
            // `cxx:208-244`: the 2D conic is centered at the origin with the same
            // radius; `BSplineCurveBuilder` moves it back through
            // `TheConic->Position()`.
            let c2d = GpCirc2d::new(GpAx22d::standard(), circ.radius());
            if parameterisation != ParameterisationType::RationalC1 || (u2 - u1) < 6.0 {
                let convert = circle_to_bspline_curve_range(&c2d, u1, u2, parameterisation)?;
                return bspline_curve_builder(&circ.position(), &convert);
            }
            // `cxx:225-242`: `U2 - U1 >= 6` splits the circle in two halves and
            // stitches them with `GeomConvert_CompCurveToBSplineCurve` to avoid
            // the numerical overflow at `U2 - U1 ~ 2*PI`.
            let u_med = (u1 + u2) * 0.5;
            let convert1 = circle_to_bspline_curve_range(&c2d, u1, u_med, parameterisation)?;
            let curve1 = bspline_curve_builder(&circ.position(), &convert1)?;
            let convert2 = circle_to_bspline_curve_range(&c2d, u_med, u2, parameterisation)?;
            let curve2 = bspline_curve_builder(&circ.position(), &convert2)?;
            // `cxx:237-241`: `Add(TheCurve2, Precision::PConfusion(), true)`
            // (After = true, WithRatio = true, MinM = 0).
            let mut cctbspl = CompCurveToBSplineCurve::from_bspline(curve1, parameterisation);
            cctbspl.add(&curve2, occt_core::precision::PCONFUSION, true, true, 0)?;
            return cctbspl
                .into_curve()
                .ok_or(ConvertError::ConstructionError);
        }

        if let Some(elips) = basis.gp_ellipse() {
            // `cxx:246-282`.
            let e2d = GpElips2d::new(
                GpAx22d::standard(),
                elips.major_radius,
                elips.minor_radius,
            );
            if parameterisation != ParameterisationType::RationalC1 || (u2 - u1) < 6.0 {
                let convert = ellipse_to_bspline_curve_range(&e2d, u1, u2, parameterisation)?;
                return bspline_curve_builder(elips.position(), &convert);
            }
            // `cxx:262-280`: the same split-and-stitch as the circle arm.
            let u_med = (u1 + u2) * 0.5;
            let convert1 = ellipse_to_bspline_curve_range(&e2d, u1, u_med, parameterisation)?;
            let curve1 = bspline_curve_builder(elips.position(), &convert1)?;
            let convert2 = ellipse_to_bspline_curve_range(&e2d, u_med, u2, parameterisation)?;
            let curve2 = bspline_curve_builder(elips.position(), &convert2)?;
            // `cxx:275-279`: `Add(TheCurve2, Precision::PConfusion(), true)`.
            let mut cctbspl = CompCurveToBSplineCurve::from_bspline(curve1, parameterisation);
            cctbspl.add(&curve2, occt_core::precision::PCONFUSION, true, true, 0)?;
            return cctbspl
                .into_curve()
                .ok_or(ConvertError::ConstructionError);
        }

        if let Some(hypr) = basis.gp_hyperbola() {
            // `cxx:284-290`.
            let h2d = GpHypr2d::new(GpAx22d::standard(), hypr.major_radius, hypr.minor_radius);
            let convert = hyperbola_to_bspline_curve(&h2d, u1, u2)?;
            return bspline_curve_builder(hypr.position(), &convert);
        }

        if let Some(parab) = basis.gp_parabola() {
            // `cxx:292-298`.
            let p2d = GpParab2d::new(GpAx22d::standard(), parab.focal);
            let convert = parabola_to_bspline_curve(&p2d, u1, u2)?;
            return bspline_curve_builder(parab.position(), &convert);
        }

        // `cxx:300-321` (`Geom_BezierCurve`: `CBez->Segment(U1, U2)`) needs
        // `Geom_BezierCurve::Segment` (`PLib::Trimming` + `PLib::CoefficientsPoles`)
        // — ported (`bezier_curve.rs::segment`, OCCT `Geom_BezierCurve.cxx:388-425`).
        // `cxx:300-321`: the copy is trimmed and re-wrapped as a B-spline over
        // `[0, 1]` with multiplicities `{Degree+1, Degree+1}`. Rational uses
        // `WeightsArray` (`GeomConvert.cxx:313-321`).
        if let Some(poles) = basis.bezier_poles() {
            let mut bez = if let Some(w) = basis.bezier_weights() {
                GeomBezierCurve::with_weights(poles.to_vec(), w.to_vec())
                    .map_err(|_| ConvertError::ConstructionError)?
            } else {
                GeomBezierCurve::new(poles.to_vec())
                    .map_err(|_| ConvertError::ConstructionError)?
            };
            bez.segment(u1, u2);
            let degree = bez.degree();
            let m = (degree + 1) as i32;
            let flat = knot_sequence(&[0.0, 1.0], &[m, m], degree as i32);
            return if let Some(w) = bez.weights.as_ref() {
                GeomBSplineCurve::rational(bez.poles.clone(), w.clone(), flat, degree)
                    .map_err(|_| ConvertError::ConstructionError)
            } else {
                GeomBSplineCurve::new(bez.poles.clone(), flat, degree)
                    .map_err(|_| ConvertError::ConstructionError)
            };
        }
        // `cxx:322-339` (`Geom_BSplineCurve`): the basis is copied, its range is
        // folded into the period (`ElCLib::AdjustPeriodic`), a full-period trim
        // drops the periodic representation (`SetNotPeriodic`) and the copy is
        // then cut with `Segment(U1, U2)`.
        if basis.bspline_poles().is_some() {
            let mut bs = GeomBSplineCurve::from_parts(
                basis.bspline_poles().ok_or(ConvertError::DomainError)?.to_vec(),
                basis.bspline_weights().map(|w| w.to_vec()),
                basis.bspline_knots().ok_or(ConvertError::DomainError)?.to_vec(),
                basis.nurbs_degree().ok_or(ConvertError::DomainError)?,
                basis.is_periodic(),
            );
            if bs.is_periodic() {
                let (uf, ul) = (bs.first_parameter(), bs.last_parameter());
                occt_core::elib::clib2d::adjust_periodic(
                    uf,
                    ul,
                    occt_core::precision::CONFUSION,
                    &mut u1,
                    &mut u2,
                );
                if (u1 - uf).abs() <= occt_core::precision::CONFUSION
                    && (u2 - ul).abs() <= occt_core::precision::CONFUSION
                {
                    bs.set_not_periodic();
                }
            }
            // `theTolerance` defaults to `Precision::PConfusion()`
            // (`Geom_BSplineCurve.hxx:322-324`); the throw is
            // `Standard_DomainError`.
            bs.segment(u1, u2, occt_core::precision::PCONFUSION)
                .map_err(|_| ConvertError::DomainError)?;
            return Ok(bs);
        }
        // `cxx:340-354`: `GeomConvert_ApproxCurve(C, 1e-4, C2, 16, 14)` on the
        // original (possibly trimmed) curve `C`; a missing result raises
        // `Standard_ConstructionError` (`:352`).
        if basis.offset_curve().is_some() {
            let appr = GeomConvertApproxCurve::new(c, 1.0e-4, Shape::C2, 16, 14);
            return match appr.curve() {
                Some(curve) => Ok(curve.clone()),
                None => Err(ConvertError::ConstructionError),
            };
        }
        // `throw Standard_DomainError("No such curve")` (`cxx:355-358`).
        return Err(ConvertError::DomainError);
    }

    // Non-trimmed arm (`cxx:361-455`).
    if let Some(elips) = c.gp_ellipse() {
        // `cxx:363-385`: `Convert_EllipseToBSplineCurve(E2d, Parameterisation)`
        // followed by `TheCurve->SetPeriodic()`.
        let e2d = GpElips2d::new(GpAx22d::standard(), elips.major_radius, elips.minor_radius);
        let convert = ellipse_to_bspline_curve(&e2d, parameterisation)?;
        let mut curve = bspline_curve_builder(elips.position(), &convert)?;
        curve.set_periodic();
        return Ok(curve);
    }

    if let Some(circ) = c.gp_circ() {
        // `cxx:387-408`.
        let c2d = GpCirc2d::new(GpAx22d::standard(), circ.radius());
        let convert = circle_to_bspline_curve(&c2d, parameterisation)?;
        let mut curve = bspline_curve_builder(&circ.position(), &convert)?;
        curve.set_periodic();
        return Ok(curve);
    }

    if let Some(poles) = c.bezier_poles() {
        // `cxx:410-429`: a Bezier is the clamped B-spline with the same poles,
        // knots `{0, 1}` of multiplicity `degree + 1`; a rational Bezier passes
        // `WeightsArray()` (`cxx:421-424`).
        let degree = poles
            .len()
            .checked_sub(1)
            .ok_or(ConvertError::ConstructionError)?;
        let mut knots = vec![0.0f64; degree + 1];
        knots.extend(std::iter::repeat(1.0).take(degree + 1));
        return if let Some(w) = c.bezier_weights() {
            GeomBSplineCurve::rational(poles.to_vec(), w.to_vec(), knots, degree)
                .map_err(|_| ConvertError::ConstructionError)
        } else {
            GeomBSplineCurve::new(poles.to_vec(), knots, degree)
                .map_err(|_| ConvertError::ConstructionError)
        };
    }

    if is_bspline(c) {
        // `cxx:431-434`: `TheCurve = C->Copy()`.
        return bspline_copy(c);
    }

    // `cxx:436-450`: `GeomConvert_ApproxCurve(C, 1e-4, C2, 16, 14)`; a missing
    // result raises `Standard_ConstructionError` (`:448`).
    if c.offset_curve().is_some() {
        let appr = GeomConvertApproxCurve::new(c, 1.0e-4, Shape::C2, 16, 14);
        return match appr.curve() {
            Some(curve) => Ok(curve.clone()),
            None => Err(ConvertError::ConstructionError),
        };
    }

    // `throw Standard_DomainError("No such curve")` (`cxx:451-454`).
    Err(ConvertError::DomainError)
}
