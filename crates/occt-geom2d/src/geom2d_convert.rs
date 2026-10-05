//! `Geom2dConvert` subset used by `GeomLib::SameRange` and
//! `ShapeBuild_Edge::TransformPCurve` (`Geom2dConvert.cxx`).
//!
//! [`curve_to_bspline_curve`] is the parameterisation-less entry the port
//! already used for the Bezier branch of `TransformPCurve`.
//!
//! [`curve_to_bspline_curve_bspl`] is the full port of
//! `Geom2dConvert::CurveToBSplineCurve(C, Parameterisation)`
//! (`Geom2dConvert.cxx:181-449`):
//!
//! * trimmed input (`:189-374`): take `BasisCurve()` and `[First, Last]`, clamp
//!   the range for a non-periodic basis (`:199-209`), then dispatch on the basis
//!   kind through [`trimmed_curve_to_bspline_curve`];
//! * plain ellipse / circle (`:379-397`): the periodic
//!   `Convert_{Ellipse,Circle}ToBSplineCurve` + `SetPeriodic()`;
//! * plain Bezier / B-spline (`:399-423`).
//!
//! The conic arms build the B-spline through
//! `occt-core::convert::{circle,ellipse,hyperbola,parabola}_to_bspline_curve`
//! (`Convert_*ToBSplineCurve`) exactly as the 3D port
//! (`occt-geom/src/convert_bspl.rs`) does, so the 2D and 3D conversions stay
//! isomorphic. The conic pole coordinates are already in the conic's own frame
//! (`Convert_*ToBSplineCurve` applies `SetTransformation(conic.XAxis(), OX2d)`),
//! so no extra transform is needed.
//!
//! UNPORTED: `Geom2dConvert_ApproxCurve` (the `Geom2d_OffsetCurve` arms
//! `:353-368`, `:425-440`). Callers hitting them get `None`.

use crate::bezier_curve::Geom2dBezierCurve;
use crate::bspline_curve::Geom2dBSplineCurve;
use crate::comp_curve_to_bspline::CompCurveToBSplineCurve;
use crate::curve::Curve2d;
use occt_core::bspl::banded_interp::knot_sequence;
use occt_core::bspl::knots::knot_sequence_periodic;
use occt_core::convert::{
    circle_to_bspline_curve, circle_to_bspline_curve_range, ellipse_to_bspline_curve,
    ellipse_to_bspline_curve_range, hyperbola_to_bspline_curve, parabola_to_bspline_curve,
    ConicToBSplineCurve, ParameterisationType,
};

/// Default parameterisation of `Geom2dConvert::CurveToBSplineCurve`
/// (`Geom2dConvert.hxx:169-171`).
pub const DEFAULT_PARAMETERISATION: ParameterisationType = ParameterisationType::TgtThetaOver2;

/// Build a 2D B-spline from a `Convert_ConicToBSplineCurve` result, the way
/// OCCT's `Geom2dConvert::BSplineCurveBuilder` (`Geom2dConvert.cxx:70-98`) does:
/// `Geom2d_BSplineCurve(Poles, Weights, Knots, Mults, Degree, IsPeriodic)`
/// followed by `updateKnots()` (the flat sequence, period-extended when
/// periodic).
fn build_conic_bspline2d(convert: &ConicToBSplineCurve) -> Option<Geom2dBSplineCurve> {
    let poles = convert.poles();
    let xs: Vec<f64> = poles.iter().map(|p| p.x()).collect();
    let ys: Vec<f64> = poles.iter().map(|p| p.y()).collect();
    let degree = convert.degree();
    let flat = if convert.is_periodic() {
        knot_sequence_periodic(convert.knots(), convert.multiplicities(), degree as i32)
    } else {
        knot_sequence(convert.knots(), convert.multiplicities(), degree as i32)
    };
    Geom2dBSplineCurve::from_flat(
        xs,
        ys,
        Some(convert.weights().to_vec()),
        flat,
        degree,
        convert.is_periodic(),
    )
    .ok()
}

/// `Geom2dConvert::CurveToBSplineCurve(C)` (`Geom2dConvert.cxx`), the exact
/// non-conic cases: Bezier -> BSpline with the identical poles, BSpline ->
/// itself. Returns `None` for any other curve type (UNPORTED).
pub fn curve_to_bspline_curve(c: &dyn Curve2d) -> Option<Box<dyn Curve2d>> {
    if c.is_bspline2d() {
        return Some(c.clone_dyn());
    }
    if c.is_bezier2d() {
        let poles = c.poles2d()?;
        let degree = poles.len().checked_sub(1)?;
        let xs: Vec<f64> = poles.iter().map(|p| p.x()).collect();
        let ys: Vec<f64> = poles.iter().map(|p| p.y()).collect();
        let mut knots = vec![0.0f64; degree + 1];
        knots.extend(std::iter::repeat(1.0f64).take(degree + 1));
        // `CBez->IsRational()` passes `WeightsArray()` (`:411-414`).
        let bs = match c.bezier_weights2d() {
            Some(w) => Geom2dBSplineCurve::rational(xs, ys, w.to_vec(), knots, degree).ok()?,
            None => Geom2dBSplineCurve::new(xs, ys, knots, degree).ok()?,
        };
        return Some(Box::new(bs));
    }
    None
}

/// `Geom2dConvert::CurveToBSplineCurve(C, Parameterisation)`
/// (`Geom2dConvert.cxx:181-449`).
pub fn curve_to_bspline_curve_bspl(
    c: &dyn Curve2d,
    parameterisation: ParameterisationType,
) -> Option<Geom2dBSplineCurve> {
    if let Some(basis) = c.trimmed_basis() {
        // `cxx:193-209`: the trim's own range, clamped for a non-periodic basis.
        let mut u1 = c.first_parameter();
        let mut u2 = c.last_parameter();
        if !basis.is_periodic() {
            if u1 < basis.first_parameter() {
                u1 = basis.first_parameter();
            }
            if u2 > basis.last_parameter() {
                u2 = basis.last_parameter();
            }
        }
        return trimmed_curve_to_bspline_curve(basis, u1, u2, parameterisation);
    }

    // Non-trimmed arm (`cxx:376-446`).
    if let Some(elips) = c.gp_elips2d() {
        // `cxx:379-387`.
        let convert = ellipse_to_bspline_curve(&elips, parameterisation).ok()?;
        let mut bs = build_conic_bspline2d(&convert)?;
        bs.set_periodic();
        return Some(bs);
    }
    if let Some(circ) = c.gp_circ2d() {
        // `cxx:389-397`.
        let convert = circle_to_bspline_curve(&circ, parameterisation).ok()?;
        let mut bs = build_conic_bspline2d(&convert)?;
        bs.set_periodic();
        return Some(bs);
    }
    if c.is_bezier2d() {
        // `cxx:399-419`: the same poles, the clamped `[0, 1]` knot vector.
        let poles = c.poles2d()?;
        let degree = poles.len().checked_sub(1)?;
        let xs: Vec<f64> = poles.iter().map(|p| p.x()).collect();
        let ys: Vec<f64> = poles.iter().map(|p| p.y()).collect();
        let mut knots = vec![0.0f64; degree + 1];
        knots.extend(std::iter::repeat(1.0f64).take(degree + 1));
        return match c.bezier_weights2d() {
            Some(w) => Geom2dBSplineCurve::rational(xs, ys, w.to_vec(), knots, degree).ok(),
            None => Geom2dBSplineCurve::new(xs, ys, knots, degree).ok(),
        };
    }
    if c.is_bspline2d() {
        // `cxx:420-423`: `TheCurve = C->Copy()`.
        return c.bspline_copy2d();
    }
    // `cxx:425-440` (Offset -> `Geom2dConvert_ApproxCurve`) is UNPORTED;
    // `cxx:442-444` else throws `Standard_DomainError`.
    None
}

/// The trimmed-curve arms of `Geom2dConvert::CurveToBSplineCurve`
/// (`Geom2dConvert.cxx:211-373`): the caller supplies the already-resolved
/// `basis` (the trim's `BasisCurve()`) and its clamped `[u1, u2]` window.
pub fn trimmed_curve_to_bspline_curve(
    basis: &dyn Curve2d,
    u1: f64,
    u2: f64,
    parameterisation: ParameterisationType,
) -> Option<Geom2dBSplineCurve> {
    // `cxx:211-226`: a 2-pole, degree-1 B-spline through the trim's endpoints.
    if basis.is_line() {
        let p0 = basis.d0(u1);
        let p1 = basis.d0(u2);
        let flat = knot_sequence(&[u1, u2], &[2, 2], 1);
        return Geom2dBSplineCurve::new(vec![p0.x(), p1.x()], vec![p0.y(), p1.y()], flat, 1).ok();
    }

    // `cxx:228-264`.
    if let Some(circ) = basis.gp_circ2d() {
        if parameterisation != ParameterisationType::RationalC1 || (u2 - u1) < 6.0 {
            let convert = circle_to_bspline_curve_range(&circ, u1, u2, parameterisation).ok()?;
            return build_conic_bspline2d(&convert);
        }
        // `cxx:244-262`: split the circle to avoid numerical overflow when
        // `U2 - U1 =~ 2*PI`.
        let u_med = (u1 + u2) * 0.5;
        let convert1 = circle_to_bspline_curve_range(&circ, u1, u_med, parameterisation).ok()?;
        let curve1 = build_conic_bspline2d(&convert1)?;
        let convert2 = circle_to_bspline_curve_range(&circ, u_med, u2, parameterisation).ok()?;
        let curve2 = build_conic_bspline2d(&convert2)?;
        let mut cctbspl = CompCurveToBSplineCurve::from_bspline(curve1, parameterisation);
        let _ = cctbspl.add(&curve2, occt_core::precision::PCONFUSION, true);
        return cctbspl.into_curve();
    }

    // `cxx:266-303`.
    if let Some(elips) = basis.gp_elips2d() {
        if parameterisation != ParameterisationType::RationalC1 || (u2 - u1) < 6.0 {
            let convert = ellipse_to_bspline_curve_range(&elips, u1, u2, parameterisation).ok()?;
            return build_conic_bspline2d(&convert);
        }
        // `cxx:283-301`: split the ellipse (`Convert_EllipseToBSplineCurve`).
        let u_med = (u1 + u2) * 0.5;
        let convert1 = ellipse_to_bspline_curve_range(&elips, u1, u_med, parameterisation).ok()?;
        let curve1 = build_conic_bspline2d(&convert1)?;
        let convert2 = ellipse_to_bspline_curve_range(&elips, u_med, u2, parameterisation).ok()?;
        let curve2 = build_conic_bspline2d(&convert2)?;
        let mut cctbspl = CompCurveToBSplineCurve::from_bspline(curve1, parameterisation);
        let _ = cctbspl.add(&curve2, occt_core::precision::PCONFUSION, true);
        return cctbspl.into_curve();
    }

    // `cxx:305-312`.
    if let Some(hypr) = basis.gp_hypr2d() {
        let convert = hyperbola_to_bspline_curve(&hypr, u1, u2).ok()?;
        return build_conic_bspline2d(&convert);
    }

    // `cxx:314-321`.
    if let Some(parab) = basis.gp_parab2d() {
        let convert = parabola_to_bspline_curve(&parab, u1, u2).ok()?;
        return build_conic_bspline2d(&convert);
    }

    // `cxx:323-345`: copy, `Segment(U1, U2)`, then a clamped `[0,1]` B-spline.
    if basis.is_bezier2d() {
        let poles = basis.poles2d()?;
        let mut cbez = match basis.bezier_weights2d() {
            Some(w) => Geom2dBezierCurve::rational(poles, w.to_vec()).ok()?,
            None => Geom2dBezierCurve::new(poles).ok()?,
        };
        cbez.segment(u1, u2);
        let degree = cbez.degree();
        let mut knots = vec![0.0f64; degree + 1];
        knots.extend(std::iter::repeat(1.0f64).take(degree + 1));
        let xs: Vec<f64> = cbez.poles.iter().map(|p| p.x()).collect();
        let ys: Vec<f64> = cbez.poles.iter().map(|p| p.y()).collect();
        return match cbez.weights {
            Some(w) => Geom2dBSplineCurve::rational(xs, ys, w, knots, degree).ok(),
            None => Geom2dBSplineCurve::new(xs, ys, knots, degree).ok(),
        };
    }

    // `cxx:347-351`: copy and `Segment(U1, U2)`.
    if basis.is_bspline2d() {
        let mut bs = basis.bspline_copy2d()?;
        bs.segment(u1, u2, occt_core::precision::PCONFUSION).ok()?;
        return Some(bs);
    }

    // `cxx:353-368` (Offset -> `Geom2dConvert_ApproxCurve`) is UNPORTED.
    None
}
