//! `GeomLib::SameRange` for 2D curves (`GeomLib.cxx:842-970`).
//!
//! The 3D overload is reached by `crates/occt-geom/src/curve_reparam.rs`; this
//! file ports the 2D one only, because `BRepFill_Sweep` needs it on pcurves
//! (`BRepFill_Sweep.cxx:1034-1042`, `cxx:1764-1769`).
//!
//! Translation notes:
//! - `occ::handle<Geom2d_Curve>&` out-parameter is the `Result` payload; the
//!   null handle OCCT leaves on the empty `cxx:906-922` branch is `Err`.
//! - `Geom2d_Line` the reference returns un-trimmed (only its location is
//!   shifted, `cxx:866-871`); so does the port.
//! - `Geom2dConvert::CurveToBSplineCurve` is
//!   [`curve_to_bspline_curve_bspl`] with the default `QuasiAngular`
//!   parameterisation (`Geom2dConvert.cxx:181-183`).
//! - `BSplCLib::Reparametrize` is `occt_core::bspl::knots::reparametrize`.

use std::sync::Arc;

use occt_core::convert::ParameterisationType;
use occt_core::gp::{GpPnt2d, GpTrsf2d, GpVec2d};
use occt_core::precision::PCONFUSION;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::geom2d_convert::curve_to_bspline_curve_bspl;
use occt_geom2d::line::Geom2dLine;
use occt_geom2d::trimmed::Geom2dTrimmedCurve;

/// `GeomLib::SameRange(Tolerance, CurvePtr, FirstOnCurve, LastOnCurve,
/// RequestedFirst, RequestedLast, NewCurvePtr)` (`GeomLib.cxx:842-970`).
pub fn same_range_2d(
    tolerance: f64,
    curve: &Arc<dyn Curve2d>,
    first_on_curve: f64,
    last_on_curve: f64,
    requested_first: f64,
    requested_last: f64,
) -> Result<Arc<dyn Curve2d>, String> {
    // `cxx:855-859`: the range already matches.
    if (last_on_curve - requested_last).abs() <= tolerance
        && (first_on_curve - requested_first).abs() <= tolerance
    {
        return Ok(curve.clone());
    }

    // `cxx:862-967`: the parametrisation length must at least be the same.
    if (last_on_curve - first_on_curve - requested_last + requested_first).abs() <= tolerance {
        // `cxx:864-871`.
        if let Some(line) = curve.gp_lin2d() {
            let du = first_on_curve - requested_first;
            let dir = *line.direction();
            let shift = GpVec2d::new(du * dir.x(), du * dir.y());
            let moved = line.translated_vec(&shift);
            return Ok(Arc::new(Geom2dLine::new(*moved.position())));
        }

        // `cxx:872-889`.
        if let Some(circ) = curve.gp_circ2d() {
            let p = circ.location();
            let du = if circ.is_direct() {
                first_on_curve - requested_first
            } else {
                requested_first - first_on_curve
            };
            let mut trsf = GpTrsf2d::identity();
            trsf.set_rotation(&p, du);
            let mut copy: Box<dyn Curve2d> = curve.clone_dyn();
            copy.transform(&trsf);
            return Ok(Arc::from(copy));
        }

        // `cxx:890-901`.
        if let Some(basis) = curve.trimmed_basis() {
            let new_basis = same_range_2d(
                tolerance,
                &Arc::from(basis.clone_dyn()),
                first_on_curve,
                last_on_curve,
                requested_first,
                requested_last,
            )?;
            return Ok(Arc::new(Geom2dTrimmedCurve::new(
                new_basis,
                requested_first,
                requested_last,
            )));
        }

        // `cxx:906-922`: the b-spline remake.
        if (last_on_curve - first_on_curve).abs() > PCONFUSION
            || (requested_last + requested_first).abs() > PCONFUSION
        {
            return remake_as_bspline_2d(
                curve,
                first_on_curve,
                last_on_curve,
                requested_first,
                requested_last,
            );
        }

        // `cxx:906-922` falls through with a null `NewCurvePtr` when both
        // guards fail.
        return Err("GeomLib::SameRange: null result curve".into());
    }

    // `cxx:924-969`: the result is segmented.
    let mut basis: Arc<dyn Curve2d> = curve.clone();
    let mut a_check: &dyn Curve2d = basis.as_ref();
    if let Some(trim_basis) = a_check.trimmed_basis() {
        basis = Arc::from(trim_basis.clone_dyn());
        a_check = basis.as_ref();
    }

    let tc: Arc<dyn Curve2d> = if a_check.is_periodic() {
        // `cxx:932-944`.
        if (last_on_curve - first_on_curve).abs() > PCONFUSION {
            Arc::new(Geom2dTrimmedCurve::new(
                curve.clone(),
                first_on_curve,
                last_on_curve,
            ))
        } else {
            Arc::new(Geom2dTrimmedCurve::new(
                curve.clone(),
                curve.first_parameter(),
                curve.last_parameter(),
            ))
        }
    } else {
        // `cxx:946-958`.
        let u_deb = curve.first_parameter().max(first_on_curve);
        let u_fin = curve.last_parameter().min(last_on_curve);
        if (u_fin - u_deb).abs() > PCONFUSION {
            Arc::new(Geom2dTrimmedCurve::new(curve.clone(), u_deb, u_fin))
        } else {
            Arc::new(Geom2dTrimmedCurve::new(
                curve.clone(),
                curve.first_parameter(),
                curve.last_parameter(),
            ))
        }
    };

    remake_trimmed_as_bspline_2d(&tc, requested_first, requested_last)
}

/// `cxx:909-921`: `Geom2d_TrimmedCurve` + `CurveToBSplineCurve` +
/// `BSplCLib::Reparametrize` + `SetKnots`.
fn remake_as_bspline_2d(
    curve: &Arc<dyn Curve2d>,
    first_on_curve: f64,
    last_on_curve: f64,
    requested_first: f64,
    requested_last: f64,
) -> Result<Arc<dyn Curve2d>, String> {
    let tc: Arc<dyn Curve2d> = Arc::new(Geom2dTrimmedCurve::new(
        curve.clone(),
        first_on_curve,
        last_on_curve,
    ));
    remake_trimmed_as_bspline_2d(&tc, requested_first, requested_last)
}

/// `cxx:915-921` and `cxx:961-968`: the shared b-spline remake.
fn remake_trimmed_as_bspline_2d(
    tc: &Arc<dyn Curve2d>,
    requested_first: f64,
    requested_last: f64,
) -> Result<Arc<dyn Curve2d>, String> {
    let mut bs = curve_to_bspline_curve_bspl(tc.as_ref(), ParameterisationType::QuasiAngular)
        .ok_or_else(|| "Geom2dConvert::CurveToBSplineCurve failed".to_string())?;
    let (mut knots, _mults) = bs.distinct_knots_and_mults();
    occt_core::bspl::knots::reparametrize(requested_first, requested_last, &mut knots);
    bs.set_knots(&knots)
        .map_err(|e| format!("Geom2d_BSplineCurve::SetKnots: {e}"))?;
    Ok(Arc::new(bs))
}

/// `gp_Pnt2d` of `curve` at `u`, used by the callers that re-trim.
pub fn curve2d_point(curve: &dyn Curve2d, u: f64) -> GpPnt2d {
    curve.d0(u)
}
