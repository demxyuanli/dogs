//! `ChFi3d_ReparamPcurv` — reparametrize a B-spline p-curve onto a requested
//! range.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_0.cxx:1467-1507`
//! (OCCT 8.0.0):
//!
//! ```text
//! void ChFi3d_ReparamPcurv(const double Uf, const double Ul, handle<Geom2d_Curve>& Pcurv)
//! ```
//!
//! The function is the `PerformSetOfKGen` fix-up for the p-curves that
//! `GeomFill_ConstrainedFilling` reparametrized onto `[0, 1]`: the B-spline is
//! first restricted to the p-curve's *own* range (`Geom2d_BSplineCurve::
//! Segment`), then its knots are rescaled onto `[Uf, Ul]` and a new
//! `Geom2d_BSplineCurve` is built from the old poles and multiplicities
//! (`cxx:1497-1504`).
//!
//! Translation notes:
//! - `occ::down_cast<Geom2d_BSplineCurve>(basis)`, where `basis` is the
//!   `BasisCurve()` of a `Geom2d_TrimmedCurve` and the curve itself otherwise,
//!   is [`Curve2d::bspline_copy2d`] — it already unwraps exactly one trimmed
//!   layer and answers `None` for every other type, which is OCCT's null
//!   `pc`,
//! - OCCT writes `Pcurv = pc` on every path, so a trimmed B-spline p-curve is
//!   replaced by its (possibly modified) basis, exactly as `cxx:1505` does,
//! - `pc->Knots()` is the flat knot sequence (`myKnots`, not the period
//!   extended `myFlatKnots`), so `BSplCLib::Reparametrize` runs over
//!   `Geom2dBSplineCurve::knots` and the rebuild hands the same flat sequence
//!   to `Geom2dBSplineCurve::from_flat` with `Periodic = false` — the
//!   `Geom2d_BSplineCurve(Poles, Knots, Mults, Degree)` constructor
//!   (`Geom2d_BSplineCurve.cxx:136-166`) leaves `myPeriodic` at its `false`
//!   default,
//! - `Mults()` is only read to build the constructor argument; the port's
//!   `from_flat` takes the already-flat sequence, so it is not needed. Its
//!   validation counterpart is `CheckCurveData(Knots, Mults, Degree,
//!   Periodic)`, which `from_flat` reproduces through `knots::check_degree`,
//! - `Segment` / the constructor raise `Standard_ConstructionError` where the
//!   port returns `Err`.

use std::sync::Arc;

use occt_core::bspl::knots::reparametrize;
use occt_core::precision::PCONFUSION;
use occt_geom2d::bspline_curve::Geom2dBSplineCurve;
use occt_geom2d::curve::Curve2d;

/// `ChFi3d_ReparamPcurv(Uf, Ul, Pcurv)` (`cxx:1471-1507`).
///
/// `pcurv` is OCCT's in-out handle: `None` answers immediately (`cxx:1473-1476`).
pub fn ch_fi3d_reparam_pcurv(
    uf: f64,
    ul: f64,
    pcurv: &mut Option<Arc<dyn Curve2d>>,
) -> Result<(), String> {
    let Some(curve) = pcurv.as_ref() else {
        return Ok(());
    };
    // `cxx:1478-1479`.
    let upcf = curve.first_parameter();
    let upcl = curve.last_parameter();

    // `cxx:1480-1489`: `basis` then `down_cast<Geom2d_BSplineCurve>`, and a
    // null cast returns.
    let Some(mut pc) = curve.bspline_copy2d() else {
        return Ok(());
    };

    // `cxx:1491-1495`.
    if (upcf - pc.first_parameter()).abs() > PCONFUSION
        || (upcl - pc.last_parameter()).abs() > PCONFUSION
    {
        pc.segment(upcf, upcl, PCONFUSION)
            .map_err(|_| "Geom2d_BSplineCurve::Segment (ChFi3d_Builder_0.cxx:1494)".to_string())?;
    }

    // `cxx:1496-1504`.
    if (uf - pc.first_parameter()).abs() > PCONFUSION
        || (ul - pc.last_parameter()).abs() > PCONFUSION
    {
        let mut kn = pc.knots.clone();
        reparametrize(uf, ul, &mut kn);
        let degree = pc.degree();
        pc = Geom2dBSplineCurve::from_flat(
            pc.xs.clone(),
            pc.ys.clone(),
            pc.weights.clone(),
            kn,
            degree,
            false,
        )
        .map_err(|e| {
            format!("Geom2d_BSplineCurve pole/knot mismatch (ChFi3d_Builder_0.cxx:1503): {e}")
        })?;
    }

    // `cxx:1505`: `Pcurv = pc;`
    *pcurv = Some(Arc::new(pc));
    Ok(())
}
