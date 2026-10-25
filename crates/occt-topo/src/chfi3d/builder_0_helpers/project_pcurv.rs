//! `ChFi3d_ProjectPCurv` — the p-curve of a 3D intersection curve, projected
//! onto an analytic surface.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_0.cxx:1513-1559`
//! (OCCT 8.0.0).
//!
//! `ProjLib_ProjectedCurve` answers a projector result rather than a
//! `Geom2d_Curve`; OCCT switches on `GetType()` and rebuilds the curve from the
//! matching `gp` object (`cxx:1527-1555`). The port reads the same state from
//! [`ProjectedCurve::result`] and rebuilds the same curve types.
//!
//! Notes:
//! - OCCT answers nothing (leaves `Pcurv` and `tolreached` untouched) on a
//!   Bezier / B-spline surface (`cxx:1520-1522`); the port states that as
//!   `Ok(None)`, and the caller keeps its `tolreached`,
//! - `default:` throws `Standard_NotImplemented`, which is the `Err`,
//! - the `GeomAdaptor_Curve` / `GeomAdaptor_Surface` arguments are the port's
//!   [`GeomAdaptorCurve`](crate::topopebrep::tool::adaptor::GeomAdaptorCurve)
//!   and `Arc<dyn Surface>`: `ProjLib_ProjectedCurve` only reads the surface.

use std::sync::Arc;

use occt_core::kernel::geomabs::CurveType;
use occt_geom::curve::Curve;
use occt_geom::projlib::ProjectedCurve;
use occt_geom::surface::Surface;
use occt_geom2d::circle::Geom2dCircle;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::ellipse::Geom2dEllipse;
use occt_geom2d::hyperbola::Geom2dHyperbola;
use occt_geom2d::line::Geom2dLine;
use occt_geom2d::parabola::Geom2dParabola;

use crate::geom_bnd_lib_surface3d::{surface_kind, GeomAbsSurfaceType};

/// `ChFi3d_ProjectPCurv(HCg, HSg, Pcurv, tol, tolreached)` (`cxx:1513-1559`).
///
/// Answers `Ok(Some(tolreached))` when a p-curve was built, `Ok(None)` when
/// OCCT leaves both out-parameters untouched (Bezier / B-spline surface).
pub fn ch_fi3d_project_pcurv(
    c3d: &Arc<dyn Curve>,
    s: &Arc<dyn Surface>,
    pcurv: &mut Option<Arc<dyn Curve2d>>,
    tol: f64,
) -> Result<Option<f64>, String> {
    // `if (HSg->GetType() != GeomAbs_BezierSurface &&
    // HSg->GetType() != GeomAbs_BSplineSurface)` (`cxx:1520-1522`).
    let kind = surface_kind(s.as_ref());
    if kind == GeomAbsSurfaceType::BezierSurface || kind == GeomAbsSurfaceType::BSplineSurface {
        return Ok(None);
    }

    // `ProjLib_ProjectedCurve Projc(HSg, HCg, tol);` (`cxx:1524`).
    let projc = ProjectedCurve::with_surface_curve_tol(Arc::clone(s), Arc::clone(c3d), tol);
    // `tolreached = Projc.GetTolerance();` (`cxx:1525`).
    let tolreached = projc.get_tolerance();

    // `switch (Projc.GetType())` (`cxx:1526-1556`).
    let curve: Arc<dyn Curve2d> = match projc.get_type() {
        // `Pcurv = new Geom2d_Line(Projc.Line());`
        CurveType::Line => Arc::new(Geom2dLine::new(*projc.result.lin.position())),
        // `Pcurv = new Geom2d_Circle(Projc.Circle());`
        CurveType::Circle => Arc::new(Geom2dCircle::new(projc.result.circ)),
        // `Pcurv = new Geom2d_Ellipse(Projc.Ellipse());`
        CurveType::Ellipse => Arc::new(Geom2dEllipse::new(projc.result.elips)),
        // `Pcurv = new Geom2d_Hyperbola(Projc.Hyperbola());`
        CurveType::Hyperbola => Arc::new(Geom2dHyperbola::new(projc.result.hypr)),
        // `Pcurv = new Geom2d_Parabola(Projc.Parabola());`
        CurveType::Parabola => Arc::new(Geom2dParabola::new(projc.result.parab)),
        // `Pcurv = Projc.Bezier();`
        CurveType::BezierCurve => match projc.result.bezier.as_ref() {
            Some(b) => b.clone() as Arc<dyn Curve2d>,
            None => {
                return Err("ChFi3d_ProjectPCurv: ProjLib result has no Bezier".to_string());
            }
        },
        // `Pcurv = Projc.BSpline();`
        CurveType::BSplineCurve => match projc.result.bspline.as_ref() {
            Some(b) => b.clone() as Arc<dyn Curve2d>,
            None => {
                return Err("ChFi3d_ProjectPCurv: ProjLib result has no BSpline".to_string());
            }
        },
        // `default: throw Standard_NotImplemented("echec approximation de la
        // pcurve ");` (`cxx:1553-1555`).
        other => {
            return Err(format!(
                "ChFi3d_ProjectPCurv: echec approximation de la pcurve ({other:?})"
            ));
        }
    };
    *pcurv = Some(curve);
    Ok(Some(tolreached))
}
