//! `ChFiKPart_ComputeData_Fcts` — the file statics every `ChFiKPart` surface
//! builder shares.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFiKPart/ChFiKPart_ComputeData_Fcts.cxx:34-155`
//! (OCCT 8.0.0).
//!
//! The ported functions are `ChFiKPart_InPeriod`, `ChFiKPart_PCurve`,
//! `ChFiKPart_IndexCurveInDS`, `ChFiKPart_IndexSurfaceInDS` and
//! `ChFiKPart_ProjPC`. The last one keeps OCCT's three-way switch on the
//! projected-curve type; its `GeomAbs_BezierCurve` / `GeomAbs_BSplineCurve`
//! arms have no producer in the port yet, because the `ProjLib` projectors
//! that call `SetBezier` / `SetBSpline` are not ported either. When the
//! projector answers not-done, OCCT falls back to `ProjLib_ComputeApprox`
//! (`ProjLib_ProjectedCurve.cxx:717-785`, 1411 lines) — **UNPORTED**, which is
//! the one answer of this function that is not a translation.

use std::sync::Arc;

use occt_core::gp::{GpLin2d, GpPnt2d};
use occt_core::kernel::geomabs::{CurveType, SurfaceType};
use occt_geom::curve::Curve as GeomCurve;
use occt_geom::projlib::{surface_type, ProjectedCurve};
use occt_geom::surface::Surface as GeomSurface;
use occt_geom2d::bezier_curve::Geom2dBezierCurve;
use occt_geom2d::bspline_curve::Geom2dBSplineCurve;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::line::Geom2dLine;

use crate::topopebrep::ds::{Curve, DataStructure, Surface};

/// `ChFiKPart_InPeriod(U, UFirst, ULast, Eps)` (`cxx:36-52`): shift `U` by
/// whole periods until it lies in `[UFirst, ULast]`, then clamp it up to
/// `UFirst` when it fell below.
pub fn ch_fi_kpart_in_period(u: f64, u_first: f64, u_last: f64, eps: f64) -> f64 {
    let period = u_last - u_first;
    let mut u = u;
    while eps < u_first - u {
        u += period;
    }
    while eps > u_last - u {
        u -= period;
    }
    if u < u_first {
        u = u_first;
    }
    u
}

/// `ChFiKPart_PCurve(UV1, UV2, Pardeb, Parfin)` (`cxx:60-75`): the straight
/// line between two `(u, v)` points as a degree-1 B-spline with knots
/// `{Pardeb, Pardeb, Parfin, Parfin}` — OCCT builds it as
/// `Geom2d_BSplineCurve(p, k, m, 1)` with `k = {Pardeb, Parfin}` and
/// `m = {2, 2}` (`cxx:66-74`), which `updateKnots` expands to exactly that
/// flat vector.
pub fn ch_fi_kpart_p_curve(
    uv1: GpPnt2d,
    uv2: GpPnt2d,
    pardeb: f64,
    parfin: f64,
) -> Result<Geom2dBSplineCurve, String> {
    Geom2dBSplineCurve::new(
        vec![uv1.x(), uv2.x()],
        vec![uv1.y(), uv2.y()],
        vec![pardeb, pardeb, parfin, parfin],
        1,
    )
    .map_err(|e| e.to_string())
}

/// `ChFiKPart_IndexCurveInDS(C, DStr)` (`cxx:142-145`):
/// `DStr.AddCurve(TopOpeBRepDS_Curve(C, 0.))`. `C` is a handle, so OCCT also
/// accepts a null one — `ChFiKPart_ComputeData_FilPlnCyl.cxx:492` /
/// `_FilPlnCon.cxx:240` pass the still-null `GCircPln` of their sphere arm.
pub fn ch_fi_kpart_index_curve_in_ds(
    c: Option<Arc<dyn GeomCurve>>,
    ds: &mut DataStructure,
) -> i32 {
    let mut curve = Curve::with_curve(c, 0.0, false);
    ds.add_curve(&mut curve)
}

/// `ChFiKPart_IndexSurfaceInDS(S, DStr)` (`cxx:152-155`):
/// `DStr.AddSurface(TopOpeBRepDS_Surface(S, 0.))`.
pub fn ch_fi_kpart_index_surface_in_ds(s: Arc<dyn GeomSurface>, ds: &mut DataStructure) -> i32 {
    ds.add_surface(&Surface::with_tolerance(s, 0.0))
}

/// `ChFiKPart_ProjPC(const GeomAdaptor_Curve& Cg, const GeomAdaptor_Surface&
/// Sg, occ::handle<Geom2d_Curve>& Pcurv)` (`cxx:83-135`).
///
/// The `Pcurv` in/out handle becomes the answer: OCCT overwrites it in all
/// three arms of the switch and never reads the incoming value (its one caller
/// passes a line it has just built, `ChFiKPart_ComputeData_Sphere.cxx:141-147`).
///
/// `Cg` is `GeomAdaptor_Curve(C, 0., ang)` (`Sphere.cxx:163`), i.e. the circle
/// restricted to `[0, ang]` in the circle's own parameterisation; the port
/// passes the equivalent `GeomTrimmedCurveBasis`.
pub fn ch_fi_kpart_proj_pc(
    surface: Arc<dyn GeomSurface>,
    curve: Arc<dyn GeomCurve>,
) -> Result<Arc<dyn Curve2d>, String> {
    // `if (Sg.GetType() < GeomAbs_BezierSurface) { ... } else { throw
    // Standard_NotImplemented("approximate pcurve on the left surface"); }`
    // (`cxx:88-135`). `GeomAbs_SurfaceType` orders the five elementary surfaces
    // before `GeomAbs_BezierSurface` (`GeomAbs_SurfaceType.hxx:21-33`), which
    // the port's `SurfaceType` mirrors.
    let elementary = matches!(
        surface_type(surface.as_ref()),
        SurfaceType::Plane
            | SurfaceType::Cylinder
            | SurfaceType::Cone
            | SurfaceType::Sphere
            | SurfaceType::Torus
    );
    if !elementary {
        return Err("Standard_NotImplemented: approximate pcurve on the left surface".into());
    }

    // `ProjLib_ProjectedCurve Projc(HSg, HCg);` (`cxx:91-92`).
    let projc = ProjectedCurve::with_surface_curve(surface, curve);

    // OCCT reaches the switch with `myResult` always done: an undone projector
    // is replaced by the approximation of `ProjLib_ComputeApprox`
    // (`ProjLib_ProjectedCurve.cxx:717-785`), which is **UNPORTED** (1411
    // lines: the advanced `ProjLib_ComputeApproxOnPolarSurface` fit plus
    // `GeomLib::SameRange`). The port therefore cannot produce the B-spline
    // OCCT would hand to the `GeomAbs_BSplineCurve` arm.
    if !projc.is_done() {
        return Err(
            "UNPORTED ProjLib_ComputeApprox (ProjLib_ComputeApprox.cxx): \
             projection of the corner circle is not analytic"
                .into(),
        );
    }

    // `switch (Projc.GetType()) { case GeomAbs_Line: Pcurv = new
    // Geom2d_Line(Projc.Line()); break; case GeomAbs_BezierCurve: { ... }
    // break; case GeomAbs_BSplineCurve: { ... } break; default: throw
    // Standard_NotImplemented("failed approximation of the pcurve "); }`
    // (`cxx:93-128`).
    match projc.get_type() {
        CurveType::Line => {
            let l: GpLin2d = projc.line()?;
            Ok(Arc::new(Geom2dLine::from_pnt_dir(
                l.location(),
                *l.direction(),
            )))
        }
        CurveType::BezierCurve => {
            // `const NCollection_Array1<gp_Pnt2d>& TP = BezProjc->Poles();
            // if (BezProjc->IsRational()) { Pcurv = new
            // Geom2d_BezierCurve(TP, BezProjc->WeightsArray()); } else {
            // Pcurv = new Geom2d_BezierCurve(TP); }` (`cxx:99-109`).
            let bez = projc
                .bezier()
                .ok_or_else(|| "Standard_NullObject: ProjLib_Projector::Bezier".to_string())?;
            let out = match &bez.weights {
                Some(w) => Geom2dBezierCurve::rational(bez.poles.clone(), w.clone()),
                None => Geom2dBezierCurve::new(bez.poles.clone()),
            }
            .map_err(|e| e.to_string())?;
            Ok(Arc::new(out))
        }
        CurveType::BSplineCurve => {
            // `const NCollection_Array1<gp_Pnt2d>& TP = BspProjc->Poles(); const
            // NCollection_Array1<double>& TK = BspProjc->Knots(); const
            // NCollection_Array1<int>& TM = BspProjc->Multiplicities(); if
            // (BspProjc->IsRational()) { Pcurv = new Geom2d_BSplineCurve(TP,
            // BspProjc->WeightsArray(), TK, TM, BspProjc->Degree()); } else {
            // Pcurv = new Geom2d_BSplineCurve(TP, TK, TM, BspProjc->Degree());
            // }` (`cxx:111-124`). The OCCT constructors expand `TK`/`TM` into
            // the flat knot sequence and reject a pole/degree mismatch, which is
            // what `Geom2dBSplineCurve::new` / `::rational` do.
            let bs = projc
                .bspline()
                .ok_or_else(|| "Standard_NullObject: ProjLib_Projector::BSpline".to_string())?;
            let (knots, mults) = bs.distinct_knots_and_mults();
            let flat = occt_core::bspl::banded_interp::knot_sequence(&knots, &mults, bs.degree() as i32);
            let out = match &bs.weights {
                Some(w) => Geom2dBSplineCurve::rational(
                    bs.xs.clone(),
                    bs.ys.clone(),
                    w.clone(),
                    flat,
                    bs.degree(),
                ),
                None => Geom2dBSplineCurve::new(bs.xs.clone(), bs.ys.clone(), flat, bs.degree()),
            }
            .map_err(|e| e.to_string())?;
            Ok(Arc::new(out))
        }
        _ => Err("Standard_NotImplemented: failed approximation of the pcurve ".into()),
    }
}
