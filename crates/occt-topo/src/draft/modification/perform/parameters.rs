//! `Parameter` (`Draft_Modification_1.cxx:2180-2270`) and `SmartParameter`
//! (`:2274-2379`): the parameter of a point on a rebuilt curve, and the
//! re-projection of a vertex point onto the pair of new pcurves.
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | `Parameter(C, P, done)` | `_1.cxx:2180-2270` | [`parameter`] |
//! | `SmartParameter(Einf, EdgeTol, Pnt, sign, S1, S2)` | `_1.cxx:2274-2379` | [`smart_parameter`] |
//!
//! # Translation notes
//!
//! * `Parameter`'s `done` out-parameter keeps OCCT's meaning: `0` for an
//!   interior minimum, `-1` / `1` for the "no extrema in range, take the first
//!   / last parameter of the trimmed curve" arm (`cxx:2243-2257`).
//! * `Extrema_ExtPC` over `GeomAdaptor_Curve(C)` (`cxx:2221-2222`) is
//!   `extrema_pc::extrema_ext_pc_range_done` over the curve's own range. Its
//!   `TrimmedSquareDistances` (`cxx:2246`) are recomputed from the two
//!   endpoint evaluations, which is what that call answers when no extremum
//!   lies inside the range.
//! * `GeomAdaptor_Curve` reparameterizes a `Geom_TrimmedCurve` onto `[0, 1]`
//!   in this port (`trimmed.rs:17-38`) where OCCT keeps the basis parameters;
//!   the range handed to the extrema engine is therefore the trimmed curve's
//!   own, which is the same geometric domain.
//! * `SmartParameter`'s `Geom2dAPI_ProjectPointOnCurve` (`cxx:2320-2322`)
//!   is `geom2d_api::project_point_on_curve_2d`, whose `None` stands for
//!   `NbPoints() == 0`.
//! * `Approx_CurveOnSurface::Perform(MaxSeg, 10, GeomAbs_C1, false, false)`
//!   (`cxx:2366`) is `ApproxCurveOnSurface::perform(max_seg, 10, 2, false,
//!   false)`; `2` is this port's `GeomAbs_C1` encoding.

use std::sync::Arc;

use occt_core::convert::ParameterisationType;
use occt_core::elib::{clib, slib};
use occt_core::gp::{GpPnt, GpPnt2d};
use occt_core::precision::{epsilon, Precision, CONFUSION};
use occt_geom::approx::curve_on_surface::ApproxCurveOnSurface;
use occt_geom::curve::Curve;
use occt_geom::curve_on_surface::CurveOnSurface;
use occt_geom::extrema_pc::extrema_ext_pc_range_done;
use occt_geom::geom_api::project_point_on_surface;
use occt_geom::geom_proj_lib;
use occt_geom::projlib::{CompProjectedCurve, ProjectedCurve2d};
use occt_geom::Surface;
use occt_geom2d::bezier_curve::Geom2dBezierCurve;
use occt_geom2d::comp_curve_to_bspline::CompCurveToBSplineCurve2d;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::geom2d_api::project_point_on_curve_2d;
use occt_geom2d::geom2d_convert::curve_to_bspline_curve_bspl;

use crate::draft::edge_info::DraftEdgeInfo;

/// `static double Parameter(C, P, done)` (`_1.cxx:2180-2270`).
pub(super) fn parameter(c: &Arc<dyn Curve>, p: &GpPnt, done: &mut i32) -> f64 {
    *done = 0;
    // `cxx:2184-2189`: the analytic arms look at the basis of a trimmed curve.
    let mut cbase = c.clone();
    if c.is_geom_trimmed() {
        if let Some(basis) = c.trim_basis() {
            cbase = basis;
        }
    }

    // `cxx:2191-2199`
    if let Some(l) = cbase.gp_line() {
        return clib::line_parameter(&l, p);
    }
    if let Some(circ) = cbase.gp_circ() {
        let mut param = clib::circle_parameter(&circ.position(), p);
        if (2.0 * std::f64::consts::PI - param).abs() <= epsilon(2.0 * std::f64::consts::PI) {
            param = 0.0;
        }
        return param;
    }
    // `cxx:2202-2212`
    if let Some(elips) = cbase.gp_ellipse() {
        let mut param = clib::parameter_elips(&elips, p);
        if (2.0 * std::f64::consts::PI - param).abs() <= epsilon(2.0 * std::f64::consts::PI) {
            param = 0.0;
        }
        return param;
    }
    if let Some(parab) = cbase.gp_parabola() {
        return clib::parameter_parab(&parab, p);
    }
    if let Some(hypr) = cbase.gp_hyperbola() {
        return clib::parameter_hypr(&hypr, p);
    }

    // `cxx:2221-2258`: `GeomAdaptor_Curve TheCurve(C); Extrema_ExtPC myExtPC(P, TheCurve)`.
    let (sols, is_done) = extrema_ext_pc_range_done(
        c.as_ref(),
        p,
        c.first_parameter(),
        c.last_parameter(),
    );
    if !is_done {
        // `cxx:2224-2227`
        panic!("Standard_Failure: Draft_Modification_1::Parameter: ExtremaPC not done.");
    }
    let mut param;
    if !sols.is_empty() {
        // `cxx:2229-2241`
        let mut jmin = 0usize;
        let mut dist2_min = sols[0].sq_dist;
        for (j, sol) in sols.iter().enumerate().skip(1) {
            if sol.sq_dist < dist2_min {
                dist2_min = sol.sq_dist;
                jmin = j;
            }
        }
        param = sols[jmin].u;
    } else {
        // `cxx:2243-2257`: `TrimmedSquareDistances`.
        let first = c.first_parameter();
        let last = c.last_parameter();
        let dist1_2 = p.square_distance(&c.d0(first));
        let dist2_2 = p.square_distance(&c.d0(last));
        if dist1_2 < dist2_2 {
            *done = -1;
            param = first;
        } else {
            *done = 1;
            param = last;
        }
    }

    // `cxx:2259-2266`
    if cbase.is_periodic() {
        let per = cbase.period();
        let tolp = Precision::parametric(CONFUSION, 1.0e2);
        if (per - param).abs() <= tolp {
            param = 0.0;
        }
    }
    param
}

/// `static double SmartParameter(Einf, EdgeTol, Pnt, sign, S1, S2)`
/// (`_1.cxx:2274-2379`).
///
/// `sign` is `done` as `parameter` reported it: `-1` prepends the vertex point
/// to the pcurve, every other value appends it (`cxx:2333-2351`).
pub(super) fn smart_parameter(
    einf: &mut DraftEdgeInfo,
    edge_tol: f64,
    pnt: &GpPnt,
    sign: i32,
    s1: &Arc<dyn Surface>,
    s2: &Arc<dyn Surface>,
) -> f64 {
    // `cxx:2281-2283`
    let tol = CONFUSION;
    let mut etol = edge_tol;

    // `cxx:2285-2307`: fill the pcurves that are still missing.
    if einf.first_pc().is_none() {
        let mut projected = None;
        if let Some(curve) = einf.geometry().cloned() {
            let (f, l) = (curve.first_parameter(), curve.last_parameter());
            projected = geom_proj_lib::curve_2d_range(&curve, f, l, s1, &mut etol);
        }
        *einf.change_first_pc() = projected;
    }
    if einf.second_pc().is_none() {
        let mut projected = None;
        if let Some(curve) = einf.geometry().cloned() {
            let (f, l) = (curve.first_parameter(), curve.last_parameter());
            projected = geom_proj_lib::curve_2d_range(&curve, f, l, s2, &mut etol);
        }
        *einf.change_second_pc() = projected;
    }

    // `cxx:2309-2311`: `GeomAPI_ProjectPointOnSurf Projector(Pnt, S1)`.
    let Some(on_s1) = project_point_on_surface(s1.as_ref(), pnt, tol) else {
        panic!("Standard_Failure: Draft_Modification_1::SmartParameter: LowerDistanceParameters.");
    };
    let p2d = GpPnt2d::new(on_s1.u, on_s1.v);

    // `cxx:2313-2318`: `NewC2d = Einf.FirstPC()`, trimmed curves unwrapped.
    let mut new_c2d = einf.first_pc().cloned();
    if let Some(basis) = new_c2d.as_ref().and_then(|c| c.trimmed_basis()) {
        new_c2d = Some(Arc::from(basis.clone_dyn()));
    }

    // `cxx:2320-2323`
    let projector = new_c2d
        .as_ref()
        .and_then(|c| project_point_on_curve_2d(c.as_ref(), &p2d, tol));
    let needs_patch = match &projector {
        None => true,
        Some(p) => p.distance > tol,
    };

    if needs_patch {
        // `cxx:2325-2353`
        let Some(mut bcurve) = new_c2d
            .as_ref()
            .and_then(|c| c.bspline_copy2d())
            .or_else(|| {
                new_c2d
                    .as_ref()
                    .and_then(|c| curve_to_bspline_curve_bspl(c.as_ref(), ParameterisationType::QuasiAngular))
            })
        else {
            panic!("Standard_Failure: Draft_Modification_1::SmartParameter: CurveToBSplineCurve.");
        };
        let last_pole = bcurve.nb_poles().saturating_sub(1);
        let patch_poles = if sign == -1 {
            // `cxx:2334-2339`: `PntArray(1) = P2d; PntArray(2) = BCurve->Pole(1)`.
            let pole1 = bcurve.pole(0).unwrap_or(p2d);
            vec![p2d, pole1]
        } else {
            // `cxx:2344-2349`: `PntArray(1) = BCurve->Pole(BCurve->NbPoles())`.
            let pole_last = bcurve.pole(last_pole).unwrap_or(p2d);
            vec![pole_last, p2d]
        };
        let Ok(patch) = Geom2dBezierCurve::new(patch_poles) else {
            panic!("Standard_Failure: Draft_Modification_1::SmartParameter: Geom2d_BezierCurve.");
        };
        let mut concat =
            CompCurveToBSplineCurve2d::from_bspline(bcurve.clone(), ParameterisationType::QuasiAngular);
        let patch_curve: Arc<dyn Curve2d> = Arc::new(patch);
        let after = sign != -1;
        let _ = concat.add(patch_curve.as_ref(), tol, after);
        if let Some(bs) = concat.into_curve() {
            bcurve = bs;
        }
        new_c2d = Some(Arc::new(bcurve));
    }

    // `cxx:2355`
    *einf.change_first_pc() = new_c2d.clone();

    // `cxx:2356-2365`: `Adaptor3d_CurveOnSurface` over `(NewC2d, S1)`, then
    // `ProjLib_HCompProjectedCurve(Hsur2, hcons, Tol, Tol)`.
    let Some(pcurve) = new_c2d else {
        panic!("Standard_Failure: Draft_Modification_1::SmartParameter: null pcurve.");
    };
    let cons: Arc<dyn Curve> = Arc::new(CurveOnSurface::new(pcurve, s1.clone()));
    let mut hproj = CompProjectedCurve::with_tolerances(s2.clone(), cons, tol, tol);

    // `cxx:2363-2366`
    let Some((udeb, ufin)) = hproj.bounds(1) else {
        panic!("Standard_Failure: Draft_Modification_1::SmartParameter: no projected part.");
    };
    let max_seg = 20 + hproj.nb_intervals(5).unwrap_or(0);
    let projector2d: Arc<dyn Curve2d> = Arc::new(ProjectedCurve2d::new(hproj));
    let mut appr = ApproxCurveOnSurface::new(projector2d, s2.clone(), udeb, ufin, tol);
    // `cxx:2366`: `Perform(MaxSeg, 10, GeomAbs_C1, false, false)`.
    let _ = appr.perform(max_seg, 10, 2, false, false);

    // `cxx:2367-2369`
    *einf.change_second_pc() = appr
        .curve2d()
        .cloned()
        .map(|c| Arc::new(c) as Arc<dyn Curve2d>);
    *einf.change_geometry() = appr.curve3d().cloned().map(|c| Arc::new(c) as Arc<dyn Curve>);
    einf.set_new_geometry(true);

    // `cxx:2371-2378`
    match einf.geometry() {
        Some(g) => {
            if sign == -1 {
                g.first_parameter()
            } else {
                g.last_parameter()
            }
        }
        None => {
            panic!("Standard_Failure: Draft_Modification_1::SmartParameter: no 3d curve.");
        }
    }
}

/// `slib::plane_parameters` / `slib::cylinder_parameters` behind the
/// `ElSLib::Parameters` switch of `Perform` (`_1.cxx:1050-1067`).
pub(super) fn surface_parameters_plane(pl: &occt_core::gp::GpPln, p: &GpPnt) -> (f64, f64) {
    slib::plane_parameters(&pl.position(), p)
}

/// `ElSLib::Parameters(Cylinder, P)` (`_1.cxx:1060-1065`).
pub(super) fn surface_parameters_cylinder(cy: &occt_core::gp::GpCylinder, p: &GpPnt) -> (f64, f64) {
    slib::cylinder_parameters(&cy.position(), p)
}
