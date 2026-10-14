//! `IntTools_Tools::SegPln` and `IntTools_Tools::ComputeTolerance`.
//!
//! Source: `IntTools_Tools.cxx:670` (SegPln) and `737` (ComputeTolerance).
//!
//! `SegPln` classifies a line segment against a plane: common block (both
//! ends inside the summed tolerance), one-sided miss, out of range, or a
//! crossing with the intersection parameter and a tolerance interval.
//!
//! `ComputeTolerance` in OCCT uses `GeomLib_CheckCurveOnSurface`, which
//! maximises the distance between the 3D curve and the surface image of the
//! pcurve with `math_PSO` plus `math_NewtonMinimum` and inflates the result by
//! the `anEps` margin (`IntTools_Tools.cxx:767-775`).
//! T-97: items below are faithful ports of the named OCCT source, but their
//! OCCT-side consumers are not all ported yet, so parts are not called from this
//! crate. The `dead_code` allowance is deliberate: **pending wiring**, not dead
//! code. Do not delete them to silence warnings (see
//! specs/_a3n00_gap_analysis.md §9.309/§9.310); wire the consumer instead.
#![allow(dead_code)]

use std::sync::Arc;

use occt_core::gp::{GpLin, GpPln, GpPnt};
use occt_core::precision::PCONFUSION;
use occt_geom::{CheckCurveOnSurface, Curve, Surface};
use occt_geom2d::curve::Curve2d;

use crate::meshing::edge_discret::CurveOnSurface;

/// Result of [`seg_pln`].
#[derive(Debug, Clone, Copy)]
pub struct SegPlnResult {
    /// OCCT `iRet`: 0 intersection, 1 common block, 2 one side, 3 out of range.
    pub status: i32,
    /// Intersection point (status 0).
    pub point: GpPnt,
    /// Parameter of the intersection on the line (status 0).
    pub t: f64,
    /// `theTolP` = `tol_lin + tol_pln` (status 0).
    pub tol_p: f64,
    /// `theTPmin` = `t - tol_pln` (status 0).
    pub t_min: f64,
    /// `theTPmax` = `t + tol_pln` (status 0).
    pub t_max: f64,
}

/// `IntTools_Tools::SegPln` (`IntTools_Tools.cxx:670`).
pub fn seg_pln(
    lin: &GpLin,
    t_lin1: f64,
    t_lin2: f64,
    tol_lin: f64,
    pln: &GpPln,
    tol_pln: f64,
) -> SegPlnResult {
    let a_tol = tol_lin + tol_pln;
    let dir_pln = pln.position().direction();
    let loc_pln = pln.position().location();
    let dir_lin = lin.direction();
    let loc_lin = lin.location();

    let p1 = GpPnt::new(
        loc_lin.x() + t_lin1 * dir_lin.x(),
        loc_lin.y() + t_lin1 * dir_lin.y(),
        loc_lin.z() + t_lin1 * dir_lin.z(),
    );
    let dist1 = dir_pln.x() * (p1.x() - loc_pln.x())
        + dir_pln.y() * (p1.y() - loc_pln.y())
        + dir_pln.z() * (p1.z() - loc_pln.z());
    let p2 = GpPnt::new(
        loc_lin.x() + t_lin2 * dir_lin.x(),
        loc_lin.y() + t_lin2 * dir_lin.y(),
        loc_lin.z() + t_lin2 * dir_lin.z(),
    );
    let dist2 = dir_pln.x() * (p2.x() - loc_pln.x())
        + dir_pln.y() * (p2.y() - loc_pln.y())
        + dir_pln.z() * (p2.z() - loc_pln.z());

    if dist1.abs() < a_tol && dist2.abs() < a_tol {
        return SegPlnResult {
            status: 1,
            point: p1,
            t: t_lin1,
            tol_p: a_tol,
            t_min: t_lin1,
            t_max: t_lin2,
        };
    }
    if dist1 * dist2 > 0.0 {
        return SegPlnResult {
            status: 2,
            point: p1,
            t: t_lin1,
            tol_p: a_tol,
            t_min: t_lin1,
            t_max: t_lin2,
        };
    }

    let (a, b, c, d) = {
        let n = dir_pln;
        let loc = loc_pln;
        (n.x(), n.y(), n.z(), -(n.x() * loc.x() + n.y() * loc.y() + n.z() * loc.z()))
    };
    let e = a * loc_lin.x() + b * loc_lin.y() + c * loc_lin.z() + d;
    let h = a * dir_lin.x() + b * dir_lin.y() + c * dir_lin.z();
    if h.abs() <= 1e-30 {
        return SegPlnResult {
            status: 2,
            point: p1,
            t: t_lin1,
            tol_p: a_tol,
            t_min: t_lin1,
            t_max: t_lin2,
        };
    }
    let tp = -e / h;
    if tp < t_lin1 - a_tol || tp > t_lin2 + a_tol {
        return SegPlnResult {
            status: 3,
            point: p1,
            t: tp,
            tol_p: a_tol,
            t_min: tp,
            t_max: tp,
        };
    }
    let p = GpPnt::new(
        loc_lin.x() + tp * dir_lin.x(),
        loc_lin.y() + tp * dir_lin.y(),
        loc_lin.z() + tp * dir_lin.z(),
    );
    SegPlnResult {
        status: 0,
        point: p,
        t: tp,
        tol_p: a_tol,
        t_min: tp - tol_pln,
        t_max: tp + tol_pln,
    }
}

/// `IntTools_Tools::ComputeTolerance` (`IntTools_Tools.cxx:737-779`).
///
/// Returns `(theMaxDist, theMaxPar)`; `None` is OCCT's `false` return, i.e. the
/// per-curve search failed. `theCurve3D` / `theCurve2D` / `theSurf` are the
/// `Handle`s of `cxx:737-739` (the port passes them by `Arc` so that no
/// geometry is copied), and `theFirst` / `theLast` are the
/// `GeomAdaptor_Curve(theCurve3D, theFirst, theLast)` /
/// `Geom2dAdaptor_Curve(theCurve2D, theFirst, theLast)` load range
/// (`cxx:749-753`) — every OCCT caller passes the curve's own
/// `FirstParameter()` / `LastParameter()` (`IntTools_FaceFace.cxx:636-637`).
///
/// UNPORTED: `theToRunParallel` (`cxx:745`, `cxx:760`). Both OCCT callers keep
/// the `false` default, and the port always runs the per-interval search
/// serially.
pub fn compute_tolerance(
    curve3d: &Arc<dyn Curve>,
    curve2d: &Arc<dyn Curve2d>,
    surf: &Arc<dyn Surface>,
    first: f64,
    last: f64,
) -> Option<(f64, f64)> {
    // A degenerate range makes `math_PSO` bail out (`cxx:654-657`) and OCCT
    // would then report `sqrt(RealLast())` as the max distance, which the
    // callers compare against `aTolC` as a real value. The port keeps its
    // pre-existing guard and reports "not computed" instead.
    if last - first < PCONFUSION {
        return None;
    }

    // `Adaptor3d_CurveOnSurface(aGeom2dAdaptorCurve, aGeomAdaptorSurface)`
    // (`cxx:756-757`).
    let an_adaptor3d_curve_on_surface =
        CurveOnSurface::new(curve2d.clone(), surf.clone(), first, last);

    // `aCS.Init(aGeomAdaptorCurve, theTolRange)` (`cxx:759`).
    let mut a_cs = CheckCurveOnSurface::with_curve(curve3d.as_ref(), PCONFUSION);
    // `aCS.SetParallel(theToRunParallel)` (`cxx:760`) is left at its default
    // `false` (`cxx:293`).
    a_cs.perform(&an_adaptor3d_curve_on_surface, curve2d.as_ref());

    if !a_cs.is_done() {
        return None;
    }

    // Obtaining a precise result is impossible with numeric methods, so a
    // margin is applied; otherwise a later, more precise minimum point would
    // invalidate the tolerance computed here (`cxx:767-775`).
    let an_eps = 1.0 + 1.0e-5;
    Some((an_eps * a_cs.max_distance(), a_cs.max_parameter()))
}

/// Correct `IsOnPave1` (`IntTools_Tools.cxx:627`): inside the range, or
/// within `aTolerance` of either bound.
pub fn is_on_pave1_occt(a_tr: f64, first: f64, last: f64, a_tolerance: f64) -> bool {
    if a_tr >= first && a_tr <= last {
        return true;
    }
    (a_tr - first).abs() <= a_tolerance || (a_tr - last).abs() <= a_tolerance
}

/// Correct `IsInRange` (`IntTools_Tools.cxx:650`).
pub fn is_in_range_occt(
    ref_first: f64,
    ref_last: f64,
    first: f64,
    last: f64,
    a_tolerance: f64,
) -> bool {
    let a_t_ref1 = ref_first - a_tolerance;
    let a_t_ref2 = ref_last + a_tolerance;
    (first >= a_t_ref1 && first <= a_t_ref2) || (last >= a_t_ref1 && last <= a_t_ref2)
}
