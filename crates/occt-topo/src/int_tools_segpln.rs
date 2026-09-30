//! `IntTools_Tools::SegPln` and `IntTools_Tools::ComputeTolerance`.
//!
//! Source: `IntTools_Tools.cxx:670` (SegPln) and `737` (ComputeTolerance).
//!
//! `SegPln` classifies a line segment against a plane: common block (both
//! ends inside the summed tolerance), one-sided miss, out of range, or a
//! crossing with the intersection parameter and a tolerance interval.
//!
//! `ComputeTolerance` in OCCT uses `GeomLib_CheckCurveOnSurface`. The port
//! samples the 3D curve against the surface image of the pcurve and returns
//! the max distance inflated by `1 + 1e-5`, matching the OCCT `anEps` margin.
//! T-97: items below are faithful ports of the named OCCT source, but their
//! OCCT-side consumers are not all ported yet, so parts are not called from this
//! crate. The `dead_code` allowance is deliberate: **pending wiring**, not dead
//! code. Do not delete them to silence warnings (see
//! specs/_a3n00_gap_analysis.md §9.309/§9.310); wire the consumer instead.
#![allow(dead_code)]

use occt_core::gp::{GpLin, GpPln, GpPnt};
use occt_geom::{Curve, Surface};
use occt_geom2d::curve::Curve2d;

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

/// Sampled `IntTools_Tools::ComputeTolerance` (`IntTools_Tools.cxx:737`).
///
/// Returns `(max_dist * (1 + 1e-5), parameter_of_max)` when both curves and
/// the surface evaluate; `None` when the range is empty.
pub fn compute_tolerance(
    curve3d: &dyn Curve,
    curve2d: &dyn Curve2d,
    surf: &dyn Surface,
    first: f64,
    last: f64,
) -> Option<(f64, f64)> {
    if last - first < occt_core::precision::PCONFUSION {
        return None;
    }
    const N: usize = 64;
    let mut max_dist = 0.0f64;
    let mut max_par = first;
    for i in 0..=N {
        let t = first + (last - first) * (i as f64 / N as f64);
        let p3 = curve3d.d0(t);
        let p2 = curve2d.d0(t);
        let ps = surf.d0(p2.x(), p2.y());
        let d = p3.distance(&ps);
        if d > max_dist {
            max_dist = d;
            max_par = t;
        }
    }
    let an_eps = 1.0 + 1.0e-5;
    Some((an_eps * max_dist, max_par))
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
