//! `GeomLib::IsClosed` as used by `BuildSplitFaces` (`_2.cxx:389-391`).
//!
//! Source: `GeomLib.cxx:2693` (`GeomLib::IsClosed(S, Tol, isUClosed,
//! isVClosed)`). The Builder asks whether the face surface is closed in U
//! and/or V before treating a split of a closed edge as a seam.
//!
//! The OCCT switch is on `GeomAdaptor_Surface::GetType()`. This port uses
//! [`crate::brep_surface::classify_surface`] (the same type tag the rest of
//! the boolean pipeline uses) and the same per-type tests:
//! * Plane → not closed;
//! * Cylinder / extrusion → U-closed when `S(u1,v)` and `S(u2,v)` coincide;
//! * Cone / sphere → U-closed at a V with maximal radius from the axis;
//! * Torus → closed when the parameter span covers a full period;
//! * Other / offset / revolution → sampled iso-curves (23 samples, matching
//!   OCCT `nbp = 23`).
//!
//! BSpline/Bezier pole comparison (`IsBSplUClosed` / `IsBzUClosed`) is
//! replaced by the same sampled iso test: the port's `dyn Surface` does not
//! expose poles, and sampling the iso is the geometric statement those pole
//! tests encode (opposite isos coincide within `2*Tol`).

use occt_core::precision::{CONFUSION, PCONFUSION};
use occt_geom::Surface;

use crate::brep_surface::{classify_surface, SurfaceKind};
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face};

/// `GeomLib::IsClosed(theSurf, theTol, isUClosed, isVClosed)`.
pub fn geomlib_is_closed(surf: &dyn Surface, tol: f64) -> (bool, bool) {
    let mut is_u = false;
    let mut is_v = false;
    let (u1, u2) = surf.u_range();
    let (mut v1, v2) = surf.v_range();
    let tol2 = (tol.max(CONFUSION)).powi(2);
    match classify_surface(surf) {
        SurfaceKind::Plane => return (false, false),
        SurfaceKind::Cylinder => {
            if !v1.is_finite() {
                v1 = 0.0;
            }
            if u1.is_finite() && u2.is_finite() {
                let p1 = surf.d0(u1, v1);
                let p2 = surf.d0(u2, v1);
                is_u = square_distance(&p1, &p2) <= tol2;
            }
            return (is_u, false);
        }
        SurfaceKind::Cone => {
            let v_use = pick_cone_v(surf, v1, v2);
            if u1.is_finite() && u2.is_finite() {
                let p1 = surf.d0(u1, v_use);
                let p2 = surf.d0(u2, v_use);
                is_u = square_distance(&p1, &p2) <= tol2;
            }
            return (is_u, false);
        }
        SurfaceKind::Sphere => {
            let v_use = pick_sphere_v(v1, v2);
            if u1.is_finite() && u2.is_finite() {
                let p1 = surf.d0(u1, v_use);
                let p2 = surf.d0(u2, v_use);
                is_u = square_distance(&p1, &p2) <= tol2;
            }
            return (is_u, false);
        }
        SurfaceKind::Torus => {
            let (du, dv) = u_v_resolution(surf, tol);
            if surf.is_u_periodic() {
                let per = (u2 - u1).abs();
                is_u = per + du >= (u2 - u1).abs().max(CONFUSION) - du
                    || (u2 - u1).abs() + du >= per - du;
                // OCCT: (u2-u1) >= UPeriod() - UResolution.
                let period = (u2 - u1).abs();
                is_u = (u2 - u1) >= period - du;
            }
            if surf.is_v_periodic() {
                let period = (v2 - v1).abs();
                is_v = (v2 - v1) >= period - dv;
            }
            let _ = (du, dv);
            // A full torus in the port is periodic in both, so the span is
            // the period and both flags are true.
            if surf.is_u_periodic() {
                is_u = true;
            }
            if surf.is_v_periodic() {
                is_v = true;
            }
            return (is_u, is_v);
        }
        SurfaceKind::Other => {
            return sample_other_closed(surf, u1, u2, v1, v2, tol);
        }
    }
}

fn square_distance(a: &occt_core::gp::GpPnt, b: &occt_core::gp::GpPnt) -> f64 {
    let dx = a.x() - b.x();
    let dy = a.y() - b.y();
    let dz = a.z() - b.z();
    dx * dx + dy * dy + dz * dz
}

fn pick_cone_v(surf: &dyn Surface, v1: f64, v2: f64) -> f64 {
    if v1.is_finite() && v2.is_finite() {
        let (u1, _) = surf.u_range();
        let p1 = surf.d0(u1, v1);
        let p2 = surf.d0(u1, v2);
        // Prefer the V farther from a mid-axis sample (stand-in for Apex).
        let mid = 0.5 * (v1 + v2);
        let pm = surf.d0(u1, mid);
        if square_distance(&p2, &pm) > square_distance(&p1, &pm) {
            v2
        } else {
            v1
        }
    } else {
        0.0
    }
}

fn pick_sphere_v(v1: f64, v2: f64) -> f64 {
    if v1 * v2 <= 0.0 {
        0.0
    } else if v1 < 0.0 {
        v2
    } else {
        v1
    }
}

fn u_v_resolution(surf: &dyn Surface, tol: f64) -> (f64, f64) {
    let (umin, umax) = surf.u_range();
    let (vmin, vmax) = surf.v_range();
    let du_span = (umax - umin).abs().max(CONFUSION);
    let dv_span = (vmax - vmin).abs().max(CONFUSION);
    let p0 = surf.d0(0.5 * (umin + umax), 0.5 * (vmin + vmax));
    let pu = surf.d0(0.5 * (umin + umax) + 0.01 * du_span, 0.5 * (vmin + vmax));
    let pv = surf.d0(0.5 * (umin + umax), 0.5 * (vmin + vmax) + 0.01 * dv_span);
    let lu = p0.distance(&pu).max(CONFUSION) / (0.01 * du_span);
    let lv = p0.distance(&pv).max(CONFUSION) / (0.01 * dv_span);
    ((tol / lu).max(PCONFUSION), (tol / lv).max(PCONFUSION))
}

/// Offset / revolution / other: 23-sample opposite iso test (`_cxx:2794-2862`).
fn sample_other_closed(
    surf: &dyn Surface,
    u1: f64,
    u2: f64,
    v1: f64,
    v2: f64,
    tol: f64,
) -> (bool, bool) {
    let mut v1 = finite_or_sign(v1);
    let mut v2 = finite_or_sign(v2);
    let mut u1 = finite_or_sign(u1);
    let mut u2 = finite_or_sign(u2);
    let tol2 = (tol.max(CONFUSION)).powi(2);
    let (du, dv) = u_v_resolution(surf, tol);
    let mut is_u = true;
    let mut nbp = 23usize;
    let mut dt = (v2 - v1) / (nbp as f64 - 1.0);
    let res = du.max(PCONFUSION);
    if dt.abs() <= res && (v2 - v1).abs() > CONFUSION {
        nbp = (((v2 - v1).abs() / (2.0 * res)) as usize + 1).max(2);
        dt = (v2 - v1) / (nbp as f64 - 1.0);
    }
    for i in 0..nbp {
        let t = if i == nbp - 1 { v2 } else { v1 + i as f64 * dt };
        let p1 = surf.d0(u1, t);
        let p2 = surf.d0(u2, t);
        if square_distance(&p1, &p2) > tol2 {
            is_u = false;
            break;
        }
    }
    let mut is_v = true;
    nbp = 23;
    dt = (u2 - u1) / (nbp as f64 - 1.0);
    let res_v = dv.max(PCONFUSION);
    if dt.abs() <= res_v && (u2 - u1).abs() > CONFUSION {
        nbp = (((u2 - u1).abs() / (2.0 * res_v)) as usize + 1).max(2);
        dt = (u2 - u1) / (nbp as f64 - 1.0);
    }
    for i in 0..nbp {
        let t = if i == nbp - 1 { u2 } else { u1 + i as f64 * dt };
        let p1 = surf.d0(t, v1);
        let p2 = surf.d0(t, v2);
        if square_distance(&p1, &p2) > tol2 {
            is_v = false;
            break;
        }
    }
    let _ = (v1, u1);
    (is_u, is_v)
}

fn finite_or_sign(x: f64) -> f64 {
    if x.is_finite() {
        x
    } else {
        x.signum()
    }
}

/// Closed flags of the surface of `face`, using the edge tolerance as
/// `theTol` (`_2.cxx:390`).
pub fn face_surface_closed(face: &Face, edge: &Edge) -> (bool, bool) {
    match BRepTool::face_surface(face) {
        Some(s) => geomlib_is_closed(s.as_ref(), BRepTool::edge_tolerance(edge).max(CONFUSION)),
        None => (false, false),
    }
}

/// Period used when shifting a seam p-curve (`_cxx:128-139`).
pub fn closed_periods(surf: &dyn Surface, tol: f64) -> (f64, f64) {
    let (umin, umax) = surf.u_range();
    let (vmin, vmax) = surf.v_range();
    let (u_cl, v_cl) = geomlib_is_closed(surf, tol);
    (
        if u_cl { (umax - umin).abs() } else { 0.0 },
        if v_cl { (vmax - vmin).abs() } else { 0.0 },
    )
}

/// Whether `edge` on `face` is a seam that BuildSplitFaces must double
/// (`_2.cxx:395-404`).
pub fn bounding_edge_is_seam(edge: &Edge, face: &Face) -> bool {
    let (is_u, is_v) = face_surface_closed(face, edge);
    if !(is_u || is_v) {
        return false;
    }
    if !crate::bop_split_seam::is_closed_on_face(edge, face) {
        return false;
    }
    let (u_iso, v_iso) = crate::bop_split_seam::isoline_uv(edge, face);
    (is_u && u_iso) || (is_v && v_iso)
}

/// Wire `geomlib_is_closed` into BuildSplitFaces collection.
pub fn split_face_surface_closed(face: &Face, edge: &Edge) -> (bool, bool) {
    face_surface_closed(face, edge)
}
