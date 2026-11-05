//! `BRepFill_CurveConstraint` — `GeomPlate_CurveConstraint` that also accepts a
//! `BRepAdaptor_Surface`.
//!
//! Source: `TKBool/BRepFill/BRepFill_CurveConstraint.cxx:32-84` and `.hxx:29-55`
//! (OCCT 8.0.0).
//!
//! The reference class adds no state to its base: both constructors fill
//! `myFrontiere` / `myTolDist` / `myTolAng` / `myTolCurv` / `myLProp` / `myOrder`
//! / `myNbPoints` / `myConstG0..2` and nothing else, so the port is two
//! constructor functions that answer the base
//! [`GeomPlateCurveConstraint`](crate::..::geom_plate::CurveConstraint) — there
//! is no separate type to carry around, and `BuildPlateSurface::Add` takes the
//! base type anyway.
//!
//! Differences from the base constructor that are reproduced here:
//!
//! * the continuity gate throws `"BRepFill : The continuity is not G0 G1 or
//!   G2"` (`cxx:48-51`) instead of the `GeomPlate` wording, and it runs
//!   *before* the frontier test;
//! * the 3D-curve constructor delegates to
//!   `GeomPlate_CurveConstraint(Boundary, Tang, NPt, TolDist)` (`cxx:79-84`),
//!   i.e. the header defaults `TolAng = 0.01`, `TolCurv = 0.1`
//!   (`GeomPlate_CurveConstraint.hxx:39-45`).
//!
//! The surface arm of the reference (`cxx:59-73`) distinguishes a
//! `GeomAdaptor_Surface` from a `BRepAdaptor_Surface` and only the latter goes
//! through `BRep_Tool::Surface(BS1->Face())`. The port's
//! [`CurveOnSurface`] always carries the `Geom_Surface` itself, so both arms
//! end in the same `myLProp.SetSurface(Surf)` and the port keeps the single
//! `BoundaryCurve::OnSurface` path.

use std::sync::Arc;

use occt_geom::curve::Curve;
use occt_geom::curve_on_surface::CurveOnSurface;
use occt_geom::geom_plate::{BoundaryCurve, CurveConstraint};

/// `BRepFill_CurveConstraint(Boundary, Order, NPt, TolDist, TolAng, TolCurv)`
/// (`cxx:32-73`).
///
/// Raises when `Order < -1 || Order > 2` (`cxx:48-51`) and — in the reference —
/// when `Boundary` is null (`cxx:56-59`); an `Arc<CurveOnSurface>` cannot be
/// null, so that arm cannot fire here.
pub fn from_curve_on_surface(
    boundary: Arc<CurveOnSurface>,
    order: i32,
    nb_points: i32,
    tol_dist: f64,
    tol_ang: f64,
    tol_curv: f64,
) -> Result<CurveConstraint, String> {
    if !(-1..=2).contains(&order) {
        return Err("BRepFill : The continuity is not G0 G1 or G2".to_string());
    }
    CurveConstraint::from_boundary(
        BoundaryCurve::OnSurface(boundary),
        order,
        nb_points,
        tol_dist,
        tol_ang,
        tol_curv,
    )
}

/// `BRepFill_CurveConstraint(Boundary, Tang, NPt, TolDist)` (`cxx:79-84`): the
/// free 3D-curve arm, which forwards to
/// `GeomPlate_CurveConstraint(Boundary, Tang, NPt, TolDist)`.
pub fn from_curve(
    boundary: Arc<dyn Curve>,
    order: i32,
    nb_points: i32,
    tol_dist: f64,
) -> Result<CurveConstraint, String> {
    // `GeomPlate_CurveConstraint.hxx:39-45`: TolAng = 0.01, TolCurv = 0.1.
    CurveConstraint::from_boundary(
        BoundaryCurve::Curve3d(boundary),
        order,
        nb_points,
        tol_dist,
        0.01,
        0.1,
    )
}
