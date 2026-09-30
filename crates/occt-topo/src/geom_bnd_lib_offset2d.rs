//! Bounding box of a 2D offset curve.
//!
//! Source: `GeomBndLib_OffsetCurve2d.cxx`. Fast paths:
//! * basis line → translate location by `offset * Normal(Direction)` then
//!   `GeomBndLib_Line2d::Box`;
//! * basis circle → `NewRadius = R + offset * sign(Normal·Radial)` then
//!   `GeomBndLib_Circle2d::Box` (or the centre if the radius vanishes).
//! Generic fallback: `GeomBndLib_OtherCurve2d::Box` on the offset curve
//! itself (not `basis_box.Enlarge(|offset|)`).
//! T-97: items below are faithful ports of the named OCCT source, but their
//! OCCT-side consumers are not all ported yet, so parts are not called from this
//! crate. The `dead_code` allowance is deliberate: **pending wiring**, not dead
//! code. Do not delete them to silence warnings (see
//! specs/_a3n00_gap_analysis.md §9.309/§9.310); wire the consumer instead.
#![allow(dead_code)]

use occt_core::bnd::BndBox2d;
use occt_core::gp::{GpCirc2d, GpLin2d, GpPnt2d, GpVec2d};
use occt_core::precision::{CONFUSION, SQUARE_CONFUSION};
use occt_geom2d::curve::Curve2d;
use occt_geom2d::offset::Geom2dOffsetCurve;

use crate::geom_bnd_lib_circle2d::box_circ_range;
use crate::geom_bnd_lib_line2d::box_lin;
use crate::geom_bnd_lib_other2d::box_other;
use crate::pcurve::{pc_curve_kind, CurveKind};

fn try_offset_line2d(
    basis: &dyn Curve2d,
    offset: f64,
    the_u1: f64,
    the_u2: f64,
    the_tol: f64,
) -> Option<BndBox2d> {
    if pc_curve_kind(basis) != CurveKind::Line {
        return None;
    }
    let loc = basis.d0(0.0);
    let (_, tan) = basis.d1(0.0);
    let dir = occt_core::gp::GpDir2d::from_vec2d(&tan).ok()?;
    let mut a_normal = GpVec2d::new(dir.y(), -dir.x());
    if a_normal.square_magnitude() <= SQUARE_CONFUSION {
        return None;
    }
    let _ = a_normal.normalize();
    a_normal.multiply_scalar(offset);
    let a_loc = GpPnt2d::new(loc.x() + a_normal.x(), loc.y() + a_normal.y());
    let a_lin = GpLin2d::from_pnt_dir(a_loc, dir);
    let mut a_local = box_lin(&a_lin, the_u1, the_u2, 0.0).ok()?;
    a_local.enlarge(the_tol);
    Some(a_local)
}

fn try_offset_circle2d(
    basis: &dyn Curve2d,
    offset: f64,
    the_u1: f64,
    the_u2: f64,
    the_tol: f64,
) -> Option<BndBox2d> {
    if pc_curve_kind(basis) != CurveKind::Circle {
        return None;
    }
    let a_u_mid = 0.5 * (the_u1 + the_u2);
    let (a_p, a_v1) = basis.d1(a_u_mid);
    if a_v1.square_magnitude() <= SQUARE_CONFUSION {
        return None;
    }
    let mut a_normal = GpVec2d::new(a_v1.y(), -a_v1.x());
    if a_normal.square_magnitude() <= SQUARE_CONFUSION {
        return None;
    }
    let _ = a_normal.normalize();
    // Recover centre as the unique point at equal distance from the circle
    // samples: for a Geom2dCircle, d0(0) and d0(pi) are opposite, midpoint
    // is the centre. That matches Circ2d::Location for the standard parametrisation.
    let p0 = basis.d0(0.0);
    let ppi = basis.d0(std::f64::consts::PI);
    let centre = GpPnt2d::new(0.5 * (p0.x() + ppi.x()), 0.5 * (p0.y() + ppi.y()));
    let a_radius = centre.distance(&p0);
    let mut a_radial = GpVec2d::new(a_p.x() - centre.x(), a_p.y() - centre.y());
    if a_radial.square_magnitude() <= SQUARE_CONFUSION {
        return None;
    }
    let _ = a_radial.normalize();
    let a_sign = a_normal.dot(&a_radial);
    if (a_sign.abs() - 1.0).abs() > 1e-9 {
        return None;
    }
    let a_new_radius = a_radius + offset * a_sign;
    let mut a_local = BndBox2d::new();
    if a_new_radius > CONFUSION {
        let mut a_offset_circ = GpCirc2d::default();
        a_offset_circ.set_location(centre);
        a_offset_circ.set_radius(a_new_radius);
        a_local = box_circ_range(&a_offset_circ, the_u1, the_u2, 0.0);
    } else if a_new_radius.abs() <= CONFUSION {
        a_local.add_point(&centre);
    } else {
        return None;
    }
    a_local.enlarge(the_tol);
    Some(a_local)
}

/// `GeomBndLib_OffsetCurve2d::Box(theU1, theU2, theTol)`.
pub fn box_offset(off: &Geom2dOffsetCurve, the_u1: f64, the_u2: f64, the_tol: f64) -> BndBox2d {
    let basis = off.basis_curve().as_ref();
    let offset = off.offset_value();
    if let Some(b) = try_offset_line2d(basis, offset, the_u1, the_u2, the_tol) {
        if !b.is_void() {
            return b;
        }
    }
    if let Some(b) = try_offset_circle2d(basis, offset, the_u1, the_u2, the_tol) {
        if !b.is_void() {
            return b;
        }
    }
    box_other(off, the_u1, the_u2, the_tol)
}

/// Offset box from the `Curve2d` trait (no downcast): OtherCurve2d Box.
pub fn box_offset_as_curve(curve: &dyn Curve2d, the_u1: f64, the_u2: f64, the_tol: f64) -> BndBox2d {
    box_other(curve, the_u1, the_u2, the_tol)
}
