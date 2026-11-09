//! `Draft_Modification::NewCurve` - the private overload
//! (`Draft_Modification_1.cxx:1963-2044`).
//!
//! # Translation notes
//!
//! * `gp_Dir Norm;` (`cxx:2005`) is a default-constructed direction, i.e.
//!   `(0, 0, 1)`, and keeps that value when the support is neither a cylinder
//!   nor a cone; the reference then feeds it to `Norm.Crossed(Direction)`
//!   unchanged, so the port does the same.
//! * The commented-out `testdir` guard (`cxx:2001-2004`) is dead code in the
//!   reference and is not resurrected here.

use std::sync::Arc;

use occt_core::elib::{slib, surface_eval};
use occt_core::gp::{GpAx1, GpDir, GpLin, GpPln};
use occt_core::precision::ANGULAR;
use occt_geom::curve::Curve;
use occt_geom::intana::IntAnaIntConicQuad;
use occt_geom::line::GeomLine;
use occt_geom::Surface;

use crate::abs::Orientation;

use super::new_surface::find_rotation;

/// `Draft_Modification::NewCurve(C, S, Oris, Direction, Angle, NeutralPlane,
/// Flag)` (`Draft_Modification_1.cxx:1963-2044`): the drafted image of `c`, or
/// `None` (the reference's null handle) when the construction is not defined.
pub fn new_curve(
    c: &Arc<dyn Curve>,
    s: &Arc<dyn Surface>,
    oris: Orientation,
    direction: &GpDir,
    angle: f64,
    neutral_plane: &GpPln,
) -> Option<Arc<dyn Curve>> {
    // `cxx:1976-1994`: on a plane the curve just rotates with the surface.
    if let Some(pl) = s.gp_pln() {
        let mut axe = GpAx1::new(pl.location(), *pl.axis().direction());
        let mut theta = 0.0;
        if find_rotation(&pl, oris, direction, angle, neutral_plane, &mut axe, &mut theta) {
            if theta.abs() > ANGULAR {
                let mut trsf = occt_core::gp::GpTrsf::identity();
                if trsf.set_rotation_ax1(&axe, theta).is_err() {
                    return None;
                }
                return Some(Arc::from(c.transformed(&trsf)));
            }
            return Some(c.clone());
        }
        return None;
    }

    // `cxx:1996-1998`: only lines are drafted analytically.
    let lin = c.gp_line()?;

    // `cxx:2005-2027`: the surface normal at the line's origin.
    let mut norm = GpDir::default_dir();
    if let Some(cy) = s.gp_cylinder() {
        let (u, v) = slib::cylinder_parameters(&cy.position(), &lin.location());
        let (_p, d1u, d1v) = surface_eval::cylinder_d1(&cy, u, v);
        norm = GpDir::from_vec(&d1u.crossed(&d1v)).ok()?;
    } else if let Some(co) = s.gp_cone() {
        let (u, v) =
            slib::cone_parameters(&co.position(), co.radius(), co.semi_angle(), &lin.location());
        let (_p, d1u, d1v) = surface_eval::cone_d1(&co, u, v);
        norm = GpDir::from_vec(&d1u.crossed(&d1v)).ok()?;
    }

    // `cxx:2029-2042`
    let inter = IntAnaIntConicQuad::line_plane(&lin, neutral_plane, ANGULAR, 0.0, 0.0);
    if inter.is_done() && inter.nb_points() != 0 {
        if oris == Orientation::Reversed {
            norm.reverse();
        }
        let pt = *inter.point(1);
        let axis = norm.crossed(direction).ok()?;
        let axrot = GpAx1::new(pt, axis);
        let mut lires = GpLin::from_pnt_dir(pt, *direction);
        lires.rotate(&axrot, angle);
        if lires.direction().dot(&lin.direction()) < 0.0 {
            lires.set_direction(lires.direction().reversed());
        }
        return Some(Arc::new(GeomLine::new(lires)));
    }
    None
}
