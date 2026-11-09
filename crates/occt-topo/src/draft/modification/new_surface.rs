//! `Draft_Modification::NewSurface` - the private overload
//! (`Draft_Modification_1.cxx:1780-1960`) - and its helper `FindRotation`
//! (`Draft_Modification_1.cxx:2403-2453`).
//!
//! # Translation notes
//!
//! * `IntAna_QuadQuadGeo::Perform(gp_Pln, gp_Cylinder, Tolang, Tol)` defaults
//!   `H` to `0` (`IntAna_QuadQuadGeo.hxx:91-94`), so the port passes `0.0` as
//!   the face height to [`quadric_quadric_plane_cylinder`].
//! * The ellipse guard of the cylinder arm (`cxx:1831-1838`) computes
//!   `isIntDone` and then the very next test rejects anything that is not a
//!   circle, so the whole arm answers `None` for an ellipse in the reference
//!   as well; the expression is kept because it is the reference's control
//!   flow, not a predicate invented here.
//! * `gp_Cone(Ax3, SemiAngle, Radius)` maps to [`GpCone::new`] `(pos, radius,
//!   semi_angle)` - the argument order of the port's constructor is the
//!   reverse of OCCT's.

use std::sync::Arc;

use occt_core::elib::clib;
use occt_core::gp::{GpAx1, GpCone, GpCylinder, GpDir, GpPln};
use occt_core::precision::{ANGULAR, CONFUSION};
use occt_geom::cone::GeomCone;
use occt_geom::cylinder::GeomCylinder;
use occt_geom::intana::{
    quadric_quadric_plane_cone, quadric_quadric_plane_cylinder, quadric_quadric_planes,
    QuadricIntersection,
};
use occt_geom::Surface;

use crate::abs::Orientation;

/// `static bool FindRotation(Pl, Oris, Direction, Angle, NeutralPlane, Axe,
/// theta)` (`Draft_Modification_1.cxx:2403-2453`): the rotation axis and angle
/// that tilt the plane `Pl` by `Angle` around its intersection line with the
/// neutral plane.
///
/// The `while (abs(theta) > M_PI)` loop folds `theta` into `[-pi, pi]`
/// (`cxx:2442-2445`). The reference answers `false` - leaving `Axe` / `theta`
/// untouched - whenever the two planes do not meet in a line or the direction
/// is too close to the intersection line.
pub fn find_rotation(
    pl: &GpPln,
    oris: Orientation,
    direction: &GpDir,
    angle: f64,
    neutral_plane: &GpPln,
    axe: &mut GpAx1,
    theta: &mut f64,
) -> bool {
    // `cxx:2410`
    let inter = quadric_quadric_planes(pl, neutral_plane, ANGULAR, CONFUSION);
    // `cxx:2412`
    let QuadricIntersection::Line(li) = inter else {
        return false;
    };
    // `cxx:2415-2416`
    let nx = li.direction();
    let Ok(ny) = pl.axis().direction().crossed(&nx) else {
        return false;
    };
    // `cxx:2418`
    let a = direction.dot(&nx);
    if a.abs() <= 1.0 - ANGULAR {
        // `cxx:2422-2423`
        let mut b = direction.dot(&ny);
        let mut c = direction.dot(&pl.axis().direction());
        let direct = pl.position().is_direct();
        if (direct && oris == Orientation::Reversed)
            || (!direct && oris == Orientation::Forward)
        {
            b = -b;
            c = -c;
        }
        // `cxx:2430-2432`
        let denom = (1.0 - a * a).sqrt();
        let sina = angle.sin();
        if denom > sina.abs() {
            let phi = (b / denom).atan2(c / denom);
            let theta0 = (sina / denom).acos();
            let mut th = theta0 - phi;
            if th.cos() < 0.0 {
                th = -theta0 - phi;
            }
            // `cxx:2442-2445`
            while th.abs() > std::f64::consts::PI {
                th += std::f64::consts::PI * if th < 0.0 { 1.0 } else { -1.0 };
            }
            // `cxx:2446-2447`
            *axe = *li.position();
            *theta = th;
            return true;
        }
    }
    false
}

/// `Draft_Modification::NewSurface(S, Oris, Direction, Angle, NeutralPlane)`
/// (`Draft_Modification_1.cxx:1780-1960`): the draft surface of `s`, or `None`
/// (the reference's null handle) when the support is not a plane, a cylinder or
/// a cone, or when the direction / neutral plane are too close to
/// perpendicular for the construction to be defined.
pub fn new_surface(
    s: &Arc<dyn Surface>,
    oris: Orientation,
    direction: &GpDir,
    angle: f64,
    neutral_plane: &GpPln,
) -> Option<Arc<dyn Surface>> {
    // `cxx:1791-1807`: plane - rotate about the neutral-plane intersection line.
    if let Some(pl) = s.gp_pln() {
        let mut axe = GpAx1::new(pl.location(), *pl.axis().direction());
        let mut theta = 0.0;
        if find_rotation(&pl, oris, direction, angle, neutral_plane, &mut axe, &mut theta) {
            if theta.abs() > ANGULAR {
                // `S->Rotated(Axe, Theta)`.
                let mut trsf = occt_core::gp::GpTrsf::identity();
                if trsf.set_rotation_ax1(&axe, theta).is_err() {
                    return None;
                }
                return Some(Arc::from(s.transformed(&trsf)));
            }
            return Some(s.clone());
        }
        return None;
    }

    // `cxx:1808-1877`: cylinder.
    if let Some(cy) = s.gp_cylinder() {
        // `cxx:1809-1818`
        let testdir = direction.dot(&neutral_plane.axis().direction());
        if testdir.abs() <= 1.0 - ANGULAR {
            return None;
        }
        // `cxx:1819-1826`
        let testdir = direction.dot(&cy.axis().direction());
        if testdir.abs() <= 1.0 - ANGULAR {
            return None;
        }
        // `cxx:1827-1852`
        if angle.abs() > ANGULAR {
            let inter = quadric_quadric_plane_cylinder(neutral_plane, &cy, ANGULAR, CONFUSION, 0.0);
            let mut is_int_done = !matches!(inter, QuadricIntersection::NotDone);
            if let QuadricIntersection::Ellipse(el) = &inter {
                is_int_done = el.major_radius() < 100000.0 * el.minor_radius();
            }
            if !is_int_done {
                return None;
            }
            let QuadricIntersection::Circle(circ) = &inter else {
                return None;
            };
            let axcone = cy.position();
            let mut alpha = angle;
            let direct = axcone.is_direct();
            if (direct && oris == Orientation::Reversed)
                || (!direct && oris == Orientation::Forward)
            {
                alpha = -alpha;
            }
            let center = circ.location();
            if testdir < 0.0 {
                alpha = -alpha;
            }
            let z = clib::line_parameter_ax1(&cy.axis(), &center);
            let mut rad = cy.radius() + z * alpha.tan();
            if rad < 0.0 {
                rad = -rad;
            } else {
                alpha = -alpha;
            }
            // `cxx:1871-1873`: `gp_Cone(axcone, alpha, Rad)`.
            let co = GpCone::new(axcone, rad, alpha).ok()?;
            return Some(Arc::new(GeomCone::new(co)));
        }
        return Some(s.clone());
    }

    // `cxx:1879-1951`: cone.
    if let Some(co1) = s.gp_cone() {
        // `cxx:1882-1888`
        let testdir = direction.dot(&neutral_plane.axis().direction());
        if testdir.abs() <= 1.0 - ANGULAR {
            return None;
        }
        // `cxx:1894-1901`
        let testdir = direction.dot(&co1.axis().direction());
        if testdir.abs() <= 1.0 - ANGULAR {
            return None;
        }
        // `cxx:1903-1911`
        let inter = quadric_quadric_plane_cone(neutral_plane, &co1, ANGULAR, CONFUSION);
        if matches!(inter, QuadricIntersection::NotDone) {
            return None;
        }
        let QuadricIntersection::Circle(circ) = &inter else {
            return None;
        };
        // `cxx:1912-1919`
        let axcone = co1.position();
        let mut alpha = angle;
        let direct = axcone.is_direct();
        if (direct && oris == Orientation::Reversed)
            || (!direct && oris == Orientation::Forward)
        {
            alpha = -alpha;
        }
        // `cxx:1920-1950`
        let center = circ.location();
        if angle.abs() > ANGULAR {
            if testdir < 0.0 {
                alpha = -alpha;
            }
            let z = clib::line_parameter_ax1(&co1.axis(), &center);
            let mut rad = circ.radius() + z * alpha.tan();
            if rad < 0.0 {
                rad = -rad;
            } else {
                alpha = -alpha;
            }
            if (alpha - co1.semi_angle()).abs() < ANGULAR {
                return Some(s.clone());
            }
            let co = GpCone::new(axcone, rad, alpha).ok()?;
            return Some(Arc::new(GeomCone::new(co)));
        }
        let cy = GpCylinder::new(axcone, circ.radius()).ok()?;
        return Some(Arc::new(GeomCylinder::new(cy)));
    }

    // `cxx:1952-1958`: "Draft_SurfNotYetImplemented".
    None
}
