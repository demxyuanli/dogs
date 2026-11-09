//! File-local `IsOnSurface` (`BRepOffset_Tool.cxx:642-743`), the curve
//! predicate `TryProject` (`:1945-2006`) uses.
//!
//! `IsOnSurface` samples the curve at `n = 5` parameters spaced by
//! `du = (f - l) / (n - 1)` - note the *negative* step, which walks from the
//! first parameter down past the last one; that is the reference behaviour and
//! is reproduced here. For every sample it evaluates the surface's own
//! parameter inversion and re-evaluates the point, so `TolReached` ends up
//! holding the last sample's distance, not the maximum (`cxx:666`, `:681`,
//! `:696`, `:711`, `:726`).
//!
//! `ElSLib::PlaneParameters` / `CylinderParameters` / ... and the matching
//! `*Value` inversions are `occt_core::elib::slib`; a surface type outside the
//! five elementary ones answers `None` (`cxx:737-739`).
//!
//! `ToSmall` (`cxx:622-638`): the "curve is a point" predicate `PipeInter`
//! (`cxx:769`) and `InterOrExtent` (`cxx:2010-2068`) use.

use std::sync::Arc;

use occt_core::elib::slib;
use occt_geom::{Curve, Surface};

use crate::geom_bnd_lib_surface3d::GeomAbsSurfaceType;
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

/// `IsOnSurface(C, S, TolConf, TolReached)` (`cxx:642-743`): `Some(TolReached)`
/// for `true`, `None` for `false`.
pub(super) fn is_on_surface(c: &dyn Curve, s: Arc<dyn Surface>, tol_conf: f64) -> Option<f64> {
    let f = c.first_parameter();
    let l = c.last_parameter();
    let n: i32 = 5;
    let du = (f - l) / (n - 1) as f64;
    let mut tol_reached = 0.0;

    // `GeomAdaptor_Surface AS(S)` (`cxx:657`).
    let a_s = SurfaceAdaptor::of_surface(s);
    match a_s.get_type() {
        GeomAbsSurfaceType::Plane => {
            // `cxx:660-671`.
            let pln = a_s.plane()?;
            let ax = pln.position();
            for i in 0..n {
                let p = c.d0(f + f64::from(i) * du);
                let (u, v) = slib::plane_parameters(&ax, &p);
                tol_reached = p.distance(&slib::plane_value(&pln, u, v));
                if tol_reached > tol_conf {
                    return None;
                }
            }
        }
        GeomAbsSurfaceType::Cylinder => {
            // `cxx:673-687`. `ElSLib::CylinderParameters` reads only `Ax`; the
            // radius argument of the OCCT call shape is unused there.
            let cyl = a_s.cylinder()?;
            let ax = cyl.position();
            for i in 0..n {
                let p = c.d0(f + f64::from(i) * du);
                let (u, v) = slib::cylinder_parameters(&ax, &p);
                tol_reached = p.distance(&slib::cylinder_value(&cyl, u, v));
                if tol_reached > tol_conf {
                    return None;
                }
            }
        }
        GeomAbsSurfaceType::Cone => {
            // `cxx:689-703`.
            let cone = a_s.cone()?;
            let ax = cone.position();
            let rad = cone.radius();
            let alp = cone.semi_angle();
            for i in 0..n {
                let p = c.d0(f + f64::from(i) * du);
                let (u, v) = slib::cone_parameters(&ax, rad, alp, &p);
                tol_reached = p.distance(&slib::cone_value(&cone, u, v));
                if tol_reached > tol_conf {
                    return None;
                }
            }
        }
        GeomAbsSurfaceType::Sphere => {
            // `cxx:705-719`.
            let sph = a_s.sphere()?;
            let ax = *sph.position();
            for i in 0..n {
                let p = c.d0(f + f64::from(i) * du);
                let (u, v) = slib::sphere_parameters(&ax, &p);
                tol_reached = p.distance(&slib::sphere_value(&sph, u, v));
                if tol_reached > tol_conf {
                    return None;
                }
            }
        }
        GeomAbsSurfaceType::Torus => {
            // `cxx:721-733`.
            let tor = a_s.torus()?;
            let ax = *tor.position();
            let r1 = tor.major_radius();
            let r2 = tor.minor_radius();
            for i in 0..n {
                let p = c.d0(f + f64::from(i) * du);
                let (u, v) = slib::torus_parameters(&ax, r1, r2, &p);
                tol_reached = p.distance(&slib::torus_value(&tor, u, v));
                if tol_reached > tol_conf {
                    return None;
                }
            }
        }
        // `cxx:736-739`.
        _ => return None,
    }
    Some(tol_reached)
}

/// `ToSmall(C)` (`cxx:622-638`): true when the three samples at the two end
/// parameters and at `f * 0.668 + l * 0.332` all stay within
/// `10 * Precision::Confusion()` of each other.
pub(super) fn to_small(c: &dyn Curve) -> bool {
    const TOL: f64 = 10.0 * occt_core::precision::CONFUSION;
    let f = c.first_parameter();
    let l = c.last_parameter();
    let m = f * 0.668 + l * 0.332;
    let p1 = c.d0(f);
    let p2 = c.d0(l);
    let p3 = c.d0(m);
    if p1.distance(&p2) > TOL {
        return false;
    }
    if p2.distance(&p3) > TOL {
        return false;
    }
    true
}
