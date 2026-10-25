//! `ChFi3d_Parameters` — analytic (u, v) of a 3D point on a surface.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_0.cxx:4671-4703`.

use occt_core::elib::slib;
use occt_core::gp::GpPnt;
use occt_core::precision;
use occt_geom::geom_api::project_point_on_surface;
use occt_geom::Surface;

use crate::geom_bnd_lib_surface3d::{surface_kind, GeomAbsSurfaceType};

/// `ChFi3d_Parameters(S, p3d, u, v)` (`cxx:4671-4703`): `(u, v)` of `p3d` on
/// `S`, exactly on the five analytic surface types and through
/// `GeomAPI_ProjectPointOnSurf` otherwise.
///
/// OCCT throws `StdFail_NotDone` when `NbPoints() != 1` (`cxx:4696-4700`); the
/// port's [`project_point_on_surface`] answers only the nearest solution, so
/// `Ok` also covers a multi-solution projection and `Err` is the no-solution
/// case.
pub fn parameters(s: &dyn Surface, p3d: &GpPnt) -> Result<(f64, f64), String> {
    match surface_kind(s) {
        GeomAbsSurfaceType::Plane => {
            let pln = s.gp_pln().expect("plane surface");
            Ok(slib::plane_parameters(&pln.position(), p3d))
        }
        GeomAbsSurfaceType::Cylinder => {
            let cy = s.gp_cylinder().expect("cylinder surface");
            Ok(slib::cylinder_parameters(&cy.position(), p3d))
        }
        GeomAbsSurfaceType::Cone => {
            let co = s.gp_cone().expect("cone surface");
            Ok(slib::cone_parameters(
                &co.position(),
                co.radius(),
                co.semi_angle(),
                p3d,
            ))
        }
        GeomAbsSurfaceType::Sphere => {
            let sp = s.gp_sphere().expect("sphere surface");
            Ok(slib::sphere_parameters(&sp.position(), p3d))
        }
        GeomAbsSurfaceType::Torus => {
            let to = s.gp_torus().expect("torus surface");
            Ok(slib::torus_parameters(
                &to.position(),
                to.major_radius(),
                to.minor_radius(),
                p3d,
            ))
        }
        _ => {
            let ps = project_point_on_surface(s, p3d, precision::CONFUSION).ok_or_else(|| {
                "StdFail_NotDone: ChFi3d_Parameters - no projection results".to_string()
            })?;
            Ok((ps.u, ps.v))
        }
    }
}
