//! `ChFi3d_Parameters` - the `(u, v)` of a 3D point on a surface.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_0.cxx:4671-4705`
//! (OCCT 8.0.0). The `GeomAbs_*` arms are `ElSLib::Parameters`; `Bezier`,
//! `BSpline` and every other type fall to `GeomAPI_ProjectPointOnSurf tool(p3d,
//! S)`, which OCCT requires to answer **exactly one** solution
//! (`cxx:4694-4701`, `StdFail_NotDone`) - the `Err` of this port. The tool's
//! constructor is `Init(P, Surface)` (`GeomAPI_ProjectPointOnSurf.cxx:106-131`):
//! `Tolerance = Precision::Confusion()` for both parameters, window = the
//! surface's own bounds, which is `Extrema_ExtPS` over that same range
//! ([`ExtPs::with_surface`]).

use std::sync::Arc;

use occt_core::elib::slib;
use occt_core::gp::GpPnt;
use occt_core::precision::CONFUSION;
use occt_geom::extrema_surf::ExtPs;
use occt_geom::Surface;

use crate::geom_bnd_lib_surface3d::GeomAbsSurfaceType;
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

/// `void ChFi3d_Parameters(const occ::handle<Geom_Surface>& S, const gp_Pnt&
/// p3d, double& u, double& v)` (`cxx:4671-4705`).
pub fn ch_fi3d_parameters(s: &Arc<dyn Surface>, p3d: &GpPnt) -> Result<(f64, f64), String> {
    // `GeomAdaptor_Surface gas(S); switch (gas.GetType())` (`cxx:4673-4675`).
    let gas = SurfaceAdaptor::of_surface(Arc::clone(s));
    match gas.get_type() {
        GeomAbsSurfaceType::Plane => {
            let pln = gas
                .plane()
                .ok_or_else(|| "ChFi3d_Parameters: the plane surface has no gp_Pln".to_string())?;
            Ok(slib::plane_parameters(&pln.pos, p3d))
        }
        GeomAbsSurfaceType::Cylinder => {
            let cyl = gas.cylinder().ok_or_else(|| {
                "ChFi3d_Parameters: the cylinder surface has no gp_Cylinder".to_string()
            })?;
            Ok(slib::cylinder_parameters(&cyl.pos, p3d))
        }
        GeomAbsSurfaceType::Cone => {
            let cone = gas.cone().ok_or_else(|| {
                "ChFi3d_Parameters: the cone surface has no gp_Cone".to_string()
            })?;
            Ok(slib::cone_parameters(
                &cone.pos,
                cone.radius,
                cone.semi_angle,
                p3d,
            ))
        }
        GeomAbsSurfaceType::Sphere => {
            let sph = s.gp_sphere().ok_or_else(|| {
                "ChFi3d_Parameters: the sphere surface has no gp_Sphere".to_string()
            })?;
            Ok(slib::sphere_parameters(&sph.pos, p3d))
        }
        GeomAbsSurfaceType::Torus => {
            let tor = s.gp_torus().ok_or_else(|| {
                "ChFi3d_Parameters: the torus surface has no gp_Torus".to_string()
            })?;
            Ok(slib::torus_parameters(
                &tor.pos,
                tor.major_radius,
                tor.minor_radius,
                p3d,
            ))
        }
        // `case GeomAbs_BezierSurface: case GeomAbs_BSplineSurface: default:`
        // (`cxx:4692-4702`).
        _ => {
            let tool = ExtPs::with_surface(p3d, s.as_ref(), CONFUSION, CONFUSION);
            if tool.nb_ext() != 1 {
                return Err("ChFi3d_Parameters() - no projection results".to_string());
            }
            let (u, v, _) = tool.point(1);
            Ok((u, v))
        }
    }
}
