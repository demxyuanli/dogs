//! `Function_SetUVBounds` — the parameter window a 3d curve spans on an
//! analytic surface, together with the two seam flags.
//!
//! Source: `ModelingData/TKGeomBase/ProjLib/ProjLib_ComputeApprox.cxx:262-1019`
//! (OCCT 8.0.0). One arm per surface type, each in its own file:
//! [`cone`] (`cxx:287-461`), [`cylinder`] (`cxx:463-591`), [`sphere`]
//! (`cxx:593-889`), [`torus`] (`cxx:890-1012`). The `default:` arm
//! (`cxx:1013-1018`) only clears the two flags.
//!
//! `UCouture` / `VCouture` mean "the curve runs around the period, so the 2d
//! result has to be opened at the seam". The windows are the ones
//! `Function_Value` (`ProjLib_ComputeApprox.cxx:103-180`) wraps into.

mod cone;
mod cylinder;
mod sphere;
mod torus;

use crate::curve::Curve;
use crate::projlib::projected_curve::surface_type;
use crate::surface::Surface;
use occt_core::kernel::geomabs::SurfaceType;

/// The five out-parameters of `Function_SetUVBounds` (`cxx:262-269`) plus the
/// two flags: `myU1`, `myU2`, `myV1`, `myV2`, `UCouture`, `VCouture`.
#[derive(Debug, Clone, Copy)]
pub struct UvBounds {
    /// `myU1` (`cxx:263`).
    pub u1: f64,
    /// `myU2` (`cxx:264`).
    pub u2: f64,
    /// `myV1` (`cxx:265`).
    pub v1: f64,
    /// `myV2` (`cxx:266`).
    pub v2: f64,
    /// `UCouture` (`cxx:267`).
    pub u_couture: bool,
    /// `VCouture` (`cxx:268`).
    pub v_couture: bool,
}

impl UvBounds {
    /// The seed the `ProjLib_Function` constructor hands in
    /// (`cxx:1040-1048`): every value zero, both flags `false`.
    pub fn new() -> Self {
        Self {
            u1: 0.0,
            u2: 0.0,
            v1: 0.0,
            v2: 0.0,
            u_couture: false,
            v_couture: false,
        }
    }
}

impl Default for UvBounds {
    fn default() -> Self {
        Self::new()
    }
}

/// `Function_SetUVBounds(...)` (`cxx:262-1019`).
pub(crate) fn function_set_uv_bounds(curve: &dyn Curve, surface: &dyn Surface) -> UvBounds {
    let mut bounds = UvBounds::new();
    let mut w1 = curve.first_parameter();
    let mut w2 = curve.last_parameter();
    let mut w = 0.5 * (w1 + w2);
    let mut p1 = curve.d0(w1);
    let mut p2 = curve.d0(w2);
    let mut p = curve.d0(w);

    match surface_type(surface) {
        SurfaceType::Cone => cone::set_uv_bounds(
            &mut bounds,
            curve,
            surface,
            &mut w1,
            &mut w2,
            &mut w,
            &mut p1,
            &mut p2,
            &mut p,
        ),
        SurfaceType::Cylinder => cylinder::set_uv_bounds(
            &mut bounds,
            curve,
            surface,
            w1,
            w2,
            &p1,
            &p2,
            &p,
        ),
        SurfaceType::Sphere => sphere::set_uv_bounds(
            &mut bounds,
            curve,
            surface,
            w1,
            w2,
            &p1,
            &mut p2,
            &p,
        ),
        SurfaceType::Torus => torus::set_uv_bounds(&mut bounds, curve, surface, w1, w2, &p1),
        _ => {
            bounds.u_couture = false;
            bounds.v_couture = false;
        }
    }

    bounds
}
