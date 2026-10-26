//! `ProjLib_ComputeApprox` file statics: the `(u, v)` evaluation of a 3d point
//! on a surface, its first derivative along the curve, and the sampling step.
//!
//! Source: `ModelingData/TKGeomBase/ProjLib/ProjLib_ComputeApprox.cxx:51-258`
//! (OCCT 8.0.0): `PlaneProjector`, `IsEqual`, `Function_Value`,
//! `Function_D1`, `Function_ComputeStep`.

use occt_core::elib::{clib, slib};
use occt_core::gp::{GpAx3, GpPnt, GpPnt2d, GpVec2d};
use occt_core::kernel::geomabs::SurfaceType;
use occt_core::precision;

use crate::curve::Curve;
use crate::gcpnts::curve_length;
use crate::projlib::projected_curve::surface_type;
use crate::surface::Surface;

/// `PlaneProjector` (`cxx:54-89`): the plane basis pre-computed once for a
/// batch of point projections.
pub(super) struct PlaneProjector {
    ox: f64,
    oy: f64,
    oz: f64,
    dxx: f64,
    dxy: f64,
    dxz: f64,
    dyx: f64,
    dyy: f64,
    dyz: f64,
}

impl PlaneProjector {
    /// `PlaneProjector(const gp_Ax3&)` (`cxx:60-77`).
    pub(super) fn new(the_pos: &GpAx3) -> Self {
        let loc = the_pos.location();
        let x_dir = the_pos.x_direction();
        let y_dir = the_pos.y_direction();
        Self {
            ox: loc.x(),
            oy: loc.y(),
            oz: loc.z(),
            dxx: x_dir.x(),
            dxy: x_dir.y(),
            dxz: x_dir.z(),
            dyx: y_dir.x(),
            dyy: y_dir.y(),
            dyz: y_dir.z(),
        }
    }

    /// `Project(const gp_Pnt&)` (`cxx:81-88`).
    pub(super) fn project(&self, the_p: &GpPnt) -> GpPnt2d {
        let dx = the_p.x() - self.ox;
        let dy = the_p.y() - self.oy;
        let dz = the_p.z() - self.oz;
        GpPnt2d::new(
            dx * self.dxx + dy * self.dxy + dz * self.dxz,
            dx * self.dyx + dy * self.dyy + dz * self.dyz,
        )
    }
}

/// `IsEqual(Check, With, Toler)` (`cxx:96-99`).
pub(super) fn is_equal(check: f64, with: f64, toler: f64) -> bool {
    (check - with).abs() < toler
}

/// `Function_Value(U, myCurve, mySurface, U1, U2, V1, V2, UCouture, VCouture)`
/// (`cxx:103-180`).
pub(super) fn function_value(
    u: f64,
    curve: &dyn Curve,
    surface: &dyn Surface,
    u1: f64,
    u2: f64,
    v1: f64,
    v2: f64,
    u_couture: bool,
    v_couture: bool,
) -> GpPnt2d {
    let p3d = curve.d0(u);
    let stype = surface_type(surface);

    let (mut s, mut t) = match stype {
        SurfaceType::Plane => slib::plane_parameters(&surface.gp_pln().expect("Plane").pos, &p3d),
        SurfaceType::Cylinder => {
            let cy = surface.gp_cylinder().expect("Cylinder");
            slib::cylinder_parameters(&cy.pos, &p3d)
        }
        SurfaceType::Cone => {
            let co = surface.gp_cone().expect("Cone");
            slib::cone_parameters(&co.pos, co.radius, co.semi_angle, &p3d)
        }
        SurfaceType::Sphere => slib::sphere_parameters(&surface.gp_sphere().expect("Sphere").pos, &p3d),
        SurfaceType::Torus => {
            let to = surface.gp_torus().expect("Torus");
            slib::torus_parameters(&to.pos, to.major_radius, to.minor_radius, &p3d)
        }
        _ => panic!("Standard_NoSuchObject: ProjLib_ComputeApprox::Value"),
    };

    if u_couture && (s < u1 || s > u2) {
        s = clib::in_period(s, u1, u2);
    }

    if v_couture {
        if stype == SurfaceType::Sphere && (s - u1).abs() > std::f64::consts::PI {
            t = std::f64::consts::PI - t;
            s = std::f64::consts::PI + s;
        }
        if t < v1 || t > v2 {
            t = clib::in_period(t, v1, v2);
        }
    }

    GpPnt2d::new(s, t)
}

/// `Function_D1(U, P, D, ...)` (`cxx:183-238`).
pub(super) fn function_d1(
    u: f64,
    curve: &dyn Curve,
    surface: &dyn Surface,
    u1: f64,
    u2: f64,
    v1: f64,
    v2: f64,
    u_couture: bool,
    v_couture: bool,
) -> Option<(GpPnt2d, GpVec2d)> {
    let p = function_value(u, curve, surface, u1, u2, v1, v2, u_couture, v_couture);

    let stype = surface_type(surface);
    match stype {
        SurfaceType::Plane
        | SurfaceType::Cone
        | SurfaceType::Cylinder
        | SurfaceType::Sphere
        | SurfaceType::Torus => {
            let (_, t) = curve.d1(u);
            let (_, d1u, d1v) = surface.d1(p.x(), p.y());

            let mut du = t.dot(&d1u);
            let mut dv = t.dot(&d1v);
            let nu = d1u.square_magnitude();
            let nv = d1v.square_magnitude();

            if nu < precision::epsilon(1.0) || nv < precision::epsilon(1.0) {
                return None;
            }

            du /= nu;
            dv /= nv;
            Some((p, GpVec2d::new(du, dv)))
        }
        _ => None,
    }
}

/// `Function_ComputeStep(myCurve, R)` (`cxx:240-258`).
pub(super) fn function_compute_step(curve: &dyn Curve, r: f64) -> f64 {
    let step0 = 0.1;
    let w1 = curve.first_parameter();
    let w2 = curve.last_parameter();
    let l = curve_length(curve);
    let mut nbp = (l / (r * std::f64::consts::FRAC_PI_4)) as i32 + 1;
    nbp = nbp.max(3);
    let mut step = (w2 - w1) / (nbp - 1) as f64;
    if step > step0 {
        step = step0;
        nbp = ((w2 - w1) / step) as i32 + 1;
        nbp = nbp.max(3);
        step = (w2 - w1) / (nbp - 1) as f64;
    }

    step
}
