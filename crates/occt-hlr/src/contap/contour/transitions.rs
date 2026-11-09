//! The three transition helpers `Contap_Contour` uses on its lines.
//!
//! Sources: `ComputeTransitionOnLine` (`Contap_Contour.cxx:893-956`),
//! `ComputeTransitionOngpLine` (`cxx:2100-2126`),
//! `ComputeTransitionOngpCircle` (`cxx:2128-2154`).

use occt_core::elib::clib;
use occt_core::elib::slib;
use occt_core::gp::{GpCirc, GpLin, GpVec};
use occt_core::precision::REAL_EPSILON;
use occt_math::function_set_root::MathFunctionSetWithDerivatives;
use occt_math::{MathMatrix, MathVector};

use occt_topo::geom_bnd_lib_surface3d::GeomAbsSurfaceType;
use occt_topo::intsurf_transition::IntSurfTypeTrans;

use super::super::surf_function::ContapSurfFunction;

/// `ComputeTransitionOnLine(SFunc, u, v, tgline)` (`cxx:893-956`).
///
/// The tangent is expressed in the `(u, v)` basis through the squared system of
/// the two first derivatives, then the gradient of the contour function is
/// combined with it.
pub(crate) fn compute_transition_on_line(
    sfunc: &mut ContapSurfFunction,
    u: f64,
    v: f64,
    tgline: &GpVec,
) -> IntSurfTypeTrans {
    let Some(surf) = sfunc.surface().cloned() else {
        return IntSurfTypeTrans::Undecided;
    };
    let Some((_pntbid, d1u, d1v)) = surf.d1(u, v) else {
        return IntSurfTypeTrans::Undecided;
    };

    let d1u_t = d1u.dot(tgline);
    let d1v_t = d1v.dot(tgline);
    let normu2 = d1u.dot(&d1u);
    let normv2 = d1v.dot(&d1v);
    let d1ud1v = d1u.dot(&d1v);
    let det = normu2 * normv2 - d1ud1v * d1ud1v;
    if det < REAL_EPSILON {
        // "On ne doit pas passer ici !!" (`cxx:916-920`).
        return IntSurfTypeTrans::Undecided;
    }

    let alpha = (d1u_t * normv2 - d1v_t * d1ud1v) / det;
    let beta = (normu2 * d1v_t - d1ud1v * d1u_t) / det;

    // `cxx:929-936`: the gradient of the contour function.
    let mut x = MathVector::new(1, 2);
    let mut df = MathMatrix::new(1, 1, 1, 2);
    x.set_value(1, u);
    x.set_value(2, v);
    sfunc.derivatives(&x, &mut df);
    let v1 = df.value(1, 1);
    let v2 = df.value(1, 2);

    let det = -v1 * beta + v2 * alpha;
    if det < REAL_EPSILON {
        return IntSurfTypeTrans::Undecided;
    }
    if det > 0.0 {
        return IntSurfTypeTrans::Out;
    }
    IntSurfTypeTrans::In
}

/// `ComputeTransitionOngpLine(SFunc, L)` (`cxx:2100-2126`).
pub(crate) fn compute_transition_on_gp_line(
    sfunc: &mut ContapSurfFunction,
    l: &GpLin,
) -> IntSurfTypeTrans {
    let Some(surf) = sfunc.surface().cloned() else {
        return IntSurfTypeTrans::Undecided;
    };
    let typ_s = surf.get_type();
    let (p, t) = clib::line_d1(l, 0.0);
    let (u, v) = match typ_s {
        GeomAbsSurfaceType::Cylinder => match surf.cylinder() {
            Some(cy) => slib::cylinder_parameters(&cy.position(), &p),
            None => (0.0, 0.0),
        },
        GeomAbsSurfaceType::Cone => match surf.cone() {
            Some(co) => slib::cone_parameters(&co.position(), co.radius(), co.semi_angle(), &p),
            None => (0.0, 0.0),
        },
        GeomAbsSurfaceType::Sphere => match surf.sphere() {
            Some(sp) => slib::sphere_parameters(sp.position(), &p),
            None => (0.0, 0.0),
        },
        _ => (0.0, 0.0),
    };
    compute_transition_on_line(sfunc, u, v, &t)
}

/// `ComputeTransitionOngpCircle(SFunc, C)` (`cxx:2128-2154`).
pub(crate) fn compute_transition_on_gp_circle(
    sfunc: &mut ContapSurfFunction,
    c: &GpCirc,
) -> IntSurfTypeTrans {
    let Some(surf) = sfunc.surface().cloned() else {
        return IntSurfTypeTrans::Undecided;
    };
    let typ_s = surf.get_type();
    let (p, t) = clib::circle_d1(c, 0.0);
    let (u, v) = match typ_s {
        GeomAbsSurfaceType::Cylinder => match surf.cylinder() {
            Some(cy) => slib::cylinder_parameters(&cy.position(), &p),
            None => (0.0, 0.0),
        },
        GeomAbsSurfaceType::Cone => match surf.cone() {
            Some(co) => slib::cone_parameters(&co.position(), co.radius(), co.semi_angle(), &p),
            None => (0.0, 0.0),
        },
        GeomAbsSurfaceType::Sphere => match surf.sphere() {
            Some(sp) => slib::sphere_parameters(sp.position(), &p),
            None => (0.0, 0.0),
        },
        _ => (0.0, 0.0),
    };
    compute_transition_on_line(sfunc, u, v, &t)
}
