//! Port of the `GeomBndLib_Surface` dispatch (`GeomBndLib_Surface.cxx:44-108`,
//! `116-265`, `320-364`) plus the concrete evaluators the analytic modules do
//! not cover: `GeomBndLib_OtherSurface.cxx`, `GeomBndLib_BSplineSurface.cxx`,
//! `GeomBndLib_BezierSurface.cxx`, `GeomBndLib_SurfaceOfRevolution.cxx`,
//! `GeomBndLib_SurfaceOfExtrusion.cxx`, `GeomBndLib_OffsetSurface.cxx`, and
//! `GeomBndLib_SamplingHelpers.pxx`.
//!
//! PARKED (no geometry type in this port - reported, not fudged):
//! - `GeomBndLib_BSplineSurface::BoxOptimal` / `OtherSurface::BoxOptimal`
//!   (`OptimizationHelpers.pxx` PSO+ Powell): `BRepBndLib::Add` only calls
//!   `Add`/`Box`, so the `AddOptimal` path is out of scope here. The
//!   `GeomBndLib_SurfaceOfExtrusion::BoxOptimal` overload
//!   (`SurfaceOfExtrusion.cxx:224-283`) belongs to that same path and is not
//!   ported either.
//! - The `catch (Standard_Failure)` fallbacks of
//!   `SurfaceOfExtrusion.cxx:163-167` / `:212-217` have no counterpart: this
//!   port's curve evaluators return empty results instead of throwing.
//! - The `GeomAdaptor_Surface` arms that rebuild a primitive from the adaptor
//!   when `initFromSurfaceType` fails (`Surface.cxx:74-107`) have no separate
//!   adaptor object in this port; the registered surface is used directly.

use std::sync::Arc;

use occt_core::bnd::BndBox;
use occt_core::gp::{GpAx2, GpCirc, GpDir, GpPnt};
use occt_core::precision::{Precision, CONFUSION};

use occt_geom::bspline_surface::GeomBSplineSurface;
use occt_geom::surface::Surface;

use crate::geom_bnd_lib_analytic3d as analytic;
use crate::geom_bnd_lib_circle3d::box_circ_range;
use crate::geom_bnd_lib_curve3d::box_curve;
use crate::geom_bnd_lib_inf3d::{open_max, open_min, open_min_max};
use crate::geom_bnd_lib_spline_helpers::compute_pole_index_range;

/// `Precision::Parametric(Precision::Confusion())` = `Confusion() * 0.01`
/// (`Precision.hxx:328-334`).
const PARAMETRIC_CONFUSION: f64 = CONFUSION * 0.01;

/// `GeomAbs_SurfaceType` as far as the dispatcher needs it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SurfaceKind {
    Plane,
    Cylinder,
    Cone,
    Sphere,
    Torus,
    BezierSurface,
    BSplineSurface,
    SurfaceOfRevolution,
    SurfaceOfExtrusion,
    OffsetSurface,
    OtherSurface,
}

/// `GeomBndLib_Surface::initFromSurface` (`Surface.cxx:136-190`): the concrete
/// geometry type, most specific first. A `Geom_RectangularTrimmedSurface`
/// reports its basis type (`GeomAdaptor_Surface::load`, `Surface.cxx:423-430`
/// recurses into the basis and dispatches there).
pub fn surface_kind(s: &dyn Surface) -> SurfaceKind {
    if let Some(basis) = s.rectangular_trimmed_basis() {
        return surface_kind(basis.as_ref());
    }
    if s.gp_pln().is_some() {
        SurfaceKind::Plane
    } else if s.gp_cylinder().is_some() {
        SurfaceKind::Cylinder
    } else if s.gp_cone().is_some() {
        SurfaceKind::Cone
    } else if s.gp_sphere().is_some() {
        SurfaceKind::Sphere
    } else if s.gp_torus().is_some() {
        SurfaceKind::Torus
    } else if s.is_bspline_surface() {
        SurfaceKind::BSplineSurface
    } else if s.osculating_bspline().is_some() {
        // Only `Geom_BSplineSurface` and `Geom_BezierSurface` answer
        // `osculating_bspline`; the former is taken above.
        SurfaceKind::BezierSurface
    } else if s.is_surface_of_revolution() {
        SurfaceKind::SurfaceOfRevolution
    } else if s.is_surface_of_linear_extrusion() {
        SurfaceKind::SurfaceOfExtrusion
    } else if s.is_offset_surface() {
        SurfaceKind::OffsetSurface
    } else {
        SurfaceKind::OtherSurface
    }
}

/// `unwrapSurface` (`GeomBndLib_Surface.cxx:35-42`): drop nested
/// `Geom_RectangularTrimmedSurface` wrappers.
pub fn unwrap_surface(s: &Arc<dyn Surface>) -> Arc<dyn Surface> {
    match s.rectangular_trimmed_basis() {
        Some(basis) => unwrap_surface(&basis),
        None => s.clone(),
    }
}

/// `Geom_BSplineSurface` view of a spline surface (`Adaptor3d_Surface::BSpline`
/// / `Bezier`). `None` for anything else.
fn spline_of(s: &dyn Surface) -> Option<GeomBSplineSurface> {
    if s.is_bspline_surface() || s.osculating_bspline().is_some() {
        s.osculating_bspline()
    } else {
        None
    }
}

/// `unique knots / mults` of a flat knot vector, plus its unique-knot count.

// ---------------------------------------------------------------------------
// Sampling helpers (`GeomBndLib_SamplingHelpers.pxx:109-222`)
// ---------------------------------------------------------------------------

/// C++ `RealToInt`: `(int)` truncation toward zero.
fn real_to_int(v: f64) -> i32 {
    v as i32
}

/// `ComputeNbUSamples(theSurf, theUMin, theUMax)` (`SamplingHelpers.pxx:109-143`).
fn nb_u_samples(s: &dyn Surface, u_min: f64, u_max: f64) -> i32 {
    let n = match surface_kind(s) {
        SurfaceKind::BezierSurface => {
            let mut n = 2 * s.nb_u_poles();
            let du = u_max - u_min;
            if du < 0.9 {
                n = real_to_int(du * n as f64) + 1;
                n = n.max(5);
            }
            n
        }
        SurfaceKind::BSplineSurface => match spline_of(s) {
            Some(b) => {
                let knots = b.distinct_knots_and_mults_u().0;
                let mut n = 2 * (b.deg_u as i32 + 1) * (knots.len() as i32 - 1);
                let (u_min_g, u_max_g) = (s.u_range().0, s.u_range().1);
                let du = (u_max - u_min) / (u_max_g - u_min_g);
                if du < 0.9 {
                    n = real_to_int(du * n as f64) + 1;
                    n = n.max(5);
                }
                n
            }
            None => 33,
        },
        _ => 33,
    };
    n.min(50)
}

/// `ComputeNbVSamples(theSurf, theVMin, theVMax)` (`SamplingHelpers.pxx:146-180`).
fn nb_v_samples(s: &dyn Surface, v_min: f64, v_max: f64) -> i32 {
    let n = match surface_kind(s) {
        SurfaceKind::BezierSurface => {
            let mut n = 2 * s.nb_v_poles();
            let dv = v_max - v_min;
            if dv < 0.9 {
                n = real_to_int(dv * n as f64) + 1;
                n = n.max(5);
            }
            n
        }
        SurfaceKind::BSplineSurface => match spline_of(s) {
            Some(b) => {
                let knots = b.distinct_knots_and_mults_v().0;
                let mut n = 2 * (b.deg_v as i32 + 1) * (knots.len() as i32 - 1);
                let (v_min_g, v_max_g) = (s.v_range().0, s.v_range().1);
                let dv = (v_max - v_min) / (v_max_g - v_min_g);
                if dv < 0.9 {
                    n = real_to_int(dv * n as f64) + 1;
                    n = n.max(5);
                }
                n
            }
            None => 33,
        },
        _ => 33,
    };
    n.min(50)
}

/// Full-range `ComputeNbUSamples(theSurf)` (`SamplingHelpers.pxx:184-203`).
fn nb_u_samples_full(s: &dyn Surface) -> i32 {
    let n: i32 = match surface_kind(s) {
        SurfaceKind::BezierSurface => 2 * s.nb_u_poles(),
        SurfaceKind::BSplineSurface => match spline_of(s) {
            Some(b) => {
                let knots = b.distinct_knots_and_mults_u().0;
                2 * (b.deg_u as i32 + 1) * (knots.len() as i32 - 1)
            }
            None => 33,
        },
        _ => 33,
    };
    n.min(50)
}

/// Full-range `ComputeNbVSamples(theSurf)` (`SamplingHelpers.pxx:206-222`).
fn nb_v_samples_full(s: &dyn Surface) -> i32 {
    let n: i32 = match surface_kind(s) {
        SurfaceKind::BezierSurface => 2 * s.nb_v_poles(),
        SurfaceKind::BSplineSurface => match spline_of(s) {
            Some(b) => {
                let knots = b.distinct_knots_and_mults_v().0;
                2 * (b.deg_v as i32 + 1) * (knots.len() as i32 - 1)
            }
            None => 33,
        },
        _ => 33,
    };
    n.min(50)
}

/// Grid sampling over `[u1,u2] x [v1,v2]` with `nu x nv` evaluations
/// (`GeomBndLib_OtherSurface.cxx:41-71`, `GeomGridEval_Surface` row-major
/// grid then `Bnd_Box::Add`).
fn grid_box(s: &dyn Surface, u1: f64, u2: f64, v1: f64, v2: f64, nu: i32, nv: i32, tol: f64) -> BndBox {
    let mut b = BndBox::new();
    if nu < 1 || nv < 1 {
        return b;
    }
    for i in 1..=nu {
        let u = u1 + (u2 - u1) * (i - 1) as f64 / (nu - 1) as f64;
        for j in 1..=nv {
            let v = v1 + (v2 - v1) * (j - 1) as f64 / (nv - 1) as f64;
            b.add_point(&s.d0(u, v));
        }
    }
    b.enlarge(tol);
    b
}

/// `GeomBndLib_OtherSurface::Box` (`OtherSurface.cxx:34-71`).
pub fn box_other_surface(s: &dyn Surface, u1: f64, u2: f64, v1: f64, v2: f64, tol: f64) -> BndBox {
    let nu = nb_u_samples(s, u1, u2);
    let nv = nb_v_samples(s, v1, v2);
    grid_box(s, u1, u2, v1, v2, nu, nv, tol)
}

// ---------------------------------------------------------------------------
// Spline surfaces
// ---------------------------------------------------------------------------

/// `ComputePolesIndexes` (`GeomBndLib_BSplineSurface.cxx:30-57`).
fn compute_poles_indexes(
    knots: &[f64],
    mults: &[i32],
    degree: i32,
    min: f64,
    max: f64,
    max_pole_idx: i32,
    is_periodic: bool,
) -> (i32, i32) {
    compute_pole_index_range(knots, mults, degree, min, max, max_pole_idx, is_periodic)
}

/// `GeomBndLib_BSplineSurface::Box` (`BSplineSurface.cxx:70-182`).
pub fn box_bspline_surface(s: &dyn Surface, u1: f64, u2: f64, v1: f64, v2: f64, tol: f64) -> BndBox {
    let Some(b) = spline_of(s) else {
        return box_other_surface(s, u1, u2, v1, v2, tol);
    };
    let (u_min_param, u_max_param) = s.u_range();
    let (v_min_param, v_max_param) = s.v_range();

    let is_out_of_bounds = (u1 - u_min_param) < -PARAMETRIC_CONFUSION
        || (v1 - v_min_param) < -PARAMETRIC_CONFUSION
        || (u2 - u_max_param) > PARAMETRIC_CONFUSION
        || (v2 - v_max_param) > PARAMETRIC_CONFUSION;

    if !is_out_of_bounds {
        let nb_u_poles = b.nb_poles_u() as i32;
        let nb_v_poles = b.nb_poles_v() as i32;
        let is_u_periodic = b.u_periodic;
        let is_v_periodic = b.v_periodic;

        let mut u_min_idx = 1;
        let mut u_max_idx = nb_u_poles;
        let mut v_min_idx = 1;
        let mut v_max_idx = nb_v_poles;

        if u1 > u_min_param || u2 < u_max_param {
            let (knots, mults) = b.distinct_knots_and_mults_u();
            let (mn, mx) = compute_poles_indexes(
                &knots,
                &mults,
                b.deg_u as i32,
                u1,
                u2,
                nb_u_poles,
                is_u_periodic,
            );
            u_min_idx = mn;
            u_max_idx = mx;
        }
        if v1 > v_min_param || v2 < v_max_param {
            let (knots, mults) = b.distinct_knots_and_mults_v();
            let (mn, mx) = compute_poles_indexes(
                &knots,
                &mults,
                b.deg_v as i32,
                v1,
                v2,
                nb_v_poles,
                is_v_periodic,
            );
            v_min_idx = mn;
            v_max_idx = mx;
        }

        let mut a_box = BndBox::new();
        for i in u_min_idx..=u_max_idx {
            let mut ip = i;
            if is_u_periodic && ip > nb_u_poles {
                ip -= nb_u_poles;
            }
            for j in v_min_idx..=v_max_idx {
                let mut jp = j;
                if is_v_periodic && jp > nb_v_poles {
                    jp -= nb_v_poles;
                }
                if ip >= 1
                    && ip <= nb_u_poles
                    && jp >= 1
                    && jp <= nb_v_poles
                {
                    a_box.add_point(&b.poles[(ip - 1) as usize][(jp - 1) as usize]);
                }
            }
        }
        a_box.enlarge(tol);
        return a_box;
    }

    let nu = nb_u_samples_full(s);
    let nv = nb_v_samples_full(s);
    grid_box(s, u1, u2, v1, v2, nu, nv, tol)
}

/// `GeomBndLib_BezierSurface::Box` (`BezierSurface.cxx:43-91`).
pub fn box_bezier_surface(s: &dyn Surface, u1: f64, u2: f64, v1: f64, v2: f64, tol: f64) -> BndBox {
    let (a_u1, a_u2) = s.u_range();
    let (a_v1, a_v2) = s.v_range();
    if (u1 - a_u1).abs() <= PARAMETRIC_CONFUSION
        && (v1 - a_v1).abs() <= PARAMETRIC_CONFUSION
        && (u2 - a_u2).abs() <= PARAMETRIC_CONFUSION
        && (v2 - a_v2).abs() <= PARAMETRIC_CONFUSION
    {
        // Full surface: convex hull of every pole (`BezierSurface.cxx:26-40`).
        let mut a_box = BndBox::new();
        if let Some(b) = spline_of(s) {
            for row in &b.poles {
                for p in row {
                    a_box.add_point(p);
                }
            }
        }
        a_box.enlarge(tol);
        return a_box;
    }
    let nu = nb_u_samples_full(s);
    let nv = nb_v_samples_full(s);
    grid_box(s, u1, u2, v1, v2, nu, nv, tol)
}

// ---------------------------------------------------------------------------
// Surface of revolution (`GeomBndLib_SurfaceOfRevolution.cxx`)
// ---------------------------------------------------------------------------

/// `addRevolutionCircle` (`SurfaceOfRevolution.cxx:35-62`).
fn add_revolution_circle(
    origin: &GpPnt,
    axis_dir: &GpDir,
    basis_pt: &GpPnt,
    u1: f64,
    u2: f64,
    b: &mut BndBox,
) {
    let delta = basis_pt.xyz().subtracted(origin.xyz());
    let h = delta.dot(axis_dir.xyz());
    let center = origin.xyz().added(&axis_dir.xyz().multiply_scalar(h));
    let radial = basis_pt.xyz().subtracted(&center);
    let radius = radial.modulus();
    if radius < CONFUSION {
        // Point on the axis - the revolution is a single point.
        b.add_point(basis_pt);
        return;
    }
    let Ok(x_dir) = GpDir::from_xyz(&radial) else {
        b.add_point(basis_pt);
        return;
    };
    let Ok(ax2) = GpAx2::new(GpPnt::from_xyz(&center), *axis_dir, x_dir) else {
        b.add_point(basis_pt);
        return;
    };
    let circ = GpCirc::new(ax2, radius);
    b.add_box(&box_circ_range(&circ, u1, u2, 0.0));
}

/// `buildRevolutionBox` (`SurfaceOfRevolution.cxx:66-161`).
fn build_revolution_box(
    curve_box: &BndBox,
    origin: &GpPnt,
    axis_dir: &GpDir,
    u1: f64,
    u2: f64,
    tol: f64,
) -> BndBox {
    let mut a_box = BndBox::new();

    if curve_box.is_void() {
        a_box.enlarge(tol);
        return a_box;
    }

    if curve_box.is_open() {
        if curve_box.has_finite_part() {
            if let Some((x_min, x_max, y_min, y_max, z_min, z_max)) = curve_box.finite_part().get() {
                for x in [x_min, x_max] {
                    for y in [y_min, y_max] {
                        for z in [z_min, z_max] {
                            add_revolution_circle(origin, axis_dir, &GpPnt::new(x, y, z), u1, u2, &mut a_box);
                        }
                    }
                }
            }
        }
        if curve_box.is_open_xmin() {
            a_box.open_xmin();
        }
        if curve_box.is_open_xmax() {
            a_box.open_xmax();
        }
        if curve_box.is_open_ymin() {
            a_box.open_ymin();
        }
        if curve_box.is_open_ymax() {
            a_box.open_ymax();
        }
        if curve_box.is_open_zmin() {
            a_box.open_zmin();
        }
        if curve_box.is_open_zmax() {
            a_box.open_zmax();
        }
        a_box.enlarge(tol);
        return a_box;
    }

    if let Some((x_min, x_max, y_min, y_max, z_min, z_max)) = curve_box.get() {
        for x in [x_min, x_max] {
            for y in [y_min, y_max] {
                for z in [z_min, z_max] {
                    add_revolution_circle(origin, axis_dir, &GpPnt::new(x, y, z), u1, u2, &mut a_box);
                }
            }
        }
    }
    a_box.enlarge(tol);
    a_box
}

/// `GeomBndLib_SurfaceOfRevolution::Box` (`SurfaceOfRevolution.cxx:176-233`).
pub fn box_revolution(s: &dyn Surface, u1: f64, u2: f64, v1: f64, v2: f64, tol: f64) -> BndBox {
    let (Some(basis_curve), Some(axis)) = (s.revolution_basis_curve(), s.revolution_axis()) else {
        return box_other_surface(s, u1, u2, v1, v2, tol);
    };
    let origin = *axis.location();
    let axis_dir = *axis.direction();

    let a_v_first = basis_curve.first_parameter();
    let a_v_last = basis_curve.last_parameter();
    let mut a_v_min = if Precision::is_negative_infinite(v1) { a_v_first } else { v1 };
    let mut a_v_max = if Precision::is_positive_infinite(v2) { a_v_last } else { v2 };
    if !basis_curve.is_periodic() {
        if a_v_min < a_v_first {
            a_v_min = a_v_first;
        } else if a_v_min > a_v_last {
            a_v_min = a_v_last;
        }
        if a_v_max < a_v_first {
            a_v_max = a_v_first;
        } else if a_v_max > a_v_last {
            a_v_max = a_v_last;
        }
        if a_v_min > a_v_max {
            std::mem::swap(&mut a_v_min, &mut a_v_max);
        }
    }

    let a_curve_box = box_curve(basis_curve.as_ref(), a_v_min, a_v_max, 0.0);
    if a_curve_box.is_void() {
        return box_other_surface(s, u1, u2, v1, v2, tol);
    }
    build_revolution_box(&a_curve_box, &origin, &axis_dir, u1, u2, tol)
}

// ---------------------------------------------------------------------------
// Surface of extrusion (`GeomBndLib_SurfaceOfExtrusion.cxx`)
// ---------------------------------------------------------------------------

/// `buildExtrusionBox` (`SurfaceOfExtrusion.cxx:29-145`): the basis curve box
/// swept by `vMin..vMax` along `dir`.
fn build_extrusion_box(
    curve_box: &BndBox,
    dir: &GpDir,
    v_min: f64,
    v_max: f64,
    tol: f64,
) -> BndBox {
    let mut a_box = BndBox::new();

    if curve_box.is_void() {
        a_box.enlarge(tol);
        return a_box;
    }

    let dir_xyz = dir.xyz();
    let is_v_min_inf = Precision::is_negative_infinite(v_min);
    let is_v_max_inf = Precision::is_positive_infinite(v_max);

    if curve_box.is_open() {
        // Propagate curve openness through the extrusion rather than returning
        // whole space (`cxx:49-73`).
        if curve_box.is_open_xmin() {
            a_box.open_xmin();
        }
        if curve_box.is_open_xmax() {
            a_box.open_xmax();
        }
        if curve_box.is_open_ymin() {
            a_box.open_ymin();
        }
        if curve_box.is_open_ymax() {
            a_box.open_ymax();
        }
        if curve_box.is_open_zmin() {
            a_box.open_zmin();
        }
        if curve_box.is_open_zmax() {
            a_box.open_zmax();
        }
        if is_v_min_inf && is_v_max_inf {
            open_min_max(dir, &mut a_box);
        } else if is_v_min_inf {
            open_min(dir, &mut a_box);
        } else if is_v_max_inf {
            open_max(dir, &mut a_box);
        }
        if curve_box.has_finite_part() {
            if let Some((x_min, x_max, y_min, y_max, z_min, z_max)) = curve_box.finite_part().get() {
                if !is_v_min_inf {
                    let shift = dir_xyz.multiplied(v_min);
                    a_box.add_point(&GpPnt::new(
                        x_min + shift.x(),
                        y_min + shift.y(),
                        z_min + shift.z(),
                    ));
                    a_box.add_point(&GpPnt::new(
                        x_max + shift.x(),
                        y_max + shift.y(),
                        z_max + shift.z(),
                    ));
                }
                if !is_v_max_inf {
                    let shift = dir_xyz.multiplied(v_max);
                    a_box.add_point(&GpPnt::new(
                        x_min + shift.x(),
                        y_min + shift.y(),
                        z_min + shift.z(),
                    ));
                    a_box.add_point(&GpPnt::new(
                        x_max + shift.x(),
                        y_max + shift.y(),
                        z_max + shift.z(),
                    ));
                }
            }
        }
        a_box.enlarge(tol);
        return a_box;
    }

    if is_v_min_inf && is_v_max_inf {
        a_box.add_box(curve_box);
        open_min_max(dir, &mut a_box);
        a_box.enlarge(tol);
        return a_box;
    }

    // The box is not void and not open here (both arms returned above), so
    // `Get()` (`cxx:115`) has a finite part.
    if let Some((x_min, x_max, y_min, y_max, z_min, z_max)) = curve_box.get() {
        if is_v_min_inf {
            open_min(dir, &mut a_box);
            let shift = dir_xyz.multiplied(v_max);
            a_box.add_point(&GpPnt::new(
                x_min + shift.x(),
                y_min + shift.y(),
                z_min + shift.z(),
            ));
            a_box.add_point(&GpPnt::new(
                x_max + shift.x(),
                y_max + shift.y(),
                z_max + shift.z(),
            ));
        } else if is_v_max_inf {
            open_max(dir, &mut a_box);
            let shift = dir_xyz.multiplied(v_min);
            a_box.add_point(&GpPnt::new(
                x_min + shift.x(),
                y_min + shift.y(),
                z_min + shift.z(),
            ));
            a_box.add_point(&GpPnt::new(
                x_max + shift.x(),
                y_max + shift.y(),
                z_max + shift.z(),
            ));
        } else {
            let shift_min = dir_xyz.multiplied(v_min);
            a_box.add_point(&GpPnt::new(
                x_min + shift_min.x(),
                y_min + shift_min.y(),
                z_min + shift_min.z(),
            ));
            a_box.add_point(&GpPnt::new(
                x_max + shift_min.x(),
                y_max + shift_min.y(),
                z_max + shift_min.z(),
            ));

            let shift_max = dir_xyz.multiplied(v_max);
            a_box.add_point(&GpPnt::new(
                x_min + shift_max.x(),
                y_min + shift_max.y(),
                z_min + shift_max.z(),
            ));
            a_box.add_point(&GpPnt::new(
                x_max + shift_max.x(),
                y_max + shift_max.y(),
                z_max + shift_max.z(),
            ));
        }
    }

    a_box.enlarge(tol);
    a_box
}

/// `GeomBndLib_SurfaceOfExtrusion::Box` (`SurfaceOfExtrusion.cxx:159-219`):
/// `P(U, V) = BasisCurve(U) + V * Direction`.
pub fn box_extrusion(s: &dyn Surface, u1: f64, u2: f64, v1: f64, v2: f64, tol: f64) -> BndBox {
    let (Some(basis_curve), Some(dir)) = (s.extrusion_basis_curve(), s.extrusion_direction())
    else {
        return box_other_surface(s, u1, u2, v1, v2, tol);
    };

    // `cxx:172-199`: clamp the U window onto the basis curve range unless the
    // basis is periodic.
    let mut a_curve_u1 = u1;
    let mut a_curve_u2 = u2;
    if !basis_curve.is_periodic() {
        let a_first = basis_curve.first_parameter();
        let a_last = basis_curve.last_parameter();
        if a_curve_u1 < a_first {
            a_curve_u1 = a_first;
        } else if a_curve_u1 > a_last {
            a_curve_u1 = a_last;
        }
        if a_curve_u2 < a_first {
            a_curve_u2 = a_first;
        } else if a_curve_u2 > a_last {
            a_curve_u2 = a_last;
        }
        if a_curve_u1 > a_curve_u2 {
            std::mem::swap(&mut a_curve_u1, &mut a_curve_u2);
        }
    }

    let a_curve_box = box_curve(basis_curve.as_ref(), a_curve_u1, a_curve_u2, 0.0);
    if a_curve_box.is_void() {
        return box_other_surface(s, u1, u2, v1, v2, tol);
    }
    build_extrusion_box(&a_curve_box, &dir, v1, v2, tol)
}

// ---------------------------------------------------------------------------
// Offset surface (`GeomBndLib_OffsetSurface.cxx`)
// ---------------------------------------------------------------------------

/// `GeomBndLib_OffsetSurface::Box` (`OffsetSurface.cxx:33-58`).
pub fn box_offset_surface(s: &dyn Surface, u1: f64, u2: f64, v1: f64, v2: f64, tol: f64) -> BndBox {
    // PARK: the analytic-equivalent arm (`OffsetSurface.cxx:38-46`) needs
    // `Geom_OffsetSurface::Surface()`, which `GeomOffsetSurface` does not
    // expose, so the conservative arm below is always taken.
    let (Some(basis), Some(offset)) = (s.offset_basis_surface(), s.offset_distance()) else {
        return box_other_surface(s, u1, u2, v1, v2, tol);
    };
    let mut a_local = box_surface(basis.as_ref(), u1, u2, v1, v2, 0.0);
    a_local.enlarge(offset.abs());
    a_local.enlarge(tol);
    a_local
}

// ---------------------------------------------------------------------------
// Dispatch
// ---------------------------------------------------------------------------

/// `GeomBndLib_Surface::Box(theUMin, theUMax, theVMin, theVMax, theTol)`
/// (`GeomBndLib_Surface.cxx:344-364`).
pub fn box_surface(s: &dyn Surface, u1: f64, u2: f64, v1: f64, v2: f64, tol: f64) -> BndBox {
    match surface_kind(s) {
        SurfaceKind::Plane => {
            let pln = s.gp_pln().expect("plane");
            let pos = pln.position();
            analytic::box_plane(&pos, u1, u2, v1, v2, tol)
        }
        SurfaceKind::Cylinder => {
            let c = s.gp_cylinder().expect("cylinder");
            let pos = c.position();
            analytic::box_cylinder(&pos, c.radius(), u1, u2, v1, v2, tol)
        }
        SurfaceKind::Cone => {
            let c = s.gp_cone().expect("cone");
            let pos = c.position();
            analytic::box_cone(&pos, c.radius(), c.semi_angle(), u1, u2, v1, v2, tol)
        }
        SurfaceKind::Sphere => {
            let sp = s.gp_sphere().expect("sphere");
            let pos = *sp.position();
            let loc = sp.location();
            analytic::box_sphere(&pos, &loc, sp.radius(), u1, u2, v1, v2, tol)
        }
        SurfaceKind::Torus => {
            let t = s.gp_torus().expect("torus");
            let pos = *t.position();
            let loc = t.location();
            analytic::box_torus(&pos, &loc, t.major_radius(), t.minor_radius(), u1, u2, v1, v2, tol)
        }
        SurfaceKind::BezierSurface => box_bezier_surface(s, u1, u2, v1, v2, tol),
        SurfaceKind::BSplineSurface => box_bspline_surface(s, u1, u2, v1, v2, tol),
        SurfaceKind::SurfaceOfRevolution => box_revolution(s, u1, u2, v1, v2, tol),
        SurfaceKind::SurfaceOfExtrusion => box_extrusion(s, u1, u2, v1, v2, tol),
        SurfaceKind::OffsetSurface => box_offset_surface(s, u1, u2, v1, v2, tol),
        SurfaceKind::OtherSurface => box_other_surface(s, u1, u2, v1, v2, tol),
    }
}

/// `GeomBndLib_Surface::Box(theUMin, ...)` preceded by `unwrapSurface`
/// (`GeomBndLib_Surface.cxx:114-131`): the `Geom_Surface` handle constructor
/// peels rectangular trims before dispatching, but keeps the *passed* range.
pub fn box_surface_handle(s: &Arc<dyn Surface>, u1: f64, u2: f64, v1: f64, v2: f64, tol: f64) -> BndBox {
    let basis = unwrap_surface(s);
    box_surface(basis.as_ref(), u1, u2, v1, v2, tol)
}

/// `BndLib_AddSurface::Add(S, UMin, UMax, VMin, VMax, Tol, B)`
/// (`BndLib_AddSurface.cxx:35-45`).
pub fn add_surface(
    s: &Arc<dyn Surface>,
    u1: f64,
    u2: f64,
    v1: f64,
    v2: f64,
    tol: f64,
    b: &mut BndBox,
) {
    b.add_box(&box_surface_handle(s, u1, u2, v1, v2, tol));
}

/// `BndLib_AddSurface::Add(S, Tol, B)` with the surface's own bounds
/// (`BndLib_AddSurface.cxx:24-32`). No caller inside `BRepBndLib::Add` (every
/// face passes its `BRepTools::UVBounds`); kept for other `BndLib_AddSurface`
/// ports.
#[allow(dead_code)]
pub fn add_surface_with_bounds(s: &Arc<dyn Surface>, tol: f64, b: &mut BndBox) {
    let (u1, u2) = s.u_range();
    let (v1, v2) = s.v_range();
    add_surface(s, u1, u2, v1, v2, tol, b);
}
