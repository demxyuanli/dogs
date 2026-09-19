//! Ports of the analytic surface boxes called by `GeomBndLib_Surface.cxx`:
//! `GeomBndLib_Plane.hxx`, `GeomBndLib_Cylinder.cxx`, `GeomBndLib_Cone.cxx`,
//! `GeomBndLib_Sphere.cxx`, and the torus arm of `BndLib.cxx` (used by
//! `GeomBndLib_Torus.cxx`).
//!
//! All entry points take the face's UV bounds, mirroring
//! `BndLib_AddSurface` / `GeomBndLib_Surface::Add(theU1,theU2,theV1,theV2,theTol)`.

use occt_core::bnd::BndBox;
use occt_core::elib::slib::{
    cone_v_iso, cylinder_v_iso, plane_value, sphere_parameters, sphere_u_iso, sphere_v_iso, torus_parameters,
    torus_u_iso, torus_v_iso,
};
use occt_core::gp::{GpAx3, GpDir, GpPnt};
use occt_core::precision::{Precision, ANGULAR, CONFUSION, PCONFUSION};

use crate::geom_bnd_lib_circle3d::box_circ_range;
use crate::geom_bnd_lib_elclib2d::in_period;
use crate::geom_bnd_lib_inf3d::{open_max, open_min, open_min_max};

const REAL_EPSILON: f64 = 2.220446049250313e-16;
const COS_PI4: f64 = 0.70710678118654746;
const COS_PI8: f64 = 0.92387953251128674;

// ---------------------------------------------------------------------------
// shared helpers
// ---------------------------------------------------------------------------

/// `ElCLib::InPeriod` (`ElCLib.cxx`) restricted to `[theU1, theU1+2*PI]`.
fn in_period_2pi(the_u: f64, the_u1: f64) -> f64 {
    in_period(the_u, the_u1, the_u1 + 2.0 * std::f64::consts::PI)
}

/// `ElSLib`-style point from two direction carriers (OCCT's templated
/// `Compute` builds `PointType(theO.Coord() + Ra*Cn*theXd.Coord() + ...)`).
fn combo(o: &GpPnt, xd: &GpPnt, yd: &GpPnt, a: f64, b: f64) -> GpPnt {
    GpPnt::new(
        o.x() + a * xd.x() + b * yd.x(),
        o.y() + a * xd.y() + b * yd.y(),
        o.z() + a * xd.z() + b * yd.z(),
    )
}

/// `BndLib.cxx:48-142` templated `Compute` for 3D.
fn compute_range(p1: f64, p2: f64, ra: f64, rb: f64, xd: &GpPnt, yd: &GpPnt, o: &GpPnt, b: &mut BndBox) {
    let (mut a_teta1, mut a_teta2) = if p2 < p1 { (p2, p1) } else { (p1, p2) };

    let a_delta = (a_teta2 - a_teta1).abs();
    if a_delta > 2.0 * std::f64::consts::PI {
        a_teta1 = 0.0;
        a_teta2 = 2.0 * std::f64::consts::PI;
    } else {
        a_teta1 %= 2.0 * std::f64::consts::PI;
        if a_teta1 < 0.0 {
            a_teta1 += 2.0 * std::f64::consts::PI;
        }
        a_teta2 = a_teta1 + a_delta;
    }

    let a_cn1 = a_teta1.cos();
    let a_sn1 = a_teta1.sin();
    let a_cn2 = a_teta2.cos();
    let a_sn2 = a_teta2.sin();
    b.add_point(&combo(o, xd, yd, ra * a_cn1, rb * a_sn1));
    b.add_point(&combo(o, xd, yd, ra * a_cn2, rb * a_sn2));

    let (a_ram, a_rbm) = if a_delta > std::f64::consts::PI / 8.0 {
        (ra / COS_PI8, rb / COS_PI8)
    } else {
        let a_tc = (a_delta / 2.0).cos();
        (ra / a_tc, rb / a_tc)
    };
    b.add_point(&combo(o, xd, yd, a_ram * a_cn1, a_rbm * a_sn1));
    b.add_point(&combo(o, xd, yd, a_ram * a_cn2, a_rbm * a_sn2));

    const A_X_MULT: [f64; 8] = [1.0, COS_PI4, 0.0, -COS_PI4, -1.0, -COS_PI4, 0.0, COS_PI4];
    const A_Y_MULT: [f64; 8] = [0.0, COS_PI4, 1.0, COS_PI4, 0.0, -COS_PI4, -1.0, -COS_PI4];

    let a_deb = (a_teta1 / (std::f64::consts::PI / 4.0)) as i32 + 1;
    let a_fin = (a_teta2 / (std::f64::consts::PI / 4.0)) as i32;
    if a_deb > a_fin {
        return;
    }
    for i in a_deb..=a_fin {
        let idx = (i % 8) as usize;
        b.add_point(&combo(o, xd, yd, a_ram * A_X_MULT[idx], a_rbm * A_Y_MULT[idx]));
    }
}

// ---------------------------------------------------------------------------
// Plane - `GeomBndLib_Plane.hxx:53-168`
// ---------------------------------------------------------------------------

/// `GeomBndLib_Plane::BaryCenter` (`Plane.hxx:91-122`).
fn plane_bary_center(pos: &GpAx3, umin: f64, umax: f64, vmin: f64, vmax: f64) -> GpPnt {
    let is_u1_inf = !umin.is_finite();
    let is_u2_inf = !umax.is_finite();
    let is_v1_inf = !vmin.is_finite();
    let is_v2_inf = !vmax.is_finite();

    let a_u = if is_u1_inf && is_u2_inf {
        0.0
    } else if is_u1_inf {
        umax - 10.0
    } else if is_u2_inf {
        umin + 10.0
    } else {
        (umin + umax) / 2.0
    };
    let a_v = if is_v1_inf && is_v2_inf {
        0.0
    } else if is_v1_inf {
        vmax - 10.0
    } else if is_v2_inf {
        vmin + 10.0
    } else {
        (vmin + vmax) / 2.0
    };
    let pln = occt_core::gp::GpPln::new(*pos);
    plane_value(&pln, a_u, a_v)
}

/// `GeomBndLib_Plane::TreatInfinitePlane` (`Plane.hxx:125-168`).
fn treat_infinite_plane(pos: &GpAx3, umin: f64, umax: f64, vmin: f64, vmax: f64, tol: f64, b: &mut BndBox) {
    let a_norm = *pos.axis().direction();
    let a_location = plane_bary_center(pos, umin, umax, vmin, vmax);

    let dx = GpDir::new(1.0, 0.0, 0.0).unwrap();
    let dy = GpDir::new(0.0, 1.0, 0.0).unwrap();
    let dz = GpDir::new(0.0, 0.0, 1.0).unwrap();
    if a_norm.is_parallel_tol(&dx, REAL_EPSILON) {
        b.add_point(&a_location);
        b.open_ymin();
        b.open_ymax();
        b.open_zmin();
        b.open_zmax();
    } else if a_norm.is_parallel_tol(&dy, REAL_EPSILON) {
        b.add_point(&a_location);
        b.open_xmin();
        b.open_xmax();
        b.open_zmin();
        b.open_zmax();
    } else if a_norm.is_parallel_tol(&dz, REAL_EPSILON) {
        b.add_point(&a_location);
        b.open_xmin();
        b.open_xmax();
        b.open_ymin();
        b.open_ymax();
    } else {
        b.set_whole();
        return;
    }
    b.enlarge(tol);
}

/// `GeomBndLib_Plane::Box(theUMin,theUMax,theVMin,theVMax,theTol)`
/// (`GeomBndLib_Plane.hxx:53-74`).
pub fn box_plane(pos: &GpAx3, umin: f64, umax: f64, vmin: f64, vmax: f64, tol: f64) -> BndBox {
    let mut a_box = BndBox::new();
    if !(umin.is_finite() && umax.is_finite() && vmin.is_finite() && vmax.is_finite()) {
        treat_infinite_plane(pos, umin, umax, vmin, vmax, tol, &mut a_box);
        return a_box;
    }
    let pln = occt_core::gp::GpPln::new(*pos);
    a_box.add_point(&plane_value(&pln, umin, vmin));
    a_box.add_point(&plane_value(&pln, umin, vmax));
    a_box.add_point(&plane_value(&pln, umax, vmin));
    a_box.add_point(&plane_value(&pln, umax, vmax));
    a_box.enlarge(tol);
    a_box
}

// ---------------------------------------------------------------------------
// Cylinder - `GeomBndLib_Cylinder.cxx:30-122`
// ---------------------------------------------------------------------------

fn compute_cylinder(pos: &GpAx3, radius: f64, umin: f64, umax: f64, vmin: f64, vmax: f64, b: &mut BndBox) {
    let a_c = cylinder_v_iso(pos, radius, vmin);
    b.add_box(&box_circ_range(&a_c, umin, umax, 0.0));
    let a_t = vmax - vmin;
    let dir = *pos.axis().direction();
    let mut a_c2 = cylinder_v_iso(pos, radius, vmax);
    let a_loc2 = GpPnt::new(
        a_c.location().x() + a_t * dir.x(),
        a_c.location().y() + a_t * dir.y(),
        a_c.location().z() + a_t * dir.z(),
    );
    a_c2.set_location(&a_loc2);
    b.add_box(&box_circ_range(&a_c2, umin, umax, 0.0));
}

/// `GeomBndLib_Cylinder::Box(theUMin,theUMax,theVMin,theVMax,theTol)`
/// (`GeomBndLib_Cylinder.cxx:53-122`).
pub fn box_cylinder(
    pos: &GpAx3,
    radius: f64,
    umin: f64,
    umax: f64,
    vmin: f64,
    vmax: f64,
    tol: f64,
) -> BndBox {
    let mut a_box = BndBox::new();
    let a_dir = *pos.axis().direction();
    if Precision::is_negative_infinite(vmin) {
        if Precision::is_positive_infinite(vmax) {
            let a_c = cylinder_v_iso(pos, radius, 0.0);
            a_box.add_box(&box_circ_range(&a_c, umin, umax, 0.0));
            open_min_max(&a_dir, &mut a_box);
        } else {
            compute_cylinder(pos, radius, umin, umax, 0.0, vmax, &mut a_box);
            open_min(&a_dir, &mut a_box);
        }
    } else if Precision::is_positive_infinite(vmin) {
        if Precision::is_positive_infinite(vmax) {
            // OCCT throws here for (vmin=+inf, vmax=+inf); the box stays empty.
        } else {
            compute_cylinder(pos, radius, umin, umax, 0.0, vmax, &mut a_box);
            open_max(&a_dir, &mut a_box);
        }
    } else if Precision::is_negative_infinite(vmax) {
        compute_cylinder(pos, radius, umin, umax, vmin, 0.0, &mut a_box);
        open_min(&a_dir, &mut a_box);
    } else if Precision::is_positive_infinite(vmax) {
        compute_cylinder(pos, radius, umin, umax, vmin, 0.0, &mut a_box);
        open_max(&a_dir, &mut a_box);
    } else {
        compute_cylinder(pos, radius, umin, umax, vmin, vmax, &mut a_box);
    }
    a_box.enlarge(tol);
    a_box
}

// ---------------------------------------------------------------------------
// Cone - `GeomBndLib_Cone.cxx:26-154`
// ---------------------------------------------------------------------------

fn compute_cone(
    pos: &GpAx3,
    radius: f64,
    s_ang: f64,
    umin: f64,
    umax: f64,
    vmin: f64,
    vmax: f64,
    b: &mut BndBox,
) {
    let a_c = cone_v_iso(pos, radius, s_ang, vmin);
    if a_c.radius() > CONFUSION {
        b.add_box(&box_circ_range(&a_c, umin, umax, 0.0));
    } else {
        b.add_point(&a_c.location());
    }

    let a_c2 = cone_v_iso(pos, radius, s_ang, vmax);
    if a_c2.radius() > CONFUSION {
        b.add_box(&box_circ_range(&a_c2, umin, umax, 0.0));
    } else {
        b.add_point(&a_c2.location());
    }
}

/// `GeomBndLib_Cone::Box(theUMin,theUMax,theVMin,theVMax,theTol)`
/// (`GeomBndLib_Cone.cxx:72-154`).
pub fn box_cone(
    pos: &GpAx3,
    radius: f64,
    s_ang: f64,
    umin: f64,
    umax: f64,
    vmin: f64,
    vmax: f64,
    tol: f64,
) -> BndBox {
    let mut a_box = BndBox::new();
    let a_dir = *pos.axis().direction();
    if Precision::is_negative_infinite(vmin) {
        if Precision::is_positive_infinite(vmax) {
            let a_c = cone_v_iso(pos, radius, s_ang, 0.0);
            if a_c.radius() > CONFUSION {
                a_box.add_box(&box_circ_range(&a_c, umin, umax, 0.0));
            } else {
                a_box.add_point(&a_c.location());
            }
            open_min_max(&a_dir, &mut a_box);
        } else {
            compute_cone(pos, radius, s_ang, umin, umax, 0.0, vmax, &mut a_box);
            open_min(&a_dir, &mut a_box);
        }
    } else if Precision::is_positive_infinite(vmin) {
        if Precision::is_positive_infinite(vmax) {
            // OCCT throws here; the box stays empty.
        } else {
            compute_cone(pos, radius, s_ang, umin, umax, 0.0, vmax, &mut a_box);
            open_max(&a_dir, &mut a_box);
        }
    } else if Precision::is_negative_infinite(vmax) {
        compute_cone(pos, radius, s_ang, umin, umax, vmin, 0.0, &mut a_box);
        open_min(&a_dir, &mut a_box);
    } else if Precision::is_positive_infinite(vmax) {
        compute_cone(pos, radius, s_ang, umin, umax, vmin, 0.0, &mut a_box);
        open_max(&a_dir, &mut a_box);
    } else {
        compute_cone(pos, radius, s_ang, umin, umax, vmin, vmax, &mut a_box);
    }
    a_box.enlarge(tol);
    a_box
}

// ---------------------------------------------------------------------------
// Sphere - `GeomBndLib_Sphere.cxx:22-158`
// ---------------------------------------------------------------------------

/// `GeomBndLib_Sphere::Box(theTol)` (`Sphere.cxx:22-32`).
pub fn box_sphere_full(loc: &GpPnt, radius: f64, tol: f64) -> BndBox {
    let mut a_box = BndBox::new();
    a_box.update(
        loc.x() - radius,
        loc.y() - radius,
        loc.z() - radius,
        loc.x() + radius,
        loc.y() + radius,
        loc.z() + radius,
    );
    a_box.enlarge(tol);
    a_box
}

/// `computeSphere` (`GeomBndLib_Sphere.cxx:38-135`).
fn compute_sphere(pos: &GpAx3, loc: &GpPnt, radius: f64, umin: f64, umax: f64, vmin: f64, vmax: f64, b: &mut BndBox) {
    let a_xmin = loc.x() - radius;
    let a_xmax = loc.x() + radius;
    let a_ymin = loc.y() - radius;
    let a_ymax = loc.y() + radius;
    let a_zmin = loc.z() - radius;
    let a_zmax = loc.z() + radius;

    let an_uper = 2.0 * std::f64::consts::PI - PCONFUSION;
    let a_vper = std::f64::consts::PI - PCONFUSION;
    if umax - umin >= an_uper && vmax - vmin >= a_vper {
        b.update(a_xmin, a_ymin, a_zmin, a_xmax, a_ymax, a_zmax);
        return;
    }

    let an_umax = umin + 2.0 * std::f64::consts::PI;

    // Six axis extrema, tested against the patch bounds.
    let extrema: [(usize, f64); 6] = [
        (0, a_xmin),
        (0, a_xmax),
        (1, a_ymin),
        (1, a_ymax),
        (2, a_zmin),
        (2, a_zmax),
    ];
    for (idx, val) in extrema {
        let mut a_pext = *loc;
        match idx {
            0 => a_pext.set_x(val),
            1 => a_pext.set_y(val),
            _ => a_pext.set_z(val),
        }
        let (mut an_u, a_v) = sphere_parameters(pos, &a_pext);
        an_u = in_period_2pi(an_u, umin);
        // Reset the modified coordinate is unnecessary: the loop re-derives
        // a_pext from `loc` each iteration.
        if an_u >= umin && an_u <= umax && a_v >= vmin && a_v <= vmax {
            b.add_point(&a_pext);
        }
        let _ = an_umax;
    }

    // Boundary iso-curves.
    let a_c = sphere_u_iso(pos, radius, umin);
    b.add_box(&box_circ_range(&a_c, vmin, vmax, 0.0));
    let a_c = sphere_u_iso(pos, radius, umax);
    b.add_box(&box_circ_range(&a_c, vmin, vmax, 0.0));
    let a_c = sphere_v_iso(pos, radius, vmin);
    b.add_box(&box_circ_range(&a_c, umin, umax, 0.0));
    let a_c = sphere_v_iso(pos, radius, vmax);
    b.add_box(&box_circ_range(&a_c, umin, umax, 0.0));
}

/// `GeomBndLib_Sphere::Box(theUMin,theUMax,theVMin,theVMax,theTol)`
/// (`GeomBndLib_Sphere.cxx:137-158`).
pub fn box_sphere(
    pos: &GpAx3,
    loc: &GpPnt,
    radius: f64,
    umin: f64,
    umax: f64,
    vmin: f64,
    vmax: f64,
    tol: f64,
) -> BndBox {
    if umin.abs() < ANGULAR
        && (umax - 2.0 * std::f64::consts::PI).abs() < ANGULAR
        && (vmin + std::f64::consts::PI / 2.0).abs() < ANGULAR
        && (vmax - std::f64::consts::PI / 2.0).abs() < ANGULAR
    {
        return box_sphere_full(loc, radius, tol);
    }
    let mut a_box = BndBox::new();
    compute_sphere(pos, loc, radius, umin, umax, vmin, vmax, &mut a_box);
    a_box.enlarge(tol);
    a_box
}

// ---------------------------------------------------------------------------
// Torus - `BndLib.cxx:1414-1610` (`BndLib::Add(const gp_Torus&, ...)`)
// ---------------------------------------------------------------------------

/// `computeDegeneratedTorus` (`BndLib.cxx:1418-1511`).
fn compute_degenerated_torus(
    pos: &GpAx3,
    loc: &GpPnt,
    major: f64,
    minor: f64,
    umin: f64,
    umax: f64,
    vmin: f64,
    vmax: f64,
    b: &mut BndBox,
) {
    let a_xmin = loc.x() - major - minor;
    let a_xmax = loc.x() + major + minor;
    let a_ymin = loc.y() - major - minor;
    let a_ymax = loc.y() + major + minor;
    let a_zmin = loc.z() - minor;
    let a_zmax = loc.z() + minor;

    let a_phi = (-major / minor).acos();

    let an_uper = 2.0 * std::f64::consts::PI - PCONFUSION;
    let a_vper = 2.0 * a_phi - PCONFUSION;
    if umax - umin >= an_uper && vmax - vmin >= a_vper {
        b.update(a_xmin, a_ymin, a_zmin, a_xmax, a_ymax, a_zmax);
        return;
    }

    let an_umax = umin + 2.0 * std::f64::consts::PI;
    let extrema: [(usize, f64); 6] = [
        (0, a_xmin),
        (0, a_xmax),
        (1, a_ymin),
        (1, a_ymax),
        (2, a_zmin),
        (2, a_zmax),
    ];
    for (idx, val) in extrema {
        let mut a_pext = *loc;
        match idx {
            0 => a_pext.set_x(val),
            1 => a_pext.set_y(val),
            _ => a_pext.set_z(val),
        }
        let (mut an_u, a_v) = torus_parameters(pos, major, minor, &a_pext);
        an_u = in_period_2pi(an_u, umin);
        if an_u >= umin && an_u <= umax && a_v >= vmin && a_v <= vmax {
            b.add_point(&a_pext);
        }
        let _ = an_umax;
    }

    let a_c = torus_u_iso(pos, major, minor, umin);
    b.add_box(&box_circ_range(&a_c, vmin, vmax, 0.0));
    let a_c = torus_u_iso(pos, major, minor, umax);
    b.add_box(&box_circ_range(&a_c, vmin, vmax, 0.0));
    let a_c = torus_v_iso(pos, major, minor, vmin);
    b.add_box(&box_circ_range(&a_c, umin, umax, 0.0));
    let a_c = torus_v_iso(pos, major, minor, vmax);
    b.add_box(&box_circ_range(&a_c, umin, umax, 0.0));
}

/// `BndLib::Add(const gp_Torus&, UMin, UMax, VMin, VMax, Tol, B)`
/// (`BndLib.cxx:1541-1610`).
pub fn box_torus(
    pos: &GpAx3,
    loc: &GpPnt,
    major: f64,
    minor: f64,
    umin: f64,
    umax: f64,
    vmin: f64,
    vmax: f64,
    tol: f64,
) -> BndBox {
    let (fi1, fi2) = if vmax < vmin {
        (
            (vmax / (std::f64::consts::PI / 4.0)).floor() as i32,
            (vmin / (std::f64::consts::PI / 4.0)).floor() as i32,
        )
    } else {
        (
            (vmin / (std::f64::consts::PI / 4.0)).floor() as i32,
            (vmax / (std::f64::consts::PI / 4.0)).floor() as i32,
        )
    };
    let fi2 = fi2 + 1;

    if fi2 < fi1 {
        return BndBox::new();
    }

    let mut b = BndBox::new();
    if major < minor {
        compute_degenerated_torus(pos, loc, major, minor, umin, umax, vmin, vmax, &mut b);
        b.enlarge(tol);
        return b;
    }

    let a_z_dir = pos.axis().direction().xyz();
    let a_loc_xyz = loc.xyz();
    let a_xd = GpPnt::new(pos.x_direction().x(), pos.x_direction().y(), pos.x_direction().z());
    let a_yd = GpPnt::new(pos.y_direction().x(), pos.y_direction().y(), pos.y_direction().z());

    let a_radius_mult: [f64; 8] = [1.0, COS_PI4, 0.0, -COS_PI4, -1.0, -COS_PI4, 0.0, COS_PI4];
    let a_z_mult: [f64; 8] = [0.0, COS_PI4, 1.0, COS_PI4, 0.0, -COS_PI4, -1.0, -COS_PI4];

    for i in fi1..=fi2 {
        let idx = (((i % 8) + 8) % 8) as usize;
        let a_radius = major + minor * a_radius_mult[idx];
        let center = GpPnt::from_xyz(&a_loc_xyz.add(&a_z_dir.multiply_scalar(minor * a_z_mult[idx])));
        compute_range(umin, umax, a_radius, a_radius, &a_xd, &a_yd, &center, &mut b);
    }

    b.enlarge(tol);
    b
}
