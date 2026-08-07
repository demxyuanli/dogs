//! Point-surface / curve-surface extrema. Port of `Extrema_ExtPS`,
//! `Extrema_ExtPElS`, `Extrema_ExtPRevS`, `Extrema_ExtPExtS`,
//! `Extrema_ExtCS`, `Extrema_ExtElCS`, `Extrema_FuncExtCS` (TKGeomBase).
//!
//! The analytic point-to-surface solvers (`Extrema_ExtPElS`: plane, sphere,
//! cylinder, cone, torus) are ported exactly from the OCCT `.cxx`. Revolved
//! and extruded surfaces (`ExtPRevS` / `ExtPExtS`) are heavy (they reduce to a
//! point-curve extrema on the generating curve) and are routed through the
//! general Newton path, exactly as OCCT does for the non-analytically-
//! computable cases.
//!
//! `dyn Surface` cannot be downcast, so analytic dispatch classifies by
//! geometric invariants (mirroring `brep_surface::classify_surface`): only
//! planes (constant normal) and spheres (equidistant samples from a solved
//! center) are classified; everything else goes through the Newton path.
//!
//! The Newton systems are solved with a NUMERIC Jacobian: the `Surface` trait
//! exposes only `d0`/`d1` (no second derivatives), so the Jacobian of the
//! orthogonality conditions `F = ((S-P)·Su, (S-P)·Sv)` is computed by central
//! finite differences of `d1` with `eps ≈ 1e-6` relative to the parameter
//! range. // ponytail: numeric Jacobian, Surface trait lacks d2

use std::cmp::Ordering;

use occt_core::elib::{clib, slib};
use occt_core::gp::{GpAx3, GpCone, GpCylinder, GpDir, GpLin, GpPln, GpPnt, GpSphere, GpTorus, GpVec};
use occt_core::precision::CONFUSION;

use crate::curve::Curve;
use crate::extrema::ExtremaPair;
use crate::surface::Surface;

/// Angle snap: OCCT's `ExtPElS_MyEps = Epsilon(2π)`; `ANGULAR` (1e-12) is the
/// crate's angular resolution and is used here for the same purpose.
const SNAP_EPS: f64 = 1e-12;

// ---------------------------------------------------------------------------
// Pair helpers.
// ---------------------------------------------------------------------------

/// Point-surface pair: `p1` is the given point, `p2` the surface point,
/// surface parameters go in `u2`/`v2` (`u1` mirrors `u2`, per the existing
/// `extrema::ExtremaPair` convention for point–surface results).
fn ps_pair(p: &GpPnt, u: f64, v: f64, q: GpPnt) -> ExtremaPair {
    ExtremaPair {
        p1: *p,
        p2: q,
        distance: p.distance(&q),
        u1: u,
        v1: None,
        u2: u,
        v2: Some(v),
    }
}

/// Curve-surface pair: `p1` on the curve (param `t`), `p2` on the surface
/// (params `u`,`v`).
fn cs_pair(pc: GpPnt, t: f64, ps: GpPnt, u: f64, v: f64, d: f64) -> ExtremaPair {
    ExtremaPair {
        p1: pc,
        p2: ps,
        distance: d,
        u1: t,
        v1: None,
        u2: u,
        v2: Some(v),
    }
}

/// `gp_Vec::AngleWithRef` — signed angle from `a` to `b` about `vref`, using
/// the direction-based `GpDir::angle_with_ref` (returns 0 on any degenerate
/// input, which OCCT then snaps to 0 anyway).
fn vec_angle_with_ref(a: &GpVec, b: &GpVec, vref: &GpVec) -> f64 {
    match (GpDir::from_vec(a), GpDir::from_vec(b), GpDir::from_vec(vref)) {
        (Ok(da), Ok(db), Ok(dr)) => da.angle_with_ref(&db, &dr),
        _ => 0.0,
    }
}

/// Clamp an unbounded parameter range to a finite sampling window.
pub(crate) fn surf_bound_u(s: &dyn Surface) -> (f64, f64) {
    let (a, b) = s.u_range();
    if a.is_finite() && b.is_finite() && b > a {
        (a, b)
    } else {
        (-10.0, 10.0)
    }
}

/// Clamp an unbounded parameter range to a finite sampling window.
pub(crate) fn surf_bound_v(s: &dyn Surface) -> (f64, f64) {
    let (a, b) = s.v_range();
    if a.is_finite() && b.is_finite() && b > a {
        (a, b)
    } else {
        (-10.0, 10.0)
    }
}

fn curve_bound(c: &dyn Curve) -> (f64, f64) {
    let (a, b) = (c.first_parameter(), c.last_parameter());
    if a.is_finite() && b.is_finite() && b > a {
        (a, b)
    } else {
        (-10.0, 10.0)
    }
}

/// Finite, sane sampling bounds for a surface (unbounded ranges clamp to ±1).
fn sample_bounds(s: &dyn Surface) -> (f64, f64, f64, f64) {
    let clamp = |a: f64, b: f64| {
        if a.is_finite() && b.is_finite() && b > a {
            (a, b)
        } else {
            (-1.0, 1.0)
        }
    };
    let (u0, u1) = clamp(s.u_range().0, s.u_range().1);
    let (v0, v1) = clamp(s.v_range().0, s.v_range().1);
    (u0, u1, v0, v1)
}

// ---------------------------------------------------------------------------
// Analytic point-to-analytic-surface solvers (port of `Extrema_ExtPElS`).
// ---------------------------------------------------------------------------

/// Closest point on a plane. Port of `Extrema_ExtPElS::Perform(gp_Pln)`:
/// project `p` along the plane normal; the in-plane parameters are the
/// coordinates of the projection in the plane's own frame.
pub fn point_plane_extrema(pl: &GpPln, p: &GpPnt) -> ExtremaPair {
    let pos = pl.position();
    let o = pl.location();
    let x = GpVec::from_xyz(pos.x_direction().xyz());
    let y = GpVec::from_xyz(pos.y_direction().xyz());
    let u = GpVec::from_pnts(&o, p).dot(&x);
    let vv = GpVec::from_pnts(&o, p).dot(&y);
    let q = slib::plane_value(pl, u, vv);
    ExtremaPair {
        p1: *p,
        p2: q,
        distance: p.distance(&q),
        u1: u,
        v1: None,
        u2: u,
        v2: Some(vv),
    }
}

/// Point to sphere: min and max. Port of
/// `Extrema_ExtPElS::Perform(gp_Sphere)`: the closest point is on the ray
/// `O→P`, the farthest on the opposite ray. Degenerate (P on the sphere axis)
/// yields the two poles.
pub fn point_sphere_extrema(sp: &GpSphere, p: &GpPnt) -> Vec<ExtremaPair> {
    let pos = sp.position();
    let o = sp.location();
    let tol = CONFUSION;
    let op = GpVec::from_pnts(&o, p);
    if op.square_magnitude() < tol * tol {
        return Vec::new(); // P == O: infinite solutions
    }
    let oz = GpVec::from_xyz(pos.direction().xyz());
    let zp = op.dot(&oz);
    let pp = p.translated_vec(&oz.multiplied_scalar(-zp));
    let opp = GpVec::from_pnts(&o, &pp);
    let my_z = GpVec::from_xyz(pos.x_direction().xyz())
        .crossed(&GpVec::from_xyz(pos.y_direction().xyz()));
    let (u1, u2, v) = if opp.square_magnitude() < tol * tol {
        let v = if zp < 0.0 {
            -std::f64::consts::FRAC_PI_2
        } else {
            std::f64::consts::FRAC_PI_2
        };
        (0.0, 0.0, v)
    } else {
        let mut u1 = vec_angle_with_ref(&GpVec::from_xyz(pos.x_direction().xyz()), &opp, &my_z);
        if u1.abs() < SNAP_EPS {
            u1 = 0.0;
        }
        let u2 = u1 + std::f64::consts::PI;
        let u1 = if u1 < 0.0 {
            u1 + 2.0 * std::f64::consts::PI
        } else {
            u1
        };
        let mut v = op.angle(&opp);
        if zp < 0.0 {
            v = -v;
        }
        (u1, u2, v)
    };
    let q1 = slib::sphere_value(sp, u1, v);
    let q2 = slib::sphere_value(sp, u2, -v);
    vec![ps_pair(p, u1, v, q1), ps_pair(p, u2, -v, q2)]
}

/// Point to cylinder: min and max. Port of
/// `Extrema_ExtPElS::Perform(gp_Cylinder)`: project onto the plane `XOY`,
/// the two extrema lie on the radial line through the projection.
pub fn point_cylinder_extrema(cy: &GpCylinder, p: &GpPnt) -> Vec<ExtremaPair> {
    let pos = cy.position();
    let o = cy.location();
    let tol = CONFUSION;
    let oz = GpVec::from_xyz(pos.direction().xyz());
    let v = GpVec::from_pnts(&o, p).dot(&oz);
    let pp = p.translated_vec(&oz.multiplied_scalar(-v));
    let opp = GpVec::from_pnts(&o, &pp);
    if opp.magnitude() < tol {
        return Vec::new(); // point on the axis: infinite solutions
    }
    let my_z = GpVec::from_xyz(pos.x_direction().xyz())
        .crossed(&GpVec::from_xyz(pos.y_direction().xyz()));
    let mut u1 = vec_angle_with_ref(&GpVec::from_xyz(pos.x_direction().xyz()), &opp, &my_z);
    if u1.abs() < SNAP_EPS {
        u1 = 0.0;
    }
    let u2 = u1 + std::f64::consts::PI;
    let u1 = if u1 < 0.0 {
        u1 + 2.0 * std::f64::consts::PI
    } else {
        u1
    };
    let q1 = slib::cylinder_value(cy, u1, v);
    let q2 = slib::cylinder_value(cy, u2, v);
    vec![ps_pair(p, u1, v, q1), ps_pair(p, u2, v, q2)]
}

/// Point to cone: min and max. Port of
/// `Extrema_ExtPElS::Perform(gp_Cone)`.
///
/// This crate parameterizes the cone so the VERTEX sits at `location()`
/// (radius-0 point); OCCT's `gp_Cone::Apex()` is that vertex, so the port maps
/// `M = S.location()`, `Vm = -RefRadius/sin(A)` and `DirZ = +OZ` (the cone
/// opens along +Z).
pub fn point_cone_extrema(co: &GpCone, p: &GpPnt) -> Vec<ExtremaPair> {
    let pos = co.position();
    let a = co.semi_angle();
    let tol = CONFUSION;
    let oz = GpVec::from_xyz(pos.direction().xyz());
    let r = co.radius();
    // True cone vertex (radius-0 point): `slib::cone_value` is parameterized
    // from the placement, so the vertex is placement − (RefRadius/tan α)·axis
    // (at surface parameter v = −RefRadius/sin α).
    let o = GpPnt::from_xyz(&pos.location().coord.subtracted(&pos.direction().xyz().multiplied(r / a.tan())));
    let mp = GpVec::from_pnts(&o, p);
    let l2 = mp.square_magnitude();
    let vm = -(r / a.sin());
    if l2 < tol * tol {
        // P coincides with the vertex: a single minimum.
        return vec![ps_pair(p, 0.0, vm, o)];
    }
    let dirz = oz; // semi-angle > 0 → cone opens along +Z
    let zp = GpVec::from_pnts(&o, p).dot(&oz);
    let pp = p.translated_vec(&oz.multiplied_scalar(-zp));
    let opp = GpVec::from_pnts(&o, &pp);
    if opp.square_magnitude() < tol * tol {
        return Vec::new(); // P on the axis: infinite solutions
    }
    let my_z = GpVec::from_xyz(pos.x_direction().xyz())
        .crossed(&GpVec::from_xyz(pos.y_direction().xyz()));
    let same = dirz.dot(&mp) >= 0.0;
    let mut u1 = vec_angle_with_ref(&GpVec::from_xyz(pos.x_direction().xyz()), &opp, &my_z);
    if u1.abs() < SNAP_EPS {
        u1 = 0.0;
    }
    if !same {
        u1 += std::f64::consts::PI;
    }
    let u2 = u1 + std::f64::consts::PI;
    let u1 = if u1 < 0.0 {
        u1 + 2.0 * std::f64::consts::PI
    } else {
        u1
    };
    let u2 = if u2 > 2.0 * std::f64::consts::PI {
        u2 - 2.0 * std::f64::consts::PI
    } else {
        u2
    };
    let b = mp.angle(&dirz);
    let l = l2.sqrt();
    let (v1, v2) = if !same {
        let b = std::f64::consts::PI - b;
        (-l * (b - a).cos(), -l * (b + a).cos())
    } else {
        (l * (b - a).cos(), l * (b + a).cos())
    };
    let sense = oz.dot(&dirz.normalized());
    let v1 = v1 * sense + vm;
    let v2 = v2 * sense + vm;
    let q1 = slib::cone_value(co, u1, v1);
    let q2 = slib::cone_value(co, u2, v2);
    vec![ps_pair(p, u1, v1, q1), ps_pair(p, u2, v2, q2)]
}

/// Point to torus: up to four extrema (near/far on both the outer and inner
/// equators). Port of `Extrema_ExtPElS::Perform(gp_Torus)`.
pub fn point_torus_extrema(to: &GpTorus, p: &GpPnt) -> Vec<ExtremaPair> {
    let pos = to.position();
    let o = to.location();
    let tol = CONFUSION;
    let tol2 = tol * tol;
    let oz = GpVec::from_xyz(pos.direction().xyz());
    let pp = p.translated_vec(&oz.multiplied_scalar(-GpVec::from_pnts(&o, p).dot(&oz)));
    let opp = GpVec::from_pnts(&o, &pp);
    let r2 = opp.square_magnitude();
    if r2 < tol2 {
        return Vec::new(); // P on the torus axis: infinite solutions
    }
    let my_z = GpVec::from_xyz(pos.x_direction().xyz())
        .crossed(&GpVec::from_xyz(pos.y_direction().xyz()));
    let mut u1 = vec_angle_with_ref(&GpVec::from_xyz(pos.x_direction().xyz()), &opp, &my_z);
    if u1.abs() < SNAP_EPS {
        u1 = 0.0;
    }
    let u2 = u1 + std::f64::consts::PI;
    let u1 = if u1 < 0.0 {
        u1 + 2.0 * std::f64::consts::PI
    } else {
        u1
    };
    let rr = r2.sqrt();
    let maj = to.major_radius();
    let oo1 = opp.divided(rr).multiplied_scalar(maj);
    let o1 = o.translated_vec(&oo1);
    let o2 = o.translated_vec(&oo1.multiplied_scalar(-1.0));
    if o1.distance(p) < tol || o2.distance(p) < tol {
        return Vec::new(); // P on a spine circle: infinite solutions
    }
    let mut v1 = vec_angle_with_ref(&opp, &GpVec::from_pnts(&o1, p), &opp.crossed(&oz));
    if v1.abs() < SNAP_EPS {
        v1 = 0.0;
    }
    let opp_rev = opp.reversed();
    let mut v2 = vec_angle_with_ref(&opp_rev, &GpVec::from_pnts(p, &o2), &opp_rev.crossed(&oz));
    if v2.abs() < SNAP_EPS {
        v2 = 0.0;
    }
    let v1 = if v1 < 0.0 {
        v1 + 2.0 * std::f64::consts::PI
    } else {
        v1
    };
    let v2 = if v2 < 0.0 {
        v2 + 2.0 * std::f64::consts::PI
    } else {
        v2
    };
    vec![
        ps_pair(p, u1, v1, slib::torus_value(to, u1, v1)),
        ps_pair(p, u1, v1 + std::f64::consts::PI, slib::torus_value(to, u1, v1 + std::f64::consts::PI)),
        ps_pair(p, u2, v2, slib::torus_value(to, u2, v2)),
        ps_pair(p, u2, v2 + std::f64::consts::PI, slib::torus_value(to, u2, v2 + std::f64::consts::PI)),
    ]
}

// ---------------------------------------------------------------------------
// Parameter recovery (ElSLib::Parameters ports).
// ---------------------------------------------------------------------------

/// (u, v) of a point on a plane (in the plane's own frame).
pub(crate) fn plane_parameters(pl: &GpPln, q: &GpPnt) -> (f64, f64) {
    let pos = pl.position();
    let o = pl.location();
    let v = GpVec::from_pnts(&o, q);
    let u = v.dot(&GpVec::from_xyz(pos.x_direction().xyz()));
    let vv = v.dot(&GpVec::from_xyz(pos.y_direction().xyz()));
    (u, vv)
}

/// (u, v) of a point on a sphere. Port of `ElSLib::Parameters(gp_Sphere)`
/// against this crate's `slib::sphere_value` parameterization.
pub(crate) fn sphere_parameters(sp: &GpSphere, q: &GpPnt) -> (f64, f64) {
    let pos = sp.position();
    let w = GpVec::from_pnts(&sp.location(), q);
    let x = w.dot(&GpVec::from_xyz(pos.x_direction().xyz()));
    let y = w.dot(&GpVec::from_xyz(pos.y_direction().xyz()));
    let z = w.dot(&GpVec::from_xyz(pos.direction().xyz()));
    let mut u = y.atan2(x);
    if u < 0.0 {
        u += 2.0 * std::f64::consts::PI;
    }
    let v = (z / sp.radius()).clamp(-1.0, 1.0).asin();
    (u, v)
}

// ---------------------------------------------------------------------------
// Analytic curve-surface solvers (port of `Extrema_ExtElCS`).
// ---------------------------------------------------------------------------

/// Line vs sphere. Port of `Extrema_ExtElCS::Perform(gp_Lin, gp_Sphere)`:
/// project the sphere center onto the line; intersecting lines contribute the
/// intersection points (distance 0), plus the far point; non-intersecting
/// lines contribute the near and far points on the sphere.
pub fn line_sphere_extrema(l: &GpLin, sp: &GpSphere) -> Vec<ExtremaPair> {
    let c = sp.location();
    let r = sp.radius();
    let d = GpVec::from_xyz(l.direction().xyz());
    let a = l.location();
    let ac = GpVec::from_pnts(&a, &c);
    let t0 = ac.dot(&d);
    let f = a.translated_vec(&d.multiplied_scalar(t0));
    let h = f.distance(&c);
    let p1 = clib::line_value(l, t0); // the perpendicular foot on the line
    let mut out: Vec<ExtremaPair> = Vec::new();
    if h <= r {
        let s = (r * r - h * h).sqrt();
        let ts = if s > 1e-12 { vec![t0 - s, t0 + s] } else { vec![t0] };
        for t in ts {
            let q = clib::line_value(l, t);
            let (u, v) = sphere_parameters(sp, &q);
            out.push(cs_pair(q, t, q, u, v, 0.0));
        }
        if h > 1e-12 {
            let dir = GpVec::from_pnts(&c, &f).divided(h);
            let far = c.translated_vec(&dir.multiplied_scalar(-r));
            let (u, v) = sphere_parameters(sp, &far);
            out.push(cs_pair(p1, t0, far, u, v, h + r));
        }
    } else if h > 1e-12 {
        let dir = GpVec::from_pnts(&c, &f).divided(h);
        let near = c.translated_vec(&dir.multiplied_scalar(r));
        let far = c.translated_vec(&dir.multiplied_scalar(-r));
        for (q, dist) in [(near, h - r), (far, h + r)] {
            let (u, v) = sphere_parameters(sp, &q);
            out.push(cs_pair(p1, t0, q, u, v, dist));
        }
    }
    out
}

/// Line vs plane. Port of `Extrema_ExtElCS::Perform(gp_Lin, gp_Pln)`: a
/// non-parallel line meets the plane (distance 0); a parallel line has a
/// constant distance.
pub fn line_plane_extrema(l: &GpLin, pl: &GpPln) -> Vec<ExtremaPair> {
    let n = GpVec::from_xyz(pl.axis().direction().xyz());
    let d = GpVec::from_xyz(l.direction().xyz());
    let dn = d.dot(&n);
    if dn.abs() > 1e-12 {
        // P(t) = A + t·d meets (X-O)·n = 0 at t = (O-A)·n / (d·n).
        let t = GpVec::from_pnts(&l.location(), &pl.location()).dot(&n) / dn;
        let q = clib::line_value(l, t);
        let (u, v) = plane_parameters(pl, &q);
        vec![cs_pair(q, t, q, u, v, 0.0)]
    } else {
        let d0 = GpVec::from_pnts(&pl.location(), &l.location()).dot(&n).abs();
        let q = clib::line_value(l, 0.0);
        let (u, v) = plane_parameters(pl, &q);
        vec![cs_pair(q, 0.0, q, u, v, d0)]
    }
}

// ---------------------------------------------------------------------------
// Surface classification by geometric invariants (brep_surface pattern,
// inlined here — occt-geom does not depend on occt-topo).
// ---------------------------------------------------------------------------

/// Unit normal of a surface at (u, v), robust to surfaces whose `d1` returns
/// zero vectors (some ported surfaces only implement `d0`).
fn surface_normal(s: &dyn Surface, u: f64, v: f64) -> GpVec {
    let (_, du, dv) = s.d1(u, v);
    let n = du.xyz().crossed(dv.xyz());
    if n.square_modulus() > 1e-30 {
        let m = n.modulus();
        return GpVec::new(n.x / m, n.y / m, n.z / m);
    }
    let eps = 1e-6;
    let (u0, u1) = s.u_range();
    let (v0, v1) = s.v_range();
    let hu = if u1 > u0 { (u1 - u0) * 1e-4 } else { eps };
    let hv = if v1 > v0 { (v1 - v0) * 1e-4 } else { eps };
    let p0 = s.d0(u, v);
    let duv = GpVec::from_pnts(&p0, &s.d0(u + hu, v));
    let dvv = GpVec::from_pnts(&p0, &s.d0(u, v + hv));
    let n = duv.xyz().crossed(dvv.xyz());
    let m = n.modulus();
    if m > 1e-30 {
        GpVec::new(n.x / m, n.y / m, n.z / m)
    } else {
        GpVec::zero()
    }
}

/// Whether all sampled surface normals are parallel (within `tol`).
fn is_planar(s: &dyn Surface, nu: usize, nv: usize, tol: f64) -> bool {
    let (u0, u1, v0, v1) = sample_bounds(s);
    let (mut first, mut first_ok) = (GpVec::zero(), false);
    for i in 0..nu {
        for j in 0..nv {
            let u = u0 + (u1 - u0) * i as f64 / (nu.max(1) - 1) as f64;
            let v = v0 + (v1 - v0) * j as f64 / (nv.max(1) - 1) as f64;
            let n = surface_normal(s, u, v);
            if n.xyz().square_modulus() < 1e-30 {
                continue;
            }
            if !first_ok {
                first = n;
                first_ok = true;
            } else if n.xyz().crossed(first.xyz()).modulus() > tol {
                return false;
            }
        }
    }
    first_ok
}

/// 3×3 determinant.
fn det3(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

/// Solve a 3×3 linear system via Cramer's rule.
fn solve3(a: &[[f64; 3]; 3], rhs: &[f64; 3]) -> Option<[f64; 3]> {
    let d = det3(a);
    if d.abs() < 1e-300 {
        return None;
    }
    let mut x = [0.0; 3];
    for k in 0..3 {
        let mut m = *a;
        for i in 0..3 {
            m[i][k] = rhs[i];
        }
        x[k] = det3(&m) / d;
    }
    Some(x)
}

/// Circumcenter of three non-collinear 3D points, if it exists.
fn circumcenter(a: &GpPnt, b: &GpPnt, c: &GpPnt) -> Option<GpPnt> {
    let d1 = GpVec::from_pnts(a, b);
    let d2 = GpVec::from_pnts(a, c);
    let n = d1.xyz().crossed(d2.xyz());
    if n.modulus() < 1e-20 {
        return None;
    }
    let n2 = |p: &GpPnt| p.coord.dot(&p.coord);
    let mat = [
        [d1.xyz().x, d1.xyz().y, d1.xyz().z],
        [d2.xyz().x, d2.xyz().y, d2.xyz().z],
        [n.x, n.y, n.z],
    ];
    let rhs = [0.5 * (n2(b) - n2(a)), 0.5 * (n2(c) - n2(a)), a.coord.dot(&n)];
    let o = solve3(&mat, &rhs)?;
    Some(GpPnt::new(o[0], o[1], o[2]))
}

/// If the surface is spherical, return its center (mid-latitude samples, away
/// from the poles where the longitude parameter collapses).
fn sphere_center(s: &dyn Surface) -> Option<GpPnt> {
    let (u0, u1, v0, v1) = sample_bounds(s);
    let u_mid = 0.5 * (u0 + u1);
    let v_mid = 0.5 * (v0 + v1);
    let span = (v1 - v0).abs().max(1e-6);
    let a = s.d0(u0, v_mid);
    let b = s.d0(u_mid, v_mid);
    let c = s.d0(u0, v_mid + 0.25 * span);
    circumcenter(&a, &b, &c)
}

/// Classify a planar surface into a `GpPln` (reconstructed normal + location;
/// the in-plane axes are chosen arbitrarily, which does not affect distances).
pub(crate) fn classify_plane(s: &dyn Surface) -> Option<GpPln> {
    if !is_planar(s, 8, 8, 1e-6) {
        return None;
    }
    let (u0, _, v0, _) = sample_bounds(s);
    let o = s.d0(u0, v0);
    let n = surface_normal(s, u0, v0);
    let d = GpDir::from_vec(&n).ok()?;
    let z_axis = GpDir::new(0.0, 0.0, 1.0).ok()?;
    let x_dir = if d.is_normal(&z_axis) {
        z_axis
    } else {
        GpDir::new(1.0, 0.0, 0.0).ok()?
    };
    GpAx3::new(o, d, &x_dir).ok().map(GpPln::new)
}

/// Classify a spherical surface into a `GpSphere`. The reconstructed axes are
/// standard at the solved center; distances are axis-independent.
pub(crate) fn classify_sphere(s: &dyn Surface) -> Option<GpSphere> {
    let center = sphere_center(s)?;
    let (u0, u1, v0, v1) = sample_bounds(s);
    let r0 = s.d0(u0, v0).distance(&center);
    if r0 < 1e-12 {
        return None;
    }
    let (nu, nv) = (8, 8);
    for i in 0..nu {
        for j in 0..nv {
            let u = u0 + (u1 - u0) * i as f64 / (nu - 1) as f64;
            let v = v0 + (v1 - v0) * j as f64 / (nv - 1) as f64;
            if (s.d0(u, v).distance(&center) - r0).abs() > 1e-4 * r0.abs().max(1.0) {
                return None;
            }
        }
    }
    let mut pos = GpAx3::standard();
    pos.set_location(center);
    Some(GpSphere { pos, radius: r0 })
}

/// Constant-tangent-direction detection ⇒ line.
fn is_line(c: &dyn Curve) -> bool {
    let (a, b) = (c.first_parameter(), c.last_parameter());
    let us: Vec<f64> = if a.is_finite() && b.is_finite() && b > a {
        (0..=5).map(|i| a + (b - a) * i as f64 / 5.0).collect()
    } else {
        vec![-2.0, -1.0, 0.0, 1.0, 2.0]
    };
    let t0 = c.d1(us[0]).1;
    let m0 = t0.magnitude();
    if m0 < 1e-9 {
        return false;
    }
    for u in &us[1..] {
        let t = c.d1(*u).1;
        let m = t.magnitude();
        if m < 1e-9 {
            return false;
        }
        if t0.cross_magnitude(&t) > 1e-7 * m0 * m {
            return false;
        }
    }
    true
}

/// Reconstruct a `GpLin` from a line curve.
fn reconstruct_line(c: &dyn Curve) -> Option<GpLin> {
    if !is_line(c) {
        return None;
    }
    let u0 = if c.first_parameter().is_finite() {
        c.first_parameter()
    } else {
        0.0
    };
    let p0 = c.d0(u0);
    let t = c.d1(u0).1;
    let d = GpDir::from_vec(&t).ok()?;
    Some(GpLin::from_pnt_dir(p0, d))
}

// ---------------------------------------------------------------------------
// General Newton path (point-surface, port of `Extrema_GenExtPS` +
// `Extrema_FuncPSNorm`, with a numeric Jacobian).
// ---------------------------------------------------------------------------

/// F(u, v) = ((S(u,v)-P)·Su, (S(u,v)-P)·Sv).
fn ps_f(s: &dyn Surface, p: &GpPnt, u: f64, v: f64) -> [f64; 2] {
    let (q, su, sv) = s.d1(u, v);
    let w = GpVec::from_pnts(p, &q);
    [w.dot(&su), w.dot(&sv)]
}

/// Numeric 2×2 Jacobian of `ps_f` by central finite differences of `d1`.
fn ps_jac(s: &dyn Surface, p: &GpPnt, u: f64, v: f64, hu: f64, hv: f64) -> [[f64; 2]; 2] {
    let fu_p = ps_f(s, p, u + hu, v);
    let fu_m = ps_f(s, p, u - hu, v);
    let fv_p = ps_f(s, p, u, v + hv);
    let fv_m = ps_f(s, p, u, v - hv);
    [
        [(fu_p[0] - fu_m[0]) / (2.0 * hu), (fv_p[0] - fv_m[0]) / (2.0 * hv)],
        [(fu_p[1] - fu_m[1]) / (2.0 * hu), (fv_p[1] - fv_m[1]) / (2.0 * hv)],
    ]
}

/// Newton on F = 0 from a seed, clamped to the parameter box. Returns the
/// converged stationary point, or `None` if the iterate diverged.
fn solve_point_surface(
    s: &dyn Surface,
    p: &GpPnt,
    u0: f64,
    v0: f64,
    ua: f64,
    ub: f64,
    va: f64,
    vb: f64,
) -> Option<(f64, f64)> {
    let cl = |x: f64, a: f64, b: f64| x.clamp(a.min(b), a.max(b));
    let (mut u, mut v) = (cl(u0, ua, ub), cl(v0, va, vb));
    let hu = ((ub - ua).abs() * 1e-6).max(1e-9);
    let hv = ((vb - va).abs() * 1e-6).max(1e-9);
    let base = ps_f(s, p, u, v).iter().map(|x| x.abs()).sum::<f64>();
    let mut converged = false;
    for _ in 0..40 {
        let f = ps_f(s, p, u, v);
        let j = ps_jac(s, p, u, v, hu, hv);
        let det = j[0][0] * j[1][1] - j[0][1] * j[1][0];
        if det.abs() < 1e-300 {
            break;
        }
        let du = (-f[0] * j[1][1] + f[1] * j[0][1]) / det;
        let dv = (-j[0][0] * f[1] + j[1][0] * f[0]) / det;
        let (un, vn) = (cl(u + du, ua, ub), cl(v + dv, va, vb));
        if (un - u).abs() < 1e-10 * (1.0 + u.abs())
            && (vn - v).abs() < 1e-10 * (1.0 + v.abs())
        {
            u = un;
            v = vn;
            converged = true;
            break;
        }
        if (un - u).abs() + (vn - v).abs() > 1e6 * (1.0 + u.abs() + v.abs()) {
            return None; // diverged
        }
        u = un;
        v = vn;
    }
    // A converged step locates a stationary point of the (finite-difference)
    // system even when the residual cannot reach machine zero at a boundary
    // (one-sided derivatives); accept it. Otherwise require a genuine root.
    if converged {
        return Some((u, v));
    }
    let resid = ps_f(s, p, u, v).iter().map(|x| x.abs()).sum::<f64>();
    if resid < base * 1e-2 + 1e-6 {
        Some((u, v))
    } else {
        None
    }
}

/// All local extrema of |S(u,v)-P| via grid seeding + Newton, deduplicated and
/// sorted by distance. Replaces the sampling path of `Extrema_ExtPS` for
/// non-analytic surfaces.
pub(crate) fn point_surface_newton_all(s: &dyn Surface, p: &GpPnt) -> Vec<ExtremaPair> {
    let (u0, u1) = surf_bound_u(s);
    let (v0, v1) = surf_bound_v(s);
    let (nu, nv) = (24, 24);
    let mut d2 = vec![vec![0.0; nv + 1]; nu + 1];
    for i in 0..=nu {
        for j in 0..=nv {
            let u = u0 + (u1 - u0) * i as f64 / nu as f64;
            let v = v0 + (v1 - v0) * j as f64 / nv as f64;
            let d = s.d0(u, v).square_distance(p);
            d2[i][j] = if d.is_finite() { d } else { f64::INFINITY };
        }
    }
    let mut seeds: Vec<(f64, f64)> = Vec::new();
    // Local extrema of the sampled squared distance (both minima and maxima of
    // the distance are stationary points of F = 0).
    for i in 1..nu {
        for j in 1..nv {
            let d = d2[i][j];
            let (du_m, du_p, dv_m, dv_p) = (d2[i - 1][j], d2[i + 1][j], d2[i][j - 1], d2[i][j + 1]);
            if d <= du_m && d <= du_p && d <= dv_m && d <= dv_p {
                seeds.push((u0 + (u1 - u0) * i as f64 / nu as f64, v0 + (v1 - v0) * j as f64 / nv as f64));
            } else if d >= du_m && d >= du_p && d >= dv_m && d >= dv_p {
                seeds.push((u0 + (u1 - u0) * i as f64 / nu as f64, v0 + (v1 - v0) * j as f64 / nv as f64));
            }
        }
    }
    // Global grid min/max + corners (boundary extrema for trimmed surfaces).
    let (mut gmin, mut gmax) = ((0usize, 0usize), (0usize, 0usize));
    for i in 0..=nu {
        for j in 0..=nv {
            if d2[i][j] < d2[gmin.0][gmin.1] {
                gmin = (i, j);
            }
            if d2[i][j] > d2[gmax.0][gmax.1] {
                gmax = (i, j);
            }
        }
    }
    seeds.push((u0 + (u1 - u0) * gmin.0 as f64 / nu as f64, v0 + (v1 - v0) * gmin.1 as f64 / nv as f64));
    seeds.push((u0 + (u1 - u0) * gmax.0 as f64 / nu as f64, v0 + (v1 - v0) * gmax.1 as f64 / nv as f64));
    seeds.push((u0, v0));
    seeds.push((u0, v1));
    seeds.push((u1, v0));
    seeds.push((u1, v1));

    let mut out: Vec<ExtremaPair> = Vec::new();
    for (u, v) in seeds {
        if let Some((uu, vv)) = solve_point_surface(s, p, u, v, u0, u1, v0, v1) {
            let q = s.d0(uu, vv);
            if !q.x().is_finite() {
                continue;
            }
            if !out.iter().any(|e: &ExtremaPair| e.p2.distance(&q) < 1e-6) {
                out.push(ps_pair(p, uu, vv, q));
            }
        }
    }
    out.sort_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap_or(Ordering::Equal));
    out
}

/// Fallback: coarse grid min + coordinate descent (degenerate surfaces).
fn fallback_point_surface(s: &dyn Surface, p: &GpPnt) -> ExtremaPair {
    let (u0, u1) = surf_bound_u(s);
    let (v0, v1) = surf_bound_v(s);
    let n = 40;
    let mut best = (0.0, 0.0, f64::INFINITY);
    for i in 0..=n {
        for j in 0..=n {
            let u = u0 + (u1 - u0) * i as f64 / n as f64;
            let v = v0 + (v1 - v0) * j as f64 / n as f64;
            let d = s.d0(u, v).square_distance(p);
            if d < best.2 {
                best = (u, v, d);
            }
        }
    }
    let (mut u, mut v) = (best.0, best.1);
    let dist2 = |u: f64, v: f64| s.d0(u, v).square_distance(p);
    let mut step = (v1 - v0).abs().max(1.0) / 4.0;
    for _ in 0..256 {
        let cur = dist2(u, v);
        let cand = [(u + step, v), (u - step, v), (u, v + step), (u, v - step)];
        let mut bi = usize::MAX;
        let mut bb = cur;
        for (k, (cu, cv)) in cand.iter().enumerate() {
            let dd = dist2(*cu, *cv);
            if dd < bb {
                bi = k;
                bb = dd;
            }
        }
        if bi == usize::MAX {
            step *= 0.5;
            if step < 1e-12 {
                break;
            }
        } else {
            match bi {
                0 => u += step,
                1 => u -= step,
                2 => v += step,
                _ => v -= step,
            }
            u = u.clamp(u0.min(u1), u0.max(u1));
            v = v.clamp(v0.min(v1), v0.max(v1));
        }
    }
    let q = s.d0(u, v);
    ps_pair(p, u, v, q)
}

// ---------------------------------------------------------------------------
// General Newton path (curve-surface, port of `Extrema_GenExtCS` +
// `Extrema_FuncExtCS`, numeric 3×3 Jacobian).
// ---------------------------------------------------------------------------

/// F(t, u, v) = ((C-S)·C′, (C-S)·Su, (C-S)·Sv).
fn cs_f(c: &dyn Curve, s: &dyn Surface, t: f64, u: f64, v: f64) -> Option<[f64; 3]> {
    let (pc, dc) = c.d1(t);
    let (ps, su, sv) = s.d1(u, v);
    if !(pc.x().is_finite() && ps.x().is_finite()) {
        return None;
    }
    let w = GpVec::from_pnts(&ps, &pc);
    Some([w.dot(&dc), w.dot(&su), w.dot(&sv)])
}

/// Newton on F = 0 (3×3) from a seed, clamped to the parameter box.
fn solve_curve_surface(
    c: &dyn Curve,
    s: &dyn Surface,
    t0: f64,
    u0: f64,
    v0: f64,
    ta: f64,
    tb: f64,
    ua: f64,
    ub: f64,
    va: f64,
    vb: f64,
) -> Option<(f64, f64, f64)> {
    let cl = |x: f64, a: f64, b: f64| x.clamp(a.min(b), a.max(b));
    let mut x = [cl(t0, ta, tb), cl(u0, ua, ub), cl(v0, va, vb)];
    let hs = [
        ((tb - ta).abs() * 1e-6).max(1e-9),
        ((ub - ua).abs() * 1e-6).max(1e-9),
        ((vb - va).abs() * 1e-6).max(1e-9),
    ];
    let base = cs_f(c, s, x[0], x[1], x[2])?.iter().map(|x| x.abs()).sum::<f64>();
    let mut converged = false;
    for _ in 0..40 {
        let f = cs_f(c, s, x[0], x[1], x[2])?;
        let mut j = [[0.0; 3]; 3];
        for k in 0..3 {
            let h = hs[k];
            let mut xp = x;
            xp[k] += h;
            let mut xm = x;
            xm[k] -= h;
            let fp = cs_f(c, s, xp[0], xp[1], xp[2])?;
            let fm = cs_f(c, s, xm[0], xm[1], xm[2])?;
            for r in 0..3 {
                j[r][k] = (fp[r] - fm[r]) / (2.0 * h);
            }
        }
        let d = match solve3(&j, &f) {
            Some(d) => d,
            None => break,
        };
        let xn = [
            cl(x[0] - d[0], ta, tb),
            cl(x[1] - d[1], ua, ub),
            cl(x[2] - d[2], va, vb),
        ];
        if (xn[0] - x[0]).abs() < 1e-10 * (1.0 + x[0].abs())
            && (xn[1] - x[1]).abs() < 1e-10 * (1.0 + x[1].abs())
            && (xn[2] - x[2]).abs() < 1e-10 * (1.0 + x[2].abs())
        {
            x = xn;
            converged = true;
            break;
        }
        if (xn[0] - x[0]).abs() + (xn[1] - x[1]).abs() + (xn[2] - x[2]).abs()
            > 1e6 * (1.0 + x[0].abs() + x[1].abs() + x[2].abs())
        {
            return None;
        }
        x = xn;
    }
    if converged {
        return Some((x[0], x[1], x[2]));
    }
    let resid = cs_f(c, s, x[0], x[1], x[2])?.iter().map(|x| x.abs()).sum::<f64>();
    if resid < base * 1e-2 + 1e-6 {
        Some((x[0], x[1], x[2]))
    } else {
        None
    }
}

/// All local extrema of the curve-surface distance via grid seeding + Newton.
pub(crate) fn curve_surface_newton_all(c: &dyn Curve, s: &dyn Surface) -> Vec<ExtremaPair> {
    let (t0, t1) = curve_bound(c);
    let (u0, u1) = surf_bound_u(s);
    let (v0, v1) = surf_bound_v(s);
    let (nt, nu, nv) = (15, 11, 11);
    let mut d2 = vec![vec![vec![0.0; nv + 1]; nu + 1]; nt + 1];
    for i in 0..=nt {
        for j in 0..=nu {
            for k in 0..=nv {
                let t = t0 + (t1 - t0) * i as f64 / nt as f64;
                let u = u0 + (u1 - u0) * j as f64 / nu as f64;
                let v = v0 + (v1 - v0) * k as f64 / nv as f64;
                let d = c.d0(t).square_distance(&s.d0(u, v));
                d2[i][j][k] = if d.is_finite() { d } else { f64::INFINITY };
            }
        }
    }
    let at = |i: usize, j: usize, k: usize| (t0 + (t1 - t0) * i as f64 / nt as f64, u0 + (u1 - u0) * j as f64 / nu as f64, v0 + (v1 - v0) * k as f64 / nv as f64);
    let mut seeds: Vec<(f64, f64, f64)> = Vec::new();
    for i in 1..nt {
        for j in 1..nu {
            for k in 1..nv {
                let d = d2[i][j][k];
                let lo = d2[i - 1][j][k] >= d
                    && d2[i + 1][j][k] >= d
                    && d2[i][j - 1][k] >= d
                    && d2[i][j + 1][k] >= d
                    && d2[i][j][k - 1] >= d
                    && d2[i][j][k + 1] >= d;
                let hi = d2[i - 1][j][k] <= d
                    && d2[i + 1][j][k] <= d
                    && d2[i][j - 1][k] <= d
                    && d2[i][j + 1][k] <= d
                    && d2[i][j][k - 1] <= d
                    && d2[i][j][k + 1] <= d;
                if lo || hi {
                    seeds.push(at(i, j, k));
                }
            }
        }
    }
    let (mut gmin, mut gmax) = ((0usize, 0usize, 0usize), (0usize, 0usize, 0usize));
    for i in 0..=nt {
        for j in 0..=nu {
            for k in 0..=nv {
                if d2[i][j][k] < d2[gmin.0][gmin.1][gmin.2] {
                    gmin = (i, j, k);
                }
                if d2[i][j][k] > d2[gmax.0][gmax.1][gmax.2] {
                    gmax = (i, j, k);
                }
            }
        }
    }
    seeds.push(at(gmin.0, gmin.1, gmin.2));
    seeds.push(at(gmax.0, gmax.1, gmax.2));

    let mut out: Vec<ExtremaPair> = Vec::new();
    for (t, u, v) in seeds {
        if let Some((tt, uu, vv)) =
            solve_curve_surface(c, s, t, u, v, t0, t1, u0, u1, v0, v1)
        {
            let pc = c.d0(tt);
            let ps = s.d0(uu, vv);
            if !pc.x().is_finite() || !ps.x().is_finite() {
                continue;
            }
            if !out
                .iter()
                .any(|e: &ExtremaPair| e.p1.distance(&pc) < 1e-6 && e.p2.distance(&ps) < 1e-6)
            {
                out.push(cs_pair(pc, tt, ps, uu, vv, pc.distance(&ps)));
            }
        }
    }
    out.sort_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap_or(Ordering::Equal));
    out
}

/// Fallback: grid min (degenerate curve-surface pairs).
fn fallback_curve_surface(c: &dyn Curve, s: &dyn Surface) -> ExtremaPair {
    let (t0, t1) = curve_bound(c);
    let (u0, u1) = surf_bound_u(s);
    let (v0, v1) = surf_bound_v(s);
    let (nt, nu, nv) = (20, 16, 16);
    let mut best = (0.0, 0.0, 0.0, f64::INFINITY);
    for i in 0..=nt {
        for j in 0..=nu {
            for k in 0..=nv {
                let t = t0 + (t1 - t0) * i as f64 / nt as f64;
                let u = u0 + (u1 - u0) * j as f64 / nu as f64;
                let v = v0 + (v1 - v0) * k as f64 / nv as f64;
                let d = c.d0(t).square_distance(&s.d0(u, v));
                if d < best.3 {
                    best = (t, u, v, d);
                }
            }
        }
    }
    let pc = c.d0(best.0);
    let ps = s.d0(best.1, best.2);
    cs_pair(pc, best.0, ps, best.1, best.2, pc.distance(&ps))
}

// ---------------------------------------------------------------------------
// Public dispatch.
// ---------------------------------------------------------------------------

/// All local extrema of the point-surface distance, deduplicated and sorted.
/// Planes and spheres are classified and solved analytically (exact); every
/// other surface goes through the grid + Newton path.
pub fn point_surface_extrema_all(s: &dyn Surface, p: &GpPnt) -> Vec<ExtremaPair> {
    if let Some(pl) = classify_plane(s) {
        return vec![point_plane_extrema(&pl, p)];
    }
    if let Some(sp) = classify_sphere(s) {
        let mut v = point_sphere_extrema(&sp, p);
        v.sort_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap_or(Ordering::Equal));
        return v;
    }
    point_surface_newton_all(s, p)
}

/// Minimum distance from `p` to `s` (with the closest point and parameters).
pub fn point_surface_extrema(s: &dyn Surface, p: &GpPnt) -> ExtremaPair {
    match point_surface_extrema_all(s, p).into_iter().next() {
        Some(e) => e,
        None => fallback_point_surface(s, p),
    }
}

/// All local extrema of the curve-surface distance, deduplicated and sorted.
/// Lines against planes/spheres are classified and solved analytically; all
/// other pairs go through the grid + Newton path.
pub fn curve_surface_extrema_all(c: &dyn Curve, s: &dyn Surface) -> Vec<ExtremaPair> {
    if let Some(l) = reconstruct_line(c) {
        if let Some(sp) = classify_sphere(s) {
            return line_sphere_extrema(&l, &sp);
        }
        if let Some(pl) = classify_plane(s) {
            return line_plane_extrema(&l, &pl);
        }
    }
    curve_surface_newton_all(c, s)
}

/// Minimum distance between curve `c` and surface `s`.
pub fn curve_surface_extrema(c: &dyn Curve, s: &dyn Surface) -> ExtremaPair {
    match curve_surface_extrema_all(c, s).into_iter().next() {
        Some(e) => e,
        None => fallback_curve_surface(c, s),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bspline_surface::{bspline_surface_uniform_knots, fit_surface_grid, GeomBSplineSurface};
    use crate::{GeomBSplineCurve, GeomLine, GeomPlane, GeomSphere};
    use occt_core::gp::{GpAx3, GpCone, GpCylinder, GpDir, GpLin, GpPnt, GpSphere as GpSphereT, GpTorus};

    const PI: f64 = std::f64::consts::PI;

    fn unit_sphere() -> GpSphere {
        GpSphereT::new(GpAx3::standard(), 1.0).unwrap()
    }

    #[test]
    fn point_sphere_extrema_min_and_max() {
        let sp = unit_sphere();
        let all = point_sphere_extrema(&sp, &GpPnt::new(3.0, 0.0, 0.0));
        assert_eq!(all.len(), 2, "min+max: {all:?}");
        let min = all.iter().min_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap()).unwrap();
        let max = all.iter().max_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap()).unwrap();
        assert!((min.distance - 2.0).abs() < 1e-9, "min {}", min.distance);
        assert!((min.p2.x() - 1.0).abs() < 1e-9, "closest {min:?}");
        assert!((max.distance - 4.0).abs() < 1e-9, "max {}", max.distance);
        assert!((max.p2.x() + 1.0).abs() < 1e-9, "farthest {max:?}");
    }

    #[test]
    fn point_plane_extrema_distance() {
        let mut pl = GpPln::new(GpAx3::standard());
        pl.set_location(&GpPnt::new(0.0, 0.0, 5.0));
        let e = point_plane_extrema(&pl, &GpPnt::new(0.0, 0.0, 0.0));
        assert!((e.distance - 5.0).abs() < 1e-9, "dist {}", e.distance);
        assert!((e.p2.z() - 5.0).abs() < 1e-9, "closest z {}", e.p2.z());
    }

    #[test]
    fn point_cylinder_extrema_min() {
        let cy = GpCylinder::new(GpAx3::standard(), 1.0).unwrap();
        let all = point_cylinder_extrema(&cy, &GpPnt::new(3.0, 0.0, 0.0));
        assert_eq!(all.len(), 2);
        let min = all.iter().min_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap()).unwrap();
        assert!((min.distance - 2.0).abs() < 1e-9, "min {}", min.distance);
        assert!((min.p2.x() - 1.0).abs() < 1e-9, "closest {min:?}");
    }

    #[test]
    fn point_cone_extrema_min() {
        // Cone placement at origin, RefRadius 1, semi-angle 45° → vertex at
        // (0,0,-1); in the XZ plane the generatrix is ρ = 1 + z. Point
        // (0.5, 0, 0.1): closest point at ρ = 0.8, z = -0.2, distance² = 0.18.
        let co = GpCone::new(GpAx3::standard(), 1.0, PI / 4.0).unwrap();
        let all = point_cone_extrema(&co, &GpPnt::new(0.5, 0.0, 0.1));
        assert_eq!(all.len(), 2, "cone extrema {all:?}");
        let min = all.iter().min_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap()).unwrap();
        assert!((min.distance - (0.18f64).sqrt()).abs() < 1e-7, "min {}", min.distance);
        assert!((min.p2.x() - 0.8).abs() < 1e-6 && (min.p2.z() - (-0.2)).abs() < 1e-6, "closest {min:?}");
    }

    #[test]
    fn point_torus_extrema_min() {
        // Major 3, minor 1; point (6,0,0): torus crosses the X axis at
        // x = 4, 2, -2, -4 → min distance 2.
        let to = GpTorus::new(GpAx3::standard(), 3.0, 1.0).unwrap();
        let all = point_torus_extrema(&to, &GpPnt::new(6.0, 0.0, 0.0));
        assert_eq!(all.len(), 4, "torus extrema {all:?}");
        let min = all.iter().min_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap()).unwrap();
        assert!((min.distance - 2.0).abs() < 1e-9, "min {}", min.distance);
        assert!((min.p2.x() - 4.0).abs() < 1e-9, "closest {min:?}");
    }

    #[test]
    fn sphere_parameters_roundtrip() {
        let sp = unit_sphere();
        for (u, v) in [(0.0, 0.0), (PI / 2.0, PI / 4.0), (PI, -PI / 3.0), (3.0 * PI / 2.0, 0.2)] {
            let q = slib::sphere_value(&sp, u, v);
            let (u2, v2) = sphere_parameters(&sp, &q);
            let q2 = slib::sphere_value(&sp, u2, v2);
            assert!(q.distance(&q2) < 1e-9, "params roundtrip (u={u},v={v})");
        }
    }

    #[test]
    fn curve_surface_line_sphere_min() {
        let line = GeomLine::new(GpLin::from_pnt_dir(
            GpPnt::new(0.0, 3.0, 0.0),
            GpDir::new(1.0, 0.0, 0.0).unwrap(),
        ));
        let sphere = GeomSphere::new(unit_sphere());
        let e = curve_surface_extrema(&line, &sphere);
        assert!((e.distance - 2.0).abs() < 1e-9, "min {}", e.distance);
        assert!((e.p1.y() - 3.0).abs() < 1e-9 && e.p1.x().abs() < 1e-9, "curve point {e:?}");
        assert!((e.p2.x() - 0.0).abs() < 1e-9 && (e.p2.y() - 1.0).abs() < 1e-9, "surf point {e:?}");
    }

    #[test]
    fn curve_surface_line_plane_intersect() {
        // Line (0,0,3) + t·(0,0,-1) pierces the plane z=0 at t=3 → distance 0.
        let line = GeomLine::new(GpLin::from_pnt_dir(
            GpPnt::new(0.0, 0.0, 3.0),
            GpDir::new(0.0, 0.0, -1.0).unwrap(),
        ));
        let pl = GpPln::new(GpAx3::standard());
        let e = curve_surface_extrema(&line, &GeomPlane::new(pl));
        assert!(e.distance.abs() < 1e-9, "line meets plane: {}", e.distance);
    }

    fn paraboloid() -> GeomBSplineSurface {
        let (nu, nv) = (3, 3);
        let points: Vec<Vec<GpPnt>> = (0..nu)
            .map(|i| {
                (0..nv)
                    .map(|j| {
                        let u = i as f64 / (nu - 1) as f64;
                        let v = j as f64 / (nv - 1) as f64;
                        GpPnt::new(u, v, u * u + v * v)
                    })
                    .collect()
            })
            .collect();
        fit_surface_grid(&points, 2, 2).unwrap()
    }

    #[test]
    fn newton_path_bspline_paraboloid_min() {
        // S(u,v) = (u, v, u²+v²) exactly (degree-2 Bernstein reproduces the
        // quadratic). Point P = (0.5, 0.5, -1). The closest point solves
        // F = 0; by symmetry u = v = t with 4t³ + 3t - 0.5 = 0 → t ≈ 0.16071.
        let s = paraboloid();
        let p = GpPnt::new(0.5, 0.5, -1.0);
        let all = point_surface_newton_all(&s, &p);
        assert!(!all.is_empty(), "no extrema found");
        let min = all.iter().min_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap()).unwrap();
        // Independent reference: dense scan over the same domain.
        let mut ref_min = f64::INFINITY;
        for i in 0..=200 {
            for j in 0..=200 {
                let u = i as f64 / 200.0;
                let v = j as f64 / 200.0;
                let d = s.d0(u, v).distance(&p);
                if d < ref_min {
                    ref_min = d;
                }
            }
        }
        assert!(
            (min.distance - ref_min).abs() < 1e-4,
            "newton {} vs reference {}",
            min.distance,
            ref_min
        );
        // Orthogonality condition holds at the closest point.
        let f = ps_f(&s, &p, min.u1, min.v2.unwrap());
        assert!(
            f[0].abs() + f[1].abs() < 1e-6,
            "orthogonality F={f:?} at u={},v={}",
            min.u1,
            min.v2.unwrap()
        );
        // Local min: neighbours are farther.
        let (u, v) = (min.u1, min.v2.unwrap());
        let eps = 1e-3;
        let du = p.distance(&s.d0(u + eps, v));
        let dv = p.distance(&s.d0(u, v + eps));
        assert!(min.distance <= du + 1e-9 && min.distance <= dv + 1e-9, "not a local min");
    }

    #[test]
    fn point_surface_extrema_dispatch_sphere() {
        let sphere = GeomSphere::new(unit_sphere());
        let e = point_surface_extrema(&sphere, &GpPnt::new(3.0, 0.0, 0.0));
        assert!((e.distance - 2.0).abs() < 1e-9, "min {}", e.distance);
    }

    #[test]
    fn classification_plane_and_sphere() {
        let sphere = GeomSphere::new(unit_sphere());
        assert!(classify_sphere(&sphere).is_some());
        assert!(classify_plane(&sphere).is_none());
        let plane = GeomPlane::new(GpPln::new(GpAx3::standard()));
        assert!(classify_plane(&plane).is_some());
        assert!(classify_sphere(&plane).is_none());
    }

    #[test]
    fn curve_surface_newton_parabola_plane_min() {
        // Curve c(t) = (t, 0, 1+t²) over [0,1] (degree-2 B-spline), plane
        // z = 0. Minimum distance 1 at t = 0, with orthogonality F ≈ 0.
        // Degree-2 B-spline through control points (0,0,1), (0.5,0,1),
        // (1,0,2) reproduces c(t) = (t, 0, 1+t²) exactly.
        let c = GeomBSplineCurve::new(
            vec![
                GpPnt::new(0.0, 0.0, 1.0),
                GpPnt::new(0.5, 0.0, 1.0),
                GpPnt::new(1.0, 0.0, 2.0),
            ],
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            2,
        )
        .unwrap();
        let s = GeomBSplineSurface::new(
            vec![
                vec![GpPnt::new(0.0, 0.0, 0.0), GpPnt::new(0.0, 1.0, 0.0)],
                vec![GpPnt::new(1.0, 0.0, 0.0), GpPnt::new(1.0, 1.0, 0.0)],
            ],
            bspline_surface_uniform_knots(2, 2, 1, 1).0,
            bspline_surface_uniform_knots(2, 2, 1, 1).1,
            1,
            1,
        )
        .unwrap();
        let all = curve_surface_newton_all(&c, &s);
        assert!(!all.is_empty(), "no curve-surface extrema");
        let min = all.iter().min_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap()).unwrap();
        assert!((min.distance - 1.0).abs() < 1e-6, "min {}", min.distance);
        let f = cs_f(&c, &s, min.u1, min.u2, min.v2.unwrap()).unwrap();
        assert!(
            f[0].abs() + f[1].abs() + f[2].abs() < 1e-6,
            "orthogonality F={f:?}"
        );
    }
}
