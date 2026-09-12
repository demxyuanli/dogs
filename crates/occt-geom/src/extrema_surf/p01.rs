use super::prelude::*;
use super::*;

/// Angle snap: OCCT's `ExtPElS_MyEps = Epsilon(2π)`; `ANGULAR` (1e-12) is the
/// crate's angular resolution and is used here for the same purpose.

pub(super) const SNAP_EPS: f64 = 1e-12;

// ---------------------------------------------------------------------------
// Pair helpers.
// ---------------------------------------------------------------------------

/// Point-surface pair: `p1` is the given point, `p2` the surface point,
/// surface parameters go in `u2`/`v2` (`u1` mirrors `u2`, per the existing
/// `extrema::ExtremaPair` convention for point–surface results).
pub(super) fn ps_pair(p: &GpPnt, u: f64, v: f64, q: GpPnt) -> ExtremaPair {
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
pub(super) fn cs_pair(pc: GpPnt, t: f64, ps: GpPnt, u: f64, v: f64, d: f64) -> ExtremaPair {
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
pub(super) fn vec_angle_with_ref(a: &GpVec, b: &GpVec, vref: &GpVec) -> f64 {
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

pub(super) fn curve_bound(c: &dyn Curve) -> (f64, f64) {
    let (a, b) = (c.first_parameter(), c.last_parameter());
    if a.is_finite() && b.is_finite() && b > a {
        (a, b)
    } else {
        (-10.0, 10.0)
    }
}

/// Finite, sane sampling bounds for a surface (unbounded ranges clamp to ±1).
pub(super) fn sample_bounds(s: &dyn Surface) -> (f64, f64, f64, f64) {
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
pub(super) fn surface_normal(s: &dyn Surface, u: f64, v: f64) -> GpVec {
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
pub(super) fn is_planar(s: &dyn Surface, nu: usize, nv: usize, tol: f64) -> bool {
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
pub(super) fn det3(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

/// Solve a 3×3 linear system via Cramer's rule.
pub(super) fn solve3(a: &[[f64; 3]; 3], rhs: &[f64; 3]) -> Option<[f64; 3]> {
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
pub(super) fn circumcenter(a: &GpPnt, b: &GpPnt, c: &GpPnt) -> Option<GpPnt> {
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
pub(super) fn sphere_center(s: &dyn Surface) -> Option<GpPnt> {
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
pub(super) fn is_line(c: &dyn Curve) -> bool {
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
pub(super) fn reconstruct_line(c: &dyn Curve) -> Option<GpLin> {
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
pub(super) fn ps_f(s: &dyn Surface, p: &GpPnt, u: f64, v: f64) -> [f64; 2] {
    let (q, su, sv) = s.d1(u, v);
    let w = GpVec::from_pnts(p, &q);
    [w.dot(&su), w.dot(&sv)]
}

/// Numeric 2×2 Jacobian of `ps_f` by central finite differences of `d1`.
pub(super) fn ps_jac(s: &dyn Surface, p: &GpPnt, u: f64, v: f64, hu: f64, hv: f64) -> [[f64; 2]; 2] {
    let fu_p = ps_f(s, p, u + hu, v);
    let fu_m = ps_f(s, p, u - hu, v);
    let fv_p = ps_f(s, p, u, v + hv);
    let fv_m = ps_f(s, p, u, v - hv);
    [
        [(fu_p[0] - fu_m[0]) / (2.0 * hu), (fv_p[0] - fv_m[0]) / (2.0 * hv)],
        [(fu_p[1] - fu_m[1]) / (2.0 * hu), (fv_p[1] - fv_m[1]) / (2.0 * hv)],
    ]
}
