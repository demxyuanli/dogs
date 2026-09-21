use super::prelude::*;
use super::*;

/// Whether the curve's tangent direction is constant over `[a, b]` (a straight
/// line): every sampled second derivative is ~0.
pub(super) fn is_line_curve(curve: &dyn occt_geom::Curve, a: f64, b: f64) -> bool {
    for i in 0..=4 {
        let t = a + (b - a) * i as f64 / 4.0;
        let (_, _, d2) = curve.d2(t);
        if d2.square_magnitude() > 1e-18 {
            return false;
        }
    }
    true
}

/// Radius of a circular arc over `[a, b]`, `None` when the sampled points are
/// not on a circle (constant tangent speed and curvature, planar binormal).
/// The tangent speed of a `GeomCircle` equals its radius, so the arc length is
/// `radius · |b − a|`.
pub(super) fn circle_radius(curve: &dyn occt_geom::Curve, a: f64, b: f64) -> Option<f64> {
    let mut speeds = Vec::with_capacity(6);
    let mut accs = Vec::with_capacity(6);
    let mut binorms: Vec<GpVec> = Vec::with_capacity(6);
    for i in 0..=5 {
        let t = a + (b - a) * i as f64 / 5.0;
        let (_, d1, d2) = curve.d2(t);
        let s = d1.magnitude();
        if s < 1e-12 {
            return None;
        }
        speeds.push(s);
        accs.push(d2.magnitude());
        binorms.push(d1.crossed(&d2).normalized());
    }
    let s0 = speeds[0];
    if speeds.iter().any(|&s| (s - s0).abs() > 1e-6 * s0.max(1.0)) {
        return None;
    }
    let a0 = accs[0];
    if a0 < 1e-12 || accs.iter().any(|&a| (a - a0).abs() > 1e-6 * a0.max(1.0)) {
        return None;
    }
    for bv in binorms.iter().skip(1) {
        if bv.cross_magnitude(&binorms[0]) > 1e-6 {
            return None;
        }
    }
    Some(s0)
}

/// Whether sampled points of the curve over `[a, b]` lie in a common plane.
pub(super) fn curve_is_planar(curve: &dyn occt_geom::Curve, a: f64, b: f64) -> bool {
    let n = 6;
    let pts: Vec<GpPnt> = (0..=n).map(|i| curve.d0(a + (b - a) * i as f64 / n as f64)).collect();
    if pts.len() < 3 {
        return false;
    }
    let v1 = GpVec::from_pnts(&pts[0], &pts[1]);
    let v2 = GpVec::from_pnts(&pts[0], &pts[pts.len() / 2]);
    let nrm = v1.crossed(&v2);
    let nm = nrm.magnitude();
    if nm < 1e-12 {
        return false;
    }
    pts.iter().all(|p| {
        GpVec::from_pnts(&pts[0], p).dot(&nrm).abs() <= 1e-6 * nm.max(1.0)
    })
}

/// Perimeter of a full ellipse with semi-axes `a` and `b` by the Ramanujan
/// approximation (the documented closed-form alternative to the arc-length
/// quadrature used for partial arcs):
/// `π·(3(a+b) − √((3a+b)·(a+3b)))`.
///
/// Exact for circles and accurate to roughly `1e-7` of the true perimeter for
/// any aspect ratio.
pub fn ellipse_length_ramanujan(a: f64, b: f64) -> f64 {
    let s = 3.0 * (a + b) - ((3.0 * a + b) * (a + 3.0 * b)).sqrt();
    PI * s
}

/// Exact length of one edge's curve.
///
/// * line — distance between the endpoints;
/// * circle — `radius · |b − a|` (the `GeomCircle` parameter is in radians);
/// * ellipse / other bounded planar curve — Gauss-Legendre quadrature of the
///   tangent speed (the elliptic integral of the second kind);
/// * anything else (B-spline, helix, …) — chord-length sampling via
///   [`crate::brep_measure::edge_length`].
pub(super) fn exact_edge_length(e: &Edge) -> f64 {
    let Some(curve) = BRepTool::edge_curve(e) else { return 0.0 };
    let (a, b) = BRepTool::edge_parameters(e);
    if !(a.is_finite() && b.is_finite() && b > a) {
        return 0.0;
    }
    match classify_curve(curve.as_ref(), a, b) {
        CurveKind::Line => curve.d0(a).distance(&curve.d0(b)),
        CurveKind::Circle => match circle_radius(curve.as_ref(), a, b) {
            Some(r) => r * (b - a),
            None => crate::brep_measure::edge_length(e, 32),
        },
        CurveKind::Ellipse => {
            let speed = |t: f64| curve.d1(t).1.magnitude();
            integrate(&speed, a, b, 24)
        }
        CurveKind::Other => crate::brep_measure::edge_length(e, 32),
    }
}

/// Exact arc length of a single edge: line/circle are computed in closed form,
/// an ellipse by Gauss-Legendre quadrature, and anything else by chord-length
/// sampling. See [`classify_curve`] for the classification rules.
pub fn edge_length_exact(e: &Edge) -> f64 {
    exact_edge_length(e)
}

/// Total length of all `shape`'s edges, exact for analytic curves (a box's 12
/// edges sum to `4·(dx+dy+dz)`, a cylinder's cap circles to `2πr` each, …).
pub fn analytic_curve_length(shape: &TopoShape) -> Result<f64, String> {
    let edges = edges_of(shape);
    if edges.is_empty() {
        return Err("analytic_curve_length: shape has no edges".into());
    }
    Ok(edges.iter().map(edge_length_exact).sum())
}

/// Total length of all `shape`'s edges; `0.0` when the shape has no edges
/// (unlike [`analytic_curve_length`], which errors).
pub fn analytic_perimeter(shape: &TopoShape) -> f64 {
    edges_of(shape).iter().map(edge_length_exact).sum()
}

/// Exact perimeter of a face's outer wire (sum of its edges' exact lengths).
pub fn face_perimeter_exact(face: &Face) -> f64 {
    let Some(w) = wires_of_face(face).first().cloned() else { return 0.0 };
    edges_of_wire(&w).iter().map(edge_length_exact).sum()
}

// ---------------------------------------------------------------------------
// Combined analytic mass properties
// ---------------------------------------------------------------------------

/// Full analytic mass properties of a solid: mass, volume, surface area,
/// centroid, inertia tensor, principal moments and principal axes.
#[derive(Debug, Clone)]
pub struct FullAnalyticProps {
    pub mass: f64,
    pub volume: f64,
    pub surface_area: f64,
    pub centroid: GpPnt,
    pub inertia: InertiaTensor,
    pub principal_moments: Vec<f64>,
    pub principal_axes: Vec<GpVec>,
    /// `true` when every face was handled analytically (no mesh fallback).
    pub exact: bool,
}

impl FullAnalyticProps {
    /// The inertia tensor about an arbitrary point `p`, shifted from the
    /// centroid by the parallel-axis theorem.
    pub fn inertia_about(&self, p: &GpPnt) -> InertiaTensor {
        let offset = GpVec::from_pnts(&self.centroid, p);
        self.inertia.translated(self.mass, &offset)
    }

    /// Radii of gyration about the principal axes: `sqrt(I/m)` per principal
    /// moment (descending order matches [`Self::principal_moments`]).
    pub fn radius_of_gyration(&self) -> Vec<f64> {
        let m = self.mass.max(1e-30);
        self.principal_moments.iter().map(|i| (i / m).sqrt()).collect()
    }
}

/// One-stop analytic mass properties: exact area/volume/centroid plus the
/// inertia tensor and principal moments/axes about the centroid.
pub fn full_analytic_properties(shape: &TopoShape, density: f64) -> Result<FullAnalyticProps, String> {
    let props = analytic_properties(shape)?;
    let inertia = inertia_tensor(shape, density)?;
    let (principal_moments, principal_axes) = principal_inertia(shape, density)?;
    Ok(FullAnalyticProps {
        mass: props.volume * density,
        volume: props.volume,
        surface_area: props.surface_area,
        centroid: props.centroid,
        inertia,
        principal_moments,
        principal_axes,
        exact: props.exact,
    })
}

// ---------------------------------------------------------------------------
// Analytic extrema
// ---------------------------------------------------------------------------

/// Exact closest point from `p` to an infinite line.
pub fn extrema_point_line(p: GpPnt, line: &GeomLine) -> ExtremaPair {
    let lin = line.lin();
    let p0 = lin.location();
    let d = *lin.direction().xyz();
    let v = GpVec::from_pnts(&p0, &p);
    let t = v.dot(&GpVec::from_xyz(&d));
    let q = GpPnt::from_xyz(&p0.coord.added(&d.multiplied(t)));
    ExtremaPair {
        p1: p,
        p2: q,
        distance: p.distance(&q),
        u1: t,
        v1: None,
        u2: t,
        v2: None,
    }
}

/// Exact closest point from `p` to a circle.
///
/// Projects `p` onto the circle's plane and solves the in-plane angle. When `p`
/// lies on the axis every circle point is equidistant; angle 0 is chosen.
pub fn extrema_point_circle(p: GpPnt, circle: &GeomCircle) -> ExtremaPair {
    let c = circle.circ();
    let center = c.location();
    let ax2 = c.position();
    let n = GpVec::from_xyz(ax2.direction().xyz());
    let xd = GpVec::from_xyz(ax2.x_direction().xyz());
    let yd = GpVec::from_xyz(ax2.y_direction().xyz());
    let r = c.radius();
    let rel = GpVec::from_pnts(&center, &p);
    let proj = rel.subtracted(&n.multiplied_scalar(rel.dot(&n)));
    let in_r = proj.magnitude();
    let angle = if in_r < 1e-12 {
        0.0
    } else {
        proj.dot(&yd).atan2(proj.dot(&xd))
    };
    let q = pnt_add_vec(
        center,
        &xd.multiplied_scalar(r * angle.cos()).add(&yd.multiplied_scalar(r * angle.sin())),
    );
    ExtremaPair {
        p1: p,
        p2: q,
        distance: p.distance(&q),
        u1: angle,
        v1: None,
        u2: angle,
        v2: None,
    }
}

/// Exact signed distance from `p` to a plane (closest point along the normal).
pub fn extrema_point_plane(p: GpPnt, plane: &GeomPlane) -> ExtremaPair {
    let pl = plane.pln();
    let loc = pl.location();
    let n = GpVec::from_xyz(pl.axis().direction().xyz());
    let xd = GpVec::from_xyz(pl.pos.x_direction().xyz());
    let yd = GpVec::from_xyz(pl.pos.y_direction().xyz());
    let rel = GpVec::from_pnts(&loc, &p);
    let d = rel.dot(&n);
    let q = GpPnt::from_xyz(&p.coord.subtracted(&n.coord.multiplied(d)));
    let u = rel.dot(&xd);
    let v = rel.dot(&yd);
    ExtremaPair {
        p1: p,
        p2: q,
        distance: d.abs(),
        u1: u,
        v1: Some(v),
        u2: u,
        v2: Some(v),
    }
}

/// Exact closest point from `p` to a sphere surface.
pub fn extrema_point_sphere(p: GpPnt, sphere: &GeomSphere) -> ExtremaPair {
    let center = sphere_center(sphere).unwrap_or_else(|| sphere.d0(0.0, 0.0));
    let r = sphere.d0(0.0, 0.0).distance(&center);
    let v = GpVec::from_pnts(&center, &p);
    let q = if v.magnitude() > 1e-12 {
        pnt_add_vec(center, &v.normalized().multiplied_scalar(r))
    } else {
        GpPnt::new(center.x() + r, center.y(), center.z())
    };
    let dist = (v.magnitude() - r).abs();
    let qrel = GpVec::from_pnts(&center, &q);
    let u = qrel.y().atan2(qrel.x());
    let vv = (qrel.z() / r.max(1e-30)).clamp(-1.0, 1.0).asin();
    ExtremaPair {
        p1: p,
        p2: q,
        distance: dist,
        u1: u,
        v1: Some(vv),
        u2: u,
        v2: Some(vv),
    }
}

/// Exact closest points between two (possibly skew) lines. Parallel lines
/// return any point pair at the perpendicular distance.
pub fn extrema_line_line(a: &GeomLine, b: &GeomLine) -> ExtremaPair {
    let la = a.lin();
    let lb = b.lin();
    let p1 = la.location();
    let p2 = lb.location();
    let d1 = GpVec::from_xyz(la.direction().xyz());
    let d2 = GpVec::from_xyz(lb.direction().xyz());
    let w0 = GpVec::from_pnts(&p2, &p1); // p1 − p2
    let aa = d1.dot(&d1);
    let bb = d1.dot(&d2);
    let cc = d2.dot(&d2);
    let dd = d1.dot(&w0);
    let ee = d2.dot(&w0);
    let denom = aa * cc - bb * bb;
    let (t, s) = if denom.abs() < 1e-12 {
        // parallel: any t; s projects p1 onto line b.
        (0.0, if cc > 1e-12 { ee / cc } else { 0.0 })
    } else {
        ((bb * ee - cc * dd) / denom, (aa * ee - bb * dd) / denom)
    };
    let q1 = GpPnt::from_xyz(&p1.coord.added(&d1.coord.multiplied(t)));
    let q2 = GpPnt::from_xyz(&p2.coord.added(&d2.coord.multiplied(s)));
    ExtremaPair {
        p1: q1,
        p2: q2,
        distance: q1.distance(&q2),
        u1: t,
        v1: None,
        u2: s,
        v2: None,
    }
}

/// Exact minimum distance between a line and a sphere:
/// `max(0, dist(line, center) − r)`. `None` when the sphere geometry is not
/// resolvable.
pub fn extrema_line_sphere(l: &GeomLine, s: &GeomSphere) -> Option<ExtremaPair> {
    let center = sphere_center(s)?;
    let r = s.d0(0.0, 0.0).distance(&center);
    let lin = l.lin();
    let p0 = lin.location();
    let d = GpVec::from_xyz(lin.direction().xyz());
    let v = GpVec::from_pnts(&p0, &center);
    let t = v.dot(&d);
    let q_line = pnt_add_vec(p0, &d.multiplied_scalar(t));
    let w = GpVec::from_pnts(&q_line, &center); // center − q_line
    let wm = w.magnitude();
    let q_sphere = if wm > 1e-12 {
        pnt_add_vec(center, &w.divided(wm).multiplied_scalar(-r))
    } else {
        GpPnt::new(center.x() + r, center.y(), center.z())
    };
    let dist = (wm - r).max(0.0);
    Some(ExtremaPair {
        p1: q_line,
        p2: q_sphere,
        distance: dist,
        u1: t,
        v1: None,
        u2: 0.0,
        v2: None,
    })
}

/// Exact nearest point on a cylinder surface to `p`.
pub(super) fn extrema_point_cylinder(s: &dyn Surface, p: GpPnt) -> Result<ExtremaPair, String> {
    let (center, ax, r) = cylinder_params(s).ok_or("extrema: cylinder params")?;
    let rel = GpVec::from_pnts(&center, &p);
    let rel_axis = rel.dot(&ax);
    let rel_plane = rel.subtracted(&ax.multiplied_scalar(rel_axis));
    let in_r = rel_plane.magnitude();
    let s0 = s.d0(0.0, 0.0);
    let xd_v = GpVec::from_pnts(&center, &s0);
    let xd = if xd_v.magnitude() > 1e-12 {
        xd_v.divided(xd_v.magnitude())
    } else {
        perpendicular(&ax)
    };
    let yd = ax.crossed(&xd).normalized();
    let angle = if in_r > 1e-12 {
        rel_plane.dot(&yd).atan2(rel_plane.dot(&xd))
    } else {
        0.0
    };
    let q = pnt_add_vec(
        center,
        &xd.multiplied_scalar(r * angle.cos())
            .add(&yd.multiplied_scalar(r * angle.sin()))
            .add(&ax.multiplied_scalar(rel_axis)),
    );
    let u = angle.rem_euclid(2.0 * PI);
    Ok(ExtremaPair {
        p1: p,
        p2: q,
        distance: p.distance(&q),
        u1: u,
        v1: Some(rel_axis),
        u2: u,
        v2: Some(rel_axis),
    })
}

/// Exact nearest point on a cone surface to `p`, solved in the axial plane.
pub(super) fn extrema_point_cone(s: &dyn Surface, p: GpPnt) -> Result<ExtremaPair, String> {
    let (apex, ax, alpha) = cone_params(s).ok_or("extrema: cone params")?;
    let rel = GpVec::from_pnts(&apex, &p);
    let z_p = rel.dot(&ax);
    let radial = rel.subtracted(&ax.multiplied_scalar(z_p));
    let r_p = radial.magnitude();
    let e_radial = if r_p > 1e-12 {
        radial.divided(r_p)
    } else {
        perpendicular(&ax)
    };
    let (sa, ca) = alpha.sin_cos();
    let t = r_p * sa + z_p * ca;
    let q = pnt_add_vec(
        apex,
        &ax.multiplied_scalar(t * ca).add(&e_radial.multiplied_scalar(t * sa)),
    );
    let u = if r_p > 1e-12 {
        let xd = e_radial;
        let yd = ax.crossed(&xd).normalized();
        let ux = e_radial.dot(&xd);
        let uy = e_radial.dot(&yd);
        uy.atan2(ux)
    } else {
        0.0
    };
    Ok(ExtremaPair {
        p1: p,
        p2: q,
        distance: p.distance(&q),
        u1: u,
        v1: Some(t),
        u2: u,
        v2: Some(t),
    })
}

/// Exact nearest point on a torus surface to `p`: project onto the center
/// circle, then onto the tube circle in the radial/axial plane.
pub(super) fn extrema_point_torus(s: &dyn Surface, p: GpPnt) -> Result<ExtremaPair, String> {
    let (center, ax, big_r, small_r) = torus_params(s).ok_or("extrema: torus params")?;
    let rel = GpVec::from_pnts(&center, &p);
    let rel_axis = rel.dot(&ax);
    let rel_plane = rel.subtracted(&ax.multiplied_scalar(rel_axis));
    let in_r = rel_plane.magnitude();
    let e_radial = if in_r > 1e-12 {
        rel_plane.divided(in_r)
    } else {
        perpendicular(&ax)
    };
    let q_center = pnt_add_vec(center, &e_radial.multiplied_scalar(big_r));
    let w = GpVec::from_pnts(&q_center, &p);
    let wm = w.magnitude();
    let q = if wm > 1e-12 {
        pnt_add_vec(q_center, &w.divided(wm).multiplied_scalar(small_r))
    } else {
        pnt_add_vec(q_center, &e_radial.multiplied_scalar(small_r))
    };
    let xd = {
        let s0 = s.d0(0.0, 0.0);
        let v = GpVec::from_pnts(&center, &s0);
        let vp = v.subtracted(&ax.multiplied_scalar(v.dot(&ax)));
        if vp.magnitude() > 1e-12 {
            vp.divided(vp.magnitude())
        } else {
            perpendicular(&ax)
        }
    };
    let yd = ax.crossed(&xd).normalized();
    let u = e_radial.dot(&yd).atan2(e_radial.dot(&xd));
    let w_dir = if wm > 1e-12 { w.divided(wm) } else { e_radial };
    let vv = w_dir.dot(&ax).atan2(w_dir.dot(&e_radial));
    Ok(ExtremaPair {
        p1: p,
        p2: q,
        distance: p.distance(&q),
        u1: u,
        v1: Some(vv),
        u2: u,
        v2: Some(vv),
    })
}

/// Exact nearest point from `p` to a surface, dispatching on the analytic
/// surface type. Unsupported surfaces fall back to the grid+descent
/// [`occt_geom::extrema::point_surface_extrema`].
pub fn extrema_point_surface_exact(s: &dyn Surface, p: GpPnt) -> Result<ExtremaPair, String> {
    match classify_surface_full(s) {
        SurfaceKind::Plane => {
            let origin = s.d0(0.0, 0.0);
            let n = surface_normal(s, 0.0, 0.0);
            if n.xyz().square_modulus() < 1e-30 {
                return Err("extrema: degenerate plane normal".into());
            }
            let z_axis = GpDir::new(0.0, 0.0, 1.0).map_err(|e| e.to_string())?;
            let xd_v = if n.is_parallel(&GpVec::from_xyz(z_axis.xyz())) {
                GpVec::new(1.0, 0.0, 0.0)
            } else {
                GpVec::from_xyz(z_axis.xyz())
            };
            let yd = n.crossed(&xd_v.normalized()).normalized();
            let rel = GpVec::from_pnts(&origin, &p);
            let d = rel.dot(&n);
            let q = GpPnt::from_xyz(&p.coord.subtracted(&n.coord.multiplied(d)));
            let u = rel.dot(&xd_v.normalized());
            let v = rel.dot(&yd);
            Ok(ExtremaPair {
                p1: p,
                p2: q,
                distance: d.abs(),
                u1: u,
                v1: Some(v),
                u2: u,
                v2: Some(v),
            })
        }
        SurfaceKind::Sphere => {
            let center = sphere_center(s).ok_or("extrema: sphere center")?;
            let r = s.d0(0.0, 0.0).distance(&center);
            let v = GpVec::from_pnts(&center, &p);
            let q = if v.magnitude() > 1e-12 {
                pnt_add_vec(center, &v.normalized().multiplied_scalar(r))
            } else {
                GpPnt::new(center.x() + r, center.y(), center.z())
            };
            let qrel = GpVec::from_pnts(&center, &q);
            let u = qrel.y().atan2(qrel.x());
            let vv = (qrel.z() / r.max(1e-30)).clamp(-1.0, 1.0).asin();
            Ok(ExtremaPair {
                p1: p,
                p2: q,
                distance: p.distance(&q),
                u1: u,
                v1: Some(vv),
                u2: u,
                v2: Some(vv),
            })
        }
        SurfaceKind::Cylinder => extrema_point_cylinder(s, p),
        SurfaceKind::Cone => extrema_point_cone(s, p),
        SurfaceKind::Torus => extrema_point_torus(s, p),
        SurfaceKind::Other => Ok(point_surface_extrema(s, &p)),
    }
}
