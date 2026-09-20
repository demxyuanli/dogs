use super::prelude::*;
use super::*;


/// Relation between two axes — port of `AxeOperator`
/// (IntAna_QuadQuadGeo.cxx). Returns `(parallel, distance, coplanar,
/// intersection point)`; `coplanar` needs the axes to (nearly) meet
/// (`distance < eps_dist` and the triple product within tolerance), and the
/// point is only computed for concurrent non-parallel axes.
pub(super) fn axe_operator(a1: &GpAx1, a2: &GpAx1, eps_dist: f64, eps_para: f64) -> (bool, f64, bool, Option<GpPnt>) {
    let mut v1 = *a1.direction();
    let mut v2 = *a2.direction();
    refine_dir(&mut v1);
    refine_dir(&mut v2);
    let (p1, p2) = (*a1.location(), *a2.location());
    let (w1, w2) = (GpVec::from_xyz(v1.xyz()), GpVec::from_xyz(v2.xyz()));
    let parallel = w1.cross_magnitude(&w2) <= eps_para;
    let distance = if parallel {
        let rel = GpVec::from_pnts(&p1, &p2);
        rel.subtracted(&w1.multiplied_scalar(rel.dot(&w1))).magnitude()
    } else {
        w1.crossed(&w2).normalized().dot(&GpVec::from_pnts(&p1, &p2)).abs()
    };
    let mut coplanar = false;
    let mut pt_intersect = None;
    if distance < eps_dist {
        let det = det33(&w1, &w2, &GpVec::from_pnts(&p2, &p1)); // rows V1, V2, P1−P2
        if det.abs() <= eps_dist {
            coplanar = true;
            if !parallel {
                // Concurrent axes: intersection point P1 + A·V1.
                let sm = GpVec::from_pnts(&p1, &p2);
                let d1 = w1.y() * w2.x() - w1.x() * w2.y();
                let d2 = w1.z() * w2.y() - w1.y() * w2.z();
                let d3 = w1.z() * w2.x() - w1.x() * w2.z();
                let a = if d1 != 0.0 && d1.abs() >= d2.abs() && d1.abs() >= d3.abs() {
                    (sm.y() * w2.x() - sm.x() * w2.y()) / d1
                } else if d2 != 0.0 && d2.abs() >= d1.abs() && d2.abs() >= d3.abs() {
                    (sm.z() * w2.y() - sm.y() * w2.z()) / d2
                } else {
                    (sm.z() * w2.x() - sm.x() * w2.z()) / d3
                };
                pt_intersect = Some(p1.translated_vec(&w1.multiplied_scalar(a)));
            }
        }
    }
    (parallel, distance, coplanar, pt_intersect)
}

/// The perpendicular segment between two non-parallel axes: signed distance and
/// the parameters of the closest points on each axis. Port of
/// `AxeOperator::Distance`.
pub(super) fn axe_distance(a1: &GpAx1, a2: &GpAx1) -> (f64, f64, f64) {
    let w1 = GpVec::from_xyz(a1.direction().xyz());
    let w2 = GpVec::from_xyz(a2.direction().xyz());
    let o1o2 = GpVec::from_pnts(a1.location(), a2.location());
    let n = w1.crossed(&w2);
    if n.magnitude() < 1e-12 {
        return (0.0, 0.0, 0.0);
    }
    let n = n.normalized();
    let d = det33(&w1, &w2, &n);
    if d != 0.0 {
        let dist = det33(&w1, &w2, &o1o2) / d;
        let p1 = det33(&o1o2, &w2, &n) / (-d);
        let p2 = det33(&w1, &o1o2, &n) / d;
        (dist, p1, p2)
    } else {
        (0.0, 0.0, 0.0)
    }
}

/// Distance from `p` to the axis line.
pub(super) fn dist_point_axis(p: &GpPnt, ax: &GpAx1) -> f64 {
    let dir = GpVec::from_xyz(ax.direction().xyz());
    let rel = GpVec::from_pnts(ax.location(), p);
    rel.subtracted(&dir.multiplied_scalar(rel.dot(&dir))).magnitude()
}

pub(super) fn mid_pnt(a: &GpPnt, b: &GpPnt) -> GpPnt {
    GpPnt::new(0.5 * (a.x() + b.x()), 0.5 * (a.y() + b.y()), 0.5 * (a.z() + b.z()))
}

/// A circle with the given center / plane normal / radius (the frame X
/// direction is arbitrary — it does not affect the curve).
pub(super) fn circle_with_normal(center: GpPnt, normal: GpDir, radius: f64) -> GpCirc {
    GpCirc::new(ax2_from_dirs(center, normal, perp_x_dir(&normal)), radius)
}

/// A plane through `pt` with the given normal.
pub(super) fn plane_normal_at(pt: GpPnt, normal: GpDir) -> Option<GpPln> {
    let ax3 = GpAx3::new(pt, normal, &perp_x_dir(&normal)).ok()?;
    Some(GpPln::new(ax3))
}

/// Cylinder ∩ cylinder. Port of `IntAna_QuadQuadGeo::Perform(gp_Cylinder,
/// gp_Cylinder)`: parallel axes give two generatrix lines (or one when tangent,
/// `Same`/`None` for coincident/disjoint), intersecting equal-radius axes give
/// the two bisector-plane ellipses, external tangency a point, and everything
/// else (`NoGeometricSolution`) → `None` for the numeric walker.
pub fn quadric_quadric_cylinder_cylinder(
    c1: &GpCylinder,
    c2: &GpCylinder,
    tol: f64,
) -> QuadricIntersection {
    let (parallel, dist, coplanar, pt_inter) =
        axe_operator(&c1.axis(), &c2.axis(), CONFUSION, ANGULAR);
    let r1 = c1.radius();
    let r2 = c2.radius();
    let rmr = (r1 - r2).abs();
    let rmr_relative = rmr / r1.max(r2);
    let dir_cyl = c1.position().direction();
    let wdir = GpVec::from_xyz(dir_cyl.xyz());

    if parallel {
        if dist <= tol {
            return if rmr <= tol {
                QuadricIntersection::Same
            } else {
                QuadricIntersection::None
            };
        }
        // Parallel axes, strictly separated. Project the 2nd location onto the
        // 1st cylinder base plane and intersect the two base circles.
        let p1 = c1.location();
        let p2t = c2.location();
        let proj = wdir.dot(&GpVec::from_pnts(&p1, &p2t));
        let p2 = p2t.translated_vec(&wdir.multiplied_scalar(-proj));
        let r1p2 = r1 + r2;
        if dist > r1p2 + tol {
            QuadricIntersection::None
        } else if (r1p2 - dist) <= f64::EPSILON {
            // External tangency: one generatrix line.
            let pt1 = p1.translated_vec(&GpVec::from_pnts(&p1, &p2).multiplied_scalar(r1 / r1p2));
            QuadricIntersection::Line(GpLin::from_pnt_dir(pt1, dir_cyl))
        } else if dist > rmr {
            // Two generatrix lines (or one when the base circles are tangent).
            let a_r1r1 = r1 * r1;
            let a_cos = 0.5 * (a_r1r1 - r2 * r2 + dist * dist) / (r1 * dist);
            let a_sin2 = 1.0 - a_cos * a_cos;
            let is_tangent = 4.0 * a_r1r1 * a_sin2 < tol * tol;
            let dir_a1a2 = GpVec::from_pnts(&p1, &p2).divided(dist);
            if is_tangent {
                let pt1 = p1.translated_vec(&dir_a1a2.multiplied_scalar(r1 * a_cos));
                QuadricIntersection::Line(GpLin::from_pnt_dir(pt1, dir_cyl))
            } else {
                let a_sin = a_sin2.sqrt();
                let axd = *c1.position().x_direction();
                let ayd = *c1.position().y_direction();
                let r1x = GpVec::from_xyz(axd.xyz()).multiplied_scalar(r1);
                let r1y = GpVec::from_xyz(ayd.xyz()).multiplied_scalar(r1);
                let adx = dir_a1a2.dot(&GpVec::from_xyz(axd.xyz()));
                let ady = dir_a1a2.dot(&GpVec::from_xyz(ayd.xyz()));
                let (ndx, ndy) = (adx * a_cos - ady * a_sin, ady * a_cos + adx * a_sin);
                let pt1 = p1
                    .translated_vec(&r1x.multiplied_scalar(ndx))
                    .translated_vec(&r1y.multiplied_scalar(ndy));
                let (ndx, ndy) = (adx * a_cos + ady * a_sin, ady * a_cos - adx * a_sin);
                let pt2 = p1
                    .translated_vec(&r1x.multiplied_scalar(ndx))
                    .translated_vec(&r1y.multiplied_scalar(ndy));
                QuadricIntersection::TwoLines(
                    GpLin::from_pnt_dir(pt1, dir_cyl),
                    GpLin::from_pnt_dir(pt2, dir_cyl),
                )
            }
        } else if dist > rmr - tol {
            // Internal tangency: one generatrix line.
            let mut r1_rmr = r1 / rmr;
            if r1 < r2 {
                r1_rmr = -r1_rmr;
            }
            let pt1 = p1.translated_vec(&GpVec::from_pnts(&p1, &p2).multiplied_scalar(r1_rmr));
            QuadricIntersection::Line(GpLin::from_pnt_dir(pt1, dir_cyl))
        } else {
            QuadricIntersection::None
        }
    } else if rmr_relative <= 1e-13 && coplanar {
        // Equal-radius cylinders with intersecting axes → two ellipses in the
        // bisector planes. Frame: `dir1`/`dir2` are the bisectors (perpendicular),
        // each used as the plane normal of one ellipse.
        let Some(pt) = pt_inter else { return QuadricIntersection::None };
        let wd2 = GpVec::from_xyz(c2.position().direction().xyz());
        let ang = wdir.angle(&wd2);
        let b = (0.5 * (PI - ang)).sin().abs();
        let a = (0.5 * ang).sin().abs();
        if a == 0.0 || b == 0.0 {
            return QuadricIntersection::Same;
        }
        let Ok(d1) = GpDir::from_vec(&wdir.added(&wd2)) else { return QuadricIntersection::None };
        let Ok(d2) = GpDir::from_vec(&wdir.subtracted(&wd2)) else { return QuadricIntersection::None };
        let (mut p1, mut p1bis) = (r1 / b, r1);
        let (mut p2, mut p2bis) = (r1 / a, r1);
        if p1 < p1bis {
            std::mem::swap(&mut p1, &mut p1bis);
        }
        if p2 < p2bis {
            std::mem::swap(&mut p2, &mut p2bis);
        }
        let Ok(ax2_1) = GpAx2::new(pt, d1, d2) else { return QuadricIntersection::None };
        let Ok(ax2_2) = GpAx2::new(pt, d2, d1) else { return QuadricIntersection::None };
        QuadricIntersection::TwoEllipses(
            GpElips::new(ax2_1, p1, p1bis),
            GpElips::new(ax2_2, p2, p2bis),
        )
    } else if (dist - r1 - r2).abs() < tol {
        // External tangency with intersecting (non-parallel) axes: a point on
        // the common perpendicular.
        let d1 = *c1.axis().direction();
        let d2 = *c2.axis().direction();
        let (_, p1p, p2p) = axe_distance(&c1.axis(), &c2.axis());
        let p1 = c1
            .axis()
            .location()
            .translated_vec(&GpVec::from_xyz(d1.xyz()).multiplied_scalar(-p1p));
        let p2 = c2
            .axis()
            .location()
            .translated_vec(&GpVec::from_xyz(d2.xyz()).multiplied_scalar(-p2p));
        let Ok(dir) = GpDir::from_vec(&GpVec::from_pnts(&p1, &p2)) else { return QuadricIntersection::None };
        let pt = p1.translated_vec(&GpVec::from_xyz(dir.xyz()).multiplied_scalar(r1));
        QuadricIntersection::Point(pt)
    } else {
        QuadricIntersection::None
    }
}

/// Cylinder ∩ sphere. Port of `IntAna_QuadQuadGeo::Perform(gp_Cylinder,
/// gp_Sphere)`: when the sphere center lies on the cylinder axis the section is
/// one or two circles (radius = cylinder radius, centered at
/// `center ± √(r_sph² − r_cyl²)·axis`); otherwise `NoGeometricSolution` → `None`.
pub fn quadric_quadric_cylinder_sphere(
    cyl: &GpCylinder,
    sph: &GpSphere,
    _tol: f64,
) -> QuadricIntersection {
    let pt = sph.location();
    // OCCT tests the axes to intersect at the sphere center exactly; a small
    // tolerance keeps the closed form on near-axis configurations.
    if dist_point_axis(&pt, &cyl.axis()) > 1e-9 {
        return QuadricIntersection::None;
    }
    let r_cyl = cyl.radius();
    let r_sph = sph.radius();
    if r_sph < r_cyl {
        return QuadricIntersection::None; // IntAna_Empty
    }
    let dist = (r_sph * r_sph - r_cyl * r_cyl).sqrt();
    let dir = cyl.position().direction();
    let w = GpVec::from_xyz(dir.xyz());
    let c1 = pt.translated_vec(&w.multiplied_scalar(dist));
    let circ1 = circle_with_normal(c1, dir, r_cyl);
    if dist > f64::EPSILON {
        let c2 = pt.translated_vec(&w.multiplied_scalar(-dist));
        QuadricIntersection::TwoCircles(circ1, circle_with_normal(c2, dir, r_cyl))
    } else {
        QuadricIntersection::Circle(circ1)
    }
}

/// Sphere ∩ cone. Port of `IntAna_QuadQuadGeo::Perform(gp_Sphere, gp_Cone)`:
/// when the sphere center lies on the cone axis the section is one or two
/// circles — the roots of the 2D cross-section quadratic
/// `(1+tg²)x² + 2·tg²·d·x + tg²·d² − r² = 0`, with `d` the apex→center
/// distance. Otherwise `NoGeometricSolution` → `None`.
pub fn quadric_quadric_sphere_cone(
    sph: &GpSphere,
    cone: &GpCone,
    _tol: f64,
) -> QuadricIntersection {
    let pt = sph.location();
    if dist_point_axis(&pt, &cone.axis()) > 1e-9 {
        return QuadricIntersection::None;
    }
    let apex = cone.apex();
    let d = pt.distance(&apex);
    let condir = if d > f64::EPSILON {
        let Ok(c) = GpDir::from_vec(&GpVec::from_pnts(&apex, &pt)) else {
            return QuadricIntersection::None;
        };
        c
    } else {
        cone.position().direction()
    };
    let rad = sph.radius();
    let tga = cone.semi_angle().tan();
    let tgatga = tga * tga;
    let roots = quadratic_roots(1.0 + tgatga, 2.0 * tgatga * d, -rad * rad + d * d * tgatga);
    if roots.is_empty() {
        return QuadricIntersection::None; // IntAna_Empty
    }
    let w = GpVec::from_xyz(condir.xyz());
    let mut circles: Vec<GpCirc> = Vec::new();
    for x in roots {
        let dpx = d + x;
        let center = apex.translated_vec(&w.multiplied_scalar(dpx));
        let r = (tga * dpx).abs();
        if r <= 0.01 * CONFUSION {
            continue; // IntAna_PointAndCircle: degenerate radius → point
        }
        circles.push(circle_with_normal(center, condir, r));
    }
    match circles.len() {
        1 => QuadricIntersection::Circle(circles.pop().unwrap()),
        2 => {
            let c2 = circles.pop().unwrap();
            let c1 = circles.pop().unwrap();
            QuadricIntersection::TwoCircles(c1, c2)
        }
        _ => QuadricIntersection::None,
    }
}

/// Cone ∩ cone. Port of `IntAna_QuadQuadGeo::Perform(gp_Cone, gp_Cone)` for the
/// tractable branches: coincident axes (two circles / a point / `Same`),
/// parallel axes with equal semi-angle (a conic in the plane through the two
/// apexes, from the plane∩cone closed form), and coincident apexes (one or two
/// generatrix lines). The common-generatrix case and everything else return
/// `None` for the numeric walker.
pub fn quadric_quadric_cone_cone(
    c1: &GpCone,
    c2: &GpCone,
    tol_ang: f64,
    tol: f64,
) -> QuadricIntersection {
    let tg1 = c1.semi_angle().tan();
    let mut tg2 = c2.semi_angle().tan();
    if tg1 * tg2 < 0.0 {
        tg2 = -tg2;
    }
    let tol2 = tol * tol;
    let ap1 = c1.apex();
    let ap2 = c2.apex();
    let d_a1a2 = ap1.square_distance(&ap2);
    let (parallel, dist_axes, _coplanar, _pt_inter) =
        axe_operator(&c1.axis(), &c2.axis(), 1e-14, ANGULAR);

    // 1 — coincident axes: two circles where the cone radii match (or the two
    // cones coincide / touch at the apex).
    if parallel && dist_axes < 1e-14 {
        let p = c1.apex();
        let d = c1.position().direction();
        let w = GpVec::from_xyz(d.xyz());
        let offset = w.dot(&GpVec::from_pnts(&p, &ap2));
        if (tg1 - tg2).abs() > ANGULAR {
            if offset.abs() < 1e-10 {
                return QuadricIntersection::Point(p);
            }
            let x1 = offset * tg2 / (tg1 + tg2);
            let x2 = offset * tg2 / (tg2 - tg1);
            let c1c = p.translated_vec(&w.multiplied_scalar(x1));
            let c2c = p.translated_vec(&w.multiplied_scalar(x2));
            QuadricIntersection::TwoCircles(
                circle_with_normal(c1c, d, (x1 * tg1).abs()),
                circle_with_normal(c2c, d, (x2 * tg1).abs()),
            )
        } else if offset.abs() < 1e-10 {
            QuadricIntersection::Same
        } else {
            let x = 0.5 * offset;
            QuadricIntersection::Circle(circle_with_normal(
                p.translated_vec(&w.multiplied_scalar(x)),
                d,
                (x * tg1).abs(),
            ))
        }
    }
    // 2 — parallel axes with (nearly) equal semi-angles: the intersection lies
    // in the plane through the two apexes; reduce to the plane∩cone conic.
    else if (tg1 - tg2).abs() < tol_ang && parallel {
        let da1 = c1.position().direction();
        let o1o2 = GpVec::from_pnts(&ap1, &ap2);
        let o1o2n = o1o2.normalized();
        let o1o2_da1 = GpVec::from_xyz(da1.xyz()).dot(&o1o2n);
        let o1_proj = o1o2n.subtracted(&GpVec::from_xyz(da1.xyz()).multiplied_scalar(o1o2_da1));
        let Ok(db1) = GpDir::from_vec(&o1_proj) else { return QuadricIntersection::None };
        let y_o1o2 = o1o2.dot(&GpVec::from_xyz(da1.xyz()));
        let abstg1 = tg1.abs();
        let x2 = (dist_axes / abstg1 - y_o1o2) * 0.5;
        let x1 = x2 + y_o1o2;
        let p1 = ap1
            .translated_vec(&GpVec::from_xyz(da1.xyz()).multiplied_scalar(x1))
            .translated_vec(&GpVec::from_xyz(db1.xyz()).multiplied_scalar(x1 * abstg1));
        let p1_m = GpVec::from_pnts(&p1, &mid_pnt(&ap1, &ap2));
        let da1_x_db1 = GpVec::from_xyz(da1.xyz()).crossed(&GpVec::from_xyz(db1.xyz()));
        let ortho = da1_x_db1.crossed(&p1_m);
        let Ok(n) = GpDir::from_vec(&ortho) else { return QuadricIntersection::None };
        let Some(pln) = plane_normal_at(p1, n) else { return QuadricIntersection::None };
        quadric_quadric_plane_cone(&pln, c1, tol_ang, tol)
    }
    // 3 — coincident apexes: one or two generatrix lines (or `None`).
    else if d_a1a2 < tol2 {
        cone_cone_common_apex(c1, c2, tg1, tg2, tol)
    } else {
        // 4/5 — common generatrix / general: no analytic closed form.
        QuadricIntersection::None
    }
}

/// Coincident-apex cones: `IntAna_QuadQuadGeo::Perform` branch 3 — a 2D
/// section analysis determines touch/intersection, then the one or two
/// generatrix lines are built through the shared apex.
pub(super) fn cone_cone_common_apex(c1: &GpCone, c2: &GpCone, tg1: f64, tg2: f64, tol: f64) -> QuadricIntersection {
    let half_pi = 0.5 * PI;
    let d1 = 1.0;
    let p0 = GpPnt2d::new(0.0, 0.0);
    let ax1 = c1.axis();
    let ax2 = c2.axis();
    let mut gamma = ax1.direction().angle(ax2.direction());
    if gamma > half_pi {
        gamma = PI - gamma;
    }
    let (cos_g, sin_g) = (gamma.cos(), gamma.sin());
    let tg_beta1 = tg1.abs();
    let tg_beta2 = tg2.abs();
    let r1 = d1 * tg_beta1;
    let p1 = GpPnt2d::new(d1, r1);
    // Project P1 onto the 2nd axis line (in the plane of the two axes) to find
    // whether the section circles overlap, touch or miss.
    let v_ax2 = GpVec2d::new(cos_g, sin_g);
    let Ok(_) = GpDir2d::from_vec2d(&v_ax2) else { return QuadricIntersection::None };
    let v = GpVec2d::new(p1.x() - p0.x(), p1.y() - p0.y());
    let mut dx = v_ax2.dot(&v);
    let pa2 = p0.translated_vec(&v_ax2.multiplied_scalar(dx));
    dx = pa2.distance(&p0);
    let r2 = dx * tg_beta2;
    let rd2 = pa2.distance(&p1);
    if rd2 > r2 + tol {
        return QuadricIntersection::None; // IntAna_Empty
    }
    let i_ret = if rd2 < r2 - tol { 2 } else { 1 };
    // 3D construction: two planes perpendicular to the axes through the ring
    // points Q1/Q2 intersect in the line through the section mid-point QX.
    let q_apex1 = c1.apex();
    let d3_ax1 = *ax1.direction();
    let w1 = GpVec::from_xyz(d3_ax1.xyz());
    let qa1 = q_apex1.translated_vec(&w1.multiplied_scalar(d1));
    let dx = w1.dot(&GpVec::from_xyz(ax2.direction().xyz()));
    let d3_ax2 = if dx < 0.0 {
        ax2.direction().reversed()
    } else {
        *ax2.direction()
    };
    let w2 = GpVec::from_xyz(d3_ax2.xyz());
    let d2 = d1 * ((1.0 + tg_beta1 * tg_beta1) / (1.0 + tg_beta2 * tg_beta2)).sqrt();
    let qa2 = q_apex1.translated_vec(&w2.multiplied_scalar(d2));
    let Some(pln1) = plane_normal_at(qa1, d3_ax1) else { return QuadricIntersection::None };
    let Some(pln2) = plane_normal_at(qa2, d3_ax2) else { return QuadricIntersection::None };
    let Some(lin) = plane_plane_line(&pln1, &pln2) else { return QuadricIntersection::None };
    let wl = GpVec::from_xyz(lin.direction().xyz());
    let orig = lin.location();
    let vr = GpVec::from_pnts(&qa1, &orig);
    let dx = wl.dot(&vr);
    let qx = orig.translated_vec(&wl.multiplied_scalar(dx));
    if i_ret == 1 {
        // One tangency line.
        let Ok(dir) = GpDir::from_vec(&GpVec::from_pnts(&q_apex1, &qx)) else {
            return QuadricIntersection::None;
        };
        QuadricIntersection::Line(GpLin::from_pnt_dir(q_apex1, dir))
    } else {
        // Two intersection lines.
        let da = qa1.distance(&qx);
        let ddx = (r1 * r1 - da * da).sqrt();
        let qx1 = qx.translated_vec(&wl.multiplied_scalar(ddx));
        let qx2 = qx.translated_vec(&wl.multiplied_scalar(-ddx));
        let Ok(dir1) = GpDir::from_vec(&GpVec::from_pnts(&q_apex1, &qx1)) else {
            return QuadricIntersection::None;
        };
        let Ok(dir2) = GpDir::from_vec(&GpVec::from_pnts(&q_apex1, &qx2)) else {
            return QuadricIntersection::None;
        };
        QuadricIntersection::TwoLines(
            GpLin::from_pnt_dir(q_apex1, dir1),
            GpLin::from_pnt_dir(q_apex1, dir2),
        )
    }
}

/// Cylinder ∩ cone. Port of `IntAna_QuadQuadGeo::Perform(gp_Cylinder,
/// gp_Cone)`: only the coincident-axis case has a closed form (two circles at
/// `apex ± r_cyl/tan(angle)` along the axis); otherwise `NoGeometricSolution`.
pub fn quadric_quadric_cylinder_cone(
    cyl: &GpCylinder,
    cone: &GpCone,
    _tol: f64,
) -> QuadricIntersection {
    let (parallel, dist, _, _) = axe_operator(&cyl.axis(), &cone.axis(), 1e-14, ANGULAR);
    if !(parallel && dist < 1e-14) {
        return QuadricIntersection::None;
    }
    let pt = cone.apex();
    let dist = cyl.radius() / cone.semi_angle().tan();
    let dir = cyl.position().direction();
    let w = GpVec::from_xyz(dir.xyz());
    let r = cyl.radius();
    QuadricIntersection::TwoCircles(
        circle_with_normal(pt.translated_vec(&w.multiplied_scalar(dist)), dir, r),
        circle_with_normal(pt.translated_vec(&w.multiplied_scalar(-dist)), dir, r),
    )
}

/// Analytic quadric-quadric intersection dispatcher for the tractable pairs.
pub fn quadric_quadric(q1: &Quadric, q2: &Quadric, tol_ang: f64, tol: f64) -> QuadricIntersection {
    use Quadric::*;
    match (q1, q2) {
        (Plane(a), Plane(b)) => quadric_quadric_planes(a, b, tol_ang, tol),
        (Plane(p), Sphere(s)) | (Sphere(s), Plane(p)) => quadric_quadric_plane_sphere(p, s),
        (Sphere(a), Sphere(b)) => quadric_quadric_sphere_sphere(a, b, tol),
        (Plane(p), Cylinder(c)) | (Cylinder(c), Plane(p)) => {
            quadric_quadric_plane_cylinder(p, c, tol_ang, tol)
        }
        (Plane(p), Cone(c)) | (Cone(c), Plane(p)) => quadric_quadric_plane_cone(p, c, tol_ang, tol),
        (Cylinder(a), Cylinder(b)) => quadric_quadric_cylinder_cylinder(a, b, tol),
        (Cylinder(c), Sphere(s)) | (Sphere(s), Cylinder(c)) => quadric_quadric_cylinder_sphere(c, s, tol),
        (Sphere(s), Cone(c)) | (Cone(c), Sphere(s)) => quadric_quadric_sphere_cone(s, c, tol),
        (Cone(a), Cone(b)) => quadric_quadric_cone_cone(a, b, tol_ang, tol),
        (Cylinder(c), Cone(k)) | (Cone(k), Cylinder(c)) => quadric_quadric_cylinder_cone(c, k, tol),
    }
}

// ---------------------------------------------------------------------------
// IntAna_IntLinTorus
// ---------------------------------------------------------------------------

/// Intersection of a line with a torus. Port of `IntAna_IntLinTorus::Perform`:
/// a quartic in the line parameter in the torus reference frame, with each
/// root verified by re-evaluation on the torus.
pub fn line_torus_intersect(l: &GpLin, t: &GpTorus) -> Vec<GpPnt> {
    let pl = l.location();
    let dl = l.direction();
    let tor_loc = t.location();
    // Reparametrize so the line location is nearest the torus location.
    let param_of_new_pl = GpVec::from_pnts(&pl, &tor_loc).dot(&GpVec::from_xyz(dl.xyz()));
    let new_pl = pl.translated_vec(&GpVec::from_xyz(dl.xyz()).multiplied_scalar(param_of_new_pl));

    // Express the line in the torus reference frame.
    let pos = t.position();
    let (xd, yd, zd) = (
        GpVec::from_xyz(pos.x_direction().xyz()),
        GpVec::from_xyz(pos.y_direction().xyz()),
        GpVec::from_xyz(pos.direction().xyz()),
    );
    let v = GpVec::from_pnts(&pos.location(), &new_pl);
    let (x0, y0, z0) = (v.dot(&xd), v.dot(&yd), v.dot(&zd));
    let dv = GpVec::from_xyz(dl.xyz());
    let (x1, y1, z1) = (dv.dot(&xd), dv.dot(&yd), dv.dot(&zd));

    let r = t.major_radius();
    let r2 = r * r;
    let rr = t.minor_radius();
    let rr2 = rr * rr;

    let a = x1 * x1 + y1 * y1 + z1 * z1;
    let b = 2.0 * (x1 * x0 + y1 * y0 + z1 * z0);
    let c = x0 * x0 + y0 * y0 + z0 * z0 - (r2 + rr2);

    let a4 = a * a;
    let a3 = 2.0 * a * b;
    let a2 = 2.0 * a * c + 4.0 * r2 * z1 * z1 + b * b;
    let a1 = 2.0 * b * c + 8.0 * r2 * z1 * z0;
    let a0 = c * c + 4.0 * r2 * (z0 * z0 - rr2);

    let mut out = Vec::new();
    for mut tt in quartic_roots(a4, a3, a2, a1, a0) {
        tt += param_of_new_pl;
        let p = clib::line_value(l, tt);
        // `IntAna_IntLinTorus.cxx:98-106`: re-parameterise the candidate on the
        // torus (`ElSLib::Parameters`) and re-evaluate it (`ElSLib::Value`); a
        // root is kept only when the **square distance** between the line point
        // and the torus point is `<= 1e-10` (OCCT counts the others as bad
        // solutions). The previous body tested an invented implicit-equation
        // residual `|(rho-R)^2 + z^2 - r^2| < 1e-7` and deduplicated points with
        // a further `1e-7` (audit A15) — OCCT stores every valid root as is.
        let (u, v) = slib::torus_parameters(&t.position(), r, rr, &p);
        let p_sol_t = slib::torus_value(t, u, v);
        if p_sol_t.square_distance(&p) > 1.0e-10 {
            continue;
        }
        out.push(p);
    }
    out.sort_by(|a, b| a.x().partial_cmp(&b.x()).unwrap_or(Ordering::Equal));
    out
}
