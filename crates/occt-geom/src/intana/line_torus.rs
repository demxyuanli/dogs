use super::prelude::*;
use super::*;

use occt_core::math_direct_poly_roots::DirectPolynomialRoots;

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

/// `(A1A2.Intersect() && Pt.Distance(A1A2.PtIntersect()) == 0.0) || A1A2.Same()`
/// (`IntAna_QuadQuadGeo.cxx:1378-1379` for cylinder/sphere, `cxx:1924-1925` for
/// sphere/cone): the sphere center must lie exactly on the point where the two
/// axes meet, or the axes must be parallel and coincident within
/// `myEPSILON_DISTANCE`. `AxeOperator`'s defaults are `theEpsDistance = 1.e-14`
/// and `theEpsAxesPara = Precision::Angular()` (`cxx:72-73`). OCCT adds no
/// tolerance of its own, so neither does this predicate.
pub fn axes_meet_at_point(a1: &GpAx1, a2: &GpAx1, pt: &GpPnt) -> bool {
    let (parallel, distance, coplanar, pt_intersect) = axe_operator(a1, a2, 1.0e-14, ANGULAR);
    // `A1A2.Intersect()` is `thecoplanar && !theparallel` (`cxx:86`).
    if coplanar && !parallel {
        if let Some(pi) = pt_intersect {
            if pt.distance(&pi) == 0.0 {
                return true;
            }
        }
    }
    // `A1A2.Same()` is `theparallel && (thedistance < myEPSILON_DISTANCE)` (`cxx:83`).
    parallel && distance < 1.0e-14
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
        } else if (r1p2 - dist) <= REAL_SMALL {
            // External tangency: one generatrix line.
            let pt1 = p1.translated_vec(&GpVec::from_pnts(&p1, &p2).multiplied_scalar(r1 / r1p2));
            QuadricIntersection::Line(GpLin::from_pnt_dir(pt1, dir_cyl))
        } else if dist > rmr {
            // Two generatrix lines (or one when the base circles are tangent).
            let r1r1 = r1 * r1;
            let cos = 0.5 * (r1r1 - r2 * r2 + dist * dist) / (r1 * dist);
            let sin2 = 1.0 - cos * cos;
            let is_tangent = 4.0 * r1r1 * sin2 < tol * tol;
            let dir_a1a2 = GpVec::from_pnts(&p1, &p2).divided(dist);
            if is_tangent {
                let pt1 = p1.translated_vec(&dir_a1a2.multiplied_scalar(r1 * cos));
                QuadricIntersection::Line(GpLin::from_pnt_dir(pt1, dir_cyl))
            } else {
                let sin_val = sin2.sqrt();
                let axd = *c1.position().x_direction();
                let ayd = *c1.position().y_direction();
                let r1x = GpVec::from_xyz(axd.xyz()).multiplied_scalar(r1);
                let r1y = GpVec::from_xyz(ayd.xyz()).multiplied_scalar(r1);
                let adx = dir_a1a2.dot(&GpVec::from_xyz(axd.xyz()));
                let ady = dir_a1a2.dot(&GpVec::from_xyz(ayd.xyz()));
                let (ndx, ndy) = (adx * cos - ady * sin_val, ady * cos + adx * sin_val);
                let pt1 = p1
                    .translated_vec(&r1x.multiplied_scalar(ndx))
                    .translated_vec(&r1y.multiplied_scalar(ndy));
                let (ndx, ndy) = (adx * cos + ady * sin_val, ady * cos - adx * sin_val);
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
        // `cxx:1224` guards this arm with `RmR_Relative <=
        // myEPSILON_CYLINDER_DELTA_RADIUS` (`InitTolerances`, `cxx:356`) and
        // `A1A2.Intersect()`; in the non-parallel branch that is exactly
        // `coplanar`.
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
    // `cxx:1378-1379`: the axes must meet exactly at the sphere center, or be
    // coincident parallel axes. The exact `AxeOperator` predicate is the OCCT
    // control flow -- no extra tolerance is added here.
    if !axes_meet_at_point(&cyl.axis(), sph.position().axis(), &pt) {
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
    // `cxx:1924-1925`: same exact-axis predicate as the cylinder/sphere arm.
    if !axes_meet_at_point(&cone.axis(), sph.position().axis(), &pt) {
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

/// `gce_MakePln(P1, P2, P3)` (`gce_MakePln.cxx:76-92`): the plane through three
/// points, X direction along `P1P2`. `None` for collinear points
/// (`gce_ColinearPoints`).
fn plane_through_three(p1: &GpPnt, p2: &GpPnt, p3: &GpPnt) -> Option<GpPln> {
    let v1 = GpVec::from_pnts(p1, p2);
    let v2 = GpVec::from_pnts(p1, p3);
    let normal = v1.crossed(&v2);
    if normal.magnitude() < REAL_SMALL {
        return None;
    }
    let n = GpDir::from_vec(&normal).ok()?;
    let x = GpDir::from_vec(&v1).ok()?;
    let ax3 = GpAx3::new(*p1, n, &x).ok()?;
    Some(GpPln::new(ax3))
}

/// `EstimDist(Con1, Con2)` (`IntAna_QuadQuadGeo.cxx:260-326`): the minimal
/// distance from an apex to any solution, from the four pairwise intersections
/// of the two cones' generatrix lines projected into the plane through both
/// apexes and the direction of the first cone. `Precision::Infinite()` when
/// that plane cannot be built.
fn estim_dist(c1: &GpCone, c2: &GpCone) -> f64 {
    let ap1 = c1.apex();
    let ap2 = c2.apex();
    let d1 = c1.position().direction();
    let ap3 = ap1.translated_vec(&GpVec::from_xyz(d1.xyz()));
    let Some(pln) = plane_through_three(&ap1, &ap2, &ap3) else {
        return INFINITE;
    };
    // `cxx:279-291`: the two axes projected, then each rotated by ±SemiAngle.
    // OCCT's `gp_Dir2d` would raise `Standard_ConstructionError` where the
    // projection vanishes; the axes are inside the plane here, so `None` is
    // unreachable for the parallel-axes caller.
    let lin1 = GpLin::from_pnt_dir(ap1, d1);
    let lin2 = GpLin::from_pnt_dir(ap2, c2.position().direction());
    let Some(ax1_2d) = project_pln_lin2d(&pln, &lin1) else {
        return INFINITE;
    };
    let Some(ax2_2d) = project_pln_lin2d(&pln, &lin2) else {
        return INFINITE;
    };
    let an_ang1 = c1.semi_angle();
    let lines1 = [
        ax1_2d.rotated(&ax1_2d.location(), an_ang1),
        ax1_2d.rotated(&ax1_2d.location(), -an_ang1),
    ];
    let an_ang2 = c2.semi_angle();
    let lines2 = [
        ax2_2d.rotated(&ax2_2d.location(), an_ang2),
        ax2_2d.rotated(&ax2_2d.location(), -an_ang2),
    ];
    let mut min_dist = [INFINITE, INFINITE];
    for l1 in &lines1 {
        for l2 in &lines2 {
            let mut an_inter = IntAna2dAnaIntersection::new();
            an_inter.perform_lin_lin(l1, l2);
            if an_inter.is_done() {
                for k in 1..=an_inter.nb_points() {
                    let p = an_inter.point(k);
                    min_dist[0] = min_dist[0].min(p.param_on_first().abs());
                    min_dist[1] = min_dist[1].min(p.param_on_second().abs());
                }
            }
        }
    }
    min_dist[0].max(min_dist[1])
}

/// Cone ∩ cone. Port of `IntAna_QuadQuadGeo::Perform(gp_Cone, gp_Cone)` for the
/// tractable branches: coincident axes (two circles / a point / `Same`),
/// parallel axes with equal semi-angle (a conic in the plane through the two
/// apexes, from the plane∩cone closed form), coincident apexes (one or two
/// generatrix lines), and intersecting axes (the shared generatrix plus the
/// nested plane∩cone conic). Everything else, and the branch-4 guards that
/// raise `Standard_ConstructionError` in OCCT, return `None`.
///
/// The angular tolerance is not an input: `IntAna_QuadQuadGeo::Perform`
/// (`IntAna_QuadQuadGeo.cxx:1433-1476`) starts from `myEPSILON_ANGLE_CONE`
/// (`= Precision::Angular()`, set by `InitTolerances`, `cxx:351-359`) and only
/// relaxes it up to `Tol` when the axes are parallel, the semi-angles are
/// within `Tol`, and the axis distance exceeds `100 * Tol` — the relaxation is
/// `Tol / EstimDist(Con1, Con2)` (clamped to `[myEPSILON_ANGLE_CONE, Tol]`),
/// which keeps a far-away apex pair from dropping to the coincident-apex
/// branch.
pub fn quadric_quadric_cone_cone(c1: &GpCone, c2: &GpCone, tol: f64) -> QuadricIntersection {
    cone_cone_analytic(c1, c2, tol).inter
}

/// The full `IntAna_QuadQuadGeo::Perform(gp_Cone, gp_Cone)` outcome. Only
/// branch 4 sets `myCommonGen` / `myPChar` (`IntAna_QuadQuadGeo.hxx:250-291`);
/// `IntPatch_ImpImpIntersection::IntCoCo` reads them to add the shared
/// generatrix to the solution (`IntPatch_ImpImpIntersection.cxx:8893-8905`,
/// `:9149-9171`).
pub struct ConeConeAnalytic {
    /// `typeres` plus the conic, reduced to this port's line list.
    pub inter: QuadricIntersection,
    /// Branch 4 only: the generatrix through both apexes and the
    /// characteristic point (`cxx:1782-1794`).
    pub common_gen: Option<(GpLin, GpPnt)>,
}

/// `Perform(gp_Cone, gp_Cone)` (`IntAna_QuadQuadGeo.cxx:1433-1890`).
pub fn cone_cone_analytic(c1: &GpCone, c2: &GpCone, tol: f64) -> ConeConeAnalytic {
    let tg1 = c1.semi_angle().tan();
    let mut tg2 = c2.semi_angle().tan();
    if tg1 * tg2 < 0.0 {
        tg2 = -tg2;
    }
    let tol2 = tol * tol;
    let ap1 = c1.apex();
    let ap2 = c2.apex();
    let d_a1a2 = ap1.square_distance(&ap2);
    let (parallel, dist_axes, coplanar, pt_inter) =
        axe_operator(&c1.axis(), &c2.axis(), 1e-14, ANGULAR);

    // `cxx:1458-1476`: relax the angular tolerance for parallel axes.
    let mut tol_ang = ANGULAR;
    if (tg1 - tg2).abs() < tol && parallel {
        if dist_axes > 100.0 * tol {
            let min_sol_dist = estim_dist(c1, c2);
            if min_sol_dist < REAL_EPSILON {
                tol_ang = tol;
            } else {
                tol_ang = ANGULAR.max(tol / min_sol_dist);
                tol_ang = tol_ang.min(tol);
            }
        }
    }

    // 1 — coincident axes: two circles where the cone radii match (or the two
    // cones coincide / touch at the apex).
    if parallel && dist_axes < 1e-14 {
        let p = c1.apex();
        let d = c1.position().direction();
        let w = GpVec::from_xyz(d.xyz());
        let offset = w.dot(&GpVec::from_pnts(&p, &ap2));
        let inter = if (tg1 - tg2).abs() > ANGULAR {
            if offset.abs() < 1e-10 {
                return ConeConeAnalytic { inter: QuadricIntersection::Point(p), common_gen: None };
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
        };
        ConeConeAnalytic { inter, common_gen: None }
    }
    // 2 — parallel axes with (nearly) equal semi-angles: the intersection lies
    // in the plane through the two apexes; reduce to the plane∩cone conic.
    else if (tg1 - tg2).abs() < tol_ang && parallel {
        let da1 = c1.position().direction();
        let o1o2 = GpVec::from_pnts(&ap1, &ap2);
        let o1o2n = o1o2.normalized();
        let o1o2_da1 = GpVec::from_xyz(da1.xyz()).dot(&o1o2n);
        let o1_proj = o1o2n.subtracted(&GpVec::from_xyz(da1.xyz()).multiplied_scalar(o1o2_da1));
        let none = ConeConeAnalytic { inter: QuadricIntersection::None, common_gen: None };
        let Ok(db1) = GpDir::from_vec(&o1_proj) else { return none };
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
        let Ok(n) = GpDir::from_vec(&ortho) else { return none };
        let Some(pln) = plane_normal_at(p1, n) else { return none };
        // `cxx:1555`: the nested plane∩cone reduce passes `Tol` for *both*
        // the angular and the linear tolerance. `cxx:1557-1603` then re-maps
        // its result: the conics are kept as they are, an `IntAna_Line`
        // (one or two nested generatrices) is reported as its first line
        // twice (`pt1 = pt2 = H.Location()`, `dir1 = dir2 = H.Direction()`,
        // `nbint = 2`), and every other nested type (parabola, point,
        // empty) becomes `IntAna_NoGeometricSolution`.
        ConeConeAnalytic {
            inter: map_plane_cone_common_axis(quadric_quadric_plane_cone(&pln, c1, tol, tol)),
            common_gen: None,
        }
    }
    // 3 — coincident apexes: one or two generatrix lines (or `None`).
    else if d_a1a2 < tol2 {
        ConeConeAnalytic { inter: cone_cone_common_apex(c1, c2, tg1, tg2, tol), common_gen: None }
    }
    // 4 — the axes intersect, so the cones share the generatrix through their
    // apexes (`cxx:1759-1886`).
    else if coplanar && !parallel {
        // `A1A2.Intersect()` (`cxx:85-88`) is `coplanar && !parallel` and is the
        // only path on which `PtIntersect()` is set (`cxx:179-207`).
        let Some(pt_axe_int) = pt_inter else {
            return ConeConeAnalytic { inter: QuadricIntersection::None, common_gen: None };
        };
        match cone_cone_common_generatrix(c1, c2, tol, &pt_axe_int) {
            Some((inter, gen)) => ConeConeAnalytic { inter, common_gen: Some(gen) },
            None => ConeConeAnalytic { inter: QuadricIntersection::None, common_gen: None },
        }
    }
    // 5 — no geometric solution.
    else {
        ConeConeAnalytic { inter: QuadricIntersection::None, common_gen: None }
    }
}

/// `cxx:1557-1603` / `cxx:1836-1881`: the nested plane∩cone reduce of the
/// parallel-axes branches is re-mapped by hand. Branch 2 keeps the four conics
/// and reports a nested `IntAna_Line` as that line twice; branch 4 (see
/// [`cone_cone_common_generatrix`]) keeps the four conics and drops the line.
fn map_plane_cone_common_axis(nested: QuadricIntersection) -> QuadricIntersection {
    match nested {
        QuadricIntersection::Ellipse(e) => QuadricIntersection::Ellipse(e),
        QuadricIntersection::Circle(c) => QuadricIntersection::Circle(c),
        QuadricIntersection::Hyperbola(h) => QuadricIntersection::Hyperbola(h),
        QuadricIntersection::Line(l) | QuadricIntersection::TwoLines(l, _) => {
            QuadricIntersection::TwoLines(l, l)
        }
        // `cxx:1556` guards the re-mapping switch with
        // `if (INTER_QUAD_PLN.IsDone())`; a declined nested result leaves
        // `typeres` at its initial `IntAna_Empty` (`cxx:960`), i.e. this port's
        // `None`.
        QuadricIntersection::NotDone => QuadricIntersection::None,
        _ => QuadricIntersection::None,
    }
}

/// Branch 4 of `Perform(gp_Cone, gp_Cone)` (`cxx:1759-1886`): the axes
/// intersect, so the two cones share the generatrix through their apexes.
///
/// Returns the analytic section reduced to the nested plane∩cone conic and the
/// `(myCommonGen, myPChar)` pair, or `None` for `IntAna_NoGeometricSolution`.
/// `None` covers the two apex-membership guards (`cxx:1768`, `cxx:1776`, where
/// `myCommonGen` is still false) and the "Something wrong ???" return
/// (`cxx:1809`, where `myCommonGen` is already true but the consumer arm that
/// precedes the generatrix append is unported). `None` also covers the two
/// `gp_Dir`/`gp_Lin` degeneracies on which OCCT raises
/// `Standard_ConstructionError` (a zero-length direction): parallel
/// `aD1`/generatrix (`cxx:1801`) or a characteristic point coinciding with the
/// generatrices' intersection (`cxx:1825`).
fn cone_cone_common_generatrix(
    c1: &GpCone,
    c2: &GpCone,
    tol: f64,
    pt_axe_int: &GpPnt,
) -> Option<(QuadricIntersection, (GpLin, GpPnt))> {
    let tol2 = tol * tol;
    // `cxx:1762-1778`: each apex must belong to the other cone (`ElSLib` round
    // trip), otherwise `IntAna_NoGeometricSolution`.
    let (u, v) = slib::cone_parameters(&c2.position(), c2.radius(), c2.semi_angle(), &c1.apex());
    let p = slib::cone_value(c2, u, v);
    if c1.apex().square_distance(&p) > tol2 {
        return None;
    }
    let (u, v) = slib::cone_parameters(&c1.position(), c1.radius(), c1.semi_angle(), &c2.apex());
    let p = slib::cone_value(c1, u, v);
    if c2.apex().square_distance(&p) > tol2 {
        return None;
    }
    // `cxx:1782-1794`: the common generatrix and the characteristic point (the
    // axes intersection projected onto that generatrix).
    let ap1 = c1.apex();
    let ap2 = c2.apex();
    let gen_line = GpLin::from_pnt_dir(ap1, GpDir::from_vec(&GpVec::from_pnts(&ap1, &ap2)).ok()?);
    let u = clib::line_parameter(&gen_line, pt_axe_int);
    let pchar = clib::line_value(&gen_line, u);
    // `cxx:1795-1830`: the other generatrices of the two cones laying in the
    // maximal plane determine that plane's normal.
    let gen1 = gen_line.rotated(&c1.axis(), PI);
    let gen2 = gen_line.rotated(&c2.axis(), PI);
    let d1_val = gen1.direction();
    let d2_val = GpDir::from_vec(
        &GpVec::from_xyz(d1_val.xyz()).crossed(&GpVec::from_xyz(gen_line.direction().xyz())),
    )
    .ok()?;
    let n = if d1_val.is_parallel_tol(&gen2.direction(), ANGULAR) {
        d1_val.crossed(&d2_val).ok()?
    } else if lin_lin_square_distance(&gen1, &gen2) > tol2 {
        // `cxx:1805-1811`: "Something wrong ???" -> `typeres =
        // IntAna_NoGeometricSolution`. Note this `return` happens *after*
        // `cxx:1781` set `myCommonGen = true`, so OCCT's `HasCommonGen()`
        // stays true: `IntCoCo` would first run its nested
        // `IntAna_NoGeometricSolution` arm (`IntPatch_ImpImpIntersection.cxx
        // :9022-9140`, unported) and then still append the generatrix
        // (`cxx:9149-9171`). That arm is not ported and the port cannot honour
        // the ordering, so this path reports no solution and no generatrix.
        return None;
    } else {
        let d1 = gen1.position().direction();
        let d2 = gen2.position().direction();
        let o1 = gen1.location();
        let o2 = gen2.location();
        let d1_dot_d2 = d1.dot(&d2);
        let sin = 1.0 - d1_dot_d2 * d1_dot_d2;
        let o1o2 = GpVec::from_pnts(&o1, &o2);
        let mut u2 = GpVec::from_xyz(d1.xyz())
            .multiplied_scalar(o1o2.dot(&GpVec::from_xyz(d1.xyz())))
            .subtracted(&o1o2)
            .dot(&GpVec::from_xyz(d2.xyz()));
        u2 /= sin;
        let p_gint = clib::line_value(&gen2, u2);
        let d1 = GpDir::from_vec(&GpVec::from_pnts(&p_gint, &pchar)).ok()?;
        d1.crossed(&d2_val).ok()?
    };
    let pln = plane_normal_at(pchar, n)?;
    // `cxx:1833-1881`: the section lies in `anIntPln`; its nested plane∩cone
    // reduce keeps only the four conics — this switch has no `IntAna_Line` arm,
    // so a nested line (or point/empty/declined) is `IntAna_NoGeometricSolution`.
    // `cxx:1835` guards the switch with `if (INTER_QUAD_PLN.IsDone())`, so a
    // declined nested result leaves `typeres` at its initial `IntAna_Empty` —
    // also `None` here.
    let inter = match quadric_quadric_plane_cone(&pln, c1, tol, tol) {
        QuadricIntersection::Ellipse(e) => QuadricIntersection::Ellipse(e),
        QuadricIntersection::Circle(c) => QuadricIntersection::Circle(c),
        QuadricIntersection::Parabola(p) => QuadricIntersection::Parabola(p),
        QuadricIntersection::Hyperbola(h) => QuadricIntersection::Hyperbola(h),
        _ => QuadricIntersection::None,
    };
    Some((inter, (gen_line, pchar)))
}

/// `gp_Lin::SquareDistance(const gp_Lin&)` (`gp_Lin.cxx:28-39` +
/// `gp_Lin.hxx:120-124`): `Distance` squared, where `Distance` is the point-line
/// distance for exactly parallel lines and the perpendicular offset otherwise.
fn lin_lin_square_distance(a: &GpLin, b: &GpLin) -> f64 {
    let d = if a.position().direction().is_parallel_tol(&b.direction(), REAL_SMALL) {
        b.distance(&a.location())
    } else {
        match GpDir::from_vec(
            &GpVec::from_xyz(a.direction().xyz()).crossed(&GpVec::from_xyz(b.direction().xyz())),
        ) {
            Ok(dir) => GpVec::from_pnts(&a.location(), &b.location())
                .dot(&GpVec::from_xyz(dir.xyz()))
                .abs(),
            Err(_) => b.distance(&a.location()),
        }
    };
    d * d
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
    // `IntAna_QuadQuadGeo.cxx:1706-1713`: the nested `IntAna_QuadQuadGeo`
    // solve of the two planes is constructed with `(Tol, Tol)` (`cxx:1709`) and
    // a failure (`!IsDone() || 0 == NbSolutions()`, `cxx:1710`) is
    // `IntAna_NoGeometricSolution` (`cxx:1712`).
    let Some(lin) = plane_plane_line(&pln1, &pln2, tol, tol) else {
        return QuadricIntersection::None;
    };
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
            // `H` defaults to 0 (`IntAna_QuadQuadGeo.hxx:94`): only `IntPCy`
            // knows the face height (`IntPatch_ImpImpIntersection.cxx:2584-2588`),
            // and it calls the arm directly.
            quadric_quadric_plane_cylinder(p, c, tol_ang, tol, 0.0)
        }
        (Plane(p), Cone(c)) | (Cone(c), Plane(p)) => quadric_quadric_plane_cone(p, c, tol_ang, tol),
        (Cylinder(a), Cylinder(b)) => quadric_quadric_cylinder_cylinder(a, b, tol),
        (Cylinder(c), Sphere(s)) | (Sphere(s), Cylinder(c)) => quadric_quadric_cylinder_sphere(c, s, tol),
        (Sphere(s), Cone(c)) | (Cone(c), Sphere(s)) => quadric_quadric_sphere_cone(s, c, tol),
        (Cone(a), Cone(b)) => quadric_quadric_cone_cone(a, b, tol),
        (Cylinder(c), Cone(k)) | (Cone(k), Cylinder(c)) => quadric_quadric_cylinder_cone(c, k, tol),
    }
}

// ---------------------------------------------------------------------------
// IntAna_IntLinTorus
// ---------------------------------------------------------------------------

/// Result of [`line_torus_intersect`] — port of `IntAna_IntLinTorus`
/// (`IntAna_IntLinTorus.hxx:29-66`), which carries a `done` flag beside the
/// points and throws `StdFail_NotDone` from every accessor when it is false
/// (`IntAna_IntLinTorus.lxx:23-44`).
#[derive(Debug, Clone)]
pub struct LineTorusIntersection {
    done: bool,
    points: Vec<GpPnt>,
}

impl LineTorusIntersection {
    /// `IsDone` (`IntAna_IntLinTorus.lxx:18-21`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// The `NbPoints` / `Value(i)` payload (`IntAna_IntLinTorus.lxx:23-44`),
    /// in the order the quartic roots were accepted (`cxx:109-115`); empty when
    /// `!is_done()`.
    pub fn points(&self) -> &[GpPnt] {
        &self.points
    }
}

/// Intersection of a line with a torus. Port of `IntAna_IntLinTorus::Perform`
/// (`IntAna_IntLinTorus.cxx:43-133`): a quartic in the line parameter in the
/// torus reference frame, with each root verified by re-evaluation on the
/// torus. `done` is the OCCT flag: `false` when the quartic solver declined
/// (`cxx:126-131`) or when it returned roots that were all rejected on
/// back-substitution (`cxx:116-120`).
pub fn line_torus_intersect(l: &GpLin, t: &GpTorus) -> LineTorusIntersection {
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
    // `IntAna_IntLinTorus.cxx:88-133`. The `Value(1 .. NbSolutions)` list is
    // walked through the ported solver itself, not through the
    // `quartic_roots` adapter, because that adapter folds `!IsDone()` into an
    // empty root list while OCCT distinguishes the two: a declined solver
    // (`cxx:126-131`) still sets `done = false` where a completed solver with
    // zero roots (`cxx:121-125`) sets `done = true`.
    let mdpr = DirectPolynomialRoots::new5(a4, a3, a2, a1, a0);
    if !mdpr.is_done() {
        return LineTorusIntersection { done: false, points: Vec::new() };
    }
    let n = mdpr.nb_solutions();
    let mut nb_bad_sol = 0;
    for i in 1..=n {
        let tt = mdpr.value(i) + param_of_new_pl;
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
            nb_bad_sol += 1;
            continue;
        }
        out.push(p);
    }
    if n > 0 && out.is_empty() && nb_bad_sol == n {
        // `cxx:116-120`: every root failed the back-substitution check, so the
        // intersection is reported as not done rather than as an empty set.
        return LineTorusIntersection { done: false, points: Vec::new() };
    }
    // `cxx:104-108` stores each accepted root at `thePoint[nbsolvalid]`, i.e. in
    // root-index order; OCCT never reorders the point list here. An earlier
    // `out.sort_by(|a, b| a.x() ...)` had no counterpart in the `.cxx` (audit
    // `occt-geom.md` section 9 item (b)); the consumers that need an order sort
    // on their own (`intcurvesurface/solvers.rs:210` by curve parameter).
    LineTorusIntersection { done: true, points: out }
}
