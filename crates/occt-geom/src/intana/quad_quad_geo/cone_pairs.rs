//! `EstimDist` and the `IntAna_QuadQuadGeo::Perform(gp_Cone, gp_Cone)` branches.

use super::super::prelude::*;
use super::*;

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
