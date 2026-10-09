//! Port of the `IntCurve_IntConicConic` Circle/Circle kernel set of
//! `IntCurve_IntConicConic_1.cxx`:
//! `CircleCircleGeometricIntersection` (`:154-355`),
//! `ProjectOnC2AndIntersectWithC2Domain` (`:58-151`) and
//! `Perform(const gp_Circ2d&, ..., const gp_Circ2d&, ...)` (`:807-1206`).

use std::f64::consts::PI;

use occt_core::elib::clib2d;
use occt_core::gp::{GpCirc2d, GpVec2d};
use occt_core::intimpargen::gen::determine_position;
use occt_core::intres2d::{
    IntRes2dDomain, IntRes2dIntersectionPoint, IntRes2dIntersectionSegment, IntRes2dTransition,
};
use occt_core::precision::RESOLUTION;

use super::int_conic_conic::IntCurveIntConicConic;
use super::int_conic_conic_tool::{
    determine_transition_lc, normalize_on_circle_domain, pi_p_pi, PeriodicInterval,
};

/// `ProjectOnC2AndIntersectWithC2Domain` (`_1.cxx:58-151`). Projects
/// `C1DomainAndRes` onto `Circle2`, intersects it with `DomainC2`, projects the
/// surviving boxes back onto `Circle1` and appends them to `SolutionC1` /
/// `SolutionC2`.
#[allow(clippy::too_many_arguments)]
fn project_on_c2_and_intersect_with_c2_domain(
    circle1: &GpCirc2d,
    circle2: &GpCirc2d,
    c1_domain_and_res: &PeriodicInterval,
    domain_c2: &PeriodicInterval,
    solution_c1: &mut [PeriodicInterval; 4],
    solution_c2: &mut [PeriodicInterval; 4],
    nb_sol_total: &mut usize,
    ident_circles: bool,
) {
    if c1_domain_and_res.is_null() {
        return; // `_1.cxx:68-71`
    }

    let c2inf = clib2d::circle_parameter_ax22d(
        circle2.position(),
        &clib2d::circle_value_ax22d(c1_domain_and_res.binf, circle1.position(), circle1.radius()),
    );
    let c2sup = clib2d::circle_parameter_ax22d(
        circle2.position(),
        &clib2d::circle_value_ax22d(c1_domain_and_res.bsup, circle1.position(), circle1.radius()),
    );

    // `PeriodicInterval C2Inter(C2inf, C2sup)` (`_1.cxx:84`) is built before
    // the `IdentCircles` adjustment below, so it uses the raw ends.
    let mut c2_inter = PeriodicInterval::from_bounds(c2inf, c2sup);

    if !ident_circles {
        if c2_inter.length() > PI {
            c2_inter.complement(); // `_1.cxx:88-92`
        }
    } else {
        let mut c2inf = c2inf;
        let mut c2sup = c2sup;
        if c2sup <= c2inf {
            c2sup += pi_p_pi();
        }
        if c2inf >= pi_p_pi() {
            c2sup -= pi_p_pi();
            c2inf -= pi_p_pi();
        }
        c2_inter.binf = c2inf;
        c2_inter.bsup = c2sup; // `_1.cxx:103-104`
        c2_inter.bsup = c2inf + c1_domain_and_res.bsup - c1_domain_and_res.binf; // `_1.cxx:105`
    }

    for i in 0..2 {
        // `i == 0` -> `FirstIntersection`, `i == 1` -> `SecondIntersection`
        // (`_1.cxx:113-115`); both shift `C2Inter` in place.
        let c2_inter_and_domain = if i == 0 {
            domain_c2.first_intersection(&mut c2_inter)
        } else {
            domain_c2.second_intersection(&mut c2_inter)
        };

        if c2_inter_and_domain.is_null() {
            continue;
        }

        let c1inf = clib2d::circle_parameter_ax22d(
            circle1.position(),
            &clib2d::circle_value_ax22d(
                c2_inter_and_domain.binf,
                circle2.position(),
                circle2.radius(),
            ),
        );
        let c1sup = clib2d::circle_parameter_ax22d(
            circle1.position(),
            &clib2d::circle_value_ax22d(
                c2_inter_and_domain.bsup,
                circle2.position(),
                circle2.radius(),
            ),
        );

        solution_c1[*nb_sol_total] = PeriodicInterval::from_bounds(c1inf, c1sup);
        if !ident_circles {
            if solution_c1[*nb_sol_total].length() > PI {
                solution_c1[*nb_sol_total].complement(); // `_1.cxx:128-133`
            }
        } else {
            if solution_c1[*nb_sol_total].bsup <= solution_c1[*nb_sol_total].binf {
                solution_c1[*nb_sol_total].bsup += pi_p_pi(); // `_1.cxx:137-140`
            }
            if solution_c1[*nb_sol_total].binf >= pi_p_pi() {
                solution_c1[*nb_sol_total].binf -= pi_p_pi(); // `_1.cxx:141-145`
                solution_c1[*nb_sol_total].bsup -= pi_p_pi();
            }
        }
        solution_c2[*nb_sol_total] = c2_inter_and_domain;
        *nb_sol_total += 1;
    }
}

/// `CircleCircleGeometricIntersection` (`_1.cxx:154-355`). Returns `nbsol`
/// (`0`, `1`, `2` or `3`) and writes the two solution intervals on `C1`.
pub(crate) fn circle_circle_geometric_intersection(
    c1: &GpCirc2d,
    c2: &GpCirc2d,
    tol: f64,
    tol_tang: f64,
    c1_res1: &mut PeriodicInterval,
    c1_res2: &mut PeriodicInterval,
) -> i32 {
    let mut c1_binf1 = 0.0f64;
    let mut c1_binf2 = 0.0f64;
    let mut c1_bsup1 = 0.0f64;
    let mut c1_bsup2 = 0.0f64;
    let d_o1_o2 = c1.location().distance(&c2.location()); // `_1.cxx:166`
    let r1 = c1.radius();
    let r2 = c2.radius();
    let abs_r1m_r2 = (r1 - r2).abs(); // `_1.cxx:168`
    let nbsol: i32;

    if d_o1_o2 > (r1 + r2 + tol) {
        // `_1.cxx:169-182`
        if d_o1_o2 > (r1 + r2 + tol_tang) {
            return 0; // `_1.cxx:171-174`
        } else {
            c1_binf1 = 0.0;
            c1_bsup1 = 0.0;
            nbsol = 1;
        }
    } else if d_o1_o2 <= tol && abs_r1m_r2 <= tol {
        return 3; // `_1.cxx:184-188`
    } else {
        let r1p_r2 = r1 + r2; // `_1.cxx:192`
        let r1p_tol = r1 + tol;
        let r1m_tol = r1 - tol;
        let r2r2 = r2 * r2;
        let r1p_tol_r1p_tol = r1p_tol * r1p_tol;
        let r1m_tol_r1m_tol = r1m_tol * r1m_tol;
        let d_o1_o2d_o1_o2 = d_o1_o2 * d_o1_o2;
        let mut d_alpha1;

        if d_o1_o2 > r1p_r2 - tol {
            // `_1.cxx:206-215`: C2 cuts the circle C1+ only.
            let dx = (r1p_tol_r1p_tol + d_o1_o2d_o1_o2 - r2r2) / (d_o1_o2 + d_o1_o2);
            let mut dy = r1p_tol_r1p_tol - dx * dx;
            dy = if dy >= 0.0 { dy.sqrt() } else { 0.0 };
            d_alpha1 = dy.atan2(dx);

            c1_binf1 = -d_alpha1;
            c1_bsup1 = d_alpha1;
            nbsol = 1;
        } else if d_o1_o2 > abs_r1m_r2 - tol {
            // `_1.cxx:220-284`: C2 cuts C1- and C1+.
            let dx = (r1p_tol_r1p_tol + d_o1_o2d_o1_o2 - r2r2) / (d_o1_o2 + d_o1_o2);
            let mut dy = r1p_tol_r1p_tol - dx * dx;
            dy = if dy >= 0.0 { dy.sqrt() } else { 0.0 };
            d_alpha1 = dy.atan2(dx);
            c1_binf1 = -d_alpha1;
            c1_bsup2 = d_alpha1;

            // `_1.cxx:231-238`: intersection C2 / C1-.
            let dx = (r1m_tol_r1m_tol + d_o1_o2d_o1_o2 - r2r2) / (d_o1_o2 + d_o1_o2);
            let mut dy = r1m_tol_r1m_tol - dx * dx;
            dy = if dy >= 0.0 { dy.sqrt() } else { 0.0 };
            d_alpha1 = dy.atan2(dx);

            c1_binf2 = d_alpha1;
            c1_bsup1 = -d_alpha1;
            let mut nb = 2;

            if dy == 0.0 {
                // `_1.cxx:245-248`: the two inner bounds coincide.
                c1_bsup1 = c1_bsup2;
                nb = 1;
            } else {
                if c1_binf1 > c1_bsup1 {
                    std::mem::swap(&mut c1_binf1, &mut c1_bsup1); // `_1.cxx:252-257`
                }
                if c1_binf2 > c1_bsup2 {
                    std::mem::swap(&mut c1_binf2, &mut c1_bsup2); // `_1.cxx:258-263`
                }
                if ((c1_binf1 <= c1_bsup2) && (c1_binf1 >= c1_binf2))
                    || ((c1_bsup1 <= c1_bsup2) && (c1_bsup1 >= c1_binf2))
                {
                    // `_1.cxx:264-283`: the two intervals are really one.
                    if c1_binf1 > c1_binf2 {
                        c1_binf1 = c1_binf2;
                    }
                    if c1_binf1 > c1_bsup2 {
                        c1_binf1 = c1_bsup2;
                    }
                    if c1_bsup1 < c1_binf2 {
                        c1_bsup1 = c1_binf2;
                    }
                    if c1_bsup1 < c1_bsup2 {
                        c1_bsup1 = c1_bsup2;
                    }
                    nb = 1;
                }
            }
            nbsol = nb;
        } else {
            // `_1.cxx:287-300`.
            if (d_o1_o2 > abs_r1m_r2 - tol_tang) && (abs_r1m_r2 - tol_tang) > 0.0 {
                c1_binf1 = 0.0;
                c1_bsup1 = 0.0;
                nbsol = 1;
            } else {
                return 0;
            }
        }
    }

    // `_1.cxx:308-329`: back into the frame of C1.
    let axe1 = GpVec2d::from_dir2d(&c1.x_axis().vdir); // `_1.cxx:312`
    let axe_o1_o2 = GpVec2d::new(
        c2.location().x() - c1.location().x(),
        c2.location().y() - c1.location().y(),
    ); // `_1.cxx:313` `gp_Vec2d(C1.Location(), C2.Location())`

    let mut d_angle1;
    if axe_o1_o2.magnitude() <= RESOLUTION {
        d_angle1 = axe1.angle(&GpVec2d::from_dir2d(&c2.x_axis().vdir)); // `_1.cxx:317-319`
    } else {
        d_angle1 = axe1.angle(&axe_o1_o2); // `_1.cxx:320-323`
    }

    if !c1.is_direct() {
        d_angle1 = -d_angle1; // `_1.cxx:325-328`
    }

    c1_binf1 += d_angle1;
    c1_bsup1 += d_angle1;

    c1_res1.set_values(c1_binf1, c1_bsup1); // `_1.cxx:336`
    if c1_res1.length() > PI {
        c1_res1.complement(); // `_1.cxx:337-340`
    }

    if nbsol == 2 {
        c1_binf2 += d_angle1;
        c1_bsup2 += d_angle1;
        c1_res2.set_values(c1_binf2, c1_bsup2); // `_1.cxx:345-347`
        if c1_res2.length() > PI {
            c1_res2.complement(); // `_1.cxx:348-350`
        }
    } else {
        c1_res2.set_null(); // `_1.cxx:352-354`
    }

    nbsol
}

impl IntCurveIntConicConic {
    /// `Perform(const gp_Circ2d&, D1, const gp_Circ2d&, D2, TolConf, Tol)`
    /// (`IntCurve_IntConicConic_1.cxx:807-1206`).
    #[allow(clippy::too_many_arguments)]
    pub fn perform_circle_circle(
        &mut self,
        circle1: &GpCirc2d,
        domain_circ1: &IntRes2dDomain,
        circle2_in: &GpCirc2d,
        domain_circ2_in: &IntRes2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        let mut circle2 = *circle2_in; // `_1.cxx:815`
        let mut domain_circ2 = *domain_circ2_in; // `_1.cxx:816`
        let indirect_circles;

        if circle1.is_direct() != circle2_in.is_direct() {
            indirect_circles = true; // `_1.cxx:819-831`
            circle2 = circle2_in.reversed();
            domain_circ2.set_bounded(
                domain_circ2_in.last_point(),
                pi_p_pi() - domain_circ2_in.last_parameter(),
                domain_circ2_in.last_tolerance(),
                domain_circ2_in.first_point(),
                pi_p_pi() - domain_circ2_in.first_parameter(),
                domain_circ2_in.first_tolerance(),
            );
            domain_circ2.set_equivalent_parameters(0.0, pi_p_pi());
        } else {
            indirect_circles = false;
        }

        self.result.reset_fields(); // `_1.cxx:833` `this->ResetFields()`

        let mut c1_int1 = PeriodicInterval::default();
        let mut c1_int2 = PeriodicInterval::default();

        // `_1.cxx:836-837`: geometric intersection, domain-agnostic.
        let nbsol =
            circle_circle_geometric_intersection(circle1, &circle2, tol_conf, tol, &mut c1_int1, &mut c1_int2);
        self.result.done = true; // `_1.cxx:838`

        if nbsol == 0 {
            return; // `_1.cxx:840-843`
        }

        // `_1.cxx:845-862`: bring the C1 domain back between 0 and 2*PI.
        let mut c1_domain = PeriodicInterval::from_domain(domain_circ1);
        let mut deltat = c1_domain.bsup - c1_domain.binf;
        if deltat >= pi_p_pi() {
            deltat = next_after_pi_p_pi(); // `_1.cxx:848-852`
        }
        while c1_domain.binf >= pi_p_pi() {
            c1_domain.binf -= pi_p_pi();
        }
        while c1_domain.binf < 0.0 {
            c1_domain.binf += pi_p_pi();
        }
        c1_domain.bsup = c1_domain.binf + deltat;

        // `_1.cxx:864-880`: same on C2.
        let mut c2_domain = PeriodicInterval::from_domain(&domain_circ2);
        deltat = c2_domain.bsup - c2_domain.binf;
        if deltat >= pi_p_pi() {
            deltat = next_after_pi_p_pi(); // `_1.cxx:868-870`
        }
        while c2_domain.binf >= pi_p_pi() {
            c2_domain.binf -= pi_p_pi();
        }
        while c2_domain.binf < 0.0 {
            c2_domain.binf += pi_p_pi();
        }
        c2_domain.bsup = c2_domain.binf + deltat;

        // `_1.cxx:882-901`.
        let ident_circles;
        if nbsol > 2 {
            // The two circles coincide within `Tol`.
            c1_int1.set_values(0.0, pi_p_pi());
            c1_int2.set_null();
            ident_circles = true;
        } else {
            ident_circles = false;
        }

        let mut nb_sol_total = 0usize;
        let mut solution_c1 = [PeriodicInterval::default(); 4];
        let mut solution_c2 = [PeriodicInterval::default(); 4];

        // `_1.cxx:915-940`: first and second intersection of C1_Int1 with the
        // C1 domain.
        let mut c1_domain_and_res = c1_domain.first_intersection(&mut c1_int1);
        project_on_c2_and_intersect_with_c2_domain(
            circle1,
            &circle2,
            &c1_domain_and_res,
            &c2_domain,
            &mut solution_c1,
            &mut solution_c2,
            &mut nb_sol_total,
            ident_circles,
        );

        c1_domain_and_res = c1_domain.second_intersection(&mut c1_int1);
        project_on_c2_and_intersect_with_c2_domain(
            circle1,
            &circle2,
            &c1_domain_and_res,
            &c2_domain,
            &mut solution_c1,
            &mut solution_c2,
            &mut nb_sol_total,
            ident_circles,
        );

        if nbsol == 2 {
            // `_1.cxx:944-964`: same for the second geometric interval.
            c1_domain_and_res = c1_domain.first_intersection(&mut c1_int2);
            project_on_c2_and_intersect_with_c2_domain(
                circle1,
                &circle2,
                &c1_domain_and_res,
                &c2_domain,
                &mut solution_c1,
                &mut solution_c2,
                &mut nb_sol_total,
                ident_circles,
            );

            c1_domain_and_res = c1_domain.second_intersection(&mut c1_int2);
            project_on_c2_and_intersect_with_c2_domain(
                circle1,
                &circle2,
                &c1_domain_and_res,
                &c2_domain,
                &mut solution_c1,
                &mut solution_c2,
                &mut nb_sol_total,
                ident_circles,
            );
        }

        // `_1.cxx:975-996`: intervals below `Tol` collapse to a point.
        let r1 = circle1.radius();
        let r2 = circle2.radius();
        let mut tol2 = tol + tol;
        if tol < 1.0e-10 {
            tol2 = 1.0e-10;
        }

        for i in 0..nb_sol_total {
            if ((r1 * solution_c1[i].length()) <= tol2) && ((r2 * solution_c2[i].length()) <= tol2) {
                let t = (solution_c1[i].binf + solution_c1[i].bsup) * 0.5;
                solution_c1[i].binf = t;
                solution_c1[i].bsup = t;

                let t = (solution_c2[i].binf + solution_c2[i].bsup) * 0.5;
                solution_c2[i].binf = t;
                solution_c2[i].bsup = t;
            }
        }

        let axis2_c1 = *circle1.position(); // `_1.cxx:1000` `Circle1.Axis()`
        let axis2_c2 = *circle2.position(); // `_1.cxx:1001` `Circle2.Axis()`

        let is_opposite = circle1.location().square_distance(&circle2.location())
            > (r1 * r1 + r2 * r2); // `_1.cxx:1007`

        for i in 0..nb_sol_total {
            // `_1.cxx:1014-1021`.
            let mut c2inf = if is_opposite { solution_c2[i].bsup } else { solution_c2[i].binf };
            let mut c2sup = if is_opposite { solution_c2[i].binf } else { solution_c2[i].bsup };
            let c1tinf = solution_c1[i].binf;
            let c2tinf = c2inf;
            let mut c1inf = normalize_on_circle_domain(c1tinf, domain_circ1);
            c2inf = normalize_on_circle_domain(c2tinf, &domain_circ2);

            // `_1.cxx:1022-1074`: clip onto the domains.
            let mut is_out_of_range = false;
            if c1inf < domain_circ1.first_parameter() {
                if c1tinf < domain_circ1.first_parameter() {
                    c1inf = domain_circ1.first_parameter();
                    is_out_of_range = true;
                } else {
                    c1inf = c1tinf;
                }
            }
            if c1inf > domain_circ1.last_parameter() {
                if c1tinf > domain_circ1.last_parameter() {
                    c1inf = domain_circ1.last_parameter();
                    is_out_of_range = true;
                } else {
                    c1inf = c1tinf;
                }
            }
            if c2inf < domain_circ2.first_parameter() {
                if c2tinf < domain_circ2.first_parameter() {
                    c2inf = domain_circ2.first_parameter();
                    is_out_of_range = true;
                } else {
                    c2inf = c2tinf;
                }
            }
            if c2inf > domain_circ2.last_parameter() {
                if c2tinf > domain_circ2.last_parameter() {
                    c2inf = domain_circ2.last_parameter();
                    is_out_of_range = true;
                } else {
                    c2inf = c2tinf;
                }
            }

            // `_1.cxx:1076-1087`: no solution inside the parametric range.
            if is_out_of_range {
                let (a_p1, _v11, _v12) = clib2d::circle_d2_ax22d(c1inf, &axis2_c1, r1);
                let (a_p2, _v21, _v22) = clib2d::circle_d2_ax22d(c2inf, &axis2_c2, r2);
                if a_p1.square_distance(&a_p2) > tol2 * tol2 {
                    continue;
                }
            }

            if indirect_circles {
                // `_1.cxx:1089-1150`.
                let (p1a, mut tan1, norm1) = clib2d::circle_d2_ax22d(c1inf, &axis2_c1, r1);
                let (p2a, mut tan2, norm2) = clib2d::circle_d2_ax22d(c2inf, &axis2_c2, r2);
                tan2.reverse();

                let pos1a = determine_position(domain_circ1, &p1a, c1inf);
                let pos2a = determine_position(domain_circ2_in, &p2a, pi_p_pi() - c2inf);
                let mut t1a = IntRes2dTransition::new();
                let mut t2a = IntRes2dTransition::new();
                determine_transition_lc(
                    pos1a,
                    &mut tan1,
                    &norm1,
                    &mut t1a,
                    pos2a,
                    &mut tan2,
                    &norm2,
                    &mut t2a,
                    tol,
                );

                let new_point1 = IntRes2dIntersectionPoint::with_transitions(
                    &p1a,
                    c1inf,
                    pi_p_pi() - c2inf,
                    &t1a,
                    &t2a,
                    false,
                );

                if (solution_c1[i].length() > 0.0) || (solution_c2[i].length() > 0.0) {
                    // `_1.cxx:1103-1148`: a non-degenerate interval.
                    let mut c1sup = normalize_on_circle_domain(solution_c1[i].bsup, domain_circ1);
                    if c1sup < c1inf {
                        c1sup += pi_p_pi();
                    }
                    c2sup = normalize_on_circle_domain(c2sup, &domain_circ2);

                    let (p1b, mut tan1, norm1) = clib2d::circle_d2_ax22d(c1sup, &axis2_c1, r1);
                    let (p2b, mut tan2, norm2) = clib2d::circle_d2_ax22d(c2sup, &axis2_c2, r2);
                    tan2.reverse();

                    let pos1b = determine_position(domain_circ1, &p1b, c1sup);
                    let pos2b = determine_position(domain_circ2_in, &p2b, pi_p_pi() - c2sup);
                    let mut t1b = IntRes2dTransition::new();
                    let mut t2b = IntRes2dTransition::new();
                    determine_transition_lc(
                        pos1b,
                        &mut tan1,
                        &norm1,
                        &mut t1b,
                        pos2b,
                        &mut tan2,
                        &norm2,
                        &mut t2b,
                        tol,
                    );

                    if is_opposite {
                        if nbsol != 3 {
                            if c2inf < c2sup {
                                c2inf += pi_p_pi();
                            }
                        }
                    } else {
                        if nbsol != 3 {
                            if c2sup < c2inf {
                                c2sup += pi_p_pi();
                            }
                        }
                    }

                    let new_point2 = IntRes2dIntersectionPoint::with_transitions(
                        &p1b,
                        c1sup,
                        pi_p_pi() - c2sup,
                        &t1b,
                        &t2b,
                        false,
                    );
                    let new_seg = IntRes2dIntersectionSegment::from_two_points(
                        &new_point1,
                        &new_point2,
                        !is_opposite,
                        false,
                    );
                    self.result.append_segment(&new_seg);
                } else {
                    self.result.append_point(&new_point1);
                }
            } else {
                // `_1.cxx:1152-1205`.
                let (p1a, mut tan1, norm1) = clib2d::circle_d2_ax22d(c1inf, &axis2_c1, r1);
                let (p2a, mut tan2, norm2) = clib2d::circle_d2_ax22d(c2inf, &axis2_c2, r2);

                let pos1a = determine_position(domain_circ1, &p1a, c1inf);
                let pos2a = determine_position(&domain_circ2, &p2a, c2inf);
                let mut t1a = IntRes2dTransition::new();
                let mut t2a = IntRes2dTransition::new();
                determine_transition_lc(
                    pos1a,
                    &mut tan1,
                    &norm1,
                    &mut t1a,
                    pos2a,
                    &mut tan2,
                    &norm2,
                    &mut t2a,
                    tol,
                );

                let new_point1 = IntRes2dIntersectionPoint::with_transitions(
                    &p1a, c1inf, c2inf, &t1a, &t2a, false,
                );

                if (solution_c1[i].length() > 0.0) || (solution_c2[i].length() > 0.0) {
                    // `_1.cxx:1162-1198`.
                    let mut c1sup = normalize_on_circle_domain(solution_c1[i].bsup, domain_circ1);
                    if c1sup < c1inf {
                        c1sup += pi_p_pi();
                    }
                    c2sup = normalize_on_circle_domain(c2sup, &domain_circ2);

                    let (p1b, mut tan1, norm1) = clib2d::circle_d2_ax22d(c1sup, &axis2_c1, r1);
                    let (p2b, mut tan2, norm2) = clib2d::circle_d2_ax22d(c2sup, &axis2_c2, r2);

                    let pos1b = determine_position(domain_circ1, &p1b, c1sup);
                    let pos2b = determine_position(&domain_circ2, &p2b, c2sup);
                    let mut t1b = IntRes2dTransition::new();
                    let mut t2b = IntRes2dTransition::new();
                    determine_transition_lc(
                        pos1b,
                        &mut tan1,
                        &norm1,
                        &mut t1b,
                        pos2b,
                        &mut tan2,
                        &norm2,
                        &mut t2b,
                        tol,
                    );

                    if is_opposite {
                        if c2inf < c2sup {
                            c2inf += pi_p_pi();
                        }
                    } else {
                        if c2sup < c2inf {
                            c2sup += pi_p_pi();
                        }
                    }

                    let new_point2 = IntRes2dIntersectionPoint::with_transitions(
                        &p1b, c1sup, c2sup, &t1b, &t2b, false,
                    );
                    let new_seg = IntRes2dIntersectionSegment::from_two_points(
                        &new_point1,
                        &new_point2,
                        is_opposite,
                        false,
                    );
                    self.result.append_segment(&new_seg);
                } else {
                    self.result.append_point(&new_point1);
                }
            }
        }
    }
}

/// `std::nextafter(PIpPI, 0.)` (`_1.cxx:850`, `:869`).
fn next_after_pi_p_pi() -> f64 {
    f64::from_bits(pi_p_pi().to_bits() - 1)
}
