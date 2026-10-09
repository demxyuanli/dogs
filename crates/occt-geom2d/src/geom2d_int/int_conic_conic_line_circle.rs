//! Port of the `IntCurve_IntConicConic` Line/Circle kernel set of
//! `IntCurve_IntConicConic_1.cxx`:
//! `LineCircleGeometricIntersection` (`:451-636`),
//! `ProjectOnLAndIntersectWithLDomain(const gp_Circ2d&, ...)` (`:360-442`) and
//! `Perform(const gp_Lin2d&, ..., const gp_Circ2d&, ...)` (`:2236-2618`).

use std::f64::consts::PI;

use occt_core::elib::clib2d;
use occt_core::gp::{GpCirc2d, GpLin2d, GpPnt2d, GpVec2d};
use occt_core::intimpargen::gen::determine_position;
use occt_core::intres2d::{
    IntRes2dDomain, IntRes2dIntersectionPoint, IntRes2dIntersectionSegment, IntRes2dPosition,
    IntRes2dTransition,
};

use super::int_conic_conic::IntCurveIntConicConic;
use super::int_conic_conic_tool::{
    determine_transition_lc, normalize_on_circle_domain, pi_p_pi, Interval, PeriodicInterval,
};

/// `const double PIsur2 = 0.5 * M_PI` (`_1.cxx:46`).
fn pi_sur_2() -> f64 {
    0.5 * PI
}

/// `LineCircleGeometricIntersection` (`_1.cxx:451-636`). Returns `nbsol` and
/// writes the two solution intervals on the circle.
pub(crate) fn line_circle_geometric_intersection(
    line: &GpLin2d,
    circle: &GpCirc2d,
    tol: f64,
    tol_tang: f64,
    c_int1: &mut PeriodicInterval,
    c_int2: &mut PeriodicInterval,
) -> i32 {
    let d_o1_o2 = line.distance(&circle.location());
    let r = circle.radius();
    let r_m_tol = r - tol;
    let mut binf1 = 0.0;
    let mut bsup1 = 0.0;
    let mut binf2 = 0.0;
    let mut bsup2 = 0.0;
    let mut nbsol;

    if d_o1_o2 > (r + tol) {
        // no intersection with the tolerance tube
        if d_o1_o2 > (r + tol_tang) {
            return 0;
        }
        binf1 = 0.0;
        bsup1 = 0.0;
        nbsol = 1;
    } else {
        let mut b2_sol = false;
        let mut d_alpha1;
        // Line cuts Circle+ (= C(x1, y1, R1 + Tol))
        if r > d_o1_o2 + tol_tang {
            let a_tol2 = tol * tol;
            let a_x2 = 4.0 * (r * r - d_o1_o2 * d_o1_o2);
            if a_x2 > a_tol2 {
                b2_sol = !b2_sol;
            }
        }
        if d_o1_o2 > r_m_tol && !b2_sol {
            // `_1.cxx:502-504`: `dy` is the constant 0.0 (the squared term is
            // commented out as a patch in OCCT), so this is `atan2(0, dx)`.
            let dx = d_o1_o2;
            let dy: f64 = 0.0;
            let dy = if dy >= 0.0 { dy.sqrt() } else { 0.0 };
            d_alpha1 = dy.atan2(dx);

            binf1 = -d_alpha1;
            bsup1 = d_alpha1;
            nbsol = 1;
        } else {
            // Intersection Line / Circle+
            let dx = d_o1_o2;
            let mut dy = r * r - dx * dx;
            dy = if dy >= 0.0 { dy.sqrt() } else { 0.0 };
            d_alpha1 = dy.atan2(dx);
            binf1 = -d_alpha1;
            bsup2 = d_alpha1;

            // Intersection Line / Circle-
            dy = r * r - dx * dx;
            dy = if dy >= 0.0 { dy.sqrt() } else { 0.0 };
            d_alpha1 = dy.atan2(dx);

            binf2 = d_alpha1;
            bsup1 = -d_alpha1;

            if (d_alpha1 * r) < tol.max(tol_tang) {
                bsup1 = bsup2;
                nbsol = 1;
            } else {
                nbsol = 2;
            }
        }
    }

    // Back to the frame of C1.
    let mut d_angle1 = circle.x_axis().vdir.angle(line.direction());
    let (a, b, c) = line.coefficients();
    let d = a * circle.location().x() + b * circle.location().y() + c;

    if d > 0.0 {
        d_angle1 += pi_sur_2();
    } else {
        d_angle1 -= pi_sur_2();
    }

    if d_angle1 < 0.0 {
        d_angle1 += pi_p_pi();
    } else if d_angle1 > pi_p_pi() {
        d_angle1 -= pi_p_pi();
    }

    binf1 += d_angle1;
    bsup1 += d_angle1;

    if !circle.is_direct() {
        let t = binf1;
        binf1 = bsup1;
        bsup1 = t;
        binf1 = -binf1;
        bsup1 = -bsup1;
    }

    c_int1.set_values(binf1, bsup1);
    if c_int1.length() > PI {
        c_int1.complement();
    }

    if nbsol == 2 {
        binf2 += d_angle1;
        bsup2 += d_angle1;

        if !circle.is_direct() {
            let t = binf2;
            binf2 = bsup2;
            bsup2 = t;
            binf2 = -binf2;
            bsup2 = -bsup2;
        }

        c_int2.set_values(binf2, bsup2);
        if c_int2.length() > PI {
            c_int2.complement();
        }
    } else {
        // Modified by Sergey KHROMOV - Thu Oct 26 17:51:05 2000
        if c_int1.bsup > pi_p_pi() && c_int1.binf < pi_p_pi() {
            nbsol = 2;
            binf2 = c_int1.binf;
            bsup2 = pi_p_pi();
            binf1 = 0.0;
            c_int1.set_values(binf1, c_int1.bsup - pi_p_pi());
            if c_int1.length() > PI {
                c_int1.complement();
            }
            c_int2.set_values(binf2, bsup2);
            if c_int2.length() > PI {
                c_int2.complement();
            }
        }
    }

    nbsol
}

/// `ProjectOnLAndIntersectWithLDomain(const gp_Circ2d&, ...)`
/// (`_1.cxx:360-442`).
#[allow(clippy::too_many_arguments)]
fn project_on_l_and_intersect_with_l_domain(
    circle: &GpCirc2d,
    line: &GpLin2d,
    c_domain_and_res: &PeriodicInterval,
    l_domain: &Interval,
    circle_solution: &mut [PeriodicInterval; 4],
    line_solution: &mut [Interval; 4],
    nb_sol_total: &mut usize,
    ref_line_domain: &IntRes2dDomain,
) {
    if c_domain_and_res.is_null() {
        return;
    }

    let linf = clib2d::parameter_lin2d(
        line,
        &clib2d::circle_value_ax22d(c_domain_and_res.binf, circle.position(), circle.radius()),
    );
    let lsup = clib2d::parameter_lin2d(
        line,
        &clib2d::circle_value_ax22d(c_domain_and_res.bsup, circle.position(), circle.radius()),
    );

    // Necessarily bounded.
    let l_inter = Interval::from_bounds(linf, lsup);
    let mut l_inter_and_domain = l_domain.intersection_with_bounded(&l_inter);

    if !l_inter_and_domain.is_null {
        let dom_linf = if ref_line_domain.has_first_point() {
            ref_line_domain.first_parameter()
        } else {
            -occt_core::precision::INFINITE
        };
        let dom_lsup = if ref_line_domain.has_last_point() {
            ref_line_domain.last_parameter()
        } else {
            occt_core::precision::INFINITE
        };

        let mut linf = l_inter_and_domain.binf;
        let mut lsup = l_inter_and_domain.bsup;

        if linf < dom_linf {
            linf = dom_linf;
        }
        if lsup < dom_linf {
            lsup = dom_linf;
        }
        if linf > dom_lsup {
            linf = dom_lsup;
        }
        if lsup > dom_lsup {
            lsup = dom_lsup;
        }

        l_inter_and_domain.binf = linf;
        l_inter_and_domain.bsup = lsup;

        let mut cinf = c_domain_and_res.binf;
        let mut csup = c_domain_and_res.bsup;
        if cinf >= csup {
            cinf = c_domain_and_res.binf;
            csup = c_domain_and_res.bsup;
        }
        circle_solution[*nb_sol_total] = PeriodicInterval::from_bounds(cinf, csup);
        if circle_solution[*nb_sol_total].length() > PI {
            circle_solution[*nb_sol_total].complement();
        }

        line_solution[*nb_sol_total] = l_inter_and_domain;
        *nb_sol_total += 1;
    }
}

impl IntCurveIntConicConic {
    /// `Perform(const gp_Lin2d& Line, const IntRes2d_Domain& LIG_Domain,
    /// const gp_Circ2d& Circle, const IntRes2d_Domain& CIRC_Domain, TolConf,
    /// Tol)` (`_1.cxx:2236-2618`).
    pub fn perform_line_circle(
        &mut self,
        line: &GpLin2d,
        lig_domain: &IntRes2dDomain,
        circle: &GpCirc2d,
        circ_domain: &IntRes2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        let the_reversed_parameters = self.result.reversed_parameters();
        self.result.reset_fields();
        self.result.set_reversed_parameters(the_reversed_parameters);

        let mut c_int1 = PeriodicInterval::default();
        let mut c_int2 = PeriodicInterval::default();

        let mut nbsol = line_circle_geometric_intersection(
            line, circle, tol_conf, tol, &mut c_int1, &mut c_int2,
        );

        self.result.done = true;

        if nbsol == 0 {
            return;
        }

        // Modified by Sergey KHROMOV - Mon Dec 18 11:13:18 2000
        if nbsol == 2 && c_int2.bsup == c_int1.binf + pi_p_pi() {
            let first_bound = circ_domain.first_parameter();
            let last_bound = circ_domain.last_parameter();
            let first_tol = circ_domain.first_tolerance();
            let last_tol = circ_domain.last_tolerance();
            if c_int1.binf == 0.0 && first_bound - first_tol > c_int1.bsup {
                nbsol = 1;
                c_int1.set_values(c_int2.binf, c_int2.bsup);
            } else if c_int2.bsup == pi_p_pi() && last_bound + last_tol < c_int2.binf {
                nbsol = 1;
            }
        }
        // Modified by Sergey KHROMOV - Mon Dec 18 11:13:20 2000 End

        let mut c_domain = PeriodicInterval::from_domain(circ_domain);
        let mut deltat = c_domain.bsup - c_domain.binf;
        while c_domain.binf >= pi_p_pi() {
            c_domain.binf -= pi_p_pi();
        }
        while c_domain.binf < 0.0 {
            c_domain.binf += pi_p_pi();
        }
        c_domain.bsup = c_domain.binf + deltat;

        // Ajout: Jeudi 28 mars 96 -- the domains are artificially enlarged.
        let mut binf_modif = c_domain.binf;
        let mut bsup_modif = c_domain.bsup;
        binf_modif -= circ_domain.first_tolerance() / circle.radius();
        bsup_modif += circ_domain.last_tolerance() / circle.radius();
        deltat = bsup_modif - binf_modif;
        if deltat <= pi_p_pi() {
            c_domain.binf = binf_modif;
            c_domain.bsup = bsup_modif;
        } else {
            let mut t = pi_p_pi() - deltat;
            t *= 0.5;
            c_domain.binf = binf_modif + t;
            c_domain.bsup = bsup_modif - t;
        }
        deltat = c_domain.bsup - c_domain.binf;
        while c_domain.binf >= pi_p_pi() {
            c_domain.binf -= pi_p_pi();
        }
        while c_domain.binf < 0.0 {
            c_domain.binf += pi_p_pi();
        }
        c_domain.bsup = c_domain.binf + deltat;

        let l_domain = Interval::from_domain(lig_domain);

        let mut nb_sol_total: usize = 0;
        let mut solution_circle = [PeriodicInterval::default(); 4];
        let mut solution_line = [Interval::default(); 4];

        let mut c_domain_and_res = c_domain.first_intersection(&mut c_int1);
        project_on_l_and_intersect_with_l_domain(
            circle,
            line,
            &c_domain_and_res,
            &l_domain,
            &mut solution_circle,
            &mut solution_line,
            &mut nb_sol_total,
            lig_domain,
        );

        c_domain_and_res = c_domain.second_intersection(&mut c_int1);
        project_on_l_and_intersect_with_l_domain(
            circle,
            line,
            &c_domain_and_res,
            &l_domain,
            &mut solution_circle,
            &mut solution_line,
            &mut nb_sol_total,
            lig_domain,
        );

        if nbsol == 2 {
            c_domain_and_res = c_domain.first_intersection(&mut c_int2);
            project_on_l_and_intersect_with_l_domain(
                circle,
                line,
                &c_domain_and_res,
                &l_domain,
                &mut solution_circle,
                &mut solution_line,
                &mut nb_sol_total,
                lig_domain,
            );

            c_domain_and_res = c_domain.second_intersection(&mut c_int2);
            project_on_l_and_intersect_with_l_domain(
                circle,
                line,
                &c_domain_and_res,
                &l_domain,
                &mut solution_circle,
                &mut solution_line,
                &mut nb_sol_total,
                lig_domain,
            );
        }

        // Computation of all transitions and positions. Intervals whose
        // `Radius * Length` is below the tolerance are reduced to points.
        let r = circle.radius();
        let mut max_tol = tol_conf;
        if max_tol < tol {
            max_tol = tol;
        }
        if max_tol < 1.0e-10 {
            max_tol = 1.0e-10;
        }

        for i in 0..nb_sol_total {
            if (r * solution_circle[i].length()) < max_tol && solution_line[i].length() < max_tol {
                let t = (solution_circle[i].binf + solution_circle[i].bsup) * 0.5;
                solution_circle[i].binf = t;
                solution_circle[i].bsup = t;

                let t = (solution_line[i].binf + solution_line[i].bsup) * 0.5;
                solution_line[i].binf = t;
                solution_line[i].bsup = t;
            }
        }

        if nb_sol_total == 0 {
            return;
        }

        let circle_axis = *circle.position();
        let line_axis = *line.position();
        let norm2 = GpVec2d::new(0.0, 0.0);
        let mut t1a = IntRes2dTransition::new();
        let mut t2a = IntRes2dTransition::new();
        let mut t1b = IntRes2dTransition::new();
        let mut t2b = IntRes2dTransition::new();
        let mut pos1a = IntRes2dPosition::Middle;
        let mut pos1b = IntRes2dPosition::Middle;
        let mut pos2a = IntRes2dPosition::Middle;
        let mut pos2b = IntRes2dPosition::Middle;

        let (mut p1a, mut tan1) =
            clib2d::circle_d1_ax22d(solution_circle[0].binf, &circle_axis, r);
        let (mut p2a, mut tan2) = clib2d::line_d1_ax2d(solution_line[0].binf, &line_axis);

        let is_opposite = tan1.dot(&tan2) < 0.0;

        for i in 0..nb_sol_total {
            // 7 aout 97: recentre Binf/Bsup so that a common portion with
            // CIRC_Domain exists.
            let mut p1 = solution_circle[i].binf;
            let mut p2 = solution_circle[i].bsup;
            let q1 = circ_domain.first_parameter();
            let q2 = circ_domain.last_parameter();
            if p1 > q2 {
                loop {
                    p1 -= pi_p_pi();
                    p2 -= pi_p_pi();
                    if !(p1 > q2) {
                        break;
                    }
                }
            } else if p2 < q1 {
                loop {
                    p1 += pi_p_pi();
                    p2 += pi_p_pi();
                    if !(p2 < q1) {
                        break;
                    }
                }
            }
            if p1 < q1 && p2 > q1 {
                p1 = q1;
            }
            if p1 < q2 && p2 > q2 {
                p2 = q2;
            }

            solution_circle[i].binf = p1;
            solution_circle[i].bsup = p2;

            let mut linf = if is_opposite { solution_line[i].bsup } else { solution_line[i].binf };
            let mut lsup = if is_opposite { solution_line[i].binf } else { solution_line[i].bsup };

            // If the parameters on the circle come first, they must be
            // returned in increasing order.
            if linf > lsup {
                let t = solution_circle[i].binf;
                solution_circle[i].binf = solution_circle[i].bsup;
                solution_circle[i].bsup = t;

                let t = linf;
                linf = lsup;
                lsup = t;
            }

            let (pa, ta, na) =
                clib2d::circle_d2_ax22d(solution_circle[i].binf, &circle_axis, r);
            p1a = pa;
            tan1 = ta;
            let norm1 = na;
            let (pb, tb) = clib2d::line_d1_ax2d(linf, &line_axis);
            p2a = pb;
            tan2 = tb;

            pos1a = determine_position(circ_domain, &p1a, solution_circle[i].binf);
            pos2a = determine_position(lig_domain, &p2a, linf);
            determine_transition_lc(
                pos1a, &mut tan1, &norm1, &mut t1a, pos2a, &mut tan2, &norm2, &mut t2a, tol,
            );

            let mut cinf;
            if pos1a == IntRes2dPosition::End {
                cinf = circ_domain.last_parameter();
                p1a = *circ_domain.last_point();
                linf = clib2d::parameter_lin2d(line, &p1a);

                let (pa, ta, na) = clib2d::circle_d2_ax22d(cinf, &circle_axis, r);
                p1a = pa;
                tan1 = ta;
                let norm1 = na;
                let (pb, tb) = clib2d::line_d1_ax2d(linf, &line_axis);
                p2a = pb;
                tan2 = tb;
                pos1a = determine_position(circ_domain, &p1a, cinf);
                pos2a = determine_position(lig_domain, &p2a, linf);
                determine_transition_lc(
                    pos1a, &mut tan1, &norm1, &mut t1a, pos2a, &mut tan2, &norm2, &mut t2a, tol,
                );
            } else if pos1a == IntRes2dPosition::Head {
                cinf = circ_domain.first_parameter();
                p1a = *circ_domain.first_point();
                linf = clib2d::parameter_lin2d(line, &p1a);

                let (pa, ta, na) = clib2d::circle_d2_ax22d(cinf, &circle_axis, r);
                p1a = pa;
                tan1 = ta;
                let norm1 = na;
                let (pb, tb) = clib2d::line_d1_ax2d(linf, &line_axis);
                p2a = pb;
                tan2 = tb;
                pos1a = determine_position(circ_domain, &p1a, cinf);
                pos2a = determine_position(lig_domain, &p2a, linf);
                determine_transition_lc(
                    pos1a, &mut tan1, &norm1, &mut t1a, pos2a, &mut tan2, &norm2, &mut t2a, tol,
                );
            } else {
                cinf = normalize_on_circle_domain(solution_circle[i].binf, circ_domain);
            }

            let new_point1 = IntRes2dIntersectionPoint::with_transitions(
                &p1a,
                linf,
                cinf,
                &t2a,
                &t1a,
                self.result.reversed_parameters(),
            );

            if (solution_line[i].length() + solution_circle[i].length()) > 0.0 {
                let (pb, tb, nb) = clib2d::circle_d2_ax22d(solution_circle[i].bsup, &circle_axis, r);
                let p1b_pt = pb;
                let mut tan1b = tb;
                let norm1b = nb;
                let (pb2, tb2) = clib2d::line_d1_ax2d(lsup, &line_axis);
                let mut p2b = pb2;
                tan2 = tb2;

                pos1b = determine_position(
                    circ_domain,
                    &p1b_pt,
                    solution_circle[i].bsup,
                );
                pos2b = determine_position(lig_domain, &p2b, lsup);
                determine_transition_lc(
                    pos1b,
                    &mut tan1b,
                    &norm1b,
                    &mut t1b,
                    pos2b,
                    &mut tan2,
                    &norm2,
                    &mut t2b,
                    tol,
                );

                let mut p1b = p1b_pt;
                let mut csup;
                if pos1b == IntRes2dPosition::End {
                    // Faithful to OCCT: `Csup = CIRC_Domain.LastParameter()`
                    // (`_1.cxx:2555-2556`); the `Perform(Line, Ellipse)` twin at
                    // `:3159` uses `DL` instead, which is reproduced there.
                    csup = circ_domain.last_parameter();
                    p1b = *circ_domain.last_point();
                    lsup = clib2d::parameter_lin2d(line, &p1b);
                    let (pb, tb, nb) = clib2d::circle_d2_ax22d(csup, &circle_axis, r);
                    p1b = pb;
                    tan1b = tb;
                    let norm1b = nb;
                    let (pb2, tb2) = clib2d::line_d1_ax2d(lsup, &line_axis);
                    p2b = pb2;
                    tan2 = tb2;

                    pos1b = determine_position(circ_domain, &p1b, csup);
                    pos2b = determine_position(lig_domain, &p2b, lsup);
                    determine_transition_lc(
                        pos1b,
                        &mut tan1b,
                        &norm1b,
                        &mut t1b,
                        pos2b,
                        &mut tan2,
                        &norm2,
                        &mut t2b,
                        tol,
                    );
                } else if pos1b == IntRes2dPosition::Head {
                    csup = circ_domain.first_parameter();
                    p1b = *circ_domain.first_point();
                    lsup = clib2d::parameter_lin2d(line, &p1b);
                    let (pb, tb, nb) = clib2d::circle_d2_ax22d(csup, &circle_axis, r);
                    p1b = pb;
                    tan1b = tb;
                    let norm1b = nb;
                    let (pb2, tb2) = clib2d::line_d1_ax2d(lsup, &line_axis);
                    p2b = pb2;
                    tan2 = tb2;

                    pos1b = determine_position(circ_domain, &p1b, csup);
                    pos2b = determine_position(lig_domain, &p2b, lsup);
                    determine_transition_lc(
                        pos1b,
                        &mut tan1b,
                        &norm1b,
                        &mut t1b,
                        pos2b,
                        &mut tan2,
                        &norm2,
                        &mut t2b,
                        tol,
                    );
                } else {
                    csup = normalize_on_circle_domain(solution_circle[i].bsup, circ_domain);
                }

                let new_point2 = IntRes2dIntersectionPoint::with_transitions(
                    &p1b,
                    lsup,
                    csup,
                    &t2b,
                    &t1b,
                    self.result.reversed_parameters(),
                );

                if (((csup - cinf).abs() * r > max_tol) && ((lsup - linf).abs() > max_tol))
                    || (t1a.transition_type() != t2a.transition_type())
                {
                    let new_seg = IntRes2dIntersectionSegment::from_two_points(
                        &new_point1,
                        &new_point2,
                        is_opposite,
                        self.result.reversed_parameters(),
                    );
                    self.result.append_segment(&new_seg);
                } else {
                    if pos1a != IntRes2dPosition::Middle || pos2a != IntRes2dPosition::Middle {
                        self.result.insert(&new_point1);
                    }
                    if pos1b != IntRes2dPosition::Middle || pos2b != IntRes2dPosition::Middle {
                        self.result.insert(&new_point2);
                    }
                }
            } else {
                self.result.insert(&new_point1);
            }
        }
    }
}
