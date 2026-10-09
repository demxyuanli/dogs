//! Port of the `IntCurve_IntConicConic` Line/Ellipse kernel set of
//! `IntCurve_IntConicConic_1.cxx`:
//! `LineEllipseGeometricIntersection` (`:2657-2775`),
//! `ProjectOnLAndIntersectWithLDomain(const gp_Elips2d&, ...)` (`:2780-2857`)
//! and `Perform(const gp_Lin2d&, ..., const gp_Elips2d&, ...)` (`:2861-3214`).

use std::f64::consts::PI;

use occt_core::elib::clib2d;
use occt_core::gp::{GpElips2d, GpLin2d, GpPnt2d, GpTrsf2d, GpVec2d};
use occt_core::intimpargen::gen::determine_position;
use occt_core::intres2d::{
    IntRes2dDomain, IntRes2dIntersectionPoint, IntRes2dIntersectionSegment, IntRes2dPosition,
    IntRes2dTransition,
};
use occt_core::precision::{epsilon, INFINITE};

use super::int_conic_conic::IntCurveIntConicConic;
use super::int_conic_conic_tool::{
    determine_transition_lc, normalize_on_circle_domain, pi_p_pi, Interval, PeriodicInterval,
};
use crate::extrema2d::ExtremaExtElC2d;

/// `LineEllipseGeometricIntersection` (`_1.cxx:2657-2775`). Returns `nbsol` and
/// writes the two solution intervals on the ellipse.
pub(crate) fn line_ellipse_geometric_intersection(
    line: &GpLin2d,
    ellipse: &GpElips2d,
    _tol_conf: f64,
    tol_tang: f64,
    e_int1: &mut PeriodicInterval,
    e_int2: &mut PeriodicInterval,
) -> i32 {
    let an_el_axis = ellipse.pos;
    let mut a_tr = GpTrsf2d::identity();
    a_tr.set_transformation_ax2d(&an_el_axis.x_axis());
    let a_t_ellipse = ellipse.transformed(&a_tr);
    let a_t_line = line.transformed(&a_tr);
    let a_dy = a_t_line.pos.vdir.y;
    let is_vert = a_dy.abs() > 1.0 - 2.0 * epsilon(1.0);

    let a = a_t_ellipse.major_radius;
    let b = a_t_ellipse.minor_radius;
    let a2 = a * a;
    let b2 = b * b;

    let mut eps0 = 1.0e-12;
    if b / a < 1.0e-5 {
        eps0 = 1.0e-6;
    }

    let (an_a, mut a_b, mut a_c) = a_t_line.coefficients();
    if is_vert {
        a_c += a_b * a_t_line.pos.loc.y();
        a_b = 0.0;
    }

    let mut x1 = 0.0;
    let mut y1 = 0.0;
    let mut x2 = 0.0;
    let mut y2 = 0.0;
    let nbsol;

    if a_b.abs() > eps0 {
        let m = -an_a / a_b;
        let m2 = m * m;
        let c = -a_c / a_b;
        let c2 = c * c;
        let mut d = a2 * m2 + b2 - c2;
        if d < 0.0 {
            // `_1.cxx:2702-2727`: the line misses the ellipse. `Extrema_ExtElC2d`
            // on the transformed pair gives the two stationary pairs; when the
            // closest one is within `TolTang` of the line the common point is
            // reported as a single tangent solution on the ellipse.
            let an_ext = ExtremaExtElC2d::new_line_ellipse(&a_t_line, &a_t_ellipse);
            let mut imin = 0usize;
            let mut dmin = f64::MAX;
            for i in 1..=an_ext.nb_ext() {
                if an_ext.square_distance(i) < dmin {
                    dmin = an_ext.square_distance(i);
                    imin = i;
                }
            }
            if imin > 0 && dmin <= tol_tang * tol_tang {
                let (_, a_p2) = an_ext.points(imin);
                let pe1 = a_p2.parameter();
                e_int1.set_values(pe1, pe1);
                return 1;
            }
            return 0;
        }
        d = d.sqrt();
        let n = a2 * m2 + b2;
        let k = a * b * d / n;
        let l = -a2 * m * c / n;
        x1 = l + k;
        y1 = m * x1 + c;
        x2 = l - k;
        y2 = m * x2 + c;
        nbsol = 2;
    } else {
        x1 = -a_c / an_a;
        if x1.abs() > a + tol_tang {
            return 0;
        } else if x1.abs() >= a - epsilon(1.0 + a) {
            nbsol = 1;
            y1 = 0.0;
        } else {
            y1 = b * (1.0 - x1 * x1 / a2).sqrt();
            x2 = x1;
            y2 = -y1;
            nbsol = 2;
        }
    }

    let a_p1 = GpPnt2d::new(x1, y1);
    let a_p2 = GpPnt2d::new(x2, y2);
    let mut pe1 = clib2d::parameter_elips2d(&a_t_ellipse, &a_p1);
    if nbsol > 1 {
        let mut pe2 = clib2d::parameter_elips2d(&a_t_ellipse, &a_p2);
        if pe2 < pe1 {
            let t = pe1;
            pe1 = pe2;
            pe2 = t;
        }
        e_int2.set_values(pe2, pe2);
    }
    e_int1.set_values(pe1, pe1);

    nbsol
}

/// `ProjectOnLAndIntersectWithLDomain(const gp_Elips2d&, ...)`
/// (`_1.cxx:2780-2857`).
#[allow(clippy::too_many_arguments)]
fn project_on_l_and_intersect_with_l_domain(
    ellipse: &GpElips2d,
    line: &GpLin2d,
    e_domain_and_res: &PeriodicInterval,
    l_domain: &Interval,
    ellipse_solution: &mut [PeriodicInterval; 4],
    line_solution: &mut [Interval; 4],
    nb_sol_total: &mut usize,
    ref_line_domain: &IntRes2dDomain,
) {
    if e_domain_and_res.is_null() {
        return;
    }

    let linf = clib2d::parameter_lin2d(
        line,
        &clib2d::ellipse_value_ax22d(
            e_domain_and_res.binf,
            &ellipse.pos,
            ellipse.major_radius,
            ellipse.minor_radius,
        ),
    );
    let lsup = clib2d::parameter_lin2d(
        line,
        &clib2d::ellipse_value_ax22d(
            e_domain_and_res.bsup,
            &ellipse.pos,
            ellipse.major_radius,
            ellipse.minor_radius,
        ),
    );

    // Necessarily bounded.
    let l_inter = Interval::from_bounds(linf, lsup);
    let mut l_inter_and_domain = l_domain.intersection_with_bounded(&l_inter);

    if !l_inter_and_domain.is_null {
        let dom_linf = if ref_line_domain.has_first_point() {
            ref_line_domain.first_parameter()
        } else {
            -INFINITE
        };
        let dom_lsup = if ref_line_domain.has_last_point() {
            ref_line_domain.last_parameter()
        } else {
            INFINITE
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

        let mut einf = e_domain_and_res.binf;
        let mut esup = e_domain_and_res.bsup;
        if einf >= esup {
            einf = e_domain_and_res.binf;
            esup = e_domain_and_res.bsup;
        }
        ellipse_solution[*nb_sol_total] = PeriodicInterval::from_bounds(einf, esup);
        if ellipse_solution[*nb_sol_total].length() > PI {
            ellipse_solution[*nb_sol_total].complement();
        }

        line_solution[*nb_sol_total] = l_inter_and_domain;
        *nb_sol_total += 1;
    }
}

impl IntCurveIntConicConic {
    /// `Perform(const gp_Lin2d& L, const IntRes2d_Domain& DL,
    /// const gp_Elips2d& E, const IntRes2d_Domain& DE, TolConf, Tol)`
    /// (`_1.cxx:2861-3210`).
    pub fn perform_line_ellipse(
        &mut self,
        l: &GpLin2d,
        dl: &IntRes2dDomain,
        e: &GpElips2d,
        de: &IntRes2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        let the_reversed_parameters = self.result.reversed_parameters();
        self.result.reset_fields();
        self.result.set_reversed_parameters(the_reversed_parameters);

        let mut e_int1 = PeriodicInterval::default();
        let mut e_int2 = PeriodicInterval::default();

        let mut nbsol =
            line_ellipse_geometric_intersection(l, e, tol_conf, tol, &mut e_int1, &mut e_int2);

        self.result.done = true;

        if nbsol == 0 {
            return;
        }

        if nbsol == 2 && e_int2.bsup == e_int1.binf + pi_p_pi() {
            let first_bound = de.first_parameter();
            let last_bound = de.last_parameter();
            let first_tol = de.first_tolerance();
            let last_tol = de.last_tolerance();
            if e_int1.binf == 0.0 && first_bound - first_tol > e_int1.bsup {
                nbsol = 1;
                e_int1.set_values(e_int2.binf, e_int2.bsup);
            } else if e_int2.bsup == pi_p_pi() && last_bound + last_tol < e_int2.binf {
                nbsol = 1;
            }
        }

        let mut e_domain = PeriodicInterval::from_domain(de);
        let mut deltat = e_domain.bsup - e_domain.binf;
        while e_domain.binf >= pi_p_pi() {
            e_domain.binf -= pi_p_pi();
        }
        while e_domain.binf < 0.0 {
            e_domain.binf += pi_p_pi();
        }
        e_domain.bsup = e_domain.binf + deltat;

        let mut binf_modif = e_domain.binf;
        let mut bsup_modif = e_domain.bsup;
        binf_modif -= de.first_tolerance() / e.minor_radius;
        bsup_modif += de.last_tolerance() / e.minor_radius;
        deltat = bsup_modif - binf_modif;
        if deltat <= pi_p_pi() {
            e_domain.binf = binf_modif;
            e_domain.bsup = bsup_modif;
        } else {
            let mut t = pi_p_pi() - deltat;
            t *= 0.5;
            e_domain.binf = binf_modif + t;
            e_domain.bsup = bsup_modif - t;
        }
        deltat = e_domain.bsup - e_domain.binf;
        while e_domain.binf >= pi_p_pi() {
            e_domain.binf -= pi_p_pi();
        }
        while e_domain.binf < 0.0 {
            e_domain.binf += pi_p_pi();
        }
        e_domain.bsup = e_domain.binf + deltat;

        let l_domain = Interval::from_domain(dl);

        let mut nb_sol_total: usize = 0;
        let mut solution_ellipse = [PeriodicInterval::default(); 4];
        let mut solution_line = [Interval::default(); 4];

        let mut e_domain_and_res = e_domain.first_intersection(&mut e_int1);
        project_on_l_and_intersect_with_l_domain(
            e,
            l,
            &e_domain_and_res,
            &l_domain,
            &mut solution_ellipse,
            &mut solution_line,
            &mut nb_sol_total,
            dl,
        );

        e_domain_and_res = e_domain.second_intersection(&mut e_int1);
        project_on_l_and_intersect_with_l_domain(
            e,
            l,
            &e_domain_and_res,
            &l_domain,
            &mut solution_ellipse,
            &mut solution_line,
            &mut nb_sol_total,
            dl,
        );

        if nbsol == 2 {
            e_domain_and_res = e_domain.first_intersection(&mut e_int2);
            project_on_l_and_intersect_with_l_domain(
                e,
                l,
                &e_domain_and_res,
                &l_domain,
                &mut solution_ellipse,
                &mut solution_line,
                &mut nb_sol_total,
                dl,
            );

            e_domain_and_res = e_domain.second_intersection(&mut e_int2);
            project_on_l_and_intersect_with_l_domain(
                e,
                l,
                &e_domain_and_res,
                &l_domain,
                &mut solution_ellipse,
                &mut solution_line,
                &mut nb_sol_total,
                dl,
            );
        }

        // Calculation of transitions at positions.
        let r = e.minor_radius;
        let mut max_tol = tol_conf;
        if max_tol < tol {
            max_tol = tol;
        }
        if max_tol < 1.0e-10 {
            max_tol = 1.0e-10;
        }

        for i in 0..nb_sol_total {
            if (r * solution_ellipse[i].length()) < max_tol && solution_line[i].length() < max_tol {
                let t = (solution_ellipse[i].binf + solution_ellipse[i].bsup) * 0.5;
                solution_ellipse[i].binf = t;
                solution_ellipse[i].bsup = t;

                let t = (solution_line[i].binf + solution_line[i].bsup) * 0.5;
                solution_line[i].binf = t;
                solution_line[i].bsup = t;
            }
        }

        if nb_sol_total == 0 {
            return;
        }

        let ellipse_axis = e.pos;
        let line_axis = l.pos;
        let norm2 = GpVec2d::new(0.0, 0.0);
        let mut t1a = IntRes2dTransition::new();
        let mut t2a = IntRes2dTransition::new();
        let mut t1b = IntRes2dTransition::new();
        let mut t2b = IntRes2dTransition::new();
        let mut pos1a;
        let mut pos1b;
        let mut pos2a;
        let mut pos2b;

        let (mut p1a, mut tan1) = clib2d::ellipse_d1_ax22d(
            solution_ellipse[0].binf,
            &ellipse_axis,
            e.major_radius,
            e.minor_radius,
        );
        let (mut p2a, mut tan2) = clib2d::line_d1_ax2d(solution_line[0].binf, &line_axis);

        let is_opposite = tan1.dot(&tan2) < 0.0;

        for i in 0..nb_sol_total {
            let mut p1 = solution_ellipse[i].binf;
            let mut p2 = solution_ellipse[i].bsup;
            let q1 = de.first_parameter();
            let q2 = de.last_parameter();

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

            solution_ellipse[i].binf = p1;
            solution_ellipse[i].bsup = p2;

            let mut linf = if is_opposite { solution_line[i].bsup } else { solution_line[i].binf };
            let mut lsup = if is_opposite { solution_line[i].binf } else { solution_line[i].bsup };

            if linf > lsup {
                let t = solution_ellipse[i].binf;
                solution_ellipse[i].binf = solution_ellipse[i].bsup;
                solution_ellipse[i].bsup = t;
                let t = linf;
                linf = lsup;
                lsup = t;
            }

            let (pa, ta, na) = clib2d::ellipse_d2_ax22d(
                solution_ellipse[i].binf,
                &ellipse_axis,
                e.major_radius,
                e.minor_radius,
            );
            p1a = pa;
            tan1 = ta;
            let norm1 = na;
            let (pb, tb) = clib2d::line_d1_ax2d(linf, &line_axis);
            p2a = pb;
            tan2 = tb;

            pos1a = determine_position(de, &p1a, solution_ellipse[i].binf);
            pos2a = determine_position(dl, &p2a, linf);
            determine_transition_lc(
                pos1a, &mut tan1, &norm1, &mut t1a, pos2a, &mut tan2, &norm2, &mut t2a, tol,
            );

            let mut einf;
            if pos1a == IntRes2dPosition::End {
                einf = de.last_parameter();
                p1a = *de.last_point();
                linf = clib2d::parameter_lin2d(l, &p1a);

                let (pa, ta, na) = clib2d::ellipse_d2_ax22d(
                    einf,
                    &ellipse_axis,
                    e.major_radius,
                    e.minor_radius,
                );
                p1a = pa;
                tan1 = ta;
                let norm1 = na;
                let (pb, tb) = clib2d::line_d1_ax2d(linf, &line_axis);
                p2a = pb;
                tan2 = tb;
                pos1a = determine_position(de, &p1a, einf);
                pos2a = determine_position(dl, &p2a, linf);
                determine_transition_lc(
                    pos1a, &mut tan1, &norm1, &mut t1a, pos2a, &mut tan2, &norm2, &mut t2a, tol,
                );
            } else if pos1a == IntRes2dPosition::Head {
                einf = de.first_parameter();
                p1a = *de.first_point();
                linf = clib2d::parameter_lin2d(l, &p1a);

                let (pa, ta, na) = clib2d::ellipse_d2_ax22d(
                    einf,
                    &ellipse_axis,
                    e.major_radius,
                    e.minor_radius,
                );
                p1a = pa;
                tan1 = ta;
                let norm1 = na;
                let (pb, tb) = clib2d::line_d1_ax2d(linf, &line_axis);
                p2a = pb;
                tan2 = tb;
                pos1a = determine_position(de, &p1a, einf);
                pos2a = determine_position(dl, &p2a, linf);
                determine_transition_lc(
                    pos1a, &mut tan1, &norm1, &mut t1a, pos2a, &mut tan2, &norm2, &mut t2a, tol,
                );
            } else {
                einf = normalize_on_circle_domain(solution_ellipse[i].binf, de);
            }

            let new_point1 = IntRes2dIntersectionPoint::with_transitions(
                &p1a,
                linf,
                einf,
                &t2a,
                &t1a,
                self.result.reversed_parameters(),
            );

            if (solution_line[i].length() + solution_ellipse[i].length()) > 0.0 {
                // Faithful to OCCT `_1.cxx:3144`, which evaluates `EllipseD2` at
                // `SolutionEllipse[i].Binf` (the Line/Circle twin at `:2546`
                // uses `Bsup`), while `DeterminePosition` below uses `Bsup`.
                let (pb, tb, nb) = clib2d::ellipse_d2_ax22d(
                    solution_ellipse[i].binf,
                    &ellipse_axis,
                    e.major_radius,
                    e.minor_radius,
                );
                let p1b_pt = pb;
                let mut tan1b = tb;
                let norm1b = nb;
                let (pb2, tb2) = clib2d::line_d1_ax2d(lsup, &line_axis);
                let mut p2b = pb2;
                tan2 = tb2;

                pos1b = determine_position(
                    de,
                    &p1b_pt,
                    solution_ellipse[i].bsup,
                );
                pos2b = determine_position(dl, &p2b, lsup);
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
                let mut esup;
                if pos1b == IntRes2dPosition::End {
                    // Faithful to OCCT `_1.cxx:3159`, which reads `DL` rather
                    // than `DE` here (unlike the Line/Circle twin at `:2556`).
                    esup = dl.last_parameter();
                    p1b = *de.last_point();
                    lsup = clib2d::parameter_lin2d(l, &p1b);
                    let (pb, tb, nb) = clib2d::ellipse_d2_ax22d(
                        esup,
                        &ellipse_axis,
                        e.major_radius,
                        e.minor_radius,
                    );
                    p1b = pb;
                    tan1b = tb;
                    let norm1b = nb;
                    let (pb2, tb2) = clib2d::line_d1_ax2d(lsup, &line_axis);
                    p2b = pb2;
                    tan2 = tb2;

                    pos1b = determine_position(de, &p1b, esup);
                    pos2b = determine_position(dl, &p2b, lsup);
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
                    esup = de.first_parameter();
                    p1b = *de.first_point();
                    lsup = clib2d::parameter_lin2d(l, &p1b);
                    let (pb, tb, nb) = clib2d::ellipse_d2_ax22d(
                        esup,
                        &ellipse_axis,
                        e.major_radius,
                        e.minor_radius,
                    );
                    p1b = pb;
                    tan1b = tb;
                    let norm1b = nb;
                    let (pb2, tb2) = clib2d::line_d1_ax2d(lsup, &line_axis);
                    p2b = pb2;
                    tan2 = tb2;

                    pos1b = determine_position(de, &p1b, esup);
                    pos2b = determine_position(dl, &p2b, lsup);
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
                    esup = normalize_on_circle_domain(solution_ellipse[i].bsup, de);
                }

                let new_point2 = IntRes2dIntersectionPoint::with_transitions(
                    &p1b,
                    lsup,
                    esup,
                    &t2b,
                    &t1b,
                    self.result.reversed_parameters(),
                );

                if (((esup - einf).abs() * r > max_tol) && ((lsup - linf).abs() > max_tol))
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
