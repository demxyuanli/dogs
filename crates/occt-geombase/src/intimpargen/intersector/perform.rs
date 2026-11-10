use crate::intcurve::iconic_tool::IntCurveIConicTool;
use crate::intcurve::pconic::IntCurvePConic;
use crate::intcurve::pconic_tool;
use occt_gp::gp::GpPnt2d;
use occt_gp::intres2d::{
    IntRes2dDomain, IntRes2dIntersectionPoint, IntRes2dIntersectionSegment, IntRes2dPosition,
};
use occt_gp::math_function_all_roots::FunctionAllRoots;
use occt_gp::math_function_sample::FunctionSample;

use super::super::gen;
use super::imp_par_tool::MyImpParTool;
use super::structure::IntImpParGenIntersector;

impl IntImpParGenIntersector {    /// `Perform` (`gxx:245-779`).
    pub fn perform(
        &mut self,
        imp_tool: &IntCurveIConicTool,
        imp_curve_domain: &IntRes2dDomain,
        par_curve: &IntCurvePConic,
        par_curve_domain: &IntRes2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        let mut head_on_imp = false;
        let mut head_on_par = false;
        let mut end_on_imp = false;
        let mut end_on_par = false;

        self.base.reset_fields();

        let mut imp_par_tool = MyImpParTool::new(imp_tool, par_curve);

        if !(par_curve_domain.has_first_point() && par_curve_domain.has_last_point()) {
            panic!("Standard_ConstructionError: Domaine sur courbe incorrect");
        }

        let nb_echantillons = pconic_tool::nb_samples_in_range(
            par_curve,
            par_curve_domain.first_parameter(),
            par_curve_domain.last_parameter(),
        );

        let mut eps_x = pconic_tool::pconic_eps_x(par_curve);
        if eps_x > 1.0e-10 {
            eps_x = 1.0e-10;
        }
        let eps_nul = if tol_conf <= 1.0e-10 { 1.0e-10 } else { tol_conf };
        let eps_dist = if tol <= 1.0e-10 { 1.0e-10 } else { tol };

        let tolerance_angulaire = eps_dist;

        if (par_curve_domain.last_parameter() - par_curve_domain.first_parameter())
            < 100.0 * eps_x
        {
            eps_x = (par_curve_domain.last_parameter()
                - par_curve_domain.first_parameter())
                * 0.01;
        }

        let sample2 = FunctionSample::new(
            par_curve_domain.first_parameter(),
            par_curve_domain.last_parameter(),
            nb_echantillons,
        );

        let sol = FunctionAllRoots::new(&mut imp_par_tool, &sample2, eps_x, eps_dist, eps_nul);

        if !sol.is_done() {
            self.base.done = false;
            return;
        }

        let nb_segments_solution = sol.nb_intervals();
        let nb_points_solution = sol.nb_points();

        // Traitement des Points Solutions
        for i in 1..=nb_points_solution {
            let param2 = sol.point(i);
            let mut pt = GpPnt2d::new(0.0, 0.0);
            let mut param1 = Self::find_u(param2, &mut pt, par_curve, imp_tool);

            if imp_curve_domain.is_closed() {
                param1 = gen::normalize_on_domain(param1, imp_curve_domain);
            }

            let mut is_on_the_imp_curve_domain = true;
            if imp_curve_domain.has_first_point() {
                if param1 < imp_curve_domain.first_parameter() {
                    if pt.distance(imp_curve_domain.first_point())
                        > imp_curve_domain.first_tolerance()
                    {
                        is_on_the_imp_curve_domain = false;
                    }
                }
            }
            if is_on_the_imp_curve_domain && imp_curve_domain.has_last_point() {
                if param1 > imp_curve_domain.last_parameter() {
                    if pt.distance(imp_curve_domain.last_point())
                        > imp_curve_domain.last_tolerance()
                    {
                        is_on_the_imp_curve_domain = false;
                    }
                }
            }

            if is_on_the_imp_curve_domain {
                let (pt1, mut tan1, norm1) = imp_tool.d2(param1);
                let (pt2, mut tan2, norm2) = pconic_tool::pconic_d2(par_curve, param2);

                let pos1 = gen::determine_position(imp_curve_domain, &pt1, param1);
                let pos2 = gen::determine_position(par_curve_domain, &pt2, param2);

                if pos1 == IntRes2dPosition::End {
                    end_on_imp = true;
                } else if pos1 == IntRes2dPosition::Head {
                    head_on_imp = true;
                }
                if pos2 == IntRes2dPosition::End {
                    end_on_par = true;
                } else if pos2 == IntRes2dPosition::Head {
                    head_on_par = true;
                }

                let mut trans1 = occt_gp::intres2d::IntRes2dTransition::new();
                let mut trans2 = occt_gp::intres2d::IntRes2dTransition::new();
                gen::determine_transition_touch(
                    pos1,
                    &mut tan1,
                    &norm1,
                    &mut trans1,
                    pos2,
                    &mut tan2,
                    &norm2,
                    &mut trans2,
                    tolerance_angulaire,
                );

                let ip = IntRes2dIntersectionPoint::with_transitions(
                    &pt1,
                    param1,
                    param2,
                    &trans1,
                    &trans2,
                    self.base.reversed_parameters(),
                );
                self.base.insert(&ip);
            }
        }

        // Traitement des Segments
        let mut inter2_and_domaine2 = vec![0.0f64; (2 + 8 * nb_segments_solution) as usize + 2];
        let mut inter1 = vec![0.0f64; (2 + 8 * nb_segments_solution) as usize + 2];

        let mut nb_segments_crees = 0;

        let mut j2 = 1;
        for j in 1..=nb_segments_solution {
            let mut ptemp = GpPnt2d::new(0.0, 0.0);
            let (param2_inf, param2_sup) = sol.interval(j);
            let mut param1_inf = Self::find_u(param2_inf, &mut ptemp, par_curve, imp_tool);
            let mut param1_sup = Self::find_u(param2_sup, &mut ptemp, par_curve, imp_tool);

            if imp_curve_domain.is_closed() {
                let (param1_origine, param1_fin) = imp_curve_domain.equivalent_parameters();
                let periode = param1_fin - param1_origine;

                while param1_inf < param1_origine {
                    param1_inf += periode;
                }
                while param1_sup < param1_origine {
                    param1_sup += periode;
                }

                let (_p2, mut t2, n2) = pconic_tool::pconic_d2(par_curve, param2_inf);
                let (_p1, mut t1, n1) = imp_tool.d2(param1_inf);
                if t1.magnitude() <= occt_foundation::precision::RESOLUTION {
                    t1 = n1;
                }
                if t2.magnitude() <= occt_foundation::precision::RESOLUTION {
                    t2 = n2;
                }

                if t1.dot(&t2) >= 0.0 {
                    // param1_inf designe un point entrant
                    if param1_inf >= param1_sup {
                        param1_sup += periode;
                    }
                } else {
                    // param1_inf : point sortant
                    if param1_inf <= param1_sup {
                        param1_inf += periode;
                    }
                }

                if imp_curve_domain.last_parameter()
                    > (if param1_inf > param1_sup {
                        param1_sup + periode
                    } else {
                        param1_inf + periode
                    })
                {
                    inter2_and_domaine2[j2 as usize] = param2_inf;
                    inter1[j2 as usize] = param1_inf + periode;
                    inter2_and_domaine2[(j2 + 1) as usize] = param2_sup;
                    inter1[(j2 + 1) as usize] = param1_sup + periode;
                    j2 += 2;
                    nb_segments_crees += 1;
                }

                if imp_curve_domain.first_parameter()
                    < (if param1_inf < param1_sup {
                        param1_sup - periode
                    } else {
                        param1_inf - periode
                    })
                {
                    inter2_and_domaine2[j2 as usize] = param2_inf;
                    inter1[j2 as usize] = param1_inf - periode;
                    inter2_and_domaine2[(j2 + 1) as usize] = param2_sup;
                    inter1[(j2 + 1) as usize] = param1_sup - periode;
                    j2 += 2;
                    nb_segments_crees += 1;
                }
            }

            inter2_and_domaine2[j2 as usize] = param2_inf;
            inter1[j2 as usize] = param1_inf;
            inter2_and_domaine2[(j2 + 1) as usize] = param2_sup;
            inter1[(j2 + 1) as usize] = param1_sup;
            j2 += 2;
        }

        let mut resultat1: Vec<f64> = Vec::new();
        let mut resultat2: Vec<f64> = Vec::new();
        let nb_segments_solution = nb_segments_solution + nb_segments_crees;
        let mut nb_resultats = nb_segments_solution * 2;

        Self::and_domaine_objet1_intersections(
            imp_tool,
            imp_curve_domain,
            par_curve,
            par_curve_domain,
            &mut nb_resultats,
            &inter2_and_domaine2,
            &inter1,
            &mut resultat1,
            &mut resultat2,
            eps_nul,
        );

        // Calcule_Toutes_Transitions (inlined in OCCT, `gxx:522-658`)
        {
            let dist_mini_imp_curve = eps_nul;
            let tolerance_angulaire_dist_mini = dist_mini_imp_curve;

            let mut k = 1;
            while k <= nb_resultats {
                let ip1 = k + 1;
                let mut only_one_point = false;

                let mut param1_on1 = resultat1[(k - 1) as usize];
                let mut param1_on2 = resultat2[(k - 1) as usize];
                let param2_on1 = resultat1[(ip1 - 1) as usize];
                let param2_on2 = resultat2[(ip1 - 1) as usize];

                let pt1_on1 = imp_tool.value(param1_on1);
                let pt2_on1 = imp_tool.value(param2_on1);
                let pt1_on2 = pconic_tool::pconic_value(par_curve, param1_on2);
                let pt2_on2 = pconic_tool::pconic_value(par_curve, param2_on2);

                if !imp_curve_domain.is_closed() {
                    if pt1_on1.distance(&pt2_on1) <= dist_mini_imp_curve
                        && pt1_on2.distance(&pt2_on2) <= dist_mini_imp_curve
                    {
                        only_one_point = true;
                    }
                }

                param1_on1 = gen::normalize_on_domain(param1_on1, imp_curve_domain);
                param1_on2 = gen::normalize_on_domain(param1_on2, par_curve_domain);

                let (mut pt1_on1, mut tan1, norm1) = imp_tool.d2(param1_on1);
                let (_pt1_on2, mut tan2, norm2) = pconic_tool::pconic_d2(par_curve, param1_on2);

                let pos1 = gen::determine_position(imp_curve_domain, &pt1_on1, param1_on1);
                let pos2 = gen::determine_position(par_curve_domain, &pt1_on2, param1_on2);

                if pos1 == IntRes2dPosition::End {
                    end_on_imp = true;
                } else if pos1 == IntRes2dPosition::Head {
                    head_on_imp = true;
                }
                if pos2 == IntRes2dPosition::End {
                    end_on_par = true;
                } else if pos2 == IntRes2dPosition::Head {
                    head_on_par = true;
                }

                let mut trans1 = occt_gp::intres2d::IntRes2dTransition::new();
                let mut trans2 = occt_gp::intres2d::IntRes2dTransition::new();
                gen::determine_transition_touch(
                    pos1,
                    &mut tan1,
                    &norm1,
                    &mut trans1,
                    pos2,
                    &mut tan2,
                    &norm2,
                    &mut trans2,
                    tolerance_angulaire_dist_mini,
                );

                // Detection du cas : l'intersection est en bout sur les 2 domaines
                if pos1 != IntRes2dPosition::Middle && pos2 != IntRes2dPosition::Middle {
                    let m = 0.5 * (pt1_on1.x() + pt1_on2.x());
                    pt1_on1.set_x(m);
                    let m = 0.5 * (pt1_on1.y() + pt1_on2.y());
                    pt1_on1.set_y(m);
                }

                let new_p1 = IntRes2dIntersectionPoint::with_transitions(
                    &pt1_on1,
                    param1_on1,
                    param1_on2,
                    &trans1,
                    &trans2,
                    self.base.reversed_parameters(),
                );

                if !only_one_point {
                    let mut new_p2 = IntRes2dIntersectionPoint::new();

                    let param2_on1 = gen::normalize_on_domain(param2_on1, imp_curve_domain);
                    let param2_on2 = gen::normalize_on_domain(param2_on2, par_curve_domain);

                    let (mut pt2_on1, mut tan1, norm1) = imp_tool.d2(param2_on1);
                    let (pt2_on2, mut tan2, norm2) = pconic_tool::pconic_d2(par_curve, param2_on2);

                    let pos1 = gen::determine_position(imp_curve_domain, &pt2_on1, param2_on1);
                    let pos2 = gen::determine_position(par_curve_domain, &pt2_on2, param2_on2);

                    if pos1 == IntRes2dPosition::End {
                        end_on_imp = true;
                    } else if pos1 == IntRes2dPosition::Head {
                        head_on_imp = true;
                    }
                    if pos2 == IntRes2dPosition::End {
                        end_on_par = true;
                    } else if pos2 == IntRes2dPosition::Head {
                        head_on_par = true;
                    }

                    let mut trans1 = occt_gp::intres2d::IntRes2dTransition::new();
                    let mut trans2 = occt_gp::intres2d::IntRes2dTransition::new();
                    gen::determine_transition_touch(
                        pos1,
                        &mut tan1,
                        &norm1,
                        &mut trans1,
                        pos2,
                        &mut tan2,
                        &norm2,
                        &mut trans2,
                        tolerance_angulaire_dist_mini,
                    );

                    if pos1 != IntRes2dPosition::Middle && pos2 != IntRes2dPosition::Middle {
                        let m = 0.5 * (pt2_on1.x() + pt2_on2.x());
                        pt2_on1.set_x(m);
                        let m = 0.5 * (pt2_on1.y() + pt2_on2.y());
                        pt2_on1.set_y(m);
                    }

                    new_p2.set_values(
                        &pt2_on1,
                        param2_on1,
                        param2_on2,
                        &trans1,
                        &trans2,
                        self.base.reversed_parameters(),
                    );

                    let segopposite = tan1.dot(&tan2) < 0.0;

                    let new_seg = IntRes2dIntersectionSegment::from_two_points(
                        &new_p1,
                        &new_p2,
                        segopposite,
                        self.base.reversed_parameters(),
                    );
                    self.base.append_segment(&new_seg);
                } else {
                    self.base.insert(&new_p1);
                }

                k += 2;
            }
        }

        // On teste les points en bouts solutions
        if !head_on_imp && imp_curve_domain.has_first_point() {
            if !head_on_par
                && imp_curve_domain
                    .first_point()
                    .distance(par_curve_domain.first_point())
                    <= imp_curve_domain
                        .first_tolerance()
                        .max(par_curve_domain.first_tolerance())
            {
                let param1 = imp_curve_domain.first_parameter();
                let param2 = par_curve_domain.first_parameter();
                let (_pt1, mut tan1, norm1) = imp_tool.d2(param1);
                let (_pt2, mut tan2, norm2) = pconic_tool::pconic_d2(par_curve, param2);
                let mut trans1 = occt_gp::intres2d::IntRes2dTransition::new();
                let mut trans2 = occt_gp::intres2d::IntRes2dTransition::new();
                gen::determine_transition_touch(
                    IntRes2dPosition::Head,
                    &mut tan1,
                    &norm1,
                    &mut trans1,
                    IntRes2dPosition::Head,
                    &mut tan2,
                    &norm2,
                    &mut trans2,
                    tolerance_angulaire,
                );
                let ip = IntRes2dIntersectionPoint::with_transitions(
                    imp_curve_domain.first_point(),
                    param1,
                    param2,
                    &trans1,
                    &trans2,
                    self.base.reversed_parameters(),
                );
                self.base.insert(&ip);
            }
            if !end_on_par
                && imp_curve_domain
                    .first_point()
                    .distance(par_curve_domain.last_point())
                    <= imp_curve_domain
                        .first_tolerance()
                        .max(par_curve_domain.last_tolerance())
            {
                let param1 = imp_curve_domain.first_parameter();
                let param2 = par_curve_domain.last_parameter();
                let (_pt1, mut tan1, norm1) = imp_tool.d2(param1);
                let (_pt2, mut tan2, norm2) = pconic_tool::pconic_d2(par_curve, param2);
                let mut trans1 = occt_gp::intres2d::IntRes2dTransition::new();
                let mut trans2 = occt_gp::intres2d::IntRes2dTransition::new();
                gen::determine_transition_touch(
                    IntRes2dPosition::Head,
                    &mut tan1,
                    &norm1,
                    &mut trans1,
                    IntRes2dPosition::End,
                    &mut tan2,
                    &norm2,
                    &mut trans2,
                    tolerance_angulaire,
                );
                let ip = IntRes2dIntersectionPoint::with_transitions(
                    imp_curve_domain.first_point(),
                    param1,
                    param2,
                    &trans1,
                    &trans2,
                    self.base.reversed_parameters(),
                );
                self.base.insert(&ip);
            }
        }

        if !end_on_imp && imp_curve_domain.has_last_point() {
            if !head_on_par
                && imp_curve_domain
                    .last_point()
                    .distance(par_curve_domain.first_point())
                    <= imp_curve_domain
                        .last_tolerance()
                        .max(par_curve_domain.first_tolerance())
            {
                let param1 = imp_curve_domain.last_parameter();
                let param2 = par_curve_domain.first_parameter();
                let (_pt1, mut tan1, norm1) = imp_tool.d2(param1);
                let (_pt2, mut tan2, norm2) = pconic_tool::pconic_d2(par_curve, param2);
                let mut trans1 = occt_gp::intres2d::IntRes2dTransition::new();
                let mut trans2 = occt_gp::intres2d::IntRes2dTransition::new();
                gen::determine_transition_touch(
                    IntRes2dPosition::End,
                    &mut tan1,
                    &norm1,
                    &mut trans1,
                    IntRes2dPosition::Head,
                    &mut tan2,
                    &norm2,
                    &mut trans2,
                    tolerance_angulaire,
                );
                let ip = IntRes2dIntersectionPoint::with_transitions(
                    imp_curve_domain.last_point(),
                    param1,
                    param2,
                    &trans1,
                    &trans2,
                    self.base.reversed_parameters(),
                );
                self.base.insert(&ip);
            }
            if !end_on_par
                && imp_curve_domain
                    .last_point()
                    .distance(par_curve_domain.last_point())
                    <= imp_curve_domain
                        .last_tolerance()
                        .max(par_curve_domain.last_tolerance())
            {
                let param1 = imp_curve_domain.last_parameter();
                let param2 = par_curve_domain.last_parameter();
                let (_pt1, mut tan1, norm1) = imp_tool.d2(param1);
                let (_pt2, mut tan2, norm2) = pconic_tool::pconic_d2(par_curve, param2);
                let mut trans1 = occt_gp::intres2d::IntRes2dTransition::new();
                let mut trans2 = occt_gp::intres2d::IntRes2dTransition::new();
                gen::determine_transition_touch(
                    IntRes2dPosition::End,
                    &mut tan1,
                    &norm1,
                    &mut trans1,
                    IntRes2dPosition::End,
                    &mut tan2,
                    &norm2,
                    &mut trans2,
                    tolerance_angulaire,
                );
                let ip = IntRes2dIntersectionPoint::with_transitions(
                    imp_curve_domain.last_point(),
                    param1,
                    param2,
                    &trans1,
                    &trans2,
                    self.base.reversed_parameters(),
                );
                self.base.insert(&ip);
            }
        }

        self.base.done = true;
    }
}
