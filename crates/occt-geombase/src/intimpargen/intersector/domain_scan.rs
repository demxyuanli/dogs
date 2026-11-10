use crate::intcurve::iconic_tool::IntCurveIConicTool;
use crate::intcurve::pconic::IntCurvePConic;
use occt_gp::gp::GpPnt2d;
use occt_gp::intres2d::IntRes2dDomain;

use super::structure::IntImpParGenIntersector;

impl IntImpParGenIntersector {
    /// `And_Domaine_Objet1_Intersections` (`gxx:42-222`).
    #[allow(clippy::too_many_arguments)]
    pub(super) fn and_domaine_objet1_intersections(
        imp_tool: &IntCurveIConicTool,
        imp_curve_domain: &IntRes2dDomain,
        par_curve: &IntCurvePConic,
        par_curve_domain: &IntRes2dDomain,
        nb_resultats: &mut i32,
        inter2_and_domain2: &[f64],
        inter1: &[f64],
        resultat1: &mut Vec<f64>,
        resultat2: &mut Vec<f64>,
        eps_nul: f64,
    ) {
        let nb_bornes_intersection = *nb_resultats;
        *nb_resultats = 0;

        let mut i = 1;
        while i <= nb_bornes_intersection {
            let mut param1 = inter1[i as usize];
            let mut param2 = inter1[(i + 1) as usize];

            let mut indice_1 = i;
            let mut indice_2 = i + 1;

            if param1 > param2 {
                std::mem::swap(&mut param1, &mut param2);
                indice_1 = i + 1;
                indice_2 = i;
            }

            let pt1 = imp_tool.value(param1);
            let pt2 = imp_tool.value(param2);
            let mut pt = GpPnt2d::new(0.0, 0.0);

            let mut is_on_the_imp_curve_domain1 = true;
            let mut is_on_the_imp_curve_domain2 = true;

            if imp_curve_domain.has_first_point() {
                if param1 < imp_curve_domain.first_parameter() {
                    if pt1.distance(imp_curve_domain.first_point())
                        > imp_curve_domain.first_tolerance()
                    {
                        is_on_the_imp_curve_domain1 = false;
                    }
                }
            }
            if is_on_the_imp_curve_domain1 && imp_curve_domain.has_last_point() {
                if param1 > imp_curve_domain.last_parameter() {
                    if pt1.distance(imp_curve_domain.last_point())
                        > imp_curve_domain.last_tolerance()
                    {
                        is_on_the_imp_curve_domain1 = false;
                    }
                }
            }

            if imp_curve_domain.has_first_point() {
                if param2 < imp_curve_domain.first_parameter() {
                    if pt2.distance(imp_curve_domain.first_point())
                        > imp_curve_domain.first_tolerance()
                    {
                        is_on_the_imp_curve_domain2 = false;
                    }
                }
            }
            if is_on_the_imp_curve_domain2 && imp_curve_domain.has_last_point() {
                if param2 > imp_curve_domain.last_parameter() {
                    if pt2.distance(imp_curve_domain.last_point())
                        > imp_curve_domain.last_tolerance()
                    {
                        is_on_the_imp_curve_domain2 = false;
                    }
                }
            }

            if is_on_the_imp_curve_domain1 {
                *nb_resultats += 1;
                resultat1.push(inter1[indice_1 as usize]);
                resultat2.push(inter2_and_domain2[indice_1 as usize]);

                if is_on_the_imp_curve_domain2 {
                    *nb_resultats += 1;
                    resultat1.push(inter1[indice_2 as usize]);
                    resultat2.push(inter2_and_domain2[indice_2 as usize]);
                } else {
                    *nb_resultats += 1;
                    let t = imp_curve_domain.last_parameter();
                    resultat1.push(t);
                    resultat2.push(Self::find_v(
                        t,
                        &mut pt,
                        imp_tool,
                        par_curve,
                        par_curve_domain,
                        inter2_and_domain2[indice_1 as usize],
                        inter2_and_domain2[indice_2 as usize],
                        eps_nul,
                    ));
                }
            } else if is_on_the_imp_curve_domain2 {
                *nb_resultats += 1;
                let t = imp_curve_domain.first_parameter();
                resultat1.push(t);
                resultat2.push(Self::find_v(
                    t,
                    &mut pt,
                    imp_tool,
                    par_curve,
                    par_curve_domain,
                    inter2_and_domain2[indice_1 as usize],
                    inter2_and_domain2[indice_2 as usize],
                    eps_nul,
                ));

                *nb_resultats += 1;
                resultat1.push(inter1[indice_2 as usize]);
                resultat2.push(inter2_and_domain2[indice_2 as usize]);
            } else if param1 < imp_curve_domain.first_parameter()
                && param2 > imp_curve_domain.last_parameter()
            {
                *nb_resultats += 1;
                let t = imp_curve_domain.first_parameter();
                resultat1.push(t);
                resultat2.push(Self::find_v(
                    t,
                    &mut pt,
                    imp_tool,
                    par_curve,
                    par_curve_domain,
                    inter2_and_domain2[indice_1 as usize],
                    inter2_and_domain2[indice_2 as usize],
                    eps_nul,
                ));

                *nb_resultats += 1;
                let t = imp_curve_domain.last_parameter();
                resultat1.push(t);
                resultat2.push(Self::find_v(
                    t,
                    &mut pt,
                    imp_tool,
                    par_curve,
                    par_curve_domain,
                    inter2_and_domain2[indice_1 as usize],
                    inter2_and_domain2[indice_2 as usize],
                    eps_nul,
                ));
            }

            i += 2;
        }
    }
}
