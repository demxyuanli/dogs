//! `IntImpParGen_Intersector` (`IntImpParGen_Intersector.gxx:42-822`) as
//! instantiated by `Geom2dInt_TheIntConicCurveOfGInter`
//! (`Geom2dInt_TheIntConicCurveOfGInter_0.cxx:32-52`): `ImpTool =
//! IntCurve_IConicTool`, `ParCurve = Adaptor2d_Curve2d` (port `Curve2d`),
//! `ParTool = Geom2dInt_Geom2dCurveTool`, `ProjectOnPCurveTool =
//! Geom2dInt_TheProjPCurOfGInter`. Monomorphised onto exactly that
//! instantiation, the way the C++ template expands. The `occt-core`
//! `intimpargen::intersector` is the same template expanded for
//! `IntCurve_IntImpConicParConic`; this is the `Geom2dInt` expansion.
use occt_core::gp::GpPnt2d;
use occt_core::intcurve::iconic_tool::IntCurveIConicTool;
use crate::curve::Curve2d;
use super::{curve_tool, proj_p_cur};
use occt_core::intres2d::{
    IntRes2dDomain, IntRes2dIntersection, IntRes2dIntersectionPoint,
    IntRes2dIntersectionSegment, IntRes2dPosition,
};
use occt_core::math_fn::{MathFunction, MathFunctionWithDerivative};
use occt_core::math_function_all_roots::FunctionAllRoots;
use occt_core::math_function_sample::FunctionSample;

use occt_core::intimpargen::gen;

/// `IntCurve_MyImpParToolOfIntImpConicParConic`
/// (`IntCurve_MyImpParToolOfIntImpConicParConic.hxx:31-57`, `_0.cxx:26-59`):
/// the `math_FunctionWithDerivative` handed to `math_FunctionAllRoots`. It
/// evaluates the signed distance from a point of the parametric curve to the
/// implicit curve.
pub struct MyImpParTool<'a> {
    the_imp_tool: &'a IntCurveIConicTool,
    the_par_curve: &'a dyn Curve2d,
}

impl<'a> MyImpParTool<'a> {
    /// `IntCurve_MyImpParToolOfIntImpConicParConic(ITool, PC)` (`_0.cxx:26-31`).
    pub fn new(the_imp_tool: &'a IntCurveIConicTool, the_par_curve: &'a dyn Curve2d) -> Self {
        Self {
            the_imp_tool,
            the_par_curve,
        }
    }
}

impl MathFunction for MyImpParTool<'_> {
    /// `Value(Param, ApproxDistance)` (`_0.cxx:33-39`).
    fn value(&mut self, param: f64, approx_distance: &mut f64) -> bool {
        *approx_distance = self
            .the_imp_tool
            .distance(&curve_tool::value(self.the_par_curve, param));
        true
    }
}

impl MathFunctionWithDerivative for MyImpParTool<'_> {
    /// `Derivative(Param, D)` (`_0.cxx:41-50`).
    fn derivative(&mut self, param: f64, d: &mut f64) -> bool {
        let grad = self
            .the_imp_tool
            .grad_distance(&curve_tool::value(self.the_par_curve, param));
        let (_pt, tan) = curve_tool::d1(self.the_par_curve, param);
        *d = grad.dot(&tan);
        true
    }

    /// `Values(Param, F, D)` (`_0.cxx:52-59`).
    fn values(&mut self, param: f64, approx_distance: &mut f64, deriv: &mut f64) -> bool {
        self.value(param, approx_distance);
        self.derivative(param, deriv);
        true
    }
}

/// `IntImpParGen_Intersector` (`IntImpParGen_Intersector.gxx`). The OCCT base
/// `IntRes2d_Intersection` is held in [`Self::base`].
pub struct Geom2dIntConicCurve {
    /// The `IntRes2d_Intersection` base.
    pub base: IntRes2dIntersection,
}

impl Default for Geom2dIntConicCurve {
    fn default() -> Self {
        Self::new()
    }
}

impl Geom2dIntConicCurve {
    /// `IntImpParGen_Intersector()` (`gxx:224-228`).
    pub fn new() -> Self {
        let mut base = IntRes2dIntersection::new();
        base.done = false;
        Self { base }
    }

    /// `IntImpParGen_Intersector(TheImpTool, TheImpCurveDomain, TheParCurve,
    /// TheParCurveDomain, TolConf, Tol)` (`gxx:233-242`).
    #[allow(clippy::too_many_arguments)]
    pub fn with_params(
        the_imp_tool: &IntCurveIConicTool,
        the_imp_curve_domain: &IntRes2dDomain,
        the_par_curve: &dyn Curve2d,
        the_par_curve_domain: &IntRes2dDomain,
        tol_conf: f64,
        tol: f64,
    ) -> Self {
        let mut s = Self::new();
        s.perform(
            the_imp_tool,
            the_imp_curve_domain,
            the_par_curve,
            the_par_curve_domain,
            tol_conf,
            tol,
        );
        s
    }

    /// `And_Domaine_Objet1_Intersections` (`gxx:42-222`).
    #[allow(clippy::too_many_arguments)]
    fn and_domaine_objet1_intersections(
        the_imp_tool: &IntCurveIConicTool,
        the_imp_curve_domain: &IntRes2dDomain,
        the_par_curve: &dyn Curve2d,
        the_par_curve_domain: &IntRes2dDomain,
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

            let pt1 = the_imp_tool.value(param1);
            let pt2 = the_imp_tool.value(param2);
            let mut pt = GpPnt2d::new(0.0, 0.0);

            let mut is_on_the_imp_curve_domain1 = true;
            let mut is_on_the_imp_curve_domain2 = true;

            if the_imp_curve_domain.has_first_point() {
                if param1 < the_imp_curve_domain.first_parameter() {
                    if pt1.distance(the_imp_curve_domain.first_point())
                        > the_imp_curve_domain.first_tolerance()
                    {
                        is_on_the_imp_curve_domain1 = false;
                    }
                }
            }
            if is_on_the_imp_curve_domain1 && the_imp_curve_domain.has_last_point() {
                if param1 > the_imp_curve_domain.last_parameter() {
                    if pt1.distance(the_imp_curve_domain.last_point())
                        > the_imp_curve_domain.last_tolerance()
                    {
                        is_on_the_imp_curve_domain1 = false;
                    }
                }
            }

            if the_imp_curve_domain.has_first_point() {
                if param2 < the_imp_curve_domain.first_parameter() {
                    if pt2.distance(the_imp_curve_domain.first_point())
                        > the_imp_curve_domain.first_tolerance()
                    {
                        is_on_the_imp_curve_domain2 = false;
                    }
                }
            }
            if is_on_the_imp_curve_domain2 && the_imp_curve_domain.has_last_point() {
                if param2 > the_imp_curve_domain.last_parameter() {
                    if pt2.distance(the_imp_curve_domain.last_point())
                        > the_imp_curve_domain.last_tolerance()
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
                    let t = the_imp_curve_domain.last_parameter();
                    resultat1.push(t);
                    resultat2.push(Self::find_v(
                        t,
                        &mut pt,
                        the_imp_tool,
                        the_par_curve,
                        the_par_curve_domain,
                        inter2_and_domain2[indice_1 as usize],
                        inter2_and_domain2[indice_2 as usize],
                        eps_nul,
                    ));
                }
            } else if is_on_the_imp_curve_domain2 {
                *nb_resultats += 1;
                let t = the_imp_curve_domain.first_parameter();
                resultat1.push(t);
                resultat2.push(Self::find_v(
                    t,
                    &mut pt,
                    the_imp_tool,
                    the_par_curve,
                    the_par_curve_domain,
                    inter2_and_domain2[indice_1 as usize],
                    inter2_and_domain2[indice_2 as usize],
                    eps_nul,
                ));

                *nb_resultats += 1;
                resultat1.push(inter1[indice_2 as usize]);
                resultat2.push(inter2_and_domain2[indice_2 as usize]);
            } else if param1 < the_imp_curve_domain.first_parameter()
                && param2 > the_imp_curve_domain.last_parameter()
            {
                *nb_resultats += 1;
                let t = the_imp_curve_domain.first_parameter();
                resultat1.push(t);
                resultat2.push(Self::find_v(
                    t,
                    &mut pt,
                    the_imp_tool,
                    the_par_curve,
                    the_par_curve_domain,
                    inter2_and_domain2[indice_1 as usize],
                    inter2_and_domain2[indice_2 as usize],
                    eps_nul,
                ));

                *nb_resultats += 1;
                let t = the_imp_curve_domain.last_parameter();
                resultat1.push(t);
                resultat2.push(Self::find_v(
                    t,
                    &mut pt,
                    the_imp_tool,
                    the_par_curve,
                    the_par_curve_domain,
                    inter2_and_domain2[indice_1 as usize],
                    inter2_and_domain2[indice_2 as usize],
                    eps_nul,
                ));
            }

            i += 2;
        }
    }

    /// `FindU` (`gxx:781-788`).
    fn find_u(
        parameter: f64,
        point: &mut GpPnt2d,
        the_par_curve: &dyn Curve2d,
        the_imp_tool: &IntCurveIConicTool,
    ) -> f64 {
        *point = curve_tool::value(the_par_curve, parameter);
        the_imp_tool.find_parameter(point)
    }

    /// `FindV` (`gxx:790-822`).
    #[allow(clippy::too_many_arguments)]
    fn find_v(
        parameter: f64,
        point: &mut GpPnt2d,
        the_imp_tool: &IntCurveIConicTool,
        the_par_curve: &dyn Curve2d,
        the_par_curve_domain: &IntRes2dDomain,
        v0: f64,
        v1: f64,
        tolerance: f64,
    ) -> f64 {
        *point = the_imp_tool.value(parameter);
        if the_par_curve_domain.is_closed() {
            let v = proj_p_cur::find_parameter(the_par_curve, point, tolerance);
            gen::normalize_on_domain(v, the_par_curve_domain)
        } else {
            let mut vv0 = v0;
            let mut vv1 = v1;
            if v1 < v0 {
                vv0 = v1;
                vv1 = v0;
            }
            let mut x = proj_p_cur::find_parameter_range(
                the_par_curve,
                point,
                vv0,
                vv1,
                tolerance,
            );
            if x > vv1 {
                x = vv1;
            } else if x < vv0 {
                x = vv0;
            }
            x
        }
    }

    /// `Perform` (`gxx:245-779`).
    pub fn perform(
        &mut self,
        the_imp_tool: &IntCurveIConicTool,
        the_imp_curve_domain: &IntRes2dDomain,
        the_par_curve: &dyn Curve2d,
        the_par_curve_domain: &IntRes2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        let mut head_on_imp = false;
        let mut head_on_par = false;
        let mut end_on_imp = false;
        let mut end_on_par = false;

        self.base.reset_fields();

        let mut the_imp_par_tool = MyImpParTool::new(the_imp_tool, the_par_curve);

        if !(the_par_curve_domain.has_first_point() && the_par_curve_domain.has_last_point()) {
            panic!("Standard_ConstructionError: Domaine sur courbe incorrect");
        }

        let nb_echantillons = curve_tool::nb_samples_curve_range(
            the_par_curve,
            the_par_curve_domain.first_parameter(),
            the_par_curve_domain.last_parameter(),
        ) as i32;

        let mut eps_x = curve_tool::eps_x();
        if eps_x > 1.0e-10 {
            eps_x = 1.0e-10;
        }
        let eps_nul = if tol_conf <= 1.0e-10 { 1.0e-10 } else { tol_conf };
        let eps_dist = if tol <= 1.0e-10 { 1.0e-10 } else { tol };

        let tolerance_angulaire = eps_dist;

        if (the_par_curve_domain.last_parameter() - the_par_curve_domain.first_parameter())
            < 100.0 * eps_x
        {
            eps_x = (the_par_curve_domain.last_parameter()
                - the_par_curve_domain.first_parameter())
                * 0.01;
        }

        let sample2 = FunctionSample::new(
            the_par_curve_domain.first_parameter(),
            the_par_curve_domain.last_parameter(),
            nb_echantillons,
        );

        let sol = FunctionAllRoots::new(&mut the_imp_par_tool, &sample2, eps_x, eps_dist, eps_nul);

        if !sol.is_done() {
            self.base.done = false;
            return;
        }

        let nb_segments_solution = sol.nb_intervals();
        let nb_points_solution = sol.nb_points();

        // Traitement des Points Solutions
        for i in 1..=nb_points_solution {
            let param2 = sol.get_point(i);
            let mut pt = GpPnt2d::new(0.0, 0.0);
            let mut param1 = Self::find_u(param2, &mut pt, the_par_curve, the_imp_tool);

            if the_imp_curve_domain.is_closed() {
                param1 = gen::normalize_on_domain(param1, the_imp_curve_domain);
            }

            let mut is_on_the_imp_curve_domain = true;
            if the_imp_curve_domain.has_first_point() {
                if param1 < the_imp_curve_domain.first_parameter() {
                    if pt.distance(the_imp_curve_domain.first_point())
                        > the_imp_curve_domain.first_tolerance()
                    {
                        is_on_the_imp_curve_domain = false;
                    }
                }
            }
            if is_on_the_imp_curve_domain && the_imp_curve_domain.has_last_point() {
                if param1 > the_imp_curve_domain.last_parameter() {
                    if pt.distance(the_imp_curve_domain.last_point())
                        > the_imp_curve_domain.last_tolerance()
                    {
                        is_on_the_imp_curve_domain = false;
                    }
                }
            }

            if is_on_the_imp_curve_domain {
                let (pt1, mut tan1, norm1) = the_imp_tool.d2(param1);
                let (pt2, mut tan2, norm2) = curve_tool::d2(the_par_curve, param2);

                let pos1 = gen::determine_position(the_imp_curve_domain, &pt1, param1);
                let pos2 = gen::determine_position(the_par_curve_domain, &pt2, param2);

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

                let mut trans1 = occt_core::intres2d::IntRes2dTransition::new();
                let mut trans2 = occt_core::intres2d::IntRes2dTransition::new();
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
            let (param2_inf, param2_sup) = sol.get_interval(j);
            let mut param1_inf = Self::find_u(param2_inf, &mut ptemp, the_par_curve, the_imp_tool);
            let mut param1_sup = Self::find_u(param2_sup, &mut ptemp, the_par_curve, the_imp_tool);

            if the_imp_curve_domain.is_closed() {
                let (param1_origine, param1_fin) = the_imp_curve_domain.equivalent_parameters();
                let periode = param1_fin - param1_origine;

                while param1_inf < param1_origine {
                    param1_inf += periode;
                }
                while param1_sup < param1_origine {
                    param1_sup += periode;
                }

                let (_p2, mut t2, n2) = curve_tool::d2(the_par_curve, param2_inf);
                let (_p1, mut t1, n1) = the_imp_tool.d2(param1_inf);
                if t1.magnitude() <= occt_core::precision::RESOLUTION {
                    t1 = n1;
                }
                if t2.magnitude() <= occt_core::precision::RESOLUTION {
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

                if the_imp_curve_domain.last_parameter()
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

                if the_imp_curve_domain.first_parameter()
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
            the_imp_tool,
            the_imp_curve_domain,
            the_par_curve,
            the_par_curve_domain,
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

                let pt1_on1 = the_imp_tool.value(param1_on1);
                let pt2_on1 = the_imp_tool.value(param2_on1);
                let pt1_on2 = curve_tool::value(the_par_curve, param1_on2);
                let pt2_on2 = curve_tool::value(the_par_curve, param2_on2);

                if !the_imp_curve_domain.is_closed() {
                    if pt1_on1.distance(&pt2_on1) <= dist_mini_imp_curve
                        && pt1_on2.distance(&pt2_on2) <= dist_mini_imp_curve
                    {
                        only_one_point = true;
                    }
                }

                param1_on1 = gen::normalize_on_domain(param1_on1, the_imp_curve_domain);
                param1_on2 = gen::normalize_on_domain(param1_on2, the_par_curve_domain);

                let (mut pt1_on1, mut tan1, norm1) = the_imp_tool.d2(param1_on1);
                let (_pt1_on2, mut tan2, norm2) = curve_tool::d2(the_par_curve, param1_on2);

                let pos1 = gen::determine_position(the_imp_curve_domain, &pt1_on1, param1_on1);
                let pos2 = gen::determine_position(the_par_curve_domain, &pt1_on2, param1_on2);

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

                let mut trans1 = occt_core::intres2d::IntRes2dTransition::new();
                let mut trans2 = occt_core::intres2d::IntRes2dTransition::new();
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

                    let param2_on1 = gen::normalize_on_domain(param2_on1, the_imp_curve_domain);
                    let param2_on2 = gen::normalize_on_domain(param2_on2, the_par_curve_domain);

                    let (mut pt2_on1, mut tan1, norm1) = the_imp_tool.d2(param2_on1);
                    let (pt2_on2, mut tan2, norm2) = curve_tool::d2(the_par_curve, param2_on2);

                    let pos1 = gen::determine_position(the_imp_curve_domain, &pt2_on1, param2_on1);
                    let pos2 = gen::determine_position(the_par_curve_domain, &pt2_on2, param2_on2);

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

                    let mut trans1 = occt_core::intres2d::IntRes2dTransition::new();
                    let mut trans2 = occt_core::intres2d::IntRes2dTransition::new();
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
        if !head_on_imp && the_imp_curve_domain.has_first_point() {
            if !head_on_par
                && the_imp_curve_domain
                    .first_point()
                    .distance(the_par_curve_domain.first_point())
                    <= the_imp_curve_domain
                        .first_tolerance()
                        .max(the_par_curve_domain.first_tolerance())
            {
                let param1 = the_imp_curve_domain.first_parameter();
                let param2 = the_par_curve_domain.first_parameter();
                let (_pt1, mut tan1, norm1) = the_imp_tool.d2(param1);
                let (_pt2, mut tan2, norm2) = curve_tool::d2(the_par_curve, param2);
                let mut trans1 = occt_core::intres2d::IntRes2dTransition::new();
                let mut trans2 = occt_core::intres2d::IntRes2dTransition::new();
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
                    the_imp_curve_domain.first_point(),
                    param1,
                    param2,
                    &trans1,
                    &trans2,
                    self.base.reversed_parameters(),
                );
                self.base.insert(&ip);
            }
            if !end_on_par
                && the_imp_curve_domain
                    .first_point()
                    .distance(the_par_curve_domain.last_point())
                    <= the_imp_curve_domain
                        .first_tolerance()
                        .max(the_par_curve_domain.last_tolerance())
            {
                let param1 = the_imp_curve_domain.first_parameter();
                let param2 = the_par_curve_domain.last_parameter();
                let (_pt1, mut tan1, norm1) = the_imp_tool.d2(param1);
                let (_pt2, mut tan2, norm2) = curve_tool::d2(the_par_curve, param2);
                let mut trans1 = occt_core::intres2d::IntRes2dTransition::new();
                let mut trans2 = occt_core::intres2d::IntRes2dTransition::new();
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
                    the_imp_curve_domain.first_point(),
                    param1,
                    param2,
                    &trans1,
                    &trans2,
                    self.base.reversed_parameters(),
                );
                self.base.insert(&ip);
            }
        }

        if !end_on_imp && the_imp_curve_domain.has_last_point() {
            if !head_on_par
                && the_imp_curve_domain
                    .last_point()
                    .distance(the_par_curve_domain.first_point())
                    <= the_imp_curve_domain
                        .last_tolerance()
                        .max(the_par_curve_domain.first_tolerance())
            {
                let param1 = the_imp_curve_domain.last_parameter();
                let param2 = the_par_curve_domain.first_parameter();
                let (_pt1, mut tan1, norm1) = the_imp_tool.d2(param1);
                let (_pt2, mut tan2, norm2) = curve_tool::d2(the_par_curve, param2);
                let mut trans1 = occt_core::intres2d::IntRes2dTransition::new();
                let mut trans2 = occt_core::intres2d::IntRes2dTransition::new();
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
                    the_imp_curve_domain.last_point(),
                    param1,
                    param2,
                    &trans1,
                    &trans2,
                    self.base.reversed_parameters(),
                );
                self.base.insert(&ip);
            }
            if !end_on_par
                && the_imp_curve_domain
                    .last_point()
                    .distance(the_par_curve_domain.last_point())
                    <= the_imp_curve_domain
                        .last_tolerance()
                        .max(the_par_curve_domain.last_tolerance())
            {
                let param1 = the_imp_curve_domain.last_parameter();
                let param2 = the_par_curve_domain.last_parameter();
                let (_pt1, mut tan1, norm1) = the_imp_tool.d2(param1);
                let (_pt2, mut tan2, norm2) = curve_tool::d2(the_par_curve, param2);
                let mut trans1 = occt_core::intres2d::IntRes2dTransition::new();
                let mut trans2 = occt_core::intres2d::IntRes2dTransition::new();
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
                    the_imp_curve_domain.last_point(),
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