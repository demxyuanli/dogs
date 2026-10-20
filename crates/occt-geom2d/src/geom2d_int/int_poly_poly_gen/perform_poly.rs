//! `Perform` with polygon-approximation refinement (`gxx:390-870`).

use super::prelude::*;
use super::super::curve_tool;
use super::super::proj_p_cur;
use super::*;

impl Geom2dIntIntPolyPolyGen {
    /// `Perform(C1, D1, C2, D2, TolConf, Tol, NbIter, DeltaU, DeltaV)`
    /// (`gxx:1040-1146`), the polygon engine of the two-curve case.
    pub(super) fn perform_poly(
        &mut self,
        c1: &dyn Curve2d,
        d1: &IntRes2dDomain,
        c2: &dyn Curve2d,
        d2: &IntRes2dDomain,
        tol_conf: f64,
        tol: f64,
        nb_iter: i32,
        delta_u: f64,
        delta_v: f64,
    ) {
        self.base.done = false;

        if nb_iter > NBITER_MAX_POLYGON {
            return;
        }

        // Number of samples running.
        let mut nbsamples_on_c1 =
            curve_tool::nb_samples_curve_range(c1, d1.first_parameter(), d1.last_parameter());
        let mut nbsamples_on_c2 =
            curve_tool::nb_samples_curve_range(c2, d2.first_parameter(), d2.last_parameter());

        if nb_iter == 0 {
            // Minimal number of points.
            nbsamples_on_c1 = nbsamples_on_c1.max(self.min_pnt_nb);
            nbsamples_on_c2 = nbsamples_on_c2.max(self.min_pnt_nb);
        } else {
            // Increase number of samples in second and next iterations.
            nbsamples_on_c1 = (5 * (nbsamples_on_c1 * nb_iter as usize)) / 4;
            nbsamples_on_c2 = (5 * (nbsamples_on_c2 * nb_iter as usize)) / 4;
        }

        let mut poly1 = Geom2dIntPolygon2d::new(c1, nbsamples_on_c1, d1, tol);
        let mut poly2 = Geom2dIntPolygon2d::new(c2, nbsamples_on_c2, d2, tol);

        if (poly1.deflection_over_estimation() > tol_conf)
            && (poly2.deflection_over_estimation() > tol_conf)
        {
            let deflection_sum = poly1.deflection_over_estimation().max(tol_conf)
                + poly2.deflection_over_estimation().max(tol_conf);

            if nbsamples_on_c2 > nbsamples_on_c1 {
                poly2.compute_with_box(c2, poly1.bounding());
                poly1.set_deflection_over_estimation(deflection_sum);
                poly1.compute_with_box(c1, poly2.bounding());
            } else {
                poly1.compute_with_box(c1, poly2.bounding());
                poly2.set_deflection_over_estimation(deflection_sum);
                poly2.compute_with_box(c2, poly1.bounding());
            }
        }

        //----------------------------------------------------------------------
        //-- if the deflection less then the Tolerance of Confusion
        //-- Then the deflection of the polygon is set in TolConf
        //-- (Detection of Tangency Zones)
        //----------------------------------------------------------------------
        if poly1.deflection_over_estimation() < tol_conf {
            poly1.set_deflection_over_estimation(tol_conf);
        }
        if poly2.deflection_over_estimation() < tol_conf {
            poly2.set_deflection_over_estimation(tol_conf);
        }

        // for case when a few polygon points were replaced by line
        // if exact solution was not found
        // then search of precise solution will be repeated
        // for polygon contains all initial points
        // secondary search will be performed only for case when initial points
        // were dropped
        let is_full_representation = (poly1.nb_segments() == nbsamples_on_c1
            && poly2.nb_segments() == nbsamples_on_c2);

        if !self.find_intersect(
            c1,
            d1,
            c2,
            d2,
            tol_conf,
            tol,
            nb_iter,
            delta_u,
            delta_v,
            &poly1,
            &poly2,
            is_full_representation,
        ) && !is_full_representation
        {
            if poly1.nb_segments() < nbsamples_on_c1 {
                poly1 = Geom2dIntPolygon2d::new(c1, nbsamples_on_c1, d1, tol);
            }
            if poly2.nb_segments() < nbsamples_on_c2 {
                poly2 = Geom2dIntPolygon2d::new(c2, nbsamples_on_c2, d2, tol);
            }

            self.find_intersect(
                c1, d1, c2, d2, tol_conf, tol, nb_iter, delta_u, delta_v, &poly1, &poly2, true,
            );
        }

        self.base.done = true;
    }

    /// `Perform(C1, D1, TolConf, Tol, NbIter, DeltaU, DeltaV)` (`gxx:390-870`),
    /// the polygon engine of the auto-intersection case.
    pub(super) fn perform_self_poly(
        &mut self,
        c1: &dyn Curve2d,
        d1: &IntRes2dDomain,
        tol_conf: f64,
        tol: f64,
        nb_iter: i32,
        delta_u: f64,
        delta_v: f64,
    ) {
        let _ = (delta_u, delta_v);
        self.base.done = false;

        let mut nbsamples =
            curve_tool::nb_samples_curve_range(c1, d1.first_parameter(), d1.last_parameter());

        if nb_iter > 3 || (nb_iter > 2 && nbsamples > 100) {
            return;
        }

        nbsamples *= 2; //---  We take systematically two times more points
                        //--   than on a normal curve.
                        //--   Auto-intersecting curves often produce
                        //--   polygons rather far from the curve with parameter ct.

        if nb_iter > 0 {
            nbsamples = (3 * (nbsamples * nb_iter as usize)) / 2;
        }
        let poly1 = Geom2dIntPolygon2d::new(c1, nbsamples, d1, tol);
        if !poly1.auto_intersection_is_possible() {
            self.base.done = true;
            return;
        }
        //----------------------------------------------------------------------
        //-- If the deflection is less than the Tolerance of Confusion
        //-- then the deflection of the polygon is set in TolConf
        //-- (Detection of Tangency Zones)
        //----------------------------------------------------------------------
        let mut poly1 = poly1;
        if poly1.deflection_over_estimation() < tol_conf {
            poly1.set_deflection_over_estimation(tol_conf);
        }

        let inter_pp = IntfInterferencePolygon2d::self_intersection(&poly1);
        let mut eip = Geom2dIntExactIntersectionPoint::new(c1, c1, tol_conf);

        //----------------------------------------------------------------------
        //-- Processing of SectionPoint
        //----------------------------------------------------------------------
        let nbsp = inter_pp.base.nb_section_points();
        if nbsp >= 1 {
            //-- filtering, filtering, filtering ...
            let mut tri_index: Vec<i32> = vec![0; nbsp + 1];
            let mut ptr_seg_index1: Vec<i32> = vec![0; nbsp + 1];
            let mut ptr_seg_index2: Vec<i32> = vec![0; nbsp + 1];
            for i in 1..=nbsp {
                tri_index[i] = i as i32;
                let spnt1 = inter_pp.base.pnt_value(i - 1);
                let (_, si1, _p1) = spnt1.info_first();
                ptr_seg_index1[i] = si1;
                let (_, si2, _p2) = spnt1.info_second();
                ptr_seg_index2[i] = si2;
            }

            loop {
                let mut triok = true;
                for tr in 1..nbsp {
                    let seg_index1 = ptr_seg_index1[tri_index[tr] as usize];
                    let seg_index_1 = ptr_seg_index1[tri_index[tr + 1] as usize];
                    let seg_index2 = ptr_seg_index2[tri_index[tr] as usize];
                    let seg_index_2 = ptr_seg_index2[tri_index[tr + 1] as usize];

                    if seg_index1 > seg_index_1 {
                        tri_index.swap(tr, tr + 1);
                        triok = false;
                    } else if seg_index1 == seg_index_1 && seg_index2 > seg_index_2 {
                        tri_index.swap(tr, tr + 1);
                        triok = false;
                    }
                }
                if triok {
                    break;
                }
            }

            //-- supression des doublons Si Si !
            for i in 1..nbsp {
                let a = (ptr_seg_index1[tri_index[i] as usize] == ptr_seg_index1[tri_index[i + 1] as usize])
                    && (ptr_seg_index2[tri_index[i] as usize]
                        == ptr_seg_index2[tri_index[i + 1] as usize]);
                if a {
                    tri_index[i] = -(i as i32);
                }
            }

            for sp in 1..=nbsp {
                if tri_index[sp] > 0 {
                    // `InterPP.PntValue(TriIndex[sp])` (`gxx:513`): OCCT's
                    // `PntValue` is 1-based; this port's is 0-based.
                    let spnt = inter_pp.base.pnt_value(tri_index[sp] as usize - 1);
                    let (_, mut seg_index1, mut param_on1) = spnt.info_first();
                    let (_, mut seg_index2, mut param_on2) = spnt.info_second();

                    if (seg_index1 - seg_index2).abs() > 1 {
                        eip.perform(
                            &poly1,
                            &poly1,
                            &mut seg_index1,
                            &mut seg_index2,
                            &mut param_on1,
                            &mut param_on2,
                        );
                        if eip.nb_roots() >= 1 {
                            //--------------------------------------------------------------------
                            //-- It is checked if the found point is a root
                            //--------------------------------------------------------------------
                            let (u, v) = eip.roots();

                            let (mut p1, tan1) = curve_tool::geom2d_d1(c1, u);
                            let (mut p2, tan2) = curve_tool::geom2d_d1(c1, v);
                            let mut dist = p1.distance(&p2);
                            let eps_x1 = 10.0 * curve_tool::geom2d_eps_x();

                            if (u - v).abs() <= eps_x1 {
                                //-----------------------------------------
                                //-- Solution not valid
                                //-----------------------------------------
                                dist = tol_conf + 1.0;
                            }

                            //-----------------------------------------------------------------
                            //-- It is checked if the point (u,v) already exists
                            //--
                            self.base.done = true;
                            let nbp = self.base.nb_points();
                            let mut p = 1usize;
                            while p <= nbp {
                                let ppt = self.base.point(p);
                                if (u - ppt.param_on_first()).abs() <= eps_x1
                                    && (v - ppt.param_on_second()).abs() <= eps_x1
                                {
                                    dist = tol_conf + 1.0;
                                    break;
                                }
                                p += 1;
                            }

                            if dist <= tol_conf {
                                //-- Or the point is already present
                                let mut pos1 = IntRes2dPosition::Middle;
                                let mut pos2 = IntRes2dPosition::Middle;
                                let mut trans1 = IntRes2dTransition::new();
                                let mut trans2 = IntRes2dTransition::new();
                                //-----------------------------------------------------------------
                                //-- Calculate Positions of Points on the curve
                                //-----------------------------------------------------------------
                                if p1.distance(self.domain_on_curve1.first_point())
                                    <= self.domain_on_curve1.first_tolerance()
                                {
                                    pos1 = IntRes2dPosition::Head;
                                } else if p1.distance(self.domain_on_curve1.last_point())
                                    <= self.domain_on_curve1.last_tolerance()
                                {
                                    pos1 = IntRes2dPosition::End;
                                }

                                if p2.distance(self.domain_on_curve2.first_point())
                                    <= self.domain_on_curve2.first_tolerance()
                                {
                                    pos2 = IntRes2dPosition::Head;
                                } else if p2.distance(self.domain_on_curve2.last_point())
                                    <= self.domain_on_curve2.last_tolerance()
                                {
                                    pos2 = IntRes2dPosition::End;
                                }
                                //-----------------------------------------------------------------
                                let mut tan1 = tan1;
                                let mut tan2 = tan2;
                                if !determine_transition_simple(
                                    pos1, &tan1, &mut trans1, pos2, &tan2, &mut trans2, tol_conf,
                                ) {
                                    let (pp1, t1, n1) = curve_tool::geom2d_d2(c1, u);
                                    let (pp2, t2, n2) = curve_tool::geom2d_d2(c1, v);
                                    p1 = pp1;
                                    p2 = pp2;
                                    tan1 = t1;
                                    tan2 = t2;
                                    determine_transition_touch(
                                        pos1, &mut tan1, &n1, &mut trans1, pos2, &mut tan2, &n2,
                                        &mut trans2, tol_conf,
                                    );
                                }
                                let ip = IntRes2dIntersectionPoint::with_transitions(
                                    &p1, u, v, &trans1, &trans2, false,
                                );
                                let _ = (tan1, tan2);
                                self.base.insert(&ip);
                            }
                        }
                    }
                }
            }
        }

        //----------------------------------------------------------------------
        //-- Processing of TangentZone
        //----------------------------------------------------------------------
        let nbtz = inter_pp.base.nb_tangent_zones();
        for tz in 1..=nbtz {
            let zone = inter_pp.base.zone_value(tz - 1);
            let nb_pnts = zone.number_of_points();
            //====================================================================
            //== Find the first and the last point in the tangency zone.
            //====================================================================
            let mut param_sup_on_curve2;
            let mut param_inf_on_curve2;
            let mut param_sup_on_curve1;
            let mut param_inf_on_curve1;
            let mut poly_u_inf = f64::MAX;
            let mut poly_v_inf = f64::MAX;
            let mut poly_u_sup = -f64::MAX;
            let mut poly_v_sup = -f64::MAX;
            param_sup_on_curve2 = -f64::MAX;
            param_sup_on_curve1 = -f64::MAX;
            param_inf_on_curve2 = f64::MAX;
            param_inf_on_curve1 = f64::MAX;
            for qq in 1..=nb_pnts {
                let spnt1 = zone.point(qq);
                let (_, mut seg_index1on_p1, mut param_on_line) = spnt1.info_first();
                if seg_index1on_p1 > poly1.nb_segments() as i32 {
                    seg_index1on_p1 -= 1;
                    param_on_line = 1.0;
                }
                if seg_index1on_p1 <= 0 {
                    seg_index1on_p1 = 1;
                    param_on_line = 0.0;
                }
                let poly_u_inf_seg = poly1.approx_param_on_curve(seg_index1on_p1 as usize, param_on_line);

                let (_, mut seg_index1on_p2, mut param_on_line) = spnt1.info_second();
                if seg_index1on_p2 > poly1.nb_segments() as i32 {
                    seg_index1on_p2 -= 1;
                    param_on_line = 1.0;
                }
                if seg_index1on_p2 <= 0 {
                    seg_index1on_p2 = 1;
                    param_on_line = 0.0;
                }
                let poly_v_inf_seg = poly1.approx_param_on_curve(seg_index1on_p2 as usize, param_on_line);

                if param_inf_on_curve1 > poly_u_inf_seg {
                    param_inf_on_curve1 = poly_u_inf_seg;
                }
                if param_inf_on_curve2 > poly_v_inf_seg {
                    param_inf_on_curve2 = poly_v_inf_seg;
                }
                if param_sup_on_curve1 < poly_u_inf_seg {
                    param_sup_on_curve1 = poly_u_inf_seg;
                }
                if param_sup_on_curve2 < poly_v_inf_seg {
                    param_sup_on_curve2 = poly_v_inf_seg;
                }
            }

            poly_u_inf = param_inf_on_curve1;
            poly_u_sup = param_sup_on_curve1;
            poly_v_inf = param_inf_on_curve2;
            poly_v_sup = param_sup_on_curve2;

            let p1a = curve_tool::geom2d_d0(c1, poly_u_inf);
            let p2a = curve_tool::geom2d_d0(c1, poly_v_inf);
            let distmemesens = p1a.square_distance(&p2a);
            let p2b = curve_tool::geom2d_d0(c1, poly_v_sup);
            let distdiffsens = p1a.square_distance(&p2b);
            if distmemesens > distdiffsens {
                std::mem::swap(&mut poly_v_inf, &mut poly_v_sup);
            }

            //-----------------------------------------------------------------
            //-- Calculate Positions of Points on the curve and
            //-- Transitions on each limit of the segment
            //-----------------------------------------------------------------
            let mut pos1 = IntRes2dPosition::Middle;
            let mut pos2 = IntRes2dPosition::Middle;
            let mut trans1 = IntRes2dTransition::new();
            let mut trans2 = IntRes2dTransition::new();

            let (mut p1, mut tan1) = curve_tool::geom2d_d1(c1, poly_u_inf);
            let (mut p2, mut tan2) = curve_tool::geom2d_d1(c1, poly_v_inf);

            if p1.distance(self.domain_on_curve1.first_point())
                <= self.domain_on_curve1.first_tolerance()
            {
                pos1 = IntRes2dPosition::Head;
            } else if p1.distance(self.domain_on_curve1.last_point())
                <= self.domain_on_curve1.last_tolerance()
            {
                pos1 = IntRes2dPosition::End;
            }
            if p2.distance(self.domain_on_curve2.first_point())
                <= self.domain_on_curve2.first_tolerance()
            {
                pos2 = IntRes2dPosition::Head;
            } else if p2.distance(self.domain_on_curve2.last_point())
                <= self.domain_on_curve2.last_tolerance()
            {
                pos2 = IntRes2dPosition::End;
            }

            if pos1 == IntRes2dPosition::Middle && pos2 != IntRes2dPosition::Middle {
                poly_u_inf = proj_p_cur::find_parameter_range(
                    c1,
                    &p2,
                    d1.first_parameter(),
                    d1.last_parameter(),
                    curve_tool::geom2d_eps_x(),
                );
            } else if pos1 != IntRes2dPosition::Middle && pos2 == IntRes2dPosition::Middle {
                poly_v_inf = proj_p_cur::find_parameter_range(
                    c1,
                    &p1,
                    d1.first_parameter(),
                    d1.last_parameter(),
                    curve_tool::geom2d_eps_x(),
                );
            } else if (param_inf_on_curve1 - param_sup_on_curve1).abs()
                > (param_inf_on_curve2 - param_sup_on_curve2).abs()
            {
                poly_v_inf = proj_p_cur::find_parameter_range(
                    c1,
                    &p1,
                    d1.first_parameter(),
                    d1.last_parameter(),
                    curve_tool::geom2d_eps_x(),
                );
            } else {
                poly_u_inf = proj_p_cur::find_parameter_range(
                    c1,
                    &p2,
                    d1.first_parameter(),
                    d1.last_parameter(),
                    curve_tool::geom2d_eps_x(),
                );
            }

            if !determine_transition_simple(pos1, &tan1, &mut trans1, pos2, &tan2, &mut trans2, tol_conf)
            {
                let (pp1, t1, n1) = curve_tool::geom2d_d2(c1, poly_u_inf);
                let (pp2, t2, n2) = curve_tool::geom2d_d2(c1, poly_v_inf);
                p1 = pp1;
                p2 = pp2;
                tan1 = t1;
                tan2 = t2;
                determine_transition_touch(
                    pos1, &mut tan1, &n1, &mut trans1, pos2, &mut tan2, &n2, &mut trans2, tol_conf,
                );
            }
            let pt_seg1 =
                IntRes2dIntersectionPoint::with_transitions(&p1, poly_u_inf, poly_v_inf, &trans1, &trans2, false);
            //----------------------------------------------------------------------

            if (poly_u_inf - poly_u_sup).abs() <= curve_tool::geom2d_eps_x()
                || (poly_v_inf - poly_v_sup).abs() <= curve_tool::geom2d_eps_x()
            {
                // bad segment
            } else {
                let (mut p1, mut tan1) = curve_tool::geom2d_d1(c1, poly_u_sup);
                let (mut p2, mut tan2) = curve_tool::geom2d_d1(c1, poly_v_sup);
                pos1 = IntRes2dPosition::Middle;
                pos2 = IntRes2dPosition::Middle;

                if p1.distance(self.domain_on_curve1.first_point())
                    <= self.domain_on_curve1.first_tolerance()
                {
                    pos1 = IntRes2dPosition::Head;
                } else if p1.distance(self.domain_on_curve1.last_point())
                    <= self.domain_on_curve1.last_tolerance()
                {
                    pos1 = IntRes2dPosition::End;
                }
                if p2.distance(self.domain_on_curve2.first_point())
                    <= self.domain_on_curve2.first_tolerance()
                {
                    pos2 = IntRes2dPosition::Head;
                } else if p2.distance(self.domain_on_curve2.last_point())
                    <= self.domain_on_curve2.last_tolerance()
                {
                    pos2 = IntRes2dPosition::End;
                }

                if pos1 == IntRes2dPosition::Middle && pos2 != IntRes2dPosition::Middle {
                    poly_u_sup = proj_p_cur::find_parameter_range(
                        c1,
                        &p2,
                        d1.first_parameter(),
                        d1.last_parameter(),
                        curve_tool::geom2d_eps_x(),
                    );
                } else if pos1 != IntRes2dPosition::Middle && pos2 == IntRes2dPosition::Middle {
                    poly_v_sup = proj_p_cur::find_parameter_range(
                        c1,
                        &p1,
                        d1.first_parameter(),
                        d1.last_parameter(),
                        curve_tool::geom2d_eps_x(),
                    );
                } else if (param_inf_on_curve1 - param_sup_on_curve1).abs()
                    > (param_inf_on_curve2 - param_sup_on_curve2).abs()
                {
                    poly_v_sup = proj_p_cur::find_parameter_range(
                        c1,
                        &p1,
                        d1.first_parameter(),
                        d1.last_parameter(),
                        curve_tool::geom2d_eps_x(),
                    );
                } else {
                    poly_u_sup = proj_p_cur::find_parameter_range(
                        c1,
                        &p2,
                        d1.first_parameter(),
                        d1.last_parameter(),
                        curve_tool::geom2d_eps_x(),
                    );
                }

                if !determine_transition_simple(
                    pos1, &tan1, &mut trans1, pos2, &tan2, &mut trans2, tol_conf,
                ) {
                    let (_, t1, n1) = curve_tool::geom2d_d2(c1, poly_u_sup);
                    let (_, t2, n2) = curve_tool::geom2d_d2(c1, poly_v_sup);
                    tan1 = t1;
                    tan2 = t2;
                    determine_transition_touch(
                        pos1, &mut tan1, &n1, &mut trans1, pos2, &mut tan2, &n2, &mut trans2,
                        tol_conf,
                    );
                }
                let pt_seg2 = IntRes2dIntersectionPoint::with_transitions(
                    &p1, poly_u_sup, poly_v_sup, &trans1, &trans2, false,
                );

                let oppos = !(tan1.dot(&tan2) > 0.0);
                if param_inf_on_curve1 > param_sup_on_curve1 {
                    self.base
                        .append_segment(&IntRes2dIntersectionSegment::from_two_points(
                            &pt_seg2, &pt_seg1, oppos, false,
                        ));
                } else {
                    self.base
                        .append_segment(&IntRes2dIntersectionSegment::from_two_points(
                            &pt_seg1, &pt_seg2, oppos, false,
                        ));
                }
            }
        } // end of processing of TangentZone

        self.base.done = true;
    }

}
