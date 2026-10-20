//! `findIntersect` (`gxx:1152-1568`).

use super::prelude::*;
use super::super::curve_tool;
use super::super::proj_p_cur;
use super::get_intersection::intersection;
use super::*;

impl Geom2dIntIntPolyPolyGen {
    /// `findIntersect(...)` (`gxx:1152-1568`).
    #[allow(clippy::too_many_arguments)]
    pub(super) fn find_intersect(
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
        poly1: &Geom2dIntPolygon2d,
        poly2: &Geom2dIntPolygon2d,
        is_full_polygon: bool,
    ) -> bool {
        let inter_pp = IntfInterferencePolygon2d::between(poly1, poly2);
        let mut eip = Geom2dIntExactIntersectionPoint::new(c1, c2, tol_conf);
        let mut u = 0.0;
        let mut v = 0.0;
        self.base.done = true; // To prevent exception in nbp=NbPoints();
        //----------------------------------------------------------------------
        //-- Processing of SectionPoint
        //----------------------------------------------------------------------
        let nbsp = inter_pp.base.nb_section_points();
        for sp in 1..=nbsp {
            let spnt = inter_pp.base.pnt_value(sp - 1);
            let (_, mut seg_index1, mut param_on1) = spnt.info_first();
            let (_, mut seg_index2, mut param_on2) = spnt.info_second();

            eip.perform(
                poly1,
                poly2,
                &mut seg_index1,
                &mut seg_index2,
                &mut param_on1,
                &mut param_on2,
            );
            let an_error_occurred = eip.an_error_occurred();

            if eip.nb_roots() == 0 && !is_full_polygon {
                return false;
            }

            if an_error_occurred {
                continue;
            }

            //--------------------------------------------------------------------
            //-- It is checked if the found point is really a root
            //--------------------------------------------------------------------
            let (uu, vv) = eip.roots();
            u = uu;
            v = vv;
            let (mut p1, mut tan1) = curve_tool::geom2d_d1(c1, u);
            let (mut p2, mut tan2) = curve_tool::geom2d_d1(c2, v);
            let mut dist = p1.distance(&p2);
            if eip.nb_roots() == 0 && dist > tol_conf {
                let mut trans = IntRes2dTransition::new();
                let mut pint =
                    IntRes2dIntersectionPoint::with_transitions(&p1, u, v, &trans, &trans, false);
                let t1f = poly1.approx_param_on_curve(seg_index1 as usize, 0.0);
                let t1l = poly1.approx_param_on_curve(seg_index1 as usize, 1.0);
                let t2f = poly2.approx_param_on_curve(seg_index2 as usize, 0.0);
                let t2l = poly2.approx_param_on_curve(seg_index2 as usize, 1.0);
                //
                let max_count = 16i32;
                let mut count = 0i32;
                intersection(
                    c1,
                    t1f,
                    t1l,
                    c2,
                    t2f,
                    t2l,
                    tol_conf,
                    max_count,
                    &mut pint,
                    &mut dist,
                    &mut count,
                );
                u = pint.param_on_first();
                v = pint.param_on_second();
                let r1 = curve_tool::geom2d_d1(c1, u);
                p1 = r1.0;
                tan1 = r1.1;
                let r2 = curve_tool::geom2d_d1(c2, v);
                p2 = r2.0;
                tan2 = r2.1;
                dist = p1.distance(&p2);
                let _ = (&mut trans, tan1, tan2);
            }
            //-----------------------------------------------------------------
            //-- It is checked if the point (u,v) does not exist already
            //--
            let nbp = self.base.nb_points();
            let eps_x1 = 10.0 * curve_tool::geom2d_eps_x();
            let eps_x2 = 10.0 * curve_tool::geom2d_eps_x();
            let mut p = 1usize;
            while p <= nbp {
                let ppt = self.base.point(p);
                if (u - ppt.param_on_first()).abs() <= eps_x1
                    && (v - ppt.param_on_second()).abs() <= eps_x2
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
                //-- Calculate the Positions of Points on the curve
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
                //-- Calculate the Transitions (see IntImpParGen.cxx)
                //-----------------------------------------------------------------
                if !determine_transition_simple(pos1, &tan1, &mut trans1, pos2, &tan2, &mut trans2, tol_conf)
                {
                    let (_, t1, n1) = curve_tool::geom2d_d2(c1, u);
                    let (_, t2, n2) = curve_tool::geom2d_d2(c2, v);
                    tan1 = t1;
                    tan2 = t2;
                    determine_transition_touch(
                        pos1, &mut tan1, &n1, &mut trans1, pos2, &mut tan2, &n2, &mut trans2,
                        tol_conf,
                    );
                }
                let ip =
                    IntRes2dIntersectionPoint::with_transitions(&p1, u, v, &trans1, &trans2, false);
                self.base.insert(&ip);
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
                let poly_u_inf_seg =
                    poly1.approx_param_on_curve(seg_index1on_p1 as usize, param_on_line);

                let (_, mut seg_index1on_p2, mut param_on_line) = spnt1.info_second();
                if seg_index1on_p2 > poly2.nb_segments() as i32 {
                    seg_index1on_p2 -= 1;
                    param_on_line = 1.0;
                }
                if seg_index1on_p2 <= 0 {
                    seg_index1on_p2 = 1;
                    param_on_line = 0.0;
                }
                let poly_v_inf_seg =
                    poly2.approx_param_on_curve(seg_index1on_p2 as usize, param_on_line);

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
            let p2a = curve_tool::geom2d_d0(c2, poly_v_inf);
            let distmemesens = p1a.square_distance(&p2a);
            let p2b = curve_tool::geom2d_d0(c2, poly_v_sup);
            let distdiffsens = p1a.square_distance(&p2b);
            if distmemesens > distdiffsens {
                std::mem::swap(&mut poly_v_inf, &mut poly_v_sup);
            }

            if ((poly1.deflection_over_estimation() > tol_conf)
                || (poly2.deflection_over_estimation() > tol_conf))
                && (nb_iter < NBITER_MAX_POLYGON)
            {
                let recurs_d1 = IntRes2dDomain::bounded(
                    &curve_tool::geom2d_value(c1, param_inf_on_curve1),
                    param_inf_on_curve1,
                    tol_conf,
                    &curve_tool::geom2d_value(c1, param_sup_on_curve1),
                    param_sup_on_curve1,
                    tol_conf,
                );
                let recurs_d2 = IntRes2dDomain::bounded(
                    &curve_tool::geom2d_value(c2, param_inf_on_curve2),
                    param_inf_on_curve2,
                    tol_conf,
                    &curve_tool::geom2d_value(c2, param_sup_on_curve2),
                    param_sup_on_curve2,
                    tol_conf,
                );
                //-- thePoly1(2) are not deleted,
                //-- finally they are destroyed.
                //-- !! No untimely return !!
                self.perform_poly(
                    c1,
                    &recurs_d1,
                    c2,
                    &recurs_d2,
                    tol,
                    tol_conf,
                    nb_iter + 1,
                    delta_u,
                    delta_v,
                );
            } else {
                //-----------------------------------------------------------------
                //-- Calculate Positions of Points on the curve and
                //-- Transitions on each limit of the segment
                //-----------------------------------------------------------------
                let mut pos1 = IntRes2dPosition::Middle;
                let mut pos2 = IntRes2dPosition::Middle;
                let mut trans1 = IntRes2dTransition::new();
                let mut trans2 = IntRes2dTransition::new();

                let (mut p1, mut tan1) = curve_tool::geom2d_d1(c1, poly_u_inf);
                let (mut p2, mut tan2) = curve_tool::geom2d_d1(c2, poly_v_inf);

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
                        c2,
                        &p1,
                        d2.first_parameter(),
                        d2.last_parameter(),
                        curve_tool::geom2d_eps_x(),
                    );
                } else if (param_inf_on_curve1 - param_sup_on_curve1).abs()
                    > (param_inf_on_curve2 - param_sup_on_curve2).abs()
                {
                    poly_v_inf = proj_p_cur::find_parameter_range(
                        c2,
                        &p1,
                        d2.first_parameter(),
                        d2.last_parameter(),
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

                if !determine_transition_simple(
                    pos1, &tan1, &mut trans1, pos2, &tan2, &mut trans2, tol_conf,
                ) {
                    let (_, t1, n1) = curve_tool::geom2d_d2(c1, poly_u_inf);
                    let (_, t2, n2) = curve_tool::geom2d_d2(c2, poly_v_inf);
                    tan1 = t1;
                    tan2 = t2;
                    determine_transition_touch(
                        pos1, &mut tan1, &n1, &mut trans1, pos2, &mut tan2, &n2, &mut trans2,
                        tol_conf,
                    );
                    let _ = (&mut p1, &mut p2);
                }
                let pt_seg1 = IntRes2dIntersectionPoint::with_transitions(
                    &p1, poly_u_inf, poly_v_inf, &trans1, &trans2, false,
                );
                //----------------------------------------------------------------------

                if (poly_u_inf - poly_u_sup).abs() <= curve_tool::geom2d_eps_x()
                    || (poly_v_inf - poly_v_sup).abs() <= curve_tool::geom2d_eps_x()
                {
                    self.base.insert(&pt_seg1);
                } else {
                    let (mut p1, mut tan1) = curve_tool::geom2d_d1(c1, poly_u_sup);
                    let (mut p2, mut tan2) = curve_tool::geom2d_d1(c2, poly_v_sup);
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
                            c2,
                            &p1,
                            d2.first_parameter(),
                            d2.last_parameter(),
                            curve_tool::geom2d_eps_x(),
                        );
                    } else if (param_inf_on_curve1 - param_sup_on_curve1).abs()
                        > (param_inf_on_curve2 - param_sup_on_curve2).abs()
                    {
                        poly_v_sup = proj_p_cur::find_parameter_range(
                            c2,
                            &p1,
                            d2.first_parameter(),
                            d2.last_parameter(),
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
                        let (_, t2, n2) = curve_tool::geom2d_d2(c2, poly_v_sup);
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
            }
        }
        true
    }

}
