//! Port of `IntCurve_IntPolyPolyGen` instantiated as
//! `Geom2dInt_TheIntPCurvePCurveOfGInter`
//! (`IntCurve_IntPolyPolyGen.gxx:1-1797`,
//! `Geom2dInt_TheIntPCurvePCurveOfGInter.hxx:35-102`,
//! `Geom2dInt_TheIntPCurvePCurveOfGInter_0.cxx:27-...`).
//!
//! The general polygon-approximation curve/curve intersector: it discretizes
//! both curves into `Geom2dIntPolygon2d`, runs `Intf_InterferencePolygon2d` on
//! them, refines every interference point with
//! `Geom2dIntExactIntersectionPoint` and turns tangency zones into
//! intersection segments (recursively refining the zone's sub-domain while the
//! polygons are still coarser than `TolConf`).
//!
//! `Perform(C1, D1, C2, D2, TolConf, Tol)` (`gxx:94-292`),
//! `Perform(C1, D1, TolConf, Tol)` (`gxx:296-387`),
//! `Perform(..., NbIter, DeltaU, DeltaV)` (`gxx:390-870` and `:1040-1146`),
//! `findIntersect` (`gxx:1152-1568`), `GetIntersection` (`gxx:1573-1783`),
//! `GetMinNbSamples` / `SetMinNbSamples` (`gxx:1787-1797`).

use occt_core::bnd::box2d::BndBox2d;
use occt_core::gp::GpPnt2d;
use occt_core::intf::{IntfInterferencePolygon2d, IntfPolygon2d};
use occt_core::intimpargen::gen::{determine_transition_simple, determine_transition_touch};
use occt_core::intres2d::{
    IntRes2dDomain, IntRes2dIntersection, IntRes2dIntersectionPoint,
    IntRes2dIntersectionSegment, IntRes2dPosition, IntRes2dTransition,
};
use occt_core::precision::{epsilon, PCONFUSION};

use crate::curve::Curve2d;
use super::curve_tool;
use super::exact_intersection_point::Geom2dIntExactIntersectionPoint;
use super::polygon2d::Geom2dIntPolygon2d;
use super::proj_p_cur;

/// `NBITER_MAX_POLYGON` (`gxx:50`).
const NBITER_MAX_POLYGON: i32 = 10;
/// `TOL_CONF_MINI` (`gxx:51`).
const TOL_CONF_MINI: f64 = 0.0000000001;
/// `TOL_MINI` (`gxx:52`).
const TOL_MINI: f64 = 0.0000000001;

/// `IntCurve_IntPolyPolyGen` (`Geom2dInt_TheIntPCurvePCurveOfGInter.hxx:35-102`),
/// carrying the `IntRes2d_Intersection` result.
#[derive(Clone, Debug)]
pub struct Geom2dIntIntPolyPolyGen {
    /// The `IntRes2d_Intersection` base class.
    pub(crate) base: IntRes2dIntersection,
    /// `DomainOnCurve1` (`hxx:96`).
    domain_on_curve1: IntRes2dDomain,
    /// `DomainOnCurve2` (`hxx:97`).
    domain_on_curve2: IntRes2dDomain,
    /// `myMinPntNb` (`hxx:100`).
    my_min_pnt_nb: usize,
}

impl Default for Geom2dIntIntPolyPolyGen {
    fn default() -> Self {
        Self::new()
    }
}

impl Geom2dIntIntPolyPolyGen {
    /// `IntCurve_IntPolyPolyGen()` (`gxx:86-90`).
    pub fn new() -> Self {
        // Minimum number of samples (`gxx:87`).
        let a_min_pnt_nb = 20;
        Self {
            base: IntRes2dIntersection::new(),
            domain_on_curve1: IntRes2dDomain::new(),
            domain_on_curve2: IntRes2dDomain::new(),
            my_min_pnt_nb: a_min_pnt_nb,
        }
    }

    /// The `IntRes2d_Intersection` base (`Geom2dInt_TheIntPCurvePCurveOfGInter`
    /// derives from it).
    pub fn result(&self) -> &IntRes2dIntersection {
        &self.base
    }

    /// `IntRes2d_Intersection::ResetFields` (inherited, `gxx:100`).
    fn reset_fields(&mut self) {
        self.base.reset_fields();
    }

    /// `GetMinNbSamples()` (`gxx:1787-1790`).
    pub fn get_min_nb_samples(&self) -> usize {
        self.my_min_pnt_nb
    }

    /// `SetMinNbSamples(theMinNbSamples)` (`gxx:1794-1797`).
    pub fn set_min_nb_samples(&mut self, the_min_nb_samples: usize) {
        self.my_min_pnt_nb = the_min_nb_samples;
    }

    /// `Perform(C1, D1, C2, D2, TolConf, Tol)` (`gxx:94-292`).
    pub fn perform(
        &mut self,
        c1: &dyn Curve2d,
        d1: &IntRes2dDomain,
        c2: &dyn Curve2d,
        d2: &IntRes2dDomain,
        the_tol_conf: f64,
        the_tol: f64,
    ) {
        self.reset_fields();
        self.domain_on_curve1 = *d1;
        self.domain_on_curve2 = *d2;
        let du = d1.last_parameter() - d1.first_parameter();
        let dv = d2.last_parameter() - d2.first_parameter();
        let tl = if the_tol < TOL_MINI { TOL_MINI } else { the_tol };
        let tl_conf = if the_tol_conf < TOL_CONF_MINI { TOL_CONF_MINI } else { the_tol_conf };
        self.perform_poly(c1, d1, c2, d2, tl_conf, tl, 0, du, dv);
        //----------------------------------------------------------------------
        //-- Processing of end points
        //----------------------------------------------------------------------
        let mut head_on1 = false;
        let mut head_on2 = false;
        let mut end_on1 = false;
        let mut end_on2 = false;

        //--------------------------------------------------------------------
        //-- The points Head Head ... End End are not rejected if
        //-- they are already present at the end of segment
        //-- PosSegment =            1    if Head Head
        //--                       2      if Head End
        //--                     4        if End  Head
        //--                   8          if End  End
        //--------------------------------------------------------------------
        let mut pos_segment: i32 = 0;

        let n = self.base.nb_points();
        for i in 1..=n {
            let pos1 = self.base.point(i).transition_of_first().position_on_curve();
            if pos1 == IntRes2dPosition::Head {
                head_on1 = true;
            } else if pos1 == IntRes2dPosition::End {
                end_on1 = true;
            }

            let pos2 = self.base.point(i).transition_of_second().position_on_curve();
            if pos2 == IntRes2dPosition::Head {
                head_on2 = true;
            } else if pos2 == IntRes2dPosition::End {
                end_on2 = true;
            }

            pos_segment |= position_bits(pos1, pos2);
        }

        let n = self.base.nb_segments();
        for i in 1..=n {
            let seg = self.base.segment(i);
            let pos1 = seg.first_point().transition_of_first().position_on_curve();
            if pos1 == IntRes2dPosition::Head {
                head_on1 = true;
            } else if pos1 == IntRes2dPosition::End {
                end_on1 = true;
            }
            let pos2 = seg.first_point().transition_of_second().position_on_curve();
            if pos2 == IntRes2dPosition::Head {
                head_on2 = true;
            } else if pos2 == IntRes2dPosition::End {
                end_on2 = true;
            }
            pos_segment |= position_bits(pos1, pos2);

            let pos1 = seg.last_point().transition_of_first().position_on_curve();
            if pos1 == IntRes2dPosition::Head {
                head_on1 = true;
            } else if pos1 == IntRes2dPosition::End {
                end_on1 = true;
            }
            let pos2 = seg.last_point().transition_of_second().position_on_curve();
            if pos2 == IntRes2dPosition::Head {
                head_on2 = true;
            } else if pos2 == IntRes2dPosition::End {
                end_on2 = true;
            }
            pos_segment |= position_bits(pos1, pos2);
        }

        let u0 = d1.first_parameter();
        let u1 = d1.last_parameter();
        let v0 = d2.first_parameter();
        let v1 = d2.last_parameter();
        let mut int_pt = IntRes2dIntersectionPoint::new();

        if d1.first_tolerance() != 0.0 || d2.first_tolerance() != 0.0 {
            if head_or_end_point(
                d1, c1, u0, d2, c2, v0, the_tol_conf, &mut int_pt, &mut head_on1, &mut head_on2,
                &mut end_on1, &mut end_on2, pos_segment,
            ) {
                self.base.insert(&int_pt);
            }
        }
        if d1.first_tolerance() != 0.0 || d2.last_tolerance() != 0.0 {
            if head_or_end_point(
                d1, c1, u0, d2, c2, v1, the_tol_conf, &mut int_pt, &mut head_on1, &mut head_on2,
                &mut end_on1, &mut end_on2, pos_segment,
            ) {
                self.base.insert(&int_pt);
            }
        }
        if d1.last_tolerance() != 0.0 || d2.first_tolerance() != 0.0 {
            if head_or_end_point(
                d1, c1, u1, d2, c2, v0, the_tol_conf, &mut int_pt, &mut head_on1, &mut head_on2,
                &mut end_on1, &mut end_on2, pos_segment,
            ) {
                self.base.insert(&int_pt);
            }
        }
        if d1.last_tolerance() != 0.0 || d2.last_tolerance() != 0.0 {
            if head_or_end_point(
                d1, c1, u1, d2, c2, v1, the_tol_conf, &mut int_pt, &mut head_on1, &mut head_on2,
                &mut end_on1, &mut end_on2, pos_segment,
            ) {
                self.base.insert(&int_pt);
            }
        }
    }

    /// `Perform(C1, D1, TolConf, Tol)` (`gxx:296-387`), the auto-intersection
    /// overload.
    pub fn perform_self(
        &mut self,
        c1: &dyn Curve2d,
        d1: &IntRes2dDomain,
        the_tol_conf: f64,
        the_tol: f64,
    ) {
        self.reset_fields();
        self.domain_on_curve1 = *d1;
        self.domain_on_curve2 = *d1;
        let du = d1.last_parameter() - d1.first_parameter();
        let tl = if the_tol < TOL_MINI { TOL_MINI } else { the_tol };
        let tl_conf = if the_tol_conf < TOL_CONF_MINI { TOL_CONF_MINI } else { the_tol_conf };
        self.perform_self_poly(c1, d1, tl_conf, tl, 0, du, du);

        //--------------------------------------------------------------------
        //-- PosSegment =            1    if Head Head
        //--                       2      if Head End
        //--                     4        if End  Head
        //--                   8          if End  End
        //--------------------------------------------------------------------
        let mut pos_segment: i32 = 0;

        let n = self.base.nb_points();
        for i in 1..=n {
            let pos1 = self.base.point(i).transition_of_first().position_on_curve();
            let pos2 = self.base.point(i).transition_of_second().position_on_curve();
            pos_segment |= position_bits(pos1, pos2);
        }

        let n = self.base.nb_segments();
        for i in 1..=n {
            let seg = self.base.segment(i);
            let pos1 = seg.first_point().transition_of_first().position_on_curve();
            let pos2 = seg.first_point().transition_of_second().position_on_curve();
            pos_segment |= position_bits(pos1, pos2);

            let pos1 = seg.last_point().transition_of_first().position_on_curve();
            let pos2 = seg.last_point().transition_of_second().position_on_curve();
            pos_segment |= position_bits(pos1, pos2);
        }
        let _ = pos_segment;
    }

    /// `Perform(C1, D1, C2, D2, TolConf, Tol, NbIter, DeltaU, DeltaV)`
    /// (`gxx:1040-1146`), the polygon engine of the two-curve case.
    fn perform_poly(
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
            nbsamples_on_c1 = nbsamples_on_c1.max(self.my_min_pnt_nb);
            nbsamples_on_c2 = nbsamples_on_c2.max(self.my_min_pnt_nb);
        } else {
            // Increase number of samples in second and next iterations.
            nbsamples_on_c1 = (5 * (nbsamples_on_c1 * nb_iter as usize)) / 4;
            nbsamples_on_c2 = (5 * (nbsamples_on_c2 * nb_iter as usize)) / 4;
        }

        let mut a_poly1 = Geom2dIntPolygon2d::new(c1, nbsamples_on_c1, d1, tol);
        let mut a_poly2 = Geom2dIntPolygon2d::new(c2, nbsamples_on_c2, d2, tol);

        if (a_poly1.deflection_over_estimation() > tol_conf)
            && (a_poly2.deflection_over_estimation() > tol_conf)
        {
            let a_deflection_sum = a_poly1.deflection_over_estimation().max(tol_conf)
                + a_poly2.deflection_over_estimation().max(tol_conf);

            if nbsamples_on_c2 > nbsamples_on_c1 {
                a_poly2.compute_with_box(c2, a_poly1.bounding());
                a_poly1.set_deflection_over_estimation(a_deflection_sum);
                a_poly1.compute_with_box(c1, a_poly2.bounding());
            } else {
                a_poly1.compute_with_box(c1, a_poly2.bounding());
                a_poly2.set_deflection_over_estimation(a_deflection_sum);
                a_poly2.compute_with_box(c2, a_poly1.bounding());
            }
        }

        //----------------------------------------------------------------------
        //-- if the deflection less then the Tolerance of Confusion
        //-- Then the deflection of the polygon is set in TolConf
        //-- (Detection of Tangency Zones)
        //----------------------------------------------------------------------
        if a_poly1.deflection_over_estimation() < tol_conf {
            a_poly1.set_deflection_over_estimation(tol_conf);
        }
        if a_poly2.deflection_over_estimation() < tol_conf {
            a_poly2.set_deflection_over_estimation(tol_conf);
        }

        // for case when a few polygon points were replaced by line
        // if exact solution was not found
        // then search of precise solution will be repeated
        // for polygon contains all initial points
        // secondary search will be performed only for case when initial points
        // were dropped
        let is_full_representation = (a_poly1.nb_segments() == nbsamples_on_c1
            && a_poly2.nb_segments() == nbsamples_on_c2);

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
            &a_poly1,
            &a_poly2,
            is_full_representation,
        ) && !is_full_representation
        {
            if a_poly1.nb_segments() < nbsamples_on_c1 {
                a_poly1 = Geom2dIntPolygon2d::new(c1, nbsamples_on_c1, d1, tol);
            }
            if a_poly2.nb_segments() < nbsamples_on_c2 {
                a_poly2 = Geom2dIntPolygon2d::new(c2, nbsamples_on_c2, d2, tol);
            }

            self.find_intersect(
                c1, d1, c2, d2, tol_conf, tol, nb_iter, delta_u, delta_v, &a_poly1, &a_poly2, true,
            );
        }

        self.base.done = true;
    }

    /// `Perform(C1, D1, TolConf, Tol, NbIter, DeltaU, DeltaV)` (`gxx:390-870`),
    /// the polygon engine of the auto-intersection case.
    fn perform_self_poly(
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

                            let (mut p1, tan1) = curve_tool::d1(c1, u);
                            let (mut p2, tan2) = curve_tool::d1(c1, v);
                            let mut dist = p1.distance(&p2);
                            let eps_x1 = 10.0 * curve_tool::eps_x();

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
                                    let (pp1, t1, n1) = curve_tool::d2(c1, u);
                                    let (pp2, t2, n2) = curve_tool::d2(c1, v);
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
                let spnt1 = zone.get_point(qq);
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

            let p1a = curve_tool::d0(c1, poly_u_inf);
            let p2a = curve_tool::d0(c1, poly_v_inf);
            let distmemesens = p1a.square_distance(&p2a);
            let p2b = curve_tool::d0(c1, poly_v_sup);
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

            let (mut p1, mut tan1) = curve_tool::d1(c1, poly_u_inf);
            let (mut p2, mut tan2) = curve_tool::d1(c1, poly_v_inf);

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
                    curve_tool::eps_x(),
                );
            } else if pos1 != IntRes2dPosition::Middle && pos2 == IntRes2dPosition::Middle {
                poly_v_inf = proj_p_cur::find_parameter_range(
                    c1,
                    &p1,
                    d1.first_parameter(),
                    d1.last_parameter(),
                    curve_tool::eps_x(),
                );
            } else if (param_inf_on_curve1 - param_sup_on_curve1).abs()
                > (param_inf_on_curve2 - param_sup_on_curve2).abs()
            {
                poly_v_inf = proj_p_cur::find_parameter_range(
                    c1,
                    &p1,
                    d1.first_parameter(),
                    d1.last_parameter(),
                    curve_tool::eps_x(),
                );
            } else {
                poly_u_inf = proj_p_cur::find_parameter_range(
                    c1,
                    &p2,
                    d1.first_parameter(),
                    d1.last_parameter(),
                    curve_tool::eps_x(),
                );
            }

            if !determine_transition_simple(pos1, &tan1, &mut trans1, pos2, &tan2, &mut trans2, tol_conf)
            {
                let (pp1, t1, n1) = curve_tool::d2(c1, poly_u_inf);
                let (pp2, t2, n2) = curve_tool::d2(c1, poly_v_inf);
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

            if (poly_u_inf - poly_u_sup).abs() <= curve_tool::eps_x()
                || (poly_v_inf - poly_v_sup).abs() <= curve_tool::eps_x()
            {
                // bad segment
            } else {
                let (mut p1, mut tan1) = curve_tool::d1(c1, poly_u_sup);
                let (mut p2, mut tan2) = curve_tool::d1(c1, poly_v_sup);
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
                        curve_tool::eps_x(),
                    );
                } else if pos1 != IntRes2dPosition::Middle && pos2 == IntRes2dPosition::Middle {
                    poly_v_sup = proj_p_cur::find_parameter_range(
                        c1,
                        &p1,
                        d1.first_parameter(),
                        d1.last_parameter(),
                        curve_tool::eps_x(),
                    );
                } else if (param_inf_on_curve1 - param_sup_on_curve1).abs()
                    > (param_inf_on_curve2 - param_sup_on_curve2).abs()
                {
                    poly_v_sup = proj_p_cur::find_parameter_range(
                        c1,
                        &p1,
                        d1.first_parameter(),
                        d1.last_parameter(),
                        curve_tool::eps_x(),
                    );
                } else {
                    poly_u_sup = proj_p_cur::find_parameter_range(
                        c1,
                        &p2,
                        d1.first_parameter(),
                        d1.last_parameter(),
                        curve_tool::eps_x(),
                    );
                }

                if !determine_transition_simple(
                    pos1, &tan1, &mut trans1, pos2, &tan2, &mut trans2, tol_conf,
                ) {
                    let (_, t1, n1) = curve_tool::d2(c1, poly_u_sup);
                    let (_, t2, n2) = curve_tool::d2(c1, poly_v_sup);
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

    /// `findIntersect(...)` (`gxx:1152-1568`).
    #[allow(clippy::too_many_arguments)]
    fn find_intersect(
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
        the_poly1: &Geom2dIntPolygon2d,
        the_poly2: &Geom2dIntPolygon2d,
        is_full_polygon: bool,
    ) -> bool {
        let inter_pp = IntfInterferencePolygon2d::between(the_poly1, the_poly2);
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
                the_poly1,
                the_poly2,
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
            let (mut p1, mut tan1) = curve_tool::d1(c1, u);
            let (mut p2, mut tan2) = curve_tool::d1(c2, v);
            let mut dist = p1.distance(&p2);
            if eip.nb_roots() == 0 && dist > tol_conf {
                let mut trans = IntRes2dTransition::new();
                let mut a_pint =
                    IntRes2dIntersectionPoint::with_transitions(&p1, u, v, &trans, &trans, false);
                let a_t1f = the_poly1.approx_param_on_curve(seg_index1 as usize, 0.0);
                let a_t1l = the_poly1.approx_param_on_curve(seg_index1 as usize, 1.0);
                let a_t2f = the_poly2.approx_param_on_curve(seg_index2 as usize, 0.0);
                let a_t2l = the_poly2.approx_param_on_curve(seg_index2 as usize, 1.0);
                //
                let a_max_count = 16i32;
                let mut a_count = 0i32;
                get_intersection(
                    c1,
                    a_t1f,
                    a_t1l,
                    c2,
                    a_t2f,
                    a_t2l,
                    tol_conf,
                    a_max_count,
                    &mut a_pint,
                    &mut dist,
                    &mut a_count,
                );
                u = a_pint.param_on_first();
                v = a_pint.param_on_second();
                let r1 = curve_tool::d1(c1, u);
                p1 = r1.0;
                tan1 = r1.1;
                let r2 = curve_tool::d1(c2, v);
                p2 = r2.0;
                tan2 = r2.1;
                dist = p1.distance(&p2);
                let _ = (&mut trans, tan1, tan2);
            }
            //-----------------------------------------------------------------
            //-- It is checked if the point (u,v) does not exist already
            //--
            let nbp = self.base.nb_points();
            let eps_x1 = 10.0 * curve_tool::eps_x();
            let eps_x2 = 10.0 * curve_tool::eps_x();
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
                    let (_, t1, n1) = curve_tool::d2(c1, u);
                    let (_, t2, n2) = curve_tool::d2(c2, v);
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
                let spnt1 = zone.get_point(qq);
                let (_, mut seg_index1on_p1, mut param_on_line) = spnt1.info_first();
                if seg_index1on_p1 > the_poly1.nb_segments() as i32 {
                    seg_index1on_p1 -= 1;
                    param_on_line = 1.0;
                }
                if seg_index1on_p1 <= 0 {
                    seg_index1on_p1 = 1;
                    param_on_line = 0.0;
                }
                let poly_u_inf_seg =
                    the_poly1.approx_param_on_curve(seg_index1on_p1 as usize, param_on_line);

                let (_, mut seg_index1on_p2, mut param_on_line) = spnt1.info_second();
                if seg_index1on_p2 > the_poly2.nb_segments() as i32 {
                    seg_index1on_p2 -= 1;
                    param_on_line = 1.0;
                }
                if seg_index1on_p2 <= 0 {
                    seg_index1on_p2 = 1;
                    param_on_line = 0.0;
                }
                let poly_v_inf_seg =
                    the_poly2.approx_param_on_curve(seg_index1on_p2 as usize, param_on_line);

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

            let p1a = curve_tool::d0(c1, poly_u_inf);
            let p2a = curve_tool::d0(c2, poly_v_inf);
            let distmemesens = p1a.square_distance(&p2a);
            let p2b = curve_tool::d0(c2, poly_v_sup);
            let distdiffsens = p1a.square_distance(&p2b);
            if distmemesens > distdiffsens {
                std::mem::swap(&mut poly_v_inf, &mut poly_v_sup);
            }

            if ((the_poly1.deflection_over_estimation() > tol_conf)
                || (the_poly2.deflection_over_estimation() > tol_conf))
                && (nb_iter < NBITER_MAX_POLYGON)
            {
                let recurs_d1 = IntRes2dDomain::bounded(
                    &curve_tool::value(c1, param_inf_on_curve1),
                    param_inf_on_curve1,
                    tol_conf,
                    &curve_tool::value(c1, param_sup_on_curve1),
                    param_sup_on_curve1,
                    tol_conf,
                );
                let recurs_d2 = IntRes2dDomain::bounded(
                    &curve_tool::value(c2, param_inf_on_curve2),
                    param_inf_on_curve2,
                    tol_conf,
                    &curve_tool::value(c2, param_sup_on_curve2),
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

                let (mut p1, mut tan1) = curve_tool::d1(c1, poly_u_inf);
                let (mut p2, mut tan2) = curve_tool::d1(c2, poly_v_inf);

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
                        curve_tool::eps_x(),
                    );
                } else if pos1 != IntRes2dPosition::Middle && pos2 == IntRes2dPosition::Middle {
                    poly_v_inf = proj_p_cur::find_parameter_range(
                        c2,
                        &p1,
                        d2.first_parameter(),
                        d2.last_parameter(),
                        curve_tool::eps_x(),
                    );
                } else if (param_inf_on_curve1 - param_sup_on_curve1).abs()
                    > (param_inf_on_curve2 - param_sup_on_curve2).abs()
                {
                    poly_v_inf = proj_p_cur::find_parameter_range(
                        c2,
                        &p1,
                        d2.first_parameter(),
                        d2.last_parameter(),
                        curve_tool::eps_x(),
                    );
                } else {
                    poly_u_inf = proj_p_cur::find_parameter_range(
                        c1,
                        &p2,
                        d1.first_parameter(),
                        d1.last_parameter(),
                        curve_tool::eps_x(),
                    );
                }

                if !determine_transition_simple(
                    pos1, &tan1, &mut trans1, pos2, &tan2, &mut trans2, tol_conf,
                ) {
                    let (_, t1, n1) = curve_tool::d2(c1, poly_u_inf);
                    let (_, t2, n2) = curve_tool::d2(c2, poly_v_inf);
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

                if (poly_u_inf - poly_u_sup).abs() <= curve_tool::eps_x()
                    || (poly_v_inf - poly_v_sup).abs() <= curve_tool::eps_x()
                {
                    self.base.insert(&pt_seg1);
                } else {
                    let (mut p1, mut tan1) = curve_tool::d1(c1, poly_u_sup);
                    let (mut p2, mut tan2) = curve_tool::d1(c2, poly_v_sup);
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
                            curve_tool::eps_x(),
                        );
                    } else if pos1 != IntRes2dPosition::Middle && pos2 == IntRes2dPosition::Middle {
                        poly_v_sup = proj_p_cur::find_parameter_range(
                            c2,
                            &p1,
                            d2.first_parameter(),
                            d2.last_parameter(),
                            curve_tool::eps_x(),
                        );
                    } else if (param_inf_on_curve1 - param_sup_on_curve1).abs()
                        > (param_inf_on_curve2 - param_sup_on_curve2).abs()
                    {
                        poly_v_sup = proj_p_cur::find_parameter_range(
                            c2,
                            &p1,
                            d2.first_parameter(),
                            d2.last_parameter(),
                            curve_tool::eps_x(),
                        );
                    } else {
                        poly_u_sup = proj_p_cur::find_parameter_range(
                            c1,
                            &p2,
                            d1.first_parameter(),
                            d1.last_parameter(),
                            curve_tool::eps_x(),
                        );
                    }

                    if !determine_transition_simple(
                        pos1, &tan1, &mut trans1, pos2, &tan2, &mut trans2, tol_conf,
                    ) {
                        let (_, t1, n1) = curve_tool::d2(c1, poly_u_sup);
                        let (_, t2, n2) = curve_tool::d2(c2, poly_v_sup);
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

/// The `PosSegment |= ...` table shared by both `Perform` wrappers
/// (`gxx:145-157` / `:176-188` / `:201-213`).
#[inline]
fn position_bits(pos1: IntRes2dPosition, pos2: IntRes2dPosition) -> i32 {
    let mut bits = 0i32;
    if pos1 == IntRes2dPosition::Head {
        if pos2 == IntRes2dPosition::Head {
            bits |= 1;
        } else if pos2 == IntRes2dPosition::End {
            bits |= 2;
        }
    } else if pos1 == IntRes2dPosition::End {
        if pos2 == IntRes2dPosition::Head {
            bits |= 4;
        } else if pos2 == IntRes2dPosition::End {
            bits |= 8;
        }
    }
    bits
}

/// `HeadOrEndPoint(...)` (`gxx:871-1032`).
#[allow(clippy::too_many_arguments)]
fn head_or_end_point(
    d1: &IntRes2dDomain,
    c1: &dyn Curve2d,
    tu: f64,
    d2: &IntRes2dDomain,
    c2: &dyn Curve2d,
    tv: f64,
    tol_conf: f64,
    int_pt: &mut IntRes2dIntersectionPoint,
    head_on1: &mut bool,
    head_on2: &mut bool,
    end_on1: &mut bool,
    end_on2: &mut bool,
    pos_segment: i32,
) -> bool {
    let mut u = tu;
    let mut v = tv;
    let svu = u;
    let svv = v;

    let (mut p1, mut t1) = curve_tool::d1(c1, u);
    let (mut p2, mut t2) = curve_tool::d1(c2, v);

    let mut pos1 = IntRes2dPosition::Middle;
    let mut pos2 = IntRes2dPosition::Middle;
    let mut trans1 = IntRes2dTransition::new();
    let mut trans2 = IntRes2dTransition::new();
    let mut sp1 = GpPnt2d::zero();
    let mut sp2 = GpPnt2d::zero();

    //----------------------------------------------------------------------
    //-- Head On 1   :        Head1 <-> P2
    if p2.distance(d1.first_point()) <= d1.first_tolerance() {
        pos1 = IntRes2dPosition::Head;
        *head_on1 = true;
        sp1 = *d1.first_point();
        u = d1.first_parameter();
    }
    //----------------------------------------------------------------------
    //-- End On 1   :         End1 <-> P2
    else if p2.distance(d1.last_point()) <= d1.last_tolerance() {
        pos1 = IntRes2dPosition::End;
        *end_on1 = true;
        sp1 = *d1.last_point();
        u = d1.last_parameter();
    }
    //----------------------------------------------------------------------
    //-- Head On 2   :        Head2 <-> P1
    else if p1.distance(d2.first_point()) <= d2.first_tolerance() {
        pos2 = IntRes2dPosition::Head;
        *head_on2 = true;
        sp2 = *d2.first_point();
        v = d2.first_parameter();
    }
    //----------------------------------------------------------------------
    //-- End On 2   :         End2 <-> P1
    else if p1.distance(d2.last_point()) <= d2.last_tolerance() {
        pos2 = IntRes2dPosition::End;
        *end_on2 = true;
        sp2 = *d2.last_point();
        v = d2.last_parameter();
    }

    let eps_x1 = curve_tool::eps_x();
    let eps_x2 = curve_tool::eps_x();

    if (pos1 != IntRes2dPosition::Middle) || (pos2 != IntRes2dPosition::Middle) {
        if pos1 == IntRes2dPosition::Middle {
            if (u - d1.first_parameter()).abs() <= eps_x1 {
                pos1 = IntRes2dPosition::Head;
                p1 = *d1.first_point();
                *head_on1 = true;
            } else if (u - d1.last_parameter()).abs() <= eps_x1 {
                pos1 = IntRes2dPosition::End;
                p1 = *d1.last_point();
                *end_on1 = true;
            }
        } else if u != tu {
            p1 = sp1;
        }

        if pos2 == IntRes2dPosition::Middle {
            if (v - d2.first_parameter()).abs() <= eps_x2 {
                pos2 = IntRes2dPosition::Head;
                *head_on2 = true;
                p2 = *d2.first_point();
                if pos1 != IntRes2dPosition::Middle {
                    p1.set_coord(0.5 * (p1.x() + p2.x()), 0.5 * (p1.y() + p2.y()));
                } else {
                    p2 = p1;
                }
            } else if (v - d2.last_parameter()).abs() <= eps_x2 {
                pos2 = IntRes2dPosition::End;
                *end_on2 = true;
                p2 = *d2.last_point();
                if pos1 != IntRes2dPosition::Middle {
                    p1.set_coord(0.5 * (p1.x() + p2.x()), 0.5 * (p1.y() + p2.y()));
                } else {
                    p2 = p1;
                }
            }
        }

        //--------------------------------------------------------------------
        //-- It is tested if a point at the end of segment already has its transitions
        //-- If Yes, the new point is not created
        //--
        //-- PosSegment =            1    if Head Head
        //--                       2      if Head End
        //--                     4        if End  Head
        //--                   8          if End  End
        //--------------------------------------------------------------------
        if pos1 == IntRes2dPosition::Head {
            if (pos2 == IntRes2dPosition::Head) && (pos_segment & 1) != 0 {
                return false;
            }
            if (pos2 == IntRes2dPosition::End) && (pos_segment & 2) != 0 {
                return false;
            }
        } else if pos1 == IntRes2dPosition::End {
            if (pos2 == IntRes2dPosition::Head) && (pos_segment & 4) != 0 {
                return false;
            }
            if (pos2 == IntRes2dPosition::End) && (pos_segment & 8) != 0 {
                return false;
            }
        }

        if !determine_transition_simple(pos1, &t1, &mut trans1, pos2, &t2, &mut trans2, tol_conf) {
            let (_, tt1, n1) = curve_tool::d2(c1, svu);
            let (_, tt2, n2) = curve_tool::d2(c2, svv);
            t1 = tt1;
            t2 = tt2;
            determine_transition_touch(
                pos1, &mut t1, &n1, &mut trans1, pos2, &mut t2, &n2, &mut trans2, tol_conf,
            );
        }
        int_pt.set_values(&p1, u, v, &trans1, &trans2, false);
        true
    } else {
        let _ = (&mut sp1, &mut sp2, &mut t1, &mut t2);
        false
    }
}

/// `GetIntersection(...)` (`gxx:1573-1783`).
#[allow(clippy::too_many_arguments)]
fn get_intersection(
    the_c1: &dyn Curve2d,
    the_t1f: f64,
    the_t1l: f64,
    the_c2: &dyn Curve2d,
    the_t2f: f64,
    the_t2l: f64,
    the_tol_conf: f64,
    the_max_count: i32,
    the_p_int: &mut IntRes2dIntersectionPoint,
    the_dist: &mut f64,
    the_count: &mut i32,
) {
    *the_count += 1;
    //
    let a_tol2 = the_tol_conf * the_tol_conf;
    let a_p_tol1 = (100.0 * epsilon(the_t1f.abs().max(the_t1l.abs()))).max(PCONFUSION);
    let a_p_tol2 = (100.0 * epsilon(the_t2f.abs().max(the_t2l.abs()))).max(PCONFUSION);
    //
    let a_p1f = curve_tool::d0(the_c1, the_t1f);
    let a_p1l = curve_tool::d0(the_c1, the_t1l);
    let mut a_b1 = BndBox2d::new();
    a_b1.add_point(&a_p1f);
    a_b1.add_point(&a_p1l);
    a_b1.enlarge(the_tol_conf);
    //
    let a_p2f = curve_tool::d0(the_c2, the_t2f);
    let a_p2l = curve_tool::d0(the_c2, the_t2l);
    let mut a_b2 = BndBox2d::new();
    a_b2.add_point(&a_p2f);
    a_b2.add_point(&a_p2l);
    a_b2.enlarge(the_tol_conf);
    //
    if a_b1.is_out_box(&a_b2) {
        *the_count -= 1;
        return;
    }
    //
    let is_small1 = (the_t1l - the_t1f) <= a_p_tol1
        || a_p1f.square_distance(&a_p1l) / 4.0 <= a_tol2;
    let is_small2 = (the_t2l - the_t2f) <= a_p_tol2
        || a_p2f.square_distance(&a_p2l) / 4.0 <= a_tol2;

    if (is_small1 && is_small2) || (*the_count > the_max_count) {
        // Seems to be intersection
        // Simple treatment of segment intersection
        let a_pnts1 = [
            *a_p1f.xy(),
            a_p1f.xy().added(&a_p1l.xy()).divided(2.0),
            *a_p1l.xy(),
        ];
        let a_pnts2 = [
            *a_p2f.xy(),
            a_p2f.xy().added(&a_p2l.xy()).divided(2.0),
            *a_p2l.xy(),
        ];
        let mut imin: i32 = -1;
        let mut jmin: i32 = -1;
        let mut dmin = f64::MAX;
        for i in 0..3 {
            for j in 0..3 {
                let d: f64 = a_pnts1[i].subtracted(&a_pnts2[j]).square_modulus();
                if d < dmin {
                    dmin = d;
                    imin = i as i32;
                    jmin = j as i32;
                }
            }
        }
        //
        dmin = dmin.sqrt();
        if *the_dist > dmin {
            *the_dist = dmin;
            //
            let t1 = match imin {
                0 => the_t1f,
                1 => (the_t1f + the_t1l) / 2.0,
                _ => the_t1l,
            };
            //
            let t2 = match jmin {
                0 => the_t2f,
                1 => (the_t2f + the_t2l) / 2.0,
                _ => the_t2l,
            };
            //
            let a_pint = GpPnt2d::from_xy(a_pnts1[imin as usize].added(&a_pnts2[jmin as usize]).divided(2.0));
            //
            let a_trans1 = IntRes2dTransition::new();
            let a_trans2 = IntRes2dTransition::new();
            the_p_int.set_values(&a_pint, t1, t2, &a_trans1, &a_trans2, false);
        }
        *the_count -= 1;
        return;
    }

    if is_small1 {
        let a_t2m = (the_t2l + the_t2f) / 2.0;
        get_intersection(
            the_c1, the_t1f, the_t1l, the_c2, the_t2f, a_t2m, the_tol_conf, the_max_count, the_p_int,
            the_dist, the_count,
        );
        get_intersection(
            the_c1, the_t1f, the_t1l, the_c2, a_t2m, the_t2l, the_tol_conf, the_max_count, the_p_int,
            the_dist, the_count,
        );
    } else if is_small2 {
        let a_t1m = (the_t1l + the_t1f) / 2.0;
        get_intersection(
            the_c1, the_t1f, a_t1m, the_c2, the_t2f, the_t2l, the_tol_conf, the_max_count, the_p_int,
            the_dist, the_count,
        );
        get_intersection(
            the_c1, a_t1m, the_t1l, the_c2, the_t2f, the_t2l, the_tol_conf, the_max_count, the_p_int,
            the_dist, the_count,
        );
    } else {
        let a_t1m = (the_t1l + the_t1f) / 2.0;
        let a_t2m = (the_t2l + the_t2f) / 2.0;
        get_intersection(
            the_c1, the_t1f, a_t1m, the_c2, the_t2f, a_t2m, the_tol_conf, the_max_count, the_p_int,
            the_dist, the_count,
        );
        get_intersection(
            the_c1, the_t1f, a_t1m, the_c2, a_t2m, the_t2l, the_tol_conf, the_max_count, the_p_int,
            the_dist, the_count,
        );
        get_intersection(
            the_c1, a_t1m, the_t1l, the_c2, the_t2f, a_t2m, the_tol_conf, the_max_count, the_p_int,
            the_dist, the_count,
        );
        get_intersection(
            the_c1, a_t1m, the_t1l, the_c2, a_t2m, the_t2l, the_tol_conf, the_max_count, the_p_int,
            the_dist, the_count,
        );
    }
}
