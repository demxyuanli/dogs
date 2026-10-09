//! Port of `IntCurve_IntConicConic::Perform(const gp_Lin2d&, ...)`
//! (`IntCurve_IntConicConic_1.cxx:1381-2232`), the Line/Line kernel, on top of
//! the `occt-core` `intres2d` result types.

use occt_core::elib::clib2d;
use occt_core::gp::{GpLin2d, GpVec2d};
use occt_core::precision::PCONFUSION;
use occt_core::intres2d::{
    IntRes2dDomain, IntRes2dIntersection, IntRes2dIntersectionPoint,
    IntRes2dIntersectionSegment, IntRes2dPosition, IntRes2dSituation,
    IntRes2dTransition, IntRes2dTypeTrans,
};

use super::conic_conic::*;

/// `IntRes2d_Transition::SetValue(false, Pos, In/Out)`.
fn set_inout(pos: IntRes2dPosition, ty: IntRes2dTypeTrans) -> IntRes2dTransition {
    IntRes2dTransition::in_out(false, pos, ty)
}

/// `IntRes2d_Transition::SetValue(false, Pos, IntRes2d_Unknown, Oppos)`.
fn set_touch(pos: IntRes2dPosition, oppos: bool) -> IntRes2dTransition {
    IntRes2dTransition::touch(false, pos, IntRes2dSituation::Unknown, oppos)
}

/// `SegmentToPoint` (`_1.cxx:2620-2653`).
fn segment_to_point(
    pa: &IntRes2dIntersectionPoint,
    t1a: &IntRes2dTransition,
    t2a: &IntRes2dTransition,
    pb: &IntRes2dIntersectionPoint,
    t1b: &IntRes2dTransition,
    t2b: &IntRes2dTransition,
) -> IntRes2dIntersectionPoint {
    if t1b.position_on_curve() == IntRes2dPosition::Middle
        && t2b.position_on_curve() == IntRes2dPosition::Middle
    {
        return *pa;
    }
    if t1a.position_on_curve() == IntRes2dPosition::Middle
        && t2a.position_on_curve() == IntRes2dPosition::Middle
    {
        return *pb;
    }
    let mut t1 = *t1a;
    let mut t2 = *t2a;
    let mut u1 = pa.param_on_first();
    let mut u2 = pa.param_on_second();
    if t1.position_on_curve() == IntRes2dPosition::Middle {
        t1.set_position(t1b.position_on_curve());
        u1 = pb.param_on_first();
    }
    if t2.position_on_curve() == IntRes2dPosition::Middle {
        t2.set_position(t2b.position_on_curve());
        u2 = pb.param_on_second();
    }
    IntRes2dIntersectionPoint::with_transitions(pa.value(), u1, u2, &t1, &t2, false)
}

/// `IntCurve_IntConicConic` (`IntCurve_IntConicConic.hxx:39-325`): carries the
/// `IntRes2d_Intersection` result of the last `Perform`.
#[derive(Clone, Debug)]
pub struct IntCurveIntConicConic {
    pub(crate) result: IntRes2dIntersection,
}

impl Default for IntCurveIntConicConic {
    fn default() -> Self {
        Self { result: IntRes2dIntersection::new() }
    }
}

impl IntCurveIntConicConic {
    pub fn new() -> Self {
        Self::default()
    }

    /// `Append(const IntRes2d_IntersectionPoint&)` (`lxx:96-99`).
    fn append_pt(&mut self, p: IntRes2dIntersectionPoint) {
        self.result.append_point(&p);
    }

    /// `Append(const IntRes2d_IntersectionSegment&)` (`lxx:90-93`).
    fn append_seg(&mut self, s: IntRes2dIntersectionSegment) {
        self.result.append_segment(&s);
    }

    /// `IntRes2d_Intersection::SetReversedParameters` (inherited).
    pub fn set_reversed_parameters(&mut self, flag: bool) {
        self.result.set_reversed_parameters(flag);
    }

    pub fn is_done(&self) -> bool {
        self.result.is_done()
    }

    pub fn result(&self) -> &IntRes2dIntersection {
        &self.result
    }

    /// `Perform(const gp_Lin2d&, D1, const gp_Lin2d&, D2, TolConf, Tol)`
    /// (`IntCurve_IntConicConic_1.cxx:1381-2232`).
    pub fn perform_line_line(
        &mut self,
        l1: &GpLin2d,
        domain1: &IntRes2dDomain,
        l2: &GpLin2d,
        domain2: &IntRes2dDomain,
        _tol_conf: f64,
        tol_r: f64,
    ) {
        self.result.reset_fields();
        let mut tol = tol_r;
        if tol < PCONFUSION {
            tol = PCONFUSION;
        }
        let (mut u1, mut u2, a_half_sin, nbsol0) = line_line_geometric_intersection(l1, l2, tol);
        let mut nbsol = nbsol0;
        let mut pt_seg1 = IntRes2dIntersectionPoint::new();
        let mut pt_seg2 = IntRes2dIntersectionPoint::new();
        let tan1 = GpVec2d::new(l1.position().vdir.x, l1.position().vdir.y);
        let tan2 = GpVec2d::new(l2.position().vdir.x, l2.position().vdir.y);
        let a_cos_t1t2 = tan1.dot(&tan2);
        let is_opposite = a_cos_t1t2 < 0.0;
        self.result.done = true;
        if nbsol == 1 && check_ll_coincidence(l1, l2, domain1, domain2, tol) {
            nbsol = 2;
        }

        if nbsol == 1 {
            let d = 0.5 * tol / a_half_sin;
            let u1m_u2 = u1 - u2;
            let u1p_u2 = u1 + u2;
            let mut u1inf = u1 - d;
            let mut u1sup = u1 + d;
            if domain1.has_first_point() && l2.distance(domain1.first_point()) < domain1.first_tolerance() {
                if u1inf > domain1.first_parameter() { u1inf = domain1.first_parameter(); }
                if u1sup < domain1.first_parameter() { u1sup = domain1.first_parameter(); }
            }
            if domain1.has_last_point() && l2.distance(domain1.last_point()) < domain1.last_tolerance() {
                if u1inf > domain1.last_parameter() { u1inf = domain1.last_parameter(); }
                if u1sup < domain1.last_parameter() { u1sup = domain1.last_parameter(); }
            }
            if domain2.has_first_point() && l1.distance(domain2.first_point()) < domain2.first_tolerance() {
                let p = clib2d::line_parameter_ax2d(l1.position(), domain2.first_point());
                if u1inf > p { u1inf = p; }
                if u1sup < p { u1sup = p; }
            }
            if domain2.has_last_point() && l1.distance(domain2.last_point()) < domain2.last_tolerance() {
                let p = clib2d::line_parameter_ax2d(l1.position(), domain2.last_point());
                if u1inf > p { u1inf = p; }
                if u1sup < p { u1sup = p; }
            }
            let (mut res1inf, mut res1sup, p1a0, p1b0) = domain_intersection(domain1, u1inf, u1sup);
            let mut pos1a = p1a0;
            let mut pos1b = p1b0;
            if (res1sup - res1inf) < 0.0 {
                // empty
            } else {
                let prod_vect_tan = tan1.x() * tan2.y() - tan1.y() * tan2.x();
                let long_mini_seg = tol;
                if ((res1sup - res1inf) <= long_mini_seg)
                    || ((pos1a == pos1b) && (pos1a != IntRes2dPosition::Middle))
                {
                    let a_cur_trans = if prod_vect_tan >= TOLERANCE_ANGULAIRE {
                        IntRes2dTypeTrans::Out
                    } else if prod_vect_tan <= -TOLERANCE_ANGULAIRE {
                        IntRes2dTypeTrans::In
                    } else {
                        IntRes2dTypeTrans::Undecided
                    };
                    let mut res_sup = res1sup;
                    if let Some(np) = compute_int_point(
                        domain1, domain2, l1, l2, a_cos_t1t2, u1, u2, res1inf, &mut res_sup, 1, a_cur_trans,
                    ) {
                        self.append_pt(np);
                    }
                } else {
                    let (u2inf, u2sup) = if is_opposite {
                        (u1p_u2 - res1sup, u1p_u2 - res1inf)
                    } else {
                        (res1inf - u1m_u2, res1sup - u1m_u2)
                    };
                    let (r2inf, r2sup, p2a0, p2b0) = domain_intersection(domain2, u2inf, u2sup);
                    let mut res2inf = r2inf;
                    let mut res2sup = r2sup;
                    let mut pos2a = p2a0;
                    let mut pos2b = p2b0;
                    let res2sup_m_res2inf = res2sup - res2inf;
                    if res2sup_m_res2inf < 0.0 {
                        // no solution
                    } else if (res2sup_m_res2inf > long_mini_seg)
                        || ((pos2a == pos2b) && (pos2a != IntRes2dPosition::Middle))
                    {
                        if is_opposite {
                            res1inf = u1p_u2 - res2sup;
                            res1sup = u1p_u2 - res2inf;
                            std::mem::swap(&mut res2inf, &mut res2sup);
                            std::mem::swap(&mut pos2a, &mut pos2b);
                        } else {
                            res1inf = u1m_u2 + res2inf;
                            res1sup = u1m_u2 + res2sup;
                        }
                        pos1a = find_position_ll(&mut res1inf, domain1);
                        pos1b = find_position_ll(&mut res1sup, domain1);
                        let t1a;
                        let t2a;
                        if prod_vect_tan >= TOLERANCE_ANGULAIRE {
                            t1a = set_inout(pos1a, IntRes2dTypeTrans::Out);
                            t2a = set_inout(pos2a, IntRes2dTypeTrans::In);
                        } else if prod_vect_tan <= -TOLERANCE_ANGULAIRE {
                            t1a = set_inout(pos1a, IntRes2dTypeTrans::In);
                            t2a = set_inout(pos2a, IntRes2dTypeTrans::Out);
                        } else {
                            t1a = set_touch(pos1a, is_opposite);
                            t2a = set_touch(pos2a, is_opposite);
                        }
                        let mut result_is_a_point = false;
                        if ((res1sup - res1inf) <= long_mini_seg)
                            || ((res2sup - res2inf).abs() <= long_mini_seg)
                        {
                            result_is_a_point = true;
                        } else {
                            if pos1a == IntRes2dPosition::Head && pos1b != IntRes2dPosition::End && u1 < res1inf {
                                result_is_a_point = true;
                                u1 = res1inf;
                                u2 = res2inf;
                            }
                            if pos1b == IntRes2dPosition::End && pos1a != IntRes2dPosition::Head && u1 > res1sup {
                                result_is_a_point = true;
                                u1 = res1sup;
                                u2 = res2sup;
                            }
                            if pos2a == IntRes2dPosition::Head {
                                if pos2b != IntRes2dPosition::End && u2 < res2inf {
                                    result_is_a_point = true;
                                    u2 = res2inf;
                                    u1 = res1inf;
                                }
                            } else if pos2a == IntRes2dPosition::End && pos2b != IntRes2dPosition::Head && u2 > res2inf {
                                result_is_a_point = true;
                                u2 = res2inf;
                                u1 = res1inf;
                            }
                            if pos2b == IntRes2dPosition::Head {
                                if pos2a != IntRes2dPosition::End && u2 < res2sup {
                                    result_is_a_point = true;
                                    u2 = res2sup;
                                    u1 = res1sup;
                                }
                            } else if pos2b == IntRes2dPosition::End && pos2a != IntRes2dPosition::Head && u2 > res2sup {
                                result_is_a_point = true;
                                u2 = res2sup;
                                u1 = res1sup;
                            }
                        }
                        if !result_is_a_point
                            && (pos1a != IntRes2dPosition::Middle || pos2a != IntRes2dPosition::Middle)
                        {
                            let (t1b, t2b) = if prod_vect_tan >= TOLERANCE_ANGULAIRE {
                                (set_inout(pos1b, IntRes2dTypeTrans::Out), set_inout(pos2b, IntRes2dTypeTrans::In))
                            } else if prod_vect_tan <= -TOLERANCE_ANGULAIRE {
                                (set_inout(pos1b, IntRes2dTypeTrans::In), set_inout(pos2b, IntRes2dTypeTrans::Out))
                            } else {
                                (set_touch(pos1b, is_opposite), set_touch(pos2b, is_opposite))
                            };
                            let pt_debut;
                            if pos1a == IntRes2dPosition::Middle {
                                let t3 = if is_opposite {
                                    if pos2a == IntRes2dPosition::Head { res2sup } else { res2inf }
                                } else if pos2a == IntRes2dPosition::Head {
                                    res2inf
                                } else {
                                    res2sup
                                };
                                pt_debut = clib2d::line_value_ax2d(t3, l2.position());
                                res1inf = clib2d::line_parameter_ax2d(l1.position(), &pt_debut);
                            } else {
                                let t4 = if pos1a == IntRes2dPosition::Head { res1inf } else { res1sup };
                                pt_debut = clib2d::line_value_ax2d(t4, l1.position());
                                res2inf = clib2d::line_parameter_ax2d(l2.position(), &pt_debut);
                            }
                            pt_seg1.set_values(&pt_debut, res1inf, res2inf, &t1a, &t2a, false);
                            if pos1b != IntRes2dPosition::Middle || pos2b != IntRes2dPosition::Middle {
                                let pt_fin;
                                if pos1b == IntRes2dPosition::Middle {
                                    pt_fin = clib2d::line_value_ax2d(res2sup, l2.position());
                                    res1sup = clib2d::line_parameter_ax2d(l1.position(), &pt_fin);
                                } else {
                                    pt_fin = clib2d::line_value_ax2d(res1sup, l1.position());
                                    res2sup = clib2d::line_parameter_ax2d(l2.position(), &pt_fin);
                                }
                                pt_seg2.set_values(&pt_fin, res1sup, res2sup, &t1b, &t2b, false);
                                self.append_seg(IntRes2dIntersectionSegment::from_two_points(
                                    &pt_seg1, &pt_seg2, is_opposite, false,
                                ));
                            } else {
                                let pos1b2 = find_position_ll(&mut u1, domain1);
                                let pos2b2 = find_position_ll(&mut u2, domain2);
                                let (t1b2, t2b2) = if prod_vect_tan >= TOLERANCE_ANGULAIRE {
                                    (set_inout(pos1b2, IntRes2dTypeTrans::Out), set_inout(pos2b2, IntRes2dTypeTrans::In))
                                } else if prod_vect_tan <= -TOLERANCE_ANGULAIRE {
                                    (set_inout(pos1b2, IntRes2dTypeTrans::In), set_inout(pos2b2, IntRes2dTypeTrans::Out))
                                } else {
                                    (set_touch(pos1b2, is_opposite), set_touch(pos2b2, is_opposite))
                                };
                                pt_seg2.set_values(
                                    &clib2d::line_value_ax2d(u2, l2.position()), u1, u2, &t1b2, &t2b2, false,
                                );
                                if (res1inf - u1).abs() > long_mini_seg && (res2inf - u2).abs() > long_mini_seg {
                                    self.append_seg(IntRes2dIntersectionSegment::from_two_points(
                                        &pt_seg1, &pt_seg2, is_opposite, false,
                                    ));
                                } else {
                                    self.append_pt(segment_to_point(
                                        &pt_seg1, &t1a, &t2a, &pt_seg2, &t1b2, &t2b2,
                                    ));
                                }
                            }
                        } else {
                            if pos1b == IntRes2dPosition::Middle { pos1b = pos1a; }
                            if pos2b == IntRes2dPosition::Middle { pos2b = pos2a; }
                            if result_is_a_point {
                                if pos1b != IntRes2dPosition::Middle || pos2b != IntRes2dPosition::Middle {
                                    let (t1b, t2b);
                                    let pt_fin;
                                    if pos1b == IntRes2dPosition::Middle {
                                        let t2 = if is_opposite {
                                            if pos2b == IntRes2dPosition::Head { res2sup } else { res2inf }
                                        } else if pos2b == IntRes2dPosition::Head {
                                            res2inf
                                        } else {
                                            res2sup
                                        };
                                        pt_fin = clib2d::line_value_ax2d(t2, l2.position());
                                        res1sup = clib2d::line_parameter_ax2d(l1.position(), &pt_fin);
                                        pos1b = find_position_ll(&mut res1sup, domain1);
                                    } else {
                                        let t1 = if pos1b == IntRes2dPosition::Head { res1inf } else { res1sup };
                                        pt_fin = clib2d::line_value_ax2d(t1, l1.position());
                                        res2sup = clib2d::line_parameter_ax2d(l2.position(), &pt_fin);
                                        pos2b = find_position_ll(&mut res2sup, domain2);
                                    }
                                    if prod_vect_tan >= TOLERANCE_ANGULAIRE {
                                        t1b = set_inout(pos1b, IntRes2dTypeTrans::Out);
                                        t2b = set_inout(pos2b, IntRes2dTypeTrans::In);
                                    } else if prod_vect_tan <= -TOLERANCE_ANGULAIRE {
                                        t1b = set_inout(pos1b, IntRes2dTypeTrans::In);
                                        t2b = set_inout(pos2b, IntRes2dTypeTrans::Out);
                                    } else {
                                        t1b = set_touch(pos1b, is_opposite);
                                        t2b = set_touch(pos2b, is_opposite);
                                    }
                                    pt_seg2.set_values(&pt_fin, res1sup, res2sup, &t1b, &t2b, false);
                                    self.append_pt(pt_seg2);
                                } else {
                                    let pos1b2 = find_position_ll(&mut u1, domain1);
                                    let pos2b2 = find_position_ll(&mut u2, domain2);
                                    let (t1b2, t2b2) = if prod_vect_tan >= TOLERANCE_ANGULAIRE {
                                        (set_inout(pos1b2, IntRes2dTypeTrans::Out), set_inout(pos2b2, IntRes2dTypeTrans::In))
                                    } else if prod_vect_tan <= -TOLERANCE_ANGULAIRE {
                                        (set_inout(pos1b2, IntRes2dTypeTrans::In), set_inout(pos2b2, IntRes2dTypeTrans::Out))
                                    } else {
                                        (set_touch(pos1b2, is_opposite), set_touch(pos2b2, is_opposite))
                                    };
                                    pt_seg1.set_values(
                                        &clib2d::line_value_ax2d(u2, l2.position()), u1, u2, &t1b2, &t2b2, false,
                                    );
                                    self.append_pt(pt_seg1);
                                }
                            } else {
                                pt_seg1.set_values(
                                    &clib2d::line_value_ax2d(u2, l2.position()), u1, u2, &t1a, &t2a, false,
                                );
                                if pos1b != IntRes2dPosition::Middle || pos2b != IntRes2dPosition::Middle {
                                    let (t1b, t2b) = if prod_vect_tan >= TOLERANCE_ANGULAIRE {
                                        (set_inout(pos1b, IntRes2dTypeTrans::Out), set_inout(pos2b, IntRes2dTypeTrans::In))
                                    } else if prod_vect_tan <= -TOLERANCE_ANGULAIRE {
                                        (set_inout(pos1b, IntRes2dTypeTrans::In), set_inout(pos2b, IntRes2dTypeTrans::Out))
                                    } else {
                                        (set_touch(pos1b, is_opposite), set_touch(pos2b, is_opposite))
                                    };
                                    let pt_fin;
                                    if pos1b == IntRes2dPosition::Middle {
                                        pt_fin = clib2d::line_value_ax2d(res2sup, l2.position());
                                        res1sup = clib2d::line_parameter_ax2d(l1.position(), &pt_fin);
                                    } else {
                                        pt_fin = clib2d::line_value_ax2d(res1sup, l1.position());
                                        res2sup = clib2d::line_parameter_ax2d(l2.position(), &pt_fin);
                                    }
                                    pt_seg2.set_values(&pt_fin, res1sup, res2sup, &t1b, &t2b, false);
                                    if (u1 - res1sup).abs() > long_mini_seg || (u2 - res2sup).abs() > long_mini_seg {
                                        self.append_seg(IntRes2dIntersectionSegment::from_two_points(
                                            &pt_seg1, &pt_seg2, is_opposite, false,
                                        ));
                                    } else {
                                        self.append_pt(segment_to_point(
                                            &pt_seg1, &t1a, &t2a, &pt_seg2, &t1b, &t2b,
                                        ));
                                    }
                                } else {
                                    self.append_pt(pt_seg1);
                                }
                            }
                        }
                    } else {
                        let a_cur_trans = if prod_vect_tan >= TOLERANCE_ANGULAIRE {
                            IntRes2dTypeTrans::In
                        } else if prod_vect_tan <= -TOLERANCE_ANGULAIRE {
                            IntRes2dTypeTrans::Out
                        } else {
                            IntRes2dTypeTrans::Undecided
                        };
                        let mut res2sup_mut = res2sup;
                        if let Some(np) = compute_int_point(
                            domain2, domain1, l2, l1, a_cos_t1t2, u2, u1, res2inf, &mut res2sup_mut, 2,
                            a_cur_trans,
                        ) {
                            self.append_pt(np);
                        }
                    }
                }
            }
        } else if nbsol == 2 {
            let org2_sur_l1 = clib2d::line_parameter_ax2d(l1.position(), &l2.location());
            let mut res_has_first_point = 0;
            let mut res_has_last_point = 0;
            if domain1.has_first_point() { res_has_first_point = 1; }
            if domain1.has_last_point() { res_has_last_point = 1; }
            if is_opposite {
                if domain2.has_last_point() { res_has_first_point += 2; }
                if domain2.has_first_point() { res_has_last_point += 2; }
            } else {
                if domain2.has_last_point() { res_has_last_point += 2; }
                if domain2.has_first_point() { res_has_first_point += 2; }
            }
            if res_has_first_point == 0 && res_has_last_point == 0 {
                self.append_seg(IntRes2dIntersectionSegment::infinite(is_opposite));
            } else {
                let mut param_start = 0.0f64;
                let mut param_start2 = 0.0f64;
                let mut param_end = 0.0f64;
                let mut param_end2 = 0.0f64;
                match res_has_first_point {
                    1 => {
                        param_start = domain1.first_parameter();
                        param_start2 = if is_opposite { org2_sur_l1 - param_start } else { param_start - org2_sur_l1 };
                    }
                    2 => {
                        if is_opposite {
                            param_start2 = domain2.last_parameter();
                            param_start = org2_sur_l1 - param_start2;
                        } else {
                            param_start2 = domain2.first_parameter();
                            param_start = org2_sur_l1 + param_start2;
                        }
                    }
                    3 => {
                        if is_opposite {
                            param_start2 = domain2.last_parameter();
                            param_start = org2_sur_l1 - param_start2;
                            if param_start < domain1.first_parameter() {
                                param_start = domain1.first_parameter();
                                param_start2 = org2_sur_l1 - param_start;
                            }
                        } else {
                            param_start2 = domain2.first_parameter();
                            param_start = org2_sur_l1 + param_start2;
                            if param_start < domain1.first_parameter() {
                                param_start = domain1.first_parameter();
                                param_start2 = param_start - org2_sur_l1;
                            }
                        }
                    }
                    _ => {}
                }
                match res_has_last_point {
                    1 => {
                        param_end = domain1.last_parameter();
                        param_end2 = if is_opposite { org2_sur_l1 - param_end } else { param_end - org2_sur_l1 };
                    }
                    2 => {
                        if is_opposite {
                            param_end2 = domain2.first_parameter();
                            param_end = org2_sur_l1 - param_end2;
                        } else {
                            param_end2 = domain2.last_parameter();
                            param_end = org2_sur_l1 + param_end2;
                        }
                    }
                    3 => {
                        if is_opposite {
                            param_end2 = domain2.first_parameter();
                            param_end = org2_sur_l1 - param_end2;
                            if param_end > domain1.last_parameter() {
                                param_end = domain1.last_parameter();
                                param_end2 = org2_sur_l1 - param_end;
                            }
                        } else {
                            param_end2 = domain2.last_parameter();
                            param_end = org2_sur_l1 + param_end2;
                            if param_end > domain1.last_parameter() {
                                param_end = domain1.last_parameter();
                                param_end2 = param_end - org2_sur_l1;
                            }
                        }
                    }
                    _ => {}
                }
                if res_has_first_point != 0 {
                    if res_has_last_point != 0 {
                        if param_end >= (param_start - tol) {
                            let mut ps = param_start;
                            let mut ps2 = param_start2;
                            let pos1 = find_position_ll(&mut ps, domain1);
                            let pos2 = find_position_ll(&mut ps2, domain2);
                            let tinf = IntRes2dTransition::touch(true, pos1, IntRes2dSituation::Unknown, is_opposite);
                            let tsup = IntRes2dTransition::touch(true, pos2, IntRes2dSituation::Unknown, is_opposite);
                            let p1 = IntRes2dIntersectionPoint::with_transitions(
                                &clib2d::line_value_ax2d(ps, l1.position()), ps, ps2, &tinf, &tsup, false,
                            );
                            if param_end > (param_start + tol) {
                                let mut pe = param_end;
                                let mut pe2 = param_end2;
                                let pos1b = find_position_ll(&mut pe, domain1);
                                let pos2b = find_position_ll(&mut pe2, domain2);
                                let tinf2 = IntRes2dTransition::touch(true, pos1b, IntRes2dSituation::Unknown, is_opposite);
                                let tsup2 = IntRes2dTransition::touch(true, pos2b, IntRes2dSituation::Unknown, is_opposite);
                                let p2 = IntRes2dIntersectionPoint::with_transitions(
                                    &clib2d::line_value_ax2d(pe, l1.position()), pe, pe2, &tinf2, &tsup2, false,
                                );
                                self.append_seg(IntRes2dIntersectionSegment::from_two_points(
                                    &p1, &p2, is_opposite, false,
                                ));
                            } else {
                                self.append_pt(p1);
                            }
                        }
                    } else {
                        let mut ps = param_start;
                        let mut ps2 = param_start2;
                        let pos1 = find_position_ll(&mut ps, domain1);
                        let pos2 = find_position_ll(&mut ps2, domain2);
                        let tinf = IntRes2dTransition::touch(true, pos1, IntRes2dSituation::Unknown, is_opposite);
                        let tsup = IntRes2dTransition::touch(true, pos2, IntRes2dSituation::Unknown, is_opposite);
                        let p = IntRes2dIntersectionPoint::with_transitions(
                            &clib2d::line_value_ax2d(ps, l1.position()), ps, ps2, &tinf, &tsup, false,
                        );
                        self.append_seg(IntRes2dIntersectionSegment::from_one_point(
                            &p, true, is_opposite, false,
                        ));
                    }
                } else {
                    let mut pe = param_end;
                    let mut pe2 = param_end2;
                    let pos1 = find_position_ll(&mut pe, domain1);
                    let pos2 = find_position_ll(&mut pe2, domain2);
                    let tinf = IntRes2dTransition::touch(true, pos1, IntRes2dSituation::Unknown, is_opposite);
                    let tsup = IntRes2dTransition::touch(true, pos2, IntRes2dSituation::Unknown, is_opposite);
                    let p2 = IntRes2dIntersectionPoint::with_transitions(
                        &clib2d::line_value_ax2d(pe, l1.position()), pe, pe2, &tinf, &tsup, false,
                    );
                    self.append_seg(IntRes2dIntersectionSegment::from_one_point(
                        &p2, false, is_opposite, false,
                    ));
                }
            }
        }
    }
}
