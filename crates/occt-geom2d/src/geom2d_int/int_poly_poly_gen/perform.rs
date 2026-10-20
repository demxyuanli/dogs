//! `Perform` entry points (`gxx:94-387`).

use super::prelude::*;
use super::geometry::{head_or_end_point, position_bits};
use super::*;

impl Geom2dIntIntPolyPolyGen {
    /// `Perform(C1, D1, C2, D2, TolConf, Tol)` (`gxx:94-292`).
    pub fn perform(
        &mut self,
        c1: &dyn Curve2d,
        d1: &IntRes2dDomain,
        c2: &dyn Curve2d,
        d2: &IntRes2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        self.reset_fields();
        self.domain_on_curve1 = *d1;
        self.domain_on_curve2 = *d2;
        let du = d1.last_parameter() - d1.first_parameter();
        let dv = d2.last_parameter() - d2.first_parameter();
        let tl = if tol < TOL_MINI { TOL_MINI } else { tol };
        let tl_conf = if tol_conf < TOL_CONF_MINI { TOL_CONF_MINI } else { tol_conf };
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
                d1, c1, u0, d2, c2, v0, tol_conf, &mut int_pt, &mut head_on1, &mut head_on2,
                &mut end_on1, &mut end_on2, pos_segment,
            ) {
                self.base.insert(&int_pt);
            }
        }
        if d1.first_tolerance() != 0.0 || d2.last_tolerance() != 0.0 {
            if head_or_end_point(
                d1, c1, u0, d2, c2, v1, tol_conf, &mut int_pt, &mut head_on1, &mut head_on2,
                &mut end_on1, &mut end_on2, pos_segment,
            ) {
                self.base.insert(&int_pt);
            }
        }
        if d1.last_tolerance() != 0.0 || d2.first_tolerance() != 0.0 {
            if head_or_end_point(
                d1, c1, u1, d2, c2, v0, tol_conf, &mut int_pt, &mut head_on1, &mut head_on2,
                &mut end_on1, &mut end_on2, pos_segment,
            ) {
                self.base.insert(&int_pt);
            }
        }
        if d1.last_tolerance() != 0.0 || d2.last_tolerance() != 0.0 {
            if head_or_end_point(
                d1, c1, u1, d2, c2, v1, tol_conf, &mut int_pt, &mut head_on1, &mut head_on2,
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
        tol_conf: f64,
        tol: f64,
    ) {
        self.reset_fields();
        self.domain_on_curve1 = *d1;
        self.domain_on_curve2 = *d1;
        let du = d1.last_parameter() - d1.first_parameter();
        let tl = if tol < TOL_MINI { TOL_MINI } else { tol };
        let tl_conf = if tol_conf < TOL_CONF_MINI { TOL_CONF_MINI } else { tol_conf };
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

}
