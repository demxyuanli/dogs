//! Position bits and end-point mapping helpers.

use super::prelude::*;
use super::super::curve_tool;

/// The `PosSegment |= ...` table shared by both `Perform` wrappers
/// (`gxx:145-157` / `:176-188` / `:201-213`).
#[inline]
pub(super) fn position_bits(pos1: IntRes2dPosition, pos2: IntRes2dPosition) -> i32 {
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
pub(super) fn head_or_end_point(
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

    let (mut p1, mut t1) = curve_tool::geom2d_d1(c1, u);
    let (mut p2, mut t2) = curve_tool::geom2d_d1(c2, v);

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

    let eps_x1 = curve_tool::geom2d_eps_x();
    let eps_x2 = curve_tool::geom2d_eps_x();

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
            let (_, tt1, n1) = curve_tool::geom2d_d2(c1, svu);
            let (_, tt2, n2) = curve_tool::geom2d_d2(c2, svv);
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
