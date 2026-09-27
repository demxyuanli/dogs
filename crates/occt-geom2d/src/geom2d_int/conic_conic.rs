//! Port of the file-static helpers shared by `IntCurve_IntConicConic::Perform`
//! (`IntCurve_IntConicConic_1.cxx`), needed by the Line/Line branch. The
//! `Perform` overloads themselves are ported on top of these.

use occt_core::elib::clib2d;
use occt_core::gp::{GpLin2d, GpPnt2d, GpVec2d};
use occt_core::precision::PCONFUSION;

use occt_core::intres2d::{IntRes2dDomain, IntRes2dIntersectionPoint, IntRes2dPosition, IntRes2dSituation, IntRes2dTransition, IntRes2dTypeTrans};

/// `TOLERANCE_ANGULAIRE` (`IntCurve_IntConicConic_Tool.cxx:20`).
pub(crate) const TOLERANCE_ANGULAIRE: f64 = 0.00000001;

/// `LineLineGeometricIntersection` (`_1.cxx:730-774`). Returns
/// `(U1, U2, SinDemiAngle, nbsol)`; `U1`/`U2` are meaningless when `nbsol != 1`.
pub(crate) fn line_line_geometric_intersection(
    l1: &GpLin2d,
    l2: &GpLin2d,
    tol: f64,
) -> (f64, f64, f64, i32) {
    let u1x = l1.position().vdir.x;
    let u1y = l1.position().vdir.y;
    let u2x = l2.position().vdir.x;
    let u2y = l2.position().vdir.y;
    let uo21x = l2.position().loc.x() - l1.position().loc.x();
    let uo21y = l2.position().loc.y() - l1.position().loc.y();
    let mut d = u1y * u2x - u1x * u2y;
    if d.abs() < TOLERANCE_ANGULAIRE {
        d = u1y * uo21x - u1x * uo21y;
        let nbsol = if d.abs() <= tol { 2 } else { 0 };
        (0.0, 0.0, 0.0, nbsol)
    } else {
        let u1 = (uo21y * u2x - uo21x * u2y) / d;
        let u2 = (uo21y * u1x - uo21x * u1y) / d;
        if d < 0.0 {
            d = -d;
        }
        if d > 1.0 {
            d = 1.0;
        }
        let sin_demi_angle = (0.5 * d.asin()).sin();
        (u1, u2, sin_demi_angle, 1)
    }
}

/// `DomainIntersection` (`_1.cxx:641-726`). Returns
/// `(Res1inf, Res1sup, PosInf, PosSup)`; the empty marker is `(1, -1, ...)`.
pub(crate) fn domain_intersection(
    domain: &IntRes2dDomain,
    u1inf: f64,
    u1sup: f64,
) -> (f64, f64, IntRes2dPosition, IntRes2dPosition) {
    let mut pos_inf = IntRes2dPosition::Middle;
    let mut pos_sup = IntRes2dPosition::Middle;
    let res1inf;
    let res1sup;
    if domain.has_first_point() {
        if u1sup < (domain.first_parameter() - domain.first_tolerance()) {
            return (1.0, -1.0, pos_inf, pos_sup);
        }
        if u1inf > (domain.first_parameter() + domain.first_tolerance()) {
            res1inf = u1inf;
            pos_inf = IntRes2dPosition::Middle;
        } else {
            res1inf = domain.first_parameter();
            pos_inf = IntRes2dPosition::Head;
        }
    } else {
        res1inf = u1inf;
        pos_inf = IntRes2dPosition::Middle;
    }
    if domain.has_last_point() {
        if u1inf > (domain.last_parameter() + domain.last_tolerance()) {
            return (1.0, -1.0, pos_inf, pos_sup);
        }
        if u1sup < (domain.last_parameter() - domain.last_tolerance()) {
            res1sup = u1sup;
            pos_sup = IntRes2dPosition::Middle;
        } else {
            res1sup = domain.last_parameter();
            pos_sup = IntRes2dPosition::End;
        }
    } else {
        res1sup = u1sup;
        pos_sup = IntRes2dPosition::Middle;
    }
    let mut res1inf = res1inf;
    let res1sup = res1sup;
    // `_1.cxx:701-711`: keep the parameters ordered.
    if res1inf > res1sup {
        if pos_sup == IntRes2dPosition::Middle {
            return (res1inf, res1inf, pos_inf, pos_sup);
        } else {
            res1inf = res1sup;
        }
    }
    (res1inf, res1sup, pos_inf, pos_sup)
}

/// `FindPositionLL` (`_1.cxx:1209-1234`): snap `param` onto the domain bounds
/// within their tolerance and report the position.
pub(crate) fn find_position_ll(param: &mut f64, domain: &IntRes2dDomain) -> IntRes2dPosition {
    let mut d_par = f64::INFINITY; // Precision::Infinite() approximates
    let mut pos = IntRes2dPosition::Middle;
    let mut res_par = *param;
    if domain.has_first_point() {
        d_par = (*param - domain.first_parameter()).abs();
        if d_par <= domain.first_tolerance() {
            res_par = domain.first_parameter();
            pos = IntRes2dPosition::Head;
        }
    }
    if domain.has_last_point() {
        let d2 = (*param - domain.last_parameter()).abs();
        if d2 <= domain.last_tolerance() && (pos == IntRes2dPosition::Middle || d2 < d_par) {
            res_par = domain.last_parameter();
            pos = IntRes2dPosition::End;
        }
    }
    *param = res_par;
    pos
}

/// `getDomainParametrs` (`_1.cxx:1240-1250`).
pub(crate) fn get_domain_parameters(domain: &IntRes2dDomain) -> (f64, f64, f64, f64) {
    let first = if domain.has_first_point() { domain.first_parameter() } else { -f64::INFINITY };
    let last = if domain.has_last_point() { domain.last_parameter() } else { f64::INFINITY };
    let tol1 = if domain.has_first_point() { domain.first_tolerance() } else { 0.0 };
    let tol2 = if domain.has_last_point() { domain.last_tolerance() } else { 0.0 };
    (first, last, tol1, tol2)
}

/// `CheckLLCoincidence` (`_1.cxx:1363-1378`): true when the two trimmed lines
/// coincide within `tol`.
pub(crate) fn check_ll_coincidence(
    l1: &GpLin2d,
    l2: &GpLin2d,
    domain1: &IntRes2dDomain,
    domain2: &IntRes2dDomain,
    tol: f64,
) -> bool {
    let is_first1 = domain1.has_first_point() && l2.distance(domain1.first_point()) < tol;
    let is_last1 = domain1.has_last_point() && l2.distance(domain1.last_point()) < tol;
    if is_first1 && is_last1 {
        return true;
    }
    let is_first2 = domain2.has_first_point() && l1.distance(domain2.first_point()) < tol;
    let is_last2 = domain2.has_last_point() && l1.distance(domain2.last_point()) < tol;
    is_first2 && is_last2
}

/// `computeIntPoint` (`_1.cxx:1254-1356`). Returns the intersection point, or
/// `None` when it lies outside both domains.
#[allow(clippy::too_many_arguments)]
pub(crate) fn compute_int_point(
    cur_domain: &IntRes2dDomain,
    other_domain: &IntRes2dDomain,
    cur_lin: &GpLin2d,
    other_lin: &GpLin2d,
    cos_t1t2: f64,
    par_cur: f64,
    par_other: f64,
    res_inf: f64,
    res_sup: &mut f64,
    num: i32,
    cur_trans: IntRes2dTypeTrans,
) -> Option<IntRes2dIntersectionPoint> {
    // `_1.cxx:1267-1270`.
    if (*res_sup - par_cur).abs() > (res_inf - par_cur).abs() {
        *res_sup = res_inf;
    }
    let mut a_res2 = par_other + (*res_sup - par_cur) * cos_t1t2;
    let (a_first2, a_last2, a_tol21, a_tol22) = get_domain_parameters(other_domain);
    if a_res2 < a_first2 - a_tol21 || a_res2 > a_last2 + a_tol22 {
        return None;
    }
    let (_, _, a_tol11, a_tol12) = get_domain_parameters(cur_domain);
    let (a_first1, a_last1, _, _) = get_domain_parameters(cur_domain);
    let a_pos1a = find_position_ll(res_sup, cur_domain);
    let a_pos2a = find_position_ll(&mut a_res2, other_domain);
    let an_other_trans = if cur_trans == IntRes2dTypeTrans::Out {
        IntRes2dTypeTrans::In
    } else if cur_trans == IntRes2dTypeTrans::In {
        IntRes2dTypeTrans::Out
    } else {
        IntRes2dTypeTrans::Undecided
    };
    let (a_t1, a_t2) = if cur_trans != IntRes2dTypeTrans::Undecided {
        (
            IntRes2dTransition::in_out(false, a_pos1a, cur_trans),
            IntRes2dTransition::in_out(false, a_pos2a, an_other_trans),
        )
    } else {
        let an_opposite = cos_t1t2 < 0.0;
        (
            IntRes2dTransition::touch(false, a_pos1a, IntRes2dSituation::Unknown, an_opposite),
            IntRes2dTransition::touch(false, a_pos2a, IntRes2dSituation::Unknown, an_opposite),
        )
    };
    let mut a_res_u1 = par_cur;
    let mut a_res_u2 = par_other;
    let is_inside1 = par_cur >= a_first1 && par_cur <= a_last1;
    let is_inside2 = par_other >= a_first2 && par_other <= a_last2;
    if !is_inside1 || !is_inside2 {
        if is_inside1 {
            let pt1 = clib2d::line_value_ax2d(a_res2, other_lin.position());
            a_res_u2 = a_res2;
            let a_par1 = clib2d::line_parameter_ax2d(cur_lin.position(), &pt1);
            a_res_u1 = if a_par1 >= a_first1 && a_par1 <= a_last1 { a_par1 } else { *res_sup };
        } else if is_inside2 {
            let a_pt1 = clib2d::line_value_ax2d(*res_sup, cur_lin.position());
            a_res_u1 = *res_sup;
            let a_par2 = clib2d::line_parameter_ax2d(other_lin.position(), &a_pt1);
            a_res_u2 = if a_par2 >= a_first2 && a_par2 <= a_last2 { a_par2 } else { a_res2 };
        } else {
            if par_cur < a_first1 - a_tol11
                || par_cur > a_last1 + a_tol12
                || par_other < a_first2 - a_tol21
                || par_other > a_last2 + a_tol22
            {
                return None;
            }
            a_res_u1 = *res_sup;
            a_res_u2 = a_res2;
        }
    }
    let p1 = clib2d::line_value_ax2d(a_res_u1, cur_lin.position());
    let p2 = clib2d::line_value_ax2d(a_res_u2, other_lin.position());
    let a_pres = GpPnt2d::new(0.5 * (p1.x() + p2.x()), 0.5 * (p1.y() + p2.y()));
    let pt = if num == 1 {
        IntRes2dIntersectionPoint::with_transitions(&a_pres, a_res_u1, a_res_u2, &a_t1, &a_t2, false)
    } else {
        IntRes2dIntersectionPoint::with_transitions(&a_pres, a_res_u2, a_res_u1, &a_t2, &a_t1, false)
    };
    Some(pt)
}

/// `Precision::PConfusion()`, re-exported for the `Perform` parameters.
pub(crate) fn pconfusion() -> f64 {
    PCONFUSION
}

/// `gp_Vec2d` dot product helper (the OCCT `Tan1.Dot(Tan2)`).
pub(crate) fn dir_dot(a: &GpVec2d, b: &GpVec2d) -> f64 {
    a.dot(b)
}

/// `gp_Vec2d` cross product helper (the OCCT `Tan1.Crossed(Tan2)`).
pub(crate) fn dir_cross(a: &GpVec2d, b: &GpVec2d) -> f64 {
    a.x() * b.y() - a.y() * b.x()
}
