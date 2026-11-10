
use super::super::interference::Interference;
use super::super::kind::Kind;
use super::super::transition::Transition;


/// `FDS_Tdata(I, SB, IB, SA, IA)` (`cxx:144-162`): the four transition sides.
pub fn fds_tdata(t: &Transition) -> (crate::abs::ShapeType, i32, crate::abs::ShapeType, i32) {
    (
        t.shape_before(),
        t.index_before(),
        t.shape_after(),
        t.index_after(),
    )
}

/// `FDS_data(I, GT1, G1, ST1, S1)` (`cxx:164-181`).
pub fn fds_data(i: &Interference) -> (Kind, i32, Kind, i32) {
    (i.geometry_type(), i.geometry(), i.support_type(), i.support())
}

/// `FDS_data(it1, I1, GT1, G1, ST1, S1)` (`cxx:184-203`): fills the four
/// outputs and reports whether the transition is `(UNKNOWN, UNKNOWN)`.
pub fn fds_data_is_unknown(i: &Interference) -> bool {
    i.transition().is_unknown()
}

/// `FDS_Idata(I, SB, IB, SA, IA, GT1, G1, ST1, S1)` (`cxx:205-223`): the
/// transition and data quadruples in one call.
pub fn fds_idata(
    i: &Interference,
) -> (
    crate::abs::ShapeType,
    i32,
    crate::abs::ShapeType,
    i32,
    Kind,
    i32,
    Kind,
    i32,
) {
    let (sb, ib, sa, ia) = fds_tdata(&i.transition());
    let (gt, g, st, s) = fds_data(i);
    (sb, ib, sa, ia, gt, g, st, s)
}

