//! The two `ChFi3d_Boite` overloads — UV boxes of a fillet to intersect.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_0.cxx:242-284`.

use occt_core::gp::GpPnt2d;

/// `ChFi3d_Boite(p1, p2, mu, Mu, mv, Mv)` (`cxx:242-257`): the UV box of two
/// points, returned as `(mu, Mu, mv, Mv)`.
pub fn boite(p1: &GpPnt2d, p2: &GpPnt2d) -> (f64, f64, f64, f64) {
    let mu = p1.x().min(p2.x());
    let big_mu = p1.x().max(p2.x());
    let mv = p1.y().min(p2.y());
    let big_mv = p1.y().max(p2.y());
    (mu, big_mu, mv, big_mv)
}

/// `ChFi3d_Boite(p1, p2, p3, p4, Du, Dv, mu, Mu, mv, Mv)` (`cxx:259-284`): the
/// UV box of four points together with its extents, returned as
/// `(Du, Dv, mu, Mu, mv, Mv)`.
pub fn boite_of_four(
    p1: &GpPnt2d,
    p2: &GpPnt2d,
    p3: &GpPnt2d,
    p4: &GpPnt2d,
) -> (f64, f64, f64, f64, f64, f64) {
    let a = p1.x().min(p2.x());
    let b = p3.x().min(p4.x());
    let mu = a.min(b);
    let a = p1.x().max(p2.x());
    let b = p3.x().max(p4.x());
    let big_mu = a.max(b);
    let a = p1.y().min(p2.y());
    let b = p3.y().min(p4.y());
    let mv = a.min(b);
    let a = p1.y().max(p2.y());
    let b = p3.y().max(p4.y());
    let big_mv = a.max(b);
    (big_mu - mu, big_mv - mv, mu, big_mu, mv, big_mv)
}
