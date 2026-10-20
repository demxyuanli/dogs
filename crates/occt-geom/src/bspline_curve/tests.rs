use super::*;
use crate::curve::Curve;
use occt_core::gp::GpPnt;

#[test]
fn linear_midpoint() {
    let c = GeomBSplineCurve::new(
        vec![GpPnt::new(0., 0., 0.), GpPnt::new(1., 0., 0.)],
        vec![0., 0., 1., 1.],
        1,
    ).unwrap();
    let p = c.d0(0.5);
    assert!((p.x() - 0.5).abs() < 1e-12);
    assert!((p.y() - 0.0).abs() < 1e-12);
    assert!((p.z() - 0.0).abs() < 1e-12);
}
