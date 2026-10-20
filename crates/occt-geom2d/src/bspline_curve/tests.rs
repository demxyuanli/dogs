use super::*;
use crate::curve::Curve2d;

#[test]
fn linear_midpoint() {
    let c = Geom2dBSplineCurve::new(
        vec![0.0, 1.0],
        vec![0.0, 0.0],
        vec![0.0, 0.0, 1.0, 1.0],
        1,
    ).unwrap();
    let p = c.d0(0.5);
    assert!((p.x() - 0.5).abs() < 1e-12);
    assert!((p.y() - 0.0).abs() < 1e-12);
}
