use super::*;
use crate::curve::Curve;
use crate::{GeomCircle, GeomLine, GeomTrimmedCurve};
use occt_core::elib::clib;
use occt_core::gp::{GpAx2, GpCirc, GpDir, GpPnt};
use std::sync::Arc;

const PI: f64 = std::f64::consts::PI;

fn unit_circle() -> GeomCircle {
    GeomCircle::new(GpCirc::new(GpAx2::standard(), 1.0))
}

fn unit_circle_arc() -> GeomTrimmedCurve {
    GeomTrimmedCurve::new(Arc::new(unit_circle()), 0.0, PI / 2.0)
}

#[test]
fn length_unit_circle() {
    let c = unit_circle();
    let l = curve_length(&c);
    assert!((l - 2.0 * PI).abs() < 1e-7, "length {l}");
}

#[test]
fn length_quarter_arc() {
    let c = unit_circle_arc();
    let l = curve_length(&c);
    assert!((l - PI / 2.0).abs() < 1e-7, "length {l}");
}

#[test]
fn abscissa_point_half_circle() {
    let c = unit_circle();
    // Arc length π from u=0 lands at u=π.
    let u = abscissa_point(&c, PI, 0.0).unwrap();
    assert!((u - PI).abs() < 1e-6, "u {u}");
    let p = c.d0(u);
    assert!(p.distance(&GpPnt::new(-1.0, 0.0, 0.0)) < 1e-6, "point {p:?}");
}

#[test]
fn abscissa_point_total_length_returns_last_parameter() {
    let c = unit_circle();
    let u = abscissa_point(&c, 2.0 * PI, 0.0).unwrap();
    assert!((u - 2.0 * PI).abs() < 1e-6, "u {u}");
    // `CPnts_AbscissaPoint.hxx:78`: "The computed point can be outside
    // of the curve 's bounds". `CPnts_AbscissaPoint::Perform` widens the
    // root bracket to `[U1 - DU, U2 + DU]` (`CPnts_AbscissaPoint.cxx:310-312`),
    // so an abscissa past the total length is not rejected: on the periodic
    // circle the integral `Length(0, X)` still matches it at `X = 3*pi`.
    let u = abscissa_point(&c, 3.0 * PI, 0.0).unwrap();
    assert!((u - 3.0 * PI).abs() < 1e-6, "u {u}");
}

#[test]
fn abscissa_point_negative_goes_backward() {
    let c = unit_circle();
    // Backward by π/2 from u=π lands at u=π/2 (the point (0,1,0)).
    let u = abscissa_point(&c, -PI / 2.0, PI).unwrap();
    assert!((u - PI / 2.0).abs() < 1e-6, "u {u}");
    let p = c.d0(u);
    assert!(p.distance(&GpPnt::new(0.0, 1.0, 0.0)) < 1e-6, "point {p:?} at u {u}");
}

#[test]
fn uniform_abscissa_unit_circle_quarter_steps() {
    let c = unit_circle();
    let params = uniform_abscissa(&c, 4).unwrap();
    assert_eq!(params.len(), 5);
    let expected = [0.0, PI / 2.0, PI, 3.0 * PI / 2.0, 2.0 * PI];
    for (g, e) in params.iter().zip(expected.iter()) {
        assert!((g - e).abs() < 1e-6, "got {g} expected {e}: {params:?}");
    }
    // Points are equally spaced on the circle.
    let pts: Vec<GpPnt> = params.iter().map(|&u| c.d0(u)).collect();
    for w in pts.windows(2) {
        let chord = w[0].distance(&w[1]);
        assert!((chord - (2.0f64).sqrt()).abs() < 1e-5, "chord {chord}");
    }
}

#[test]
fn uniform_abscissa_line_is_uniform() {
    // A 10-unit trimmed line has parameter range [0, 1]; equal arc-length
    // spacing of 10/4 = 2.5 lands at params 0, 0.25, 0.5, 0.75, 1.
    let line = GeomLine::from_pnt_dir(GpPnt::new(0., 0., 0.), GpDir::new(1., 0., 0.).unwrap());
    let c = GeomTrimmedCurve::new(Arc::new(line), 0.0, 10.0);
    let params = uniform_abscissa(&c, 4).unwrap();
    for (i, u) in params.iter().enumerate() {
        assert!((u - 0.25 * i as f64).abs() < 1e-9, "params {params:?}");
    }
}

#[test]
fn uniform_abscissa_degenerate_zero_length() {
    let line = GeomLine::from_pnt_dir(GpPnt::new(1., 0., 0.), GpDir::new(1., 0., 0.).unwrap());
    let c = GeomTrimmedCurve::new(Arc::new(line), 0.0, 0.0);
    let params = uniform_abscissa(&c, 3).unwrap();
    assert_eq!(params, vec![0.0, 0.0, 0.0, 0.0]);
}

#[test]
fn quasi_uniform_abscissa_line() {
    // 10-unit trimmed line over [0, 1]: 5 points at equal chord intervals
    // 2.5 → params 0, 0.25, 0.5, 0.75, 1.
    let line = GeomLine::from_pnt_dir(GpPnt::new(0., 0., 0.), GpDir::new(1., 0., 0.).unwrap());
    let c = GeomTrimmedCurve::new(Arc::new(line), 0.0, 10.0);
    let params = quasi_uniform_abscissa(&c, 5).unwrap();
    assert_eq!(params.len(), 5);
    assert!((params[0] - 0.0).abs() < 1e-9);
    assert!((params[4] - 1.0).abs() < 1e-9);
    assert!((params[1] - 0.25).abs() < 1e-6, "params {params:?}");
    assert!((params[2] - 0.5).abs() < 1e-6, "params {params:?}");
    assert!((params[3] - 0.75).abs() < 1e-6, "params {params:?}");
}

#[test]
fn sampled_circle_points_consistent() {
    let c = unit_circle();
    let params = uniform_abscissa(&c, 8).unwrap();
    for &u in &params {
        let p = c.d0(u);
        let q = clib::circle_value(c.circ(), u);
        assert!(p.distance(&q) < 1e-12);
    }
}
