//! `ChFi3d_ApproxByC2` — re-fit a curve that is not `C2` by a `C2` B-spline.
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_0.cxx:5835-5851`
//! (OCCT 8.0.0), declaration `ChFi3d_Builder_0.hxx:525-526`.
//!
//! OCCT samples the curve at 101 points (`cxx:5841-5847`: the first point is
//! `First`, then `First + i * delta` for `i = 0..99`, and the last slot is
//! `Last` again, so `Last` is present twice in the sequence) and runs
//! `GeomAPI_PointsToBSpline(Points, Approx_ChordLength, 3, 8, GeomAbs_C2,
//! 1.000001e-3)` (`cxx:5849-5850`).
//!
//! The `OccHandle<Geom_BSplineCurve>` return becomes an `Option`: every OCCT
//! caller (`cxx:4988-4991`, `:5075-5078`) assigns the result without a null
//! check, so an `Err` here is a hard failure at the call site rather than a
//! silent fallback.

use occt_core::approx::ApproxParametrizationType;
use occt_core::kernel::geomabs::Shape;
use occt_geom::bspline_curve::GeomBSplineCurve;
use occt_geom::geom_api_points_to_bspline::PointsToBSpline;
use occt_geom::Curve;

/// `ChFi3d_ApproxByC2(C)` (`cxx:5835-5851`).
pub fn approx_by_c2(c: &dyn Curve) -> Option<GeomBSplineCurve> {
    let first = c.first_parameter();
    let last = c.last_parameter();
    let nb_points = 101usize;
    let delta = (last - first) / (nb_points - 1) as f64;

    let mut points = Vec::with_capacity(nb_points);
    for i in 1..=nb_points - 1 {
        points.push(c.d0(first + (i - 1) as f64 * delta));
    }
    points.push(c.d0(last));

    let approx = PointsToBSpline::with_params_type(
        &points,
        ApproxParametrizationType::ChordLength,
        3,
        8,
        Shape::C2,
        1.000001e-3,
    );
    approx.curve().ok().cloned()
}
