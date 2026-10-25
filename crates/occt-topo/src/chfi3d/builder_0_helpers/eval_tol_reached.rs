//! `ChFi3d_EvalTolReached` — the largest 3D gap between an intersection curve
//! and its two p-curves, sampled at 45 parameters.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_0.cxx:3555-3594`
//! (OCCT 8.0.0).
//!
//! Translation notes:
//! - `Adaptor3d_Surface` becomes a pair of `(Surface, SurfaceAdaptor)` views;
//!   only `Value(u, v)` and the surface handed to the p-curve's own domain is
//!   read, so the plain `Surface` trait is enough,
//! - `pc1->Value(t).Coord(u, v)` and `S1->Value(u, v)` are `Curve2d::value` and
//!   `Surface::d0` in the port.

use occt_core::precision::CONFUSION;
use occt_geom::curve::Curve;
use occt_geom::surface::Surface;
use occt_geom2d::curve::Curve2d;

/// `ChFi3d_EvalTolReached(S1, pc1, S2, pc2, C)` (`cxx:3555-3594`).
pub fn ch_fi3d_eval_tol_reached(
    s1: &dyn Surface,
    pc1: &dyn Curve2d,
    s2: &dyn Surface,
    pc2: &dyn Curve2d,
    c: &dyn Curve,
) -> f64 {
    let mut distmax: f64 = 0.0;

    let f = c.first_parameter();
    let l = c.last_parameter();
    let nbp = 45;
    let step = 1.0 / (nbp - 1) as f64;
    for i in 0..nbp {
        // `t = step * i; t = (1 - t) * f + t * l;` (`cxx:3570-3572`).
        let t = step * i as f64;
        let t = (1.0 - t) * f + t * l;
        // `pc1->Value(t).Coord(u, v); gp_Pnt pS1 = S1->Value(u, v);`
        let uv = pc1.value(t);
        let p_s1 = s1.d0(uv.x(), uv.y());
        // `pc2->Value(t).Coord(u, v); gp_Pnt pS2 = S2->Value(u, v);`
        let uv = pc2.value(t);
        let p_s2 = s2.d0(uv.x(), uv.y());
        // `gp_Pnt pC = C->Value(t);`
        let p_c = c.d0(t);
        let mut d = p_s1.square_distance(&p_c);
        if d > distmax {
            distmax = d;
        }
        d = p_s2.square_distance(&p_c);
        if d > distmax {
            distmax = d;
        }
        d = p_s1.square_distance(&p_s2);
        if d > distmax {
            distmax = d;
        }
    }
    // `distmax = 1.5 * sqrt(distmax); distmax = max(distmax,
    // Precision::Confusion());` (`cxx:3592-3593`).
    let distmax = 1.5 * distmax.sqrt();
    distmax.max(CONFUSION)
}
