//! Closed-form quadric pair intersections (plane/sphere/cylinder/cone/torus).
//! Source: analytic helpers previously in `intpatch.rs`.

use std::sync::Arc;

use occt_core::gp::{
    GpAx1, GpAx2, GpAx3, GpCirc, GpCone, GpCylinder, GpDir, GpElips, GpPln, GpPnt, GpSphere, GpTorus,
    GpVec,
};
use occt_geom::{
    Curve, GeomCircle, GeomCylinder, GeomEllipse, GeomLine, GeomPlane, GeomSphere, GeomTorus,
    Surface,
};

use super::geom::project_params;
use super::IntersectionCurve;

// ---------------------------------------------------------------------------
// Curve construction helpers
// ---------------------------------------------------------------------------

/// A `GeomCircle` curve with the given center, plane normal and radius.
fn circle_curve(center: GpPnt, normal: GpDir, radius: f64) -> Result<Arc<dyn Curve>, String> {
    let z = GpDir::new(0.0, 0.0, 1.0).map_err(|_| "z dir")?;
    let x_dir = if normal.is_normal(&z) { z } else { GpDir::new(1.0, 0.0, 0.0).map_err(|_| "x dir")? };
    let ax2 = GpAx2::new(center, normal, x_dir).map_err(|e| e.to_string())?;
    Ok(Arc::new(GeomCircle::new(GpCirc::new(ax2, radius))))
}

/// An in-plane unit direction perpendicular to `normal`, for building frames.
fn perpendicular_dir(normal: &GpDir) -> Result<GpDir, String> {
    let z = GpDir::new(0.0, 0.0, 1.0).map_err(|_| "z dir")?;
    Ok(if normal.is_normal(&z) { z } else { GpDir::new(1.0, 0.0, 0.0).map_err(|_| "x dir")? })
}

/// Build an `IntersectionCurve` by sampling `curve` on a `n`-point grid and
/// projecting every sample onto both surfaces.
pub(crate) fn sample_curve_on(
    curve: Arc<dyn Curve>,
    a: &dyn Surface,
    b: &dyn Surface,
    n: usize,
) -> IntersectionCurve {
    let (f0, f1) = (curve.first_parameter(), curve.last_parameter());
    let n = n.max(2);
    let mut points = Vec::with_capacity(n);
    let mut on_a = Vec::with_capacity(n);
    let mut on_b = Vec::with_capacity(n);
    for i in 0..n {
        let t = if (f1 - f0).is_finite() {
            f0 + (f1 - f0) * i as f64 / (n - 1) as f64
        } else {
            -5.0 + 10.0 * i as f64 / (n - 1) as f64
        };
        let p = curve.d0(t);
        points.push(p);
        on_a.push(project_params(a, &p));
        on_b.push(project_params(b, &p));
    }
    IntersectionCurve { curve, points, on_a, on_b }
}

/// T-28 step 4 adapter: `intana`'s analytic conic -> the sampled
/// `IntersectionCurve`s `surface_surface_intersection` returns. Same job as
/// `FaceFace::conics_to_curves` (`int_face_face_analytic.rs:408`) but at the
/// `IntersectionCurve` level. Sampling densities mirror the closed forms this
/// replaces (48 for circles/ellipses, 32 for lines).
pub fn ic_list_from_quadric(
    qi: occt_geom::intana::QuadricIntersection,
    a: &dyn Surface,
    b: &dyn Surface,
) -> Vec<IntersectionCurve> {
    use occt_geom::intana::QuadricIntersection::*;
    let line = |l: occt_core::gp::GpLin| sample_curve_on(Arc::new(GeomLine::new(l)) as Arc<dyn Curve>, a, b, 32);
    let circ = |c: occt_core::gp::GpCirc| sample_curve_on(Arc::new(GeomCircle::new(c)) as Arc<dyn Curve>, a, b, 48);
    let elps = |e: occt_core::gp::GpElips| sample_curve_on(Arc::new(GeomEllipse::new(e)) as Arc<dyn Curve>, a, b, 48);
    match qi {
        Line(l) => vec![line(l)],
        TwoLines(l1, l2) => vec![line(l1), line(l2)],
        Circle(c) => vec![circ(c)],
        TwoCircles(c1, c2) => vec![circ(c1), circ(c2)],
        Ellipse(e) => vec![elps(e)],
        TwoEllipses(e1, e2) => vec![elps(e1), elps(e2)],
        Parabola(_) | Hyperbola(_) | Point(_) | Same | None => Vec::new(),
    }
}

/// T-28 step 4 adapter: `intana`'s torus result -> the sampled
/// `IntersectionCurve`s `surface_surface_intersection` returns. Sampling
/// density (32) mirrors the closed form this replaces.
pub fn ic_list_from_torus(
    qi: occt_geom::intana::TorusIntersection,
    a: &dyn Surface,
    b: &dyn Surface,
) -> Vec<IntersectionCurve> {
    match qi {
        occt_geom::intana::TorusIntersection::Circles(cs) => cs
            .into_iter()
            .map(|c| sample_curve_on(Arc::new(GeomCircle::new(c)) as Arc<dyn Curve>, a, b, 32))
            .collect(),
        _ => Vec::new(),
    }
}
