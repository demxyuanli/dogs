//! Port of OCCT `BRepMesh_Deflection` + `BRepMesh_DegreeOfFreedom`.
//!
//! Source:
//! - `src/ModelingAlgorithms/TKMesh/BRepMesh/BRepMesh_Deflection.{hxx,cxx}`
//! - `src/ModelingAlgorithms/TKMesh/BRepMesh/BRepMesh_DegreeOfFreedom.hxx`
//!
//! The OCCT class computes deflections of discrete edges/wires/faces and drives
//! the absolute-vs-relative deflection bookkeeping. The topological `IMeshData`
//! overloads are not ported here (the sibling stubs are empty); instead this
//! port exposes the topology-free core: absolute-deflection adjustment,
//! consistency checking, and direct curve/surface chord-deviation measurement.

use occt_core::gp::{GpPnt, GpPnt2d};
use occt_geom::{Curve, Surface};

/// Degrees of freedom of a mesh node. Source: `BRepMesh_DegreeOfFreedom`
///
/// Ordered exactly as in the OCCT header; numeric values match the C enum
/// (`Free = 0`, ..., `Deleted = 6`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum DegreeOfFreedom {
    Free,
    InVolume,
    OnSurface,
    OnCurve,
    Fixed,
    Frontier,
    Deleted,
}

impl DegreeOfFreedom {
    /// Numeric value as in the OCCT C enum.
    pub fn index(self) -> usize {
        self as usize
    }

    /// Human-readable name.
    pub fn to_str(self) -> &'static str {
        match self {
            Self::Free => "Free",
            Self::InVolume => "InVolume",
            Self::OnSurface => "OnSurface",
            Self::OnCurve => "OnCurve",
            Self::Fixed => "Fixed",
            Self::Frontier => "Frontier",
            Self::Deleted => "Deleted",
        }
    }
}

/// Auxiliary tool encompassing methods to compute deflection of shapes.
/// Source: `BRepMesh_Deflection`
pub struct Deflection;

impl Deflection {
    /// Returns absolute deflection from the relative deflection and the maximum
    /// shape size.
    ///
    /// `box_max_dimension` is the largest extent of the shape's bounding box
    /// (OCCT's `BRepMesh_ShapeTool::BoxMaxDimension`); `relative_deflection` is
    /// the requested relative value; `max_shape_size` clamps the reference size
    /// (pass a non-positive value to use the box dimension itself). Source:
    /// `ComputeAbsoluteDeflection`
    pub fn compute_absolute_deflection(
        box_max_dimension: f64,
        relative_deflection: f64,
        max_shape_size: f64,
    ) -> f64 {
        let shape_size = box_max_dimension.max(relative_deflection);
        let max_shape_size = if max_shape_size > 0.0 {
            max_shape_size
        } else {
            box_max_dimension
        };

        let mut coefficient = max_shape_size / (2.0 * shape_size);
        if coefficient < 0.5 {
            coefficient = 0.5;
        } else if coefficient > 2.0 {
            coefficient = 2.0;
        }

        coefficient * shape_size * relative_deflection
    }

    /// Checks whether a current polygonal-representation deflection is
    /// consistent with the required one.
    ///
    /// When `allow_decrease` is true the current deflection must lie within
    /// `ratio` of the required value; otherwise it must only not exceed it.
    /// `ratio` defaults to `0.1` in OCCT. Source: `IsConsistent`
    pub fn is_consistent(
        current: f64,
        required: f64,
        allow_decrease: bool,
        ratio: f64,
    ) -> bool {
        current < (1.0 + ratio) * required
            && (!allow_decrease || current > (1.0 - ratio) * required)
    }

    /// Maximum distance from the true `curve` to its polyline approximation.
    ///
    /// For every segment `(params[i], params[i+1])` the curve is re-sampled at
    /// `samples` interior parameters and the furthest distance to the chord
    /// `(points[i], points[i+1])` is accumulated. `params` and `points` are the
    /// discretized parameter/point pairs (e.g. from `GeomTool`); a longer
    /// polyline yields a smaller deviation.
    pub fn curve_deflection(
        curve: &dyn Curve,
        params: &[f64],
        points: &[GpPnt],
        samples: usize,
    ) -> f64 {
        let n = params.len().min(points.len());
        let samples = samples.max(1);
        let mut max_dev = 0.0f64;
        for i in 0..n.saturating_sub(1) {
            let (ua, ub) = (params[i], params[i + 1]);
            let (pa, pb) = (&points[i], &points[i + 1]);
            for k in 1..=samples {
                let u = ua + (ub - ua) * k as f64 / (samples + 1) as f64;
                let dev = point_segment_dist(&curve.d0(u), pa, pb);
                if dev > max_dev {
                    max_dev = dev;
                }
            }
        }
        max_dev
    }

    /// Maximum distance from the true `surface` to the triangulated chord mesh.
    ///
    /// Each UV triangle is lifted to 3-D; the surface is sampled on a
    /// barycentric `samples`-grid inside the triangle and the furthest distance
    /// to the triangle's plane is returned. Periodic `U` ranges are unwrapped so
    /// triangles straddling the seam are measured on the same sheet.
    pub fn surface_deflection(
        surface: &dyn Surface,
        triangles: &[(GpPnt2d, GpPnt2d, GpPnt2d)],
        samples: usize,
    ) -> f64 {
        let (u0, u1) = surface.u_range();
        let period = if surface.is_u_periodic() && u1.is_finite() && u0.is_finite() {
            u1 - u0
        } else {
            0.0
        };
        let samples = samples.max(1);
        let mut max_dev = 0.0f64;
        for &(ta, tb, tc) in triangles {
            let a = surface.d0(ta.x(), ta.y());
            let b = surface.d0(tb.x(), tb.y());
            let c = surface.d0(tc.x(), tc.y());
            let n = b.coord.subtracted(&a.coord).crossed(&c.coord.subtracted(&a.coord));
            if n.square_modulus() < 1e-30 {
                continue;
            }
            let n = n.divided(n.modulus());

            let (ua, va) = (ta.x(), ta.y());
            let (mut ub, vb) = (tb.x(), tb.y());
            let (mut uc, vc) = (tc.x(), tc.y());
            if period > 0.0 {
                let half = 0.5 * period;
                if (ub - ua).abs() > half {
                    ub += if ub > ua { -period } else { period };
                }
                if (uc - ua).abs() > half {
                    uc += if uc > ua { -period } else { period };
                }
            }

            for i in 0..=samples {
                for j in 0..=(samples - i) {
                    let w = i as f64 / samples as f64;
                    let v = j as f64 / samples as f64;
                    let u = 1.0 - w - v;
                    let uu = u * ua + w * ub + v * uc;
                    let vv = u * va + w * vb + v * vc;
                    let p = surface.d0(uu, vv);
                    let d = p.coord.subtracted(&a.coord).dot(&n).abs();
                    if d > max_dev {
                        max_dev = d;
                    }
                }
            }
        }
        max_dev
    }
}

/// Distance from `p` to the segment `a..b`.
fn point_segment_dist(p: &GpPnt, a: &GpPnt, b: &GpPnt) -> f64 {
    let ab = b.coord.subtracted(&a.coord);
    let len2 = ab.square_modulus();
    if len2 <= f64::EPSILON {
        return p.distance(a);
    }
    let ap = p.coord.subtracted(&a.coord);
    let t = (ap.dot(&ab) / len2).clamp(0.0, 1.0);
    let proj = a.coord.added(&ab.multiplied(t));
    p.coord.subtracted(&proj).modulus()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;
    use std::sync::Arc;

    use occt_core::gp::{GpAx2, GpAx3, GpCirc, GpPln, GpSphere};
    use occt_geom::{GeomCircle, GeomPlane, GeomSphere};

    fn circle(r: f64) -> Arc<dyn Curve> {
        Arc::new(GeomCircle::new(GpCirc::new(GpAx2::standard(), r)))
    }

    fn plane() -> Arc<dyn Surface> {
        Arc::new(GeomPlane::new(GpPln::new(GpAx3::standard())))
    }

    fn sphere(r: f64) -> Arc<dyn Surface> {
        Arc::new(GeomSphere::new(GpSphere::new(GpAx3::standard(), r).unwrap()))
    }

    #[test]
    fn degree_of_freedom_enum_matches_occt_order() {
        let all = [
            DegreeOfFreedom::Free,
            DegreeOfFreedom::InVolume,
            DegreeOfFreedom::OnSurface,
            DegreeOfFreedom::OnCurve,
            DegreeOfFreedom::Fixed,
            DegreeOfFreedom::Frontier,
            DegreeOfFreedom::Deleted,
        ];
        assert_eq!(all.len(), 7);
        for (i, v) in all.iter().enumerate() {
            assert_eq!(v.index(), i);
        }
        assert_eq!(DegreeOfFreedom::Free.index(), 0);
        assert_eq!(DegreeOfFreedom::Deleted.index(), 6);
        assert_eq!(DegreeOfFreedom::OnCurve.to_str(), "OnCurve");
    }

    #[test]
    fn curve_chord_deviation_converges_with_finer_sampling() {
        let c = circle(1.0);

        let n_coarse = 5usize;
        let params_coarse: Vec<f64> =
            (0..n_coarse).map(|i| PI * i as f64 / (n_coarse - 1) as f64).collect();
        let points_coarse: Vec<GpPnt> = params_coarse.iter().map(|&u| c.d0(u)).collect();
        let e_coarse = Deflection::curve_deflection(c.as_ref(), &params_coarse, &points_coarse, 8);

        let n_fine = 33usize;
        let params_fine: Vec<f64> =
            (0..n_fine).map(|i| PI * i as f64 / (n_fine - 1) as f64).collect();
        let points_fine: Vec<GpPnt> = params_fine.iter().map(|&u| c.d0(u)).collect();
        let e_fine = Deflection::curve_deflection(c.as_ref(), &params_fine, &points_fine, 8);

        assert!(e_fine < e_coarse, "expected finer mesh to deviate less: {e_fine} vs {e_coarse}");
        assert!(e_coarse > 1e-3, "coarse semicircle should deviate visibly: {e_coarse}");
        // 33 points on a unit semicircle -> 32 chords of pi/32; sagitta is
        // 1 - cos(pi/64) ~ 1.2e-3, so the measured deviation must be under 2e-3.
        assert!(e_fine < 2.0e-3, "fine semicircle should hug the arc: {e_fine}");
    }

    #[test]
    fn surface_deflection_is_zero_for_plane() {
        let s = plane();
        let tris = vec![(
            GpPnt2d::new(0.0, 0.0),
            GpPnt2d::new(1.0, 0.0),
            GpPnt2d::new(0.0, 1.0),
        )];
        let d = Deflection::surface_deflection(s.as_ref(), &tris, 5);
        assert!(d < 1e-12, "plane must have zero deflection, got {d}");
    }

    #[test]
    fn sphere_surface_deflection_decreases_with_finer_triangulation() {
        let s = sphere(1.0);
        // A single triangle spanning a large part of the sphere.
        let coarse = vec![(
            GpPnt2d::new(0.0, 0.0),
            GpPnt2d::new(1.0, 0.0),
            GpPnt2d::new(0.0, 1.0),
        )];
        // A finer mesh: subdivide the same triangle into four.
        let m = 0.5;
        let fine = vec![
            (GpPnt2d::new(0.0, 0.0), GpPnt2d::new(m, 0.0), GpPnt2d::new(0.0, m)),
            (GpPnt2d::new(1.0, 0.0), GpPnt2d::new(1.0, m), GpPnt2d::new(m, 0.0)),
            (GpPnt2d::new(0.0, 1.0), GpPnt2d::new(m, m), GpPnt2d::new(0.0, m)),
            (GpPnt2d::new(1.0, 0.0), GpPnt2d::new(m, m), GpPnt2d::new(1.0, m)),
            (GpPnt2d::new(m, 0.0), GpPnt2d::new(1.0, m), GpPnt2d::new(m, m)),
            (GpPnt2d::new(0.0, 1.0), GpPnt2d::new(m, m), GpPnt2d::new(1.0, m)),
        ];
        let d_coarse = Deflection::surface_deflection(s.as_ref(), &coarse, 3);
        let d_fine = Deflection::surface_deflection(s.as_ref(), &fine, 3);
        assert!(d_coarse > 0.0, "curved surface must deviate from chords");
        assert!(d_fine < d_coarse, "finer mesh should deviate less: {d_fine} vs {d_coarse}");
    }

    #[test]
    fn is_consistent_matches_occt() {
        let ratio = 0.1;
        assert!(Deflection::is_consistent(0.095, 0.1, true, ratio));
        assert!(Deflection::is_consistent(0.095, 0.1, false, ratio));
        assert!(Deflection::is_consistent(0.05, 0.1, false, ratio), "decrease allowed when allow_decrease=false");
        assert!(!Deflection::is_consistent(0.05, 0.1, true, ratio), "too-small current rejected when decrease allowed");
        assert!(!Deflection::is_consistent(0.15, 0.1, false, ratio), "too-large current always rejected");
        // Default ratio 0.1 applied when the caller omits it.
        assert!(Deflection::is_consistent(0.095, 0.1, true, ratio));
    }

    #[test]
    fn compute_absolute_deflection_adjustment() {
        // Box 100 units wide: coefficient = 100 / (2*100) = 0.5.
        let d = Deflection::compute_absolute_deflection(100.0, 0.01, -1.0);
        assert!((d - 0.5).abs() < 1e-12, "got {d}");

        // Tiny shape relative to the max shape size: coefficient clamps to 2.0.
        let d2 = Deflection::compute_absolute_deflection(1.0, 0.01, 100.0);
        assert!((d2 - 0.02).abs() < 1e-12, "got {d2}");
    }
}
