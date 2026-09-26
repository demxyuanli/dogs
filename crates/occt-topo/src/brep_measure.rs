//! Analytic measurements of a BRep model — edge lengths, wire perimeters,
//! face areas, shape totals.
//! Source: `BRep_Tool`, `GProp`-style integrals over the boundary.

use occt_core::gp::GpPnt;

use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, TopoShape, Wire};
use crate::topo_tools_full::{edges_of, edges_of_wire, faces_of, wires_of_face};

/// Length of a single edge: polyline length of its curve sampled uniformly.
pub fn edge_length(e: &Edge, samples: usize) -> f64 {
    let Some(curve) = BRepTool::edge_curve(e) else { return 0.0 };
    let (a, b) = BRepTool::edge_parameters(e);
    if !a.is_finite() || !b.is_finite() || b <= a {
        return 0.0;
    }
    let n = samples.max(2);
    let mut prev = curve.d0(a);
    let mut len = 0.0;
    for i in 1..=n {
        let u = a + (b - a) * i as f64 / n as f64;
        let p = curve.d0(u);
        len += prev.distance(&p);
        prev = p;
    }
    len
}

/// Total length of a wire's boundary edges.
pub fn wire_length(wire: &Wire, samples: usize) -> f64 {
    edges_of_wire(wire).iter().map(|e| edge_length(e, samples)).sum()
}

/// Surface area of a planar-ish face via a `nu × nv` grid triangulation.
/// More accurate than the naive surface patch integral when the face is
/// trimmed, but here it samples the full natural surface.
pub fn face_area(f: &Face, nu: usize, nv: usize) -> f64 {
    let Some(s) = BRepTool::face_surface(f) else { return 0.0 };
    let (u0, u1) = s.u_range();
    let (v0, v1) = s.v_range();
    if !(u0.is_finite() && u1.is_finite() && v0.is_finite() && v1.is_finite()) {
        // Unbounded surface (a plane) has infinite area.
        return f64::INFINITY;
    }
    let nu = nu.max(2);
    let nv = nv.max(2);
    let mut area = 0.0;
    for i in 0..nu {
        for j in 0..nv {
            let p00 = s.d0(u0 + (u1 - u0) * i as f64 / nu as f64, v0 + (v1 - v0) * j as f64 / nv as f64);
            let p10 = s.d0(u0 + (u1 - u0) * (i + 1) as f64 / nu as f64, v0 + (v1 - v0) * j as f64 / nv as f64);
            let p01 = s.d0(u0 + (u1 - u0) * i as f64 / nu as f64, v0 + (v1 - v0) * (j + 1) as f64 / nv as f64);
            let p11 = s.d0(u0 + (u1 - u0) * (i + 1) as f64 / nu as f64, v0 + (v1 - v0) * (j + 1) as f64 / nv as f64);
            let tri_area = |a: &GpPnt, b: &GpPnt, c: &GpPnt| {
                0.5 * occt_core::gp::GpVec::from_pnts(a, b)
                    .xyz()
                    .crossed(&occt_core::gp::GpVec::from_pnts(a, c).xyz())
                    .modulus()
            };
            area += tri_area(&p00, &p10, &p11) + tri_area(&p00, &p11, &p01);
        }
    }
    area
}

/// Perimeter of a shape: the sum of all its edges' lengths.
pub fn shape_perimeter(shape: &TopoShape, samples: usize) -> f64 {
    edges_of(shape).iter().map(|e| edge_length(e, samples)).sum()
}

/// Total surface area of a shape's faces (finite surfaces only; unbounded
/// planar faces contribute `INFINITY`).
pub fn shape_analytic_area(shape: &TopoShape, nu: usize, nv: usize) -> f64 {
    faces_of(shape).iter().map(|f| face_area(f, nu, nv)).sum()
}

/// Number of boundary wires in a shape.
pub fn wire_count(shape: &TopoShape) -> usize {
    faces_of(shape).iter().flat_map(wires_of_face).count()
}

/// Total boundary edge count (distinct edges, not occurrences).
pub fn boundary_edge_count(shape: &TopoShape) -> usize {
    edges_of(shape).len()
}

/// Length of the longest edge in a shape.
pub fn max_edge_length(shape: &TopoShape, samples: usize) -> f64 {
    edges_of(shape)
        .iter()
        .map(|e| edge_length(e, samples))
        .fold(0.0, f64::max)
}

/// Length of the shortest edge in a shape.
pub fn min_edge_length(shape: &TopoShape, samples: usize) -> f64 {
    edges_of(shape)
        .iter()
        .map(|e| edge_length(e, samples))
        .fold(f64::INFINITY, f64::min)
}

/// Whether the shape has any unbounded (infinite) edge or face — a hint that
/// it represents an infinite analytic surface rather than a closed solid.
pub fn has_infinite_geometry(shape: &TopoShape) -> bool {
    edges_of(shape)
        .iter()
        .any(|e| {
            let (a, b) = BRepTool::edge_parameters(e);
            !a.is_finite() || !b.is_finite()
        }) || faces_of(shape).iter().any(|f| {
            match BRepTool::face_surface(f) {
                Some(s) => {
                    let (u0, u1) = s.u_range();
                    let (v0, v1) = s.v_range();
                    !(u0.is_finite() && u1.is_finite() && v0.is_finite() && v1.is_finite())
                }
                None => true,
            }
        })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitives::BRepPrimBox;
    use crate::shape::Face;


    fn approx(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-3 * b.abs().max(1.0)
    }

    #[test]
    fn box_edges_and_perimeter() {
        let b = BRepPrimBox::make_box(2.0, 3.0, 4.0);
        // 4 edges of length 2, 4 of 3, 4 of 4.
        let all: Vec<f64> = edges_of(&b.solid.0).iter().map(|e| edge_length(e, 16)).collect();
        assert_eq!(all.len(), 12);
        for &l in &all {
            assert!(approx(l, 2.0) || approx(l, 3.0) || approx(l, 4.0), "edge length {l}");
        }
        let perim = shape_perimeter(&b.solid.0, 16);
        assert!(approx(perim, 4.0 * (2.0 + 3.0 + 4.0)));
        assert_eq!(boundary_edge_count(&b.solid.0), 12);
    }

    #[test]
    fn sphere_face_area_finite_and_plane_infinite() {
        use occt_core::gp::GpAx3;
        use occt_geom::GeomSphere;
        use std::sync::Arc;
        use crate::builder::TopoBuilder;

        let b = TopoBuilder::new();
        // Bounded surface: sphere → finite area ≈ 4πr².
        let sphere = b.make_face(Arc::new(GeomSphere::new(
            occt_core::gp::GpSphere::new(GpAx3::standard(), 2.0).unwrap(),
        )), &[]);
        let a = face_area(&sphere, 24, 24);
        assert!(a > 45.0 && a < 55.0, "sphere area {a} (expect ≈ 50.27)");

        // Unbounded surface: plane → infinite analytic area.
        let pln = b.make_face_plane(&occt_core::gp::GpPln::new(GpAx3::standard()));
        assert!(face_area(&pln, 8, 8).is_infinite());
        assert!(has_infinite_geometry(&pln.0));
    }

    #[test]
    fn min_max_edge() {
        let b = BRepPrimBox::make_box(1.0, 2.0, 5.0);
        assert!(approx(min_edge_length(&b.solid.0, 16), 1.0));
        assert!(approx(max_edge_length(&b.solid.0, 16), 5.0));
    }

    #[test]
    fn infinite_geometry_detection() {
        use occt_core::gp::GpAx3;
        use occt_geom::GeomSphere;
        use std::sync::Arc;
        use crate::builder::TopoBuilder;

        // Box edges are finite, but its faces are unbounded planes (natural
        // surface ranges) → the shape reports infinite geometry.
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        assert!(has_infinite_geometry(&b.solid.0));

        // A bounded surface (sphere) face reports finite.
        let tb = TopoBuilder::new();
        let sphere = tb.make_face(Arc::new(GeomSphere::new(
            occt_core::gp::GpSphere::new(GpAx3::standard(), 1.0).unwrap(),
        )), &[]);
        assert!(!has_infinite_geometry(&sphere.0));

        // An unregistered empty face has no surface → treated as infinite.
        let face = Face::new();
        assert!(has_infinite_geometry(&face.0));
    }
}
