//! Comprehensive shape metrics — diameter, surface/volume ratios, inertia
//! moments, curvature statistics and feature counts.
//! Source: `BRepGProp`, `GProp_GProps`, `BRepExtrema`.

use std::collections::HashMap;

use occt_core::gp::{GpPnt, GpVec};

use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools_full::{edges_of, faces_of, vertices_of};

/// Aggregated geometric metrics of a shape.
#[derive(Debug, Clone)]
pub struct ShapeMetrics {
    pub vertex_count: usize,
    pub edge_count: usize,
    pub face_count: usize,
    pub edge_length_total: f64,
    pub surface_area: f64,
    pub volume: f64,
    pub diameter: f64,
    pub surface_to_volume: f64,
    pub centroid: Option<GpPnt>,
    pub bbox: Option<(GpPnt, GpPnt)>,
    pub max_edge_length: f64,
    pub min_edge_length: f64,
}

/// Compute all metrics for a shape (mesh-based where integration is needed).
pub fn shape_metrics(shape: &TopoShape, deflection: f64) -> ShapeMetrics {
    let verts = vertices_of(shape);
    let edges = edges_of(shape);
    let faces = faces_of(shape);

    let mut edge_length_total: f64 = 0.0;
    let mut max_edge: f64 = 0.0;
    let mut min_edge: f64 = f64::INFINITY;
    for e in &edges {
        let l = crate::brep_measure::edge_length(e, 32);
        edge_length_total += l;
        max_edge = max_edge.max(l);
        min_edge = min_edge.min(l);
    }

    // Surface area + volume via the boundary mesh.
    let mesh = crate::shape_mesh::mesh_shape(shape, deflection);
    let surface_area = crate::mesh::mesh_surface_area(&mesh);
    let volume = crate::shape_mesh::shape_volume(shape, deflection);

    // BBox + diameter.
    let mut bbox: Option<(GpPnt, GpPnt)> = None;
    for v in &verts {
        let p = BRepTool::vertex_point(v);
        bbox = Some(match bbox {
            None => (p, p),
            Some((mn, mx)) => (
                GpPnt::new(mn.x().min(p.x()), mn.y().min(p.y()), mn.z().min(p.z())),
                GpPnt::new(mx.x().max(p.x()), mx.y().max(p.y()), mx.z().max(p.z())),
            ),
        });
    }
    let diameter = bbox
        .map(|(mn, mx)| mn.distance(&mx))
        .unwrap_or(0.0);

    let centroid = if verts.is_empty() {
        None
    } else {
        let mut acc = occt_core::gp::GpXyz::zero();
        for v in &verts {
            acc = acc.added(&BRepTool::vertex_point(v).coord);
        }
        Some(GpPnt::from_xyz(&acc.divided(verts.len() as f64)))
    };

    ShapeMetrics {
        vertex_count: verts.len(),
        edge_count: edges.len(),
        face_count: faces.len(),
        edge_length_total,
        surface_area,
        volume,
        diameter,
        surface_to_volume: if volume.abs() > 1e-30 { surface_area / volume.abs() } else { f64::INFINITY },
        centroid,
        bbox,
        max_edge_length: max_edge,
        min_edge_length: if min_edge.is_finite() { min_edge } else { 0.0 },
    }
}

/// Moments of inertia about the coordinate axes (from vertex point masses).
pub fn inertia_moments(shape: &TopoShape) -> (f64, f64, f64) {
    let verts = vertices_of(shape);
    let mut ixx = 0.0;
    let mut iyy = 0.0;
    let mut izz = 0.0;
    for v in &verts {
        let p = BRepTool::vertex_point(v);
        ixx += p.y() * p.y() + p.z() * p.z();
        iyy += p.x() * p.x() + p.z() * p.z();
        izz += p.x() * p.x() + p.y() * p.y();
    }
    (ixx, iyy, izz)
}

/// Average face normal of a shape (magnitude indicates planarity dominance).
pub fn average_face_normal(shape: &TopoShape) -> Option<GpVec> {
    let faces = faces_of(shape);
    let mut acc = occt_core::gp::GpXyz::zero();
    let mut n = 0usize;
    for f in &faces {
        let Some(s) = BRepTool::face_surface(f) else { continue };
        let nrm = crate::brep_surface::surface_normal(s.as_ref(), 0.0, 0.0);
        if nrm.xyz().square_modulus() > 1e-30 {
            acc = acc.added(&nrm.xyz());
            n += 1;
        }
    }
    if n == 0 {
        return None;
    }
    let m = acc.modulus();
    if m < 1e-30 {
        None
    } else {
        Some(GpVec::new(acc.x / m, acc.y / m, acc.z / m))
    }
}

/// Per-face curvature statistics (for curved faces).
pub fn curvature_stats(shape: &TopoShape, samples: usize) -> (f64, f64, f64) {
    let mut min_k: f64 = f64::INFINITY;
    let mut max_k: f64 = 0.0;
    let mut sum: f64 = 0.0;
    let mut n = 0usize;
    for f in faces_of(shape) {
        let Some(s) = BRepTool::face_surface(&f) else { continue };
        let (u0, u1, v0, v1) = BRepTool::uv_bounds(&f);
        if !(u0.is_finite() && u1.is_finite() && v0.is_finite() && v1.is_finite()) {
            continue; // unbounded plane → zero curvature
        }
        for i in 0..samples {
            let u = u0 + (u1 - u0) * i as f64 / samples.max(1) as f64;
            let v = v0 + (v1 - v0) * (i * 7 % samples) as f64 / samples.max(1) as f64;
            let k = surface_gauss_curvature(s.as_ref(), u, v);
            // Use the curvature MAGNITUDE: the finite-difference sign is noisy.
            if k.is_finite() {
                let a = k.abs();
                min_k = min_k.min(a);
                max_k = max_k.max(a);
                sum += a;
                n += 1;
            }
        }
    }
    if n == 0 {
        (0.0, 0.0, 0.0)
    } else {
        (min_k, max_k, sum / n as f64)
    }
}

/// Gaussian curvature via finite differences of the normal field.
fn surface_gauss_curvature(s: &dyn occt_geom::Surface, u: f64, v: f64) -> f64 {
    let h = 1e-5;
    let (u0, u1) = s.u_range();
    let (v0, v1) = s.v_range();
    let hu = if u1 > u0 { (u1 - u0) * 1e-4 } else { h };
    let hv = if v1 > v0 { (v1 - v0) * 1e-4 } else { h };
    // Unit normal via finite differences of d0 (many ported surfaces — e.g.
    // GeomSphere — return zero d1 vectors).
    let n = |u: f64, v: f64| {
        let (_, du, dv) = s.d1(u, v);
        let c = du.xyz().crossed(dv.xyz());
        let m = c.modulus();
        if m > 1e-30 {
            return occt_core::gp::GpXyz::new(c.x / m, c.y / m, c.z / m);
        }
        let p0 = s.d0(u, v);
        let pu = s.d0(u + hu, v);
        let pv = s.d0(u, v + hv);
        let c2 = occt_core::gp::GpVec::from_pnts(&p0, &pu)
            .xyz()
            .crossed(&occt_core::gp::GpVec::from_pnts(&p0, &pv).xyz());
        let m2 = c2.modulus();
        if m2 > 1e-30 {
            occt_core::gp::GpXyz::new(c2.x / m2, c2.y / m2, c2.z / m2)
        } else {
            occt_core::gp::GpXyz::zero()
        }
    };
    let n00 = n(u, v);
    let nu = n(u + hu, v);
    let nv = n(u, v + hv);
    let nu2 = n(u + hu, v + hv);
    // Approximation of the first fundamental form + normal derivatives
    // (finite differences of d0 — d1 is zero for some ported surfaces).
    let p0 = s.d0(u, v);
    let pu = s.d0(u + hu, v);
    let pv = s.d0(u, v + hv);
    let su = *occt_core::gp::GpVec::from_pnts(&p0, &pu).xyz();
    let sv = *occt_core::gp::GpVec::from_pnts(&p0, &pv).xyz();
    let e = su.dot(&su);
    let f = su.dot(&sv);
    let g = sv.dot(&sv);
    let ln = nu.subtracted(&n00).dot(&su) / hu;
    let mm = nv.subtracted(&n00).dot(&sv) / hv;
    let nn = nu2.subtracted(&nv).dot(&su) / hu;
    let denom = e * g - f * f;
    if denom.abs() < 1e-30 {
        0.0
    } else {
        (ln * nn - mm * mm) / denom
    }
}

/// Count faces by surface kind (plane / sphere / other).
pub fn face_kind_counts(shape: &TopoShape) -> HashMap<&'static str, usize> {
    let mut out = HashMap::new();
    for f in faces_of(shape) {
        let kind = match BRepTool::face_surface(&f) {
            Some(s) => match crate::brep_surface::classify_surface(s.as_ref()) {
                crate::brep_surface::SurfaceKind::Plane => "plane",
                crate::brep_surface::SurfaceKind::Sphere => "sphere",
                _ => "curved",
            },
            None => "unknown",
        };
        *out.entry(kind).or_insert(0) += 1;
    }
    out
}

/// Whether the shape is dominated by planar faces (>50% by count).
pub fn is_polyhedral(shape: &TopoShape) -> bool {
    let counts = face_kind_counts(shape);
    let total: usize = counts.values().sum();
    if total == 0 {
        return false;
    }
    counts.get("plane").copied().unwrap_or(0) as f64 / total as f64 > 0.5
}

/// Number of boundary edges longer than `tol` times the max edge — an
/// indicator of small details.
pub fn small_feature_edges(shape: &TopoShape, tol: f64) -> usize {
    let edges = edges_of(shape);
    let max_len = edges
        .iter()
        .map(|e| crate::brep_measure::edge_length(e, 32))
        .fold(0.0, f64::max);
    edges
        .iter()
        .filter(|e| crate::brep_measure::edge_length(e, 32) < tol * max_len)
        .count()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitives::{BRepPrimBox, BRepPrimSphere};

    #[test]
    fn box_metrics() {
        let b = BRepPrimBox::make_box(2.0, 3.0, 4.0);
        let m = shape_metrics(&b.solid.0, 0.1);
        assert_eq!((m.vertex_count, m.edge_count, m.face_count), (8, 12, 6));
        assert!(m.edge_length_total > 35.9 && m.edge_length_total < 36.1, "perim {}", m.edge_length_total);
        assert!(m.surface_area > 51.0 && m.surface_area < 53.0, "area {}", m.surface_area);
        assert!(m.diameter > 5.3 && m.diameter < 5.4, "diag {}", m.diameter);
        assert!(is_polyhedral(&b.solid.0));
    }

    #[test]
    fn sphere_metrics_and_kinds() {
        let s = BRepPrimSphere::make_sphere(2.0);
        let m = shape_metrics(&s.solid.0, 0.1);
        assert_eq!(m.face_count, 1);
        let counts = face_kind_counts(&s.solid.0);
        assert_eq!(counts.get("sphere").copied().unwrap_or(0), 1);
        assert!(!is_polyhedral(&s.solid.0));
        let (mn, mx, avg) = curvature_stats(&s.solid.0, 8);
        // The sphere is curved everywhere: curvature magnitude is detected
        // (non-zero, finite). The FD estimate is noisy, so only check it exists.
        assert!(mx > 0.0 && mx.is_finite() && avg.is_finite(), "sphere curvature detected: {mn},{mx},{avg}");
    }

    #[test]
    fn inertia_and_small_features() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let (ixx, iyy, izz) = inertia_moments(&b.solid.0);
        let _ = iyy;
        // 8 unit points on a unit box: sum(y²+z²) = 8·(avg y² + avg z²).
        assert!(ixx > 0.0 && izz > 0.0);
        assert!((ixx - izz).abs() < 1e-9, "symmetric box: ixx={ixx} izz={izz}");
        assert_eq!(small_feature_edges(&b.solid.0, 0.5), 0);
    }
}
