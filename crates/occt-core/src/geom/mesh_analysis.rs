//! Triangle mesh quality and topology analysis.
//! Source: `BRepMesh`, `Poly_MeshPurpose`, mesh validation utilities.

use crate::gp::{GpPnt, GpVec};

/// Triangle quality metrics for a single triangle.
#[derive(Debug, Clone, Copy)]
pub struct TriangleQuality {
    /// Ratio of twice the inradius to the circumradius (equilateral = 1).
    pub aspect_ratio: f64,
    /// Minimum angle in degrees.
    pub min_angle_deg: f64,
    /// Maximum angle in degrees.
    pub max_angle_deg: f64,
    /// Area of the triangle.
    pub area: f64,
    /// Whether the triangle is degenerate (near-zero area).
    pub degenerate: bool,
}

/// Compute quality metrics for one triangle.
pub fn triangle_quality(a: &GpPnt, b: &GpPnt, c: &GpPnt) -> TriangleQuality {
    let ab = b.coord.subtracted(&a.coord);
    let ac = c.coord.subtracted(&a.coord);
    let bc = c.coord.subtracted(&b.coord);
    let la = ab.modulus();
    let lb = ac.modulus();
    let lc = bc.modulus();
    let area = 0.5 * ab.crossed(&ac).modulus();
    let degenerate = area < 1e-15 * (la * lb).max(1.0);
    // Angles via the law of cosines.
    let ang = |x: f64, y: f64, z: f64| {
        let denom = 2.0 * x * y;
        if denom.abs() < 1e-30 {
            std::f64::consts::FRAC_PI_2
        } else {
            ((x * x + y * y - z * z) / denom).clamp(-1.0, 1.0).acos()
        }
    };
    let angle_a = ang(lb, lc, la); // at vertex A (sides AC, AB, BC)
    let angle_b = ang(la, lc, lb);
    let angle_c = ang(la, lb, lc);
    let (min_a, max_a) = (angle_a.min(angle_b).min(angle_c), angle_a.max(angle_b).max(angle_c));
    // Aspect ratio: for a triangle, 1 for equilateral, large for skinny.
    let semiperim = (la + lb + lc) * 0.5;
    let aspect_ratio = if degenerate {
        f64::INFINITY
    } else {
        let inradius = area / semiperim.max(1e-30);
        let circumradius = la * lb * lc / (4.0 * area.max(1e-30));
        if circumradius < 1e-30 {
            f64::INFINITY
        } else {
            (inradius / circumradius * 2.0).recip()
        }
    };
    TriangleQuality {
        aspect_ratio,
        min_angle_deg: min_a.to_degrees(),
        max_angle_deg: max_a.to_degrees(),
        area,
        degenerate,
    }
}

/// Aggregate mesh quality over all triangles.
#[derive(Debug, Clone)]
pub struct MeshQualitySummary {
    pub triangle_count: usize,
    pub degenerate_count: usize,
    pub min_aspect: f64,
    pub avg_aspect: f64,
    pub max_aspect: f64,
    pub min_angle_deg: f64,
    pub max_angle_deg: f64,
    pub total_area: f64,
}

/// Analyze a triangle mesh.
pub fn analyze_mesh(vertices: &[GpPnt], triangles: &[(usize, usize, usize)]) -> MeshQualitySummary {
    let mut deg = 0usize;
    let mut total_area = 0.0;
    let mut min_asp = f64::INFINITY;
    let mut max_asp: f64 = 0.0;
    let mut asp_sum: f64 = 0.0;
    let mut min_ang: f64 = 90.0;
    let mut max_ang: f64 = 0.0;
    for &(i, j, k) in triangles {
        if i >= vertices.len() || j >= vertices.len() || k >= vertices.len() {
            continue;
        }
        let q = triangle_quality(&vertices[i], &vertices[j], &vertices[k]);
        if q.degenerate {
            deg += 1;
        }
        total_area += q.area;
        min_asp = min_asp.min(q.aspect_ratio);
        max_asp = max_asp.max(q.aspect_ratio);
        asp_sum += q.aspect_ratio;
        min_ang = min_ang.min(q.min_angle_deg);
        max_ang = max_ang.max(q.max_angle_deg);
    }
    let n = triangles.len().max(1) as f64;
    MeshQualitySummary {
        triangle_count: triangles.len(),
        degenerate_count: deg,
        min_aspect: min_asp,
        avg_aspect: asp_sum / n,
        max_aspect: max_asp,
        min_angle_deg: min_ang,
        max_angle_deg: max_ang,
        total_area,
    }
}

/// Fraction of degenerate triangles in [0, 1].
pub fn degenerate_fraction(summary: &MeshQualitySummary) -> f64 {
    if summary.triangle_count == 0 {
        0.0
    } else {
        summary.degenerate_count as f64 / summary.triangle_count as f64
    }
}

/// Compute per-vertex averaged normals (area-weighted).
pub fn vertex_normals(vertices: &[GpPnt], triangles: &[(usize, usize, usize)]) -> Vec<GpVec> {
    let mut acc = vec![(0.0f64, 0.0, 0.0); vertices.len()];
    for &(i, j, k) in triangles {
        if i >= vertices.len() || j >= vertices.len() || k >= vertices.len() {
            continue;
        }
        let n = vertices[j]
            .coord
            .subtracted(&vertices[i].coord)
            .crossed(&vertices[k].coord.subtracted(&vertices[i].coord));
        acc[i].0 += n.x;
        acc[i].1 += n.y;
        acc[i].2 += n.z;
        acc[j].0 += n.x;
        acc[j].1 += n.y;
        acc[j].2 += n.z;
        acc[k].0 += n.x;
        acc[k].1 += n.y;
        acc[k].2 += n.z;
    }
    acc.into_iter()
        .map(|(x, y, z)| {
            let m = (x * x + y * y + z * z).sqrt();
            if m > 1e-30 {
                GpVec::new(x / m, y / m, z / m)
            } else {
                GpVec::new(0.0, 0.0, 1.0)
            }
        })
        .collect()
}

/// Whether the mesh forms a closed 2-manifold: every edge used by exactly two
/// triangles (or one for a boundary). Returns (open_edge_count, non_manifold).
pub fn mesh_edge_topology(triangles: &[(usize, usize, usize)]) -> (usize, usize) {
    use std::collections::HashMap;
    let mut counts: HashMap<(usize, usize), usize> = HashMap::new();
    for &(i, j, k) in triangles {
        for (a, b) in [(i, j), (j, k), (k, i)] {
            let key = if a < b { (a, b) } else { (b, a) };
            *counts.entry(key).or_insert(0) += 1;
        }
    }
    let open = counts.values().filter(|&&c| c == 1).count();
    let non_manifold = counts.values().filter(|&&c| c > 2).count();
    (open, non_manifold)
}

/// Signed volume of a closed mesh via the divergence theorem.
pub fn mesh_signed_volume(vertices: &[GpPnt], triangles: &[(usize, usize, usize)]) -> f64 {
    let mut vol = 0.0;
    for &(i, j, k) in triangles {
        if i >= vertices.len() || j >= vertices.len() || k >= vertices.len() {
            continue;
        }
        let a = vertices[i].coord;
        let b = vertices[j].coord;
        let c = vertices[k].coord;
        vol += a.dot(&b.crossed(&c)) / 6.0;
    }
    vol
}

#[cfg(test)]
mod tests {
    use super::*;

    fn box_mesh() -> (Vec<GpPnt>, Vec<(usize, usize, usize)>) {
        // 8 corners, 12 triangles (2 per face), outward CCW.
        let v: Vec<GpPnt> = vec![
            GpPnt::new(0.,0.,0.), GpPnt::new(1.,0.,0.), GpPnt::new(1.,1.,0.), GpPnt::new(0.,1.,0.),
            GpPnt::new(0.,0.,1.), GpPnt::new(1.,0.,1.), GpPnt::new(1.,1.,1.), GpPnt::new(0.,1.,1.),
        ];
        let t = vec![
            (0,3,2),(0,2,1), // -Z
            (4,5,6),(4,6,7), // +Z
            (0,1,5),(0,5,4), // -Y
            (3,7,6),(3,6,2), // +Y
            (0,4,7),(0,7,3), // -X
            (1,2,6),(1,6,5), // +X
        ];
        (v, t)
    }

    #[test]
    fn box_mesh_quality() {
        let (v, t) = box_mesh();
        let s = analyze_mesh(&v, &t);
        assert_eq!(s.triangle_count, 12);
        assert_eq!(s.degenerate_count, 0);
        assert!(s.total_area > 5.9 && s.total_area < 6.1);
        assert!(s.min_angle_deg > 40.0, "min angle {}", s.min_angle_deg);
    }

    #[test]
    fn box_mesh_closed() {
        let (_, t) = box_mesh();
        let (open, nonman) = mesh_edge_topology(&t);
        assert_eq!(open, 0, "closed box");
        assert_eq!(nonman, 0);
    }

    #[test]
    fn box_signed_volume() {
        let (v, t) = box_mesh();
        let vol = mesh_signed_volume(&v, &t);
        assert!(vol > 0.99 && vol < 1.01, "volume {vol}");
    }

    #[test]
    fn normals_unit() {
        let (v, t) = box_mesh();
        let n = vertex_normals(&v, &t);
        assert_eq!(n.len(), 8);
        for nn in &n {
            let m = (nn.x() * nn.x() + nn.y() * nn.y() + nn.z() * nn.z()).sqrt();
            assert!((m - 1.0).abs() < 1e-9);
        }
    }

    #[test]
    fn degenerate_triangle_detected() {
        let q = triangle_quality(&GpPnt::new(0.,0.,0.), &GpPnt::new(1.,0.,0.), &GpPnt::new(2.,0.,0.));
        assert!(q.degenerate);
        assert!(q.area < 1e-12);
    }
}
