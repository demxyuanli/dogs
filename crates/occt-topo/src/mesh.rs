//! Shape → mesh conversion (deflection-based tessellation).
//! Source: `BRepMesh_IncrementalMesh` (simplified wireframe/box meshing).
use crate::shape::{Solid, TopoShape, Face};
use crate::abs::ShapeType;
use occt_core::gp::GpPnt;
use occt_core::bnd::BndBox;
use occt_core::poly::Triangulation;
use occt_core::poly::triangulation::Triangle;

/// Simple mesh representation of a solid (points + triangles + quads).
#[derive(Debug, Clone)]
pub struct ShapeMesh {
    pub vertices: Vec<GpPnt>,
    pub triangles: Vec<Triangle>,
    pub source_shape: ShapeType,
}

/// Mesh a box primitive (axis-aligned) into 12 triangles (6 faces × 2).
/// This produces the exact triangulation of the box surface.
pub fn mesh_box(corners: (GpPnt, GpPnt)) -> ShapeMesh {
    let (lo, hi) = corners;
    // 8 vertices
    let v = [
        lo,
        GpPnt::new(hi.x(), lo.y(), lo.z()),
        GpPnt::new(hi.x(), hi.y(), lo.z()),
        GpPnt::new(lo.x(), hi.y(), lo.z()),
        GpPnt::new(lo.x(), lo.y(), hi.z()),
        GpPnt::new(hi.x(), lo.y(), hi.z()),
        GpPnt::new(hi.x(), hi.y(), hi.z()),
        GpPnt::new(lo.x(), hi.y(), hi.z()),
    ];
    // 6 faces, 2 triangles each (outward normals)
    let tris = [
        (0,1,5),(0,5,4), // -Y
        (1,2,6),(1,6,5), // +X
        (2,3,7),(2,7,6), // +Y
        (3,0,4),(3,4,7), // -X
        (0,3,2),(0,2,1), // -Z
        (4,5,6),(4,6,7), // +Z
    ];
    ShapeMesh {
        vertices: v.to_vec(),
        triangles: tris.iter().map(|&(a,b,c)| Triangle::new(a,b,c)).collect(),
        source_shape: ShapeType::Solid,
    }
}

/// Mesh a sphere (lat/long grid) with given stacks/slices. Radius r at origin.
pub fn mesh_sphere(radius: f64, stacks: usize, slices: usize) -> ShapeMesh {
    let mut vertices = Vec::new();
    let mut triangles = Vec::new();
    for i in 0..=stacks {
        let phi = std::f64::consts::PI * i as f64 / stacks as f64;
        let y = radius * phi.cos();
        let r = radius * phi.sin();
        for j in 0..=slices {
            let theta = 2.0 * std::f64::consts::PI * j as f64 / slices as f64;
            vertices.push(GpPnt::new(r * theta.cos(), y, r * theta.sin()));
        }
    }
    let stride = slices + 1;
    for i in 0..stacks {
        for j in 0..slices {
            let a = i * stride + j;
            let b = a + 1;
            let c = a + stride;
            let d = c + 1;
            triangles.push(Triangle::new(a, b, c));
            triangles.push(Triangle::new(b, d, c));
        }
    }
    ShapeMesh { vertices, triangles, source_shape: ShapeType::Solid }
}

/// Mesh a cylinder (side + caps) with given slices, radius r, height h along Y.
pub fn mesh_cylinder(radius: f64, height: f64, slices: usize) -> ShapeMesh {
    let half = height * 0.5;
    let mut vertices = Vec::new();
    let mut triangles = Vec::new();
    // Side
    for i in 0..=slices {
        let theta = 2.0 * std::f64::consts::PI * i as f64 / slices as f64;
        let x = radius * theta.cos();
        let z = radius * theta.sin();
        vertices.push(GpPnt::new(x, -half, z));
        vertices.push(GpPnt::new(x, half, z));
    }
    let n_side = (slices + 1) * 2;
    for i in 0..slices {
        let a = 2*i; let b = 2*i+1; let c = 2*i+2; let d = 2*i+3;
        triangles.push(Triangle::new(a, b, c));
        triangles.push(Triangle::new(b, d, c));
    }
    // Caps (fan)
    let top_idx = vertices.len();
    let bot_idx = vertices.len() + 1;
    vertices.push(GpPnt::new(0.0, half, 0.0));
    vertices.push(GpPnt::new(0.0, -half, 0.0));
    for i in 0..slices {
        let tb = 2*i+1; let tt = 2*i+3; // top ring
        let bb = 2*i; let bt = 2*i+2;   // bottom ring
        triangles.push(Triangle::new(top_idx, tt, tb)); // top cap
        triangles.push(Triangle::new(bot_idx, bb, bt)); // bottom cap
    }
    ShapeMesh { vertices, triangles, source_shape: ShapeType::Solid }
}

/// Compute a bounding box from a triangulated mesh.
pub fn mesh_bbox(mesh: &ShapeMesh) -> BndBox {
    let mut b = BndBox::new();
    for v in &mesh.vertices { b.add_point(v); }
    b
}

/// Total surface area of a triangle mesh.
pub fn mesh_surface_area(mesh: &ShapeMesh) -> f64 {
    let mut area = 0.0;
    for t in &mesh.triangles {
        let a = &mesh.vertices[t.n0]; let b = &mesh.vertices[t.n1]; let c = &mesh.vertices[t.n2];
        let ab = b.coord.subtracted(&a.coord);
        let ac = c.coord.subtracted(&a.coord);
        area += 0.5 * ab.crossed(&ac).modulus();
    }
    area
}

/// Convert a ShapeMesh to a Poly Triangulation.
pub fn to_triangulation(mesh: &ShapeMesh) -> Triangulation {
    Triangulation::new(mesh.vertices.clone(), mesh.triangles.clone())
}

/// Compute outward-facing vertex normals for a mesh.
pub fn compute_vertex_normals(mesh: &ShapeMesh) -> Vec<GpPnt> {
    let n = mesh.vertices.len();
    let mut accum = vec![occt_core::gp::GpXyz::zero(); n];
    for t in &mesh.triangles {
        let a = &mesh.vertices[t.n0]; let b = &mesh.vertices[t.n1]; let c = &mesh.vertices[t.n2];
        let nrm = b.coord.subtracted(&a.coord).crossed(&c.coord.subtracted(&a.coord));
        accum[t.n0] = accum[t.n0].added(&nrm);
        accum[t.n1] = accum[t.n1].added(&nrm);
        accum[t.n2] = accum[t.n2].added(&nrm);
    }
    accum.into_iter().map(|v| {
        let m = v.modulus();
        if m > 1e-30 { GpPnt::from_xyz(&v.divided(m)) } else { GpPnt::zero() }
    }).collect()
}

/// Mesh quality: fraction of degenerate (zero-area) triangles.
pub fn degenerate_fraction(mesh: &ShapeMesh) -> f64 {
    let total = mesh.triangles.len();
    if total == 0 { return 0.0; }
    let mut deg = 0;
    for t in &mesh.triangles {
        let a = &mesh.vertices[t.n0]; let b = &mesh.vertices[t.n1]; let c = &mesh.vertices[t.n2];
        let ab = b.coord.subtracted(&a.coord);
        let ac = c.coord.subtracted(&a.coord);
        if ab.crossed(&ac).square_modulus() < 1e-24 { deg += 1; }
    }
    deg as f64 / total as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn box_has_12_triangles() {
        let mesh = mesh_box((GpPnt::new(0.,0.,0.), GpPnt::new(1.,1.,1.)));
        assert_eq!(mesh.triangles.len(), 12);
        assert_eq!(mesh.vertices.len(), 8);
        assert!((mesh_surface_area(&mesh) - 6.0).abs() < 1e-12);
    }

    #[test]
    fn sphere_area() {
        let mesh = mesh_sphere(1.0, 20, 20);
        let area = mesh_surface_area(&mesh);
        let expect = 4.0*std::f64::consts::PI;
        assert!((area - expect).abs() < 0.5, "sphere area {area} vs {expect}");
    }

    #[test]
    fn cylinder_side_area() {
        let mesh = mesh_cylinder(1.0, 2.0, 32);
        let area = mesh_surface_area(&mesh);
        // side = 2πrh = 4π, caps = 2π → total 6π
        let expect2 = 6.0*std::f64::consts::PI;
        assert!((area - expect2).abs() < 1.0, "cyl area {area} vs {expect2}");
    }

    #[test]
    fn no_degenerate() {
        let mesh = mesh_box((GpPnt::new(0.,0.,0.), GpPnt::new(1.,1.,1.)));
        assert!(degenerate_fraction(&mesh) < 1e-12);
    }
}
