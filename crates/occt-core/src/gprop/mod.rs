//! Geometric properties (centroid, area, volume, inertia). Source: `GProp/`
//! Computes global properties of geometric primitives and meshes.
use crate::gp::GpPnt;

/// Mass/area/volume properties for a triangulated mesh.
#[derive(Debug, Clone)]
pub struct GProperties {
    pub mass: f64,
    pub center: GpPnt,
    pub inertia: [f64; 9], // 3x3 inertia matrix (row-major)
}

impl GProperties {
    pub fn new() -> Self { Self { mass:0., center:GpPnt::zero(), inertia:[0.;9] } }

    /// Add contribution from a point mass.
    pub fn add_point_mass(&mut self, p: &GpPnt, mass: f64) {
        if mass.abs() < 1e-30 { return; }
        self.mass += mass;
        let cx = self.center.x() + mass * (p.x() - self.center.x()) / self.mass;
        let cy = self.center.y() + mass * (p.y() - self.center.y()) / self.mass;
        let cz = self.center.z() + mass * (p.z() - self.center.z()) / self.mass;
        self.center = GpPnt::new(cx, cy, cz);
    }

    /// Add contribution from a triangle (uniform density).
    pub fn add_triangle(&mut self, a: &GpPnt, b: &GpPnt, c: &GpPnt, density: f64) {
        let ab = b.coord.subtracted(&a.coord);
        let ac = c.coord.subtracted(&a.coord);
        let area = 0.5 * ab.crossed(&ac).modulus();
        let mass = area * density;
        if mass.abs() < 1e-30 { return; }
        let centroid = GpPnt::from_xyz(&a.coord.added(&ab.added(&ac)).divided(3.0));
        self.add_point_mass(&centroid, mass);

        // Inertia contribution (parallel axis theorem)
        let idx = |p: &GpPnt, i, j| -> f64 { let v=[p.x(),p.y(),p.z()]; (if i==j {1.0} else {0.0}) * (v[0]*v[0]+v[1]*v[1]+v[2]*v[2]) - v[i]*v[j] };
        for i in 0..3 { for j in 0..3 {
            let ij = idx(a,i,j) + idx(b,i,j) + idx(c,i,j);
            self.inertia[i*3+j] += mass * ij / 6.0;
        }}
    }
}

/// Compute centroid of a point set.
pub fn centroid_of_points(points: &[GpPnt]) -> GpPnt {
    let n = points.len();
    if n == 0 { return GpPnt::zero(); }
    let mut c = GpXyz::zero();
    for p in points { c = c.added(&p.coord); }
    GpPnt::from_xyz(&c.divided(n as f64))
}

use crate::gp::GpXyz;

/// Compute area of a polygon in 3D (via Newell's method).
pub fn polygon_area(vertices: &[GpPnt]) -> f64 {
    let n = vertices.len();
    if n < 3 { return 0.0; }
    let mut normal = GpXyz::zero();
    for i in 0..n {
        let j = (i + 1) % n;
        normal.x += (vertices[i].y() - vertices[j].y()) * (vertices[i].z() + vertices[j].z());
        normal.y += (vertices[i].z() - vertices[j].z()) * (vertices[i].x() + vertices[j].x());
        normal.z += (vertices[i].x() - vertices[j].x()) * (vertices[i].y() + vertices[j].y());
    }
    0.5 * normal.modulus()
}

/// Signed volume of a closed triangulated mesh (tetrahedron sum).
pub fn mesh_volume(vertices: &[GpPnt], triangles: &[(usize, usize, usize)]) -> f64 {
    let mut vol = 0.0;
    for &(i, j, k) in triangles {
        let a = &vertices[i]; let b = &vertices[j]; let c = &vertices[k];
        vol += a.coord.dot_cross(&b.coord, &c.coord);
    }
    vol.abs() / 6.0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn centroid_square() {
        let pts = vec![GpPnt::new(0.,0.,0.),GpPnt::new(1.,0.,0.),GpPnt::new(1.,1.,0.),GpPnt::new(0.,1.,0.)];
        let c = centroid_of_points(&pts);
        assert!((c.x()-0.5).abs()<1e-14);
        assert!((c.y()-0.5).abs()<1e-14);
    }

    #[test]
    fn triangle_volume() {
        let v = vec![GpPnt::new(0.,0.,0.),GpPnt::new(1.,0.,0.),GpPnt::new(0.,1.,0.),GpPnt::new(0.,0.,1.)];
        let vol = mesh_volume(&v, &[(0,1,2),(0,1,3),(0,2,3),(1,2,3)]);
        assert!((vol - 1.0/6.0).abs() < 1e-14);
    }
}
pub mod inertia;
pub mod gprops;
