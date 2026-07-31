//! Full 3D polygon implementation. Source: `Poly_Polygon3D.hxx`
use crate::gp::GpPnt;

/// 3D polygon with optional per-node normals and parameters.
#[derive(Debug, Clone)]
pub struct Polygon3D {
    pub nodes: Vec<GpPnt>,
    pub normals: Option<Vec<GpPnt>>,
    pub params: Option<Vec<f64>>,
    pub deflection: f64,
    pub is_closed: bool,
}

impl Polygon3D {
    pub fn new(nodes: Vec<GpPnt>) -> Self { Self { nodes, normals: None, params: None, deflection: 0.0, is_closed: false } }
    pub fn with_normals(mut self, n: Vec<GpPnt>) -> Self { self.normals = Some(n); self }
    pub fn with_params(mut self, p: Vec<f64>) -> Self { self.params = Some(p); self }
    pub fn closed(mut self) -> Self { self.is_closed = true; self }

    pub fn nb_nodes(&self) -> usize { self.nodes.len() }
    pub fn node(&self, i: usize) -> &GpPnt { &self.nodes[i] }
    pub fn normal(&self, i: usize) -> Option<&GpPnt> { self.normals.as_ref().and_then(|n| n.get(i)) }
    pub fn parameter(&self, i: usize) -> Option<f64> { self.params.as_ref().and_then(|p| p.get(i).copied()) }

    /// Chord length of polygon.
    pub fn length(&self) -> f64 {
        let n = self.nodes.len(); if n < 2 { return 0.0; }
        let mut total = 0.0;
        for i in 0..n-1 { total += self.nodes[i].coord.subtracted(&self.nodes[i+1].coord).modulus(); }
        if self.is_closed && n > 2 { total += self.nodes[n-1].coord.subtracted(&self.nodes[0].coord).modulus(); }
        total
    }

    /// Centroid of polygon nodes.
    pub fn centroid(&self) -> GpPnt {
        let n = self.nodes.len(); if n == 0 { return GpPnt::zero(); }
        let mut c = crate::gp::GpXyz::zero();
        for p in &self.nodes { c = c.added(&p.coord); }
        GpPnt::from_xyz(&c.divided(n as f64))
    }

    /// Plane normal via Newell's method.
    pub fn plane_normal(&self) -> GpPnt {
        let n = self.nodes.len(); if n < 3 { return GpPnt::zero(); }
        let mut nx = 0.0; let mut ny = 0.0; let mut nz = 0.0;
        for i in 0..n {
            let j = (i+1) % n;
            nx += (self.nodes[i].y() - self.nodes[j].y()) * (self.nodes[i].z() + self.nodes[j].z());
            ny += (self.nodes[i].z() - self.nodes[j].z()) * (self.nodes[i].x() + self.nodes[j].x());
            nz += (self.nodes[i].x() - self.nodes[j].x()) * (self.nodes[i].y() + self.nodes[j].y());
        }
        let len = (nx*nx + ny*ny + nz*nz).sqrt();
        if len < 1e-30 { GpPnt::zero() } else { GpPnt::new(nx/len, ny/len, nz/len) }
    }

    /// Compute bounding box.
    pub fn bounding_box(&self) -> crate::bnd::BndBox {
        let mut b = crate::bnd::BndBox::new();
        for n in &self.nodes { b.add_point(n); }
        b
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn polygon_length() {
        let p = Polygon3D::new(vec![GpPnt::new(0.,0.,0.), GpPnt::new(3.,0.,0.), GpPnt::new(3.,4.,0.)]);
        assert!((p.length() - 7.0).abs() < 1e-14);
    }
    #[test]
    fn polygon_normal() {
        let p = Polygon3D::new(vec![GpPnt::new(0.,0.,0.), GpPnt::new(1.,0.,0.), GpPnt::new(0.,1.,0.)]);
        let n = p.plane_normal();
        assert!((n.z() - 1.0).abs() < 1e-14);
    }
}
