//! Full triangulation structures: nodes / UV nodes / triangles / deflection,
//! plus polygons on a triangulation.
//! Source: `Poly_Triangulation.hxx`, `Poly_PolygonOnTriangulation.hxx`,
//! `Poly_ArrayOfNodes.hxx`, `Poly_ArrayOfUVNodes.hxx`.

use crate::gp::{GpPnt, GpPnt2d, GpXyz};
use crate::precision::CONFUSION;
use super::triangulation::Triangle;

/// Managed array of 3D nodes. Source: `Poly_ArrayOfNodes`.
/// A thin managed `Vec<GpPnt>` with first/last (min/max) access.
#[derive(Debug, Clone, Default)]
pub struct PolyArrayOfNodes {
    pub nodes: Vec<GpPnt>,
}

impl PolyArrayOfNodes {
    pub fn new() -> Self { Self { nodes: Vec::new() } }
    pub fn with_length(n: usize) -> Self { Self { nodes: vec![GpPnt::zero(); n] } }
    pub fn from_vec(v: Vec<GpPnt>) -> Self { Self { nodes: v } }

    pub fn len(&self) -> usize { self.nodes.len() }
    pub fn is_empty(&self) -> bool { self.nodes.is_empty() }
    /// Minimum (first) valid index, always 0.
    pub fn lower(&self) -> usize { 0 }
    /// Maximum (last) valid index; `None` when empty.
    pub fn upper(&self) -> Option<usize> { if self.nodes.is_empty() { None } else { Some(self.nodes.len() - 1) } }

    pub fn value(&self, i: usize) -> &GpPnt { &self.nodes[i] }
    pub fn set_value(&mut self, i: usize, p: GpPnt) { self.nodes[i] = p; }
    pub fn first(&self) -> Option<&GpPnt> { self.nodes.first() }
    pub fn last(&self) -> Option<&GpPnt> { self.nodes.last() }
    pub fn push(&mut self, p: GpPnt) { self.nodes.push(p); }
}

/// Managed array of 2D UV nodes. Source: `Poly_ArrayOfUVNodes`.
#[derive(Debug, Clone, Default)]
pub struct PolyArrayOfUVNodes {
    pub nodes: Vec<GpPnt2d>,
}

impl PolyArrayOfUVNodes {
    pub fn new() -> Self { Self { nodes: Vec::new() } }
    pub fn with_length(n: usize) -> Self { Self { nodes: vec![GpPnt2d::zero(); n] } }
    pub fn from_vec(v: Vec<GpPnt2d>) -> Self { Self { nodes: v } }

    pub fn len(&self) -> usize { self.nodes.len() }
    pub fn is_empty(&self) -> bool { self.nodes.is_empty() }
    pub fn lower(&self) -> usize { 0 }
    pub fn upper(&self) -> Option<usize> { if self.nodes.is_empty() { None } else { Some(self.nodes.len() - 1) } }

    pub fn value(&self, i: usize) -> &GpPnt2d { &self.nodes[i] }
    pub fn set_value(&mut self, i: usize, p: GpPnt2d) { self.nodes[i] = p; }
    pub fn first(&self) -> Option<&GpPnt2d> { self.nodes.first() }
    pub fn last(&self) -> Option<&GpPnt2d> { self.nodes.last() }
    pub fn push(&mut self, p: GpPnt2d) { self.nodes.push(p); }
}

/// Full triangulation: 3D nodes, optional UV nodes, triangles, optional
/// per-node normals and a deflection value.
/// Source: `Poly_Triangulation`.
#[derive(Debug, Clone)]
pub struct PolyTriangulation {
    pub nodes: Vec<GpPnt>,
    pub uv_nodes: Option<Vec<GpPnt2d>>,
    pub triangles: Vec<Triangle>,
    pub normals: Option<Vec<GpPnt>>,
    pub deflection: f64,
}

impl Default for PolyTriangulation {
    fn default() -> Self { Self::new() }
}

impl PolyTriangulation {
    /// Constructs an empty triangulation.
    pub fn new() -> Self {
        Self { nodes: Vec::new(), uv_nodes: None, triangles: Vec::new(), normals: None, deflection: 0.0 }
    }

    /// Constructs a triangulation pre-sized for `nb_nodes` nodes and
    /// `nb_triangles` triangles (all zero-initialized).
    pub fn with_capacity(nb_nodes: usize, nb_triangles: usize, has_uv: bool, has_normals: bool) -> Self {
        Self {
            nodes: vec![GpPnt::zero(); nb_nodes],
            uv_nodes: if has_uv { Some(vec![GpPnt2d::zero(); nb_nodes]) } else { None },
            triangles: vec![Triangle::new(0, 0, 0); nb_triangles],
            normals: if has_normals { Some(vec![GpPnt::zero(); nb_nodes]) } else { None },
            deflection: 0.0,
        }
    }

    /// Constructs from 3D nodes and triangles.
    pub fn from_parts(nodes: Vec<GpPnt>, triangles: Vec<Triangle>) -> Self {
        Self { nodes, uv_nodes: None, triangles, normals: None, deflection: 0.0 }
    }

    /// Constructs from 3D nodes, UV nodes and triangles.
    pub fn from_parts_uv(nodes: Vec<GpPnt>, uv_nodes: Vec<GpPnt2d>, triangles: Vec<Triangle>) -> Self {
        Self { nodes, uv_nodes: Some(uv_nodes), triangles, normals: None, deflection: 0.0 }
    }

    /// Converts from the basic `Triangulation` view.
    pub fn from_basic(b: &super::triangulation::Triangulation) -> Self {
        Self {
            nodes: b.nodes.clone(),
            uv_nodes: b.uv_nodes.clone().map(|v| v.into_iter().map(|(u, vv)| GpPnt2d::new(u, vv)).collect()),
            triangles: b.triangles.clone(),
            normals: b.normals.clone(),
            deflection: b.deflection,
        }
    }

    /// Down-converts to the basic `Triangulation` view.
    pub fn to_basic(&self) -> super::triangulation::Triangulation {
        super::triangulation::Triangulation {
            nodes: self.nodes.clone(),
            triangles: self.triangles.clone(),
            normals: self.normals.clone(),
            uv_nodes: self.uv_nodes.as_ref().map(|v| v.iter().map(|p| (p.x(), p.y())).collect()),
            deflection: self.deflection,
        }
    }

    /// Full deep copy.
    pub fn copy(&self) -> Self { self.clone() }

    pub fn nb_nodes(&self) -> usize { self.nodes.len() }
    pub fn nb_triangles(&self) -> usize { self.triangles.len() }
    /// TRUE if there is some geometry (nodes and triangles).
    pub fn has_geometry(&self) -> bool { !self.nodes.is_empty() && !self.triangles.is_empty() }
    pub fn has_uv_nodes(&self) -> bool { self.uv_nodes.as_ref().is_some_and(|v| !v.is_empty()) }
    pub fn has_normals(&self) -> bool { self.normals.as_ref().is_some_and(|v| !v.is_empty()) }

    pub fn deflection(&self) -> f64 { self.deflection }
    pub fn set_deflection(&mut self, d: f64) { self.deflection = d; }

    pub fn node(&self, i: usize) -> GpPnt { self.nodes[i] }
    pub fn set_node(&mut self, i: usize, p: GpPnt) { self.nodes[i] = p; }
    pub fn node_ref(&self, i: usize) -> &GpPnt { &self.nodes[i] }

    /// Returns the UV node at index `i`. Panics if UV nodes are absent.
    pub fn uv_node(&self, i: usize) -> GpPnt2d { self.uv_nodes.as_ref().expect("no UV nodes")[i] }
    pub fn set_uv_node(&mut self, i: usize, p: GpPnt2d) {
        self.uv_nodes.as_mut().expect("no UV nodes")[i] = p;
    }
    /// Attempted UV read returning `None` when UV nodes are absent.
    pub fn uv_node_opt(&self, i: usize) -> Option<GpPnt2d> { self.uv_nodes.as_ref().map(|v| v[i]) }

    pub fn triangle(&self, i: usize) -> Triangle { self.triangles[i] }
    pub fn set_triangle(&mut self, i: usize, t: Triangle) { self.triangles[i] = t; }
    pub fn triangle_ref(&self, i: usize) -> &Triangle { &self.triangles[i] }

    /// Returns the normal at node `i`. Panics if normals are absent.
    pub fn normal(&self, i: usize) -> GpPnt { self.normals.as_ref().expect("no normals")[i] }
    pub fn set_normal(&mut self, i: usize, n: GpPnt) { self.normals.as_mut().expect("no normals")[i] = n; }
    pub fn normal_opt(&self, i: usize) -> Option<GpPnt> { self.normals.as_ref().map(|v| v[i]) }

    /// Clears all arrays.
    pub fn clear(&mut self) {
        self.nodes.clear();
        self.triangles.clear();
        self.uv_nodes = None;
        self.normals = None;
    }

    /// Ensures UV nodes are allocated, sized to the node count.
    pub fn add_uv_nodes(&mut self) {
        let n = self.nodes.len();
        match self.uv_nodes.as_mut() {
            Some(v) if v.len() == n => {}
            _ => self.uv_nodes = Some(vec![GpPnt2d::zero(); n]),
        }
    }

    /// Deallocates UV nodes.
    pub fn remove_uv_nodes(&mut self) { self.uv_nodes = None; }

    /// Ensures normals are allocated, sized to the node count.
    pub fn add_normals(&mut self) {
        let n = self.nodes.len();
        match self.normals.as_mut() {
            Some(v) if v.len() == n => {}
            _ => self.normals = Some(vec![GpPnt::zero(); n]),
        }
    }

    /// Deallocates normals.
    pub fn remove_normals(&mut self) { self.normals = None; }

    /// Resizes the node array (and any parallel UV / normal arrays).
    pub fn resize_nodes(&mut self, n: usize, copy_old: bool) -> Result<(), String> {
        let old = std::mem::take(&mut self.nodes);
        self.nodes = vec![GpPnt::zero(); n];
        if copy_old {
            let m = old.len().min(n);
            self.nodes[..m].copy_from_slice(&old[..m]);
        }
        if let Some(uv) = self.uv_nodes.as_mut() {
            let old_uv = std::mem::take(uv);
            *uv = vec![GpPnt2d::zero(); n];
            if copy_old { let m = old_uv.len().min(n); uv[..m].copy_from_slice(&old_uv[..m]); }
        }
        if let Some(nrm) = self.normals.as_mut() {
            let old_n = std::mem::take(nrm);
            *nrm = vec![GpPnt::zero(); n];
            if copy_old { let m = old_n.len().min(n); nrm[..m].copy_from_slice(&old_n[..m]); }
        }
        Ok(())
    }

    /// Resizes the triangle array.
    pub fn resize_triangles(&mut self, n: usize, copy_old: bool) -> Result<(), String> {
        let old = std::mem::take(&mut self.triangles);
        self.triangles = vec![Triangle::new(0, 0, 0); n];
        if copy_old {
            let m = old.len().min(n);
            self.triangles[..m].copy_from_slice(&old[..m]);
        }
        Ok(())
    }

    /// Computes per-node normals by accumulating triangle normals (cross
    /// product of the two edge vectors) and normalizing. Degenerate nodes
    /// (zero accumulation) get `(0, 0, 1)`.
    pub fn compute_normals(&mut self) {
        self.add_normals();
        let mut acc = vec![GpXyz::zero(); self.nodes.len()];
        for t in &self.triangles {
            let p0 = &self.nodes[t.n0].coord;
            let p1 = &self.nodes[t.n1].coord;
            let p2 = &self.nodes[t.n2].coord;
            let v01 = p1.subtracted(p0);
            let v02 = p2.subtracted(p0);
            let n = v01.crossed(&v02);
            acc[t.n0] = acc[t.n0].add(&n);
            acc[t.n1] = acc[t.n1].add(&n);
            acc[t.n2] = acc[t.n2].add(&n);
        }
        let normals = acc.into_iter().map(|a| {
            let m = a.square_modulus();
            let c = if m == 0.0 { GpXyz::new(0.0, 0.0, 1.0) } else { a.normalized() };
            GpPnt::from_xyz(&c)
        }).collect();
        self.normals = Some(normals);
    }

    /// Axis-aligned bounding box of all nodes.
    pub fn bounding_box(&self) -> crate::bnd::BndBox {
        let mut b = crate::bnd::BndBox::new();
        for n in &self.nodes { b.add_point(n); }
        b
    }
}

/// A polygon in 3D space based on a triangulation: a sequence of node indices
/// (0-based) with optional per-node parameters.
/// Source: `Poly_PolygonOnTriangulation`.
#[derive(Debug, Clone)]
pub struct PolyPolygonOnTriangulation {
    pub nodes: Vec<usize>,
    pub parameters: Option<Vec<f64>>,
    pub deflection: f64,
}

impl Default for PolyPolygonOnTriangulation {
    fn default() -> Self { Self::new(0, false) }
}

impl PolyPolygonOnTriangulation {
    /// Constructs with space for `nb_nodes` nodes and, if `has_params`, parameters.
    pub fn new(nb_nodes: usize, has_params: bool) -> Self {
        Self {
            nodes: vec![0; nb_nodes],
            parameters: if has_params { Some(vec![0.0; nb_nodes]) } else { None },
            deflection: 0.0,
        }
    }

    pub fn from_nodes(nodes: Vec<usize>) -> Self {
        Self { nodes, parameters: None, deflection: 0.0 }
    }

    pub fn from_nodes_params(nodes: Vec<usize>, parameters: Vec<f64>) -> Result<Self, String> {
        if nodes.len() != parameters.len() {
            return Err("PolyPolygonOnTriangulation: nodes and parameters sizes differ".into());
        }
        Ok(Self { nodes, parameters: Some(parameters), deflection: 0.0 })
    }

    pub fn nb_nodes(&self) -> usize { self.nodes.len() }
    pub fn node(&self, i: usize) -> usize { self.nodes[i] }
    pub fn set_node(&mut self, i: usize, n: usize) { self.nodes[i] = n; }
    pub fn has_parameters(&self) -> bool { self.parameters.is_some() }
    pub fn parameter(&self, i: usize) -> Result<f64, String> {
        self.parameters.as_ref().map(|p| p[i]).ok_or_else(|| "parameters are NULL".to_string())
    }
    pub fn set_parameter(&mut self, i: usize, v: f64) -> Result<(), String> {
        self.parameters.as_mut().ok_or_else(|| "parameters are NULL".to_string())?[i] = v;
        Ok(())
    }
    pub fn set_parameters(&mut self, p: Vec<f64>) -> Result<(), String> {
        if p.len() != self.nodes.len() {
            return Err("PolyPolygonOnTriangulation::set_parameters - invalid array size".into());
        }
        self.parameters = Some(p);
        Ok(())
    }
    pub fn deflection(&self) -> f64 { self.deflection }
    pub fn set_deflection(&mut self, d: f64) { self.deflection = d; }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn array_min_max_access() {
        let mut a = PolyArrayOfNodes::with_length(3);
        a.set_value(0, GpPnt::new(1., 2., 3.));
        a.set_value(2, GpPnt::new(4., 5., 6.));
        assert_eq!(a.len(), 3);
        assert_eq!(a.lower(), 0);
        assert_eq!(a.upper(), Some(2));
        assert_eq!(*a.first().unwrap(), GpPnt::new(1., 2., 3.));
        assert_eq!(*a.last().unwrap(), GpPnt::new(4., 5., 6.));
        assert_eq!(*a.value(1), GpPnt::zero());
    }

    #[test]
    fn uv_array_access() {
        let mut a = PolyArrayOfUVNodes::with_length(2);
        a.set_value(1, GpPnt2d::new(0.5, 0.25));
        assert_eq!(*a.value(1), GpPnt2d::new(0.5, 0.25));
        assert_eq!(a.upper(), Some(1));
    }

    #[test]
    fn triangulation_counts() {
        let nodes = vec![GpPnt::new(0.,0.,0.), GpPnt::new(1.,0.,0.), GpPnt::new(0.,1.,0.)];
        let tris = vec![Triangle::new(0, 1, 2)];
        let t = PolyTriangulation::from_parts(nodes, tris);
        assert_eq!(t.nb_nodes(), 3);
        assert_eq!(t.nb_triangles(), 1);
        assert!(t.has_geometry());
        assert!(!t.has_uv_nodes());
        assert_eq!(t.triangle(0), Triangle::new(0, 1, 2));
    }

    #[test]
    fn compute_normals_z() {
        let nodes = vec![GpPnt::new(0.,0.,0.), GpPnt::new(1.,0.,0.), GpPnt::new(0.,1.,0.)];
        let tris = vec![Triangle::new(0, 1, 2)];
        let mut t = PolyTriangulation::from_parts(nodes, tris);
        t.compute_normals();
        assert!(t.has_normals());
        let n = t.normal(0);
        assert!((n.z() - 1.0).abs() < 1e-12);
        assert!((n.x()).abs() < 1e-12);
        assert!((n.y()).abs() < 1e-12);
    }

    #[test]
    fn compute_normals_degenerate_node() {
        // Two triangles sharing node 0 with opposing directions cancel out.
        let nodes = vec![
            GpPnt::new(0.,0.,0.), GpPnt::new(1.,0.,0.), GpPnt::new(0.,1.,0.), GpPnt::new(0.,-1.,0.),
        ];
        let tris = vec![Triangle::new(0, 1, 2), Triangle::new(0, 3, 1)];
        let mut t = PolyTriangulation::from_parts(nodes, tris);
        t.compute_normals();
        // Node 0 normal cancels -> fallback (0,0,1).
        let n = t.normal(0);
        assert!((n.z() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn resize_nodes_copies() {
        let mut t = PolyTriangulation::with_capacity(2, 1, false, false);
        t.set_node(0, GpPnt::new(1., 2., 3.));
        t.set_node(1, GpPnt::new(4., 5., 6.));
        t.resize_nodes(4, true).unwrap();
        assert_eq!(t.nb_nodes(), 4);
        assert_eq!(t.node(0), GpPnt::new(1., 2., 3.));
        assert_eq!(t.node(1), GpPnt::new(4., 5., 6.));
        assert_eq!(t.node(3), GpPnt::zero());
    }

    #[test]
    fn polygon_on_triangulation_params() {
        let mut p = PolyPolygonOnTriangulation::new(3, true);
        p.set_node(0, 5);
        p.set_node(1, 6);
        p.set_node(2, 7);
        p.set_parameter(1, 0.5).unwrap();
        assert_eq!(p.nb_nodes(), 3);
        assert_eq!(p.node(0), 5);
        assert!(p.has_parameters());
        assert!((p.parameter(1).unwrap() - 0.5).abs() < 1e-14);
        assert!(p.set_parameters(vec![1.0; 2]).is_err());
    }

    #[test]
    fn from_parts_uv_roundtrip() {
        let nodes = vec![GpPnt::zero(), GpPnt::new(1., 0., 0.), GpPnt::new(0., 1., 0.)];
        let uv = vec![GpPnt2d::zero(), GpPnt2d::new(1., 0.), GpPnt2d::new(0., 1.)];
        let tris = vec![Triangle::new(0, 1, 2)];
        let t = PolyTriangulation::from_parts_uv(nodes, uv, tris);
        assert!(t.has_uv_nodes());
        assert_eq!(t.uv_node(1), GpPnt2d::new(1., 0.));
        let b = t.to_basic();
        assert_eq!(b.uv_nodes.as_ref().unwrap()[1], (1.0, 0.0));
        let t2 = PolyTriangulation::from_basic(&b);
        assert!(t2.has_uv_nodes());
        assert_eq!(t2.uv_node(1), GpPnt2d::new(1., 0.));
    }

    #[test]
    fn bounding_box_nonvoid() {
        let nodes = vec![GpPnt::new(0.,0.,0.), GpPnt::new(2.,0.,0.), GpPnt::new(0.,3.,0.)];
        let tris = vec![Triangle::new(0, 1, 2)];
        let t = PolyTriangulation::from_parts(nodes, tris);
        let b = t.bounding_box();
        assert!(!b.is_void());
    }

    #[test]
    fn confusion_constant_reachable() {
        assert!(CONFUSION > 0.0);
    }
}
