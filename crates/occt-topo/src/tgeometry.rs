//! Geometry attachment — side-table storing analytic geometry for TShapes.
//!
//! OCCT stores geometry *inside* the TShape sub-classes (`BRep_TVertex`,
//! `BRep_TEdge`, `BRep_TFace`). This port keeps `TShape` itself light
//! (Clone + Debug, no trait objects) and attaches geometry via a process-wide
//! registry keyed by TShape identity (pointer address). `BRepBuilder` writes
//! geometry here; `BRepTool` reads it back. This mirrors `BRep_Tool`'s role of
//! separating topological structure from geometric content.
//!
//! The `TShape` address is stable for the lifetime of the `Arc`, so pointer
//! keys are valid until the shape is dropped. `clear_shape` (called on Drop
//! hooks, or manually) releases the entry.
//!
//! Source: `BRep_TVertex` / `BRep_TEdge` / `BRep_TFace` (TKBRep).

use std::collections::HashMap;
use std::sync::{Arc, RwLock, OnceLock};

use occt_core::gp::GpPnt;
use occt_geom::{Curve, Surface};
use occt_geom2d::curve::Curve2d;

use crate::shape::TopoShape;

/// Vertex geometry — a 3D point and a tolerance. (BRep_TVertex)
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VertexGeom {
    pub point: GpPnt,
    pub tolerance: f64,
}

/// Edge geometry — an underlying curve and a parameter range. (BRep_TEdge)
pub struct EdgeGeom {
    pub curve: Arc<dyn Curve>,
    pub first: f64,
    pub last: f64,
    pub tolerance: f64,
    pub same_parameter: bool,
    pub same_range: bool,
    pub degenerated: bool,
    /// Per-face 2D pcurves, keyed by face pointer identity (see `shape_key`).
    /// Mirrors `BRep_TEdge`'s list of `(face -> Geom2d_Curve)`: a seam edge of a
    /// periodic surface carries *two* pcurves on the same face (one per side of
    /// the seam, e.g. the cone's `u = 0` and `u = 2π` sides), stored in
    /// forward-then-reversed order.
    pub pcurves: HashMap<usize, Vec<Arc<dyn Curve2d>>>,
}

impl EdgeGeom {
    pub fn new(curve: Arc<dyn Curve>, first: f64, last: f64) -> Self {
        Self {
            curve,
            first,
            last,
            tolerance: 0.0,
            same_parameter: true,
            same_range: true,
            degenerated: false,
            pcurves: HashMap::new(),
        }
    }
    pub fn curve(&self) -> Arc<dyn Curve> { self.curve.clone() }
    pub fn parameters(&self) -> (f64, f64) { (self.first, self.last) }

    /// Attach the (single) pcurve of this edge on the face identified by
    /// `face_key`, replacing any previously attached pcurves.
    pub fn set_pcurve(&mut self, face_key: usize, c: Arc<dyn Curve2d>) {
        self.pcurves.insert(face_key, vec![c]);
    }

    /// Replace the pcurves of this edge on the face (one for a normal edge, two
    /// in forward-then-reversed order for a seam edge).
    pub fn set_pcurves(&mut self, face_key: usize, cs: Vec<Arc<dyn Curve2d>>) {
        self.pcurves.insert(face_key, cs);
    }

    /// The first pcurve on the face (the single pcurve of a normal edge).
    pub fn get_pcurve(&self, face_key: usize) -> Option<Arc<dyn Curve2d>> {
        self.pcurves.get(&face_key).and_then(|v| v.first().cloned())
    }

    /// All pcurves on the face (two, in forward-then-reversed order, for a seam
    /// edge).
    pub fn get_pcurves(&self, face_key: usize) -> Vec<Arc<dyn Curve2d>> {
        self.pcurves.get(&face_key).cloned().unwrap_or_default()
    }
}

/// Face geometry — an underlying surface and a tolerance. (BRep_TFace)
pub struct FaceGeom {
    pub surface: Arc<dyn Surface>,
    pub tolerance: f64,
    pub natural_restriction: bool,
}

impl FaceGeom {
    pub fn new(surface: Arc<dyn Surface>) -> Self {
        Self { surface, tolerance: 0.0, natural_restriction: true }
    }
    pub fn surface(&self) -> Arc<dyn Surface> { self.surface.clone() }
}

/// Process-wide geometry side-table keyed by `TShape` pointer identity.
pub struct GeometryRegistry {
    vertices: RwLock<HashMap<usize, VertexGeom>>,
    edges: RwLock<HashMap<usize, EdgeGeom>>,
    faces: RwLock<HashMap<usize, FaceGeom>>,
}

/// Registry key: the address of the `TShape` stored inside the shared
/// `RwLock`. This is stable for the lifetime of the `Arc` (the `TShape` never
/// moves) and matches what `TShape::drop` can recover from `&self`, so entries
/// are reliably removed when the last handle to a shape is dropped.
fn key(s: &TopoShape) -> usize {
    let lock = s.tshape.read().expect("poisoned TShape lock");
    std::ptr::addr_of!(*lock) as usize
}

impl GeometryRegistry {
    /// The shared registry. Shapes and geometry live for the whole process,
    /// matching OCCT's reference-counted Handle model.
    pub fn global() -> &'static GeometryRegistry {
        static REGISTRY: OnceLock<GeometryRegistry> = OnceLock::new();
        REGISTRY.get_or_init(|| GeometryRegistry {
            vertices: RwLock::new(HashMap::new()),
            edges: RwLock::new(HashMap::new()),
            faces: RwLock::new(HashMap::new()),
        })
    }

    // ---- vertices ----

    pub fn set_vertex(&self, s: &TopoShape, geom: VertexGeom) {
        self.vertices.write().unwrap().insert(key(s), geom);
    }

    pub fn vertex_geom(&self, s: &TopoShape) -> Option<VertexGeom> {
        self.vertices.read().unwrap().get(&key(s)).copied()
    }

    /// The vertex point, or the origin when unregistered (placeholder fallback).
    pub fn vertex_point(&self, s: &TopoShape) -> GpPnt {
        self.vertex_geom(s).map(|g| g.point).unwrap_or_else(GpPnt::zero)
    }

    pub fn vertex_tolerance(&self, s: &TopoShape) -> f64 {
        self.vertex_geom(s).map(|g| g.tolerance).unwrap_or(0.0)
    }

    // ---- edges ----

    pub fn set_edge(&self, s: &TopoShape, geom: EdgeGeom) {
        self.edges.write().unwrap().insert(key(s), geom);
    }

    pub fn edge_geom(&self, s: &TopoShape) -> Option<EdgeGeom> {
        self.edges.read().unwrap().get(&key(s)).map(|g| {
            // EdgeGeom is not Clone (Arc<dyn Curve> is Clone, but we rebuild a
            // fresh struct to avoid needing Clone on the whole thing).
            EdgeGeom {
                curve: g.curve.clone(),
                first: g.first,
                last: g.last,
                tolerance: g.tolerance,
                same_parameter: g.same_parameter,
                same_range: g.same_range,
                degenerated: g.degenerated,
                pcurves: g.pcurves.clone(),
            }
        })
    }

    /// The underlying curve handle (clone of the Arc).
    pub fn edge_curve(&self, s: &TopoShape) -> Option<Arc<dyn Curve>> {
        self.edge_geom(s).map(|g| g.curve)
    }

    pub fn edge_parameters(&self, s: &TopoShape) -> (f64, f64) {
        match self.edge_geom(s) {
            Some(g) => (g.first, g.last),
            None => (f64::NEG_INFINITY, f64::INFINITY),
        }
    }

    pub fn edge_tolerance(&self, s: &TopoShape) -> f64 {
        self.edge_geom(s).map(|g| g.tolerance).unwrap_or(0.0)
    }

    pub fn same_parameter(&self, s: &TopoShape) -> bool {
        self.edge_geom(s).map(|g| g.same_parameter).unwrap_or(false)
    }

    pub fn is_degenerated_edge(&self, s: &TopoShape) -> bool {
        self.edge_geom(s).map(|g| g.degenerated).unwrap_or(false)
    }

    // ---- edge p-curves ----

    /// The pcurve of edge `s` on the face identified by `face_key` (see
    /// `shape_key`), if one has been attached. For a seam edge this is the
    /// forward pcurve; use [`GeometryRegistry::edge_pcurves`] to get both.
    pub fn edge_pcurve(&self, s: &TopoShape, face_key: usize) -> Option<Arc<dyn Curve2d>> {
        self.edges.read().unwrap().get(&key(s)).and_then(|g| g.get_pcurve(face_key))
    }

    /// All pcurves of edge `s` on the face identified by `face_key`
    /// (forward-then-reversed for a seam edge, one for a normal edge).
    pub fn edge_pcurves(&self, s: &TopoShape, face_key: usize) -> Vec<Arc<dyn Curve2d>> {
        self.edges.read().unwrap().get(&key(s)).map(|g| g.get_pcurves(face_key)).unwrap_or_default()
    }

    /// Attach a pcurve to edge `s` for the face identified by `face_key`.
    /// Mirrors `BRep_Builder::UpdateEdge(edge, curve2d, face, tol)`.
    pub fn set_edge_pcurve(&self, s: &TopoShape, face_key: usize, curve: Arc<dyn Curve2d>) {
        if let Some(g) = self.edges.write().unwrap().get_mut(&key(s)) {
            g.set_pcurve(face_key, curve);
        }
    }

    /// Replace the pcurves of edge `s` on the face identified by `face_key`.
    /// Mirrors the seam overload `BRep_Builder::UpdateEdge(edge, c1, c2, face)`.
    pub fn set_edge_pcurves(&self, s: &TopoShape, face_key: usize, curves: Vec<Arc<dyn Curve2d>>) {
        if let Some(g) = self.edges.write().unwrap().get_mut(&key(s)) {
            g.set_pcurves(face_key, curves);
        }
    }

    /// Stable pointer-identity key for a shape, usable as a `HashMap` key.
    /// Same value as the internal `key` used to index the side-table.
    pub fn shape_key(s: &TopoShape) -> usize {
        key(s)
    }

    // ---- faces ----

    pub fn set_face(&self, s: &TopoShape, geom: FaceGeom) {
        self.faces.write().unwrap().insert(key(s), geom);
    }

    pub fn face_geom(&self, s: &TopoShape) -> Option<FaceGeom> {
        self.faces.read().unwrap().get(&key(s)).map(|g| FaceGeom {
            surface: g.surface.clone(),
            tolerance: g.tolerance,
            natural_restriction: g.natural_restriction,
        })
    }

    /// The underlying surface handle (clone of the Arc).
    pub fn face_surface(&self, s: &TopoShape) -> Option<Arc<dyn Surface>> {
        self.face_geom(s).map(|g| g.surface)
    }

    pub fn face_tolerance(&self, s: &TopoShape) -> f64 {
        self.face_geom(s).map(|g| g.tolerance).unwrap_or(0.0)
    }

    pub fn natural_restriction(&self, s: &TopoShape) -> bool {
        self.face_geom(s).map(|g| g.natural_restriction).unwrap_or(true)
    }

    // ---- lifecycle ----

    /// Drop all geometry entries owned by `s`. Call when a shape is discarded
    /// to keep the side-table from growing without bound.
    pub fn clear_shape(&self, s: &TopoShape) {
        let k = key(s);
        self.vertices.write().unwrap().remove(&k);
        self.edges.write().unwrap().remove(&k);
        self.faces.write().unwrap().remove(&k);
    }

    /// Remove every geometry entry keyed by the raw `TShape` address. Called
    /// from `TShape::drop` so entries die with their shape — this prevents a
    /// stale entry from leaking into a future shape that reuses the address.
    pub fn remove_by_ptr(&self, ptr: usize) {
        self.vertices.write().unwrap().remove(&ptr);
        self.edges.write().unwrap().remove(&ptr);
        self.faces.write().unwrap().remove(&ptr);
    }

    /// Number of live entries (vertices + edges + faces).
    pub fn len(&self) -> usize {
        self.vertices.read().unwrap().len()
            + self.edges.read().unwrap().len()
            + self.faces.read().unwrap().len()
    }

    pub fn is_empty(&self) -> bool { self.len() == 0 }

    /// Remove every entry (for tests / teardown).
    pub fn clear_all(&self) {
        self.vertices.write().unwrap().clear();
        self.edges.write().unwrap().clear();
        self.faces.write().unwrap().clear();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::abs::ShapeType;
    use occt_geom::GeomLine;

    #[test]
    fn vertex_roundtrip() {
        let reg = GeometryRegistry::global();
        let v = TopoShape::new(ShapeType::Vertex);
        reg.set_vertex(&v, VertexGeom { point: GpPnt::new(1.0, 2.0, 3.0), tolerance: 1e-3 });
        let g = reg.vertex_geom(&v).expect("registered");
        assert!(g.point.is_equal(&GpPnt::new(1.0, 2.0, 3.0)));
        assert_eq!(g.tolerance, 1e-3);
        assert!(reg.vertex_point(&v).is_equal(&GpPnt::new(1.0, 2.0, 3.0)));
        reg.clear_shape(&v);
        assert!(reg.vertex_geom(&v).is_none());
    }

    #[test]
    fn edge_roundtrip_with_curve() {
        use occt_core::gp::{GpAx1, GpDir, GpLin};
        let reg = GeometryRegistry::global();
        let e = TopoShape::new(ShapeType::Edge);
        let lin = GpLin::new(GpAx1::new(GpPnt::zero(), GpDir::new(1.0, 0.0, 0.0).unwrap()));
        let curve = Arc::new(GeomLine::new(lin));
        reg.set_edge(&e, EdgeGeom::new(curve, 0.0, 5.0));
        assert_eq!(reg.edge_parameters(&e), (0.0, 5.0));
        let c = reg.edge_curve(&e).expect("curve present");
        let p = c.d0(2.5);
        assert!((p.x() - 2.5).abs() < 1e-12);
        reg.clear_shape(&e);
        assert_eq!(reg.edge_parameters(&e), (f64::NEG_INFINITY, f64::INFINITY));
    }

    #[test]
    fn unregistered_shapes_fall_back() {
        let reg = GeometryRegistry::global();
        let e = TopoShape::new(ShapeType::Edge);
        // Default: unbounded range, no curve, origin point for vertex.
        assert_eq!(reg.edge_parameters(&e), (f64::NEG_INFINITY, f64::INFINITY));
        assert!(reg.edge_curve(&e).is_none());
        let v = TopoShape::new(ShapeType::Vertex);
        assert!(reg.vertex_point(&v).is_equal(&GpPnt::zero()));
    }
}
