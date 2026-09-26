//! TShape — underlying shape data (geometric + topological content).
//! Source: `TopoDS_TShape`
use std::collections::HashMap;
use std::sync::{Arc, RwLock};
use std::fmt;
use occt_geom2d::curve::Curve2d;
use crate::abs::{ShapeType, ShapeFlags};
use crate::shape::TopoShape;
use occt_core::toploc::TopLocLocation;

/// Thread-safe shared handle to a TShape (replaces Handle(TopoDS_TShape)).
pub type HandleTShape = Arc<RwLock<TShape>>;

/// Shared shape data. OCCT's TopoDS_TShape (reference-counted, shared).
/// Use `Arc<TShape>` + `RwLock` for interior mutability like OCCT's Handle(TShape).
pub struct TShape {
    pub shape_type: ShapeType,
    pub flags: ShapeFlags,
    pub location: TopLocLocation,
    /// Process-unique id, assigned at construction. The geometry side-table is
    /// keyed by this shape's heap address; the id lets `Drop` verify that the
    /// entry still belongs to *this* shape, so a stale drop can never erase the
    /// geometry of a later shape that reused the freed address (parallel tests
    /// made that race observable — `groove_cuts_cylinder` flaked ~20% of full
    /// `--lib` runs).
    pub id: u64,
    /// Real children list (OCCT's `TopoDS_TShape::myShapes`, which stores
    /// `TopoDS_Shape` = TShape + Location + Orientation). A wire holds its
    /// edges, a face its wires, a solid its shells, a compound arbitrary
    /// shapes — each with its own orientation.
    pub children: Vec<TopoShape>,
}

/// Source of `TShape::id` (one per construction, process-wide).
fn next_shape_id() -> u64 {
    use std::sync::atomic::{AtomicU64, Ordering};
    static NEXT: AtomicU64 = AtomicU64::new(1);
    NEXT.fetch_add(1, Ordering::Relaxed)
}

impl fmt::Debug for TShape {
    // Do not recurse into children (that would loop on cyclic topology);
    // print type + child count instead.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "TShape({} children={})", self.shape_type.to_str(), self.children.len())
    }
}

impl TShape {
    pub fn new(shape_type: ShapeType) -> Self {
        Self {
            shape_type,
            flags: ShapeFlags::default(),
            location: TopLocLocation::identity(),
            id: next_shape_id(),
            children: Vec::new(),
        }
    }

    pub fn shape_type(&self) -> ShapeType { self.shape_type }
    pub fn free(&self) -> bool { self.flags.free }
    pub fn set_free(&mut self, v: bool) { self.flags.free = v; }
    pub fn closed(&self) -> bool { self.flags.closed }
    pub fn set_closed(&mut self, v: bool) { self.flags.closed = v; }
    pub fn infinite(&self) -> bool { self.flags.infinite }
    pub fn set_infinite(&mut self, v: bool) { self.flags.infinite = v; }
    pub fn modified(&self) -> bool { self.flags.modified }
    pub fn set_modified(&mut self, v: bool) { self.flags.modified = v; }
    pub fn nb_children(&self) -> usize { self.children.len() }
    pub fn set_location(&mut self, l: &TopLocLocation) { self.location = l.clone(); }
    pub fn child(&self, i: usize) -> Option<TopoShape> { self.children.get(i).cloned() }
    pub fn add_child(&mut self, c: TopoShape) { self.children.push(c); }
}

impl Drop for TShape {
    /// Release the geometry side-table entries for this shape when the last
    /// `Arc` handle is dropped. This keeps the process-wide registry from
    /// growing unbounded and, critically, prevents a stale entry (keyed by a
    /// now-reused heap address) from leaking into an unrelated later shape.
    fn drop(&mut self) {
        crate::tgeometry::GeometryRegistry::global()
            .remove_by_ptr(self as *const TShape as usize, self.id);
    }
}

/// Vertex shape data — 3D point (from BRep_TVertex).
#[derive(Debug)]
pub struct VertexShape {
    pub base: TShape,
    pub point: occt_core::gp::GpPnt,
    pub tolerance: f64,
}

/// Per-face 2D pcurves of an edge, and their `BRep_GCurve` ranges.
///
/// Lifted from `GeometryRegistry`'s `EdgeGeom` payload (`tgeometry.rs:64-83`)
/// as part of T-25: the two maps are kept under **one** lock so that a
/// `(curves, ranges)` pair can never be observed half-updated.
///
/// Keyed by `face_key` (`GeometryRegistry::shape_key(face)`).
/// `BRep_Tool::CurveOnSurface` keys a representation by the `Geom_Surface`, so
/// faces sharing one surface share these pcurves: `GeometryRegistry`'s
/// lookup falls back to another face key on the same surface (see
/// `GeometryRegistry::edge_pcurves`, `tgeometry.rs:310-337`).
/// `Debug` is manual because `dyn Curve2d` is not `Debug`.
#[derive(Default)]
pub struct EdgePcurves {
    /// Forward-then-reversed for a seam edge, one entry for a normal edge.
    pub curves: HashMap<usize, Vec<Arc<dyn Curve2d>>>,
    /// `BRep_GCurve` First/Last of the CurveOnSurface representation.
    pub ranges: HashMap<usize, (f64, f64)>,
}

impl fmt::Debug for EdgePcurves {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("EdgePcurves")
            .field("curves", &self.curves.len())
            .field("ranges", &self.ranges.len())
            .finish()
    }
}

/// Edge shape data — curve + parameter range (from BRep_TEdge).
#[derive(Debug)]
pub struct EdgeShape {
    pub base: TShape,
    pub first: f64,
    pub last: f64,
    pub tolerance: f64,
    pub same_parameter: bool,
    pub same_range: bool,
    pub degenerated: bool,
    /// Per-face pcurves (T-25: moved off the `GeometryRegistry` side table).
    pub pcurves: RwLock<EdgePcurves>,
}

/// Wire shape data.
#[derive(Debug)]
pub struct WireShape { pub base: TShape }

/// Face shape data — surface (from BRep_TFace).
#[derive(Debug)]
pub struct FaceShape {
    pub base: TShape,
    pub tolerance: f64,
    pub natural_restriction: bool,
}

/// Shell shape data.
#[derive(Debug)]
pub struct ShellShape { pub base: TShape }

/// Solid shape data.
#[derive(Debug)]
pub struct SolidShape { pub base: TShape }

/// Compound shape data — groups arbitrary sub-shapes.
#[derive(Debug)]
pub struct CompoundShape { pub base: TShape }

impl VertexShape { pub fn new(p: occt_core::gp::GpPnt) -> Self { Self { base: TShape::new(ShapeType::Vertex), point: p, tolerance: 0.0 } } }
impl EdgeShape {
    /// All pcurves of this edge on `face_key` (`GeometryRegistry::edge_pcurves`
    /// direct hit; the surface-identity fallback stays in the registry until
    /// the face surface slot moves too).
    pub fn pcurves_on(&self, face_key: usize) -> Vec<Arc<dyn Curve2d>> {
        self.pcurves.read().unwrap().curves.get(&face_key).cloned().unwrap_or_default()
    }

    /// The first pcurve on `face_key` (`GeometryRegistry::edge_pcurve`).
    pub fn pcurve_on(&self, face_key: usize) -> Option<Arc<dyn Curve2d>> {
        self.pcurves_on(face_key).into_iter().next()
    }

    /// Attach a pcurve (`BRep_Builder::UpdateEdge(edge, c2d, face, tol)`).
    pub fn set_pcurve_on(&self, face_key: usize, curve: Arc<dyn Curve2d>) {
        let mut g = self.pcurves.write().unwrap();
        let v = g.curves.entry(face_key).or_default();
        if v.is_empty() {
            v.push(curve);
        } else {
            v[0] = curve;
        }
    }

    /// Replace the pcurves of `face_key` (seam overload
    /// `BRep_Builder::UpdateEdge(edge, c1, c2, face)`).
    pub fn set_pcurves_on(&self, face_key: usize, curves: Vec<Arc<dyn Curve2d>>) {
        self.pcurves.write().unwrap().curves.insert(face_key, curves);
    }

    /// `BRep_GCurve` range of the curve-on-surface on `face_key`.
    pub fn pcurve_range_on(&self, face_key: usize) -> Option<(f64, f64)> {
        self.pcurves.read().unwrap().ranges.get(&face_key).copied()
    }

    /// `BRep_Builder::Range(edge, face, first, last)`.
    pub fn set_pcurve_range_on(&self, face_key: usize, range: (f64, f64)) {
        self.pcurves.write().unwrap().ranges.insert(face_key, range);
    }

    /// Drop the curve-on-surface representation(s) on `face_key`.
    pub fn remove_pcurves_on(&self, face_key: usize) {
        let mut g = self.pcurves.write().unwrap();
        g.curves.remove(&face_key);
        g.ranges.remove(&face_key);
    }

    pub fn new() -> Self { Self { base: TShape::new(ShapeType::Edge), first: f64::NEG_INFINITY, last: f64::INFINITY, tolerance: 0.0, same_parameter: false, same_range: false, degenerated: false, pcurves: RwLock::new(EdgePcurves::default()) } } }
impl WireShape { pub fn new() -> Self { Self { base: TShape::new(ShapeType::Wire) } } }
impl FaceShape { pub fn new() -> Self { Self { base: TShape::new(ShapeType::Face), tolerance: 0.0, natural_restriction: false } } }
impl ShellShape { pub fn new() -> Self { Self { base: TShape::new(ShapeType::Shell) } } }
impl SolidShape { pub fn new() -> Self { Self { base: TShape::new(ShapeType::Solid) } } }
impl CompoundShape { pub fn new() -> Self { Self { base: TShape::new(ShapeType::Compound) } } }

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn vertex_defaults() {
        let v = VertexShape::new(occt_core::gp::GpPnt::new(1.,2.,3.));
        assert_eq!(v.base.shape_type(), ShapeType::Vertex);
        assert_eq!(v.point.x(), 1.0);
    }

    #[test]
    fn edge_bounds() {
        let e = EdgeShape::new();
        assert!(e.first.is_infinite());
        assert!(!e.degenerated);
    }
}
