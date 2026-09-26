//! TShape — underlying shape data (geometric + topological content).
//! Source: `TopoDS_TShape`
use std::collections::HashMap;
use std::sync::{Arc, RwLock};
use std::fmt;
use occt_geom::Curve;
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
    /// Per-face 2D pcurves of an EDGE shape (T-25: lifted from the
    /// `GeometryRegistry` side table, `tgeometry.rs:64-83`). `None` on
    /// non-edge shapes and on edges that never got one.
    ///
    /// `TopoShape` already carries `Arc<RwLock<TShape>>`, so the guard
    /// provides the interior mutability: no extra lock per slot.
    pub edge_pcurves: Option<EdgePcurves>,
    /// The edge's 3D geometry (T-25 batch 2). `None` until something writes it.
    pub edge_core: Option<EdgeGeomCore>,
    /// The vertex's geometry (T-25 batch 3). `None` until something writes it.
    pub vertex_core: Option<VertexGeomCore>,
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
            edge_pcurves: None,
            edge_core: None,
            vertex_core: None,
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

    /// The edge's pcurve store, created on first use (only meaningful for
    /// `ShapeType::Edge`; other shapes are simply never asked).
    pub fn edge_pcurves_mut(&mut self) -> &mut EdgePcurves {
        self.edge_pcurves.get_or_insert_with(EdgePcurves::default)
    }

    /// The edge's pcurve store, if it has one yet.
    pub fn edge_pcurves(&self) -> Option<&EdgePcurves> { self.edge_pcurves.as_ref() }

    /// The edge's 3D geometry store, created on first use (T-25 batch 2).
    pub fn edge_core_mut(&mut self) -> &mut EdgeGeomCore {
        self.edge_core.get_or_insert_with(EdgeGeomCore::default)
    }

    /// The edge's 3D geometry, if it has one yet.
    pub fn edge_core(&self) -> Option<&EdgeGeomCore> { self.edge_core.as_ref() }

    /// The vertex's geometry store, created on first use (T-25 batch 3).
    pub fn vertex_core_mut(&mut self) -> &mut VertexGeomCore {
        self.vertex_core.get_or_insert_with(VertexGeomCore::default)
    }

    /// The vertex's geometry, if it has one yet.
    pub fn vertex_core(&self) -> Option<&VertexGeomCore> { self.vertex_core.as_ref() }
}

impl Drop for TShape {
    /// Release the geometry side-table entries for this shape when the last
    /// `Arc` handle is dropped. This keeps the process-wide registry from
    /// growing unbounded and, critically, prevents a stale entry (keyed by a
    /// now-reused heap address) from leaking into an unrelated later shape.
    fn drop(&mut self) {
        crate::tgeometry::GeometryRegistry::global()
            .remove_by_ptr(self as *const TShape as usize, self.id);
        // T-25: the geometry lives on the shape itself, so it dies with it.
        self.edge_pcurves = None;
        self.edge_core = None;
        self.vertex_core = None;
    }
}

/// Vertex shape data — 3D point (from BRep_TVertex).
#[derive(Debug)]
pub struct VertexShape {
    pub base: TShape,
    pub point: occt_core::gp::GpPnt,
    pub tolerance: f64,
}

/// The geometry of a vertex (T-25 batch 3: lifted from `VertexGeom`,
/// `tgeometry.rs:27-30`).
#[derive(Default, Clone, Copy, PartialEq)]
pub struct VertexGeomCore {
    pub point: occt_core::gp::GpPnt,
    pub tolerance: f64,
}

/// The 3D geometry of an edge (T-25 batch 2: lifted from `EdgeGeom`,
/// `tgeometry.rs:64-83` — `curve` plus the six `BRep_TEdge` scalars).
///
/// `curve` is optional because a `BRep_TEdge` may carry no 3D curve.
#[derive(Default, Clone)]
pub struct EdgeGeomCore {
    pub curve: Option<Arc<dyn Curve>>,
    pub first: f64,
    pub last: f64,
    pub tolerance: f64,
    pub same_parameter: bool,
    pub same_range: bool,
    pub degenerated: bool,
}

impl fmt::Debug for EdgeGeomCore {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("EdgeGeomCore")
            .field("has_curve", &self.curve.is_some())
            .field("first", &self.first)
            .field("last", &self.last)
            .field("tolerance", &self.tolerance)
            .finish()
    }
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
#[derive(Default, Clone)]
pub struct EdgePcurves {
    /// Forward-then-reversed for a seam edge, one entry for a normal edge.
    pub curves: HashMap<usize, Vec<Arc<dyn Curve2d>>>,
    /// `BRep_GCurve` First/Last of the CurveOnSurface representation.
    pub ranges: HashMap<usize, (f64, f64)>,
}

impl EdgePcurves {
    /// Attach the (single) pcurve on `face_key`, replacing any previous ones
    /// (`EdgeGeom::set_pcurve`, `tgeometry.rs:104-107`).
    pub fn set_pcurve(&mut self, face_key: usize, c: Arc<dyn Curve2d>, first: f64, last: f64) {
        self.curves.insert(face_key, vec![c.clone()]);
        self.init_range(face_key, Some(c.as_ref()), first, last);
    }

    /// Replace the pcurves on `face_key` (one for a normal edge, two in
    /// forward-then-reversed order for a seam edge).
    pub fn set_pcurves(&mut self, face_key: usize, cs: Vec<Arc<dyn Curve2d>>, first: f64, last: f64) {
        let head = cs.first().cloned();
        self.curves.insert(face_key, cs);
        self.init_range(face_key, head.as_deref(), first, last);
    }

    /// `UpdateCurves` (`BRep_Builder.cxx:149-164`): the new COS range is the
    /// 3D range when finite, else the pcurve's own `[First, Last]`.
    fn init_range(&mut self, face_key: usize, pc: Option<&dyn Curve2d>, first: f64, last: f64) {
        let (f, l) = if first.is_finite() && last.is_finite() {
            (first, last)
        } else if let Some(c) = pc {
            (c.first_parameter(), c.last_parameter())
        } else {
            return;
        };
        self.ranges.insert(face_key, (f, l));
    }

    /// The first pcurve on `face_key`.
    pub fn get_pcurve(&self, face_key: usize) -> Option<Arc<dyn Curve2d>> {
        self.curves.get(&face_key).and_then(|v| v.first().cloned())
    }

    /// All pcurves on `face_key`.
    pub fn get_pcurves(&self, face_key: usize) -> Vec<Arc<dyn Curve2d>> {
        self.curves.get(&face_key).cloned().unwrap_or_default()
    }
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
    pub fn new() -> Self { Self { base: TShape::new(ShapeType::Edge), first: f64::NEG_INFINITY, last: f64::INFINITY, tolerance: 0.0, same_parameter: false, same_range: false, degenerated: false } } }
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
