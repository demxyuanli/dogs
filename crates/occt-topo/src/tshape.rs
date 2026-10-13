//! TShape — underlying shape data (geometric + topological content).
//! Source: `TopoDS_TShape`
use std::collections::HashMap;
use std::ops::Index;
use std::sync::{Arc, RwLock};
use std::fmt;
use occt_geom::{Curve, Surface};
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
    /// The face's geometry (T-25 batch 4). `None` until something writes it.
    pub face_core: Option<FaceGeomCore>,
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
            face_core: None,
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

    /// The face's geometry store, created on first use (T-25 batch 4).
    pub fn face_core_mut(&mut self) -> &mut FaceGeomCore {
        self.face_core.get_or_insert_with(FaceGeomCore::default)
    }

    /// The face's geometry, if it has one yet.
    pub fn face_core(&self) -> Option<&FaceGeomCore> { self.face_core.as_ref() }
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
        self.face_core = None;
    }
}

/// Vertex shape data — 3D point (from BRep_TVertex).
#[derive(Debug)]
pub struct VertexShape {
    pub base: TShape,
    pub point: occt_core::gp::GpPnt,
    pub tolerance: f64,
}

/// The geometry of a face (T-25 batch 4: lifted from `FaceGeom`,
/// `tgeometry.rs:143-148`).
///
/// `surface` is optional because a face is built without one
/// (`BRep_Builder::MakeFace(F)` leaves it null in OCCT).
#[derive(Clone)]
pub struct FaceGeomCore {
    pub surface: Option<Arc<dyn Surface>>,
    pub tolerance: f64,
    pub natural_restriction: bool,
}

impl Default for FaceGeomCore {
    fn default() -> Self {
        // `BRep_TFace::BRep_TFace()` (`BRep_TFace.cxx:29-30`): tolerance
        // `RealEpsilon()`, no surface, `myNaturalRestriction` false.
        Self {
            surface: None,
            tolerance: occt_core::precision::REAL_EPSILON,
            natural_restriction: false,
        }
    }
}

impl fmt::Debug for FaceGeomCore {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("FaceGeomCore")
            .field("has_surface", &self.surface.is_some())
            .field("tolerance", &self.tolerance)
            .field("natural_restriction", &self.natural_restriction)
            .finish()
    }
}

/// The geometry of a vertex (T-25 batch 3: lifted from `VertexGeom`,
/// `tgeometry.rs:27-30`).
#[derive(Clone, Copy, PartialEq)]
pub struct VertexGeomCore {
    pub point: occt_core::gp::GpPnt,
    pub tolerance: f64,
}

impl Default for VertexGeomCore {
    fn default() -> Self {
        // `BRep_TVertex::BRep_TVertex()` (`BRep_TVertex.cxx:27`): tolerance
        // `RealEpsilon()`.
        Self { point: occt_core::gp::GpPnt::zero(), tolerance: occt_core::precision::REAL_EPSILON }
    }
}

/// The 3D geometry of an edge (T-25 batch 2: lifted from `EdgeGeom`,
/// `tgeometry.rs:64-83` — `curve` plus the six `BRep_TEdge` scalars).
///
/// `curve` is optional because a `BRep_TEdge` may carry no 3D curve.
#[derive(Clone)]
pub struct EdgeGeomCore {
    pub curve: Option<Arc<dyn Curve>>,
    /// `true` when the edge's `BRep_Curve3D` representation exists but holds a
    /// *null* curve (`BRep_Builder::UpdateCurves`, `BRep_Builder.cxx:58-95`:
    /// `new BRep_Curve3D(nullCurve, L)`), as installed by
    /// `BRep_Builder::Degenerated(E, true)` (`BRep_Builder.cxx:1073-1085`) or by
    /// `ShapeBuild_Edge::MakeEdge(edge, C2d, S, L, p1, p2)`
    /// (`ShapeBuild_Edge.cxx:851-870`) for an edge that has a pcurve only.
    ///
    /// `BRep_Tool::Curve` then answers a null handle (`BRep_Tool.cxx:172-196`),
    /// so `ShapeAnalysis_Edge::HasCurve3d` is false and
    /// `ShapeFix_Edge::FixAddCurve3d` /`ShapeBuild_Edge::BuildCurve3d`
    /// (`BRepLib.cxx:320-325`) rebuild the 3D curve from the pcurve. `curve`
    /// still carries what `BRepAdaptor_Curve(E, F)` builds from that pcurve
    /// (`Adaptor3d_CurveOnSurface`), because `BRepMesh` tessellates through
    /// `BRepAdaptor_Curve(E, F)` rather than `BRep_Tool::Curve`.
    pub curve3d_null: bool,
    pub first: f64,
    pub last: f64,
    pub tolerance: f64,
    pub same_parameter: bool,
    pub same_range: bool,
    pub degenerated: bool,
    /// `BRep_PointOnCurve` parameters written by `BRep_Builder::UpdateVertex`
    /// (`BRepSweep_Rotation::SetDirectingParameter` / `SetGeneratingParameter`).
    /// Keyed by `GeometryRegistry::shape_key` of the vertex.
    pub vertex_params: Vec<(usize, f64)>,
}

impl fmt::Debug for EdgeGeomCore {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("EdgeGeomCore")
            .field("has_curve", &self.curve.is_some())
            .field("first", &self.first)
            .field("last", &self.last)
            .field("tolerance", &self.tolerance)
            .finish()    }
}

impl Default for EdgeGeomCore {
    fn default() -> Self {
        // `BRep_TEdge::BRep_TEdge()` (`BRep_TEdge.cxx:33`): tolerance
        // `RealEpsilon()`; `SameParameter`/`SameRange` start `true`
        // (`cxx:31-38`), the rest false / unbounded range.
        Self {
            curve: None,
            curve3d_null: false,
            first: f64::NEG_INFINITY,
            last: f64::INFINITY,
            tolerance: occt_core::precision::REAL_EPSILON,
            same_parameter: false,
            same_range: false,
            degenerated: false,
            vertex_params: Vec::new(),
        }
    }
}

/// Keyed slots with the traversal order of an OCCT list.
///
/// `BRep_TEdge::myCurves` is an `NCollection_List` (`BRep_TEdge.hxx:81`) and
/// `BRep_Builder::UpdateCurves` (`BRep_Builder.cxx:104-167`) *removes* the
/// representation it is about to replace and *appends* the new one at the end,
/// so `BRep_Tool::CurveOnSurface(E, ..., Index)` walks the representations in
/// append order. The port needs the same order for the same decisions, which a
/// `HashMap` cannot give (its iteration order is seeded per process): every
/// caller that reads "the first representation" would otherwise see a different
/// one on every run.
///
/// Lookup stays hash-based; only iteration follows `order`.
#[derive(Clone)]
pub struct OrderedSlots<V> {
    order: Vec<usize>,
    map: HashMap<usize, V>,
}

impl<V> Default for OrderedSlots<V> {
    fn default() -> Self {
        Self { order: Vec::new(), map: HashMap::new() }
    }
}

impl<V> OrderedSlots<V> {
    /// `BRep_Builder::UpdateCurves` (`BRep_Builder.cxx:104-167`): the
    /// representation on this key is removed and the new one is **appended**,
    /// so writing a pcurve that is already there moves it to the end of the
    /// list.
    pub fn set(&mut self, key: usize, value: V) -> Option<V> {
        self.order.retain(|k| *k != key);
        self.order.push(key);
        self.map.insert(key, value)
    }

    /// `ShapeBuild_Edge::CopyPCurves` (`ShapeBuild_Edge.cxx:395-399`): a key
    /// that is already present keeps its position in the list.
    pub fn insert(&mut self, key: usize, value: V) -> Option<V> {
        if !self.map.contains_key(&key) {
            self.order.push(key);
        }
        self.map.insert(key, value)
    }

    pub fn get(&self, key: &usize) -> Option<&V> { self.map.get(key) }

    pub fn get_mut(&mut self, key: &usize) -> Option<&mut V> { self.map.get_mut(key) }

    pub fn contains_key(&self, key: &usize) -> bool { self.map.contains_key(key) }

    pub fn remove(&mut self, key: &usize) -> Option<V> {
        let out = self.map.remove(key);
        if out.is_some() {
            self.order.retain(|k| k != key);
        }
        out
    }

    /// Re-key in place (the representation keeps its position in the list).
    pub fn rename_key(&mut self, from: usize, to: usize) {
        if from == to {
            return;
        }
        let Some(value) = self.map.remove(&from) else {
            return;
        };
        match self.order.iter().position(|k| *k == from) {
            Some(pos) => self.order[pos] = to,
            None => self.order.push(to),
        }
        self.map.insert(to, value);
        let Self { order, map } = self;
        let mut seen = std::collections::HashSet::new();
        order.retain(|k| map.contains_key(k) && seen.insert(*k));
    }

    pub fn iter(&self) -> impl Iterator<Item = (&usize, &V)> {
        self.order.iter().filter_map(move |k| self.map.get(k).map(|v| (k, v)))
    }

    pub fn keys(&self) -> impl Iterator<Item = &usize> {
        self.order.iter().filter(|k| self.map.contains_key(k))
    }

    /// Same key order as `self`, with each key remapped (`BRepBuilderAPI_Copy`
    /// re-keying the representations onto the copied surfaces).
    pub fn map_keys(&self, f: impl Fn(usize) -> usize) -> Self
    where
        V: Clone,
    {
        let mut out = Self::default();
        for (k, v) in self.iter() {
            out.insert(f(*k), v.clone());
        }
        out
    }

    /// Run `f` on every value in list order. The `&mut` form of
    /// [`OrderedSlots::iter`]: an iterator of `&mut` cannot borrow the map
    /// repeatedly, so the traversal is driven by the order list instead.
    pub fn for_each_mut(&mut self, mut f: impl FnMut(&usize, &mut V)) {
        let Self { order, map } = self;
        for k in order.iter() {
            if let Some(v) = map.get_mut(k) {
                f(k, v);
            }
        }
    }

    pub fn len(&self) -> usize { self.map.len() }

    pub fn is_empty(&self) -> bool { self.map.is_empty() }
}

impl<V> Index<&usize> for OrderedSlots<V> {
    type Output = V;
    fn index(&self, key: &usize) -> &V {
        self.map.get(key).expect("ordered slot key")
    }
}

impl<'a, V> IntoIterator for &'a OrderedSlots<V> {
    type Item = (&'a usize, &'a V);
    type IntoIter = Box<dyn Iterator<Item = (&'a usize, &'a V)> + 'a>;
    fn into_iter(self) -> Self::IntoIter {
        Box::new(self.iter())
    }
}

impl<V: fmt::Debug> fmt::Debug for OrderedSlots<V> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_map().entries(self.iter()).finish()
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
    pub curves: OrderedSlots<Vec<Arc<dyn Curve2d>>>,
    /// `BRep_GCurve` First/Last of the CurveOnSurface representation.
    pub ranges: OrderedSlots<(f64, f64)>,
}

impl EdgePcurves {
    /// Attach the (single) pcurve on `face_key`, replacing any previous ones
    /// (`EdgeGeom::set_pcurve`, `tgeometry.rs:104-107`).
    pub fn set_pcurve(&mut self, face_key: usize, c: Arc<dyn Curve2d>, first: f64, last: f64) {
        self.curves.set(face_key, vec![c.clone()]);
        self.init_range(face_key, Some(c.as_ref()), first, last);
    }

    /// Replace the pcurves on `face_key` (one for a normal edge, two in
    /// forward-then-reversed order for a seam edge).
    pub fn set_pcurves(&mut self, face_key: usize, cs: Vec<Arc<dyn Curve2d>>, first: f64, last: f64) {
        let head = cs.first().cloned();
        self.curves.set(face_key, cs);
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

impl VertexShape { pub fn new(p: occt_core::gp::GpPnt) -> Self { Self { base: TShape::new(ShapeType::Vertex), point: p, tolerance: occt_core::precision::REAL_EPSILON } } }
impl EdgeShape {
    pub fn new() -> Self { Self { base: TShape::new(ShapeType::Edge), first: f64::NEG_INFINITY, last: f64::INFINITY, tolerance: occt_core::precision::REAL_EPSILON, same_parameter: false, same_range: false, degenerated: false } } }
impl WireShape { pub fn new() -> Self { Self { base: TShape::new(ShapeType::Wire) } } }
impl FaceShape { pub fn new() -> Self { Self { base: TShape::new(ShapeType::Face), tolerance: occt_core::precision::REAL_EPSILON, natural_restriction: false } } }
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
