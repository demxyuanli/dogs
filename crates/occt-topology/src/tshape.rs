//! TShape — underlying shape data (geometric + topological content).
//! Source: `TopoDS_TShape`
use crate::abs::{ShapeFlags, ShapeType};
use crate::shape::TopoShape;
use occt_core::poly::Triangulation;
use occt_core::toploc::TopLocLocation;
use occt_geom::{Curve, Surface};
use occt_geom2d::curve::Curve2d;
use std::fmt;
use std::sync::{Arc, RwLock};

mod ordered_slots;
pub use ordered_slots::OrderedSlots;

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
    /// `GeometryRegistry` side table, `tgeometry/mod.rs`). `None` on
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
        write!(
            f,
            "TShape({} children={})",
            self.shape_type.to_str(),
            self.children.len()
        )
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

    pub fn shape_type(&self) -> ShapeType {
        self.shape_type
    }
    pub fn free(&self) -> bool {
        self.flags.free
    }
    pub fn set_free(&mut self, v: bool) {
        self.flags.free = v;
    }
    pub fn is_closed(&self) -> bool {
        self.flags.closed
    }
    pub fn set_closed(&mut self, v: bool) {
        self.flags.closed = v;
    }
    pub fn infinite(&self) -> bool {
        self.flags.infinite
    }
    pub fn set_infinite(&mut self, v: bool) {
        self.flags.infinite = v;
    }
    pub fn modified(&self) -> bool {
        self.flags.modified
    }
    pub fn set_modified(&mut self, v: bool) {
        self.flags.modified = v;
    }
    pub fn nb_children(&self) -> usize {
        self.children.len()
    }
    pub fn set_location(&mut self, l: &TopLocLocation) {
        self.location = l.clone();
    }
    pub fn child(&self, i: usize) -> Option<TopoShape> {
        self.children.get(i).cloned()
    }
    pub fn add_child(&mut self, c: TopoShape) {
        self.children.push(c);
    }

    /// The edge's pcurve store, created on first use (only meaningful for
    /// `ShapeType::Edge`; other shapes are simply never asked).
    pub fn edge_pcurves_mut(&mut self) -> &mut EdgePcurves {
        self.edge_pcurves.get_or_insert_with(EdgePcurves::default)
    }

    /// The edge's pcurve store, if it has one yet.
    pub fn edge_pcurves(&self) -> Option<&EdgePcurves> {
        self.edge_pcurves.as_ref()
    }

    /// The edge's 3D geometry store, created on first use (T-25 batch 2).
    pub fn edge_core_mut(&mut self) -> &mut EdgeGeomCore {
        self.edge_core.get_or_insert_with(EdgeGeomCore::default)
    }

    /// The edge's 3D geometry, if it has one yet.
    pub fn edge_core(&self) -> Option<&EdgeGeomCore> {
        self.edge_core.as_ref()
    }

    /// The vertex's geometry store, created on first use (T-25 batch 3).
    pub fn vertex_core_mut(&mut self) -> &mut VertexGeomCore {
        self.vertex_core.get_or_insert_with(VertexGeomCore::default)
    }

    /// The vertex's geometry, if it has one yet.
    pub fn vertex_core(&self) -> Option<&VertexGeomCore> {
        self.vertex_core.as_ref()
    }

    /// The face's geometry store, created on first use (T-25 batch 4).
    pub fn face_core_mut(&mut self) -> &mut FaceGeomCore {
        self.face_core.get_or_insert_with(FaceGeomCore::default)
    }

    /// The face's geometry, if it has one yet.
    pub fn face_core(&self) -> Option<&FaceGeomCore> {
        self.face_core.as_ref()
    }

    /// `BRep_TFace::Triangulation(Poly_MeshPurpose_NONE)` (`BRep_TFace.cxx:47-53`):
    /// the active triangulation, if any.
    ///
    /// UNPORTED: the `Poly_MeshPurpose` selection (`cxx:54-70`) — this crate's
    /// `Triangulation` carries no purpose bits, so only the active one is
    /// answerable.
    pub fn face_triangulation(&self) -> Option<&Arc<Triangulation>> {
        let core = self.face_core()?;
        core.active_triangulation
            .and_then(|i| core.triangulations.get(i))
    }

    /// `BRep_TFace::Triangulations()` (`BRep_TFace.hxx:112-115`): every
    /// triangulation attached to the face.
    pub fn face_triangulations(&self) -> &[Arc<Triangulation>] {
        match self.face_core() {
            Some(c) => &c.triangulations,
            None => &[],
        }
    }

    /// `BRep_TFace::Triangulation(theTriangulation, theToReset)`
    /// (`BRep_TFace.cxx:75-136`).
    ///
    /// With `to_reset` (or a null triangulation) the list is replaced by the
    /// single input and that one becomes active. Otherwise an input already in
    /// the list is made active; an input not in the list replaces the active
    /// entry in place. The active bit OCCT keeps on `Poly_Triangulation` is the
    /// `active_triangulation` index here.
    pub fn set_face_triangulation(&mut self, tri: Option<Arc<Triangulation>>, to_reset: bool) {
        let core = self.face_core_mut();
        if to_reset || tri.is_none() {
            core.active_triangulation = None;
            core.triangulations.clear();
            if let Some(t) = tri {
                core.triangulations.push(t);
                core.active_triangulation = Some(0);
            }
            return;
        }
        let Some(t) = tri else { return };
        // `cxx:99-118`: already listed -> make it the active one.
        if let Some(i) = core.triangulations.iter().position(|e| Arc::ptr_eq(e, &t)) {
            core.active_triangulation = Some(i);
            return;
        }
        // `cxx:119-135`: replace the active entry, keeping the list length.
        if let Some(a) = core.active_triangulation {
            if let Some(slot) = core.triangulations.get_mut(a) {
                *slot = t;
            }
        }
    }

    /// `BRep_TFace::Triangulations(theTriangulations, theActiveTriangulation)`
    /// (`BRep_TFace.cxx:139-177`): replace the whole list; a null active handle
    /// selects the first entry.
    pub fn set_face_triangulations(
        &mut self,
        mut tris: Vec<Arc<Triangulation>>,
        active: Option<usize>,
    ) {
        let core = self.face_core_mut();
        if tris.is_empty() {
            core.active_triangulation = None;
            core.triangulations.clear();
            return;
        }
        core.active_triangulation = Some(match active {
            Some(i) if i < tris.len() => i,
            _ => 0,
        });
        core.triangulations.clear();
        core.triangulations.append(&mut tris);
    }
}

impl Drop for TShape {
    /// Release the geometry side-table entries for this shape when the last
    /// `Arc` handle is dropped. This keeps the process-wide registry from
    /// growing unbounded and, critically, prevents a stale entry (keyed by a
    /// now-reused heap address) from leaking into an unrelated later shape.
    fn drop(&mut self) {
        // The registry only registers a shape through `GeometryRegistry::set_*`,
        // and each of those also fills the corresponding geometry slot on this
        // same `TShape`. With every slot empty there is no entry to release, so
        // a pure-topology shape does not have to take the process-wide lock.
        if self.edge_pcurves.is_some()
            || self.edge_core.is_some()
            || self.vertex_core.is_some()
            || self.face_core.is_some()
        {
            crate::tgeometry::GeometryRegistry::global()
                .remove_by_ptr(self as *const TShape as usize, self.id);
        }
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
/// `tgeometry/mod.rs`).
///
/// `surface` is optional because a face is built without one
/// (`BRep_Builder::MakeFace(F)` leaves it null in OCCT).
#[derive(Clone)]
pub struct FaceGeomCore {
    pub surface: Option<Arc<dyn Surface>>,
    pub tolerance: f64,
    pub natural_restriction: bool,
    /// `BRep_TFace::myTriangulations` (`BRep_TFace.hxx:139`): every
    /// triangulation attached to the face, in insertion order. A face built by
    /// meshing holds one; a glTF "proxy" face also holds one while carrying no
    /// surface (`BRep_Builder::MakeFace(face, triangulation)`,
    /// `BRep_Builder.cxx:516-527`).
    pub triangulations: Vec<Arc<Triangulation>>,
    /// `BRep_TFace::myActiveTriangulation` (`BRep_TFace.hxx:140`), stored as an
    /// index into `triangulations` because this crate has no handle identity
    /// for `Poly_Triangulation`.
    pub active_triangulation: Option<usize>,
}

impl Default for FaceGeomCore {
    fn default() -> Self {
        // `BRep_TFace::BRep_TFace()` (`BRep_TFace.cxx:29-30`): tolerance
        // `RealEpsilon()`, no surface, `myNaturalRestriction` false; the
        // triangulation list and the active handle start empty.
        Self {
            surface: None,
            tolerance: occt_core::precision::REAL_EPSILON,
            natural_restriction: false,
            triangulations: Vec::new(),
            active_triangulation: None,
        }
    }
}

impl fmt::Debug for FaceGeomCore {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("FaceGeomCore")
            .field("has_surface", &self.surface.is_some())
            .field("tolerance", &self.tolerance)
            .field("natural_restriction", &self.natural_restriction)
            .field("triangulations", &self.triangulations.len())
            .field("active_triangulation", &self.active_triangulation)
            .finish()
    }
}

/// The geometry of a vertex (T-25 batch 3: lifted from `VertexGeom`,
/// `tgeometry/mod.rs`).
#[derive(Clone, PartialEq)]
pub struct VertexGeomCore {
    pub point: occt_core::gp::GpPnt,
    pub tolerance: f64,
    /// `BRep_TVertex::myPoints` (`BRep_TVertex.hxx:96`): the vertex's point
    /// representations, in `NCollection_List` append order. Read by
    /// `BRep_Tool::Parameter(V, E)` (`BRep_Tool.cxx:1401-1492`),
    /// `BRep_Tool::Parameter(V, E, S, L)` (`cxx:1615-1628`) and
    /// `BRep_Tool::Parameters(V, F)` (`cxx:1663-1677`); written by
    /// `BRep_Builder::UpdateVertex` for an `INTERNAL` / `EXTERNAL` vertex
    /// (`BRep_Builder.cxx:1287-1300`) and rebuilt by
    /// `ShapeAnalysis_TransferParametersProj::CopyNMVertex`
    /// (`ShapeAnalysis_TransferParametersProj.cxx:734-781`).
    pub points: Vec<crate::brep_point_representation::PointRepresentation>,
}

impl Default for VertexGeomCore {
    fn default() -> Self {
        // `BRep_TVertex::BRep_TVertex()` (`BRep_TVertex.cxx:27`): tolerance
        // `RealEpsilon()`, empty point representation list.
        Self {
            point: occt_core::gp::GpPnt::zero(),
            tolerance: occt_core::precision::REAL_EPSILON,
            points: Vec::new(),
        }
    }
}

impl VertexGeomCore {
    /// `UpdatePoints(lpr, p, C, L)` (`BRep_Builder.cxx:407-433`): rewrite the
    /// `Parameter` of the first representation that is a point on `c`, or
    /// append a fresh `BRep_PointOnCurve`.
    pub fn update_point_on_curve(&mut self, p: f64, c: &Arc<dyn Curve>) {
        use crate::brep_point_representation::PointRepresentation;
        for pr in self.points.iter_mut() {
            if pr.is_point_on_curve_of(c) {
                pr.set_parameter(p);
                return;
            }
        }
        self.points.push(PointRepresentation::OnCurve {
            curve: Arc::clone(c),
            parameter: p,
        });
    }

    /// `UpdatePoints(lpr, p, PC, S, L)` (`BRep_Builder.cxx:436-463`).
    pub fn update_point_on_curve_on_surface(
        &mut self,
        p: f64,
        pc: &Arc<dyn occt_geom2d::curve::Curve2d>,
        s: &Arc<dyn Surface>,
    ) {
        use crate::brep_point_representation::PointRepresentation;
        for pr in self.points.iter_mut() {
            if pr.is_point_on_curve_on_surface_of(pc, s) {
                pr.set_parameter(p);
                return;
            }
        }
        self.points.push(PointRepresentation::OnCurveOnSurface {
            pcurve: Arc::clone(pc),
            surface: Arc::clone(s),
            parameter: p,
        });
    }

    /// `UpdatePoints(lpr, p1, p2, S, L)` (`BRep_Builder.cxx:466-497`): the
    /// first parameter is written back to the matched representation and the
    /// second through `Parameter2`.
    pub fn update_point_on_surface(&mut self, p1: f64, p2: f64, s: &Arc<dyn Surface>) {
        use crate::brep_point_representation::PointRepresentation;
        for pr in self.points.iter_mut() {
            if pr.is_point_on_surface_of(s) {
                pr.set_parameter(p1);
                pr.set_parameter2(p2);
                return;
            }
        }
        self.points.push(PointRepresentation::OnSurface {
            surface: Arc::clone(s),
            parameter: p1,
            parameter2: p2,
        });
    }
}

/// The 3D geometry of an edge (T-25 batch 2: lifted from `EdgeGeom`,
/// `tgeometry/mod.rs` — `curve` plus the six `BRep_TEdge` scalars).
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
    /// The `BRep_Polygon3D` representation (`BRep_Builder::UpdateEdge(E, P, L)`,
    /// `BRep_Builder.cxx:749-785`): a polygon-backed edge with no 3D curve.
    /// glTF line primitives create exactly this kind of edge —
    /// `aBuilder.MakeEdge(anEdge, aPoly)` (`RWGltf_GltfJsonParser.cxx:1994-1996`,
    /// `.lxx:43-47`). `BRep_Tool::Curve` stays null for such an edge; the
    /// polygon is read through `BRep_Tool::Polygon3D` (`BRep_Tool.cxx:723-741`).
    pub polygon3d: Option<Arc<occt_core::poly::Polygon3D>>,
    /// The `BRep_PolygonOnTriangulation` / `BRep_PolygonOnClosedTriangulation`
    /// representations, in `BRep_TEdge::myCurves` append order
    /// (`BRep_Builder::UpdateEdge(E, P, T, L)`, `BRep_Builder.cxx:790-832`, and
    /// the two-polygon overload `:837-883`). Read through
    /// [`crate::brep_tool::BRepTool::polygon_on_triangulation`]
    /// (`BRep_Tool.cxx:683-712`).
    ///
    /// UNPORTED: nothing writes this slot yet — the only producer in OCCT is
    /// `BRepMesh_ModelPostProcessor::PolygonCommitter::commitPolygons`
    /// (`BRepMesh_ModelPostProcessor.cxx:90-155`) through
    /// `BRepMesh_ShapeTool::UpdateEdge(E, P, T, L)`, which needs the pcurve
    /// node indices `IMeshData_PCurve::Index` that
    /// `BRepMesh_BaseMeshAlgo::initDataStructure` (`BRepMesh_BaseMeshAlgo.cxx:98-108`)
    /// assigns. `meshing::mesh_algo::build_data_structure` leaves them at 0, so
    /// reading this slot answers `None` exactly as OCCT does for an edge that
    /// was never meshed.
    pub polygons_on_triangulation: Vec<PolygonOnTriangulationRepr>,
    /// `BRep_CurveOn2Surfaces` entries in `BRep_TEdge::myCurves`
    /// (`BRep_Builder::Continuity` / `UpdateCurves`, `BRep_Builder.cxx:376-405`).
    /// Matched by surface handles plus locations relative to the edge
    /// (`BRep_CurveOn2Surfaces::IsRegularity`, `.cxx:58-65`). Copied by
    /// `BRep_TEdge::EmptyCopy` (`BRep_TEdge.cxx:116-120`).
    pub regularities: Vec<CurveOn2SurfacesRepr>,
}

/// `BRep_CurveOn2Surfaces` (`BRep_CurveOn2Surfaces.hxx` / `.cxx:29-109`):
/// regularity between two surfaces on an edge. Locations are relative to the
/// edge (`L.Predivided(E.Location())` at write time).
#[derive(Clone)]
pub struct CurveOn2SurfacesRepr {
    pub surface1: Arc<dyn Surface>,
    pub surface2: Arc<dyn Surface>,
    pub location1: TopLocLocation,
    pub location2: TopLocLocation,
    /// `GeomAbs_Shape` ordinal (`GeomAbs_C0 = 0` .. `GeomAbs_CN = 6`).
    pub continuity: u8,
}

impl CurveOn2SurfacesRepr {
    /// `BRep_CurveOn2Surfaces::IsRegularity(S1,S2,L1,L2)` (`.cxx:58-65`).
    pub fn is_regularity(
        &self,
        s1: &Arc<dyn Surface>,
        s2: &Arc<dyn Surface>,
        l1: &TopLocLocation,
        l2: &TopLocLocation,
    ) -> bool {
        (Arc::ptr_eq(&self.surface1, s1)
            && Arc::ptr_eq(&self.surface2, s2)
            && self.location1.is_equal(l1)
            && self.location2.is_equal(l2))
            || (Arc::ptr_eq(&self.surface1, s2)
                && Arc::ptr_eq(&self.surface2, s1)
                && self.location1.is_equal(l2)
                && self.location2.is_equal(l1))
    }
}

/// `BRep_PolygonOnTriangulation` (`BRep_PolygonOnTriangulation.hxx:24-64`)
/// together with the `BRep_PolygonOnClosedTriangulation` subclass
/// (`BRep_PolygonOnClosedTriangulation.hxx:24-50`), which only adds
/// `myPolygon2`.
///
/// OCCT keys the representation by `(Poly_Triangulation, TopLoc_Location)`
/// (`BRep_PolygonOnTriangulation::IsPolygonOnTriangulation`, `.cxx:36-49`);
/// this crate's shapes carry no location on the representation, so the
/// triangulation handle is the whole key. `polygon2` is `Some` for the closed
/// (seam) form only.
#[derive(Clone)]
pub struct PolygonOnTriangulationRepr {
    pub triangulation: Arc<Triangulation>,
    pub polygon: Arc<occt_core::poly::PolyPolygonOnTriangulation>,
    pub polygon2: Option<Arc<occt_core::poly::PolyPolygonOnTriangulation>>,
}

impl fmt::Debug for PolygonOnTriangulationRepr {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("PolygonOnTriangulationRepr")
            .field("nb_nodes", &self.polygon.nb_nodes())
            .field("closed", &self.polygon2.is_some())
            .finish()
    }
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
            polygon3d: None,
            polygons_on_triangulation: Vec::new(),
            regularities: Vec::new(),
        }
    }
}

/// Per-face 2D pcurves of an edge, and their `BRep_GCurve` ranges.
///
/// Lifted from `GeometryRegistry`'s `EdgeGeom` payload (`tgeometry/mod.rs`)
/// as part of T-25: the two maps are kept under **one** lock so that a
/// `(curves, ranges)` pair can never be observed half-updated.
///
/// Keyed by `face_key` (`GeometryRegistry::shape_key(face)`).
/// `BRep_Tool::CurveOnSurface` keys a representation by the `Geom_Surface`, so
/// faces sharing one surface share these pcurves: `GeometryRegistry`'s
/// lookup falls back to another face key on the same surface (see
/// `GeometryRegistry::edge_pcurves`, `tgeometry/registry_pcurves.rs`).
/// `Debug` is manual because `dyn Curve2d` is not `Debug`.
#[derive(Default, Clone)]
pub struct EdgePcurves {
    /// Forward-then-reversed for a seam edge, one entry for a normal edge.
    pub curves: OrderedSlots<Vec<Arc<dyn Curve2d>>>,
    /// `BRep_GCurve` First/Last of the CurveOnSurface representation.
    pub ranges: OrderedSlots<(f64, f64)>,
    /// `BRep_CurveRepresentation::Location` relative to the owning edge.
    pub locations: OrderedSlots<TopLocLocation>,
}

impl EdgePcurves {
    /// Attach the (single) pcurve on `face_key`, replacing any previous ones
    /// (`EdgeGeom::set_pcurve`, `tgeometry/mod.rs`).
    pub fn set_pcurve(&mut self, face_key: usize, c: Arc<dyn Curve2d>, first: f64, last: f64) {
        self.set_pcurve_at(face_key, TopLocLocation::identity(), c, first, last);
    }

    pub fn set_pcurve_at(
        &mut self,
        face_key: usize,
        location: TopLocLocation,
        c: Arc<dyn Curve2d>,
        first: f64,
        last: f64,
    ) {
        self.curves.set(face_key, vec![c.clone()]);
        self.locations.set(face_key, location);
        self.init_range(face_key, Some(c.as_ref()), first, last);
    }

    /// Replace the pcurves on `face_key` (one for a normal edge, two in
    /// forward-then-reversed order for a seam edge).
    pub fn set_pcurves(
        &mut self,
        face_key: usize,
        cs: Vec<Arc<dyn Curve2d>>,
        first: f64,
        last: f64,
    ) {
        self.set_pcurves_at(face_key, TopLocLocation::identity(), cs, first, last);
    }

    pub fn set_pcurves_at(
        &mut self,
        face_key: usize,
        location: TopLocLocation,
        cs: Vec<Arc<dyn Curve2d>>,
        first: f64,
        last: f64,
    ) {
        let head = cs.first().cloned();
        self.curves.set(face_key, cs);
        self.locations.set(face_key, location);
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
    pub fn pcurve(&self, face_key: usize) -> Option<Arc<dyn Curve2d>> {
        self.curves.get(&face_key).and_then(|v| v.first().cloned())
    }

    /// All pcurves on `face_key`.
    pub fn pcurves(&self, face_key: usize) -> Vec<Arc<dyn Curve2d>> {
        self.curves.get(&face_key).cloned().unwrap_or_default()
    }

    pub fn location(&self, face_key: usize) -> TopLocLocation {
        self.locations
            .get(&face_key)
            .cloned()
            .unwrap_or_else(TopLocLocation::identity)
    }
}

impl fmt::Debug for EdgePcurves {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("EdgePcurves")
            .field("curves", &self.curves.len())
            .field("ranges", &self.ranges.len())
            .field("locations", &self.locations.len())
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
pub struct WireShape {
    pub base: TShape,
}

/// Face shape data — surface (from BRep_TFace).
#[derive(Debug)]
pub struct FaceShape {
    pub base: TShape,
    pub tolerance: f64,
    pub natural_restriction: bool,
}

/// Shell shape data.
#[derive(Debug)]
pub struct ShellShape {
    pub base: TShape,
}

/// Solid shape data.
#[derive(Debug)]
pub struct SolidShape {
    pub base: TShape,
}

/// Compound shape data — groups arbitrary sub-shapes.
#[derive(Debug)]
pub struct CompoundShape {
    pub base: TShape,
}

impl VertexShape {
    pub fn new(p: occt_core::gp::GpPnt) -> Self {
        Self {
            base: TShape::new(ShapeType::Vertex),
            point: p,
            tolerance: occt_core::precision::REAL_EPSILON,
        }
    }
}
impl EdgeShape {
    pub fn new() -> Self {
        Self {
            base: TShape::new(ShapeType::Edge),
            first: f64::NEG_INFINITY,
            last: f64::INFINITY,
            tolerance: occt_core::precision::REAL_EPSILON,
            same_parameter: false,
            same_range: false,
            degenerated: false,
        }
    }
}
impl WireShape {
    pub fn new() -> Self {
        Self {
            base: TShape::new(ShapeType::Wire),
        }
    }
}
impl FaceShape {
    pub fn new() -> Self {
        Self {
            base: TShape::new(ShapeType::Face),
            tolerance: occt_core::precision::REAL_EPSILON,
            natural_restriction: false,
        }
    }
}
impl ShellShape {
    pub fn new() -> Self {
        Self {
            base: TShape::new(ShapeType::Shell),
        }
    }
}
impl SolidShape {
    pub fn new() -> Self {
        Self {
            base: TShape::new(ShapeType::Solid),
        }
    }
}
impl CompoundShape {
    pub fn new() -> Self {
        Self {
            base: TShape::new(ShapeType::Compound),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn vertex_defaults() {
        let v = VertexShape::new(occt_core::gp::GpPnt::new(1., 2., 3.));
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
