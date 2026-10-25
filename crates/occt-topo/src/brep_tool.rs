//! Geometric access to topology — a port of `BRep_Tool`.
//!
//! Reads real geometry (curves, surfaces, tolerances) back from the
//! `GeometryRegistry` side-table that `TopoBuilder` populates. Unregistered
//! shapes fall back to OCCT-like defaults (origin point, unbounded ranges,
//! no curve/surface).

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpPnt2d};
use occt_core::poly::Triangulation;
use occt_geom::{Curve, Surface};

use crate::abs::{Orientation, ShapeType};
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{explorer_shapes, is_same};

/// Namespace for BRep_Tool-style static queries.
pub struct BRepTool;

impl BRepTool {
    /// `BRep_Tool::IsClosed(E, F)` (`BRep_Tool.cxx:795-805`): true when the
    /// edge has a CurveOnSurface representation with two pcurves on the face's
    /// surface (`BRep_Tool.cxx:814-841`).
    ///
    /// UNPORTED: the triangulation arm (`BRep_Tool.cxx:803-804`, `:849-...`);
    /// the port does not attach a `Poly_Triangulation` to faces for this query.
    pub fn is_closed_edge_face(edge: &Edge, face: &Face) -> bool {
        let Some(surf) = GeometryRegistry::global().face_surface(&face.0) else {
            return false;
        };
        // `BRep_Tool.cxx:819-822`: a plane is never a closed surface.
        if surf.gp_pln().is_some() {
            return false;
        }
        let face_key = GeometryRegistry::shape_key(&face.0);
        GeometryRegistry::global().edge_pcurves(&edge.0, face_key).len() > 1
    }

    /// `BRep_Tool::IsClosed(theShape)` (`BRep_Tool.cxx:1707-1757`).
    ///
    /// - `SHELL` (`cxx:1709-1729`): fold every non-degenerated, non
    ///   `INTERNAL`/`EXTERNAL` edge occurrence in a `TopTools_ShapeMapHasher`
    ///   map (add when absent, remove when present); closed when at least one
    ///   edge was counted and the map ends empty.
    /// - `WIRE` (`cxx:1730-1748`): same fold over the vertex occurrences.
    /// - `EDGE` (`cxx:1750-1756`): `TopExp::Vertices(E, V1, V2)` with
    ///   `V1.IsSame(V2)`.
    /// - Any other type (`cxx:1757`): the shape's own `Closed()` flag.
    ///
    /// Each arm starts from `theShape.Oriented(TopAbs_FORWARD)`, so the
    /// explored child orientations are cumulated with the root's forward
    /// orientation; the fold key is `IsSame` (`TShape` + location), which is
    /// what `NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>` tests.
    ///
    /// UNPORTED: nothing in the queried arms; the port has no `Poly_Triangulation`
    /// arm here (that belongs to `IsClosed(E, F)`).
    pub fn is_closed(shape: &TopoShape) -> bool {
        match shape.shape_type() {
            ShapeType::Shell => Self::fold_occurrences_closed(shape, ShapeType::Edge, |e| {
                !Self::is_degenerated(&Edge(e.clone()))
                    && !matches!(e.orientation(), Orientation::Internal | Orientation::External)
            }),
            ShapeType::Wire => Self::fold_occurrences_closed(shape, ShapeType::Vertex, |v| {
                !matches!(v.orientation(), Orientation::Internal | Orientation::External)
            }),
            ShapeType::Edge => match crate::topo_tools::edge_vertices(&Edge(shape.clone())) {
                (Some(a), Some(b)) => is_same(&a.0, &b.0),
                _ => false,
            },
            _ => shape.is_closed(),
        }
    }

    /// The fold shared by the `SHELL` and `WIRE` arms of [`Self::is_closed`]:
    /// walk the forward-oriented shape's `to_find` occurrences, keep those
    /// `accept` lets through, and toggle each in a `TopTools_ShapeMapHasher`
    /// map; closed when `hasBound && map.IsEmpty()`.
    fn fold_occurrences_closed(
        shape: &TopoShape,
        to_find: ShapeType,
        accept: impl Fn(&TopoShape) -> bool,
    ) -> bool {
        let mut root = shape.clone();
        root.set_orientation(Orientation::Forward);
        let mut folded: Vec<TopoShape> = Vec::new();
        let mut has_bound = false;
        for s in explorer_shapes(&root, to_find) {
            if !accept(&s) {
                continue;
            }
            has_bound = true;
            match folded.iter().position(|f| is_same(f, &s)) {
                Some(i) => {
                    folded.remove(i);
                }
                None => folded.push(s),
            }
        }
        has_bound && folded.is_empty()
    }

    /// The location of a vertex (from `BRep_TVertex`), or the origin when
    /// the vertex has no registered geometry.
    pub fn vertex_point(v: &Vertex) -> GpPnt {
        GeometryRegistry::global().vertex_point(&v.0)
    }

    /// `BRep_Tool::Triangulation(theFace, theLocation, NONE)`
    /// (`BRep_Tool.cxx:111-119`): the face's active triangulation, if any.
    ///
    /// OCCT also returns the face location through the out parameter
    /// (`cxx:116`); the port's face location lives on the `TopoShape`, and no
    /// caller of this accessor needs it, so only the triangulation is returned.
    ///
    /// UNPORTED: the `Poly_MeshPurpose` argument (`cxx:115`) — see
    /// [`crate::tshape::TShape::face_triangulation`].
    pub fn face_triangulation(face: &Face) -> Option<Arc<Triangulation>> {
        let ts = face.0.tshape.read().unwrap();
        ts.face_triangulation().cloned()
    }

    /// `BRep_Tool::Triangulations(theFace, theLocation)`
    /// (`BRep_Tool.cxx:123-130`): every triangulation attached to the face.
    pub fn face_triangulations(face: &Face) -> Vec<Arc<Triangulation>> {
        let ts = face.0.tshape.read().unwrap();
        ts.face_triangulations().to_vec()
    }

    /// `BRep_Tool::Polygon3D(theEdge, theLocation)`
    /// (`BRep_Tool.cxx:723-741`): the polygon of a polygon-backed edge, i.e. the
    /// `BRep_Polygon3D` representation `BRep_Builder::UpdateEdge(E, P, L)`
    /// installs. `None` for an edge that carries only a curve.
    pub fn edge_polygon3d(edge: &Edge) -> Option<Arc<occt_core::poly::Polygon3D>> {
        let ts = edge.0.tshape.read().unwrap();
        ts.edge_core().and_then(|c| c.polygon3d.clone())
    }

    /// Version of Tolerance with a minimal value -- a port of
    /// `BRep_Tool::Tolerance(Vertex)` (`BRep_Tool.cxx:1314-1332`):
    ///
    /// ```text
    /// double p = aTVert->Tolerance();      // BRep_TVertex::Tolerance()
    /// constexpr double pMin = Precision::Confusion();
    /// if (p > pMin) return p; else return pMin;
    /// ```
    ///
    /// The clamp is load bearing: a `BRep_TVertex` starts at `RealEpsilon()`
    /// (`BRep_TVertex.cxx:27`), so the raw value is below `Confusion` and the
    /// clamp is what makes `BRep_Tool::Tolerance` answer `Confusion` for an
    /// untouched vertex. `BRepTools_WireExplorer::Init`
    /// (`BRepTools_WireExplorer.cxx:148-164`) treats a wire whose vertex *and*
    /// edge tolerances are all below `Confusion` as empty, which can only
    /// happen through the raw `BRep_TVertex` accessor, never through this one.
    pub fn vertex_tolerance(v: &Vertex) -> f64 {
        clamp_to_confusion(GeometryRegistry::global().vertex_tolerance(&v.0))
    }

    /// The curve underlying an edge (`Handle(Geom_Curve)`), if registered.
    pub fn edge_curve(e: &Edge) -> Option<Arc<dyn Curve>> {
        GeometryRegistry::global().edge_curve(&e.0)
    }

    /// Parameter range of an edge. Unregistered edges report the unbounded
    /// range `(-inf, +inf)`.
    pub fn edge_parameters(e: &Edge) -> (f64, f64) {
        GeometryRegistry::global().edge_parameters(&e.0)
    }

    /// Edge tolerance -- `BRep_Tool::Tolerance(Edge)` (`BRep_Tool.cxx:881-894`),
    /// which returns `BRep_TEdge::Tolerance()` raised to `Precision::Confusion()`
    /// when the stored value is below it.
    pub fn edge_tolerance(e: &Edge) -> f64 {
        clamp_to_confusion(GeometryRegistry::global().edge_tolerance(&e.0))
    }

    /// The surface underlying a face (`Handle(Geom_Surface)`), if registered.
    pub fn face_surface(f: &Face) -> Option<Arc<dyn Surface>> {
        GeometryRegistry::global().face_surface(&f.0)
    }

    /// Face tolerance -- `BRep_Tool::Tolerance(Face)` (`BRep_Tool.cxx:137-150`),
    /// `BRep_TFace::Tolerance()` raised to `Precision::Confusion()`.
    pub fn face_tolerance(f: &Face) -> f64 {
        clamp_to_confusion(GeometryRegistry::global().face_tolerance(&f.0))
    }

    /// Endpoint points of the edge: the curve evaluated at the first and last
    /// parameters. `None` when the edge has no registered curve.
    pub fn edge_vertices(e: &Edge) -> Option<(GpPnt, GpPnt)> {
        let curve = Self::edge_curve(e)?;
        let (first, last) = Self::edge_parameters(e);
        Some((curve.d0(first), curve.d0(last)))
    }

    /// Whether the edge is same-parameter (`BRep_TEdge` flag).
    pub fn same_parameter(e: &Edge) -> bool {
        GeometryRegistry::global().same_parameter(&e.0)
    }

    /// Whether the edge is degenerated (`BRep_TEdge` flag).
    pub fn is_degenerated(e: &Edge) -> bool {
        GeometryRegistry::global().is_degenerated_edge(&e.0)
    }

    /// The first (`i == 0`) or last parameter of the edge.
    pub fn parameter_on_edge(e: &Edge, i: usize) -> f64 {
        let (first, last) = Self::edge_parameters(e);
        if i == 0 {
            first
        } else {
            last
        }
    }

    /// A closed edge has its first parameter equal to its last within tolerance.
    pub fn is_closed_edge(e: &Edge) -> bool {
        let (first, last) = Self::edge_parameters(e);
        (first - last).abs() <= Self::edge_tolerance(e)
    }

    /// Whether the face is bounded by its natural parametric curves
    /// (`BRep_TFace` flag).
    pub fn natural_restriction(f: &Face) -> bool {
        GeometryRegistry::global().natural_restriction(&f.0)
    }

    /// Parametric bounds of the face: `(u_min, u_max, v_min, v_max)` from the
    /// surface's `u_range`/`v_range`. Unbounded fallback for unregistered faces.
    pub fn uv_bounds(f: &Face) -> (f64, f64, f64, f64) {
        match Self::face_surface(f) {
            Some(s) => {
                let (u1, u2) = s.u_range();
                let (v1, v2) = s.v_range();
                (u1, u2, v1, v2)
            }
            None => (f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY),
        }
    }

    /// World-coordinate vertex point: the registered local point transformed
    /// by the shape's accumulated `location` (`TopoDS_Shape::Location`).
    pub fn vertex_point_world(v: &Vertex) -> GpPnt {
        let t = v.0.location().transformation();
        Self::vertex_point(v).transformed(&t)
    }

    /// World-coordinate edge curve: the registered local curve transformed by
    /// the edge's location.
    pub fn edge_curve_world(e: &Edge) -> Option<Arc<dyn Curve>> {
        let t = e.0.location().transformation();
        Self::edge_curve(e).map(|c| Arc::from(c.transformed(&t)))
    }

    /// World-coordinate edge curve as `BRepAdaptor_Curve(E)` evaluates it
    /// (`BRepAdaptor_Curve.cxx`: `myCurve` materialized from
    /// `BRep_Tool::Curve(E, L)` when the edge has a `BRep_Curve3D`, otherwise
    /// the `Adaptor3d_CurveOnSurface` built from the pcurve, and `myTrsf = L`
    /// applied to every returned point). `edge_curve_world` cannot serve here:
    /// it answers `None` on a pcurve-only edge (`curve3d_null`), which the
    /// adaptor still evaluates.
    pub fn edge_curve_adaptor_world(e: &Edge) -> Option<Arc<dyn Curve>> {
        let geom = GeometryRegistry::global().edge_geom(&e.0)?;
        let loc = e.0.location();
        if loc.is_identity() {
            return Some(geom.curve);
        }
        Some(Arc::from(geom.curve.transformed(&loc.transformation())))
    }

    /// World-coordinate face surface: the registered local surface transformed
    /// by the face's location.
    pub fn face_surface_world(f: &Face) -> Option<Arc<dyn Surface>> {
        let t = f.0.location().transformation();
        Self::face_surface(f).map(|s| Arc::from(s.transformed(&t)))
    }

    /// World-coordinate face UV bounds (from the transformed surface).
    pub fn uv_bounds_world(f: &Face) -> (f64, f64, f64, f64) {
        match Self::face_surface_world(f) {
            Some(s) => {
                let (u1, u2) = s.u_range();
                let (v1, v2) = s.v_range();
                (u1, u2, v1, v2)
            }
            None => (f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY),
        }
    }

    /// Structural integrity of a shape list. Always true in this port:
    /// TShape types are fixed at construction.
    pub fn check_shape_integrity(shapes: &[TopoShape]) -> bool {
        let _ = shapes;
        true
    }

    /// Human-readable type name.
    pub fn shape_type_str(t: ShapeType) -> &'static str {
        t.to_str()
    }

    /// `BRep_Tool::Parameter(V, E)` (`BRep_Tool.cxx:1340-1397`).
    ///
    /// The returned parameter is the edge's `First` or `Last` **chosen by the
    /// orientation the vertex carries inside the edge** (`cxx:1360-1383` walks
    /// the edge's vertex children and matches with `IsSame`), flipped when the
    /// edge handle itself is `REVERSED` (`cxx:1368`, `:1390`, `:1396`). The
    /// caller's own vertex orientation only matters for a degenerated edge that
    /// stores no vertex children (`cxx:1354-1357`).
    pub fn parameter_of_vertex(edge: &TopoShape, vertex: &TopoShape) -> f64 {
        let key = GeometryRegistry::shape_key(vertex);
        let (first, last) = Self::edge_parameters(&Edge(edge.clone()));
        let mut found: Option<Orientation> = None;
        let mut rev = false;
        let mut has_vertex_children = false;
        if let Ok(ts) = edge.tshape.read() {
            if let Some(core) = ts.edge_core() {
                if let Some((_, p)) = core.vertex_params.iter().find(|(k, _)| *k == key) {
                    return *p;
                }
            }
            for child in &ts.children {
                if child.shape_type() != ShapeType::Vertex {
                    continue;
                }
                has_vertex_children = true;
                if !is_same(child, vertex) {
                    continue;
                }
                match found {
                    None => found = Some(child.orientation()),
                    Some(cur) => {
                        rev = edge.orientation() == Orientation::Reversed;
                        if child.orientation() == vertex.orientation() {
                            found = Some(cur);
                        }
                    }
                }
            }
        }
        let orient = if let Some(o) = found {
            o
        } else if !has_vertex_children && Self::is_degenerated(&Edge(edge.clone())) {
            vertex.orientation()
        } else {
            // UNPORTED: `cxx:1401-1492` walks the vertex's point representations
            // (`BRep_TVertex::Points()`) against the 3D curve and then against
            // the p-curve, and returns `false` (`cxx:1494`) when none matches.
            // The port's `TopoDS_Vertex` stores no point list, so the branch
            // answers the edge's first parameter instead of failing.
            return first;
        };
        match orient {
            Orientation::Forward => {
                if rev {
                    last
                } else {
                    first
                }
            }
            Orientation::Reversed => {
                if rev {
                    first
                } else {
                    last
                }
            }
            _ => first,
        }
    }

    /// `BRep_Tool::Parameter(V, E, F)` (`BRep_Tool.cxx:1519-1524`) ->
    /// `Parameter(V, E, S, L)` (`cxx:1532-1620`): the parameter of `v` on the
    /// p-curve of `e` on `f`.
    ///
    /// The vertex-in-edge search and the `rev` latch are the `BRep_Tool::
    /// Parameter(V, E)` ones (`cxx:1537-1566`); the difference is the range the
    /// FORWARD / REVERSED arms read — `BRep_Tool::Range(E, S, L, f, l)`
    /// (`cxx:1572-1582`), i.e. the **p-curve** range, not the 3D one. The
    /// INTERNAL arm's `BRep_PointRepresentation` walk (`cxx:1585-1597`) is
    /// UNPORTED for the same reason as in [`Self::parameter_of_vertex`] and
    /// falls back to the 3D parameter.
    pub fn parameter_of_vertex_on_face(vertex: &TopoShape, edge: &Edge, face: &Face) -> f64 {
        let (first, last) = Self::pcurve_range(edge, face);

        let key = GeometryRegistry::shape_key(vertex);
        let mut found: Option<Orientation> = None;
        let mut rev = false;
        if let Ok(ts) = edge.0.tshape.read() {
            if let Some(core) = ts.edge_core() {
                if let Some((_, p)) = core.vertex_params.iter().find(|(k, _)| *k == key) {
                    return *p;
                }
            }
            for child in &ts.children {
                if child.shape_type() != ShapeType::Vertex || !is_same(child, vertex) {
                    continue;
                }
                match found {
                    None => found = Some(child.orientation()),
                    Some(cur) => {
                        rev = edge.0.orientation() == Orientation::Reversed;
                        if child.orientation() == vertex.orientation() {
                            found = Some(cur);
                        }
                    }
                }
            }
        }

        match found {
            Some(Orientation::Forward) => {
                if rev {
                    last
                } else {
                    first
                }
            }
            Some(Orientation::Reversed) => {
                if rev {
                    first
                } else {
                    last
                }
            }
            // `cxx:1585-1597`: the vertex's point representations; not stored by
            // the port (see `parameter_of_vertex`).
            _ => Self::parameter_of_vertex(&edge.0, vertex),
        }
    }

    /// `BRep_Tool::Range(E, S, L, f, l)` (`BRep_Tool.cxx:1004-1013`): the
    /// stored p-curve's range, or the 3D edge range when no p-curve exists.
    fn pcurve_range(e: &Edge, f: &Face) -> (f64, f64) {
        let face_key = GeometryRegistry::shape_key(&f.0);
        if let Some((t1, t2)) = GeometryRegistry::global().pcurve_range(&e.0, face_key) {
            if t1.is_finite() && t2.is_finite() {
                return (t1, t2);
            }
        }
        Self::edge_parameters(e)
    }

    /// `BRep_Tool::UVPoints(E, F, PFirst, PLast)` (`BRep_Tool.cxx:1090-1106`):    /// the UV coordinates of the edge's first and last vertices on the face's
    /// surface.
    ///
    /// `cxx:1096-1101` reverses a local copy of the edge when the face is
    /// `REVERSED`; the representation walk (`cxx:1017-1036`) then answers
    /// `UVPoints2` (the seam p-curve) when that local edge is `REVERSED` and
    /// the representation is a closed-surface one. Both p-curves are evaluated
    /// at the same `First`/`Last` (`BRep_CurveOnClosedSurface.cxx:123-134`).
    ///
    /// The `Geom_RectangularTrimmedSurface` unwrapping of `cxx:1042-1053` is
    /// folded into `gp_pln`, which already answers the basis plane of a
    /// trimmed plane.
    pub fn uv_points(edge: &Edge, face: &Face) -> (GpPnt2d, GpPnt2d) {
        use occt_core::precision::Precision;

        let face_reversed = face.0.orientation() == Orientation::Reversed;
        let edge_reversed = edge.0.orientation().is_reversed();
        // `BRep_Tool.cxx:1098-1101`: `aLocalEdge.Reverse()`.
        let local_reversed = edge_reversed != face_reversed;

        let face_key = GeometryRegistry::shape_key(&face.0);
        let pcs = GeometryRegistry::global().edge_pcurves(&edge.0, face_key);
        if !pcs.is_empty() {
            let pc = if pcs.len() >= 2 && local_reversed {
                &pcs[1]
            } else {
                &pcs[0]
            };
            let (first, last) = GeometryRegistry::global()
                .pcurve_range(&edge.0, face_key)
                .unwrap_or_else(|| Self::edge_parameters(edge));
            let mut pfirst = GpPnt2d::new(0.0, 0.0);
            let mut plast = GpPnt2d::new(0.0, 0.0);
            if !Precision::is_negative_infinite(first) {
                pfirst = pc.d0(first);
            }
            if !Precision::is_positive_infinite(last) {
                plast = pc.d0(last);
            }
            return (pfirst, plast);
        }

        // `cxx:1039-1083`: no representation -- for a plane, project the two
        // vertices through `ElSLib::Parameters`.
        if let Some(surf) = GeometryRegistry::global().face_surface(&face.0) {
            if let Some(pln) = surf.gp_pln() {
                let pos = pln.position();
                let loc = face.0.location();
                let invert = if loc.is_identity() {
                    None
                } else {
                    loc.transformation().inverted().ok()
                };
                // `cxx:1059-1061`: `Vf.Move(L.Inverted())`.
                let project = |v: &Vertex| -> GpPnt2d {
                    let p = match &invert {
                        Some(t) => Self::vertex_point_world(v).transformed(t),
                        None => Self::vertex_point(v),
                    };
                    let (u, vv) = occt_core::elib::slib::plane_parameters(&pos, &p);
                    GpPnt2d::new(u, vv)
                };
                let (vf, vl) = crate::topo_tools::edge_vertices(edge);
                let pfirst = match vf {
                    Some(v) => project(&v),
                    None => GpPnt2d::new(0.0, 0.0),
                };
                let plast = match vl {
                    Some(v) => project(&v),
                    None => GpPnt2d::new(0.0, 0.0),
                };
                return (pfirst, plast);
            }
        }
        (GpPnt2d::new(0.0, 0.0), GpPnt2d::new(0.0, 0.0))
    }
}

/// The four lines every `BRep_Tool::Tolerance` overload runs on the value it
/// reads off the `TShape` (`BRep_Tool.cxx:140-149`, `:884-893`, `:1323-1332`):
///
/// ```text
/// constexpr double pMin = Precision::Confusion();
/// if (p > pMin) return p; else return pMin;
/// ```
///
/// Note this is `BRep_Tool` behaviour only: the underlying accessors
/// (`BRep_TVertex::Tolerance`, `BRep_TEdge::Tolerance`, `BRep_TFace::Tolerance`)
/// keep answering the stored `RealEpsilon()`, which is what the persistence
/// layer (`BinTools_ShapeSet`, `XmlMDF`) writes out.
fn clamp_to_confusion(p: f64) -> f64 {
    let p_min = occt_core::precision::CONFUSION;
    if p > p_min {
        p
    } else {
        p_min
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::builder::TopoBuilder;
    use occt_core::gp::{GpAx3, GpPln};

    fn edge() -> Edge {
        Edge(TopoShape::new(ShapeType::Edge))
    }

    /// Release registry entries for a shape tree so tests don't leave stale
    /// geometry keyed by a freed Arc address in the process-wide side-table.
    fn clear_tree(s: &TopoShape) {
        GeometryRegistry::global().clear_shape(s);
        let children = s.tshape.read().unwrap().children.clone();
        for c in children {
            clear_tree(&c);
        }
    }

    #[test]
    fn edge_parameters_are_unbounded() {
        let (first, last) = BRepTool::edge_parameters(&edge());
        assert_eq!(first, f64::NEG_INFINITY);
        assert_eq!(last, f64::INFINITY);
    }

    #[test]
    fn unbounded_edge_is_not_closed() {
        assert!(!BRepTool::is_closed_edge(&edge()));
    }

    #[test]
    fn shape_type_str_names_vertex() {
        assert_eq!(BRepTool::shape_type_str(ShapeType::Vertex), "Vertex");
    }

    #[test]
    fn empty_shape_list_is_integral() {
        let shapes: Vec<TopoShape> = Vec::new();
        assert!(BRepTool::check_shape_integrity(&shapes));
    }

    #[test]
    fn make_vertex_stores_point() {
        let b = TopoBuilder::new();
        let v = b.make_vertex(GpPnt::new(1.0, 2.0, 3.0), 1e-3);
        assert!(BRepTool::vertex_point(&v).is_equal(&GpPnt::new(1.0, 2.0, 3.0)));
        assert_eq!(BRepTool::vertex_tolerance(&v), 1e-3);
        clear_tree(&v.0);
    }

    #[test]
    fn make_edge_stores_curve() {
        let b = TopoBuilder::new();
        let e = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(10.0, 0.0, 0.0));
        assert_eq!(BRepTool::edge_parameters(&e), (0.0, 10.0));
        let c = BRepTool::edge_curve(&e).expect("curve present");
        let mid = c.d0(5.0);
        assert!(mid.is_equal(&GpPnt::new(5.0, 0.0, 0.0)));
        assert!(BRepTool::same_parameter(&e));
        assert!(!BRepTool::is_degenerated(&e));
        assert_eq!(BRepTool::parameter_on_edge(&e, 0), 0.0);
        assert_eq!(BRepTool::parameter_on_edge(&e, 1), 10.0);
        clear_tree(&e.0);
    }

    #[test]
    fn segment_has_two_vertices_and_endpoints() {
        let b = TopoBuilder::new();
        let p1 = GpPnt::new(0.0, 0.0, 0.0);
        let p2 = GpPnt::new(3.0, 4.0, 0.0);
        let e = b.make_edge_segment(&p1, &p2);
        assert_eq!(e.0.tshape.read().unwrap().children.len(), 2);
        let (a, z) = BRepTool::edge_vertices(&e).expect("endpoints");
        assert!(a.is_equal(&p1));
        assert!(z.is_equal(&p2));
        clear_tree(&e.0);
    }

    #[test]
    fn make_face_plane_registers_surface() {
        let b = TopoBuilder::new();
        let f = b.make_face_plane(&GpPln::new(GpAx3::standard()));
        assert!(BRepTool::face_surface(&f).is_some());
        let (u1, u2, v1, v2) = BRepTool::uv_bounds(&f);
        assert_eq!(u1, f64::NEG_INFINITY);
        assert_eq!(u2, f64::INFINITY);
        assert_eq!(v1, f64::NEG_INFINITY);
        assert_eq!(v2, f64::INFINITY);
        assert!(BRepTool::natural_restriction(&f));
        clear_tree(&f.0);
    }

    #[test]
    fn child_trees_build() {
        let b = TopoBuilder::new();
        let e1 = b.make_edge_segment(&GpPnt::zero(), &GpPnt::new(1.0, 0.0, 0.0));
        let e2 = b.make_edge_segment(&GpPnt::new(1.0, 0.0, 0.0), &GpPnt::new(1.0, 1.0, 0.0));
        let wire = b.make_wire(&[e1, e2]);
        let face = b.make_face_plane(&GpPln::new(GpAx3::standard()));
        let shell = b.make_shell(&[face.clone()]);
        let solid = b.make_solid(&[shell.clone()]);
        assert_eq!(wire.0.tshape.read().unwrap().children.len(), 2);
        assert_eq!(face.0.tshape.read().unwrap().children.len(), 0);
        assert_eq!(shell.0.tshape.read().unwrap().children.len(), 1);
        assert_eq!(solid.0.tshape.read().unwrap().children.len(), 1);
        clear_tree(&solid.0);
        clear_tree(&wire.0);
    }

    #[test]
    fn unregistered_shapes_fall_back() {
        let v = Vertex::new();
        assert!(BRepTool::vertex_point(&v).is_equal(&GpPnt::zero()));
        assert_eq!(BRepTool::vertex_tolerance(&v), occt_core::precision::CONFUSION);
        assert!(BRepTool::edge_curve(&edge()).is_none());
        assert!(BRepTool::edge_vertices(&edge()).is_none());
        assert!(BRepTool::face_surface(&Face::new()).is_none());
        assert_eq!(
            BRepTool::uv_bounds(&Face::new()),
            (f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY)
        );
        assert!(BRepTool::natural_restriction(&Face::new()));
    }
}
