use std::sync::Arc;

use occt_core::gp::GpPnt;
use occt_core::poly::Triangulation;
use occt_geom::{Curve, Surface};

use crate::abs::ShapeType;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::explorer_shapes;

use super::helpers::clamp_to_confusion;
use super::BRepTool;

impl BRepTool {
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

    /// `BRep_Tool::PolygonOnTriangulation(E, T, L)` (`BRep_Tool.cxx:683-712`):
    /// the polygon of `theEdge` in the parametric space of the face whose
    /// triangulation is `tri`, i.e. the `BRep_PolygonOnTriangulation`
    /// representation `BRep_Builder::UpdateEdge(E, P, T, L)` installs.
    ///
    /// The OCCT loop matches `cr->IsPolygonOnTriangulation(T, l)` and then, for
    /// a closed-triangulation representation with a reversed edge, answers
    /// `PolygonOnTriangulation2()` instead of `PolygonOnTriangulation()`
    /// (`cxx:698-707`). `reversed` is that `Eisreversed`
    /// (`cxx:689`, `E.Orientation() == TopAbs_REVERSED`), passed in because
    /// `HLRBRep_PolyAlgo::InitBiPointsWithConnexity` reverses a *local copy* of
    /// the edge (`HLRBRep_PolyAlgo.cxx:1076`) before the second lookup
    /// (`:1079`).
    ///
    /// UNPORTED: `TopLoc_Location` — the representation carries no location in
    /// this crate, so `L.Predivided(E.Location())` has nothing to compare and
    /// the triangulation handle is the whole key.
    pub fn polygon_on_triangulation(
        edge: &Edge,
        tri: &Arc<occt_core::poly::Triangulation>,
        reversed: bool,
    ) -> Option<Arc<occt_core::poly::PolyPolygonOnTriangulation>> {
        let ts = edge.0.tshape.read().unwrap();
        let core = ts.edge_core()?;
        let repr = core
            .polygons_on_triangulation
            .iter()
            .find(|r| Arc::ptr_eq(&r.triangulation, tri))?;
        if reversed {
            if let Some(p2) = &repr.polygon2 {
                return Some(p2.clone());
            }
        }
        Some(repr.polygon.clone())
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

    /// `BRep_Tool::MaxTolerance(theShape, theSubShape)`
    /// (`BRep_Tool.cxx:1792-1824`): the maximum `Tolerance` over every sub-shape
    /// of `sub_shape` found by a `TopExp_Explorer`. Any other type answers the
    /// initial `0.0`, exactly as the reference's unhandled branches do.
    pub fn max_tolerance(shape: &TopoShape, sub_shape: ShapeType) -> f64 {
        let mut tol: f64 = 0.0;
        match sub_shape {
            // `cxx:1798-1805`.
            ShapeType::Face => {
                for s in explorer_shapes(shape, ShapeType::Face) {
                    if let Some(f) = Face::wrap(s) {
                        tol = tol.max(Self::face_tolerance(&f));
                    }
                }
            }
            // `cxx:1806-1813`.
            ShapeType::Edge => {
                for s in explorer_shapes(shape, ShapeType::Edge) {
                    if let Some(e) = Edge::wrap(s) {
                        tol = tol.max(Self::edge_tolerance(&e));
                    }
                }
            }
            // `cxx:1814-1821`.
            ShapeType::Vertex => {
                for s in explorer_shapes(shape, ShapeType::Vertex) {
                    if let Some(v) = Vertex::wrap(s) {
                        tol = tol.max(Self::vertex_tolerance(&v));
                    }
                }
            }
            _ => {}
        }
        tol
    }
}
