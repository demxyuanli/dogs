//! Topological shape builder — constructs shapes from sub-shapes.
//! Source: `BRep_Builder`
use std::sync::Arc;

use occt_core::gp::{GpAx2, GpCirc, GpDir, GpLin, GpPln, GpPnt, GpTrsf, GpVec};
use occt_core::poly::Triangulation;
use occt_geom::{Curve, GeomCircle, GeomLine, GeomPlane, Surface};

use crate::abs::{Orientation, ShapeType};
use crate::shape::{Compound, Edge, Face, Shell, Solid, TopoShape, Vertex, Wire};
use crate::tgeometry::{EdgeGeom, FaceGeom, GeometryRegistry, VertexGeom};

/// Builds topological shapes (mirrors BRep_Builder).
#[derive(Debug, Default)]
pub struct TopoBuilder;

impl TopoBuilder {
    pub fn new() -> Self { Self }

    /// Make an empty shape of given type.
    pub fn make_shape(&self, shape_type: ShapeType) -> TopoShape {
        TopoShape::new(shape_type)
    }

    /// Add sub-shape. Source: `TopoDS_Builder::Add` (TopoDS_Builder.cxx).
    ///
    /// Stores a full `TopoShape` (TShape + location + orientation). If the
    /// parent *view* is `REVERSED`, the stored child is reversed; if the
    /// parent location is not identity, the stored child is moved by
    /// `parent.location.inverted()`. The component TShape is frozen
    /// (`Free(false)`). Type compatibility is the OCCT `aTb` table; an
    /// incompatible pair is skipped (OCCT throws `TopoDS_UnCompatibleShapes`).
    pub fn add(&self, shape: &mut TopoShape, sub: &TopoShape) {
        // OCCT freezes the component first (also guards self-insertion).
        sub.set_free(false);
        if !shape.free() {
            return;
        }
        if !self.compatible(shape.shape_type(), sub.shape_type()) {
            return;
        }
        let mut child = sub.clone();
        if shape.orientation() == crate::abs::Orientation::Reversed {
            child.reverse();
        }
        if !shape.location().is_identity() {
            child.move_location(&shape.location().inverted());
        }
        if let Ok(mut t) = shape.tshape.write() {
            t.add_child(child);
            t.set_modified(true);
        }
    }

    /// OCCT `TopoDS_Builder::Add` compatibility bits: which parent types may
    /// contain a given component type.
    fn compatible(&self, parent: ShapeType, component: ShapeType) -> bool {
        use ShapeType::*;
        match component {
            Compound => parent == Compound,
            CompSolid => parent == Compound,
            Solid => parent == Compound || parent == CompSolid,
            Shell => parent == Compound || parent == Solid,
            Face => parent == Compound || parent == Shell,
            Wire => parent == Compound || parent == Face,
            Edge => parent == Compound || parent == Solid || parent == Wire,
            Vertex => parent == Compound || parent == Solid || parent == Face || parent == Edge,
            Shape => false,
        }
    }

    /// Add edge to wire.
    pub fn add_edge(&self, wire: &mut Wire, edge: &Edge) { self.add(&mut wire.0, &edge.0); }
    /// Add wire to face.
    pub fn add_wire(&self, face: &mut Face, wire: &Wire) { self.add(&mut face.0, &wire.0); }
    /// Add face to shell.
    pub fn add_face(&self, shell: &mut Shell, face: &Face) { self.add(&mut shell.0, &face.0); }
    /// Add shell to solid.
    pub fn add_shell(&self, solid: &mut Solid, shell: &Shell) { self.add(&mut solid.0, &shell.0); }
    /// Add any shape to compound.
    pub fn add_compound(&self, comp: &mut Compound, shape: &TopoShape) { self.add(&mut comp.0, shape); }

    /// Make a vertex with a point. The point and tolerance are registered in
    /// the geometry side-table (`BRep_TVertex`).
    pub fn make_vertex(&self, p: GpPnt, tolerance: f64) -> Vertex {
        let v = Vertex::new();
        v.set_point(p);
        v.set_tolerance(tolerance);
        v.set_free(true);
        v
    }

    /// Make an edge over `curve` on the parameter range `[first, last]`.
    /// Registers an `EdgeGeom` in the geometry side-table (`BRep_TEdge`).
    pub fn make_edge(&self, curve: Arc<dyn Curve>, first: f64, last: f64) -> Edge {
        let e = Edge::new();
        GeometryRegistry::global().set_edge(&e.0, EdgeGeom::new(curve, first, last));
        e
    }

    /// `BRep_Builder::MakeEdge(E, thePolygon)` (`BRep_Builder.lxx:43-47`):
    /// `MakeEdge(E)` then `UpdateEdge(E, P)`, which appends a `BRep_Polygon3D`
    /// representation holding `polygon` and the identity location
    /// (`BRep_Builder.cxx:749-785`). The edge has no `Geom_Curve`.
    pub fn make_edge_polygon(&self, polygon: Arc<occt_core::poly::Polygon3D>) -> Edge {
        let e = Edge::new();
        {
            let mut ts = e.0.tshape.write().unwrap();
            ts.edge_core_mut().polygon3d = Some(polygon);
            // `cxx:783`: `TE->Modified(true)`.
            ts.set_modified(true);
        }
        e
    }

    /// `BRepPrim_Builder::MakeDegeneratedEdge` (`BRepPrim_OneAxis.cxx:934`,
    /// `:1212`, `:1264`): an edge whose 3D representation is the single point
    /// `p` (a sphere/cone pole, a torus seam point). OCCT stores a `BRep_Curve3D`
    /// representation holding a **null** curve — `BRep_Builder::Degenerated(E,
    /// true)` calls `UpdateCurves(TE->ChangeCurves(), null, E.Location())`
    /// (`BRep_Builder.cxx:1073-1087`), which appends the representation with a
    /// null curve (`cxx:80-92`) — so `BRep_Tool::Curve` answers null while the
    /// representation itself is present. That presence is load-bearing for
    /// `BRepCheck_Edge::Minimum` (`BRepCheck_Edge.cxx:100-147`): a present
    /// representation with a null curve raises neither `No3DCurve`
    /// (`cxx:118-123`, `exist` is set by the representation) nor
    /// `InvalidDegeneratedFlag` (`cxx:144-147`, `myCref` stays null), and the
    /// range check (`cxx:149-252`) is skipped because `myCref` is null. The mesh
    /// pipeline then uses the edge's pcurve only (`meshing/edge_discret.rs:720-
    /// 723`, `:744`). `EdgeGeom::curve` is non-nullable here, so a constant-point
    /// curve carries `p`; every consumer must key off the `degenerated` flag and
    /// `curve3d_null`, never off this curve's shape.
    pub fn make_degenerated_edge(&self, p: &GpPnt) -> Edge {
        let e = self.make_edge(Arc::new(DegeneratePointCurve { p: *p }), 0.0, 0.0);
        GeometryRegistry::global().set_degenerated(&e.0, true);
        GeometryRegistry::global().set_curve3d_null(&e.0, true);
        e
    }

    /// Make a straight segment from `p1` to `p2`, adding registered endpoint
    /// vertices as children (`BRepBuilderAPI_MakeEdge(P1, P2)`).
    ///
    /// The endpoint vertices carry `BRepLib::Precision()` — `Init(C, S, P1, P2)`
    /// builds them with `Tol = BRepLib::Precision()`, i.e. `Precision::Confusion()`
    /// (`BRepLib_MakeEdge.cxx:826-830`). That tolerance is load-bearing: every
    /// consumer that derives a working precision from vertex tolerances
    /// (`BRepTools_WireExplorer`'s `dfVertToler`, `BRepLib::UpdateTolerances`)
    /// degenerates when a constructed vertex carries no tolerance at all.
    pub fn make_edge_segment(&self, p1: &GpPnt, p2: &GpPnt) -> Edge {
        let dir = GpDir::from_vec(&GpVec::from_pnts(p1, p2))
            .expect("make_edge_segment: p1 and p2 must be distinct");
        let lin = GpLin::from_pnt_dir(*p1, dir);
        let mut e = self.make_edge(Arc::new(GeomLine::new(lin)), 0.0, p1.distance(p2));
        let v1 = self.make_vertex(*p1, occt_core::precision::CONFUSION);
        let v2 = self.make_vertex(*p2, occt_core::precision::CONFUSION);
        self.add_edge_vertices(&mut e, &v1, &v2);
        e
    }

    /// Make a straight segment edge from `p1` to `p2`, reusing the given
    /// endpoint vertices instead of creating fresh ones — so edges sharing an
    /// endpoint reference the same vertex `TShape` (needed for a meaningful
    /// Euler characteristic / vertex-sharing across faces).
    pub fn make_edge_segment_with_vertices(
        &self,
        p1: &GpPnt,
        p2: &GpPnt,
        v1: &Vertex,
        v2: &Vertex,
    ) -> Edge {
        let dir = GpDir::from_vec(&GpVec::from_pnts(p1, p2))
            .expect("make_edge_segment_with_vertices: p1 and p2 must be distinct");
        let lin = GpLin::from_pnt_dir(*p1, dir);
        let mut e = self.make_edge(Arc::new(GeomLine::new(lin)), 0.0, p1.distance(p2));
        self.add_edge_vertices(&mut e, v1, v2);
        e
    }

    /// `BRepLib_MakeEdge`: first vertex FORWARD, last vertex REVERSED.
    pub(crate) fn add_edge_vertices(&self, e: &mut Edge, v1: &Vertex, v2: &Vertex) {
        let mut a = v1.0.clone();
        a.set_orientation(Orientation::Forward);
        let mut b = v2.0.clone();
        b.set_orientation(Orientation::Reversed);
        self.add(&mut e.0, &a);
        self.add(&mut e.0, &b);
    }

    /// Make a circular arc edge in the plane `axis` with the given radius and
    /// parameter range (`BRepBuilderAPI_MakeEdge(gp_Circ, ...)`).
    pub fn make_edge_circle(&self, axis: &GpAx2, radius: f64, first: f64, last: f64) -> Edge {
        self.make_edge(Arc::new(GeomCircle::new(GpCirc::new(*axis, radius))), first, last)
    }

    /// Make a wire containing `edges` (`BRep_Builder` level append:
    /// `BRep_Builder::MakeWire` + `BRep_Builder::Add`).
    ///
    /// **Not** `BRepBuilderAPI_MakeWire`: this only stores the edges in the given
    /// order with the orientation each occurrence already carries. OCCT's
    /// `BRepLib_MakeWire::Add` additionally reverses an edge when the accumulated
    /// wire requires it and unifies its endpoints with the wire's existing
    /// vertices by geometric proximity (`BRepLib_MakeWire.cxx:223-285`,
    /// `:304-368`); [`crate::brep_lib_make_wire::MakeWire`] is that class, and
    /// call sites that model `BRepBuilderAPI_MakeWire` must go through it
    /// (`builder_area::make_wire_from_edges`).
    pub fn make_wire(&self, edges: &[Edge]) -> Wire {
        let mut wire = Wire::new();
        for e in edges {
            self.add_edge(&mut wire, e);
        }
        if std::env::var("WDBGNEW").is_ok() {
            let oris = edges
                .iter()
                .map(|e| format!("{:?}", e.0.orientation()))
                .collect::<Vec<_>>()
                .join(",");
            eprintln!(
                "MAKEWIRE wptr={:#x} n={} [{}]",
                Arc::as_ptr(&wire.0.tshape) as usize,
                edges.len(),
                oris
            );
        }
        wire
    }

    /// Make a face over `surface`, adding `wires` as children. Registers a
    /// `FaceGeom` (`BRep_TFace`) with natural restriction (no trimming wires).
    pub fn make_face(&self, surface: Arc<dyn Surface>, wires: &[Wire]) -> Face {
        let mut face = Face::new();
        GeometryRegistry::global().set_face(&face.0, FaceGeom::new(surface));
        for w in wires {
            self.add_wire(&mut face, w);
        }
        face
    }

    /// `BRep_Builder::UpdateFace(F, S, L, Tol)` (`BRep_Builder.cxx:508-536`):
    /// `make_face` followed by the explicit tolerance the construction site
    /// passes, instead of the `BRep_TFace` default `RealEpsilon()`.
    ///
    /// The tolerance is load bearing for the mesh pipeline: `BRepBndLib::Add`
    /// hands it to `BndLib_AddSurface::Add` (`BRepBndLib.cxx:109`), which lands
    /// in `Bnd_Box::Gap`.
    pub fn make_face_with_tolerance(
        &self,
        surface: Arc<dyn Surface>,
        wires: &[Wire],
        tolerance: f64,
    ) -> Face {
        let face = self.make_face(surface, wires);
        GeometryRegistry::global().set_face_tolerance(&face.0, tolerance);
        face
    }

    /// Make an unbounded planar face (`BRepBuilderAPI_MakeFace(gp_Pln)`).
    ///
    /// The whole plane is a *natural restriction* face (`BRepPrim_FaceBuilder.cxx:168`
    /// sets `NaturalRestriction(myFace, true)` for the unbounded face builder),
    /// which is what makes `BRepToIGES_BRShell.cxx:381` write it as a whole
    /// surface. A *bounded* face — `make_face` plus wires — keeps the
    /// `BRep_TFace` default (false).
    pub fn make_face_plane(&self, pln: &GpPln) -> Face {
        let face = self.make_face(Arc::new(GeomPlane::new(pln.clone())), &[]);
        GeometryRegistry::global().set_natural_restriction(&face.0, true);
        face
    }

    /// `BRep_Builder::MakeFace(F, theTriangulation)` (`BRep_Builder.cxx:516-527`):
    /// a face whose only content is a triangulation — no surface, no wires.
    ///
    /// This is the glTF "proxy" face (`RWGltf_GltfJsonParser.cxx:2002-2006`);
    /// `BRep_TFace::Triangulation(tri)` defaults `theToReset` to TRUE
    /// (`BRep_TFace.hxx:99-100`), so the list holds exactly this triangulation
    /// and it is active.
    pub fn make_face_from_triangulation(&self, tri: Arc<Triangulation>) -> Face {
        let mut face = Face::new();
        {
            let mut ts = face.0.tshape.write().unwrap();
            ts.set_face_triangulation(Some(tri), true);
        }
        face
    }

    /// `BRep_Builder::UpdateFace(theFace, theTriangulation, theToReset)`
    /// (`BRep_Builder.cxx:582-594`): replace (or re-activate) the face's
    /// triangulation and flag the TShape as modified.
    pub fn update_face_triangulation(
        &self,
        face: &Face,
        tri: Option<Arc<Triangulation>>,
        to_reset: bool,
    ) {
        {
            let mut ts = face.0.tshape.write().unwrap();
            ts.set_face_triangulation(tri, to_reset);
            // `cxx:593`: `theFace.TShape()->Modified(true)`.
            ts.set_modified(true);
        }
    }

    /// Make a shell containing `faces` (`BRepBuilderAPI_MakeShell`).
    pub fn make_shell(&self, faces: &[Face]) -> Shell {
        let mut shell = Shell::new();
        for f in faces {
            self.add_face(&mut shell, f);
        }
        shell
    }

    /// Make a solid containing `shells` (`BRepBuilderAPI_MakeSolid`).
    pub fn make_solid(&self, shells: &[Shell]) -> Solid {
        let mut solid = Solid::new();
        for s in shells {
            self.add_shell(&mut solid, s);
        }
        solid
    }

    /// Make a compound of `shapes` (`BRepBuilderAPI_MakeCompound`).
    pub fn make_compound_of(&self, shapes: &[TopoShape]) -> Compound {
        let mut comp = Compound::new();
        for s in shapes {
            self.add_compound(&mut comp, s);
        }
        comp
    }

    /// Set the free flag.
    pub fn make_free(&self, shape: &TopoShape) { shape.set_free(true); }

    /// `TopoDS_Builder::Remove` (`TopoDS_Builder.cxx:107-133`). A frozen parent
    /// is left unchanged; OCCT throws `TopoDS_FrozenShape` in that case.
    pub fn remove(&self, shape: &mut TopoShape, sub: &TopoShape) {
        if !shape.free() {
            return;
        }
        let mut key = sub.clone();
        if shape.orientation() == crate::abs::Orientation::Reversed {
            key.reverse();
        }
        if !shape.location().is_identity() {
            key.move_location(&shape.location().inverted());
        }
        if let Ok(mut t) = shape.tshape.write() {
            if let Some(i) = t.children.iter().position(|c| {
                std::sync::Arc::ptr_eq(&c.tshape, &key.tshape)
                    && c.orientation() == key.orientation()
                    && c.location().is_equal(key.location())
            }) {
                t.children.remove(i);
                t.set_modified(true);
            }
        }
    }
}

impl Vertex {
    /// Set the vertex tolerance: sets the check flag and, when geometry is
    /// registered, keeps the side-table `VertexGeom` tolerance in sync.
    pub fn set_tolerance(&self, tol: f64) {
        self.0.tshape.write().unwrap().flags.check = tol > 1e-15;
        if let Some(mut g) = GeometryRegistry::global().vertex_geom(&self.0) {
            g.tolerance = tol;
            GeometryRegistry::global().set_vertex(&self.0, g);
        }
    }

    /// Register the vertex point in the geometry side-table, preserving any
    /// tolerance already registered (`BRep_TVertex::Pnt`).
    pub fn set_point(&self, p: GpPnt) {
        let tol = GeometryRegistry::global().vertex_geom(&self.0).map(|g| g.tolerance).unwrap_or(0.0);
        GeometryRegistry::global().set_vertex(&self.0, VertexGeom { point: p, tolerance: tol });
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use occt_core::gp::GpAx3;

    fn nb_children(s: &TopoShape) -> usize {
        s.tshape.read().unwrap().children.len()
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
    fn build_vertex() {
        let b = TopoBuilder::new();
        let v = b.make_vertex(GpPnt::new(1.,2.,3.), 0.001);
        assert!(v.is_vertex());
        assert!(v.free());
        clear_tree(&v.0);
    }

    #[test]
    fn make_compound() {
        let b = TopoBuilder::new();
        let mut c = Compound::new();
        let v = b.make_vertex(GpPnt::zero(), 0.0);
        b.add_compound(&mut c, &v.0);
        assert_eq!(c.shape_type(), ShapeType::Compound);
        clear_tree(&v.0);
        clear_tree(&c.0);
    }

    #[test]
    fn add_reversed_parent_stores_reversed_child() {
        let b = TopoBuilder::new();
        let mut wire = Wire::new();
        wire.0.set_orientation(crate::abs::Orientation::Reversed);
        let e = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(1.0, 0.0, 0.0));
        b.add_edge(&mut wire, &e);
        let stored = wire.0.tshape.read().unwrap().children[0].clone();
        assert!(stored.orientation().is_reversed());
        let viewed = crate::iterator::ShapeIterator::of_shape(&wire.0).next().unwrap();
        assert!(viewed.orientation().is_forward());
        assert!(!e.free());
        clear_tree(&wire.0);
    }

    #[test]
    fn wire_shell_solid_build_child_trees() {
        let b = TopoBuilder::new();
        let e1 = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(1.0, 0.0, 0.0));
        let e2 = b.make_edge_segment(&GpPnt::new(1.0, 0.0, 0.0), &GpPnt::new(1.0, 1.0, 0.0));
        let wire = b.make_wire(&[e1, e2]);
        assert_eq!(nb_children(&wire.0), 2);

        let face = b.make_face_plane(&GpPln::new(GpAx3::standard()));
        assert_eq!(nb_children(&face.0), 0);

        let shell = b.make_shell(&[face.clone()]);
        assert_eq!(nb_children(&shell.0), 1);

        let solid = b.make_solid(&[shell.clone()]);
        assert_eq!(nb_children(&solid.0), 1);

        let comp = b.make_compound_of(&[wire.0.clone(), face.0.clone()]);
        assert_eq!(nb_children(&comp.0), 2);
        clear_tree(&comp.0);
    }
}

/// Port-internal stand-in for OCCT's **absent** 3D curve of a degenerated edge
/// (`BRep_Builder::Degenerated` drops it, `BRep_Builder.cxx:1073-1085`).
/// `EdgeGeom::curve` is a non-nullable `Arc<dyn Curve>` in this port, so a
/// constant-point curve keeps `p` available while every derivative is zero and
/// `[First, Last] = [0, 0]`. Consumers must rely on the `degenerated` flag
/// (`GeometryRegistry::is_degenerated_edge`), never on this curve's shape.
#[derive(Debug, Clone)]
struct DegeneratePointCurve {
    p: GpPnt,
}

impl Curve for DegeneratePointCurve {
    fn d0(&self, _u: f64) -> GpPnt {
        self.p
    }
    fn d1(&self, _u: f64) -> (GpPnt, GpVec) {
        (self.p, GpVec::zero())
    }
    fn d2(&self, _u: f64) -> (GpPnt, GpVec, GpVec) {
        (self.p, GpVec::zero(), GpVec::zero())
    }
    fn first_parameter(&self) -> f64 {
        0.0
    }
    fn last_parameter(&self) -> f64 {
        0.0
    }
    fn continuity(&self) -> u8 {
        6 // GeomAbs_CN
    }
    fn transform(&mut self, t: &GpTrsf) {
        self.p = self.p.transformed(t);
    }
    fn reverse(&mut self) {}
    fn clone_dyn(&self) -> Box<dyn Curve> {
        Box::new(self.clone())
    }
}
