//! `BRepFilletAPI_MakeFillet2d` (`BRepFilletAPI_MakeFillet2d.cxx:24-171`,
//! `.lxx:18-71`, `.hxx:42-297`).

use crate::chfi2d::{ChFi2dBuilder, ConstructionError};
use crate::shape::{Edge, Face, TopoShape, Vertex};

/// `BRepFilletAPI_MakeFillet2d` (`hxx:54-297`).
///
/// `myMakeChFi2d` is the [`ChFi2dBuilder`] the class wraps; `shape` / `done` /
/// `generated` are `BRepBuilderAPI_MakeShape`'s `myShape` / `myDone` /
/// `myGenerated` (`BRepBuilderAPI_MakeShape.hxx:57-60`).
#[derive(Debug, Clone)]
pub struct MakeFillet2d {
    /// `myMakeChFi2d`.
    make_chfi2d: ChFi2dBuilder,
    /// `myShape`.
    shape: Option<TopoShape>,
    /// `myDone`.
    done: bool,
    /// `myGenerated`.
    generated: Vec<TopoShape>,
}

impl Default for MakeFillet2d {
    fn default() -> Self {
        Self::new()
    }
}

impl MakeFillet2d {
    /// `BRepFilletAPI_MakeFillet2d()` (`cxx:26`): default-constructs
    /// `myMakeChFi2d`, whose `status` starts at `ChFi2d_NotPlanar`
    /// (`ChFi2d_Builder.cxx:66-69`).
    pub fn new() -> Self {
        Self {
            make_chfi2d: ChFi2dBuilder::new(),
            shape: None,
            done: false,
            generated: Vec::new(),
        }
    }

    /// `BRepFilletAPI_MakeFillet2d(F)` (`cxx:29-32`): `myMakeChFi2d.Init(F)`.
    pub fn with_face(f: &Face) -> Self {
        Self {
            make_chfi2d: ChFi2dBuilder::with_face(f),
            shape: None,
            done: false,
            generated: Vec::new(),
        }
    }

    /// `Init(F)` (`cxx:36-39`).
    pub fn init(&mut self, f: &Face) {
        self.make_chfi2d.init(f);
    }

    /// `Init(RefFace, ModFace)` (`cxx:44-47`).
    pub fn init_from_faces(&mut self, ref_face: &Face, mod_face: &Face) {
        self.make_chfi2d.init_from_faces(ref_face, mod_face);
    }

    /// `AddFillet(V, Radius)` (`cxx:50-53`).
    ///
    /// `None` is the null edge OCCT returns whenever the construction status is
    /// not `ChFi2d_IsDone`.
    pub fn add_fillet(&mut self, v: &Vertex, radius: f64) -> Option<Edge> {
        self.make_chfi2d.add_fillet(&v.0, radius)
    }

    /// `ModifyFillet(Fillet, Radius)` (`cxx:57-60`).
    pub fn modify_fillet(&mut self, fillet: &Edge, radius: f64) -> Option<Edge> {
        self.make_chfi2d.modify_fillet(fillet, radius)
    }

    /// `RemoveFillet(Fillet)` (`cxx:64-67`): the vertex the fillet was attached
    /// to, or `None` for OCCT's null vertex.
    pub fn remove_fillet(&mut self, fillet: &Edge) -> Option<Vertex> {
        self.make_chfi2d
            .remove_fillet(fillet)
            .and_then(Vertex::wrap)
    }

    /// `AddChamfer(E1, E2, D1, D2)` (`cxx:72-78`).
    pub fn add_chamfer(&mut self, e1: &Edge, e2: &Edge, d1: f64, d2: f64) -> Option<Edge> {
        self.make_chfi2d.add_chamfer_d1d2(e1, e2, d1, d2)
    }

    /// `AddChamfer(E, V, D, Ang)` (`cxx:82-88`).
    pub fn add_chamfer_with_vertex(
        &mut self,
        e: &Edge,
        v: &Vertex,
        d: f64,
        ang: f64,
    ) -> Option<Edge> {
        self.make_chfi2d.add_chamfer_da(e, &v.0, d, ang)
    }

    /// `ModifyChamfer(Chamfer, E1, E2, D1, D2)` (`cxx:91-98`).
    ///
    /// `E1` is commented out in OCCT's definition (`ChFi2d_Builder_0.cxx:253`)
    /// and the port's `ChFi2dBuilder::modify_chamfer_d1d2` only takes `E2`.
    pub fn modify_chamfer(
        &mut self,
        chamfer: &Edge,
        e2: &Edge,
        d1: f64,
        d2: f64,
    ) -> Option<Edge> {
        self.make_chfi2d.modify_chamfer_d1d2(chamfer, e2, d1, d2)
    }

    /// `ModifyChamfer(Chamfer, E, D, Ang)` (`cxx:102-108`).
    pub fn modify_chamfer_with_vertex(
        &mut self,
        chamfer: &Edge,
        e: &Edge,
        d: f64,
        ang: f64,
    ) -> Option<Edge> {
        self.make_chfi2d.modify_chamfer_da(chamfer, e, d, ang)
    }

    /// `RemoveChamfer(Chamfer)` (`cxx:112-115`).
    pub fn remove_chamfer(&mut self, chamfer: &Edge) -> Option<Vertex> {
        self.make_chfi2d
            .remove_chamfer(chamfer)
            .and_then(Vertex::wrap)
    }

    /// `BasisEdge(E)` (`cxx:119-122` -> `ChFi2d_Builder.cxx:747-762`): the input
    /// edge when it is not a descendant.
    pub fn basis_edge(&self, e: &Edge) -> Edge {
        self.make_chfi2d.basis_edge(e)
    }

    /// `Build(theRange)` (`cxx:126-137`): `Done()` + `myShape = Result()` when
    /// the builder status is `ChFi2d_IsDone`, `NotDone()` otherwise.
    pub fn build(&mut self) {
        if self.make_chfi2d.status() == ConstructionError::IsDone {
            self.done = true;
            self.shape = self.make_chfi2d.result().map(|f| f.0);
        } else {
            self.done = false;
        }
    }

    /// `BRepBuilderAPI_MakeShape::IsDone()` (`BRepBuilderAPI_MakeShape.hxx:47`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `BRepBuilderAPI_MakeShape::Shape()` (`BRepBuilderAPI_MakeShape.hxx:50`).
    pub fn shape(&self) -> Option<&TopoShape> {
        self.shape.as_ref()
    }

    /// `Modified(S)` (`cxx:144-148`): `myGenerated` is cleared and refilled with
    /// `DescendantEdge(TopoDS::Edge(S))`.
    ///
    /// OCCT raises `Standard_TypeMismatch` when `S` is not an edge and
    /// `Standard_NoSuchObject` when it has no descendant; the port answers an
    /// empty list on both arms.
    pub fn modified(&mut self, s: &TopoShape) -> &[TopoShape] {
        self.generated.clear();
        if let Some(edge) = Edge::wrap(s.clone()) {
            if let Some(descendant) = self.make_chfi2d.descendant_edge(&edge) {
                self.generated.push(descendant.0);
            }
        }
        &self.generated
    }

    /// `NbCurves()` (`cxx:152-155`).
    pub fn nb_curves(&self) -> usize {
        self.nb_fillet() + self.nb_chamfer()
    }

    /// `NewEdges(I)` (`cxx:158-171`): 1-based over the fillets first, then the
    /// chamfers.
    ///
    /// OCCT's `NCollection_Sequence::Value` raises `Standard_OutOfRange` outside
    /// `1..NbCurves()`; the port answers an empty list.
    pub fn new_edges(&mut self, i: usize) -> &[TopoShape] {
        self.generated.clear();
        if i == 0 {
            return &self.generated;
        }
        let nf = self.nb_fillet();
        if i <= nf {
            self.generated.push(self.fillet_edges()[i - 1].clone());
        } else if i - nf <= self.nb_chamfer() {
            self.generated.push(self.chamfer_edges()[i - nf - 1].clone());
        }
        &self.generated
    }

    /// `IsModified(E)` (`.lxx:19-22`).
    pub fn is_modified(&self, e: &Edge) -> bool {
        self.make_chfi2d.is_modified(e)
    }

    /// `FilletEdges()` (`.lxx:26-29`).
    pub fn fillet_edges(&self) -> &[TopoShape] {
        self.make_chfi2d.fillet_edges()
    }

    /// `NbFillet()` (`.lxx:33-36`).
    pub fn nb_fillet(&self) -> usize {
        self.make_chfi2d.nb_fillet()
    }

    /// `ChamferEdges()` (`.lxx:40-43`).
    pub fn chamfer_edges(&self) -> &[TopoShape] {
        self.make_chfi2d.chamfer_edges()
    }

    /// `NbChamfer()` (`.lxx:47-50`).
    pub fn nb_chamfer(&self) -> usize {
        self.make_chfi2d.nb_chamfer()
    }

    /// `HasDescendant(E)` (`.lxx:54-57`).
    pub fn has_descendant(&self, e: &Edge) -> bool {
        self.make_chfi2d.has_descendant(e)
    }

    /// `DescendantEdge(E)` (`.lxx:61-64`): `None` is OCCT's
    /// `Standard_NoSuchObject` on an unbound edge.
    pub fn descendant_edge(&self, e: &Edge) -> Option<Edge> {
        self.make_chfi2d.descendant_edge(e)
    }

    /// `Status()` (`.lxx:68-71`): `myMakeChFi2d.Status()`.
    pub fn status(&self) -> ConstructionError {
        self.make_chfi2d.status()
    }

    /// The wrapped builder, for the callers that need the `ChFi2d_Builder`
    /// queries `BRepFilletAPI_MakeFillet2d` does not forward.
    pub fn builder(&self) -> &ChFi2dBuilder {
        &self.make_chfi2d
    }
}
