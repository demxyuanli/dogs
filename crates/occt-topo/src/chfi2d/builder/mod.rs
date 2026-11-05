//! `ChFi2d_Builder` (`ChFi2d_Builder.hxx:38-282`, `ChFi2d_Builder.cxx:1-1289`,
//! `ChFi2d_Builder_0.cxx:1-890`): fillets and chamfers on the wire of a planar
//! face.
//!
//! The class walks a face, replaces the two edges adjacent to a vertex by
//! trimmed copies plus the new fillet/chamfer edge, and records in `history`
//! which basis edge each modified edge was derived from.
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | ctors, `Init`, `IsIssuedFrom` | `Builder.cxx:66-239` | [`init`] |
//! | `AddFillet`, `ModifyFillet`, `RemoveFillet`, `ComputeFillet`, `IsAFillet`, `IsAChamfer` | `Builder.cxx:243-498`, `:504-530`, `:1233-1265` | [`fillet`] |
//! | `BuildFilletEdge` | `Builder.cxx:764-1230` | [`fillet_edge`] |
//! | `AddChamfer`, `ModifyChamfer`, `RemoveChamfer`, `ComputeChamfer` | `_0.cxx:72-308`, `:181-248` | [`chamfer`] |
//! | `BuildChamferEdge`, `ComputePoint`, `OrientChamfer` | `_0.cxx:524-865` | [`chamfer_edge`] |
//! | `BuildNewWire`, `BuildNewEdge` | `Builder.cxx:533-681` | [`wire_build`] |
//! | `UpDateHistory`, `BasisEdge` | `Builder.cxx:684-762` | [`history`] |
//! | `RemoveFillet` / `RemoveChamfer` shared body | `Builder.cxx:290-497`, `_0.cxx:310-519` | [`remove`] |
//!
//! Translation notes:
//! - OCCT's out-parameters (`TopoDS_Edge& TrimE1`, `TopoDS_Vertex& NewExtr1`)
//!   become `Option<...>` returns; `None` is the null handle OCCT leaves behind
//!   on an early `return`,
//! - `NCollection_DataMap<TopoDS_Shape, TopoDS_Shape>` keyed by
//!   `TopTools_ShapeMapHasher` is [`ChFi2dBuilder::history`], an insertion
//!   ordered `Vec` of `(basis, modified)` pairs compared with
//!   [`crate::topo_tools::is_same`],
//! - `BRepLib::BuildCurves3d(filletEdge)` (`Builder.cxx:1228`) is a no-op here:
//!   the port builds the fillet edge with its 3D curve already registered,
//!   which is the `cxx:321-325` early return of `BRepLib::BuildCurve3d`,
//! - `BRepLib_MakeEdge(C2d, S, V1, V2, p1, p2)` + `BRepLib::BuildCurves3d`
//!   (`Builder.cxx:1197`, `:1205`, `:1210`, `:1215`) store a p-curve on the
//!   plane only, and `BuildCurve3d` then lifts it with `GeomLib::To3d`
//!   (`BRepLib.cxx:361-365`). The port performs the same lift in one step
//!   ([`fillet_edge::make_edge_from_2d_circle`]) and stores the 3D curve, which
//!   is exactly the representation `BuildCurve3d` writes.

mod chamfer;
mod chamfer_edge;
mod fillet;
mod fillet_edge;
mod history;
mod init;
mod remove;
mod wire_build;

use std::sync::Arc;

use occt_core::gp::GpPln;
use occt_geom::Surface;
use occt_geom2d::adaptor_curve::type2d;
use occt_core::kernel::geomabs::CurveType;

use crate::abs::Orientation;
use crate::boptools_2d::curve_on_surface;
use crate::chfi2d::ConstructionError;
use crate::shape::{Edge, Face, TopoShape};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{is_same, wires_of_face};

/// `ChFi2d_Builder` (`ChFi2d_Builder.hxx:42-282`).
#[derive(Debug, Clone)]
pub struct ChFi2dBuilder {
    /// `status`.
    status: ConstructionError,
    /// `refFace`: the face handed to `Init`, whose plane every rebuilt face
    /// reuses.
    ref_face: Option<Face>,
    /// `newFace`: the face being edited.
    new_face: Option<Face>,
    /// `fillets`.
    fillets: Vec<TopoShape>,
    /// `chamfers`.
    chamfers: Vec<TopoShape>,
    /// `history` (`NCollection_DataMap<TopoDS_Shape, TopoDS_Shape>`): basis
    /// edge -> its modified descendant.
    history: Vec<(TopoShape, TopoShape)>,
}

impl ChFi2dBuilder {
    /// `Status()` (`ChFi2d_Builder.lxx:87-90`).
    pub fn status(&self) -> ConstructionError {
        self.status
    }

    /// `Result()` (`ChFi2d_Builder.lxx:25-31`): the modified face, carrying the
    /// reference face's orientation.
    pub fn result(&self) -> Option<Face> {
        let face = self.new_face.clone()?;
        let orient = self.ref_face.as_ref().map(|f| f.0.orientation()).unwrap_or(Orientation::Forward);
        Some(Face(face.0.oriented(orient)))
    }

    /// `IsModified(E)` (`ChFi2d_Builder.lxx:37-40`) / `HasDescendant(E)`
    /// (`.lxx:70-73`): both answer `history.IsBound(E)`.
    pub fn is_modified(&self, e: &Edge) -> bool {
        self.history.iter().any(|(k, _)| is_same(k, &e.0))
    }

    /// `HasDescendant(E)` (`ChFi2d_Builder.lxx:70-73`): the same
    /// `history.IsBound(E)` as [`ChFi2dBuilder::is_modified`]. OCCT keeps both
    /// names on the class.
    pub fn has_descendant(&self, e: &Edge) -> bool {
        self.is_modified(e)
    }

    /// `FilletEdges()` (`lxx:44-48`).
    pub fn fillet_edges(&self) -> &[TopoShape] {
        &self.fillets
    }

    /// `NbFillet()` (`lxx:57-60`).
    pub fn nb_fillet(&self) -> usize {
        self.fillets.len()
    }

    /// `ChamferEdges()` (`lxx:52-55`).
    pub fn chamfer_edges(&self) -> &[TopoShape] {
        &self.chamfers
    }

    /// `NbChamfer()` (`lxx:65-68`).
    pub fn nb_chamfer(&self) -> usize {
        self.chamfers.len()
    }

    /// `DescendantEdge(E)` (`lxx:78-81`): `TopoDS::Edge(history.Find(E))`.
    pub fn descendant_edge(&self, e: &Edge) -> Option<Edge> {
        self.find_history(e).and_then(|s| Edge::wrap(s.clone()))
    }

    /// The plane of `refFace`, i.e. `BRepAdaptor_Surface(refFace).Plane()`.
    ///
    /// A face over a *trimmed* plane reports the basis plane, which is what
    /// `Adaptor3d_Surface::Plane` does for `Geom_RectangularTrimmedSurface`
    /// (`Adaptor3d_Surface.cxx:71-80` forwards `Plane()` to the basis).
    fn ref_plane(&self) -> Option<GpPln> {
        let face = self.ref_face.as_ref()?;
        let surf = GeometryRegistry::global().face_surface(&face.0)?;
        surf.gp_pln()
            .or_else(|| surf.rectangular_trimmed_basis().and_then(|b| b.gp_pln()))
    }

    /// The surface handle `refFace` carries (the plane every rebuilt face is
    /// built over).
    fn ref_surface(&self) -> Option<Arc<dyn Surface>> {
        let face = self.ref_face.as_ref()?;
        GeometryRegistry::global().face_surface(&face.0)
    }

    /// `BRepLib_MakeFace(Adaptor3dSurface.Plane(), newWire)` followed by
    /// `newFace = mFace` (`Builder.cxx:597-599`, `_0.cxx:515-518`).
    ///
    /// `Init(Pl, false, Confusion)` + `Add(W)` (`BRepLib_MakeFace.cxx:265-266`)
    /// is `make_face` over the reference plane plus the wire.
    fn rebuild_face(&mut self, wire: &crate::shape::Wire) {
        let Some(surf) = self.ref_surface() else { return };
        let face = crate::builder::TopoBuilder::new().make_face(surf, std::slice::from_ref(wire));
        self.new_face = Some(face);
    }

    /// Whether the reference face's wire is closed (`Builder.cxx:542-549`):
    /// `TopExp_Explorer(refFace, TopAbs_WIRE)` then `aWire.Closed()` on the
    /// first wire found; `true` when the face has no wire, which is the value
    /// the C++ local starts from.
    fn ref_wire_closed(&self) -> bool {
        let Some(face) = self.ref_face.as_ref() else { return true };
        match wires_of_face(face).first() {
            Some(w) => w.0.is_closed(),
            None => true,
        }
    }
}

/// `IsLineOrCircle(E, F)` (`Builder.cxx:1266-1289`, `_0.cxx:869-890`): the
/// edge's p-curve on `F`, unwrapped from `Geom2d_TrimmedCurve`, is a
/// `Geom2d_Circle` or a `Geom2d_Line`.
///
/// `false` where OCCT dereferences a null `CurveOnSurface` handle and throws.
fn is_line_or_circle(edge: &Edge, face: &Face) -> bool {
    let Some(c) = curve_on_surface(edge, face) else {
        return false;
    };
    matches!(type2d(c.as_ref()), CurveType::Circle | CurveType::Line)
}
