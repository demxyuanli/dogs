//! `BRepOffset_Offset` (`BRepOffset_Offset.hxx`, `.cxx`): the elementary
//! offset of one face, kept per source face by `BRepOffset_MakeOffset`
//! (`myMapSF`) and queried through `Face()` / `Generated(E)` / `Status()`.
//!
//! Ported:
//!
//! * `Init(Face, Offset, Created, OffsetOutside, JoinType)`
//!   (`BRepOffset_Offset.cxx:466-1094`): `:473-478` source record and
//!   orientation-corrected offset, `:480-499` trimmed-surface unwrap and
//!   location bake, `:500-522` cone particular case, `:524`
//!   `BRepOffset::Surface`, `:526-807` degenerated-edge singularity handling
//!   (`singularity.rs`), `:811-820` face creation, `:822-870` `MapSS` seeding
//!   from `Created`, `:872-1089` wire rebuild, `:1091-1093` orientation and
//!   `BRepTools::Update`.
//! * `Face()` (`:1728-1731`), `InitialShape()` (`.lxx:22-25`) and `Status()`
//!   (`:1811-1814`).
//! * `Generated()` (`:1735-1807`), including its positional pairing of the
//!   source face edges with the offset face edges.
//! * `ComputeCurve3d` (`:190-385`) and the three `UpdateEdge` helpers
//!   (`:115-183`), in `compute_curve3d.rs`.
//! * `Init(Path, Edge1, Edge2, Offset, Polynomial, Tol, Conti)` (`:1098-1107`)
//!   and `Init(Path, Edge1, Edge2, Offset, FirstEdge, LastEdge, ...)`
//!   (`:1112-1553`), in `init_pipe.rs`. `BRepOffset_MakeOffset::
//!   BuildOffsetByArc` (`:1988`) is the caller.
//! * `Init(Vertex, LEdge, Offset, Polynomial, TolApp, Conti)`
//!   (`:1559-1696`), in `init_sphere.rs`; `BuildOffsetByArc`'s
//!   sphere-on-vertex step (`BRepOffset_MakeOffset.cxx:2046`) is the caller.
//!
//! Wiring: `curve_face_offset::offset_face` sends every face that is not an
//! exact plane or sphere through [`BRepOffsetOffset::init_face`] with an empty
//! `Created`, `OffsetOutside = true` and `JoinType = GeomAbs_Arc` (the
//! defaults `BRepOffset_MakeOffset::Initialize` install), which replaces the
//! arm that used to keep the source wires.
//!
//! UNPORTED, with the OCCT lines that keep it out:
//!
//! * `Init(Edge, Offset)` (`cxx:1698-1724`): the `GeomFill_Pipe` on a curve
//!   translated by its own location, in `init_edge.rs`;
//!   `BiTgte_Blend::ComputeCenters` (`BiTgte_Blend.cxx:1401`) is its caller for
//!   a free boundary the ball rolls on.
//! * `myMap` (`hxx:157`): `BRepOffset_Offset.cxx` never binds it (no `myMap`
//!   occurrence in the file), so the port keeps no map.
//!
//! Port guards, each replacing an OCCT exception or a null-handle
//! dereference, never a branch:
//!
//! * a source edge without a p-curve on `CurFace` (`cxx:888`, `:894`, `:1045`)
//!   or with an endpoint missing from `TopExp::Vertices` is dropped from the
//!   rebuilt wire: OCCT dereferences the null handle there.
//! * `BRep_Tool::Parameter` raising `Standard_NoSuchObject` (`cxx:912-913`,
//!   `:1046-1047`) maps to the same drop, through
//!   `parameter_of_vertex_on_face_checked`.
//! * `BRepOffset::Surface` raising `Standard_ConstructionError` becomes the
//!   `Err` of [`BRepOffsetOffset::init_face`].
//! * `Geom_Line(P, Dir)` with a null `Dir` (`cxx:650`, `:691`, `:732`, `:775`)
//!   throws; the port skips that singularity block.
//! * `BRepTools::Update` (`cxx:1093`) only drops cached UV points and the
//!   triangulation; the port has neither.
//!
//! Location model: the new face carries the source location `L` (`cxx:815`),
//! the new edges carry the local curve and the same `L` as their view
//! location (OCCT's `BRep_GCurve` representation location), and the new
//! vertices carry the location-baked point with an identity location
//! (`cxx:929-933`, `:977-978`), exactly as OCCT arranges it.

mod compute_curve3d;
mod init_edge;
mod init_face;
mod init_pipe;
mod init_sphere;
mod singularity;

use crate::abs::{Orientation, ShapeType};
use crate::brep_offset::tool::ShapeMap;
use crate::brep_offset::{BRepOffsetStatus, BRepOffsetTool};
use crate::shape::{Edge, Face, TopoShape, Wire};
use crate::topo_tools::{edges_of_wire, explorer_shapes, is_same, wires_of_face};

pub use compute_curve3d::{
    compute_curve3d, update_edge_curve3d, update_edge_pcurve, update_edge_pcurves,
};

/// `BRepOffset_Offset` (`BRepOffset_Offset.hxx:45-158`).
#[derive(Debug, Clone)]
pub struct BRepOffsetOffset {
    /// `myShape` (`hxx:154`): the source face (or, once the unported
    /// constructors exist, edge/vertex).
    my_shape: Option<TopoShape>,
    /// `myStatus` (`hxx:155`), filled by `BRepOffset::Surface` (`cxx:524`).
    my_status: BRepOffsetStatus,
    /// `myFace` (`hxx:156`).
    my_face: Option<Face>,
}

impl Default for BRepOffsetOffset {
    /// `BRepOffset_Offset::BRepOffset_Offset()` (`cxx:389`): every member
    /// default-initialized. `myStatus` is left indeterminate by OCCT's
    /// defaulted constructor; the port starts at `Good`.
    fn default() -> Self {
        Self {
            my_shape: None,
            my_status: BRepOffsetStatus::Good,
            my_face: None,
        }
    }
}

impl BRepOffsetOffset {
    pub fn new() -> Self {
        Self::default()
    }

    /// `BRepOffset_Offset::InitialShape()` (`.lxx:22-25`).
    pub fn initial_shape(&self) -> Option<&TopoShape> {
        self.my_shape.as_ref()
    }

    /// `BRepOffset_Offset::Face()` (`cxx:1728-1731`).
    pub fn face(&self) -> Option<&Face> {
        self.my_face.as_ref()
    }

    /// `BRepOffset_Offset::Status()` (`cxx:1811-1814`).
    pub fn status(&self) -> BRepOffsetStatus {
        self.my_status
    }

    /// `BRepOffset_Offset::Generated(const TopoDS_Shape&)`
    /// (`cxx:1735-1807`).
    ///
    /// FACE arm (`cxx:1741-1759`): walk the source face (`FORWARD`) and the
    /// offset face (`FORWARD`) edge lists in lockstep and answer the offset
    /// edge that pairs with `shape`; a `REVERSED` source face reverses the
    /// answer. EDGE arm (`cxx:1762-1801`): only the first wire of `myFace`
    /// is read, `shape` is matched against the source edge's **last**
    /// vertex (`cxx:1773`) and the third or fourth edge of that wire is
    /// answered (`cxx:1770-1787`); a `REVERSED` `myFace` reverses it
    /// (`cxx:1796-1799`). Every other shape type answers a null shape
    /// (`cxx:1802-1803`).
    pub fn generated(&self, shape: &TopoShape) -> Option<TopoShape> {
        let my_shape = self.my_shape.as_ref()?;
        match my_shape.shape_type() {
            ShapeType::Face => {
                let my_face = self.my_face.as_ref()?;
                let src = Face(my_shape.oriented(Orientation::Forward));
                let dst = Face(my_face.0.oriented(Orientation::Forward));
                let src_edges = explorer_shapes(&src.0, ShapeType::Edge);
                let dst_edges = explorer_shapes(&dst.0, ShapeType::Edge);
                for (a, b) in src_edges.iter().zip(dst_edges.iter()) {
                    if !is_same(a, shape) {
                        continue;
                    }
                    return Some(if my_shape.orientation() == Orientation::Reversed {
                        b.oriented(Orientation::Reversed)
                    } else {
                        b.clone()
                    });
                }
                None
            }
            ShapeType::Edge => {
                let my_face = self.my_face.as_ref()?;
                let (_, v2) = BRepOffsetTool::edge_vertices(&Edge(my_shape.clone()));
                let v2 = v2?;
                let fwd = Face(my_face.0.oriented(Orientation::Forward));
                let first_wire = wires_of_face(&fwd).into_iter().next()?;
                let wire = Wire(first_wire.0.oriented(Orientation::Forward));
                let edges = edges_of_wire(&wire);
                let mut it = edges.iter().skip(2);
                let answered = if is_same(&v2.0, shape) {
                    it.next()?
                } else {
                    it.next()?;
                    it.next()?
                };
                Some(if my_face.0.orientation() == Orientation::Reversed {
                    answered.0.oriented(Orientation::Reversed)
                } else {
                    answered.0.clone()
                })
            }
            _ => None,
        }
    }
}

/// `NCollection_DataMap<TopoDS_Shape, TopoDS_Shape, TopTools_ShapeMapHasher>`:
/// the `Created` map and `MapSS` (`cxx:822`).
pub(crate) type MapSs = ShapeMap;
