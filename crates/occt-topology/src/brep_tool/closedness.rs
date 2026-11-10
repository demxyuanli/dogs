use crate::abs::{Orientation, ShapeType};
use crate::shape::{Edge, Face, TopoShape};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{explorer_shapes, is_same};

use super::BRepTool;

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
}
