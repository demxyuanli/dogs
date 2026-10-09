//! Port of `BRepClass_Edge` (`BRepClass_Edge.cxx`, `BRepClass_Edge.hxx`):
//! an edge of a face, the successor edge at its last vertex, and the
//! tolerance / bounding-box switches read by `BRepClass_Intersector`.

use occt_core::precision::INFINITE;

use crate::abs::Orientation;
use crate::shape::{Edge, Face, Vertex};
use crate::topo_tools_full::edge_vertices;

/// `BRepClass_Edge`.
#[derive(Clone, Debug)]
pub struct BRepClassEdge {
    edge: Option<Edge>,
    face: Option<Face>,
    next_edge: Option<Edge>,
    max_tolerance: f64,
    use_bnd_box: bool,
}

impl Default for BRepClassEdge {
    fn default() -> Self {
        Self::new()
    }
}

impl BRepClassEdge {
    /// `BRepClass_Edge()`: null edge and face, `MaxTolerance` infinite,
    /// bounding-box prefilter off.
    pub fn new() -> Self {
        Self {
            edge: None,
            face: None,
            next_edge: None,
            max_tolerance: INFINITE,
            use_bnd_box: false,
        }
    }

    /// `BRepClass_Edge(E, F)`.
    pub fn from_edge_face(edge: Edge, face: Face) -> Self {
        Self {
            edge: Some(edge),
            face: Some(face),
            ..Self::new()
        }
    }

    /// `Edge()`; `None` stands for a null edge.
    pub fn edge(&self) -> Option<&Edge> {
        self.edge.as_ref()
    }

    /// `Face()`; `None` stands for a null face.
    pub fn face(&self) -> Option<&Face> {
        self.face.as_ref()
    }

    /// `NextEdge()`; `None` stands for a null edge.
    pub fn next_edge(&self) -> Option<&Edge> {
        self.next_edge.as_ref()
    }

    /// `MaxTolerance()`.
    pub fn max_tolerance(&self) -> f64 {
        self.max_tolerance
    }

    /// `SetMaxTolerance(theValue)`.
    pub fn set_max_tolerance(&mut self, value: f64) {
        self.max_tolerance = value;
    }

    /// `UseBndBox()`.
    pub fn use_bnd_box(&self) -> bool {
        self.use_bnd_box
    }

    /// `SetUseBndBox(theValue)`.
    pub fn set_use_bnd_box(&mut self, value: bool) {
        self.use_bnd_box = value;
    }

    /// `SetNextEdge(theMapVE)` (`BRepClass_Edge.cxx:31-61`). `map_ve` is the
    /// vertex-to-edges map of `TopExp::MapShapesAndAncestors`, given as
    /// pairs. The successor is set only when the last vertex is shared by
    /// exactly two edges; when `Seek` finds no entry OCCT dereferences a null
    /// list pointer, which is treated here as "no successor".
    pub fn set_next_edge(&mut self, map_ve: &[(Vertex, Vec<Edge>)]) {
        let Some(edge) = self.edge.as_ref() else {
            return;
        };
        if map_ve.is_empty() {
            return;
        }
        // `TopExp::Vertices(myEdge, aVF, aVL, true)`: cumulative orientation.
        let (vf, vl) = vertices_cum_ori(edge);
        let Some(vl) = vl else {
            return;
        };
        if let Some(vf) = vf.as_ref() {
            if vl.same_tshape(vf) {
                return;
            }
        }
        let Some((_, list)) = map_ve.iter().find(|(v, _)| v.same_tshape(&vl)) else {
            return;
        };
        if list.len() == 2 {
            for e in list {
                if !e.same_tshape(edge) {
                    self.next_edge = Some(e.clone());
                }
            }
        }
    }
}

/// `TopExp::Vertices(E, V1, V2, CumOri = true)`: the first and last vertex of
/// `edge` with the edge orientation composed in. A reversed edge swaps the
/// two stored endpoints. The plain (`CumOri = false`) pair comes from
/// `topo_tools_full::edge_vertices`.
pub(super) fn vertices_cum_ori(edge: &Edge) -> (Option<Vertex>, Option<Vertex>) {
    let (first, last) = edge_vertices(edge);
    if edge.orientation() == Orientation::Reversed {
        (last, first)
    } else {
        (first, last)
    }
}
