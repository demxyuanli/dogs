//! The edge-orientation helpers of `BRepOffset_Inter2d.cxx`:
//! `CommonVertex` (`:83-97`), `DefineClosedness` (`:101-124`) and
//! `GetEdgesOrientedInFace` (`:128-266`).
//!
//! `GetEdgesOrientedInFace` orders the edges of one face into a connected
//! chain: it starts at a vertex carried by a single edge, then walks by shared
//! vertices. A *closed* set of edges (every vertex has two incident edges) has
//! no such start vertex, so the chain is seeded with the vertex whose two
//! incident pcurves differ most along the `DefineClosedness` coordinate and the
//! walk goes on until it returns to that seed.
//!
//! Port note: two null `TopoDS_Shape`s are `IsSame` in OCCT, so the reference
//! can build a one-element sequence holding a null edge when no start vertex is
//! found on a non-closed set (`:184-232`). The port answers the empty `Vec`
//! there: a `Vec<Edge>` cannot hold a null edge, and every reader of the
//! sequence immediately casts the entry to `TopoDS::Edge` and dereferences it.

use std::collections::HashMap;

use crate::abs::{Orientation, ShapeType};
use crate::brep_adaptor_curve2d::Curve2dAdaptor;
use crate::brep_algo::AsDes;
use crate::brep_tool::BRepTool;
use crate::brep_tools::is_really_closed;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::topo_tools::{is_same, shape_map_key, ShapeMapKey};

/// `CommonVertex(E1, E2)` (`cxx:83-97`): `TopExp::Vertices(E1, V1, V2, true)`
/// then `V1[1]` before `V1[0]`, i.e. the *last* endpoint of `E1` wins.
///
/// Returns the shared vertex, or `None` for the null `V` of the reference
/// (`IsSame` of two null shapes never reaches the port because the endpoint
/// checks are on `Option`).
pub(crate) fn common_vertex(e1: &Edge, e2: &Edge) -> Option<Vertex> {
    let v1 = cumulative_vertices(e1);
    let v2 = cumulative_vertices(e2);
    if let Some(last1) = &v1[1] {
        if v2.iter().flatten().any(|v| is_same(&last1.0, &v.0)) {
            return Some(last1.clone());
        }
    }
    if let Some(first1) = &v1[0] {
        if v2.iter().flatten().any(|v| is_same(&first1.0, &v.0)) {
            return Some(first1.clone());
        }
    }
    None
}

/// `TopExp::Vertices(E, V1, V2, CumOri = true)` (`TopExp.cxx:168-...`): the
/// edge's two vertices, swapped when the edge itself is `REVERSED`.
pub(crate) fn cumulative_vertices(e: &Edge) -> [Option<Vertex>; 2] {
    let (first, last) = crate::topo_tools::edge_vertices(e);
    if e.0.orientation() == Orientation::Reversed {
        [last, first]
    } else {
        [first, last]
    }
}

/// `TopExp::FirstVertex(E, CumOri = true)`.
pub(crate) fn cumulative_first_vertex(e: &Edge) -> Option<Vertex> {
    cumulative_vertices(e)[0].clone()
}

/// `TopExp::LastVertex(E, CumOri = true)`.
pub(crate) fn cumulative_last_vertex(e: &Edge) -> Option<Vertex> {
    cumulative_vertices(e)[1].clone()
}

/// `DefineClosedness(theFace)` (`cxx:101-124`): `1` when the face's first
/// really-closed edge has a pcurve tangent parallel to `OY` at its first
/// parameter, `2` when parallel to `OX`, `0` when the face has no really-closed
/// edge.
///
/// The pcurve is read with `BRep_Tool::CurveOnSurface(anEdge, theFace, fpar,
/// lpar)` and the tangent with `DN(fpar, 1)`, then crossed with `gp::DX2d()`
/// and `gp::DY2d()`: a tangent parallel to `OY` makes the `DX2d` cross product
/// the LARGER of the two, hence the `1` / `2` split.
pub(crate) fn define_closedness(face: &Face) -> i32 {
    for shape in crate::topo_tools::explorer_shapes(&face.0, ShapeType::Edge) {
        let edge = Edge(shape);
        if !is_really_closed(&edge, face) {
            continue;
        }
        let Some((pc, fpar, _lpar)) = crate::boptools_2d::curve_on_surface_range(&edge, face) else {
            // Port guard: `aPCurve->DN` at `cxx:111` faults on a null handle.
            continue;
        };
        let tangent = pc.eval_dn(fpar, 1);
        let cross1 = tangent.crossed(&occt_core::gp::GpVec2d::new(1.0, 0.0));
        let cross2 = tangent.crossed(&occt_core::gp::GpVec2d::new(0.0, 1.0));
        return if cross2.abs() < cross1.abs() { 1 } else { 2 };
    }
    0
}

/// `GetEdgesOrientedInFace(theShape, theFace, theAsDes, theSeqEdges)`
/// (`cxx:128-266`): the edges of `theShape` that also lie on `theFace`,
/// ordered into a chain.
///
/// `theSeqEdges` is an out-parameter and is *appended* to, as in the reference;
/// the caller passes an empty sequence.
pub(crate) fn get_edges_oriented_in_face(
    shape: &TopoShape,
    face: &Face,
    as_des: &AsDes,
    seq_edges: &mut Vec<Edge>,
) {
    // `cxx:134-149`: keep the AsDes occurrence (orientation included) of every
    // edge of `theShape` that is `IsSame` to an edge of the face.
    let face_edges = as_des.descendant(&face.0);
    for shape_edge in crate::topo_tools::explorer_shapes(shape, ShapeType::Edge) {
        for edge_in_face in face_edges {
            if is_same(edge_in_face, &shape_edge) {
                seq_edges.push(Edge(edge_in_face.clone()));
                break;
            }
        }
    }

    // `cxx:151-154`.
    if seq_edges.len() == 1 {
        return;
    }

    // `cxx:156-161`: `TopExp::MapShapesAndAncestors(seqEdges(ii), VERTEX, EDGE,
    // aVEmap)` accumulated over every sequence entry, so the map is an
    // `IndexedDataMap` ordered by first insertion (`TopExp.cxx:80-120`).
    let mut vemap = IndexedAncestorMap::new();
    for e in seq_edges.iter() {
        vemap.add_edge(e);
    }

    // `cxx:163-181`: the first vertex carried by exactly one edge, taken when
    // it is that edge's *first* vertex.
    let mut first_vertex: Option<Vertex> = None;
    let mut first_edge: Option<Edge> = None;
    for (vertex, e_list) in vemap.iter() {
        if e_list.len() == 1 {
            let edge = Edge(e_list[0].clone());
            let v1 = cumulative_first_vertex(&edge);
            if let Some(v1) = v1 {
                if is_same(&v1.0, &vertex.0) {
                    first_vertex = Some(vertex.clone());
                    first_edge = Some(edge);
                    break;
                }
            }
        }
    }

    if first_edge.is_none() {
        // `cxx:183-232`: a closed set of edges.
        let ind_coord = define_closedness(face);
        if ind_coord != 0 {
            let mut max_delta = 0.0f64;
            for (vertex, e_list) in vemap.iter() {
                // `AFirst/Last()` of the reference only exist for a non-empty
                // list; `MapShapesAndAncestors` binds every vertex it meets to
                // at least one edge.
                let (Some(edge1), Some(edge2)) = (e_list.first(), e_list.last()) else {
                    continue;
                };
                let edge1 = Edge(edge1.clone());
                let edge2 = Edge(edge2.clone());
                // `BRep_Tool::Parameter(aVertex, anEdge)`.
                let param1 = BRepTool::parameter_of_vertex(&edge1.0, &vertex.0);
                let param2 = BRepTool::parameter_of_vertex(&edge2.0, &vertex.0);
                // `BRepAdaptor_Curve2d(anEdge, theFace).Value(aParam)`.
                let Some(p1) = Curve2dAdaptor::of_edge_face(&edge1, face).d0(param1) else {
                    continue;
                };
                let Some(p2) = Curve2dAdaptor::of_edge_face(&edge2, face).d0(param2) else {
                    continue;
                };
                let delta = if ind_coord == 1 {
                    (p1.x() - p2.x()).abs()
                } else {
                    (p1.y() - p2.y()).abs()
                };
                if delta > max_delta {
                    max_delta = delta;
                    first_vertex = Some(vertex.clone());
                }
            }
            // `cxx:219-230`: the edge of that vertex that starts there.
            if let Some(fv) = first_vertex.as_ref() {
                if let Some(e_list) = vemap.get(fv) {
                    for e in e_list {
                        let edge = Edge(e.clone());
                        if let Some(v1) = cumulative_first_vertex(&edge) {
                            if is_same(&v1.0, &fv.0) {
                                first_edge = Some(edge);
                                break;
                            }
                        }
                    }
                }
            }
        }
    }

    // `cxx:234-239`.
    let nb_edges = seq_edges.len();
    seq_edges.clear();
    let Some(mut edge) = first_edge else {
        // Port guard: the reference appends the null `aFirstEdge` here
        // (`cxx:236`) and breaks out of the walk on `null.IsSame(null)`.
        return;
    };
    seq_edges.push(edge.clone());

    // `cxx:240-268`.
    loop {
        let last_vertex = cumulative_last_vertex(&edge);
        // `aLastVertex.IsSame(aFirstVertex)`: two null vertices are `IsSame`
        // too, so a missing endpoint on either side ends the walk.
        let stop = match (last_vertex.as_ref(), first_vertex.as_ref()) {
            (Some(a), Some(b)) => is_same(&a.0, &b.0),
            (None, None) => true,
            _ => false,
        };
        if stop {
            break;
        }
        let Some(last_vertex) = last_vertex else {
            break;
        };
        // `cxx:246-249`: `FindFromKey` raises for an unbound key; the port
        // guards it.
        let Some(e_list) = vemap.get(&last_vertex) else {
            break;
        };
        if e_list.len() == 1 {
            break;
        }
        // `cxx:251-258`.
        edge = if is_same(&e_list[0], &edge.0) {
            Edge(e_list[e_list.len() - 1].clone())
        } else {
            Edge(e_list[0].clone())
        };
        // `cxx:260-265`.
        seq_edges.push(edge.clone());
        if seq_edges.len() == nb_edges {
            break;
        }
    }
}

/// `NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>>`
/// restricted to what `GetEdgesOrientedInFace` needs: insertion-ordered keys
/// and per-key append-only lists. `NCollection_List` order is preserved, hence
/// the `Vec`.
struct IndexedAncestorMap {
    keys: Vec<Vertex>,
    index: HashMap<ShapeMapKey, usize>,
    lists: Vec<Vec<TopoShape>>,
}

impl IndexedAncestorMap {
    fn new() -> Self {
        Self {
            keys: Vec::new(),
            index: HashMap::new(),
            lists: Vec::new(),
        }
    }

    /// One `TopExp::MapShapesAndAncestors(edge, VERTEX, EDGE, map)` pass: each
    /// vertex of the edge gets the edge appended as an ancestor.
    fn add_edge(&mut self, edge: &Edge) {
        for v in crate::topo_tools::vertices_of(&edge.0) {
            let slot = match self.index.get(&shape_map_key(&v.0)) {
                Some(slot) => *slot,
                None => {
                    let slot = self.keys.len();
                    self.index.insert(shape_map_key(&v.0), slot);
                    self.keys.push(v);
                    self.lists.push(Vec::new());
                    slot
                }
            };
            self.lists[slot].push(edge.0.clone());
        }
    }

    fn iter(&self) -> impl Iterator<Item = (&Vertex, &Vec<TopoShape>)> {
        self.keys.iter().zip(self.lists.iter())
    }

    fn get(&self, vertex: &Vertex) -> Option<&Vec<TopoShape>> {
        self.index
            .get(&shape_map_key(&vertex.0))
            .map(|slot| &self.lists[*slot])
    }
}

