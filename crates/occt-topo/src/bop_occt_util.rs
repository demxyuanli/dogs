//! Shared helpers for the OCCT Builder translation modules.
//!
//! The Builder stages in `BOPAlgo_Builder_2.cxx` / `_3.cxx` and
//! `BOPAlgo_Tools.cxx` share a small set of identity, fence, iterator and
//! orientation helpers. This module is the Rust equivalent of those locals:
//! `TopTools_MapOfShape` fence maps, `TopoDS_Iterator` child walks,
//! `TopExp::MapShapes` typed collectors, `BRep_Builder` add wrappers, and
//! the `myShapesSD.IsBound` / `myImages.IsBound` lookups the solid and face
//! stages both need.
//!
//! Identity is [`GeometryRegistry::shape_key`] (TShape pointer), matching
//! `TopTools_ShapeMapHasher` on the TShape handle. Do not mix this key with
//! `Arc::as_ptr` aliases used in a few older solid-stage helpers.

use std::collections::{HashMap, HashSet};

use crate::abs::{Orientation, ShapeType};
use crate::bop_hist::BopHistory;
use crate::bopds::BopdsDS;
use crate::builder::TopoBuilder;
use crate::iterator::cumulated_children;
use crate::shape::{Edge, Face, Shell, Solid, TopoShape, Vertex, Wire};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools_full::{edges_of, faces_of, vertices_of};

/// TShape identity used as a `TopTools_MapOfShape` key.
pub fn shape_key(s: &TopoShape) -> usize {
    GeometryRegistry::shape_key(s)
}

/// True when `a` and `b` share a TShape (`IsSame` without location).
pub fn same_shape(a: &TopoShape, b: &TopoShape) -> bool {
    a.same_tshape(b)
}

/// Ordered fence map: insert returns true when `s` was new.
pub fn fence_add(fence: &mut Vec<TopoShape>, s: TopoShape) -> bool {
    if fence.iter().any(|x| x.same_tshape(&s)) {
        false
    } else {
        fence.push(s);
        true
    }
}

/// Hash-set fence keyed by [`shape_key`].
pub fn fence_insert(fence: &mut HashSet<usize>, s: &TopoShape) -> bool {
    fence.insert(shape_key(s))
}

/// Direct children with orientation/location composition (`TopoDS_Iterator`).
pub fn iter_children(s: &TopoShape) -> Vec<TopoShape> {
    cumulated_children(s)
}

/// Direct children whose type equals `ty`.
pub fn iter_typed(s: &TopoShape, ty: ShapeType) -> Vec<TopoShape> {
    iter_children(s)
        .into_iter()
        .filter(|c| c.shape_type() == ty)
        .collect()
}

/// Recursive explorer (`TopExp_Explorer`) collecting every sub-shape of `ty`.
pub fn explore(s: &TopoShape, ty: ShapeType) -> Vec<TopoShape> {
    match ty {
        ShapeType::Vertex => vertices_of(s).into_iter().map(|v| v.0).collect(),
        ShapeType::Edge => edges_of(s).into_iter().map(|e| e.0).collect(),
        ShapeType::Face => faces_of(s).into_iter().map(|f| f.0).collect(),
        ShapeType::Wire => collect_type(s, ShapeType::Wire),
        ShapeType::Shell => collect_type(s, ShapeType::Shell),
        ShapeType::Solid => collect_type(s, ShapeType::Solid),
        _ => collect_type(s, ty),
    }
}

fn collect_type(s: &TopoShape, ty: ShapeType) -> Vec<TopoShape> {
    let mut out = Vec::new();
    let mut seen: HashSet<usize> = HashSet::new();
    fn walk(s: &TopoShape, ty: ShapeType, out: &mut Vec<TopoShape>, seen: &mut HashSet<usize>) {
        if s.shape_type() == ty && seen.insert(shape_key(s)) {
            out.push(s.clone());
        }
        for c in iter_children(s) {
            walk(&c, ty, out, seen);
        }
    }
    walk(s, ty, &mut out, &mut seen);
    out
}

/// `TopExp::MapShapes(s, type, map)` — unique TShape set of that type.
pub fn map_shapes(s: &TopoShape, ty: ShapeType) -> HashSet<usize> {
    explore(s, ty).into_iter().map(|x| shape_key(&x)).collect()
}

/// `TopExp::MapShapesAndAncestors(s, sub, anc, map)`.
pub fn map_shapes_and_ancestors(
    s: &TopoShape,
    sub: ShapeType,
    anc: ShapeType,
) -> HashMap<usize, Vec<TopoShape>> {
    let mut out: HashMap<usize, Vec<TopoShape>> = HashMap::new();
    for a in explore(s, anc) {
        for k in map_shapes(&a, sub) {
            let list = out.entry(k).or_default();
            if !list.iter().any(|x| x.same_tshape(&a)) {
                list.push(a.clone());
            }
        }
    }
    out
}

/// Whether the history table has an image list for `s` (`myImages.IsBound`).
pub fn images_bound(h: &BopHistory, s: &TopoShape) -> bool {
    h.has_image(s)
}

/// Image pieces of `s`, or an empty slice when unbound.
pub fn images_of<'a>(h: &'a BopHistory, s: &TopoShape) -> &'a [TopoShape] {
    h.image(s).unwrap_or(&[])
}

/// Reverse a shape in place (`TopoDS_Shape::Reverse`).
pub fn reverse_shape(s: &mut TopoShape) {
    s.reverse();
}

/// Copy of `s` with orientation `or`.
pub fn oriented(s: &TopoShape, or: Orientation) -> TopoShape {
    s.oriented(or)
}

/// Copy of `s` with `FORWARD` orientation.
pub fn forward_copy(s: &TopoShape) -> TopoShape {
    s.oriented(Orientation::Forward)
}

/// Copy of `s` with `REVERSED` orientation.
pub fn reversed_copy(s: &TopoShape) -> TopoShape {
    s.oriented(Orientation::Reversed)
}

/// Copy of `s` with `INTERNAL` orientation.
pub fn internal_copy(s: &TopoShape) -> TopoShape {
    s.oriented(Orientation::Internal)
}

/// Face-info record for DS index `i`, if the pool holds one (`HasFaceInfo`).
pub fn face_info_of(ds: &BopdsDS, i: usize) -> Option<&crate::bopds::BopdsFaceInfo> {
    ds.face_info_pool().iter().find(|fi| fi.face_index == i)
}

/// True when the DS records a face-info entry for source index `i`.
pub fn has_face_info(ds: &BopdsDS, i: usize) -> bool {
    face_info_of(ds, i).is_some()
}

/// Alone vertices of a face (`BOPDS_DS::AloneVertices`): face vertices that
/// belong to no edge of that face.
pub fn alone_vertices(ds: &BopdsDS, face_idx: usize) -> Vec<usize> {
    let Some(si) = ds.shape_info(face_idx) else {
        return Vec::new();
    };
    let mut edge_vertices: HashSet<usize> = HashSet::new();
    let mut vertex_indices: Vec<usize> = Vec::new();
    for &sub in &si.sub_indices {
        let Some(sub_info) = ds.shape_info(sub) else { continue };
        match sub_info.shape_type() {
            ShapeType::Edge => {
                for &v in &sub_info.sub_indices {
                    edge_vertices.insert(v);
                }
            }
            ShapeType::Vertex => vertex_indices.push(sub),
            _ => {}
        }
    }
    vertex_indices
        .into_iter()
        .filter(|v| !edge_vertices.contains(v))
        .collect()
}

/// Same-domain index of a DS-indexed shape (`BOPDS_DS::HasShapeSD`).
pub fn ds_shape_sd(ds: &BopdsDS, s: &TopoShape) -> Option<usize> {
    ds.index(s).and_then(|i| ds.has_shape_sd(i))
}

/// Empty solid with the given orientation (`BRep_Builder::MakeSolid`).
pub fn make_empty_solid(or: Orientation) -> Solid {
    let mut s = Solid::new();
    s.0.set_orientation(or);
    s
}

/// Empty shell with the given orientation (`BRep_Builder::MakeShell`).
pub fn make_empty_shell(or: Orientation) -> Shell {
    let mut s = Shell::new();
    s.0.set_orientation(or);
    s
}

/// Empty wire (`BRep_Builder::MakeWire`).
pub fn make_empty_wire(or: Orientation) -> Wire {
    let mut w = Wire::new();
    w.0.set_orientation(or);
    w
}

/// Add `sub` to `parent` through [`TopoBuilder::add`].
pub fn builder_add(parent: &mut TopoShape, sub: &TopoShape) {
    TopoBuilder::new().add(parent, sub);
}

/// Add a face to a shell.
pub fn add_face_to_shell(shell: &mut Shell, face: &Face) {
    TopoBuilder::new().add_face(shell, face);
}

/// Add a shell to a solid.
pub fn add_shell_to_solid(solid: &mut Solid, shell: &Shell) {
    TopoBuilder::new().add_shell(solid, shell);
}

/// Add a wire to a face.
pub fn add_wire_to_face(face: &mut Face, wire: &Wire) {
    TopoBuilder::new().add_wire(face, wire);
}

/// Add an edge to a wire.
pub fn add_edge_to_wire(wire: &mut Wire, edge: &Edge) {
    TopoBuilder::new().add_edge(wire, edge);
}

/// Typed view of a shape as a face, when the type matches.
pub fn as_face(s: &TopoShape) -> Option<Face> {
    Face::wrap(s.clone())
}

/// Typed view of a shape as an edge, when the type matches.
pub fn as_edge(s: &TopoShape) -> Option<Edge> {
    Edge::wrap(s.clone())
}

/// Typed view of a shape as a solid, when the type matches.
pub fn as_solid(s: &TopoShape) -> Option<Solid> {
    Solid::wrap(s.clone())
}

/// Typed view of a shape as a shell, when the type matches.
pub fn as_shell(s: &TopoShape) -> Option<Shell> {
    Shell::wrap(s.clone())
}

/// Typed view of a shape as a wire, when the type matches.
pub fn as_wire(s: &TopoShape) -> Option<Wire> {
    Wire::wrap(s.clone())
}

/// Typed view of a shape as a vertex, when the type matches.
pub fn as_vertex(s: &TopoShape) -> Option<Vertex> {
    Vertex::wrap(s.clone())
}

/// Unique `(edge, first, last)` tuples from a pave-block list, dropping
/// non-edge DS entries. Two blocks of the same original keep both ranges.
pub fn unique_edge_paves(ds: &BopdsDS, src: &[(usize, f64, f64)]) -> Vec<(usize, f64, f64)> {
    let mut dst = Vec::new();
    for &(e, fl, ll) in src {
        let is_edge = ds
            .shape(e)
            .map(|s| s.shape_type() == ShapeType::Edge)
            .unwrap_or(false);
        if is_edge
            && !dst.iter().any(|t: &(usize, f64, f64)| {
                t.0 == e && (t.1 - fl).abs() <= 1e-7 && (t.2 - ll).abs() <= 1e-7
            })
        {
            dst.push((e, fl, ll));
        }
    }
    dst
}

/// Split edge of a pave block: the block of `e_idx` whose range matches
/// `[fl, ll]` yields `pb.edge()`; otherwise the whole edge.
pub fn on_face_split_edge(ds: &BopdsDS, e_idx: usize, fl: f64, ll: f64) -> Option<TopoShape> {
    let tol = 1e-7;
    let blocks = ds.pave_blocks(e_idx);
    if let Some(pb) = blocks.iter().find(|pb| {
        let (a, b) = pb.range();
        ((a - fl).abs() <= tol && (b - ll).abs() <= tol)
            || ((a - ll).abs() <= tol && (b - fl).abs() <= tol)
    }) {
        if let Some(sp) = ds.shape(pb.edge()).cloned() {
            return Some(sp);
        }
    }
    ds.shape(e_idx).cloned()
}

/// Append both orientations of each on-face split edge to `le`.
pub fn append_on_face_paves(ds: &BopdsDS, paves: &[(usize, f64, f64)], le: &mut Vec<Edge>) {
    for &(e_idx, fl, ll) in paves {
        let Some(mut sp) = on_face_split_edge(ds, e_idx, fl, ll) else {
            continue;
        };
        sp.set_orientation(Orientation::Forward);
        le.push(Edge(sp.clone()));
        sp.set_orientation(Orientation::Reversed);
        le.push(Edge(sp));
    }
}

/// Source solids of the DS, in index order.
pub fn source_solids(ds: &BopdsDS) -> Vec<TopoShape> {
    let n = ds.nb_source_shapes();
    let mut out = Vec::new();
    for i in 0..n {
        let Some(si) = ds.shape_info(i) else { continue };
        if si.shape_type() == ShapeType::Solid {
            out.push(si.shape().clone());
        }
    }
    out
}

/// Source faces of the DS, in index order.
pub fn source_faces(ds: &BopdsDS) -> Vec<(usize, TopoShape)> {
    let n = ds.nb_source_shapes();
    let mut out = Vec::new();
    for i in 0..n {
        let Some(si) = ds.shape_info(i) else { continue };
        if si.shape_type() == ShapeType::Face {
            out.push((i, si.shape().clone()));
        }
    }
    out
}

/// Candidate faces of FillIn3DParts: every source FACE replaced by its image
/// splits, an un-split face kept as-is, fenced by TShape identity.
pub fn collect_candidate_faces(h: &BopHistory, ds: &BopdsDS) -> Vec<TopoShape> {
    let mut fence: HashSet<usize> = HashSet::new();
    let mut out = Vec::new();
    for (_, face) in source_faces(ds) {
        if images_bound(h, &face) {
            for im in images_of(h, &face) {
                if fence.insert(shape_key(im)) {
                    out.push(im.clone());
                }
            }
        } else if fence.insert(shape_key(&face)) {
            out.push(face);
        }
    }
    out
}

/// Own internal (non-shell) children of a solid (`OwnInternalShapes`).
pub fn own_internal_shapes(solid: &TopoShape) -> Vec<TopoShape> {
    iter_children(solid)
        .into_iter()
        .filter(|c| c.shape_type() != ShapeType::Shell)
        .collect()
}

/// Flatten compounds into a list (`BOPTools_AlgoTools::TreatCompound`).
pub fn treat_compound(s: &TopoShape, out: &mut Vec<TopoShape>, fence: &mut HashSet<usize>) {
    if s.shape_type() == ShapeType::Compound {
        for c in iter_children(s) {
            treat_compound(&c, out, fence);
        }
        return;
    }
    if fence.insert(shape_key(s)) {
        out.push(s.clone());
    }
}

/// Whether a wire's first and last vertices coincide (`BRep_Tool::IsClosed`).
pub fn wire_is_closed(wire: &Wire) -> bool {
    let edges = crate::topo_tools_full::edges_of_wire(wire);
    if edges.is_empty() {
        return false;
    }
    let (first_a, _) = crate::topo_tools_full::edge_vertices(&edges[0]);
    let last = edges.last().unwrap();
    let (_, last_b) = crate::topo_tools_full::edge_vertices(last);
    match (first_a, last_b) {
        (Some(a), Some(b)) => a.0.same_tshape(&b.0),
        _ => false,
    }
}

/// Merge two ancestor maps, appending unique ancestors per key.
pub fn merge_ancestor_maps(
    dst: &mut HashMap<usize, Vec<TopoShape>>,
    src: HashMap<usize, Vec<TopoShape>>,
) {
    for (k, vals) in src {
        let list = dst.entry(k).or_default();
        for v in vals {
            if !list.iter().any(|x| x.same_tshape(&v)) {
                list.push(v);
            }
        }
    }
}

/// Origins back-map insert (`myOrigins.Bound` / `ChangeSeek`).
pub fn origins_append(
    origins: &mut HashMap<usize, Vec<TopoShape>>,
    split: &TopoShape,
    source: TopoShape,
) {
    origins
        .entry(shape_key(split))
        .or_default()
        .push(source);
}

/// Deduplicate a shape list by TShape identity, preserving order.
pub fn unique_shapes(shapes: &[TopoShape]) -> Vec<TopoShape> {
    let mut fence: HashSet<usize> = HashSet::new();
    let mut out = Vec::new();
    for s in shapes {
        if fence.insert(shape_key(s)) {
            out.push(s.clone());
        }
    }
    out
}

/// Pair of orientations (FORWARD, REVERSED) of `s`.
pub fn both_orientations(s: &TopoShape) -> [TopoShape; 2] {
    [forward_copy(s), reversed_copy(s)]
}

/// Whether `or` is INTERNAL.
pub fn is_internal(or: Orientation) -> bool {
    or == Orientation::Internal
}

/// Whether `or` is REVERSED.
pub fn is_reversed(or: Orientation) -> bool {
    or == Orientation::Reversed
}

/// Whether `or` is FORWARD.
pub fn is_forward(or: Orientation) -> bool {
    or == Orientation::Forward
}

/// Bounding-box of a shape from the DS when indexed, otherwise computed.
pub fn shape_box_of(ds: &BopdsDS, s: &TopoShape) -> occt_core::bnd::BndBox {
    if let Some(i) = ds.index(s) {
        if let Some(b) = ds.box_of(i) {
            if !b.is_void() {
                return b.clone();
            }
        }
    }
    crate::bbox_from_geometry::shape_bbox(s)
}

/// Copy a list of shapes.
pub fn clone_list(src: &[TopoShape]) -> Vec<TopoShape> {
    src.to_vec()
}

/// True when `list` contains a shape with the same TShape as `s`.
pub fn list_contains(list: &[TopoShape], s: &TopoShape) -> bool {
    list.iter().any(|x| x.same_tshape(s))
}

/// Append `s` to `list` if not already present by TShape identity.
pub fn list_append_unique(list: &mut Vec<TopoShape>, s: TopoShape) -> bool {
    if list_contains(list, &s) {
        false
    } else {
        list.push(s);
        true
    }
}

/// Wires that are direct children of a face.
pub fn wires_of_face(face: &Face) -> Vec<Wire> {
    iter_typed(&face.0, ShapeType::Wire)
        .into_iter()
        .filter_map(|w| Wire::wrap(w))
        .collect()
}

/// Shells that are direct children of a solid.
pub fn shells_of_solid(solid: &TopoShape) -> Vec<Shell> {
    iter_typed(solid, ShapeType::Shell)
        .into_iter()
        .filter_map(|s| Shell::wrap(s))
        .collect()
}

/// Faces that are direct children of a shell.
pub fn faces_of_shell(shell: &TopoShape) -> Vec<Face> {
    iter_typed(shell, ShapeType::Face)
        .into_iter()
        .filter_map(|f| Face::wrap(f))
        .collect()
}

/// Edges that are direct children of a wire.
pub fn edges_of_wire_direct(wire: &TopoShape) -> Vec<Edge> {
    iter_typed(wire, ShapeType::Edge)
        .into_iter()
        .filter_map(|e| Edge::wrap(e))
        .collect()
}

/// Vertices that are direct children of an edge.
pub fn vertices_of_edge_direct(edge: &TopoShape) -> Vec<Vertex> {
    iter_typed(edge, ShapeType::Vertex)
        .into_iter()
        .filter_map(|v| Vertex::wrap(v))
        .collect()
}

/// Image of an edge, or the edge itself when unbound.
pub fn edge_or_images(h: &BopHistory, e: &TopoShape) -> Vec<TopoShape> {
    if images_bound(h, e) {
        images_of(h, e).to_vec()
    } else {
        vec![e.clone()]
    }
}

/// Image of a face, or the face itself when unbound.
pub fn face_or_images(h: &BopHistory, f: &TopoShape) -> Vec<TopoShape> {
    if images_bound(h, f) {
        images_of(h, f).to_vec()
    } else {
        vec![f.clone()]
    }
}

/// Image of a solid, or the solid itself when unbound.
pub fn solid_or_images(h: &BopHistory, s: &TopoShape) -> Vec<TopoShape> {
    if images_bound(h, s) {
        images_of(h, s).to_vec()
    } else {
        vec![s.clone()]
    }
}

/// Record every image of `old` from `news`.
pub fn bind_images(h: &mut BopHistory, old: &TopoShape, news: &[TopoShape]) {
    for n in news {
        h.add_image(old, n.clone());
    }
}

/// DS index of a shape, if registered.
pub fn ds_index(ds: &BopdsDS, s: &TopoShape) -> Option<usize> {
    ds.index(s)
}

/// Shape stored at DS index `i`.
pub fn ds_shape(ds: &BopdsDS, i: usize) -> Option<TopoShape> {
    ds.shape(i).cloned()
}

/// Shape type at DS index `i`.
pub fn ds_type(ds: &BopdsDS, i: usize) -> Option<ShapeType> {
    ds.shape_info(i).map(|si| si.shape_type())
}

/// Number of source shapes.
pub fn nb_source(ds: &BopdsDS) -> usize {
    ds.nb_source_shapes()
}

/// Iterate source indices of a given type.
pub fn source_indices_of(ds: &BopdsDS, ty: ShapeType) -> Vec<usize> {
    let n = ds.nb_source_shapes();
    let mut out = Vec::new();
    for i in 0..n {
        if ds_type(ds, i) == Some(ty) {
            out.push(i);
        }
    }
    out
}

/// Host the solid-image stages need (`myDS`, `myImages`, `myOrigins`,
/// `myShapesSD`, `myReport`). Implemented by `BopBuilder` so this module
/// does not depend on that type.
pub trait BopSolidHost {
    fn ds(&self) -> &BopdsDS;
    fn history(&self) -> &BopHistory;
    fn history_mut(&mut self) -> &mut BopHistory;
    fn fuzzy_value(&self) -> f64;
    fn arguments(&self) -> &[TopoShape];
    fn origins_mut(&mut self) -> &mut HashMap<usize, Vec<TopoShape>>;
    fn add_warning(&mut self, msg: String);
    fn seek_shapes_sd(&self, shape: &TopoShape) -> Option<TopoShape>;
    fn bind_shapes_sd(&mut self, shape: TopoShape, sd: TopoShape);
}
