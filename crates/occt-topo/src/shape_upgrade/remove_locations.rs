//! `ShapeUpgrade_RemoveLocations` — rebuild a shape tree with the datum
//! transforms baked into the geometry.
//!
//! Source: `ShapeUpgrade_RemoveLocations.cxx:49-314`, `.hxx:36-77`,
//! `.lxx:17-48`.
//!
//! The class walks the tree child by child. A node whose location is not
//! identity (and whose type passes the removal level) is rebuilt with the
//! location folded into its own geometry: a face gets a transformed surface, an
//! edge a transformed 3D curve, a vertex a transformed point. Faces are rebuilt
//! before their edges, so an edge's pcurve can be re-attached to the *new*
//! surface of its ancestor face. Topology (children, orientation, `Closed`) is
//! copied; only the locations disappear.
//!
//! `BRepFill_Section` (`BRepFill_Section.cxx:46-51`) is the first caller: it
//! runs `SetRemoveLevel(TopAbs_COMPOUND)` + `Remove(profile)` so a section wire
//! reaching it under any location is baked to world coordinates.

use std::collections::HashMap;
use std::sync::Arc;

use occt_core::toploc::{LocationKey, TopLocLocation};
use occt_geom::Surface;
use occt_geom2d::curve::Curve2d;

use crate::abs::{Orientation, ShapeType};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::tgeometry::GeometryRegistry;

/// `(TShape identity, Location)` — the comparison unit of
/// `TopTools_ShapeMapHasher` / `TopoDS_Shape::IsSame`
/// (`TopoDS_Shape.hxx:262-267`), i.e. of `myMapNewShapes`.
type ShapeKey = (usize, LocationKey);

fn map_key(s: &TopoShape) -> ShapeKey {
    (s.shape_addr(), s.location().key())
}

/// `TopAbs_ShapeEnum` as an integer (`TopAbs.hxx:24-35`): `COMPOUND` is 0 and
/// `SHAPE` is 8. `ShapeType` declares its variants in that order.
fn level_index(t: ShapeType) -> i32 {
    t as i32
}

/// `ShapeUpgrade_RemoveLocations` (`ShapeUpgrade_RemoveLocations.hxx:36-77`).
#[derive(Debug)]
pub struct RemoveLocations {
    /// `myLevelRemoving`, default `TopAbs_SHAPE` (`cxx:40-43`).
    level_removing: ShapeType,
    /// `myShape` — the result of the last `Remove`.
    shape: Option<TopoShape>,
    /// `myMapNewShapes`: old shape (TShape + location) -> rebuilt shape.
    new_shapes: HashMap<ShapeKey, TopoShape>,
}

impl Default for RemoveLocations {
    fn default() -> Self {
        Self::new()
    }
}

impl RemoveLocations {
    pub fn new() -> Self {
        Self { level_removing: ShapeType::Shape, shape: None, new_shapes: HashMap::new() }
    }

    /// `SetRemoveLevel` (`.lxx:24-31`): a level coarser than `FACE` is clamped
    /// to `FACE`.
    pub fn set_remove_level(&mut self, level: ShapeType) {
        self.level_removing =
            if level_index(level) > level_index(ShapeType::Face) { ShapeType::Face } else { level };
    }

    /// `RemoveLevel` (`.lxx:35-38`).
    pub fn remove_level(&self) -> ShapeType {
        self.level_removing
    }

    /// `Remove` (`cxx:49-58`). Answers whether anything was rebuilt; the
    /// result is read through [`RemoveLocations::result`].
    pub fn remove(&mut self, the_shape: &TopoShape) -> bool {
        let shtype = the_shape.shape_type();
        let is_remove_loc = (shtype != ShapeType::Compound
            && self.level_removing == ShapeType::Shape)
            || level_index(self.level_removing) <= level_index(shtype);
        // `cxx:50-51`: `myShape = aShape` before the call; `cxx:311` overwrites
        // it with the rebuilt root (`theAncShape` is a null shape, `cxx:54`).
        let (is_done, new_shape) = self.make_new_shape(the_shape, None, is_remove_loc);
        self.shape = Some(new_shape);
        is_done
    }

    /// `GetResult` (`.lxx:17-20`).
    pub fn result(&self) -> Option<&TopoShape> {
        self.shape.as_ref()
    }

    /// `ModifiedShape` (`.lxx:41-48`): the rebuilt counterpart of a sub-shape,
    /// or the sub-shape itself when it was not recorded.
    pub fn modified_shape(&self, the_init_shape: &TopoShape) -> TopoShape {
        self.new_shapes.get(&map_key(the_init_shape)).cloned().unwrap_or_else(|| the_init_shape.clone())
    }

    /// `MakeNewShape` (`cxx:175-314`).
    fn make_new_shape(
        &mut self,
        the_shape: &TopoShape,
        the_anc: Option<&TopoShape>,
        the_remove_loc: bool,
    ) -> (bool, TopoShape) {
        let shtype = the_shape.shape_type();
        let builder = TopoBuilder::new();

        // `cxx:183-188`: when this node keeps its location, the *key* shape
        // drops it (it is re-applied to the result at `cxx:306-309`).
        let mut a_shape = the_shape.clone();
        if !the_remove_loc && !the_shape.location().is_identity() {
            a_shape.set_location(&TopLocLocation::identity());
        }

        // `cxx:189-204`.
        let key = map_key(&a_shape);
        let is_bound = self.new_shapes.contains_key(&key);
        let mut a_new_shape: Option<TopoShape> = None;
        if is_bound {
            let mut found = self.new_shapes[&key].clone();
            found.set_orientation(the_shape.orientation());
            if !the_remove_loc && !the_shape.location().is_identity() {
                found.set_location(the_shape.location());
            }
            if shtype != ShapeType::Edge {
                return (true, found);
            }
            a_new_shape = Some(found);
        }

        // `cxx:206-211`.
        let mut is_remove_loc = the_remove_loc;
        if !the_remove_loc {
            is_remove_loc = (shtype != ShapeType::Compound
                && self.level_removing == ShapeType::Shape)
                || level_index(self.level_removing) <= level_index(shtype);
        }

        // `cxx:213-218`.
        let mut a_rebuild = false;
        let mut an_anc: Option<TopoShape> = the_anc.cloned();
        if shtype == ShapeType::Face {
            an_anc = Some(a_shape.clone());
        }

        if is_remove_loc
            && (!a_shape.location().is_identity()
                || shtype == ShapeType::Edge
                || shtype == ShapeType::Face)
        {
            match shtype {
                // `cxx:224-233`.
                ShapeType::Face => {
                    let Some(old_face) = Face::wrap(a_shape.clone()) else {
                        unreachable!("shape_type() is Face")
                    };
                    if let Some(new_face) = rebuild_face(&old_face) {
                        a_rebuild = true;
                        self.new_shapes.insert(map_key(&old_face.0), new_face.0.clone());
                        a_new_shape = Some(new_face.0);
                    }
                }
                // `cxx:235-256`.
                ShapeType::Edge => {
                    let Some(old_edge) = Edge::wrap(a_shape.clone()) else {
                        unreachable!("shape_type() is Edge")
                    };
                    // `cxx:239-249`: the ancestor face and its rebuilt
                    // counterpart, when this edge hangs under a face at all.
                    let ancestor_face = an_anc.as_ref().and_then(|a| Face::wrap(a.clone()));
                    let new_face = ancestor_face.as_ref().map(|f| {
                        self.new_shapes
                            .get(&map_key(&f.0))
                            .and_then(|m| Face::wrap(m.clone()))
                            .unwrap_or_else(|| f.clone())
                    });
                    let mut anew_edge = match &a_new_shape {
                        Some(s) => Edge::wrap(s.clone()).unwrap_or_else(Edge::new),
                        None => Edge::new(),
                    };
                    a_rebuild = rebuild_edge(
                        &old_edge,
                        &mut anew_edge,
                        ancestor_face.as_ref(),
                        new_face.as_ref(),
                        is_bound,
                    );
                    a_new_shape = Some(anew_edge.0);
                }
                // `cxx:257-266`.
                ShapeType::Vertex => {
                    let Some(v) = Vertex::wrap(a_shape.clone()) else {
                        unreachable!("shape_type() is Vertex")
                    };
                    a_rebuild = true;
                    a_new_shape = Some(rebuild_vertex(&v).0);
                }
                _ => {}
            }
        }
        let mut is_done = a_rebuild;

        // `cxx:270-309`: rebuild the container, removing the locations of the
        // sub-shapes in dependence of `theRemoveLoc`.
        if !is_bound {
            if !a_rebuild {
                // `cxx:277-279`: a topology-only copy; this operation never
                // changes topology, so `Closed` carries over as-is.
                let empty = the_shape.empty_copied();
                empty.set_closed(the_shape.is_closed());
                a_new_shape = Some(empty);
            }
            // `cxx:280-284`.
            let old_loc = the_shape.location().clone();
            if let Some(ns) = a_new_shape.as_mut() {
                if !old_loc.is_identity() {
                    ns.set_location(&TopLocLocation::identity());
                }
            }
            // `cxx:285-287`.
            let orient = the_shape.orientation();
            if let Some(ns) = a_new_shape.as_mut() {
                ns.set_orientation(Orientation::Forward);
            }
            // `cxx:288-295`: `TopoDS_Iterator(aShape, false, isRemoveLoc)`.
            for subshape in iterator_children(&a_shape, false, is_remove_loc) {
                let (is_done_sub, anew_subshape) =
                    self.make_new_shape(&subshape, an_anc.as_ref(), is_remove_loc);
                is_done = is_done || is_done_sub;
                if let Some(ns) = a_new_shape.as_mut() {
                    builder.add(ns, &anew_subshape);
                }
            }
            // `cxx:297-303`.
            if is_done {
                if let Some(ns) = a_new_shape.as_mut() {
                    ns.set_orientation(orient);
                }
            } else {
                a_new_shape = Some(a_shape.clone());
            }
            // `cxx:304-309`: the map keeps the location-free shape, the answer
            // gets the old location back.
            let mut bound = a_new_shape.unwrap_or_else(|| the_shape.clone());
            self.new_shapes.insert(key, bound.clone());
            if !the_remove_loc && !old_loc.is_identity() {
                bound.set_location(&old_loc);
            }
            a_new_shape = Some(bound);
        }

        // `cxx:311-313`.
        let new_shape = a_new_shape.unwrap_or_else(|| the_shape.clone());
        (is_done || is_bound, new_shape)
    }
}

/// `TopoDS_Iterator(S, cumOri, cumLoc)` (`TopoDS_Iterator.cxx:26-82`):
/// children as stored, optionally composed with the parent orientation and
/// moved by the parent location.
fn iterator_children(s: &TopoShape, cum_ori: bool, cum_loc: bool) -> Vec<TopoShape> {
    let stored = s.tshape.read().expect("poisoned TShape lock").children.clone();
    let parent_ori = s.orientation();
    let parent_loc = s.location().clone();
    let do_ori = cum_ori && !parent_ori.is_forward();
    let do_loc = cum_loc && !parent_loc.is_identity();
    if !do_ori && !do_loc {
        return stored;
    }
    stored
        .into_iter()
        .map(|mut child| {
            // `cxx:77`: always composed; `Compose(FORWARD, or)` is `or`.
            child.set_orientation(Orientation::compose(parent_ori, child.orientation()));
            if do_loc {
                // `cxx:79-80`: `Move(myLocation, false)` -> `child * parent`.
                child.move_location(&parent_loc);
            }
            child
        })
        .collect()
}

/// `RebuildShape(Face)` (`cxx:63-77`): bake the face's location into its
/// surface. `None` when the location is identity (nothing to rebuild).
fn rebuild_face(the_face: &Face) -> Option<Face> {
    let a_loc = the_face.0.location().clone();
    if a_loc.is_identity() {
        return None;
    }
    // `cxx:67`: `BRep_Tool::Surface(theFace, aLoc)`.
    let a_surf = BRepTool::face_surface(the_face)?;
    let new_surf: Arc<dyn Surface> = Arc::from(a_surf.transformed(&a_loc.transformation()));
    // `cxx:71-72`: `MakeFace(theNewFace, anewSurf, Tolerance(theFace))`.
    Some(TopoBuilder::new().make_face_with_tolerance(
        new_surf,
        &[],
        BRepTool::face_tolerance(the_face),
    ))
}

/// `RebuildShape(Edge, ..., theFace, theNewFace, isBound)` (`cxx:81-157`).
fn rebuild_edge(
    the_edge: &Edge,
    the_new_edge: &mut Edge,
    the_face: Option<&Face>,
    the_new_face: Option<&Face>,
    is_bound: bool,
) -> bool {
    let geometry = GeometryRegistry::global();
    let mut is_rebuild = false;
    // `cxx:89-113`: the 3D curve is only rebuilt for a first occurrence; a
    // bound edge keeps the one its first occurrence produced.
    if !is_bound {
        let (c3d, a_loc, first3d, last3d) = match geometry.edge_geom(&the_edge.0) {
            // `BRep_Tool::Curve` answers null for a null representation.
            Some(g) if !g.curve3d_null => (
                Some(g.curve),
                the_edge.0.location().clone(),
                g.first,
                g.last,
            ),
            // `cxx:205`: `Curve` on an edge without a 3D representation sets
            // `First = Last = 0` and leaves the location identity.
            _ => (None, TopLocLocation::identity(), 0.0, 0.0),
        };
        if let Some(c3d) = c3d {
            // `cxx:96-102`.
            let curve: Arc<dyn occt_geom::Curve> = if a_loc.is_identity() {
                c3d
            } else {
                Arc::from(c3d.transformed(&a_loc.transformation()))
            };
            // `cxx:105-106`: `UpdateEdge(theNewEdge, C3d, Tolerance(theEdge))`
            // then `Range(theNewEdge, First3d, Last3d)`.
            geometry.set_edge_curve3d(&the_new_edge.0, curve, BRepTool::edge_tolerance(the_edge));
            geometry.set_edge_range(&the_new_edge.0, first3d, last3d);
        }
        // `cxx:107`.
        the_new_edge.0.set_orientation(the_edge.0.orientation());
        // `cxx:108-111`.
        if BRepTool::is_degenerated(the_edge) {
            geometry.set_degenerated(&the_new_edge.0, true);
            geometry.set_curve3d_null(&the_new_edge.0, true);
        }
        is_rebuild = true;
    }

    // `cxx:114-156`.
    if let (Some(face), Some(new_face)) = (the_face, the_new_face) {
        // `cxx:117`: the *old* face's surface decides whether pcurves are
        // carried over; a plane is left without them.
        if let Some(a_surf) = BRepTool::face_surface(face) {
            if a_surf.gp_pln().is_none() {
                // `cxx:121`.
                let (c2d, first2d, last2d) = curve_on_surface(the_edge, face);
                let mut c2d1: Option<Arc<dyn Curve2d>> = None;
                if BRepTool::is_closed_edge_face(the_edge, face) {
                    // `cxx:123-124`.
                    if !BRepTool::is_closed_edge_face(the_new_edge, new_face) {
                        // `cxx:126-127`: the reversed edge of a closed
                        // representation carries the second pcurve.
                        let mut tmp_e = the_edge.clone();
                        tmp_e.0.reverse();
                        let (other, _, _) = curve_on_surface(&tmp_e, face);
                        c2d1 = other;
                        // `cxx:128-133`.
                        let mut or_edge = the_new_edge.0.orientation();
                        if face.0.orientation() == Orientation::Reversed {
                            or_edge = or_edge.reversed();
                        }
                        // `cxx:136-144`: `UpdateEdge(edge, c1, c2, face, 0)`.
                        let mut pair = Vec::with_capacity(2);
                        if or_edge.is_forward() {
                            pair.extend(c2d.clone());
                            pair.extend(c2d1.clone());
                        } else {
                            pair.extend(c2d1.clone());
                            pair.extend(c2d.clone());
                        }
                        match pair.len() {
                            2 => geometry.set_edge_pcurves(
                                &the_new_edge.0,
                                GeometryRegistry::shape_key(&new_face.0),
                                pair,
                            ),
                            1 => geometry.set_edge_pcurve(
                                &the_new_edge.0,
                                GeometryRegistry::shape_key(&new_face.0),
                                pair.into_iter().next().expect("one pcurve"),
                            ),
                            _ => {}
                        }
                    }
                } else if let Some(c2d) = c2d.clone() {
                    // `cxx:148`.
                    geometry.set_edge_pcurve(
                        &the_new_edge.0,
                        GeometryRegistry::shape_key(&new_face.0),
                        c2d,
                    );
                }
                // `cxx:151-154`.
                if c2d.is_some() || c2d1.is_some() {
                    geometry.set_pcurve_range(
                        &the_new_edge.0,
                        GeometryRegistry::shape_key(&new_face.0),
                        first2d,
                        last2d,
                    );
                }
            }
        }
    }
    is_rebuild
}

/// `RebuildShape(Vertex)` (`cxx:162-170`).
fn rebuild_vertex(the_vertex: &Vertex) -> Vertex {
    // `cxx:168`: `BRep_Tool::Pnt` already carries the vertex location.
    let p1 = BRepTool::vertex_point_world(the_vertex);
    let mut vertex = TopoBuilder::new().make_vertex(p1, BRepTool::vertex_tolerance(the_vertex));
    vertex.0.set_orientation(the_vertex.0.orientation());
    vertex
}

/// `BRep_Tool::CurveOnSurface(E, F, First, Last)` (`BRep_Tool.cxx:301-315`,
/// `:327-372`): the pcurve of `edge` on `face` with the representation's
/// parameter range.
///
/// A reversed face reverses the local edge (`cxx:310-313`), and a reversed
/// local edge selects the *second* pcurve of a closed-surface representation
/// (`cxx:356-361`) — the port keeps those two pcurves in
/// forward-then-reversed order, so the selection is an index.
///
/// UNPORTED: the `CurveOnPlane` fallback (`cxx:363-373`, `:377-450`) — when the
/// face carries no CurveOnSurface representation for this surface, OCCT
/// projects the 3D curve through `GeomProjLib::ProjectOnPlane` +
/// `ProjLib_ProjectedCurve`. This port answers "no pcurve" instead, which
/// leaves the COS range unset for that edge.
fn curve_on_surface(edge: &Edge, face: &Face) -> (Option<Arc<dyn Curve2d>>, f64, f64) {
    let geometry = GeometryRegistry::global();
    // `cxx:306-313`: `aLocalEdge.Reverse()` for a reversed face, then
    // `Eisreversed = (E.Orientation() == REVERSED)`.
    let local_reversed =
        edge.0.orientation().is_reversed() != face.0.orientation().is_reversed();
    let face_key = GeometryRegistry::shape_key(&face.0);
    let pcurves = geometry.edge_pcurves(&edge.0, face_key);
    if pcurves.is_empty() {
        return (None, 0.0, 0.0);
    }
    let pc = if pcurves.len() >= 2 && local_reversed {
        pcurves[1].clone()
    } else {
        pcurves[0].clone()
    };
    let (first, last) = geometry
        .pcurve_range(&edge.0, face_key)
        .unwrap_or_else(|| BRepTool::edge_parameters(edge));
    (Some(pc), first, last)
}
