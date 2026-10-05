//! ShapeBuild_ReShape equivalence: `Apply` / `applyImpl`
//! (`ShapeBuild_ReShape.cxx:191-324`) plus ShapeFix_ComposeShell::ApplyContext
//! (`ShapeFix_ComposeShell.cxx:382-446`).

use std::collections::HashMap;
use std::sync::Arc;

use crate::abs::{Orientation, ShapeType};
use crate::builder::TopoBuilder;
use crate::shape::{Edge, TopoShape};

use super::wire_segment::WireSegment;

/// ShapeBuild_ReShape substitution interface.
pub trait ReShape {
    /// `ShapeBuild_ReShape::Apply(shape)` (`ShapeBuild_ReShape.cxx:191-196` ->
    /// `applyImpl` `:200-324` with `until = TopAbs_SHAPE`): the recorded
    /// replacement of `s` when there is one, otherwise `s` rebuilt from its
    /// replaced sub-shapes.
    fn apply(&mut self, s: &TopoShape) -> TopoShape;
    /// ShapeBuild_ReShape::Replace: record a substitution. OCCT keys its
    /// TopTools_ShapeMapHasher by IsSame, i.e. by TShape identity.
    fn replace(&mut self, old: &TopoShape, new: &TopoShape);
    /// `ShapeBuild_ReShape::IsRecorded`.
    fn is_recorded(&self, _s: &TopoShape) -> bool {
        false
    }
    /// `ShapeBuild_ReShape::Value`.
    fn value(&self, _s: &TopoShape) -> Option<TopoShape> {
        None
    }
}

/// ReShape with no registered substitutions.
pub struct IdentityReShape;

impl ReShape for IdentityReShape {
    fn apply(&mut self, s: &TopoShape) -> TopoShape {
        s.clone()
    }
    fn replace(&mut self, _old: &TopoShape, _new: &TopoShape) {}
}

/// Minimal ShapeBuild_ReShape: the direct old->new bindings recorded by
/// Replace, plus `applyImpl`'s rebuild (`ShapeBuild_ReShape.cxx:200-324`) of a
/// shape whose direct sub-shapes were replaced.
///
/// The rebuild is what ComposeShell relies on: `SplitWire` empty-copies the
/// endpoint of the first edge it splits and records `Replace(prevV, fV)`
/// (`ShapeFix_ComposeShell.cxx:1262-1296`), so *every other* edge that used
/// the original vertex has to come back rebuilt onto the copy - that is what
/// keeps the split wire closed and lets `CollectWires` join the pieces
/// (`cxx:2119-2126` applies the context to every wire of the sequence right
/// after `SplitByLine`).
#[derive(Default, Clone)]
pub struct MapReShape {
    /// `BRepTools_ReShape::myMap` (`BRepTools_ReShape.hxx:239-240`) is a
    /// `NCollection_DataMap<TopoDS_Shape, TReplacement, TopTools_ShapeMapHasher>`:
    /// it stores the **key shape** and the replacement. The key shape must be
    /// kept alive, otherwise the `TopTools_ShapeMapHasher` key (the TShape
    /// address, `TopTools_ShapeMapHasher.cxx:22-26`) becomes dangling and a
    /// later shape allocated in the freed slot would be read as the replaced
    /// one (a stale edge key matching a Face defeated the edge count of
    /// `ShapeFix_ComposeShell::ApplyContext`, `cxx:415-442`). Storing the pair
    /// keeps the same lifetime as OCCT.
    map: HashMap<usize, (TopoShape, TopoShape)>,
}

impl MapReShape {
    pub fn new() -> Self {
        Self::default()
    }

    /// `ShapeBuild_ReShape::applyImpl` (`ShapeBuild_ReShape.cxx:200-324`) with
    /// `until = TopAbs_SHAPE`. `in_flight` is OCCT's DFS cycle guard
    /// (`cxx:220-226`): a shape already being descended is returned through its
    /// direct replacement only.
    fn apply_impl(&mut self, s: &TopoShape, in_flight: &mut Vec<usize>) -> TopoShape {
        if s.is_null() {
            return s.clone(); // `cxx:205-208`.
        }
        // `cxx:210-226`: direct replacement first, `NULL` when removed.
        let new_shape = match self.value(s) {
            Some(v) => v,
            None => s.clone(),
        };
        if in_flight.contains(&key(s)) {
            return new_shape;
        }
        if !crate::topo_tools_full::is_same(&new_shape, s) {
            // `cxx:228-237`: the modifications apply to the replacement as well.
            in_flight.push(key(s));
            let res = self.apply_impl(&new_shape, in_flight);
            in_flight.pop();
            return res;
        }
        let st = s.shape_type();
        if matches!(st, ShapeType::Vertex | ShapeType::Shape) {
            return s.clone(); // `cxx:244-247`.
        }
        let children: Vec<TopoShape> = {
            let ts = s.tshape.read().expect("poisoned TShape lock");
            ts.children.clone()
        };
        if children.is_empty() {
            return s.clone(); // `cxx:302-305`: nothing was modified.
        }
        let mut modified = false;
        let mut new_children: Vec<TopoShape> = Vec::with_capacity(children.len());
        for c in &children {
            let nc = self.apply_impl(c, in_flight);
            if !crate::topo_tools_full::is_same(&nc, c) {
                modified = true;
            }
            new_children.push(nc);
        }
        if !modified {
            return s.clone(); // `cxx:302-305`.
        }
        // `cxx:249-254`: `EmptyCopied()` (a new empty shape of the same
        // type) put FORWARD while the components are added, so
        // `TopoDS_Builder::Add` (`TopoDS_Builder.cxx:74-91`) stores each one
        // with the very orientation it has in `theShape`; the source
        // orientation is restored at `cxx:317`.
        let builder = TopoBuilder::new();
        let mut copy = TopoShape::new(st);
        copy.set_orientation(Orientation::Forward);
        copy.set_location(s.location());
        if st == ShapeType::Edge {
            // `BRep_TEdge::EmptyCopy` (`BRep_TEdge.cxx:104-129`) keeps the
            // tolerance, the `SameParameter` / `SameRange` / `Degenerated`
            // flags and a copy of every curve representation, and the edge core
            // carries the 3D range (`cxx:308-311` restores it through
            // `ShapeBuild_Edge::CopyRanges`).
            copy_edge_geom(&copy, s);
        } else if st == ShapeType::Face {
            // `BRep_TFace::EmptyCopy` (`BRep_TFace.cxx:37-43`) keeps the
            // surface, the location and the tolerance. Without this a rebuilt
            // face lost its `Geom_Surface` and every later query (bounding box,
            // pcurve projection, meshing) saw a null-surface face.
            copy_face_geom(&copy, s);
        }
        for (orig, nc) in children.iter().zip(new_children.iter()) {
            if nc.is_null() {
                continue; // `cxx:271-276`: removed component (DONE4).
            }
            if st == ShapeType::Compound || nc.shape_type() == orig.shape_type() {
                builder.add(&mut copy, nc); // `cxx:278-281`.
                continue;
            }
            // `cxx:282-299`: the replacement holds the components of the
            // replaced shape's type; anything else is dropped (FAIL1). The
            // components are taken through `TopoDS_Iterator aSubIt(aNewShape)`
            // (`ShapeBuild_ReShape.cxx:284`), i.e. with `cumOri = true`: each
            // component is oriented by `aNewShape`'s orientation first.
            let subs = {
                let ts = nc.tshape.read().expect("poisoned TShape lock");
                ts.children.clone()
            };
            let nc_ori = nc.orientation();
            for sub in subs {
                if sub.shape_type() == orig.shape_type() {
                    let mut c = sub;
                    c.set_orientation(Orientation::compose(nc_ori, c.orientation()));
                    builder.add(&mut copy, &c);
                }
            }
        }
        // `cxx:306-315`: the 3D range of an edge is restored by `CopyRanges`
        // (the edge core carries it, see `copy_edge_geom`), while a WIRE / SHELL
        // carries a `Closed` flag that `EmptyCopied` drops and OCCT recomputes
        // with `BRep_Tool::IsClosed(aResult)`.
        if matches!(st, ShapeType::Wire | ShapeType::Shell) {
            copy.set_closed(brep_tool_is_closed(&copy));
        }
        copy.set_orientation(s.orientation()); // `cxx:317`.
        let result = copy;
        // `cxx:320`: the rebuilt shape is recorded, so every later `Apply` of
        // the same source shape returns this very shape.
        self.replace(s, &result);
        result
    }
}

/// `TopoDS_Shape::EmptyCopied` (`TopoDS_Shape.hxx:294-302`) plus
/// `BRep_TEdge::EmptyCopy` (`BRep_TEdge.cxx:104-129`) and `CopyRanges`
/// (`ShapeBuild_Edge.cxx:125`) for an edge: the tolerance, the
/// `SameParameter` / `SameRange` / `Degenerated` flags, the 3D range and a copy
/// of every curve representation (3D curve and pcurves), with no components.
fn copy_edge_geom(dst: &TopoShape, src: &TopoShape) {
    let (core, pcs) = {
        let ts = src.tshape.read().expect("poisoned TShape lock");
        (ts.edge_core().cloned(), ts.edge_pcurves().cloned())
    };
    let mut ts = dst.tshape.write().expect("poisoned TShape lock");
    if let Some(c) = core {
        let d = ts.edge_core_mut();
        d.curve = c.curve;
        d.first = c.first;
        d.last = c.last;
        d.tolerance = c.tolerance;
        d.same_parameter = c.same_parameter;
        d.same_range = c.same_range;
        d.degenerated = c.degenerated;
    }
    if let Some(p) = pcs {
        *ts.edge_pcurves_mut() = p;
    }
}

fn key(s: &TopoShape) -> usize {
    Arc::as_ptr(&s.tshape) as usize
}

/// `BRep_TFace::EmptyCopy` (`BRep_TFace.cxx:37-43`) for a face: the surface and
/// the tolerance are carried over (the location is restored by the caller).
/// `myNaturalRestriction` deliberately stays at its constructor default
/// (`false`, `BRep_TFace.cxx:28-31`): `EmptyCopy` does not copy it either.
fn copy_face_geom(dst: &TopoShape, src: &TopoShape) {
    let core = {
        let ts = src.tshape.read().expect("poisoned TShape lock");
        ts.face_core().cloned()
    };
    if let Some(c) = &core {
        let mut ts = dst.tshape.write().expect("poisoned TShape lock");
        let d = ts.face_core_mut();
        d.surface = c.surface.clone();
        d.tolerance = c.tolerance;
    }
    // `BRep_Tool::CurveOnSurface` resolves the edge's representation through
    // the face's surface, so the rebuilt face must be registered against that
    // very surface; otherwise `repr_key` falls back to the face pointer and
    // every pcurve (keyed by the surface) becomes invisible on the copy.
    if let Some(surf) = core.as_ref().and_then(|c| c.surface.clone()) {
        crate::tgeometry::GeometryRegistry::global().register_face_surface(dst, &surf);
    }
}

/// `BRep_Tool::IsClosed(theShape)` (`BRep_Tool.cxx:1707-1755`) for a WIRE or a
/// SHELL, which is the `cxx:313` recomputation for a rebuilt composite.
///
/// OCCT walks the FORWARD-oriented shape with a `TopExp_Explorer`, skips the
/// components whose explored orientation is INTERNAL or EXTERNAL (plus
/// degenerated edges for a shell), and toggles each remaining boundary
/// component in a map keyed by `IsSame` (TShape identity, `TopTools_ShapeMapHasher`):
/// `if (!aMap.Add(x)) aMap.Remove(x);`. The shape is closed when at least one
/// boundary component exists and the map ends up empty, i.e. when every one of
/// them appears an even number of times.
///
/// This is *not* the stored `Closed` flag (dropped by `EmptyCopied`) and *not*
/// a positional chain test.
fn brep_tool_is_closed(s: &TopoShape) -> bool {
    let mut has_bound = false;
    let mut odd: Vec<usize> = Vec::new();
    match s.shape_type() {
        ShapeType::Wire => {
            for v in crate::topo_tools_full::vertices_of(s) {
                if matches!(v.0.orientation(), Orientation::Internal | Orientation::External) {
                    continue;
                }
                has_bound = true;
                toggle(&mut odd, key(&v.0));
            }
        }
        ShapeType::Shell => {
            for e in crate::topo_tools_full::edges_of(s) {
                if crate::brep_tool::BRepTool::is_degenerated(&e)
                    || matches!(e.0.orientation(), Orientation::Internal | Orientation::External)
                {
                    continue;
                }
                has_bound = true;
                toggle(&mut odd, key(&e.0));
            }
        }
        // `BRep_Tool.cxx:1755`: `return theShape.Closed();`
        _ => return s.closed(),
    }
    has_bound && odd.is_empty()
}

/// `NCollection_Map::Add` + `Remove` pair (`BRep_Tool.cxx:1722-1724`): keep the
/// members with an odd number of occurrences.
fn toggle(set: &mut Vec<usize>, k: usize) {
    match set.iter().position(|x| *x == k) {
        Some(i) => {
            set.remove(i);
        }
        None => set.push(k),
    }
}

impl ReShape for MapReShape {
    fn apply(&mut self, s: &TopoShape) -> TopoShape {
        self.apply_impl(s, &mut Vec::new())
    }
    /// `BRepTools_ReShape::replace` (`BRepTools_ReShape.cxx:164-209`): the
    /// replacement is stored **orientation-normalised** — when the key shape is
    /// REVERSED, both it and the replacement are reversed before the bind; for
    /// INTERNAL/EXTERNAL the replacement keeps the relative orientation and the
    /// key becomes FORWARD. `Value` re-applies the query's orientation.
    fn replace(&mut self, old: &TopoShape, new: &TopoShape) {
        let mut shape = old.clone();
        let mut newshape = new.clone();
        if shape.orientation() == Orientation::Reversed {
            shape.reverse();
            newshape.reverse();
        } else if matches!(shape.orientation(), Orientation::Internal | Orientation::External) {
            newshape.set_orientation(if newshape.orientation() == shape.orientation() {
                Orientation::Forward
            } else {
                Orientation::Reversed
            });
            shape.set_orientation(Orientation::Forward);
        }
        self.map.insert(key(&shape), (shape, newshape));
    }
    fn is_recorded(&self, s: &TopoShape) -> bool {
        self.map.contains_key(&key(s))
    }
    /// `BRepTools_ReShape::Value` (`BRepTools_ReShape.cxx:230-279`): the stored
    /// replacement, REVERSED again when the queried shape is REVERSED.
    ///
    /// Without this a seam edge that appears twice in a wire (Forward and
    /// Reversed) got the *same* stored orientation for both occurrences, which
    /// collapsed the periodic face's boundary (a3n00 f173/f176).
    fn value(&self, s: &TopoShape) -> Option<TopoShape> {
        let mut res = self.map.get(&key(s))?.1.clone();
        if s.orientation() == Orientation::Reversed {
            res.reverse();
        }
        if matches!(s.orientation(), Orientation::Internal | Orientation::External) {
            res.set_orientation(s.orientation());
        }
        Some(res)
    }
}
/// `Handle(ShapeBuild_ReShape)`: the *same* instance reaches every fix tool of
/// one `ShapeFix_Shape` pass.
///
/// `ShapeFix_Shape.cxx:75` creates the context once, and
/// `ShapeFix_Solid.cxx:468` / `ShapeFix_Shell.cxx:108` /
/// `ShapeFix_Face.cxx:379` hand that same handle down to every face. So a
/// substitution recorded by one face's `ShapeFix_ComposeShell`
/// (`ShapeFix_ComposeShell.cxx:3456-3457`, `:1397`, `:1412`) is visible to
/// every later `Context()->Apply` of the pass - in particular to the final
/// `myResult = Context()->Apply(S)` (`ShapeFix_Shape.cxx:257`,
/// `ShapeFix_Shell.cxx:139`), which is what rewrites a face whose *neighbour*
/// was split: T0M face 1693 (torus, R=12.5/r=3) gets its bottom R=9.5 circle
/// cut at U=4.97868 by the seam the neighbouring cylinder
/// (face 1695, R=9.5) inserts at its own existing vertex azimuth, and has 5
/// edges only because of that shared pass.
#[derive(Clone, Default)]
pub struct SharedReShape(std::rc::Rc<std::cell::RefCell<MapReShape>>);

impl SharedReShape {
    pub fn new() -> Self {
        Self::default()
    }
}

impl ReShape for SharedReShape {
    fn apply(&mut self, s: &TopoShape) -> TopoShape {
        self.0.borrow_mut().apply(s)
    }
    fn replace(&mut self, old: &TopoShape, new: &TopoShape) {
        self.0.borrow_mut().replace(old, new);
    }
    fn is_recorded(&self, s: &TopoShape) -> bool {
        self.0.borrow().is_recorded(s)
    }
    fn value(&self, s: &TopoShape) -> Option<TopoShape> {
        self.0.borrow().value(s)
    }
}

/// ApplyContext (cxx:382-446): substitute edge iedge through context and
/// return how many edges took its place (1 when nothing changed).
pub fn apply_context(wire: &mut WireSegment, iedge: usize, context: &mut dyn ReShape) -> i32 {
    let edge = match wire.edge(iedge) {
        Some(e) => e.clone(),
        None => return 1,
    };
    let res = context.apply(&edge.0);
    if crate::topo_tools_full::is_same(&res, &edge.0) {
        return 1;
    }
    if res.shape_type() == ShapeType::Edge {
        wire.set_edge(iedge, Edge(res));
        return 1;
    }
    let mut index = iedge as i32;
    let segw: Vec<Edge> = {
        let ts = res.tshape.read().expect("poisoned TShape lock");
        ts.children
            .iter()
            .filter(|c| c.is_edge())
            .map(|c| Edge(c.clone()))
            .collect()
    };
    if !segw.is_empty() {
        let (iumin, iumax, ivmin, ivmax) =
            wire.get_patch_index(iedge).unwrap_or((0, 0, 0, 0));
        let nb_edges = segw.len() as i32;
        let eo = edge.0.orientation();
        for i in 1..=nb_edges {
            let ind = if eo == Orientation::Forward || eo == Orientation::Internal {
                i
            } else {
                nb_edges - i + 1
            };
            let ae = segw[(ind - 1) as usize].clone();
            if i == 1 {
                wire.set_edge(index as usize, ae);
            } else {
                wire.add_edge_patch(index as usize, ae, iumin, iumax, ivmin, ivmax);
            }
            index += 1;
        }
    }
    index - iedge as i32
}
