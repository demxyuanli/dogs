//! Minimal ShapeBuild_ReShape equivalent (the goal allows a minimal stand-in)
//! plus ShapeFix_ComposeShell::ApplyContext (ShapeFix_ComposeShell.cxx:382-446).

use std::collections::HashMap;
use std::sync::Arc;

use crate::abs::{Orientation, ShapeType};
use crate::shape::{Edge, TopoShape};

use super::wire_segment::WireSegment;

/// ShapeBuild_ReShape substitution interface.
pub trait ReShape {
    /// ShapeBuild_ReShape::Apply: the registered replacement of s, or s itself
    /// when nothing is bound to it.
    fn apply(&self, s: &TopoShape) -> TopoShape;
    /// ShapeBuild_ReShape::Replace: record a substitution. OCCT keys its
    /// TopTools_ShapeMapHasher by IsSame, i.e. by TShape identity.
    fn replace(&mut self, old: &TopoShape, new: &TopoShape);

    /// `ShapeBuild_ReShape::IsRecorded` (`ShapeBuild_ReShape.cxx:...`).
    fn is_recorded(&self, _s: &TopoShape) -> bool {
        false
    }

    /// `ShapeBuild_ReShape::Value` (`ShapeBuild_ReShape.cxx:...`).
    fn value(&self, _s: &TopoShape) -> Option<TopoShape> {
        None
    }

}

/// ReShape with no registered substitutions.
pub struct IdentityReShape;

impl ReShape for IdentityReShape {
    fn apply(&self, s: &TopoShape) -> TopoShape {
        s.clone()
    }
    fn replace(&mut self, _old: &TopoShape, _new: &TopoShape) {}
}

/// Minimal ShapeBuild_ReShape: the direct old->new bindings recorded by
/// Replace. Unlike OCCT this does not rebuild composite shapes from per-sub
/// replacements, which is all ComposeShell needs: it binds whole edges to the
/// wire of their sub-edges and vertices to vertices, then reads those bindings
/// back through Apply.
#[derive(Default, Clone)]
pub struct MapReShape {
    map: HashMap<usize, TopoShape>,
}

impl MapReShape {
    pub fn new() -> Self {
        Self::default()
    }
}

fn key(s: &TopoShape) -> usize {
    Arc::as_ptr(&s.tshape) as usize
}

impl ReShape for MapReShape {
    fn apply(&self, s: &TopoShape) -> TopoShape {
        // `BRepTools_ReShape::Apply` (`BRepTools_ReShape.cxx:...`) goes through
        // `Status` -> `Value`.
        self.value(s).unwrap_or_else(|| s.clone())
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
        self.map.insert(key(&shape), newshape);
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
        let mut res = self.map.get(&key(s))?.clone();
        if s.orientation() == Orientation::Reversed {
            res.reverse();
        }
        if matches!(s.orientation(), Orientation::Internal | Orientation::External) {
            res.set_orientation(s.orientation());
        }
        Some(res)
    }
}

/// ApplyContext (cxx:382-446): substitute edge iedge through context and
/// return how many edges took its place (1 when nothing changed).
pub fn apply_context(wire: &mut WireSegment, iedge: usize, context: &dyn ReShape) -> i32 {
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
