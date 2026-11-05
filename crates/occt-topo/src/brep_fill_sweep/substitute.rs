//! `BRepFill_Sweep::Substitute` (`BRepFill_Sweep.cxx:1341-1367`): register the
//! replacement of an edge and of both its vertices in a
//! [`BRepTools_Substitution`].
//!
//! Translation note: `TopoDS_Vertex` / `TopoDS_Edge` can be null in the
//! reference (`TopExp::Vertices` leaves them null), and a null shape is not
//! representable as a [`TopoShape`]; the port skips those entries instead of
//! inserting a null key. Every edge built by `BuildEdge` has both vertices, so
//! the arms are unreachable from the `BRepFill_Sweep` call sites.

use crate::abs::Orientation;
use crate::brep_tools_substitution::BRepToolsSubstitution;
use crate::shape::Edge;
use crate::topo_tools::edge_vertices;

/// `Substitute(aSubstitute, Old, New)` (`cxx:1341-1367`).
pub fn substitute(a_substitute: &mut BRepToolsSubstitution, old: &Edge, new: &Edge) {
    let (old_v1, old_v2) = edge_vertices(old);
    let (new_v1, new_v2) = edge_vertices(new);

    // `cxx:1352-1363`.
    if let (Some(ov), Some(nv)) = (&old_v1, &new_v1) {
        if !a_substitute.is_copied(&ov.0) {
            a_substitute.substitute(&ov.0, &[nv.0.oriented(Orientation::Forward)]);
        }
    }
    if let (Some(ov), Some(nv)) = (&old_v2, &new_v2) {
        if !a_substitute.is_copied(&ov.0) {
            a_substitute.substitute(&ov.0, &[nv.0.oriented(Orientation::Forward)]);
        }
    }

    // `cxx:1364-1367`.
    if !a_substitute.is_copied(&old.0) {
        a_substitute.substitute(&old.0, &[new.0.oriented(Orientation::Forward)]);
    }
}
