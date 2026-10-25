//! `BuildListConnexEdge` (`TopOpeBRepTool_FuseEdges.cxx:480-537`) and
//! `NextConnexEdge` (`cxx:549-618`).
//!
//! The chain of fusible edges around one vertex is walked in both directions
//! from the seed edge: forward from its last vertex (`cxx:495-513`) and
//! backward from its first vertex with `NCollection_List::Prepend`
//! (`cxx:519-536`), which is why the result list is ordered along the chain.
//!
//! `NextConnexEdge` is the three-condition gate of the class:
//!
//! 1. `cxx:559` — the vertex has exactly two connex edges;
//! 2. `cxx:579-611` — both edges have the same `EDGE -> FACE` ancestor list and
//!    that list is shorter than three;
//! 3. `cxx:603-606` — [`same_support`] holds.
//!
//! `myAvoidEdg` is tested before the other two (`cxx:569-572`).
//!
//! The `while` loops also stop as soon as the candidate is already in
//! `mapUniqEdg` (`cxx:498-501`, `cxx:521-524`); the walk deliberately keeps the
//! candidate that closes a loop out of the list, so `NextConnexEdge`'s own
//! `IsSame` probe (`cxx:562-567`) is what refuses to step back onto the seed.

use crate::abs::Orientation;
use crate::shape::TopoShape;
use crate::topo_tools::{is_same, shape_map_key, ShapeMapKey};

use super::same_support::same_support;
use super::{cumulated_vertex, FuseEdges};

impl FuseEdges {
    /// `BuildListConnexEdge(theEdge, theMapUniq, theLstEdg)` (`cxx:480-537`).
    ///
    /// `theLstEdg` is cleared first (`cxx:492`); `theMapUniq` is the caller's
    /// set of already-visited edges, appended to as the walk proceeds.
    pub(super) fn build_list_connex_edge(
        &self,
        the_edge: &TopoShape,
        the_map_uniq: &mut Vec<ShapeMapKey>,
        the_lst_edg: &mut Vec<TopoShape>,
    ) {
        // `cxx:487`: `VL = TopExp::LastVertex(TopoDS::Edge(theEdge), true)`.
        let mut vl = cumulated_vertex(the_edge, true);

        the_lst_edg.clear();
        the_lst_edg.push(the_edge.clone());
        the_map_uniq.push(shape_map_key(the_edge));

        // `cxx:495-513`: walk forward from the last vertex.
        let mut edgecur = the_edge.clone();
        while let Some(v) = vl.as_ref() {
            let Some(edgeconnex) = self.next_connex_edge(v, &edgecur) else {
                break;
            };
            // `cxx:498-501`
            if the_map_uniq.contains(&shape_map_key(&edgeconnex)) {
                break;
            }
            the_lst_edg.push(edgeconnex.clone());
            edgecur = edgeconnex;
            // `cxx:505-510`: an INTERNAL / EXTERNAL edge is never walked
            // through.
            let ori2 = edgecur.orientation();
            if ori2 == Orientation::External || ori2 == Orientation::Internal {
                break;
            }
            vl = cumulated_vertex(&edgecur, true);
            the_map_uniq.push(shape_map_key(&edgecur));
        }

        // `cxx:515-536`: walk backward from the first vertex, prepending.
        edgecur = the_edge.clone();
        let mut vf = cumulated_vertex(the_edge, false);
        while let Some(v) = vf.as_ref() {
            let Some(edgeconnex) = self.next_connex_edge(v, &edgecur) else {
                break;
            };
            // `cxx:521-524`
            if the_map_uniq.contains(&shape_map_key(&edgeconnex)) {
                break;
            }
            the_lst_edg.insert(0, edgeconnex.clone());
            edgecur = edgeconnex;
            // `cxx:528-533`
            let ori2 = edgecur.orientation();
            if ori2 == Orientation::External || ori2 == Orientation::Internal {
                break;
            }
            vf = cumulated_vertex(&edgecur, false);
            the_map_uniq.push(shape_map_key(&edgecur));
        }
    }

    /// `NextConnexEdge(theVertex, theEdge, theEdgeConnex)` (`cxx:549-618`).
    ///
    /// `None` is the C++ `false`.
    pub(super) fn next_connex_edge(
        &self,
        the_vertex: &TopoShape,
        the_edge: &TopoShape,
    ) -> Option<TopoShape> {
        // `cxx:555`: `myMapVerLstEdg.FindFromKey(theVertex)`.
        let lmap_edg = self.my_map_ver_lst_edg.find_from_key(the_vertex)?;

        // 1st condition (`cxx:559`)
        if lmap_edg.len() != 2 {
            return None;
        }
        let mut it = lmap_edg.iter();
        let mut cand = it.next()?.clone();
        if is_same(the_edge, &cand) {
            cand = it.next()?.clone();
        }

        // `cxx:569-572`
        if self.is_avoided(&cand) {
            return None;
        }

        // 2nd condition (`cxx:576-600`)
        // `FindFromKey` raises `Standard_NoSuchObject` in OCCT when the edge is
        // under no face; the port reads that as an empty ancestor list, which
        // fails the extent comparison unless the other side is empty too.
        let lmap_fac1 = self
            .my_map_edg_lst_fac
            .find_from_key(the_edge)
            .unwrap_or(&[]);
        let lmap_fac2 = self.my_map_edg_lst_fac.find_from_key(&cand).unwrap_or(&[]);
        if lmap_fac1.len() != lmap_fac2.len() || lmap_fac1.len() >= 3 {
            return None;
        }
        // `cxx:585-599`: every face of the first list must also be in the
        // second (`face1.IsSame(face2)`).
        for face1 in lmap_fac1 {
            if !lmap_fac2.iter().any(|face2| is_same(face1, face2)) {
                return None;
            }
        }

        // 3rd condition (`cxx:603-606`)
        let (e1, e2) = (
            crate::shape::Edge::wrap(the_edge.clone())?,
            crate::shape::Edge::wrap(cand.clone())?,
        );
        if !same_support(&e1, &e2) {
            return None;
        }
        Some(cand)
    }
}
