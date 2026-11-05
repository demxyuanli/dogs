//! `TrimEdge` (`BRepFill_Evolved.cxx:2770-2866`).
//!
//! Given the control edges of a bisecting-locus interval and the parameters of
//! the vertices on them, `TrimEdge` sorts the pairs by parameter, drops the
//! vertices that are not proof points of the control edges, merges the
//! duplicates, and rebuilds one edge per (start, end) pair.

use crate::abs::Orientation;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, TopoShape};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::is_same;

use super::side::double_or_not_in_face;

/// `BRepLib::UpdateTolerances(theEdge, false)` (`BRepLib.cxx:1966-1970`) for a
/// bare edge: `InternalUpdateTolerances` (`cxx:1744-1965`) finds no face
/// ancestor, so only the "VERTEX tolerance >= EDGE tolerance" half of the rule
/// can change anything (`cxx:1859-1965`).
fn update_vertex_tolerances(edge: &Edge) {
    let registry = GeometryRegistry::global();
    let tol = BRepTool::edge_tolerance(edge);
    for v in crate::topo_tools::vertices_of(&edge.0) {
        let current = registry
            .vertex_geom(&v.0)
            .map(|g| g.tolerance)
            .unwrap_or(0.0);
        if current < tol {
            registry.set_vertex_tolerance(&v.0, tol);
        }
    }
}

/// `TrimEdge(Edge, TheEdgesControle, TheVer, ThePar, S)`
/// (`BRepFill_Evolved.cxx:2770-2866`).
///
/// `the_ver` / `the_par` are sorted in place, exactly like the reference's
/// `NCollection_Sequence::Exchange` walk.
pub fn trim_edge(
    edge: &Edge,
    the_edges_controle: &[TopoShape],
    the_ver: &mut Vec<TopoShape>,
    the_par: &mut Vec<f64>,
) -> Vec<TopoShape> {
    let mut out: Vec<TopoShape> = Vec::new();
    let builder = TopoBuilder::new();
    let registry = GeometryRegistry::global();

    // `cxx:2783-2793`: the sorting pass of the whole sequence.
    let mut change = true;
    while change {
        change = false;
        for i in 0..the_par.len().saturating_sub(1) {
            if the_par[i] > the_par[i + 1] {
                the_par.swap(i, i + 1);
                the_ver.swap(i, i + 1);
                change = true;
            }
        }
    }

    // `cxx:2799-2811`: a vertex that is not a proof point is removed.
    if !BRepTool::is_degenerated(edge) {
        let mut k = 0;
        while k < the_ver.len() {
            let Some(v) = crate::shape::Vertex::wrap(the_ver[k].clone()) else {
                k += 1;
                continue;
            };
            if double_or_not_in_face(the_edges_controle, &v) {
                the_ver.remove(k);
                the_par.remove(k);
            } else {
                k += 1;
            }
        }
    }

    // `cxx:2818-2836`: merge the doubled vertices of a non-degenerated edge.
    if !BRepTool::is_degenerated(edge) {
        let mut k = 0;
        while k + 1 < the_ver.len() {
            if is_same(&the_ver[k], &the_ver[k + 1]) {
                the_ver.remove(k + 1);
                the_par.remove(k + 1);
                let doubled = match crate::shape::Vertex::wrap(the_ver[k].clone()) {
                    Some(v) => double_or_not_in_face(the_edges_controle, &v),
                    None => false,
                };
                if doubled {
                    the_ver.remove(k);
                    the_par.remove(k);
                }
                // `cxx:2833-2834`: the reference's `k--` puts the next `k++`
                // back on the same index, so the pair is compared again.
                continue;
            }
            k += 1;
        }
    }

    // `cxx:2841-2865`: one edge per uneven/even pair.
    let mut k = 0;
    while k + 1 < the_ver.len() {
        // `Edge.EmptyCopied()` keeps the edge's curve representations and its
        // orientation, which is what `B.Range` re-trims below.
        let mut new_edge = edge.0.empty_copied();
        if let Some(geom) = registry.edge_geom(edge) {
            registry.set_edge(&new_edge, geom);
        }
        registry.copy_edge_pcurve_slots(&new_edge, &edge.0);

        if new_edge.orientation() == Orientation::Reversed {
            builder.add(&mut new_edge, &the_ver[k].oriented(Orientation::Reversed));
            builder.add(&mut new_edge, &the_ver[k + 1].oriented(Orientation::Forward));
        } else {
            builder.add(&mut new_edge, &the_ver[k].oriented(Orientation::Forward));
            builder.add(&mut new_edge, &the_ver[k + 1].oriented(Orientation::Reversed));
        }
        registry.set_ranges_all(&new_edge, the_par[k], the_par[k + 1]);

        // `cxx:2860`: `BRepLib::UpdateTolerances(NewEdge, false)`.
        update_vertex_tolerances(&Edge(new_edge.clone()));

        out.push(new_edge);
        k += 2;
    }

    out
}
