//! `GroupShape` (`TopOpeBRepBuild_FuseFace.cxx:845-971`) — group the shapes of
//! one list into families of shapes that share at least one edge, and write
//! `representative -> family` into `mymapShLSh`.
//!
//! `mylist` is either the reconstructed faces of `myLRF` (`cxx:220`), the wires
//! collected from one face family (`cxx:273`) or the edges collected from one
//! wire family (`cxx:314`, with `Keep_Edge = true`).
//!
//! With `Keep_Edge == false` the shapes that share an edge are merged into one
//! family. With `Keep_Edge == true` no family is merged: the result is
//! `edge -> the shapes carrying that edge`, and — because the map key ignores
//! orientation — that list holds every *occurrence* of the edge, which is what
//! `PerformFace` reads to collect the orientations of a shared edge
//! (`cxx:330-425`).

use crate::abs::{Orientation, ShapeType};
use crate::shape::TopoShape;
use crate::topo_tools::{explorer_shapes, is_same};

use super::shape_list_map::ShapeListMap;

/// `GroupShape(mylist, Keep_Edge, mymapShLSh)` (`cxx:845-971`).
pub(crate) fn group_shape(list: &[TopoShape], keep_edge: bool, out: &mut ShapeListMap) {
    // `locmapEdgLSh` / `locmapShLSh` (`cxx:851-852`).
    let mut map_edg_l_sh = ShapeListMap::new();
    let mut map_sh_l_sh = ShapeListMap::new();

    // `cxx:860-941`: index every shape under its own edges.
    for shap1 in list {
        map_sh_l_sh.bind(shap1, vec![shap1.clone()]);

        for edg1 in explorer_shapes(shap1, ShapeType::Edge) {
            // `cxx:872-879`: an INTERNAL / EXTERNAL edge only counts when the
            // caller asked to keep edges.
            let ori = edg1.orientation();
            let mut edge_ok = true;
            if ori == Orientation::Internal || ori == Orientation::External {
                edge_ok = false;
            }
            if !(edge_ok || keep_edge) {
                continue;
            }

            if !map_edg_l_sh.is_bound(&edg1) {
                map_edg_l_sh.bind(&edg1, vec![shap1.clone()]);
                continue;
            }

            // `cxx:890-891`: append the shape to the edge's list.
            map_edg_l_sh.append(&edg1, shap1);
            if keep_edge {
                continue;
            }

            // `cxx:893-898`: first shape of the list bound to this edge.
            let shap2 = map_edg_l_sh
                .find(&edg1)
                .expect("edge bound above")
                .first()
                .cloned()
                .expect("list appended above");

            // `cxx:900-909`: are the first shape and the current shape
            // already in the same family?
            let already_linked = map_sh_l_sh
                .find(&shap2)
                .map(|l| l.iter().any(|s| is_same(s, shap1)))
                .unwrap_or(false);
            if already_linked {
                continue;
            }

            // `cxx:911-936`: merge the two families by appending each family's
            // member list to every member of the other one. `Lmap1` / `Lmap2`
            // are copies, so the appends do not extend the walked lists.
            let lmap1: Vec<TopoShape> = map_sh_l_sh.find(shap1).unwrap_or(&[]).to_vec();
            let lmap2: Vec<TopoShape> = map_sh_l_sh.find(&shap2).unwrap_or(&[]).to_vec();
            for shap in &lmap1 {
                map_sh_l_sh.change_find(shap).extend(lmap2.iter().cloned());
            }
            for shap in &lmap2 {
                map_sh_l_sh.change_find(shap).extend(lmap1.iter().cloned());
            }
        }
    }

    // `cxx:943-970`: drop the families already reported, or — with
    // `Keep_Edge` — report the edge map itself.
    out.clear();
    if !keep_edge {
        let mut seen: Vec<_> = Vec::new();
        for (shap1, lmap_sh) in map_sh_l_sh.pairs() {
            let key = crate::topo_tools::shape_map_key(shap1);
            if seen.contains(&key) {
                continue;
            }
            seen.push(key);
            out.bind(shap1, lmap_sh.to_vec());
            for shap2 in lmap_sh {
                let k2 = crate::topo_tools::shape_map_key(shap2);
                if !seen.contains(&k2) {
                    seen.push(k2);
                }
            }
        }
    } else {
        for (edg, lmap_sh) in map_edg_l_sh.pairs() {
            out.bind(edg, lmap_sh.to_vec());
        }
    }
}
