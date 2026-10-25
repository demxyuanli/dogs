//! `GroupEdge` (`TopOpeBRepBuild_FuseFace.cxx:975-1074`) — chain the edges
//! that meet at a vertex the way `PerformEdge` collected them.
//!
//! Source table: the vertex map (`myVerLEdg`) is `vertex -> the two edges that
//! meet there` (`cxx:640-661`). This function walks those pairs and builds
//! `edge -> the chain of edges reachable from it`, with the OCCT rule that the
//! two edges of a vertex are put in *one* family when neither of them is bound
//! yet, and that the new edge absorbs the old edge's chain when exactly one of
//! the two already belongs to a chain.

use super::shape_list_map::ShapeListMap;

/// `GroupEdge(mymapVerLEdg, mymapEdgLEdg)` (`cxx:975-1074`).
pub fn group_edge(map_ver_l_edg: &ShapeListMap, out: &mut ShapeListMap) {
    // `mapEdgLEdg` — the chain map built by this call (`cxx:984`).
    let mut map_edg_l_edg = ShapeListMap::new();

    for (_, lmap_edg) in map_ver_l_edg.pairs() {
        // `cxx:996-997`: the vertex holds exactly two edges (filtered at
        // `cxx:653-657`).
        let (Some(edg1), Some(edg2)) = (lmap_edg.first().cloned(), lmap_edg.get(1).cloned())
        else {
            continue;
        };

        // `cxx:1000-1022`.
        let edge1_add = !map_edg_l_edg.is_bound(&edg1);
        let edge2_add = !map_edg_l_edg.is_bound(&edg2);
        let edgold;
        let edgnew;
        if !edge1_add {
            edgold = edg1.clone();
            edgnew = edg2.clone();
        } else {
            edgold = edg2.clone();
            edgnew = edg1.clone();
        }

        if !(edge1_add || edge2_add) {
            // `cxx:1024-1026`: both edges already have a chain.
            continue;
        }
        if edge1_add && edge2_add {
            // `cxx:1027-1032`: neither is chained yet — one chain for both.
            map_edg_l_edg.bind(&edg1, lmap_edg.to_vec());
            map_edg_l_edg.bind(&edg2, lmap_edg.to_vec());
            continue;
        }

        // `cxx:1033-1051`: the bound edge is `edgold`, the free one `edgnew`;
        // the free edge inherits the whole chain.
        map_edg_l_edg.bind(&edgnew, vec![edgnew.clone()]);
        let lmap_edg1 = map_edg_l_edg
            .find(&edgold)
            .expect("edgold is bound")
            .to_vec();
        for edg22 in &lmap_edg1 {
            map_edg_l_edg.append(&edgnew, edg22);
            map_edg_l_edg.append(edg22, &edgnew);
        }
    }

    // `cxx:1053-1074`: keep one chain per set of already reported edges.
    let mut seen: Vec<_> = Vec::new();
    for (edg1, lmap_edg) in map_edg_l_edg.pairs() {
        let key = crate::topo_tools::shape_map_key(edg1);
        if seen.contains(&key) {
            continue;
        }
        seen.push(key);
        out.bind(edg1, lmap_edg.to_vec());
        for edg2 in lmap_edg {
            let k2 = crate::topo_tools::shape_map_key(edg2);
            if !seen.contains(&k2) {
                seen.push(k2);
            }
        }
    }
}