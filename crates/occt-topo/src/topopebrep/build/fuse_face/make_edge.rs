//! `MakeEdge` (`TopOpeBRepBuild_FuseFace.cxx:1078-1177`) — rebuild each chain
//! of same-support edges of `GroupEdge` as one edge spanning the two extreme
//! vertices.
//!
//! `PerformEdge` calls it as `MakeEdge(mapEdgLEdg)` (`cxx:672`) on a **local**
//! map, so the rebuilt chains never leave the function: OCCT's last statement
//! `mymapEdgLEdg = mapEdgLEdg` (`cxx:1176`) writes into that local copy. What
//! this function does change is the control flow of `PerformEdge` — the early
//! `return` at `cxx:1147` abandons the remaining chains — and the edges it
//! creates in the shared shape registry, exactly as the C++ creates
//! `TopoDS_Edge`s that nothing references afterwards.
//!
//! Per chain (`cxx:1095-1176`):
//!
//! 1. `mapVerInt` (`cxx:1100-1119`): a vertex seen once is an end point
//!    (`VertexExtrem = 1`), a vertex seen twice is interior (`0`). `Bind`
//!    overwrites the value and keeps the insertion position, so the second
//!    sighting flips `1` to `0`.
//! 2. `myEdgeLV` / `myEdgeLMV` (`cxx:1121-1136`): ends and interior vertices,
//!    in `mapVerInt` order. `myEdgeLMV` is filled with `INTERNAL`-oriented
//!    copies and never read — its only consumer, `cxx:1166-1172`, is commented
//!    out in OCCT.
//! 3. `number != 2` returns from the **whole** function (`cxx:1138-1148`).
//! 4. `verf` / `verl` are the two ends oriented `FORWARD` (`cxx:1149-1157`),
//!    `curv = BRep_Tool::Curve(edg1, loc, first, last)` (`cxx:1158-1161`) — the
//!    curve **without** `loc` applied, while the vertices are in world
//!    coordinates; `BRepLib_MakeEdge ME(curv, verf, verl)` (`cxx:1163`) is
//!    `Init(C, V1, V2)` (`BRepLib_MakeEdge.cxx:542-572`), i.e.
//!    [`init_3d_projecting`], and `edgnew = ME.Edge()` is
//!    `BRepLib_MakeEdge.cxx:780-795`: `MakeEdge(E, C, preci)`, `Add(E, V1)`,
//!    `Add(E, V2)`, `Range`, `Degenerated`, `Closed`.
//! 5. `mapEdgLEdg.Bind(edgnew, LmapEdg)` (`cxx:1174`).
//!
//! `BRep_Tool::Curve` answers a null handle only for an edge without a
//! `BRep_Curve3D`; `Project` would then dereference it (`cxx:60-62`), so the
//! port skips that chain instead of creating a degenerate edge.

use crate::abs::{Orientation, ShapeType};
use crate::brep_lib_make_edge::{add_endpoint_vertices, init_3d_projecting};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Vertex};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{explorer_shapes, shape_map_key, ShapeMapKey};

use super::shape_list_map::ShapeListMap;

/// `MakeEdge(mymapEdgLEdg)` (`cxx:1078-1177`).
pub fn make_edge(map_edg_l_edg: &ShapeListMap) {
    for (edg1, lmap_edg) in map_edg_l_edg.pairs() {
        // `mapVerInt` (`cxx:1100`): vertex -> 1 for an extreme vertex (seen
        // once), 0 for an interior one (seen twice). `Bind` overwrites, so a
        // repeated vertex flips 1 to 0 (`cxx:1115-1120`).
        let mut map_ver_int: Vec<(ShapeMapKey, crate::shape::TopoShape, bool)> = Vec::new();
        for edg2 in lmap_edg {
            let Some(edge) = Edge::wrap(edg2.clone()) else {
                continue;
            };
            for ver in explorer_shapes(&edge.0, ShapeType::Vertex) {
                let key = shape_map_key(&ver);
                if let Some(e) = map_ver_int.iter_mut().find(|(k, _, _)| *k == key) {
                    e.2 = false;
                } else {
                    map_ver_int.push((key, ver, true));
                }
            }
        }

        // `cxx:1122-1136`: the extreme vertices go to `myEdgeLV`; the interior
        // ones to `myEdgeLMV`, which OCCT builds and then never reads.
        let mut my_edge_lv: Vec<crate::shape::TopoShape> = Vec::new();
        for (_, ver, extrem) in map_ver_int.iter() {
            if *extrem {
                my_edge_lv.push(ver.clone());
            } else {
                // `cxx:1134`: `ver.Oriented(TopAbs_INTERNAL)`.
                let _ = ver.oriented(Orientation::Internal);
            }
        }

        // `cxx:1138-1147`: a chain that does not end on exactly two extreme
        // vertices aborts the whole reconstruction.
        if my_edge_lv.len() != 2 {
            return;
        }

        // `cxx:1149-1157`: both ends oriented FORWARD.
        let verf = Vertex(my_edge_lv[0].oriented(Orientation::Forward));
        let verl = Vertex(my_edge_lv[1].oriented(Orientation::Forward));

        // `cxx:1158-1161`: `BRep_Tool::Curve(edg1, loc, first, last)` without
        // `loc`; `edg1` is the map key, the chain's representative edge.
        let Some(old) = Edge::wrap(edg1.clone()) else {
            continue;
        };
        let Some(curv) = BRepTool::edge_curve(&old) else {
            continue;
        };

        // `cxx:1163-1164`: `BRepLib_MakeEdge ME(curv, verf, verl)`;
        // `edgnew = ME.Edge()`. `!IsDone()` leaves `myShape` null in OCCT, and a
        // null edge bound into a map that is then discarded changes nothing.
        let Some(build) = init_3d_projecting(&curv, Some(&verf), Some(&verl)) else {
            continue;
        };

        let b = TopoBuilder::new();
        let mut new_edge = b.make_edge(curv.clone(), build.first, build.last);
        add_endpoint_vertices(&b, &mut new_edge, &build);
        // `cxx:792-793`: `B.Degenerated(E, degenerated); E.Closed(closed)`.
        if build.closed {
            new_edge.0.set_closed(true);
        }
        if build.degenerated {
            GeometryRegistry::global().set_degenerated(&new_edge.0, true);
        }

        // `cxx:1174`: `mapEdgLEdg.Bind(edgnew, LmapEdg)` — the map is the local
        // of `PerformEdge` (`cxx:667-672`), so the entry is dropped with it.
    }
}
