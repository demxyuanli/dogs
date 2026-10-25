//! `TopOpeBRepBuild_FuseFace::PerformEdge` (`cxx:590-681`), `ClearEdge`
//! (`cxx:688-822`) and `ClearVertex` (`cxx:829-841`).
//!
//! `PerformEdge` fuses the co-surfacic edges of the already fused faces: it
//! collects `vertex -> the edges meeting there` (`cxx:604-640`), keeps the
//! vertices that belong to exactly two same-support edges (`cxx:642-662`),
//! chains them (`GroupEdge`) and rebuilds the chained edges (`MakeEdge`).
//!
//! `ClearEdge` drops the INTERNAL / EXTERNAL edges of every face and rebuilds
//! each face from its remaining wires with `BRepLib_MakeFace(wire, false)`
//! (`cxx:773-795`). `ClearVertex` is an empty body in OCCT — its only content is
//! `OCCT_DEBUG` tracing (`cxx:830-842`) — so the port keeps only the `Done`
//! flag side effects of the caller.

use crate::abs::{ShapeType, Orientation};
use crate::brep_lib_make_face::make_face_from_wire;
use crate::builder::TopoBuilder;
use crate::shape::{TopoShape, Wire};
use crate::topo_tools::explorer_shapes;

use super::group_edge::group_edge;
use super::make_edge::make_edge;
use super::perform_face::make_wire_from_list;
use super::same_support::same_support;
use super::shape_list_map::ShapeListMap;
use super::FuseFace;

impl FuseFace {
    /// `PerformEdge()` (`cxx:590-681`).
    pub fn perform_edge(&mut self) {
        // `cxx:604-640`: `mapVerLEdg` — one edge occurrence per vertex
        // occurrence, so an edge whose two vertices repeat is listed twice.
        let mut map_ver_l_edg = ShapeListMap::new();
        let lff = self.my_lff.clone();
        for fac in &lff {
            for wir in explorer_shapes(fac, ShapeType::Wire) {
                for edg in explorer_shapes(&wir, ShapeType::Edge) {
                    for ver in explorer_shapes(&edg, ShapeType::Vertex) {
                        if map_ver_l_edg.is_bound(&ver) {
                            map_ver_l_edg.append(&ver, &edg);
                        } else {
                            map_ver_l_edg.bind(&ver, vec![edg.clone()]);
                        }
                    }
                }
            }
        }

        // `cxx:642-663`: only the vertices with exactly two same-support edges.
        let map_tampon = map_ver_l_edg.clone();
        map_ver_l_edg.clear();
        for (ver, lmap_edg) in map_tampon.pairs() {
            if lmap_edg.len() != 2 {
                continue;
            }
            if same_support(&lmap_edg[0], &lmap_edg[1]) {
                map_ver_l_edg.bind(ver, lmap_edg.to_vec());
            }
        }

        // `cxx:665-671`.
        let mut map_edg_l_edg = ShapeListMap::new();
        group_edge(&map_ver_l_edg, &mut map_edg_l_edg);
        make_edge(&map_edg_l_edg);

        // `cxx:673-674`.
        self.my_modified = true;
        self.my_done = true;
    }

    /// `ClearEdge()` (`cxx:688-822`).
    pub fn clear_edge(&mut self) {
        // `cxx:699`: the rebuilt faces replace `myLFF` only when all of them
        // were built.
        let mut my_lff_new: Vec<TopoShape> = Vec::new();

        let lff = self.my_lff.clone();
        for fac in &lff {
            // `cxx:706-760`: rebuild the wires without the INTERNAL /
            // EXTERNAL edges, which are collected on the side.
            let mut my_face_lw: Vec<Wire> = Vec::new();
            for wir in explorer_shapes(fac, ShapeType::Wire) {
                let mut my_wire_le: Vec<TopoShape> = Vec::new();
                for edg in explorer_shapes(&wir, ShapeType::Edge) {
                    match edg.orientation() {
                        Orientation::Internal => self.my_lie.push(edg),
                        Orientation::External => self.my_lee.push(edg),
                        _ => my_wire_le.push(edg),
                    }
                }
                if my_wire_le.is_empty() {
                    continue;
                }
                let Some(w) = make_wire_from_list(&my_wire_le) else {
                    self.no_fuse();
                    return;
                };
                my_face_lw.push(w);
            }

            // `cxx:764-776`.
            if my_face_lw.is_empty() {
                self.no_fuse();
                return;
            }

            // `cxx:778-795`: `BRepLib_MakeFace MF(wir, OnlyPlane = false)` and
            // `MF.Add(wir1)` for the remaining wires.
            let Some(mut face) = make_face_from_wire(&my_face_lw[0], false) else {
                self.no_fuse();
                return;
            };
            let b = TopoBuilder::new();
            for w in &my_face_lw[1..] {
                b.add(&mut face.0, &w.0);
            }
            my_lff_new.push(face.0);
        }

        // `cxx:800-813`.
        if my_lff_new.is_empty() {
            self.no_fuse();
            return;
        }
        self.my_lff = my_lff_new;

        // `cxx:815-816`.
        self.my_modified = true;
        self.my_done = true;
    }

    /// `ClearVertex()` (`cxx:826-843`): an empty body in OCCT — the only
    /// statements are the `OCCT_DEBUG` traces at `cxx:830-842`, and the two
    /// `myLIV` / `myLEV` / `myLMV` lists it is named after are never written
    /// anywhere in the class (`cxx:110-114` only clears them in `Init`).
    pub fn clear_vertex(&mut self) {}
}
