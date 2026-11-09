//! `LocOpe_Spliter::Perform`, stages 1-6
//! (`LocOpe_Spliter.cxx:63-266`).

use crate::abs::{Orientation, ShapeType};
use crate::brep_tool::BRepTool;
use crate::brep_tools_substitution::BRepToolsSubstitution;
use crate::builder::TopoBuilder;
use crate::locope::shape_maps::{ShapeDataMap, ShapeIndexedMap, ShapeSet};
use crate::locope::split_shape::SplitShape;
use crate::locope::wires_on_shape::WiresOnShape;
use crate::shape::{Edge, Face, Shell, TopoShape, Vertex, Wire};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{explorer_shapes, is_same};

use super::statics::{put, rebuild_wires};
use super::Spliter;

impl Spliter {
    /// `Perform(PW)` (`cxx:63-266`): vertex substitution, section bookkeeping,
    /// wire rebuild, descendant table and the resulting shape. Stages 7-9 live
    /// in [`super::finish`] and run only when this stage produced a shape.
    pub fn perform(&mut self, pw: &mut WiresOnShape) -> Result<(), String> {
        // `cxx:66-70`.
        let Some(my_shape) = self.my_shape.clone() else {
            return Err("LocOpe_Spliter: null shape".to_string());
        };
        self.my_done = false;
        self.my_map.clear();
        self.my_res = None;

        // `cxx:72`: `Put(myShape, myMap)`.
        put(&my_shape, &mut self.my_map);

        // `cxx:74-78`.
        let mut map_v = ShapeSet::new();
        let mut map_e = ShapeSet::new();
        let mut edg_on_edg = ShapeDataMap::new();
        let mut map_fe: ShapeIndexedMap<Vec<TopoShape>> = ShapeIndexedMap::new();
        let mut the_subs = BRepToolsSubstitution::new();
        let mut lsubs: Vec<TopoShape> = Vec::new();
        let registry = GeometryRegistry::global();

        // The `PW` pairs the iterator hands out. `myMapEF` is filled by
        // `InitEdgeIterator` and never changes afterwards, so the two reference
        // loops (`cxx:88-115` and `cxx:153-198`) read the same snapshot.
        let mut pw_items: Vec<(Edge, Face, Option<Edge>)> = Vec::new();
        pw.init_edge_iterator()?;
        while pw.more_edge() {
            pw_items.push((pw.edge(), pw.on_face(), pw.on_edge()));
            pw.next_edge();
        }

        // ---- stage 1: vertex substitution (`cxx:82-120`) -------------------
        for (edg, _, _) in &pw_items {
            // `cxx:91`: `mapE.Add(edg)`.
            map_e.add(&edg.0);
            for vtx in explorer_shapes(&edg.0, ShapeType::Vertex) {
                let vtx = Vertex(vtx);
                // `cxx:95`.
                if map_v.contains(&vtx.0) {
                    continue;
                }
                // `cxx:97`.
                let Some(vb) = pw.on_vertex(&vtx) else {
                    continue;
                };
                // `cxx:99`.
                map_v.add(&vtx.0);
                // `cxx:100-103`.
                if is_same(&vtx.0, &vb.0) {
                    continue;
                }
                lsubs.clear();
                // `cxx:105-111`.
                let vsub = vtx.0.oriented(Orientation::Forward);
                let p1 = BRepTool::vertex_point_world(&Vertex(vsub.clone()));
                let p2 = BRepTool::vertex_point_world(&vb);
                let d = p1.distance(&p2) + BRepTool::vertex_tolerance(&vb);
                // `BB.UpdateVertex(vsub, d)` (`cxx:109`).
                registry.set_vertex_tolerance(&vsub, d);
                lsubs.push(vsub);
                // `cxx:111`: `theSubs.Substitute(Vb.Oriented(FORWARD), lsubs)`.
                let _ = the_subs.substitute(&vb.0.oriented(Orientation::Forward), &lsubs);
            }
        }

        // ---- `cxx:118-141`: fold the substitution into `myMap` ------------
        the_subs.build(&my_shape);
        let keys: Vec<TopoShape> = self.my_map.iter().map(|(k, _)| k.clone()).collect();
        for k in keys {
            if !the_subs.is_copied(&k) {
                continue;
            }
            // `cxx:135-137`: `Clear()` then `Append(lsub.First())`. An empty
            // `lsub` would raise in `NCollection_List::First`; it cannot happen
            // for the one-element substitutions above, so it is skipped.
            let Some(first) = the_subs.copy(&k).and_then(|l| l.first().cloned()) else {
                continue;
            };
            self.my_map.bind(&k, vec![first]);
        }

        // `cxx:143`: `myRes = myMap(myShape).First()`.
        let mut res = self
            .my_map
            .find(&my_shape)
            .and_then(|l| l.first().cloned())
            .ok_or_else(|| "LocOpe_Spliter: the shape has no descendant".to_string())?;
        // `cxx:144`: `LocOpe_SplitShape theCFace(myRes)`.
        let mut cface = SplitShape::with_shape(&res);

        // ---- stage 3: sections on faces and edges (`cxx:150-199`) -----------
        let mut the_faces_with_section = ShapeSet::new();
        for (edg, fac_on, ebis) in &pw_items {
            // `cxx:155-175`: vertices on edges.
            for vtx in explorer_shapes(&edg.0, ShapeType::Vertex) {
                let vtx = Vertex(vtx);
                // `cxx:158-159`.
                if map_v.contains(&vtx.0) {
                    continue;
                }
                map_v.add(&vtx.0);
                let Some((ed, prm)) = pw.on_edge_of_vertex_from(&vtx, edg) else {
                    continue;
                };
                // `cxx:166-168`.
                if !self.my_map.is_bound(&ed.0) {
                    continue;
                }
                let Some(ed) = self
                    .my_map
                    .find(&ed.0)
                    .and_then(|l| l.first().cloned())
                    .and_then(Edge::wrap)
                else {
                    continue;
                };
                // `cxx:169`.
                cface.add_vertex(&vtx, prm, &ed);
            }

            // `cxx:176-198`.
            if let Some(ebis) = ebis {
                edg_on_edg.bind(&edg.0, &ebis.0);
                continue;
            }
            let fac = fac_on;
            // `cxx:181-184`.
            if !self.my_map.is_bound(&fac.0) {
                continue;
            }
            let is_face_with_sec = pw.is_face_with_section(&fac.0);
            let Some(fac) = self
                .my_map
                .find(&fac.0)
                .and_then(|l| l.first().cloned())
                .and_then(Face::wrap)
            else {
                continue;
            };
            // `cxx:187-190`.
            if is_face_with_sec {
                the_faces_with_section.add(&fac.0);
            }
            // `cxx:191-198`: `mapFE(fac).Append(edg)`.
            if !map_fe.contains(&fac.0) {
                map_fe.add(&fac.0, Vec::new());
            }
            if let Some(l) = map_fe.change_from_key(&fac.0) {
                l.push(edg.0.clone());
            }
        }

        // ---- stage 4: rebuild the wires of every touched face (`cxx:202-222`)
        for i in 1..=map_fe.extent() {
            let Some(fac) = Face::wrap(map_fe.find_key(i).clone()) else {
                continue;
            };
            let ledges = map_fe.find_from_index(i).clone();
            let new_list = rebuild_wires(&ledges, pw)?;
            if the_faces_with_section.contains(&fac.0) {
                // `cxx:212-214`: `theCFace.Add(ledges, fac)`.
                cface.add_wires(&new_list, &fac);
            } else {
                // `cxx:216-221`.
                for w in &new_list {
                    if let Some(wire) = Wire::wrap(w.clone()) {
                        cface.add_wire(&wire, &fac);
                    }
                }
            }
        }

        // ---- stage 5: descendants (`cxx:224-234`) --------------------------
        let keys: Vec<TopoShape> = self.my_map.iter().map(|(k, _)| k.clone()).collect();
        for k in keys {
            let Some(first) = self.my_map.find(&k).and_then(|l| l.first().cloned()) else {
                continue;
            };
            let d = cface.descendant_shapes(&first);
            self.my_map.bind(&k, d);
        }

        // ---- stage 6: the resulting shape (`cxx:236-267`) ------------------
        // `cxx:236`: `myMap(myShape)`.
        let lres = self.my_map.find(&my_shape).cloned().unwrap_or_default();
        let typ = my_shape.shape_type();
        let builder = TopoBuilder::new();
        if typ == ShapeType::Face && lres.len() >= 2 {
            // `cxx:240-250`: a shell of the pieces, in place of the face.
            let mut shell = Shell::new();
            for s in &lres {
                builder.add(&mut shell.0, &s.oriented(my_shape.orientation()));
            }
            res = shell.0;
        } else if typ == ShapeType::Edge && lres.len() >= 2 {
            // `cxx:251-262`: a wire of the pieces, in place of the edge.
            let mut wire = Wire::new();
            for s in &lres {
                builder.add(&mut wire.0, &s.oriented(my_shape.orientation()));
            }
            res = wire.0;
        } else {
            // `cxx:263-267`: exactly one descendant, else the reference returns
            // with `myDone` still `false`.
            if lres.len() != 1 {
                return Ok(());
            }
            res = lres[0].clone();
        }
        self.my_res = Some(res.clone());

        // ---- stages 7-9 (`cxx:268-577`) -----------------------------------
        self.finish(&mut res, &edg_on_edg, &map_e)
    }
}
