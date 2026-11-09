//! `LocOpe_Spliter::Perform`, stages 7-9
//! (`LocOpe_Spliter.cxx:268-577`).
//!
//! - stage 7 (`cxx:268-378`): re-orient the section edges that were projected
//!   onto an edge of the shape, then fold that substitution through the whole
//!   descendant table;
//! - stage 8 (`cxx:380-452`): chain the degenerated edges into wires and bind
//!   every vertex of such a wire to the first one, so a degenerated edge keeps a
//!   single vertex;
//! - stage 9 (`cxx:454-577`): `DirectLeft` / `Left`, the faces on the left of
//!   the sections, grown across the faces that carry no section.

use crate::abs::{Orientation, ShapeType};
use crate::brep_tool::BRepTool;
use crate::brep_tools_substitution::BRepToolsSubstitution;
use crate::builder::TopoBuilder;
use crate::locope::shape_maps::{ShapeDataMap, ShapeSet};
use crate::shape::{Edge, TopoShape, Wire};
use crate::topo_tools::{edge_vertices, explorer_shapes, is_same, map_shapes, wire_vertices};

use super::statics::select;
use super::Spliter;

impl Spliter {
    /// Stages 7-9 of `Perform`, run once the resulting shape is known.
    pub(super) fn finish(
        &mut self,
        res: &mut TopoShape,
        edg_on_edg: &ShapeDataMap,
        map_e: &ShapeSet,
    ) -> Result<(), String> {
        // ---- stage 7: edge-on-edge substitution (`cxx:268-346`) ------------
        let mut the_subs = BRepToolsSubstitution::new();
        let mut lsubs: Vec<TopoShape> = Vec::new();

        for (e1, e2) in edg_on_edg.iter() {
            let Some(e1) = Edge::wrap(e1.clone()) else {
                continue;
            };
            // `cxx:280-281`.
            let (Some(vf1), Some(vl1)) = edge_vertices(&e1) else {
                continue;
            };
            lsubs.clear();
            // `cxx:283`: the descendants of `e2`.
            for cand in self.my_map.find(e2).cloned().unwrap_or_default() {
                let Some(cand) = Edge::wrap(cand) else {
                    continue;
                };
                // `cxx:286`.
                let (Some(vf2), Some(vl2)) = edge_vertices(&cand) else {
                    continue;
                };
                // `cxx:288-303`.
                if !is_same(&vl1.0, &vf1.0) {
                    if is_same(&vf1.0, &vf2.0) && is_same(&vl1.0, &vl2.0) {
                        lsubs.push(cand.0.oriented(Orientation::Forward));
                    } else if is_same(&vf1.0, &vl2.0) && is_same(&vl1.0, &vf2.0) {
                        lsubs.push(cand.0.oriented(Orientation::Reversed));
                    }
                } else if is_same(&vf2.0, &vl2.0) && is_same(&vl2.0, &vl1.0) {
                    // `cxx:304-325`: a closed edge, decided by the tangent at the
                    // start parameter. `BRep_Tool::Curve` + `Loc.Transformation()`
                    // is `edge_curve_world`; `C->D1(f, ...)` is `d1(f)`.
                    let (Some(c1), Some(c2)) = (
                        BRepTool::edge_curve_world(&e1),
                        BRepTool::edge_curve_world(&cand),
                    ) else {
                        continue;
                    };
                    let (f1, _) = BRepTool::edge_parameters(&e1);
                    let (f2, _) = BRepTool::edge_parameters(&cand);
                    let v1 = c1.d1(f1).1;
                    let v2 = c2.d1(f2).1;
                    if v1.dot(&v2) > 0.0 {
                        lsubs.push(cand.0.oriented(Orientation::Forward));
                    } else {
                        lsubs.push(cand.0.oriented(Orientation::Reversed));
                    }
                }
            }
            // `cxx:327-330`.
            if lsubs.len() >= 2 {
                select(&e1, &mut lsubs);
            }
            // `cxx:332-338`: `e1` takes the orientation of the chosen section.
            if lsubs.len() == 1 {
                let ebase = lsubs[0].clone();
                let replacement = e1.0.oriented(ebase.orientation());
                let _ = the_subs.substitute(&ebase, &[replacement]);
            }
            // `cxx:339-345`: otherwise the reference only prints under
            // `OCCT_DEBUG` and leaves the edge alone.
        }

        // `cxx:347`: `theSubs.Build(myRes)`.
        the_subs.build(res);

        // `cxx:348-375`: rewrite every descendant list with the substitution.
        let keys: Vec<TopoShape> = self.my_map.iter().map(|(k, _)| k.clone()).collect();
        for k in keys {
            let ldesc = self.my_map.find(&k).cloned().unwrap_or_default();
            let mut newdesc: Vec<TopoShape> = Vec::with_capacity(ldesc.len());
            for item in &ldesc {
                if the_subs.is_copied(item) {
                    // `cxx:358-366`: `lsub.First()`; an empty answer would raise
                    // in the reference, so it is dropped here.
                    if let Some(first) = the_subs.copy(item).and_then(|l| l.first().cloned()) {
                        newdesc.push(first);
                    }
                } else {
                    // `cxx:369-371`.
                    newdesc.push(item.clone());
                }
            }
            self.my_map.bind(&k, newdesc);
        }

        // `cxx:377-380`.
        if the_subs.is_copied(res) {
            if let Some(first) = the_subs.copy(res).and_then(|l| l.first().cloned()) {
                *res = first;
            }
        }

        // ---- stage 8: degenerated edges (`cxx:380-452`) --------------------
        the_subs.clear();
        let builder = TopoBuilder::new();
        let emap = map_shapes(res, &[ShapeType::Edge]);
        let mut deg_edges: Vec<TopoShape> = Vec::new();
        for e in emap.get(&ShapeType::Edge).cloned().unwrap_or_default() {
            if let Some(ed) = Edge::wrap(e.clone()) {
                if BRepTool::is_degenerated(&ed) {
                    deg_edges.push(e);
                }
            }
        }

        // `cxx:396-431`: chain the degenerated edges into wire pieces.
        let mut deg_wires: Vec<TopoShape> = Vec::new();
        while !deg_edges.is_empty() {
            let mut deg_wire = Wire::new();
            builder.add(&mut deg_wire.0, &deg_edges[0]);
            deg_edges.remove(0);
            loop {
                // `cxx:410`: `TopExp::Vertices(aDegWire, Vfirst, Vlast)`.
                let (Some(vfirst), Some(vlast)) = wire_vertices(&deg_wire) else {
                    break;
                };
                let mut found = false;
                for i in 0..deg_edges.len() {
                    let Some(an_edge) = Edge::wrap(deg_edges[i].clone()) else {
                        continue;
                    };
                    let (Some(v1), Some(v2)) = edge_vertices(&an_edge) else {
                        continue;
                    };
                    // `cxx:417`.
                    if is_same(&v1.0, &vfirst.0)
                        || is_same(&v1.0, &vlast.0)
                        || is_same(&v2.0, &vfirst.0)
                        || is_same(&v2.0, &vlast.0)
                    {
                        builder.add(&mut deg_wire.0, &deg_edges[i]);
                        deg_edges.remove(i);
                        found = true;
                        break;
                    }
                }
                // `cxx:424-427`.
                if !found {
                    break;
                }
            }
            deg_wires.push(deg_wire.0);
        }

        // `cxx:433-449`: every vertex of a degenerated wire becomes the first.
        for w in &deg_wires {
            let vmap = map_shapes(w, &[ShapeType::Vertex]);
            let verts = vmap.get(&ShapeType::Vertex).cloned().unwrap_or_default();
            let Some(first) = verts.first() else {
                continue;
            };
            // `cxx:438`.
            let lv = vec![first.oriented(Orientation::Forward)];
            for v in verts.iter().skip(1) {
                // `cxx:441`.
                if !is_same(v, first) {
                    let _ = the_subs.substitute(v, &lv);
                }
            }
        }

        // `cxx:447-452`.
        the_subs.build(res);
        if the_subs.is_copied(res) {
            if let Some(first) = the_subs.copy(res).and_then(|l| l.first().cloned()) {
                *res = first;
            }
        }

        // ---- stage 9: DirectLeft / Left (`cxx:454-577`) --------------------
        // `cxx:454-456`.
        self.my_d_left.clear();
        self.my_left.clear();
        let mut map_v = ShapeSet::new();

        // `cxx:459-487`: a face touches a section when it carries one of the
        // `mapE` edges *with the same orientation*.
        for fac in explorer_shapes(res, ShapeType::Face) {
            let mut touches = false;
            for edg in explorer_shapes(&fac, ShapeType::Edge) {
                let mut hit = false;
                for key in map_e.iter() {
                    // `cxx:469`.
                    if is_same(key, &edg) && edg.orientation() == key.orientation() {
                        hit = true;
                        break;
                    }
                }
                if hit {
                    touches = true;
                    break;
                }
            }
            // `cxx:478-486`.
            if touches {
                self.my_d_left.push(fac.clone());
                self.my_left.push(fac);
            } else {
                map_v.add(&fac);
            }
        }

        // `cxx:522-533`: toggle the sections of the left faces.
        let mut mapebord = ShapeSet::new();
        for fac in &self.my_left {
            for edg in explorer_shapes(fac, ShapeType::Edge) {
                if !map_e.contains(&edg) && !mapebord.add(&edg) {
                    mapebord.remove(&edg);
                }
            }
        }

        // `cxx:535-573`: grow `myLeft` across the faces free of sections.
        while !mapebord.is_empty() {
            // `cxx:537-539`: the map's first key.
            let edg = mapebord.iter()[0].clone();
            // `cxx:541-556`: the first face of `mapV` that carries `edg`.
            let mut left_face: Option<TopoShape> = None;
            for fac in map_v.iter() {
                if explorer_shapes(fac, ShapeType::Edge)
                    .iter()
                    .any(|e| is_same(e, &edg))
                {
                    left_face = Some(fac.clone());
                    break;
                }
            }
            if let Some(fac) = left_face {
                // `cxx:558-563`.
                for e in explorer_shapes(&fac, ShapeType::Edge) {
                    if !mapebord.add(&e) {
                        mapebord.remove(&e);
                    }
                }
                // `cxx:564-566`.
                map_v.remove(&fac);
                self.my_left.push(fac);
            } else {
                // `cxx:568-571`.
                mapebord.remove(&edg);
            }
        }

        // `cxx:575-576`.
        self.my_res = Some(res.clone());
        self.my_done = true;
        Ok(())
    }
}
