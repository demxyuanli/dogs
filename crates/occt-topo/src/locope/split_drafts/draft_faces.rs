//! Rebuilding the split pieces of `F` on the draft planes
//! (`LocOpe_SplitDrafts.cxx:1102-1358`).
//!
//! Every face the split produced under `F` that lies on the left of the wire is
//! recreated on `NewSg`, the others on `NewSd`, wire by wire. Edges that sit on
//! the original wire take their image from `MapSg` / `MapSd`; the rest get a
//! fresh section edge unless the `MapEV` vertices still need the recovered edge
//! appended (`MapSonS`), which is what `BRepTools_Substitution` then swaps in the
//! rest of `Res`.

use crate::abs::{Orientation, ShapeType};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::locope::shape_maps::ShapeSet;
use crate::locope::split_drafts::perform::reverse_orientation;
use crate::locope::split_drafts::state::PerformState;
use crate::locope::split_drafts::statics::{contains, new_edge};
use crate::locope::split_drafts::SplitDrafts;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::topo_tools::{edge_vertices, explorer_shapes, is_same};

impl SplitDrafts {
    /// `cxx:1102-1358`. `Ok(false)` is the reference's bare `return`.
    pub(crate) fn build_draft_faces(
        &mut self,
        f: &Face,
        st: &mut PerformState,
    ) -> Result<bool, String> {
        let b = TopoBuilder::new();
        let mut mapedgadded = ShapeSet::new();
        let split_faces: Vec<TopoShape> = self.my_map.find(&f.0).cloned().unwrap_or_default();
        let split_faces_key = split_faces.clone();

        for fac in split_faces {
            // `cxx:1108`: `theMap` is reused as a per-face "seen" set.
            st.the_map.clear();

            // `cxx:1109-1121`.
            let is_left = contains(&st.the_left, &fac);
            let new_s = if is_left { st.new_sg.clone() } else { st.new_sd.clone() };
            let mut drft_face = b
                .make_face_with_tolerance(
                    new_s.clone(),
                    &[],
                    BRepTool::face_tolerance(&Face(fac.clone())),
                )
                .0;

            // `cxx:1123-1124`.
            let forward_fac = fac.oriented(Orientation::Forward);
            for wir in explorer_shapes(&forward_fac, ShapeType::Wire) {
                let mut new_wire_edges: Vec<Edge> = Vec::new();
                let forward_wir = wir.oriented(Orientation::Forward);
                for edg in explorer_shapes(&forward_wir, ShapeType::Edge) {
                    // `cxx:1131-1134`.
                    if !st.the_map.add(&edg) {
                        continue;
                    }
                    if st.map_w.is_bound(&edg) {
                        // `cxx:1136-1155`: an edge of the original wire.
                        let mut ored = edg.orientation();
                        let list = if is_left {
                            st.map_sg.find(&edg).cloned().unwrap_or_default()
                        } else {
                            st.map_sd.find(&edg).cloned().unwrap_or_default()
                        };
                        for v in list {
                            if v.orientation() == Orientation::Reversed {
                                ored = reverse_orientation(ored);
                            }
                            new_wire_edges.push(Edge(v.oriented(ored)));
                        }
                        continue;
                    }

                    // `cxx:1158-1168`.
                    let Some(indedg) = st.map_ev.find_index(&edg) else {
                        return Ok(false);
                    };
                    let gen_f: TopoShape = st.map_ev.find_from_index(indedg).clone();
                    let forward = edg.oriented(Orientation::Forward);
                    let (vf, vl) = edge_vertices(&Edge(forward));
                    let (Some(vf), Some(vl)) = (vf, vl) else {
                        return Ok(false);
                    };
                    let (Some(gvf), Some(gvl)) = (
                        st.map_ev.find_from_key(&vf.0).cloned(),
                        st.map_ev.find_from_key(&vl.0).cloned(),
                    ) else {
                        return Ok(false);
                    };
                    let gvf_is_vertex = gvf.shape_type() == ShapeType::Vertex;
                    let gvl_is_vertex = gvl.shape_type() == ShapeType::Vertex;

                    let new_edg: Edge;
                    if gvf_is_vertex && gvl_is_vertex {
                        // `cxx:1178-1201`.
                        if !is_same(&gvf, &vf.0) || !is_same(&gvl, &vl.0) {
                            let Some(e) = new_edge(
                                &Edge(edg.clone()),
                                &gen_f,
                                &new_s,
                                &Vertex(gvf.clone()),
                                &Vertex(gvl.clone()),
                            ) else {
                                return Ok(false);
                            };
                            st.map_son_s.bind(&edg, e.0.clone());
                            new_edg = e;
                        } else {
                            // `cxx:1200-1202`: frozen, the edge is reused as is.
                            new_wire_edges.push(Edge(edg.clone()));
                            continue;
                        }
                    } else {
                        // `cxx:1206-1238`.
                        let vff = if gvf_is_vertex {
                            Vertex(gvf.clone())
                        } else {
                            let map = if is_left { &st.map_sg } else { &st.map_sd };
                            match map.find(&vf.0).and_then(|l| l.first()) {
                                Some(s) => Vertex(s.clone()),
                                None => return Ok(false),
                            }
                        };
                        let vll = if gvl_is_vertex {
                            Vertex(gvl.clone())
                        } else {
                            let map = if is_left { &st.map_sg } else { &st.map_sd };
                            match map.find(&vl.0).and_then(|l| l.first()) {
                                Some(s) => Vertex(s.clone()),
                                None => return Ok(false),
                            }
                        };
                        let Some(e) = new_edge(
                            &Edge(edg.clone()),
                            &gen_f,
                            &new_s,
                            &vff,
                            &vll,
                        ) else {
                            return Ok(false);
                        };
                        new_edg = e;

                        // `cxx:1245-1320`.
                        let wvf = st.map_w.is_bound(&vf.0);
                        let wvl = st.map_w.is_bound(&vl.0);
                        if !wvf && !wvl {
                            st.map_son_s.bind(&edg, new_edg.0.clone());
                        } else if wvf {
                            if gvf.shape_type() != ShapeType::Edge
                                || mapedgadded.contains(&gvf)
                            {
                                st.map_son_s.bind(&edg, new_edg.0.clone());
                            } else {
                                // `cxx:1262-1287`: the recovered edge closes the
                                // wire with the new section edge.
                                let mut new_wir = b.make_wire(&[new_edg.clone()]);
                                let (_vf2, vl2) = match edge_vertices(&Edge(gvf.clone())) {
                                    (Some(a), Some(c)) => (a, c),
                                    _ => return Ok(false),
                                };
                                let oriented = if is_same(&vl2, &vff.0) {
                                    gvf.oriented(Orientation::Forward)
                                } else {
                                    gvf.oriented(Orientation::Reversed)
                                };
                                b.add_edge(&mut new_wir, &Edge(oriented));
                                mapedgadded.add(&gvf);
                                st.map_son_s.bind(&edg, new_wir.0.clone());
                            }
                        } else if gvl.shape_type() != ShapeType::Edge
                            || mapedgadded.contains(&gvl)
                        {
                            st.map_son_s.bind(&edg, new_edg.0.clone());
                        } else {
                            // `cxx:1294-1317`.
                            let mut new_wir = b.make_wire(&[new_edg.clone()]);
                            let (vf2, _vl2) = match edge_vertices(&Edge(gvl.clone())) {
                                (Some(a), Some(c)) => (a, c),
                                _ => return Ok(false),
                            };
                            let oriented = if is_same(&vf2, &vll.0) {
                                gvl.oriented(Orientation::Forward)
                            } else {
                                gvl.oriented(Orientation::Reversed)
                            };
                            b.add_edge(&mut new_wir, &Edge(oriented));
                            mapedgadded.add(&gvl);
                            st.map_son_s.bind(&edg, new_wir.0.clone());
                        }
                    }

                    // `cxx:1321-1328`.
                    let oriented = if new_edg.0.orientation() == Orientation::Reversed {
                        new_edg.0.oriented(reverse_orientation(edg.orientation()))
                    } else {
                        new_edg.0.oriented(edg.orientation())
                    };
                    new_wire_edges.push(Edge(oriented));
                }

                // `cxx:1332`.
                let new_wir = b.make_wire(&new_wire_edges);
                let oriented = new_wir.0.oriented(wir.orientation());
                b.add(&mut drft_face, &oriented);
            }

            // `cxx:1335`.
            st.the_faces.push(drft_face);
        }

        // `cxx:1337-1348`: `BRepTools_Substitution::Substitute`.
        for (key, value) in st.map_son_s.iter() {
            let lsubs: Vec<TopoShape> = explorer_shapes(value, ShapeType::Edge);
            st.subs.substitute(key, &lsubs);
        }

        // `cxx:1350-1357`.
        let res = match &st.res {
            Some(r) => r.clone(),
            None => return Ok(false),
        };
        for face in explorer_shapes(&res, ShapeType::Face) {
            if contains(&split_faces_key, &face) {
                continue;
            }
            st.subs.build(&face);
        }
        Ok(true)
    }
}
