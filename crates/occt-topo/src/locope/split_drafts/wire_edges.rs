//! Reconstruction of the boundary edges of `theW`
//! (`LocOpe_SplitDrafts.cxx:868-1100`).
//!
//! Each edge of the split wire `theW` gets one image on each draft plane
//! (`NewEdge`). When the two images collapse onto the same curve the pair
//! becomes a single edge (or even the original edge, when nothing moved);
//! otherwise the two images plus the `MapEV` faces hanging off the endpoints
//! close into a 2/3/4-sided face through `MakeFace`.

use crate::abs::{Orientation, ShapeType};
use crate::brep_tool::BRepTool;
use crate::locope::split_drafts::perform::edges_of_shape;
use crate::locope::split_drafts::state::PerformState;
use crate::locope::split_drafts::statics::{make_face, new_edge};
use crate::locope::split_drafts::SplitDrafts;
use crate::shape::{Edge, TopoShape, Vertex};
use crate::topo_tools::{edge_vertices, is_same};

impl SplitDrafts {
    /// `cxx:868-1100`. `Ok(false)` is the reference's bare `return`.
    pub(crate) fn build_wire_edges(&mut self, st: &mut PerformState) -> Result<bool, String> {
        // `cxx:868-870`.
        st.the_map.clear();
        for edg in edges_of_shape(&st.the_w.0) {
            if !st.the_map.add(&edg) {
                continue;
            }
            // `cxx:876-878`: `MapEV(indedg)`; the reference raises when the edge
            // is absent.
            let Some(indedg) = st.map_ev.find_index(&edg) else {
                return Ok(false);
            };
            let gen_f0: TopoShape = st.map_ev.find_from_index(indedg).clone();
            // `cxx:879-880`.
            st.map_sg.bind(&edg, Vec::new());
            st.map_sd.bind(&edg, Vec::new());

            // `cxx:881-888`.
            let forward = edg.oriented(Orientation::Forward);
            let (vf, vl) = edge_vertices(&Edge(forward));
            let (Some(vf), Some(vl)) = (vf, vl) else {
                return Ok(false);
            };
            let gvf = st.map_ev.find_from_key(&vf.0).cloned();
            let gvl = st.map_ev.find_from_key(&vl.0).cloned();
            let (Some(gvf), Some(gvl)) = (gvf, gvl) else {
                return Ok(false);
            };
            let gvf_is_vertex = gvf.shape_type() == ShapeType::Vertex;
            let gvl_is_vertex = gvl.shape_type() == ShapeType::Vertex;

            // `cxx:999-1019`.
            let first_of = |map: &crate::locope::shape_maps::ShapeValueMap<Vec<TopoShape>>,
                            key: &TopoShape|
             -> Option<Vertex> { map.find(key)?.first().map(|s| Vertex(s.clone())) };

            let (vfg, vfd) = if gvf_is_vertex {
                let v = Vertex(gvf.clone());
                (v.clone(), v)
            } else {
                match (first_of(&st.map_sg, &vf.0), first_of(&st.map_sd, &vf.0)) {
                    (Some(a), Some(b)) => (a, b),
                    _ => return Ok(false),
                }
            };
            let (vlg, vld) = if gvl_is_vertex {
                let v = Vertex(gvl.clone());
                (v.clone(), v)
            } else {
                match (first_of(&st.map_sg, &vl.0), first_of(&st.map_sd, &vl.0)) {
                    (Some(a), Some(b)) => (a, b),
                    _ => return Ok(false),
                }
            };

            // `cxx:1021-1031`.
            let Some(new_edgg) = new_edge(&Edge(edg.clone()), &gen_f0, &st.new_sg, &vfg, &vlg)
            else {
                return Ok(false);
            };
            let Some(new_edgd) = new_edge(&Edge(edg.clone()), &gen_f0, &st.new_sd, &vfd, &vld)
            else {
                return Ok(false);
            };

            // `cxx:1032-1078`.
            let mut isedg = false;
            if gvf_is_vertex && gvl_is_vertex {
                let (fg, lg) = BRepTool::edge_parameters(&new_edgg);
                let (fd, ld) = BRepTool::edge_parameters(&new_edgd);
                let (Some(cg), Some(cd)) = (
                    BRepTool::edge_curve(&new_edgg),
                    BRepTool::edge_curve(&new_edgd),
                ) else {
                    return Ok(false);
                };
                let pg = cg.d0((fg + lg) / 2.0);
                let pd = cd.d0((fd + ld) / 2.0);
                // `cxx:1047`: the reference takes `NewEdgg`'s tolerance twice.
                let tol = BRepTool::edge_tolerance(&new_edgg)
                    .max(BRepTool::edge_tolerance(&new_edgg));
                if pg.square_distance(&pd) <= tol * tol {
                    isedg = true;
                    // `cxx:1050-1062`.
                    let mut modified = true;
                    if is_same(&gvf, &vf.0) && is_same(&gvl, &vl.0) {
                        let (ed, el_d) = BRepTool::edge_parameters(&Edge(edg.clone()));
                        let Some(orig) = BRepTool::edge_curve(&Edge(edg.clone())) else {
                            return Ok(false);
                        };
                        let pd2 = orig.d0((ed + el_d) / 2.0);
                        let tol2 = BRepTool::edge_tolerance(&new_edgg)
                            .max(BRepTool::edge_tolerance(&Edge(edg.clone())));
                        if pg.square_distance(&pd2) <= tol2 * tol2 {
                            modified = false;
                        }
                    }
                    // `cxx:1063-1078`.
                    if modified {
                        st.map_w.bind(&edg, Some(new_edgg.0.clone()));
                        st.map_sg.bind(&edg, vec![new_edgg.0.clone()]);
                        st.map_sd.bind(&edg, vec![new_edgg.0.clone()]);
                    } else {
                        st.map_w.bind(&edg, Some(edg.clone()));
                        st.map_sg.bind(&edg, vec![edg.clone()]);
                        st.map_sd.bind(&edg, vec![edg.clone()]);
                    }
                }
            }

            // `cxx:1080-1099`.
            if !isedg {
                st.map_sg.bind(&edg, vec![new_edgg.0.clone()]);
                st.map_sd.bind(&edg, vec![new_edgd.0.clone()]);
                let mut theedges: Vec<TopoShape> = vec![new_edgg.0.clone(), new_edgd.0.clone()];
                if gvf.shape_type() == ShapeType::Edge {
                    theedges.push(gvf.clone());
                }
                if gvl.shape_type() == ShapeType::Edge {
                    theedges.push(gvl.clone());
                }
                let mut gen_f = gen_f0.clone();
                make_face(&mut gen_f, &mut theedges);
                // `MapEV(indedg)` is the same slot `MakeFace` filled.
                *st.map_ev.change_from_index(indedg) = gen_f.clone();
                st.map_w.bind(&edg, Some(gen_f));
            }
        }
        Ok(true)
    }
}
