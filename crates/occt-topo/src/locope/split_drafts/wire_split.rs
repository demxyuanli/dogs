//! The `theW` split stage of `LocOpe_SplitDrafts::Perform`
//! (`LocOpe_SplitDrafts.cxx:177-351`).
//!
//! The two draft planes `NewSg` / `NewSd` meet in a line. Sweeping the wire
//! `W` along that line (`GeomFill_Pipe`) gives a surface that every edge of `W`
//! crosses exactly where it has to be cut for the drafted walls to meet; the
//! crossings are inserted into a `LocOpe_SplitShape` and the resulting pieces
//! replace `W`'s edges in `myMap`.

use occt_core::precision::CONFUSION;

use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::locope::split_drafts::perform::{intcs_perform, pipe_surface};
use crate::locope::split_drafts::perform::edges_of_shape;
use crate::locope::split_drafts::state::PerformState;
use crate::locope::split_drafts::SplitDrafts;
use crate::locope::SplitShape;
use crate::shape::{Edge, Face, TopoShape, Vertex, Wire};
use crate::topo_tools::{self, edge_vertices, is_same};

impl SplitDrafts {
    /// `cxx:177-351`. `Ok(false)` is one of the reference's bare `return`s,
    /// which leave `myResult` null.
    pub(crate) fn split_wire(
        &mut self,
        w: &Wire,
        _f: &Face,
        st: &mut PerformState,
    ) -> Result<bool, String> {
        if st.i2s.is_done() && st.i2s.nb_lines() > 0 {
            // `cxx:186-195`: the pipe surface along the plane/plane line.
            let Some(line1) = st.i2s.line(1).map(|l| l.curve.clone()) else {
                return Ok(true);
            };
            let spl = pipe_surface(&st.line_pipe, line1)?;

            // `cxx:197`.
            let mut splw = SplitShape::with_shape(&w.0);

            // `cxx:200-296`.
            let b = TopoBuilder::new();
            for edg in edges_of_shape(&w.0) {
                // `cxx:202`: `if (theMap.Add(edg))`.
                if !st.the_map.add(&edg) {
                    continue;
                }
                let edge = Edge(edg.clone());
                let Some(curve) = BRepTool::edge_curve(&edge) else {
                    continue;
                };
                let (f, l) = BRepTool::edge_parameters(&edge);
                // `cxx:208-213`: `intcs.Perform(HAC, HAS)`.
                let Some(intcs) = intcs_perform(curve.as_ref(), spl.as_ref()) else {
                    continue;
                };

                // `cxx:216-218`: two or more segments is not implemented.
                if intcs.nb_segments() >= 2 {
                    continue;
                }

                if intcs.nb_segments() == 1 {
                    // `cxx:221-260`.
                    let p1 = intcs.segment(0).first_point();
                    let p2 = intcs.segment(0).second_point();
                    let pf = p1.pnt();
                    let pl = p2.pnt();
                    let (v_f, v_l) = edge_vertices(&edge);
                    let (Some(v_f), Some(v_l)) = (v_f, v_l) else {
                        continue;
                    };
                    let pf_e = BRepTool::vertex_point(&v_f);
                    let pl_e = BRepTool::vertex_point(&v_l);
                    let mut tol_f = BRepTool::vertex_tolerance(&v_f);
                    let mut tol_l = BRepTool::vertex_tolerance(&v_l);
                    tol_f *= tol_f;
                    tol_l *= tol_l;

                    let dff = pf.square_distance(&pf_e);
                    let dfl = pf.square_distance(&pl_e);
                    let dlf = pl.square_distance(&pf_e);
                    let dll = pl.square_distance(&pl_e);

                    if (dff <= tol_f && dll <= tol_l) || (dlf <= tol_f && dfl <= tol_l) {
                        continue;
                    }
                    // `cxx:246-259`.
                    let v_new_f = b.make_vertex(pf, CONFUSION);
                    let v_new_l = b.make_vertex(pl, CONFUSION);
                    if p1.param() >= f
                        && p1.param() <= l
                        && p2.param() >= f
                        && p2.param() <= l
                    {
                        splw.add_vertex(&v_new_f, p1.param(), &edge);
                        splw.add_vertex(&v_new_l, p2.param(), &edge);
                    } else {
                        continue;
                    }
                } else if intcs.nb_points() != 0 {
                    // `cxx:261-294`.
                    let (v_f, v_l) = edge_vertices(&edge);
                    let (Some(v_f), Some(v_l)) = (v_f, v_l) else {
                        continue;
                    };
                    let pf = BRepTool::vertex_point(&v_f);
                    let pl = BRepTool::vertex_point(&v_l);
                    let mut tol_f = BRepTool::vertex_tolerance(&v_f);
                    let mut tol_l = BRepTool::vertex_tolerance(&v_l);
                    tol_f *= tol_f;
                    tol_l *= tol_l;

                    for i in 0..intcs.nb_points() {
                        let pi = intcs.point(i);
                        let p = pi.pnt();
                        let dif = p.square_distance(&pf);
                        let dil = p.square_distance(&pl);
                        if dif <= tol_f {
                            // `cxx:279-281`: an empty arm in the reference.
                        } else if dil <= tol_l {
                            // `cxx:281-283`: an empty arm in the reference.
                        } else if pi.param() >= f && pi.param() <= l {
                            let v_new: Vertex = b.make_vertex(p, CONFUSION);
                            splw.add_vertex(&v_new, pi.param(), &edge);
                        }
                    }
                }
            }

            // `cxx:298-306`.
            let lres = splw.descendant_shapes(&w.0);
            if lres.len() != 1 {
                return Ok(false);
            }
            if !is_same(&w.0, &lres[0]) {
                st.the_w = Wire(lres[0].clone());
            }

            // `cxx:310-330`: every edge of `W` (and its vertices) that the
            // first loop did not bind gets the split's descendants.
            for edg in edges_of_shape(&w.0) {
                if self.my_map.is_bound(&edg) {
                    continue;
                }
                let mut list: Vec<TopoShape> = Vec::new();
                for d in splw.descendant_shapes(&edg) {
                    list.push(d);
                }
                self.my_map.bind(&edg, list);
                for v in topo_tools::vertices_of(&edg) {
                    if self.my_map.is_bound(&v.0) {
                        continue;
                    }
                    self.my_map.bind(&v.0, Vec::new());
                    if let Some(slot) = self.my_map.change(&v.0) {
                        slot.push(v.0.clone());
                    }
                }
            }
        } else {
            // `cxx:332-351`.
            for edg in edges_of_shape(&w.0) {
                if self.my_map.is_bound(&edg) {
                    continue;
                }
                self.my_map.bind(&edg, vec![edg.clone()]);
                for v in topo_tools::vertices_of(&edg) {
                    if self.my_map.is_bound(&v.0) {
                        continue;
                    }
                    self.my_map.bind(&v.0, vec![v.0.clone()]);
                }
            }
        }
        Ok(true)
    }
}
