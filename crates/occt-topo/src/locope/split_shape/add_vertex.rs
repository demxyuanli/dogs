//! `LocOpe_SplitShape::Add(const TopoDS_Vertex& V, const double P, const TopoDS_Edge& E)`
//! (`LocOpe_SplitShape.cxx:143-249`).
//!
//! Splits one descendant edge at parameter `P` by inserting vertex `V`, and
//! records the two halves in `myMap` in place of the original edge.

use crate::abs::{Orientation, ShapeType};
use crate::brep_builder_update_vertex::update_vertex_on_edge;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::locope::split_shape::{empty_copied_edge, SplitShape};
use crate::shape::{Edge, TopoShape, Vertex};
use crate::topo_tools::edge_vertices;

impl SplitShape {
    /// `Add(V, P, E)` (`cxx:143-249`).
    ///
    /// The reference writes straight into `myMap`, and so does this method: the
    /// split only becomes visible once `Rebuild` runs.
    pub fn add_vertex(&mut self, v: &Vertex, p: f64, e: &Edge) {
        // `cxx:145-148`.
        if !self.can_split(e) {
            return;
        }

        let b = TopoBuilder::new();

        // `cxx:150-152`: `NCollection_List<TopoDS_Shape>& le = myMap(E)`.
        let mut le = self.map_of(&e.0);

        // `cxx:153-156`.
        if le.is_empty() {
            le.push(e.0.clone());
        }

        // `cxx:158-168`: the first descendant whose range strictly contains `P`.
        let tol = occt_core::precision::PCONFUSION;
        let mut idx = 0usize;
        let mut found = false;
        while idx < le.len() {
            let Some(edg) = Edge::wrap(le[idx].clone()) else {
                break;
            };
            let (f, l) = BRepTool::edge_parameters(&edg);
            if p > f + tol && p < l - tol {
                found = true;
                break;
            }
            idx += 1;
        }
        // `cxx:169-172`.
        if !found {
            self.set_map_of(&e.0, le);
            return;
        }

        // `cxx:173-175`: `le.Remove(itl)`.
        let edg = Edge::wrap(le[idx].clone()).expect("the element was just wrapped");
        let head: Vec<TopoShape> = le[..idx].to_vec();
        let rest: Vec<TopoShape> = le[idx + 1..].to_vec();

        // `cxx:176-229`.
        if matches!(
            v.0.orientation(),
            Orientation::Forward | Orientation::Reversed
        ) {
            let mut edg = edg;
            edg.0.set_orientation(Orientation::Forward);
            // `cxx:179-183`.
            let (cur_v1, cur_v2) = edge_vertices(&edg);
            let a_par1 = cur_v1
                .as_ref()
                .map(|v1| BRepTool::parameter_of_vertex(&edg.0, &v1.0))
                .unwrap_or(0.0);
            let a_par2 = cur_v2
                .as_ref()
                .map(|v2| BRepTool::parameter_of_vertex(&edg.0, &v2.0))
                .unwrap_or(0.0);

            // `cxx:185-193`: two `EmptyCopied` edges, both `FORWARD`.
            let mut e1 = empty_copied_edge(&edg);
            e1.0.set_orientation(Orientation::Forward);
            let mut e2 = empty_copied_edge(&edg);
            e2.0.set_orientation(Orientation::Forward);

            // `cxx:194-196`.
            let new_vtx = v.clone();
            let a_tol_split_v = BRepTool::vertex_tolerance(v);

            // `cxx:199-203`: a degenerated edge grows the start vertex's
            // tolerance to the splitting vertex's; `aTolV2` reuses it
            // (`cxx:217-219`).
            let a_tol_v1 = cur_v1
                .as_ref()
                .map(|v1| {
                    if BRepTool::is_degenerated(&edg) {
                        BRepTool::vertex_tolerance(v1).max(a_tol_split_v)
                    } else {
                        BRepTool::vertex_tolerance(v1)
                    }
                })
                .unwrap_or(0.0);

            let mut a_new_list = head;

            if let Some(mut cur_v1) = cur_v1 {
                cur_v1.0.set_orientation(Orientation::Forward);
                // `cxx:198`.
                b.add(&mut e1.0, &cur_v1.0);
                // `cxx:204`.
                let _ = update_vertex_on_edge(&cur_v1.0, a_par1, &e1.0, a_tol_v1);

                // `cxx:205-208`.
                let mut new_vtx_rev = new_vtx.clone();
                new_vtx_rev.0.set_orientation(Orientation::Reversed);
                b.add(&mut e1.0, &new_vtx_rev.0);
                let _ = update_vertex_on_edge(
                    &new_vtx_rev.0,
                    p,
                    &e1.0,
                    BRepTool::vertex_tolerance(v),
                );

                // `cxx:209-212`.
                let mut new_vtx_fwd = new_vtx.clone();
                new_vtx_fwd.0.set_orientation(Orientation::Forward);
                b.add(&mut e2.0, &new_vtx_fwd.0);
                let _ = update_vertex_on_edge(
                    &new_vtx_fwd.0,
                    p,
                    &e2.0,
                    BRepTool::vertex_tolerance(v),
                );
            }

            // `cxx:214-219`.
            if let Some(mut cur_v2) = cur_v2 {
                cur_v2.0.set_orientation(Orientation::Reversed);
                b.add(&mut e2.0, &cur_v2.0);
                let a_tol_v2 = if BRepTool::is_degenerated(&edg) {
                    a_tol_v1
                } else {
                    BRepTool::vertex_tolerance(&cur_v2)
                };
                let _ = update_vertex_on_edge(&cur_v2.0, a_par2, &e2.0, a_tol_v2);
            }

            // `cxx:221-228`: the two halves, then the untouched remainder.
            a_new_list.push(e1.0.clone());
            a_new_list.push(e2.0.clone());
            a_new_list.extend(rest);

            // `cxx:229-230`: `myMap.UnBind(E); myMap.Bind(E, aNewList);`.
            self.map_remove(&e.0);
            let _ = self.map_bind(&e.0, a_new_list);
        } else {
            // `cxx:231-248`.
            let mut e1 = empty_copied_edge(&edg);
            for child in
                crate::topo_tools::explorer_shapes_avoiding(&edg.0, ShapeType::Vertex, ShapeType::Edge)
            {
                let Some(vtx) = Vertex::wrap(child) else {
                    continue;
                };
                let f = BRepTool::parameter_of_vertex(&edg.0, &vtx.0);
                b.add(&mut e1.0, &vtx.0);
                let _ = update_vertex_on_edge(&vtx.0, f, &e1.0, BRepTool::vertex_tolerance(&vtx));
            }
            b.add(&mut e1.0, &v.0);
            let _ = update_vertex_on_edge(&v.0, p, &e1.0, BRepTool::vertex_tolerance(v));

            // `cxx:248`: `le.Append(E1)`.
            let mut le_rest = head;
            le_rest.extend(rest);
            le_rest.push(e1.0.clone());
            self.set_map_of(&e.0, le_rest);
        }
    }
}
