//! `BRepFill_Pipe::BuildHistory` (`BRepFill_Pipe.cxx:1100-1220`).

use crate::abs::{Orientation, ShapeType};
use crate::brep_fill_pipe::make_shape::update_map;
use crate::brep_fill_pipe::Pipe;
use crate::brep_fill_sweep::BrepFillSweep;
use crate::brep_tools_wire_explorer::WireExplorer;
use crate::brep_tool::BRepTool;
use crate::iterator::ShapeIterator;
use crate::shape::{Shell, TopoShape, Wire};
use crate::topo_tools::{edge_vertices, shape_map_key};

impl Pipe {
    /// `BuildHistory(theSweep, theSection)` (`cxx:1100-1220`): record the
    /// generated shapes of the section edges, of the section vertices and of the
    /// spine edges / vertices.
    ///
    /// Null values are not appended: OCCT's `UpdateMap` stores them, but a null
    /// shape has no use in `Generated`.
    pub(crate) fn build_history(
        &mut self,
        the_sweep: &BrepFillSweep,
        the_section: &TopoShape,
    ) -> Result<(), String> {
        // `cxx:1102`: the U isos of the sweep.
        let an_u_edges = the_sweep.inter_faces().clone();

        // `cxx:1105-1145`: the wire section.
        if the_section.shape_type() == ShapeType::Wire {
            let Some(a_section) = Wire::wrap(the_section.clone()) else {
                return Err("Standard_TypeMismatch: BRepFill_Pipe::BuildHistory".to_string());
            };
            let mut wexp = WireExplorer::of_wire(&a_section);
            let mut inde = 0usize;
            // The reference's `for` runs `wexp.Next()` in its increment, so
            // every `continue` below still advances the explorer.
            loop {
                if !wexp.more() {
                    break;
                }
                inde += 1;
                let current = wexp.current_edge();
                wexp.next();
                let Some(an_edge) = current else {
                    continue;
                };
                if BRepTool::is_degenerated(&an_edge) {
                    continue;
                }
                if self.gen_map().contains_key(&shape_map_key(&an_edge.0)) {
                    continue;
                }
                // `cxx:1118-1124`.
                let (v0, v1) = edge_vertices(&an_edge);
                let a_tape = the_sweep.tape(inde).cloned().and_then(Shell::wrap);
                let mut u_index = [inde, inde + 1];
                let vertices = [v0, v1];
                if an_edge.0.orientation() == Orientation::Reversed {
                    u_index.swap(0, 1);
                }
                // `cxx:1126-1139`.
                for kk in 0..2 {
                    let Some(vertex) = vertices[kk].as_ref() else {
                        continue;
                    };
                    if self.gen_map().contains_key(&shape_map_key(&vertex.0)) {
                        continue;
                    }
                    let list = self.gen_map_mut().entry(shape_map_key(&vertex.0)).or_default();
                    for jj in 1..=an_u_edges.nb_cols() {
                        if let Some(u_edge) = an_u_edges.value(u_index[kk], jj) {
                            list.push(u_edge.clone());
                        }
                    }
                }
                // `cxx:1141-1145`.
                if let Some(tape) = a_tape {
                    let list = self.gen_map_mut().entry(shape_map_key(&an_edge.0)).or_default();
                    for child in ShapeIterator::of_shape(&tape.0) {
                        list.push(child);
                    }
                }
            }
        }

        // `cxx:1147-1220`: the spine edges and vertices.
        let a_faces = the_sweep.sub_shape().clone();
        let a_v_edges = the_sweep.sections().clone();
        let Some(spine) = self.spine().cloned() else {
            return Ok(());
        };
        let mut wexp = WireExplorer::of_wire(&spine);
        let mut inde = 0usize;
        loop {
            let to_exit = !wexp.more();
            inde += 1;
            if !to_exit {
                if let Some(an_edge_of_spine) = wexp.current_edge() {
                    for i in 1..=a_faces.nb_rows() {
                        if let Some(a_face) = a_faces.value(i, inde) {
                            update_map(self.gen_map_mut(), &an_edge_of_spine.0, a_face);
                        }
                    }
                }
            }
            // `cxx:1200-1214`: the vertex of the current step, which is the last
            // vertex of the spine on the exit iteration.
            if let Some(a_vertex_of_spine) = wexp.current_vertex().cloned() {
                for i in 1..=a_v_edges.nb_rows() {
                    if let Some(a_v_edge) = a_v_edges.value(i, inde) {
                        update_map(self.gen_map_mut(), &a_vertex_of_spine.0, a_v_edge);
                    }
                }
            }
            if to_exit {
                break;
            }
            if wexp.more() {
                wexp.next();
            }
        }
        Ok(())
    }
}
