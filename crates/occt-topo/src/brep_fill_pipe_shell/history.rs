//! `BRepFill_PipeShell::BuildHistory` (`BRepFill_PipeShell.cxx:1282-1615`).

use std::collections::HashMap;

use crate::abs::{Orientation, ShapeType};
use crate::brep_fill_pipe_shell::PipeShell;
use crate::brep_fill_sweep::BrepFillSweep;
use crate::brep_tool::BRepTool;
use crate::brep_tools_wire_explorer::WireExplorer;
use crate::iterator::ShapeIterator;
use crate::shape::{Edge, Shell, TopoShape, Wire};
use crate::builder::TopoBuilder;
use crate::topo_tools::{
    edge_vertices, is_same, map_shapes_and_ancestors, shape_map_key, ShapeMapKey,
};

impl PipeShell {
    /// `BuildHistory(theSweep)` (`cxx:1282-1615`): fill `myGenMap` with the
    /// images of the section edges / vertices and of the spine edges / vertices.
    pub fn build_history(&mut self, the_sweep: &BrepFillSweep) -> Result<(), String> {
        let builder = TopoBuilder::new();
        let an_u_edges = the_sweep.inter_faces();
        let nb_rows = an_u_edges.nb_rows();
        let nb_cols = an_u_edges.nb_cols();
        let mut ind_wire_map: HashMap<i32, TopoShape> = HashMap::new();
        let mut inde = 1usize;

        for indw in 1..=self.my_seq.len() {
            let Some(section_shape) = self.my_seq[indw - 1].original_shape().cloned() else {
                continue;
            };
            if self.my_seq[indw - 1].is_punctual() {
                // `cxx:1298-1312`: a punctual end section generates every wire
                // found along the path.
                let mut generated: Vec<TopoShape> = Vec::new();
                for i in 1..=nb_rows {
                    for j in 1..=nb_cols {
                        if let Some(shape) = an_u_edges.value(i, j) {
                            generated.push(shape.clone());
                        }
                    }
                }
                self.my_gen_map
                    .entry(shape_map_key(&section_shape))
                    .or_default()
                    .extend(generated);
                continue;
            }

            let a_section = Wire::wrap(section_shape).ok_or_else(|| {
                "BRepFill_PipeShell::BuildHistory: the section is not a wire".to_string()
            })?;
            // `cxx:1318-1322`: the real index of the section on the path.
            let ind_of_w = *self.my_ind_of_sec.get(indw - 1).ok_or_else(|| {
                "BRepFill_PipeShell::BuildHistory: no section index".to_string()
            })?;
            let the_wire = self
                .w_seq
                .get(ind_of_w - 1)
                .ok_or_else(|| "BRepFill_PipeShell::BuildHistory: no working wire".to_string())?
                .clone();

            let mut wexp_sec = WireExplorer::of_wire(&a_section);
            while wexp_sec.more() {
                let an_original_edge = wexp_sec.current().cloned();
                wexp_sec.next();
                let Some(an_original_edge) = an_original_edge else {
                    inde += 1;
                    continue;
                };

                let Some(an_edge) = self.my_seq[indw - 1]
                    .modified_shape(&an_original_edge.0)
                    .and_then(Edge::wrap)
                else {
                    inde += 1;
                    continue;
                };
                if BRepTool::is_degenerated(&an_edge) {
                    inde += 1;
                    continue;
                }

                let mut a_shell = Shell::new();
                let (a_vertex0, a_vertex1) = edge_vertices(&an_original_edge);
                let a_vertex = [a_vertex0.map(|v| v.0), a_vertex1.map(|v| v.0)];
                let sign_of_an_edge = if an_original_edge.orientation() == Orientation::Forward {
                    1i32
                } else {
                    -1
                };

                // `cxx:1339-1353`: the `inde`-th non-degenerated edge of the
                // working wire.
                let mut the_edge: Option<Edge> = None;
                let mut i = 1usize;
                let mut wexp = WireExplorer::of_wire(&the_wire);
                while wexp.more() {
                    if let Some(current) = wexp.current_edge() {
                        the_edge = Some(current);
                    }
                    if !BRepTool::is_degenerated(&an_edge) {
                        if i == inde {
                            break;
                        }
                        i += 1;
                    }
                    wexp.next();
                }

                // `cxx:1357-1377`: every split of `theEdge` contributes its
                // "tape" (a shell running from one end of the path to the
                // other) to `aShell`.
                let new_edges: Vec<TopoShape> = the_edge
                    .as_ref()
                    .and_then(|edge| self.my_edge_new_edges.get(&shape_map_key(&edge.0)))
                    .cloned()
                    .unwrap_or_default();
                let section_arc = self
                    .section()
                    .cloned()
                    .ok_or_else(|| "BRepFill_PipeShell::BuildHistory: no section law".to_string())?;
                let mut sign_of_a_new_edge = 0i32;
                let mut sign_of_index = 0i32;
                for a_new_edge in &new_edges {
                    sign_of_a_new_edge = if a_new_edge.orientation() == Orientation::Forward {
                        1
                    } else {
                        -1
                    };
                    // `IndexOfEdge` throws for an unknown edge. The index is
                    // signed: `-ii` for a `REVERSED` section edge
                    // (`BRepFill_NSections.cxx:489-495`).
                    let an_ind_e = section_arc.base().index_of_edge(a_new_edge)?;
                    sign_of_index = if an_ind_e > 0 { 1 } else { -1 };
                    if let Some(a_tape) = the_sweep.tape(an_ind_e.unsigned_abs() as usize) {
                        for shape in ShapeIterator::of_shape(a_tape) {
                            add_if_not_null(&builder, &mut a_shell.0, &shape);
                        }
                    }
                }

                // `cxx:1380-1406`: the two u-edges the section vertices lie on.
                let to_reverse = sign_of_an_edge * sign_of_a_new_edge * sign_of_index;
                let nb_law = section_arc.base().nb_law() as i32;
                let first_index = new_edges
                    .first()
                    .map(|shape| section_arc.base().index_of_edge(shape))
                    .transpose()?
                    .map(|index| index.abs())
                    .unwrap_or(0);
                let last_index = new_edges
                    .last()
                    .map(|shape| section_arc.base().index_of_edge(shape))
                    .transpose()?
                    .map(|index| index.abs())
                    .unwrap_or(0);
                let mut u_index = [first_index, last_index + to_reverse];
                if to_reverse == -1 {
                    u_index[0] += 1;
                    u_index[1] += 1;
                }
                if section_arc.base().is_u_closed() {
                    if u_index[0] > nb_law {
                        u_index[0] = 1;
                    }
                    if u_index[1] > nb_law {
                        u_index[1] = 1;
                    }
                }
                if sign_of_an_edge == -1 || sign_of_a_new_edge == -1 {
                    u_index.swap(0, 1);
                }

                let vemap =
                    map_shapes_and_ancestors(&a_shell.0, ShapeType::Vertex, ShapeType::Edge);

                for kk in 0..2 {
                    let Some(vertex) = a_vertex[kk].as_ref() else {
                        continue;
                    };
                    let vertex_key = shape_map_key(vertex);
                    if self.my_gen_map.contains_key(&vertex_key) {
                        continue;
                    }
                    // `cxx:1419-1431`: the rail was already built for this
                    // u-edge.
                    if let Some(a_wire) = ind_wire_map.get(&u_index[kk]).cloned() {
                        let generated: Vec<TopoShape> = ShapeIterator::of_shape(&a_wire)
                            .collect();
                        self.my_gen_map
                            .entry(vertex_key)
                            .or_default()
                            .extend(generated);
                        continue;
                    }

                    // `cxx:1434-1442`: the u-edges of this row, first column to
                    // last.
                    let mut seq_edges: Vec<Option<TopoShape>> = Vec::new();
                    for jj in 1..=nb_cols {
                        seq_edges.push(an_u_edges.value(u_index[kk] as usize, jj).cloned());
                    }

                    // `cxx:1444-1508`: assemble the rail along the path,
                    // filling the holes with the shell edges.
                    let location_closed = self
                        .location()
                        .map(|loc| loc.base().is_closed())
                        .unwrap_or(false);
                    let Some(a_wire) =
                        assemble_rail(&builder, &seq_edges, &vemap, location_closed)
                    else {
                        continue;
                    };

                    let generated: Vec<TopoShape> = ShapeIterator::of_shape(&a_wire.0).collect();
                    self.my_gen_map
                        .entry(vertex_key)
                        .or_default()
                        .extend(generated);
                    ind_wire_map.insert(u_index[kk], a_wire.0.clone());
                }

                // `cxx:1524-1530`: the original section edge generates the
                // shell's faces.
                let generated: Vec<TopoShape> = ShapeIterator::of_shape(&a_shell.0).collect();
                self.my_gen_map
                    .entry(shape_map_key(&an_original_edge.0))
                    .or_default()
                    .extend(generated);

                inde += 1;
            }
        }

        // `cxx:1536-1613`: the spine sub-shapes.
        let a_faces = the_sweep.sub_shape();
        let a_v_edges = the_sweep.sections();
        let nb_faces = a_faces.nb_rows();
        let nb_v_edges = a_v_edges.nb_rows();

        let mut wexp = WireExplorer::of_wire(&self.my_spine.clone());
        let mut inde = 0usize;
        let mut to_exit = false;
        loop {
            if !wexp.more() {
                to_exit = true;
            }
            inde += 1;

            if !to_exit {
                if let Some(an_edge_of_spine) = wexp.current().cloned() {
                    let mut generated: Vec<TopoShape> = Vec::new();
                    for i in 1..=nb_faces {
                        if let Some(face) = a_faces.value(i, inde) {
                            if face.shape_type() == ShapeType::Face {
                                generated.push(face.clone());
                            }
                        }
                    }
                    self.my_gen_map
                        .entry(shape_map_key(&an_edge_of_spine.0))
                        .or_default()
                        .extend(generated);
                }
            }

            if let Some(a_vertex_of_spine) = wexp.current_vertex().map(|v| v.0.clone()) {
                let key = shape_map_key(&a_vertex_of_spine);
                let mut generated: Vec<TopoShape> = Vec::new();
                for i in 1..=nb_v_edges {
                    let Some(vshape) = a_v_edges.value(i, inde) else {
                        continue;
                    };
                    if vshape.shape_type() == ShapeType::Edge
                        || vshape.shape_type() == ShapeType::Face
                    {
                        generated.push(vshape.clone());
                    } else {
                        for sub in ShapeIterator::of_shape(vshape) {
                            if sub.shape_type() == ShapeType::Edge
                                || sub.shape_type() == ShapeType::Face
                            {
                                generated.push(sub);
                            } else {
                                // it is a wire
                                for child in ShapeIterator::of_shape(&sub) {
                                    generated.push(child);
                                }
                            }
                        }
                    }
                }
                self.my_gen_map.entry(key).or_default().extend(generated);
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

/// `cxx:1444-1508`: the rail wire of one section vertex, holes filled from the
/// shell. `None` is the `continue` of `cxx:1449-1451`.
fn assemble_rail(
    builder: &TopoBuilder,
    seq_edges: &[Option<TopoShape>],
    vemap: &HashMap<ShapeMapKey, Vec<TopoShape>>,
    location_closed: bool,
) -> Option<Wire> {
    let first_edge = seq_edges.first()?.clone()?;
    let mut a_wire = Wire::new();
    add_if_not_null(builder, &mut a_wire.0, &first_edge);

    let (first_vertex, mut cur_vertex) = edge_ends(Some(&first_edge));
    let mut cur_edge: Option<TopoShape> = None;

    for slot in seq_edges.iter().skip(1) {
        cur_edge = slot.clone();
        let (vfirst, vlast) = edge_ends(cur_edge.as_ref());

        if same_end(&cur_vertex, &vfirst) {
            // `CurVertex = Vlast` of `cxx:1462` is overwritten at the end of the
            // iteration, so only the hole case reads `CurVertex` here.
        } else {
            // A hole: take the shell edge that joins the two vertices.
            for candidate in candidates_at(vemap, vfirst.as_ref()) {
                if let Some(current) = cur_edge.as_ref() {
                    if is_same(&candidate, current) {
                        continue;
                    }
                }
                let (v1, v2) = edge_ends(Some(&candidate));
                if same_end(&v1, &cur_vertex) || same_end(&v2, &cur_vertex) {
                    add_if_not_null(builder, &mut a_wire.0, &candidate);
                    break;
                }
            }
        }

        cur_vertex = vlast.clone();
        if let Some(current) = cur_edge.as_ref() {
            add_if_not_null(builder, &mut a_wire.0, current);
        }
    }

    // Case of a closed wire (`cxx:1489-1508`).
    if location_closed && !same_end(&cur_vertex, &first_vertex) {
        for candidate in candidates_at(vemap, cur_vertex.as_ref()) {
            if let Some(current) = cur_edge.as_ref() {
                if is_same(&candidate, current) {
                    continue;
                }
            }
            let (v1, v2) = edge_ends(Some(&candidate));
            if same_end(&v1, &first_vertex) || same_end(&v2, &first_vertex) {
                add_if_not_null(builder, &mut a_wire.0, &candidate);
                break;
            }
        }
    }

    Some(a_wire)
}

/// `VEmap.FindFromKey(V)`; an unknown or null vertex has no ancestor.
fn candidates_at(
    vemap: &HashMap<ShapeMapKey, Vec<TopoShape>>,
    vertex: Option<&TopoShape>,
) -> Vec<TopoShape> {
    match vertex {
        Some(vertex) => vemap.get(&shape_map_key(vertex)).cloned().unwrap_or_default(),
        None => Vec::new(),
    }
}

/// `TopExp::Vertices(E, V1, V2)` (`TopExp.cxx:182-210`) with `CumOri` false; a
/// null edge has two null vertices.
fn edge_ends(edge: Option<&TopoShape>) -> (Option<TopoShape>, Option<TopoShape>) {
    match edge.and_then(|shape| Edge::wrap(shape.clone())) {
        Some(edge) => {
            let (v1, v2) = edge_vertices(&edge);
            (v1.map(|v| v.0), v2.map(|v| v.0))
        }
        None => (None, None),
    }
}

/// `TopoDS_Shape::IsSame` where a null shape is the same as another null shape
/// (both hold a null `TShape` handle in OCCT).
fn same_end(a: &Option<TopoShape>, b: &Option<TopoShape>) -> bool {
    match (a, b) {
        (None, None) => true,
        (Some(x), Some(y)) => is_same(x, y),
        _ => false,
    }
}

/// `BRep_Builder::Add` skips a null component.
fn add_if_not_null(builder: &TopoBuilder, parent: &mut TopoShape, sub: &TopoShape) {
    if sub.is_null() {
        return;
    }
    builder.add(parent, sub);
}
