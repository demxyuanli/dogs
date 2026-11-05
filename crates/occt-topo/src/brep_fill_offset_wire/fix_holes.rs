//! `BRepFill_OffsetWire::FixHoles` (`cxx:1583-1920`): removes duplicated edges
//! from the offset wires and stitches the dangling ends of the open ones, so a
//! closed offset comes out without holes.

use std::collections::HashMap;


use occt_core::precision::CONFUSION;

use crate::abs::ShapeType;
use crate::brep_fill_offset_wire::vertices::update_vertex_tolerance;
use crate::brep_fill_offset_wire::BRepFillOffsetWire;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, TopoShape, Vertex, Wire};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{
    edges_of_wire, is_same, map_shapes_and_ancestors, shape_map_key, vertices_of, wire_vertices,
    wires_of, ShapeMapKey,
};

impl BRepFillOffsetWire {
    /// `FixHoles()` (`cxx:1583-1920`).
    pub(crate) fn fix_holes(&mut self) {
        let builder = TopoBuilder::new();

        // `cxx:1589-1600`: the tolerance of the spine's vertices, x100.
        let mut max_tol = 0.0f64;
        if let Some(spine) = &self.my_spine {
            for v in vertices_of(&spine.0) {
                let tol = BRepTool::vertex_tolerance(&v);
                if tol > max_tol {
                    max_tol = tol;
                }
            }
        }
        max_tol *= 100.0;

        let mut closed_wires: Vec<TopoShape> = Vec::new();
        let mut unclosed_wires: Vec<TopoShape> = Vec::new();
        let mut isolated_wires: Vec<TopoShape> = Vec::new();

        // `cxx:1602-1644`.
        let shape = self.my_shape.clone().expect("FixHoles: null myShape");
        for w in wires_of(&shape) {
            let mut ee_map: HashMap<ShapeMapKey, Vec<TopoShape>> = HashMap::new();
            for e in edges_of_wire(&w) {
                // `cxx:1612-1621`: the first occurrence binds the key, later ones
                // are the duplicates to remove.
                ee_map.entry(shape_map_key(&e.0)).or_default().push(e.0.clone());
            }
            let mut a_wire = w.0.clone();
            a_wire.set_free(true);
            for list in ee_map.values() {
                for e in list.iter().skip(1) {
                    builder.remove(&mut a_wire, e);
                }
            }
            if a_wire.is_closed() {
                closed_wires.push(a_wire);
            } else {
                unclosed_wires.push(a_wire);
            }
        }

        // `cxx:1646-1892`.
        while !unclosed_wires.is_empty() {
            let mut base = unclosed_wires[0].clone();
            let (vf, vl) = wire_vertices(&Wire(base.clone()));
            let (Some(vf), Some(vl)) = (vf, vl) else {
                panic!("BRepFill_OffsetWire::FixHoles(): Wrong wire.");
            };
            let pf = GeometryRegistry::global().vertex_point(&vf.0);
            let pl = GeometryRegistry::global().vertex_point(&vl.0);
            let mut dist_f = f64::MAX;
            let mut dist_l = f64::MAX;
            let mut index_f = 0usize;
            let mut index_l = 0usize;
            let mut is_first_f = false;
            let mut is_first_l = false;

            for i in 1..unclosed_wires.len() {
                let (v1, v2) = wire_vertices(&Wire(unclosed_wires[i].clone()));
                let (Some(v1), Some(v2)) = (v1, v2) else {
                    panic!("BRepFill_OffsetWire::FixHoles(): Wrong wire.");
                };
                let p1 = GeometryRegistry::global().vertex_point(&v1.0);
                let p2 = GeometryRegistry::global().vertex_point(&v2.0);
                let mut dist = pf.distance(&p1);
                if dist < dist_f {
                    dist_f = dist;
                    index_f = i + 1;
                    is_first_f = true;
                }
                dist = pf.distance(&p2);
                if dist < dist_f {
                    dist_f = dist;
                    index_f = i + 1;
                    is_first_f = false;
                }
                dist = pl.distance(&p1);
                if dist < dist_l {
                    dist_l = dist;
                    index_l = i + 1;
                    is_first_l = true;
                }
                dist = pl.distance(&p2);
                if dist < dist_l {
                    dist_l = dist;
                    index_l = i + 1;
                    is_first_l = false;
                }
            }

            if dist_f > max_tol {
                index_f = 0;
            }
            if dist_l > max_tol {
                index_l = 0;
            }

            let mut try_to_close = true;
            if dist_f <= max_tol && dist_l <= max_tol && index_f == index_l && is_first_f == is_first_l {
                if dist_f < dist_l {
                    dist_l = f64::MAX;
                    index_l += 1;
                } else {
                    dist_f = f64::MAX;
                    index_f += 1;
                }
                try_to_close = false;
            }

            // `cxx:1731-1784`.
            if dist_f <= max_tol {
                let the_wire = Wire(unclosed_wires[index_f - 1].clone());
                let (v1, v2) = wire_vertices(&the_wire);
                let vemap =
                    map_shapes_and_ancestors(&the_wire.0, ShapeType::Vertex, ShapeType::Edge);
                let anchor = if is_first_f { v1.clone() } else { v2.clone() };
                let the_edge = vemap
                    .get(&shape_map_key(
                        &anchor
                            .expect("FixHoles: wire without vertices")
                            .0,
                    ))
                    .and_then(|l| l.first())
                    .map(|e| Edge(e.clone()))
                    .expect("FixHoles: vertex without an incident edge");
                for item in edges_of_wire(&the_wire) {
                    let mut an_edge = item.0.clone();
                    if is_first_f {
                        an_edge.reverse();
                    }
                    if !is_same(&an_edge, &the_edge.0) {
                        builder.add(&mut base, &an_edge);
                    }
                }
                let the_vertex = if is_first_f {
                    v1.expect("FixHoles: wire without vertices")
                } else {
                    v2.expect("FixHoles: wire without vertices")
                };
                let common_tol = BRepTool::vertex_tolerance(&vf)
                    .max(BRepTool::vertex_tolerance(&the_vertex));
                if dist_f <= common_tol {
                    let_edge_stitch(
                        &builder,
                        &mut base,
                        the_edge,
                        the_vertex,
                        &vf,
                        is_first_f,
                        common_tol,
                    );
                } else {
                    let mut the_edge = the_edge.0.clone();
                    if is_first_f {
                        the_edge.reverse();
                    }
                    builder.add(&mut base, &the_edge);
                    let new_edge = make_edge_between(&builder, &the_vertex, &vf);
                    builder.add(&mut base, &new_edge.0);
                }
            }

            // `cxx:1785-1836`.
            if dist_l <= max_tol && index_l != index_f {
                let the_wire = Wire(unclosed_wires[index_l - 1].clone());
                let (v1, v2) = wire_vertices(&the_wire);
                let vemap =
                    map_shapes_and_ancestors(&the_wire.0, ShapeType::Vertex, ShapeType::Edge);
                let anchor = if is_first_l { v1.clone() } else { v2.clone() };
                let the_edge = vemap
                    .get(&shape_map_key(
                        &anchor
                            .expect("FixHoles: wire without vertices")
                            .0,
                    ))
                    .and_then(|l| l.first())
                    .map(|e| Edge(e.clone()))
                    .expect("FixHoles: vertex without an incident edge");
                for item in edges_of_wire(&the_wire) {
                    let mut an_edge = item.0.clone();
                    if !is_first_l {
                        an_edge.reverse();
                    }
                    if !is_same(&an_edge, &the_edge.0) {
                        builder.add(&mut base, &an_edge);
                    }
                }
                let the_vertex = if is_first_l {
                    v1.expect("FixHoles: wire without vertices")
                } else {
                    v2.expect("FixHoles: wire without vertices")
                };
                let common_tol = BRepTool::vertex_tolerance(&vl)
                    .max(BRepTool::vertex_tolerance(&the_vertex));
                if dist_l <= common_tol {
                    let_edge_stitch(
                        &builder,
                        &mut base,
                        the_edge,
                        the_vertex,
                        &vl,
                        !is_first_l,
                        common_tol,
                    );
                } else {
                    let mut the_edge = the_edge.0.clone();
                    if !is_first_l {
                        the_edge.reverse();
                    }
                    builder.add(&mut base, &the_edge);
                    let new_edge = make_edge_between(&builder, &vl, &the_vertex);
                    builder.add(&mut base, &new_edge.0);
                }
            }

            // `cxx:1837-1866`: try to close the resulting wire.
            if try_to_close {
                let (vf2, vl2) = wire_vertices(&Wire(base.clone()));
                let (Some(vf2), Some(vl2)) = (vf2, vl2) else {
                    panic!("BRepFill_OffsetWire::FixHoles(): Wrong wire.");
                };
                let common_tol = BRepTool::vertex_tolerance(&vf2)
                    .max(BRepTool::vertex_tolerance(&vl2));
                let vemap = map_shapes_and_ancestors(&base, ShapeType::Vertex, ShapeType::Edge);
                let e_last = vemap
                    .get(&shape_map_key(&vl2.0))
                    .and_then(|l| l.first())
                    .map(|e| e.clone())
                    .expect("FixHoles: vertex without an incident edge");
                let p_f = GeometryRegistry::global().vertex_point(&vf2.0);
                let p_l = GeometryRegistry::global().vertex_point(&vl2.0);
                let dist = p_f.distance(&p_l);
                if dist <= common_tol {
                    let mut e_last = e_last;
                    e_last.set_free(true);
                    let mut vf_new = vf2.0.clone();
                    vf_new.set_orientation(vl2.0.orientation());
                    builder.remove(&mut e_last, &vl2.0);
                    builder.add(&mut e_last, &vf_new);
                    update_vertex_tolerance(&Vertex(vf_new.clone()), common_tol);
                    base.set_closed(true);
                } else if dist <= max_tol {
                    let new_edge = make_edge_between(&builder, &vf2, &vl2);
                    builder.add(&mut base, &new_edge.0);
                    base.set_closed(true);
                }
            }

            // `cxx:1868-1891`.
            if dist_f <= max_tol {
                unclosed_wires.remove(index_f - 1);
            }
            if dist_l <= max_tol && index_l != index_f {
                if dist_f <= max_tol && index_l > index_f {
                    index_l -= 1;
                }
                unclosed_wires.remove(index_l - 1);
            }
            if base.is_closed() {
                closed_wires.push(base);
                unclosed_wires.remove(0);
            } else if dist_f > max_tol && dist_l > max_tol {
                isolated_wires.push(base);
                unclosed_wires.remove(0);
            }
        }

        // `cxx:1894-1919`.
        if closed_wires.len() + isolated_wires.len() == 1 {
            self.my_shape = Some(if !closed_wires.is_empty() {
                closed_wires.remove(0)
            } else {
                isolated_wires.remove(0)
            });
        } else {
            let mut all = closed_wires;
            all.extend(isolated_wires);
            let compound = builder.make_compound_of(&all);
            self.my_shape = Some(compound.0);
        }
    }
}

/// `cxx:1758-1768` / `cxx:1808-1818`: `Free`, `Remove` the old vertex, `Add`
/// the kept one, raise its tolerance, then restore the edge orientation.
fn let_edge_stitch(
    builder: &TopoBuilder,
    base: &mut TopoShape,
    mut the_edge: Edge,
    the_vertex: Vertex,
    kept: &Vertex,
    reverse_after: bool,
    common_tol: f64,
) {
    the_edge.0.set_free(true);
    let mut kept = kept.0.clone();
    kept.set_orientation(the_vertex.0.orientation());
    builder.remove(&mut the_edge.0, &the_vertex.0);
    builder.add(&mut the_edge.0, &kept);
    update_vertex_tolerance(&Vertex(kept.clone()), common_tol);
    if reverse_after {
        the_edge.0.reverse();
    }
    builder.add(base, &the_edge.0);
}

/// `BRepLib_MakeEdge(V1, V2)` between two existing vertices. OCCT raises
/// `BRepLib_LineThroughIdenticPoints` for coincident points; the port panics.
fn make_edge_between(builder: &TopoBuilder, v1: &Vertex, v2: &Vertex) -> Edge {
    let reg = GeometryRegistry::global();
    let p1 = reg.vertex_point(&v1.0);
    let p2 = reg.vertex_point(&v2.0);
    if p1.distance(&p2) <= CONFUSION {
        panic!("BRepLib_MakeEdge: BRepLib_LineThroughIdenticPoints");
    }
    builder.make_edge_segment_with_vertices(&p1, &p2, v1, v2)
}
