//! `BRepSweep_Trsf::SetContinuity` (`BRepSweep_Trsf.cxx:97-201`).

use occt_core::precision::CONFUSION;
use occt_geom::curve::Curve;

use crate::abs::{Orientation, ShapeType};
use crate::brep_lib_continuity::set_stored_continuity;
use crate::brep_lprop;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, TopoShape};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools_full::edges_of;

use super::grid::{shape_of, SweepOps};
use super::num_shape::NumShape;

const TA: f64 = 0.00175;

pub fn set_continuity(host: &dyn SweepOps, gen: &TopoShape, dir: NumShape) {
    if gen.shape_type() == ShapeType::Edge {
        set_edge(host, gen, dir);
    } else if gen.shape_type() == ShapeType::Wire {
        set_wire(host, gen, dir);
    }
}

fn set_edge(host: &dyn SweepOps, gen: &TopoShape, dir: NumShape) {
    if !host.has_shape(gen, dir) {
        return;
    }
    let (d, f) = edge_ends(gen);
    if let (Some(d), Some(f)) = (d, f) {
        if same_located(&d, &f) {
            let tol3d = CONFUSION.max(2.0 * BRepTool::vertex_tolerance(&crate::shape::Vertex(d.clone())));
            let (ud, uf) = end_parameters(gen, &d, &f);
            match pair_continuity(gen, &d, ud, gen, &f, uf, tol3d) {
                Ok(cont) if cont >= 1 => {
                    write_same_face(host, &d, dir, gen, dir, cont);
                }
                Err(()) => {
                    host.tables().fail("BRepLProp::Continuity");
                    return;
                }
                _ => {}
            }
        }
    }
    if dir.closed {
        let dirv = host.tables().dir_tool.shape(2);
        let Some(face_shape) = shape_of(host, gen, dir) else {
            return;
        };
        let Some(sub) = shape_of(host, gen, dirv) else {
            return;
        };
        if !host.gdd_is_to_add(&face_shape, &sub, gen, dir, dirv) {
            return;
        }
        let Some(curve) = BRepTool::edge_curve(&Edge(gen.clone())) else {
            host.tables().fail("BRep_Tool::Curve");
            return;
        };
        let cont = curve.continuity();
        if cont >= 1 {
            write_same_face(host, gen, dirv, gen, dir, cont);
        }
    }
}

fn set_wire(host: &dyn SweepOps, gen: &TopoShape, dir: NumShape) {
    let pairs = vertex_edges(gen);
    for (vertex, edges) in pairs {
        let n = edges.len();
        for j in 0..n {
            for k in (j + 1)..n {
                let e1 = &edges[j];
                let e2 = &edges[k];
                if same_located(e1, e2) || !host.has_shape(e1, dir) || !host.has_shape(e2, dir) {
                    continue;
                }
                let u1 = vertex_parameter(&vertex, e1);
                let u2 = vertex_parameter(&vertex, e2);
                let tol3d =
                    CONFUSION.max(2.0 * BRepTool::vertex_tolerance(&crate::shape::Vertex(vertex.clone())));
                match pair_continuity(e1, &vertex, u1, e2, &vertex, u2, tol3d) {
                    Ok(cont) if cont >= 1 => {
                        let Some(swept_v) = shape_of(host, &vertex, dir) else {
                            continue;
                        };
                        let Some(f1) = shape_of(host, e1, dir) else {
                            continue;
                        };
                        let Some(f2) = shape_of(host, e2, dir) else {
                            continue;
                        };
                        if swept_v.shape_type() != ShapeType::Edge
                            || f1.shape_type() != ShapeType::Face
                            || f2.shape_type() != ShapeType::Face
                        {
                            host.tables().fail("BRepSweep_Trsf::SetContinuity");
                            return;
                        }
                        set_stored_continuity(
                            &Edge(swept_v),
                            &Face(f1),
                            &Face(f2),
                            cont,
                        );
                    }
                    Err(()) => {
                        host.tables().fail("BRepLProp::Continuity");
                        return;
                    }
                    _ => {}
                }
            }
        }
    }
}

fn write_same_face(
    host: &dyn SweepOps,
    edge_gen: &TopoShape,
    edge_dir: NumShape,
    face_gen: &TopoShape,
    face_dir: NumShape,
    cont: u8,
) {
    let Some(edge) = shape_of(host, edge_gen, edge_dir) else {
        return;
    };
    let Some(face) = shape_of(host, face_gen, face_dir) else {
        return;
    };
    if edge.shape_type() != ShapeType::Edge || face.shape_type() != ShapeType::Face {
        host.tables().fail("BRepSweep_Trsf::SetContinuity");
        return;
    }
    set_stored_continuity(&Edge(edge), &Face(face.clone()), &Face(face), cont);
}

fn pair_continuity(
    e1: &TopoShape,
    _v1: &TopoShape,
    u1: f64,
    e2: &TopoShape,
    _v2: &TopoShape,
    u2: f64,
    tl: f64,
) -> Result<u8, ()> {
    let c1 = BRepTool::edge_curve_world(&Edge(e1.clone())).ok_or(())?;
    let c2 = BRepTool::edge_curve_world(&Edge(e2.clone())).ok_or(())?;
    let range1 = BRepTool::edge_parameters(&Edge(e1.clone()));
    let range2 = BRepTool::edge_parameters(&Edge(e2.clone()));
    brep_lprop::continuity(
        c1.as_ref(),
        e1.orientation(),
        u1,
        range1,
        c2.as_ref(),
        e2.orientation(),
        u2,
        range2,
        same_located(e1, e2),
        tl,
        TA,
    )
}

/// `BRep_Tool.cxx:1432-1445`: a closed 3D curve reports `First` or `Last` from
/// the vertex orientation. One parameter slot cannot hold both ends.
fn end_parameters(edge: &TopoShape, d: &TopoShape, f: &TopoShape) -> (f64, f64) {
    let (first, last) = BRepTool::edge_parameters(&Edge(edge.clone()));
    let ud = if d.orientation() == Orientation::Reversed {
        last
    } else {
        first
    };
    let uf = if f.orientation() == Orientation::Reversed {
        last
    } else {
        first
    };
    (ud, uf)
}

fn vertex_parameter(vertex: &TopoShape, edge: &TopoShape) -> f64 {
    let key = GeometryRegistry::shape_key(vertex);
    if let Ok(ts) = edge.tshape.read() {
        if let Some(core) = ts.edge_core() {
            if let Some((_, p)) = core.vertex_params.iter().find(|(k, _)| *k == key) {
                return *p;
            }
        }
    }
    let (first, last) = BRepTool::edge_parameters(&Edge(edge.clone()));
    if vertex.orientation() == Orientation::Reversed {
        last
    } else {
        first
    }
}

fn edge_ends(edge: &TopoShape) -> (Option<TopoShape>, Option<TopoShape>) {
    let stored = edge
        .tshape
        .read()
        .expect("poisoned TShape lock")
        .children
        .clone();
    let parent_loc = edge.location().clone();
    let move_loc = !parent_loc.is_identity();
    let mut first = None;
    let mut last = None;
    for mut child in stored {
        if child.shape_type() != ShapeType::Vertex {
            continue;
        }
        if move_loc {
            child.move_location(&parent_loc);
        }
        if child.orientation() == Orientation::Forward {
            first = Some(child);
        } else if child.orientation() == Orientation::Reversed {
            last = Some(child);
        }
    }
    (first, last)
}

fn vertex_edges(wire: &TopoShape) -> Vec<(TopoShape, Vec<TopoShape>)> {
    let mut keys = Vec::new();
    let mut lists: Vec<Vec<TopoShape>> = Vec::new();
    for edge in edges_of(wire) {
        let (a, b) = edge_ends(&edge.0);
        for v in [a, b].into_iter().flatten() {
            if let Some(i) = keys.iter().position(|k: &TopoShape| same_located(k, &v)) {
                if !lists[i].iter().any(|e| same_located(e, &edge.0)) {
                    lists[i].push(edge.0.clone());
                }
            } else {
                keys.push(v);
                lists.push(vec![edge.0.clone()]);
            }
        }
    }
    keys.into_iter().zip(lists).collect()
}

fn same_located(a: &TopoShape, b: &TopoShape) -> bool {
    a.same_tshape(b) && a.location().is_equal(b.location())
}
