//! `VerticesAssembling` (`cxx:3548-3602`) with `GlueVertices`
//! (`cxx:3302-3547`), `CreateNewNodes` (`cxx:3085-3232`), `IsMergedVertices`
//! (`cxx:3234-3300`), plus the two reshaping statics `replaceNMVertices`
//! (`cxx:3608-3700`) and `ReplaceEdge` (`cxx:3707-3780`).

use occt_core::gp::GpPnt;

use crate::abs::{Orientation, ShapeType};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, TopoShape, Vertex};
use crate::shape_fix_compose_shell::MapReShape;
use crate::tgeometry::GeometryRegistry;

use super::closed_surface::is_closed_shape;
use super::shape_maps::{ShapeKeyedMap, ShapeSet, VertexInspector};
use super::Sewing;

/// `CreateNewNodes(NodeNearestNode, NodeVertices, aVertexNode, aNodeEdges)`
/// (`cxx:3085-3232`).
fn create_new_nodes(
    node_nearest_node: &ShapeKeyedMap<TopoShape>,
    node_vertices: &ShapeKeyedMap<Vec<TopoShape>>,
    vertex_node: &mut ShapeKeyedMap<TopoShape>,
    node_edges: &mut ShapeKeyedMap<Vec<TopoShape>>,
) -> bool {
    let builder = TopoBuilder::new();
    let mut old_node_new_node: ShapeKeyedMap<TopoShape> = ShapeKeyedMap::new();
    let mut new_node_old_nodes: ShapeKeyedMap<Vec<TopoShape>> = ShapeKeyedMap::new();
    let entries: Vec<(TopoShape, TopoShape)> = node_nearest_node
        .iter()
        .map(|(k, v)| (k.clone(), v.clone()))
        .collect();
    for (old_node1, old_node2) in entries {
        // The second node should also be in the map (`cxx:2810-2811`).
        if !node_nearest_node.contains(&old_node2) {
            continue;
        }
        if let Some(newnode1) = old_node_new_node.get(&old_node1) {
            if let Some(newnode2) = old_node_new_node.get(&old_node2) {
                if !crate::topo_tools::is_same(&newnode1, &newnode2) {
                    // Change the data for new node #2 (`cxx:2821-2831`).
                    let moved: Vec<TopoShape> = node_old_nodes(new_node_old_nodes.find(&newnode2));
                    {
                        let lnode1 = new_node_old_nodes.change(&newnode1);
                        for node2 in &moved {
                            lnode1.push(node2.clone());
                            old_node_new_node.bind(node2, newnode1.clone());
                        }
                    }
                    let _ = node_old_nodes;
                    new_node_old_nodes.remove(&newnode2);
                }
            } else {
                // Old node #2 is not bound: add it to old node #1.
                old_node_new_node.bind(&old_node2, newnode1.clone());
                new_node_old_nodes.change(&newnode1).push(old_node2.clone());
            }
        } else if let Some(newnode2) = old_node_new_node.get(&old_node2) {
            // Old node #1 is not bound: add it to old node #2.
            old_node_new_node.bind(&old_node1, newnode2.clone());
            new_node_old_nodes.change(&newnode2).push(old_node1.clone());
        } else {
            // The nodes are not bound: create a new node.
            let newnode = builder.make_vertex(GpPnt::zero(), 0.0);
            old_node_new_node.bind(&old_node1, newnode.0.clone());
            old_node_new_node.bind(&old_node2, newnode.0.clone());
            new_node_old_nodes.add(&newnode.0, vec![old_node1.clone(), old_node2.clone()]);
        }
    }
    if new_node_old_nodes.extent() == 0 {
        return false;
    }
    let created: Vec<(TopoShape, Vec<TopoShape>)> = new_node_old_nodes
        .iter()
        .map(|(k, v)| (k.clone(), v.clone()))
        .collect();
    for (newnode, old_nodes) in created {
        let mut x = 0.0;
        let mut y = 0.0;
        let mut z = 0.0;
        let mut lvert: Vec<TopoShape> = Vec::new();
        let mut medge = ShapeSet::new();
        let mut ledge: Vec<TopoShape> = Vec::new();
        for oldnode in &old_nodes {
            if let Some(vertices) = node_vertices.get(oldnode) {
                for v in &vertices {
                    vertex_node.bind(v, newnode.clone());
                    let p = BRepTool::vertex_point_world(&Vertex(v.clone()));
                    x += p.x();
                    y += p.y();
                    z += p.z();
                    lvert.push(v.clone());
                }
            }
            if let Some(edges) = node_edges.get(oldnode) {
                for e in &edges {
                    if medge.add(e) {
                        ledge.push(e.clone());
                    }
                }
            }
            node_edges.remove(oldnode);
        }
        node_edges.bind(&newnode, ledge);
        let center = GpPnt::new(
            x / lvert.len() as f64,
            y / lvert.len() as f64,
            z / lvert.len() as f64,
        );
        // Calculate the new node tolerance (`cxx:2912-2921`).
        let mut toler = 0.0f64;
        for v in &lvert {
            let t = center.distance(&BRepTool::vertex_point_world(&Vertex(v.clone())))
                + BRepTool::vertex_tolerance(&Vertex(v.clone()));
            if toler < t {
                toler = t;
            }
        }
        GeometryRegistry::global().set_vertex(
            &newnode,
            crate::tgeometry::VertexGeom {
                point: center,
                tolerance: toler,
            },
        );
    }
    true
}

/// `TopTools_ListOfShape` copy of a lookup, empty when unbound.
fn node_old_nodes(v: Option<&Vec<TopoShape>>) -> Vec<TopoShape> {
    v.cloned().unwrap_or_default()
}

/// `IsMergedVertices(face1, e1, e2, vtx1, vtx2)` (`cxx:3234-3300`).
fn is_merged_vertices(
    face1: Option<&TopoShape>,
    e1: &TopoShape,
    e2: &TopoShape,
    vtx1: &TopoShape,
    vtx2: &TopoShape,
) -> i32 {
    let Some(face1) = face1 else {
        // Case of floating edges.
        return if is_closed_shape(e1, vtx1, vtx2) { 0 } else { 1 };
    };
    // Find the wires containing the given edges (`cxx:2938-2951`).
    let mut wire1: Option<TopoShape> = None;
    let mut wire2: Option<TopoShape> = None;
    for wire in crate::topo_tools::explorer_shapes(face1, ShapeType::Wire) {
        if wire1.is_some() && wire2.is_some() {
            break;
        }
        for child in wire.tshape.read().expect("poisoned TShape lock").children.clone() {
            if wire1.is_none() && crate::topo_tools::is_same(e1, &child) {
                wire1 = Some(wire.clone());
            }
            if wire2.is_none() && crate::topo_tools::is_same(e2, &child) {
                wire2 = Some(wire.clone());
            }
        }
    }
    let mut status = 0;
    if let (Some(wire1), Some(wire2)) = (wire1, wire2) {
        if crate::topo_tools::is_same(&wire1, &wire2) {
            for child in wire1.tshape.read().expect("poisoned TShape lock").children.clone() {
                let (ve1, ve2) = crate::topo_tools::edge_vertices(&Edge(child.clone()));
                let (Some(ve1), Some(ve2)) = (ve1, ve2) else {
                    continue;
                };
                let a = (crate::topo_tools::is_same(&ve1.0, vtx1)
                    && crate::topo_tools::is_same(&ve2.0, vtx2))
                    || (crate::topo_tools::is_same(&ve2.0, vtx1)
                        && crate::topo_tools::is_same(&ve1.0, vtx2));
                if a {
                    return if is_closed_shape(&child, vtx1, vtx2) { 0 } else { 1 };
                }
            }
            if is_closed_shape(&wire1, vtx1, vtx2) {
                let wire = crate::shape::Wire(wire1.clone());
                let (v1, v2) = match crate::topo_tools::wire_vertices(&wire) {
                    (Some(a), Some(b)) => (a.0, b.0),
                    _ => return 0,
                };
                let is_end_vertex = (crate::topo_tools::is_same(&v1, vtx1)
                    && crate::topo_tools::is_same(&v2, vtx2))
                    || (crate::topo_tools::is_same(&v2, vtx1)
                        && crate::topo_tools::is_same(&v1, vtx2));
                if !is_end_vertex {
                    status = 1;
                }
            } else {
                status = 1;
            }
        } else {
            status = -1;
        }
    }
    status
}

/// `GlueVertices(aVertexNode, aNodeEdges, aBoundFaces, Tolerance)`
/// (`cxx:3302-3547`).
fn glue_vertices(
    vertex_node: &mut ShapeKeyedMap<TopoShape>,
    node_edges: &mut ShapeKeyedMap<Vec<TopoShape>>,
    bound_faces: &ShapeKeyedMap<Vec<TopoShape>>,
    tolerance: f64,
) -> bool {
    // Create the map node -> vertices (`cxx:2988-3011`).
    let mut node_vertices: ShapeKeyedMap<Vec<TopoShape>> = ShapeKeyedMap::new();
    let mut inspector = VertexInspector::new(tolerance);
    let entries: Vec<(TopoShape, TopoShape)> =
        vertex_node.iter().map(|(k, v)| (k.clone(), v.clone())).collect();
    for (vertex, node) in entries {
        if node_vertices.contains(&node) {
            node_vertices.change(&node).push(vertex);
        } else {
            let p = BRepTool::vertex_point_world(&Vertex(node.clone()));
            node_vertices.add(&node, vec![vertex]);
            inspector.add(p.xyz());
        }
    }
    let nb_nodes = node_vertices.extent();
    // Merge the nearest nodes (`cxx:3015-3181`).
    let mut node_nearest_node: ShapeKeyedMap<TopoShape> = ShapeKeyedMap::new();
    for i in 1..=nb_nodes {
        let Some(node1) = node_vertices.key_at(i).cloned() else {
            continue;
        };
        let pt1 = BRepTool::vertex_point_world(&Vertex(node1.clone()));
        inspector.set_current(pt1.xyz());
        let p_min = occt_core::gp::GpXyz::new(
            pt1.x() - tolerance,
            pt1.y() - tolerance,
            pt1.z() - tolerance,
        );
        let p_max = occt_core::gp::GpXyz::new(
            pt1.x() + tolerance,
            pt1.y() + tolerance,
            pt1.z() + tolerance,
        );
        inspector.clear_res_list();
        inspector.inspect_range(&p_min, &p_max);
        if inspector.res_ind().is_empty() {
            continue;
        }
        let ledges1 = node_edges.get(&node1).unwrap_or_default();
        let mut seq_nodes: Vec<TopoShape> = Vec::new();
        let mut list_nodes_same_edge: Vec<TopoShape> = Vec::new();
        let candidates: Vec<usize> = inspector.res_ind().to_vec();
        for index in candidates {
            let Some(node2) = node_vertices.key_at(index).cloned() else {
                continue;
            };
            if crate::topo_tools::is_same(&node1, &node2) {
                continue;
            }
            let ledges2 = node_edges.get(&node2).unwrap_or_default();
            // Check the merging condition for the pair of nodes
            // (`cxx:3044-3126`).
            let mut status = 0i32;
            let mut is_same_edge = false;
            'outer: for e1 in &ledges1 {
                // Obtain the real vertex from the edge (`cxx:3050-3062`).
                let mut v1 = node1.clone();
                let (ov1, ov2) = crate::topo_tools::edge_vertices(&Edge(e1.clone()));
                if let (Some(ov1), Some(ov2)) = (ov1, ov2) {
                    if let Some(n) = vertex_node.get(&ov1.0) {
                        if crate::topo_tools::is_same(&node1, &n) {
                            v1 = ov1.0.clone();
                        }
                    }
                    if let Some(n) = vertex_node.get(&ov2.0) {
                        if crate::topo_tools::is_same(&node1, &n) {
                            v1 = ov2.0.clone();
                        }
                    }
                }
                // Create the map of faces for e1 (`cxx:3066-3075`).
                let mut faces1 = ShapeSet::new();
                if let Some(lfac1) = bound_faces.get(e1) {
                    for f in &lfac1 {
                        if !crate::topo_tools::is_same(f, &TopoShape::new(ShapeType::Shape)) {
                            faces1.add(f);
                        }
                    }
                }
                for e2 in &ledges2 {
                    if status != 0 || is_same_edge {
                        break;
                    }
                    let mut v2 = node2.clone();
                    let (ov1, ov2) = crate::topo_tools::edge_vertices(&Edge(e2.clone()));
                    if let (Some(ov1), Some(ov2)) = (ov1, ov2) {
                        if let Some(n) = vertex_node.get(&ov1.0) {
                            if crate::topo_tools::is_same(&node2, &n) {
                                v2 = ov1.0.clone();
                            }
                        }
                        if let Some(n) = vertex_node.get(&ov2.0) {
                            if crate::topo_tools::is_same(&node2, &n) {
                                v2 = ov2.0.clone();
                            }
                        }
                    }
                    if let Some(lfac2) = bound_faces.get(e2) {
                        for f in &lfac2 {
                            if status != 0 || is_same_edge {
                                break;
                            }
                            if faces1.contains(f) {
                                let stat = is_merged_vertices(Some(f), e1, e2, &v1, &v2);
                                if stat == 1 {
                                    is_same_edge = true;
                                } else {
                                    status = stat;
                                }
                            }
                        }
                    } else if faces1.is_empty() && crate::topo_tools::is_same(e1, e2) {
                        let stat = is_merged_vertices(None, e1, e1, &v1, &v2);
                        if stat == 1 {
                            is_same_edge = true;
                        } else {
                            status = stat;
                        }
                        break;
                    }
                }
                if status != 0 || is_same_edge {
                    break 'outer;
                }
            }
            if status != 0 {
                continue;
            }
            if is_same_edge {
                list_nodes_same_edge.push(node2.clone());
            }
            // Append the near node to the sequence (`cxx:3131-3149`).
            let pt2 = BRepTool::vertex_point_world(&Vertex(node2.clone()));
            let dist = pt1.distance(&pt2);
            if dist < tolerance {
                let mut is_ins = false;
                for kk in 1..=seq_nodes.len() {
                    let pt = BRepTool::vertex_point_world(&Vertex(seq_nodes[kk - 1].clone()));
                    if dist < pt1.distance(&pt) {
                        seq_nodes.insert(kk - 1, node2.clone());
                        is_ins = true;
                        break;
                    }
                }
                if !is_ins {
                    seq_nodes.push(node2.clone());
                }
            }
        }
        if !seq_nodes.is_empty() {
            // Remove the nodes near to some other from the same edge
            // (`cxx:3151-3175`).
            for n2 in &list_nodes_same_edge {
                let p2 = BRepTool::vertex_point_world(&Vertex(n2.clone()));
                let mut k = 1usize;
                while k <= seq_nodes.len() {
                    let n1 = &seq_nodes[k - 1];
                    let keep = if !crate::topo_tools::is_same(n1, n2) {
                        let p1 = BRepTool::vertex_point_world(&Vertex(n1.clone()));
                        p2.distance(&p1) >= pt1.distance(&p1)
                    } else {
                        false
                    };
                    if keep {
                        k += 1;
                        continue;
                    }
                    seq_nodes.remove(k - 1);
                }
            }
            if !seq_nodes.is_empty() {
                node_nearest_node.add(&node1, seq_nodes[0].clone());
            }
        }
        inspector.clear_res_list();
    }
    if node_nearest_node.is_empty() {
        return false;
    }
    create_new_nodes(&node_nearest_node, &node_vertices, vertex_node, node_edges)
}

impl Sewing {
    /// `VerticesAssembling` (`cxx:3548-3602`).
    pub(crate) fn vertices_assembling(&mut self) {
        let nb_vert = self.vertex_node.extent();
        let nb_vert_free = self.vertex_node_free.extent();
        if nb_vert == 0 && nb_vert_free == 0 {
            return;
        }
        // Fill the map node -> sections (`cxx:3196-3213`).
        let bounds = self.bound_faces.keys();
        for bound in bounds {
            for node in Sewing::raw_children(&bound) {
                self.node_sections.change(&node).push(bound.clone());
            }
        }
        // Glue the vertices (`cxx:3215-3234`).
        if nb_vert > 0 {
            while glue_vertices(
                &mut self.vertex_node,
                &mut self.node_sections,
                &self.bound_faces,
                self.tolerance,
            ) {}
        }
        if nb_vert_free > 0 {
            while glue_vertices(
                &mut self.vertex_node_free,
                &mut self.node_sections,
                &self.bound_faces,
                self.tolerance,
            ) {}
        }
    }
}

/// `replaceNMVertices(edge, V1, V2, theReShape)` (`cxx:3608-3700`).
pub(crate) fn replace_nm_vertices(
    edge: &TopoShape,
    v1: &TopoShape,
    v2: &TopoShape,
    reshape: &mut MapReShape,
) {
    let (nm_vert, nm_pars) = super::same_parameter::find_nm_vertices(edge);
    if nm_vert.is_empty() {
        return;
    }
    let (first, last) = GeometryRegistry::global().edge_parameters(edge);
    let Some(c3d) = GeometryRegistry::global().edge_curve(edge) else {
        return;
    };
    let _ = c3d;
    let mut ed_vert: Vec<TopoShape> = Vec::new();
    let mut ed_pars: Vec<f64> = Vec::new();
    let nb = nm_pars.len();
    for i in 0..nb {
        let apar = nm_pars[i];
        if (apar - first).abs() <= occt_core::precision::PCONFUSION {
            reshape.replace(&nm_vert[i], v1);
            continue;
        }
        if (apar - last).abs() <= occt_core::precision::PCONFUSION {
            reshape.replace(&nm_vert[i], v2);
            continue;
        }
        let a_v = nm_vert[i].clone();
        let mut inserted = false;
        let mut j = 0usize;
        while j < ed_pars.len() {
            let apar2 = ed_pars[j];
            if (apar - apar2).abs() <= occt_core::precision::PCONFUSION {
                reshape.replace(&a_v, &ed_vert[j]);
                inserted = true;
                break;
            } else if apar < apar2 {
                // `BRep_PointOnCurve(apar, c3d, aLoc)` on a copy of the vertex.
                let anew_v = a_v.empty_copied();
                GeometryRegistry::global().set_vertex(
                    &anew_v,
                    crate::tgeometry::VertexGeom {
                        point: BRepTool::vertex_point_world(&Vertex(a_v.clone())),
                        tolerance: BRepTool::vertex_tolerance(&Vertex(a_v.clone())),
                    },
                );
                GeometryRegistry::global().set_vertex_parameter_on_edge(edge, &anew_v, apar);
                ed_vert.insert(j, anew_v.clone());
                ed_pars.insert(j, apar);
                reshape.replace(&a_v, &anew_v);
                inserted = true;
                break;
            }
            j += 1;
        }
        if !inserted {
            let anew_v = a_v.empty_copied();
            GeometryRegistry::global().set_vertex(
                &anew_v,
                crate::tgeometry::VertexGeom {
                    point: BRepTool::vertex_point_world(&Vertex(a_v.clone())),
                    tolerance: BRepTool::vertex_tolerance(&Vertex(a_v.clone())),
                },
            );
            GeometryRegistry::global().set_vertex_parameter_on_edge(edge, &anew_v, apar);
            ed_vert.push(anew_v.clone());
            ed_pars.push(apar);
            reshape.replace(&a_v, &anew_v);
        }
    }
    let newnb = ed_pars.len();
    if newnb < nb {
        // `cxx:3313-3323`: the edge is rebuilt with the surviving vertices.
        let mut anew_edge = crate::shape_fix_compose_shell::empty_copied_edge(edge);
        let an_ori = edge.orientation();
        anew_edge.set_orientation(Orientation::Forward);
        let builder = TopoBuilder::new();
        builder.add(&mut anew_edge, v1);
        builder.add(&mut anew_edge, v2);
        for v in &ed_vert {
            builder.add(&mut anew_edge, v);
        }
        anew_edge.set_orientation(an_ori);
        reshape.replace(edge, &anew_edge);
    }
}

/// `ReplaceEdge(oldEdge, theNewShape, aReShape)` (`cxx:3707-3780`).
pub(crate) fn replace_edge(
    old_edge: &TopoShape,
    the_new_shape: &TopoShape,
    reshape: &mut MapReShape,
) {
    let old_shape = reshape.apply_option(old_edge).unwrap_or_else(|| old_edge.clone());
    let new_shape = reshape
        .apply_option(the_new_shape)
        .unwrap_or_else(|| the_new_shape.clone());
    if crate::topo_tools::is_same(&old_shape, &new_shape) || reshape.is_recorded(&new_shape) {
        return;
    }
    reshape.replace(&old_shape, &new_shape);
    let (v1old, v2old) = match crate::topo_tools::edge_vertices(&Edge(old_shape.clone())) {
        (Some(a), Some(b)) => (a.0, b.0),
        _ => return,
    };
    let orold = old_shape.orientation();
    let ornew = orold;
    let mut v1new: Option<TopoShape> = None;
    let mut v2new: Option<TopoShape> = None;
    let mut ornew = ornew;
    if new_shape.shape_type() == ShapeType::Edge {
        let a_en = Edge(new_shape.clone());
        if let (Some(a), Some(b)) = crate::topo_tools::edge_vertices(&a_en) {
            v1new = Some(a.0.clone());
            v2new = Some(b.0.clone());
            ornew = a_en.0.orientation();
            replace_nm_vertices(&a_en.0, &a.0, &b.0, reshape);
        }
    } else if new_shape.shape_type() == ShapeType::Wire {
        for ed in crate::topo_tools::explorer_shapes(&new_shape, ShapeType::Edge) {
            let ed = Edge(ed);
            ornew = ed.0.orientation();
            if let (Some(a), Some(b)) = crate::topo_tools::edge_vertices(&ed) {
                replace_nm_vertices(&ed.0, &a.0, &b.0, reshape);
                if v1new.is_none() {
                    v1new = Some(a.0.clone());
                }
                v2new = Some(b.0.clone());
            }
        }
    }
    let (Some(mut v1new), Some(mut v2new)) = (v1new, v2new) else {
        return;
    };
    v1new.set_orientation(v1old.orientation());
    v2new.set_orientation(v2old.orientation());
    if crate::topo_tools::is_same(&v1old, &v2old)
        && !crate::topo_tools::is_same(&v1old, &v1new)
        && !reshape.is_recorded(&v1new)
    {
        reshape.replace(&v1old, &v1new);
        return;
    }
    if orold == ornew {
        v1new.set_orientation(v1old.orientation());
        v2new.set_orientation(v2old.orientation());
        if !crate::topo_tools::is_same(&v1old, &v1new)
            && !crate::topo_tools::is_same(&v1old, &v2new)
            && !reshape.is_recorded(&v1new)
        {
            reshape.replace(&v1old, &v1new);
        }
        if !crate::topo_tools::is_same(&v2old, &v2new)
            && !crate::topo_tools::is_same(&v2old, &v1new)
            && !reshape.is_recorded(&v2new)
        {
            reshape.replace(&v2old, &v2new);
        }
    } else {
        v1new.set_orientation(v2old.orientation());
        v2new.set_orientation(v1old.orientation());
        if !crate::topo_tools::is_same(&v1old, &v2new)
            && !crate::topo_tools::is_same(&v1old, &v1new)
            && !reshape.is_recorded(&v2new)
        {
            reshape.replace(&v1old, &v2new);
        }
        if !crate::topo_tools::is_same(&v2old, &v2new)
            && !crate::topo_tools::is_same(&v2old, &v1new)
            && !reshape.is_recorded(&v1new)
        {
            reshape.replace(&v2old, &v1new);
        }
    }
}




