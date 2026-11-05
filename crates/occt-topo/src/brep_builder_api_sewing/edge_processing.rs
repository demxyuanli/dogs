//! `GetSeqEdges` (`cxx:4604-4650`), `GetFreeWires` (`cxx:4653-4705`),
//! `IsDegeneratedWire` (`cxx:4709-4792`), `DegeneratedSection`
//! (`cxx:4800-4898`), `EdgeProcessing` (`cxx:4911-4996`), `EdgeRegularity`
//! (`cxx:5002-5026`) and `SameParameterShape` (`cxx:5881-5908`).

use occt_core::precision::CONFUSION;

use crate::abs::{Orientation, ShapeType};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::tgeometry::GeometryRegistry;

use super::shape_maps::ShapeKeyedMap;
use super::vertices::replace_edge;
use super::Sewing;

/// `GetSeqEdges(edge, seqEdges, VertEdge)` (`cxx:4604-4650`): walk the edges
/// chained to `edge` through their vertices.
fn get_seq_edges(
    edge: &TopoShape,
    seq_edges: &mut Vec<TopoShape>,
    vert_edge: &ShapeKeyedMap<Vec<TopoShape>>,
) {
    let mut num_v = 0usize;
    for v1 in Sewing::raw_children(edge) {
        if v1.shape_type() != ShapeType::Vertex {
            continue;
        }
        num_v += 1;
        if let Some(list_edges) = vert_edge.get(&v1) {
            for edge1 in &list_edges {
                if crate::topo_tools::is_same(edge1, edge) {
                    continue;
                }
                let mut is_contained = false;
                let mut index = 1usize;
                for (i, value) in seq_edges.iter().enumerate() {
                    is_contained = crate::topo_tools::is_same(value, edge1);
                    if !is_contained && crate::topo_tools::is_same(value, edge) {
                        index = i + 1;
                    }
                }
                if !is_contained {
                    if num_v == 1 {
                        seq_edges.insert(index - 1, edge1.clone());
                    } else {
                        seq_edges.insert(index, edge1.clone());
                    }
                    get_seq_edges(edge1, seq_edges, vert_edge);
                }
            }
        }
    }
}

/// `IsDegeneratedWire(wire)` (`cxx:4709-4792`).
fn is_degenerated_wire(wire: &TopoShape) -> bool {
    if wire.shape_type() != ShapeType::Wire {
        return false;
    }
    let mut wire_length = 0.0;
    let mut nume = 0usize;
    let mut is_small = 0usize;
    let mut v1: Option<TopoShape> = None;
    let mut v2: Option<TopoShape> = None;
    for edge in crated_children(wire) {
        nume += 1;
        let (ve1, ve2) = crate::topo_tools::edge_vertices(&Edge(edge.clone()));
        let (Some(ve1), Some(ve2)) = (ve1, ve2) else {
            continue;
        };
        let ve1 = ve1.0;
        let ve2 = ve2.0;
        if nume == 1 {
            v1 = Some(ve1.clone());
            v2 = Some(ve2.clone());
        } else {
            if let Some(cur) = &v1 {
                if crate::topo_tools::is_same(&ve1, cur) {
                    v1 = Some(ve2.clone());
                }
            }
            if let Some(cur) = &v2 {
                if crate::topo_tools::is_same(&ve1, cur) {
                    v2 = Some(ve2.clone());
                }
            }
            if let Some(cur) = &v1 {
                if crate::topo_tools::is_same(&ve2, cur) {
                    v1 = Some(ve1.clone());
                }
            }
            if let Some(cur) = &v2 {
                if crate::topo_tools::is_same(&ve2, cur) {
                    v2 = Some(ve1.clone());
                }
            }
        }
        if let Some(c3d) = BRepTool::edge_curve_adaptor_world(&Edge(edge.clone())) {
            let (first, last) = BRepTool::edge_parameters(&Edge(edge.clone()));
            let pfirst = c3d.d0(first);
            let plast = c3d.d0(last);
            let pmid = c3d.d0((first + last) * 0.5);
            let length = if pfirst.distance(&plast) > pfirst.distance(&pmid) {
                pfirst.distance(&plast)
            } else {
                occt_geom::gcpnts::curve_length_range(c3d.as_ref(), first, last, CONFUSION)
            };
            let tole = BRepTool::vertex_tolerance(&Vertex(ve1))
                + BRepTool::vertex_tolerance(&Vertex(ve2));
            if length <= tole {
                is_small += 1;
            }
            wire_length += length;
        }
    }
    if is_small == nume {
        return true;
    }
    let tol = match (&v1, &v2) {
        (Some(a), Some(b)) => {
            BRepTool::vertex_tolerance(&Vertex(a.clone()))
                + BRepTool::vertex_tolerance(&Vertex(b.clone()))
        }
        _ => 0.0,
    };
    if wire_length > tol {
        return false;
    }
    true
}

/// `TopoDS_Iterator(shape, Standard_False)`.
fn crated_children(shape: &TopoShape) -> Vec<TopoShape> {
    shape.tshape.read().expect("poisoned TShape lock").children.clone()
}

/// `DegeneratedSection(section, face)` (`cxx:4800-4898`).
fn degenerated_section(section: &TopoShape, face: Option<&TopoShape>) -> Option<TopoShape> {
    let sec = Edge(section.clone());
    // Return if the section is already degenerated (`cxx:4296-4298`).
    if BRepTool::is_degenerated(&sec) {
        return Some(section.clone());
    }
    let (first, last) = BRepTool::edge_parameters(&sec);
    let Some(c3d) = BRepTool::edge_curve_adaptor_world(&sec) else {
        // `cxx:4306-4311`: no 3D curve - mark it degenerated in place.
        let edge1 = section.clone();
        GeometryRegistry::global().set_degenerated(&edge1, true);
        return Some(edge1);
    };
    // Test if the new edge is degenerated (`cxx:4320-4338`).
    let (v1, v2) = match crate::topo_tools::edge_vertices(&sec) {
        (Some(a), Some(b)) => (a, b),
        _ => return None,
    };
    let p1 = BRepTool::vertex_point_world(&v1);
    let p3 = BRepTool::vertex_point_world(&v2);
    let p2 = c3d.d0(0.5 * (first + last));
    // Processing (`cxx:4340-4367`).
    let edge = crate::shape_fix_compose_shell::empty_copied_edge(section);
    let builder = TopoBuilder::new();
    if crate::topo_tools::is_same(&v1.0, &v2.0) {
        let mut an_edge = edge.oriented(Orientation::Forward);
        builder.add(&mut an_edge, &v1.0.oriented(Orientation::Forward));
        builder.add(&mut an_edge, &v2.0.oriented(Orientation::Reversed));
    } else {
        let new_vertex = if p1.distance(&p3) < BRepTool::vertex_tolerance(&v1) {
            v1.0.clone()
        } else if p1.distance(&p3) < BRepTool::vertex_tolerance(&v2) {
            v2.0.clone()
        } else {
            let d1 = BRepTool::vertex_tolerance(&v1) + p2.distance(&p1);
            let d2 = BRepTool::vertex_tolerance(&v2) + p2.distance(&p3);
            let new_tolerance = d1.max(d2);
            builder.make_vertex(p2, new_tolerance).0
        };
        let mut an_edge = edge.oriented(Orientation::Forward);
        builder.add(&mut an_edge, &new_vertex.oriented(Orientation::Forward));
        builder.add(&mut an_edge, &new_vertex.oriented(Orientation::Reversed));
    }
    let (rfirst, rlast) = BRepTool::edge_parameters(&sec);
    GeometryRegistry::global().set_edge_range(&edge, rfirst, rlast);
    GeometryRegistry::global().set_degenerated(&edge, true);
    if let Some(face) = face {
        let fac = Face(face.clone());
        if let Some((a_c2dt, af, al)) = crate::boptools_2d::curve_on_surface_range(&sec, &fac) {
            // `aBuilder.UpdateEdge(edge, null, 0)` then, when the new edge has
            // no pcurve, `UpdateEdge(edge, aC2dt, face, 0)` (`cxx:4378-4383`).
            if crate::boptools_2d::curve_on_surface(&Edge(edge.clone()), &fac).is_none() {
                GeometryRegistry::global().set_edge_pcurve(
                    &edge,
                    GeometryRegistry::shape_key(&fac.0),
                    a_c2dt,
                );
                GeometryRegistry::global().set_pcurve_range(
                    &edge,
                    GeometryRegistry::shape_key(&fac.0),
                    af,
                    al,
                );
            }
        }
    }
    Some(edge)
}

impl Sewing {
    /// `GetFreeWires(MapFreeEdges, seqWires)` (`cxx:4653-4705`).
    fn get_free_wires(&self, map_free_edges: &mut ShapeVecSet) -> Vec<TopoShape> {
        let mut vert_edge: ShapeKeyedMap<Vec<TopoShape>> = ShapeKeyedMap::new();
        let mut seq_free_edges: Vec<TopoShape> = Vec::new();
        for edge in map_free_edges.keys() {
            seq_free_edges.push(edge.clone());
            for v in crated_children(&edge) {
                if v.shape_type() != ShapeType::Vertex {
                    continue;
                }
                vert_edge.change(&v).push(edge.clone());
            }
        }
        let builder = TopoBuilder::new();
        let mut seq_wires: Vec<TopoShape> = Vec::new();
        let mut i = 0usize;
        while i < seq_free_edges.len() {
            let edge = seq_free_edges[i].clone();
            i += 1;
            if !map_free_edges.contains(&edge) {
                continue;
            }
            let mut seq_edges: Vec<TopoShape> = vec![edge.clone()];
            get_seq_edges(&edge, &mut seq_edges, &vert_edge);
            let mut wire = TopoShape::new(ShapeType::Wire);
            for e in &seq_edges {
                builder.add(&mut wire, e);
                map_free_edges.remove(e);
            }
            seq_wires.push(wire);
            if map_free_edges.is_empty() {
                break;
            }
        }
        seq_wires
    }

    /// `EdgeProcessing` (`cxx:4911-4996`): after the merging, every free bound
    /// that ended up in a degenerated wire is turned into a degenerated
    /// section.
    pub(crate) fn edge_processing(&mut self) {
        let mut map_free_edges = ShapeVecSet::new();
        let mut edge_face: ShapeKeyedMap<TopoShape> = ShapeKeyedMap::new();
        let bounds = self.bound_faces.keys();
        for bound in bounds {
            let list_faces = self.bound_faces.get(&bound).unwrap_or_default();
            if list_faces.len() != 1 {
                continue;
            }
            if self.bound_sections.contains(&bound) {
                let sections = self.bound_sections.get(&bound).unwrap_or_default();
                for section in &sections {
                    if !self.merged_edges.contains(section) {
                        let edge = self
                            .reshape
                            .apply_option(section)
                            .unwrap_or_else(|| section.clone());
                        if !map_free_edges.contains(&edge) {
                            let face = list_faces[0].clone();
                            edge_face.bind(&edge, face);
                            map_free_edges.push(edge);
                        }
                    }
                }
            } else if !self.merged_edges.contains(&bound) {
                let edge = self
                    .reshape
                    .apply_option(&bound)
                    .unwrap_or_else(|| bound.clone());
                if !map_free_edges.contains(&edge) {
                    let face = list_faces[0].clone();
                    edge_face.bind(&edge, face);
                    map_free_edges.push(edge);
                }
            }
        }
        if !map_free_edges.is_empty() {
            let seq_wires = self.get_free_wires(&mut map_free_edges);
            for wire in seq_wires {
                if !is_degenerated_wire(&wire) {
                    continue;
                }
                for edge_in in crated_children(&wire) {
                    let edge = self
                        .reshape
                        .apply_option(&edge_in)
                        .unwrap_or_else(|| edge_in.clone());
                    let face = edge_face.get(&edge);
                    let Some(degedge) = degenerated_section(&edge, face.as_ref()) else {
                        continue;
                    };
                    if !crate::topo_tools::is_same(&degedge, &edge) {
                        let mut reshape = std::mem::take(&mut self.reshape);
                        replace_edge(&edge, &degedge, &mut reshape);
                        self.reshape = reshape;
                    }
                    if BRepTool::is_degenerated(&Edge(degedge.clone())) {
                        self.degenerated.add(&degedge);
                    }
                }
            }
        }
    }

    /// `EdgeRegularity` (`cxx:5002-5026`): re-encode the regularity of the
    /// merged edges that end up shared by exactly two faces.
    pub(crate) fn edge_regularity(&mut self) {
        let Some(sewed) = self.sewed_shape.clone() else {
            return;
        };
        let map_ef = crate::topo_tools::map_shapes_and_ancestors(
            &sewed,
            ShapeType::Edge,
            ShapeType::Face,
        );
        let merged: Vec<TopoShape> = self.merged_edges.keys();
        for edge in merged {
            let Some(an_edge) = self.reshape.apply_option(&edge) else {
                continue;
            };
            let key = crate::topo_tools::shape_map_key(&an_edge);
            if let Some(faces) = map_ef.get(&key) {
                if faces.len() == 2 {
                    crate::brep_lib_continuity::encode_regularity_edge(
                        &Edge(an_edge),
                        &Face(faces[0].clone()),
                        &Face(faces[1].clone()),
                        occt_core::precision::ANGULAR,
                    );
                }
            }
        }
        self.merged_edges.clear();
    }

    /// `SameParameterShape` (`cxx:5881-5908`).
    pub(crate) fn same_parameter_shape(&mut self) {
        if !self.same_parameter_mode {
            return;
        }
        let Some(sewed) = self.sewed_shape.clone() else {
            return;
        };
        for edge in crate::topo_tools::explorer_shapes(&sewed, ShapeType::Edge) {
            let sec = Edge(edge);
            let tol = BRepTool::edge_tolerance(&sec);
            crate::brep_lib_same_parameter::same_parameter_edge_inplace(&sec, tol);
        }
    }
}

/// A `TopTools_IndexedMapOfShape` used as an ordered set by
/// `GetFreeWires` (`cxx:4653-4705`).
pub(crate) struct ShapeVecSet {
    items: Vec<TopoShape>,
}

impl ShapeVecSet {
    pub(crate) fn new() -> Self {
        Self { items: Vec::new() }
    }
    pub(crate) fn push(&mut self, s: TopoShape) {
        self.items.push(s);
    }
    pub(crate) fn contains(&self, s: &TopoShape) -> bool {
        self.items.iter().any(|x| crate::topo_tools::is_same(x, s))
    }
    pub(crate) fn remove(&mut self, s: &TopoShape) {
        self.items.retain(|x| !crate::topo_tools::is_same(x, s));
    }
    pub(crate) fn is_empty(&self) -> bool {
        self.items.is_empty()
    }
    pub(crate) fn keys(&self) -> Vec<TopoShape> {
        self.items.clone()
    }
}




