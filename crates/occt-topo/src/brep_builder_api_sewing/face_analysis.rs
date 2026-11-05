//! `FaceAnalysis` (`cxx:2596-2867`) and `FindFreeBoundaries`
//! (`cxx:2877-3074`).

use occt_core::gp::GpPnt;
use occt_core::precision::CONFUSION;

use crate::abs::{Orientation, ShapeType};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::shape_fix_compose_shell::{empty_copied_edge, empty_copied_face};
use crate::tgeometry::GeometryRegistry;

use super::shape_maps::{ShapeKeyedMap, ShapeSet};
use super::Sewing;

impl Sewing {
    /// `FaceAnalysis` (`cxx:2596-2867`): remove the too-small edges (they are
    /// made degenerated), glue their endpoints and drop the faces whose every
    /// edge is small.
    pub(crate) fn face_analysis(&mut self) {
        if self.shape.is_some() && self.old_shapes.is_empty() {
            let shape = self.shape.clone().unwrap();
            self.add(&shape);
            self.shape = None;
        }
        let builder = TopoBuilder::new();
        let mut small_edges = ShapeSet::new();
        let mut glued_vertices: ShapeKeyedMap<Vec<TopoShape>> = ShapeKeyedMap::new();
        let nb_shapes = self.old_shapes.extent();
        for i in 1..=nb_shapes {
            let Some(shape) = self.old_shapes.key_at(i).cloned() else {
                continue;
            };
            for face_shape in crate::topo_tools::explorer_shapes(&shape, ShapeType::Face) {
                let face = Face(face_shape);
                let mut nb_edges = 0usize;
                let mut nb_small = 0usize;
                // Build the replacing face (`cxx:2376-2379`).
                let mut nface = empty_copied_face(&face.0).oriented(Orientation::Forward);
                let mut is_face_changed = false;
                let forward_face = face.0.oriented(Orientation::Forward);
                for wire_shape in crate::iterator::cumulated_children(&forward_face) {
                    if wire_shape.shape_type() != ShapeType::Wire {
                        continue;
                    }
                    let wire = wire_shape;
                    // Build the replacing wire (`cxx:2392-2394`).
                    let mut nwire = wire.empty_copied().oriented(Orientation::Forward);
                    let mut is_wire_changed = false;
                    let forward_wire = wire.oriented(Orientation::Forward);
                    for edge_shape in crate::iterator::cumulated_children(&forward_wire) {
                        let edge = Edge(edge_shape);
                        nb_edges += 1;
                        // Process a degenerated edge (`cxx:2406-2411`).
                        if BRepTool::is_degenerated(&edge) {
                            builder.add(&mut nwire, &edge.0); // Old edge kept.
                            self.degenerated.add(&edge.0);
                            nb_small += 1;
                            continue;
                        }
                        let mut is_small = small_edges.contains(&edge.0);
                        if !is_small {
                            // Check for a small edge (`cxx:2419-2457`).
                            if let Some(c3d) = BRepTool::edge_curve(&edge) {
                                let (first, last) = BRepTool::edge_parameters(&edge);
                                let npt = 5;
                                let pf = c3d.d0(first);
                                let pl = c3d.d0(last);
                                let cp = GpPnt::new(
                                    (pf.x() + pl.x()) * 0.5,
                                    (pf.y() + pl.y()) * 0.5,
                                    (pf.z() + pl.z()) * 0.5,
                                );
                                let mut maxdist = 0.0f64;
                                let delta = (last - first) / (npt - 1) as f64;
                                for idx in 0..npt {
                                    let dist = cp.distance(&c3d.d0(first + idx as f64 * delta));
                                    if maxdist < dist {
                                        maxdist = dist;
                                    }
                                }
                                is_small = 2.0 * maxdist <= self.min_tolerance;
                            }
                            if is_small {
                                // Store the small edge in the map (`cxx:2461`).
                                small_edges.add(&edge.0);
                                let (v1, v2) =
                                    match crate::topo_tools::edge_vertices(&edge) {
                                        (Some(a), Some(b)) => (a.0, b.0),
                                        _ => (TopoShape::new(ShapeType::Vertex), TopoShape::new(ShapeType::Vertex)),
                                    };
                                let nv1 = self.reshape.apply_option(&v1).unwrap_or_else(|| v1.clone());
                                let nv2 = self.reshape.apply_option(&v2).unwrap_or_else(|| v2.clone());
                                // Store the glued vertices (`cxx:2464-2509`).
                                if !crate::topo_tools::is_same(&nv1, &v1) {
                                    let mut vlist1 = glued_vertices.change(&nv1).clone();
                                    if !crate::topo_tools::is_same(&nv2, &v2) {
                                        // Merge the lists of glued vertices.
                                        if !crate::topo_tools::is_same(&nv1, &nv2) {
                                            if let Some(list) = glued_vertices.get(&nv2) {
                                                for v in list {
                                                    self.reshape.replace(
                                                        &v,
                                                        &nv1.oriented(v.orientation()),
                                                    );
                                                    vlist1.push(v);
                                                }
                                            }
                                            glued_vertices.remove(&nv2);
                                        }
                                    } else {
                                        vlist1.push(v2.clone());
                                        self.reshape.replace(
                                            &v2,
                                            &nv1.oriented(v2.orientation()),
                                        );
                                    }
                                    glued_vertices.bind(&nv1, vlist1);
                                } else if !crate::topo_tools::is_same(&nv2, &v2) {
                                    // Add the first vertex to the existing list.
                                    glued_vertices.change(&nv2).push(v1.clone());
                                    self.reshape.replace(&v1, &nv2.oriented(v1.orientation()));
                                } else if !crate::topo_tools::is_same(&v1, &v2) {
                                    // Record the new glued vertices.
                                    let nv = builder.make_vertex(GpPnt::zero(), 0.0);
                                    let vlist = vec![v1.clone(), v2.clone()];
                                    glued_vertices.add(&nv.0, vlist);
                                    self.reshape.replace(&v1, &nv.0.oriented(v1.orientation()));
                                    self.reshape.replace(&v2, &nv.0.oriented(v2.orientation()));
                                }
                            }
                        }
                        if is_small {
                            // Replace the small edge (`cxx:2514-2540`).
                            nb_small += 1;
                            let fedge = edge.0.oriented(Orientation::Forward);
                            let fe = Edge(fedge);
                            if let Some((c2d, pfirst, plast)) =
                                crate::boptools_2d::curve_on_surface_range(&fe, &face)
                            {
                                let nedge = TopoShape::new(ShapeType::Edge);
                                let fk = GeometryRegistry::shape_key(&face.0);
                                GeometryRegistry::global().set_edge_pcurve(&nedge, fk, c2d);
                                GeometryRegistry::global().set_edge_range(&nedge, pfirst, plast);
                                GeometryRegistry::global().set_degenerated(&nedge, true);
                                let (v1, v2) = match crate::topo_tools::edge_vertices(&fe) {
                                    (Some(a), Some(b)) => (a.0, b.0),
                                    _ => continue,
                                };
                                let a1 = self.reshape.apply_option(&v1).unwrap_or(v1.clone());
                                let a2 = self.reshape.apply_option(&v2).unwrap_or(v2.clone());
                                builder.add(&mut nedge.clone(), &a1.oriented(v1.orientation()));
                                builder.add(&mut nedge.clone(), &a2.oriented(v2.orientation()));
                                builder.add(&mut nwire, &nedge.oriented(edge.0.orientation()));
                                self.degenerated.add(&nedge);
                            }
                            is_wire_changed = true;
                        } else {
                            builder.add(&mut nwire, &edge.0); // Old edge kept.
                        }
                    }
                    if is_wire_changed {
                        builder.add(&mut nface, &nwire.oriented(wire.orientation()));
                        is_face_changed = true;
                    } else {
                        builder.add(&mut nface, &wire);
                    }
                }
                // Remove the small face (`cxx:2556-2564`).
                if nb_small == nb_edges {
                    self.little_face.add(&face.0);
                    self.reshape.remove(&face.0);
                } else if is_face_changed {
                    self.reshape
                        .replace(&face.0, &nface.oriented(face.0.orientation()));
                }
            }
        }

        // Update the glued vertices (`cxx:2572-2602`).
        let glued: Vec<(TopoShape, Vec<TopoShape>)> = glued_vertices
            .iter()
            .map(|(k, v)| (k.clone(), v.clone()))
            .collect();
        for (vnew, vlist) in glued {
            let mut x = 0.0;
            let mut y = 0.0;
            let mut z = 0.0;
            let mut nb_points = 0usize;
            for v in &vlist {
                let p = BRepTool::vertex_point_world(&Vertex(v.clone()));
                x += p.x();
                y += p.y();
                z += p.z();
                nb_points += 1;
            }
            if nb_points > 0 {
                let vp = GpPnt::new(x / nb_points as f64, y / nb_points as f64, z / nb_points as f64);
                let mut tol = 0.0f64;
                let mut mtol = 0.0f64;
                for v in &vlist {
                    let vtol = BRepTool::vertex_tolerance(&Vertex(v.clone()));
                    if mtol < vtol {
                        mtol = vtol;
                    }
                    let d = vp.distance(&BRepTool::vertex_point_world(&Vertex(v.clone())));
                    if tol < d {
                        tol = d;
                    }
                }
                GeometryRegistry::global().set_vertex(
                    &vnew,
                    crate::tgeometry::VertexGeom {
                        point: vp,
                        tolerance: tol + mtol,
                    },
                );
            }
        }

        // Update the input shapes (`cxx:2604-2605`).
        let keys: Vec<TopoShape> = self.old_shapes.keys();
        for key in keys {
            let applied = self.reshape.apply_option(&key).unwrap_or_else(|| key.clone());
            self.old_shapes.bind(&key, applied);
        }
    }

    /// `FindFreeBoundaries` (`cxx:2877-3074`): the edge -> faces map, the
    /// bound map and the two vertex/node maps.
    pub(crate) fn find_free_boundaries(&mut self) {
        // Take the context shape into account (`cxx:2619-2637`).
        let mut new_shapes: Vec<TopoShape> = Vec::new();
        if let Some(shape) = self.shape.clone() {
            if self.old_shapes.is_empty() {
                self.add(&shape);
                self.shape = None;
            } else if let Some(applied) = self.reshape.apply_option(&shape) {
                if !new_shapes.iter().any(|s| crate::topo_tools::is_same(s, &applied)) {
                    new_shapes.push(applied);
                }
            }
        }
        // Create the map edge -> faces (`cxx:2638-2657`).
        let mut edge_faces: ShapeKeyedMap<Vec<TopoShape>> = ShapeKeyedMap::new();
        let nb_shapes = self.old_shapes.extent();
        for i in 1..=nb_shapes {
            let Some(shape) = self.old_shapes.key_at(i).cloned() else {
                continue;
            };
            if !new_shapes.iter().any(|s| crate::topo_tools::is_same(s, &shape)) {
                new_shapes.push(shape.clone());
            }
            for edge in crate::topo_tools::explorer_shapes(&shape, ShapeType::Edge) {
                if !edge_faces.contains(&edge) {
                    edge_faces.add(&edge, Vec::new());
                }
            }
        }
        // Fill the map edge -> faces (`cxx:2658-2694`).
        let mut map_faces = ShapeSet::new();
        for shape in &new_shapes {
            for face in crate::topo_tools::explorer_shapes(shape, ShapeType::Face) {
                if !map_faces.add(&face) {
                    continue;
                }
                for wire in crate::iterator::cumulated_children(&face) {
                    if wire.shape_type() != ShapeType::Wire {
                        continue;
                    }
                    for edge in crate::iterator::cumulated_children(&wire) {
                        if edge_faces.contains(&edge) {
                            edge_faces.change(&edge).push(face.clone());
                        }
                    }
                }
            }
        }
        // Find the free boundaries (`cxx:2695-2780`).
        let entries: Vec<(TopoShape, Vec<TopoShape>)> = edge_faces
            .iter()
            .map(|(k, v)| (k.clone(), v.clone()))
            .collect();
        for (edge_in, list_faces) in entries {
            let mut edge = edge_in;
            let nb_faces = list_faces.len();
            if edge.orientation() == Orientation::Internal {
                continue;
            }
            let mut is_seam = false;
            if nb_faces == 1 {
                let face = Face(list_faces[0].clone());
                is_seam = BRepTool::is_closed_edge_face(&Edge(edge.clone()), &face);
                if is_seam {
                    // `cxx:2712-2739`: a seam read as a single-face bound is
                    // re-cut into a copy carrying one pcurve.
                    let mut anew_edge = empty_copied_edge(&edge);
                    anew_edge.set_orientation(Orientation::Forward);
                    let builder = TopoBuilder::new();
                    for v in crate::iterator::cumulated_children(&edge) {
                        builder.add(&mut anew_edge, &v);
                    }
                    if let Some((c2dold, first2d, last2d)) =
                        crate::boptools_2d::curve_on_surface_range(&Edge(edge.clone()), &face)
                    {
                        let fk = GeometryRegistry::shape_key(&face.0);
                        GeometryRegistry::global().set_edge_pcurve(&anew_edge, fk, c2dold);
                        let (a_first, a_last) = BRepTool::edge_parameters(&Edge(edge.clone()));
                        GeometryRegistry::global().set_edge_range(&anew_edge, a_first, a_last);
                        GeometryRegistry::global().set_pcurve_range(
                            &anew_edge, fk, first2d, last2d,
                        );
                    }
                    self.reshape.replace(&edge, &anew_edge);
                    edge = anew_edge;
                    is_seam = false;
                }
            }
            let is_bound_float = self.floating_edges_mode && nb_faces == 0;
            let is_bound = self.face_mode
                && ((self.nonmanifold && nb_faces > 0) || (nb_faces == 1 && !is_seam));
            if is_bound || is_bound_float {
                if BRepTool::is_degenerated(&Edge(edge.clone())) {
                    continue;
                }
                self.bound_faces.bind(&edge, list_faces.clone());
                let (v_first, v_last) = match crate::topo_tools::edge_vertices(&Edge(edge.clone())) {
                    (Some(a), Some(b)) => (a, b),
                    _ => continue,
                };
                if v_first.0.orientation() == Orientation::Internal
                    || v_last.0.orientation() == Orientation::Internal
                {
                    continue;
                }
                if is_bound {
                    if !self.vertex_node.contains(&v_first.0) {
                        self.vertex_node.add(&v_first.0, v_first.0.clone());
                    }
                    if !self.vertex_node.contains(&v_last.0) {
                        self.vertex_node.add(&v_last.0, v_last.0.clone());
                    }
                } else {
                    if !self.vertex_node_free.contains(&v_first.0) {
                        self.vertex_node_free.add(&v_first.0, v_first.0.clone());
                    }
                    if !self.vertex_node_free.contains(&v_last.0) {
                        self.vertex_node_free.add(&v_last.0, v_last.0.clone());
                    }
                }
            }
        }
        let _ = CONFUSION;
    }
}

