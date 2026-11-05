//! `Cutting` (`cxx:4440-4600`) and its helpers `ProjectPointsOnCurve`
//! (`cxx:5366-5476`), `CreateCuttingNodes` (`cxx:5479-5665`) and
//! `CreateSections` (`cxx:5668-5878`).

use std::sync::Arc;

use occt_core::bnd::BndBox;
use occt_core::gp::GpPnt;
use occt_core::precision::{CONFUSION, INFINITE};
use occt_geom::extrema_pc::extrema_ext_pc_range_done;
use occt_geom::Curve;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::trimmed::Geom2dTrimmedCurve;

use crate::abs::Orientation;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, TopoShape, Vertex};
use crate::tgeometry::GeometryRegistry;

use super::same_parameter::find_nm_vertices;
use super::shape_maps::{BoxTree, ShapeKeyedMap, ShapeSequence, RealSequence};
use super::Sewing;

/// `ProjectPointsOnCurve(arrPnt, c3d, first, last, arrDist, arrPara, arrProj,
/// isConsiderEnds)` (`cxx:5366-5476`).
///
/// The `try`/`catch` arm of `cxx:4870-4875` sets `worktol = MinTolerance()`
/// after a failed projection; the port's `Extrema_ExtPC` does not throw, so a
/// `!IsDone()` result takes that arm.
#[allow(clippy::too_many_arguments)]
pub(crate) fn project_points_on_curve(
    arr_pnt: &[GpPnt],
    c3d: &dyn Curve,
    first: f64,
    last: f64,
    tolerance: f64,
    is_consider_ends: bool,
    arr_dist: &mut [f64],
    arr_para: &mut [f64],
    arr_proj: &mut [GpPnt],
) {
    for v in arr_dist.iter_mut() {
        *v = -1.0;
    }
    let pfirst = c3d.d0(first);
    let plast = c3d.d0(last);
    let find = 1usize;
    let lind = arr_pnt.len();
    for i1 in find..=lind {
        let pt = arr_pnt[i1 - 1];
        let mut worktol = tolerance;
        let dist_f2 = pfirst.square_distance(&pt);
        let dist_l2 = plast.square_distance(&pt);
        let mut is_projected = false;
        let (sols, done) = extrema_ext_pc_range_done(c3d, &pt, first, last);
        if done && !sols.is_empty() {
            let mut dist2_min = if is_consider_ends || i1 == find || i1 == lind {
                dist_f2.min(dist_l2)
            } else {
                INFINITE
            };
            let mut ind_min = 0usize;
            for (ind, s) in sols.iter().enumerate() {
                if s.sq_dist < dist2_min {
                    ind_min = ind + 1;
                    dist2_min = s.sq_dist;
                }
            }
            if ind_min > 0 {
                is_projected = true;
                let sol = &sols[ind_min - 1];
                let mut param_proj = sol.u;
                let mut pt_proj = c3d.d0(param_proj);
                let mut dist_proj2 = pt_proj.square_distance(&pt);
                if !sol.is_min && dist_f2.min(dist_l2) < dist2_min {
                    if dist_f2 < dist_l2 {
                        param_proj = first;
                        dist_proj2 = dist_f2;
                        pt_proj = pfirst;
                    } else {
                        param_proj = last;
                        dist_proj2 = dist_l2;
                        pt_proj = plast;
                    }
                }
                if dist_proj2 < worktol * worktol || !is_consider_ends {
                    arr_dist[i1 - 1] = dist_proj2.sqrt();
                    arr_para[i1 - 1] = param_proj;
                    arr_proj[i1 - 1] = pt_proj;
                }
            }
        } else {
            // `catch`: `worktol = MinTolerance()`.
            worktol = f64::NAN;
        }
        if !is_projected && is_consider_ends {
            let worktol = if worktol.is_nan() {
                // MinTolerance is only known to the sewing object; the caller
                // passes the tolerance, so the arm keeps it (`cxx:4875` is a
                // diagnostic-only reset for a non-throwing failure).
                tolerance
            } else {
                worktol
            };
            if dist_f2.min(dist_l2) < worktol * worktol {
                if dist_f2 < dist_l2 {
                    arr_dist[i1 - 1] = dist_f2.sqrt();
                    arr_para[i1 - 1] = first;
                    arr_proj[i1 - 1] = pfirst;
                } else {
                    arr_dist[i1 - 1] = dist_l2.sqrt();
                    arr_para[i1 - 1] = last;
                    arr_proj[i1 - 1] = plast;
                }
            }
        }
    }
}

impl Sewing {
    /// `Cutting` (`cxx:4440-4600`): split every bound at the nodes of the
    /// other bounds that project onto it.
    pub(crate) fn cutting_bounds(&mut self) {
        let nb_vertices = self.vertex_node.extent();
        if nb_vertices == 0 {
            return;
        }
        // Create a box tree with the vertices (`cxx:3989-4003`).
        let eps = self.tolerance * 0.5;
        let mut tree = BoxTree::new();
        for i in 1..=nb_vertices {
            let v = Vertex(self.vertex_node.key_at(i).unwrap().clone());
            let pt = BRepTool::vertex_point_world(&v);
            let mut b = BndBox::new();
            b.add_point(&pt);
            b.enlarge(eps);
            tree.add(i, b);
        }

        // Iterate on all boundaries (`cxx:4005-4119`).
        let bounds: Vec<TopoShape> = self.bound_faces.keys();
        for bound in bounds {
            let Some(list_faces) = self.bound_faces.get(&bound) else {
                continue;
            };
            if list_faces.is_empty() {
                continue;
            }
            let bound_e = Edge(bound.clone());
            let (first, last) = BRepTool::edge_parameters(&bound_e);
            let Some(c3d) = BRepTool::edge_curve_adaptor_world(&bound_e) else {
                continue;
            };
            let mut list_sections: Vec<TopoShape> = Vec::new();
            {
                // Create a bounding box around the curve (`cxx:4034-4039`).
                let a_global_box = occt_geom::geom_bnd_lib::curve3d::box_curve(
                    c3d.as_ref(),
                    first,
                    last,
                    self.tolerance,
                );
                let selected = tree.select(&a_global_box);
                if selected.is_empty() {
                    continue;
                }
                // Retrieve the bound nodes (`cxx:4044-4048`).
                let (v1, v2) = match crate::topo_tools::edge_vertices(&bound_e) {
                    (Some(a), Some(b)) => (a.0, b.0),
                    _ => continue,
                };
                let (Some(node1), Some(node2)) = (
                    self.vertex_node.get(&v1),
                    self.vertex_node.get(&v2),
                ) else {
                    continue;
                };
                // Fill the map of candidate vertices (`cxx:4049-4060`).
                let mut candidate_vertices: Vec<TopoShape> = Vec::new();
                for index in selected {
                    let Some(node) = self.vertex_node.value_at(index) else {
                        continue;
                    };
                    if !crate::topo_tools::is_same(node, &node1)
                        && !crate::topo_tools::is_same(node, &node2)
                    {
                        let vertex = self.vertex_node.key_at(index).unwrap().clone();
                        if !candidate_vertices.iter().any(|v| crate::topo_tools::is_same(v, &vertex)) {
                            candidate_vertices.push(vertex);
                        }
                    }
                }
                let nb_candidates = candidate_vertices.len();
                if nb_candidates == 0 {
                    continue;
                }
                // Project the vertices on the curve (`cxx:4063-4069`).
                let mut arr_para = vec![0.0f64; nb_candidates];
                let mut arr_dist = vec![-1.0f64; nb_candidates];
                let mut arr_pnt = vec![GpPnt::zero(); nb_candidates];
                let mut arr_proj = vec![GpPnt::zero(); nb_candidates];
                for j in 0..nb_candidates {
                    arr_pnt[j] = BRepTool::vertex_point_world(&Vertex(candidate_vertices[j].clone()));
                }
                project_points_on_curve(
                    &arr_pnt,
                    c3d.as_ref(),
                    first,
                    last,
                    self.tolerance,
                    true,
                    &mut arr_dist,
                    &mut arr_para,
                    &mut arr_proj,
                );
                // Create the cutting nodes (`cxx:4070-4085`).
                let mut seq_node: ShapeSequence<bool> = ShapeSequence::new();
                let mut seq_para = RealSequence::new();
                self.create_cutting_nodes(
                    &candidate_vertices,
                    &bound,
                    &v1,
                    &v2,
                    &arr_dist,
                    &arr_para,
                    &arr_pnt,
                    &mut seq_node,
                    &mut seq_para,
                );
                if seq_para.len() == 0 {
                    continue;
                }
                // Create the cutting sections (`cxx:4086-4087`).
                self.create_sections(&bound, &seq_node, &seq_para, &mut list_sections);
            }
            if list_sections.len() > 1 {
                // Modification of the maps (`cxx:4089-4117`).
                for section in &list_sections {
                    let mut vertices: Vec<TopoShape> = Vec::new();
                    for child in crate::iterator::cumulated_children(section) {
                        vertices.push(child);
                    }
                    for vertex in vertices {
                        let mut vertex = vertex;
                        if let Some(node) = self.vertex_node.get(&vertex) {
                            vertex = node;
                        }
                        self.node_sections
                            .change(&vertex)
                            .push(section.clone());
                    }
                    self.section_bound.bind(section, bound.clone());
                }
                self.bound_sections.bind(&bound, list_sections);
            }
        }
    }

    /// `CreateCuttingNodes(MapVert, bound, vfirst, vlast, arrDist, arrPara,
    /// arrPnt, seqVert, seqPara)` (`cxx:5479-5665`).
    #[allow(clippy::too_many_arguments)]
    fn create_cutting_nodes(
        &mut self,
        map_vert: &[TopoShape],
        bound: &TopoShape,
        vfirst: &TopoShape,
        vlast: &TopoShape,
        arr_dist: &[f64],
        arr_para: &[f64],
        arr_pnt: &[GpPnt],
        seq_vert: &mut ShapeSequence<bool>,
        seq_para: &mut RealSequence,
    ) {
        let nb_proj = map_vert.len();
        // Reorder the projections by distance (`cxx:4913-4937`).
        let mut ordered_index: Vec<usize> = Vec::new();
        let mut ordered_distance: Vec<f64> = Vec::new();
        for i in 1..=nb_proj {
            let dist_proj = arr_dist[i - 1];
            if dist_proj < 0.0 {
                continue;
            }
            let mut is_inserted = false;
            for j in 0..ordered_index.len() {
                if dist_proj < ordered_distance[j] {
                    ordered_index.insert(j, i);
                    ordered_distance.insert(j, dist_proj);
                    is_inserted = true;
                    break;
                }
            }
            if !is_inserted {
                ordered_index.push(i);
                ordered_distance.push(dist_proj);
            }
        }
        let nb_proj = ordered_index.len();
        if nb_proj == 0 {
            return;
        }
        let builder = TopoBuilder::new();
        // Insert the two initial vertices (`cxx:4942-4965`).
        let mut seq_dist: Vec<f64> = Vec::new();
        let mut seq_pnt: Vec<GpPnt> = Vec::new();
        {
            let bound_e = Edge(bound.clone());
            let (first, last) = BRepTool::edge_parameters(&bound_e);
            seq_vert.prepend(vfirst, true);
            seq_vert.append(vlast, true);
            seq_para.prepend(first);
            seq_para.append(last);
            seq_dist.insert(0, -1.0);
            seq_dist.push(-1.0);
            seq_pnt.insert(0, gac_value(bound, first));
            seq_pnt.push(gac_value(bound, last));
        }
        // `TopTools_IndexedDataMapOfShapeShape` with a null replacement.
        let mut node_cutting_vertex: ShapeKeyedMap<Option<TopoShape>> = ShapeKeyedMap::new();
        for i in 1..=nb_proj {
            let index = ordered_index[i - 1];
            let dis_proj = arr_dist[index - 1];
            let pnt_proj = arr_pnt[index - 1];
            // Skip a node already bound to a cutting vertex (`cxx:4975-4978`).
            let Some(node) = self.vertex_node.get(&map_vert[index - 1]) else {
                continue;
            };
            if node_cutting_vertex.contains(&node) {
                continue;
            }
            // Find the closest vertex (`cxx:4980-4994`).
            let mut index_min = 1usize;
            let mut dist_min = pnt_proj.distance(&seq_pnt[0]);
            for j in 2..=seq_pnt.len() {
                let dist = pnt_proj.distance(&seq_pnt[j - 1]);
                if dist < dist_min {
                    dist_min = dist;
                    index_min = j;
                }
            }
            if dist_min <= (dis_proj * 0.1).max(self.min_tolerance) {
                let jdist = seq_dist[index_min - 1];
                if jdist < 0.0 {
                    // Bind a new cutting node (end vertex only).
                    seq_dist[index_min - 1] = dis_proj;
                    let cvertex = seq_vert.value(index_min).unwrap().clone();
                    node_cutting_vertex.bind(&node, Some(cvertex));
                } else {
                    // Bind a secondary cutting node: a null replacement.
                    node_cutting_vertex.bind(&node, None);
                }
            } else {
                // Build a new cutting vertex (`cxx:5012-5031`).
                let cvertex = builder.make_vertex(pnt_proj, CONFUSION);
                node_cutting_vertex.bind(&node, Some(cvertex.0.clone()));
                let par_proj = arr_para[index - 1];
                let mut inserted = false;
                for j in 2..=seq_para.len() {
                    if par_proj <= seq_para.value(j).unwrap() {
                        seq_vert.insert_before(j, &cvertex.0, true);
                        seq_para.insert_before(j, par_proj);
                        seq_dist.insert(j - 1, dis_proj);
                        seq_pnt.insert(j - 1, pnt_proj);
                        inserted = true;
                        break;
                    }
                }
                let _ = inserted;
            }
        }
        // Fill the map for the cutting nodes (`cxx:5034-5068`).
        let entries: Vec<(TopoShape, Option<TopoShape>)> = node_cutting_vertex
            .iter()
            .map(|(k, v)| (k.clone(), v.clone()))
            .collect();
        for (node, cnode) in entries {
            let Some(cnode) = cnode else {
                continue; // Skip secondary nodes.
            };
            let cnode = match self.vertex_node.get(&cnode) {
                Some(end_node) => end_node,
                None => {
                    // Create the link: cutting vertex -> node.
                    self.cutting_node.bind(&cnode, vec![node.clone()]);
                    cnode
                }
            };
            // Create the link: node -> cutting vertex.
            self.cutting_node.change(&node).push(cnode);
        }
        // Remove the two initial vertices (`cxx:5070-5075`).
        seq_vert.remove(1);
        seq_vert.remove(seq_vert.len());
        seq_para.remove(1);
        seq_para.remove(seq_para.len());
    }

    /// `CreateSections(section, seqNode, seqPara, listEdge)`
    /// (`cxx:5668-5878`).
    fn create_sections(
        &mut self,
        section: &TopoShape,
        seq_node: &ShapeSequence<bool>,
        seq_para: &RealSequence,
        list_edge: &mut Vec<TopoShape>,
    ) {
        let sec = Edge(section.clone());
        // To keep the NM vertices on the edge (`cxx:5088-5090`).
        let (mut nm_vert, mut nm_pars) = find_nm_vertices(section);
        let builder = TopoBuilder::new();
        let (first, last) = BRepTool::edge_parameters(&sec);
        let len = seq_para.len() + 1;
        for i in 1..=len {
            let edge = crate::shape_fix_compose_shell::empty_copied_edge(section);
            let (par1, par2, v1, v2);
            if i == 1 {
                par1 = first;
                par2 = seq_para.value(i).unwrap_or(first);
                v1 = first_vertex(&sec);
                v2 = seq_node.value(i).cloned();
            } else if i == len {
                par1 = seq_para.value(i - 1).unwrap_or(first);
                par2 = last;
                v1 = seq_node.value(i - 1).cloned();
                v2 = last_vertex(&sec);
            } else {
                par1 = seq_para.value(i - 1).unwrap_or(first);
                par2 = seq_para.value(i).unwrap_or(last);
                v1 = seq_node.value(i - 1).cloned();
                v2 = seq_node.value(i).cloned();
            }
            let mut an_edge = edge.oriented(Orientation::Forward);
            if let Some(v1) = v1 {
                builder.add(&mut an_edge, &v1.oriented(Orientation::Forward));
            }
            if let Some(v2) = v2 {
                builder.add(&mut an_edge, &v2.oriented(Orientation::Reversed));
            }
            GeometryRegistry::global().set_edge_range(&edge, par1, par2);
            // Keep the NM vertices that fall in this parameter range.
            let mut k = 1usize;
            while k <= nm_pars.len() {
                let apar = nm_pars[k - 1];
                if apar >= par1 && apar <= par2 {
                    builder.add(&mut an_edge, &nm_vert[k - 1]);
                    nm_vert.remove(k - 1);
                    nm_pars.remove(k - 1);
                } else {
                    k += 1;
                }
            }
            list_edge.push(edge);
        }
        let list_faces = match self.bound_faces.get(section) {
            Some(f) => f,
            None => return,
        };
        if list_faces.is_empty() {
            return;
        }
        let tol_edge = BRepTool::edge_tolerance(&sec);
        // Add the cutting pcurves (`cxx:5163-5274`).
        for fac_shape in &list_faces {
            let fac = crate::shape::Face(fac_shape.clone());
            let Some((c2d, _, _)) = crate::boptools_2d::curve_on_surface_range(&sec, &fac) else {
                continue;
            };
            let is_seam = BRepTool::is_closed_edge_face(&sec, &fac);
            let mut c2d1: Option<Arc<dyn Curve2d>> = None;
            if is_seam {
                let sec_rev = sec.oriented(Orientation::Reversed);
                match crate::boptools_2d::curve_on_surface_range(&Edge(sec_rev), &fac) {
                    Some((c, _, _)) => c2d1 = Some(c),
                    None => continue,
                }
            }
            for edge in list_edge.iter() {
                let (par1, par2) = BRepTool::edge_parameters(&Edge(edge.clone()));
                let c2d_new = copy_curve2d(&c2d);
                let c2d1_new = c2d1.as_ref().map(copy_curve2d);
                let _ = (par1, par2);
                if !is_seam && c2d1_new.is_none() {
                    GeometryRegistry::global().set_edge_pcurve(
                        edge,
                        GeometryRegistry::shape_key(&fac.0),
                        c2d_new,
                    );
                    GeometryRegistry::global().set_edge_tolerance(edge, tol_edge);
                } else {
                    let mut ori = edge.orientation();
                    if fac.0.orientation() == Orientation::Reversed {
                        ori = ori.reversed();
                    }
                    let fk = GeometryRegistry::shape_key(&fac.0);
                    let second = c2d1_new.unwrap_or_else(|| c2d_new.clone());
                    if ori == Orientation::Forward {
                        GeometryRegistry::global().set_edge_pcurves(
                            edge,
                            fk,
                            vec![c2d_new, second],
                        );
                    } else {
                        GeometryRegistry::global().set_edge_pcurves(
                            edge,
                            fk,
                            vec![second, c2d_new],
                        );
                    }
                    GeometryRegistry::global().set_edge_tolerance(edge, tol_edge);
                }
            }
        }
    }
}

/// `GeomAdaptor_Curve::Value(u)` on the world curve of `shape`.
fn gac_value(shape: &TopoShape, u: f64) -> GpPnt {
    match GeometryRegistry::global().edge_curve(shape) {
        Some(c) => c.d0(u),
        None => GpPnt::zero(),
    }
}

/// `Geom2d_Curve::Copy()`.
fn copy_curve2d(c: &Arc<dyn Curve2d>) -> Arc<dyn Curve2d> {
    Arc::from(c.clone_dyn())
}

/// `TopExp::FirstVertex(E)` (`TopExp.cxx:118-131`).
fn first_vertex(e: &Edge) -> Option<TopoShape> {
    let (v1, v2) = crate::topo_tools::edge_vertices(e);
    let v1 = v1?.0;
    let v2 = v2.map(|v| v.0);
    if e.0.orientation() == Orientation::Reversed {
        v2
    } else {
        Some(v1)
    }
}

/// `TopExp::LastVertex(E)` (`TopExp.cxx:147-160`).
fn last_vertex(e: &Edge) -> Option<TopoShape> {
    let (v1, v2) = crate::topo_tools::edge_vertices(e);
    let v1 = v1?.0;
    let v2 = v2.map(|v| v.0);
    if e.0.orientation() == Orientation::Reversed {
        Some(v1)
    } else {
        v2
    }
}

#[allow(dead_code)]
fn unused(_: &Geom2dTrimmedCurve) {}


