//! `TopOpeBRepBuild_Builder` — the `Make*` family and the two wire repair
//! helpers of `TopOpeBRepBuild_Merge.cxx`.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_Merge.cxx:48-170`
//! (`CorrectEdgeOrientation`, `CorrectUnclosedWire`), `:374-522`
//! (`MakeSolids`, `MakeShells`, `MakeFaces`, `MakeEdges`) (OCCT 8.0.0).
//!
//! The `Make*` members turn one `AreaBuilder` iteration into the new
//! sub-shapes of the result: `MakeEdges` walks `EdgeBuilder` vertices,
//! `MakeFaces` walks `FaceBuilder` wires and edges, `MakeSolids` /
//! `MakeShells` walk `SolidBuilder` shells and faces. Every shape goes through
//! `TopOpeBRepDS_BuildTool`, so the port keeps the `BRep_Builder` behaviour in
//! one place.
//!
//! Deviation, both wire helpers: `BRepCheck_Analyzer(newWire, false)`
//! (`Merge.cxx:489`) disables the geometric controls, while
//! [`crate::brep_check::Analyzer`] has no such flag yet, so the status list is
//! the one of `GeomControls = true`. Only the `IsValid()` gate and the
//! `BRepCheck_BadOrientationOfSubshape` test read that list, and both are
//! topological (`BRepCheck_Wire`), so the repaired wires are the same.

use crate::abs::{Orientation, ShapeType};
use crate::bop_occt_util::iter_children;
use crate::brep_check::{Analyzer, Status};
use crate::builder::TopoBuilder;
use crate::shape::{Edge, TopoShape};
use crate::topo_tools::{
    edge_vertices, explorer_shapes, is_equal, is_same, map_shapes_and_ancestors,
};

use super::super::edge_builder::EdgeBuilder;
use super::super::face_builder::FaceBuilder;
use super::super::solid_builder::SolidBuilder;
use super::Builder;

impl Builder {
    /// `MakeSolids(SOBU, L)` (`Merge.cxx:374-406`).
    pub fn make_solids(&mut self, sobu: &mut SolidBuilder, l: &mut Vec<TopoShape>) {
        sobu.init_solid();
        while sobu.more_solid() {
            let mut new_solid = TopoShape::new(ShapeType::Solid);
            self.my_build_tool.make_solid(&mut new_solid);
            sobu.init_shell();
            while sobu.more_shell() {
                let new_shell = if sobu.is_old_shell() {
                    // cxx:386-389
                    sobu.old_shell()
                        .unwrap_or_else(|| TopoShape::new(ShapeType::Shell))
                } else {
                    // cxx:390-397
                    let mut s = TopoShape::new(ShapeType::Shell);
                    self.my_build_tool.make_shell(&mut s);
                    sobu.init_face();
                    while sobu.more_face() {
                        if let Some(f) = sobu.face().cloned() {
                            self.my_build_tool.add_shell_face(&mut s, &f);
                        }
                        sobu.next_face();
                    }
                    s
                };
                // cxx:399: `myBuildTool.Closed(newShell, true)` (the exact
                // closed character of the shell is not checked, NYI in OCCT).
                self.my_build_tool.closed(&new_shell, true);
                self.my_build_tool.add_solid_shell(&mut new_solid, &new_shell);
                sobu.next_shell();
            }
            l.push(new_solid);
            sobu.next_solid();
        }
    }

    /// `MakeShells(SOBU, L)` (`Merge.cxx:409-432`).
    pub fn make_shells(&mut self, sobu: &mut SolidBuilder, l: &mut Vec<TopoShape>) {
        sobu.init_shell();
        while sobu.more_shell() {
            let new_shell = if sobu.is_old_shell() {
                // cxx:416-419
                sobu.old_shell()
                    .unwrap_or_else(|| TopoShape::new(ShapeType::Shell))
            } else {
                // cxx:420-427
                let mut s = TopoShape::new(ShapeType::Shell);
                self.my_build_tool.make_shell(&mut s);
                sobu.init_face();
                while sobu.more_face() {
                    if let Some(f) = sobu.face().cloned() {
                        self.my_build_tool.add_shell_face(&mut s, &f);
                    }
                    sobu.next_face();
                }
                s
            };
            l.push(new_shell);
            sobu.next_shell();
        }
    }

    /// `MakeFaces(aFace, FABU, L)` (`Merge.cxx:435-513`).
    pub fn make_faces(
        &mut self,
        a_face: &TopoShape,
        fabu: &mut FaceBuilder,
        l: &mut Vec<TopoShape>,
    ) {
        let hashds = self.data_structure().is_some();

        fabu.init_face();
        while fabu.more_face() {
            // cxx:454
            let mut new_face = self.my_build_tool.copy_face(a_face);

            // cxx:455-464: a new surface replaces the one copied from `aFace`.
            let mut hns = false;
            let mut new_surface = None;
            if hashds {
                if let Some(hds) = self.data_structure() {
                    let bds = hds.ds();
                    hns = bds.has_new_surface(a_face);
                    if hns {
                        new_surface = bds.new_surface(a_face).cloned();
                    }
                }
            }
            if let Some(su) = new_surface {
                self.my_build_tool.update_surface(&new_face, su);
            }

            fabu.init_wire();
            while fabu.more_wire() {
                let is_old = fabu.is_old_wire();
                let mut new_wire = if is_old {
                    // cxx:469-473
                    fabu.old_wire()
                        .unwrap_or_else(|| TopoShape::new(ShapeType::Wire))
                } else {
                    // cxx:474-485
                    let mut w = TopoShape::new(ShapeType::Wire);
                    self.my_build_tool.make_wire(&mut w);
                    fabu.init_edge();
                    while fabu.more_edge() {
                        let e = fabu.edge().clone();
                        if hns {
                            self.my_build_tool.update_surface_edge(&e, a_face, &new_face);
                        }
                        self.my_build_tool.add_wire_edge(&mut w, &e);
                        fabu.next_edge();
                    }
                    w
                };

                // cxx:487-505
                if !is_old {
                    let bca = Analyzer::new(new_wire.clone());
                    if !bca.is_valid() {
                        new_wire.set_free(true);
                        correct_unclosed_wire(&mut new_wire);
                        // `bcr->Status()` is the list recorded before the
                        // repair, which is what OCCT iterates (`cxx:495-503`).
                        let bad_orientation = bca
                            .result(&new_wire)
                            .map(|r| {
                                r.statuses
                                    .iter()
                                    .any(|s| *s == Status::BadOrientationOfSubshape)
                            })
                            .unwrap_or(false);
                        if bad_orientation {
                            correct_edge_orientation(&mut new_wire);
                        }
                    }
                }

                // cxx:504: `myBuildTool.Closed(newWire, true)`.
                self.my_build_tool.closed(&new_wire, true);
                self.my_build_tool.add_face_wire(&mut new_face, &new_wire);
                fabu.next_wire();
            }

            l.push(new_face);
            fabu.next_face();
        }
    }

    /// `MakeEdges(anEdge, EDBU, L)` (`Merge.cxx:516-624`).
    pub fn make_edges(
        &mut self,
        an_edge: &TopoShape,
        edbu: &mut EdgeBuilder,
        l: &mut Vec<TopoShape>,
    ) {
        edbu.init_edge();
        while edbu.more_edge() {
            // cxx:530-538: an edge carrying at most one vertex is dropped.
            let mut nloop = 0;
            edbu.init_vertex();
            while edbu.more_vertex() {
                nloop += 1;
                edbu.next_vertex();
            }
            if nloop <= 1 {
                edbu.next_edge();
                continue;
            }

            // cxx:541
            let mut new_edge = self.my_build_tool.copy_edge(an_edge);

            let mut has_vertex = false;
            edbu.init_vertex();
            while edbu.more_vertex() {
                let Some(mut v) = edbu.vertex() else {
                    edbu.next_vertex();
                    continue;
                };
                let vori = v.orientation();

                // cxx:549-555: a same-domain vertex is replaced by its
                // reference vertex, with the local orientation kept.
                let hassd = self
                    .data_structure()
                    .map(|hds| hds.has_same_domain(&v, false))
                    .unwrap_or(false);
                if hassd {
                    if let Some(hds) = self.data_structure() {
                        let iref = hds.same_domain_reference(&v);
                        if let Some(vref) = hds.shape(iref, false) {
                            v = vref;
                        }
                    }
                    v.set_orientation(vori);
                }

                let ori_v = v.orientation();
                if ori_v != Orientation::External {
                    // cxx:560-591 "betonnage": look for an equivalent vertex
                    // already on the new edge.
                    let mut equa_found = false;
                    for ve in explorer_shapes(&new_edge, ShapeType::Vertex) {
                        let ori_ve = ve.orientation();
                        if is_equal(&v, &ve) {
                            equa_found = true;
                            break;
                        } else if ori_ve == Orientation::Forward || ori_ve == Orientation::Reversed
                        {
                            if ori_v == ori_ve {
                                equa_found = true;
                                break;
                            }
                        } else if ori_ve == Orientation::Internal
                            || ori_ve == Orientation::External
                        {
                            let par_v = edbu.parameter();
                            let par_ve = crate::brep_tool::BRepTool::parameter_of_vertex(
                                &new_edge, &ve,
                            );
                            if par_v == par_ve {
                                equa_found = true;
                                break;
                            }
                        }
                    }
                    if !equa_found {
                        // cxx:592-598
                        has_vertex = true;
                        let par_v = edbu.parameter();
                        self.my_build_tool.add_edge_vertex(&mut new_edge, &v);
                        self.my_build_tool.parameter(&new_edge, &v, par_v);
                    }
                }
                edbu.next_vertex();
            }

            if has_vertex {
                l.push(new_edge);
            }
            edbu.next_edge();
        }
    }
}

/// `CorrectEdgeOrientation(aWire)` (`Merge.cxx:48-136`).
///
/// Re-orders — and where needed reverses — the edges of a wire so that each
/// edge continues the previous one. Edges that connect at neither end are
/// retried after the ones already placed; what stays unmatched is appended to
/// the result unchanged (`cxx:123-127`).
pub(crate) fn correct_edge_orientation(a_wire: &mut TopoShape) {
    let mut an_edge_list: Vec<TopoShape> = iter_children(a_wire);
    let mut n = an_edge_list.len();
    if n <= 1 {
        return;
    }

    let mut an_aux_list: Vec<TopoShape> = Vec::new();
    let mut a_true_edge_list: Vec<TopoShape> = Vec::new();

    let mut it = 0usize;
    let mut an_cur_edge = an_edge_list[0].clone();
    let (mut vf, mut vl) = edge_end_vertices(&an_cur_edge);
    a_true_edge_list.push(an_cur_edge.clone());
    it += 1;

    let mut append = true;
    while n > 0 && append {
        while it < an_edge_list.len() {
            append = false;
            an_cur_edge = an_edge_list[it].clone();
            let (v1f, v1l) = edge_end_vertices(&an_cur_edge);
            if same_opt(&v1f, &vl) {
                a_true_edge_list.push(an_cur_edge.clone());
                vl = v1l;
                append = true;
                it += 1;
                continue;
            }
            if same_opt(&v1l, &vf) {
                a_true_edge_list.push(an_cur_edge.clone());
                vf = v1f;
                append = true;
                it += 1;
                continue;
            }
            if same_opt(&v1l, &vl) {
                let mut an_rev_edge = an_cur_edge.clone();
                an_rev_edge.reverse();
                a_true_edge_list.push(an_rev_edge);
                vl = v1f;
                append = true;
                it += 1;
                continue;
            }
            if same_opt(&v1f, &vf) {
                let mut an_rev_edge = an_cur_edge.clone();
                an_rev_edge.reverse();
                a_true_edge_list.push(an_rev_edge);
                vf = v1l;
                append = true;
                it += 1;
                continue;
            }
            an_aux_list.push(an_cur_edge.clone());
            it += 1;
        }

        // cxx:117-121
        an_edge_list = std::mem::take(&mut an_aux_list);
        n = an_edge_list.len();
        it = 0;
    }

    if n > 0 {
        a_true_edge_list.extend(an_edge_list);
    }

    // cxx:129-135: `aWire.Nullify(); BB.MakeWire(...)` then add the edges.
    *a_wire = TopoShape::new(ShapeType::Wire);
    let builder = TopoBuilder::new();
    for e in a_true_edge_list {
        builder.add(a_wire, &e);
    }
}

/// `CorrectUnclosedWire(aWire)` (`Merge.cxx:138-171`).
///
/// Drops the edges that carry at most one vertex (`cxx:144-153`) and the edges
/// one of whose vertices is not shared with any other edge of the wire
/// (`cxx:155-169`).
pub(crate) fn correct_unclosed_wire(a_wire: &mut TopoShape) {
    let builder = TopoBuilder::new();

    for ed in iter_children(a_wire) {
        // `ed.NbChildren()`: the direct sub-shapes of the edge (its vertices).
        let nbv = ed.tshape.read().map(|t| t.children.len()).unwrap_or(0);
        if nbv <= 1 {
            builder.remove(a_wire, &ed);
        }
    }

    // `TopExp::MapShapesAndAncestors(aWire, VERTEX, EDGE, VElists)`.
    let vlists = map_shapes_and_ancestors(a_wire, ShapeType::Vertex, ShapeType::Edge);
    for elist in vlists.values() {
        if elist.len() == 1 {
            builder.remove(a_wire, &elist[0]);
        }
    }
}

/// `TopExp::Vertices(E, vf, vl, true)` (`TopExp.cxx:166-190`): the two
/// bounding vertices of `E`.
fn edge_end_vertices(e: &TopoShape) -> (Option<TopoShape>, Option<TopoShape>) {
    let Some(edge) = Edge::wrap(e.clone()) else {
        return (None, None);
    };
    let (a, b) = edge_vertices(&edge);
    (a.map(|v| v.0), b.map(|v| v.0))
}

/// `TopoDS_Shape::IsSame` on possibly-null shapes: two null shapes are the
/// same shape in OCCT (`TopoDS_Shape.lxx:43-52`).
fn same_opt(a: &Option<TopoShape>, b: &Option<TopoShape>) -> bool {
    match (a, b) {
        (None, None) => true,
        (Some(a), Some(b)) => is_same(a, b),
        _ => false,
    }
}

