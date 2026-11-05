//! The edge sections of `BRepFill_Sweep::BuildShell` (`cxx:2665-2989`): the
//! singular cases, the iso-u edges (`cxx:2845-2911`), the iso-v edges
//! (`cxx:2913-2986`) and the `RebuildTopOrBottomEdge` helper
//! (`cxx:4047-4151`).
//!
//! Translation notes:
//! - `BuildEdge(S, isUiso, val, VF, VL, Tol)` is
//!   [`crate::brep_fill_sweep::build_edge::build_edge_iso`] and `UpdateEdge(E, S, isUiso, val)`
//!   is [`crate::brep_fill_sweep::update_edge_iso::update_edge_iso`].
//! - `RebuildTopOrBottomEdge(aNewEdge, anEdge, ReversedEdges)` takes `anEdge`
//!   by reference and reverses it in place; the port works on a copy of the
//!   table cell and writes it back through
//!   [`ShellBuild::rebuild_top_or_bottom_edge`].
//! - A null `TopoDS_Shape` is `None`, so `TopoDS::Edge(UEdge(...))` on a null
//!   cell is reported instead of dereferenced.

use std::sync::Arc;

use occt_core::precision::CONFUSION;
use occt_core::toploc::TopLocLocation;
use occt_geom::curve::Curve;

use crate::abs::Orientation;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, TopoShape, Vertex};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{is_same, shape_map_key};

use super::{oriented_edge_vertices, ShellBuild};
use crate::brep_fill_sweep::build_edge::build_edge_iso;
use crate::brep_fill_sweep::update_edge_iso::update_edge_iso;

/// `TopoDS_Shape::IsSame` with the null-shape arm of OCCT where two null shapes
/// are the same.
fn same_shape(a: &Option<TopoShape>, b: &Option<TopoShape>) -> bool {
    match (a, b) {
        (Some(x), Some(y)) => is_same(x, y),
        (None, None) => true,
        _ => false,
    }
}

impl ShellBuild<'_> {
    /// `Vertex(row, col)` as a `TopoDS_Vertex`.
    fn vertex_cell(&self, row: usize, col: usize) -> Option<Vertex> {
        self.tables
            .vertex
            .value(row, col)
            .filter(|s| s.shape_type() == crate::abs::ShapeType::Vertex)
            .map(|s| Vertex(s.clone()))
    }

    /// `TopoDS::Edge(UEdge(row, col))`.
    fn u_edge_cell(&self, row: usize, col: usize) -> Option<Edge> {
        self.tables.u_edge.value(row, col).map(|s| Edge(s.clone()))
    }

    /// `TopoDS::Edge(VEdge(row, col))`.
    fn v_edge_cell(&self, row: usize, col: usize) -> Option<Edge> {
        self.tables.v_edge.value(row, col).map(|s| Edge(s.clone()))
    }

    /// (2) `cxx:2665-2989`: the vertices, the singular cases and the iso edges
    /// of every `(isec, ipath)` cell. `Ok(false)` is the reference's `return
    /// false` for a failed sweep.
    pub(super) fn construct_vertices_and_edges(&mut self) -> Result<bool, String> {
        for ipath in 1..=self.nb_path {
            for isec in 1..=self.nb_law {
                if self.locals.is_built[isec] {
                    continue;
                }
                // (2.1) `cxx:2692-2798`.
                self.construct_vertices(ipath, isec)?;

                let Some(s) = self.tables.tab_s.value(isec, ipath).cloned() else {
                    return Err(format!(
                        "BuildShell: null sweep surface ({isec}, {ipath}) at BRepFill_Sweep.cxx:2676"
                    ));
                };
                let exuv = self.tables.exch_uv.value(isec, ipath);
                let (mut u_first, mut u_last) = s.u_range();
                let (mut v_first, mut v_last) = s.v_range();
                if self.tables.u_reverse.value(isec, ipath) {
                    if exuv {
                        std::mem::swap(&mut v_first, &mut v_last);
                    } else {
                        std::mem::swap(&mut u_first, &mut u_last);
                    }
                }

                // Singular cases (`cxx:2800-2814`).
                let singv = self.merge_vertex_cells((isec, ipath + 1), (isec + 1, ipath + 1));
                let singu = self.merge_vertex_cells((isec + 1, ipath), (isec + 1, ipath + 1));
                if singu || singv {
                    let tol = self
                        .sweep
                        .my_tol_3d
                        .max(self.tables.tab_err.value(isec, ipath));
                    self.tables
                        .degenerated
                        .set_value(isec, ipath, crate::brep_fill_sweep::helpers::is_degen(&s, tol));
                }

                if self.tables.degenerated.value(isec, ipath) {
                    // `cxx:2815-2839`: particular construction of edges.
                    self.hasdegen = true;
                    if self.tables.u_edge.is_null(isec + 1, ipath) {
                        if singu {
                            if let Some(cell) = self.tables.vertex.value(isec + 1, ipath).cloned() {
                                let mut v = Vertex(cell);
                                let e = crate::brep_fill_sweep::helpers::null_edge(&mut v);
                                self.tables.u_edge.set_value(isec + 1, ipath, Some(e.0));
                            }
                        } else {
                            let previous = self.tables.u_edge.value(isec, ipath).cloned();
                            self.tables.u_edge.set_value(isec + 1, ipath, previous);
                        }
                    }
                    if self.tables.v_edge.is_null(isec, ipath + 1) {
                        if singv {
                            if let Some(cell) = self.tables.vertex.value(isec, ipath + 1).cloned() {
                                let mut v = Vertex(cell);
                                let e = crate::brep_fill_sweep::helpers::null_edge(&mut v);
                                self.tables.v_edge.set_value(isec, ipath + 1, Some(e.0));
                            }
                        } else {
                            let previous = self.tables.v_edge.value(isec, ipath).cloned();
                            self.tables.v_edge.set_value(isec, ipath + 1, previous);
                        }
                    }
                } else {
                    // `cxx:2842-2849`: swap the iso directions on an exchanged
                    // sweep.
                    if exuv {
                        std::mem::swap(&mut u_first, &mut v_first);
                        std::mem::swap(&mut u_last, &mut v_last);
                    }

                    // (2.2) Iso-u (`cxx:2851-2911`).
                    if isec == 1 && self.tables.u_edge.is_null(1, ipath) {
                        let v_at = self.tables.vertex.value(1, ipath).cloned();
                        let v_next = self.tables.vertex.value(1, ipath + 1).cloned();
                        if !same_shape(&v_at, &v_next) {
                            // `cxx:2859-2866`.
                            if let (Some(a), Some(b)) = (&v_at, &v_next) {
                                if a.shape_type() == crate::abs::ShapeType::Vertex
                                    && b.shape_type() == crate::abs::ShapeType::Vertex
                                {
                                    let p1 = BRepTool::vertex_point(&Vertex(a.clone()));
                                    let p2 = BRepTool::vertex_point(&Vertex(b.clone()));
                                    if p1.distance(&p2) <= self.sweep.my_tol_3d {
                                        self.tables.vertex.set_value(1, ipath + 1, Some(a.clone()));
                                    }
                                }
                            }
                        }
                        let vf = self
                            .vertex_cell(1, ipath)
                            .ok_or_else(|| "BuildShell: null vertex (1, ipath)".to_string())?;
                        let vl = self
                            .vertex_cell(1, ipath + 1)
                            .ok_or_else(|| "BuildShell: null vertex (1, ipath+1)".to_string())?;
                        let e = build_edge_iso(&s, !exuv, u_first, &vf, &vl, self.sweep.my_tol_3d)?;
                        self.tables.u_edge.set_value(1, ipath, Some(e.0));
                    } else {
                        let Some(existing) = self.u_edge_cell(isec, ipath) else {
                            // `cxx:2872-2875`: the sweep failed.
                            return Ok(false);
                        };
                        let mut e = existing;
                        update_edge_iso(&mut e, &s, !exuv, u_first)?;
                    }

                    // `cxx:2879-2911`.
                    if self.u_close && isec == self.nb_law {
                        if self.tables.u_edge.is_null(1, ipath) {
                            // Degenerated case.
                            let vf = self.vertex_cell(isec + 1, ipath).ok_or_else(|| {
                                "BuildShell: null vertex (isec+1, ipath)".to_string()
                            })?;
                            let vl = self.vertex_cell(isec + 1, ipath + 1).ok_or_else(|| {
                                "BuildShell: null vertex (isec+1, ipath+1)".to_string()
                            })?;
                            let e =
                                build_edge_iso(&s, !exuv, u_last, &vf, &vl, self.sweep.my_tol_3d)?;
                            self.tables.u_edge.set_value(isec + 1, ipath, Some(e.0));
                        } else if let Some(shared) = self.u_edge_cell(1, ipath) {
                            let mut e = shared.clone();
                            update_edge_iso(&mut e, &s, !exuv, u_last)?;
                            self.tables.u_edge.set_value(isec + 1, ipath, Some(shared.0));
                        }
                    } else if self.tables.u_edge.is_null(isec + 1, ipath) {
                        let vf = self.vertex_cell(isec + 1, ipath).ok_or_else(|| {
                            "BuildShell: null vertex (isec+1, ipath)".to_string()
                        })?;
                        let vl = self.vertex_cell(isec + 1, ipath + 1).ok_or_else(|| {
                            "BuildShell: null vertex (isec+1, ipath+1)".to_string()
                        })?;
                        let e = build_edge_iso(&s, !exuv, u_last, &vf, &vl, self.sweep.my_tol_3d)?;
                        self.tables.u_edge.set_value(isec + 1, ipath, Some(e.0));
                    } else if let Some(existing) = self.u_edge_cell(isec + 1, ipath) {
                        let mut e = existing;
                        update_edge_iso(&mut e, &s, !exuv, u_last)?;
                    }

                    // (2.3) Iso-v (`cxx:2913-2986`).
                    if ipath == 1 {
                        let vf = self
                            .vertex_cell(isec, 1)
                            .ok_or_else(|| "BuildShell: null vertex (isec, 1)".to_string())?;
                        let vl = self
                            .vertex_cell(isec + 1, 1)
                            .ok_or_else(|| "BuildShell: null vertex (isec+1, 1)".to_string())?;
                        let new_first =
                            build_edge_iso(&s, exuv, v_first, &vf, &vl, self.sweep.my_tol_3d)?;
                        if self.tables.v_edge.is_null(isec, ipath) {
                            self.tables.v_edge.set_value(isec, ipath, Some(new_first.0));
                        } else if let Some(mut old) = self.v_edge_cell(isec, ipath) {
                            // `cxx:2919-2929`: rebuild the first edge.
                            self.rebuild_top_or_bottom_edge(&new_first, &mut old);
                            self.tables.v_edge.set_value(isec, ipath, Some(old.0.clone()));
                            if self.reversed_edges.contains(&shape_map_key(&old.0)) {
                                if let Some(w) = self.sweep.first_shape.as_mut() {
                                    crate::brep_fill_sweep::helpers::reverse_edge_in_first_or_last_wire(
                                        &mut w.0, &old.0,
                                    );
                                }
                                if let Some(se) = self.locals.start_edges[isec].as_mut() {
                                    se.reverse();
                                }
                            }
                        }
                    } else {
                        let Some(mut existing) = self.v_edge_cell(isec, ipath) else {
                            return Ok(false);
                        };
                        update_edge_iso(&mut existing, &s, exuv, v_first)?;
                    }

                    // `cxx:2938-2986`.
                    if self.v_close && ipath == self.nb_path {
                        if self.tables.v_edge.is_null(isec, 1) {
                            // Degenerated case.
                            let vf = self.vertex_cell(isec, ipath + 1).ok_or_else(|| {
                                "BuildShell: null vertex (isec, ipath+1)".to_string()
                            })?;
                            let vl = self.vertex_cell(isec + 1, ipath + 1).ok_or_else(|| {
                                "BuildShell: null vertex (isec+1, ipath+1)".to_string()
                            })?;
                            let e = build_edge_iso(&s, exuv, v_last, &vf, &vl, self.sweep.my_tol_3d)?;
                            self.tables.v_edge.set_value(isec, ipath + 1, Some(e.0));
                        } else if let Some(shared) = self.v_edge_cell(isec, 1) {
                            let mut e = shared.clone();
                            update_edge_iso(&mut e, &s, exuv, v_last)?;
                            self.tables.v_edge.set_value(isec, ipath + 1, Some(shared.0));
                        }
                    } else if self.tables.v_edge.is_null(isec, ipath + 1) {
                        let vf = self.vertex_cell(isec, ipath + 1).ok_or_else(|| {
                            "BuildShell: null vertex (isec, ipath+1)".to_string()
                        })?;
                        let vl = self.vertex_cell(isec + 1, ipath + 1).ok_or_else(|| {
                            "BuildShell: null vertex (isec+1, ipath+1)".to_string()
                        })?;
                        let e = build_edge_iso(&s, exuv, v_last, &vf, &vl, self.sweep.my_tol_3d)?;
                        self.tables.v_edge.set_value(isec, ipath + 1, Some(e.0));
                    } else if ipath != self.nb_path
                        || self.v_close
                        || (self.global_vclose && self.i_last == self.nb_path_law + 1)
                    {
                        // `cxx:2966-2969`.
                        if let Some(mut existing) = self.v_edge_cell(isec, ipath + 1) {
                            update_edge_iso(&mut existing, &s, exuv, v_last)?;
                        }
                    } else {
                        // `cxx:2971-2985`: rebuild the last edge.
                        let vf = self.vertex_cell(isec, ipath + 1).ok_or_else(|| {
                            "BuildShell: null vertex (isec, ipath+1)".to_string()
                        })?;
                        let vl = self.vertex_cell(isec + 1, ipath + 1).ok_or_else(|| {
                            "BuildShell: null vertex (isec+1, ipath+1)".to_string()
                        })?;
                        let new_last =
                            build_edge_iso(&s, exuv, v_last, &vf, &vl, self.sweep.my_tol_3d)?;
                        if let Some(mut old) = self.v_edge_cell(isec, ipath + 1) {
                            self.rebuild_top_or_bottom_edge(&new_last, &mut old);
                            self.tables.v_edge.set_value(isec, ipath + 1, Some(old.0.clone()));
                            if self.reversed_edges.contains(&shape_map_key(&old.0)) {
                                if let Some(w) = self.sweep.last_shape.as_mut() {
                                    crate::brep_fill_sweep::helpers::reverse_edge_in_first_or_last_wire(
                                        &mut w.0, &old.0,
                                    );
                                }
                            }
                        }
                    }
                }
            }
        }
        Ok(true)
    }

    /// `RebuildTopOrBottomEdge(aNewEdge, anEdge, ReversedEdges)`
    /// (`cxx:4047-4151`): replace the top or bottom section edge by the sweep's
    /// iso curve, reversing it when the two disagree in direction.
    fn rebuild_top_or_bottom_edge(&mut self, a_new_edge: &Edge, an_edge: &mut Edge) {
        let reg = GeometryRegistry::global();
        let a_new_curve: Option<Arc<dyn Curve>> = reg.edge_curve(&a_new_edge.0);
        let (fpar, lpar) = BRepTool::edge_parameters(a_new_edge);

        let mut to_reverse = false;
        let is_degen = BRepTool::is_degenerated(a_new_edge);
        if !is_degen {
            if let Some(new_curve) = &a_new_curve {
                if let (Some(v1), Some(v2)) = oriented_edge_vertices(an_edge) {
                    if !is_same(&v1.0, &v2.0) {
                        let (new_v1, new_v2) = oriented_edge_vertices(a_new_edge);
                        if let (Some(new_v1), Some(new_v2)) = (new_v1, new_v2) {
                            // `V1.Location(Identity)` (`cxx:4069`).
                            let mut v1 = v1;
                            v1.0.set_location(&TopLocLocation::identity());
                            if !is_same(&v1.0, &new_v1.0) {
                                if is_same(&v1.0, &new_v2.0) {
                                    to_reverse = true;
                                } else {
                                    // `gp_Pnt::IsEqual(NewPnt1, TolSum)`.
                                    let pnt1 = BRepTool::vertex_point(&v1);
                                    let new_pnt1 = BRepTool::vertex_point(&new_v1);
                                    let tol_sum = BRepTool::vertex_tolerance(&v1)
                                        + BRepTool::vertex_tolerance(&new_v1);
                                    if !(pnt1.distance(&new_pnt1) <= tol_sum) {
                                        to_reverse = true;
                                    }
                                }
                            }
                        }
                    } else {
                        // `cxx:4089-4100`: a closed section edge, compare the
                        // mid tangents.
                        let (old_first, old_last) = BRepTool::edge_parameters(an_edge);
                        if let Some(old_curve) = reg.edge_curve(&an_edge.0) {
                            let (_, old_d1) = old_curve.d1(0.5 * (old_first + old_last));
                            let (_, new_d1) = new_curve.d1(0.5 * (fpar + lpar));
                            if old_d1.dot(&new_d1) < 0.0 {
                                to_reverse = true;
                            }
                        }
                    }
                }
            }
        }

        an_edge.0.set_location(&TopLocLocation::identity());
        // `TEdge->Tolerance(BRep_Tool::Tolerance(aNewEdge))` (`cxx:4103-4105`).
        reg.set_edge_tolerance(&an_edge.0, BRepTool::edge_tolerance(a_new_edge));
        reg.set_edge_range(&an_edge.0, fpar, lpar);
        if let Some(curve) = &a_new_curve {
            reg.set_edge_curve3d(&an_edge.0, curve.clone(), CONFUSION);
        }
        // `cxx:4110-4124`: copy every CurveOnSurface representation.
        for (surf, pcurve, _f, _l) in reg.edge_pcurve_reps(&a_new_edge.0) {
            reg.update_edge_pcurve_on_surface(&an_edge.0, &surf, Some(pcurve), CONFUSION);
        }

        // `cxx:4126-4148`.
        an_edge.0.set_free(true);
        let (v1, v2) = oriented_edge_vertices(an_edge);
        let mut forward = an_edge.0.oriented(Orientation::Forward);
        forward.set_free(true);
        let b = TopoBuilder::new();
        if let Some(v) = &v1 {
            b.remove(&mut forward, &v.0);
        }
        if let Some(v) = &v2 {
            b.remove(&mut forward, &v.0);
        }
        let mut v1 = v1;
        let mut v2 = v2;
        if let Some(v) = v1.as_mut() {
            v.0.set_location(&TopLocLocation::identity());
        }
        if let Some(v) = v2.as_mut() {
            v.0.set_location(&TopLocLocation::identity());
        }
        if to_reverse {
            if let Some(v) = v2.as_mut() {
                v.0.set_orientation(Orientation::Forward);
            }
            if let Some(v) = v1.as_mut() {
                v.0.set_orientation(Orientation::Reversed);
            }
        }
        if let Some(v) = &v1 {
            b.add(&mut forward, &v.0);
        }
        if let Some(v) = &v2 {
            b.add(&mut forward, &v.0);
        }

        if to_reverse {
            an_edge.0.reverse();
            self.reversed_edges.insert(shape_map_key(&an_edge.0));
        }
        reg.set_degenerated(&an_edge.0, is_degen);
    }
}
