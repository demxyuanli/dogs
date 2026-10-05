//! Port of ShapeFix_ComposeShell::SplitWire (ShapeFix_ComposeShell.cxx:942-1429).
//!
//! myLoc is the identity here: this port's shapes carry no location, so the
//! gp_Trsf branch (cxx:959-962), the two Pnt transforms (cxx:1025-1029), the
//! MakeVertex transform (cxx:1240) and the T argument of CheckByCurve3d all
//! collapse to the identity.

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpPnt2d, GpTrsf};
use occt_core::precision::{CONFUSION, PCONFUSION};
use occt_geom::Curve;
use occt_geom2d::curve::Curve2d;

use crate::abs::Orientation;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Vertex};
use crate::tgeometry::GeometryRegistry;

use super::helpers::*;
use super::reshape::apply_context;
use super::reshape::ReShape;
use super::shell::ComposeShell;
use super::wire_segment::WireSegment;

/// ShapeAnalysis_Edge::Curve3d(edge, C, First, Last) with the default
/// orient=true (ShapeAnalysis_Edge.cxx:100-125): the stored range, swapped
/// for a REVERSED edge.
fn curve3d_oriented(edge: &Edge) -> Option<(Arc<dyn Curve>, f64, f64)> {
    let curve = BRepTool::edge_curve(edge)?;
    let (mut first, mut last) = BRepTool::edge_parameters(edge);
    if edge.0.orientation().is_reversed() {
        std::mem::swap(&mut first, &mut last);
    }
    Some((curve, first, last))
}

/// TopoDS_Shape::EmptyCopied for a vertex: a fresh TShape carrying the same
/// point and tolerance (BRep_TVertex::EmptyCopy) and the same orientation
/// (TopoDS_Shape.hxx:294-302). Used by cxx:1266-1282 so that fixing
/// SameParameter later cannot widen the tolerance of the original vertex.
fn empty_copied_vertex(v: &Vertex) -> Vertex {
    let mut nv = Vertex::new();
    nv.set_point(BRepTool::vertex_point(v));
    nv.set_tolerance(BRepTool::vertex_tolerance(v));
    nv.0.set_orientation(v.0.orientation());
    nv
}

impl ComposeShell {
    /// ShapeFix_ComposeShell::SplitWire (ShapeFix_ComposeShell.cxx:942-1429).
    ///
    /// indexes is the 1-based edge index of every split point, values the
    /// matching pcurve parameter and segment_codes the position code of each
    /// segment (all 1-based in OCCT; the slices here keep OCCT's 1-based
    /// meaning). vertices receives the cut vertices in walk order.
    pub fn split_wire(
        &mut self,
        wire: &mut WireSegment,
        indexes: &mut Vec<i32>,
        values: &[f64],
        vertices: &mut Vec<Vertex>,
        segment_codes: &[i32],
        is_cut_by_u: bool,
        cut_index: i32,
    ) -> WireSegment {
        let builder = TopoBuilder::new();
        let reg = GeometryRegistry::global();
        let mut result = WireSegment::new();
        let Some(face) = self.face.clone() else {
            return result;
        };
        let nb_splits = indexes.len();
        let mut start = 1usize;
        let an_wire_orient = wire.orientation();
        // cxx:958-962: T is myLoc.Inverted().Transformation(), the identity here.
        let t = GpTrsf::identity();
        let code_at = |k: usize| -> i32 {
            if k >= 1 && k <= segment_codes.len() {
                segment_codes[k - 1]
            } else {
                0
            }
        };

        let mut i = 1usize;
        while i <= wire.nb_edges() {
            // cxx:969: for an already split seam edge, redistribute its split points.
            let nsplit = apply_context(wire, i, self.context_mut());
            let mut redo = false;
            if nsplit != 1 {
                distribute_split_points(wire.edges(), i, nsplit.max(0) as usize, indexes, values);
                if nsplit <= 0 {
                    // cxx:974-980: edge dismissed; i-- then i++ re-runs the same i.
                    redo = true;
                }
            }
            if !redo {
                let edge = wire.edge(i).cloned();
                let patch = wire.get_patch_index(i).unwrap_or((0, 0, 0, 0));
                if let Some(edge) = edge {
                    let (iumin, iumax, ivmin, ivmax) = patch;
                    // cxx:988: position code of the first segment of the edge.
                    let mut code = code_at(if start > 1 { start - 1 } else { segment_codes.len() });
                    // cxx:991-995.
                    let mut stop = start;
                    while stop <= nb_splits && indexes[stop - 1] == i as i32 {
                        stop += 1;
                    }
                    if stop == start {
                        // cxx:996-1004.
                        result.add_edge_patch(0, edge.clone(), iumin, iumax, ivmin, ivmax);
                        if code != 0 || wire.orientation() != Orientation::External {
                            define_patch(&mut result, code, is_cut_by_u, cut_index, 0);
                        }
                    } else {
                        // cxx:1006-1015: non-manifold vertices of the edge.
                        let a_nm_vertices: Vec<Vertex> = {
                            let ts = edge.0.tshape.read().expect("poisoned TShape lock");
                            ts.children
                                .iter()
                                .filter(|c| c.is_vertex())
                                .filter(|c| {
                                    let o = c.orientation();
                                    o != Orientation::Forward && o != Orientation::Reversed
                                })
                                .map(|c| Vertex(c.clone()))
                                .collect()
                        };
                        let c3d_data = curve3d_oriented(&edge);
                        let (c3d, f3d, l3d) = match &c3d_data {
                            Some((c, f, l)) => (Some(c.clone()), *f, *l),
                            None => (None, 0.0, 0.0),
                        };
                        let pc = crate::boptools_2d::curve_on_surface_oriented(&edge, &face, true);
                        if pc.is_none() {
                            // cxx:1041-1044: sae.PCurve failed -> FAIL2. OCCT then
                            // dereferences a null C2d; the port cannot, so the edge is
                            // appended unsplit (the failure branch is unreachable for the
                            // STEP data this goal targets).
                            self.status |= super::shell::SHAPEEXTEND_FAIL2;
                        }
                        let prev_v0 = crate::shhealing::first_vertex(&edge);
                        let last_v0 = crate::shhealing::last_vertex(&edge);
                        match (prev_v0, last_v0, pc) {
                            (Some(prev_v0), Some(last_v0), Some((c2d, first_par, last_par))) => {
                                let mut prev_v = prev_v0;
                                let mut last_v = last_v0;
                                let mut prev_v_tol = self.limit_tolerance(BRepTool::vertex_tolerance(&prev_v));
                                let last_v_tol = self.limit_tolerance(BRepTool::vertex_tolerance(&last_v));
                                let mut prev_v_pnt = BRepTool::vertex_point(&prev_v);
                                let last_v_pnt = BRepTool::vertex_point(&last_v);
                                let tol_edge = BRepTool::edge_tolerance(&edge);
                                // cxx:1046-1087: sequence of non-manifold parameters.
                                let mut a_nm_vertices = a_nm_vertices;
                                let mut a_nm_params: Vec<f64> = Vec::new();
                                if !a_nm_vertices.is_empty() {
                                    for v in &a_nm_vertices {
                                        let ap_v = BRepTool::vertex_point(v);
                                        let apar = if let Some(c) = &c3d {
                                            crate::shhealing::shape_analysis_curve::project_adaptor(
                                                c.as_ref(), &ap_v, CONFUSION, true,
                                            ).param
                                        } else {
                                            let surf = BRepTool::face_surface(&face)
                                                .expect("face without surface");
                                            let a_p2d = crate::pcurve_full::surface_value_of_uv(
                                                surf.as_ref(), &ap_v, CONFUSION,
                                            );
                                            occt_geom2d::extrema2d::point_curve_extrema2d(
                                                c2d.as_ref(), &a_p2d,
                                            ).u1
                                        };
                                        a_nm_params.push(apar);
                                    }
                                }
                                // cxx:1089-1102.
                                let sp = f3d == first_par && l3d == last_par;
                                let span2d = last_par - first_par;
                                let mut prev_par = first_par;
                                let mut prev_pnt2d = c2d.d0(prev_par);
                                let last_pnt2d = c2d.d0(last_par);
                                let mut prev_pnt = self.grid.value_pnt(&prev_pnt2d);
                                let last_pnt = self.grid.value_pnt(&last_pnt2d);
                                let is_periodic = c2d.is_periodic();
                                let a_period = if is_periodic { c2d.period() } else { 0.0 };

                                // cxx:1104-1108.
                                let nb_edges_start = result.nb_edges();
                                let mut splitted = false;
                                let mut curr_par = last_par;
                                let mut j = start;
                                while j <= stop {
                                    if !splitted && j >= stop {
                                        break;
                                    }
                                    curr_par = if j < stop { values[j - 1] } else { last_par };
                                    // cxx:1116-1128: shift a periodic pcurve parameter into range.
                                    if is_periodic {
                                        let hi = last_par.max(first_par) + PCONFUSION;
                                        let lo = last_par.min(first_par) - PCONFUSION;
                                        if curr_par > hi || curr_par < lo {
                                            let shift = crate::shhealing::adjust_by_period(
                                                curr_par, (first_par + last_par) * 0.5, a_period,
                                            );
                                            curr_par += shift;
                                        }
                                    }
                                    let mut curr_pnt2d = GpPnt2d::new(0.0, 0.0);
                                    let mut curr_pnt = GpPnt::zero();
                                    let mut do_cut = true;
                                    let mut v_opt: Option<Vertex> = None;
                                    if (curr_par - last_par).abs() < PCONFUSION {
                                        // cxx:1136-1140.
                                        v_opt = Some(last_v.clone());
                                        do_cut = false;
                                    } else if (curr_par - prev_par).abs() < PCONFUSION {
                                        // cxx:1141-1146.
                                        vertices.push(prev_v.clone());
                                        code = code_at(j);
                                        prev_par = curr_par;
                                        j += 1;
                                        continue;
                                    } else {
                                        curr_pnt2d = c2d.d0(curr_par);
                                        curr_pnt = self.grid.value_pnt(&curr_pnt2d);
                                        let mid_last = c2d.d0(0.5 * (curr_par + last_par));
                                        let last_ok = curr_pnt.distance(&last_v_pnt) <= last_v_tol
                                            && check_by_curve_3d(
                                                &last_v_pnt,
                                                c3d.as_ref(),
                                                f3d + (curr_par - first_par) * (l3d - f3d) / span2d,
                                                &t,
                                                last_v_tol + 2.0 * CONFUSION,
                                            )
                                            && last_pnt.distance(&self.grid.value_pnt(&mid_last)) <= last_v_tol;
                                        if last_ok {
                                            // cxx:1151-1184.
                                            v_opt = Some(last_v.clone());
                                            let (u_res, v_res) =
                                                self.split_res(is_cut_by_u, cut_index, last_v_tol);
                                            if is_coincided(&last_pnt2d, &curr_pnt2d, u_res, v_res, last_v_tol)
                                                && is_coincided(
                                                    &last_pnt2d, &mid_last, u_res, v_res, last_v_tol,
                                                )
                                            {
                                                do_cut = false;
                                            }
                                        } else {
                                            let mid_prev = c2d.d0(0.5 * (curr_par + prev_par));
                                            let prev_ok = curr_pnt.distance(&prev_v_pnt) <= prev_v_tol
                                                && check_by_curve_3d(
                                                    &prev_v_pnt,
                                                    c3d.as_ref(),
                                                    f3d + (curr_par - first_par) * (l3d - f3d) / span2d,
                                                    &t,
                                                    prev_v_tol + 2.0 * CONFUSION,
                                                )
                                                && prev_pnt.distance(&self.grid.value_pnt(&mid_prev)) <= prev_v_tol;
                                            if prev_ok {
                                                // cxx:1185-1222.
                                                v_opt = Some(prev_v.clone());
                                                let (u_res, v_res) =
                                                    self.split_res(is_cut_by_u, cut_index, prev_v_tol);
                                                if is_coincided(
                                                    &prev_pnt2d, &curr_pnt2d, u_res, v_res, prev_v_tol,
                                                ) && is_coincided(
                                                    &prev_pnt2d, &mid_prev, u_res, v_res, prev_v_tol,
                                                ) {
                                                    vertices.push(prev_v.clone());
                                                    code = code_at(j);
                                                    prev_par = curr_par;
                                                    j += 1;
                                                    continue;
                                                }
                                            } else if BRepTool::is_degenerated(&edge)
                                                && crate::topo_tools_full::is_same(&prev_v.0, &last_v.0)
                                            {
                                                // cxx:1224-1229.
                                                v_opt = Some(prev_v.clone());
                                            }
                                        }
                                    }
                                    // cxx:1231-1235: classification code for this segment.
                                    if j > start {
                                        code = code_at(if j > 1 { j - 1 } else { segment_codes.len() });
                                    }
                                    let v;
                                    match v_opt {
                                        None => {
                                            // cxx:1238-1242.
                                            let nv = builder.make_vertex(curr_pnt, tol_edge);
                                            vertices.push(nv.clone());
                                            v = nv;
                                        }
                                        Some(vv) if !do_cut => {
                                            // cxx:1243-1255: adjusted to the end of the edge.
                                            while j < stop {
                                                vertices.push(last_v.clone());
                                                j += 1;
                                            }
                                            if !splitted {
                                                break;
                                            }
                                            curr_par = last_par;
                                            v = vv;
                                        }
                                        Some(vv) => {
                                            vertices.push(vv.clone());
                                            v = vv;
                                        }
                                    }
                                    let mut v = v;
                                    // cxx:1261-1293: protect the original end vertices.
                                    if !splitted {
                                        let f_v = empty_copied_vertex(&prev_v);
                                        self.context_mut().replace(&prev_v.0, &f_v.0);
                                        let l_v = if crate::topo_tools_full::is_same(&prev_v.0, &last_v.0) {
                                            let mut tmp = f_v.clone();
                                            tmp.0.set_orientation(last_v.0.orientation());
                                            tmp
                                        } else {
                                            let nv = empty_copied_vertex(&last_v);
                                            self.context_mut().replace(&last_v.0, &nv.0);
                                            nv
                                        };
                                        if crate::topo_tools_full::is_same(&v.0, &last_v.0) {
                                            v = l_v.clone();
                                        } else if crate::topo_tools_full::is_same(&v.0, &prev_v.0) {
                                            v = f_v.clone();
                                        }
                                        last_v = l_v;
                                        prev_v = f_v;
                                    }
                                    // cxx:1295-1307.
                                    splitted = true;
                                    prev_v.0.set_orientation(Orientation::Forward);
                                    v.0.set_orientation(Orientation::Reversed);
                                    let mut an_init_edge = edge.clone();
                                    let is_manifold = matches!(
                                        edge.0.orientation(),
                                        Orientation::Forward | Orientation::Reversed
                                    );
                                    if !is_manifold {
                                        an_init_edge.0.set_orientation(Orientation::Forward);
                                    }
                                    let mut new_edge = crate::shhealing::copy_replace_vertices_with(
                                        &an_init_edge, Some(&prev_v), Some(&v),
                                    );
                                    // cxx:1309-1337: internal vertices of the edge.
                                    let mut n = 0usize;
                                    while n < a_nm_params.len() {
                                        let apar = a_nm_params[n];
                                        let atmp_v = self.context_mut().apply(&a_nm_vertices[n].0);
                                        let mut removed = false;
                                        if (apar - prev_par).abs() <= PCONFUSION {
                                            self.context_mut().replace(&atmp_v, &prev_v.0);
                                            a_nm_params.remove(n);
                                            a_nm_vertices.remove(n);
                                            removed = true;
                                        } else if (apar - curr_par).abs() <= PCONFUSION {
                                            self.context_mut().replace(&atmp_v, &v.0);
                                            a_nm_params.remove(n);
                                            a_nm_vertices.remove(n);
                                            removed = true;
                                        }
                                        if apar > prev_par && apar < curr_par {
                                            builder.add(&mut new_edge.0, &atmp_v);
                                            if n < a_nm_params.len() {
                                                a_nm_params.remove(n);
                                                a_nm_vertices.remove(n);
                                                removed = true;
                                            }
                                        }
                                        if !removed {
                                            n += 1;
                                        }
                                    }
                                    // cxx:1339-1344.
                                    crate::shhealing::copy_pcurves(&new_edge, &an_init_edge);
                                    let mut tool = crate::shhealing::transfer_params::TransferParametersProj::new();
                                    tool.set_max_tolerance(self.max_tolerance());
                                    tool.init(&an_init_edge, &face);
                                    tool.transfer_range(&mut new_edge, prev_par, curr_par, true);
                                    // cxx:1346-1356.
                                    if !is_manifold {
                                        if code == IOR_UNDEF {
                                            new_edge.0.set_orientation(Orientation::External);
                                        } else {
                                            new_edge.0.set_orientation(edge.0.orientation());
                                        }
                                    }
                                    // cxx:1358-1361.
                                    if !sp && !BRepTool::is_degenerated(&new_edge) {
                                        reg.set_same_range(&new_edge.0, false);
                                    }
                                    // cxx:1362-1366.
                                    if code == 0 && wire.orientation() == Orientation::External {
                                        code = if is_cut_by_u == (j == 1) { 1 } else { 2 };
                                    }
                                    result.add_edge_patch(0, new_edge, iumin, iumax, ivmin, ivmax);
                                    define_patch(&mut result, code, is_cut_by_u, cut_index, 0);
                                    // cxx:1371-1376.
                                    prev_v = v.clone();
                                    prev_v_tol = self.limit_tolerance(BRepTool::vertex_tolerance(&v));
                                    prev_v_pnt = BRepTool::vertex_point(&v);
                                    prev_pnt = curr_pnt;
                                    prev_pnt2d = curr_pnt2d;
                                    // cxx:1108 for-increment.
                                    prev_par = curr_par;
                                    j += 1;
                                }
                                // cxx:1380-1398: record the replacement in the context.
                                if splitted {
                                    let n = result.nb_edges();
                                    let mut wire_edges: Vec<Edge> = Vec::new();
                                    for k in nb_edges_start..n {
                                        let one_based = if matches!(
                                            edge.0.orientation(),
                                            Orientation::Forward | Orientation::Internal
                                        ) {
                                            k + 1
                                        } else {
                                            n - k + nb_edges_start
                                        };
                                        if let Some(e) = result.edge(one_based) {
                                            wire_edges.push(e.clone());
                                        }
                                    }
                                    let res_wire = builder.make_wire(&wire_edges);
                                    self.context_mut().replace(&edge.0, &res_wire.0);
                                } else {
                                    // cxx:1399-1425.
                                    if an_wire_orient == Orientation::Internal && code == 0 {
                                        let mut e_edge = edge.clone();
                                        if e_edge.0.orientation() == Orientation::Internal {
                                            e_edge.0.set_orientation(Orientation::Forward);
                                        }
                                        let mut e1 = crate::shhealing::copy_replace_vertices(&e_edge);
                                        let face_key = GeometryRegistry::shape_key(&face.0);
                                        let c2d2: Arc<dyn Curve2d> = Arc::from(c2d.clone_dyn());
                                        reg.set_edge_pcurves(&e1.0, face_key, vec![c2d.clone(), c2d2]);
                                        e1.0.set_orientation(Orientation::External);
                                        self.context_mut().replace(&edge.0, &e1.0);
                                        result.add_edge_patch(0, e1, iumin, iumax, ivmin, ivmax);
                                    } else {
                                        result.add_edge_patch(0, edge.clone(), iumin, iumax, ivmin, ivmax);
                                    }
                                    if code == 0 && wire.orientation() == Orientation::External {
                                        code = if is_cut_by_u
                                            == ((first_par - curr_par).abs() < (last_par - curr_par).abs())
                                        {
                                            2
                                        } else {
                                            1
                                        };
                                    }
                                    define_patch(&mut result, code, is_cut_by_u, cut_index, 0);
                                }
                            }
                            _ => {
                                // cxx:996-1004 fallback: no endpoints/pcurve means the
                                // edge cannot be split. OCCT would dereference null; the
                                // port appends it unsplit.
                                result.add_edge_patch(0, edge.clone(), iumin, iumax, ivmin, ivmax);
                                define_patch(&mut result, code, is_cut_by_u, cut_index, 0);
                            }
                        }
                    }
                    start = stop;
                }
            }
            if !redo {
                i += 1;
            }
        }
        // cxx:1427.
        result.set_orientation(an_wire_orient);
        result
    }

    /// cxx:1163-1174 / 1199-1210: the u/v resolution used by IsCoincided.
    fn split_res(&self, is_cut_by_u: bool, cut_index: i32, vtol: f64) -> (f64, f64) {
        let (mut u_res, mut v_res) = (self.u_resolution, self.v_resolution);
        if is_cut_by_u {
            let grid_res = get_grid_resolution(self.grid.u_joint_values(), cut_index) / vtol;
            u_res = self.u_resolution.min(grid_res);
        } else {
            let grid_res = get_grid_resolution(self.grid.v_joint_values(), cut_index) / vtol;
            v_res = self.v_resolution.min(grid_res);
        }
        (u_res, v_res)
    }
}
