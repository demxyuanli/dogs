//! `perform` (node insertion for a model face) and the structure initialisation.

use std::collections::HashMap;

use occt_core::gp::{GpPnt, GpPnt2d, GpXY};
use occt_core::precision::SQUARE_CONFUSION;
use occt_geom::geom_api::project_point_on_surface;
use occt_geom::Surface;

use crate::abs::Orientation;

use super::super::data_model::*;
use super::super::delaun::Delaun;
use super::super::delaun_data::DelaunDataStructure;
use super::super::delaun_types::{DelaunLink, DelaunTriangle, DelaunVertex, VertexState};
use super::super::face_discret::Classifier;
#[allow(unused_imports)]
use super::super::mesh_algo::*;
use super::super::face_discret::PointState;
use super::super::deflection_control::{
    DelaunayDeflectionControlMeshAlgo, DeflectionMesh,
};
use super::super::parameters::MeshParameters;
use super::super::range_splitter::{
    mesh_surface_type, create_range_splitter, factory_uses_deflection_control,
    grab_params_of_internal_edges, RangeSplitter, MeshSurfaceType,
};
use super::*;
use super::parametric::uv_bounds;

use super::*;

impl DelaunayNodeInsertionMeshAlgo {
    /// Performs node insertion for a model face: extracts the boundary UV points
    /// and constraint edges from the face's wires → edges → pcurves (populated by
    /// `edge_discret`/`face_discret`), registers them together with any internal
    /// face points, triangulates, and returns the triangles.
    ///
    /// ponytail: `MeshModel` has no triangle slot yet — the committing
    /// `collectTriangles`/`commitSurfaceTriangulation` belong to the stub
    /// `mesh_algo.rs` module — so the triangulation is returned and stored in
    /// [`Self::last_result`] for the caller to commit.
    pub fn perform(
        &mut self,
        model: &mut MeshModel,
        face_index: usize,
        params: &MeshParameters,
    ) -> Result<TriangulationResult, String> {
        let (uv, _p3d, _constraints, wires_uv) = Self::collect_boundary_uv(model, face_index)?;
        if crate::debug_flags::mbdiag() && uv.len() <= 1 {
            eprintln!(
                "MBDIAG thin face={face_index} uv={} wires={} edges={:?}",
                uv.len(),
                wires_uv.len(),
                wires_uv.iter().map(|w| w.len()).collect::<Vec<_>>()
            );
        }
        if uv.is_empty() {
            return Err(format!(
                "DelaunayNodeInsertionMeshAlgo::perform: face {face_index} has no boundary UV points"
            ));
        }
        self.boundary_uv = uv.clone();

        // TEMP diagnostic (TORSBD=<x1>,<y1>,<z1>,<x2>,<y2>,<z2>): per-wire-edge
        // pcurve point arrays of the matching face, so they can be diffed against
        // `occt_probe --healerdbg` (the same six bbox coordinates).
        if let Ok(spec) = std::env::var("TORSBD") {
            let t: Vec<f64> = spec
                .split(',')
                .filter_map(|s| s.trim().parse::<f64>().ok())
                .collect();
            if t.len() == 6 {
                let bb = crate::brep_bnd_lib::shape_bnd_box(&model.face(face_index)?.face().0);
                let (mn, mx) = (bb.corner_min(), bb.corner_max());
                {
                    let hit = (mn.x() - t[0]).abs() < 1e-6
                        && (mn.y() - t[1]).abs() < 1e-6
                        && (mn.z() - t[2]).abs() < 1e-6
                        && (mx.x() - t[3]).abs() < 1e-6
                        && (mx.y() - t[4]).abs() < 1e-6
                        && (mx.z() - t[5]).abs() < 1e-6;
                    if hit {                        let face = model.face(face_index)?;
                        for (wi, &wire_index) in face.wires().iter().enumerate() {
                            let wire = model.wire(wire_index)?;
                            for j in 0..wire.edges_nb() {
                                let edge_index = wire.edge(j)?;
                                let orientation = wire.edge_orientation(j)?;
                                let rev = wire.edge_reverse_walk(j);
                                let edge = model.edge(edge_index)?;
                                let Some(pc) = edge.pcurve_for(face_index, orientation) else {
                                    eprintln!("TORSBD face={face_index} w={wi} e={j} pc=none");
                                    continue;
                                };
                                let pts = pc.points();
                                let n = pts.len();
                                let f = pts.first().copied().unwrap_or(GpPnt2d::new(0.0, 0.0));
                                let l = pts.last().copied().unwrap_or(GpPnt2d::new(0.0, 0.0));
                                eprintln!(
                                    "TORSBD face={face_index} w={wi} e={j} nb={n} fwd={} rev={rev} \
same_param={} same_range={} free={} degen={} defl={:.15} ang={:.15} \
first=({:.15},{:.15}) last=({:.15},{:.15})",
                                    pc.is_forward() as u8,
                                    edge.same_param() as u8,
                                    edge.same_range() as u8,
                                    edge.is_free() as u8,
                                    edge.is_degenerated() as u8,
                                    edge.deflection(),
                                    edge.angular_deflection(),
                                    f.x(),
                                    f.y(),
                                    l.x(),
                                    l.y()
                                );
                                let mut s = format!("TORSBDV face={face_index} w={wi} e={j} v=");
                                for p in pts {
                                    s += &format!("{:.15},", p.y());
                                }
                                eprintln!("{s}");
                            }
                        }
                    }
                }
            }
        }

        // `BRepMesh_NodeInsertionMeshAlgo::initDataStructure` (`hxx:47-112`):
        // Reset the splitter, feed the discrete wire points (`collectWirePoints`
        // calls `myRangeSplitter.AddPoint` per pcurve sample, `hxx:179`),
        // AdjustRange, SetCellSize/SetTolerance, then RegisterWire with the
        // splitter range. Without AddPoint a plane's range stays the adaptor's
        // infinite bounds and the cell size collapses the whole face.
        self.splitter = None;
        self.classifier = Classifier::new();
        let mut range_invalid = false;
        {
            let face = model.face(face_index)?;
            if let Some(surface) = face.surface() {
                let mut sp = create_range_splitter(surface.as_ref());
                sp.reset(face, params);
                for w in &wires_uv {
                    for &p in w {
                        sp.add_point(p);
                    }
                }
                sp.adjust_range();
                if sp.is_valid() {
                    let (du, dv) = sp.delta();
                    let (tu, tv) = sp.tolerance_uv();
                    if du.abs() > 1e-16 && dv.abs() > 1e-16 {
                        self.structure
                            .set_cell_size(14.0 * tu / du, 14.0 * tv / dv);
                        self.structure.set_tolerance(tu / du, tv / dv);
                    }
                    for w in &wires_uv {
                        self.classifier
                            .register_wire(w, sp.tolerance_uv(), sp.range_u(), sp.range_v());
                    }
                } else {
                    if crate::debug_flags::mbdiag() {
                        let (su0, su1) = surface.as_ref().u_range();
                        let (sv0, sv1) = surface.as_ref().v_range();
                        let (ru0, ru1) = sp.range_u();
                        let (rv0, rv1) = sp.range_v();
                        let (du, dv) = sp.delta();
                        let (bu0, bu1, bv0, bv1) = uv_bounds(&uv);
                        eprintln!(
                            "MBDIAG rng face={face_index} surf={:?} sup={} svp={} su=[{su0:.6},{su1:.6}] sv=[{sv0:.6},{sv1:.6}] du=[{ru0:.6},{ru1:.6}] dv=[{rv0:.6},{rv1:.6}] delta=({du:.6},{dv:.6}) buv=[{bu0:.6},{bu1:.6}]x[{bv0:.6},{bv1:.6}] nuv={}",
                            mesh_surface_type(surface.as_ref()),
                            surface.as_ref().is_u_periodic(),
                            surface.as_ref().is_v_periodic(),
                            uv.len()
                        );
                    }
                    range_invalid = true;
                }
                if !range_invalid {
                    self.splitter = Some(sp);
                }
            }
        }
        // `BRepMesh_NodeInsertionMeshAlgo::initDataStructure`
        // (`BRepMesh_NodeInsertionMeshAlgo.hxx:79-83`): when `AdjustRange` leaves
        // the range splitter invalid, OCCT marks the face failed and returns
        // false, so `BRepMesh_BaseMeshAlgo::process` (`BRepMesh_BaseMeshAlgo.cxx:52-59`)
        // skips `generateMesh`/`commitSurfaceTriangulation` for it. That is what
        // keeps `BRepMesh_GeomTool::CellsCount` (`BRepMesh_GeomTool.cxx:465-511`)
        // from ever seeing a degenerate discrete range: the cylinder branch
        // divides by the V extent (`BRepMesh_GeomTool.cxx:495-501`), and an
        // invalid splitter implies that extent is zero (a zero-extent range
        // samples to a zero-length `computeLengthV`, hence `IsValid` == false).
        // Without this guard a zero V extent reaches `initCirclesTool`
        // (`BRepMesh_Delaun.cxx:280-304`) as a huge cell count and a microscopic
        // cell size, and binding the circumcircles grows the cell grid without
        // bound.
        if range_invalid {
            model
                .face_mut(face_index)
                .map_err(|e| format!("DelaunayNodeInsertionMeshAlgo::perform: {e}"))?
                .set_status(MeshStatus::FAILURE);
            return Err(format!(
                "DelaunayNodeInsertionMeshAlgo::perform: face {face_index} has an invalid discrete range"
            ));
        }
        if self.classifier.wires_nb() == 0 {
            let (umin, umax, vmin, vmax) = uv_bounds(&uv);
            for w in &wires_uv {
                self.classifier
                    .register_wire(w, (1e-9, 1e-9), (umin, umax), (vmin, vmax));
            }
        }
        self.init_data_structure(model, face_index)?;
        self.diag_dump_structure(face_index);
        if let Ok(spec) = std::env::var("DPROBE") {
            if spec.trim().parse::<i32>().ok() == Some(face_index as i32) {
                let nb = self.structure.nb_nodes() as i32;
                let (mut nfree, mut nfront, mut nfixed) = (0usize, 0usize, 0usize);
                for i in 1..=nb {
                    match self.structure.node(i).state {
                        VertexState::Free => nfree += 1,
                        VertexState::Frontier => nfront += 1,
                        VertexState::Fixed => nfixed += 1,
                        _ => {}
                    }
                }
                eprintln!(
                    "DPROBE face={face_index} nb_nodes={nb} nb_links={} uv={} free={nfree} frontier={nfront} fixed={nfixed}",
                    self.structure.nb_links(),
                    uv.len()
                );
                for i in 1..=nb {
                    let n = self.structure.node(i);
                    let p3d = self
                        .nodes_map
                        .get(i as usize - 1)
                        .copied()
                        .unwrap_or_else(GpPnt::zero);
                    eprintln!(
                        "DPROBE n i={i} st={} uv=({:.6},{:.6}) p=({:.6},{:.6},{:.6})",
                        n.state.to_str(),
                        n.location.x(),
                        n.location.y(),
                        p3d.x(),
                        p3d.y(),
                        p3d.z()
                    );
                }
            }
        }

        // Internal (in-face) 3D points -> UV via the face surface, dropped when
        // outside the face (`insertInternalVertex` classifier check).
        let mut internal_uv: Vec<GpPnt2d> = Vec::new();
        let face = model.face(face_index)?;
        if let Some(surface) = face.surface() {
            for p in face.points() {
                if let Some(proj) = project_point_on_surface(surface.as_ref(), p, 1e-7) {
                    let uv_p = GpPnt2d::new(proj.u, proj.v);
                    if self.classifier.is_inside(&uv_p) {
                        internal_uv.push(uv_p);
                    }
                }
            }
        }
        self.insert_internal_nodes(&internal_uv)?;

        let face_deflection = model.face(face_index)?.deflection();
        self.optimized_face_deflection = None;
        if self.pre_process_surface_nodes {
            self.generate_surface_nodes(model, face_index, params)?;
            let result = self.finish_mesh(&[], params, face_deflection, face_index);
            if self.delaun_failed {
                model
                    .face_mut(face_index)
                    .map_err(|e| format!("DelaunayNodeInsertionMeshAlgo::perform: {e}"))?
                    .set_status(MeshStatus::FAILURE);
                return Err(format!(
                    "DelaunayNodeInsertionMeshAlgo::perform: face {face_index}: BRepMesh_Delaun::addTriangle failed, no triangulation committed"
                ));
            }
            if let Some(d) = self.optimized_face_deflection.take() {
                model
                    .face_mut(face_index)
                    .map_err(|e| format!("DelaunayNodeInsertionMeshAlgo::perform: {e}"))?
                    .set_deflection(d);
            }
            return result;
        }
        let insert = self.list_surface_nodes(model, face_index, params)?;
        let result = self.finish_mesh(&insert, params, face_deflection, face_index);
        if std::env::var("DPROBE").ok().and_then(|v| v.trim().parse::<i32>().ok())
            == Some(face_index as i32)
        {
            if let Ok(res) = &result {
                eprintln!(
                    "DPROBE result face={face_index} nodes={} tris={} insert={}",
                    res.nodes.len(),
                    res.triangles.len(),
                    insert.len()
                );
                for (i, n) in res.nodes.iter().enumerate() {
                    eprintln!(
                        "DPROBE r i={} st={} uv=({:.6},{:.6}) p=({:.6},{:.6},{:.6})",
                        i + 1,
                        n.state.to_str(),
                        n.location.x(),
                        n.location.y(),
                        n.p3d.x(),
                        n.p3d.y(),
                        n.p3d.z()
                    );
                }
            }
        }
        if let Some(d) = self.optimized_face_deflection.take() {
            model
                .face_mut(face_index)
                .map_err(|e| format!("DelaunayNodeInsertionMeshAlgo::perform: {e}"))?
                .set_deflection(d);
        }
        if self.delaun_failed {
            // `BRepMesh_BaseMeshAlgo::Perform` (`BRepMesh_BaseMeshAlgo.cxx:40-62`):
            // the swallowed `Standard_OutOfRange` leaves the face with
            // `IMeshData_Failure` and no triangulation, and the remaining faces of
            // the shape are still meshed.
            model
                .face_mut(face_index)
                .map_err(|e| format!("DelaunayNodeInsertionMeshAlgo::perform: {e}"))?
                .set_status(MeshStatus::FAILURE);
            return Err(format!(
                "DelaunayNodeInsertionMeshAlgo::perform: face {face_index}: BRepMesh_Delaun::addTriangle failed, no triangulation committed"
            ));
        }
        result
    }

    /// Collects boundary UV for two OCCT paths that must not be mixed:
    ///
    /// * Delaunay frontier: every pcurve sample `0..=last` with consecutive
    ///   links only (`BRepMesh_BaseMeshAlgo.cxx:99-123`). No last-to-first
    ///   wrap of the concatenated wire — a periodic seam's ends differ in UV
    ///   and are the same 3D vertex; OCCT never adds that diagonal.
    /// * Classifier wire: skip last sample in traversal order
    ///   (`BRepMesh_NodeInsertionMeshAlgo.hxx:152-180` `collectWirePoints`).
    pub(super) fn collect_boundary_uv(
        model: &MeshModel,
        face_index: usize,
    ) -> Result<(Vec<GpPnt2d>, Vec<GpPnt>, Vec<(i32, i32)>, Vec<Vec<GpPnt2d>>), String> {
        let face = model.face(face_index)?;
        let mut uv: Vec<GpPnt2d> = Vec::new();
        let mut p3d: Vec<GpPnt> = Vec::new();
        let mut constraints: Vec<(i32, i32)> = Vec::new();
        let mut wires_uv: Vec<Vec<GpPnt2d>> = Vec::new();

        for (wire_it, &wire_index) in face.wires().iter().enumerate() {
            let wire = model.wire(wire_index)?;
            // `NodeInsertionMeshAlgo::initDataStructure` (`hxx:70-74`):
            // skip a self-intersecting wire, and skip an OpenWire unless it is
            // the outer wire (`aWireIt == 0`). Those wires are not AddPoint'd
            // and not RegisterWire'd, so they cannot flip registerSurfaceNodes
            // IN to OUT (`BRepMesh_Classifier.cxx:35-59`).
            let skip_classifier = wire.is_status(MeshStatus::SELF_INTERSECTING_WIRE)
                || (wire.is_status(MeshStatus::OPEN_WIRE) && wire_it != 0);
            let mut classifier_wire: Vec<GpPnt2d> = Vec::new();
            for j in 0..wire.edges_nb() {
                let edge_index = wire.edge(j)?;
                let orientation = wire.edge_orientation(j)?;
                let reverse_walk = wire.edge_reverse_walk(j);
                let edge = model.edge(edge_index)?;
                let Some(pcurve) = edge.pcurve_for(face_index, orientation) else {
                    continue;
                };
                let pts = pcurve.points();
                let pts3d = edge.discretization().points();
                let n = pts.len();
                if n == 0 {
                    continue;
                }
                if std::env::var("BPUV").is_ok() {
                    let (mut minv, mut maxv) = (f64::INFINITY, f64::NEG_INFINITY);
                    let (mut minu, mut maxu) = (f64::INFINITY, f64::NEG_INFINITY);
                    for p in pts {
                        minv = minv.min(p.y());
                        maxv = maxv.max(p.y());
                        minu = minu.min(p.x());
                        maxu = maxu.max(p.x());
                    }
                    eprintln!(
                        "BPUV face={face_index} wire={wire_it} edge={edge_index} n={n} \
first=({:.17e},{:.17e}) last=({:.17e},{:.17e}) \
u=[{:.17e},{:.17e}] v=[{:.17e},{:.17e}]",
                        pts[0].x(),
                        pts[0].y(),
                        pts[n - 1].x(),
                        pts[n - 1].y(),
                        minu,
                        maxu,
                        minv,
                        maxv
                    );
                }
                // Frontier links: every sample in the same traversal order as
                // `collectWirePoints`, including the last (`initDataStructure`
                // registers 0..=last). `Ordered < 0` walks the CheckOrder
                // pcurve backwards (same iso, not PCurve2).
                let start = uv.len();
                let reverse = !pcurve.is_forward() ^ reverse_walk;
                if !reverse {
                    for i in 0..n {
                        uv.push(pts[i]);
                        p3d.push(pts3d[i]);
                    }
                } else {
                    for i in (0..n).rev() {
                        uv.push(pts[i]);
                        p3d.push(pts3d[i]);
                    }
                }
                for t in start..uv.len() - 1 {
                    if uv[t] != uv[t + 1] {
                        constraints.push((t as i32, (t + 1) as i32));
                    }
                }

                // `collectWirePoints` (`hxx:153-180`): walk the pcurve in the
                // effective edge orientation. `reverse_walk` models
                // `ShapeExtend_WireData::Edge(signed)` = `Edge.Reverse()` on a
                // seam; `BRepMesh_ShapeVisitor::addWire` stores that reversed
                // orientation into `IMeshData_Wire`, so `GetEdgeOrientation`
                // (and thus `IMeshData_PCurve::IsForward`) already includes it.
                // The port keeps the original orientation plus `reverse_walk`,
                // so the classifier polygon must XOR it exactly like the
                // frontier links above.
                if !reverse {
                    classifier_wire.extend_from_slice(&pts[..n.saturating_sub(1)]);
                } else {
                    for i in (1..n).rev() {
                        classifier_wire.push(pts[i]);
                    }
                }
            }
            if !skip_classifier && classifier_wire.len() >= 2 {
                wires_uv.push(classifier_wire);
            }
        }

        Ok((uv, p3d, constraints, wires_uv))
    }

    /// Frontier/Fixed links from each pcurve in parameter order.
    /// Source: `BRepMesh_BaseMeshAlgo::initDataStructure` (`BRepMesh_BaseMeshAlgo.cxx:75-128`)
    /// via `NodeInsertionMeshAlgo::initDataStructure` returning `BaseAlgo::initDataStructure()`.
    pub(super) fn init_data_structure(
        &mut self,
        model: &MeshModel,
        face_index: usize,
    ) -> Result<(), String> {
        let wire_indices: Vec<usize> = model.face(face_index)?.wires().to_vec();
        for (wire_it, &wire_index) in wire_indices.iter().enumerate() {
            let (skip, slot_edges) = {
                let wire = model.wire(wire_index)?;
                // `NodeInsertionMeshAlgo::initDataStructure` (`hxx:70-74`):
                // skip a self-intersecting wire, and skip an OpenWire unless
                // it is the outer wire (`aWireIt == 0`).
                if wire.is_status(MeshStatus::SELF_INTERSECTING_WIRE)
                    || (wire.is_status(MeshStatus::OPEN_WIRE) && wire_it != 0)
                {
                    (true, Vec::new())
                } else {
                    let mut edges = Vec::new();
                    for j in 0..wire.edges_nb() {
                        edges.push(wire.edge(j)?);
                    }
                    (false, edges)
                }
            };
            if skip {
                continue;
            }
            // `BRepMesh_BaseMeshAlgo::initDataStructure` (`cxx:87-125`): for every
            // wire slot, *every* pcurve of that edge on this face is registered,
            // in `aDEdge->GetPCurves(myDFace)` order (the faces' pcurve index
            // list). A seam edge the wire walks twice therefore contributes both
            // of its period copies, and the walk order is the pcurve's own point
            // order - OCCT never consults the CheckOrder slot orientation here.
            // The link orientation is `fixSeamEdgeOrientation` (`cxx:197-225`).
            for edge_index in slot_edges {
                // `model` is borrowed independently of `self`, so the pcurve
                // and discretization point slices can be read in place while
                // nodes and links are registered; copying them per pcurve
                // (a seam edge contributes two) only churned memory.
                let edge = model.edge(edge_index)?;
                let pts3d = edge.discretization().points();
                for pc_index in edge.pcurves_for(face_index) {
                    let pc = edge.pcurve(pc_index)?;
                    let ori = Self::fix_seam_edge_orientation(&edge, pc_index, face_index)?;
                    let pts = pc.points();
                    // TEMPORARY diagnostic (DELREGSKEL): the registration skeleton in the
                    // `occt_probe --structdump` `OCCTSTRUCT E` line format, so the pcurve
                    // visit order and the effective link orientation can be diffed against
                    // `BRepMesh_BaseMeshAlgo::initDataStructure`.
                    if crate::debug_flags::del_reg_skel() {
                        let first = pts.first().copied().unwrap_or_else(GpPnt2d::zero);
                        let last = pts.last().copied().unwrap_or_else(GpPnt2d::zero);
                        eprintln!(
                            "REG f={} w={} edge={} pc={} ori={:?} n={} first=({:.9},{:.9}) last=({:.9},{:.9})",
                            face_index,
                            wire_it,
                            edge_index,
                            pc_index,
                            ori,
                            pts.len(),
                            first.x(),
                            first.y(),
                            last.x(),
                            last.y()
                        );
                    }
                    let mut prev = -1i32;
                    for (i, &uv) in pts.iter().enumerate() {
                        let p3d = pts3d.get(i).copied().unwrap_or_else(GpPnt::zero);
                        let node = self.register_node(p3d, uv, VertexState::Frontier, false);
                        if prev != -1 && prev != node {
                            let links_nb = self.structure.nb_links() as i32;
                            let link_index = self.add_link_to_mesh(prev, node, ori);
                            if wire_it != 0 && link_index <= links_nb {
                                self.structure
                                    .set_link_movability(link_index, VertexState::Fixed);
                            }
                        }
                        prev = node;
                    }
                }
            }
        }
        Ok(())
    }

}
