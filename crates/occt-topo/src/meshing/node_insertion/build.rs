//! `perform` (node insertion for a model face) and the structure initialisation.


use occt_core::gp::{GpPnt, GpPnt2d};
use occt_geom::geom_api::project_point_on_surface;


use super::super::delaun_types::VertexState;
use super::super::face_discret::Classifier;
#[allow(unused_imports)]
use super::super::mesh_algo::*;
use super::super::parameters::MeshParameters;
use super::super::range_splitter::{mesh_surface_type, create_range_splitter};
use super::*;
use super::parametric::uv_bounds;
use crate::shape::Edge;

/// Map a local `MeshEdge` polyline point into `occurrence`'s world frame.
///
/// The discrete edge stores the registry curve without the first-seen
/// occurrence's location (`MeshEdge::new`); each pcurve keeps the wire-edge
/// TopoDS location. World = `T_occ * P_local` (`BRepAdaptor_Curve::D0`).
fn occurrence_world_point(_mesh_edge: &Edge, occurrence: &Edge, p: &GpPnt) -> GpPnt {
    let occ_loc = occurrence.0.location();
    if occ_loc.is_identity() {
        return *p;
    }
    p.transformed(&occ_loc.transformation())
}


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
                if crate::debug_flags::bpuv() {
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

    /// Frontier/Fixed links from the wire-occurrence pcurve only.
    /// Source: `BRepMesh_BaseMeshAlgo::initDataStructure` (`BRepMesh_BaseMeshAlgo.cxx:75-145`,
    /// OCCT 8.0.1 #1338): a wire defines one topological occurrence of an edge;
    /// only that occurrence's pcurve creates mesh links. Other face pcurves still
    /// receive node indices for `Poly_PolygonOnTriangulation`.
    pub(super) fn init_data_structure(
        &mut self,
        model: &mut MeshModel,
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
                        edges.push((wire.edge(j)?, wire.edge_orientation(j)?));
                    }
                    (false, edges)
                }
            };
            if skip {
                continue;
            }
            // OCCT 8.0.1 (`cxx:89-145`): register the wire-occurrence pcurve with
            // links, then copy node indices onto every other face pcurve without
            // adding duplicate seam constraints.
            for (edge_index, wire_ori) in slot_edges {
                let edge = model.edge(edge_index)?;
                let pts3d = edge.discretization().points();
                let all_pc = edge.pcurves_for(face_index);
                // `GetPCurve(face, wireOrientation)` (`cxx:92-93`): prefer the
                // pcurve whose orientation matches the wire occurrence.
                let primary_index = all_pc
                    .iter()
                    .copied()
                    .find(|&i| {
                        edge.pcurve(i)
                            .map(|p| p.orientation() == wire_ori)
                            .unwrap_or(false)
                    })
                    .or_else(|| all_pc.last().copied());
                let Some(primary_index) = primary_index else {
                    continue;
                };

                let mut pc_registrations: Vec<(usize, Vec<i32>)> = Vec::new();
                {
                    let pc = edge.pcurve(primary_index)?;
                    let ori = Self::fix_seam_edge_orientation(&edge, primary_index, face_index)?;
                    let pts = pc.points();
                    if crate::debug_flags::del_reg_skel() {
                        let first = pts.first().copied().unwrap_or_else(GpPnt2d::zero);
                        let last = pts.last().copied().unwrap_or_else(GpPnt2d::zero);
                        eprintln!(
                            "REG f={} w={} edge={} pc={} ori={:?} n={} first=({:.9},{:.9}) last=({:.9},{:.9})",
                            face_index,
                            wire_it,
                            edge_index,
                            primary_index,
                            ori,
                            pts.len(),
                            first.x(),
                            first.y(),
                            last.x(),
                            last.y()
                        );
                    }
                    let mut prev = -1i32;
                    let mut pc_nodes: Vec<i32> = Vec::with_capacity(pts.len());
                    for (i, &uv) in pts.iter().enumerate() {
                        let p3d_raw = pts3d.get(i).copied().unwrap_or_else(GpPnt::zero);
                        let p3d = occurrence_world_point(edge.edge(), pc.occurrence(), &p3d_raw);
                        let node = self.register_node(p3d, uv, VertexState::Frontier, false);
                        pc_nodes.push(node);
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
                    pc_registrations.push((primary_index, pc_nodes));
                }

                // Other pcurves: node indices only, no mesh links (`cxx:126-144`).
                for pc_index in all_pc {
                    if pc_index == primary_index {
                        continue;
                    }
                    let pc = edge.pcurve(pc_index)?;
                    let pts = pc.points();
                    let mut pc_nodes: Vec<i32> = Vec::with_capacity(pts.len());
                    for (i, &uv) in pts.iter().enumerate() {
                        let p3d_raw = pts3d.get(i).copied().unwrap_or_else(GpPnt::zero);
                        let p3d = occurrence_world_point(edge.edge(), pc.occurrence(), &p3d_raw);
                        let node = self.register_node(p3d, uv, VertexState::Frontier, false);
                        pc_nodes.push(node);
                    }
                    pc_registrations.push((pc_index, pc_nodes));
                }

                for (pc_index, pc_nodes) in pc_registrations {
                    let pcurve = model.edge_mut(edge_index)?.pcurve_mut(pc_index)?;
                    for (i, &node) in pc_nodes.iter().enumerate() {
                        *pcurve.index_mut(i)? = node;
                    }
                }
            }
        }
        Ok(())
    }

}
