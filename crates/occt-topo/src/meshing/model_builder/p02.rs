use super::prelude::*;
use super::*;

/// Approximate length of an edge from its curve endpoints.
fn edge_length(e: &Edge) -> f64 {
    match BRepTool::edge_vertices(e) {
        Some((a, b)) => a.distance(&b),
        None => 0.0,
    }
}

/// `BRepMesh_Deflection::ComputeAbsoluteDeflection`: convert a relative
/// deflection into an absolute one for a shape of bounding-box max dimension
/// `shape_size` inside a shape of max dimension `max_shape_size`. The coefficient
/// `max_size / (2·shape_size)` is clamped to `[0.5, 2.0]`.
fn compute_absolute_deflection(shape_size: f64, max_shape_size: f64, relative: f64) -> f64 {
    if relative <= 0.0 {
        return 0.0;
    }
    // OCCT seeds aShapeSize with the relative value and overwrites it with the
    // shape's bbox max dimension when non-void.
    let shape_size = if shape_size > 0.0 { shape_size } else { relative };
    let max_size = if max_shape_size > 0.0 { max_shape_size } else { shape_size };
    let coeff = (max_size / (2.0 * shape_size)).clamp(0.5, 2.0);
    coeff * shape_size * relative
}

/// Distance from each edge endpoint's vertex point to the curve at the matching
/// parameter (`-1.0` when a vertex or curve is missing) — the `aDistF`/`aDistL`
/// of `BRepMesh_Deflection::ComputeDeflection`.
fn edge_vertex_adjust(e: &MeshEdge) -> (f64, f64) {
    let Some(curve) = e.curve() else {
        return (-1.0, -1.0);
    };
    let (v1, v2) = edge_vertices(e.edge());
    let (f, l) = (e.first_parameter(), e.last_parameter());
    let df = v1
        .map(|v| BRepTool::vertex_point(&v).distance(&curve.d0(f)))
        .unwrap_or(-1.0);
    let dl = v2
        .map(|v| BRepTool::vertex_point(&v).distance(&curve.d0(l)))
        .unwrap_or(-1.0);
    (df, dl)
}

/// Tool for pre-processing a discrete model before meshing.
/// Port of `BRepMesh_ModelPreProcessor`: initializes per-entity deflection
/// limits and marks all entities as `Outdated` (nothing to reuse in a fresh
/// model — the real triangulation-consistency reuse check lands in Wave 2/3).
pub struct ModelPreProcessor;

impl ModelPreProcessor {
    /// `BRepMesh_Deflection::ComputeDeflection(edge)` — set the linear and angular
    /// deflection of one edge (absolute deflection in relative mode, plus the
    /// vertex-adjustment floor).
    pub fn compute_edge_deflection(
        model: &mut MeshModel,
        edge_index: usize,
        params: &MeshParameters,
    ) -> Result<(), String> {
        let max_size = model.max_size();
        let edge = model.edge(edge_index)?;
        let size = edge_length(edge.edge());
        let mut lin = if params.relative {
            compute_absolute_deflection(size, max_size, params.deflection)
        } else {
            params.deflection
        };
        let (df, dl) = edge_vertex_adjust(edge);
        lin = lin.max(df.max(dl));
        let e = model.edge_mut(edge_index)?;
        e.set_deflection(lin.max(1e-7));
        e.set_angular_deflection(params.angle);
        e.set_status(MeshStatus::OUTDATED);
        Ok(())
    }

    /// `BRepMesh_Deflection::ComputeDeflection(wire)` — mean of the wire's edge
    /// deflections (or `params.deflection` for an empty wire).
    pub fn compute_wire_deflection(
        model: &mut MeshModel,
        wire_index: usize,
        params: &MeshParameters,
    ) -> Result<(), String> {
        let w = model.wire(wire_index)?;
        let def = if w.edges_nb() > 0 {
            let sum: f64 = (0..w.edges_nb())
                .map(|j| model.edge(w.edge(j)?).map(|e| e.deflection()))
                .collect::<Result<Vec<_>, _>>()?
                .into_iter()
                .sum();
            sum / w.edges_nb() as f64
        } else {
            params.deflection
        };
        model.wire_mut(wire_index)?.set_deflection(def);
        Ok(())
    }

    /// `BRepMesh_Deflection::ComputeDeflection(face)` — the face interior
    /// deflection (or its absolute form in relative mode) floored by the mean of
    /// the wire deflections and `2·MaxFaceTolerance`.
    pub fn compute_face_deflection(
        model: &mut MeshModel,
        face_index: usize,
        params: &MeshParameters,
    ) -> Result<(), String> {
        let (face_shape, wires, force) = {
            let f = model.face(face_index)?;
            (f.face().clone(), f.wires().to_vec(), params.force_face_deflection)
        };
        let interior = if params.relative {
            let size = box_max_dimension(&shape_bbox(&face_shape)).unwrap_or(0.0);
            compute_absolute_deflection(size, -1.0, params.deflection_interior)
        } else {
            params.deflection_interior
        };
        let mut face_def = 0.0;
        if !force {
            if !wires.is_empty() {
                let sum: f64 = wires
                    .iter()
                    .map(|&wi| model.wire(wi).map(|w| w.deflection()))
                    .collect::<Result<Vec<_>, _>>()?
                    .into_iter()
                    .sum();
                face_def = sum / wires.len() as f64;
            }
            face_def =
                (2.0 * ShapeTool::max_face_tolerance(&face_shape)).max(face_def);
        }
        face_def = interior.max(face_def);
        let fm = model.face_mut(face_index)?;
        fm.set_deflection(face_def.max(1e-7));
        fm.set_status(MeshStatus::OUTDATED);
        Ok(())
    }

    /// `BRepMesh_ModelPreProcessor::performInternal` — seam-edge amplification on
    /// cone faces, triangulation-consistency (reuse) check, and cleanup of
    /// outdated polygons. A fresh model has no stored triangulation, so the
    /// consistency/cleanup are no-ops.
    pub fn perform(model: &mut MeshModel, params: &MeshParameters) -> bool {
        amplify_cone_seams(model, params);
        model.faces_nb() != 0 || model.edges_nb() != 0
    }
}

/// `BRepMesh_ModelPreProcessor::SeamEdgeAmplifier` — split 2-point cone seam
/// edges so `collectWirePoints` (which drops each pcurve's last sample) still
/// leaves a vertical UV side. Source: `BRepMesh_ModelPreProcessor.cxx:118-151`.
fn amplify_cone_seams(model: &mut MeshModel, params: &MeshParameters) {
    let faces_nb = model.faces_nb();
    for face_index in 0..faces_nb {
        amplify_cone_face(model, face_index, params);
    }
}

fn amplify_cone_face(model: &mut MeshModel, face_index: usize, params: &MeshParameters) {
    let (wire_index, n_edges) = {
        let Ok(face) = model.face(face_index) else {
            return;
        };
        if face.is_status(MeshStatus::FAILURE) {
            return;
        }
        let Some(surface) = face.surface() else {
            return;
        };
        if classify_surface(surface.as_ref()) != SurfaceType::Cone {
            return;
        }
        if face.wires_nb() == 0 {
            return;
        }
        let Ok(wire_index) = face.wire(0) else {
            return;
        };
        let Ok(wire) = model.wire(wire_index) else {
            return;
        };
        (wire_index, wire.edges_nb())
    };
    if n_edges < 2 {
        return;
    }
    for edge_slot in 0..n_edges - 1 {
        let Ok(edge_index) = model.wire(wire_index).and_then(|w| w.edge(edge_slot)) else {
            continue;
        };
        let is_seam = {
            let Ok(edge) = model.edge(edge_index) else {
                continue;
            };
            let fwd = pcurve_index_for(edge, face_index, Orientation::Forward);
            let rev = pcurve_index_for(edge, face_index, Orientation::Reversed);
            match (fwd, rev) {
                (Some(a), Some(b)) => a != b,
                _ => false,
            }
        };
        if !is_seam {
            continue;
        }
        let n3d = model
            .edge(edge_index)
            .map(|e| e.discretization().parameters_nb())
            .unwrap_or(0);
        if n3d == 2 {
            let du = cone_seam_step(model, face_index, params).abs();
            split_cone_seam(model, edge_index, face_index, du);
        }
        return;
    }
}

/// `SeamEdgeAmplifier::getConeStep`.
fn cone_seam_step(model: &MeshModel, face_index: usize, params: &MeshParameters) -> f64 {
    let Ok(face) = model.face(face_index) else {
        return 0.0;
    };
    let Ok(wire_index) = face.wire(0) else {
        return 0.0;
    };
    let Ok(wire) = model.wire(wire_index) else {
        return 0.0;
    };
    let mut splitter = ConeRangeSplitter::new();
    splitter.reset(face, params);
    for j in 0..wire.edges_nb() {
        let Ok(ei) = wire.edge(j) else {
            continue;
        };
        let Ok(ori) = wire.edge_orientation(j) else {
            continue;
        };
        let Ok(edge) = model.edge(ei) else {
            continue;
        };
        let Some(pc) = edge.pcurve_for(face_index, ori) else {
            continue;
        };
        for k in 0..pc.parameters_nb() {
            if let Ok(uv) = pc.get_point(k) {
                splitter.add_point(uv);
            }
        }
    }
    let mut steps_nb = (0i32, 0i32);
    splitter.get_split_steps(params, &mut steps_nb).1
}

/// `SeamEdgeAmplifier::splitEdge`.
fn split_cone_seam(model: &mut MeshModel, edge_index: usize, face_index: usize, du: f64) -> bool {
    let (first, last, n_pc, y0, y1, last_u) = {
        let Ok(edge) = model.edge(edge_index) else {
            return false;
        };
        if edge.pcurves_nb() < 2 {
            return false;
        }
        let pc0 = match edge.pcurve(0) {
            Ok(p) => p,
            Err(_) => return false,
        };
        let n = pc0.parameters_nb();
        if n < 2 {
            return false;
        }
        let Ok(p0) = pc0.get_point(0) else {
            return false;
        };
        let Ok(p1) = pc0.get_point(n - 1) else {
            return false;
        };
        (
            edge.first_parameter(),
            edge.last_parameter(),
            edge.pcurves_nb(),
            p0.y(),
            p1.y(),
            p1.x(),
        )
    };
    let a_mod = (y0 - y1).abs();
    if a_mod < RESOLUTION {
        return false;
    }
    let dt = (last - first).abs() / a_mod * du;
    if !(dt.is_finite() && dt > PCONFUSION) {
        return false;
    }
    let (curve3d, topo_edge, topo_face) = {
        let Ok(edge) = model.edge(edge_index) else {
            return false;
        };
        let Some(curve3d) = edge.curve() else {
            return false;
        };
        let topo_edge = edge.edge().clone();
        let Ok(topo_face) = model.face(face_index).map(|f| f.face().clone()) else {
            return false;
        };
        (curve3d, topo_edge, topo_face)
    };
    {
        let Ok(em) = model.edge_mut(edge_index) else {
            return false;
        };
        if !split_curve_3d(curve3d.as_ref(), em.discretization_mut(), dt) {
            return false;
        }
    }
    let mut pc1 = match make_pcurve_full(
        &Edge(topo_edge.0.oriented(Orientation::Forward)),
        &topo_face,
    ) {
        Ok(c) => c,
        Err(_) => return false,
    };
    let mut pc2 = match make_pcurve_full(
        &Edge(topo_edge.0.oriented(Orientation::Reversed)),
        &topo_face,
    ) {
        Ok(c) => c,
        Err(_) => return false,
    };
    let geom_first = pc1.d0(pc1.first_parameter());
    if (last_u - geom_first.x()).abs() > CONFUSION {
        std::mem::swap(&mut pc1, &mut pc2);
    }
    if n_pc < 2 {
        return false;
    }
    if let Ok(em) = model.edge_mut(edge_index) {
        if let Ok(p0) = em.pcurve_mut(0) {
            split_curve_2d(pc1.as_ref(), p0, dt);
        }
        if let Ok(p1) = em.pcurve_mut(1) {
            split_curve_2d(pc2.as_ref(), p1, dt);
        }
    }
    true
}

/// `SeamEdgeAmplifier::splitCurve` for the 3D polyline.
fn split_curve_3d(geom: &dyn Curve, disc: &mut MeshCurve, dt: f64) -> bool {
    let n = disc.parameters_nb();
    if n < 2 {
        return false;
    }
    let Ok(first) = disc.get_parameter(0) else {
        return false;
    };
    let Ok(last) = disc.get_parameter(n - 1) else {
        return false;
    };
    let reversed = first > last;
    let mut updated = false;
    let mut k = 1i32;
    loop {
        let curr = first + f64::from(k) * dt * if reversed { -1.0 } else { 1.0 };
        if seam_split_done(curr, last, reversed) {
            break;
        }
        let pos = disc.parameters_nb() - 1;
        if disc.insert_point(pos, geom.d0(curr), curr).is_err() {
            break;
        }
        updated = true;
        k += 1;
    }
    updated
}

/// `SeamEdgeAmplifier::splitCurve` for a pcurve polyline.
fn split_curve_2d(geom: &dyn Curve2d, disc: &mut MeshPCurve, dt: f64) -> bool {
    let n = disc.parameters_nb();
    if n < 2 {
        return false;
    }
    let Ok(first) = disc.get_parameter(0) else {
        return false;
    };
    let Ok(last) = disc.get_parameter(n - 1) else {
        return false;
    };
    let reversed = first > last;
    let mut updated = false;
    let mut k = 1i32;
    loop {
        let curr = first + f64::from(k) * dt * if reversed { -1.0 } else { 1.0 };
        if seam_split_done(curr, last, reversed) {
            break;
        }
        let pos = disc.parameters_nb() - 1;
        if disc.insert_point(pos, geom.d0(curr), curr).is_err() {
            break;
        }
        updated = true;
        k += 1;
    }
    updated
}

fn seam_split_done(curr: f64, last: f64, reversed: bool) -> bool {
    if reversed {
        curr - last < PCONFUSION
    } else {
        !(curr - last < -PCONFUSION)
    }
}

/// `IMeshData_Edge::GetPCurve(face, orientation)` — matching handle index.
fn pcurve_index_for(edge: &MeshEdge, face: usize, orientation: Orientation) -> Option<usize> {
    let ids = edge.pcurves_for(face);
    ids.iter()
        .copied()
        .find(|&i| {
            edge.pcurve(i)
                .map(|p| p.orientation() == orientation)
                .unwrap_or(false)
        })
        .or_else(|| ids.last().copied())
}
