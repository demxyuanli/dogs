//! Port of OCCT model_builder — Wave 1 BRepMesh.
//!
//! `BRepMesh_ModelBuilder` turns a `TopoDS_Shape` into a discrete `MeshModel`
//! (faces, wires, edges, 3D curves and 2D pcurves) and
//! `BRepMesh_ModelPreProcessor` initializes per-entity deflection/status before
//! the edge/face discretizers run.
//!
//! Reference: `BRepMesh_ModelBuilder.hxx/.cxx`, `BRepMesh_ModelPreProcessor.hxx/.cxx`,
//! `BRepMeshData_Model.hxx`, `BRepMesh_ShapeVisitor.cxx`, `IMeshData_*.hxx`.
//!
//! Data types (`MeshModel`/`MeshEdge`/`MeshFace`/`MeshWire`/`MeshCurve`/
//! `MeshPCurve`/`MeshStatus`) come from `data_model` and `MeshParameters` from
//! `parameters` — the Wave 1 sibling modules.

use std::collections::HashMap;
use std::sync::Arc;

use occt_core::bnd::BndBox;
use occt_core::gp::GpPnt2d;
use occt_geom2d::curve::Curve2d;

use crate::abs::{Orientation, ShapeType};
use crate::bbox_from_geometry::shape_bbox;
use crate::brep_tool::BRepTool;
use crate::pcurve_full::make_pcurve_full;
use crate::shape::{Edge, Face, TopoShape, Wire};
use crate::topo_tools_full::{edge_vertices, edges_of, edges_of_wire, faces_of, wires_of_face};

use super::data_model::{MeshEdge, MeshModel, MeshStatus};
use super::parameters::MeshParameters;
use super::wire_order::{WireOrder, WireOrderStatus};

/// Maximum dimension of a bounding box (`BRepMesh_ShapeTool::BoxMaxDimension`).
fn box_max_dimension(b: &BndBox) -> Option<f64> {
    if b.is_void() {
        return None;
    }
    let (xmin, xmax, ymin, ymax, zmin, zmax) = b.get()?;
    Some((xmax - xmin).max(ymax - ymin).max(zmax - zmin))
}

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

/// Port of `BRepMesh_ShapeVisitor::addWire` (2D pcurve mode).
///
/// Orders the wire's edges into one connected chain
/// (`ShapeAnalysis_Wire::CheckOrder` → `ShapeAnalysis_WireOrder::Perform`) and,
/// for every non-EXTERNAL edge, registers its pcurve on the face and its
/// position/orientation in the wire chain. Returns `false` when the wire cannot
/// be ordered (a missing pcurve → `ShapeExtend_FAIL`); a wire whose edges needed
/// reversing sets the face `UNORIENTED_WIRE` status (`ShapeExtend_DONE3`).
fn add_wire(
    model: &mut MeshModel,
    face_index: usize,
    wire: &Wire,
    edge_index: &mut HashMap<usize, usize>,
) -> bool {
    let face = match model.face(face_index) {
        Ok(f) => f.face().clone(),
        Err(_) => return false,
    };
    // `ShapeExtend_WireData(theWire, chained=true, manifold=false)` keeps only
    // FORWARD/REVERSED edges in the main list; INTERNAL/EXTERNAL are non-manifold.
    let stored: Vec<Edge> = edges_of_wire(wire)
        .into_iter()
        .filter(|e| {
            let o = e.0.orientation();
            o == Orientation::Forward || o == Orientation::Reversed
        })
        .collect();
    if stored.is_empty() {
        return false;
    }

    // 2D endpoints in each edge's *traversal* direction, mirroring
    // `ShapeAnalysis_Edge::PCurve(..., orient=true)`: a reversed edge toggles
    // `cf`/`cl`, so its start point is the pcurve's natural end.
    let mut order = WireOrder::new();
    for e in &stored {
        // Pass the edge's *actual* orientation: a seam edge has two pcurves
        // (one per side) and must hand each traversal its own side, not the
        // forward pcurve twice.
        let pc = match make_pcurve_full(e, &face) {
            Ok(pc) => pc,
            Err(_) => return false,
        };
        let fwd = Edge(e.0.oriented(Orientation::Forward));
        let (a, b) = BRepTool::edge_parameters(&fwd);
        if !a.is_finite() || !b.is_finite() || b - a < 1e-15 {
            return false;
        }
        let (p_a, p_b) = (pc.d0(a), pc.d0(b));
        let (begin, end) = if e.0.orientation().is_reversed() {
            (p_b, p_a)
        } else {
            (p_a, p_b)
        };
        order.add_edge(begin, end);
    }
    order.perform();

    if order.status() == WireOrderStatus::Reversed {
        if let Ok(f) = model.face_mut(face_index) {
            f.set_status(MeshStatus::UNORIENTED_WIRE);
        }
    }
    if order.nb_edges() != stored.len() {
        return false;
    }

    let wire_index = model.add_wire(wire.clone());
    for i in 1..=stored.len() {
        let signed = order.ordered(i);
        let e = &stored[signed.unsigned_abs() as usize - 1];
        let orientation = if signed < 0 {
            e.0.orientation().reversed()
        } else {
            e.0.orientation()
        };
        if orientation == Orientation::External {
            continue;
        }
        let key = Arc::as_ptr(&e.0.tshape) as usize;
        let eidx = match edge_index.get(&key) {
            Some(&i) => i,
            None => {
                let i = model.add_edge(e.clone());
                edge_index.insert(key, i);
                i
            }
        };
        if let Ok(edge) = model.edge_mut(eidx) {
            edge.add_pcurve(face_index, orientation);
        }
        if let Ok(w) = model.wire_mut(wire_index) {
            w.add_edge(eidx, orientation);
        }
    }
    if let Ok(f) = model.face_mut(face_index) {
        f.add_wire(wire_index);
    }
    true
}

/// Port of `ShapeAnalysis::OuterWire`: the first wire whose 2D signed area
/// (`ShapeAnalysis::TotCross2D`) is non-negative, else the last wire.
fn outer_wire(face: &Face) -> Option<Wire> {
    let wires = wires_of_face(face);
    if wires.is_empty() {
        return None;
    }
    for (i, w) in wires.iter().enumerate() {
        if i == wires.len() - 1 {
            return Some(w.clone());
        }
        if wire_area_2d(w, face) >= 0.0 {
            return Some(w.clone());
        }
    }
    None
}

/// `ShapeAnalysis::TotCross2D` — signed 2D area of a wire's pcurves (trapezoid
/// rule over sampled pcurve points, sequence reversed per REVERSED edge).
fn wire_area_2d(wire: &Wire, face: &Face) -> f64 {
    let mut totcross = 0.0;
    let mut uv0: Option<GpPnt2d> = None;
    let mut fuv = GpPnt2d::new(0.0, 0.0);
    let mut nbc = 0usize;
    for e in edges_of_wire(wire) {
        let fwd = Edge(e.0.oriented(Orientation::Forward));
        let Ok(pc) = make_pcurve_full(&fwd, face) else { continue };
        let (a, b) = BRepTool::edge_parameters(&fwd);
        if !a.is_finite() || !b.is_finite() {
            continue;
        }
        let mut pts = sample_pcurve(pc.as_ref(), a, b);
        if e.0.orientation().is_reversed() {
            pts.reverse();
        }
        nbc += 1;
        if nbc == 1 {
            fuv = pts[0];
            uv0 = Some(pts[0]);
        }
        for p in &pts {
            totcross += (fuv.x() - p.x()) * (fuv.y() + p.y()) / 2.0;
            fuv = *p;
        }
    }
    if let Some(u0) = uv0 {
        totcross += (fuv.x() - u0.x()) * (fuv.y() + u0.y()) / 2.0;
    }
    totcross
}

/// Uniform sample of a pcurve over `[a, b]` (endpoints + interior), standing in
/// for `ShapeAnalysis_Curve::GetSamplePoints`.
fn sample_pcurve(pc: &dyn Curve2d, a: f64, b: f64) -> Vec<GpPnt2d> {
    const N: usize = 8;
    let mut pts = Vec::with_capacity(N);
    for i in 0..N {
        let t = a + (b - a) * i as f64 / (N - 1) as f64;
        pts.push(pc.d0(t));
    }
    pts
}

/// Tool for building a discrete model from a topological shape.
/// Port of `BRepMesh_ModelBuilder` (+ the visiting logic of
/// `BRepMesh_ShapeVisitor`).
pub struct ModelBuilder;

impl ModelBuilder {
    /// Build the discrete model of `shape` under `params`.
    ///
    /// Mirrors `BRepMesh_ModelBuilder::performInternal` plus the
    /// `BRepMesh_ShapeVisitor` walk: every face gets its outer wire first
    /// (`ShapeAnalysis::OuterWire`), then its inner wires, each wire reordered
    /// into a connected chain (`ShapeAnalysis_Wire::CheckOrder`), with EXTERNAL
    /// edges skipped and shared edges deduplicated by TShape identity. Free edges
    /// (not bounding any face) are added as well.
    ///
    /// Returns `Err` when the shape is empty (void bounding box → `Message_Fail1`).
    pub fn build_model(shape: &TopoShape, params: &MeshParameters) -> Result<MeshModel, String> {
        let mut model = MeshModel::new(shape.clone());

        // Maximum size of the shape's bounding box (used by relative deflection).
        let bbox = shape_bbox(shape);
        if bbox.is_void() {
            return Err("BRepMesh_ModelBuilder::build_model: empty shape".to_string());
        }
        let max_size = if params.relative {
            box_max_dimension(&bbox).unwrap_or(0.0)
        } else {
            params.deflection.max(params.deflection_interior.max(0.0))
        };
        model.set_max_size(max_size);

        // Visit(Face): outer wire first; a failure on the outer wire fails the
        // face, a failure on an inner wire only marks it unoriented.
        let mut edge_index: HashMap<usize, usize> = HashMap::new();
        for f in faces_of(shape) {
            let face_index = model.add_face(f.clone());

            let outer = outer_wire(&f);
            if let Some(outer_wire) = &outer {
                if !add_wire(&mut model, face_index, outer_wire, &mut edge_index) {
                    if let Ok(fm) = model.face_mut(face_index) {
                        fm.set_status(MeshStatus::FAILURE);
                    }
                    continue;
                }
            }

            for w in wires_of_face(&f) {
                if let Some(outer_wire) = &outer {
                    if Arc::ptr_eq(&w.0.tshape, &outer_wire.0.tshape) {
                        continue;
                    }
                }
                if !add_wire(&mut model, face_index, &w, &mut edge_index) {
                    if let Ok(fm) = model.face_mut(face_index) {
                        fm.set_status(MeshStatus::UNORIENTED_WIRE);
                    }
                }
            }
        }

        // Free edges that do not bound any face.
        for e in edges_of(shape) {
            let key = Arc::as_ptr(&e.0.tshape) as usize;
            if !edge_index.contains_key(&key) {
                model.add_edge(e.clone());
            }
        }

        Ok(model)
    }
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
                (2.0 * super::shape_tool::ShapeTool::max_face_tolerance(&face_shape)).max(face_def);
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
    ///
    /// ponytail: `SeamEdgeAmplifier` (cone seam-edge splitting) is not ported — a
    /// density refinement for cone seam edges.
    pub fn perform(model: &mut MeshModel, _params: &MeshParameters) -> bool {
        model.faces_nb() != 0 || model.edges_nb() != 0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitives::BRepPrimBox;

    fn unit_box() -> TopoShape {
        BRepPrimBox::make_box(1.0, 1.0, 1.0).solid.0.clone()
    }

    #[test]
    fn build_model_box_counts() {
        let shape = unit_box();
        let params = MeshParameters::default();
        let model = ModelBuilder::build_model(&shape, &params).expect("model built");
        assert_eq!(model.faces_nb(), 6, "box has 6 faces");
        assert_eq!(model.edges_nb(), 12, "box has 12 edges");
        for i in 0..model.faces_nb() {
            let f = model.face(i).expect("face index");
            assert_eq!(f.wires_nb(), 1, "each box face has one wire");
            let wi = f.wire(0).expect("wire index");
            assert_eq!(
                model.wire(wi).expect("wire in model").edges_nb(),
                4,
                "each wire has 4 edges"
            );
        }
        // Every edge is shared by exactly two faces → two pcurves.
        for i in 0..model.edges_nb() {
            assert_eq!(model.edge(i).expect("edge index").pcurves_nb(), 2, "interior edge bounds two faces");
        }
    }

    #[test]
    fn build_model_sets_max_size() {
        let shape = unit_box();
        // Non-relative mode: max_size is max(deflection, deflection_interior),
        // matching BRepMesh_ModelBuilder::performInternal.
        let params = MeshParameters::default();
        let model = ModelBuilder::build_model(&shape, &params).expect("model built");
        assert!(
            (model.max_size() - params.deflection).abs() < 1e-12,
            "max_size {}",
            model.max_size()
        );

        // Relative mode: max_size is the bounding-box maximum dimension.
        let rel_params = MeshParameters { relative: true, ..MeshParameters::default() };
        let model = ModelBuilder::build_model(&shape, &rel_params).expect("model built");
        assert!((model.max_size() - 1.0).abs() < 1e-9, "rel max_size {}", model.max_size());
    }

    #[test]
    fn build_model_empty_shape_fails() {
        let shape = TopoShape::new(ShapeType::Compound);
        let params = MeshParameters::default();
        assert!(ModelBuilder::build_model(&shape, &params).is_err());
    }

    #[test]
    fn preprocessor_sets_deflections_and_status() {
        let shape = unit_box();
        let params = MeshParameters::default();
        let mut model = ModelBuilder::build_model(&shape, &params).expect("model built");
        for i in 0..model.edges_nb() {
            ModelPreProcessor::compute_edge_deflection(&mut model, i, &params).expect("edge deflection");
        }
        for i in 0..model.wires_nb() {
            ModelPreProcessor::compute_wire_deflection(&mut model, i, &params).expect("wire deflection");
        }
        for i in 0..model.faces_nb() {
            ModelPreProcessor::compute_face_deflection(&mut model, i, &params).expect("face deflection");
        }
        for i in 0..model.edges_nb() {
            let e = model.edge(i).expect("edge index");
            assert!(e.deflection() >= params.deflection);
            assert!(e.is_status(MeshStatus::OUTDATED));
        }
        for i in 0..model.faces_nb() {
            let f = model.face(i).expect("face index");
            assert!(f.deflection() >= params.deflection);
            assert!(f.is_status(MeshStatus::OUTDATED));
        }
    }

    #[test]
    fn preprocessor_relative_mode_scales_by_edge_size() {
        let shape = unit_box();
        let params = MeshParameters { relative: true, ..MeshParameters::default() };
        let mut model = ModelBuilder::build_model(&shape, &params).expect("model built");
        for i in 0..model.edges_nb() {
            ModelPreProcessor::compute_edge_deflection(&mut model, i, &params).expect("edge deflection");
        }
        for i in 0..model.edges_nb() {
            let d = model.edge(i).expect("edge index").deflection();
            // Unit box: edge size 1, model max size 1 → coefficient 1/(2·1)=0.5,
            // so the relative deflection becomes 0.5·deflection
            // (`BRepMesh_Deflection::ComputeAbsoluteDeflection`).
            assert!((d - 0.5 * params.deflection).abs() < 1e-9, "rel deflection {d}");
        }
    }
}
