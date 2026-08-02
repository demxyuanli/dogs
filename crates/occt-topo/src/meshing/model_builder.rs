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

use crate::abs::ShapeType;
use crate::bbox_from_geometry::shape_bbox;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, TopoShape};
use crate::topo_tools_full::{edges_of, edges_of_wire, faces_of, wires_of_face};

use super::data_model::{MeshModel, MeshStatus};
use super::parameters::MeshParameters;

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

/// Tool for building a discrete model from a topological shape.
/// Port of `BRepMesh_ModelBuilder` (+ the visiting logic of
/// `BRepMesh_ShapeVisitor`).
pub struct ModelBuilder;

impl ModelBuilder {
    /// Build the discrete model of `shape` under `params`.
    ///
    /// Mirrors `BRepMesh_ModelBuilder::performInternal` plus the
    /// `BRepMesh_ShapeVisitor` walk: distinct edges are collected once and
    /// shared between faces; every face gets its wires, every wire its ordered
    /// edges, and each edge gets one pcurve per adjacent face. Free edges (not
    /// bounding any face) are added as well.
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

        // Visit faces → wires → edges, deduplicating edges by TShape identity
        // so a shared edge yields a single MeshEdge (ShapeVisitor::Visit(Edge)).
        let mut edge_index: HashMap<usize, usize> = HashMap::new();
        for f in faces_of(shape) {
            let face_index = model.add_face(f.clone());
            let wires = wires_of_face(&f);
            for w in wires {
                let wire_index = model.add_wire(w.clone());
                for e in edges_of_wire(&w) {
                    let key = Arc::as_ptr(&e.0.tshape) as usize;
                    let eidx = match edge_index.get(&key) {
                        Some(&i) => i,
                        None => {
                            let i = model.add_edge(e.clone());
                            edge_index.insert(key, i);
                            i
                        }
                    };
                    let orientation = e.0.orientation();
                    model
                        .edge_mut(eidx)
                        .expect("edge index just added")
                        .add_pcurve(face_index, orientation);
                    model
                        .wire_mut(wire_index)
                        .expect("wire index just added")
                        .add_edge(eidx, orientation);
                }
                model
                    .face_mut(face_index)
                    .expect("face index just added")
                    .add_wire(wire_index);
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
    /// Initialize deflections and statuses of every face/edge in `model`.
    ///
    /// In relative mode each edge's deflection is `params.deflection * edge length`
    /// (OCCT: `<deflection> * size of edge`); otherwise it is
    /// `max(deflection, deflection_interior)`. Face deflection is the maximum
    /// deflection of the face's boundary edges. Returns `false` for an empty
    /// model.
    pub fn perform(model: &mut MeshModel, params: &MeshParameters) -> bool {
        if model.faces_nb() == 0 && model.edges_nb() == 0 {
            return false;
        }

        for i in 0..model.edges_nb() {
            let len = edge_length(model.edge(i).expect("edge index in range").edge());
            let d = if params.relative {
                if len > 0.0 {
                    params.deflection * len
                } else {
                    params.deflection
                }
            } else {
                params.deflection.max(params.deflection_interior.max(0.0))
            };
            let e = model.edge_mut(i).expect("edge index in range");
            e.set_deflection(d.max(1e-7));
            e.set_angular_deflection(params.angle);
            // Fresh model: no prior triangulation, so nothing is reusable.
            e.set_status(MeshStatus::OUTDATED);
        }

        for i in 0..model.faces_nb() {
            // Face deflection is the maximum deflection of its boundary edges.
            let mut d = 0.0f64;
            let wire_indices: Vec<usize> = model.face(i).expect("face index in range").wires().to_vec();
            for wi in wire_indices {
                let edge_indices: Vec<usize> =
                    model.wire(wi).expect("wire index in range").edges().to_vec();
                for ei in edge_indices {
                    d = d.max(model.edge(ei).expect("edge index in range").deflection());
                }
            }
            let f = model.face_mut(i).expect("face index in range");
            f.set_deflection(d.max(params.deflection));
            f.set_status(MeshStatus::OUTDATED);
        }

        true
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
        assert!(ModelPreProcessor::perform(&mut model, &params));
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
        assert!(ModelPreProcessor::perform(&mut model, &params));
        for i in 0..model.edges_nb() {
            let d = model.edge(i).expect("edge index").deflection();
            // unit box edges are length 1 → relative deflection == base deflection.
            assert!((d - params.deflection).abs() < 1e-9, "rel deflection {d}");
        }
    }
}
