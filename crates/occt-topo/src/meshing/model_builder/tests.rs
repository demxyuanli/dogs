use super::prelude::*;
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
