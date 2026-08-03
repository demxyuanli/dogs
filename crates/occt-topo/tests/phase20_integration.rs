//! Phase 20 integration — BOPAlgo_Builder end-to-end boolean over boxes.
//!
//! The precise boolean pipeline (PaveFiller intersection → BOPAlgo_Builder
//! rebuild) is exercised through the public `bop_builder2::{fuse, cut, common}`
//! entry points on real box geometry, and the rebuilt solids' volumes are
//! checked against the analytic expectations.

use occt_topo::bop_builder2::{common, cut, fuse};
use occt_topo::shape::TopoShape;

/// Axis-aligned unit box translated by `(dx, dy, dz)`.
fn box_at(dx: f64, dy: f64, dz: f64) -> TopoShape {
    let b = occt_topo::primitives::BRepPrimBox::make_box(1.0, 1.0, 1.0);
    let solid = b.solid.0;
    let v = occt_core::gp::GpVec::new(dx, dy, dz);
    occt_topo::transform::translated(&solid, &v)
}

/// Volume of a solid via `brep_gprop` (mesh integration). Boxes are exact to
/// mesh resolution; tolerance 0.05 on a unit box.
fn volume(s: &TopoShape) -> f64 {
    if s.shape_type() != occt_topo::abs::ShapeType::Solid {
        // Non-solid (face/shell/wire) has no volume.
        return 0.0;
    }
    let mp = occt_topo::brep_gprop::shape_mass_properties(s, 0.05);
    mp.volume
}

/// Sum of volumes of the shapes in a compound (or single shape).
fn total_volume(s: &TopoShape) -> f64 {
    let shapes =
        occt_topo::topo_tools_full::shapes_of(s, occt_topo::abs::ShapeType::Solid);
    if shapes.is_empty() {
        return volume(s);
    }
    shapes.iter().map(volume).sum()
}

#[test]
fn fuse_two_overlapping_boxes() {
    // Boxes [0,1]³ and [0.5,1.5]×[0,1]² overlap on [0.5,1]×[0,1]² = volume 0.5.
    // Fuse union volume = 1 + 1 − 0.5 = 1.5.
    let a = box_at(0.0, 0.0, 0.0);
    let b = box_at(0.5, 0.0, 0.0);
    let r = fuse(&[a], &[b]).expect("fuse");
    let v = total_volume(&r);
    assert!((v - 1.5).abs() < 0.15, "fuse volume {v}, expected 1.5");
}

#[test]
fn fuse_two_disjoint_boxes() {
    // Boxes [0,1]³ and [3,4]×[0,1]² — no overlap, union volume 2.0.
    let a = box_at(0.0, 0.0, 0.0);
    let b = box_at(3.0, 0.0, 0.0);
    let r = fuse(&[a], &[b]).expect("fuse");
    let v = total_volume(&r);
    assert!((v - 2.0).abs() < 0.15, "fuse volume {v}, expected 2.0");
}

#[test]
fn cut_box_removes_overlap() {
    // Cut box [0,1]³ by [0.5,1.5]×[0,1]² removes [0.5,1]×[0,1]² (volume 0.5).
    // Remaining = 1.0 − 0.5 = 0.5.
    let a = box_at(0.0, 0.0, 0.0);
    let b = box_at(0.5, 0.0, 0.0);
    let r = cut(&[a], &[b]).expect("cut");
    let v = total_volume(&r);
    assert!((v - 0.5).abs() < 0.15, "cut volume {v}, expected 0.5");
}

#[test]
fn common_two_overlapping_boxes() {
    // Common of [0,1]³ and [0.5,1.5]×[0,1]² = [0.5,1]×[0,1]² = volume 0.5.
    let a = box_at(0.0, 0.0, 0.0);
    let b = box_at(0.5, 0.0, 0.0);
    let r = common(&[a], &[b]).expect("common");
    let v = total_volume(&r);
    assert!((v - 0.5).abs() < 0.15, "common volume {v}, expected 0.5");
}

#[test]
fn full_pipeline_runs_without_errors() {
    // Smoke test: the whole pipeline on two overlapping boxes completes and
    // returns a non-empty shape.
    let a = box_at(0.0, 0.0, 0.0);
    let b = box_at(0.5, 0.0, 0.0);
    let r = fuse(&[a], &[b]).expect("fuse");
    let n = occt_topo::topo_tools_full::shapes_of(&r, occt_topo::abs::ShapeType::Solid)
        .len();
    assert!(n > 0, "result should contain shapes");
}
